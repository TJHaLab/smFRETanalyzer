import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

/**
 * The measurement, scored against a field whose PSF is known exactly.
 *
 * Two things here are the reason the method looks the way it does, and both are checked by
 * breaking them rather than only by asserting the good case works:
 *
 *   - **contaminated pixels are dropped, not contaminated spots.** A FRET field is crowded, so
 *     rejecting every spot with a neighbour throws most of the field away and the survivors still
 *     have neighbour light in their wings;
 *   - **the patch border is not background.** An aberrated Airy still holds around a percent of
 *     its peak at the edge of a ten pixel patch, so subtracting the border median takes real
 *     light off every point of the profile - which lands almost entirely on the wings, the one
 *     part this model exists to measure.
 *
 * The fields are rendered by *point sampling* the model, not by integrating it over each pixel.
 * That is the opposite of what SyntheticField does and it is deliberate: what is under test here
 * is the extraction, masking, pooling and fitting, so the field is made to be exactly what the
 * model says and any recovered error belongs to the pipeline rather than to binning. Spots sit on
 * integer coordinates for the same reason - it removes the centring error, which is measured
 * separately below.
 */
class PSFMeasurementTest {

    private static final double TRUE_SIGMA = 1.4;
    private static final double TRUE_WAVES = 0.4;
    private static final double BACKGROUND = 20.0;
    private static final double PEAK = 100.0;

    private static final int FIELD = 600;
    private static final int PATCH = 10;
    private static final double NEIGHBOUR_MASK = 6.0;
    private static final int BINS = 20;

    // Spots on the crowded field, giving about the density of the real data. Note that a
    // neighbour further than patch + mask - 16 pixels at these settings - cannot reach a patch
    // pixel at all, so a fixture whose spots are all further apart than that makes masking a
    // no-op and the two paths return bit-identical answers.
    private static final int CROWDED_SPOTS = 1000;

    // How far each spot is rendered. Far past the patch and its border, so nothing the
    // measurement looks at sees the truncation.
    private static final double RENDER_REACH = 30.0;

    /** The model on a fine radial grid, so a field can be rendered by interpolation. */
    private static final double GRID_STEP = 0.01;
    private static final double[] GRID = buildGrid();

    private static double[] buildGrid() {
        int n = (int) Math.ceil(RENDER_REACH / GRID_STEP) + 2;
        double[] radii = new double[n];
        for (int i = 0; i < n; i++) {
            radii[i] = i * GRID_STEP;
        }
        return smFRETPSF.airyProfile(radii, TRUE_SIGMA, TRUE_WAVES);
    }

    private static double shapeAt(double radius) {
        if (radius >= RENDER_REACH) {
            return 0.0;
        }
        double position = radius / GRID_STEP;
        int index = (int) position;
        double fraction = position - index;
        return GRID[index] + fraction * (GRID[index + 1] - GRID[index]);
    }

    /** A field of identical spots on a grid, each peaking PEAK above the background. */
    private static float[] field(List<double[]> spots) {
        float[] pixels = new float[FIELD * FIELD];
        java.util.Arrays.fill(pixels, (float) BACKGROUND);

        int reach = (int) Math.ceil(RENDER_REACH);
        for (double[] spot : spots) {
            int cx = (int) Math.round(spot[0]);
            int cy = (int) Math.round(spot[1]);
            for (int y = Math.max(0, cy - reach); y <= Math.min(FIELD - 1, cy + reach); y++) {
                for (int x = Math.max(0, cx - reach); x <= Math.min(FIELD - 1, cx + reach); x++) {
                    double r = Math.hypot(x - spot[0], y - spot[1]);
                    pixels[y * FIELD + x] += (float) (PEAK * shapeAt(r));
                }
            }
        }
        return pixels;
    }

    /**
     * A randomly placed field at roughly the density of the real data.
     *
     * Random rather than a grid, and that is not a detail. On a regular grid every neighbour sits
     * at the identical distance, so either none of them reaches into the patch or all of them
     * strip the same slice off every spot's border at once - at spacing 13 that rejects 1600 of
     * 1764 spots, a failure mode no real field produces. A random field has the distribution the
     * method was designed against: most spots comfortable, some with a close neighbour, a few
     * unusable.
     *
     * hel1 carries about 430 spots in a 256 by 512 half, one per 300 pixels; 1000 spots in
     * 600 by 600 is one per 360.
     */
    private static List<double[]> randomField(int count, long seed) {
        java.util.Random random = new java.util.Random(seed);
        List<double[]> spots = new ArrayList<>();

        // A minimum separation, because MaximumFinder cannot report two maxima in adjacent
        // pixels either - without it the fixture is full of pairs that are really one spot.
        double minimumSeparation = 4.0;
        int attempts = 0;
        while ((spots.size() < count) && (attempts < (200 * count))) {
            attempts++;
            double x = 30 + random.nextInt(FIELD - 60);
            double y = 30 + random.nextInt(FIELD - 60);

            boolean clear = true;
            for (double[] other : spots) {
                if (Math.hypot(other[0] - x, other[1] - y) < minimumSeparation) {
                    clear = false;
                    break;
                }
            }
            if (clear) {
                spots.add(new double[] {x, y});
            }
        }
        return spots;
    }

    /** Spot positions on a square grid, in the reloaded layout: x at column 0. */
    private static List<double[]> grid(int spacing, int margin) {
        List<double[]> spots = new ArrayList<>();
        for (int y = margin; y <= (FIELD - 1 - margin); y += spacing) {
            for (int x = margin; x <= (FIELD - 1 - margin); x += spacing) {
                spots.add(new double[] {x, y});
            }
        }
        return spots;
    }

    private static double[][] asArray(List<double[]> spots) {
        return spots.toArray(new double[0][]);
    }

    private static smFRETPSF.Measurement measure(List<double[]> spots, double neighbourMask) {
        return smFRETPSF.analyse(
                smFRETPSF.extract(field(spots), FIELD, FIELD, asArray(spots), PATCH, neighbourMask),
                BINS);
    }

    /**
     * The baseline: a field sparse enough that neighbours barely matter, so what is being scored
     * is the pooling, the border correction and the fit rather than the masking.
     */
    @Test
    @DisplayName("a sparse field recovers the PSF it was built from")
    void sparseFieldRecoversTheTruth() {
        List<double[]> spots = grid(50, 40);
        smFRETPSF.Measurement measurement = measure(spots, NEIGHBOUR_MASK);

        assertEquals(spots.size(), measurement.samples.spotsUsed,
                "every spot on a sparse field should be usable");

        smFRETPSF.Fit fit = measurement.withPedestal;
        assertEquals(TRUE_SIGMA, fit.sigma, 0.08, "sigma");
        assertEquals(TRUE_WAVES, fit.waves, 0.08, "waves");

        // The field has no halo in it, so there is nothing for a pedestal to absorb - which is
        // what makes the 3.3% found on the real data worth believing.
        assertTrue(fit.pedestal < 0.01,
                "a pedestal of " + fit.pedestal + " on a field that has none");
    }

    /**
     * The case the method is for. At this spacing every spot has neighbours inside its patch, and
     * masking their pixels has to recover the same answer the sparse field gave.
     */
    @Test
    @DisplayName("a crowded field recovers it too, once neighbours are masked")
    void crowdedFieldRecoversTheTruth() {

        List<double[]> spots = randomField(CROWDED_SPOTS, 20260810L);
        smFRETPSF.Measurement measurement = measure(spots, NEIGHBOUR_MASK);

        // The claim masking exists to support, scored against the thing it replaced rather than
        // against a threshold: how many spots would be left if every spot with a neighbour
        // reaching into its patch were thrown away instead.
        int isolated = 0;
        for (double[] spot : spots) {
            boolean alone = true;
            for (double[] other : spots) {
                double gap = Math.hypot(other[0] - spot[0], other[1] - spot[1]);
                if ((gap > 1.0e-6) && (gap < Math.hypot(PATCH, PATCH))) {
                    alone = false;
                    break;
                }
            }
            if (alone) {
                isolated++;
            }
        }

        assertTrue(measurement.samples.spotsUsed > (3 * isolated),
                "masking kept " + measurement.samples.spotsUsed + " of " + spots.size()
                        + " where rejecting whole spots would keep " + isolated
                        + " - not enough of a difference to be worth the method");

        smFRETPSF.Fit fit = measurement.withPedestal;
        assertEquals(TRUE_SIGMA, fit.sigma, 0.12, "sigma");
        assertEquals(TRUE_WAVES, fit.waves, 0.12, "waves");
    }

    /**
     * The same crowded field with the masking switched off, which is what the measurement would
     * be if neighbour light were simply tolerated. The neighbours' wings pile into every patch,
     * so the profile reads too bright out at radius and the fit answers with a PSF that is wrong
     * in a specific direction.
     */
    @Test
    @DisplayName("without masking, the same field gives a worse answer")
    void maskingIsDoingWork() {
        List<double[]> spots = randomField(CROWDED_SPOTS, 20260810L);

        smFRETPSF.Measurement masked = measure(spots, NEIGHBOUR_MASK);
        smFRETPSF.Measurement unmasked = measure(spots, 0.0);

        double maskedError = Math.abs(masked.withPedestal.waves - TRUE_WAVES)
                + Math.abs(masked.withPedestal.sigma - TRUE_SIGMA);
        double unmaskedError = Math.abs(unmasked.withPedestal.waves - TRUE_WAVES)
                + Math.abs(unmasked.withPedestal.sigma - TRUE_SIGMA);

        assertTrue(maskedError < unmaskedError,
                "masked error " + maskedError + " against unmasked " + unmaskedError
                        + " - if masking is not helping, the fixture is not crowded enough to"
                        + " be testing anything");
    }

    /**
     * The border correction, checked by looking at what the profile would have been without it.
     * The uncorrected profile has had a real slice of the PSF subtracted off every point, and
     * because that slice is a constant it is negligible against the core and large against the
     * wings - so a fit to it comes back with too little aberration.
     */
    @Test
    @DisplayName("the border correction is what keeps the wings")
    void borderCorrectionMatters() {
        List<double[]> spots = grid(50, 40);
        smFRETPSF.Measurement measurement = measure(spots, NEIGHBOUR_MASK);

        assertTrue(measurement.borderSettled, "the border correction did not settle");
        assertTrue(measurement.borderLevel > 0.001,
                "border level " + measurement.borderLevel
                        + " - an aberrated Airy should have left more than this on the border");

        // Refit the uncorrected profile over the same bins the corrected one used.
        boolean[] usable = measurement.usable();
        List<Double> radii = new ArrayList<>();
        List<Double> raw = new ArrayList<>();
        for (int i = 0; i < usable.length; i++) {
            if (usable[i]) {
                radii.add(measurement.binCentre[i]);
                raw.add(measurement.binRaw[i]);
            }
        }
        double[] r = new double[radii.size()];
        double[] v = new double[raw.size()];
        for (int i = 0; i < r.length; i++) {
            r[i] = radii.get(i);
            v[i] = raw.get(i);
        }
        smFRETPSF.Fit uncorrected = smFRETPSF.fit(r, v, true);

        // The direction first: subtracting the border took wing light off, so an uncorrected fit
        // has to report *less* aberration than a corrected one.
        assertTrue(uncorrected.waves < measurement.withPedestal.waves,
                "uncorrected fit gave " + uncorrected.waves + " waves against the corrected "
                        + measurement.withPedestal.waves + " - the correction should be adding"
                        + " aberration back, not taking it away");

        // And the size of it. Stated as a ratio rather than an absolute margin because how much
        // the border holds depends on the PSF and the patch: on this fixture the corrected fit
        // lands about twenty times closer to the truth, and a third of that is still decisive.
        double correctedError = Math.abs(measurement.withPedestal.waves - TRUE_WAVES);
        double uncorrectedError = Math.abs(uncorrected.waves - TRUE_WAVES);
        assertTrue(uncorrectedError > (3.0 * correctedError),
                "uncorrected was " + uncorrectedError + " from the truth against the corrected "
                        + correctedError + ", which is not a large enough difference to say the"
                        + " correction is doing the work");
    }

    /**
     * A spot with a neighbour close enough to mask its own centre is two spots, not one, and the
     * peak the whole patch is normalised by would not be either of them.
     */
    @Test
    @DisplayName("a spot whose own centre is masked is skipped")
    void unresolvableSpotsAreSkipped() {
        List<double[]> spots = new ArrayList<>();
        spots.add(new double[] {100, 100});
        spots.add(new double[] {103, 100});        // Inside the 6 px mask of the first.
        spots.add(new double[] {300, 300});        // Alone.

        smFRETPSF.Samples samples = smFRETPSF.extract(field(spots), FIELD, FIELD,
                asArray(spots), PATCH, NEIGHBOUR_MASK);

        assertEquals(1, samples.spotsUsed, "only the isolated spot should have been usable");
        assertEquals(3, samples.spotsTotal);
    }

    @Test
    @DisplayName("spots too near the frame edge are skipped")
    void edgeSpotsAreSkipped() {
        List<double[]> spots = new ArrayList<>();
        spots.add(new double[] {5, 300});                    // Patch would run off the left.
        spots.add(new double[] {300, FIELD - 6});            // And off the bottom.
        spots.add(new double[] {300, 300});

        smFRETPSF.Samples samples = smFRETPSF.extract(field(spots), FIELD, FIELD,
                asArray(spots), PATCH, NEIGHBOUR_MASK);

        assertEquals(1, samples.spotsUsed);
    }

    /**
     * The 2D image is pooled the same way the radial profile is - a ratio of sums per pixel - so
     * the two are the same measurement seen two ways and the image's own radial average has to
     * agree with the profile.
     */
    @Test
    @DisplayName("the 2D image agrees with the radial profile")
    void imageAndProfileAgree() {
        List<double[]> spots = grid(50, 40);
        smFRETPSF.Measurement measurement = measure(spots, NEIGHBOUR_MASK);
        double[] image = measurement.samples.image();
        int size = measurement.samples.size;

        assertEquals(1.0, image[PATCH * size + PATCH], 1.0e-9,
                "the centre pixel is the peak everything was normalised by");

        // Average the image over the annulus each bin covers, and compare against the bin.
        double width = (double) PATCH / BINS;
        for (int bin = 0; bin < BINS; bin++) {
            if (measurement.binCount[bin] <= 50) {
                continue;
            }
            double total = 0.0;
            int count = 0;
            for (int dy = 0; dy < size; dy++) {
                for (int dx = 0; dx < size; dx++) {
                    double r = Math.hypot(dx - PATCH, dy - PATCH);
                    if (((int) (r / width)) != bin) {
                        continue;
                    }
                    double value = image[dy * size + dx];
                    if (!Double.isNaN(value)) {
                        total += value;
                        count++;
                    }
                }
            }
            if (count == 0) {
                continue;
            }

            // Not identical: the bin pools every contributing pixel of every spot at once, while
            // this averages per-pixel ratios that were each pooled over spots. On a uniform field
            // they agree closely, which is all this is checking.
            assertEquals(measurement.binRaw[bin], total / count, 0.01,
                    "bin " + bin + " at radius " + measurement.binCentre[bin]);
        }
    }

    /**
     * The image the panel draws is border corrected, the same as the profile beside it.
     *
     * Samples.image() is the raw pooled ratio with the border median still out of it, so its
     * outer pixels sit near zero by construction and some go negative - on the example data the
     * donor bottoms out at -0.0024, which clamps to black and reads as "no light here" when the
     * profile below it is saying 0.024.
     */
    @Test
    @DisplayName("the drawn image has the border put back, like the profile")
    void correctedImageMatchesTheProfileFooting() {
        List<double[]> spots = grid(50, 40);
        smFRETPSF.Measurement measurement = measure(spots, NEIGHBOUR_MASK);

        double[] raw = measurement.samples.image();
        double[] corrected = measurement.correctedImage();
        assertTrue(measurement.borderLevel > 0.0, "nothing to correct on this fixture");

        for (int i = 0; i < raw.length; i++) {
            if (Double.isNaN(raw[i])) {
                assertTrue(Double.isNaN(corrected[i]), "an unmeasured pixel should stay unmeasured");
                continue;
            }
            assertEquals((raw[i] * (1.0 - measurement.borderLevel)) + measurement.borderLevel,
                    corrected[i], 1.0e-12, "at index " + i);

            // Never downward. Not strictly upward either: the peak is 1 by construction and
            // v(1-L)+L has a fixed point there, which is the property that keeps the corrected
            // image still normalised to its own peak.
            assertTrue(corrected[i] >= (raw[i] - 1.0e-12),
                    "the correction should not lower anything, at index " + i);
        }

        int size = measurement.samples.size;
        assertTrue(corrected[0] > raw[0], "the corner should have been lifted");
        assertEquals(1.0, corrected[(size / 2) * size + (size / 2)], 1.0e-12,
                "the peak is a fixed point of the correction");

        // And it is the same lift the profile got, so the two agree at the centre.
        assertEquals(measurement.binProfile[0], corrected[PATCH * size + PATCH], 1.0e-9,
                "the centre pixel and the first bin are both the peak");
    }

    /**
     * Pixels no spot could contribute to are NaN rather than zero, so a display can tell "nothing
     * measured here" from "no light here" - on a crowded field at a large mask those are very
     * different statements.
     */
    @Test
    @DisplayName("unmeasured pixels are marked, not zeroed")
    void unmeasuredPixelsAreNaN() {

        // Exactly one spot may contribute, or the question is meaningless: a second usable spot
        // covers the hole in the first one's patch from its own vantage point, since the image is
        // an average over spots rather than a picture of one. So A is isolated enough to be
        // usable, and B and C sit two pixels apart, which masks each other's centres and takes
        // both of them out - while still masking part of A's patch.
        List<double[]> spots = new ArrayList<>();
        spots.add(new double[] {300, 300});
        spots.add(new double[] {312, 300});
        spots.add(new double[] {314, 300});

        smFRETPSF.Samples samples = smFRETPSF.extract(field(spots), FIELD, FIELD,
                asArray(spots), PATCH, NEIGHBOUR_MASK);
        assertEquals(1, samples.spotsUsed, "only the first spot should be usable");

        double[] image = samples.image();
        int size = samples.size;

        // Absolute (308, 300), four pixels from B, so inside its mask.
        assertTrue(Double.isNaN(image[PATCH * size + (PATCH + 8)]),
                "the pixel under the neighbour should be unmeasured");
        assertTrue(!Double.isNaN(image[PATCH * size + PATCH]), "the centre should be measured");
        assertEquals(0, samples.imageCount[PATCH * size + (PATCH + 8)], "contributors");
    }
}
