import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import ij.ImagePlus;
import ij.process.FloatProcessor;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

/**
 * Sigma clipped normalized convolution, and the bias correction that makes it unbiased.
 *
 * The estimate is what gets subtracted from every trace, and a trace is 4 pi sigma^2 times the
 * difference - so one ADU of error here is 50 units of trace. Two failure modes are pinned:
 *
 *   - reading *high* where the spots are, which is what happens when the estimator trusts the
 *     mask completely. A real PSF's wings reach past any masking radius a crowded field can
 *     afford, so the light in them lands in pixels the estimator was told to trust, and the
 *     background is then subtracted twice. The clip exists to find those pixels itself;
 *   - reading *low* everywhere, which is what one sided clipping does on its own. Keeping only
 *     the pixels below kappa spread leaves a truncated sample whose mean sits below the true
 *     level - about 1 ADU, since four rounds compound it - and that is put back analytically.
 *
 * The second is the one worth a test, because it is invisible: it happens on a field with nothing
 * in it but noise, where there is no contamination to blame and the answer simply reads low.
 */
class BackgroundEstimateTest {

    private static final int WIDTH = 160;
    private static final int HEIGHT = 160;
    private static final double LEVEL = 20.0;

    private static ImagePlus ones(int width, int height) {
        float[] pixels = new float[width * height];
        java.util.Arrays.fill(pixels, 1.0f);
        return new ImagePlus("mask", new FloatProcessor(width, height, pixels, null));
    }

    /**
     * A spot finder wired up to estimate a background and nothing else. Both masks are open, so
     * every pixel is trusted to start with and anything excluded is excluded by the clip.
     */
    private static smFRETSpotFinder finder(double kappa) {
        smFRETSpotFinder finder = new smFRETSpotFinder();
        finder.backgroundKappa = kappa;
        finder.spotSigma = 2.0;
        finder.overlapMask = ones(WIDTH, HEIGHT);
        finder.backgroundMask = ones(WIDTH, HEIGHT);
        return finder;
    }

    private static double meanOf(ImagePlus image, int margin) {
        double total = 0.0;
        int count = 0;
        for (int y = margin; y < (HEIGHT - margin); y++) {
            for (int x = margin; x < (WIDTH - margin); x++) {
                total += image.getProcessor().getf(x, y);
                count++;
            }
        }
        return total / count;
    }

    private static ImagePlus imageOf(float[] pixels) {
        return new ImagePlus("field", new FloatProcessor(WIDTH, HEIGHT, pixels, null));
    }

    @Test
    @DisplayName("a flat noiseless field estimates to its own level")
    void flatFieldIsRecoveredExactly() {
        float[] pixels = new float[WIDTH * HEIGHT];
        java.util.Arrays.fill(pixels, (float) LEVEL);

        ImagePlus estimate = finder(1.3).backgroundEstimate(imageOf(pixels));

        for (int i = 0; i < pixels.length; i++) {
            assertEquals(LEVEL, estimate.getProcessor().getf(i % WIDTH, i / WIDTH), 1.0e-3,
                    "at index " + i);
        }
    }

    /**
     * Noise and nothing else. There is no contamination here at all, so the only thing that can
     * move the answer off the true level is the truncation the one sided clip introduces - which
     * is exactly the bias truncationConstants exists to put back. Without that correction this
     * reads about an ADU low.
     */
    @Test
    @DisplayName("a noisy flat field is not biased low by the one sided clip")
    void theClipDoesNotBiasTheLevelLow() {
        Random random = new Random(20260810L);
        float[] pixels = new float[WIDTH * HEIGHT];
        for (int i = 0; i < pixels.length; i++) {
            pixels[i] = (float) (LEVEL + Math.sqrt(LEVEL) * random.nextGaussian());
        }

        ImagePlus estimate = finder(1.3).backgroundEstimate(imageOf(pixels));

        // Away from the edges, where a normalized convolution genuinely has less to work with.
        double mean = meanOf(estimate, (int) Math.ceil(3.0 * 14.0));
        assertEquals(LEVEL, mean, 0.25,
                "the estimate read " + mean + " against a true " + LEVEL
                        + " - an uncorrected one sided clip reads about 1 ADU low");
    }

    /**
     * The same, at a range of clip thresholds. The correction is a function of kappa, so if it
     * were wrong the error would grow as the clip tightens - a tight clip removes more of the
     * sample and so needs more put back.
     */
    @Test
    @DisplayName("the level is unbiased at every clip threshold")
    void theCorrectionHoldsAcrossKappa() {
        Random random = new Random(99L);
        float[] pixels = new float[WIDTH * HEIGHT];
        for (int i = 0; i < pixels.length; i++) {
            pixels[i] = (float) (LEVEL + Math.sqrt(LEVEL) * random.nextGaussian());
        }

        for (double kappa : new double[] {0.8, 1.0, 1.3, 1.8, 3.0}) {
            double mean = meanOf(finder(kappa).backgroundEstimate(imageOf(pixels)),
                    (int) Math.ceil(3.0 * 14.0));
            assertEquals(LEVEL, mean, 0.3, "at kappa " + kappa + " the estimate read " + mean);
        }
    }

    /**
     * The case the clip is really for: spots present and *not* masked out. The estimator has to
     * find them itself and still land on the background, rather than averaging the spot light in.
     */
    @Test
    @DisplayName("unmasked spots are clipped out rather than averaged in")
    void spotsAreFoundByTheClip() {
        Random random = new Random(7L);

        // Roughly the density of the real data - hel1 carries about 430 spots in a 256 by 512
        // half, one per 300 pixels, and a 16 pixel grid here is one per 356. Density is the whole
        // question for a clip: at one spot per 121 pixels the field is *half* spot light by area
        // and no clip can recover the background, because most pixels are contaminated.
        List<SyntheticField.Spot> spots = new ArrayList<>();
        for (int y = 12; y < (HEIGHT - 12); y += 16) {
            for (int x = 12; x < (WIDTH - 12); x += 16) {
                spots.add(new SyntheticField.Spot(x, y, 1000.0));
            }
        }

        float[] pixels = SyntheticField.half(WIDTH, HEIGHT, spots, 2.0, LEVEL);
        for (int i = 0; i < pixels.length; i++) {
            pixels[i] += (float) (Math.sqrt(LEVEL) * random.nextGaussian());
        }

        // What a plain average would give: the background plus all the spot light spread over the
        // field. Measuring it rather than deriving it keeps the comparison honest.
        double raw = 0.0;
        for (float value : pixels) {
            raw += value;
        }
        raw /= pixels.length;

        double mean = meanOf(finder(1.3).backgroundEstimate(imageOf(pixels)),
                (int) Math.ceil(3.0 * 14.0));

        // Told nothing about where the spots are, so this is the clip working alone. It does not
        // remove all of the contamination - a PSF's wings are spread too thin to clip - but it
        // has to remove most of it, which is the claim: telling the estimator where the spots are
        // is not what makes it work.
        assertTrue(Math.abs(mean - LEVEL) < 0.5 * Math.abs(raw - LEVEL),
                "the estimate read " + mean + " where a plain average of the same field reads "
                        + raw + " and the true background is " + LEVEL);
    }

    /**
     * The realistic arrangement: the clip *and* the mask, which is what the pipeline runs. The
     * mask takes out the cores at 2.5 sigma and the clip cleans up the wings that reach past it,
     * and between them the level comes back.
     */
    @Test
    @DisplayName("with the spot mask as well, the level comes back")
    void maskedSpotsLeaveTheLevelAlone() {
        Random random = new Random(7L);
        List<SyntheticField.Spot> spots = new ArrayList<>();
        double[][] flagged = new double[81][3];
        int index = 0;
        for (int y = 12; y < (HEIGHT - 12); y += 16) {
            for (int x = 12; x < (WIDTH - 12); x += 16) {
                spots.add(new SyntheticField.Spot(x, y, 1000.0));
                flagged[index][0] = 1.0;
                flagged[index][1] = x;
                flagged[index][2] = y;
                index++;
            }
        }
        double[][] present = java.util.Arrays.copyOf(flagged, index);

        float[] pixels = SyntheticField.half(WIDTH, HEIGHT, spots, 2.0, LEVEL);
        for (int i = 0; i < pixels.length; i++) {
            pixels[i] += (float) (Math.sqrt(LEVEL) * random.nextGaussian());
        }

        smFRETSpotFinder finder = finder(1.3);
        finder.backgroundMask = finder.neighborhoodMaskToBackgroundMask(
                finder.createSpotsNeighborhoodMask(present, WIDTH, HEIGHT,
                        smFRETSpotFinder.spotMarginFor(2.0)));

        double mean = meanOf(finder.backgroundEstimate(imageOf(pixels)),
                (int) Math.ceil(3.0 * 14.0));
        assertEquals(LEVEL, mean, 0.3,
                "the estimate read " + mean + " on a field whose background is " + LEVEL);
    }

    /**
     * A slowly varying illumination profile is the thing the sigma 14 smoothing exists to follow,
     * so the estimate has to track it rather than flattening it to a single number.
     *
     * With shot noise on it, deliberately. On a *noiseless* structured field the clip degenerates:
     * the residual it measures is then the estimator's own edge behaviour rather than noise, the
     * robust spread comes out near zero, and clipping at 1.3 of near zero removes a systematic
     * half of the field. That cannot happen on real data - a background always carries noise, and
     * the flat noiseless case is covered above where the residual is identically zero and the
     * loop exits at once.
     */
    @Test
    @DisplayName("a smooth illumination gradient is followed")
    void anIlluminationGradientIsTracked() {
        Random random = new Random(31337L);
        float[] pixels = new float[WIDTH * HEIGHT];
        for (int y = 0; y < HEIGHT; y++) {
            for (int x = 0; x < WIDTH; x++) {
                double level = LEVEL + 10.0 * x / WIDTH;
                pixels[y * WIDTH + x] = (float) (level + Math.sqrt(level) * random.nextGaussian());
            }
        }

        ImagePlus estimate = finder(1.3).backgroundEstimate(imageOf(pixels));

        int margin = (int) Math.ceil(3.0 * 14.0);
        for (int y = margin; y < (HEIGHT - margin); y += 7) {
            for (int x = margin; x < (WIDTH - margin); x += 7) {
                assertEquals(LEVEL + 10.0 * x / WIDTH, estimate.getProcessor().getf(x, y), 0.4,
                        "at " + x + "," + y);
            }
        }
    }

    /**
     * backgroundKappa above zero overrides the default, which is the escape hatch for a real
     * illumination profile with structure a hard clip would eat - the example long movie wanted
     * 1.8. Zero means "use the derived default", which is 1.3 now rather than a function of the
     * spot size.
     */
    @Test
    @DisplayName("kappa is a constant default that an explicit setting overrides")
    void kappaOverride() {
        smFRETSpotFinder finder = new smFRETSpotFinder();
        finder.spotSigma = 2.0;

        finder.backgroundKappa = 0.0;
        assertEquals(1.3, finder.clippingThreshold(), 1.0e-12, "the derived default");

        finder.spotSigma = 3.0;
        assertEquals(1.3, finder.clippingThreshold(), 1.0e-12,
                "the default no longer depends on the spot size");

        finder.backgroundKappa = 1.8;
        assertEquals(1.8, finder.clippingThreshold(), 1.0e-12, "an explicit override");
    }
}
