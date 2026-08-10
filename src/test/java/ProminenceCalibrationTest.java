import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import ij.ImagePlus;
import ij.process.FloatProcessor;
import java.util.Arrays;
import java.util.Collections;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

/**
 * What spotProminence means: a fraction of what a clean lone spot scores.
 *
 * That reading is the whole of the current design. The parameter used to be an unnormalized ratio
 * whose ideal ran from 1.649 at spotSigma 1 to 4.010 at spotSigma 3, so the old default of 1.8
 * rejected *every* real spot at spotSigma 1. The fix was a thin annulus at two sigma plus an
 * explicit normalization by what a noiseless isolated spot of that size gives, and these pin the
 * three properties that buys - a lone spot scores 1.0, it scores 1.0 at any size, and it scores
 * 1.0 at any brightness.
 *
 * The brightness one is not decoration. The old ring clamped each pixel at 1 ADU, an absolute
 * floor, which made a dim spot's noise ring score enormously and a bright spot's real PSF wings
 * score honestly - so the filter preferentially threw away the *good* spots, 52% of those at SNR
 * 15-25 against 16% below SNR 10 on the example data.
 */
class ProminenceCalibrationTest {

    private static final int SIZE = 96;
    private static final int CENTRE = SIZE / 2;

    /**
     * A spot finder set up to measure prominence and nothing else.
     *
     * cameraBlackLevel and the background image are both zero, which drives the ring's noise bias
     * term to exactly zero - the term is 1.2816 * sqrt(level / (gain * frames)) - so what these
     * measure is the geometry rather than the geometry plus a noise estimate.
     */
    private static smFRETSpotFinder finderAt(double sigma) {
        smFRETSpotFinder finder = new smFRETSpotFinder();
        finder.spotSigma = sigma;
        finder.cameraGain = 1.0;
        finder.cameraBlackLevel = 0;
        finder.cachedFrames = 20;
        finder.spotProminence = 0.4;
        return finder;
    }

    /**
     * A *point sampled* Gaussian, which is what idealProminence's reference profile is.
     *
     * Deliberately not SyntheticField's pixel integrated spot. This is the exact case the
     * normalization is defined against, so it must come out at exactly 1.0; the realistic binned
     * case is a separate test below with its own, looser, expectation.
     */
    private static ImagePlus pointSampledSpot(double sigma, double amplitude, double... extra) {
        float[] pixels = new float[SIZE * SIZE];
        for (int y = 0; y < SIZE; y++) {
            for (int x = 0; x < SIZE; x++) {
                double dx = x - CENTRE;
                double dy = y - CENTRE;
                pixels[y * SIZE + x] =
                        (float) (amplitude * Math.exp(-(dx * dx + dy * dy) / (2.0 * sigma * sigma)));
            }
        }

        // Optional second spot, as {offsetX, offsetY, amplitude}, for the doublet case.
        if (extra.length == 3) {
            for (int y = 0; y < SIZE; y++) {
                for (int x = 0; x < SIZE; x++) {
                    double dx = x - CENTRE - extra[0];
                    double dy = y - CENTRE - extra[1];
                    pixels[y * SIZE + x] +=
                            (float) (extra[2] * Math.exp(-(dx * dx + dy * dy) / (2.0 * sigma * sigma)));
                }
            }
        }
        return new ImagePlus("spot", new FloatProcessor(SIZE, SIZE, pixels, null));
    }

    private static ImagePlus zeroBackground() {
        return new ImagePlus("background", new FloatProcessor(SIZE, SIZE, new float[SIZE * SIZE], null));
    }

    /** The flag-prefixed in-memory layout: good, x, y. */
    private static double[][] oneSpotAtCentre() {
        return new double[][] {{1.0, CENTRE, CENTRE}};
    }

    private static double prominenceOf(smFRETSpotFinder finder, ImagePlus image) {
        double[][] filtered = finder.spotFilterProminence(oneSpotAtCentre(), image, zeroBackground());
        return filtered[0][filtered[0].length - 1];
    }

    @ParameterizedTest(name = "spotSigma={0}")
    @ValueSource(doubles = {1.0, 1.5, 2.0, 2.5, 3.0})
    @DisplayName("a noiseless lone spot scores 1.0 at every spot size")
    void loneSpotScoresOne(double sigma) {
        smFRETSpotFinder finder = finderAt(sigma);
        assertEquals(1.0, prominenceOf(finder, pointSampledSpot(sigma, 100.0)), 1.0e-5,
                "an isolated spot at spotSigma " + sigma);
    }

    /**
     * The property the old absolute ADU floor destroyed. A ten thousandfold range of brightness
     * has to give the same number, because the measurement is a ratio and the reference it is
     * divided by does not depend on the spot at all.
     */
    @ParameterizedTest(name = "amplitude={0}")
    @ValueSource(doubles = {0.5, 5.0, 100.0, 5000.0})
    @DisplayName("prominence does not depend on how bright the spot is")
    void prominenceIsBrightnessIndependent(double amplitude) {
        smFRETSpotFinder finder = finderAt(2.0);
        assertEquals(1.0, prominenceOf(finder, pointSampledSpot(2.0, amplitude)), 1.0e-4,
                "an isolated spot of amplitude " + amplitude);
    }

    /**
     * A real camera bins the PSF over its pixels, so the measured ratio is not exactly the
     * continuous one the reference is computed from. It should still be close, which is what
     * makes 1.0 a usable reading on real data rather than only on paper - on hel1 and hel7 the
     * median lands at 0.996 and 0.975.
     */
    @ParameterizedTest(name = "spotSigma={0}")
    @ValueSource(doubles = {1.0, 1.5, 2.0, 2.5, 3.0})
    @DisplayName("a pixel integrated lone spot still scores about 1.0")
    void binnedLoneSpotScoresAboutOne(double sigma) {
        smFRETSpotFinder finder = finderAt(sigma);

        float[] pixels = SyntheticField.half(SIZE, SIZE,
                Collections.singletonList(new SyntheticField.Spot(CENTRE, CENTRE, 4000.0)),
                sigma, 0.0);
        ImagePlus image = new ImagePlus("spot", new FloatProcessor(SIZE, SIZE, pixels, null));

        double prominence = prominenceOf(finder, image);
        assertTrue((prominence > 0.9) && (prominence < 1.1),
                "pixel integrated lone spot at spotSigma " + sigma + " scored " + prominence);
    }

    /**
     * A pair straddling the measured pixel, which is where MaximumFinder actually puts it.
     *
     * Centring matters and getting it wrong makes the filter look better than it is: measuring at
     * one member of a pair puts the whole of the other member on one side of the ring, so even a
     * tight pair scores low. What the pipeline really sees is a single maximum at the combined
     * peak, between the two, and a ring that is lifted symmetrically - a much weaker signal.
     */
    private static ImagePlus centredPair(double sigma, double separation) {
        float[] pixels = new float[SIZE * SIZE];
        for (double cx : new double[] {CENTRE - separation / 2.0, CENTRE + separation / 2.0}) {
            for (int y = 0; y < SIZE; y++) {
                for (int x = 0; x < SIZE; x++) {
                    double dx = x - cx;
                    double dy = y - CENTRE;
                    pixels[y * SIZE + x] +=
                            (float) (100.0 * Math.exp(-(dx * dx + dy * dy) / (2.0 * sigma * sigma)));
                }
            }
        }
        return new ImagePlus("pair", new FloatProcessor(SIZE, SIZE, pixels, null));
    }

    /**
     * What the filter is for: the score has to fall as a pair opens up, over the whole range, so
     * that a single threshold means the same ordering everywhere.
     */
    @ParameterizedTest(name = "spotSigma={0}")
    @ValueSource(doubles = {1.0, 2.0, 3.0})
    @DisplayName("prominence falls monotonically as a pair separates")
    void prominenceFallsWithSeparation(double sigma) {
        smFRETSpotFinder finder = finderAt(sigma);

        double previous = Double.MAX_VALUE;
        for (double separation : new double[] {0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0}) {
            double prominence = prominenceOf(finder, centredPair(sigma, separation * sigma));
            assertTrue(prominence < previous,
                    "separation " + separation + " sigma scored " + prominence
                            + ", not below the previous " + previous);
            previous = prominence;
        }
    }

    /**
     * Two coincident molecules are one molecule of twice the brightness, and prominence is
     * brightness independent, so this is exactly 1.0 rather than approximately. It is the
     * cleanest statement of what ADR-0001 is about: at zero separation there is nothing in the
     * image to see, so no statistic computed from the image can see it.
     */
    @Test
    @DisplayName("a coincident pair is indistinguishable from one brighter spot")
    void aCoincidentPairScoresExactlyOne() {
        assertEquals(1.0, prominenceOf(finderAt(2.0), centredPair(2.0, 0.0)), 1.0e-5);
    }

    /**
     * The blind spot, pinned so it is not rediscovered as a bug. A pair at a quarter of a sigma
     * still scores near the isolated value and sails through the default 0.4 - which is the
     * premise ADR-0001 rests on, and the reason the prominence filter is not a doublet detector.
     * The band the filter really works in is the 1.5 to 2.5 sigma one below.
     */
    @ParameterizedTest(name = "spotSigma={0}")
    @ValueSource(doubles = {1.0, 2.0, 3.0})
    @DisplayName("a pair at a quarter sigma passes, a pair at two sigma does not survive as well")
    void theFilterSeesOnlyWellSeparatedPairs(double sigma) {
        smFRETSpotFinder finder = finderAt(sigma);

        double tight = prominenceOf(finder, centredPair(sigma, 0.25 * sigma));
        assertTrue(tight > 0.97,
                "a pair at 0.25 sigma scored " + tight + " - if the ring can now see that,"
                        + " ADR-0001 needs revisiting with new numbers");

        double open = prominenceOf(finder, centredPair(sigma, 2.0 * sigma));
        assertTrue(open < 0.7,
                "a pair at 2 sigma scored " + open + ", which the filter should be catching");
    }

    /**
     * A spot whose ring reads at or below zero has no measurable neighbour, so it passes rather
     * than dividing by nothing. Capped so the CSV cannot end up holding an infinity.
     */
    @Test
    @DisplayName("an empty ring is capped rather than infinite")
    void anEmptyRingIsCapped() {
        smFRETSpotFinder finder = finderAt(2.0);

        // A single bright pixel: the ring at two sigma sees only zeros.
        float[] pixels = new float[SIZE * SIZE];
        pixels[CENTRE * SIZE + CENTRE] = 500.0f;
        ImagePlus spike = new ImagePlus("spike", new FloatProcessor(SIZE, SIZE, pixels, null));

        double prominence = prominenceOf(finder, spike);
        assertTrue(Double.isFinite(prominence), "prominence was " + prominence);
        assertEquals(10.0, prominence, 0.0, "the documented cap");
    }

    /**
     * The ring geometry itself: a thin annulus at two sigma, never empty, and never reaching back
     * in toward the spot the way the old (2s-1, 2s+1) ring did at large sigma.
     */
    @ParameterizedTest(name = "spotSigma={0}")
    @ValueSource(doubles = {1.0, 1.5, 2.0, 2.5, 3.0})
    @DisplayName("the ring is a thin annulus at two sigma")
    void ringGeometry(double sigma) {
        int[][] ring = finderAt(sigma).prominenceRing();

        assertTrue(ring.length >= 16,
                "the ring holds only " + ring.length + " pixels at spotSigma " + sigma);
        for (int[] offset : ring) {
            double radius = Math.hypot(offset[0], offset[1]);
            assertEquals(2.0 * sigma, radius, 0.6 + 1.0e-9,
                    "ring pixel " + Arrays.toString(offset) + " at spotSigma " + sigma);
        }
    }

    /**
     * The filter marks a spot bad by clearing column 0, and records the score whether or not it
     * passed - the score is a diagnostic column in the CSV, so it has to survive rejection.
     */
    @Test
    @DisplayName("the flag is cleared but the score is still recorded")
    void rejectionKeepsTheScore() {
        smFRETSpotFinder finder = finderAt(2.0);
        finder.spotProminence = 2.0;      // Above what any real spot can reach.

        double[][] filtered = finder.spotFilterProminence(oneSpotAtCentre(),
                pointSampledSpot(2.0, 100.0), zeroBackground());

        assertEquals(0.0, filtered[0][0], 0.0, "the good flag");
        assertEquals(1.0, filtered[0][filtered[0].length - 1], 1.0e-5, "the recorded prominence");
    }
}
