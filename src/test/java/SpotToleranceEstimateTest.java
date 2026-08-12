import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.scijava.log.StderrLogService;

/**
 * Estimating SpotTolerance from the image, which is issue #10.
 *
 * The parameter is MaximumFinder's flood fill depth, in raw image units. That is what made it
 * confusing to set: SpotThreshold sits next to it and reads in sigma, so the two looked like a
 * pair of brightness cuts, but only one of them means the same thing on two microscopes. What
 * counts as deep enough to separate two peaks depends on where the noise floor of the averaged
 * image sits, and that moves with the camera, the gain and the number of frames averaged.
 *
 * So the number is measured rather than set. These check the measurement - that it recovers a
 * known noise level, that neither spots nor an illumination gradient move it, and that the
 * plumbing around it honours an explicitly set value. How well the estimate scores against a
 * fixed 5 is a simulation question, in SIMULATION.md and re-derivable with
 * tools/tolerance_sweep.py.
 */
class SpotToleranceEstimateTest {

    private static final int WIDTH = 256;
    private static final int HEIGHT = 256;

    /** Flat background plus Gaussian noise of a known sigma. */
    private static float[] noisyField(double background, double sigma, long seed) {
        Random random = new Random(seed);
        float[] pixels = new float[WIDTH * HEIGHT];
        for (int i = 0; i < pixels.length; i++) {
            pixels[i] = (float) (background + sigma * random.nextGaussian());
        }
        return pixels;
    }

    /**
     * The measurement itself, against a known answer.
     *
     * Over two orders of magnitude, because the whole point is that it tracks a noise floor
     * that moves - a estimator that were only right near 1.4 would reproduce the old fixed
     * default and fix nothing.
     */
    @Test
    @DisplayName("it recovers a known noise level")
    void recoversKnownNoise() {
        for (double sigma : new double[] {0.05, 0.5, 1.4, 5.0}) {
            double measured = smFRETSpotFinder.pixelNoise(
                    noisyField(100.0, sigma, 7L), WIDTH, HEIGHT);
            assertEquals(sigma, measured, 0.03 * sigma,
                    "measured " + measured + " for a true sigma of " + sigma);
        }
    }

    /**
     * Spots must not raise it, which is why this is a MAD of differences rather than a standard
     * deviation. A plain sigma over a field of spots measures the spots, and the tolerance
     * derived from it would grow with the density of the sample.
     */
    @Test
    @DisplayName("spots do not raise the estimate")
    void spotsDoNotCount() {
        float[] clean = noisyField(100.0, 1.0, 11L);
        double before = smFRETSpotFinder.pixelNoise(clean, WIDTH, HEIGHT);

        // 400 bright spots, far above anything the noise does.
        float[] withSpots = clean.clone();
        Random random = new Random(3L);
        for (int i = 0; i < 400; i++) {
            int x = 2 + random.nextInt(WIDTH - 4);
            int y = 2 + random.nextInt(HEIGHT - 4);
            for (int dy = -1; dy <= 1; dy++) {
                for (int dx = -1; dx <= 1; dx++) {
                    withSpots[(y + dy) * WIDTH + (x + dx)] += 200.0f;
                }
            }
        }
        double after = smFRETSpotFinder.pixelNoise(withSpots, WIDTH, HEIGHT);

        assertEquals(before, after, 0.05 * before,
                "spots moved the noise estimate from " + before + " to " + after);
    }

    /**
     * Nor may the illumination profile. The fields these run on are lit by a Gaussian beam wide
     * enough that the background smoothing has something to do, and differencing is what removes
     * it - a estimator taken on the pixel values themselves would read the beam, not the noise.
     */
    @Test
    @DisplayName("an illumination gradient does not raise the estimate")
    void gradientDoesNotCount() {
        float[] flat = noisyField(100.0, 1.0, 13L);
        double before = smFRETSpotFinder.pixelNoise(flat, WIDTH, HEIGHT);

        float[] lit = flat.clone();
        for (int y = 0; y < HEIGHT; y++) {
            for (int x = 0; x < WIDTH; x++) {
                double r2 = Math.pow(x - WIDTH / 2.0, 2) + Math.pow(y - HEIGHT / 2.0, 2);
                lit[y * WIDTH + x] += (float) (80.0 * Math.exp(-r2 / (2 * 180.0 * 180.0)));
            }
        }
        double after = smFRETSpotFinder.pixelNoise(lit, WIDTH, HEIGHT);

        assertEquals(before, after, 0.05 * before,
                "the beam moved the noise estimate from " + before + " to " + after);
    }

    /** An image with nothing to measure returns 0, which the caller turns into the fallback. */
    @Test
    @DisplayName("a uniform or degenerate image reports no noise")
    void nothingToMeasure() {
        assertEquals(0.0, smFRETSpotFinder.pixelNoise(new float[WIDTH * HEIGHT], WIDTH, HEIGHT),
                1.0e-12, "a uniform image");
        assertEquals(0.0, smFRETSpotFinder.pixelNoise(new float[4], 1, 4), 1.0e-12,
                "one pixel wide, nothing to difference");
    }

    /** The tolerance actually used: estimated at 0, and the set value otherwise. */
    @Test
    @DisplayName("zero estimates, anything else is taken as given")
    void zeroMeansEstimate() {
        smFRETSpotFinder finder = new smFRETSpotFinder();
        finder.log = new StderrLogService();

        smFRETSpotFinder.Shared image = new smFRETSpotFinder.Shared(WIDTH, HEIGHT);
        System.arraycopy(noisyField(100.0, 2.0, 17L), 0, image.pixels, 0, image.pixels.length);

        finder.spotTolerance = 0.0;
        double estimated = finder.effectiveTolerance(image);
        double wanted = smFRETSpotFinder.TOLERANCE_SIGMAS * 2.0;
        assertEquals(wanted, estimated, 0.05 * wanted,
                "estimated " + estimated + " where the noise is 2.0");

        finder.spotTolerance = 8.0;
        assertEquals(8.0, finder.effectiveTolerance(image), 1.0e-9, "a set value is used as is");

        // A flat image has no noise to measure, and handing MaximumFinder a zero would return
        // every maximum in the field.
        finder.spotTolerance = 0.0;
        assertEquals(smFRETSpotFinder.TOLERANCE_FALLBACK,
                finder.effectiveTolerance(new smFRETSpotFinder.Shared(WIDTH, HEIGHT)), 1.0e-9,
                "an unmeasurable image falls back");
    }

    /**
     * The old preference key is ignored, and that is the whole point of renaming it.
     *
     * Everyone who has ever run this plugin has 5.0 saved under "spotTolerance" - the default
     * while it was set by hand. Restoring it would pin every existing user to a hand-set value
     * they never chose and hide the estimate behind it, so the estimate would ship and nobody
     * would get it.
     */
    @Test
    @DisplayName("a preference file from before the change does not pin the old default")
    void theOldKeyIsNotRestored() {
        Map<String, String> old = new HashMap<>();
        old.put("spotTolerance", "5.0");

        smFRETSpotFinder fresh = new smFRETSpotFinder();
        fresh.applySettings(old);
        assertEquals(0.0, fresh.spotTolerance, 1.0e-9, "the stale key was restored");

        // The new key still round trips, so an explicit override survives a restart.
        smFRETSpotFinder set = new smFRETSpotFinder();
        set.spotTolerance = 7.5;
        smFRETSpotFinder restored = new smFRETSpotFinder();
        restored.applySettings(set.currentSettings());
        assertEquals(7.5, restored.spotTolerance, 1.0e-9, "an explicit value did not survive");
        assertTrue(set.currentSettings().containsKey("spotToleranceAuto"),
                "saved under the new key");
    }
}
