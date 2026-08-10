import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import ij.ImagePlus;
import ij.process.FloatProcessor;
import java.util.Collections;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

/**
 * What spotThreshold means: real standard deviations, at any spot size.
 *
 * The number spotFilterSNR reports is meant to be the matched filter significance of the
 * detection, N / (2 sigma sqrt(pi B)) for a spot of N photoelectrons on a background of B per
 * pixel. Two bugs had to be fixed before it was:
 *
 *   - the norm was 2 pi sigma^2 rather than 4 pi sigma^2, leaving the spot's own width out of the
 *     quadrature sum. That recovered N/2 and took the noise term from half the right area, so the
 *     ratio was short of the true significance by exactly sqrt(2) - a threshold of 6 behaved like
 *     8.49;
 *   - the clip bias made the reported value climb with the spot size, +15% from sigma 1 to 3 at
 *     fixed true significance and +46% at fixed molecular brightness, so the same threshold meant
 *     something different on every movie.
 *
 * Both are pinned here against the closed form, which is the point: this test knows what the
 * answer should be rather than what the code last said.
 */
class SnrCalibrationTest {

    private static final int SIZE = 128;
    private static final int CENTRE = SIZE / 2;
    private static final double BACKGROUND = 20.0;

    /** N / (2 sigma sqrt(pi B)) - the matched filter significance the report is meant to be. */
    private static double trueSignificance(double photons, double sigma, double background) {
        return photons / (2.0 * sigma * Math.sqrt(Math.PI * background));
    }

    /**
     * What a camera's pixels cost, and it is not a calibration error.
     *
     * A real detector integrates the PSF over each pixel rather than sampling it, which lowers
     * the peak of the binned spot below the continuous Gaussian's by 1/(1 + 1/24 sigma^2). The
     * spot finder measures that peak, so its report is short of the continuous significance by
     * this factor - 4% at spotSigma 1, 0.5% at spotSigma 3. It is the same factor that keeps
     * trace recovery just under 1.
     *
     * Dividing it out is what turns "the report drifts by 4% across the spot size range" into the
     * question actually worth asking, which is whether anything *else* drifts.
     */
    private static double pixelIntegration(double sigma) {
        return 1.0 / (1.0 + 1.0 / (24.0 * sigma * sigma));
    }

    /** The photon count that gives a particular true significance at this spot size. */
    private static double photonsFor(double significance, double sigma, double background) {
        return significance * 2.0 * sigma * Math.sqrt(Math.PI * background);
    }

    private static smFRETSpotFinder finderAt(double sigma) {
        smFRETSpotFinder finder = new smFRETSpotFinder();
        finder.spotSigma = sigma;
        finder.cameraGain = 1.0;
        finder.cameraBlackLevel = 0;
        finder.spotThreshold = 6.0;
        return finder;
    }

    /**
     * The measured SNR of one spot of a known photon count on a known flat background.
     *
     * Noiseless and with the background handed in exactly rather than estimated, so what is being
     * measured is the arithmetic of the report and not the estimator underneath it.
     */
    private static double reportedSnr(double sigma, double photons) {
        float[] pixels = SyntheticField.half(SIZE, SIZE,
                Collections.singletonList(new SyntheticField.Spot(CENTRE, CENTRE, photons)),
                sigma, BACKGROUND);
        ImagePlus image = new ImagePlus("field", new FloatProcessor(SIZE, SIZE, pixels, null));

        float[] flat = new float[SIZE * SIZE];
        java.util.Arrays.fill(flat, (float) BACKGROUND);
        ImagePlus background = new ImagePlus("bg", new FloatProcessor(SIZE, SIZE, flat, null));

        double[][] filtered = finderAt(sigma).spotFilterSNR(
                new double[][] {{1.0, CENTRE, CENTRE}}, image, background);
        return filtered[0][filtered[0].length - 1];
    }

    /**
     * The closed form, at the default spot size. The shortfall against it is pixel integration -
     * the same 1/(1 + 1/24 sigma^2) that keeps trace recovery just under 1 - so the tolerance is
     * a couple of percent rather than exact.
     */
    @Test
    @DisplayName("the reported SNR is the matched filter significance")
    void reportedSnrIsTheTrueSignificance() {
        double sigma = 2.0;
        double photons = 2000.0;

        double expected = trueSignificance(photons, sigma, BACKGROUND);
        double measured = reportedSnr(sigma, photons);

        assertEquals(expected * pixelIntegration(sigma), measured, 0.01 * expected,
                "reported " + measured + " against a true significance of " + expected);

        // The pre-fix norm would have reported this over sqrt(2). Stated separately because it is
        // the specific regression, and it is far outside the tolerance above.
        assertTrue(measured > 1.2 * (expected / Math.sqrt(2.0)),
                "reported " + measured + " is close to the pre norm fix value of "
                        + (expected / Math.sqrt(2.0)));
    }

    /**
     * Size independence, which is what lets spotThreshold be one number for every movie. A spot
     * at a *fixed true significance* has to report that significance whatever its width - the
     * drift was 15% before the truncation correction and is 0.5% after it.
     */
    @ParameterizedTest(name = "spotSigma={0}")
    @ValueSource(doubles = {1.0, 1.5, 2.0, 2.5, 3.0})
    @DisplayName("a spot of fixed significance reports it at any spot size")
    void reportIsSizeIndependent(double sigma) {
        double target = 8.0;
        double measured = reportedSnr(sigma, photonsFor(target, sigma, BACKGROUND));
        double expected = target * pixelIntegration(sigma);

        assertEquals(expected, measured, 0.01 * target,
                "at spotSigma " + sigma + " a true 8 sigma spot reported " + measured
                        + " against " + expected + " predicted by pixel integration alone");
    }

    /**
     * The claim in its strongest form: once the known pixel integration factor is divided out,
     * nothing else varies with the spot size. That is what makes spotThreshold one number for
     * every movie - before the truncation correction this drifted 15%, and the residual here is
     * what is left.
     */
    @Test
    @DisplayName("no drift with spot size beyond pixel integration")
    void thereIsNoDriftWithSpotSize() {
        double target = 8.0;
        double lowest = Double.MAX_VALUE;
        double highest = -Double.MAX_VALUE;

        for (double sigma = 1.0; sigma <= 3.0; sigma += 0.25) {
            double corrected = reportedSnr(sigma, photonsFor(target, sigma, BACKGROUND))
                    / pixelIntegration(sigma);
            lowest = Math.min(lowest, corrected);
            highest = Math.max(highest, corrected);
        }

        double drift = (highest - lowest) / target;
        assertTrue(drift < 0.02,
                "corrected SNR spread " + lowest + " to " + highest + " over spotSigma 1 to 3,"
                        + " a residual drift of " + (100.0 * drift) + "%");
    }

    @Test
    @DisplayName("SNR scales with the photon count and falls with the background")
    void snrScalesTheRightWay() {
        double sigma = 2.0;

        assertEquals(2.0, reportedSnr(sigma, 4000.0) / reportedSnr(sigma, 2000.0), 0.02,
                "twice the photons should be twice the significance");

        // The filter records the score whether or not the spot passed, and clears the flag when
        // it did not - the score is a diagnostic column in the CSV and has to survive rejection.
        smFRETSpotFinder finder = finderAt(sigma);
        finder.spotThreshold = 1000.0;

        float[] pixels = SyntheticField.half(SIZE, SIZE,
                Collections.singletonList(new SyntheticField.Spot(CENTRE, CENTRE, 2000.0)),
                sigma, BACKGROUND);
        float[] flat = new float[SIZE * SIZE];
        java.util.Arrays.fill(flat, (float) BACKGROUND);

        double[][] filtered = finder.spotFilterSNR(new double[][] {{1.0, CENTRE, CENTRE}},
                new ImagePlus("field", new FloatProcessor(SIZE, SIZE, pixels, null)),
                new ImagePlus("bg", new FloatProcessor(SIZE, SIZE, flat, null)));

        assertEquals(0.0, filtered[0][0], 0.0, "the good flag should be cleared");
        assertTrue(filtered[0][filtered[0].length - 1] > 0.0, "the score should still be recorded");
    }

    /**
     * The black level is multiplied by how many channels went into the image spots are found in -
     * two for the sum, one otherwise. On a single channel the old fixed factor of two
     * over-subtracted a whole black level and shifted every SNR.
     */
    @Test
    @DisplayName("the black level is counted once per channel in the analysis image")
    void blackLevelFollowsTheChannelCount() {
        smFRETSpotFinder finder = new smFRETSpotFinder();

        finder.spotChannel = "sum";
        assertEquals(2, finder.channelCount());

        finder.spotChannel = "donor";
        assertEquals(1, finder.channelCount());

        finder.spotChannel = "acceptor";
        assertEquals(1, finder.channelCount());
    }
}
