import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

/**
 * The background smoothing, and the coarse grid it is really done on.
 *
 * Smoothing at sigma 14 was ~70% of stage 3, so maskedSmooth bins by 4, blurs at 3.5 on the
 * coarse grid and interpolates back - worth 209 s to 103 s on a 1295 frame movie. What has to
 * survive that is the thing the estimator is for: a normalized convolution returns the mean of
 * the pixels it was told to trust, and over a field where every trusted pixel holds the same
 * value that mean is that value, everywhere, including at the edges where a plain blur would
 * droop toward zero.
 *
 * That invariant is exact, which makes it a much sharper test than comparing against a
 * full resolution blur: binning, blurring, interpolating and dividing all have to be right for a
 * constant to come back a constant.
 */
class DecimatedSmoothingTest {

    private static final int WIDTH = 200;
    private static final int HEIGHT = 160;
    private static final double SIGMA = 14.0;

    private static smFRETSpotFinder.Shared filled(int width, int height, float value) {
        smFRETSpotFinder.Shared image = new smFRETSpotFinder.Shared(width, height);
        java.util.Arrays.fill(image.pixels, value);
        return image;
    }

    private static boolean[] allTrusted(int count) {
        boolean[] keep = new boolean[count];
        java.util.Arrays.fill(keep, true);
        return keep;
    }

    @Test
    @DisplayName("a constant field smooths to that constant, edges included")
    void constantFieldIsPreserved() {
        smFRETSpotFinder finder = new smFRETSpotFinder();
        smFRETSpotFinder.Shared image = filled(WIDTH, HEIGHT, 17.5f);

        smFRETSpotFinder.Shared smoothed = finder.maskedSmooth(image, allTrusted(image.pixels.length),
                SIGMA, new float[image.pixels.length]);

        for (int y = 0; y < HEIGHT; y++) {
            for (int x = 0; x < WIDTH; x++) {
                assertEquals(17.5, smoothed.pixels[y * WIDTH + x], 1.0e-3,
                        "at " + x + "," + y);
            }
        }
    }

    /**
     * The same, with most of the field masked out. Every *trusted* pixel still holds 17.5, so the
     * weighted mean of the trusted pixels is still 17.5 - the masked pixels contribute to neither
     * the numerator nor the denominator, which is the whole point of dividing the two blurs.
     */
    @Test
    @DisplayName("masked pixels do not pull the estimate toward zero")
    void maskedPixelsAreIgnoredRatherThanCountedAsZero() {
        smFRETSpotFinder finder = new smFRETSpotFinder();
        smFRETSpotFinder.Shared image = filled(WIDTH, HEIGHT, 17.5f);

        // Everything a spot would have masked out, plus arbitrary contamination in the gaps, so
        // that counting a masked pixel as a zero would be obvious.
        boolean[] keep = allTrusted(image.pixels.length);
        for (int y = 0; y < HEIGHT; y += 7) {
            for (int x = 0; x < WIDTH; x += 5) {
                keep[y * WIDTH + x] = false;
                image.pixels[y * WIDTH + x] = 900.0f;
            }
        }

        smFRETSpotFinder.Shared smoothed = finder.maskedSmooth(image, keep, SIGMA,
                new float[image.pixels.length]);

        for (int i = 0; i < smoothed.pixels.length; i++) {
            assertEquals(17.5, smoothed.pixels[i], 1.0e-3, "at index " + i);
        }
    }

    /**
     * A linear ramp is what catches a centring error, and this is the reason unbin offsets by
     * half a block. The binned sample stands for the *centre* of the block it summed, half a
     * block in from its corner; interpolating as though it stood for the corner shifts the whole
     * estimate by (decimation - 1) / 2 pixels, which on a ramp is a constant offset of that many
     * units and on real data is a background subtracted from the wrong place.
     *
     * Only the interior is checked. Near an edge the mean of the available neighbours genuinely
     * is not the centre value for a ramp, so the estimator is right to differ there.
     */
    @Test
    @DisplayName("a ramp comes back unshifted")
    void aRampIsNotShifted() {
        smFRETSpotFinder finder = new smFRETSpotFinder();
        smFRETSpotFinder.Shared image = new smFRETSpotFinder.Shared(WIDTH, HEIGHT);
        for (int y = 0; y < HEIGHT; y++) {
            for (int x = 0; x < WIDTH; x++) {
                image.pixels[y * WIDTH + x] = x;
            }
        }

        smFRETSpotFinder.Shared smoothed = finder.maskedSmooth(image, allTrusted(image.pixels.length),
                SIGMA, new float[image.pixels.length]);

        int margin = (int) Math.ceil(3.0 * SIGMA);
        for (int y = margin; y < (HEIGHT - margin); y++) {
            for (int x = margin; x < (WIDTH - margin); x++) {

                // A shift of (4-1)/2 = 1.5 would show here as 1.5, so 0.05 is a tenth of the
                // smallest error this is meant to catch.
                assertEquals(x, smoothed.pixels[y * WIDTH + x], 0.05, "at " + x + "," + y);
            }
        }
    }

    @ParameterizedTest(name = "decimation={0}")
    @ValueSource(ints = {2, 3, 4, 8})
    @DisplayName("bin sums each block, short blocks included")
    void binSumsBlocks(int decimation) {

        // Deliberately not a multiple of the decimation, so the last row and column of blocks
        // hang off the edge and are short. The weight image is what accounts for that downstream.
        int width = 4 * decimation + 1;
        int height = 3 * decimation + 2;
        smFRETSpotFinder.Shared image = filled(width, height, 1.0f);

        smFRETSpotFinder.Shared binned = smFRETSpotFinder.bin(image, decimation);

        assertEquals((width + decimation - 1) / decimation, binned.width, "binned width");
        assertEquals((height + decimation - 1) / decimation, binned.height, "binned height");

        // Every pixel of the input lands in exactly one block, so the totals must agree.
        double inputTotal = 0.0;
        for (float value : image.pixels) {
            inputTotal += value;
        }
        double binnedTotal = 0.0;
        for (float value : binned.pixels) {
            binnedTotal += value;
        }
        assertEquals(inputTotal, binnedTotal, 1.0e-6, "no pixel lost or double counted");

        // A full interior block holds decimation^2 ones.
        assertEquals(decimation * decimation, binned.pixels[0], 1.0e-6, "the first block");
    }

    @ParameterizedTest(name = "decimation={0}")
    @ValueSource(ints = {2, 3, 4, 8})
    @DisplayName("unbin of a constant is that constant")
    void unbinOfConstantIsConstant(int decimation) {
        smFRETSpotFinder.Shared binned = filled(11, 9, 5.0f);
        smFRETSpotFinder.Shared target =
                new smFRETSpotFinder.Shared(11 * decimation, 9 * decimation);

        smFRETSpotFinder.unbin(binned, decimation, target);

        for (int i = 0; i < target.pixels.length; i++) {
            assertEquals(5.0, target.pixels[i], 1.0e-6, "at index " + i);
        }
    }

    /**
     * Positions off either end clamp to the last sample rather than extrapolating, matching
     * extendBorder - so the weights are always a convex combination and unbin can never overshoot
     * the range of its input.
     */
    @Test
    @DisplayName("interpolation weights stay inside the binned grid")
    void interpolationWeightsClampAtTheEnds() {
        int size = 41;
        int decimation = 4;
        int binnedSize = (size + decimation - 1) / decimation;

        int[] lower = new int[size];
        int[] upper = new int[size];
        float[] weight = new float[size];
        smFRETSpotFinder.interpolationWeights(size, binnedSize, decimation, lower, upper, weight);

        for (int p = 0; p < size; p++) {
            assertTrue((lower[p] >= 0) && (lower[p] < binnedSize), "lower at " + p + ": " + lower[p]);
            assertTrue((upper[p] >= 0) && (upper[p] < binnedSize), "upper at " + p + ": " + upper[p]);
            assertTrue((weight[p] >= 0.0f) && (weight[p] <= 1.0f), "weight at " + p + ": " + weight[p]);
            assertTrue(upper[p] >= lower[p], "inverted at " + p);
        }

        // The centre of the first block is at (decimation - 1) / 2 = 1.5, so pixel 0 sits before
        // the first sample and clamps onto it exactly.
        assertEquals(0, lower[0]);
        assertEquals(0.0f, weight[0], 0.0f);
    }

    /**
     * The decimation is an approximation, and this is how good it has to stay. Binning is a box
     * prefilter that adds decimation^2 / 12 to the kernel variance - 1.3 against sigma 14's 196 -
     * so the coarse answer and the full resolution one must agree to far less than an ADU on a
     * background of order 20.
     *
     * Note the decimation^2: bin *sums* each block rather than averaging it, so smoothDecimated
     * comes back that many times larger than a blur. maskedSmooth never notices, because it
     * decimates the weighted sum and the weight the same way and the factor cancels when it
     * divides one by the other - but a test comparing this against a blur directly has to put it
     * back.
     */
    @Test
    @DisplayName("the coarse grid agrees with a full resolution blur")
    void decimatedMatchesFullResolution() {
        int decimation = 4;
        double blockScale = decimation * decimation;
        smFRETSpotFinder.Shared image = new smFRETSpotFinder.Shared(WIDTH, HEIGHT);
        for (int y = 0; y < HEIGHT; y++) {
            for (int x = 0; x < WIDTH; x++) {

                // A smoothly varying illumination profile, which is what the background is.
                image.pixels[y * WIDTH + x] = (float) (20.0
                        + 4.0 * Math.sin(x / 60.0) * Math.cos(y / 50.0));
            }
        }

        smFRETSpotFinder.Shared coarse = new smFRETSpotFinder.Shared(WIDTH, HEIGHT);
        smFRETSpotFinder.smoothDecimated(image, SIGMA, decimation, coarse);

        smFRETSpotFinder.Shared full = new smFRETSpotFinder.Shared(WIDTH, HEIGHT);
        smFRETSpotFinder.gauss(SIGMA, net.imglib2.view.Views.extendZero(image.img), full.img);

        double sumSquares = 0.0;
        int count = 0;
        int margin = (int) Math.ceil(3.0 * SIGMA);
        for (int y = margin; y < (HEIGHT - margin); y++) {
            for (int x = margin; x < (WIDTH - margin); x++) {
                double difference = (coarse.pixels[y * WIDTH + x] / blockScale)
                        - full.pixels[y * WIDTH + x];
                sumSquares += difference * difference;
                count++;
            }
        }
        double rms = Math.sqrt(sumSquares / count);
        assertTrue(rms < 0.05, "rms difference against a full resolution blur was " + rms);
    }
}
