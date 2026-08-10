import static org.junit.jupiter.api.Assertions.assertEquals;

import ij.ImageStack;
import ij.process.FloatProcessor;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

/**
 * The temporal window behind the background estimate.
 *
 * Two things here were wrong before and are worth holding down. It used to call Filters3D with
 * the parameter as a z *radius*, so backgroundAverageNFrames of 30 averaged 61 frames - and on a
 * 30 frame movie that covered the whole movie for every frame, making the background estimate
 * identical on every frame and unable to track anything time varying at all. And its x/y radius
 * of 1 meant it averaged spatially as well, which the parameter never claimed.
 *
 * So: the number is the whole window, the window shrinks at the ends rather than inventing
 * frames, and nothing spatial happens.
 */
class TemporalMeanTest {

    /** A stack of 1x1 frames holding 0, 1, 2, ... so a frame's value is its own index. */
    private static ImageStack counting(int frames) {
        ImageStack stack = new ImageStack(1, 1);
        for (int i = 0; i < frames; i++) {
            stack.addSlice(new FloatProcessor(1, 1, new float[] {i}, null));
        }
        return stack;
    }

    private static double[] meanOf(ImageStack stack, int window) {
        ImageStack averaged = new smFRETAnalyzer().temporalMean(stack, window);
        double[] values = new double[averaged.size()];
        for (int i = 0; i < values.length; i++) {
            values[i] = averaged.getProcessor(i + 1).getf(0, 0);
        }
        return values;
    }

    @Test
    @DisplayName("a window of one leaves the movie alone")
    void windowOfOneIsIdentity() {
        double[] averaged = meanOf(counting(6), 1);
        for (int i = 0; i < 6; i++) {
            assertEquals(i, averaged[i], 1.0e-6, "frame " + i);
        }
    }

    /**
     * Three means three, not seven. This is the radius-versus-width bug stated as a number: under
     * the old behaviour frame 3 would have been the mean of frames 0 to 6.
     */
    @Test
    @DisplayName("the window is a width, not a radius")
    void windowIsAWidth() {
        double[] averaged = meanOf(counting(9), 3);

        // Interior frames: the mean of themselves and one either side, which for a counting
        // sequence is the frame's own index.
        for (int i = 1; i <= 7; i++) {
            assertEquals(i, averaged[i], 1.0e-6, "frame " + i);
        }

        // The ends, where the window has fewer frames to work with rather than invented ones.
        assertEquals(0.5, averaged[0], 1.0e-6, "first frame is the mean of 0 and 1");
        assertEquals(7.5, averaged[8], 1.0e-6, "last frame is the mean of 7 and 8");
    }

    /**
     * An even window cannot sit symmetrically on a frame, so it leans one frame forward: before =
     * (n-1)/2 and after = n-1-before. For a width of 4 that is one before and two after.
     */
    @Test
    @DisplayName("an even window leans forward")
    void evenWindowLeansForward() {
        double[] averaged = meanOf(counting(10), 4);

        // Frame 4 averages 3, 4, 5, 6.
        assertEquals(4.5, averaged[4], 1.0e-6, "frame 4");

        // Frame 0 has nothing before it: 0, 1, 2.
        assertEquals(1.0, averaged[0], 1.0e-6, "frame 0");

        // Frame 9 is the last: 8, 9.
        assertEquals(8.5, averaged[9], 1.0e-6, "frame 9");
    }

    /**
     * The case the old bug was worst in: a window as wide as the movie. It must still vary from
     * frame to frame, because the window shrinks at the ends rather than covering everything
     * everywhere - a background estimate identical on every frame cannot track lamp drift.
     */
    @Test
    @DisplayName("a window as wide as the movie still varies frame to frame")
    void fullWidthWindowStillVaries() {
        int frames = 30;
        double[] averaged = meanOf(counting(frames), frames);

        assertEquals(frames, averaged.length, "one output frame per input frame");

        // before = 14, after = 15, so frame 0 averages 0..15 and the last averages 15..29.
        assertEquals(7.5, averaged[0], 1.0e-6, "frame 0");
        assertEquals(22.0, averaged[frames - 1], 1.0e-6, "the last frame");

        double spread = 0.0;
        for (int i = 1; i < frames; i++) {
            spread += Math.abs(averaged[i] - averaged[i - 1]);
        }
        org.junit.jupiter.api.Assertions.assertTrue(spread > 0.0,
                "the estimate was identical on every frame, which is the Filters3D bug");
    }

    /**
     * Nothing spatial. Each pixel is averaged along time only, so a frame with structure in it
     * keeps that structure exactly when the window is one.
     */
    @Test
    @DisplayName("neighbouring pixels are not averaged together")
    void nothingSpatialHappens() {
        int width = 5;
        int height = 4;
        ImageStack stack = new ImageStack(width, height);
        for (int frame = 0; frame < 3; frame++) {
            float[] pixels = new float[width * height];

            // One bright pixel in a field of zeros. Any spatial averaging spreads it.
            pixels[2 * width + 3] = 100.0f;
            stack.addSlice(new FloatProcessor(width, height, pixels, null));
        }

        ImageStack averaged = new smFRETAnalyzer().temporalMean(stack, 3);
        FloatProcessor middle = (FloatProcessor) averaged.getProcessor(2);

        assertEquals(100.0, middle.getf(3, 2), 1.0e-6, "the bright pixel");
        assertEquals(0.0, middle.getf(2, 2), 1.0e-6, "its left neighbour");
        assertEquals(0.0, middle.getf(4, 2), 1.0e-6, "its right neighbour");
        assertEquals(0.0, middle.getf(3, 1), 1.0e-6, "the pixel above it");
        assertEquals(0.0, middle.getf(3, 3), 1.0e-6, "the pixel below it");
    }
}
