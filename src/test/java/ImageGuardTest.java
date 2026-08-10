import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ShortProcessor;
import java.awt.image.IndexColorModel;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

/**
 * What the three analysis stages refuse to measure, and why the image type alone is not the test.
 *
 * Every stage assumes one number per pixel. An RGB image packs three into an int and an indexed
 * colour image stores palette entries, so in both cases the arithmetic would run and produce
 * numbers that mean nothing - which is worse than an error.
 *
 * The subtle case is the one these exist for: a colour indexed PNG or GIF opens through ImagePlus
 * as GRAY8 rather than COLOR_256, so the type says nothing and only the processor's lookup table
 * gives it away. That test is deliberately *not* applied above 8 bits, where a colour LUT is a
 * display choice over real intensities.
 */
class ImageGuardTest {

    /** A palette that is not a grey ramp, which is what isColorLut keys off. */
    private static IndexColorModel colourPalette() {
        byte[] reds = new byte[256];
        byte[] greens = new byte[256];
        byte[] blues = new byte[256];
        for (int i = 0; i < 256; i++) {
            reds[i] = (byte) i;
            greens[i] = (byte) (255 - i);
            blues[i] = (byte) ((i * 7) % 256);
        }
        return new IndexColorModel(8, 256, reds, greens, blues);
    }

    @Test
    @DisplayName("8, 16 and 32 bit grayscale are accepted")
    void grayscaleIsAccepted() {
        assertDoesNotThrow(() -> smFRETChannelMapper.requireGrayscale(
                new ImagePlus("8", new ByteProcessor(4, 4)), "the image"));
        assertDoesNotThrow(() -> smFRETChannelMapper.requireGrayscale(
                new ImagePlus("16", new ShortProcessor(4, 4)), "the image"));
        assertDoesNotThrow(() -> smFRETChannelMapper.requireGrayscale(
                new ImagePlus("32", new FloatProcessor(4, 4)), "the image"));
    }

    @Test
    @DisplayName("an RGB image is refused")
    void rgbIsRefused() {
        smFRETAnalysisException thrown = assertThrows(smFRETAnalysisException.class,
                () -> smFRETChannelMapper.requireGrayscale(
                        new ImagePlus("rgb", new ColorProcessor(4, 4)), "the movie"));
        assertTrue(thrown.getMessage().contains("RGB"), thrown.getMessage());
        assertTrue(thrown.getMessage().contains("the movie"),
                "the message should name what was being opened: " + thrown.getMessage());
    }

    /**
     * The case the type check misses. The message has to say the values may be palette indices
     * and how to clear the table, because the alternative - refusing silently or measuring the
     * indices - are both worse than a recoverable error.
     */
    @Test
    @DisplayName("8 bit with a colour lookup table is refused, with the remedy")
    void colourIndexedEightBitIsRefused() {
        ByteProcessor indexed = new ByteProcessor(4, 4);
        indexed.setColorModel(colourPalette());
        ImagePlus image = new ImagePlus("indexed", indexed);

        smFRETAnalysisException thrown = assertThrows(smFRETAnalysisException.class,
                () -> smFRETChannelMapper.requireGrayscale(image, "the movie"));
        assertTrue(thrown.getMessage().contains("palette indices"), thrown.getMessage());
        assertTrue(thrown.getMessage().contains("Lookup Tables"),
                "the message should say how to clear the table: " + thrown.getMessage());
    }

    /**
     * The deliberate asymmetry. A 16 bit image cannot be palette indexed, so a colour LUT there
     * is display only and no reason to refuse real intensities.
     */
    @Test
    @DisplayName("a colour lookup table above 8 bits is display only and allowed")
    void colourLutAboveEightBitsIsAllowed() {
        ShortProcessor sixteenBit = new ShortProcessor(4, 4);
        sixteenBit.setColorModel(colourPalette());
        assertDoesNotThrow(() -> smFRETChannelMapper.requireGrayscale(
                new ImagePlus("16 bit with a LUT", sixteenBit), "the image"));
    }

    @Test
    @DisplayName("an unreadable image is reported rather than dereferenced")
    void nullImageIsReported() {
        assertThrows(smFRETAnalysisException.class,
                () -> smFRETChannelMapper.requireGrayscale(null, "the image"));
        assertThrows(smFRETAnalysisException.class,
                () -> smFRETChannelMapper.requireGrayscale(new ImagePlus(), "the image"));
    }

    /**
     * A stack with one non-singleton axis is taken as time whichever way ImageJ labelled it,
     * because ImageJ cannot tell a movie from a depth stack saved as a plain TIFF - the example
     * data reports z=30, t=1 and so does a genuine 30 slice volume. Demanding t>1 would refuse
     * every movie this plugin has ever run on.
     */
    @Test
    @DisplayName("a plain stack is accepted whichever axis it claims")
    void singleAxisStacksAreAccepted() {
        ImageStack stack = new ImageStack(4, 4);
        for (int i = 0; i < 6; i++) {
            stack.addSlice(new ByteProcessor(4, 4));
        }

        ImagePlus asDepth = new ImagePlus("z=6", stack);
        assertDoesNotThrow(() -> smFRETChannelMapper.requireTimeStack(asDepth, "the movie"));

        ImagePlus asTime = new ImagePlus("t=6", stack);
        asTime.setDimensions(1, 1, 6);
        assertDoesNotThrow(() -> smFRETChannelMapper.requireTimeStack(asTime, "the movie"));
    }

    @Test
    @DisplayName("a stack with both depth and time is refused")
    void fourDimensionalStacksAreRefused() {
        ImageStack stack = new ImageStack(4, 4);
        for (int i = 0; i < 6; i++) {
            stack.addSlice(new ByteProcessor(4, 4));
        }
        ImagePlus fourD = new ImagePlus("zt", stack);
        fourD.setDimensions(1, 2, 3);

        smFRETAnalysisException thrown = assertThrows(smFRETAnalysisException.class,
                () -> smFRETChannelMapper.requireTimeStack(fourD, "the movie"));
        assertTrue(thrown.getMessage().contains("Hyperstack to Stack"), thrown.getMessage());
    }

    /**
     * More than one channel means the two FRET channels are stored separately rather than side by
     * side within each frame, so the movie would be measured as one twice as long as it is.
     */
    @Test
    @DisplayName("a multi channel stack is refused, and told how to split")
    void multiChannelStacksAreRefused() {
        ImageStack stack = new ImageStack(4, 4);
        for (int i = 0; i < 6; i++) {
            stack.addSlice(new ByteProcessor(4, 4));
        }
        ImagePlus twoChannel = new ImagePlus("c", stack);
        twoChannel.setDimensions(2, 1, 3);

        smFRETAnalysisException thrown = assertThrows(smFRETAnalysisException.class,
                () -> smFRETChannelMapper.requireTimeStack(twoChannel, "the movie"));
        assertTrue(thrown.getMessage().contains("Split Channels"), thrown.getMessage());
    }
}
