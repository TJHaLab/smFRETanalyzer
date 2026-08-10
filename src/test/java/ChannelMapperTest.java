import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import java.io.File;
import net.imglib2.img.display.imagej.ImageJFunctions;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.scijava.log.StderrLogService;

/**
 * The mapping: solving an affine from three landmark pairs, and applying it.
 *
 * TurboReg measures the mapping and imglib2 applies it, which is only safe because the three
 * pairs determine the affine *exactly* - it is a solve by Cramer's rule, not a fit - so the
 * transform reconstructed from a mapping file is the one TurboReg found, to 2e-13 px. These pin
 * that solve, and the one behaviour of the warp that everything downstream depends on: a pixel
 * whose inverse mapped coordinate falls outside the source is left at zero, which is what
 * createOverlapMask reads to find the region the two channels share.
 */
class ChannelMapperTest {

    private static smFRETChannelMapper mapper() {
        smFRETChannelMapper mapper = new smFRETChannelMapper();
        mapper.log = new StderrLogService();
        return mapper;
    }

    /** Apply a solved affine the way AffineTransform2D.set reads it: two rows of three. */
    private static double[] apply(double[] affine, double x, double y) {
        return new double[] {
            affine[0] * x + affine[1] * y + affine[2],
            affine[3] * x + affine[4] * y + affine[5],
        };
    }

    @Test
    @DisplayName("coincident landmarks solve to the identity")
    void coincidentLandmarksGiveIdentity() {
        double[][] points = {{10.0, 20.0}, {80.0, 30.0}, {40.0, 90.0}};
        double[] affine = smFRETChannelMapper.solveAffine(points, points);

        assertEquals(1.0, affine[0], 1.0e-12);
        assertEquals(0.0, affine[1], 1.0e-12);
        assertEquals(0.0, affine[2], 1.0e-12);
        assertEquals(0.0, affine[3], 1.0e-12);
        assertEquals(1.0, affine[4], 1.0e-12);
        assertEquals(0.0, affine[5], 1.0e-12);
    }

    /**
     * A solve, not a fit: an affine put through three pairs must come back exactly, with no
     * residual to speak of. 1e-10 on coordinates of order 100 is 1e-12 relative, which is
     * double precision arithmetic rather than an approximation.
     */
    @Test
    @DisplayName("an arbitrary affine round trips through its landmarks")
    void arbitraryAffineRoundTrips() {
        double angle = 0.13;
        double scale = 1.04;
        double[] truth = {
            scale * Math.cos(angle), -scale * Math.sin(angle), 7.25,
            scale * Math.sin(angle), scale * Math.cos(angle), -3.5,
        };

        double[][] source = {{10.0, 20.0}, {80.0, 30.0}, {40.0, 90.0}};
        double[][] target = new double[3][];
        for (int i = 0; i < 3; i++) {
            target[i] = apply(truth, source[i][0], source[i][1]);
        }

        double[] solved = smFRETChannelMapper.solveAffine(source, target);
        for (int i = 0; i < 6; i++) {
            assertEquals(truth[i], solved[i], 1.0e-10, "coefficient " + i);
        }
    }

    @Test
    @DisplayName("collinear landmarks are refused rather than solved")
    void collinearLandmarksThrow() {
        double[][] collinear = {{0.0, 0.0}, {10.0, 10.0}, {20.0, 20.0}};
        smFRETAnalysisException thrown = assertThrows(smFRETAnalysisException.class,
                () -> smFRETChannelMapper.solveAffine(collinear, collinear));
        assertTrue(thrown.getMessage().contains("collinear"), thrown.getMessage());
    }

    @Test
    @DisplayName("an identity mapping warps an image to itself")
    void identityWarpIsANoOp(@TempDir File directory) {
        int width = 64;
        int height = 48;
        File mappingFile = SyntheticField.writeIdentityMapping(directory, "identity.json",
                2 * width, height);

        smFRETChannelMapper mapper = mapper();
        mapper.loadMappingJSON(mappingFile.toString());

        float[] pixels = new float[width * height];
        for (int i = 0; i < pixels.length; i++) {
            pixels[i] = (i % 97) + 0.25f;
        }
        ImagePlus source = new ImagePlus("source", new FloatProcessor(width, height, pixels, null));

        ImagePlus warped = mapper.transformImagePlus(
                ImageJFunctions.convertFloat(source), width, height);

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                assertEquals(pixels[y * width + x], warped.getProcessor().getf(x, y), 1.0e-4,
                        "at " + x + "," + y);
            }
        }
    }

    /**
     * The out of range rule, and it is on the inverse mapped *coordinate* rather than on the
     * sampled value. A target pixel is written only when it lands strictly inside the source and
     * is left at zero otherwise, which is exactly what createOverlapMask keys off - it adds the
     * two warped halves and keeps only where both contributed. Extending with zero and letting
     * the interpolator blend across the boundary instead moved a three pixel band by up to
     * 19 ADU.
     */
    @Test
    @DisplayName("pixels mapping off the source are zero, not extrapolated")
    void outOfRangePixelsAreZero(@TempDir File directory) throws Exception {
        int width = 40;
        int height = 32;
        int shift = 5;

        // target = source + shift in x, so a target pixel below x = shift inverse maps to a
        // negative source x and has nothing to read.
        double[][] source = {{8.0, 6.0}, {30.0, 9.0}, {15.0, 25.0}};
        double[][] target = new double[3][];
        for (int i = 0; i < 3; i++) {
            target[i] = new double[] {source[i][0] + shift, source[i][1]};
        }

        java.util.Map<String, Object> mapping = new java.util.HashMap<>();
        mapping.put("source points", source);
        mapping.put("target points", target);
        mapping.put("image width", 2 * width);
        mapping.put("image height", height);
        File mappingFile = new File(directory, "shifted.json");
        new com.fasterxml.jackson.databind.ObjectMapper().writeValue(mappingFile, mapping);

        smFRETChannelMapper mapper = mapper();
        mapper.loadMappingJSON(mappingFile.toString());

        float[] pixels = new float[width * height];
        java.util.Arrays.fill(pixels, 100.0f);
        ImagePlus warped = mapper.transformImagePlus(
                ImageJFunctions.convertFloat(
                        new ImagePlus("source", new FloatProcessor(width, height, pixels, null))),
                width, height);

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < shift; x++) {
                assertEquals(0.0, warped.getProcessor().getf(x, y), 0.0,
                        "off the source at " + x + "," + y);
            }
            for (int x = shift; x < width; x++) {
                assertEquals(100.0, warped.getProcessor().getf(x, y), 1.0e-4,
                        "inside the source at " + x + "," + y);
            }
        }
    }

    @Test
    @DisplayName("a mapping file with the key but not the contents says so")
    void malformedMappingIsReported(@TempDir File directory) throws Exception {
        java.util.Map<String, Object> mapping = new java.util.HashMap<>();
        mapping.put("source points", new double[][] {{1.0, 2.0}});     // one pair, not three
        mapping.put("target points", new double[][] {{1.0, 2.0}});
        mapping.put("image width", 100);
        mapping.put("image height", 50);
        File file = new File(directory, "short.json");
        new com.fasterxml.jackson.databind.ObjectMapper().writeValue(file, mapping);

        smFRETAnalysisException thrown = assertThrows(smFRETAnalysisException.class,
                () -> mapper().loadMappingJSON(file.toString()));
        assertTrue(thrown.getMessage().contains("contents could not be read"), thrown.getMessage());
    }

    @Test
    @DisplayName("the average is the mean over the requested range")
    void averageIsTheMeanOverTheRange() {
        ImageStack stack = new ImageStack(2, 2);
        for (int i = 0; i < 10; i++) {
            float[] pixels = new float[4];
            java.util.Arrays.fill(pixels, i);
            stack.addSlice(new FloatProcessor(2, 2, pixels, null));
        }
        ImagePlus movie = new ImagePlus("movie", stack);

        // Slices 3 to 6, one based and inclusive, so frames 2..5 hold 2,3,4,5 and mean 3.5.
        ImagePlus averaged = mapper().averageImagePlus(movie, 3, 6);
        assertEquals(3.5, averaged.getProcessor().getf(0, 0), 1.0e-6);

        // Always float out, whatever went in - the quantization back into the input's type is
        // what used to put a rounding step between every stage.
        assertTrue(averaged.getProcessor() instanceof FloatProcessor,
                "averaged image was " + averaged.getProcessor().getClass().getSimpleName());
    }

    @Test
    @DisplayName("a single frame movie averages to itself")
    void singleFrameIsReturnedUnchanged() {
        ImagePlus single = new ImagePlus("single",
                new FloatProcessor(2, 2, new float[] {1.0f, 2.0f, 3.0f, 4.0f}, null));
        ImagePlus averaged = mapper().averageImagePlus(single, 1, 30);
        assertEquals(1.0, averaged.getProcessor().getf(0, 0), 0.0);
        assertEquals(4.0, averaged.getProcessor().getf(1, 1), 0.0);
    }

    /**
     * frameCount reads the stack size rather than getNSlices, which is what makes a stack saved
     * with t=30, z=1 work - that used to report one frame, average nothing, and analyse frame 1
     * alone.
     */
    @Test
    @DisplayName("the frame count follows whichever axis ImageJ labelled")
    void frameCountReadsTheStackSize() {
        ImageStack stack = new ImageStack(2, 2);
        for (int i = 0; i < 7; i++) {
            stack.addSlice(new ByteProcessor(2, 2));
        }

        ImagePlus asDepth = new ImagePlus("z", stack);
        assertEquals(7, smFRETChannelMapper.frameCount(asDepth), "labelled as slices");

        ImagePlus asTime = new ImagePlus("t", stack);
        asTime.setDimensions(1, 1, 7);
        assertEquals(7, smFRETChannelMapper.frameCount(asTime), "labelled as frames");
        assertNotEquals(7, asTime.getNSlices(), "getNSlices is what this must not be reading");
    }

    @Test
    @DisplayName("toFloat converts once and preserves the values")
    void toFloatPreservesValues() {
        ImageStack stack = new ImageStack(2, 2);
        stack.addSlice(new ByteProcessor(2, 2, new byte[] {0, 17, (byte) 200, (byte) 255}, null));
        ImagePlus eightBit = new ImagePlus("8 bit", stack);

        ImagePlus asFloat = smFRETChannelMapper.toFloat(eightBit);
        assertEquals(32, asFloat.getBitDepth());
        assertEquals(0.0, asFloat.getProcessor().getf(0, 0), 0.0);
        assertEquals(17.0, asFloat.getProcessor().getf(1, 0), 0.0);
        assertEquals(200.0, asFloat.getProcessor().getf(0, 1), 0.0);
        assertEquals(255.0, asFloat.getProcessor().getf(1, 1), 0.0);

        // Already float, so handed straight back rather than copied.
        ImagePlus alreadyFloat = new ImagePlus("32 bit", new FloatProcessor(2, 2));
        assertTrue(smFRETChannelMapper.toFloat(alreadyFloat) == alreadyFloat);
    }
}
