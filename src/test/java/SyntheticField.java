import com.fasterxml.jackson.databind.ObjectMapper;

import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileSaver;
import ij.process.FloatProcessor;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * Two channel fields with a known answer, and the identity mapping that goes with them.
 *
 * The point of generating rather than checking in data is that the integrated intensity of every
 * spot is known exactly, which is what lets the trace measurement be scored rather than merely
 * compared against itself. Everything here is deterministic - noise comes from a seeded Random -
 * so a failure is reproducible.
 *
 * **Spots are pixel *integrated*, not point sampled, and that is deliberate.** A camera collects
 * the light falling across the whole of a pixel, so the value at pixel p of a spot carrying N
 * photons is N times the Gaussian's mass over that pixel, not N times its height at the centre.
 * The difference is the entire reason measureTimeTraces recovers 0.96 N rather than N at
 * spotSigma 1: it reads the peak of a smoothed image and divides by 4 pi sigma^2, which is the
 * right factor for a continuous Gaussian and short of it for a binned one by 1/(1 + 1/24 sigma^2).
 * Point sampling here would hide that and let the tests claim a recovery the pipeline does not
 * actually have.
 *
 * The normal CDF used for the integration is the production one. It is not independently
 * implemented here because it is itself checked against scipy in TruncationConstantsTest, and a
 * second approximation of the same function would only be a second thing to get wrong.
 */
final class SyntheticField {

    private SyntheticField() {
    }

    /** Where a spot is in its half frame, and how many photons it carries in total. */
    static final class Spot {
        final double intensity;
        final double x;
        final double y;

        Spot(double x, double y, double intensity) {
            this.x = x;
            this.y = y;
            this.intensity = intensity;
        }
    }

    /**
     * One channel's pixels: a flat background with pixel integrated Gaussians added.
     *
     * Each spot is only evaluated within five sigma of its centre. Past that the Gaussian's mass
     * in a pixel is below 1e-6 of the peak's, which is far under the noise of any field worth
     * measuring, and evaluating the whole frame per spot makes a crowded field quadratic.
     */
    static float[] half(int width, int height, List<Spot> spots, double sigma, double background) {
        float[] pixels = new float[width * height];
        java.util.Arrays.fill(pixels, (float) background);

        for (Spot spot : spots) {
            int reach = (int) Math.ceil(5.0 * sigma);
            int xLow = Math.max(0, (int) Math.floor(spot.x) - reach);
            int xHigh = Math.min(width - 1, (int) Math.ceil(spot.x) + reach);
            int yLow = Math.max(0, (int) Math.floor(spot.y) - reach);
            int yHigh = Math.min(height - 1, (int) Math.ceil(spot.y) + reach);

            for (int y = yLow; y <= yHigh; y++) {
                double massY = smFRETSpotFinder.normalCdf((y + 0.5 - spot.y) / sigma)
                        - smFRETSpotFinder.normalCdf((y - 0.5 - spot.y) / sigma);
                for (int x = xLow; x <= xHigh; x++) {
                    double massX = smFRETSpotFinder.normalCdf((x + 0.5 - spot.x) / sigma)
                            - smFRETSpotFinder.normalCdf((x - 0.5 - spot.x) / sigma);
                    pixels[y * width + x] += (float) (spot.intensity * massX * massY);
                }
            }
        }
        return pixels;
    }

    /**
     * The two halves side by side, which is the layout every stage here assumes: donor left,
     * acceptor right.
     */
    static FloatProcessor frame(int halfWidth, int height, List<Spot> donor, List<Spot> acceptor,
                                double sigma, double background) {
        float[] left = half(halfWidth, height, donor, sigma, background);
        float[] right = half(halfWidth, height, acceptor, sigma, background);

        float[] pixels = new float[2 * halfWidth * height];
        for (int y = 0; y < height; y++) {
            System.arraycopy(left, y * halfWidth, pixels, y * 2 * halfWidth, halfWidth);
            System.arraycopy(right, y * halfWidth, pixels, y * 2 * halfWidth + halfWidth, halfWidth);
        }
        return new FloatProcessor(2 * halfWidth, height, pixels, null);
    }

    /**
     * A movie of identical frames, optionally with independent noise on each.
     *
     * Identical rather than bleaching, because what these tests measure is a level rather than a
     * time course, and a constant movie makes the expected trace a single number for every frame.
     * Pass noise 0 for an exactly reproducible field with no statistics in it at all.
     */
    static ImagePlus movie(int halfWidth, int height, List<Spot> donor, List<Spot> acceptor,
                           double sigma, double background, int frames, double noise, long seed) {
        Random random = new Random(seed);
        ImageStack stack = new ImageStack(2 * halfWidth, height);

        for (int frame = 0; frame < frames; frame++) {
            FloatProcessor processor = frame(halfWidth, height, donor, acceptor, sigma, background);
            if (noise > 0.0) {
                float[] pixels = (float[]) processor.getPixels();
                for (int i = 0; i < pixels.length; i++) {
                    pixels[i] += (float) (noise * random.nextGaussian());
                }
            }
            stack.addSlice("frame " + (frame + 1), processor);
        }
        return new ImagePlus("synthetic", stack);
    }

    /**
     * A grid of spots far enough apart that none is any other's neighbour.
     *
     * Spots land on integer coordinates on purpose. MaximumFinder reports the integer position of
     * the brightest pixel, and measureTimeTraces reads the smoothed image at exactly that pixel,
     * so a spot centred between pixels would be measured off its own peak and the recovery would
     * carry that miss as well as the effect under test.
     */
    static List<Spot> grid(int width, int height, int spacing, int margin, double intensity) {
        List<Spot> spots = new ArrayList<>();
        for (int y = margin; y <= (height - 1 - margin); y += spacing) {
            for (int x = margin; x <= (width - 1 - margin); x += spacing) {
                spots.add(new Spot(x, y, intensity));
            }
        }
        return spots;
    }

    static File writeMovie(File directory, String name, ImagePlus movie) {
        File file = new File(directory, name);
        boolean saved = new FileSaver(movie).saveAsTiff(file.toString());
        if (!saved) {
            throw new IllegalStateException("could not write " + file);
        }
        return file;
    }

    /**
     * A mapping whose affine is the identity, so the acceptor half lands on the donor's frame
     * unchanged.
     *
     * Three coincident landmark pairs give the identity exactly - solveAffine is a solve rather
     * than a fit - so nothing is approximated here. The recorded width and height are the *whole*
     * two channel frame, which is what splitImagePlus checks the input against.
     *
     * The landmarks are placed in the same triangle smFRETChannelMapper uses, purely so that a
     * mapping written here looks like one written by the real thing.
     */
    static File writeIdentityMapping(File directory, String name, int fullWidth, int height) {
        int halfWidth = fullWidth / 2;
        int iw = halfWidth / 4;
        int ih = height / 4;
        double[][] points = {
            {2 * iw, ih},
            {iw, 3 * ih},
            {3 * iw, 3 * ih},
        };

        Map<String, Object> mapping = new HashMap<>();
        mapping.put("source points", points);
        mapping.put("target points", points);
        mapping.put("image width", fullWidth);
        mapping.put("image height", height);

        File file = new File(directory, name);
        try {
            new ObjectMapper().writeValue(file, mapping);
        } catch (Exception e) {
            throw new IllegalStateException("could not write " + file, e);
        }
        return file;
    }
}
