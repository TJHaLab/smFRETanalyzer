import static org.junit.jupiter.api.Assertions.assertTrue;

import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.scijava.log.StderrLogService;

/**
 * The predicted best spotSigma, checked against what the pipeline actually recovers.
 *
 * filterResponse is a continuous integral. The real thing is not: the camera bins the PSF into
 * pixels, MaximumFinder reports the integer pixel rather than the centroid so the filter is
 * applied off centre by up to half a pixel, and the background subtracted is an estimate rather
 * than the truth. Each of those could move the optimum, and all of them in the same direction -
 * a wider filter is more forgiving of a centring error and less sensitive to a background
 * estimate's own structure. So the prediction is worth nothing until it has been run.
 *
 * This measures trace SNR the way somebody looking at a trace would: a field of identical
 * non-bleaching molecules, so a trace should be a flat line, and its SNR is its mean over its
 * standard deviation. Spots are placed at **random sub-pixel positions**, which is what puts the
 * centring error into the answer.
 */
class PSFFilterSweepTest {

    private static final double TRUE_SIGMA = 1.5;
    private static final double TRUE_WAVES = 0.4;

    private static final int HALF_WIDTH = 160;
    private static final int HEIGHT = 160;
    // A per trace standard deviation over this many frames is uncertain by about 13%, and the
    // median over the spots on the field brings that to a few percent - which is what it takes to
    // see a curve whose whole span is 10%. Fewer frames was tried at 22 and the peak became
    // noise: the settings at 1.8, 2.0 and 2.4 came out within 1% of each other.
    private static final int FRAMES = 30;
    private static final double BACKGROUND = 20.0;
    private static final double AMPLITUDE = 900.0;

    // Sub-samples per pixel per axis when rendering. The camera integrates the PSF over each
    // pixel and 3x3 midpoints is enough to stand in for that - the difference against 5x5 is far
    // below the noise this field carries.
    private static final int SUPERSAMPLE = 3;

    private static final double RENDER_REACH = 25.0;
    private static final double GRID_STEP = 0.02;
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
        return GRID[index] + ((position - index) * (GRID[index + 1] - GRID[index]));
    }

    /** One pixel's worth of the PSF, integrated rather than sampled at its centre. */
    private static double pixel(double dx, double dy) {
        double total = 0.0;
        for (int sy = 0; sy < SUPERSAMPLE; sy++) {
            for (int sx = 0; sx < SUPERSAMPLE; sx++) {
                double ox = ((sx + 0.5) / SUPERSAMPLE) - 0.5;
                double oy = ((sy + 0.5) / SUPERSAMPLE) - 0.5;
                total += shapeAt(Math.hypot(dx + ox, dy + oy));
            }
        }
        return total / (SUPERSAMPLE * SUPERSAMPLE);
    }

    /** A two channel movie of identical non-bleaching molecules on a noisy background. */
    private static ImagePlus movie(List<double[]> spots, long seed) {
        Random random = new Random(seed);
        int width = 2 * HALF_WIDTH;
        ImageStack stack = new ImageStack(width, HEIGHT);

        float[] clean = new float[width * HEIGHT];
        java.util.Arrays.fill(clean, (float) BACKGROUND);
        int reach = (int) Math.ceil(RENDER_REACH);
        for (double[] spot : spots) {
            for (int half = 0; half < 2; half++) {
                double cx = spot[0] + (half * HALF_WIDTH);
                double cy = spot[1];
                for (int y = Math.max(0, (int) cy - reach);
                        y <= Math.min(HEIGHT - 1, (int) cy + reach); y++) {
                    for (int x = Math.max(0, (int) cx - reach);
                            x <= Math.min(width - 1, (int) cx + reach); x++) {
                        clean[y * width + x] += (float) (AMPLITUDE * pixel(x - cx, y - cy));
                    }
                }
            }
        }

        // Background dominated noise, which is the regime filterResponse assumes and the one a
        // dim single molecule is actually in.
        double noise = Math.sqrt(BACKGROUND);
        for (int frame = 0; frame < FRAMES; frame++) {
            float[] pixels = clean.clone();
            for (int i = 0; i < pixels.length; i++) {
                pixels[i] += (float) (noise * random.nextGaussian());
            }
            stack.addSlice(new FloatProcessor(width, HEIGHT, pixels, null));
        }
        return new ImagePlus("sweep", stack);
    }

    /** Well separated spots at random sub-pixel positions. */
    private static List<double[]> spots(long seed) {
        Random random = new Random(seed);
        List<double[]> spots = new ArrayList<>();
        for (int y = 30; y <= (HEIGHT - 31); y += 26) {
            for (int x = 30; x <= (HALF_WIDTH - 31); x += 26) {
                spots.add(new double[] {x + random.nextDouble() - 0.5,
                        y + random.nextDouble() - 0.5});
            }
        }
        return spots;
    }

    /** Median over spots of a trace's mean over its standard deviation. */
    private static double medianSnr(double[][] traces, int channel) {
        List<Double> snr = new ArrayList<>();
        for (int spot = 0; spot < (traces.length / 2); spot++) {
            double[] trace = traces[(2 * spot) + channel];

            double mean = 0.0;
            for (double v : trace) {
                mean += v;
            }
            mean /= trace.length;

            double sum = 0.0;
            for (double v : trace) {
                sum += (v - mean) * (v - mean);
            }
            double deviation = Math.sqrt(sum / (trace.length - 1));
            if (deviation > 0.0) {
                snr.add(mean / deviation);
            }
        }
        java.util.Collections.sort(snr);
        return snr.get(snr.size() / 2);
    }

    @Test
    @DisplayName("the predicted best spotSigma is the one the pipeline actually prefers")
    void predictionMatchesTheSweep(@TempDir File directory) {
        StderrLogService log = new StderrLogService();

        List<double[]> placed = spots(11L);
        ImagePlus movie = movie(placed, 22L);
        File movieFile = SyntheticField.writeMovie(directory, "sweep.tif", movie);
        File mappingFile = SyntheticField.writeIdentityMapping(directory, "sweep_mapping.json",
                2 * HALF_WIDTH, HEIGHT);

        // Stage 2 once, at the fitted core width. The spot set and the masks are held fixed
        // across the sweep on purpose: what is being isolated is the filter, and re-finding spots
        // at every setting would change which molecules are being measured as well.
        smFRETSpotFinder finder = new smFRETSpotFinder();
        finder.log = log;
        finder.inputImageName = movieFile;
        finder.mappingFile = mappingFile;
        finder.startSlice = 1;
        finder.endSlice = FRAMES;
        finder.spotSigma = TRUE_SIGMA;
        finder.spotThreshold = 6.0;
        finder.spotTolerance = 5.0;
        finder.spotProminence = 0.4;
        finder.cameraBlackLevel = 0;
        finder.cameraGain = 1.0;
        finder.spotSpacing = 3;
        finder.edgeMargin = 5;
        finder.backgroundKappa = 0.0;
        finder.findSpots();

        String root = new File(directory, "sweep").toString();
        smFRETAnalyzer analyzer = new smFRETAnalyzer();
        analyzer.log = log;
        analyzer.backgroundAverageNFrames = FRAMES;
        analyzer.smfsf.log = log;
        analyzer.smfsf.loadMappingJSON(mappingFile.toString());
        analyzer.smfsf.loadMasks(root + "_spotf_masks.tif");
        double[][] found = analyzer.smfsf.loadSpotLocations(root + "_spotf_spots.csv");
        assertTrue(found.length > (0.8 * placed.size()),
                "only " + found.length + " of " + placed.size() + " spots were found");

        ImagePlus image = smFRETChannelMapper.toFloat(
                smFRETFiles.openImage(movieFile, "the image"));

        // The background estimate does not depend on spotSigma, so it is computed once - which is
        // also most of the cost of a stage 3 run.
        List<ImagePlus> background = analyzer.backGroundEstimation(image);

        // Fine where the peak is, coarse at the ends where the only job is to show the curve
        // really does fall away.
        double[] settings = {1.0, 1.4, 1.6, 1.8, 2.0, 2.4, 3.0};
        double[] measured = new double[settings.length];
        int bestAt = 0;
        for (int i = 0; i < settings.length; i++) {
            double[][] traces = analyzer.measureTimeTraces(image, background, found,
                    settings[i], 1.0);
            measured[i] = medianSnr(traces, 0);
            if (measured[i] > measured[bestAt]) {
                bestAt = i;
            }
            log.info(String.format("spotSigma %.2f -> median trace SNR %.2f",
                    settings[i], measured[i]));
        }

        smFRETPSF.FilterOptimum predicted = smFRETPSF.optimalFilter(
                new double[][] {smFRETPSF.airyProfile(smFRETPSF.SNR_RADII, TRUE_SIGMA, TRUE_WAVES)},
                0.02);
        log.info(String.format("predicted %.2f, band %.2f to %.2f; swept best %.2f",
                predicted.best, predicted.low, predicted.high, settings[bestAt]));

        // The sweep's own maximum has to land inside the predicted band. Comparing the two
        // maxima directly would be asking too much of a curve this flat: the SNR at neighbouring
        // settings differs by well under the scatter of a median over a few dozen traces.
        assertTrue((settings[bestAt] >= (predicted.low - 0.2))
                        && (settings[bestAt] <= (predicted.high + 0.2)),
                "the sweep peaked at spotSigma " + settings[bestAt]
                        + ", outside the predicted band " + predicted.low + " to "
                        + predicted.high);

        // And the shape agrees where it matters: the prediction says the ends of the range are
        // clearly worse, so the measurement has to see that too. If it does not, the curve is
        // being swamped by noise and the agreement above means nothing.
        assertTrue(measured[0] < (0.97 * measured[bestAt]),
                "spotSigma 1.0 measured " + measured[0] + " against a best of "
                        + measured[bestAt] + " - too flat for this test to be measuring anything");
        assertTrue(measured[settings.length - 1] < (0.97 * measured[bestAt]),
                "spotSigma 3.0 measured " + measured[settings.length - 1]
                        + " against a best of " + measured[bestAt]);
    }
}
