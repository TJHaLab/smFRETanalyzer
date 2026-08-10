import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ij.ImagePlus;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.scijava.log.StderrLogService;

/**
 * Stages 2 and 3 end to end on a field whose answer is known.
 *
 * This is the test the repository most needed. A trace is the integrated intensity of a spot, and
 * for the whole life of the plugin it was exactly **half** of that: measureTimeTraces divided by
 * 2 pi sigma^2, the factor that turns a normalized Gaussian into a unit height one, where the
 * right factor is 4 pi sigma^2 because the spot's own width and the kernel's add in quadrature.
 * Nothing caught it, because a trace only looks wrong if you know what it should be - FRET
 * efficiencies were unaffected, since the factor cancels out of A/(D+A).
 *
 * So the field is generated with a known photon count per spot and the recovered trace is scored
 * against it. The expected answer is not 1: a camera integrates the PSF over its pixels and the
 * measurement reads the peak, which costs 1/(1 + 1/24 sigma^2) - 4% at spotSigma 1 and 0.5% at
 * spotSigma 3 - plus a little more from spot light leaking into the background estimate. What
 * must never come back is 0.48.
 *
 * Stage 1 is not exercised: it is the one part that needs TurboReg installed, and what it
 * produces is a mapping file, which SyntheticField writes directly as the identity.
 */
class TraceRecoveryTest {

    private static final int HALF_WIDTH = 160;
    private static final int HEIGHT = 160;
    private static final int FRAMES = 5;
    private static final double BACKGROUND = 20.0;
    private static final double DONOR_PHOTONS = 2000.0;

    // Half the donor, so a swapped channel is a factor of two rather than a subtlety.
    private static final double ACCEPTOR_FRACTION = 0.5;

    /** What pixel integration alone costs, which is most of the shortfall below 1. */
    private static double pixelIntegration(double sigma) {
        return 1.0 / (1.0 + 1.0 / (24.0 * sigma * sigma));
    }

    /**
     * Everything stage 3 needs, having actually run stage 2 to produce it.
     */
    private static final class Run {
        double[][] spots;           // The reloaded layout: x, y, snr, prominence.
        double[][] traces;          // [2 * spot][frame], donor at 2j and acceptor at 2j+1.
        List<SyntheticField.Spot> placed;
        String root;
    }

    /**
     * A sparse field of equal spots, run through spot finding and trace measurement.
     *
     * Sparse on purpose. Spot light that reaches past the mask lands in the background estimate
     * and is then subtracted from the traces, so a crowded field would fold a density dependent
     * bias into a number that is supposed to be measuring the norm.
     */
    private static Run measure(File directory, double sigma, int spacing, double donorPhotons) {
        StderrLogService log = new StderrLogService();

        List<SyntheticField.Spot> donor = new ArrayList<>();
        List<SyntheticField.Spot> acceptor = new ArrayList<>();
        for (int y = 30; y <= (HEIGHT - 31); y += spacing) {
            for (int x = 30; x <= (HALF_WIDTH - 31); x += spacing) {
                donor.add(new SyntheticField.Spot(x, y, donorPhotons));
                acceptor.add(new SyntheticField.Spot(x, y, ACCEPTOR_FRACTION * donorPhotons));
            }
        }

        ImagePlus movie = SyntheticField.movie(HALF_WIDTH, HEIGHT, donor, acceptor, sigma,
                BACKGROUND, FRAMES, 0.0, 1L);
        File movieFile = SyntheticField.writeMovie(directory, "sim.tif", movie);
        File mappingFile = SyntheticField.writeIdentityMapping(directory, "sim_mapping.json",
                2 * HALF_WIDTH, HEIGHT);

        smFRETSpotFinder stage2 = new smFRETSpotFinder();
        stage2.log = log;
        stage2.inputImageName = movieFile;
        stage2.mappingFile = mappingFile;
        stage2.startSlice = 1;
        stage2.endSlice = FRAMES;
        stage2.spotSigma = sigma;
        stage2.spotThreshold = 6.0;
        stage2.spotTolerance = 5.0;
        stage2.spotProminence = 0.4;
        stage2.cameraBlackLevel = 0;
        stage2.cameraGain = 1.0;
        stage2.spotSpacing = 3;
        stage2.edgeMargin = 5;
        stage2.backgroundKappa = 0.0;
        stage2.findSpots();

        String root = new File(directory, "sim").toString();

        smFRETAnalyzer stage3 = new smFRETAnalyzer();
        stage3.log = log;
        stage3.backgroundAverageNFrames = FRAMES;
        stage3.smfsf.log = log;
        stage3.smfsf.loadMappingJSON(mappingFile.toString());
        stage3.smfsf.loadMasks(root + "_spotf_masks.tif");

        Run result = new Run();
        result.placed = donor;
        result.root = root;
        result.spots = stage3.smfsf.loadSpotLocations(root + "_spotf_spots.csv");

        ImagePlus image = smFRETChannelMapper.toFloat(
                smFRETFiles.openImage(movieFile, "the image"));
        List<ImagePlus> background = stage3.backGroundEstimation(image);
        result.traces = stage3.measureTimeTraces(image, background, result.spots, sigma, 1.0);
        return result;
    }

    /** The donor trace of spot j, averaged over the movie. */
    private static double donorLevel(Run run, int spot) {
        double total = 0.0;
        for (double value : run.traces[2 * spot]) {
            total += value;
        }
        return total / run.traces[2 * spot].length;
    }

    private static double acceptorLevel(Run run, int spot) {
        double total = 0.0;
        for (double value : run.traces[2 * spot + 1]) {
            total += value;
        }
        return total / run.traces[2 * spot + 1].length;
    }

    private static double median(List<Double> values) {
        List<Double> sorted = new ArrayList<>(values);
        java.util.Collections.sort(sorted);
        int middle = sorted.size() / 2;
        return ((sorted.size() % 2) != 0) ? sorted.get(middle)
                : 0.5 * (sorted.get(middle - 1) + sorted.get(middle));
    }

    /**
     * The headline number: a spot carrying a known count of photoelectrons has to measure as that
     * count, to within the pixel integration the measurement cannot see.
     */
    @ParameterizedTest(name = "spotSigma={0}")
    @ValueSource(doubles = {1.0, 2.0, 3.0})
    @DisplayName("a spot of known integrated intensity measures as that intensity")
    void recoversKnownIntegratedIntensity(double sigma, @TempDir File directory) {
        Run run = measure(directory, sigma, 40, DONOR_PHOTONS);
        assertTrue(run.spots.length > 0, "spot finding produced nothing to measure");

        List<Double> recoveries = new ArrayList<>();
        for (int j = 0; j < run.spots.length; j++) {
            recoveries.add(donorLevel(run, j) / DONOR_PHOTONS);
        }
        double recovery = median(recoveries);

        // The specific regression, stated on its own because it is what this test is for: the old
        // norm recovered almost exactly half.
        assertTrue(recovery > 0.8,
                "recovered " + recovery + " of the true intensity at spotSigma " + sigma
                        + " - the pre-fix norm recovered about 0.48");

        double predicted = pixelIntegration(sigma);
        assertTrue(recovery < (predicted + 0.03),
                "recovered " + recovery + " at spotSigma " + sigma
                        + ", above the " + predicted + " pixel integration allows");
        assertTrue(recovery > (predicted - 0.08),
                "recovered " + recovery + " at spotSigma " + sigma + " against " + predicted
                        + " predicted - more spot light in the background estimate than expected");
    }

    /**
     * The link the whole viewer rests on: measureTimeTraces walks the spot table in load order,
     * so trace row j belongs to spot table row j. Nothing checks it at runtime beyond the lengths
     * matching, and if it ever slipped every trace would be attributed to the wrong molecule -
     * which would look like noisy data rather than like a bug.
     *
     * Giving each spot its own brightness is what makes that checkable: the trace has to match
     * the spot the table says it is, at the position the table gives.
     */
    @Test
    @DisplayName("trace row j belongs to spot table row j")
    void traceRowsFollowTheSpotTable(@TempDir File directory) {
        StderrLogService log = new StderrLogService();

        // A distinct brightness per spot, spread widely enough that a mix-up cannot be mistaken
        // for measurement error.
        List<SyntheticField.Spot> donor = new ArrayList<>();
        List<SyntheticField.Spot> acceptor = new ArrayList<>();
        double photons = 1200.0;
        for (int y = 30; y <= (HEIGHT - 31); y += 40) {
            for (int x = 30; x <= (HALF_WIDTH - 31); x += 40) {
                donor.add(new SyntheticField.Spot(x, y, photons));
                acceptor.add(new SyntheticField.Spot(x, y, ACCEPTOR_FRACTION * photons));
                photons += 900.0;
            }
        }

        ImagePlus movie = SyntheticField.movie(HALF_WIDTH, HEIGHT, donor, acceptor, 2.0,
                BACKGROUND, FRAMES, 0.0, 1L);
        File movieFile = SyntheticField.writeMovie(directory, "sim.tif", movie);
        File mappingFile = SyntheticField.writeIdentityMapping(directory, "sim_mapping.json",
                2 * HALF_WIDTH, HEIGHT);

        smFRETSpotFinder stage2 = new smFRETSpotFinder();
        stage2.log = log;
        stage2.inputImageName = movieFile;
        stage2.mappingFile = mappingFile;
        stage2.startSlice = 1;
        stage2.endSlice = FRAMES;
        stage2.spotSigma = 2.0;
        stage2.spotThreshold = 6.0;
        stage2.spotTolerance = 5.0;
        stage2.spotProminence = 0.4;
        stage2.cameraBlackLevel = 0;
        stage2.cameraGain = 1.0;
        stage2.spotSpacing = 3;
        stage2.edgeMargin = 5;
        stage2.backgroundKappa = 0.0;
        stage2.findSpots();

        String root = new File(directory, "sim").toString();
        smFRETAnalyzer stage3 = new smFRETAnalyzer();
        stage3.log = log;
        stage3.backgroundAverageNFrames = FRAMES;
        stage3.smfsf.log = log;
        stage3.smfsf.loadMappingJSON(mappingFile.toString());
        stage3.smfsf.loadMasks(root + "_spotf_masks.tif");

        Run run = new Run();
        run.spots = stage3.smfsf.loadSpotLocations(root + "_spotf_spots.csv");
        ImagePlus image = smFRETChannelMapper.toFloat(smFRETFiles.openImage(movieFile, "the image"));
        run.traces = stage3.measureTimeTraces(image,
                stage3.backGroundEstimation(image), run.spots, 2.0, 1.0);

        assertEquals(run.spots.length, run.traces.length / 2,
                "one pair of traces per spot in the table");
        assertTrue(run.spots.length >= 6, "only " + run.spots.length + " spots to check with");

        for (int j = 0; j < run.spots.length; j++) {
            double x = run.spots[j][0];
            double y = run.spots[j][1];

            // The spot the table row is standing at, found by position rather than by order.
            SyntheticField.Spot nearest = null;
            double best = Double.MAX_VALUE;
            for (SyntheticField.Spot candidate : donor) {
                double distance = Math.hypot(candidate.x - x, candidate.y - y);
                if (distance < best) {
                    best = distance;
                    nearest = candidate;
                }
            }
            assertTrue(best < 2.0, "table row " + j + " at " + x + "," + y
                    + " is not near any spot that was placed");

            double recovery = donorLevel(run, j) / nearest.intensity;
            assertTrue((recovery > 0.85) && (recovery < 1.05),
                    "row " + j + " at " + x + "," + y + " measured "
                            + donorLevel(run, j) + " where the spot there carries "
                            + nearest.intensity + " - traces and the spot table are out of step");
        }
    }

    /**
     * Donor is the left half and acceptor is the right, in that order, all the way through to the
     * trace matrix. A swap would leave every FRET efficiency reflected about 0.5 and nothing else
     * would complain.
     */
    @Test
    @DisplayName("the donor is the left half and the acceptor the right")
    void channelsAreNotSwapped(@TempDir File directory) {
        Run run = measure(directory, 2.0, 40, DONOR_PHOTONS);

        List<Double> ratios = new ArrayList<>();
        for (int j = 0; j < run.spots.length; j++) {
            ratios.add(acceptorLevel(run, j) / donorLevel(run, j));
        }

        assertEquals(ACCEPTOR_FRACTION, median(ratios), 0.03,
                "the acceptor should be half the donor, as generated");
    }

    /**
     * The whole of stage 3 through its own run(), which is the path a macro takes - so this is
     * also the only test that exercises the HDF5 and .traces writers and the spot finder JSON
     * contract between the two stages.
     */
    @Test
    @DisplayName("run() writes an h5 whose traces are the measured ones")
    void runWritesTheExpectedOutputs(@TempDir File directory) throws Exception {
        StderrLogService log = new StderrLogService();
        Run reference = measure(directory, 2.0, 40, DONOR_PHOTONS);

        smFRETAnalyzer stage3 = new smFRETAnalyzer();
        stage3.log = log;
        stage3.spotJSONFile = new File(reference.root + "_spotf_finding.json");
        stage3.backgroundAverageNFrames = FRAMES;
        stage3.run();

        File h5 = new File(reference.root + ".h5");
        File traces = new File(reference.root + ".traces");
        assertTrue(h5.exists(), "no .h5 at " + h5);
        assertTrue(traces.exists(), "no .traces at " + traces);

        try (IHDF5Reader reader = HDF5Factory.openForReading(h5)) {
            float[][] target = reader.readFloatMatrix("target-traces");
            float[][] source = reader.readFloatMatrix("source-traces");

            assertEquals(reference.spots.length, target.length, "one donor trace per spot");
            assertEquals(reference.spots.length, source.length, "one acceptor trace per spot");
            assertEquals(FRAMES, target[0].length, "one point per frame");

            for (int j = 0; j < target.length; j++) {
                for (int t = 0; t < FRAMES; t++) {
                    assertEquals(reference.traces[2 * j][t], target[j][t], 1.0e-2,
                            "donor spot " + j + " frame " + t);
                    assertEquals(reference.traces[2 * j + 1][t], source[j][t], 1.0e-2,
                            "acceptor spot " + j + " frame " + t);
                }
            }

            // The spot table travels with the traces, and in the reloaded layout - x at column 0.
            double[][] spots = reader.readDoubleMatrix("spots");
            assertEquals(reference.spots.length, spots.length);
            assertEquals(reference.spots[0][0], spots[0][0], 1.0e-9, "x");
            assertEquals(reference.spots[0][1], spots[0][1], 1.0e-9, "y");
        }
    }

    /**
     * Stage 2's four outputs, which are the contract with stage 3. Named individually because a
     * missing one of them is the failure that shows up several minutes into a batch run.
     */
    @Test
    @DisplayName("spot finding writes all four of its outputs")
    void spotFindingWritesItsOutputs(@TempDir File directory) {
        Run run = measure(directory, 2.0, 40, DONOR_PHOTONS);

        for (String suffix : new String[] {"_spotf_finding.json", "_spotf_spots.csv",
                                           "_spotf_masks.tif", "_spotf_qc_image.tif"}) {
            File file = new File(run.root + suffix);
            assertTrue(file.exists() && (file.length() > 0), "missing or empty: " + file);
        }

        // Every spot placed on a clean sparse field should survive every filter.
        assertEquals(run.placed.size(), run.spots.length,
                "placed " + run.placed.size() + " spots and recovered " + run.spots.length);
    }
}
