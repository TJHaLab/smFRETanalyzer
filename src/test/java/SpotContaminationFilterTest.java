import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import ij.ImagePlus;
import ij.measure.ResultsTable;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.scijava.log.StderrLogService;

/**
 * The contamination filter, driven through the spot finder rather than called directly.
 *
 * <p>{@code SpotQualityForestTest} proves the model evaluates like its trainer and
 * {@code SpotQualityFeatureTest} proves it is fed the right numbers. Neither says the plugin
 * is wired to them: that the score reaches the spot table, that the threshold actually drops
 * spots, or that the column contract survived a filter being added. That is what this covers.
 *
 * <p><b>The scores on this field are not meaningful as contamination.</b> {@code
 * SyntheticField} paints Gaussian spots and the model was trained on aberrated Airy patterns,
 * whose wings are the thing it mostly reads, so everything here scores between 0.109 and 0.461
 * with a median of 0.112 - a floor far above the near zero a clean spot earns on real data.
 * The thresholds below are picked from that measured range. What is being tested is the
 * plumbing; whether the number is right is settled by the fixture tests and by SIMULATION.md.
 */
class SpotContaminationFilterTest {

    private static final int HALF_WIDTH = 128;
    private static final int HEIGHT = 128;
    private static final int FRAMES = 12;
    private static final double SIGMA = 1.5;

    /** A field with well separated spots and, deliberately, some very close pairs. */
    private static List<SyntheticField.Spot> field() {
        List<SyntheticField.Spot> spots = new ArrayList<>();
        for (int y = 20; y <= (HEIGHT - 20); y += 20) {
            for (int x = 20; x <= (HALF_WIDTH - 20); x += 20) {
                spots.add(new SyntheticField.Spot(x, y, 900.0));
                // Every third site gets a partner three pixels away - far enough that both
                // are found, close enough that each is badly contaminated by the other.
                if (((x + y) % 60) == 0) {
                    spots.add(new SyntheticField.Spot(x + 3.0, y + 1.0, 800.0));
                }
            }
        }
        return spots;
    }

    private static smFRETSpotFinder finder(File directory, double contamination) {
        List<SyntheticField.Spot> spots = field();
        // Real shot noise rather than a noiseless field: the model reads a radial maximum
        // per annulus, which on a perfectly smooth image is a quantity it never saw.
        ImagePlus movie = SyntheticField.movie(HALF_WIDTH, HEIGHT, spots, spots, SIGMA,
                20.0, FRAMES, Math.sqrt(20.0), 7L);
        File movieFile = SyntheticField.writeMovie(directory, "field.tif", movie);
        File mappingFile = SyntheticField.writeIdentityMapping(directory, "field_mapping.json",
                2 * HALF_WIDTH, HEIGHT);

        smFRETSpotFinder finder = new smFRETSpotFinder();
        finder.log = new StderrLogService();
        finder.inputImageName = movieFile;
        finder.mappingFile = mappingFile;
        finder.startSlice = 1;
        finder.endSlice = FRAMES;
        finder.spotSigma = SIGMA;
        finder.spotThreshold = 6.0;
        finder.spotTolerance = 5.0;
        finder.spotContamination = contamination;
        finder.cameraBlackLevel = 0;
        finder.cameraGain = 1.0;
        finder.spotSpacing = 3;
        finder.edgeMargin = 5;
        finder.backgroundKappa = 0.0;
        return finder;
    }

    private static ResultsTable table(File directory) {
        return ResultsTable.open2(new File(directory, "field_spotf_spots.csv").toString());
    }

    @Test
    @DisplayName("the score reaches the spot table as a column of fractions")
    void theScoreIsWritten(@TempDir File directory) {
        finder(directory, 1.0).findSpots();
        ResultsTable rt = table(directory);

        assertTrue(rt.getCounter() > 0, "no spots found");
        assertTrue(rt.getColumnIndex("contamination") != ResultsTable.COLUMN_NOT_FOUND,
                "no contamination column in the spot table");

        double largest = 0.0;
        for (int i = 0; i < rt.getCounter(); i++) {
            double value = rt.getValue("contamination", i);
            assertTrue((value >= 0.0) && (value <= 1.0), "not a fraction: " + value);
            largest = Math.max(largest, value);
        }
        // A field built with close pairs in it has to produce some contaminated spots, or
        // the score is being computed on something other than this image.
        // The close pairs have to move the score, or it is being computed on something
        // other than this image. The spread is what matters, not the level.
        assertTrue(largest > 0.3, "nothing scored above 0.3, largest was " + largest);
    }

    /**
     * The score is a diagnostic column, so switching the filter off must still record it.
     * An earlier version returned early on a threshold of 1 and wrote a column of zeros.
     */
    @Test
    @DisplayName("the score is recorded even with the filter switched off")
    void theScoreSurvivesTheFilterBeingOff(@TempDir File directory) {
        finder(directory, 1.0).findSpots();
        ResultsTable rt = table(directory);

        double largest = 0.0;
        for (int i = 0; i < rt.getCounter(); i++) {
            largest = Math.max(largest, rt.getValue("contamination", i));
        }
        assertTrue(largest > 0.0, "the column is all zeros with the filter off");
    }

    /**
     * The threshold has to do something, and it has to do it in the right direction.
     */
    @Test
    @DisplayName("a stricter threshold keeps fewer spots")
    void theThresholdBites(@TempDir File directory) {
        finder(directory, 1.0).findSpots();
        int unfiltered = table(directory).getCounter();

        finder(directory, 0.12).findSpots();   // just above this field's median
        int strict = table(directory).getCounter();

        assertTrue(strict < unfiltered,
                "strict kept " + strict + " of " + unfiltered + ", so the threshold did nothing");
        assertTrue(strict > 0, "the strict threshold emptied the field");
    }

    /**
     * What survives a threshold must be what scored below it. Checking the surviving scores
     * rather than only the count is what catches a filter that drops the right *number* of
     * spots for the wrong reason.
     */
    @Test
    @DisplayName("every surviving spot scored below the threshold")
    void survivorsAreUnderTheThreshold(@TempDir File directory) {
        double threshold = 0.12;
        finder(directory, threshold).findSpots();
        ResultsTable rt = table(directory);

        assertTrue(rt.getCounter() > 0, "no spots survived");
        for (int i = 0; i < rt.getCounter(); i++) {
            assertTrue(rt.getValue("contamination", i) <= threshold,
                    "kept a spot scoring " + rt.getValue("contamination", i));
        }
    }

    /**
     * Setting it to 1 switches the filter off rather than rejecting everything, which is what
     * the documentation promises and the opposite of what a naive {@code >=} would do.
     */
    @Test
    @DisplayName("a threshold of 1 keeps everything")
    void oneKeepsEverything(@TempDir File directory) {
        finder(directory, 1.0).findSpots();
        int kept = table(directory).getCounter();

        finder(directory, 0.99).findSpots();
        assertEquals(kept, table(directory).getCounter(),
                "1.0 and 0.99 should keep the same spots, since nothing here scores that high");
    }

    /**
     * The round trip: what the finder wrote is what the loader reads back, in the reloaded
     * layout with the flag gone and contamination last.
     */
    @Test
    @DisplayName("the contamination column survives the save and load round trip")
    void theColumnRoundTrips(@TempDir File directory) {
        smFRETSpotFinder finder = finder(directory, 1.0);
        finder.findSpots();

        double[][] reloaded = finder.loadSpotLocations(
                new File(directory, "field_spotf_spots.csv").toString());
        ResultsTable rt = table(directory);

        assertEquals(rt.getCounter(), reloaded.length, "row count");
        assertEquals(finder.columnHeaders.size(), reloaded[0].length, "column count");
        for (int i = 0; i < reloaded.length; i++) {
            assertEquals(rt.getValue("contamination", i), reloaded[i][4], 1.0e-9,
                    "contamination in the last column, row " + i);
            assertEquals(rt.getValue("x", i), reloaded[i][0], 1.0e-9, "x in column 0");
        }
    }
}
