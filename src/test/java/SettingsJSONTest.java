import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertInstanceOf;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import ij.ImagePlus;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.scijava.log.StderrLogService;

/**
 * What a run records about itself, when its File objects are not exactly java.io.File.
 *
 * <p>A file chooser does not promise to hand back a `java.io.File`. On Windows JFileChooser
 * returns `sun.awt.shell.Win32ShellFolder2`, a live handle on a shell item that extends File.
 * That distinction is invisible everywhere in this pipeline except one place: Jackson registers
 * its File serializer against `java.io.File` **exactly**, so a subclass misses the lookup and is
 * written as a bean - an object of paths and disk sizes where every reader expects a string.
 *
 * <p>Under Java 16 and later it is worse than a wrong shape. `sun.awt.shell` is no longer open
 * for reflection, so Jackson finds nothing it can read, throws part way through the file, and
 * leaves a truncated JSON. That file is newer than it was, which is all the batch's "did spot
 * finding write its output" check looks at, so the movie was reported as having got past a stage
 * it died in. Fiji's Java 21 is what turned a latent bug into a broken batch.
 *
 * <p>The subclass below is a plain one rather than a shell folder, because no portable class
 * reproduces the reflection failure and the tests have to run everywhere. It does reproduce the
 * cause - dispatch on the exact class - and so fails on the unfixed writer on every platform.
 */
class SettingsJSONTest {

    private static final int FRAMES = 6;
    private static final int HALF_WIDTH = 128;
    private static final int HEIGHT = 128;

    /** A File that is not exactly java.io.File, standing in for what a chooser returns. */
    private static final class ChooserFile extends File {
        ChooserFile(File file) {
            super(file.getPath());
        }
    }

    private static File writeMovie(File directory, String name, long seed) {
        List<SyntheticField.Spot> donor = new ArrayList<>();
        List<SyntheticField.Spot> acceptor = new ArrayList<>();
        for (int y = 30; y <= (HEIGHT - 31); y += 30) {
            for (int x = 30; x <= (HALF_WIDTH - 31); x += 30) {
                donor.add(new SyntheticField.Spot(x, y, 3000.0));
                acceptor.add(new SyntheticField.Spot(x, y, 1500.0));
            }
        }
        ImagePlus movie = SyntheticField.movie(HALF_WIDTH, HEIGHT, donor, acceptor, 2.0,
                20.0, FRAMES, 0.0, seed);
        return SyntheticField.writeMovie(directory, name, movie);
    }

    private static smFRETSpotFinder finder(File movie, File mapping) {
        smFRETSpotFinder finder = new smFRETSpotFinder();
        finder.log = new StderrLogService();
        finder.showQCImage = false;
        finder.inputImageName = movie;
        finder.mappingFile = mapping;
        finder.startSlice = 1;
        finder.endSlice = FRAMES;
        finder.spotSigma = 2.0;
        finder.spotThreshold = 6.0;
        finder.spotTolerance = 5.0;
        finder.spotContamination = 1.0;
        finder.cameraBlackLevel = 0;
        finder.cameraGain = 1.0;
        finder.spotSpacing = 3;
        finder.edgeMargin = 5;
        finder.backgroundKappa = 0.0;
        finder.spotChannel = "donor";
        return finder;
    }

    /**
     * Both paths come back as strings, whatever kind of File went in.
     *
     * Asserted as String and not merely as equal, because the failure this guards writes a
     * *different shape* rather than a different value, and every reader of these two keys casts
     * them to String.
     */
    @Test
    @DisplayName("paths are recorded as strings, whatever kind of File was chosen")
    void pathsAreWrittenAsStrings(@TempDir File directory) {
        File movie = writeMovie(directory, "chooser.tif", 5L);
        File mapping = SyntheticField.writeIdentityMapping(directory, "chooser_mapping.json",
                2 * HALF_WIDTH, HEIGHT);

        finder(new ChooserFile(movie), new ChooserFile(mapping)).findSpots();

        File json = new File(directory, "chooser_spotf_finding.json");
        assertTrue(json.isFile(), "the run wrote no " + json);

        Map<String, Object> written = smFRETFiles.readSpotFinderJSON(json);
        assertInstanceOf(String.class, written.get("image name"), "image name");
        assertInstanceOf(String.class, written.get("mapping file"), "mapping file");
        assertEquals(movie.toString(), written.get("image name"));
        assertEquals(mapping.toString(), written.get("mapping file"));
    }

    /**
     * The whole batch stage on a movie added through the file chooser rather than dropped.
     *
     * The two ways of adding movies hand over different kinds of File - drag and drop gives plain
     * ones, the chooser does not - and only one of them was ever exercised, which is how this
     * shipped. Trace measurement reads back what spot finding wrote, so it fails here too, and on
     * every Java version rather than only the one that throws.
     */
    @Test
    @DisplayName("a movie added through the file chooser analyses end to end")
    void aChooserMovieIsAnalysed(@TempDir File directory) {
        File example = writeMovie(directory, "example.tif", 5L);
        File mapping = SyntheticField.writeIdentityMapping(directory, "example_mapping.json",
                2 * HALF_WIDTH, HEIGHT);
        finder(example, mapping).findSpots();
        File exampleJSON = new File(directory, "example_spotf_finding.json");
        assertTrue(exampleJSON.isFile(), "the example run wrote no " + exampleJSON);

        File movie = writeMovie(directory, "batch01.tif", 11L);
        smFRETBatchAnalyzer plugin = new smFRETBatchAnalyzer();
        plugin.log = new StderrLogService();
        String problem = plugin.analyze(new ChooserFile(movie),
                smFRETBatchAnalyzer.Settings.read(exampleJSON), FRAMES, (done, what) -> { });

        assertNull(problem, "the movie should have analysed: " + problem);
        assertTrue(new File(directory, "batch01.h5").isFile(), "no traces beside the movie");
    }
}
