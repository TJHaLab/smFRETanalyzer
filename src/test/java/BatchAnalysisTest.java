import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import ij.ImagePlus;
import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.scijava.log.StderrLogService;

/**
 * Batch analysis (issue #21): copying one run's settings onto a folder of movies.
 *
 * Two things here are not obvious from the plugin's description and are the reason it has tests.
 *
 * The first is that the *mapping is not asked for*. A spot finder JSON records the mapping its
 * own run used, so an example run already names one - which makes a separate mapping input a
 * question whose answer is written in the file the user just chose. What has to work instead is
 * saying so clearly when that path has gone stale.
 *
 * The second is that an estimated spot tolerance is copied as **the decision to estimate**, not
 * as the number the example run arrived at. The JSON records the effective value, so copying it
 * literally would pin forty movies to one movie's noise floor - which is the exact thing issue
 * #10 introduced the estimate to stop.
 */
class BatchAnalysisTest {

    private static final int FRAMES = 6;
    private static final int HALF_WIDTH = 128;
    private static final int HEIGHT = 128;

    /** A two channel movie with a sparse field of equal spots. */
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

    /**
     * An example run: one movie through spot finding, leaving the JSON a batch copies.
     *
     * Returns the JSON. tolerance 0 means the run estimated it, which is the shipped default and
     * so the case a real template will usually be in.
     */
    private static File exampleRun(File directory, double tolerance) {
        File movie = writeMovie(directory, "example.tif", 5L);
        File mapping = SyntheticField.writeIdentityMapping(directory, "example_mapping.json",
                2 * HALF_WIDTH, HEIGHT);

        smFRETSpotFinder finder = new smFRETSpotFinder();
        finder.log = new StderrLogService();
        finder.showQCImage = false;
        finder.inputImageName = movie;
        finder.mappingFile = mapping;
        finder.startSlice = 1;
        finder.endSlice = FRAMES;
        finder.spotSigma = 2.0;
        finder.spotThreshold = 6.0;
        finder.spotTolerance = tolerance;
        finder.spotContamination = 1.0;
        finder.cameraBlackLevel = 0;
        finder.cameraGain = 1.0;
        finder.spotSpacing = 3;
        finder.edgeMargin = 5;
        finder.backgroundKappa = 0.0;
        finder.spotChannel = "donor";
        finder.findSpots();

        File json = new File(directory, "example_spotf_finding.json");
        assertTrue(json.isFile(), "the example run wrote no " + json);
        return json;
    }

    private static smFRETBatchAnalyzer batch() {
        smFRETBatchAnalyzer plugin = new smFRETBatchAnalyzer();
        plugin.log = new StderrLogService();
        return plugin;
    }

    /** Discards the progress reports, for the tests that only care about the outcome. */
    private static final smFRETBatchAnalyzer.Progress QUIET = (stages, what) -> { };

    /**
     * The settings come across, including the ones a batch would otherwise silently default.
     */
    @Test
    @DisplayName("an example run's settings are copied")
    void settingsAreCopied(@TempDir File directory) {
        File json = exampleRun(directory, 8.0);
        smFRETBatchAnalyzer.Settings settings = smFRETBatchAnalyzer.Settings.read(json);

        assertEquals(1, settings.startSlice);
        assertEquals(FRAMES, settings.endSlice);
        assertEquals(2.0, settings.spotSigma, 1.0e-9);
        assertEquals(6.0, settings.spotThreshold, 1.0e-9);
        assertEquals(1.0, settings.spotContamination, 1.0e-9);
        assertEquals(3, settings.spotSpacing);
        assertEquals(5, settings.edgeMargin);
        assertEquals("donor", settings.spotChannel, "the spot channel is a setting like any other");
        assertNull(settings.problem(), "the mapping it names is right there");
        assertEquals("example_mapping.json", settings.mappingFile.getName());
    }

    /**
     * A tolerance that was typed in is a choice, and is copied as it stands.
     */
    @Test
    @DisplayName("a set tolerance is copied as the number it was")
    void aSetToleranceIsCopied(@TempDir File directory) {
        smFRETBatchAnalyzer.Settings settings =
                smFRETBatchAnalyzer.Settings.read(exampleRun(directory, 8.0));

        assertFalse(settings.toleranceEstimated);
        assertEquals(8.0, settings.spotTolerance, 1.0e-9);
    }

    /**
     * A tolerance that was measured is copied as the decision to measure.
     *
     * The JSON records the *effective* value - the estimate, a number like 6.01 that the example
     * movie's noise floor produced - precisely so that the run can be explained. Copying that
     * number onto the batch would give every movie the example movie's noise floor, which defeats
     * the estimate. What is copied is the zero that asks for it to be measured again.
     */
    @Test
    @DisplayName("an estimated tolerance is copied as the decision to estimate")
    void anEstimatedToleranceIsRemeasured(@TempDir File directory) {
        File json = exampleRun(directory, 0.0);

        // The JSON should hold the estimate, not the zero - otherwise this test would pass for
        // the wrong reason and the settings would be right by accident.
        Object recorded = smFRETFiles.readSpotFinderJSON(json).get("spot tolerance");
        assertTrue(((Number) recorded).doubleValue() > 0.0,
                "the example JSON should record the value used, not the 0 that asked for it");

        smFRETBatchAnalyzer.Settings settings = smFRETBatchAnalyzer.Settings.read(json);
        assertTrue(settings.toleranceEstimated);
        assertEquals(0.0, settings.spotTolerance, 1.0e-9,
                "the batch should measure each movie, not inherit one movie's estimate");
        assertTrue(settings.describe().contains("estimated"), settings.describe());
    }

    /**
     * A JSON from before the current keys existed still runs, at the shipped defaults.
     *
     * Real templates will be old ones - the file on disk from the run someone did last month -
     * and every key this reads has been added at some point. An absent key has to leave its
     * default rather than fail the batch.
     */
    @Test
    @DisplayName("an old settings file falls back to the defaults")
    void anOldSettingsFileStillWorks(@TempDir File directory) throws Exception {
        File mapping = SyntheticField.writeIdentityMapping(directory, "old_mapping.json",
                2 * HALF_WIDTH, HEIGHT);

        // What a pre-1.2 run wrote: a prominence that no longer exists, no spot channel, no
        // contamination, and a tolerance with nothing saying whether it was measured.
        File json = new File(directory, "old_spotf_finding.json");
        Files.write(json.toPath(), ("{\"spot threshold\": 7.5, \"spot prominence\": 1.3,"
                + " \"spot sigma\": 1.4, \"camera gain\": 2.0, \"spot spacing\": 2,"
                + " \"spot margin\": 4, \"edge margin\": 9, \"camera black\": 100,"
                + " \"spot tolerance\": 10.0, \"spots file\": \"old_spotf_spots.csv\","
                + " \"mapping file\": " + quote(mapping) + "}").getBytes());

        smFRETBatchAnalyzer.Settings settings = smFRETBatchAnalyzer.Settings.read(json);
        smFRETSpotFinder defaults = new smFRETSpotFinder();

        assertEquals(7.5, settings.spotThreshold, 1.0e-9, "a key that is there is read");
        assertEquals(9, settings.edgeMargin);
        assertEquals(100, settings.cameraBlackLevel);
        assertEquals(10.0, settings.spotTolerance, 1.0e-9, "no 'estimated' key means it was set");
        assertFalse(settings.toleranceEstimated);

        assertEquals(defaults.spotChannel, settings.spotChannel, "absent, so the default");
        assertEquals(defaults.spotContamination, settings.spotContamination, 1.0e-9);
        assertEquals(defaults.startSlice, settings.startSlice);
        assertNull(settings.problem());
    }

    /**
     * A settings file whose mapping has moved says so, and names the file it wanted.
     *
     * This is the failure the missing mapping input trades for, so it is the one that has to be
     * good: the message has to name the path and say what to do, because the user has no field to
     * correct and would otherwise be stuck.
     */
    @Test
    @DisplayName("a moved mapping file is named, with what to do about it")
    void aMovedMappingIsNamed(@TempDir File directory) throws Exception {
        File json = exampleRun(directory, 5.0);
        File mapping = new File(directory, "example_mapping.json");
        Files.move(mapping.toPath(), new File(directory, "elsewhere.json").toPath());

        smFRETBatchAnalyzer.Settings settings = smFRETBatchAnalyzer.Settings.read(json);
        String problem = settings.problem();
        assertNotNull(problem, "a mapping that is not there should be a problem");
        assertTrue(problem.contains("example_mapping.json"), problem);
        assertTrue(problem.contains("smFRET Spot Finder"), "say what to do: " + problem);
    }

    /**
     * Anything that is not a spot finder JSON is refused as what it is.
     *
     * The mapping JSON is the one worth checking: it sits in the same folder, is named similarly,
     * and parses just as cleanly - so without the key check a batch would run on a file with none
     * of the settings in it and quietly use every default.
     */
    @Test
    @DisplayName("a mapping JSON is not a settings file")
    void aMappingJsonIsRefused(@TempDir File directory) {
        File mapping = SyntheticField.writeIdentityMapping(directory, "m.json",
                2 * HALF_WIDTH, HEIGHT);
        smFRETAnalysisException thrown = assertThrows(smFRETAnalysisException.class,
                () -> smFRETBatchAnalyzer.Settings.read(mapping));
        assertTrue(thrown.getMessage().contains("m.json"), thrown.getMessage());
    }

    /**
     * Both stages, over a movie that was not the one the settings came from.
     *
     * The whole plugin in one test: the outputs land beside the *batch* movie rather than beside
     * the example, and the stages that were told nothing about each other agree on where they are.
     */
    @Test
    @DisplayName("a batch movie is analysed with the example run's settings")
    void aMovieIsAnalysed(@TempDir File directory) {
        File json = exampleRun(directory, 5.0);
        File movie = writeMovie(directory, "batch01.tif", 11L);

        smFRETBatchAnalyzer plugin = batch();
        List<String> stages = new ArrayList<>();
        String problem = plugin.analyze(movie, smFRETBatchAnalyzer.Settings.read(json), FRAMES,
                (done, what) -> stages.add(what));

        assertNull(problem, "the movie should have analysed: " + problem);
        assertEquals(java.util.Arrays.asList("finding spots", "measuring traces", "done"), stages);

        File batchJson = new File(directory, "batch01_spotf_finding.json");
        File batchH5 = new File(directory, "batch01.h5");
        assertTrue(batchJson.isFile(), "no settings JSON beside the batch movie");
        assertTrue(batchH5.isFile(), "no traces beside the batch movie");
        assertTrue(new File(directory, "batch01_analysis").isDirectory(), "no analysis folder");

        // And the run really used the settings rather than the defaults.
        java.util.Map<String, Object> written = smFRETFiles.readSpotFinderJSON(batchJson);
        assertEquals("donor", written.get("spot channel"));
        assertEquals(movie.toString(), written.get("image name"));
        assertEquals(new File(directory, "example_mapping.json").toString(),
                written.get("mapping file"), "the mapping named by the example run");
    }

    /**
     * A movie that cannot be analysed is reported rather than passing silently.
     *
     * Both stages catch their own exceptions and log them, which is right for a menu item and
     * leaves a batch with nothing to test - so the outcome is judged by whether the output
     * changed. Without that, a failure looks exactly like a success.
     */
    @Test
    @DisplayName("a movie that fails is reported")
    void aFailedMovieIsReported(@TempDir File directory) throws Exception {
        File json = exampleRun(directory, 5.0);

        File notAMovie = new File(directory, "broken.tif");
        Files.write(notAMovie.toPath(), "this is not a movie".getBytes());

        smFRETBatchAnalyzer plugin = batch();
        String problem = plugin.analyze(notAMovie, smFRETBatchAnalyzer.Settings.read(json),
                FRAMES, QUIET);

        assertNotNull(problem, "a movie that cannot be read should not report success");
        assertTrue(problem.startsWith("spot finding failed"),
                "the verdict should come first, so a truncated row still says what happened: "
                        + problem);
        assertTrue(problem.contains("broken_spotf_finding.json"), problem);
    }

    /**
     * A stale output from an earlier run is not mistaken for this one's.
     *
     * The trap in judging success by the output: a movie analysed last week already has all of
     * its files, so "the file is there" would pass for every movie in a batch that failed on
     * every movie.
     */
    @Test
    @DisplayName("an earlier run's output does not count as this run's")
    void aStaleOutputIsNotSuccess(@TempDir File directory) throws Exception {
        File json = exampleRun(directory, 5.0);
        File movie = writeMovie(directory, "batch01.tif", 11L);

        smFRETBatchAnalyzer plugin = batch();
        assertNull(plugin.analyze(movie, smFRETBatchAnalyzer.Settings.read(json), FRAMES, QUIET));

        // Now break the movie, leaving every output of the good run in place, and re-run. The
        // files are all still there and all still valid; none of them is from this run.
        Files.write(movie.toPath(), "this is not a movie any more".getBytes());
        String problem = plugin.analyze(movie, smFRETBatchAnalyzer.Settings.read(json),
                FRAMES, QUIET);

        assertNotNull(problem, "the previous run's files should not count as this run's");
    }

    /**
     * The queue takes movies and leaves the files that sit beside them alone.
     */
    @Test
    @DisplayName("only movies are queued")
    void onlyMoviesAreQueued(@TempDir File directory) {
        File json = exampleRun(directory, 5.0);
        File movie = writeMovie(directory, "batch01.tif", 11L);

        assertTrue(smFRETFiles.isMovie(movie), "a movie is a movie");
        assertFalse(smFRETFiles.isMovie(json), "a settings JSON is not");
        assertFalse(smFRETFiles.isMovie(new File(directory, "example_mapping.json")),
                "nor is a mapping JSON");
        assertFalse(smFRETFiles.isMovie(new File(directory, "nothing.tif")),
                "nor is a file that is not there");
    }

    /**
     * The settings line wraps rather than clipping, and a long path can break.
     *
     * Both are one property really: a JLabel given plain text draws a single line and truncates
     * it to an ellipsis. That is survivable for a status line and not for this one, because the
     * message that matters most here is the one reporting a *path* - and the path sits in the
     * middle of the sentence, so the half a narrow window drops is the half naming the file the
     * user has to go and find. A settings file made on another machine records a path that will
     * never resolve here, so this is the ordinary case rather than a corner of it.
     *
     * The break opportunities are the other half: a path has no spaces, so with nowhere legal to
     * break it the label clips it *and* takes the whole path as its minimum width, which puts the
     * message off its own right edge instead of wrapping it.
     */
    @Test
    @DisplayName("the settings line wraps, and a long path can break")
    void theSettingsLineWraps() {
        assertTrue(smFRETBatchAnalyzer.wrapped("x").startsWith("<html>"),
                "plain text in a JLabel is one clipped line, not a wrapped paragraph");

        String windows = smFRETBatchAnalyzer.breakable(
                "C:\\Users\\xiaot\\Desktop\\test\\film1_mapping.json");
        assertTrue(windows.contains("\u200b"), "a Windows path should gain break opportunities");
        assertEquals("C:\\Users\\xiaot\\Desktop\\test\\film1_mapping.json",
                windows.replace("\u200b", ""),
                "and nothing else - the path still has to read as itself");

        assertTrue(smFRETBatchAnalyzer.breakable("/home/x/y.json").contains("\u200b"),
                "a POSIX path too");

        // Prose is left alone, so a message without a path in it is not peppered with them.
        String prose = "re-run smFRET Spot Finder on one movie";
        assertEquals(prose, smFRETBatchAnalyzer.breakable(prose));
    }

    /**
     * Names going into those labels are escaped.
     *
     * File names are user data. One containing an ampersand or an angle bracket would otherwise be
     * eaten by the HTML parser - so the message reporting a file would mangle exactly the names
     * most likely to have been chosen badly in the first place.
     */
    @Test
    @DisplayName("a file name with markup in it survives the label")
    void namesAreEscaped() {
        assertEquals("a&amp;b &lt;c&gt;.json", smFRETBatchAnalyzer.escaped("a&b <c>.json"));

        // The ampersand is replaced first, or the escapes escape each other.
        assertEquals("&amp;lt;", smFRETBatchAnalyzer.escaped("&lt;"));
    }

    private static String quote(File file) {
        return "\"" + file.toString().replace("\\", "\\\\") + "\"";
    }
}
