import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import com.fasterxml.jackson.databind.ObjectMapper;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.process.ByteProcessor;
import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.Map;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

/**
 * Saying what the file actually is.
 *
 * The Browse dialog cannot filter - SciJava's single file chooser takes no filter, and a typed
 * path or a macro bypasses the widget anyway - so the check has to run when the plugin does. The
 * two JSON files are the reason any of this exists: they sit beside each other with the same
 * extension and are told apart only by their contents, and handing the spot finder's own output
 * back to it as a mapping used to parse cleanly and fail several steps later as "Cannot transform
 * image, transform model not set", which names neither the file nor the mistake.
 *
 * So these assert on the *messages*. A check that merely throws would satisfy a test that only
 * looked for an exception, and be no better than what it replaced.
 */
class FileSniffTest {

    private static File mappingJSON(File directory) {
        return SyntheticField.writeIdentityMapping(directory, "movie_mapping.json", 128, 64);
    }

    private static File spotFinderJSON(File directory) throws Exception {
        Map<String, Object> json = new HashMap<>();
        json.put("spots file", "movie_spotf_spots.csv");
        json.put("masks file", "movie_spotf_masks.tif");
        json.put("image name", "movie.tif");
        json.put("mapping file", "movie_mapping.json");
        json.put("root name", "movie");
        json.put("spot sigma", 2.0);
        json.put("camera black", 0);
        json.put("camera gain", 1.0);

        File file = new File(directory, "movie_spotf_finding.json");
        new ObjectMapper().writeValue(file, json);
        return file;
    }

    private static File tiff(File directory) {
        ImagePlus image = new ImagePlus("movie", new ByteProcessor(8, 8));
        File file = new File(directory, "movie.tif");
        new FileSaver(image).saveAsTiff(file.toString());
        return file;
    }

    private static File hdf5(File directory) throws Exception {

        // Only the signature is looked at, so a real HDF5 file is not needed to test the sniff.
        byte[] signature = {(byte) 0x89, 'H', 'D', 'F', '\r', '\n', 0x1a, '\n', 0, 0, 0, 0};
        File file = new File(directory, "movie.h5");
        Files.write(file.toPath(), signature);
        return file;
    }

    private static File spotTable(File directory) throws Exception {
        File file = new File(directory, "movie_spotf_spots.csv");
        Files.write(file.toPath(),
                " ,x,y,snr,prominence\n1,10,20,8.5,0.94\n".getBytes(StandardCharsets.UTF_8));
        return file;
    }

    /** The mistake the whole class exists for, in one direction. */
    @Test
    @DisplayName("a spot finder JSON offered as a mapping names both")
    void spotFinderJsonIsNotAMapping(@TempDir File directory) throws Exception {
        File finding = spotFinderJSON(directory);

        smFRETAnalysisException thrown = assertThrows(smFRETAnalysisException.class,
                () -> smFRETFiles.readMappingJSON(finding.toString()));

        String message = thrown.getMessage();
        assertTrue(message.contains("is a spot finder JSON file"), message);
        assertTrue(message.contains("not a channel mapping JSON file"), message);
        assertTrue(message.contains("_mapping.json"), "the message should say which file to pick: " + message);
    }

    /** And the other. */
    @Test
    @DisplayName("a mapping JSON offered as a spot finder file names both")
    void mappingJsonIsNotASpotFinderFile(@TempDir File directory) {
        File mapping = mappingJSON(directory);

        smFRETAnalysisException thrown = assertThrows(smFRETAnalysisException.class,
                () -> smFRETFiles.readSpotFinderJSON(mapping));

        String message = thrown.getMessage();
        assertTrue(message.contains("is a channel mapping JSON file"), message);
        assertTrue(message.contains("not a spot finder JSON file"), message);
        assertTrue(message.contains("_spotf_finding.json"), message);
    }

    @Test
    @DisplayName("the right JSON is read rather than merely accepted")
    void theRightJsonParses(@TempDir File directory) throws Exception {
        Map<String, Object> mapping = smFRETFiles.readMappingJSON(mappingJSON(directory).toString());
        assertTrue(mapping.containsKey("source points"));
        assertEquals(128, ((Number) mapping.get("image width")).intValue());

        Map<String, Object> finding = smFRETFiles.readSpotFinderJSON(spotFinderJSON(directory));
        assertEquals("movie.tif", finding.get("image name"));
    }

    @Test
    @DisplayName("JSON this pipeline did not write is distinguished from the wrong one of ours")
    void foreignJsonIsDistinguished(@TempDir File directory) throws Exception {
        File file = new File(directory, "something_else.json");
        Files.write(file.toPath(), "{\"unrelated\": 1}".getBytes(StandardCharsets.UTF_8));

        smFRETAnalysisException thrown = assertThrows(smFRETAnalysisException.class,
                () -> smFRETFiles.readMappingJSON(file.toString()));
        assertTrue(thrown.getMessage().contains("not one this pipeline wrote"), thrown.getMessage());
    }

    @Test
    @DisplayName("a movie offered where a mapping was wanted is named as a TIFF")
    void tiffIsNotAMapping(@TempDir File directory) {
        smFRETAnalysisException thrown = assertThrows(smFRETAnalysisException.class,
                () -> smFRETFiles.readMappingJSON(tiff(directory).toString()));
        assertTrue(thrown.getMessage().contains("is a TIFF image"), thrown.getMessage());
    }

    @Test
    @DisplayName("a JSON offered where the movie was wanted is named as JSON")
    void jsonIsNotAMovie(@TempDir File directory) {
        File mapping = mappingJSON(directory);
        smFRETAnalysisException thrown = assertThrows(smFRETAnalysisException.class,
                () -> smFRETFiles.openImage(mapping, "the image"));

        String message = thrown.getMessage();
        assertTrue(message.contains("is a channel mapping JSON file"), message);
        assertTrue(message.contains("not an image"), message);
        assertTrue(message.contains("two channel movie"), message);
    }

    @Test
    @DisplayName("a real TIFF opens")
    void tiffOpens(@TempDir File directory) {
        assertDoesNotThrow(() -> smFRETFiles.openImage(tiff(directory), "the image"));
    }

    @Test
    @DisplayName("HDF5 is recognized by its signature, both ways round")
    void hdf5IsRecognized(@TempDir File directory) throws Exception {
        assertDoesNotThrow(() -> smFRETFiles.requireHDF5(hdf5(directory)));

        smFRETAnalysisException thrown = assertThrows(smFRETAnalysisException.class,
                () -> smFRETFiles.requireHDF5(tiff(directory)));
        assertTrue(thrown.getMessage().contains("is a TIFF image"), thrown.getMessage());
        assertTrue(thrown.getMessage().contains("not an HDF5 trace file"), thrown.getMessage());
    }

    /**
     * ResultsTable.open2 does not fail on a binary file - it returns a table of whatever it made
     * of the bytes - so without this check the caller's missing-column test fires and blames an
     * older spot finder, which is a different problem with a different fix.
     */
    @Test
    @DisplayName("a binary file is refused as a spot table before it can be parsed")
    void binaryIsNotASpotTable(@TempDir File directory) throws Exception {
        assertDoesNotThrow(() -> smFRETFiles.requireSpotTable(spotTable(directory)));

        smFRETAnalysisException thrown = assertThrows(smFRETAnalysisException.class,
                () -> smFRETFiles.requireSpotTable(hdf5(directory)));
        assertTrue(thrown.getMessage().contains("not a spot table CSV"), thrown.getMessage());
        assertTrue(thrown.getMessage().contains("_spotf_spots.csv"), thrown.getMessage());
    }

    /**
     * The three states before the contents matter at all, each with its own message because each
     * is a different mistake: a moved output, a permissions problem, and a run that died partway
     * through writing.
     */
    @Test
    @DisplayName("missing, empty and directory are separate messages")
    void preContentChecksAreDistinct(@TempDir File directory) throws Exception {
        smFRETAnalysisException missing = assertThrows(smFRETAnalysisException.class,
                () -> smFRETFiles.requireReadable(new File(directory, "gone.json"), "a mapping"));
        assertTrue(missing.getMessage().contains("does not exist"), missing.getMessage());

        File empty = new File(directory, "empty.json");
        Files.write(empty.toPath(), new byte[0]);
        smFRETAnalysisException blank = assertThrows(smFRETAnalysisException.class,
                () -> smFRETFiles.requireReadable(empty, "a mapping"));
        assertTrue(blank.getMessage().contains("is empty"), blank.getMessage());

        smFRETAnalysisException folder = assertThrows(smFRETAnalysisException.class,
                () -> smFRETFiles.requireReadable(directory, "a mapping"));
        assertTrue(folder.getMessage().contains("is a directory"), folder.getMessage());

        smFRETAnalysisException none = assertThrows(smFRETAnalysisException.class,
                () -> smFRETFiles.requireReadable(null, "a mapping"));
        assertTrue(none.getMessage().contains("no file was given"), none.getMessage());
    }

    /**
     * Where the generated files go, from the movie's root name.
     *
     * The bare name is the case worth pinning: a macro that passes a relative path gets a root
     * with no directory at all, and the folder still has to come out beside it rather than in
     * whatever the working directory happens to be.
     */
    @Test
    @DisplayName("the analysis root is a folder beside the movie, named from it")
    void theAnalysisRootIsBesideTheMovie() {
        assertEquals(new File("/data/hel1_analysis/hel1").getPath(),
                smFRETFiles.analysisRoot("/data/hel1"));
        assertEquals(new File("hel1_analysis/hel1").getPath(),
                smFRETFiles.analysisRoot("hel1"));

        // A dot in the directory rather than the file - the root has already had the extension
        // stripped by the time it gets here, so nothing may be stripped again.
        assertEquals(new File("/data/run.2/hel1_analysis/hel1").getPath(),
                smFRETFiles.analysisRoot("/data/run.2/hel1"));
    }

    /**
     * The folder is made before anything is written to it. FileSaver.saveAsTiff returns a boolean
     * nobody checks, so a missing folder would otherwise be an analysis that reports success and
     * leaves nothing behind.
     */
    @Test
    @DisplayName("the analysis folder is created, and saying so beats a silent failure")
    void theAnalysisFolderIsCreated(@TempDir File directory) throws Exception {
        String root = smFRETFiles.createAnalysisRoot(new File(directory, "hel1").toString());
        File folder = new File(directory, "hel1_analysis");
        assertTrue(folder.isDirectory(), "no folder at " + folder);
        assertEquals(new File(folder, "hel1").getPath(), root);

        // Twice in a row is the ordinary case - re-running spot finding on the same movie.
        assertDoesNotThrow(() -> smFRETFiles.createAnalysisRoot(new File(directory, "hel1").toString()));

        // A file where the folder should go. Not a case a user hits often, but the alternative
        // is mkdirs returning false and every save afterwards quietly doing nothing.
        File blocked = new File(directory, "movie_analysis");
        Files.write(blocked.toPath(), "not a folder".getBytes(StandardCharsets.UTF_8));
        smFRETAnalysisException thrown = assertThrows(smFRETAnalysisException.class,
                () -> smFRETFiles.createAnalysisRoot(new File(directory, "movie").toString()));
        assertTrue(thrown.getMessage().contains("movie_analysis"), thrown.getMessage());
    }
}
