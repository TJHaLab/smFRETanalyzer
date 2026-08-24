import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.scijava.log.StderrLogService;

/**
 * Pooling several trace files into one histogram (issue #22).
 *
 * The rule that carries the risk is the frame count: one frame range slider governs every trace,
 * so the pool can only span the frames the *shortest* file has. Take the longest instead and
 * every loop reads off the end of a short file's row; take neither and the slider means a
 * different interval per file, which is worse than a crash because it produces a plausible
 * number. Half the tests here are that rule from a different angle.
 *
 * combine() is static and takes its input, so most of this runs without an H5 or a window. The
 * reading and ordering tests use real HDF5 files, written directly rather than through the
 * pipeline - what is under test is the pool, and a full stage 2 plus stage 3 run per file would
 * pay seconds for traces whose values nothing here looks at.
 */
class TracePoolTest {

    /** A file's worth of traces: every donor value i+1, every acceptor value -(i+1). */
    private static smFRETTraceHistogram.TraceFile traces(String name, int spots, int frames) {
        float[][] target = new float[spots][frames];
        float[][] source = new float[spots][frames];
        for (int i = 0; i < spots; i++) {
            Arrays.fill(target[i], i + 1.0f);
            Arrays.fill(source[i], -(i + 1.0f));
        }
        return new smFRETTraceHistogram.TraceFile(new File(name), target, source);
    }

    /** A real H5 holding just the two datasets the histogram reads. */
    private static File writeH5(File directory, String name, int spots, int frames) {
        File file = new File(directory, name);
        float[][] target = new float[spots][frames];
        float[][] source = new float[spots][frames];
        for (int i = 0; i < spots; i++) {
            Arrays.fill(target[i], i + 1.0f);
            Arrays.fill(source[i], -(i + 1.0f));
        }
        try (IHDF5Writer writer = HDF5Factory.configure(file).writer()) {
            writer.writeFloatMatrix("target-traces", target);
            writer.writeFloatMatrix("source-traces", source);
        }
        return file;
    }

    private static smFRETTraceHistogram histogram() {
        smFRETTraceHistogram plugin = new smFRETTraceHistogram();
        plugin.log = new StderrLogService();
        return plugin;
    }

    /**
     * Two files become one set of traces, in pool order and with nothing dropped.
     */
    @Test
    @DisplayName("pooled files are laid end to end")
    void filesAreLaidEndToEnd() {
        List<smFRETTraceHistogram.TraceFile> pool = new ArrayList<>();
        pool.add(traces("a.h5", 3, 10));
        pool.add(traces("b.h5", 2, 10));

        smFRETTraceHistogram.Combined combined = smFRETTraceHistogram.combine(pool);

        assertEquals(5, combined.nSpots(), "every trace from both files");
        assertEquals(10, combined.nFrames);

        // The first file's rows first, then the second's - which is what makes the pool order
        // shown in the panel the order the CSV header lists.
        for (int i = 0; i < 3; i++) {
            assertEquals(i + 1.0f, combined.target[i][0], "a.h5 trace " + i);
        }
        for (int i = 0; i < 2; i++) {
            assertEquals(i + 1.0f, combined.target[3 + i][0], "b.h5 trace " + i);
        }
    }

    /**
     * The shortest file sets the frame count. The point of the issue, and the one that would
     * otherwise index off the end of a row.
     */
    @Test
    @DisplayName("the shortest file sets the frame count")
    void theShortestFileSetsTheFrameCount() {
        List<smFRETTraceHistogram.TraceFile> pool = new ArrayList<>();
        pool.add(traces("long.h5", 2, 100));
        pool.add(traces("short.h5", 2, 40));
        pool.add(traces("middle.h5", 2, 60));

        assertEquals(40, smFRETTraceHistogram.combine(pool).nFrames);
    }

    /**
     * The shortest file is only the shortest of the *ticked* ones - untick it and the pool gets
     * its frames back. A frame count taken over the whole pool would leave the histogram
     * truncated to a file that is no longer in it, with nothing on screen saying why.
     */
    @Test
    @DisplayName("unticking the shortest file restores the frame range")
    void untickingTheShortestFileRestoresTheRange() {
        List<smFRETTraceHistogram.TraceFile> pool = new ArrayList<>();
        pool.add(traces("long.h5", 2, 100));
        smFRETTraceHistogram.TraceFile shortFile = traces("short.h5", 3, 40);
        pool.add(shortFile);

        assertEquals(40, smFRETTraceHistogram.combine(pool).nFrames);
        assertEquals(5, smFRETTraceHistogram.combine(pool).nSpots());

        shortFile.included = false;
        smFRETTraceHistogram.Combined combined = smFRETTraceHistogram.combine(pool);
        assertEquals(100, combined.nFrames, "the unticked file should not still be limiting");
        assertEquals(2, combined.nSpots(), "its traces should be gone too");
    }

    /**
     * Everything unticked is a state the panel can be put in with two clicks, so it has to be a
     * plot that says so rather than an exception.
     *
     * The filter bounds are the trap: computeFilterBounds() walks the traces to find them, and
     * with no traces to walk both ends are left at their sentinels, which the slider rescale
     * casts to an int.
     */
    @Test
    @DisplayName("an empty pool is a state, not an error")
    void anEmptyPoolIsSafe() {
        List<smFRETTraceHistogram.TraceFile> pool = new ArrayList<>();
        smFRETTraceHistogram.TraceFile only = traces("a.h5", 4, 10);
        only.included = false;
        pool.add(only);

        smFRETTraceHistogram.Combined combined = smFRETTraceHistogram.combine(pool);
        assertEquals(0, combined.nSpots());
        assertEquals(0, combined.nFrames);

        smFRETTraceHistogram plugin = histogram();
        plugin.pool.add(only);
        plugin.rebuild();
        assertEquals(0, plugin.nSpots);

        for (int f = 0; f < smFRETTraceHistogram.N_FILTERS; f++) {
            assertTrue(plugin.filterMin[f] > -1.0e9 && plugin.filterMax[f] < 1.0e9,
                    "filter " + f + " bounds left at their sentinels: "
                            + plugin.filterMin[f] + " to " + plugin.filterMax[f]);
            assertTrue(plugin.filterMax[f] > plugin.filterMin[f], "filter " + f + " is inverted");
        }

        // And it still bins, because the plot panel asks it to before it draws the message.
        smFRETTraceHistogram.Histogram result = plugin.computeHistogram(
                smFRETTraceHistogram.TYPE_FRET, 1, 1,
                smFRETTraceHistogram.Filters.none(), 10, new smFRETTraceHistogram.Corrections(0, 0, 0));
        assertEquals(0, result.nSpotsUsed);
    }

    /**
     * Files read from disk, pooled, and put in name order however they arrived.
     *
     * Name order rather than added order because a pool is assembled by dropping a folder in,
     * where the drop order is whatever the file manager felt like, and the names are the only
     * thing the user picked.
     */
    @Test
    @DisplayName("pooled files are read and held in name order")
    void filesAreReadAndOrdered(@TempDir File directory) {
        File c = writeH5(directory, "c.h5", 1, 10);
        File a = writeH5(directory, "a.h5", 2, 10);
        File b = writeH5(directory, "B.h5", 3, 10);

        smFRETTraceHistogram plugin = histogram();
        assertTrue(plugin.addTraceFiles(Arrays.asList(c, a, b)).isEmpty(), "all three should read");

        List<String> names = new ArrayList<>();
        for (smFRETTraceHistogram.TraceFile traces : plugin.pool) {
            names.add(traces.name());
        }
        assertEquals(Arrays.asList("a.h5", "B.h5", "c.h5"), names,
                "name order, and case should not put B after c");
        assertEquals(6, plugin.nSpots, "every trace from all three");
    }

    /**
     * The same file twice adds it once. Dropping a folder onto a pool that already holds some of
     * it is the ordinary way to get here, and the alternative is a file silently counted twice -
     * which does not look like a mistake in the histogram, it just weights that movie double.
     */
    @Test
    @DisplayName("the same file is not pooled twice")
    void theSameFileIsNotPooledTwice(@TempDir File directory) {
        File a = writeH5(directory, "a.h5", 4, 10);

        smFRETTraceHistogram plugin = histogram();
        plugin.addTraceFiles(Arrays.asList(a));
        plugin.addTraceFiles(Arrays.asList(a, new File(directory, "./a.h5")));

        assertEquals(1, plugin.pool.size(), "one entry per file");
        assertEquals(4, plugin.nSpots, "and its traces counted once");
    }

    /**
     * A file that will not read is reported and the rest of the drop is kept.
     *
     * The usual way to get one is sweeping up a stray while selecting the wanted files, and
     * failing the whole drop would make the user redo a selection that was almost right.
     */
    @Test
    @DisplayName("a file that will not read leaves the others pooled")
    void aBadFileLeavesTheOthersPooled(@TempDir File directory) throws Exception {
        File good = writeH5(directory, "good.h5", 3, 10);
        File bad = new File(directory, "notes.h5");
        Files.write(bad.toPath(), "this is not an H5 file".getBytes());

        smFRETTraceHistogram plugin = histogram();
        List<String> problems = plugin.addTraceFiles(Arrays.asList(good, bad));

        assertEquals(1, problems.size(), "one problem reported: " + problems);
        assertTrue(problems.get(0).contains("notes.h5"), problems.get(0));
        assertEquals(1, plugin.pool.size(), "the good file should still be pooled");
        assertEquals(3, plugin.nSpots);
    }

    /**
     * The seeding file is different: it is the one chosen in the dialog, so there is nothing left
     * to carry on with and the plugin has to say so rather than open an empty window.
     */
    @Test
    @DisplayName("a bad seeding file is an error, not a warning")
    void aBadSeedingFileThrows(@TempDir File directory) throws Exception {
        File bad = new File(directory, "notes.h5");
        Files.write(bad.toPath(), "this is not an H5 file".getBytes());

        smFRETTraceHistogram plugin = histogram();
        smFRETAnalysisException thrown = assertThrows(smFRETAnalysisException.class,
                () -> plugin.loadTraces(bad));
        assertTrue(thrown.getMessage().contains("notes.h5"), thrown.getMessage());
    }

    /**
     * A file that could not be read is reported with its path wrapped, not as one enormous line.
     *
     * JOptionPane sizes itself to its longest line, and these messages name a path because naming
     * the file is the whole of what they are for - so without a width the dialog reporting a
     * stray dropped file is wider than the screen.
     */
    @Test
    @DisplayName("the not-added message is wrapped, not one line per path")
    void theProblemDialogIsWrapped(@TempDir File directory) throws Exception {
        File bad = new File(directory, "notes.h5");
        Files.write(bad.toPath(), "this is not an H5 file".getBytes());

        smFRETTraceHistogram plugin = histogram();
        List<String> problems = plugin.addTraceFiles(Arrays.asList(
                writeH5(directory, "good.h5", 2, 10), bad));
        assertEquals(1, problems.size(), problems.toString());

        String shown = smFRETSwing.dialogMessage("Not added:\n" + problems.get(0));
        assertTrue(shown.contains("width:"), "no width, so the dialog grows to the path: " + shown);
        assertFalse(shown.contains("\n"), "newlines should be line breaks the dialog honours");
        assertTrue(shown.contains("\u200b"), "the path should be able to break");
    }

    /**
     * Ticking is a view of the pool, not a change to it - the file stays loaded so that ticking
     * it back on does not re-read it.
     */
    @Test
    @DisplayName("unticking keeps the file in the pool")
    void untickingKeepsTheFile(@TempDir File directory) {
        smFRETTraceHistogram plugin = histogram();
        plugin.addTraceFiles(Arrays.asList(writeH5(directory, "a.h5", 2, 10),
                writeH5(directory, "b.h5", 3, 10)));
        assertEquals(5, plugin.nSpots);

        plugin.pool.get(0).included = false;
        plugin.rebuild();

        assertEquals(2, plugin.pool.size(), "the file should still be pooled");
        assertEquals(1, plugin.includedFiles().size());
        assertEquals(3, plugin.nSpots, "only the ticked file's traces");
        assertFalse(plugin.poolLabel().contains("2 files"), plugin.poolLabel());
    }
}
