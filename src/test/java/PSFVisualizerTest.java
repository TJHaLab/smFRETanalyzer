import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import ij.ImagePlus;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.scijava.log.StderrLogService;

/**
 * The plugin's plumbing: what it reads, what it does not need, and what it says when a file is
 * missing.
 *
 * The measurement itself is scored in PSFMeasurementTest against a field of known PSF; this is
 * about the wiring around it. The one behaviour worth stating outright is that a PSF is measured
 * from the *field* rather than from the traces, so unlike smFRETTraceVisualizer this does not
 * need the trace H5 and runs as soon as spot finding has - which is when it is most useful, since
 * the PSF is what says whether spotSigma was a sensible choice in the first place.
 */
class PSFVisualizerTest {

    private static final int HALF_WIDTH = 160;
    private static final int HEIGHT = 160;
    private static final int FRAMES = 4;

    /**
     * A field run through spot finding, and nothing else - deliberately no stage 3.
     *
     * Returns the spot finder JSON.
     */
    private static File spotFinderOutput(File directory) {
        StderrLogService log = new StderrLogService();

        List<SyntheticField.Spot> spots = new ArrayList<>();
        for (int y = 30; y <= (HEIGHT - 31); y += 32) {
            for (int x = 30; x <= (HALF_WIDTH - 31); x += 32) {
                spots.add(new SyntheticField.Spot(x, y, 3000.0));
            }
        }

        ImagePlus movie = SyntheticField.movie(HALF_WIDTH, HEIGHT, spots, spots, 2.0,
                20.0, FRAMES, 0.0, 5L);
        File movieFile = SyntheticField.writeMovie(directory, "psf.tif", movie);
        File mappingFile = SyntheticField.writeIdentityMapping(directory, "psf_mapping.json",
                2 * HALF_WIDTH, HEIGHT);

        smFRETSpotFinder finder = new smFRETSpotFinder();
        finder.log = log;
        finder.inputImageName = movieFile;
        finder.mappingFile = mappingFile;
        finder.startSlice = 1;
        finder.endSlice = FRAMES;
        finder.spotSigma = 2.0;
        finder.spotThreshold = 6.0;
        finder.spotTolerance = 5.0;
        finder.spotContamination = 1.0;   // no quality filtering in this test
        finder.cameraBlackLevel = 0;
        finder.cameraGain = 1.0;
        finder.spotSpacing = 3;
        finder.edgeMargin = 5;
        finder.backgroundKappa = 0.0;
        finder.findSpots();

        return new File(directory, "psf_spotf_finding.json");
    }

    private static smFRETPSFVisualizer visualizer(File json) {
        smFRETPSFVisualizer plugin = new smFRETPSFVisualizer();
        plugin.log = new StderrLogService();
        plugin.spotJSONFile = json;
        return plugin;
    }

    /**
     * The difference from smFRETTraceVisualizer, stated as a test: no trace H5 anywhere, and it
     * still loads. That plugin refuses without one, and it is right to - it shows traces.
     */
    @Test
    @DisplayName("it loads without a trace h5, which stage 3 has not written")
    void loadsWithoutTraces(@TempDir File directory) {
        File json = spotFinderOutput(directory);
        assertFalse(new File(directory, "psf.h5").exists(),
                "the fixture should not have run stage 3");

        smFRETPSFVisualizer plugin = visualizer(json);
        plugin.load(json);
    }

    /**
     * A headless run measures and reports rather than trying to open a window. That is what a
     * macro driving this would be after, and it is also what makes the whole path testable.
     */
    @Test
    @DisplayName("a headless run measures both channels")
    void headlessRunMeasures(@TempDir File directory) {
        File json = spotFinderOutput(directory);
        smFRETPSFVisualizer plugin = visualizer(json);

        plugin.run();

        // run() catches its own exceptions and reports them, so the check that it worked is that
        // the state it builds is there and sane.
        for (int c = 0; c < 2; c++) {
            smFRETPSF.Measurement measured = plugin.measurement(c);
            assertTrue(measured != null, "channel " + c + " was not measured");
            assertTrue(measured.samples.spotsUsed > 0,
                    "channel " + c + " used no spots of " + measured.samples.spotsTotal);
            assertTrue(measured.withPedestal != null, "channel " + c + " did not fit");

            // A field of Gaussians is not an Airy pattern, so the fitted aberration is a
            // compromise rather than a truth to check against - what matters here is that the
            // whole path ran and produced a core width in the region of the spots it was given.
            assertTrue((measured.withPedestal.sigma > 0.8) && (measured.withPedestal.sigma < 4.0),
                    "channel " + c + " fitted sigma " + measured.withPedestal.sigma);
        }
    }

    /**
     * A movie or a spot table that has moved since spot finding ran is the usual failure, and the
     * message has to name the file rather than arriving as a null somewhere downstream.
     */
    @Test
    @DisplayName("a missing input is named")
    void missingInputIsNamed(@TempDir File directory) {
        File json = spotFinderOutput(directory);
        assertTrue(new File(directory, "psf_analysis/psf_spotf_spots.csv").delete(),
                "could not stage the fault");

        smFRETPSFVisualizer plugin = visualizer(json);
        smFRETAnalysisException thrown = assertThrows(smFRETAnalysisException.class,
                () -> plugin.load(json));
        assertTrue(thrown.getMessage().contains("psf_spotf_spots.csv"), thrown.getMessage());
    }

    /**
     * The panels lay out and paint without a display.
     *
     * They had no coverage at all until a caption was found clipped off the bottom of the image
     * panels, which no test could have caught because nothing ever painted them. buildContent()
     * hands back a plain JPanel precisely so this is possible - a JFrame cannot even be
     * constructed headless.
     *
     * This does not check that the result *looks* right; it checks that every paint path runs,
     * which covers the log floor, the legend, the NaN branch and the fit curve at sizes where
     * things are tight enough to go negative.
     */
    @Test
    @DisplayName("the panels paint at any size without throwing")
    void panelsPaint(@TempDir File directory) throws Exception {
        File json = spotFinderOutput(directory);
        smFRETPSFVisualizer plugin = visualizer(json);
        plugin.run();

        for (int[] size : new int[][] {{1000, 760}, {640, 480}, {320, 240}, {120, 90}}) {
            javax.swing.JPanel content = plugin.buildContent();

            // The same sequence show() runs. Without it the labels are still blank and the
            // panels paint an empty measurement, which exercises far less of the drawing.
            plugin.update();
            content.setSize(size[0], size[1]);
            layOut(content);

            BufferedImage image = new BufferedImage(size[0], size[1],
                    BufferedImage.TYPE_INT_RGB);
            Graphics2D g2 = image.createGraphics();
            content.paint(g2);
            g2.dispose();

            // Something was drawn - a panel that silently painted nothing would pass a
            // "did not throw" test on its own.
            if (size[0] >= 320) {
                boolean marked = false;
                for (int y = 0; (y < size[1]) && !marked; y++) {
                    for (int x = 0; x < size[0]; x++) {
                        if ((image.getRGB(x, y) & 0xffffff) != 0xffffff) {
                            marked = true;
                            break;
                        }
                    }
                }
                assertTrue(marked, "nothing was drawn at " + size[0] + "x" + size[1]);
            }
        }
    }

    private static void layOut(java.awt.Component c) {
        c.doLayout();
        if (c instanceof java.awt.Container) {
            for (java.awt.Component child : ((java.awt.Container) c).getComponents()) {
                layOut(child);
            }
        }
    }

    /**
     * The spot table is read in the reloaded layout - x at column 0 - which is the one the
     * measurement indexes. smFRETSpotFinder holds the other one, flag-prefixed, and reading a
     * patch centred on a flag would put every patch at the top left corner.
     */
    @Test
    @DisplayName("spots are read in the reloaded layout")
    void spotsAreInTheReloadedLayout(@TempDir File directory) {
        File json = spotFinderOutput(directory);
        smFRETPSFVisualizer plugin = visualizer(json);
        plugin.load(json);

        double[][] spots = plugin.spots();
        assertTrue(spots.length > 0, "no spots loaded");
        for (double[] spot : spots) {
            // Width follows smFRETSpotFinder.columnHeaders, and this assertion is here to
            // fail when a filter is added without the headers being updated with it - the
            // round trip only lines up while the two agree.
            assertEquals(5, spot.length, "x, y, snr, prominence, contamination");
            assertTrue((spot[0] >= 0) && (spot[0] < HALF_WIDTH),
                    "x out of the half frame: " + spot[0]);
            assertTrue((spot[1] >= 0) && (spot[1] < HEIGHT),
                    "y out of the half frame: " + spot[1]);
        }
    }
}
