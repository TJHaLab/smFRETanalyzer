import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.HashMap;
import java.util.Map;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

/**
 * Settings survive from one session to the next (issue #7).
 *
 * <p>SciJava's own persistence does not cover this plugin: it saves a module's inputs when
 * the module completes, and pressing Find spots invokes a callback rather than running the
 * module, so a value typed into the dialog was never saved while a value passed by a macro
 * was. The round trip is therefore done by hand, and this tests it without needing a SciJava
 * context or a dialog.
 */
class SpotFinderSettingsTest {

    private static smFRETSpotFinder configured() {
        smFRETSpotFinder finder = new smFRETSpotFinder();
        finder.startSlice = 5;
        finder.endSlice = 44;
        finder.spotThreshold = 7.5;
        finder.spotTolerance = 3.0;
        finder.spotContamination = 0.27;
        finder.spotSigma = 1.7;
        finder.cameraBlackLevel = 100;
        finder.cameraGain = 0.48;
        finder.spotSpacing = 4;
        finder.edgeMargin = 9;
        finder.backgroundKappa = 1.3;
        finder.spotChannel = "donor";
        return finder;
    }

    @Test
    @DisplayName("every saved setting comes back")
    void theRoundTripIsFaithful() {
        Map<String, String> saved = configured().currentSettings();

        smFRETSpotFinder fresh = new smFRETSpotFinder();
        fresh.applySettings(saved);

        assertEquals(5, fresh.startSlice.intValue(), "startSlice");
        assertEquals(44, fresh.endSlice.intValue(), "endSlice");
        assertEquals(7.5, fresh.spotThreshold, 1.0e-9, "spotThreshold");
        assertEquals(3.0, fresh.spotTolerance, 1.0e-9, "spotTolerance");
        assertEquals(0.27, fresh.spotContamination, 1.0e-9, "spotContamination");
        assertEquals(1.7, fresh.spotSigma, 1.0e-9, "spotSigma");
        assertEquals(100, fresh.cameraBlackLevel.intValue(), "cameraBlackLevel");
        assertEquals(0.48, fresh.cameraGain, 1.0e-9, "cameraGain");
        assertEquals(4, fresh.spotSpacing.intValue(), "spotSpacing");
        assertEquals(9, fresh.edgeMargin.intValue(), "edgeMargin");
        assertEquals(1.3, fresh.backgroundKappa, 1.0e-9, "backgroundKappa");
        assertEquals("donor", fresh.spotChannel, "spotChannel");
    }

    /**
     * The slice range is persisted on purpose - working through a folder of movies with one
     * range is the case the issue was raised about.
     */
    @Test
    @DisplayName("the slice range is among the settings")
    void theSliceRangePersists() {
        Map<String, String> saved = configured().currentSettings();
        assertEquals("5", saved.get("startSlice"));
        assertEquals("44", saved.get("endSlice"));
    }

    /**
     * The file paths are deliberately absent. Restoring them would leave a fresh session
     * pointed at the last movie, looking ready to run and silently analysing the wrong data.
     */
    @Test
    @DisplayName("file paths are never persisted")
    void filePathsAreNotSaved() {
        Map<String, String> saved = configured().currentSettings();
        assertTrue(!saved.containsKey("inputImageName"), "the movie path was saved");
        assertTrue(!saved.containsKey("mappingFile"), "the mapping path was saved");
    }

    /**
     * Stored preferences outlive the code that wrote them - the real file still holds a key
     * for spotProminence, which no longer exists - so a stale or malformed entry has to leave
     * its parameter at the default rather than stop the plugin opening.
     */
    @Test
    @DisplayName("unparseable and unknown settings are ignored, not fatal")
    void junkIsSurvivable() {
        Map<String, String> junk = new HashMap<>();
        junk.put("spotThreshold", "not a number");
        junk.put("spotSigma", "");
        junk.put("spotChannel", "ultraviolet");
        junk.put("spotProminence", "0.4");
        junk.put("somethingFromTheFuture", "17");

        smFRETSpotFinder fresh = new smFRETSpotFinder();
        double threshold = fresh.spotThreshold;
        fresh.applySettings(junk);

        assertEquals(threshold, fresh.spotThreshold, 1.0e-9, "left at the default");
        assertEquals("sum", fresh.spotChannel, "an unknown channel is refused");
    }

    @Test
    @DisplayName("no saved settings at all is not an error")
    void nothingSavedIsFine() {
        smFRETSpotFinder fresh = new smFRETSpotFinder();
        fresh.applySettings(null);
        fresh.applySettings(new HashMap<>());
        assertEquals(6.0, fresh.spotThreshold, 1.0e-9);
    }

    /**
     * SciJava's own persistence must be off for the two file parameters.
     *
     * <p>It is a separate mechanism from the settings map above: SciJava saves every
     * @Parameter when a module completes, so a macro run writes the paths and the next dialog
     * restores them. Excluding them from currentSettings() does not prevent that - only
     * persist = false does. Asserted on the annotation because the failure is invisible in
     * any headless test: the plugin works perfectly and just opens pointing at the wrong file.
     */
    @Test
    @DisplayName("SciJava is told not to persist the file parameters")
    void filePathsHavePersistDisabled() throws NoSuchFieldException {
        for (String name : new String[] {"inputImageName", "mappingFile"}) {
            org.scijava.plugin.Parameter annotation = smFRETSpotFinder.class
                    .getDeclaredField(name)
                    .getAnnotation(org.scijava.plugin.Parameter.class);
            assertTrue(annotation != null, name + " has no @Parameter");
            assertTrue(!annotation.persist(),
                    name + " would be restored into a new session by SciJava");
        }
    }
}
