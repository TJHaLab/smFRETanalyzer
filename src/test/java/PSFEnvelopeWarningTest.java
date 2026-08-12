import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

/**
 * The warning that the measured PSF is one SpotContamination was never tested on.
 *
 * <p>This is the only place in the plugin that measures a PSF, so it is the only place that
 * can notice. The failure it guards against is quiet and is not a graceful degradation:
 * outside the trained range the contamination score can rank spots *worse* than the untrained
 * prominence statistic, so silence would leave a user with an unusual PSF trusting a filter
 * that may be actively misleading.
 */
class PSFEnvelopeWarningTest {

    /**
     * Only sigma and waves matter here; the rest of a Fit is not consulted.
     *
     * <p>Note the constructor takes sigma and waves *first*, which is not the order the
     * fields are declared in - passing them in field order compiles cleanly and silently
     * hands the check an amplitude where it wanted a sigma.
     */
    private static smFRETPSF.Fit fit(double sigma, double waves) {
        return new smFRETPSF.Fit(sigma, waves, 1.0, 0.0, 0.01);
    }

    @Test
    @DisplayName("nothing is said about a PSF the model was trained on")
    void silentInsideTheEnvelope() {
        assertEquals("", smFRETPSFVisualizer.envelopeWarning(fit(1.5, 0.3), fit(1.6, 0.35)));
        // The example data's own PSF, which had better not warn.
        assertEquals("", smFRETPSFVisualizer.envelopeWarning(fit(1.34, 0.42), fit(1.55, 0.36)));
        // The corners of the range are inside it.
        assertEquals("", smFRETPSFVisualizer.envelopeWarning(fit(1.0, 0.0), fit(2.5, 0.5)));
    }

    @Test
    @DisplayName("the warning names which channel is out of range")
    void itNamesTheChannel() {
        String donorOnly = smFRETPSFVisualizer.envelopeWarning(fit(0.7, 0.3), fit(1.5, 0.3));
        String acceptorOnly = smFRETPSFVisualizer.envelopeWarning(fit(1.5, 0.3), fit(2.9, 0.3));
        String both = smFRETPSFVisualizer.envelopeWarning(fit(0.7, 0.3), fit(2.9, 0.8));

        assertTrue(donorOnly.startsWith("The donor channel is"), donorOnly);
        assertTrue(acceptorOnly.startsWith("The acceptor channel is"), acceptorOnly);
        assertTrue(both.startsWith("Both channels are"), both);
    }

    /**
     * Aberration is the axis that matters most - it is where a narrowly trained model fell
     * below prominence - so a core inside the range with too much aberration must still warn.
     */
    @Test
    @DisplayName("too much aberration warns even with a normal core")
    void aberrationAloneTriggersIt() {
        assertTrue(!smFRETPSFVisualizer.envelopeWarning(fit(1.5, 0.9), fit(1.5, 0.3)).isEmpty(),
                "0.9 waves is past where the Airy model even stays monotone");
    }

    @Test
    @DisplayName("the warning says what the range is, so it can be acted on")
    void itQuotesTheRange() {
        String warning = smFRETPSFVisualizer.envelopeWarning(fit(0.7, 0.3), fit(1.5, 0.3));
        assertTrue(warning.contains("1.0 to 2.5"), warning);
        assertTrue(warning.contains("0.0 to 0.5"), warning);
        assertTrue(warning.contains("SpotContamination"), warning);
    }
}
