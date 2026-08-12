import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

/**
 * The Java ensemble has to give the same answers as the scikit-learn one it was exported from.
 *
 * <p>A tree ensemble is easy to port and easy to port <i>subtly wrong</i>: children are stored
 * as offsets within their own tree, so forgetting the tree base reads real values out of the
 * wrong nodes and returns plausible numbers rather than crashing. Nothing but a fixture
 * catches that. {@code spot_quality_fixture.txt} is 200 feature vectors spread across the
 * range of the training predictions, with what scikit-learn returns for each.
 */
class SpotQualityForestTest {

    private static final String FIXTURE = "/spot_quality_fixture.txt";

    private static final class Fixture {
        final List<double[]> features = new ArrayList<>();
        final List<Double> expected = new ArrayList<>();
    }

    private static Fixture fixture() {
        Fixture out = new Fixture();
        BufferedReader reader = new BufferedReader(new InputStreamReader(
                SpotQualityForestTest.class.getResourceAsStream(FIXTURE),
                Charset.forName("UTF-8")));
        try {
            String line;
            while ((line = reader.readLine()) != null) {
                int arrow = line.indexOf("->");
                if (line.startsWith("#") || arrow < 0) {
                    continue;
                }
                String[] parts = line.substring(0, arrow).trim().split(" +");
                double[] values = new double[parts.length];
                for (int i = 0; i < parts.length; i++) {
                    values[i] = Double.parseDouble(parts[i]);
                }
                out.features.add(values);
                out.expected.add(Double.parseDouble(line.substring(arrow + 2).trim()));
            }
            reader.close();
        } catch (IOException e) {
            throw new IllegalStateException(e);
        }
        return out;
    }

    @Test
    @DisplayName("the shipped model loads and reports what it is")
    void theModelLoads() {
        smFRETSpotQuality.Forest forest = smFRETSpotQuality.shipped();

        assertEquals(100, forest.treeCount(), "trees");
        assertEquals(43, forest.featureCount(), "features");
        assertEquals(0.2, forest.defaultThreshold(), 1.0e-9, "default threshold");
        assertTrue(forest.nodeCount() > 500, "nodes: " + forest.nodeCount());
    }

    /**
     * The whole point of the fixture. Every prediction, to well below anything that could
     * change a filtering decision.
     */
    @Test
    @DisplayName("every fixture prediction matches scikit-learn")
    void predictionsMatchTheTrainer() {
        smFRETSpotQuality.Forest forest = smFRETSpotQuality.shipped();
        Fixture fixture = fixture();
        assertTrue(fixture.features.size() >= 100,
                "only " + fixture.features.size() + " fixture rows");

        double worst = 0.0;
        for (int i = 0; i < fixture.features.size(); i++) {
            double actual = forest.raw(fixture.features.get(i));
            worst = Math.max(worst, Math.abs(actual - fixture.expected.get(i)));
        }
        assertTrue(worst < 1.0e-9, "worst disagreement with scikit-learn was " + worst);
    }

    /**
     * The ensemble is a sum of regression trees with nothing bounding it, and its training
     * predictions really do run slightly negative, so the clamp is load bearing rather than
     * defensive.
     */
    @Test
    @DisplayName("predictions are clamped to a fraction")
    void predictionsAreAFraction() {
        smFRETSpotQuality.Forest forest = smFRETSpotQuality.shipped();
        Fixture fixture = fixture();

        boolean sawNegativeRaw = false;
        for (double[] features : fixture.features) {
            double predicted = forest.predict(features);
            assertTrue((predicted >= 0.0) && (predicted <= 1.0),
                    "prediction outside [0, 1]: " + predicted);
            sawNegativeRaw |= forest.raw(features) < 0.0;
        }
        assertTrue(sawNegativeRaw,
                "no fixture row exercised the clamp, so it is not being tested");
    }

    @Test
    @DisplayName("the PSF envelope is carried with the model and is not open ended")
    void theEnvelopeIsRecorded() {
        smFRETSpotQuality.Forest forest = smFRETSpotQuality.shipped();

        assertTrue(forest.covers(1.5, 0.3), "the PSF it was trained at");
        assertTrue(forest.covers(1.34, 0.42), "the example data's PSF");
        assertTrue(!forest.covers(0.6, 0.3), "a core far narrower than anything trained");
        assertTrue(!forest.covers(1.5, 0.9), "past where the Airy model reorganizes");
    }

    /**
     * The response table is what lets the neighbour feature be computed without any Airy
     * maths at run time. It falls away with distance and reaches zero past its end.
     *
     * <p>The sweep starts at 1 pixel rather than 0 because of the next test.
     */
    @Test
    @DisplayName("the neighbour response falls with distance and ends at zero")
    void theResponseTableFallsAway() {
        smFRETSpotQuality.Forest forest = smFRETSpotQuality.shipped();

        double previous = forest.responseAt(1.0);
        assertTrue(previous > 0.0, "a spot must respond to a neighbour one pixel away");
        for (double distance = 2.0; distance <= 20.0; distance += 1.0) {
            double here = forest.responseAt(distance);
            assertTrue(here < previous, "response did not fall by " + distance + " px");
            previous = here;
        }
        assertEquals(0.0, forest.responseAt(1000.0), 0.0, "past the table");
    }

    /**
     * The table peaks at half a pixel, not at zero, and is symmetric about it - so a
     * neighbour one pixel away weighs exactly as much as the spot itself.
     *
     * <p>This looks like a bug and is not. The plugin reads the blurred image at a *pixel*,
     * so the measurement kernel sits on that pixel's centre, while the table walks a
     * neighbour outwards from the pixel's integer coordinate. A neighbour at 0 and a
     * neighbour at 1 are then both half a pixel from the kernel centre and must weigh the
     * same. Integer coordinates are the right convention here because detections are integer
     * coordinates - MaximumFinder reports the brightest pixel, not a centroid - so the
     * separations this table is ever asked about are hypotenuses of integers.
     *
     * <p>Pinned because a future reader will otherwise 'fix' the asymmetry and silently
     * shift every neighbour feature the model was trained on.
     */
    @Test
    @DisplayName("the response peaks half a pixel out, and that is the pixel grid")
    void theResponsePeaksOffCentre() {
        smFRETSpotQuality.Forest forest = smFRETSpotQuality.shipped();

        // Equal to a few parts in 1e8 rather than exactly. The symmetry is exact in
        // principle; the table is not, because an Airy pattern has no closed form pixel
        // integral and is sub-sampled. An earlier version of this file asserted 1e-9 and
        // passed only because the table was being written with seven digits.
        assertEquals(forest.responseAt(0.0), forest.responseAt(1.0), 1.0e-7,
                "zero and one pixel are equidistant from the kernel centre");
        assertTrue(forest.responseAt(0.5) > forest.responseAt(0.0),
                "the peak sits at the pixel centre, half a pixel out");
    }

    @Test
    @DisplayName("a truncated model file is an error with a message, not a wrong answer")
    void aTruncatedModelIsRefused() {
        String truncated = "version 1\nfeatures 3\nbaseline 0.1\ntree 0 3\n0 0 1.0 1 2 0.0\n";
        smFRETAnalysisException thrown = assertThrows(smFRETAnalysisException.class,
                () -> smFRETSpotQuality.load(new ByteArrayInputStream(
                        truncated.getBytes(Charset.forName("UTF-8")))));
        assertTrue(thrown.getMessage().contains("part way through"), thrown.getMessage());
    }

    @Test
    @DisplayName("a feature vector of the wrong length is refused")
    void theFeatureCountIsChecked() {
        smFRETSpotQuality.Forest forest = smFRETSpotQuality.shipped();
        smFRETAnalysisException thrown = assertThrows(smFRETAnalysisException.class,
                () -> forest.predict(new double[] {1.0, 2.0, 3.0}));
        assertTrue(thrown.getMessage().contains("43"), thrown.getMessage());
    }
}
