import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

/**
 * The 43 numbers this computes have to be the 43 numbers the model was trained on.
 *
 * <p>The forest fixture proves the ensemble evaluates correctly; it says nothing about
 * whether it is being fed the right thing. That is where the porting risk actually is - a
 * patch off by a pixel, a median over the wrong ring, a kernel truncated at a different
 * radius. Every one of those produces plausible features and a quietly wrong score, and no
 * amount of reading the code catches it.
 *
 * <p>{@code spot_quality_feature_fixture.txt} carries a small image <b>written out in full</b>
 * rather than described by a formula, so the two sides cannot build subtly different
 * pictures, along with the features Python computes from it. The Python that produced it was
 * itself checked against the stored dataset for three real fields and agreed exactly, so this
 * is a comparison against the pipeline that trained the model rather than against a second
 * guess at it.
 */
class SpotQualityFeatureTest {

    private static final String FIXTURE = "/spot_quality_feature_fixture.txt";

    private static final class Fixture {
        int width;
        int height;
        double gain;
        double black;
        int frames;
        float[] image;
        double[] x;
        double[] y;
        double[] snr;
        double[] prominence;
        List<double[]> features = new ArrayList<>();
        double[] predictions;
    }

    private static Fixture fixture() {
        Fixture out = new Fixture();
        List<double[]> spots = new ArrayList<>();
        List<Float> pixels = new ArrayList<>();

        BufferedReader reader = new BufferedReader(new InputStreamReader(
                SpotQualityFeatureTest.class.getResourceAsStream(FIXTURE),
                Charset.forName("UTF-8")));
        try {
            String line;
            boolean readingImage = false;
            boolean readingFeatures = false;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) {
                    continue;
                }
                String[] parts = line.split(" +");
                if ("width".equals(parts[0])) {
                    out.width = Integer.parseInt(parts[1]);
                } else if ("height".equals(parts[0])) {
                    out.height = Integer.parseInt(parts[1]);
                } else if ("gain".equals(parts[0])) {
                    out.gain = Double.parseDouble(parts[1]);
                } else if ("black".equals(parts[0])) {
                    out.black = Double.parseDouble(parts[1]);
                } else if ("frames".equals(parts[0])) {
                    out.frames = Integer.parseInt(parts[1]);
                } else if ("spot".equals(parts[0])) {
                    spots.add(new double[] {Double.parseDouble(parts[1]),
                            Double.parseDouble(parts[2]), Double.parseDouble(parts[3]),
                            Double.parseDouble(parts[4])});
                } else if ("image".equals(parts[0])) {
                    readingImage = true;
                    readingFeatures = false;
                } else if ("features".equals(parts[0])) {
                    readingImage = false;
                    readingFeatures = true;
                } else if ("predictions".equals(parts[0])) {
                    readingImage = false;
                    readingFeatures = false;
                    String values = reader.readLine().trim();
                    String[] each = values.split(" +");
                    out.predictions = new double[each.length];
                    for (int i = 0; i < each.length; i++) {
                        out.predictions[i] = Double.parseDouble(each[i]);
                    }
                } else if (readingImage) {
                    for (String value : parts) {
                        pixels.add((float) Double.parseDouble(value));
                    }
                } else if (readingFeatures) {
                    double[] values = new double[parts.length];
                    for (int i = 0; i < parts.length; i++) {
                        values[i] = Double.parseDouble(parts[i]);
                    }
                    out.features.add(values);
                }
            }
            reader.close();
        } catch (IOException e) {
            throw new IllegalStateException(e);
        }

        out.image = new float[pixels.size()];
        for (int i = 0; i < out.image.length; i++) {
            out.image[i] = pixels.get(i);
        }
        out.x = new double[spots.size()];
        out.y = new double[spots.size()];
        out.snr = new double[spots.size()];
        out.prominence = new double[spots.size()];
        for (int i = 0; i < spots.size(); i++) {
            out.x[i] = spots.get(i)[0];
            out.y[i] = spots.get(i)[1];
            out.snr[i] = spots.get(i)[2];
            out.prominence[i] = spots.get(i)[3];
        }
        return out;
    }

    private static double[][] computed(Fixture fixture) {
        return smFRETSpotQuality.features(fixture.image, fixture.width, fixture.height,
                fixture.x, fixture.y, fixture.snr, fixture.prominence, fixture.gain,
                fixture.black, fixture.frames, smFRETSpotQuality.shipped());
    }

    @Test
    @DisplayName("the fixture is what it claims to be")
    void theFixtureLoads() {
        Fixture fixture = fixture();
        assertEquals(fixture.width * fixture.height, fixture.image.length, "pixels");
        assertEquals(8, fixture.x.length, "spots");
        assertEquals(8, fixture.features.size(), "feature rows");
        assertEquals(43, fixture.features.get(0).length, "features per spot");
    }

    /**
     * The test this file exists for. Every feature of every spot, named when it disagrees,
     * because "feature 27 is wrong" is a diagnosis and "the features are wrong" is not.
     */
    @Test
    @DisplayName("every feature matches the Python that trained the model")
    void featuresMatchTheTrainer() {
        Fixture fixture = fixture();
        double[][] actual = computed(fixture);
        String[] names = smFRETSpotQuality.shipped().featureNames;

        double worst = 0.0;
        String worstName = "";
        for (int spot = 0; spot < actual.length; spot++) {
            double[] expected = fixture.features.get(spot);
            for (int i = 0; i < expected.length; i++) {
                double scale = Math.max(1.0, Math.abs(expected[i]));
                double error = Math.abs(actual[spot][i] - expected[i]) / scale;
                if (error > worst) {
                    worst = error;
                    worstName = "spot " + spot + " feature " + i
                            + ((i < names.length) ? (" (" + names[i] + ")") : "")
                            + ": " + actual[spot][i] + " against " + expected[i];
                }
            }
        }
        // Not 1e-9. The flux is an integral of a Gaussian over pixels, and Commons
        // Math's erf and the C library's erf are both correct and not identical, so
        // a few parts in 1e8 is the floor here and tightening past it would only
        // pin one implementation's rounding. What that is worth is settled by
        // predictionsMatchTheTrainer below, which asserts on the number the
        // plugin actually acts on.
        assertTrue(worst < 1.0e-7, "worst relative disagreement " + worst + " at " + worstName);
    }

    /**
     * End to end against Python: the features go through the ensemble and the
     * contamination fraction has to come out the same.
     *
     * <p>This is the assertion that matters. A feature can disagree in its eighth digit
     * without anything downstream noticing; a prediction that disagrees is a spot filtered
     * differently. The tolerance here is far below the resolution of any threshold anyone
     * would set.
     */
    @Test
    @DisplayName("predicted contamination matches the Python end to end")
    void predictionsMatchTheTrainer() {
        Fixture fixture = fixture();
        double[][] values = computed(fixture);
        smFRETSpotQuality.Forest forest = smFRETSpotQuality.shipped();

        assertTrue(fixture.predictions != null, "the fixture carries no predictions");
        double worst = 0.0;
        for (int spot = 0; spot < values.length; spot++) {
            worst = Math.max(worst,
                    Math.abs(forest.raw(values[spot]) - fixture.predictions[spot]));
        }
        assertTrue(worst < 1.0e-6, "worst prediction disagreement " + worst);
    }

    /**
     * The two spots against the frame edge are in the fixture on purpose - edge replication
     * is the part of the patch cut most likely to be ported differently, and a spot at the
     * corner exercises both axes at once.
     */
    @Test
    @DisplayName("spots against the frame edge agree too")
    void edgeSpotsAgree() {
        Fixture fixture = fixture();
        double[][] actual = computed(fixture);

        for (int spot : new int[] {6, 7}) {
            double[] expected = fixture.features.get(spot);
            for (int i = 0; i < expected.length; i++) {
                assertEquals(expected[i], actual[spot][i],
                        1.0e-7 * Math.max(1.0, Math.abs(expected[i])),
                        "edge spot " + spot + ", feature " + i);
            }
        }
    }

    /**
     * End to end: features into the ensemble gives a contamination fraction. The close pair
     * in the fixture has to score above the isolated bright spot, or the whole thing is
     * measuring something other than what it claims.
     */
    @Test
    @DisplayName("a close pair scores worse than an isolated spot")
    void thePairIsFlagged() {
        Fixture fixture = fixture();
        double[][] values = computed(fixture);
        smFRETSpotQuality.Forest forest = smFRETSpotQuality.shipped();

        double bright = forest.predict(values[3]);
        double pairLeft = forest.predict(values[1]);
        double pairRight = forest.predict(values[2]);

        assertTrue((pairLeft > bright) && (pairRight > bright),
                "pair " + pairLeft + "/" + pairRight + " against isolated " + bright);
        for (double[] row : values) {
            double predicted = forest.predict(row);
            assertTrue((predicted >= 0.0) && (predicted <= 1.0), "out of range: " + predicted);
        }
    }
}
