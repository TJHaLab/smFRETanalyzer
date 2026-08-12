import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;

/**
 * Predicts how much of a spot's measured signal belongs to something other than the molecule
 * at its centre.
 *
 * <p>A trace is measured by blurring the frame and reading one pixel, so a molecule a few
 * pixels away contributes to a trace that is not its own. The number this predicts is
 *
 * <pre>contamination = 1 - (largest single molecule's contribution) / (total)</pre>
 *
 * which is 0 for a clean isolated molecule, about 0.5 for a pair too close to resolve, and
 * 0.75 and up for a spot sitting on a dye aggregate. It is <i>calibrated</i> rather than
 * merely monotone - 0.2 really does mean about 20% - which is what makes it a number worth
 * showing rather than only a ranking.
 *
 * <p><b>The model is a gradient boosted ensemble trained on simulated data</b>, shipped as
 * {@code spot_quality_model.txt}. Trees rather than a formula because the thing that has to
 * be inferred is the light of molecules the spot finder never detected: given every molecule,
 * a closed form gets this almost exactly right, but given only the reported spots the same
 * formula is worse than useless. Detection efficiency runs from about 90% on a sparse field
 * to 50% on a crowded one, and an unreported molecule contaminates exactly as much as a
 * reported one. Its light is still in the pixels, which is what the model reads. See
 * SIMULATION.md.
 *
 * <p><b>The model is specific to a range of point spread functions</b>, recorded in the
 * resource as {@code trained_sigma_range} and {@code trained_waves_range}. Inside that range
 * it beats every alternative tested; outside it, it is untested and can be worse than nothing.
 * Two attempts to remove that dependence - recomputing every feature at the measured sigma,
 * and handing the model the fitted PSF as an input - both changed nothing, so it is a property
 * of the problem rather than of this implementation.
 *
 * <p>This class holds no ImageJ and no Swing, so it can be scored against fixtures generated
 * by the Python that trained it, the same way {@link smFRETPSF} can.
 */
public class smFRETSpotQuality {

    /** Where the shipped model lives on the classpath. */
    static final String MODEL_RESOURCE = "/spot_quality_model.txt";

    /**
     * A gradient boosted ensemble, flattened into arrays.
     *
     * <p>Every tree's nodes are concatenated into one set of arrays with a per tree offset,
     * rather than kept as objects. A prediction walks 100 short paths, so the difference
     * between chasing pointers and indexing parallel arrays is the whole cost of the thing.
     *
     * <p>Node children are stored as offsets <i>within their own tree</i>, which is how
     * scikit-learn writes them; {@link #predict} adds the tree's base. Getting that wrong
     * gives a model that reads plausible values from the wrong nodes rather than crashing,
     * which is why the fixture test exists.
     */
    public static final class Forest {

        final double baseline;
        final int featureCount;
        final int[] treeStart;
        final boolean[] leaf;
        final int[] feature;
        final double[] threshold;
        final int[] left;
        final int[] right;
        final double[] value;

        /** Metadata carried alongside, so the plugin can report what it is applying. */
        final double defaultThreshold;
        final double analysisSigma;
        final double psfSigma;
        final double psfWaves;
        final double[] sigmaRange;
        final double[] wavesRange;
        final double responseStep;
        final double[] response;
        final String[] featureNames;

        Forest(double baseline, int featureCount, int[] treeStart, boolean[] leaf, int[] feature,
                double[] threshold, int[] left, int[] right, double[] value,
                double defaultThreshold, double analysisSigma, double psfSigma, double psfWaves,
                double[] sigmaRange, double[] wavesRange, double responseStep, double[] response,
                String[] featureNames) {
            this.baseline = baseline;
            this.featureCount = featureCount;
            this.treeStart = treeStart;
            this.leaf = leaf;
            this.feature = feature;
            this.threshold = threshold;
            this.left = left;
            this.right = right;
            this.value = value;
            this.defaultThreshold = defaultThreshold;
            this.analysisSigma = analysisSigma;
            this.psfSigma = psfSigma;
            this.psfWaves = psfWaves;
            this.sigmaRange = sigmaRange;
            this.wavesRange = wavesRange;
            this.responseStep = responseStep;
            this.response = response;
            this.featureNames = featureNames;
        }

        public int treeCount() {
            return treeStart.length;
        }

        public int nodeCount() {
            return leaf.length;
        }

        public int featureCount() {
            return featureCount;
        }

        public double defaultThreshold() {
            return defaultThreshold;
        }

        /**
         * The raw ensemble sum, before clamping. Exposed only so the fixture test can compare
         * against what scikit-learn produces, which is unclamped.
         */
        double raw(double[] features) {
            if (features.length != featureCount) {
                throw new smFRETAnalysisException("the spot quality model wants " + featureCount
                        + " features, was given " + features.length);
            }
            double total = baseline;
            for (int tree = 0; tree < treeStart.length; tree++) {
                int base = treeStart[tree];
                int node = base;
                while (!leaf[node]) {
                    // A NaN feature takes the left branch, matching scikit-learn's default
                    // for a model trained with no missing values. A comparison against NaN
                    // is false either way, so this is what falls out of `<=` regardless -
                    // it is written down because it is a decision, not an accident.
                    node = base + ((features[feature[node]] <= threshold[node])
                            ? left[node] : right[node]);
                }
                total += value[node];
            }
            return total;
        }

        /**
         * The predicted contamination fraction, clamped to [0, 1].
         *
         * <p>The ensemble is a sum of regression trees with nothing bounding it, so on a very
         * clean spot it can return a small negative number - the training predictions run to
         * -0.015. A negative fraction of a spot's light belonging to a neighbour is not a
         * quantity, so it is clamped here rather than shown to anyone.
         */
        public double predict(double[] features) {
            double raw = raw(features);
            return (raw < 0.0) ? 0.0 : ((raw > 1.0) ? 1.0 : raw);
        }

        /**
         * Whether a measured PSF falls inside what the model was trained on.
         *
         * <p>Outside it the prediction is not merely less accurate - at high aberration and a
         * wide core a model trained elsewhere scores below doing nothing at all - so a caller
         * that can measure the PSF should say so rather than quietly scoring anyway.
         */
        public boolean covers(double sigma, double waves) {
            return (sigma >= sigmaRange[0]) && (sigma <= sigmaRange[1])
                    && (waves >= wavesRange[0]) && (waves <= wavesRange[1]);
        }

        /** What a unit intensity neighbour {@code distance} pixels away reads at a spot. */
        public double responseAt(double distance) {
            if (distance <= 0.0) {
                return response[0];
            }
            double position = distance / responseStep;
            int index = (int) position;
            if (index >= (response.length - 1)) {
                return 0.0;
            }
            double fraction = position - index;
            return response[index] + (fraction * (response[index + 1] - response[index]));
        }
    }

    private static Forest shipped;

    /** The model that ships with the plugin, parsed once. */
    public static synchronized Forest shipped() {
        if (shipped == null) {
            InputStream stream = smFRETSpotQuality.class.getResourceAsStream(MODEL_RESOURCE);
            if (stream == null) {
                throw new smFRETAnalysisException("the spot quality model is missing from the "
                        + "jar (" + MODEL_RESOURCE + "); the plugin cannot score spots");
            }
            shipped = load(stream);
        }
        return shipped;
    }

    /**
     * Parse a model file.
     *
     * <p>Plain text rather than a serialized object, so that a model change shows up as a
     * reviewable diff and so that nothing here depends on a library version matching whatever
     * wrote it.
     */
    public static Forest load(InputStream stream) {
        double baseline = 0.0;
        double defaultThreshold = 0.2;
        double analysisSigma = 2.0;
        double psfSigma = 1.5;
        double psfWaves = 0.3;
        double responseStep = 0.25;
        double[] sigmaRange = {1.0, 2.5};
        double[] wavesRange = {0.0, 0.5};
        double[] response = new double[0];
        String[] featureNames = new String[0];
        int featureCount = -1;

        List<Integer> starts = new ArrayList<>();
        List<Boolean> leaf = new ArrayList<>();
        List<Integer> feature = new ArrayList<>();
        List<Double> threshold = new ArrayList<>();
        List<Integer> left = new ArrayList<>();
        List<Integer> right = new ArrayList<>();
        List<Double> value = new ArrayList<>();

        BufferedReader reader = new BufferedReader(
                new InputStreamReader(stream, Charset.forName("UTF-8")));
        try {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) {
                    continue;
                }
                int space = line.indexOf(' ');
                String key = (space < 0) ? line : line.substring(0, space);
                String rest = (space < 0) ? "" : line.substring(space + 1).trim();

                if ("baseline".equals(key)) {
                    baseline = Double.parseDouble(rest);
                } else if ("threshold_default".equals(key)) {
                    defaultThreshold = Double.parseDouble(rest);
                } else if ("analysis_sigma".equals(key)) {
                    analysisSigma = Double.parseDouble(rest);
                } else if ("psf_sigma".equals(key)) {
                    psfSigma = Double.parseDouble(rest);
                } else if ("psf_waves".equals(key)) {
                    psfWaves = Double.parseDouble(rest);
                } else if ("trained_sigma_range".equals(key)) {
                    sigmaRange = parseDoubles(rest.split(" +"));
                } else if ("trained_waves_range".equals(key)) {
                    wavesRange = parseDoubles(rest.split(" +"));
                } else if ("features".equals(key)) {
                    featureCount = Integer.parseInt(rest);
                } else if ("feature_names".equals(key)) {
                    featureNames = rest.split(",");
                } else if ("response_step".equals(key)) {
                    responseStep = Double.parseDouble(rest);
                } else if ("response".equals(key)) {
                    response = parseDoubles(rest.split(" +"));
                } else if ("tree".equals(key)) {
                    String[] parts = rest.split(" +");
                    int nodes = Integer.parseInt(parts[1]);
                    starts.add(leaf.size());
                    for (int i = 0; i < nodes; i++) {
                        String node = reader.readLine();
                        if (node == null) {
                            throw new smFRETAnalysisException(
                                    "the spot quality model ends part way through a tree");
                        }
                        String[] fields = node.trim().split(" +");
                        leaf.add(!"0".equals(fields[0]));
                        feature.add(Integer.parseInt(fields[1]));
                        threshold.add(Double.parseDouble(fields[2]));
                        left.add(Integer.parseInt(fields[3]));
                        right.add(Integer.parseInt(fields[4]));
                        value.add(Double.parseDouble(fields[5]));
                    }
                }
            }
        } catch (IOException e) {
            throw new smFRETAnalysisException("could not read the spot quality model: "
                    + e.getMessage());
        } finally {
            try {
                reader.close();
            } catch (IOException ignored) {
                // Nothing useful to do about a failure to close a resource stream.
            }
        }

        if ((featureCount <= 0) || starts.isEmpty()) {
            throw new smFRETAnalysisException(
                    "the spot quality model has no trees or no feature count");
        }

        return new Forest(baseline, featureCount, toIntArray(starts), toBooleanArray(leaf),
                toIntArray(feature), toDoubleArray(threshold), toIntArray(left),
                toIntArray(right), toDoubleArray(value), defaultThreshold, analysisSigma,
                psfSigma, psfWaves, sigmaRange, wavesRange, responseStep, response,
                featureNames);
    }

    private static double[] parseDoubles(String[] parts) {
        double[] out = new double[parts.length];
        for (int i = 0; i < parts.length; i++) {
            out[i] = Double.parseDouble(parts[i]);
        }
        return out;
    }

    private static int[] toIntArray(List<Integer> values) {
        int[] out = new int[values.size()];
        for (int i = 0; i < out.length; i++) {
            out[i] = values.get(i);
        }
        return out;
    }

    private static double[] toDoubleArray(List<Double> values) {
        double[] out = new double[values.size()];
        for (int i = 0; i < out.length; i++) {
            out[i] = values.get(i);
        }
        return out;
    }

    private static boolean[] toBooleanArray(List<Boolean> values) {
        boolean[] out = new boolean[values.size()];
        for (int i = 0; i < out.length; i++) {
            out[i] = values.get(i);
        }
        return out;
    }
}
