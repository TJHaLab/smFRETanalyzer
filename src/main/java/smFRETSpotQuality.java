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
     * Half width of the patch cut around each spot, in pixels.
     *
     * <p>Not a free choice: the model was trained on 25x25 patches and every radial feature
     * is an annulus of one of them, so changing this does not widen the model's view, it
     * feeds it different numbers under the same names.
     */
    static final int PATCH_RADIUS = 12;

    /** Where the Gaussian measurement kernel is truncated, matching the trainer. */
    static final double KERNEL_EXTENT_SIGMA = 4.0;

    /** Neighbours further than this contribute nothing worth summing. */
    static final double NEIGHBOUR_REACH = 45.0;

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

    /**
     * The 43 features the model wants, one row per spot.
     *
     * <p>Everything comes from the analysis image and the detection list, which is all the
     * spot finder has when it would apply a score. The two absolute quantities - the local
     * background and the spot's flux - are converted to <b>photoelectrons per frame</b> using
     * the gain, black level and frame count, because the model was trained on one camera and
     * an ADU on another camera is a different amount of light. The frame count goes in as a
     * feature too, since it sets how noisy the averaged image is.
     *
     * @param image the frame-averaged analysis image, row major
     * @param x     detection x coordinates, in pixels, as MaximumFinder reported them
     * @param frames how many frames went into the average
     */
    public static double[][] features(float[] image, int width, int height, double[] x,
            double[] y, double[] snr, double[] prominence, double gain, double blackLevel,
            int frames, Forest forest) {

        int count = x.length;
        int size = (2 * PATCH_RADIUS) + 1;
        double[][] out = new double[count][forest.featureCount];

        double norm = 2.0 * Math.PI * forest.analysisSigma * forest.analysisSigma;
        double own = forest.responseAt(0.0);
        double density = (count / (double) (width * height)) * 1000.0;
        double logFrames = Math.log10(frames);

        double[][] patch = new double[size][size];
        double[] measured = new double[count];
        double[] background = new double[count];

        // Two passes: every spot's own flux has to exist before any spot's neighbour sum can
        // be worked out, since a neighbour's contribution is scaled by what it measures.
        for (int spot = 0; spot < count; spot++) {
            cutPatch(image, width, height, (int) x[spot], (int) y[spot], patch);
            background[spot] = borderMedian(patch);
            measured[spot] = matchedFilter(image, width, height, x[spot], y[spot],
                    forest.analysisSigma) - (norm * background[spot]);
        }

        for (int spot = 0; spot < count; spot++) {
            cutPatch(image, width, height, (int) x[spot], (int) y[spot], patch);
            double[] row = out[spot];

            double nearest = 999.0;
            double within5 = 0.0;
            double within10 = 0.0;
            double within15 = 0.0;
            double overlap = 0.0;
            for (int other = 0; other < count; other++) {
                if (other == spot) {
                    continue;
                }
                double gap = Math.hypot(x[other] - x[spot], y[other] - y[spot]);
                if (gap < nearest) {
                    nearest = gap;
                }
                if (gap < 5.0) {
                    within5 += 1.0;
                }
                if (gap < 10.0) {
                    within10 += 1.0;
                }
                if (gap < 15.0) {
                    within15 += 1.0;
                }
                if (gap < NEIGHBOUR_REACH) {
                    overlap += forest.responseAt(gap) * (measured[other] / own);
                }
            }
            overlap *= gain;

            double backgroundElectrons = (background[spot] - (2.0 * blackLevel)) * gain;
            double flux = Math.max(measured[spot] * gain, 1.0);

            double[] shape = shapeStatistics(patch, background[spot], forest.psfSigma);

            row[0] = prominence[spot];
            row[1] = snr[spot];
            row[2] = Math.log10(flux);
            row[3] = backgroundElectrons;
            row[4] = overlap;
            row[5] = overlap / flux;
            row[6] = nearest;
            row[7] = within5;
            row[8] = within10;
            row[9] = within15;
            row[10] = density;
            row[11] = shape[0];
            row[12] = shape[1];
            row[13] = shape[2];
            row[14] = shape[0] * snr[spot];
            row[15] = logFrames;

            radial(patch, background[spot], row, 16);
        }
        return out;
    }

    /**
     * A square of the image around a pixel, with the frame edge replicated outwards.
     *
     * <p>Replicated rather than zero filled, and rather than dropping the spot: a spot near
     * the edge is a real molecule and zeros would read as a very dark neighbourhood, which is
     * a shape the model has never seen. This matches what the trainer did.
     */
    static void cutPatch(float[] image, int width, int height, int cx, int cy,
            double[][] patch) {
        int size = patch.length;
        for (int dy = 0; dy < size; dy++) {
            int row = clamp((cy + dy) - PATCH_RADIUS, 0, height - 1);
            for (int dx = 0; dx < size; dx++) {
                int column = clamp((cx + dx) - PATCH_RADIUS, 0, width - 1);
                patch[dy][dx] = image[(row * width) + column];
            }
        }
    }

    private static int clamp(int value, int low, int high) {
        return (value < low) ? low : ((value > high) ? high : value);
    }

    /**
     * Median of the patch's outermost ring: the background estimate the model was trained on.
     *
     * <p>Not the spot finder's own background estimate, good though that is - the model was
     * fitted against this one, and the two do not agree. The ring reads slightly high because
     * the PSF's wings sit on it, and the model has learned that.
     */
    static double borderMedian(double[][] patch) {
        int size = patch.length;
        double[] ring = new double[(4 * size) - 4];
        int at = 0;
        for (int i = 0; i < size; i++) {
            ring[at++] = patch[0][i];
        }
        for (int i = 0; i < size; i++) {
            ring[at++] = patch[size - 1][i];
        }
        for (int i = 1; i < (size - 1); i++) {
            ring[at++] = patch[i][0];
        }
        for (int i = 1; i < (size - 1); i++) {
            ring[at++] = patch[i][size - 1];
        }
        java.util.Arrays.sort(ring);
        int middle = ring.length / 2;
        return ((ring.length % 2) == 0)
                ? (0.5 * (ring[middle - 1] + ring[middle])) : ring[middle];
    }

    /**
     * The spot's flux, measured exactly the way smFRETAnalyzer measures a trace: blur with a
     * Gaussian of the analysis sigma, read the value at the spot's pixel, scale by 2 pi s^2.
     *
     * <p>The kernel is <i>integrated</i> over each pixel rather than sampled at its centre,
     * and is centred on the pixel's centre rather than on the reported coordinate, both of
     * which match the trainer. Sampling instead of integrating would bias this by several
     * percent at these sigmas.
     */
    static double matchedFilter(float[] image, int width, int height, double x, double y,
            double sigma) {
        double centreX = Math.floor(x) + 0.5;
        double centreY = Math.floor(y) + 0.5;
        double extent = KERNEL_EXTENT_SIGMA * sigma;

        int x0 = Math.max(0, (int) Math.floor(centreX - extent));
        int x1 = Math.min(width, (int) Math.ceil(centreX + extent) + 1);
        int y0 = Math.max(0, (int) Math.floor(centreY - extent));
        int y1 = Math.min(height, (int) Math.ceil(centreY + extent) + 1);
        if ((x0 >= x1) || (y0 >= y1)) {
            return 0.0;
        }

        double[] weightX = pixelWeights(centreX, sigma, x0, x1);
        double[] weightY = pixelWeights(centreY, sigma, y0, y1);

        double total = 0.0;
        for (int j = 0; j < weightY.length; j++) {
            for (int i = 0; i < weightX.length; i++) {
                total += weightY[j] * weightX[i];
            }
        }
        if (total <= 0.0) {
            return 0.0;
        }

        double sum = 0.0;
        for (int j = 0; j < weightY.length; j++) {
            int offset = (y0 + j) * width;
            for (int i = 0; i < weightX.length; i++) {
                sum += image[offset + x0 + i] * ((weightY[j] * weightX[i]) / total);
            }
        }
        return 2.0 * Math.PI * sigma * sigma * sum;
    }

    /** Integral of a unit normal over each pixel in [low, high). */
    private static double[] pixelWeights(double centre, double sigma, int low, int high) {
        double[] weights = new double[high - low];
        double previous = normalCdf((low - centre) / sigma);
        for (int i = 0; i < weights.length; i++) {
            double next = normalCdf(((low + i + 1) - centre) / sigma);
            weights[i] = next - previous;
            previous = next;
        }
        return weights;
    }

    private static double normalCdf(double value) {
        return 0.5 * (1.0 + org.apache.commons.math3.special.Erf.erf(value / Math.sqrt(2.0)));
    }

    /**
     * Ellipticity, moment size and the reduced residual of a single Gaussian fit.
     *
     * <p>All three at the PSF sigma the model was trained at rather than at whatever this
     * data's PSF is - deliberately, since recomputing them at the measured sigma was tried
     * and changed nothing.
     *
     * <p>The Gaussian fit is a closed form: with the width known and the centroid taken from
     * the moments, amplitude and a local offset are a 2x2 normal equation rather than an
     * iteration that could fail to converge.
     */
    static double[] shapeStatistics(double[][] patch, double background, double sigma) {
        int size = patch.length;
        int radius = size / 2;
        double twoSigmaSquared = 2.0 * sigma * sigma;

        double total = 0.0;
        double sumX = 0.0;
        double sumY = 0.0;
        for (int j = 0; j < size; j++) {
            double dy = j - radius;
            for (int i = 0; i < size; i++) {
                double dx = i - radius;
                double signal = patch[j][i] - background;
                double weighted = ((signal > 0.0) ? signal : 0.0)
                        * Math.exp(-(((dx * dx) + (dy * dy)) / twoSigmaSquared));
                total += weighted;
                sumX += weighted * dx;
                sumY += weighted * dy;
            }
        }
        double weight = (total > 0.0) ? total : 1.0;
        double centroidX = sumX / weight;
        double centroidY = sumY / weight;

        double mxx = 0.0;
        double myy = 0.0;
        double mxy = 0.0;
        for (int j = 0; j < size; j++) {
            double dy = j - radius;
            for (int i = 0; i < size; i++) {
                double dx = i - radius;
                double signal = patch[j][i] - background;
                double weighted = ((signal > 0.0) ? signal : 0.0)
                        * Math.exp(-(((dx * dx) + (dy * dy)) / twoSigmaSquared));
                double ox = dx - centroidX;
                double oy = dy - centroidY;
                mxx += weighted * ox * ox;
                myy += weighted * oy * oy;
                mxy += weighted * ox * oy;
            }
        }
        mxx /= weight;
        myy /= weight;
        mxy /= weight;

        double momentSize = mxx + myy;
        double ellipticity = 0.0;
        if (momentSize > 0.0) {
            double e1 = (mxx - myy) / momentSize;
            double e2 = (2.0 * mxy) / momentSize;
            ellipticity = Math.hypot(e1, e2);
        } else {
            momentSize = 0.0;
        }

        double sumModelModel = 0.0;
        double sumModel = 0.0;
        double sumModelSignal = 0.0;
        double sumSignal = 0.0;
        int pixels = size * size;
        for (int j = 0; j < size; j++) {
            double dy = (j - radius) - centroidY;
            for (int i = 0; i < size; i++) {
                double dx = (i - radius) - centroidX;
                double model = Math.exp(-(((dx * dx) + (dy * dy)) / twoSigmaSquared));
                double signal = patch[j][i] - background;
                sumModelModel += model * model;
                sumModel += model;
                sumModelSignal += model * signal;
                sumSignal += signal;
            }
        }
        double determinant = (sumModelModel * pixels) - (sumModel * sumModel);
        if (Math.abs(determinant) <= 1.0e-12) {
            determinant = 1.0;
        }
        double amplitude = ((pixels * sumModelSignal) - (sumModel * sumSignal)) / determinant;
        double offset = ((sumModelModel * sumSignal) - (sumModel * sumModelSignal))
                / determinant;

        double chiSquared = 0.0;
        for (int j = 0; j < size; j++) {
            double dy = (j - radius) - centroidY;
            for (int i = 0; i < size; i++) {
                double dx = (i - radius) - centroidX;
                double model = Math.exp(-(((dx * dx) + (dy * dy)) / twoSigmaSquared));
                double residual = (patch[j][i] - background) - (amplitude * model) - offset;
                double variance = Math.max(patch[j][i], 1.0);
                chiSquared += (residual * residual) / variance;
            }
        }
        return new double[] {ellipticity, momentSize, chiSquared / (pixels - 2)};
    }

    /**
     * The radial block: mean and max in each one pixel annulus, then the log of the centre.
     *
     * <p>The patch is divided by its own centre pixel first, so the shape is scale free and
     * the brightness goes in separately. Annuli rather than pixels because contamination does
     * not care which side the neighbour is on - a representation that has to learn that
     * spends its capacity on it. Measured: 23 radial numbers beat 24 principal components of
     * the raw patch, 0.903 against 0.735.
     */
    static void radial(double[][] patch, double background, double[] row, int at) {
        int size = patch.length;
        int radius = size / 2;
        double centre = patch[radius][radius] - background;
        double scale = (Math.abs(centre) > 1.0e-6) ? centre : 1.0;

        for (int ring = 0; ring <= radius; ring++) {
            double sum = 0.0;
            double largest = Double.NEGATIVE_INFINITY;
            int seen = 0;
            for (int j = 0; j < size; j++) {
                double dy = j - radius;
                for (int i = 0; i < size; i++) {
                    double dx = i - radius;
                    double distance = Math.hypot(dx, dy);
                    if ((distance >= (ring - 0.5)) && (distance < (ring + 0.5))) {
                        double value = (patch[j][i] - background) / scale;
                        sum += value;
                        if (value > largest) {
                            largest = value;
                        }
                        seen++;
                    }
                }
            }
            row[at + ring] = (seen > 0) ? (sum / seen) : 0.0;
            row[at + radius + 1 + ring] = (seen > 0) ? largest : 0.0;
        }
        row[at + (2 * (radius + 1))] = Math.log10(Math.max(centre, 1.0));
    }
}
