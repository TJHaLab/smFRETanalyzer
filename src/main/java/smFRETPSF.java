/*
 * Measuring the point spread function from the spots a field already contains, and fitting an
 * aberrated Airy pattern to it.
 *
 * This is a port of scripts/fit_psf.py in the smfret-simulator repository, which is where the
 * method was worked out and checked against simulated data of known PSF. The numbers it produced
 * on the Ha lab example data - sigma 1.36 px, 0.41 waves of spherical aberration, over a flat
 * pedestal of 3.3% of peak - are what the `ha-lab` simulation preset uses, so keeping this
 * faithful to it is what lets a measurement here be compared against a simulation there.
 *
 * Why an Airy pattern with spherical aberration rather than a Gaussian: real spots on this data
 * carry several percent of their peak out at 5 to 10 pixels, where a Gaussian of the same core
 * width has fallen to 1e-5. Those wings are what a background estimator has to cope with and what
 * makes neighbouring molecules contaminate each other, so they are not a detail.
 *
 * No Swing and no ImageJ here, deliberately: everything in this class takes a float[] and returns
 * numbers, so the measurement can be exercised against a field of known PSF without a display.
 */

import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.fitting.leastsquares.ParameterValidator;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.special.BesselJ;
import org.apache.commons.math3.util.Pair;

import java.util.ArrayList;
import java.util.List;

final class smFRETPSF {

    private smFRETPSF() {
    }

    /**
     * First zero of the unaberrated Airy pattern, in units of the equivalent Gaussian sigma.
     *
     * The usual approximation is sigma = 0.21 lambda/NA against a first zero at 0.61 lambda/NA.
     * Sizing the Airy pattern by a sigma rather than by its own zero is what keeps this
     * comparable with spotSigma and with the Gaussian the rest of the pipeline assumes.
     */
    static final double AIRY_ZERO_SIGMAS = 0.61 / 0.21;

    /** Argument of the first zero of J1, which is where the unaberrated pattern vanishes. */
    private static final double AIRY_FIRST_ZERO_V = 3.8317;

    /**
     * The settings the Python arrived at, and what they are trading off.
     *
     * The patch has to be wide enough to hold the wings the model exists to measure - an
     * aberrated Airy is still at a percent of its peak at ten pixels - and narrow enough that a
     * crowded field can still find border pixels for it. The mask is four sigma of the core at
     * the fitted size; radii of 5 and 6 recovered the simulated truth equally well and 6 had the
     * lower residual.
     *
     * Note the two together: a neighbour further than patch + mask cannot reach any pixel of the
     * patch, so at these defaults nothing beyond 16 pixels is masked at all.
     */
    static final int DEFAULT_BINS = 20;
    static final double DEFAULT_NEIGHBOUR_MASK = 6.0;
    static final int DEFAULT_PATCH = 10;

    // Pupil samples in the aperture integral. Kept at the 2000 the Python used rather than
    // reduced: with the J0 table below each sample is two array reads, so the whole integral is
    // microseconds and there is nothing to buy by approximating it further.
    private static final int APERTURE_SAMPLES = 2000;

    // Patch border pixels a spot must have left after neighbour masking, as a fraction of the
    // whole border. Below about a quarter the median is being taken over one side of the patch
    // and has stopped being a local background.
    private static final double MIN_BORDER_FRACTION = 0.25;

    // The border correction is solved rather than assumed, so it iterates. Four passes is what
    // the simulated data needed; the budget is generous because each pass is one more least
    // squares fit over twenty points and the whole thing runs in milliseconds.
    private static final int BORDER_ITERATIONS = 8;
    private static final double BORDER_TOLERANCE = 1.0e-5;

    // Fit bounds, from the Python. Guardrails rather than working constraints - the optimum on
    // both real and simulated data sits well inside them.
    private static final double[] FIT_LOWER = {0.5, 0.0, 0.5, 0.0};
    private static final double[] FIT_UPPER = {6.0, 3.0, 2.0, 0.2};
    private static final double[] FIT_START = {2.0, 0.3, 1.0, 0.03};

    /**
     * J0 on a lookup table, because the fit asks for millions of them.
     *
     * The aperture integral evaluates J0 once per pupil sample per radius, and a fit runs the
     * model a few hundred times - which through a general order Bessel routine is seconds per
     * channel and not something an interactive control can wait for. Linear interpolation on this
     * grid is accurate to about 1e-7, far past what a fit to twenty measured points can use.
     *
     * The range covers the largest argument reachable: v = r * 3.8317 / (AIRY_ZERO_SIGMAS *
     * sigma), which at the smallest fitted sigma of 0.5 and a patch of 20 pixels is 53.
     */
    private static final double J0_MAX = 128.0;
    private static final int J0_SAMPLES = 1 << 17;
    private static final double J0_STEP = J0_MAX / (J0_SAMPLES - 1);
    private static final double[] J0_TABLE = buildJ0Table();

    private static double[] buildJ0Table() {
        double[] table = new double[J0_SAMPLES];
        for (int i = 0; i < J0_SAMPLES; i++) {
            table[i] = BesselJ.value(0.0, i * J0_STEP);
        }
        return table;
    }

    private static double besselJ0(double x) {
        double a = Math.abs(x);
        if (a >= J0_MAX) {

            // Unreachable through airyProfile at any fitted sigma, and going to the real routine
            // rather than clamping means a caller outside those bounds gets an answer, not a lie.
            return BesselJ.value(0.0, a);
        }
        double position = a / J0_STEP;
        int index = (int) position;
        double fraction = position - index;
        return J0_TABLE[index] + fraction * (J0_TABLE[index + 1] - J0_TABLE[index]);
    }

    /**
     * The aberrated Airy shape, peak normalised.
     *
     * The diffraction pattern of a circular pupil carrying primary spherical aberration, the
     * wavefront error W(rho) = A rho^4:
     *
     *   I(r) proportional to | 2 integral_0^1 exp(2 pi i A rho^4) J0(v rho) rho drho |^2
     *
     * with v = 3.8317 r / r_zero, so v = 3.8317 - the first zero of the unaberrated pattern -
     * lands at r_zero. There is no closed form, hence the numerical integral.
     *
     * Aberration fills the wings by taking light *out of* the core, so it cannot thicken one
     * without widening the other - which is why sigma and waves have to be fitted together and
     * why fitting the wings alone gives a core visibly broader than the real spots.
     */
    static double[] airyProfile(double[] radii, double sigma, double waves) {
        double firstZero = AIRY_ZERO_SIGMAS * sigma;

        // The pupil sampling and the aberration phase depend only on waves, so they are built
        // once per call rather than once per radius. Midpoints, so the sum is a midpoint rule
        // rather than a trapezoid missing its endpoints.
        double[] rho = new double[APERTURE_SAMPLES];
        double[] weightCos = new double[APERTURE_SAMPLES];
        double[] weightSin = new double[APERTURE_SAMPLES];
        for (int i = 0; i < APERTURE_SAMPLES; i++) {
            double r = (i + 0.5) / APERTURE_SAMPLES;
            double phase = 2.0 * Math.PI * waves * r * r * r * r;
            rho[i] = r;
            weightCos[i] = r * Math.cos(phase);
            weightSin[i] = r * Math.sin(phase);
        }

        double peak = intensity(0.0, rho, weightCos, weightSin);
        double[] out = new double[radii.length];
        for (int i = 0; i < radii.length; i++) {
            double v = radii[i] * AIRY_FIRST_ZERO_V / firstZero;
            out[i] = intensity(v, rho, weightCos, weightSin) / peak;
        }
        return out;
    }

    /** The squared modulus of the pupil integral at one radius. */
    private static double intensity(double v, double[] rho, double[] weightCos, double[] weightSin) {
        double real = 0.0;
        double imaginary = 0.0;
        for (int i = 0; i < rho.length; i++) {
            double kernel = besselJ0(v * rho[i]);
            real += kernel * weightCos[i];
            imaginary += kernel * weightSin[i];
        }
        return real * real + imaginary * imaginary;
    }

    /**
     * A fitted profile: the shape, how tall it is, and what flat term it sits on.
     */
    static final class Fit {
        final double amplitude;
        final double pedestal;
        final double rms;
        final double sigma;
        final double waves;

        Fit(double sigma, double waves, double amplitude, double pedestal, double rms) {
            this.sigma = sigma;
            this.waves = waves;
            this.amplitude = amplitude;
            this.pedestal = pedestal;
            this.rms = rms;
        }

        /** Where the first zero of the unaberrated pattern of this size would fall, in pixels. */
        double firstZero() {
            return sigma * AIRY_ZERO_SIGMAS;
        }

        /** The fitted curve at a set of radii, pedestal included. */
        double[] at(double[] radii) {
            double[] shape = airyProfile(radii, sigma, waves);
            double[] out = new double[radii.length];
            for (int i = 0; i < radii.length; i++) {
                out[i] = amplitude * shape[i] + pedestal;
            }
            return out;
        }
    }

    /**
     * Every uncontaminated pixel of every usable spot, pooled, plus the same thing kept as an
     * image.
     *
     * Separate from the binning and the fit because it is the expensive half and depends on
     * different settings: the patch and the neighbour mask change it, the bin count and the
     * pedestal do not.
     */
    static final class Samples {
        final int[] imageCount;         // How many spots contributed to each pixel.
        final double[] imagePeak;       // Sum of the contributing spots' peaks, per pixel.
        final double[] imageValue;      // Sum of their background subtracted values, per pixel.
        final int invalidPixels;        // Dropped for falling outside the mapped region.
        final int patch;
        final double[] peak;            // The spot peak behind each pooled sample.
        final double[] radius;
        final int size;                 // 2 * patch + 1.
        final int spotsTotal;
        final int spotsUsed;
        final double[] value;

        Samples(int patch, double[] radius, double[] value, double[] peak,
                double[] imageValue, double[] imagePeak, int[] imageCount,
                int spotsUsed, int spotsTotal, int invalidPixels) {
            this.invalidPixels = invalidPixels;
            this.patch = patch;
            this.size = 2 * patch + 1;
            this.radius = radius;
            this.value = value;
            this.peak = peak;
            this.imageValue = imageValue;
            this.imagePeak = imagePeak;
            this.imageCount = imageCount;
            this.spotsUsed = spotsUsed;
            this.spotsTotal = spotsTotal;
        }

        /**
         * The measured PSF as an image, peak normalised, with NaN where nothing contributed.
         *
         * A ratio of sums per pixel, the same way the radial profile is pooled, so a faint spot
         * whose centre pixel happens to read three ADU cannot contribute a patch three times too
         * tall. NaN rather than zero for an empty pixel so the display can tell "no data" from
         * "no light" - at a large neighbour mask on a crowded field the two are very different.
         */
        double[] image() {
            double[] out = new double[imageValue.length];
            for (int i = 0; i < out.length; i++) {
                out[i] = (imageCount[i] > 0) ? (imageValue[i] / imagePeak[i]) : Double.NaN;
            }
            return out;
        }
    }

    /**
     * A measured radial profile, what it took to get there, and the fits to it.
     */
    static final class Measurement {
        final double[] binCentre;       // Mean radius of what landed in each bin.
        final int[] binCount;
        final double[] binProfile;      // Border corrected.
        final double[] binRaw;          // Before the border correction, for showing the difference.
        final double borderLevel;
        final boolean borderSettled;
        final Fit psfOnly;
        final Samples samples;
        final Fit withPedestal;

        Measurement(Samples samples, double[] binCentre, double[] binRaw, double[] binProfile,
                    int[] binCount, double borderLevel, boolean borderSettled,
                    Fit psfOnly, Fit withPedestal) {
            this.samples = samples;
            this.binCentre = binCentre;
            this.binRaw = binRaw;
            this.binProfile = binProfile;
            this.binCount = binCount;
            this.borderLevel = borderLevel;
            this.borderSettled = borderSettled;
            this.psfOnly = psfOnly;
            this.withPedestal = withPedestal;
        }

        /**
         * The measured PSF as an image, on the same footing as the radial profile.
         *
         * Samples.image() is the raw pooled ratio, with the border median still subtracted out of
         * it - so its outer pixels read near zero by construction and some go slightly negative.
         * The profile has the border put back before it is fitted or drawn, and the image has to
         * have it too, or the two halves of the display are showing the same measurement on
         * different footings and the image looks black where the profile says there is light.
         */
        double[] correctedImage() {
            double[] image = samples.image();
            for (int i = 0; i < image.length; i++) {
                if (!Double.isNaN(image[i])) {
                    image[i] = (image[i] * (1.0 - borderLevel)) + borderLevel;
                }
            }
            return image;
        }

        /** Bins holding enough pixels to mean anything, which is what the fit ran over. */
        boolean[] usable() {
            boolean[] good = new boolean[binCount.length];
            for (int i = 0; i < good.length; i++) {
                good[i] = binCount[i] > MIN_BIN_SAMPLES;
            }
            return good;
        }
    }

    // A bin holding fewer than this is one or two pixels of one or two spots, and its value is
    // noise rather than a measurement.
    private static final int MIN_BIN_SAMPLES = 50;

    /**
     * Cut a patch around every spot, drop the pixels a neighbour reaches into, and pool the rest.
     *
     * **Contaminated pixels are dropped, not contaminated spots**, and that is the whole reason
     * this works on a real field. Rejecting every spot with a close neighbour is the obvious
     * approach and it fails twice over: a FRET field is crowded, so on simulated data it throws
     * away 135 of 235 spots, and the 100 that survive still have neighbour light in their wings -
     * fitting those returns 1.53 sigma and 0.36 waves where the truth is 1.45 and 0.41. Masking
     * the pixels instead keeps 234 spots and returns 1.49 and 0.41.
     *
     * **Pixels outside the mapped region are dropped along with the contaminated ones.** Warping
     * the acceptor onto the donor's frame leaves a band around the edges where the inverse mapped
     * coordinate falls outside the source, and transformImagePlus writes an exact zero there -
     * the same sentinel createOverlapMask reads to find the region the two channels share. Those
     * are not measurements of anything, and a patch reaching into that band would have them in
     * its border median, in its peak and in the pool. On the example data the band is 1.5% of the
     * acceptor half and three of its patches touch it; on a movie with a larger channel offset it
     * would be far more. The donor is unwarped and has none, so this costs it nothing.
     *
     * Patches are centred on the integer pixel the spot finder reported. MaximumFinder returns
     * the brightest pixel rather than a centroid, so a spot sits up to half a pixel off centre
     * and the pooled core comes out slightly broader than the truth - 1.49 against 1.45 on the
     * simulated check. That is the method as validated; sub-pixel centring would tighten it but
     * is a different method and would need checking against known truth before being believed.
     */
    static Samples extract(float[] pixels, int width, int height, double[][] spots,
                           int patch, double neighbourMask) {
        int size = 2 * patch + 1;

        // A neighbour further than this cannot reach any pixel of the patch.
        double reach = Math.hypot(patch, patch) + neighbourMask;

        double[] patchRadius = new double[size * size];
        boolean[] edge = new boolean[size * size];
        int borderPixels = 0;
        for (int dy = 0; dy < size; dy++) {
            for (int dx = 0; dx < size; dx++) {
                int i = dy * size + dx;
                patchRadius[i] = Math.hypot(dx - patch, dy - patch);
                edge[i] = (dx == 0) || (dy == 0) || (dx == (size - 1)) || (dy == (size - 1));
                if (edge[i]) {
                    borderPixels++;
                }
            }
        }
        int minBorder = Math.max(8, (int) Math.round(MIN_BORDER_FRACTION * borderPixels));

        List<double[]> pooled = new ArrayList<>();
        double[] imageValue = new double[size * size];
        double[] imagePeak = new double[size * size];
        int[] imageCount = new int[size * size];
        boolean[] keep = new boolean[size * size];
        float[] border = new float[size * size];
        int invalidPixels = 0;
        int used = 0;

        for (double[] spot : spots) {
            int cx = (int) spot[0];
            int cy = (int) spot[1];

            // Two pixels of slack past the patch, matching the Python: a patch running to the
            // very edge of the frame has a border made of frame edge rather than of field.
            if ((cx < (patch + 2)) || (cx >= (width - patch - 2))
                    || (cy < (patch + 2)) || (cy >= (height - patch - 2))) {
                continue;
            }

            // Out of range pixels start excluded rather than being cleared afterwards, so they
            // are absent from the border median and the peak test as well as from the pool.
            int dropped = 0;
            for (int dy = 0; dy < size; dy++) {
                for (int dx = 0; dx < size; dx++) {
                    boolean mapped =
                            pixels[(cy - patch + dy) * width + (cx - patch + dx)] != 0.0f;
                    keep[dy * size + dx] = mapped;
                    if (!mapped) {
                        dropped++;
                    }
                }
            }

            for (double[] other : spots) {
                double gapX = other[0] - spot[0];
                double gapY = other[1] - spot[1];
                double gap = Math.hypot(gapX, gapY);
                if ((gap <= 1.0e-6) || (gap >= reach)) {
                    continue;
                }
                for (int dy = 0; dy < size; dy++) {
                    for (int dx = 0; dx < size; dx++) {
                        double px = cx - patch + dx;
                        double py = cy - patch + dy;
                        if (Math.hypot(px - other[0], py - other[1]) < neighbourMask) {
                            keep[dy * size + dx] = false;
                        }
                    }
                }
            }

            // A neighbour close enough to mask this spot's own centre means the two are not
            // separable, and the peak here is not one spot's.
            if (!keep[patch * size + patch]) {
                continue;
            }

            int borderCount = 0;
            for (int i = 0; i < keep.length; i++) {
                if (edge[i] && keep[i]) {
                    int dx = i % size;
                    int dy = i / size;
                    border[borderCount++] = pixels[(cy - patch + dy) * width + (cx - patch + dx)];
                }
            }
            if (borderCount < minBorder) {
                continue;
            }

            // The same quickselect the background estimator uses, for the same reason: the middle
            // element does not require ordering the rest.
            double level = smFRETSpotFinder.median(border, borderCount);
            double peak = pixels[cy * width + cx] - level;
            if (peak <= 0.0) {
                continue;
            }

            invalidPixels += dropped;
            for (int i = 0; i < keep.length; i++) {
                if (!keep[i]) {
                    continue;
                }
                int dx = i % size;
                int dy = i / size;
                double value = pixels[(cy - patch + dy) * width + (cx - patch + dx)] - level;
                pooled.add(new double[] {patchRadius[i], value, peak});
                imageValue[i] += value;
                imagePeak[i] += peak;
                imageCount[i]++;
            }
            used++;
        }

        double[] radius = new double[pooled.size()];
        double[] value = new double[pooled.size()];
        double[] peak = new double[pooled.size()];
        for (int i = 0; i < pooled.size(); i++) {
            radius[i] = pooled.get(i)[0];
            value[i] = pooled.get(i)[1];
            peak[i] = pooled.get(i)[2];
        }
        return new Samples(patch, radius, value, peak, imageValue, imagePeak, imageCount,
                used, spots.length, invalidPixels);
    }

    /**
     * Bin the pooled samples, correct for the border, and fit.
     *
     * **Each bin is a ratio of sums, not a mean of ratios.** Dividing each patch by its own centre
     * pixel divides by a small noisy number - spots in the real 8 bit data rise a median of 9 ADU
     * above background - and whatever the border median failed to remove is magnified into the
     * wings as a floor. Summing first weights every spot by its brightness and the noise averages
     * out instead. On the real data the difference is a fitted pedestal of 10% of peak against 3%.
     *
     * A bin is reported at the mean radius of what fell into it rather than at the middle of its
     * interval. Pixel radii are hypot of two integers and take a few discrete values which need
     * not sit near the middle: the first bin holds only the centre pixel, at radius 0 rather than
     * the 0.25 its interval suggests. Evaluating the model at interval centres instead doubles the
     * residual and costs about 0.01 waves.
     */
    static Measurement analyse(Samples samples, int bins) {
        double[] centre = new double[bins];
        double[] raw = new double[bins];
        int[] count = new int[bins];
        double[] valueSum = new double[bins];
        double[] peakSum = new double[bins];
        double[] radiusSum = new double[bins];

        double width = (double) samples.patch / bins;
        for (int i = 0; i < samples.radius.length; i++) {
            int bin = (int) (samples.radius[i] / width);
            if ((bin < 0) || (bin >= bins)) {
                continue;
            }
            valueSum[bin] += samples.value[i];
            peakSum[bin] += samples.peak[i];
            radiusSum[bin] += samples.radius[i];
            count[bin]++;
        }
        for (int i = 0; i < bins; i++) {
            centre[i] = (count[i] > 0) ? (radiusSum[i] / count[i]) : ((i + 0.5) * width);
            raw[i] = (count[i] > 0) ? (valueSum[i] / peakSum[i]) : 0.0;
        }

        // Only the bins with enough in them are fitted, and the border correction is solved over
        // the same set, so the two cannot disagree about which points the model is answering to.
        int usable = 0;
        for (int i = 0; i < bins; i++) {
            if (count[i] > MIN_BIN_SAMPLES) {
                usable++;
            }
        }
        double[] fitRadius = new double[usable];
        double[] fitRaw = new double[usable];
        int at = 0;
        for (int i = 0; i < bins; i++) {
            if (count[i] > MIN_BIN_SAMPLES) {
                fitRadius[at] = centre[i];
                fitRaw[at] = raw[i];
                at++;
            }
        }

        double level = 0.0;
        boolean settled = false;
        if (usable >= 4) {
            double[] solved = borderLevel(fitRadius, fitRaw, samples.patch);
            level = solved[0];
            settled = solved[1] > 0.5;
        }

        double[] profile = new double[bins];
        for (int i = 0; i < bins; i++) {
            profile[i] = raw[i] * (1.0 - level) + level;
        }
        double[] fitProfile = new double[usable];
        for (int i = 0; i < usable; i++) {
            fitProfile[i] = fitRaw[i] * (1.0 - level) + level;
        }

        Fit psfOnly = (usable >= 3) ? fit(fitRadius, fitProfile, false) : null;
        Fit withPedestal = (usable >= 4) ? fit(fitRadius, fitProfile, true) : null;

        return new Measurement(samples, centre, raw, profile, count, level, settled,
                psfOnly, withPedestal);
    }

    /**
     * How much of the peak the border median took off, and whether it settled.
     *
     * **The patch border is not background.** An aberrated Airy still holds around 1% of its peak
     * out at the edge of a ten pixel patch, so subtracting the border median takes real light off
     * every point of the profile. Against the core that is nothing; against the wings, which are
     * the entire reason for using this model, it is a large fraction - fitting the uncorrected
     * profile against simulated data of known sigma 1.45 and 0.41 waves returns 1.55 and 0.32, a
     * PSF with the right core and wings 30 to 50% too weak.
     *
     * Solved rather than assumed, because how much PSF sits on the border is a property of the
     * PSF being fitted. Each pass lifts the profile by the current estimate, refits, and reads
     * the median of the fitted PSF over the border pixels - a median, matching the median that
     * was taken over the measured ones.
     *
     * Only the PSF term is added back, never the pedestal. A flat term subtracts out of the
     * border median exactly as it subtracts out of the patch, so it never reached the profile and
     * there is nothing to restore; adding it back anyway would make each pass add a fixed amount
     * to the last and never converge.
     */
    private static double[] borderLevel(double[] radii, double[] raw, int patch) {
        double[] borderRadii = borderRadii(patch);
        double level = 0.0;

        for (int pass = 0; pass < BORDER_ITERATIONS; pass++) {
            double[] lifted = new double[raw.length];
            for (int i = 0; i < raw.length; i++) {
                lifted[i] = raw[i] * (1.0 - level) + level;
            }

            Fit fit = fit(radii, lifted, true);
            double[] shape = airyProfile(borderRadii, fit.sigma, fit.waves);
            float[] psf = new float[shape.length];
            for (int i = 0; i < shape.length; i++) {
                psf[i] = (float) (fit.amplitude * shape[i]);
            }

            double updated = smFRETSpotFinder.median(psf, psf.length);
            boolean settled = Math.abs(updated - level) < BORDER_TOLERANCE;
            level = updated;
            if (settled) {
                return new double[] {level, 1.0};
            }
        }
        return new double[] {level, 0.0};
    }

    /** Radii of the patch border pixels, which is what the border median was taken over. */
    private static double[] borderRadii(int patch) {
        int size = 2 * patch + 1;
        List<Double> radii = new ArrayList<>();
        for (int dy = 0; dy < size; dy++) {
            for (int dx = 0; dx < size; dx++) {
                if ((dx == 0) || (dy == 0) || (dx == (size - 1)) || (dy == (size - 1))) {
                    radii.add(Math.hypot(dx - patch, dy - patch));
                }
            }
        }
        double[] out = new double[radii.size()];
        for (int i = 0; i < out.length; i++) {
            out[i] = radii.get(i);
        }
        return out;
    }

    /**
     * Least squares fit of size and aberration, optionally over a flat pedestal.
     *
     * Weighted by 1/sqrt(value) so the wings, which are what a Gaussian gets wrong and therefore
     * the entire point of this model, carry real influence rather than being swamped by the core.
     *
     * The bounds are applied by clamping in a ParameterValidator rather than by a genuinely
     * bounded algorithm, which is what Commons Math offers. That is adequate here because they
     * are guardrails against a runaway rather than active constraints - on both the real and the
     * simulated data the optimum sits well inside them - but a fit that comes back sitting
     * exactly on a bound should be read as "did not converge", not as a measurement.
     */
    static Fit fit(double[] radii, double[] profile, boolean pedestal) {
        int nParameters = pedestal ? 4 : 3;

        double[] weights = new double[profile.length];
        for (int i = 0; i < weights.length; i++) {
            weights[i] = 1.0 / Math.sqrt(Math.max(profile[i], 0.01));
        }

        MultivariateJacobianFunction model = point -> {
            double[] p = point.toArray();
            double[] residual = residual(p, radii, profile, weights, pedestal);

            // Forward differences. The model is a numerical integral with no derivative to hand,
            // and at twenty points by four parameters the extra evaluations are free.
            double[][] jacobian = new double[radii.length][nParameters];
            for (int j = 0; j < nParameters; j++) {
                double step = 1.0e-6 * Math.max(1.0, Math.abs(p[j]));
                double[] shifted = p.clone();
                shifted[j] += step;
                double[] moved = residual(shifted, radii, profile, weights, pedestal);
                for (int i = 0; i < radii.length; i++) {
                    jacobian[i][j] = (moved[i] - residual[i]) / step;
                }
            }
            return new Pair<RealVector, RealMatrix>(new ArrayRealVector(residual, false),
                    new Array2DRowRealMatrix(jacobian, false));
        };

        ParameterValidator clamp = point -> {
            double[] p = point.toArray();
            for (int i = 0; i < nParameters; i++) {
                p[i] = Math.min(FIT_UPPER[i], Math.max(FIT_LOWER[i], p[i]));
            }
            return new ArrayRealVector(p, false);
        };

        LeastSquaresProblem problem = new LeastSquaresBuilder()
                .start(java.util.Arrays.copyOf(FIT_START, nParameters))
                .model(model)
                .target(new double[radii.length])
                .parameterValidator(clamp)
                .lazyEvaluation(false)
                .maxEvaluations(1000)
                .maxIterations(200)
                .build();

        LeastSquaresOptimizer.Optimum optimum = new LevenbergMarquardtOptimizer().optimize(problem);
        double[] p = optimum.getPoint().toArray();

        double sumSquares = 0.0;
        double[] finalResidual = residual(p, radii, profile, weights, pedestal);
        for (double r : finalResidual) {
            sumSquares += r * r;
        }
        return new Fit(p[0], p[1], p[2], pedestal ? p[3] : 0.0,
                Math.sqrt(sumSquares / finalResidual.length));
    }

    /**
     * The radial grid the filter calculation integrates over.
     *
     * Out to fifty pixels because an aberrated Airy still has something there and a wide filter
     * will reach for it - the integrand is the PSF times the filter, and the filter is what makes
     * it converge. A step of 0.05 px is far finer than the pixel grid the measurement came off.
     */
    static final double[] SNR_RADII = snrRadii();

    private static double[] snrRadii() {
        double[] radii = new double[1000];
        for (int i = 0; i < radii.length; i++) {
            radii[i] = i * 0.05;
        }
        return radii;
    }

    // Where the filter search looks, and how finely. The answer is reported to two decimals, so
    // a step of 0.005 is already past what is reportable.
    private static final double FILTER_LOWEST = 0.3;
    private static final double FILTER_HIGHEST = 6.0;
    private static final double FILTER_STEP = 0.005;

    /**
     * The trace SNR a Gaussian filter of this width recovers from this PSF, up to a constant.
     *
     * smFRETAnalyzer measures a trace by convolving with a *normalized* Gaussian of width
     * spotSigma, reading the peak, and multiplying by 4 pi spotSigma^2. For a PSF p(r) carrying N
     * photons on a background of variance sigma_b^2 per pixel, that gives
     *
     *   signal = 4 pi sf^2 . N . integral p(r) G_sf(r) 2 pi r dr
     *   noise  = 4 pi sf^2 . sigma_b . sqrt(integral G_sf^2 dA) = 2 sf sqrt(pi) . sigma_b
     *
     * and the ratio reduces to (1/sf) integral p(r) exp(-r^2 / 2 sf^2) r dr, which is what this
     * returns. **N and sigma_b drop out**, so where it peaks is a property of the PSF's shape
     * alone - not of how bright the molecules are or how high the background is.
     *
     * The profile need not be normalised: a constant factor scales every filter width equally and
     * moves neither the maximum nor the relative width of the region around it.
     *
     * Two assumptions worth knowing. The noise is taken as background dominated and spatially
     * uncorrelated, which is the same assumption spotFilterSNR already makes, and the background
     * estimate is taken as noiseless - it is smoothed at sigma 14, so against a filter of sigma 2
     * its own noise is smaller by more than an order of magnitude. And this is a continuous
     * integral where the real thing runs on pixels; see the sweep in PSFFilterSweepTest for how
     * far apart the two end up.
     */
    static double filterResponse(double[] profile, double filterSigma) {
        double total = 0.0;
        double step = SNR_RADII[1] - SNR_RADII[0];
        double scale = 1.0 / (2.0 * filterSigma * filterSigma);
        for (int i = 0; i < SNR_RADII.length; i++) {
            double r = SNR_RADII[i];
            total += profile[i] * Math.exp(-r * r * scale) * r;
        }
        return (total * step) / filterSigma;
    }

    /** The PSF of a fit, sampled where filterResponse wants it. */
    static double[] snrProfile(Fit fit) {
        return airyProfile(SNR_RADII, fit.sigma, fit.waves);
    }

    /**
     * The best filter width, and how far either side of it barely matters.
     *
     * The band is the point of this as much as the maximum is. On the example data the optimum
     * sits at 1.7 but anything from about 1.3 to 2.2 is within 2% of it, so quoting the maximum
     * on its own would suggest a precision the measurement does not have and send somebody
     * chasing a third decimal that is worth nothing.
     */
    static final class FilterOptimum {
        final double best;
        final double high;
        final double low;
        final double tolerance;

        FilterOptimum(double best, double low, double high, double tolerance) {
            this.best = best;
            this.low = low;
            this.high = high;
            this.tolerance = tolerance;
        }
    }

    /**
     * Where a Gaussian filter recovers the most trace SNR from these PSFs together.
     *
     * With one profile this is that channel's optimum. With two it is the geometric mean of the
     * channels' SNR, each first divided by its own best - so the objective reads as "the typical
     * fraction of what each channel could manage", and a band drawn at 2% of it means 2% per
     * channel rather than 2% of a product that would be 1% each. A FRET trace needs both channels
     * measured well, which is why this is a joint answer and not the donor's.
     */
    static FilterOptimum optimalFilter(double[][] profiles, double tolerance) {
        int steps = (int) Math.round((FILTER_HIGHEST - FILTER_LOWEST) / FILTER_STEP) + 1;

        double[][] response = new double[profiles.length][steps];
        double[] channelBest = new double[profiles.length];
        for (int c = 0; c < profiles.length; c++) {
            for (int i = 0; i < steps; i++) {
                response[c][i] = filterResponse(profiles[c], FILTER_LOWEST + (i * FILTER_STEP));
                channelBest[c] = Math.max(channelBest[c], response[c][i]);
            }
        }

        double[] objective = new double[steps];
        int at = 0;
        for (int i = 0; i < steps; i++) {
            double product = 1.0;
            for (int c = 0; c < profiles.length; c++) {
                product *= (channelBest[c] > 0.0) ? (response[c][i] / channelBest[c]) : 0.0;
            }
            objective[i] = Math.pow(product, 1.0 / profiles.length);
            if (objective[i] > objective[at]) {
                at = i;
            }
        }

        double threshold = objective[at] * (1.0 - tolerance);
        int low = at;
        while ((low > 0) && (objective[low - 1] >= threshold)) {
            low--;
        }
        int high = at;
        while ((high < (steps - 1)) && (objective[high + 1] >= threshold)) {
            high++;
        }

        return new FilterOptimum(FILTER_LOWEST + (at * FILTER_STEP),
                FILTER_LOWEST + (low * FILTER_STEP),
                FILTER_LOWEST + (high * FILTER_STEP), tolerance);
    }

    /** Weighted residual of the model against the measured profile. */
    private static double[] residual(double[] p, double[] radii, double[] profile,
                                     double[] weights, boolean pedestal) {
        double[] shape = airyProfile(radii, p[0], p[1]);
        double offset = pedestal ? p[3] : 0.0;

        double[] out = new double[radii.length];
        for (int i = 0; i < radii.length; i++) {
            out[i] = (p[2] * shape[i] + offset - profile[i]) * weights[i];
        }
        return out;
    }
}
