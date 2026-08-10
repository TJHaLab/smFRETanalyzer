import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.apache.commons.math3.special.BesselJ;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

/**
 * The aberrated Airy pattern itself, checked against what it reduces to.
 *
 * With no aberration the pupil integral has a closed form - the ordinary Airy pattern,
 * (2 J1(v) / v)^2 - so the numerical integral can be scored against an analytic answer rather
 * than against itself. That is the check worth having here, because everything downstream
 * assumes this shape is right and a fit will happily report parameters for a wrong one.
 *
 * The aberrated case has no closed form, so it is pinned by its properties instead: light leaves
 * the core for the wings, monotonically in the aberration, and the total is conserved.
 */
class PSFModelTest {

    /** The analytic unaberrated Airy pattern, peak normalised. */
    private static double analyticAiry(double v) {
        if (Math.abs(v) < 1.0e-9) {
            return 1.0;
        }
        double amplitude = 2.0 * BesselJ.value(1.0, v) / v;
        return amplitude * amplitude;
    }

    /** Radius in pixels corresponding to a dimensionless v, at this sigma. */
    private static double radiusFor(double v, double sigma) {
        return v * smFRETPSF.AIRY_ZERO_SIGMAS * sigma / 3.8317;
    }

    @ParameterizedTest(name = "sigma={0}")
    @ValueSource(doubles = {1.0, 1.36, 2.0, 3.0})
    @DisplayName("with no aberration it is the analytic Airy pattern")
    void unaberratedMatchesTheClosedForm(double sigma) {

        // Out to v = 15, which is past the fourth zero and covers everything a ten pixel patch
        // can see at any sigma the fit allows.
        double[] radii = new double[60];
        double[] expected = new double[radii.length];
        for (int i = 0; i < radii.length; i++) {
            double v = 0.25 * (i + 1);
            radii[i] = radiusFor(v, sigma);
            expected[i] = analyticAiry(v);
        }

        double[] measured = smFRETPSF.airyProfile(radii, sigma, 0.0);
        for (int i = 0; i < radii.length; i++) {
            assertEquals(expected[i], measured[i], 2.0e-6,
                    "at v = " + (0.25 * (i + 1)) + ", sigma " + sigma);
        }
    }

    @Test
    @DisplayName("the profile is normalised to a peak of one")
    void peakIsOne() {
        for (double waves : new double[] {0.0, 0.4, 1.0, 2.0}) {
            double[] value = smFRETPSF.airyProfile(new double[] {0.0}, 1.4, waves);
            assertEquals(1.0, value[0], 1.0e-12, "at " + waves + " waves");
        }
    }

    /**
     * The sizing convention: an unaberrated pattern vanishes at AIRY_ZERO_SIGMAS times the
     * equivalent Gaussian sigma. This is what makes the fitted sigma comparable with spotSigma
     * rather than being a number in its own units.
     */
    @ParameterizedTest(name = "sigma={0}")
    @ValueSource(doubles = {1.0, 1.36, 2.0, 3.0})
    @DisplayName("the first zero lands where the sizing convention says")
    void firstZeroIsWhereItShouldBe(double sigma) {
        double zero = smFRETPSF.AIRY_ZERO_SIGMAS * sigma;

        double[] atZero = smFRETPSF.airyProfile(new double[] {zero}, sigma, 0.0);
        assertTrue(atZero[0] < 1.0e-6, "should vanish at " + zero + " but was " + atZero[0]);

        // And it really is a zero rather than a flat tail: either side of it there is light.
        double[] around = smFRETPSF.airyProfile(
                new double[] {0.8 * zero, 1.2 * zero}, sigma, 0.0);
        assertTrue(around[0] > 1.0e-3, "no light inside the first zero");
        assertTrue(around[1] > 1.0e-3, "no light outside the first zero");
    }

    /**
     * What aberration does to the shape.
     *
     * Note what this can and cannot say. The physical statement is that aberration takes light
     * *out of the core* - the Strehl ratio falls - but this profile is peak normalised, so the
     * peak is 1 by construction and that loss is divided straight back out. What is visible here
     * is the consequence: everything else rises *relative to* the peak, because the peak falls
     * faster than any other part of the pattern. So the peak normalised core broadens rather than
     * narrowing, which is the opposite of what the physical description sounds like it predicts.
     *
     * Monotone only over the range that matters. Past about 0.8 waves the pattern reorganizes and
     * the trend reverses - measured, the wing at six sigma peaks near 0.151 at 0.8 waves and is
     * back to 0.067 by 1.2 - so a fit reporting more than about a wave is not on this branch and
     * should not be read as "even more aberrated".
     */
    @Test
    @DisplayName("aberration broadens the peak normalised pattern, the wings most")
    void aberrationFillsTheWings() {
        double sigma = 1.4;
        double[] core = {1.0 * sigma};
        double[] wing = {6.0 * sigma};

        double firstCore = 0.0;
        double firstWing = 0.0;
        double previousCore = 0.0;
        double previousWing = 0.0;
        for (double waves : new double[] {0.0, 0.1, 0.2, 0.4, 0.8}) {
            double atCore = smFRETPSF.airyProfile(core, sigma, waves)[0];
            double atWing = smFRETPSF.airyProfile(wing, sigma, waves)[0];

            assertTrue(atCore > previousCore,
                    "the core should broaden relative to the peak, at " + waves + " waves");
            assertTrue(atWing > previousWing,
                    "the wings should fill, at " + waves + " waves");
            if (waves == 0.0) {
                firstCore = atCore;
                firstWing = atWing;
            }
            previousCore = atCore;
            previousWing = atWing;
        }

        // And the wings gain far more than the core does - which is the whole reason a fit can
        // tell aberration apart from simply a wider spot.
        double coreGain = previousCore / firstCore;
        double wingGain = previousWing / firstWing;
        assertTrue(wingGain > 10.0 * coreGain,
                "the wings gained " + wingGain + "x against the core's " + coreGain
                        + "x, which is not enough separation to fit the two apart");
    }

    /**
     * The number that makes this model worth having. A Gaussian of the same core width has
     * essentially nothing left at five to ten pixels, where the real spots carry percents of
     * their peak - four orders of magnitude, which is what a background estimator has to cope
     * with and what makes neighbours contaminate each other.
     */
    @Test
    @DisplayName("the wings are orders of magnitude above a Gaussian's")
    void wingsAreFarAboveAGaussian() {
        double sigma = 1.36;

        // The gap is steeply radius dependent, because a Gaussian falls as exp(-r^2) while this
        // falls roughly as a power law - so a single threshold across the range would be either
        // vacuous at ten pixels or false at five. Measured ratios at 0.41 waves are 48, 6.3e5 and
        // 2.6e9; these are those with a wide margin.
        double[] radii = {5.0, 8.0, 10.0};
        double[] leastRatio = {20.0, 1.0e5, 1.0e8};

        double[] airy = smFRETPSF.airyProfile(radii, sigma, 0.41);
        for (int i = 0; i < radii.length; i++) {
            double gaussian = Math.exp(-(radii[i] * radii[i]) / (2.0 * sigma * sigma));
            assertTrue(airy[i] > (leastRatio[i] * gaussian),
                    "at " + radii[i] + " px the Airy holds " + airy[i]
                            + " against a Gaussian's " + gaussian
                            + ", a ratio of " + (airy[i] / gaussian));

            // The absolute claim, and the one that matters for contamination: percents of peak
            // out where a Gaussian has nothing at all.
            assertTrue(airy[i] > 0.002,
                    "at " + radii[i] + " px the model holds only " + airy[i] + " of its peak");
        }
    }

    /**
     * The fit has to invert the model: given a profile the model itself produced, it must return
     * the parameters that produced it. Anything less and a measurement is uninterpretable.
     */
    @ParameterizedTest(name = "sigma={0}")
    @ValueSource(doubles = {1.0, 1.36, 2.0, 2.5})
    @DisplayName("the fit recovers parameters it was given")
    void fitRoundTrips(double sigma) {
        double waves = 0.41;
        double amplitude = 1.0;
        double pedestal = 0.03;

        double[] radii = new double[20];
        for (int i = 0; i < radii.length; i++) {
            radii[i] = 0.5 * i;
        }
        double[] shape = smFRETPSF.airyProfile(radii, sigma, waves);
        double[] profile = new double[radii.length];
        for (int i = 0; i < radii.length; i++) {
            profile[i] = amplitude * shape[i] + pedestal;
        }

        smFRETPSF.Fit fit = smFRETPSF.fit(radii, profile, true);

        assertEquals(sigma, fit.sigma, 0.02, "sigma");
        assertEquals(waves, fit.waves, 0.02, "waves");
        assertEquals(amplitude, fit.amplitude, 0.02, "amplitude");
        assertEquals(pedestal, fit.pedestal, 0.005, "pedestal");
        assertTrue(fit.rms < 1.0e-3, "rms was " + fit.rms + " on a profile the model generated");
    }

    @Test
    @DisplayName("a fit without a pedestal reports none")
    void fitWithoutPedestal() {
        double[] radii = new double[20];
        for (int i = 0; i < radii.length; i++) {
            radii[i] = 0.5 * i;
        }
        double[] profile = smFRETPSF.airyProfile(radii, 1.5, 0.3);

        smFRETPSF.Fit fit = smFRETPSF.fit(radii, profile, false);
        assertEquals(0.0, fit.pedestal, 0.0);
        assertEquals(1.5, fit.sigma, 0.02);
        assertEquals(0.3, fit.waves, 0.02);
    }

    /**
     * The J0 lookup table the aperture integral runs on. It exists because a fit asks for
     * millions of Bessel evaluations, and it is only safe if it agrees with the real routine.
     */
    @Test
    @DisplayName("the tabulated J0 agrees with the library")
    void tabulatedBesselIsAccurate() {
        for (double x = 0.0; x < 60.0; x += 0.017) {

            // Reached through airyProfile, which is the only way the table is used.
            double[] one = smFRETPSF.airyProfile(new double[] {x}, 1.0, 0.0);
            double v = x * 3.8317 / smFRETPSF.AIRY_ZERO_SIGMAS;
            assertEquals(analyticAiry(v), one[0], 2.0e-6, "at radius " + x);
        }
    }
}
