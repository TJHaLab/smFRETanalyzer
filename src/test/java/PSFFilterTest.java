import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

/**
 * Which spotSigma extracts the most trace SNR, given the PSF that was measured.
 *
 * smFRETAnalyzer measures a trace with a Gaussian matched filter, but the PSF is not a Gaussian -
 * an aberrated Airy carries a real fraction of its light out where a Gaussian of the same core
 * has nothing, and a filter matched to the core alone throws that away. So the best spotSigma is
 * not the fitted core width, and how far off it is depends on the aberration.
 *
 * The case with a known answer anchors all of it: for a Gaussian PSF the best Gaussian filter is
 * the PSF's own width, exactly. That is the classic matched filter result and it is what says the
 * derivation behind filterResponse is right rather than merely plausible.
 */
class PSFFilterTest {

    private static final double TOLERANCE = 0.02;

    private static double[] gaussianProfile(double sigma) {
        double[] profile = new double[smFRETPSF.SNR_RADII.length];
        for (int i = 0; i < profile.length; i++) {
            double r = smFRETPSF.SNR_RADII[i];
            profile[i] = Math.exp(-(r * r) / (2.0 * sigma * sigma));
        }
        return profile;
    }

    private static double optimumFor(double[] profile) {
        return smFRETPSF.optimalFilter(new double[][] {profile}, TOLERANCE).best;
    }

    /**
     * The anchor. A matched filter is optimal when it matches, so a Gaussian PSF must want a
     * Gaussian filter of exactly its own width - and it does, to the resolution of the search.
     */
    @ParameterizedTest(name = "sigma={0}")
    @ValueSource(doubles = {1.0, 1.5, 2.0, 2.5, 3.0})
    @DisplayName("a Gaussian PSF wants a filter of its own width")
    void gaussianWantsItsOwnWidth(double sigma) {
        assertEquals(sigma, optimumFor(gaussianProfile(sigma)), 0.01,
                "the matched filter result should hold exactly");
    }

    /**
     * The problem has no length scale of its own, so doubling the PSF doubles the answer. That is
     * what collapses a table over sigma and aberration into one function of aberration alone.
     */
    @ParameterizedTest(name = "waves={0}")
    @ValueSource(doubles = {0.0, 0.2, 0.4, 0.5})
    @DisplayName("the optimum scales with the PSF, so only the aberration matters")
    void theOptimumIsScaleFree(double waves) {
        double reference = optimumFor(smFRETPSF.airyProfile(smFRETPSF.SNR_RADII, 1.0, waves))
                / 1.0;

        for (double sigma : new double[] {1.0, 1.5, 2.0, 2.5}) {
            double ratio = optimumFor(smFRETPSF.airyProfile(smFRETPSF.SNR_RADII, sigma, waves))
                    / sigma;
            assertEquals(reference, ratio, 0.01,
                    "optimum/sigma at sigma " + sigma + ", " + waves + " waves");
        }
    }

    /**
     * What aberration does to the answer, and the reason this feature exists: light moves out of
     * the core into the wings, so the filter that collects the most of it widens. Below about 0.3
     * waves it barely moves, which is worth knowing before anyone retunes on a clean system.
     */
    @Test
    @DisplayName("aberration widens the best filter, and only past about 0.3 waves")
    void aberrationWidensTheOptimum() {
        double sigma = 1.5;
        double previous = 0.0;
        double atZero = 0.0;
        double atThree = 0.0;

        for (double waves : new double[] {0.0, 0.1, 0.2, 0.3, 0.4, 0.5}) {
            double best = optimumFor(smFRETPSF.airyProfile(smFRETPSF.SNR_RADII, sigma, waves));
            assertTrue(best > previous,
                    "the optimum should widen with aberration, at " + waves + " waves");
            if (waves == 0.0) {
                atZero = best;
            }
            if (waves == 0.3) {
                atThree = best;
            }
            previous = best;
        }

        assertTrue((atThree / atZero) < 1.15,
                "up to 0.3 waves the optimum moved by " + (100.0 * ((atThree / atZero) - 1.0))
                        + "%, which is more than 'barely'");
        assertTrue((previous / atZero) > 1.35,
                "by 0.5 waves the optimum should have moved substantially, moved to "
                        + (previous / atZero));
    }

    /**
     * An unaberrated Airy is not quite the Gaussian it is sized against - the equivalence between
     * the two is an approximation - so its best filter sits a little under its nominal sigma.
     * Recorded because a reader expecting exactly 1.00 would otherwise read 0.98 as an error.
     */
    @Test
    @DisplayName("an unaberrated Airy wants slightly less than its nominal sigma")
    void unaberratedAiryIsNotQuiteAGaussian() {
        double ratio = optimumFor(smFRETPSF.airyProfile(smFRETPSF.SNR_RADII, 2.0, 0.0)) / 2.0;
        assertTrue((ratio > 0.95) && (ratio < 1.0),
                "optimum/sigma for an unaberrated Airy was " + ratio);
    }

    /**
     * The band, which is the more useful half of the answer. The maximum of a matched filter
     * response is quadratically flat around its peak, so a wide range of settings is nearly as
     * good and quoting the maximum alone overstates what is known.
     */
    @Test
    @DisplayName("the band brackets the optimum and is genuinely wide")
    void theBandBracketsTheOptimum() {
        double[] profile = smFRETPSF.airyProfile(smFRETPSF.SNR_RADII, 1.35, 0.42);
        smFRETPSF.FilterOptimum optimum =
                smFRETPSF.optimalFilter(new double[][] {profile}, TOLERANCE);

        assertTrue(optimum.low <= optimum.best, "band should contain the optimum");
        assertTrue(optimum.high >= optimum.best, "band should contain the optimum");
        assertTrue((optimum.high - optimum.low) > 0.5,
                "the 2% band was only " + optimum.low + " to " + optimum.high
                        + ", which is tighter than a matched filter's peak can be");

        // Everything inside the band really is within the tolerance, and just outside it is not.
        double best = smFRETPSF.filterResponse(profile, optimum.best);
        assertTrue(smFRETPSF.filterResponse(profile, optimum.low) >= (best * (1.0 - TOLERANCE)),
                "the low end of the band is outside its own tolerance");
        assertTrue(smFRETPSF.filterResponse(profile, optimum.high) >= (best * (1.0 - TOLERANCE)),
                "the high end of the band is outside its own tolerance");
        assertTrue(smFRETPSF.filterResponse(profile, optimum.low - 0.05) < best,
                "outside the band should be worse");
    }

    /**
     * A tighter tolerance can only narrow the band, never widen it.
     */
    @Test
    @DisplayName("the band narrows as the tolerance tightens")
    void tighterToleranceNarrowsTheBand() {
        double[] profile = smFRETPSF.airyProfile(smFRETPSF.SNR_RADII, 1.35, 0.42);

        double previous = Double.MAX_VALUE;
        for (double tolerance : new double[] {0.05, 0.02, 0.01, 0.002}) {
            smFRETPSF.FilterOptimum optimum =
                    smFRETPSF.optimalFilter(new double[][] {profile}, tolerance);
            double width = optimum.high - optimum.low;
            assertTrue(width <= previous,
                    "the band widened as the tolerance tightened, at " + tolerance);
            previous = width;
        }
    }

    /**
     * spotSigma is one parameter for both channels, so the recommendation is a joint one - and it
     * has to land between the two channels' own answers rather than following either.
     */
    @Test
    @DisplayName("the joint optimum sits between the two channels'")
    void theJointOptimumIsBetweenTheChannels() {
        double[] donor = smFRETPSF.airyProfile(smFRETPSF.SNR_RADII, 1.343, 0.418);
        double[] acceptor = smFRETPSF.airyProfile(smFRETPSF.SNR_RADII, 1.552, 0.362);

        double donorBest = optimumFor(donor);
        double acceptorBest = optimumFor(acceptor);
        double joint = smFRETPSF.optimalFilter(new double[][] {donor, acceptor}, TOLERANCE).best;

        assertTrue(donorBest < acceptorBest, "the fixture should have the two disagreeing");
        assertTrue((joint > donorBest) && (joint < acceptorBest),
                "joint optimum " + joint + " is not between " + donorBest + " and "
                        + acceptorBest);
    }

    /**
     * Two identical channels must give the same answer as one, or the geometric mean is doing
     * something other than averaging.
     */
    @Test
    @DisplayName("duplicating a channel changes nothing")
    void duplicateChannelsAgreeWithOne() {
        double[] profile = smFRETPSF.airyProfile(smFRETPSF.SNR_RADII, 1.5, 0.4);

        smFRETPSF.FilterOptimum single =
                smFRETPSF.optimalFilter(new double[][] {profile}, TOLERANCE);
        smFRETPSF.FilterOptimum doubled =
                smFRETPSF.optimalFilter(new double[][] {profile, profile}, TOLERANCE);

        assertEquals(single.best, doubled.best, 1.0e-9, "optimum");
        assertEquals(single.low, doubled.low, 1.0e-9, "band low");
        assertEquals(single.high, doubled.high, 1.0e-9, "band high");
    }

    /**
     * The measured answer on the example data, so a change in the model or the integration shows
     * up as a number rather than as a shrug. The fits are the ones hel1 gives.
     */
    @Test
    @DisplayName("hel1's measured PSF wants a filter near 1.7")
    void hel1LandsWhereItDid() {
        double[] donor = smFRETPSF.airyProfile(smFRETPSF.SNR_RADII, 1.343, 0.418);
        double[] acceptor = smFRETPSF.airyProfile(smFRETPSF.SNR_RADII, 1.552, 0.362);
        smFRETPSF.FilterOptimum joint =
                smFRETPSF.optimalFilter(new double[][] {donor, acceptor}, TOLERANCE);

        assertEquals(1.72, joint.best, 0.06, "joint optimum on hel1");

        // And the reason the macros' spotsigma of 2 was never costing anything measurable.
        double[][] both = {donor, acceptor};
        smFRETPSF.FilterOptimum band = smFRETPSF.optimalFilter(both, TOLERANCE);
        assertTrue((band.low < 1.5) && (band.high > 2.0),
                "the 2% band was " + band.low + " to " + band.high
                        + ", which should comfortably contain both the fitted core and 2.0");
    }
}
