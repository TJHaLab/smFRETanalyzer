import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

/**
 * The correction for what one-sided clipping does to the background level.
 *
 * Keeping only the pixels below kappa spread leaves a truncated sample, so the level reads low
 * even where there is nothing but noise to clip - about 1 ADU, which is 50 units of trace. Both
 * halves of the fix are checked: the mean offset phi(k)/Phi(k), and the factor that unbiases the
 * MAD, which robustSpread reads off the already-truncated residuals and would otherwise report
 * too small, making the correction scaled by it too small as well.
 *
 * The reference values are scipy's, to nine decimals. The production normal CDF is Abramowitz and
 * Stegun 7.1.26, good to about 1e-7, so the tolerances here are that accuracy rather than the
 * reference's.
 */
class TruncationConstantsTest {

    // A&S 7.1.26 is quoted at |error| < 1.5e-7 on the erf, so half that on the CDF. Measured
    // worst case over the points below is 6.9e-8, at the 97.5% quantile.
    private static final double CDF_ACCURACY = 2.0e-7;

    // The constants are solved from that CDF by bisection, evaluating it several times per
    // answer, so their error is a small multiple of the CDF's rather than equal to it -
    // measured at 3.9e-7 on the MAD factor at kappa 3.
    private static final double CONSTANT_ACCURACY = 1.0e-6;

    @ParameterizedTest(name = "kappa={0}")
    @CsvSource({
        "0.5, 0.509160434, 0.718043485",
        "1.0, 0.287599971, 0.838502109",
        "1.3, 0.189735035, 0.896395001",
        "2.0, 0.055247863, 0.974016520",
        "3.0, 0.004437839, 0.998425994",
    })
    @DisplayName("mean offset and MAD factor match scipy")
    void matchesScipy(double kappa, double expectedOffset, double expectedMadFactor) {
        smFRETSpotFinder.truncationConstants(kappa);
        assertEquals(expectedOffset, smFRETSpotFinder.cachedMeanOffset, CONSTANT_ACCURACY,
                "mean offset at kappa " + kappa);
        assertEquals(expectedMadFactor, smFRETSpotFinder.cachedMadFactor, CONSTANT_ACCURACY,
                "MAD factor at kappa " + kappa);
    }

    /**
     * Both limits, which is what says the correction is a correction rather than a fudge: a clip
     * far out in the tail removes nothing, so it must ask for nothing back.
     */
    @Test
    @DisplayName("a distant clip corrects by nothing")
    void bothConstantsReachTheirUntruncatedLimits() {
        smFRETSpotFinder.truncationConstants(8.0);
        assertEquals(0.0, smFRETSpotFinder.cachedMeanOffset, 1.0e-9);
        assertEquals(1.0, smFRETSpotFinder.cachedMadFactor, 1.0e-5);
    }

    /**
     * A tighter clip cuts more of the sample away, so it has to put more back and the MAD it
     * reads is further below the truth. Monotone in kappa, both of them.
     */
    @Test
    @DisplayName("a tighter clip needs a larger correction")
    void correctionsAreMonotoneInKappa() {
        double previousOffset = Double.MAX_VALUE;
        double previousMadFactor = -1.0;
        for (double kappa = 0.4; kappa <= 6.0; kappa += 0.2) {
            smFRETSpotFinder.truncationConstants(kappa);
            double offset = smFRETSpotFinder.cachedMeanOffset;
            double madFactor = smFRETSpotFinder.cachedMadFactor;

            assertTrue(offset < previousOffset,
                    "mean offset should fall with kappa, at " + kappa);
            assertTrue(madFactor > previousMadFactor,
                    "MAD factor should rise with kappa, at " + kappa);
            assertTrue((madFactor > 0.0) && (madFactor <= 1.0 + 1.0e-9),
                    "MAD factor out of range at " + kappa + ": " + madFactor);

            previousOffset = offset;
            previousMadFactor = madFactor;
        }
    }

    /**
     * The cache is keyed on kappa alone and is static, so it outlives any one estimate. A stale
     * hit would apply one kappa's correction to another's clip.
     */
    @Test
    @DisplayName("the cache re-solves when kappa changes")
    void cacheTracksKappa() {
        smFRETSpotFinder.truncationConstants(1.3);
        double atThirteen = smFRETSpotFinder.cachedMeanOffset;

        smFRETSpotFinder.truncationConstants(3.0);
        smFRETSpotFinder.truncationConstants(1.3);

        assertEquals(atThirteen, smFRETSpotFinder.cachedMeanOffset, 0.0);
    }

    @Test
    @DisplayName("the normal CDF is accurate over the range the clip uses")
    void normalCdfIsAccurate() {

        // Reference values from scipy.stats.norm.cdf.
        assertEquals(0.5, smFRETSpotFinder.normalCdf(0.0), CDF_ACCURACY);
        assertEquals(0.841344746, smFRETSpotFinder.normalCdf(1.0), CDF_ACCURACY);
        assertEquals(0.158655254, smFRETSpotFinder.normalCdf(-1.0), CDF_ACCURACY);
        assertEquals(0.975000000, smFRETSpotFinder.normalCdf(1.959963985), CDF_ACCURACY);
        assertEquals(0.903199515, smFRETSpotFinder.normalCdf(1.3), CDF_ACCURACY);
        assertEquals(0.999968329, smFRETSpotFinder.normalCdf(4.0), CDF_ACCURACY);

        // Symmetry is exact away from the origin, and it is exact rather than approximate because
        // both sides evaluate the same polynomial in |x| and differ only in the sign applied to
        // the result - so whatever the approximation's error is, it cancels in the sum.
        for (double x = 0.25; x <= 5.0; x += 0.25) {
            assertEquals(1.0, smFRETSpotFinder.normalCdf(x) + smFRETSpotFinder.normalCdf(-x),
                    1.0e-15, "symmetry at " + x);
        }

        // Zero is the one place it is not exact, and only because negative zero does not test as
        // less than zero, so both calls take the positive branch and the approximation's residual
        // at the origin - 1e-9 rather than 0 - is added twice instead of cancelling. Far below
        // anything the clip correction can notice, but it is why the loop above starts past it.
        assertEquals(0.5, smFRETSpotFinder.normalCdf(0.0), CDF_ACCURACY);
        assertEquals(1.0, smFRETSpotFinder.normalCdf(0.0) + smFRETSpotFinder.normalCdf(-0.0),
                1.0e-8);
    }

    /**
     * The 90th percentile of a standard normal is the constant the prominence ring's noise bias
     * is taken off with, so this quantile and promNoiseBias have to agree.
     */
    @Test
    @DisplayName("the quantile inverts the CDF")
    void normalQuantileInvertsTheCdf() {
        for (double p = 0.02; p < 0.99; p += 0.02) {
            double x = smFRETSpotFinder.normalQuantile(p);
            assertEquals(p, smFRETSpotFinder.normalCdf(x), CDF_ACCURACY, "round trip at p " + p);
        }
        assertEquals(1.2815515655, smFRETSpotFinder.normalQuantile(0.9), 1.0e-5);
    }
}
