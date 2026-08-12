import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

/**
 * The binning behind the trace histograms, and the corrections applied before it.
 *
 * Three properties here are exact identities rather than tolerances, which is what makes them
 * worth having: a leakage of 1 collapses the Total histogram onto the raw acceptor because the
 * donor cancels out of the sum; a leakage of f moves the acceptor by f times the donor; and
 * setting every correction back to zero restores the filter bounds exactly. Any change to the
 * correction arithmetic that is wrong will break at least one of them.
 *
 * The traces are injected rather than loaded from an H5, so these run without touching the
 * HDF5 library or the disk. PipelineOutputTest covers the reading path.
 */
class HistogramTest {

    private static smFRETTraceHistogram withTraces(float[][] donor, float[][] acceptor) {
        smFRETTraceHistogram histogram = new smFRETTraceHistogram();
        histogram.targetTraces = donor;
        histogram.sourceTraces = acceptor;
        histogram.nSpots = donor.length;
        histogram.nFrames = donor[0].length;
        return histogram;
    }

    /**
     * One trace, so the auto-ranged histogram's lower edge *is* that trace's value - which turns
     * the binning into an exact readout of the arithmetic under test.
     */
    private static double singleValue(smFRETTraceHistogram histogram, int type,
                                      smFRETTraceHistogram.Corrections corrections) {
        smFRETTraceHistogram.Histogram result = histogram.computeHistogram(type, 1,
                histogram.nFrames, smFRETTraceHistogram.Filters.none(), 10, corrections);
        assertEquals(1, result.nSpotsUsed, "the fixture should contribute exactly one point");
        return result.lo;
    }

    /**
     * The same readout for FRET, which cannot use the trick above: its range is *fixed* at -0.2
     * to 1.2 rather than auto-ranged, so the lower edge is that constant and says nothing about
     * the data. The bin the single trace landed in does, to within a bin width - 0.001 here.
     */
    private static double singleFretValue(smFRETTraceHistogram histogram,
                                          smFRETTraceHistogram.Corrections corrections) {
        int bins = 1400;
        smFRETTraceHistogram.Histogram result = histogram.computeHistogram(
                smFRETTraceHistogram.TYPE_FRET, 1, histogram.nFrames,
                smFRETTraceHistogram.Filters.none(), bins, corrections);

        assertEquals(1, result.nPoints, "the fixture should bin exactly one point");
        for (int i = 0; i < bins; i++) {
            if (result.counts[i] > 0) {
                return result.lo + (i + 0.5) * result.binWidth;
            }
        }
        throw new AssertionError("no populated bin");
    }

    private static smFRETTraceHistogram.Corrections none() {
        return new smFRETTraceHistogram.Corrections(0.0, 0.0, 0.0);
    }

    /**
     * D + (A - D) = A, exactly. The Total histogram under a leakage of 1 has to be the
     * uncorrected acceptor histogram, bin for bin - which also says the correction is applied
     * before the averaging rather than after.
     */
    @Test
    @DisplayName("a leakage of one collapses Total onto the raw acceptor")
    void leakageOfOneCancelsTheDonor() {
        float[][] donor = {{100.0f, 120.0f, 90.0f}, {40.0f, 45.0f, 50.0f}, {300.0f, 280.0f, 310.0f}};
        float[][] acceptor = {{200.0f, 180.0f, 220.0f}, {10.0f, 12.0f, 8.0f}, {60.0f, 70.0f, 65.0f}};
        smFRETTraceHistogram histogram = withTraces(donor, acceptor);

        smFRETTraceHistogram.Histogram total = histogram.computeHistogram(
                smFRETTraceHistogram.TYPE_TOTAL, 1, 3, smFRETTraceHistogram.Filters.none(), 16,
                new smFRETTraceHistogram.Corrections(0.0, 0.0, 1.0));

        smFRETTraceHistogram.Histogram rawAcceptor = histogram.computeHistogram(
                smFRETTraceHistogram.TYPE_ACCEPTOR, 1, 3, smFRETTraceHistogram.Filters.none(), 16, none());

        assertArrayEquals(rawAcceptor.counts, total.counts, "bin counts");
        assertEquals(rawAcceptor.lo, total.lo, 1.0e-12, "lower edge");
        assertEquals(rawAcceptor.binWidth, total.binWidth, 1.0e-12, "bin width");
    }

    /**
     * Leakage comes off the acceptor as a fraction of the *baseline corrected* donor, since
     * leakage is a fraction of real donor emission and a residual offset is not. Both halves of
     * that are checked: the size of the shift, and that a donor baseline changes it.
     */
    @Test
    @DisplayName("leakage moves the acceptor by that fraction of the donor")
    void leakageMovesTheAcceptorByAFractionOfTheDonor() {
        float[][] donor = {{100.0f, 200.0f, 300.0f}};       // mean 200
        float[][] acceptor = {{50.0f, 60.0f, 70.0f}};       // mean 60
        smFRETTraceHistogram histogram = withTraces(donor, acceptor);

        double uncorrected = singleValue(histogram, smFRETTraceHistogram.TYPE_ACCEPTOR, none());
        assertEquals(60.0, uncorrected, 1.0e-6);

        double leaked = singleValue(histogram, smFRETTraceHistogram.TYPE_ACCEPTOR,
                new smFRETTraceHistogram.Corrections(0.0, 0.0, 0.08));
        assertEquals(60.0 - 0.08 * 200.0, leaked, 1.0e-6, "0.08 of the donor mean");

        // With a donor baseline of 50 the donor mean is 150, so the leakage subtracted is
        // 0.08 * 150 rather than 0.08 * 200.
        double withBaseline = singleValue(histogram, smFRETTraceHistogram.TYPE_ACCEPTOR,
                new smFRETTraceHistogram.Corrections(50.0, 0.0, 0.08));
        assertEquals(60.0 - 0.08 * 150.0, withBaseline, 1.0e-6,
                "leakage should use the baseline corrected donor");
    }

    @Test
    @DisplayName("baselines are subtracted, so a negative one adds")
    void baselinesAreSubtracted() {
        float[][] donor = {{100.0f, 100.0f}};
        float[][] acceptor = {{40.0f, 40.0f}};
        smFRETTraceHistogram histogram = withTraces(donor, acceptor);

        assertEquals(70.0, singleValue(histogram, smFRETTraceHistogram.TYPE_DONOR,
                new smFRETTraceHistogram.Corrections(30.0, 0.0, 0.0)), 1.0e-6);
        assertEquals(130.0, singleValue(histogram, smFRETTraceHistogram.TYPE_DONOR,
                new smFRETTraceHistogram.Corrections(-30.0, 0.0, 0.0)), 1.0e-6);
        assertEquals(15.0, singleValue(histogram, smFRETTraceHistogram.TYPE_ACCEPTOR,
                new smFRETTraceHistogram.Corrections(0.0, 25.0, 0.0)), 1.0e-6);
    }

    /**
     * FRET is the ratio of the interval averaged intensities, not the average of the per frame
     * ratios. The two differ by about 0.2 in the mean on the example data, so a fixture where
     * they differ obviously is the right one to pin it with.
     */
    @Test
    @DisplayName("FRET is the ratio of the averages, not the average of the ratios")
    void fretIsARatioOfAverages() {
        float[][] donor = {{10.0f, 1000.0f}};
        float[][] acceptor = {{1000.0f, 10.0f}};
        smFRETTraceHistogram histogram = withTraces(donor, acceptor);

        // mean(A) / (mean(D) + mean(A)) = 505 / 1010 = 0.5. The mean of the per frame ratios is
        // also 0.5 here, since the fixture is symmetric - so this one only fixes the value.
        assertEquals(0.5, singleFretValue(histogram, none()), 1.0e-3);

        // An asymmetric fixture, where the two definitions genuinely disagree: the ratio of the
        // averages is 100 / 150.5 = 0.664, where the average of the per frame ratios 100/101 and
        // 100/200 would be 0.745.
        smFRETTraceHistogram other = withTraces(new float[][] {{1.0f, 100.0f}},
                new float[][] {{100.0f, 100.0f}});
        assertEquals(100.0 / 150.5, singleFretValue(other, none()), 1.0e-3);
        assertTrue(Math.abs(singleFretValue(other, none()) - 0.745) > 0.05,
                "this is the value averaging the per frame ratios would give");
    }

    /**
     * The intensity range is an all-frames test at *both* ends, which makes it far stricter than
     * it looks: one frame outside drops the whole trace rather than diluting its average. That is
     * what catches a molecule which bleaches partway through, and equally what makes a single
     * bright spike enough to lose a trace.
     */
    @Test
    @DisplayName("one frame outside the range drops the whole trace")
    void theRangeFilterIsAnAllFramesTest() {
        float[][] donor = {{100.0f, 100.0f, 100.0f}, {100.0f, 100.0f, 100.0f}};
        float[][] acceptor = {{100.0f, 100.0f, 100.0f}, {100.0f, 5000.0f, 100.0f}};
        smFRETTraceHistogram histogram = withTraces(donor, acceptor);

        smFRETTraceHistogram.Histogram unfiltered = histogram.computeHistogram(
                smFRETTraceHistogram.TYPE_TOTAL, 1, 3, smFRETTraceHistogram.Filters.none(), 8, none());
        assertEquals(2, unfiltered.nSpotsUsed);

        // The second trace averages 1900, well inside a maximum of 3000 - but one of its frames
        // is 5100, so it goes entirely.
        smFRETTraceHistogram.Histogram filtered = histogram.computeHistogram(
                smFRETTraceHistogram.TYPE_TOTAL, 1, 3, smFRETTraceHistogram.Filters.only(smFRETTraceHistogram.FILTER_TOTAL, -1.0e9, 3000.0), 8, none());
        assertEquals(1, filtered.nSpotsUsed, "the spiking trace should be dropped whole");
    }

    @Test
    @DisplayName("the range applies to whichever quantity is selected")
    void theRangeAppliesToTheSelectedQuantity() {
        float[][] donor = {{500.0f, 500.0f}};
        float[][] acceptor = {{10.0f, 10.0f}};
        smFRETTraceHistogram histogram = withTraces(donor, acceptor);

        // A maximum of 100 on the acceptor keeps it; the same maximum on the total does not.
        assertEquals(1, histogram.computeHistogram(smFRETTraceHistogram.TYPE_TOTAL, 1, 2,
                smFRETTraceHistogram.Filters.only(smFRETTraceHistogram.FILTER_ACCEPTOR, -1.0e9, 100.0), 8, none()).nSpotsUsed);
        assertEquals(0, histogram.computeHistogram(smFRETTraceHistogram.TYPE_TOTAL, 1, 2,
                smFRETTraceHistogram.Filters.only(smFRETTraceHistogram.FILTER_TOTAL, -1.0e9, 100.0), 8, none()).nSpotsUsed);
    }

    /**
     * The three ranges are ANDed, which is the whole of issue #8.
     *
     * They used to be a combo box beside one slider, so only one could be live and filtering on
     * donor brightness meant giving up filtering on total brightness. Nothing about the two made
     * them alternatives except the widget.
     */
    @Test
    @DisplayName("a trace has to pass all three ranges, not the selected one")
    void allThreeRangesApply() {
        // Donor 500, acceptor 10, total 510.
        float[][] donor = {{500.0f, 500.0f}};
        float[][] acceptor = {{10.0f, 10.0f}};
        smFRETTraceHistogram histogram = withTraces(donor, acceptor);

        double[] min = {-1.0e9, -1.0e9, -1.0e9};
        double[] max = {1.0e9, 1.0e9, 1.0e9};

        // Wide open on all three keeps it, which is the state the sliders start in.
        assertEquals(1, kept(histogram, new smFRETTraceHistogram.Filters(min, max)));

        // A donor range that keeps it, and an acceptor range that does not. Under the old
        // single filter this combination could not be expressed at all.
        min[smFRETTraceHistogram.FILTER_DONOR] = 100.0;
        max[smFRETTraceHistogram.FILTER_DONOR] = 900.0;
        assertEquals(1, kept(histogram, new smFRETTraceHistogram.Filters(min, max)),
                "the donor range alone should keep it");

        max[smFRETTraceHistogram.FILTER_ACCEPTOR] = 5.0;
        assertEquals(0, kept(histogram, new smFRETTraceHistogram.Filters(min, max)),
                "the acceptor range should drop it even though the donor range keeps it");
    }

    /**
     * Which filter did the rejecting, counted per filter rather than once per trace.
     *
     * A trace failing two of them counts in both, so the counts do not sum to the number
     * rejected. That is deliberate: it is what says whether widening one slider on its own
     * would bring anything back, which is the question three live sliders create.
     */
    @Test
    @DisplayName("rejections are attributed to every filter that would have caused them")
    void rejectionsAreAttributed() {
        float[][] donor = {{500.0f, 500.0f}};
        float[][] acceptor = {{10.0f, 10.0f}};
        smFRETTraceHistogram histogram = withTraces(donor, acceptor);

        double[] min = {-1.0e9, -1.0e9, -1.0e9};
        double[] max = {1.0e9, 1.0e9, 1.0e9};
        max[smFRETTraceHistogram.FILTER_DONOR] = 1.0;      // 500 fails this
        max[smFRETTraceHistogram.FILTER_TOTAL] = 1.0;      // 510 fails this too

        smFRETTraceHistogram.Histogram result = histogram.computeHistogram(
                smFRETTraceHistogram.TYPE_TOTAL, 1, 2,
                new smFRETTraceHistogram.Filters(min, max), 8, none());

        assertEquals(0, result.nSpotsUsed, "nothing should survive");
        assertEquals(1, result.nRejectedBy[smFRETTraceHistogram.FILTER_DONOR], "donor");
        assertEquals(1, result.nRejectedBy[smFRETTraceHistogram.FILTER_TOTAL], "total");
        assertEquals(0, result.nRejectedBy[smFRETTraceHistogram.FILTER_ACCEPTOR],
                "the acceptor range was wide open and must not be blamed");
    }

    /** Traces surviving the three ranges. */
    private static int kept(smFRETTraceHistogram histogram, smFRETTraceHistogram.Filters filters) {
        return histogram.computeHistogram(smFRETTraceHistogram.TYPE_TOTAL, 1, 2,
                filters, 8, none()).nSpotsUsed;
    }

    /**
     * Naming the range that emptied the histogram, which is what issue #8 left open.
     *
     * Three live ranges can be set to an empty intersection. Preventing that would mean
     * recomputing every slider's limits from the traces on every drag; saying which one did it
     * costs nothing and answers the question the user actually has.
     */
    @Test
    @DisplayName("an empty result names the range responsible")
    void anEmptyResultIsExplained() {
        // Two traces, donor 500 / acceptor 10 and donor 20 / acceptor 400.
        float[][] donor = {{500.0f, 500.0f}, {20.0f, 20.0f}};
        float[][] acceptor = {{10.0f, 10.0f}, {400.0f, 400.0f}};
        smFRETTraceHistogram histogram = withTraces(donor, acceptor);

        double[] min = {-1.0e9, -1.0e9, -1.0e9};
        double[] max = {1.0e9, 1.0e9, 1.0e9};

        // Nothing is rejected, so there is nothing to explain.
        assertNull(smFRETTraceHistogram.emptyExplanation(
                run(histogram, min, max), histogram.nSpots));

        // One range that excludes both traces on its own.
        max[smFRETTraceHistogram.FILTER_TOTAL] = 1.0;
        String one = smFRETTraceHistogram.emptyExplanation(
                run(histogram, min, max), histogram.nSpots);
        assertTrue(one.contains("Total (D+A)"), one);
        assertTrue(one.contains("on its own"), one);

        // Two of them, where widening either alone would still leave nothing.
        max[smFRETTraceHistogram.FILTER_DONOR] = 1.0;
        String two = smFRETTraceHistogram.emptyExplanation(
                run(histogram, min, max), histogram.nSpots);
        assertTrue(two.contains("Donor (target)") && two.contains("Total (D+A)"), two);
        assertTrue(two.contains("will not be enough"), two);

        // Neither range excludes everything by itself, but together they leave nothing: the
        // donor cut keeps only the dim trace and the acceptor cut keeps only the bright one.
        double[] lo = {-1.0e9, -1.0e9, -1.0e9};
        double[] hi = {1.0e9, 100.0, 100.0};
        smFRETTraceHistogram.Histogram result = run(histogram, lo, hi);
        assertEquals(0, result.nSpotsUsed, "the intersection should be empty");
        assertEquals(1, result.nRejectedBy[smFRETTraceHistogram.FILTER_DONOR], "one each");
        assertEquals(1, result.nRejectedBy[smFRETTraceHistogram.FILTER_ACCEPTOR], "one each");
        String neither = smFRETTraceHistogram.emptyExplanation(result, histogram.nSpots);
        assertTrue(neither.contains("No single one"), neither);
    }

    /** The per-range rejection counts, as the status line prints them. */
    @Test
    @DisplayName("the status breakdown lists only the ranges that rejected something")
    void theBreakdownSkipsIdleRanges() {
        float[][] donor = {{500.0f, 500.0f}};
        float[][] acceptor = {{10.0f, 10.0f}};
        smFRETTraceHistogram histogram = withTraces(donor, acceptor);

        double[] min = {-1.0e9, -1.0e9, -1.0e9};
        double[] max = {1.0e9, 1.0e9, 1.0e9};
        assertNull(smFRETTraceHistogram.describeRejections(run(histogram, min, max)),
                "an unfiltered histogram should carry no breakdown at all");

        max[smFRETTraceHistogram.FILTER_DONOR] = 1.0;
        String text = smFRETTraceHistogram.describeRejections(run(histogram, min, max));
        assertEquals("donor 1", text,
                "the two open ranges must not appear");
    }

    private static smFRETTraceHistogram.Histogram run(smFRETTraceHistogram histogram,
                                                      double[] min, double[] max) {
        return histogram.computeHistogram(smFRETTraceHistogram.TYPE_TOTAL, 1, 2,
                new smFRETTraceHistogram.Filters(min, max), 8, none());
    }

    /**
     * FRET is plotted over a fixed -0.2 to 1.2 so the noise skirts stay visible, and anything
     * past that is *counted* rather than folded into the end bin. The cast in the bin index
     * truncates toward zero, so a value just below the lower edge would otherwise land in bin 0.
     */
    @Test
    @DisplayName("FRET outside the fixed range is counted, not folded into the end bin")
    void outOfRangeFretIsCounted() {

        // Donor far larger than the acceptor in the negative direction: A/(D+A) below -0.2.
        float[][] donor = {{100.0f}, {100.0f}};
        float[][] acceptor = {{-50.0f}, {50.0f}};
        smFRETTraceHistogram histogram = withTraces(donor, acceptor);

        smFRETTraceHistogram.Histogram result = histogram.computeHistogram(
                smFRETTraceHistogram.TYPE_FRET, 1, 1, smFRETTraceHistogram.Filters.none(), 10, none());

        assertEquals(smFRETTraceHistogram.FRET_MIN, result.lo, 1.0e-12, "the fixed lower edge");
        assertEquals(2, result.nSpotsUsed, "both traces are in range for the intensity filter");
        assertEquals(1, result.nPoints, "only the one inside -0.2 to 1.2 is binned");
        assertEquals(1, result.nOutside, "the other is reported rather than dropped silently");
    }

    @Test
    @DisplayName("a trace whose channels sum to zero is skipped rather than blowing up")
    void aZeroTotalIsSkippedForFret() {
        float[][] donor = {{100.0f}, {50.0f}};
        float[][] acceptor = {{-100.0f}, {50.0f}};
        smFRETTraceHistogram histogram = withTraces(donor, acceptor);

        smFRETTraceHistogram.Histogram result = histogram.computeHistogram(
                smFRETTraceHistogram.TYPE_FRET, 1, 1, smFRETTraceHistogram.Filters.none(), 10, none());

        assertEquals(1, result.nPoints, "the divide by zero trace should not be binned");
        for (int count : result.counts) {
            assertTrue(count <= 1, "no bin should hold an infinity");
        }
    }

    /**
     * The load-bearing consequence of applying corrections inline: the filter bounds are measured
     * on the corrected traces, so they have to be recomputed whenever a correction changes.
     * Restoring the corrections has to restore the bounds *exactly*, or an end of the slider
     * parked at its extreme would start excluding traces it never used to.
     */
    @Test
    @DisplayName("zeroing the corrections restores the filter bounds exactly")
    void zeroingCorrectionsRestoresTheBounds() {
        float[][] donor = {{100.0f, 250.0f, 80.0f}, {30.0f, 45.0f, 900.0f}};
        float[][] acceptor = {{200.0f, 40.0f, 610.0f}, {12.0f, 8.0f, 77.0f}};
        smFRETTraceHistogram histogram = withTraces(donor, acceptor);

        histogram.computeFilterBounds(none());
        double[] originalMin = histogram.filterMin.clone();
        double[] originalMax = histogram.filterMax.clone();

        histogram.computeFilterBounds(new smFRETTraceHistogram.Corrections(37.0, -11.0, 0.13));
        assertTrue(histogram.filterMin[smFRETTraceHistogram.FILTER_DONOR] != originalMin[smFRETTraceHistogram.FILTER_DONOR],
                "a correction should have moved the bounds");

        histogram.computeFilterBounds(none());
        assertArrayEquals(originalMin, histogram.filterMin, 0.0, "minima");
        assertArrayEquals(originalMax, histogram.filterMax, 0.0, "maxima");
    }

    /**
     * The bounds come from the per frame extremes rather than from the trace averages, and that
     * is what guarantees an end parked at its extreme excludes nothing - since the filter itself
     * is a per frame test.
     */
    @Test
    @DisplayName("the bounds are per frame extremes, not per trace averages")
    void boundsArePerFrame() {
        float[][] donor = {{10.0f, 500.0f}};
        float[][] acceptor = {{0.0f, 0.0f}};
        smFRETTraceHistogram histogram = withTraces(donor, acceptor);

        histogram.computeFilterBounds(none());
        assertEquals(10.0, histogram.filterMin[smFRETTraceHistogram.FILTER_DONOR], 1.0e-9);
        assertEquals(500.0, histogram.filterMax[smFRETTraceHistogram.FILTER_DONOR], 1.0e-9,
                "the average of 255 would silently exclude the bright frame");
    }

    @Test
    @DisplayName("the frame interval is inclusive at both ends")
    void theFrameIntervalIsInclusive() {
        float[][] donor = {{0.0f, 10.0f, 20.0f, 30.0f, 40.0f}};
        float[][] acceptor = {{0.0f, 0.0f, 0.0f, 0.0f, 0.0f}};
        smFRETTraceHistogram histogram = withTraces(donor, acceptor);

        smFRETTraceHistogram.Histogram result = histogram.computeHistogram(
                smFRETTraceHistogram.TYPE_DONOR, 2, 4, smFRETTraceHistogram.Filters.none(), 10, none());

        // Frames 2, 3 and 4 hold 10, 20 and 30, so three frames averaging 20.
        assertEquals(20.0, result.lo, 1.0e-9);
    }
}
