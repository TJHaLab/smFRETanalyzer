import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.Arrays;
import java.util.Random;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

/**
 * The quickselect behind every background median.
 *
 * It replaced Arrays.sort because it is called twice per clipping round over every trusted pixel
 * and was 14% of the profile on a long movie. The whole argument for that swap is that selection
 * returns exactly the same value a sort would, so that is what these check - against a sorted
 * reference, over the input shapes a naive pivot goes quadratic on.
 *
 * The other property tested here is easy to lose and would be silent: median() leaves the same
 * multiset behind, only reordered. robustSpread relies on it, reusing the scratch array it just
 * took a median of to hold the absolute deviations.
 */
class SpotFinderSelectionTest {

    /** What median() has to agree with. */
    private static double sortedMedian(float[] values, int count) {
        float[] sorted = Arrays.copyOf(values, count);
        Arrays.sort(sorted);
        if ((count % 2) != 0) {
            return sorted[count / 2];
        }
        return 0.5 * ((double) sorted[count / 2 - 1] + (double) sorted[count / 2]);
    }

    @ParameterizedTest(name = "count={0}")
    @ValueSource(ints = {1, 2, 3, 4, 5, 8, 9, 100, 101, 1000, 1001})
    @DisplayName("median matches a sort, at both parities")
    void medianMatchesSort(int count) {
        Random random = new Random(count * 7919L);
        for (int trial = 0; trial < 40; trial++) {
            float[] values = new float[count];
            for (int i = 0; i < count; i++) {
                values[i] = (float) (random.nextGaussian() * 50.0);
            }
            assertEquals(sortedMedian(values, count), smFRETSpotFinder.median(values.clone(), count),
                    1.0e-9, "count " + count + " trial " + trial);
        }
    }

    /**
     * The shapes a median-of-three pivot exists to survive. A naive pivot goes quadratic on
     * sorted and reverse-sorted input, and heavy duplication is what breaks a partition that
     * does not handle equal keys.
     */
    @Test
    @DisplayName("median matches a sort on the adversarial shapes")
    void medianMatchesSortOnDegenerateInput() {
        int[] counts = {1, 2, 3, 16, 17, 999, 1000};
        for (int count : counts) {
            float[] ascending = new float[count];
            float[] descending = new float[count];
            float[] identical = new float[count];
            float[] twoValued = new float[count];
            for (int i = 0; i < count; i++) {
                ascending[i] = i;
                descending[i] = count - i;
                identical[i] = 42.0f;
                twoValued[i] = ((i % 2) == 0) ? 1.0f : 2.0f;
            }

            for (float[] values : new float[][] {ascending, descending, identical, twoValued}) {
                assertEquals(sortedMedian(values, count),
                        smFRETSpotFinder.median(values.clone(), count), 1.0e-9,
                        "count " + count + " of " + Arrays.toString(
                                Arrays.copyOf(values, Math.min(count, 6))));
            }
        }
    }

    @Test
    @DisplayName("select returns every order statistic")
    void selectReturnsEveryOrderStatistic() {
        Random random = new Random(20260810L);
        for (int trial = 0; trial < 200; trial++) {
            int count = 1 + random.nextInt(60);
            float[] values = new float[count];
            for (int i = 0; i < count; i++) {

                // A small range against the count, so duplicates are common rather than incidental.
                values[i] = random.nextInt(10);
            }
            float[] sorted = Arrays.copyOf(values, count);
            Arrays.sort(sorted);

            for (int k = 0; k < count; k++) {
                assertEquals(sorted[k], smFRETSpotFinder.select(values.clone(), count, k), 0.0f,
                        "k " + k + " of " + count);
            }
        }
    }

    /**
     * robustSpread takes a median of the residuals in scratch, then overwrites those same entries
     * with their absolute deviations. That is only correct if the median left the values there -
     * reordered is fine, lost or duplicated is not.
     */
    @Test
    @DisplayName("median preserves the multiset it was given")
    void medianPreservesTheMultiset() {
        Random random = new Random(4242L);
        for (int trial = 0; trial < 100; trial++) {
            int count = 1 + random.nextInt(200);

            // Longer than count, to check the tail past count is untouched as well.
            float[] values = new float[count + 17];
            for (int i = 0; i < values.length; i++) {
                values[i] = random.nextInt(50);
            }
            float[] original = values.clone();

            smFRETSpotFinder.median(values, count);

            float[] beforeHead = Arrays.copyOf(original, count);
            float[] afterHead = Arrays.copyOf(values, count);
            Arrays.sort(beforeHead);
            Arrays.sort(afterHead);
            assertArrayEquals(beforeHead, afterHead, 0.0f, "the first " + count + " entries");

            assertArrayEquals(Arrays.copyOfRange(original, count, original.length),
                    Arrays.copyOfRange(values, count, values.length), 0.0f,
                    "the entries past count");
        }
    }

    /**
     * The prominence ring's summary. Interpolated between neighbours, so the answer generally is
     * not one of the inputs - which is the part a naive nearest-rank implementation gets wrong.
     */
    @Test
    @DisplayName("percentile90 interpolates between the neighbouring order statistics")
    void percentile90Interpolates() {

        // 11 values, so 0.9*(11-1) = 9.0 lands exactly on an order statistic.
        float[] exact = new float[11];
        for (int i = 0; i < 11; i++) {
            exact[i] = i;
        }
        assertEquals(9.0, smFRETSpotFinder.percentile90(exact), 1.0e-9);

        // 6 values, so 0.9*5 = 4.5 lands halfway between the 5th and 6th.
        float[] halfway = {0.0f, 1.0f, 2.0f, 3.0f, 10.0f, 20.0f};
        assertEquals(15.0, smFRETSpotFinder.percentile90(halfway), 1.0e-9);

        assertEquals(7.0, smFRETSpotFinder.percentile90(new float[] {7.0f}), 1.0e-9);
    }

    @Test
    @DisplayName("percentile90 does not reorder its argument")
    void percentile90LeavesItsArgumentAlone() {
        float[] values = {5.0f, 1.0f, 9.0f, 3.0f, 7.0f};
        float[] original = values.clone();
        smFRETSpotFinder.percentile90(values);
        assertArrayEquals(original, values, 0.0f);
    }
}
