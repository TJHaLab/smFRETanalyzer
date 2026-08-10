import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.awt.event.KeyEvent;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

/**
 * The two handle slider, which Swing does not have and which is therefore hand rolled.
 *
 * The keys are not a convenience. A track a few hundred pixels wide cannot address every frame of
 * a 1295 frame movie, so the arrow keys are the only way to land on an exact interval, and PgUp
 * and PgDn are the only way to walk a movie in non overlapping windows. That last one turns on an
 * off by one that is invisible unless you check it: the interval is *inclusive*, so a five frame
 * window 10..14 has to step by five to reach 15..19, and stepping by the width alone would leave
 * frame 14 in both windows.
 *
 * The other property worth pinning is that low <= high is enforced inside the component. update()
 * in the histogram has no crossover handling at all because of it.
 */
class RangeSliderTest {

    private static smFRETTraceHistogram.RangeSlider slider(int minimum, int maximum) {
        return new smFRETTraceHistogram.RangeSlider(minimum, maximum);
    }

    @Test
    @DisplayName("a new slider is open to its full range")
    void startsFullyOpen() {
        smFRETTraceHistogram.RangeSlider slider = slider(1, 30);
        assertEquals(1, slider.getLow());
        assertEquals(30, slider.getHigh());
    }

    @Test
    @DisplayName("low never crosses high")
    void lowNeverCrossesHigh() {
        smFRETTraceHistogram.RangeSlider slider = slider(1, 30);

        slider.setValues(20, 5);
        assertTrue(slider.getLow() <= slider.getHigh(),
                "got " + slider.getLow() + ".." + slider.getHigh());

        slider.setValues(-100, 1000);
        assertEquals(1, slider.getLow(), "clamped to the minimum");
        assertEquals(30, slider.getHigh(), "clamped to the maximum");

        // Collapsing onto a single value is allowed; inverting is not.
        slider.setValues(15, 15);
        assertEquals(15, slider.getLow());
        assertEquals(15, slider.getHigh());
    }

    @Test
    @DisplayName("left and right move the interval by one, keeping its width")
    void arrowsMoveByOne() {
        smFRETTraceHistogram.RangeSlider slider = slider(1, 30);
        slider.setValues(5, 10);

        slider.onKey(KeyEvent.VK_RIGHT);
        assertEquals(6, slider.getLow());
        assertEquals(11, slider.getHigh());

        slider.onKey(KeyEvent.VK_LEFT);
        assertEquals(5, slider.getLow());
        assertEquals(10, slider.getHigh());
    }

    /**
     * The off by one. Successive presses have to tile the movie: each window starts exactly one
     * past where the last one ended, with no frame in two windows and none skipped.
     */
    @Test
    @DisplayName("page down tiles the movie without repeating a frame")
    void pageDownTiles() {
        smFRETTraceHistogram.RangeSlider slider = slider(1, 100);
        slider.setValues(1, 5);

        int previousHigh = slider.getHigh();
        for (int press = 0; press < 5; press++) {
            slider.onKey(KeyEvent.VK_PAGE_DOWN);
            assertEquals(previousHigh + 1, slider.getLow(),
                    "window " + press + " should start one past the last one's end");
            assertEquals(4, slider.getHigh() - slider.getLow(), "the width should not change");
            previousHigh = slider.getHigh();
        }

        slider.onKey(KeyEvent.VK_PAGE_UP);
        assertEquals(previousHigh - 5 - 4, slider.getLow(), "page up is the exact inverse");
    }

    @Test
    @DisplayName("home and end jump to the limits keeping the width")
    void homeAndEndKeepTheWidth() {
        smFRETTraceHistogram.RangeSlider slider = slider(1, 30);
        slider.setValues(10, 19);
        int width = slider.getHigh() - slider.getLow();

        slider.onKey(KeyEvent.VK_HOME);
        assertEquals(1, slider.getLow());
        assertEquals(1 + width, slider.getHigh());

        slider.onKey(KeyEvent.VK_END);
        assertEquals(30 - width, slider.getLow());
        assertEquals(30, slider.getHigh());
    }

    @Test
    @DisplayName("up and down resize rather than move")
    void upAndDownResize() {
        smFRETTraceHistogram.RangeSlider slider = slider(1, 30);
        slider.setValues(10, 15);

        slider.onKey(KeyEvent.VK_UP);
        assertEquals(10, slider.getLow(), "the start should not move");
        assertEquals(16, slider.getHigh());

        slider.onKey(KeyEvent.VK_DOWN);
        assertEquals(15, slider.getHigh());

        // Shrinking stops at a single frame rather than inverting.
        for (int i = 0; i < 20; i++) {
            slider.onKey(KeyEvent.VK_DOWN);
        }
        assertEquals(10, slider.getLow());
        assertEquals(10, slider.getHigh());
    }

    /**
     * Sliding stops at the ends rather than letting the width shrink against them, so walking to
     * the end of a movie and back returns the same window rather than a narrower one.
     */
    @Test
    @DisplayName("sliding into an end preserves the width")
    void slidingStopsAtTheEnds() {
        smFRETTraceHistogram.RangeSlider slider = slider(1, 30);
        slider.setValues(27, 30);

        for (int i = 0; i < 5; i++) {
            slider.onKey(KeyEvent.VK_RIGHT);
        }
        assertEquals(27, slider.getLow(), "should have stopped rather than slid off the end");
        assertEquals(30, slider.getHigh());
        assertEquals(3, slider.getHigh() - slider.getLow(), "width preserved");

        slider.setValues(1, 4);
        for (int i = 0; i < 5; i++) {
            slider.onKey(KeyEvent.VK_LEFT);
        }
        assertEquals(1, slider.getLow());
        assertEquals(4, slider.getHigh());
    }

    /**
     * setRange is what rescaling to a different quantity or a reloaded file goes through, and it
     * has to pull the interval inside the new limits rather than leaving it dangling outside.
     */
    @Test
    @DisplayName("rescaling brings the interval inside the new limits")
    void setRangeClampsTheInterval() {
        smFRETTraceHistogram.RangeSlider slider = slider(1, 1000);
        slider.setValues(400, 900);

        slider.setRange(1, 100);
        assertTrue(slider.getLow() >= 1, "low " + slider.getLow());
        assertTrue(slider.getHigh() <= 100, "high " + slider.getHigh());
        assertTrue(slider.getLow() <= slider.getHigh());
    }

    @Test
    @DisplayName("a degenerate range does not divide by zero")
    void aSingleValueRangeIsSafe() {
        smFRETTraceHistogram.RangeSlider slider = slider(7, 7);
        assertEquals(7, slider.getLow());
        assertEquals(7, slider.getHigh());

        slider.onKey(KeyEvent.VK_RIGHT);
        slider.onKey(KeyEvent.VK_PAGE_DOWN);
        slider.onKey(KeyEvent.VK_END);
        assertEquals(7, slider.getLow());
        assertEquals(7, slider.getHigh());
    }
}
