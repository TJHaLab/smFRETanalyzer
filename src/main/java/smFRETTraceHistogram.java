/*
 * This class plots histograms of the time traces measured by smFRETAnalyzer.
 *
 * Unlike the other plugins in this package it is interactive rather than batch - it opens a
 * window whose histogram is recomputed as the controls are adjusted, so that thresholds can
 * be chosen by eye. It reads the '.h5' file written by smFRETAnalyzer.
 */

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;

import ij.IJ;

import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Menu;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;


// The single-type import of org.scijava.plugin.Menu shadows java.awt.Menu from the wildcard
// import above, so Menu here is the SciJava annotation.
@Plugin(type = Command.class,
        menu = {@Menu(label = "Plugins"),
                @Menu(label = "smFRET"),
                @Menu(label = "smFRET Trace Histograms", weight = 4.0)})
public class smFRETTraceHistogram implements Command {

    @Parameter
    LogService log;

    // The file the pool starts with, not the only one it can hold. A SciJava dialog takes one
    // file, and more are added in the window itself (issue #22) - which also keeps this
    // parameter, and so every macro that sets it, exactly as it was.
    @Parameter(description = "H5 file written by smFRET Time Traces, more can be added in the window",
               label = "Trace H5 file", style = "open")
    File h5File;

    // Histogram types, indices match the order of the type radio buttons.
    static final int TYPE_FRET = 0;
    static final int TYPE_DONOR = 1;
    static final int TYPE_ACCEPTOR = 2;
    static final int TYPE_TOTAL = 3;
    private static final String[] TYPE_NAMES = {"FRET efficiency", "Donor (target)", "Acceptor (source)", "Total (D+A)"};

    // The three quantities an intensity range can be set on. One slider each, all three live at
    // once and ANDed (issue #8) - they used to be a combo box beside a single slider, which made
    // them mutually exclusive for no reason other than the widget. Filtering on donor brightness
    // and on total brightness are not alternatives, and wanting both is the ordinary case.
    static final int FILTER_TOTAL = 0;
    static final int FILTER_DONOR = 1;
    static final int FILTER_ACCEPTOR = 2;
    static final int N_FILTERS = 3;
    static final String[] FILTER_NAMES = {"Total (D+A)", "Donor (target)", "Acceptor (source)"};

    // For the status line, where all three appear at once and the parenthesised halves would
    // crowd out the counts they are there to label.
    static final String[] SHORT_FILTER_NAMES = {"total", "donor", "acceptor"};

    // FRET efficiency is plotted over a fixed range, slightly wider than [0,1] so that the
    // noise skirts either side of the physical range stay visible.
    static final double FRET_MIN = -0.2;
    static final double FRET_MAX = 1.2;

    // The traces as smFRETAnalyzer wrote them, no corrections applied.
    private static final Corrections NO_CORRECTIONS = new Corrections(0.0, 0.0, 0.0);

    // Wide enough for a typical movie name at the default font, and fixed rather than packed to
    // the longest name so that adding one long name does not shove the histogram sideways.
    private static final int POOL_WIDTH = 210;

    static final String WINDOW_TITLE = "smFRET Trace Histograms";

    // Member variables.
    private JSpinner acceptorBaselineSpinner;
    private JSlider binsSlider;
    private JSpinner donorBaselineSpinner;
    final double[] filterMax = new double[N_FILTERS];
    final double[] filterMin = new double[N_FILTERS];
    private RangeSlider frameRangeSlider;
    private final boolean isHeadless = GraphicsEnvironment.isHeadless();
    private JSpinner leakageSpinner;
    private JButton saveCsvButton;
    private JButton savePngButton;
    private final RangeSlider[] valueRangeSliders = new RangeSlider[N_FILTERS];
    int nFrames = 0;
    int nSpots = 0;
    private HistogramPanel plotPanel;

    // Every file loaded, ticked or not, in name order. The combined matrices above are
    // rebuilt from whichever of these are ticked (issue #22).
    final java.util.List<TraceFile> pool = new ArrayList<>();
    private JPanel poolPanel;
    private Histogram result;
    float[][] sourceTraces;      // [spot][frame], acceptor.
    private JLabel statusLabel;
    private boolean suspendUpdates = false;
    float[][] targetTraces;      // [spot][frame], donor.
    private JRadioButton[] typeButtons;

    /**
     * The result of binning the traces, everything the plot panel needs to draw itself.
     */
    static class Histogram {
        double binWidth;
        int[] counts;
        double lo;
        int maxCount;
        int nOutside;       // Traces dropped for falling outside [lo,hi].
        int nPoints;        // Traces actually binned.
        int nSpotsUsed;     // Traces inside every intensity range.

        // Per filter, how many traces that filter would have rejected. A trace failing two of
        // them counts in both, so these do not sum to the number rejected - that is the point,
        // since it is what says whether widening one slider alone would bring anything back.
        int[] nRejectedBy = new int[N_FILTERS];
        String valueLabel;
    }

    /**
     * The per frame corrections applied to the traces before anything is measured from them.
     *
     * Both baselines are *subtracted*, so a negative value adds one. Leakage is the fraction of
     * the donor that appears in the acceptor channel, and it is taken off the acceptor using the
     * *baseline corrected* donor: leakage is a fraction of real donor emission, and an offset
     * that survived into the trace is not donor emission.
     *
     * Grouped rather than passed as three more doubles because computeHistogram() already takes
     * seven positional arguments, and three unlabelled doubles in a row are easy to transpose.
     */
    static class Corrections {

        final double acceptorBaseline;
        final double donorBaseline;
        final double leakage;

        Corrections(double donorBaseline, double acceptorBaseline, double leakage) {
            this.acceptorBaseline = acceptorBaseline;
            this.donorBaseline = donorBaseline;
            this.leakage = leakage;
        }

        /**
         * Takes the *corrected* donor rather than the raw one, both because that is the quantity
         * the leakage is a fraction of and so that callers do not correct the donor twice.
         */
        double correctAcceptor(double rawAcceptor, double correctedDonor) {
            return rawAcceptor - acceptorBaseline - leakage * correctedDonor;
        }

        double correctDonor(double rawDonor) {
            return rawDonor - donorBaseline;
        }

        String describe() {
            return "baseline D " + compact(donorBaseline) + " / A " + compact(acceptorBaseline)
                    + ", leakage " + compact(leakage);
        }

        boolean isIdentity() {
            return (acceptorBaseline == 0.0) && (donorBaseline == 0.0) && (leakage == 0.0);
        }
    }

    /**
     * The three intensity ranges, one per quantity. A trace has to satisfy all of them.
     *
     * Grouped for the same reason Corrections is: this replaced a filter index and one pair of
     * limits, and passing three pairs positionally would take computeHistogram() to twelve
     * arguments, six of them adjacent unlabelled doubles.
     *
     * Indexed by the FILTER_* constants rather than named per channel, so that the per filter
     * accounting is an array walk and filterValue() can be reused as it stands.
     */
    static final class Filters {

        final double[] min;
        final double[] max;

        Filters(double[] min, double[] max) {
            this.min = min.clone();
            this.max = max.clone();
        }

        /** Everything open, which is what the sliders read before anyone touches them. */
        static Filters none() {
            double[] min = new double[N_FILTERS];
            double[] max = new double[N_FILTERS];
            java.util.Arrays.fill(min, -Double.MAX_VALUE);
            java.util.Arrays.fill(max, Double.MAX_VALUE);
            return new Filters(min, max);
        }

        /** One filter open and the rest closed down to it - the shape the old single one had. */
        static Filters only(int filterType, double low, double high) {
            Filters filters = none();
            filters.min[filterType] = low;
            filters.max[filterType] = high;
            return filters;
        }

        /** Whether this quantity's range excludes the value, for one frame. */
        boolean excludes(int filterType, double value) {
            return (value < min[filterType]) || (value > max[filterType]);
        }
    }

    /**
     * Bin the loaded traces. Takes its settings as arguments rather than reading the controls
     * directly so that the binning can be exercised without a GUI.
     */
    Histogram computeHistogram(int type, int firstFrame, int lastFrame,
                               Filters filters, int nBins,
                               Corrections corrections) {

        // One point per trace, the average over the selected interval. For FRET the donor and
        // acceptor are averaged first and the ratio taken from those averages - averaging the
        // per frame ratios instead would not give the same answer.
        double[] values = new double[nSpots];
        int nValues = 0;
        int nSpotsUsed = 0;
        int nIntervalFrames = lastFrame - firstFrame + 1;
        int[] nRejectedBy = new int[N_FILTERS];

        for (int i = 0; i < nSpots; i++) {
            double donorSum = 0.0;
            double acceptorSum = 0.0;

            // Which of the three ranges this trace fell out of, if any. Recorded per filter
            // rather than as one boolean because a trace that only just fails is the thing the
            // user needs to find, and with three sliders live the question is always which one.
            boolean[] rejectedBy = new boolean[N_FILTERS];

            for (int t = firstFrame - 1; t < lastFrame; t++) {
                double frameDonor = corrections.correctDonor(targetTraces[i][t]);
                double frameAcceptor = corrections.correctAcceptor(sourceTraces[i][t], frameDonor);
                donorSum += frameDonor;
                acceptorSum += frameAcceptor;

                // The whole trace goes if any single frame in the interval falls outside any of
                // the ranges, so a molecule that bleaches part way through contributes nothing
                // rather than a diluted average. The maximum works the same way by design, which
                // makes it strict: one bright frame is enough to drop a trace. That is what
                // catches an aggregate, and it also means a single spike will do it.
                //
                // Not short circuited on the first failure - every filter that would have
                // rejected this trace is recorded, so the accounting can say that widening one
                // slider alone would not bring it back.
                for (int f = 0; f < N_FILTERS; f++) {
                    double frameValue = filterValue(f, frameDonor, frameAcceptor,
                            frameDonor + frameAcceptor);
                    if (filters.excludes(f, frameValue)) {
                        rejectedBy[f] = true;
                    }
                }
            }

            boolean rejected = false;
            for (int f = 0; f < N_FILTERS; f++) {
                if (rejectedBy[f]) {
                    nRejectedBy[f] += 1;
                    rejected = true;
                }
            }
            if (rejected) {
                continue;
            }

            double donor = donorSum / nIntervalFrames;
            double acceptor = acceptorSum / nIntervalFrames;
            double total = donor + acceptor;

            double value;
            if (type == TYPE_FRET) {
                // A near zero total makes the ratio meaningless, not just noisy.
                if (Math.abs(total) < 1.0e-9) {
                    continue;
                }
                value = acceptor / total;
            } else if (type == TYPE_DONOR) {
                value = donor;
            } else if (type == TYPE_ACCEPTOR) {
                value = acceptor;
            } else {
                value = total;
            }

            values[nValues++] = value;
            nSpotsUsed += 1;
        }

        Histogram hist = new Histogram();
        hist.counts = new int[nBins];
        hist.nSpotsUsed = nSpotsUsed;
        hist.nRejectedBy = nRejectedBy;
        hist.valueLabel = TYPE_NAMES[type];

        // Fixed range for FRET efficiency, auto range for the intensity histograms.
        double lo;
        double hi;
        if (type == TYPE_FRET) {
            lo = FRET_MIN;
            hi = FRET_MAX;
        } else {
            lo = Double.MAX_VALUE;
            hi = -Double.MAX_VALUE;
            for (int i = 0; i < nValues; i++) {
                if (values[i] < lo) { lo = values[i]; }
                if (values[i] > hi) { hi = values[i]; }
            }
            if (nValues == 0) {
                lo = 0.0;
                hi = 1.0;
            }
        }
        if (hi <= lo) {
            hi = lo + 1.0;
        }

        hist.lo = lo;
        hist.binWidth = (hi - lo) / nBins;

        double hiEdge = lo + hist.binWidth * nBins;
        for (int i = 0; i < nValues; i++) {

            // Range test the value rather than the bin index. A cast truncates toward zero, so a
            // value just below lo gives bin 0 and would be silently folded into the first bin
            // instead of being counted as out of range.
            if ((values[i] < lo) || (values[i] > hiEdge)) {
                hist.nOutside += 1;
                continue;
            }

            int bin = (int) ((values[i] - lo) / hist.binWidth);

            // The largest value lands one past the last bin, keep it rather than dropping it.
            if (bin >= nBins) {
                bin = nBins - 1;
            }
            hist.counts[bin] += 1;
            hist.nPoints += 1;
        }

        for (int count : hist.counts) {
            if (count > hist.maxCount) {
                hist.maxCount = count;
            }
        }

        return hist;
    }

    /**
     * Why nothing survived the three ranges, or null when something did.
     *
     * Three live filters can be set to an empty intersection, and the issue that asked for them
     * left what to do about it open. Preventing it is the wrong answer: whether an intersection
     * is empty depends on the traces rather than on the numbers, so the limits would have to be
     * recomputed on every drag and the sliders would fight the user's hand. Explaining it costs
     * nothing and answers the question they actually have, which is which of three sliders to
     * move.
     *
     * A filter that rejected every trace is one that would still reject every trace with the
     * other two wide open, so widening it is necessary. When several did, all of them are named
     * - moving one would not be enough. When none did, the ranges are individually survivable
     * and only their intersection is empty, which is worth saying outright because no single
     * slider looks guilty.
     */
    /**
     * How many traces each range rejected, or null when none of them rejected anything.
     *
     * A trace failing two ranges is counted against both, so these do not sum to the number
     * dropped. That is the useful behaviour: it says whether opening one slider would recover
     * anything, which a single total could not.
     */
    static String describeRejections(Histogram result) {
        StringBuilder out = new StringBuilder();
        for (int f : new int[] {FILTER_DONOR, FILTER_ACCEPTOR, FILTER_TOTAL}) {
            if (result.nRejectedBy[f] <= 0) {
                continue;
            }
            if (out.length() > 0) {
                out.append(", ");
            }
            out.append(SHORT_FILTER_NAMES[f]).append(' ').append(String.format("%,d", result.nRejectedBy[f]));
        }
        return (out.length() == 0) ? null : out.toString();
    }

    static String emptyExplanation(Histogram result, int nSpots) {
        if ((nSpots <= 0) || (result.nSpotsUsed > 0)) {
            return null;
        }

        java.util.List<String> culprits = new java.util.ArrayList<>();
        for (int f : new int[] {FILTER_DONOR, FILTER_ACCEPTOR, FILTER_TOTAL}) {
            if (result.nRejectedBy[f] >= nSpots) {
                culprits.add(FILTER_NAMES[f]);
            }
        }

        if (culprits.isEmpty()) {
            return "No traces pass all three ranges. No single one excludes everything, so it is"
                    + " their overlap that is empty - widen whichever matters least.";
        }
        if (culprits.size() == 1) {
            return "No traces pass all three ranges. " + culprits.get(0)
                    + " rejects every trace on its own - widen it.";
        }
        return "No traces pass all three ranges. " + String.join(" and ", culprits)
                + " each reject every trace on their own - widening one will not be enough.";
    }

    /**
     * The intensity the range slider is currently applied to.
     */
    static double filterValue(int filterType, double donor, double acceptor, double total) {
        if (filterType == FILTER_DONOR) {
            return donor;
        }
        if (filterType == FILTER_ACCEPTOR) {
            return acceptor;
        }
        return total;
    }

    /**
     * Panel that draws the current histogram.
     */
    private class HistogramPanel extends JPanel {

        private static final int MARGIN_BOTTOM = 46;
        private static final int MARGIN_LEFT = 66;
        private static final int MARGIN_RIGHT = 18;
        private static final int MARGIN_TOP = 18;

        HistogramPanel() {
            setPreferredSize(new Dimension(660, 340));
            setBackground(Color.WHITE);
        }

        /**
         * The empty-set explanation, centred and wrapped to the plot area.
         *
         * Wrapped by hand because the sentence names filters whose labels vary in length and the
         * window is resizable, so there is no width at which one line is safe.
         */
        private void drawEmptyMessage(Graphics2D g2, String message, int plotWidth, int plotHeight) {
            g2.setColor(new Color(150, 60, 40));
            FontMetrics metrics = g2.getFontMetrics();

            java.util.List<String> lines = new java.util.ArrayList<>();
            StringBuilder line = new StringBuilder();
            for (String word : message.split(" ")) {
                String candidate = (line.length() == 0) ? word : line + " " + word;
                if ((metrics.stringWidth(candidate) > plotWidth) && (line.length() > 0)) {
                    lines.add(line.toString());
                    line = new StringBuilder(word);
                } else {
                    line = new StringBuilder(candidate);
                }
            }
            if (line.length() > 0) {
                lines.add(line.toString());
            }

            int lineHeight = metrics.getHeight();
            int y = MARGIN_TOP + (plotHeight - lines.size() * lineHeight) / 2 + metrics.getAscent();
            for (String text : lines) {
                g2.drawString(text, MARGIN_LEFT + (plotWidth - metrics.stringWidth(text)) / 2, y);
                y += lineHeight;
            }
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            if (result == null) {
                return;
            }

            Graphics2D g2 = (Graphics2D) g.create();
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            g2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

            int plotWidth = getWidth() - MARGIN_LEFT - MARGIN_RIGHT;
            int plotHeight = getHeight() - MARGIN_TOP - MARGIN_BOTTOM;
            if ((plotWidth < 10) || (plotHeight < 10)) {
                g2.dispose();
                return;
            }

            // An empty plot is where the user looks first, so the reason goes here rather than
            // only in the status line under it. Drawn instead of the axes: empty axes invite the
            // reading that the data is wrong, when what is wrong is a slider.
            if (nSpots == 0) {
                drawEmptyMessage(g2, "No files are ticked. Tick one in the list on the left to"
                        + " plot it.", plotWidth, plotHeight);
                g2.dispose();
                return;
            }

            String empty = emptyExplanation(result, nSpots);
            if (empty != null) {
                drawEmptyMessage(g2, empty, plotWidth, plotHeight);
                g2.dispose();
                return;
            }

            int nBins = result.counts.length;
            double yScale = (result.maxCount > 0) ? ((double) plotHeight / (double) result.maxCount) : 0.0;

            // Bars.
            g2.setColor(new Color(70, 115, 175));
            for (int i = 0; i < nBins; i++) {
                if (result.counts[i] == 0) {
                    continue;
                }
                int x0 = MARGIN_LEFT + (int) Math.round((double) i * plotWidth / nBins);
                int x1 = MARGIN_LEFT + (int) Math.round((double) (i + 1) * plotWidth / nBins);
                int h = (int) Math.round(result.counts[i] * yScale);
                int barWidth = Math.max(1, x1 - x0 - 1);
                g2.fillRect(x0, MARGIN_TOP + plotHeight - h, barWidth, h);
            }

            // Axes.
            g2.setColor(Color.DARK_GRAY);
            g2.drawLine(MARGIN_LEFT, MARGIN_TOP + plotHeight, MARGIN_LEFT + plotWidth, MARGIN_TOP + plotHeight);
            g2.drawLine(MARGIN_LEFT, MARGIN_TOP, MARGIN_LEFT, MARGIN_TOP + plotHeight);

            FontMetrics fm = g2.getFontMetrics();

            // X ticks, placed on round multiples of a step rather than at fixed fractions of the
            // axis. For the fixed FRET range this works out as -0.2, 0.0, 0.2 ... 1.2.
            double hi = result.lo + result.binWidth * nBins;
            double range = hi - result.lo;
            double step = niceTickStep(range);
            int decimals = tickDecimals(step);
            double eps = step * 1.0e-6;

            for (int k = (int) Math.ceil(result.lo / step - eps); ; k++) {
                double value = k * step;
                if (value > (hi + eps)) {
                    break;
                }

                // Multiplying out can leave a tiny residue where the tick should be exactly zero.
                if (Math.abs(value) < eps) {
                    value = 0.0;
                }

                int x = MARGIN_LEFT + (int) Math.round((value - result.lo) / range * plotWidth);
                g2.drawLine(x, MARGIN_TOP + plotHeight, x, MARGIN_TOP + plotHeight + 4);
                String label = formatTick(value, decimals);
                g2.drawString(label, x - fm.stringWidth(label) / 2, MARGIN_TOP + plotHeight + 18);
            }

            // Y ticks.
            for (int i = 0; i <= 4; i++) {
                double frac = i / 4.0;
                int y = MARGIN_TOP + plotHeight - (int) Math.round(frac * plotHeight);
                g2.drawLine(MARGIN_LEFT - 4, y, MARGIN_LEFT, y);
                String label = Integer.toString((int) Math.round(frac * result.maxCount));
                g2.drawString(label, MARGIN_LEFT - 8 - fm.stringWidth(label), y + fm.getAscent() / 2 - 1);
            }

            // Axis titles.
            String xTitle = result.valueLabel;
            g2.drawString(xTitle, MARGIN_LEFT + (plotWidth - fm.stringWidth(xTitle)) / 2, getHeight() - 8);

            Graphics2D g2r = (Graphics2D) g2.create();
            g2r.rotate(-Math.PI / 2.0, 16, MARGIN_TOP + plotHeight / 2.0);
            g2r.drawString("counts", 16 - fm.stringWidth("counts") / 2, MARGIN_TOP + plotHeight / 2.0f);
            g2r.dispose();

            g2.dispose();
        }
    }

    /**
     * Two handle slider for choosing the frame interval to average over.
     *
     * Swing has no range slider. One ordinary slider per end can express the same interval, but
     * makes the common operation - moving a window of fixed width through the movie - awkward,
     * because both ends have to be dragged separately and kept in step. Here dragging either
     * handle resizes the interval and dragging the bar between them slides it at constant width.
     *
     * Only low <= high states are reachable, so callers do not have to order the two values.
     */
    static class RangeSlider extends JComponent {

        private static final int BAR_HEIGHT = 6;
        private static final int THUMB_SIZE = 13;

        private static final int DRAG_NONE = 0;
        private static final int DRAG_LOW = 1;
        private static final int DRAG_HIGH = 2;
        private static final int DRAG_BAR = 3;

        private int dragMode = DRAG_NONE;
        private int grabHigh;               // Interval at the start of a bar drag, so that the
        private int grabLow;                // width is preserved exactly however far it is
        private int grabValue;              // dragged, rather than drifting by a rounding error.
        private int high;
        private final ArrayList<ChangeListener> listeners = new ArrayList<>();
        private int low;
        private int maximum;
        private int minimum;
        private boolean valueIsAdjusting;

        RangeSlider(int minimum, int maximum) {
            this.minimum = minimum;
            this.maximum = Math.max(maximum, minimum);
            low = this.minimum;
            high = this.maximum;

            setPreferredSize(new Dimension(200, THUMB_SIZE + 8));
            setFocusable(true);

            // The keys are worth more than the mouse on a long movie - a track a few hundred
            // pixels wide cannot address every frame of it - and there is no way to discover
            // them, nor that the slider has to be clicked first to take the focus, other than
            // being told here. HTML because a tooltip does not otherwise wrap.
            //
            // Worded in steps rather than frames: both the frame interval and the intensity
            // range are this same component.
            setToolTipText("<html>Drag either end to resize the interval, "
                    + "drag the middle to slide it.<br>"
                    + "Click first, then:<br>"
                    + "&nbsp;&nbsp;<b>\u2190</b> <b>\u2192</b> move the interval by one step<br>"
                    + "&nbsp;&nbsp;<b>PgUp</b> <b>PgDn</b> move it a whole interval<br>"
                    + "&nbsp;&nbsp;<b>Home</b> <b>End</b> jump to the start or the end<br>"
                    + "&nbsp;&nbsp;<b>\u2191</b> <b>\u2193</b> resize it</html>");

            MouseAdapter mouse = new MouseAdapter() {
                @Override
                public void mousePressed(MouseEvent e) {
                    requestFocusInWindow();
                    onPress(e.getX());
                }

                @Override
                public void mouseDragged(MouseEvent e) {
                    onDrag(e.getX());
                }

                @Override
                public void mouseReleased(MouseEvent e) {
                    dragMode = DRAG_NONE;

                    // Fire once more on release so listeners that only act on settled values
                    // see the final interval.
                    if (valueIsAdjusting) {
                        valueIsAdjusting = false;
                        fireChange();
                    }
                }
            };
            addMouseListener(mouse);
            addMouseMotionListener(mouse);

            addKeyListener(new KeyAdapter() {
                @Override
                public void keyPressed(KeyEvent e) {
                    onKey(e.getKeyCode());
                }
            });

            addFocusListener(new FocusAdapter() {
                @Override
                public void focusGained(FocusEvent e) {
                    repaint();
                }

                @Override
                public void focusLost(FocusEvent e) {
                    repaint();
                }
            });
        }

        void addChangeListener(ChangeListener listener) {
            listeners.add(listener);
        }

        private int clamp(int value) {
            if (value < minimum) {
                return minimum;
            }
            if (value > maximum) {
                return maximum;
            }
            return value;
        }

        private void fireChange() {
            ChangeEvent event = new ChangeEvent(this);
            for (ChangeListener listener : listeners) {
                listener.stateChanged(event);
            }
        }

        int getHigh() {
            return high;
        }

        int getLow() {
            return low;
        }

        int getMaximum() {
            return maximum;
        }

        int getMinimum() {
            return minimum;
        }

        private void onDrag(int x) {
            if (dragMode == DRAG_NONE) {
                return;
            }
            valueIsAdjusting = true;
            int value = xToValue(x);

            if (dragMode == DRAG_LOW) {
                setValues(Math.min(value, high), high);
            } else if (dragMode == DRAG_HIGH) {
                setValues(low, Math.max(value, low));
            } else {
                slide(grabLow + (value - grabValue), grabHigh - grabLow);
            }
        }

        void onKey(int keyCode) {
            int width = high - low;
            if (keyCode == KeyEvent.VK_LEFT) {
                slide(low - 1, width);
            } else if (keyCode == KeyEvent.VK_RIGHT) {
                slide(low + 1, width);
            } else if (keyCode == KeyEvent.VK_PAGE_UP) {

                // A whole interval per press rather than a fixed number of frames, so each press
                // lands on the next window that does not overlap this one. A fixed step would be
                // either useless on a 1295 frame movie or far too coarse on a 30 frame one; this
                // scales with whatever is being looked at.
                //
                // width + 1, not width: the interval is inclusive, so 10..14 is five frames and
                // stepping by four would leave frame 14 in both windows.
                slide(low - (width + 1), width);
            } else if (keyCode == KeyEvent.VK_PAGE_DOWN) {
                slide(low + (width + 1), width);
            } else if (keyCode == KeyEvent.VK_HOME) {
                slide(minimum, width);
            } else if (keyCode == KeyEvent.VK_END) {
                slide(maximum - width, width);
            } else if (keyCode == KeyEvent.VK_UP) {
                setValues(low, high + 1);
            } else if (keyCode == KeyEvent.VK_DOWN) {
                setValues(low, Math.max(high - 1, low));
            }
        }

        private void onPress(int x) {
            int lowX = valueToX(low);
            int highX = valueToX(high);
            int lowDistance = Math.abs(x - lowX);
            int highDistance = Math.abs(x - highX);

            if ((lowDistance <= THUMB_SIZE) || (highDistance <= THUMB_SIZE)) {

                // Which handle was grabbed. When the interval is empty the two coincide, and the
                // side of the handle that was clicked decides - otherwise a collapsed interval
                // could only ever be reopened in one direction.
                if ((lowDistance < highDistance) || ((lowDistance == highDistance) && (x <= lowX))) {
                    dragMode = DRAG_LOW;
                } else {
                    dragMode = DRAG_HIGH;
                }
            } else if ((x > lowX) && (x < highX)) {
                dragMode = DRAG_BAR;
                grabValue = xToValue(x);
                grabLow = low;
                grabHigh = high;
            } else {

                // On the track outside the interval, bring the nearer end out to meet the click.
                dragMode = (xToValue(x) < low) ? DRAG_LOW : DRAG_HIGH;
            }
            onDrag(x);
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);

            Graphics2D g2 = (Graphics2D) g.create();
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

            int barY = (getHeight() - BAR_HEIGHT) / 2;
            int lowX = valueToX(low);
            int highX = valueToX(high);

            g2.setColor(new Color(205, 205, 205));
            g2.fillRoundRect(trackX(), barY, trackWidth(), BAR_HEIGHT, BAR_HEIGHT, BAR_HEIGHT);

            // The selected interval, in the same blue as the histogram bars.
            g2.setColor(new Color(70, 115, 175));
            g2.fillRoundRect(lowX, barY, Math.max(1, highX - lowX), BAR_HEIGHT, BAR_HEIGHT, BAR_HEIGHT);

            paintThumb(g2, lowX);
            paintThumb(g2, highX);

            g2.dispose();
        }

        private void paintThumb(Graphics2D g2, int x) {
            int y = (getHeight() - THUMB_SIZE) / 2;
            g2.setColor(Color.WHITE);
            g2.fillOval(x - THUMB_SIZE / 2, y, THUMB_SIZE, THUMB_SIZE);
            g2.setColor(isFocusOwner() ? new Color(40, 80, 140) : new Color(90, 90, 90));
            g2.drawOval(x - THUMB_SIZE / 2, y, THUMB_SIZE - 1, THUMB_SIZE - 1);
        }

        /**
         * Set the limits, keeping the interval inside them.
         */
        void setRange(int newMinimum, int newMaximum) {
            minimum = newMinimum;
            maximum = Math.max(newMaximum, newMinimum);
            setValues(low, high);
            repaint();
        }

        /**
         * Set the interval, clamped to the limits and to low <= high.
         */
        void setValues(int newLow, int newHigh) {
            int clampedLow = clamp(newLow);
            int clampedHigh = Math.max(clamp(newHigh), clampedLow);
            if ((clampedLow == low) && (clampedHigh == high)) {
                return;
            }

            low = clampedLow;
            high = clampedHigh;
            repaint();
            fireChange();
        }

        /**
         * Move the interval to start at newLow without changing its width, stopping at the ends
         * rather than letting the width shrink against them.
         */
        private void slide(int newLow, int width) {
            int clampedLow = Math.min(Math.max(newLow, minimum), maximum - width);
            setValues(clampedLow, clampedLow + width);
        }

        // Half a handle of padding at each end, so that the handles stay inside the component
        // when the interval is at either limit.
        private int trackWidth() {
            return Math.max(1, getWidth() - THUMB_SIZE);
        }

        private int trackX() {
            return THUMB_SIZE / 2;
        }

        private int valueToX(int value) {
            if (maximum == minimum) {
                return trackX();
            }
            return trackX() + (int) Math.round((double) (value - minimum) * trackWidth() / (maximum - minimum));
        }

        private int xToValue(int x) {
            if (maximum == minimum) {
                return minimum;
            }
            double fraction = (double) (x - trackX()) / trackWidth();
            return clamp(minimum + (int) Math.round(fraction * (maximum - minimum)));
        }
    }

    /**
     * A round tick spacing (1, 2, 2.5 or 5 times a power of ten) for an axis of this range. The
     * FRET range of 1.4 gives 0.2.
     */
    private static double niceTickStep(double range) {
        if (!(range > 0.0)) {
            return 1.0;
        }

        double raw = range / 7.0;
        double magnitude = Math.pow(10.0, Math.floor(Math.log10(raw)));
        double normalized = raw / magnitude;

        double step;
        if (normalized <= 1.0) {
            step = 1.0;
        } else if (normalized <= 2.0) {
            step = 2.0;
        } else if (normalized <= 2.5) {
            step = 2.5;
        } else if (normalized <= 5.0) {
            step = 5.0;
        } else {
            step = 10.0;
        }
        return step * magnitude;
    }

    /**
     * How many decimals a tick label needs to distinguish one step from the next.
     */
    private static int tickDecimals(double step) {
        if (step >= 1.0) {
            return 0;
        }
        if (step >= 0.1) {
            return 1;
        }
        if (step >= 0.01) {
            return 2;
        }
        return 3;
    }

    /**
     * Short tick labels, the intensity histograms can run to large values.
     */
    private static String formatTick(double value, int decimals) {
        if (Math.abs(value) >= 100000.0) {
            return String.format("%.1e", value);
        }
        return String.format("%." + decimals + "f", value);
    }

    /**
     * A correction value for the status line and the CSV header, without the trailing zeroes.
     *
     * Not DecimalFormat, which would put a comma where a locale wants a decimal separator and
     * make the CSV header lie. %f always emits a decimal point, so stripping trailing zeroes
     * cannot eat a digit of an integer the way it would on a bare "100".
     */
    private static String compact(double value) {
        String text = String.format(java.util.Locale.US, "%.6f", value);
        while (text.endsWith("0")) {
            text = text.substring(0, text.length() - 1);
        }
        return text.endsWith(".") ? text.substring(0, text.length() - 1) : text;
    }

    /**
     * One trace H5 in the pool: its traces, and whether it is currently contributing.
     *
     * The traces are held per file rather than only in the combined matrices so that unticking a
     * file is instant and re-reads nothing. Ticking is the control a user reaches for to ask "is
     * this one movie dragging the distribution?", and a question asked that often should not cost
     * a disk read every time it is asked.
     */
    static final class TraceFile {

        final File file;
        boolean included = true;
        final float[][] source;      // [spot][frame], acceptor.
        final float[][] target;      // [spot][frame], donor.

        TraceFile(File file, float[][] target, float[][] source) {
            this.file = file;
            this.source = source;
            this.target = target;
        }

        /** Frames this file can offer, which is what makes it the shortest one or not. */
        int frames() {
            return Math.min(target[0].length, source[0].length);
        }

        String name() {
            return file.getName();
        }

        int spots() {
            return target.length;
        }
    }

    /**
     * The traces of several files laid end to end, as one pooled experiment.
     *
     * nFrames is the *shortest* file's rather than the longest: one frame range slider governs
     * every trace, so they all have to span the same interval for it to mean anything, and
     * truncating the long files is the only way to get that without inventing data. The surplus
     * frames of a longer movie are simply never looked at - every per spot loop is bounded by
     * nFrames, so no read runs off the end of a short file's row.
     */
    static final class Combined {

        final int nFrames;
        final float[][] source;
        final float[][] target;

        Combined(float[][] target, float[][] source, int nFrames) {
            this.nFrames = nFrames;
            this.source = source;
            this.target = target;
        }

        int nSpots() {
            return target.length;
        }
    }

    /**
     * Concatenate the traces of the files that are ticked, skipping the rest.
     *
     * Static, and given its input rather than reading the pool, so that the combining rule - which
     * is the whole of what issue #22 asked for - can be tested without an H5 file or a window.
     */
    static Combined combine(java.util.List<TraceFile> files) {
        int nSpots = 0;
        int nFrames = Integer.MAX_VALUE;
        for (TraceFile traces : files) {
            if (traces.included) {
                nSpots += traces.spots();
                nFrames = Math.min(nFrames, traces.frames());
            }
        }
        if (nSpots == 0) {
            return new Combined(new float[0][], new float[0][], 0);
        }

        float[][] target = new float[nSpots][];
        float[][] source = new float[nSpots][];
        int at = 0;
        for (TraceFile traces : files) {
            if (!traces.included) {
                continue;
            }

            // The rows are shared rather than copied. Nothing downstream writes to a trace - the
            // corrections are applied to each value as it is read - so copying would double what
            // a pool costs in memory and buy nothing.
            for (int i = 0; i < traces.spots(); i++) {
                target[at] = traces.target[i];
                source[at] = traces.source[i];
                at += 1;
            }
        }
        return new Combined(target, source, nFrames);
    }

    /**
     * Read one trace H5, with the checks that make a wrong file say so.
     */
    static TraceFile readTraceFile(File file) {

        // Checked before the HDF5 library is asked to open it, whose complaint about anything else
        // is a library level error with a stack trace and no mention of which file was wrong.
        smFRETFiles.requireHDF5(file);

        float[][] target;
        float[][] source;
        try (IHDF5Reader reader = HDF5Factory.openForReading(file)) {
            target = reader.readFloatMatrix("target-traces");
            source = reader.readFloatMatrix("source-traces");
        }

        if ((target.length == 0) || (source.length == 0)) {
            throw new smFRETAnalysisException("No traces in " + file);
        }
        if (target.length != source.length) {
            throw new smFRETAnalysisException("Target and source trace counts differ ("
                    + target.length + " vs " + source.length + ") in " + file);
        }
        return new TraceFile(file, target, source);
    }

    /**
     * Load one file as the whole pool, which is what opening the plugin does.
     */
    void loadTraces(File file) {
        pool.clear();
        java.util.List<String> problems = addTraceFiles(java.util.Collections.singletonList(file));
        if (!problems.isEmpty()) {

            // The seeding file is the one the user chose in the dialog, so a failure here is not
            // a stray in a multi-file drop and there is nothing left to carry on with.
            throw new smFRETAnalysisException(problems.get(0));
        }
        log.info("loaded " + nSpots + " traces of " + nFrames + " frames from " + file);
    }

    /**
     * Add files to the pool, keeping it in name order and free of duplicates.
     *
     * Returns a line per file that could not be read, empty when all of them were. A bad file in
     * a dropped selection leaves the rest of the selection loaded rather than failing the whole
     * drop: the usual way to get one is to sweep up a stray while selecting the wanted ones, and
     * discarding the good files along with it would be the wrong trade.
     */
    java.util.List<String> addTraceFiles(java.util.List<File> files) {
        java.util.List<String> problems = new java.util.ArrayList<>();
        for (File file : files) {
            if (inPool(file)) {
                continue;
            }
            try {
                pool.add(readTraceFile(file));
            } catch (Exception e) {
                problems.add(file.getName() + " - " + e.getMessage());
            }
        }
        sortPool();
        rebuild();
        return problems;
    }

    /** Whether this file is already pooled, compared by path rather than by name. */
    private boolean inPool(File file) {
        String wanted = pathKey(file);
        for (TraceFile traces : pool) {
            if (pathKey(traces.file).equals(wanted)) {
                return true;
            }
        }
        return false;
    }

    /**
     * The identity of a file, for telling "already pooled" from "another file of the same name".
     *
     * Canonical rather than absolute, which is the whole point: an absolute path keeps whatever
     * spelling it arrived with, so the same file dropped from two places - through a symlinked
     * directory, or with a './' in it - would compare unequal and be pooled twice. That does not
     * look like a mistake in the histogram, it just silently weights one movie double.
     *
     * Falls back to the absolute path when the file system will not answer, which is the best
     * that can be done and no worse than not trying.
     */
    private static String pathKey(File file) {
        try {
            return file.getCanonicalPath();
        } catch (java.io.IOException e) {
            return file.getAbsolutePath();
        }
    }

    /**
     * Name order, which is the order a folder of movies is named in and so the order they are
     * thought about in. Ties break on the full path, because two directories of repeats will
     * often hold the same file names.
     */
    private void sortPool() {
        pool.sort((a, b) -> {
            int byName = a.name().compareToIgnoreCase(b.name());
            return (byName != 0) ? byName : pathKey(a.file).compareTo(pathKey(b.file));
        });
    }

    /**
     * Rebuild the combined traces from whichever files are ticked.
     */
    void rebuild() {
        Combined combined = combine(pool);
        targetTraces = combined.target;
        sourceTraces = combined.source;
        nSpots = combined.nSpots();
        nFrames = combined.nFrames;
        computeFilterBounds(corrections());
    }

    /** The ticked files, in the order they were concatenated. */
    java.util.List<TraceFile> includedFiles() {
        java.util.List<TraceFile> included = new java.util.ArrayList<>();
        for (TraceFile traces : pool) {
            if (traces.included) {
                included.add(traces);
            }
        }
        return included;
    }

    /**
     * The file pool, down the left hand side of the window.
     *
     * A column rather than a row because the pool grows downwards and the histogram beside it
     * wants the width: a folder of twenty repeats is an ordinary thing to pool, and a horizontal
     * list of twenty file names is not readable at any window size.
     *
     * Each file is a tick box, so unticking one takes it out of the histogram without taking it
     * out of the pool. That is the comparison this panel exists for - whether one movie is
     * dragging the distribution - and it is a different question from "I chose the wrong file",
     * which is what the remove button beside it answers.
     */
    private JComponent buildPoolPanel() {
        poolPanel = new JPanel();
        poolPanel.setLayout(new BoxLayout(poolPanel, BoxLayout.Y_AXIS));
        poolPanel.setBackground(Color.WHITE);

        JScrollPane scroll = new JScrollPane(poolPanel,
                JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scroll.setPreferredSize(new Dimension(POOL_WIDTH, 0));
        scroll.getVerticalScrollBar().setUnitIncrement(12);

        // Dropping files anywhere over the column adds them. Once, on the scroll pane, which
        // covers the rows inside it too - see smFRETSwing for why that is enough.
        scroll.setTransferHandler(smFRETSwing.fileDropHandler(this::onFilesAdded));

        // Kept alongside the drop target rather than replacing it: dragging is quicker when a
        // file manager is already open, and unreachable when it is not.
        JButton addButton = new JButton("Add...");
        addButton.setToolTipText("Add trace H5 files to the pool. Files can also be dropped here.");
        addButton.addActionListener(e -> onAdd());

        JLabel heading = new JLabel("Trace files");
        heading.setBorder(new EmptyBorder(0, 2, 2, 2));

        JLabel hint = new JLabel("<html><i>drop files here</i></html>");
        hint.setBorder(new EmptyBorder(2, 2, 0, 2));
        hint.setForeground(Color.GRAY);

        JPanel buttons = new JPanel(new BorderLayout());
        buttons.add(hint, BorderLayout.NORTH);
        buttons.add(addButton, BorderLayout.SOUTH);

        JPanel column = new JPanel(new BorderLayout(0, 2));
        column.setBorder(new EmptyBorder(6, 8, 6, 4));
        column.add(heading, BorderLayout.NORTH);
        column.add(scroll, BorderLayout.CENTER);
        column.add(buttons, BorderLayout.SOUTH);
        return column;
    }

    /**
     * Rebuild the pool rows from the pool.
     *
     * Rebuilt wholesale rather than patched because adding a file re-sorts the list, so the row
     * order after any change is not a function of the row order before it.
     */
    void refreshPoolPanel() {
        if (poolPanel == null) {
            return;
        }
        poolPanel.removeAll();
        for (TraceFile traces : pool) {
            poolPanel.add(buildPoolRow(traces));
        }

        // Holds the rows at the top of the column. Without it BoxLayout spreads them down the
        // whole height and a pool of two files sits with a gap between them.
        poolPanel.add(Box.createVerticalGlue());
        poolPanel.revalidate();
        poolPanel.repaint();
    }

    /**
     * One row of the pool: tick box, name, trace count and a remove button.
     */
    private JComponent buildPoolRow(TraceFile traces) {
        JCheckBox tick = new JCheckBox(traces.name(), traces.included);
        tick.setBackground(Color.WHITE);

        // The name is often too long for the column, and the part that identifies a movie is
        // usually its tail rather than its head. The tooltip carries the full path and the shape
        // of the file, which is what says why one of them is limiting the frame range.
        tick.setToolTipText(smFRETSwing.dialogMessage(pathKey(traces.file) + "\n"
                + String.format("%,d traces of %,d frames", traces.spots(), traces.frames())));
        tick.addActionListener(e -> {
            traces.included = tick.isSelected();
            onPoolTicked();
        });

        JButton remove = new JButton("×");
        remove.setToolTipText("Remove " + traces.name() + " from the pool");
        remove.setMargin(new Insets(0, 4, 0, 4));
        remove.setFocusable(false);

        // The last file cannot be removed. An empty pool has no directory to open a save dialog
        // in and no name to title a plot with, and "untick it" already covers wanting it out of
        // the histogram - so the button that would empty the pool is the one that is not needed.
        remove.setEnabled(pool.size() > 1);
        remove.addActionListener(e -> {
            pool.remove(traces);
            rebuild();
            onPoolChanged();
        });

        JPanel row = new JPanel(new BorderLayout(2, 0));
        row.setBackground(Color.WHITE);
        row.add(tick, BorderLayout.CENTER);
        row.add(remove, BorderLayout.EAST);

        // BoxLayout otherwise stretches every row to the tallest one it can, which on a short
        // pool makes each name a band the height of the column.
        row.setMaximumSize(new Dimension(Integer.MAX_VALUE, row.getPreferredSize().height));
        row.setAlignmentX(Component.LEFT_ALIGNMENT);
        return row;
    }

    /**
     * Files were added or removed: the pool is materially different, so the sliders start over.
     */
    void onPoolChanged() {
        refreshPoolPanel();
        resetSliderRanges();
        retitle();
        update();
    }

    /**
     * A file was ticked or unticked: the same pool seen a different way, so the filters stay put.
     */
    private void onPoolTicked() {
        rebuild();
        rescaleSliderRanges();
        retitle();
        update();
    }

    /**
     * Prompt for files to add.
     */
    private void onAdd() {
        JFileChooser chooser = new JFileChooser(poolDirectory());
        chooser.setDialogTitle("Select smFRET trace H5 files");
        chooser.setMultiSelectionEnabled(true);
        if (chooser.showOpenDialog(ownerWindow()) != JFileChooser.APPROVE_OPTION) {
            return;
        }
        onFilesAdded(java.util.Arrays.asList(chooser.getSelectedFiles()));
    }

    /**
     * Add files from a drop or the chooser, reporting whichever of them could not be read.
     */
    private void onFilesAdded(java.util.List<File> files) {
        java.util.List<String> problems;
        try {
            problems = addTraceFiles(files);
        } catch (Exception e) {
            log.info(e);
            JOptionPane.showMessageDialog(ownerWindow(),
                    smFRETSwing.dialogMessage("Could not add files:\n" + e.getMessage()),
                    WINDOW_TITLE, JOptionPane.ERROR_MESSAGE);
            return;
        }
        onPoolChanged();

        // After the pool has been updated, not instead of updating it - the files that did read
        // are already in and the window should show them while this is on screen.
        if (!problems.isEmpty()) {
            // Through dialogMessage because these name paths: a plain string would make
            // JOptionPane as wide as the longest of them, which is wider than the screen.
            JOptionPane.showMessageDialog(ownerWindow(),
                    smFRETSwing.dialogMessage("Not added:\n" + String.join("\n", problems)),
                    WINDOW_TITLE, JOptionPane.WARNING_MESSAGE);
        }
    }

    /**
     * The window the pool panel is in, or null before it is in one.
     *
     * Looked up rather than held, which is what lets every control in this panel be built before
     * the frame exists and still address it afterwards. The alternative was threading a JFrame
     * through six methods whose only use for it was setTitle.
     */
    private java.awt.Window ownerWindow() {
        return (poolPanel == null) ? null : SwingUtilities.getWindowAncestor(poolPanel);
    }

    /** Retitle the window for the current selection, if there is a window yet. */
    private void retitle() {
        java.awt.Window window = ownerWindow();
        if (window instanceof JFrame) {
            ((JFrame) window).setTitle(WINDOW_TITLE + " - " + poolLabel());
        }
    }

    /**
     * How the ticked files are named in the window title, the saved plot and the status line.
     */
    String poolLabel() {
        java.util.List<TraceFile> included = includedFiles();
        if (included.isEmpty()) {
            return "no files selected";
        }
        if (included.size() == 1) {
            return included.get(0).name();
        }
        return included.size() + " files";
    }

    /** Where the save and add dialogs open, which is beside whichever file is first. */
    private File poolDirectory() {
        java.util.List<TraceFile> included = includedFiles();
        File first = included.isEmpty() ? pool.get(0).file : included.get(0).file;
        return first.getParentFile();
    }

    /**
     * The path a saved CSV or PNG is offered under.
     *
     * A pool is named after its first file with a suffix rather than after all of them, because
     * the alternative is a file name that grows with the pool - and the files that went into it
     * are written into the CSV header, where they can be read without being in the name.
     */
    private String poolSaveRoot() {
        java.util.List<TraceFile> included = includedFiles();
        if (included.size() == 1) {
            return stripExtension(included.get(0).file);
        }
        File first = included.isEmpty() ? pool.get(0).file : included.get(0).file;
        return stripExtension(first) + "_pool";
    }

    /**
     * Range of each quantity the intensity range slider can be applied to, which sets its limits.
     *
     * Measured on the corrected traces, so it has to be redone whenever a correction changes -
     * subtracting a baseline moves the data out from under bounds taken on the raw traces, and an
     * end of the slider parked at a stale extreme would then silently exclude traces.
     */
    void computeFilterBounds(Corrections corrections) {
        for (int f = 0; f < FILTER_NAMES.length; f++) {
            filterMin[f] = Double.MAX_VALUE;
            filterMax[f] = -Double.MAX_VALUE;
        }
        for (int i = 0; i < nSpots; i++) {
            for (int t = 0; t < nFrames; t++) {
                double donor = corrections.correctDonor(targetTraces[i][t]);
                double acceptor = corrections.correctAcceptor(sourceTraces[i][t], donor);
                for (int f = 0; f < FILTER_NAMES.length; f++) {
                    double value = filterValue(f, donor, acceptor, donor + acceptor);
                    if (value < filterMin[f]) { filterMin[f] = value; }
                    if (value > filterMax[f]) { filterMax[f] = value; }
                }
            }
        }
        for (int f = 0; f < FILTER_NAMES.length; f++) {

            // With nothing ticked the loop above ran zero times and both bounds are still their
            // sentinels, which resetFilterSliderRange would cast to an int and overflow. Nothing
            // is being filtered in that state, so any finite range will do.
            if (nSpots == 0) {
                filterMin[f] = 0.0;
                filterMax[f] = 1.0;
            } else if (filterMax[f] <= filterMin[f]) {
                filterMax[f] = filterMin[f] + 1.0;
            }
        }
    }

    /**
     * The corrections the spinboxes currently ask for.
     *
     * loadTraces() runs before the window is built, and again from Browse... after it is, so this
     * has to answer in both states - uncorrected before there are controls to read.
     */
    private Corrections corrections() {
        if (donorBaselineSpinner == null) {
            return NO_CORRECTIONS;
        }
        return new Corrections(((Number) donorBaselineSpinner.getValue()).doubleValue(),
                ((Number) acceptorBaselineSpinner.getValue()).doubleValue(),
                ((Number) leakageSpinner.getValue()).doubleValue());
    }

    /**
     * A correction changed: the traces have moved, so the intensity range slider is rescaled to
     * them and reopened, exactly as switching the quantity it applies to already does.
     */
    private void onCorrectionChanged() {
        if (suspendUpdates) {
            return;
        }
        computeFilterBounds(corrections());
        resetFilterSliderRange();
        update();
    }

    /**
     * Write the current histogram out as a CSV table.
     */
    private void onSaveCsv(JFrame frame) {
        JFileChooser chooser = new JFileChooser(poolDirectory());
        chooser.setSelectedFile(new File(poolSaveRoot() + "_histogram.csv"));
        if (chooser.showSaveDialog(frame) != JFileChooser.APPROVE_OPTION) {
            return;
        }
        try (PrintWriter writer = new PrintWriter(chooser.getSelectedFile())) {
            writer.println("# " + result.valueLabel + " from " + poolLabel());

            // Every ticked file with its own trace count, which is the only record of what went
            // into a pooled histogram - the file name carries the first one and nothing else.
            for (TraceFile traces : includedFiles()) {
                writer.println("# file: " + traces.file + " (" + traces.spots() + " traces of "
                        + traces.frames() + " frames)");
            }
            writer.println("# frames pooled: " + nFrames + ", the shortest file's");
            // All three ranges, including the ones left wide open. A saved histogram has to say
            // what was filtered, and an omitted line reads as "no filter" rather than "this one
            // was untouched" only if you already know how many there are.
            StringBuilder ranges = new StringBuilder();
            for (int f : new int[] {FILTER_DONOR, FILTER_ACCEPTOR, FILTER_TOTAL}) {
                ranges.append(", ").append(FILTER_NAMES[f]).append(' ')
                        .append(valueRangeSliders[f].getLow()).append('-')
                        .append(valueRangeSliders[f].getHigh());
            }
            writer.println("# frames " + frameRangeSlider.getLow() + "-" + frameRangeSlider.getHigh()
                    + ranges + ", " + result.nPoints + " of " + result.nSpotsUsed
                    + " traces in range");

            // Written even when they are all zero, so that a saved histogram says what was done
            // to the traces rather than leaving it to be inferred from the absence of a line.
            writer.println("# corrections: " + corrections().describe());
            writer.println("bin_center,count");
            for (int i = 0; i < result.counts.length; i++) {
                writer.println((result.lo + (i + 0.5) * result.binWidth) + "," + result.counts[i]);
            }
        } catch (Exception e) {
            log.info(e);
            IJ.handleException(e);
        }
    }

    /**
     * Write the current plot out as a PNG.
     */
    private void onSavePng(JFrame frame) {
        JFileChooser chooser = new JFileChooser(poolDirectory());
        chooser.setSelectedFile(new File(poolSaveRoot() + "_histogram.png"));
        if (chooser.showSaveDialog(frame) != JFileChooser.APPROVE_OPTION) {
            return;
        }
        try {
            ImageIO.write(renderPlotImage(), "png", chooser.getSelectedFile());
        } catch (Exception e) {
            log.info(e);
            IJ.handleException(e);
        }
    }

    /**
     * Render the current plot for saving, titled with the H5 file name.
     *
     * The window itself shows the file name in its header row, but a saved plot travels on its
     * own and would otherwise lose track of which data it came from, so the title is added here
     * rather than being drawn into the panel on screen.
     */
    BufferedImage renderPlotImage() {
        String title = poolLabel();
        int titleHeight = 30;

        BufferedImage image = new BufferedImage(plotPanel.getWidth(),
                plotPanel.getHeight() + titleHeight, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2 = image.createGraphics();
        g2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

        g2.setColor(plotPanel.getBackground());
        g2.fillRect(0, 0, image.getWidth(), titleHeight);

        g2.setColor(Color.DARK_GRAY);
        g2.setFont(g2.getFont().deriveFont(Font.BOLD, 14.0f));
        FontMetrics fm = g2.getFontMetrics();
        g2.drawString(title, (image.getWidth() - fm.stringWidth(title)) / 2,
                (titleHeight + fm.getAscent()) / 2 - 2);

        g2.translate(0, titleHeight);
        plotPanel.paint(g2);
        g2.dispose();

        return image;
    }

    /**
     * Point the frame and intensity sliders at the currently loaded traces.
     */
    private void resetSliderRanges() {
        boolean wasSuspended = suspendUpdates;
        suspendUpdates = true;
        try {
            frameRangeSlider.setRange(1, Math.max(1, nFrames));
            frameRangeSlider.setValues(1, Math.max(1, nFrames));
            resetFilterSliderRange();
        } finally {
            suspendUpdates = wasSuspended;
        }
    }

    /**
     * Re-scale the sliders to the ticked files without discarding the limits the user set.
     *
     * The counterpart to resetSliderRanges(): ticking a file in or out is the same experiment
     * seen a different way, and the question being asked is whether that one movie moves the
     * distribution - which it cannot answer if the filters reset underneath it every time.
     *
     * A handle sitting at its limit is the exception. That is not a number anyone chose, it is
     * the slider saying "no limit", so it follows the limit rather than staying behind at the old
     * value - otherwise ticking in a brighter file would leave an untouched slider quietly
     * excluding its traces. Handles that were moved keep their value, clamped to the new range.
     */
    private void rescaleSliderRanges() {
        boolean wasSuspended = suspendUpdates;
        suspendUpdates = true;
        try {
            boolean framesWereOpen = (frameRangeSlider.getLow() <= frameRangeSlider.getMinimum())
                    && (frameRangeSlider.getHigh() >= frameRangeSlider.getMaximum());
            frameRangeSlider.setRange(1, Math.max(1, nFrames));
            if (framesWereOpen) {
                frameRangeSlider.setValues(1, Math.max(1, nFrames));
            }

            for (int f = 0; f < N_FILTERS; f++) {
                RangeSlider slider = valueRangeSliders[f];
                boolean lowWasOpen = slider.getLow() <= slider.getMinimum();
                boolean highWasOpen = slider.getHigh() >= slider.getMaximum();
                int lo = (int) Math.floor(filterMin[f]);
                int hi = (int) Math.ceil(filterMax[f]);
                slider.setRange(lo, hi);
                slider.setValues(lowWasOpen ? lo : slider.getLow(),
                        highWasOpen ? hi : slider.getHigh());
            }
        } finally {
            suspendUpdates = wasSuspended;
        }
    }

    /**
     * Point the intensity range slider at the range of the currently selected filter quantity.
     * Both handles start at the extremes so that nothing is hidden until the user asks for it.
     *
     * The bounds are the per-frame minimum and maximum rather than the range of the trace
     * averages, which is what guarantees neither end can silently exclude a trace when it is
     * parked at its extreme.
     */
    private void resetFilterSliderRange() {
        boolean wasSuspended = suspendUpdates;
        suspendUpdates = true;
        try {
            for (int f = 0; f < N_FILTERS; f++) {
                int lo = (int) Math.floor(filterMin[f]);
                int hi = (int) Math.ceil(filterMax[f]);
                valueRangeSliders[f].setRange(lo, hi);
                valueRangeSliders[f].setValues(lo, hi);
            }
        } finally {
            suspendUpdates = wasSuspended;
        }
    }

    /** The three ranges as the sliders currently read them. */
    private Filters filters() {
        double[] min = new double[N_FILTERS];
        double[] max = new double[N_FILTERS];
        for (int f = 0; f < N_FILTERS; f++) {
            min[f] = valueRangeSliders[f].getLow();
            max[f] = valueRangeSliders[f].getHigh();
        }
        return new Filters(min, max);
    }

    /**
     * Which histogram type is currently selected.
     */
    private int selectedType() {
        for (int i = 0; i < typeButtons.length; i++) {
            if (typeButtons[i].isSelected()) {
                return i;
            }
        }
        return TYPE_FRET;
    }

    /**
     * File name without its extension, used to suggest save names.
     */
    private static String stripExtension(File file) {
        String name = file.toString();
        int dotIndex = name.lastIndexOf('.');
        if (dotIndex > 0) {
            name = name.substring(0, dotIndex);
        }
        return name;
    }

    /**
     * Recompute the histogram and redraw. Called whenever a control changes.
     */
    private void update() {
        if (suspendUpdates) {
            return;
        }

        Corrections corrections = corrections();
        result = computeHistogram(selectedType(),
                frameRangeSlider.getLow(),
                frameRangeSlider.getHigh(),
                filters(),
                binsSlider.getValue(),
                corrections);

        // The file count only when there is more than one, so a single file reads exactly as it
        // did before the pool existed.
        int nIncluded = includedFiles().size();
        String from = (nIncluded > 1) ? String.format(" from %d files", nIncluded) : "";

        String status = String.format("%,d of %,d traces%s · frames %d-%d (%,d wide)",
                result.nSpotsUsed, nSpots, from,
                frameRangeSlider.getLow(), frameRangeSlider.getHigh(),
                frameRangeSlider.getHigh() - frameRangeSlider.getLow() + 1);
        if (result.nOutside > 0) {
            status += String.format(" · %,d outside range", result.nOutside);
        }

        // What each range is costing, before it costs everything. Only the ones actually
        // rejecting something are listed - three zeroes on an unfiltered histogram would be
        // noise, and the point is to make a slider that is biting stand out.
        String rejected = describeRejections(result);
        if (rejected != null) {
            status += " · rejected: " + rejected;
        }

        // Only when they are doing something, so that the common uncorrected case reads as it
        // did before rather than carrying three zeroes around.
        if (!corrections.isIdentity()) {
            status += " · " + corrections.describe();
        }
        statusLabel.setText(status);

        // A full width row is enough for the usual line, but a narrow window and three biting
        // ranges will still outrun it, and a JLabel truncates to an ellipsis without saying what
        // it dropped. The tooltip is the same text, so nothing is unreachable.
        statusLabel.setToolTipText(status);

        // Nothing ticked means there is no histogram to write, and a PNG of the message that
        // says so is not a plot anybody wanted to save.
        if (saveCsvButton != null) {
            saveCsvButton.setEnabled(nSpots > 0);
            savePngButton.setEnabled(nSpots > 0);
        }

        plotPanel.repaint();
    }

    /**
     * Build the window.
     */
    private void showWindow() {
        JFrame frame = new JFrame(WINDOW_TITLE + " - " + poolLabel());
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

        // Histogram type.
        typeButtons = new JRadioButton[TYPE_NAMES.length];
        ButtonGroup typeGroup = new ButtonGroup();
        JPanel typePanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 12, 2));
        typePanel.add(new JLabel("Histogram:"));
        for (int i = 0; i < TYPE_NAMES.length; i++) {
            typeButtons[i] = new JRadioButton(TYPE_NAMES[i], i == TYPE_FRET);
            typeButtons[i].addActionListener(e -> update());
            typeGroup.add(typeButtons[i]);
            typePanel.add(typeButtons[i]);
        }

        // The pool lives in its own column down the left (issue #22), so the top of the window
        // is the histogram type on its own.
        JComponent poolColumn = buildPoolPanel();

        JPanel topPanel = new JPanel();
        topPanel.setLayout(new BoxLayout(topPanel, BoxLayout.Y_AXIS));
        topPanel.add(typePanel);

        // Plot.
        plotPanel = new HistogramPanel();
        plotPanel.setBorder(new EmptyBorder(6, 6, 6, 6));

        // Controls. The sliders are created before resetSliderRanges() fills in their limits.
        // One point per trace means far fewer points than the old per frame histogram, so the
        // default bin count is correspondingly lower.
        binsSlider = new JSlider(10, 200, 30);
        frameRangeSlider = new RangeSlider(1, Math.max(1, nFrames));

        // One slider per quantity, all three live at once. Their limits are unrelated to each
        // other - a total runs to roughly the sum of the two channels - so each is scaled to its
        // own range and opened to it, which is what makes "parked at the extreme excludes
        // nothing" true of all three independently.
        for (int f = 0; f < N_FILTERS; f++) {
            valueRangeSliders[f] = new RangeSlider(0, 1);
        }

        // Corrections. Spin boxes rather than sliders because these are set to a measured number -
        // a baseline read off a blank region, a leakage measured on a donor only sample - rather
        // than dialled in by eye, and typing the number is the point.
        donorBaselineSpinner = correctionSpinner(new SpinnerNumberModel(0.0, -1.0e9, 1.0e9, 1.0), "0.###",
                "Subtracted from every donor value. Negative adds one.");
        acceptorBaselineSpinner = correctionSpinner(new SpinnerNumberModel(0.0, -1.0e9, 1.0e9, 1.0), "0.###",
                "Subtracted from every acceptor value. Negative adds one.");
        // Four decimals rather than the baselines' three: a leakage coefficient is a small
        // fraction, and the editor commits what it displays, so the format is the precision.
        leakageSpinner = correctionSpinner(new SpinnerNumberModel(0.0, 0.0, 1.0, 0.01), "0.####",
                "Fraction of the baseline corrected donor removed from the acceptor.");

        JPanel correctionPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 6, 1));
        correctionPanel.add(new JLabel("Baseline   donor"));
        correctionPanel.add(donorBaselineSpinner);
        correctionPanel.add(new JLabel("acceptor"));
        correctionPanel.add(acceptorBaselineSpinner);
        correctionPanel.add(Box.createHorizontalStrut(14));
        correctionPanel.add(new JLabel("Donor leakage"));
        correctionPanel.add(leakageSpinner);

        JPanel controlPanel = new JPanel(new GridBagLayout());
        controlPanel.setBorder(new EmptyBorder(4, 8, 4, 8));
        addSliderRow(controlPanel, 0, new JLabel("Bins"), binsSlider, () -> Integer.toString(binsSlider.getValue()));
        addSliderRow(controlPanel, 1, new JLabel("Frames"), frameRangeSlider,
                () -> frameRangeSlider.getLow() + "-" + frameRangeSlider.getHigh());
        // Donor, acceptor and total in that order rather than the FILTER_* order, which puts
        // total first. The two channels are what the plot is of, and the total is derived from
        // them, so reading them in that order matches the rest of the window.
        int row = 2;
        for (int f : new int[] {FILTER_DONOR, FILTER_ACCEPTOR, FILTER_TOTAL}) {
            final RangeSlider slider = valueRangeSliders[f];
            addSliderRow(controlPanel, row++, new JLabel(FILTER_NAMES[f]), slider,
                    () -> slider.getLow() + "-" + slider.getHigh());
        }

        GridBagConstraints correctionConstraints = new GridBagConstraints();
        correctionConstraints.anchor = GridBagConstraints.WEST;
        correctionConstraints.gridwidth = 3;
        correctionConstraints.gridx = 0;
        correctionConstraints.gridy = row;
        correctionConstraints.insets = new Insets(1, 0, 1, 6);
        controlPanel.add(correctionPanel, correctionConstraints);

        // Status and save buttons.
        statusLabel = new JLabel(" ");
        saveCsvButton = new JButton("Save CSV...");
        saveCsvButton.addActionListener(e -> onSaveCsv(frame));
        savePngButton = new JButton("Save PNG...");
        savePngButton.addActionListener(e -> onSavePng(frame));

        // The status line gets the full width of the window, with the buttons on their own row
        // underneath. It used to share a row with them, so the buttons ate the right hand end of
        // it - which was survivable when it read "412 of 500 traces" and is not now that it
        // carries a per range rejection breakdown as well.
        JPanel bottomPanel = new JPanel(new BorderLayout(0, 4));
        bottomPanel.setBorder(new EmptyBorder(2, 10, 8, 10));
        bottomPanel.add(statusLabel, BorderLayout.NORTH);
        JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.RIGHT, 6, 0));
        buttonPanel.add(saveCsvButton);
        buttonPanel.add(savePngButton);
        bottomPanel.add(buttonPanel, BorderLayout.SOUTH);

        JPanel southPanel = new JPanel(new BorderLayout());
        southPanel.add(controlPanel, BorderLayout.CENTER);
        southPanel.add(bottomPanel, BorderLayout.SOUTH);

        frame.getContentPane().setLayout(new BorderLayout());
        frame.getContentPane().add(topPanel, BorderLayout.NORTH);
        frame.getContentPane().add(poolColumn, BorderLayout.WEST);
        frame.getContentPane().add(plotPanel, BorderLayout.CENTER);
        frame.getContentPane().add(southPanel, BorderLayout.SOUTH);

        refreshPoolPanel();
        resetSliderRanges();
        update();

        frame.pack();
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }

    /**
     * showWindow() helper, one correction spin box.
     */
    private JSpinner correctionSpinner(SpinnerNumberModel model, String format, String tip) {
        JSpinner spinner = new JSpinner(model);
        JSpinner.NumberEditor editor = new JSpinner.NumberEditor(spinner, format);
        editor.getTextField().setColumns(6);
        spinner.setEditor(editor);
        spinner.setToolTipText(tip);
        spinner.addChangeListener(e -> onCorrectionChanged());
        return spinner;
    }

    /**
     * showWindow() helper, adds one labelled slider row with a live value readout.
     *
     * JSlider and RangeSlider have the same addChangeListener() method but no common supertype
     * that declares it, hence the pair of overloads over a shared layout helper.
     */
    private void addSliderRow(JPanel parent, int row, JComponent label, JSlider slider,
                              java.util.function.Supplier<String> valueText) {
        JLabel valueLabel = addControlRow(parent, row, label, slider, valueText);
        slider.addChangeListener(e -> {
            valueLabel.setText(valueText.get());
            update();
        });
    }

    /**
     * showWindow() helper, the range slider flavour of addSliderRow().
     */
    private void addSliderRow(JPanel parent, int row, JComponent label, RangeSlider slider,
                              java.util.function.Supplier<String> valueText) {
        JLabel valueLabel = addControlRow(parent, row, label, slider, valueText);
        slider.addChangeListener(e -> {
            valueLabel.setText(valueText.get());
            update();
        });
    }

    /**
     * addSliderRow() helper, lays out one row and returns its value readout.
     */
    private JLabel addControlRow(JPanel parent, int row, JComponent label, JComponent control,
                                 java.util.function.Supplier<String> valueText) {
        GridBagConstraints c = new GridBagConstraints();
        c.gridy = row;
        c.insets = new Insets(1, 2, 1, 6);

        c.gridx = 0;
        c.anchor = GridBagConstraints.WEST;
        parent.add(label, c);

        c.gridx = 1;
        c.fill = GridBagConstraints.HORIZONTAL;
        c.weightx = 1.0;
        parent.add(control, c);

        c.gridx = 2;
        c.fill = GridBagConstraints.NONE;
        c.weightx = 0.0;
        JLabel valueLabel = new JLabel(valueText.get());

        // Wide enough for the frame row's "1234-5678", the longest readout of the three.
        valueLabel.setPreferredSize(new Dimension(78, valueLabel.getPreferredSize().height));
        parent.add(valueLabel, c);

        return valueLabel;
    }

    /**
     * Run ...
     */
    @Override
    public void run() {
        try {
            if (isHeadless) {
                log.info("smFRET Trace Histograms is interactive and cannot run headless");
                return;
            }

            log.info("loading traces from " + h5File);
            loadTraces(h5File);

            SwingUtilities.invokeLater(this::showWindow);

        } catch (smFRETAnalysisException e) {

            // This plugin's own validation, so the message is the whole of what is worth showing.
            smFRETFiles.report(log, e);
        } catch (Exception e) {
            log.info(e);
            IJ.handleException(e);
        }
    }
}
