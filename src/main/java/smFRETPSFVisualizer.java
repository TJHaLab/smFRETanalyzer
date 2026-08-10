/*
 * This class measures the point spread function in both channels and shows it four ways: the
 * donor and acceptor PSFs as images, and their radial profiles against a fitted aberrated Airy
 * pattern.
 *
 * Interactive, and writes nothing unless asked. Its input is the spot finder JSON, like
 * smFRETTraceVisualizer - but unlike that one it does not need the trace H5, because a PSF is
 * measured from the field rather than from the traces. It therefore runs straight after stage 2.
 *
 * The measurement itself is in smFRETPSF, which has no Swing and no ImageJ in it so that it can
 * be scored against a field of known PSF. What is here is loading, the panels and the controls.
 */

import ij.ImagePlus;
import ij.process.ImageProcessor;

import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Menu;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.PrintWriter;
import java.util.Map;


// The single-type import of org.scijava.plugin.Menu shadows java.awt.Menu from the wildcard
// import above, so Menu here is the SciJava annotation.
@Plugin(type = Command.class,
        menu = {@Menu(label = "Plugins"),
                @Menu(label = "smFRET"),
                @Menu(label = "smFRET PSF Visualizer", weight = 6.0)})
public class smFRETPSFVisualizer implements Command {

    @Parameter
    LogService log;

    @Parameter(description = "JSON file written by smFRET Spot Finder", label = "Spot Finder JSON file", style = "open")
    File spotJSONFile;

    // The suffix smFRETSpotFinder appends to its root name. Stripping it recovers the root, which
    // is what the other outputs are named from - the recorded "root name" is relative to whatever
    // directory stage 2 ran in and does not survive being opened from anywhere else.
    private static final String JSON_SUFFIX = "_spotf_finding.json";

    private static final Color ACCEPTOR_COLOR = new Color(200, 60, 40);
    private static final Color DONOR_COLOR = new Color(30, 140, 60);
    private static final Color FIT_COLOR = new Color(70, 115, 175);
    private static final Color RAW_COLOR = new Color(170, 170, 170);

    // Both plots are logarithmic, because on a linear axis everything past the core is a flat
    // line along zero and the entire question - how much light is out in the wings - becomes
    // invisible. The floor is chosen from the data rather than fixed: on the example data the
    // profile bottoms out around 0.024, so a fixed 1e-4 floor left half the axis empty.
    private static final double FLOOR_LOWEST = 1.0e-6;
    private static final double FLOOR_HIGHEST = 1.0e-1;

    // What the fit is, said once where it is being looked at.
    private static final String MODEL_DESCRIPTION =
            "Data fit with an Airy disk model with primary spherical aberration";

    // How much trace SNR may be given up before a spotSigma is called worse than the best one.
    // The maximum of a matched filter is quadratically flat, so quoting it alone would suggest a
    // precision that is not there - on the example data anything from 1.35 to 2.25 is inside this.
    private static final double FILTER_TOLERANCE = 0.02;

    // Member variables.
    private JSpinner binsSpinner;
    private smFRETPSF.Measurement[] channel = new smFRETPSF.Measurement[2];
    private final boolean isHeadless = GraphicsEnvironment.isHeadless();
    private int firstFrame = 1;
    private JFrame frame;
    private int lastFrame = 1;
    private ImagePlus movie;
    private JSpinner maskSpinner;
    private JSpinner patchSpinner;
    private JLabel adviceLabel;
    private JCheckBox pedestalBox;
    private ProfilePanel[] profilePanel = new ProfilePanel[2];
    private PsfImagePanel[] imagePanel = new PsfImagePanel[2];
    private String root;
    private final smFRETChannelMapper smfcm = new smFRETChannelMapper();
    private final smFRETSpotFinder smfsf = new smFRETSpotFinder();
    private double[][] spots;
    private double spotSigma = 2.0;
    private JPanel plotPanels;
    private JLabel statusLabel;
    private boolean suspendUpdates = false;

    // The averaged and split halves. Built once: the frame range is fixed at whatever spot
    // finding used, so nothing the user can change from here affects them, and re-measuring after
    // a patch or mask change costs only the extraction.
    private float[][] cachedHalf = new float[2][];
    private int cachedWidth;
    private int cachedHeight;

    static final int DONOR = 0;
    static final int ACCEPTOR = 1;
    private static final String[] CHANNEL_NAMES = {"Donor (target)", "Acceptor (source)"};

    /**
     * Everything the plugin needs, found from the spot finder JSON.
     *
     * No trace H5 is looked for. A PSF is measured from the field, so this can run as soon as
     * spot finding has - which is when it is most useful, since the PSF is what tells you whether
     * spotSigma was a sensible choice in the first place.
     */
    void load(File jsonFile) {
        Map<String, Object> mapping = smFRETFiles.readSpotFinderJSON(jsonFile);

        String path = jsonFile.getAbsolutePath();
        root = path.endsWith(JSON_SUFFIX)
                ? path.substring(0, path.length() - JSON_SUFFIX.length())
                : (String) mapping.get("root name");
        File jsonDir = jsonFile.getAbsoluteFile().getParentFile();

        spotSigma = ((Number) mapping.get("spot sigma")).doubleValue();

        // The same frames spot finding averaged. Not offered as a control: a PSF is a property of
        // the optics rather than of the movie, so averaging a different stretch of the same field
        // is a way to get a noisier version of the same answer. Taking the recorded range also
        // means the PSF is measured from the very image the spots were found in.
        firstFrame = intOr(mapping.get("start slice"), 1);
        lastFrame = intOr(mapping.get("end slice"), Integer.MAX_VALUE);

        File spotsFile = locate((String) mapping.get("spots file"), jsonDir, root + "_spotf_spots.csv");
        File imageFile = locate((String) mapping.get("image name"), jsonDir, null);
        File mappingFile = locate((String) mapping.get("mapping file"), jsonDir, null);
        for (File needed : new File[] {spotsFile, imageFile, mappingFile}) {
            if ((needed == null) || !needed.exists()) {
                throw new smFRETAnalysisException("Could not find " + needed
                        + ", named by " + jsonFile);
            }
        }

        // The mapper does the image work - averaging, splitting and the warp - and the spot
        // finder is here only to read its own spot table, which knows the column layout.
        smfsf.log = log;
        smfcm.log = log;
        smfcm.loadMappingJSON(mappingFile.toString());
        spots = smfsf.loadSpotLocations(spotsFile.toString());

        movie = smFRETFiles.openImage(imageFile, "the image");
        movie = smFRETChannelMapper.toFloat(movie);

        log.info("loaded " + spots.length + " spots from " + root);
    }

    /** A recorded integer, or a default when an older spot finder did not write it. */
    private static int intOr(Object value, int fallback) {
        return (value instanceof Number) ? ((Number) value).intValue() : fallback;
    }

    /**
     * The recorded path if it is still there, otherwise the same file beside the JSON.
     */
    private static File locate(String recorded, File jsonDir, String derived) {
        if (recorded != null) {
            File asRecorded = new File(recorded);
            if (asRecorded.exists()) {
                return asRecorded;
            }
            File beside = new File(jsonDir, new File(recorded).getName());
            if (beside.exists()) {
                return beside;
            }
        }
        return (derived == null) ? new File(String.valueOf(recorded)) : new File(derived);
    }

    /**
     * Average the chosen frames and split them, unless that has already been done for this range.
     *
     * The acceptor is measured on the *warped* half, the same one smFRETAnalyzer measures traces
     * from. That is a real choice with a cost: warping interpolates, so what comes back is the
     * acceptor PSF convolved with the interpolation kernel and its core reads slightly broader
     * than it is. Measuring on the unwarped half instead would need the spot centres inverse
     * mapped into the source frame, and would then be reporting a PSF at positions nothing else
     * in the pipeline uses. Consistency with how the acceptor is actually measured is worth more
     * here than the last few percent of core width.
     */
    private void buildHalves() {
        if (cachedHalf[DONOR] != null) {
            return;
        }

        // averageImagePlus clamps the range to the movie, so an end slice past the last frame -
        // which is what the Integer.MAX_VALUE fallback is - simply means "to the end".
        ImagePlus averaged = smfcm.averageImagePlus(movie, firstFrame, lastFrame);
        java.util.List<ImagePlus> halves = smfcm.splitImagePlus(averaged, true);

        for (int c = 0; c < 2; c++) {
            ImageProcessor processor = halves.get(c).getProcessor();
            cachedWidth = processor.getWidth();
            cachedHeight = processor.getHeight();
            cachedHalf[c] = new float[cachedWidth * cachedHeight];
            for (int y = 0; y < cachedHeight; y++) {
                for (int x = 0; x < cachedWidth; x++) {
                    cachedHalf[c][y * cachedWidth + x] = processor.getf(x, y);
                }
            }
        }
    }

    /**
     * Re-measure both channels with the current settings and redraw.
     */
    void update() {
        if (suspendUpdates) {
            return;
        }
        try {
            if (frame != null) {
                frame.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
            }

            int patch = ((Number) patchSpinner.getValue()).intValue();
            double mask = ((Number) maskSpinner.getValue()).doubleValue();
            int bins = ((Number) binsSpinner.getValue()).intValue();

            buildHalves();
            for (int c = 0; c < 2; c++) {
                channel[c] = smFRETPSF.analyse(
                        smFRETPSF.extract(cachedHalf[c], cachedWidth, cachedHeight, spots,
                                patch, mask),
                        bins);
                imagePanel[c].repaint();
                profilePanel[c].repaint();
            }
            updateAdvice();
            statusLabel.setText(describe());
        } catch (smFRETAnalysisException e) {
            smFRETFiles.report(log, e);
        } catch (Exception e) {
            log.info(e);
            statusLabel.setText("measurement failed: " + e.getMessage());
        } finally {
            if (frame != null) {
                frame.setCursor(Cursor.getDefaultCursor());
            }
        }
    }

    /** The status line: what was measured, and what the fit says. */
    private String describe() {
        smFRETPSF.Measurement donor = channel[DONOR];
        StringBuilder text = new StringBuilder();
        text.append(String.format("%,d of %,d spots used · frames %d-%d",
                donor.samples.spotsUsed, donor.samples.spotsTotal,
                firstFrame, Math.min(lastFrame, smFRETChannelMapper.frameCount(movie))));

        for (int c = 0; c < 2; c++) {
            smFRETPSF.Fit fit = fitOf(c);
            if (fit == null) {
                text.append(" · ").append(c == DONOR ? "D" : "A").append(" no fit");
                continue;
            }
            text.append(String.format(" · %s sigma %.2f px, %.2f waves",
                    (c == DONOR) ? "D" : "A", fit.sigma, fit.waves));
            if (pedestalBox.isSelected()) {
                text.append(String.format(", pedestal %.3f", fit.pedestal));
            }
            text.append(String.format(", rms %.3f", fit.rms));
        }

        if (!channel[DONOR].borderSettled) {
            text.append(" · border correction did not settle, treat the wings with care");
        }
        return text.toString();
    }

    /**
     * The bottom of a log axis, chosen so the data fills it.
     *
     * The decade at or below the smallest value being drawn, clamped so a single stray point near
     * zero cannot stretch the axis to nothing. A fixed floor is what left the plots using half
     * their height: the profile on the example data bottoms out around 0.024, two decades above
     * the 1e-4 the axis used to start at.
     */
    private static double floorFor(double smallest) {
        if (!(smallest > 0.0) || Double.isNaN(smallest)) {
            return FLOOR_HIGHEST;
        }
        double decade = Math.pow(10.0, Math.floor(Math.log10(smallest)));
        return Math.min(FLOOR_HIGHEST, Math.max(FLOOR_LOWEST, decade));
    }

    /**
     * What spotSigma would extract the most trace SNR from the PSF that was just measured.
     *
     * smFRETAnalyzer measures a trace with a Gaussian matched filter, and the PSF is not a
     * Gaussian - an aberrated Airy holds a real fraction of its light out where a Gaussian of the
     * same core has none, so a filter matched to the core alone throws that away. The best width
     * is therefore not the fitted core, and how far off it is depends on the aberration: barely
     * at all below 0.3 waves, half again as wide by 0.5.
     *
     * Only the PSF term goes into this, never the fitted pedestal. A flat halo is background, and
     * putting it in would make the objective grow with the filter area and push the answer up for
     * a reason that has nothing to do with the molecule.
     *
     * spotSigma is one parameter for both channels, so the recommendation is the joint one - the
     * geometric mean of the two channels' SNR, each over its own best. The per channel answers go
     * beside it because when they disagree that is worth seeing.
     */
    private void updateAdvice() {
        if ((adviceLabel == null) || (channel[DONOR] == null) || (channel[ACCEPTOR] == null)) {
            return;
        }

        smFRETPSF.Fit donor = fitOf(DONOR);
        smFRETPSF.Fit acceptor = fitOf(ACCEPTOR);
        if ((donor == null) || (acceptor == null)) {
            adviceLabel.setText("Not enough of a profile to recommend a spotSigma");
            return;
        }

        double[] donorProfile = smFRETPSF.snrProfile(donor);
        double[] acceptorProfile = smFRETPSF.snrProfile(acceptor);
        smFRETPSF.FilterOptimum joint = smFRETPSF.optimalFilter(
                new double[][] {donorProfile, acceptorProfile}, FILTER_TOLERANCE);
        smFRETPSF.FilterOptimum forDonor = smFRETPSF.optimalFilter(
                new double[][] {donorProfile}, FILTER_TOLERANCE);
        smFRETPSF.FilterOptimum forAcceptor = smFRETPSF.optimalFilter(
                new double[][] {acceptorProfile}, FILTER_TOLERANCE);

        adviceLabel.setText(String.format(
                "Best SpotSigma for trace SNR: %.2f \u00b7 within %.0f%% of best over %.2f to %.2f"
                        + " \u00b7 donor %.2f, acceptor %.2f",
                joint.best, 100.0 * FILTER_TOLERANCE, joint.low, joint.high,
                forDonor.best, forAcceptor.best));
    }

    /** Whichever of the two fits the pedestal checkbox is asking for. */
    private smFRETPSF.Fit fitOf(int c) {
        if (channel[c] == null) {
            return null;
        }
        return pedestalBox.isSelected() ? channel[c].withPedestal : channel[c].psfOnly;
    }

    /**
     * The measured PSF as an image, on a log scale.
     *
     * Log because the wings are the point: on a linear stretch a pattern holding a percent of its
     * peak at eight pixels is indistinguishable from one holding nothing there, which is exactly
     * the distinction this plugin exists to show.
     *
     * Pixels no spot could contribute to - masked out in every patch that covered them - are
     * drawn in a flat colour rather than as black, because "nothing measured here" and "no light
     * here" are very different statements on a crowded field.
     */
    private class PsfImagePanel extends JPanel {

        private final int index;

        PsfImagePanel(int index) {
            this.index = index;
            setPreferredSize(new Dimension(300, 300));
            setBackground(Color.WHITE);
            setBorder(new EmptyBorder(6, 6, 6, 6));
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            smFRETPSF.Measurement measurement = channel[index];
            if (measurement == null) {
                return;
            }

            Graphics2D g2 = (Graphics2D) g.create();
            g2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING,
                    RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

            int size = measurement.samples.size;

            // Border corrected, so this and the profile below it are the same measurement on the
            // same footing. Without it the outer pixels read near zero by construction - the
            // border median having been subtracted out of every patch - and a few go negative.
            double[] image = measurement.correctedImage();

            double smallest = Double.MAX_VALUE;
            int unmeasured = 0;
            for (double value : image) {
                if (Double.isNaN(value)) {
                    unmeasured++;
                } else if (value > 0.0) {
                    smallest = Math.min(smallest, value);
                }
            }
            double floor = floorFor(smallest);

            // Room for the title above and the scale caption below. Leaving the caption out of
            // this is what clipped it off the bottom of the panel.
            int titleHeight = 20;
            int captionHeight = 18;
            int side = Math.min(getWidth() - 12, getHeight() - 12 - titleHeight - captionHeight);
            if (side < 20) {
                g2.dispose();
                return;
            }
            int left = (getWidth() - side) / 2;
            int top = titleHeight + ((getHeight() - titleHeight - side) / 2);

            BufferedImage rendered = new BufferedImage(size, size, BufferedImage.TYPE_INT_RGB);
            double logFloor = Math.log10(floor);
            for (int y = 0; y < size; y++) {
                for (int x = 0; x < size; x++) {
                    double value = image[y * size + x];
                    int rgb;
                    if (Double.isNaN(value)) {

                        // A colour no part of the greyscale ramp reaches, so an unmeasured pixel
                        // cannot be mistaken for a dark one.
                        rgb = new Color(60, 70, 110).getRGB();
                    } else {
                        double scaled = (Math.log10(Math.max(value, floor)) - logFloor)
                                / (0.0 - logFloor);
                        int level = (int) Math.round(255.0 * Math.min(1.0, Math.max(0.0, scaled)));
                        rgb = new Color(level, level, level).getRGB();
                    }
                    rendered.setRGB(x, y, rgb);
                }
            }

            g2.setColor(Color.DARK_GRAY);
            g2.setFont(g2.getFont().deriveFont(Font.BOLD, 12.0f));
            g2.drawString(CHANNEL_NAMES[index] + " PSF", 8, 14);

            g2.drawImage(rendered, left, top, side, side, null);
            g2.setColor(new Color(120, 120, 120));
            g2.drawRect(left, top, side, side);

            g2.setFont(g2.getFont().deriveFont(Font.PLAIN, 10.0f));
            g2.setColor(Color.GRAY);
            // The blue is only mentioned when there is some. On a field of a few hundred spots
            // every offset in the patch picks up a contribution from something and there are
            // none at all, so a permanent note about them would be a permanent puzzle.
            String caption = String.format("%d px across, log scale %.0e to 1", size, floor);
            if (unmeasured > 0) {
                caption += String.format(" · %d px unmeasured (blue)", unmeasured);
            }
            g2.drawString(caption, left, top + side + 12);

            g2.dispose();
        }
    }

    /**
     * The radial profile against the fit, on a log axis.
     */
    private class ProfilePanel extends JPanel {

        private static final int MARGIN_BOTTOM = 34;
        private static final int MARGIN_LEFT = 52;
        private static final int MARGIN_RIGHT = 12;
        private static final int MARGIN_TOP = 24;

        private double floor = FLOOR_HIGHEST;
        private final int index;

        ProfilePanel(int index) {
            this.index = index;
            setPreferredSize(new Dimension(300, 260));
            setBackground(Color.WHITE);
        }

        private int xFor(double radius, double maxRadius, int width) {
            return MARGIN_LEFT + (int) Math.round((radius / maxRadius) * width);
        }

        private int yFor(double value, int height) {
            double logFloor = Math.log10(floor);
            double scaled = (Math.log10(Math.max(value, floor)) - logFloor) / (0.0 - logFloor);
            return MARGIN_TOP + height - (int) Math.round(scaled * height);
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            smFRETPSF.Measurement measurement = channel[index];
            if (measurement == null) {
                return;
            }

            Graphics2D g2 = (Graphics2D) g.create();
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            g2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING,
                    RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

            int width = getWidth() - MARGIN_LEFT - MARGIN_RIGHT;
            int height = getHeight() - MARGIN_TOP - MARGIN_BOTTOM;
            if ((width < 40) || (height < 40)) {
                g2.dispose();
                return;
            }

            double maxRadius = measurement.samples.patch;
            boolean[] usable = measurement.usable();

            // The axis is scaled to what is actually drawn - both the measured points and the
            // fitted curve, since a pedestal below the last point would otherwise fall off the
            // bottom of a plot that claims to show it.
            smFRETPSF.Fit fit = fitOf(index);
            double smallest = Double.MAX_VALUE;
            for (int i = 0; i < usable.length; i++) {
                if (usable[i]) {
                    smallest = Math.min(smallest, Math.min(measurement.binProfile[i],
                            measurement.binRaw[i]));
                }
            }
            if (fit != null) {
                double[] ends = fit.at(new double[] {maxRadius});
                smallest = Math.min(smallest, ends[0]);
                if (pedestalBox.isSelected() && (fit.pedestal > 0.0)) {
                    smallest = Math.min(smallest, fit.pedestal);
                }
            }
            floor = floorFor(smallest);

            // Decade gridlines, which on a log axis are the only tick positions that mean
            // anything.
            g2.setFont(g2.getFont().deriveFont(Font.PLAIN, 10.0f));
            for (double decade = 1.0; decade >= floor; decade /= 10.0) {
                int y = yFor(decade, height);
                g2.setColor(new Color(232, 232, 232));
                g2.drawLine(MARGIN_LEFT, y, MARGIN_LEFT + width, y);
                g2.setColor(Color.GRAY);
                String label = (decade >= 1.0) ? "1" : String.format("%.0e", decade);
                g2.drawString(label, 8, y + 4);
            }
            g2.setColor(new Color(120, 120, 120));
            g2.drawLine(MARGIN_LEFT, MARGIN_TOP, MARGIN_LEFT, MARGIN_TOP + height);
            g2.drawLine(MARGIN_LEFT, MARGIN_TOP + height, MARGIN_LEFT + width, MARGIN_TOP + height);

            for (int r = 0; r <= (int) maxRadius; r += 2) {
                int x = xFor(r, maxRadius, width);
                g2.setColor(Color.GRAY);
                g2.drawLine(x, MARGIN_TOP + height, x, MARGIN_TOP + height + 3);
                g2.drawString(String.valueOf(r), x - 4, MARGIN_TOP + height + 15);
            }
            g2.drawString("radius (pixels)", MARGIN_LEFT + (width / 2) - 34,
                    MARGIN_TOP + height + 29);

            // The uncorrected points behind the corrected ones, so it is visible where the border
            // correction mattered and where it did not.
            for (int i = 0; i < measurement.binCentre.length; i++) {
                if (!usable[i]) {
                    continue;
                }
                int x = xFor(measurement.binCentre[i], maxRadius, width);
                g2.setColor(RAW_COLOR);
                g2.drawOval(x - 3, yFor(measurement.binRaw[i], height) - 3, 6, 6);
            }

            if (fit != null) {
                double[] fine = new double[120];
                for (int i = 0; i < fine.length; i++) {
                    fine[i] = (maxRadius * i) / (fine.length - 1);
                }
                double[] curve = fit.at(fine);

                g2.setColor(FIT_COLOR);
                g2.setStroke(new BasicStroke(1.6f));
                for (int i = 1; i < fine.length; i++) {
                    g2.drawLine(xFor(fine[i - 1], maxRadius, width), yFor(curve[i - 1], height),
                            xFor(fine[i], maxRadius, width), yFor(curve[i], height));
                }
                g2.setStroke(new BasicStroke(1.0f));

                if (pedestalBox.isSelected() && (fit.pedestal > floor)) {
                    int y = yFor(fit.pedestal, height);
                    g2.setColor(new Color(FIT_COLOR.getRed(), FIT_COLOR.getGreen(),
                            FIT_COLOR.getBlue(), 120));
                    g2.drawLine(MARGIN_LEFT, y, MARGIN_LEFT + width, y);
                }
            }

            Color measured = (index == DONOR) ? DONOR_COLOR : ACCEPTOR_COLOR;
            for (int i = 0; i < measurement.binCentre.length; i++) {
                if (!usable[i]) {
                    continue;
                }
                int x = xFor(measurement.binCentre[i], maxRadius, width);
                int y = yFor(measurement.binProfile[i], height);
                g2.setColor(measured);
                g2.fillOval(x - 3, y - 3, 6, 6);
            }

            g2.setColor(Color.DARK_GRAY);
            g2.setFont(g2.getFont().deriveFont(Font.BOLD, 12.0f));
            g2.drawString(CHANNEL_NAMES[index] + " radial profile", 8, 14);

            if (fit != null) {
                g2.setFont(g2.getFont().deriveFont(Font.PLAIN, 11.0f));
                g2.setColor(FIT_COLOR);
                g2.drawString(String.format("sigma %.2f px, %.2f waves", fit.sigma, fit.waves),
                        MARGIN_LEFT + width - 150, MARGIN_TOP + 2);
            }

            drawLegend(g2, measured, fit, width, height);
            g2.dispose();
        }

        /**
         * What the two kinds of circle are, which is not guessable from the plot.
         *
         * Bottom left, because the profile runs from top left to bottom right and that corner is
         * the one it never reaches. Drawn last, over a panel-coloured box, so a curve that does
         * stray into it does not read through the text.
         */
        private void drawLegend(Graphics2D g2, Color measured, smFRETPSF.Fit fit,
                                int width, int height) {
            boolean showPedestal = (fit != null) && pedestalBox.isSelected()
                    && (fit.pedestal > floor);

            java.util.List<String> labels = new java.util.ArrayList<>();
            labels.add("measured");
            labels.add("before border correction");
            labels.add("fit");
            if (showPedestal) {
                labels.add("pedestal");
            }

            g2.setFont(g2.getFont().deriveFont(Font.PLAIN, 10.0f));
            FontMetrics metrics = g2.getFontMetrics();
            int textWidth = 0;
            for (String label : labels) {
                textWidth = Math.max(textWidth, metrics.stringWidth(label));
            }

            int rowHeight = metrics.getHeight();
            int boxWidth = textWidth + 30;
            int boxHeight = (rowHeight * labels.size()) + 8;
            int left = MARGIN_LEFT + 6;
            int top = MARGIN_TOP + height - boxHeight - 6;

            // Not fully opaque: the legend should never be the reason a data point is invisible.
            g2.setColor(new Color(255, 255, 255, 215));
            g2.fillRect(left, top, boxWidth, boxHeight);
            g2.setColor(new Color(205, 205, 205));
            g2.drawRect(left, top, boxWidth, boxHeight);

            int markerX = left + 12;
            int y = top + 4 + metrics.getAscent();
            for (String label : labels) {
                int markerY = y - (metrics.getAscent() / 2) + 1;
                if ("measured".equals(label)) {
                    g2.setColor(measured);
                    g2.fillOval(markerX - 3, markerY - 3, 6, 6);
                } else if ("before border correction".equals(label)) {
                    g2.setColor(RAW_COLOR);
                    g2.drawOval(markerX - 3, markerY - 3, 6, 6);
                } else if ("fit".equals(label)) {
                    g2.setColor(FIT_COLOR);
                    g2.setStroke(new BasicStroke(1.6f));
                    g2.drawLine(markerX - 6, markerY, markerX + 6, markerY);
                    g2.setStroke(new BasicStroke(1.0f));
                } else {
                    g2.setColor(new Color(FIT_COLOR.getRed(), FIT_COLOR.getGreen(),
                            FIT_COLOR.getBlue(), 120));
                    g2.drawLine(markerX - 6, markerY, markerX + 6, markerY);
                }

                g2.setColor(Color.DARK_GRAY);
                g2.drawString(label, markerX + 12, y);
                y += rowHeight;
            }
        }
    }

    /**
     * Build and show the window.
     */
    private void show() {
        frame = new JFrame("smFRET PSF - " + new File(root).getName());
        frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);

        frame.getContentPane().add(buildContent(), BorderLayout.CENTER);
        frame.pack();
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);

        update();
    }

    /**
     * Everything inside the window, as one panel, with no window around it.
     *
     * Separate from show() because a JPanel can be laid out and painted with no display, where a
     * JFrame cannot be constructed at all - so this is what lets the four panels be rendered and
     * looked at from a headless test. It is also the unit the PNG export captures.
     */
    JPanel buildContent() {
        JPanel panels = new JPanel(new GridLayout(2, 2, 4, 4));
        for (int c = 0; c < 2; c++) {
            imagePanel[c] = new PsfImagePanel(c);
            panels.add(imagePanel[c]);
        }
        for (int c = 0; c < 2; c++) {
            profilePanel[c] = new ProfilePanel(c);
            panels.add(profilePanel[c]);
        }

        patchSpinner = new JSpinner(new SpinnerNumberModel(smFRETPSF.DEFAULT_PATCH, 5, 20, 1));
        maskSpinner = new JSpinner(new SpinnerNumberModel(
                smFRETPSF.DEFAULT_NEIGHBOUR_MASK, 0.0, 15.0, 0.5));
        binsSpinner = new JSpinner(new SpinnerNumberModel(smFRETPSF.DEFAULT_BINS, 6, 40, 1));
        pedestalBox = new JCheckBox("Fit a pedestal", true);

        // None of these is guessable from its label, and two of them interact in a way that is
        // invisible from either one alone, so the tooltips say so.
        String patchTip = "<html>Half width of the square cut around each spot, in pixels.<br>"
                + "Wide enough to hold the wings the model is fitting, narrow enough that a"
                + " crowded<br>field still leaves border pixels around each spot - the border is"
                + " what the local<br>background is taken from.</html>";
        String maskTip = "<html>Pixels this close to <i>another</i> spot are discarded, in"
                + " pixels.<br>"
                + "Contaminated pixels are dropped rather than contaminated spots, which is what"
                + " keeps<br>a crowded field usable. <b>A neighbour further away than Patch +"
                + " Neighbour mask cannot<br>reach a patch pixel at all</b>, so raising this past"
                + " that does nothing until Patch goes up<br>too. Set it to 0 to switch the"
                + " masking off.</html>";
        String binsTip = "<html>How many radial bins the pooled pixels are averaged into.<br>"
                + "More bins resolve the profile more finely and put fewer pixels in each; bins"
                + " holding<br>too few pixels are dropped from the plot and from the fit.</html>";

        patchSpinner.setToolTipText(patchTip);
        maskSpinner.setToolTipText(maskTip);
        binsSpinner.setToolTipText(binsTip);
        pedestalBox.setToolTipText("<html>Fit a flat term under the PSF.<br>"
                + "On the example data this cuts the residual more than fourfold and leaves about"
                + " 3% of<br>peak, flat out where any real wing would still be decaying - which"
                + " looks like a scattered<br>light halo rather than aberration. Simulated data,"
                + " which has no halo, fits a pedestal<br>of almost zero.</html>");

        patchSpinner.addChangeListener(e -> update());
        maskSpinner.addChangeListener(e -> update());
        binsSpinner.addChangeListener(e -> update());

        // Only the drawing and the status line change, so this does not re-measure.
        pedestalBox.addActionListener(e -> {
            for (int c = 0; c < 2; c++) {
                profilePanel[c].repaint();
            }
            updateAdvice();
            statusLabel.setText(describe());
        });

        JPanel controls = new JPanel(new GridBagLayout());
        controls.setBorder(new EmptyBorder(6, 8, 4, 8));
        GridBagConstraints at = new GridBagConstraints();
        at.insets = new Insets(2, 4, 2, 4);
        at.anchor = GridBagConstraints.WEST;

        at.gridy = 0;
        at.gridx = 0;
        JLabel patchLabel = new JLabel("Patch");
        patchLabel.setToolTipText(patchTip);
        controls.add(patchLabel, at);
        at.gridx = 1;
        controls.add(patchSpinner, at);

        at.gridx = 2;
        JLabel maskLabel = new JLabel("Neighbour mask");
        maskLabel.setToolTipText(maskTip);
        controls.add(maskLabel, at);
        at.gridx = 3;
        controls.add(maskSpinner, at);

        at.gridx = 4;
        JLabel binsLabel = new JLabel("Bins");
        binsLabel.setToolTipText(binsTip);
        controls.add(binsLabel, at);
        at.gridx = 5;
        controls.add(binsSpinner, at);

        at.gridx = 6;
        controls.add(pedestalBox, at);

        JLabel modelLabel = new JLabel(MODEL_DESCRIPTION);
        modelLabel.setBorder(new EmptyBorder(2, 8, 0, 8));
        modelLabel.setForeground(FIT_COLOR);

        adviceLabel = new JLabel(" ");
        adviceLabel.setBorder(new EmptyBorder(2, 8, 0, 8));
        adviceLabel.setToolTipText("<html>The Gaussian filter width that recovers the most trace"
                + " SNR from the PSF measured here.<br>"
                + "smFRET Time Traces measures a trace by convolving with a Gaussian of"
                + " SpotSigma and reading<br>the peak, so this is the matched filter width - and"
                + " because the PSF is not a Gaussian it is<br>not the fitted core."
                + " Aberration is the only thing that moves it: barely below 0.3 waves,<br>"
                + "half again as wide by 0.5.<br><br>"
                + "<b>Two things it does not know.</b> The noise is taken as background dominated,"
                + " which is the same<br>assumption the spot finder's SNR already makes."
                + " And SpotSigma also drives spot <i>finding</i>,<br>not just trace extraction -"
                + " re-running the spot finder with a new value changes which<br>molecules are"
                + " found, the masks and the prominence, so this is advice about traces"
                + " alone.</html>");

        statusLabel = new JLabel(" ");
        statusLabel.setBorder(new EmptyBorder(2, 8, 4, 8));

        JButton saveCsv = new JButton("Save CSV");
        saveCsv.addActionListener(e -> onSaveCsv());
        JButton savePng = new JButton("Save PNG");
        savePng.addActionListener(e -> onSavePng());

        JPanel buttons = new JPanel(new FlowLayout(FlowLayout.RIGHT, 6, 2));
        buttons.add(saveCsv);
        buttons.add(savePng);

        JPanel captions = new JPanel(new GridLayout(3, 1));
        captions.add(modelLabel);
        captions.add(adviceLabel);
        captions.add(statusLabel);

        JPanel bottom = new JPanel(new BorderLayout());
        bottom.add(controls, BorderLayout.NORTH);
        bottom.add(captions, BorderLayout.CENTER);
        bottom.add(buttons, BorderLayout.SOUTH);

        plotPanels = panels;

        JPanel content = new JPanel(new BorderLayout());
        content.add(panels, BorderLayout.CENTER);
        content.add(bottom, BorderLayout.SOUTH);
        return content;
    }

    /**
     * The measured profiles and both fits, as a table.
     */
    private void onSaveCsv() {
        JFileChooser chooser = new JFileChooser(new File(root).getParentFile());
        chooser.setSelectedFile(new File(root + "_psf.csv"));
        if (chooser.showSaveDialog(frame) != JFileChooser.APPROVE_OPTION) {
            return;
        }

        try (PrintWriter writer = new PrintWriter(chooser.getSelectedFile())) {
            writer.println("# PSF measured from " + spotJSONFile);
            writer.println("# " + MODEL_DESCRIPTION);
            writer.println("# " + adviceLabel.getText());
            writer.println("# frames " + firstFrame + "-"
                    + Math.min(lastFrame, smFRETChannelMapper.frameCount(movie))
                    + " (as used for spot finding), patch " + patchSpinner.getValue()
                    + ", neighbour mask " + maskSpinner.getValue()
                    + ", bins " + binsSpinner.getValue());
            writer.println("# spot sigma used for finding: " + spotSigma);

            for (int c = 0; c < 2; c++) {
                smFRETPSF.Measurement m = channel[c];
                writer.println("# " + CHANNEL_NAMES[c] + ": " + m.samples.spotsUsed + " of "
                        + m.samples.spotsTotal + " spots, border correction "
                        + m.borderLevel + (m.borderSettled ? "" : " (did not settle)"));
                writeFit(writer, "#   PSF only      ", m.psfOnly);
                writeFit(writer, "#   PSF + pedestal", m.withPedestal);
            }

            writer.println("channel,radius,measured,uncorrected,pixels");
            for (int c = 0; c < 2; c++) {
                smFRETPSF.Measurement m = channel[c];
                boolean[] usable = m.usable();
                for (int i = 0; i < m.binCentre.length; i++) {
                    if (!usable[i]) {
                        continue;
                    }
                    writer.println(CHANNEL_NAMES[c] + "," + m.binCentre[i] + ","
                            + m.binProfile[i] + "," + m.binRaw[i] + "," + m.binCount[i]);
                }
            }
        } catch (Exception e) {
            log.info(e);
            JOptionPane.showMessageDialog(frame, "Could not write the table:\n" + e.getMessage(),
                    "smFRET PSF", JOptionPane.ERROR_MESSAGE);
        }
    }

    private static void writeFit(PrintWriter writer, String label, smFRETPSF.Fit fit) {
        if (fit == null) {
            writer.println(label + " did not fit");
            return;
        }
        writer.println(String.format("%s sigma %.4f px, %.4f waves, amplitude %.4f,"
                        + " pedestal %.4f, rms %.4f, first zero %.4f px",
                label, fit.sigma, fit.waves, fit.amplitude, fit.pedestal, fit.rms,
                fit.firstZero()));
    }

    /**
     * The four panels as one image, titled so a saved plot can be identified on its own.
     */
    private void onSavePng() {
        JFileChooser chooser = new JFileChooser(new File(root).getParentFile());
        chooser.setSelectedFile(new File(root + "_psf.png"));
        if (chooser.showSaveDialog(frame) != JFileChooser.APPROVE_OPTION) {
            return;
        }

        try {
            Container panels = plotPanels;
            int titleHeight = 28;
            BufferedImage image = new BufferedImage(panels.getWidth(),
                    panels.getHeight() + titleHeight + 50, BufferedImage.TYPE_INT_RGB);

            Graphics2D g2 = image.createGraphics();
            g2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING,
                    RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
            g2.setColor(Color.WHITE);
            g2.fillRect(0, 0, image.getWidth(), image.getHeight());

            g2.setColor(Color.DARK_GRAY);
            g2.setFont(g2.getFont().deriveFont(Font.BOLD, 14.0f));
            String title = new File(root).getName();
            FontMetrics metrics = g2.getFontMetrics();
            g2.drawString(title, (image.getWidth() - metrics.stringWidth(title)) / 2, 19);

            g2.translate(0, titleHeight);
            panels.paint(g2);

            g2.translate(0, panels.getHeight());
            g2.setFont(g2.getFont().deriveFont(Font.PLAIN, 11.0f));
            g2.setColor(FIT_COLOR);
            g2.drawString(MODEL_DESCRIPTION, 8, 13);
            g2.drawString(adviceLabel.getText(), 8, 27);
            g2.setColor(Color.DARK_GRAY);
            g2.drawString(describe(), 8, 41);
            g2.dispose();

            ImageIO.write(image, "png", chooser.getSelectedFile());
        } catch (Exception e) {
            log.info(e);
            JOptionPane.showMessageDialog(frame, "Could not write the image:\n" + e.getMessage(),
                    "smFRET PSF", JOptionPane.ERROR_MESSAGE);
        }
    }

    /** What the last measurement of a channel came to. Null before anything has been measured. */
    smFRETPSF.Measurement measurement(int which) {
        return channel[which];
    }

    /** The spot table, in the reloaded layout: x at column 0. */
    double[][] spots() {
        return spots;
    }

    @Override
    public void run() {
        try {
            load(spotJSONFile);
            if (isHeadless) {

                // Nothing to show, and nothing this plugin does is worth doing without somewhere
                // to show it - it writes only on request. Measure once so a headless run still
                // reports the fit to the log, which is what a macro would be after.
                buildHalves();
                for (int c = 0; c < 2; c++) {
                    smFRETPSF.Measurement m = smFRETPSF.analyse(
                            smFRETPSF.extract(cachedHalf[c], cachedWidth, cachedHeight, spots,
                                    smFRETPSF.DEFAULT_PATCH, smFRETPSF.DEFAULT_NEIGHBOUR_MASK),
                            smFRETPSF.DEFAULT_BINS);
                    channel[c] = m;
                    if (m.withPedestal == null) {
                        log.info(CHANNEL_NAMES[c] + ": " + m.samples.spotsUsed
                                + " spots, too few usable bins to fit");
                        continue;
                    }
                    log.info(CHANNEL_NAMES[c] + ": " + m.samples.spotsUsed + " spots, sigma "
                            + m.withPedestal.sigma + " px, " + m.withPedestal.waves + " waves");
                }

                if ((channel[DONOR].withPedestal != null)
                        && (channel[ACCEPTOR].withPedestal != null)) {
                    smFRETPSF.FilterOptimum joint = smFRETPSF.optimalFilter(new double[][] {
                        smFRETPSF.snrProfile(channel[DONOR].withPedestal),
                        smFRETPSF.snrProfile(channel[ACCEPTOR].withPedestal)}, FILTER_TOLERANCE);
                    log.info(String.format("best spotSigma for trace SNR: %.2f"
                            + " (within %.0f%% over %.2f to %.2f)",
                            joint.best, 100.0 * FILTER_TOLERANCE, joint.low, joint.high));
                }
                return;
            }
            SwingUtilities.invokeLater(this::show);
        } catch (smFRETAnalysisException e) {
            smFRETFiles.report(log, e);
        } catch (Exception e) {
            log.info(e);
            log.error(e);
        }
    }
}
