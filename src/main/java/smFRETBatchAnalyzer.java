/*
 * Run spot finding and trace measurement over a folder of movies, unattended.
 *
 * The two stages this drives are unchanged: it constructs an smFRETSpotFinder and an
 * smFRETAnalyzer per movie and runs them exactly as the menu items would. What it adds is the
 * settings, which come from an example spot finder JSON rather than from a dialog - so a batch is
 * "do to these forty movies what I just did to this one", and the run that is being copied has
 * been looked at on screen first.
 */

import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Menu;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.io.File;
import java.util.ArrayList;
import java.util.Map;


// The single-type import of org.scijava.plugin.Menu shadows java.awt.Menu from the wildcard
// import above, so Menu here is the SciJava annotation.
@Plugin(type = Command.class,
        menu = {@Menu(label = "Plugins"),
                @Menu(label = "smFRET"),
                @Menu(label = "smFRET Batch Analysis", weight = 7.0)})
public class smFRETBatchAnalyzer implements Command {

    @Parameter
    LogService log;

    static final String WINDOW_TITLE = "smFRET Batch Analysis";

    // For anything that went wrong, in the settings line and on a movie's row.
    private static final Color PROBLEM_COLOR = new Color(150, 60, 40);

    // What a movie is doing, which is both the row's text and what the counts at the end are of.
    static final int QUEUED = 0;
    static final int RUNNING = 1;
    static final int DONE = 2;
    static final int FAILED = 3;
    static final int STOPPED = 4;

    // Member variables.
    private JButton addButton;
    private JButton analyzeButton;
    private JSpinner backgroundSpinner;
    private JButton chooseButton;
    final java.util.List<Job> jobs = new ArrayList<>();
    private JPanel jobsPanel;
    private JProgressBar overallBar;
    private volatile boolean running = false;
    private JLabel settingsLabel;
    Settings settings;
    private JLabel statusLabel;
    private volatile boolean stopRequested = false;
    private File templateFile;

    /**
     * One movie and what has happened to it.
     *
     * The state is written from the worker thread and read from the event thread, hence volatile:
     * the row is repainted from an invokeLater that runs after the write, and without it there is
     * no guarantee the repaint sees the value that prompted it.
     */
    static final class Job {

        final File movie;
        volatile String detail = "";
        volatile int state = QUEUED;

        // 0 = nothing done, 1 = spot finding done, 2 = traces done. The progress bar's fill, so
        // it shows finished work rather than a guess at how far through the current stage it is.
        volatile int stagesDone = 0;

        Job(File movie) {
            this.movie = movie;
        }

        String name() {
            return movie.getName();
        }
    }

    /**
     * The spot finder settings a batch runs with, read from an example run's JSON.
     *
     * Every field is optional in the file. JSONs written before 1.2 have no spot channel and a
     * "spot prominence" that no longer exists, and the ones written before that have no slice
     * range either - so an absent key leaves the shipped default rather than failing the batch,
     * which is the same tolerance smFRETAnalyzer already applies to "background kappa".
     */
    static final class Settings {

        final double backgroundKappa;
        final int cameraBlackLevel;
        final double cameraGain;
        final int edgeMargin;
        final int endSlice;
        final File mappingFile;
        final String spotChannel;
        final double spotContamination;
        final double spotSigma;
        final int spotSpacing;
        final double spotThreshold;
        final double spotTolerance;
        final int startSlice;

        // Whether the example run measured its tolerance rather than being told one. Kept because
        // it is the difference between copying a number and copying a decision - see read().
        final boolean toleranceEstimated;

        private Settings(Map<String, Object> json) {
            smFRETSpotFinder defaults = new smFRETSpotFinder();

            backgroundKappa = asDouble(json.get("background kappa"), defaults.backgroundKappa);
            cameraBlackLevel = asInt(json.get("camera black"), defaults.cameraBlackLevel);
            cameraGain = asDouble(json.get("camera gain"), defaults.cameraGain);
            edgeMargin = asInt(json.get("edge margin"), defaults.edgeMargin);
            endSlice = asInt(json.get("end slice"), defaults.endSlice);
            spotContamination = asDouble(json.get("spot contamination"), defaults.spotContamination);
            spotSigma = asDouble(json.get("spot sigma"), defaults.spotSigma);
            spotSpacing = asInt(json.get("spot spacing"), defaults.spotSpacing);
            spotThreshold = asDouble(json.get("spot threshold"), defaults.spotThreshold);
            startSlice = asInt(json.get("start slice"), defaults.startSlice);

            String channel = (String) json.get("spot channel");
            spotChannel = ("sum".equals(channel) || "donor".equals(channel)
                    || "acceptor".equals(channel)) ? channel : defaults.spotChannel;

            // The mapping the example run used, which is why this plugin does not ask for one.
            // It is recorded as an absolute path, so it is right whenever the example JSON and
            // its mapping have not moved since - and when they have, saying so is more use than
            // a second file chooser that would usually be answered with the same file.
            Object mapping = json.get("mapping file");
            mappingFile = (mapping == null) ? null : new File(mapping.toString());

            // A run that measured its tolerance is copied as *the decision to measure*, not as
            // the number it arrived at. The JSON records the effective value - the estimate - and
            // reusing that would pin every movie in the batch to one movie's noise floor, which
            // is the thing estimating it was for. A tolerance that was actually typed in is a
            // choice, and is copied as it stands.
            toleranceEstimated = Boolean.TRUE.equals(json.get("spot tolerance estimated"));
            spotTolerance = toleranceEstimated ? 0.0
                    : asDouble(json.get("spot tolerance"), defaults.spotTolerance);
        }

        /**
         * Read an example spot finder JSON.
         *
         * Throws if it is not one - readSpotFinderJSON checks for a key only the spot finder
         * writes, so a mapping JSON or a movie chosen here is named as such rather than arriving
         * as a missing key halfway through the first movie.
         */
        static Settings read(File jsonFile) {
            return new Settings(smFRETFiles.readSpotFinderJSON(jsonFile));
        }

        /** Copy these onto a spot finder, leaving it the movie and the mapping to be told. */
        void applyTo(smFRETSpotFinder finder) {
            finder.backgroundKappa = backgroundKappa;
            finder.cameraBlackLevel = cameraBlackLevel;
            finder.cameraGain = cameraGain;
            finder.edgeMargin = edgeMargin;
            finder.endSlice = endSlice;
            finder.spotChannel = spotChannel;
            finder.spotContamination = spotContamination;
            finder.spotSigma = spotSigma;
            finder.spotSpacing = spotSpacing;
            finder.spotThreshold = spotThreshold;
            finder.spotTolerance = spotTolerance;
            finder.startSlice = startSlice;
        }

        /** The settings as one line, so what a batch is about to do is visible before it runs. */
        String describe() {
            return String.format("slices %d-%d · sigma %.2f · threshold %.1f · tolerance %s"
                            + " · contamination %.2f · %s · gain %.2f, black %d",
                    startSlice, endSlice, spotSigma, spotThreshold,
                    toleranceEstimated ? "estimated" : String.format("%.1f", spotTolerance),
                    spotContamination, spotChannel, cameraGain, cameraBlackLevel);
        }

        /**
         * Why this template cannot be used, or null when it can.
         *
         * Only the mapping is checked, because it is the only thing here that is a file - every
         * other setting is a number that is as valid for one movie as another.
         */
        String problem() {
            if (mappingFile == null) {
                return "the settings file records no mapping file - it was written by a version"
                        + " that did not save one, so re-run smFRET Spot Finder on one movie to"
                        + " make a current settings file.";
            }
            if (!mappingFile.isFile()) {
                return "the mapping file it names is not there: " + mappingFile
                        + " - move it back, or re-run smFRET Spot Finder on one movie with the"
                        + " mapping you want and use that run's settings file.";
            }
            return null;
        }

        private static double asDouble(Object value, double fallback) {
            return (value instanceof Number) ? ((Number) value).doubleValue() : fallback;
        }

        private static int asInt(Object value, int fallback) {
            return (value instanceof Number) ? ((Number) value).intValue() : fallback;
        }
    }

    /**
     * Told which stage a movie has reached, so the analysis can report progress without knowing
     * there is a window. Called from the worker thread.
     */
    interface Progress {
        void stage(int stagesDone, String what);
    }

    /**
     * Run both stages over one movie. Returns null on success, or why it failed.
     *
     * Both stages catch their own exceptions and report them to the log rather than throwing -
     * which is right for a menu item, where the user is watching, and leaves a batch with nothing
     * to test. So success is judged by the output: each stage's file has to have *changed*, which
     * is what distinguishes a run that worked from a stale file left by an earlier one.
     *
     * Package visible and taking its settings rather than reading the window, so a batch can be
     * exercised end to end without a display.
     */
    String analyze(File movie, Settings settings, int backgroundFrames, Progress progress) {
        String root = movie.toString();
        int dot = root.lastIndexOf('.');
        if (dot > 0) {
            root = root.substring(0, dot);
        }

        File json = new File(root + "_spotf_finding.json");
        File h5 = new File(root + ".h5");
        String jsonBefore = stamp(json);
        String h5Before = stamp(h5);

        progress.stage(0, "finding spots");
        smFRETSpotFinder finder = new smFRETSpotFinder();
        finder.log = log;

        // No UIService and no PrefService injected, and neither is an oversight. showQCImage
        // keeps findSpots away from the one, and a null PrefService makes its saveSettings a no
        // op - a batch should not overwrite what the interactive dialog will open with next time.
        finder.showQCImage = false;
        finder.inputImageName = movie;
        finder.mappingFile = settings.mappingFile;
        settings.applyTo(finder);
        finder.findSpots();

        if (!changed(jsonBefore, stamp(json))) {
            return "spot finding failed - no " + json.getName() + " was written, see the log";
        }

        progress.stage(1, "measuring traces");
        smFRETAnalyzer traces = new smFRETAnalyzer();
        traces.log = log;
        traces.spotJSONFile = json;
        traces.backgroundAverageNFrames = backgroundFrames;
        traces.run();

        if (!changed(h5Before, stamp(h5))) {
            return "trace measurement failed - no " + h5.getName() + " was written, see the log";
        }

        progress.stage(2, "done");
        return null;
    }

    /**
     * A file's identity for "did this run rewrite it", or null when it is not there.
     *
     * Compared rather than tested against the time the run started, which would be the obvious
     * thing: that assumes the clock the file system stamps with is the clock this JVM reads, and
     * on a network share it is not. Comparing the file against itself needs no such assumption.
     */
    private static String stamp(File file) {
        return file.exists() ? (file.lastModified() + ":" + file.length()) : null;
    }

    private static boolean changed(String before, String after) {
        return (after != null) && !after.equals(before);
    }

    /**
     * Work through the queued movies. Runs on its own thread.
     *
     * One at a time. Two movies at once would halve nothing that matters - both stages are
     * already bounded by reading the movie and by memory, and ImageJ's statics are not something
     * to share between two analyses on a promise.
     */
    private void runBatch() {
        int done = 0;
        int failed = 0;
        int total = jobs.size();

        for (Job job : jobs) {
            if (stopRequested) {
                job.state = STOPPED;
                job.detail = "stopped before this movie";
                repaintJob(job);
                continue;
            }

            job.state = RUNNING;
            job.stagesDone = 0;
            job.detail = "starting";
            repaintJob(job);
            log.info("batch: " + job.name());

            String problem;
            try {
                problem = analyze(job.movie, settings, backgroundFrames(), (stages, what) -> {
                    job.stagesDone = stages;
                    job.detail = what;
                    repaintJob(job);
                });
            } catch (Exception e) {

                // Neither stage is supposed to throw, so anything here is unexpected. It stops
                // this movie and not the batch: the usual cause is something about this one
                // movie, and thirty-nine others should not be lost to it.
                log.info(e);
                problem = String.valueOf(e.getMessage());
            }

            if (problem == null) {
                job.state = DONE;
                job.stagesDone = 2;
                job.detail = "done";
                done += 1;
            } else {
                job.state = FAILED;
                job.detail = problem;
                failed += 1;
            }
            repaintJob(job);

            final int finished = done + failed;
            SwingUtilities.invokeLater(() -> {
                overallBar.setValue(finished);
                overallBar.setString(finished + " of " + total);
            });
        }

        final int okCount = done;
        final int failCount = failed;
        SwingUtilities.invokeLater(() -> {
            running = false;
            stopRequested = false;
            setStatus(String.format("%d of %d analysed%s", okCount, total,
                    (failCount > 0) ? (", " + failCount + " failed - see the log") : ""));
            enableControls(true);
        });
    }

    /**
     * One movie's row, kept so that progress can be written into it rather than rebuilt.
     *
     * Rebuilding the list on every stage change would be simpler and is what the trace histogram's
     * pool does, but that list changes shape when it changes at all; this one changes twice per
     * movie while the user is watching, and rebuilding it would throw the scroll position away
     * each time - on a batch of forty movies, exactly when it matters.
     */
    private static final class JobRow {

        final JProgressBar bar;
        final JLabel detail;
        final JPanel panel;
        final JButton remove;

        JobRow(JPanel panel, JProgressBar bar, JLabel detail, JButton remove) {
            this.bar = bar;
            this.detail = detail;
            this.panel = panel;
            this.remove = remove;
        }
    }

    private final java.util.Map<Job, JobRow> rowFor = new java.util.LinkedHashMap<>();

    /**
     * Run ...
     */
    @Override
    public void run() {
        if (GraphicsEnvironment.isHeadless()) {
            log.info(WINDOW_TITLE + " is interactive and cannot run headless");
            return;
        }
        SwingUtilities.invokeLater(this::showWindow);
    }

    /**
     * Build the window.
     */
    private void showWindow() {
        JFrame frame = new JFrame(WINDOW_TITLE);
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

        // Closing the window during a batch stops it after the current movie. The worker is a
        // daemon thread and would otherwise carry on writing files for however long is left,
        // with nothing on screen to say so - and the obvious way to stop a run you have changed
        // your mind about is to close the window.
        frame.addWindowListener(new java.awt.event.WindowAdapter() {
            @Override
            public void windowClosing(java.awt.event.WindowEvent e) {
                if (running) {
                    stopRequested = true;
                    log.info("batch: window closed, stopping after the current movie");
                }
            }
        });

        // Settings row. The mapping is not asked for: the settings file records the one its own
        // run used, which is the mapping that goes with those settings.
        chooseButton = new JButton("Settings...");
        chooseButton.setToolTipText("Choose the '_spotf_finding.json' written by an smFRET Spot"
                + " Finder run whose settings this batch should copy.");
        chooseButton.addActionListener(e -> onChooseTemplate());

        settingsLabel = new JLabel("no settings chosen");
        settingsLabel.setForeground(Color.GRAY);

        JPanel settingsPanel = new JPanel(new BorderLayout(8, 0));
        settingsPanel.setBorder(new EmptyBorder(8, 10, 4, 10));
        settingsPanel.add(chooseButton, BorderLayout.WEST);
        settingsPanel.add(settingsLabel, BorderLayout.CENTER);

        // Movies.
        jobsPanel = new JPanel();
        jobsPanel.setLayout(new BoxLayout(jobsPanel, BoxLayout.Y_AXIS));
        jobsPanel.setBackground(Color.WHITE);

        JScrollPane scroll = new JScrollPane(jobsPanel,
                JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scroll.setPreferredSize(new Dimension(720, 260));
        scroll.getVerticalScrollBar().setUnitIncrement(14);

        // Once, on the scroll pane, which covers the rows inside it too - see smFRETSwing.
        // Refused while a batch is running, since there would be nothing to add movies to.
        scroll.setTransferHandler(smFRETSwing.fileDropHandler(() -> !running, this::onFilesAdded));

        addButton = new JButton("Add movies...");
        addButton.addActionListener(e -> onAddMovies());

        JLabel hint = new JLabel("<html><i>drop movies here</i></html>");
        hint.setForeground(Color.GRAY);

        JPanel listHeader = new JPanel(new BorderLayout());
        listHeader.setBorder(new EmptyBorder(0, 10, 2, 10));
        listHeader.add(new JLabel("Movies"), BorderLayout.WEST);
        listHeader.add(hint, BorderLayout.EAST);

        JPanel listPanel = new JPanel(new BorderLayout(0, 2));
        listPanel.setBorder(new EmptyBorder(0, 10, 4, 10));
        listPanel.add(listHeader, BorderLayout.NORTH);
        listPanel.add(scroll, BorderLayout.CENTER);

        // Trace measurement's one setting. There is no traces JSON to copy it from - stage 3 has
        // this and the spot finder JSON and nothing else - so it is typed here.
        backgroundSpinner = new JSpinner(new SpinnerNumberModel(30, 1, 100000, 1));
        ((JSpinner.NumberEditor) backgroundSpinner.getEditor()).getTextField().setColumns(5);
        backgroundSpinner.setToolTipText("Frames each background estimate is averaged over,"
                + " the one setting smFRET Time Traces takes.");

        overallBar = new JProgressBar(0, 1);
        overallBar.setStringPainted(true);

        analyzeButton = new JButton("Analyze");
        analyzeButton.addActionListener(e -> onAnalyze());

        JPanel runPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 8, 2));
        runPanel.add(new JLabel("Background frames"));
        runPanel.add(backgroundSpinner);
        runPanel.add(Box.createHorizontalStrut(12));
        runPanel.add(addButton);

        statusLabel = new JLabel(" ");

        // The status line gets the full width with the button under it, the same shape the trace
        // histogram's ended up: these messages name a file and say what went wrong with it, and
        // sharing a row with a button truncates exactly the useful end of them.
        JPanel bottomPanel = new JPanel(new BorderLayout(0, 4));
        bottomPanel.setBorder(new EmptyBorder(2, 10, 8, 10));
        bottomPanel.add(statusLabel, BorderLayout.NORTH);

        JPanel buttonRow = new JPanel(new BorderLayout(8, 0));
        buttonRow.add(overallBar, BorderLayout.CENTER);
        buttonRow.add(analyzeButton, BorderLayout.EAST);
        bottomPanel.add(buttonRow, BorderLayout.SOUTH);

        JPanel southPanel = new JPanel();
        southPanel.setLayout(new BoxLayout(southPanel, BoxLayout.Y_AXIS));
        runPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
        bottomPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
        southPanel.add(runPanel);
        southPanel.add(bottomPanel);

        frame.getContentPane().setLayout(new BorderLayout());
        frame.getContentPane().add(settingsPanel, BorderLayout.NORTH);
        frame.getContentPane().add(listPanel, BorderLayout.CENTER);
        frame.getContentPane().add(southPanel, BorderLayout.SOUTH);

        refreshJobsPanel();
        updateReadiness();

        frame.pack();
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }

    /**
     * Rebuild the movie rows.
     */
    private void refreshJobsPanel() {
        jobsPanel.removeAll();
        rowFor.clear();
        for (Job job : jobs) {
            JobRow row = buildJobRow(job);
            rowFor.put(job, row);
            jobsPanel.add(row.panel);
        }

        // Holds the rows at the top. Without it BoxLayout spreads three movies down the whole
        // height of the list.
        jobsPanel.add(Box.createVerticalGlue());
        jobsPanel.revalidate();
        jobsPanel.repaint();
    }

    /**
     * One movie's row: name, a progress bar, what it is doing, and a remove button.
     */
    private JobRow buildJobRow(Job job) {
        JLabel name = new JLabel(job.name());
        name.setPreferredSize(new Dimension(220, name.getPreferredSize().height));
        name.setToolTipText(job.movie.getAbsolutePath());

        // Two steps, and the fill is *finished* stages rather than a guess at how far through the
        // current one it is. Neither stage reports its own progress, so a bar that crept would be
        // inventing the number; the text beside it says what is happening now.
        JProgressBar bar = new JProgressBar(0, 2);
        bar.setPreferredSize(new Dimension(110, 14));

        JLabel detail = new JLabel();
        JButton remove = new JButton("×");
        remove.setMargin(new Insets(0, 4, 0, 4));
        remove.setFocusable(false);
        remove.addActionListener(e -> {
            jobs.remove(job);
            refreshJobsPanel();
            resetProgress();
        });

        JPanel left = new JPanel(new FlowLayout(FlowLayout.LEFT, 6, 1));
        left.setBackground(Color.WHITE);
        left.add(name);
        left.add(bar);

        JPanel panel = new JPanel(new BorderLayout(6, 0));
        panel.setBackground(Color.WHITE);
        panel.add(left, BorderLayout.WEST);
        panel.add(detail, BorderLayout.CENTER);
        panel.add(remove, BorderLayout.EAST);
        panel.setMaximumSize(new Dimension(Integer.MAX_VALUE, panel.getPreferredSize().height));
        panel.setAlignmentX(Component.LEFT_ALIGNMENT);

        JobRow row = new JobRow(panel, bar, detail, remove);
        paintJob(job, row);
        return row;
    }

    /** Write a job's state into its row. Event thread only. */
    private void paintJob(Job job, JobRow row) {
        row.bar.setValue(job.stagesDone);
        row.remove.setEnabled(!running);

        String text = job.detail;
        if (job.state == QUEUED) {
            text = "queued";
        }
        row.detail.setText(text);

        // The failure messages are a sentence long and name a file, so they will outrun the row.
        // The tooltip is the whole of it, and it is in the log as well.
        row.detail.setToolTipText(text);
        row.detail.setForeground((job.state == FAILED) ? PROBLEM_COLOR
                : ((job.state == DONE) ? new Color(40, 110, 60) : Color.DARK_GRAY));
    }

    /** Ask for a job's row to be repainted. Called from the worker thread. */
    private void repaintJob(Job job) {
        SwingUtilities.invokeLater(() -> {
            JobRow row = rowFor.get(job);
            if (row != null) {
                paintJob(job, row);
            }
        });
    }

    /**
     * Choose the example run whose settings the batch copies.
     */
    private void onChooseTemplate() {
        JFileChooser chooser = new JFileChooser(defaultDirectory());
        chooser.setDialogTitle("Select a spot finder settings file (_spotf_finding.json)");
        if (chooser.showOpenDialog(ownerWindow()) != JFileChooser.APPROVE_OPTION) {
            return;
        }
        File chosen = chooser.getSelectedFile();
        try {
            settings = Settings.read(chosen);
        } catch (Exception e) {
            settings = null;
            templateFile = null;
            settingsLabel.setForeground(PROBLEM_COLOR);
            settingsLabel.setText(wrapped("<b>" + escaped(chosen.getName()) + "</b><br>"
                    + breakable(escaped(String.valueOf(e.getMessage())))));
            settingsLabel.setToolTipText(e.getMessage());
            setStatus(e.getMessage());
            updateReadiness();
            return;
        }

        templateFile = chosen;
        String problem = settings.problem();
        if (problem != null) {
            settingsLabel.setForeground(PROBLEM_COLOR);

            // The name on its own line, then the sentence. The path is the middle of that
            // sentence and is the part being reported, so it is also the part that must not be
            // what a narrow window drops.
            settingsLabel.setText(wrapped("<b>" + escaped(chosen.getName()) + "</b><br>"
                    + breakable(escaped(problem))));
            settingsLabel.setToolTipText(wrapped(escaped(problem)));
        } else {
            settingsLabel.setForeground(Color.DARK_GRAY);
            settingsLabel.setText(wrapped("<b>" + escaped(chosen.getName()) + "</b> · "
                    + escaped(settings.describe()) + "<br>mapping "
                    + escaped(settings.mappingFile.getName())));
            settingsLabel.setToolTipText(wrapped(escaped(chosen.getAbsolutePath()) + "<br>"
                    + escaped(settings.describe())));
        }
        setStatus(" ");
        updateReadiness();
    }

    /**
     * Prompt for movies to add.
     */
    private void onAddMovies() {
        JFileChooser chooser = new JFileChooser(defaultDirectory());
        chooser.setDialogTitle("Select the movies to analyse");
        chooser.setMultiSelectionEnabled(true);
        if (chooser.showOpenDialog(ownerWindow()) != JFileChooser.APPROVE_OPTION) {
            return;
        }
        onFilesAdded(java.util.Arrays.asList(chooser.getSelectedFiles()));
    }

    /**
     * Queue movies, skipping ones already queued and anything that is not an image.
     *
     * The check is a sniff of the file's first bytes rather than its extension, so dropping a
     * folder's worth of files queues the movies in it and says what it left out - which is the
     * ordinary way to fill this list, and would otherwise queue the JSONs and the CSVs alongside
     * them and fail on each in turn several minutes apart.
     */
    private void onFilesAdded(java.util.List<File> files) {
        java.util.List<String> skipped = new ArrayList<>();
        for (File file : files) {
            if (queued(file)) {
                continue;
            }
            if (smFRETFiles.isMovie(file)) {
                jobs.add(new Job(file));
            } else {
                skipped.add(file.getName());
            }
        }

        // Name order, so that a batch runs in the order the movies are named rather than in
        // whatever order the file manager handed them over.
        jobs.sort((a, b) -> {
            int byName = a.name().compareToIgnoreCase(b.name());
            return (byName != 0) ? byName : pathKey(a.movie).compareTo(pathKey(b.movie));
        });

        refreshJobsPanel();
        resetProgress();
        if (!skipped.isEmpty()) {
            setStatus("Not movies, skipped: " + String.join(", ", skipped));
        }
    }

    /** Start the batch. */
    private void onAnalyze() {
        if (running) {
            stopRequested = true;
            analyzeButton.setEnabled(false);
            setStatus("stopping after the current movie...");
            return;
        }

        String problem = settings.problem();
        if (problem != null) {
            JOptionPane.showMessageDialog(ownerWindow(), problem, WINDOW_TITLE,
                    JOptionPane.ERROR_MESSAGE);
            return;
        }

        // A fresh start for every movie, including the ones a previous run finished, so that
        // pressing Analyze twice means the same thing both times.
        for (Job job : jobs) {
            job.state = QUEUED;
            job.stagesDone = 0;
            job.detail = "";
        }
        refreshJobsPanel();

        running = true;
        stopRequested = false;
        overallBar.setMaximum(jobs.size());
        overallBar.setValue(0);
        overallBar.setString("0 of " + jobs.size());
        setStatus("analysing " + jobs.size() + " movies...");
        enableControls(false);

        // Off the event thread, or the window would not repaint for the hours this can take.
        Thread worker = new Thread(this::runBatch, "smFRET batch");
        worker.setDaemon(true);
        worker.start();
    }

    /** Whether this movie is already queued, compared by canonical path. */
    private boolean queued(File file) {
        String wanted = pathKey(file);
        for (Job job : jobs) {
            if (pathKey(job.movie).equals(wanted)) {
                return true;
            }
        }
        return false;
    }

    /**
     * The identity of a file, for telling "already queued" from "another movie of that name".
     *
     * Canonical rather than absolute: an absolute path keeps whatever spelling it arrived with, so
     * the same movie dropped from two places would be queued twice and analysed twice, the second
     * run overwriting the first's output while the first was still being looked at.
     */
    private static String pathKey(File file) {
        try {
            return file.getCanonicalPath();
        } catch (java.io.IOException e) {
            return file.getAbsolutePath();
        }
    }

    /**
     * Lay text out as HTML so that a JLabel wraps it instead of clipping it.
     *
     * A JLabel given plain text draws one line and truncates it to an ellipsis. That is survivable
     * for a status line and not for this one, which reports a path: the path sits in the middle of
     * the sentence, so the clipped half is the half naming the file the user has to go and find.
     * HTML text is laid out to the width the label was given and reflows when the window resizes.
     */
    static String wrapped(String body) {
        return "<html>" + body + "</html>";
    }

    /**
     * Escape text being put into one of those labels.
     *
     * File names and paths are user data, and one containing an ampersand or an angle bracket
     * would otherwise be swallowed by the HTML parser - so the message reporting a file would
     * mangle exactly the names most likely to have been chosen badly in the first place.
     */
    static String escaped(String text) {
        return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;");
    }

    /**
     * Give a long path somewhere to break.
     *
     * A path contains no spaces, so the line breaker has nowhere legal to break it: the label
     * clips it instead, and its minimum width becomes the width of the whole path - so a narrow
     * window gets a message running off its own right edge rather than a wrapped one. A zero width
     * space after each separator is a break opportunity that adds no visible character, so the
     * path still reads as a path and can still be copied out of the tooltip unchanged.
     *
     * Only separators, so ordinary prose is untouched.
     */
    static String breakable(String text) {
        return text.replace("/", "/\u200b").replace("\\", "\\\u200b");
    }

    private int backgroundFrames() {
        return ((Number) backgroundSpinner.getValue()).intValue();
    }

    /** Where the file choosers open: beside the movies if there are any, else the settings. */
    private File defaultDirectory() {
        if (!jobs.isEmpty()) {
            return jobs.get(0).movie.getParentFile();
        }
        return (templateFile == null) ? null : templateFile.getParentFile();
    }

    private void enableControls(boolean enabled) {
        addButton.setEnabled(enabled);
        backgroundSpinner.setEnabled(enabled);
        chooseButton.setEnabled(enabled);
        for (JobRow row : rowFor.values()) {
            row.remove.setEnabled(enabled);
        }
        analyzeButton.setEnabled(true);
        analyzeButton.setText(enabled ? "Analyze" : "Stop");
        if (enabled) {
            updateReadiness();
        }
    }

    /** The window the controls are in, or null before they are in one. */
    private java.awt.Window ownerWindow() {
        return (jobsPanel == null) ? null : SwingUtilities.getWindowAncestor(jobsPanel);
    }

    private void setStatus(String status) {
        statusLabel.setText(status);
        statusLabel.setToolTipText(status);
    }

    /**
     * Analyze is available only once there is something to run and settings to run it with.
     *
     * Disabled rather than pressable-and-complaining, because the two things missing are named on
     * the button's own tooltip and in the status line, and both are one click away.
     */
    /** The queue changed, so a previous batch's count no longer describes it. */
    private void resetProgress() {
        overallBar.setValue(0);
        updateReadiness();
    }

    private void updateReadiness() {
        boolean haveSettings = (settings != null) && (settings.problem() == null);
        boolean ready = haveSettings && !jobs.isEmpty();
        analyzeButton.setEnabled(ready);

        if (ready) {
            analyzeButton.setToolTipText("Analyse " + jobs.size() + " movies with these settings.");
        } else if (!haveSettings) {
            analyzeButton.setToolTipText("Choose a spot finder settings file first.");
        } else {
            analyzeButton.setToolTipText("Add the movies to analyse first.");
        }

        overallBar.setMaximum(Math.max(1, jobs.size()));

        // Left alone once a batch has run, so that its result is still on screen afterwards -
        // this is called again when the controls come back, and "4 movies" over a full bar reads
        // as if nothing had happened.
        if (overallBar.getValue() == 0) {
            overallBar.setString(jobs.size() + " movies");
        }
    }
}
