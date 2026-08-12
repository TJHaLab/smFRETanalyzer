import org.scijava.Context;
import org.scijava.log.LogService;
import java.io.File;

/**
 * Stage 2 on a simulated movie with the tolerance and contamination forced, for the sweeps.
 *
 * Pass a contamination of 1 to disable that filter, which is what the sweeps do: it is the last
 * stage and affects nothing upstream, so one run with it off gives the score of every surviving
 * spot and any threshold can then be applied offline.
 *
 * **1 is off, not 0, and the two are opposite ends.** This argument used to be spotProminence,
 * which rejected *below* its value and so was disabled with a large negative number;
 * spotContamination rejects *above* its value, and the prediction is clamped to [0, 1], so 1
 * keeps everything and 0 keeps almost nothing. A sweep that carried the old sentinel over would
 * reject every spot and score an empty field.
 *
 * args: &lt;dir&gt; &lt;sigma&gt; &lt;tolerance&gt; &lt;contamination&gt; &lt;frames&gt; [threshold]
 */
public class RunTune {
    public static void main(String[] args) throws Exception {
        String dir = args[0];

        Context context = new Context(LogService.class);
        LogService log = context.getService(LogService.class);

        smFRETSpotFinder s = new smFRETSpotFinder();
        s.log = log;
        s.inputImageName = new File(dir, "sim.tif");
        s.mappingFile = new File(dir, "sim_mapping.json");
        s.startSlice = 1;
        s.endSlice = Integer.parseInt(args[4]);
        s.spotSigma = Double.parseDouble(args[1]);
        s.spotTolerance = Double.parseDouble(args[2]);
        s.spotContamination = Double.parseDouble(args[3]);
        s.spotThreshold = (args.length > 5) ? Double.parseDouble(args[5]) : 6.0;
        s.cameraBlackLevel = 5; s.cameraGain = 1.0;
        s.spotSpacing = 3; s.edgeMargin = 5; s.backgroundKappa = 0.0;
        s.run();
        System.exit(0);
    }
}
