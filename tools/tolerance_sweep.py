"""Is a spotTolerance estimated from the image better than a fixed one? Behind issue #10.

MaximumFinder's tolerance is a flood fill depth in raw image units, so what counts as deep
enough to separate two peaks depends on where the noise floor of the averaged image sits.
That moves with the camera, the gain and above all the number of frames averaged - which is
why the fixed default of 5, swept at 20 frames, needed revisiting whenever any of those
changed. See the SpotTolerance section of SIMULATION.md.

This sweeps the two against each other across spot size, density and frame count, scoring the
way tune_sweep.py does: recall of truth spots that should be findable (true SNR >= the
threshold), plus false positives, with the match radius scaled by spot size because
MaximumFinder returns an integer coordinate.

The K values are swept as well as compared. The point is not only which wins but how sharply
the answer depends on the constant: a flat response is what makes the estimate safe to ship,
since the acceptor half is warped before the two are added and interpolation costs about 10%
of the measured noise on real data.

    python3 ../tools/tolerance_sweep.py            # writes tolerance_sweep.csv

Run it from a scratch directory. It writes 40 simulated movies.
"""
import csv
import os
import subprocess
import sys

import numpy as np
import tifffile

HERE = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.join(HERE, "gen_calib.py")
CLASSES = os.path.join(os.path.dirname(HERE), "target", "classes")

SIGMAS = [1.0, 1.5, 2.0, 2.5, 3.0]
DENSITIES = [400, 150]
FRAME_COUNTS = [5, 10, 20, 40]
KS = [2.0, 2.5, 3.0, 3.5, 4.0]
FIXED = 5.0                # what the parameter defaulted to while it was set by hand
MEDIAN_N = 400.0
THRESHOLD = 6.0
SEED = 11


def _tool(name, default):
    return os.environ.get(name, default)


JAVA = _tool("SMFRET_JAVA", "/usr/lib/jvm/java-11-openjdk-amd64/bin/java")
PYTHON = _tool("SMFRET_PYTHON", sys.executable)


def classpath():
    cp = os.environ.get("SMFRET_CP")
    if cp:
        return cp
    with open(os.path.join(HERE, "cp.txt")) as fh:
        return fh.read().strip()


def match_radius(sigma):
    """MaximumFinder returns an integer coordinate, so localisation error grows with the PSF."""
    return min(3.0, max(1.5, 1.2 * sigma))


def spots_csv(directory):
    """The spot table, in the analysis folder the plugin writes it to (issue #9)."""
    return os.path.join(directory, "sim_analysis", "sim_spotf_spots.csv")


def pixel_noise(image):
    """The same estimator the plugin uses: MAD of horizontal differences, undone for the pair."""
    d = np.diff(image, axis=1).ravel()
    return 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2.0)


def summed_average(path, frames):
    with tifffile.TiffFile(path) as tf:
        movie = np.stack([p.asarray() for p in tf.pages[:frames]]).astype(np.float64)
    average = movie.mean(axis=0)
    width = average.shape[1]
    return average[:, : width // 2] + average[:, width // 2 :]


def score(directory, truth, findable, match):
    table = np.genfromtxt(spots_csv(directory), delimiter=",", names=True, dtype=None)
    if table.size == 0:
        return 0.0, 0, 0
    px = np.atleast_1d(table["x"]).astype(float)
    py = np.atleast_1d(table["y"]).astype(float)
    dist = np.hypot(px[:, None] - truth["cx"][None, :], py[:, None] - truth["cy"][None, :])
    nearest = dist.argmin(axis=1)
    good = dist[np.arange(len(px)), nearest] <= match
    found = np.zeros(len(truth["cx"]), dtype=bool)
    found[nearest[good]] = True
    recall = float(found[findable].mean()) if findable.any() else 0.0
    return recall, int((~good).sum()), len(px)


def main():
    cp = os.pathsep.join([".", CLASSES, classpath()])
    rows = []

    for nspots in DENSITIES:
        for sigma in SIGMAS:
            for frames in FRAME_COUNTS:
                directory = f"tol10_{nspots}_{sigma}_{frames}"
                if not os.path.exists(os.path.join(directory, "sim.tif")):
                    subprocess.run(
                        [PYTHON, GEN, str(sigma), str(MEDIAN_N), directory, str(SEED)],
                        check=True, stdout=subprocess.DEVNULL,
                        env=dict(os.environ, CALIB_SPOTS=str(nspots),
                                 CALIB_FRAMES=str(frames)))

                truth = np.genfromtxt(os.path.join(directory, "truth.csv"),
                                      delimiter=",", names=True)
                findable = truth["snr"] >= THRESHOLD
                match = match_radius(sigma)

                noise = pixel_noise(summed_average(os.path.join(directory, "sim.tif"), frames))

                # 0 asks the plugin to estimate; the K sweep forces the value instead, so that
                # one run per K measures the constant rather than only the shipped choice.
                trials = [("fixed", FIXED)] + [(f"k{k}", k * noise) for k in KS]
                for label, tol in trials:
                    subprocess.run(
                        [JAVA, "-Djava.awt.headless=true", "-Xmx4g", "-cp", cp, "RunTune",
                         directory, str(sigma), f"{tol:.4f}", "1.0", str(frames),
                         str(THRESHOLD)],
                        check=True, capture_output=True, text=True)
                    recall, fp, n = score(directory, truth, findable, match)
                    rows.append({"spots": nspots, "sigma": sigma, "frames": frames,
                                 "noise": round(noise, 4), "which": label,
                                 "tol": round(tol, 3), "recall": round(recall, 4),
                                 "fp": fp, "found": n})
                    print(f"{nspots:4d} sigma={sigma} frames={frames:2d} noise={noise:6.3f} "
                          f"{label:>6} tol={tol:6.2f} recall={recall:.3f} fp={fp}",
                          flush=True)

    with open("tolerance_sweep.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote tolerance_sweep.csv, {len(rows)} rows")


if __name__ == "__main__":
    sys.exit(main())
