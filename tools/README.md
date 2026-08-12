# Calibration tools

The simulations and sweeps that set the spot finding and background constants. They are here so
the numbers quoted in `CLAUDE.md` and `docs/adr/` can be re-derived rather than taken on trust —
every default they produced is stated there as fact, and none of it is reproducible without
these.

Nothing here is part of the plugin. The Maven build ignores this directory.

These are not the test suite, and the two do different jobs. `src/test/java` pins the constants
and the arithmetic against known answers on generated fields, in ten seconds and with no data —
it tells you something *broke*. The sweeps here tell you what a constant should *be*, which needs
hours and large simulated movies. Re-derive with these; guard with those.

## Requirements

Python 3 with `numpy`, `scipy` and `tifffile`, and a built plugin plus its dependency classpath:

```bash
JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64 mvn -o package
mvn -o dependency:build-classpath -Dmdep.outputFile=tools/cp.txt
```

`tools/cp.txt` is absolute paths into `~/.m2`, so it is machine specific and gitignored. Set
`SMFRET_CP` instead if you prefer, and `SMFRET_JAVA` / `SMFRET_JAVAC` if your JDK 11 is
elsewhere. Run the sweeps from a scratch directory — they write simulated movies, which are
large.

## The generators

| | |
|---|---|
| `gen_calib.py` | Random field at real density, log-normal brightness, Gaussian beam. Writes the true background, so the estimator can be scored against a known answer. Behind `backgroundKappa`, `spotMargin`, `spotTolerance` and `spotThreshold`. |
| `gen_shape.py` | Adds a **labelled** population — singles, doublets at 0.25–3σ, aggregates — so a filter can be scored on what it rejects. Behind `spotProminence` (since retired) and ADR-0001. |

**Two methodology traps, both of which produced wrong answers first:**

1. **Hold the SNR constant, not the brightness.** Fixing `median N` while varying σ puts the
   large-σ fields at the detection limit (detections fell 386 → 121 from σ=1 to 3), so a filter
   looks useless there when the field is simply too dim. Pass `N` proportional to σ.
2. **Scale the match radius with σ.** `MaximumFinder` returns the integer coordinate of the
   brightest pixel, so localization error grows with the PSF — median 0.82 px at σ=1 against
   1.26 px at σ=3. A fixed 1.5 px radius scored 20% of real detections as false positives at
   σ=2.5. `tune_sweep.py` uses `min(3, max(1.5, 1.2σ))`.

## The sweeps

```bash
python3 ../tools/tune_sweep.py && python3 ../tools/analyze_tune.py   # spotTolerance -> 5
python3 ../tools/prom_design.py                                      # spotProminence -> 0.4 (retired)
python3 ../tools/shape_compare.py                                    # ADR-0001
```

`spotProminence` no longer exists as a parameter — the contamination model replaced it as the
last-stage filter, and `SIMULATION.md` is where that default comes from. `prom_design.py` still
runs and still produces 0.4: prominence is measured on every spot and written to the spot table
as a column, it just no longer rejects anything on its own. It is kept because that number is
quoted as fact in `docs/adr/` and this is what re-derives it.

`prom_design.py` and `shape_compare.py` expect `shape_<sigma>` directories from `gen_shape.py`;
`tune_sweep.py` generates its own. The last-stage filter is swept offline from one run per field
with it disabled, since it affects nothing upstream: one run gives the score of every surviving
spot and any threshold can then be applied to the table.

**Disabling it means `RunTune ... 1.0`, not a negative number.** That argument was
`spotProminence` until the contamination model replaced it, and the two run in opposite
directions - prominence rejected *below* its value, contamination rejects *above* it, clamped to
[0, 1]. Carrying the old `-1000` sentinel over rejects every spot and scores an empty field,
which looks like a broken sweep rather than a wrong flag.

## The Java drivers

`RunTune.java` runs stage 2 on a simulated movie with the tolerance and contamination forced.
`Rig.java` runs stages 2 and 3 on a real movie, which is how the `.h5` the viewers need gets
rebuilt. Both are compiled into the working directory by `rig.py`, or by hand:

```bash
javac -cp target/classes:$(cat tools/cp.txt) -d . tools/RunTune.java
```

**Run the plugins with `-Dsmfret.diagnostics=true`.** The sweeps score the background estimate
by reading `_spotf_bg_smooth.tif` and `_spotf_fg_smooth.tif`, which the plugins only write when
that property is set - it is off for normal use. `tune_sweep.py` passes it already; anything
driving `RunTune` or `Rig` by hand has to.

**Everything the plugins generate now lands in `<movie>_analysis/`,** not beside the movie -
only the two JSONs and the trace `.h5` stay put. So a sweep run in `tol_400_2.0/` reads its spot
table from `tol_400_2.0/sim_analysis/sim_spotf_spots.csv`. The scripts here were updated with
the layout; a scratch directory from an older run has the files one level up and will not score.

## A caveat on the numbers

The example movies are 30 frames and their end-of-movie drop is **the illumination laser being
switched off**, not photobleaching. That makes the dark-frame level a valid zero reference but
means the short movies must not be used to judge background estimation — the temporal window
spans the switch-off. Simulation with a known answer is the primary evidence; the real movies
are a sanity check.
