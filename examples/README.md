# Example data

`hel1.tif` is a real two-channel smFRET movie of a DNA-only sample: 512 x 512, 8 bit, 30 frames,
donor in the left half and acceptor in the right, which is the layout every plugin here assumes.

It is enough on its own to run the whole pipeline. Each stage writes its outputs next to the
input, so work on a copy if you want to keep the folder clean.

1. **Plugins > smFRET > smFRET Channel Mapping** on `hel1.tif`, frames 1 to 30. Writes
   `hel1_mapping.json`. This is the one stage that needs [TurboReg](https://imagej.net/plugins/turboreg)
   installed.
2. **smFRET Spot Finder** on `hel1.tif` plus the mapping file from step 1. Writes
   `hel1_spotf_finding.json` and three companion files. The defaults are tuned; see the
   parameter reference in [../README.md](../README.md) before changing them.
3. **smFRET Time Traces** on `hel1_spotf_finding.json`. Writes `hel1.h5` and `hel1.traces`.

Only the two JSONs and `hel1.h5` land beside the movie; the rest go in a `hel1_analysis` folder
created alongside it. See [../README.md](../README.md) for the full layout.

The two viewers and the PSF visualizer then read those outputs. smFRETPSFVisualizer only needs
the spot-finder JSON, so it can be run straight after step 2.

One thing to know before you look at the traces: intensities collapse over the last few frames
because **the illumination laser is switched off**, not because the molecules photobleached.
