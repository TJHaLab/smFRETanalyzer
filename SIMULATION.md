# How the defaults were chosen

Almost every default in this plugin was set by running the analysis against
simulated data whose right answer is known, rather than by tuning it until
real data looked reasonable. This file records what was measured and what the
numbers were. Nothing here is needed to use the plugin - the
[README](README.md) is the reference for that - but if you are wondering
whether a default suits your data, this is where to find out what it assumes.

The simulator is a separate project,
[smfret-simulator](https://github.com/HazenBabcock/smfret-simulator). It
generates two-channel movies with a known PSF, a known background and known
per-molecule intensities, writes them in the format this plugin reads, and can
run the real plugin over them headless.

## Why simulation rather than real data

Real data has no ground truth. You cannot ask a real movie which of its
detections are single molecules, what the background truly was under a spot,
or how much of a trace belongs to the molecule next to it. Every one of those
is exactly the quantity a default has to be chosen against.

Two methodology points cost real mistakes before they were understood, and
they apply to any study of this kind:

* **Hold signal-to-noise constant, not brightness.** Varying spot size at
  fixed molecular brightness makes the large-sigma fields sit at the detection
  limit - detections fell from 386 to 121 going from SpotSigma 1 to 3 - so a
  filter looks useless there when the field is simply too dim to see.
* **Scale the match radius with spot size.** MaximumFinder reports the integer
  coordinate of the brightest pixel, so localisation error grows with the PSF.
  A fixed 1.5 pixel radius for calling a detection "correct" scored 20% of the
  real detections as false positives at SpotSigma 2.5.

## SpotThreshold - default 6

A signal-to-noise ratio, in real standard deviations, so 6 means a 6 sigma
detection.

Detection is not what limits it. In simulation at 30 averaged frames there are
no false positives at any threshold at or above 3, and about half the spots at
a given true SNR pass a threshold set to that value, which is what should
happen. Treat it as a trace quality bar rather than a detection cut.

It is meant to mean the same thing at every spot size, and that was checked:
on simulated movies from SpotSigma 1.0 to 3.0, analysed at the matching
SpotSigma, the reported SNR of a spot at fixed true significance varies by
**0.5%**. Getting there required correcting a bias in the background estimator
- before that it drifted 15% across the same range, and 46% at fixed molecular
brightness, which is the physically real case.

A fixed threshold still selects fewer molecules as spots get larger. That is
real rather than an artifact: a wider PSF spreads the same photons over more
background, so the same molecule genuinely is harder to detect.

## SpotTolerance - default 5

The flood fill depth MaximumFinder uses, which behaves as a minimum peak
amplitude: a spot whose peak sits less than this above the plain connecting it
to a brighter neighbour is merged into that neighbour and lost.

Swept from 1 to 20 against ground truth at spot sizes 1.0 to 3.0, two
densities and 20 averaged frames, scored on how many real spots survive every
filter. **5 recovers 98.3% of the best achievable on average and 96.0% in the
worst condition; 10 manages 79.1% and 43.6%.**

Both directions do damage, and the low end is worse than it looks. Every
maximum builds the masks, including the noise ones - they trip the spot
proximity filter and cull their real neighbours, and they mask background
pixels out of the estimate. At tolerance 1 with 400 real spots, 6208 maxima
were found and **3 spots survived**, with the background error at 44.5 ADU.
Too high and real spots are merged away, and the background degrades again
because the missed spots go unmasked.

It is not truly size independent - the useful upper limit falls from about 20
at SpotSigma 1 to about 4 at SpotSigma 3, because peak amplitude goes as
1/SpotSigma^2 - but 5 is within 4% of the best at every size tested. If you
average many more or many fewer than 20 frames, or your background is very
different, the noise floor moves and this is the parameter to revisit.

## SpotContamination - default 0.20

This is the newest default and the one with the most work behind it. The full
study is in
[spot-quality-ml](https://github.com/HazenBabcock/spot-quality-ml); what
follows is the summary.

### What the score predicts

The pipeline measures a trace by blurring the frame and reading one pixel, so
a molecule a few pixels away contributes to a trace that is not its own. That
measurement is linear in the image, and the image is a sum of per-molecule
contributions, so in simulation the signal at a spot splits **exactly** into
one term per molecule. The quantity being predicted is

```
contamination = 1 - (largest single molecule's contribution) / (total)
```

"what fraction of the signal measured here does not belong to its largest
contributor". A clean isolated molecule scores 0. Two molecules too close to
tell apart score about 0.5 each. A spot sitting on a dye aggregate scores
0.75 and up.

### Why a model rather than a formula

Given every molecule's position and brightness, a formula reproduces this
almost perfectly - R^2 of 0.98. Given only the spots the finder actually
reported, the same formula gets **R^2 of -0.18**. The difference is entirely
molecules the spot finder never detected: detection efficiency runs from about
90% on a sparse field down to 50% on a crowded one, and an unreported molecule
contaminates its neighbours exactly as much as a reported one does.

Their light is still in the image even though they are missing from the spot
list, which is why the score reads the pixels around each spot rather than
working from the spot list alone.

### What it uses

42 numbers per spot, all available the moment spot finding finishes: the
spot's own SNR and prominence, its flux and local background in
photoelectrons, the shape of its radial profile out to 12 pixels, second
moment and single-Gaussian-fit statistics, how many other spots sit within 5,
10 and 15 pixels, and what the detected spots alone predict should be landing
on it. These feed a small gradient boosted ensemble - 100 trees of depth 3,
about 1200 nodes - trained on 55308 simulated spots across 240 fields.

### How well it works

Held out by density, so every test spot is at a crowding level the model never
trained on:

| | rank correlation with the truth, within a field | AUC for finding a spot over 10% contaminated |
|---|---|---|
| prominence alone | 0.666 | 0.896 |
| **the model** | **0.831** | **0.971** |

The within-field number is the honest one - it is computed inside each field
and then averaged, so a model that merely learned which fields are crowded
would score nothing for it.

It is also **calibrated**, which is why the number is worth showing rather
than just ranking on. Spots predicted at 0.07 average 0.071 in truth; at 0.14,
0.140; at 0.29, 0.296. And it identifies dye aggregates at an AUC of 0.991.

### Choosing 0.20

The default follows the same principle as the prominence filter it replaces:
favour keeping spots, at the cost of letting more contaminated ones through.

| threshold | clean spots lost | contaminated spots caught |
|---|---|---|
| 0.15 | 0.78% | 73% |
| **0.20** | **0.27%** | **60%** |
| 0.25 | 0.10% | 50% |
| 0.30 | 0.05% | 39% |
| *old prominence filter* | *0.32%* | *25%* |

0.20 loses slightly **fewer** clean spots than the prominence default did
while catching two and a half times as many contaminated ones. On a field at
the density of the example data it rejects about 10% of spots and loses no
clean ones at all.

### What it does not cover

* **The PSF.** A model trained at one PSF does not transfer to another - at
  0.5 waves of aberration and a wide core it falls *below* untrained
  prominence. The shipped model is trained across sigma 1.0 to 2.5 and 0 to
  0.5 waves and beats prominence everywhere in that range, but a badly
  aberrated system outside it is untested. Two attempts to fix this by telling
  the model what the PSF is - recomputing every feature at the measured sigma,
  and handing it the fitted sigma and aberration as inputs - both changed
  nothing. The PSF dependence is in the relationship itself, not in how the
  numbers are scaled.
* **Anything below the resolution limit.** Where a spot's contamination comes
  from molecules the finder reported, the model is accurate to about half a
  percent. Where it comes from molecules nothing reported, it is off by four
  or five percent. A molecule sitting exactly on top of another is absent from
  the image as well as from the spot list, and no model recovers it. That is a
  floor, not a shortfall.
* **Real data.** The quantity being predicted cannot be measured on a real
  movie, so the score has never been validated against one directly. The
  nearest indirect check is that a flagged trace should be more likely to show
  two-step photobleaching.

## BackgroundKappa - default 1.3, set 0 to use it

How far above the background estimate a pixel may sit, in robust standard
deviations, before it is treated as spot light and excluded.

Swept from 0.4 to 10 against a known background at spot sizes 1 to 3 and at
the spot density of real data, scored on the bias of the estimate **at spot
locations**, since that is what gets subtracted from a trace. It does not
depend on spot size: aggregated over sigma the best value is flat from about
1.0 to 1.6, and within that range it is not even monotone. On the example data
the change from the old size-dependent rule is invisible - background 17.441
against 17.466, 421 spots against 418.

An earlier version used `1.5 - 0.3 * SpotSigma`. That line was fitted while
the estimator still had a bias of its own, so it was partly compensating for
the estimator rather than describing the spots.

Clipping only from above leaves the estimate reading low by about an ADU,
which is enough to matter - one ADU of background is roughly 50 units of
trace. The estimator puts that bias back analytically, so this parameter
controls which pixels are excluded without also shifting the level.

## Spot masking radius and background smoothing

Not user-facing, but derived the same way and worth recording.

The **masking radius** does scale with spot size: `max(2, round(2.5 *
SpotSigma))`, swept from 2 to 10 pixels at sigma 1.0 to 3.0 against a known
background. Against the old fixed 4 this cuts the bias at spot locations from
0.225 to 0.098 ADU at 400 spots per channel - worth about 15 units of trace at
SpotSigma 3, systematically.

It scales at 2.5 sigma rather than the 3 sigma a PSF actually reaches, because
masking is not free: every masked pixel is one the estimator cannot use. At
sigma 3 the error bottoms at a radius of 8 and gets *worse* by 10 even as the
bias keeps falling.

The **smoothing scale** showed no trend against spot size at all and is a
constant 14 pixels.

## The SpotSigma that smFRETPSFVisualizer recommends

smFRETAnalyzer measures a trace with a **Gaussian** matched filter, while a
real PSF is not a Gaussian. An aberrated Airy pattern carries a real fraction
of its light out where a Gaussian of the same core has none, so a filter
matched to the core alone throws that away and the best width is wider than
the fitted core. Working the matched filter through gives a signal-to-noise
proportional to

```
(1 / s) * integral p(r) exp(-r^2 / 2 s^2) r dr
```

for a filter of width `s` and a PSF `p(r)`. Brightness and background cancel,
so **where this peaks is a property of the PSF's shape alone**. For a Gaussian
PSF it peaks at exactly the PSF's own width, which is the standard matched
filter result and is what says the derivation is right rather than merely
plausible.

Aberration is the only thing that moves it, and the answer scales with the
PSF, so the whole two-dimensional question collapses to one multiplier:

| spherical aberration | 0 | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 waves |
| --- | --- | --- | --- | --- | --- | --- |
| best SpotSigma / fitted sigma | 0.98 | 0.99 | 1.02 | 1.08 | 1.21 | 1.49 |

Below about 0.3 waves the fitted core is already the right filter.

**The band matters more than the number**, because a matched filter's response
is flat around its peak. On the example data anything from 1.4 to 2.2 is
within 2% of the best, which is why `spotsigma=2` in the example macros was
never costing anything measurable even though the fitted core is 1.34.

The prediction was checked against the pipeline rather than only derived: on a
simulated aberrated field of known PSF, sweeping the real trace extraction
from SpotSigma 1.0 to 3.0 puts the measured peak at **1.80 where the formula
predicts 1.81** - so pixel binning, the sub-pixel centring error and an
estimated background together move it by less than the sweep can resolve.

Two things it does not account for: the noise is taken as background
dominated, which is the same assumption the spot finder's SNR makes, and
SpotSigma also drives spot *finding*, so re-running with a new value changes
which molecules are found.

## The PSF measurement

Two choices in smFRETPSFVisualizer were made against simulated data.

**Contaminated pixels are masked, not contaminated spots discarded.** Throwing
away every molecule with a close neighbour is the obvious approach and it does
not work on a crowded field: on a simulated field at real density it keeps 100
of 235 spots and still returns the wrong answer, 1.53 sigma and 0.36 waves
against a truth of 1.45 and 0.41. Masking the pixels keeps 234 spots and
returns 1.49 and 0.41.

**The fitted pedestal is real.** On the example data, fitting a flat term
under the PSF cuts the residual more than fourfold and leaves about 3% of
peak, flat from 5 to 10 pixels where any real PSF wing would still be
decaying. Run the same measurement on simulated data, which has no scattered
light halo, and the pedestal comes out near zero - which is what makes the 3%
on real data worth believing rather than an artifact of the method.

## Things that were tried and rejected

Recorded so they are not re-tried without new evidence.

* **A shape-based filter for close pairs.** SNR-normalised ellipticity is a
  genuinely better statistic than prominence for separating doublets - 0.473
  against 0.538 on the scoring used - but at a single constant threshold it
  buys 5% on a clean field against 38% on a dirty one, and the two optima are
  far enough apart that the dirty-field setting does real damage to clean
  data. See [ADR-0001](docs/adr/0001-no-shape-based-doublet-rejection.md). The
  contamination score is in a sense the resolution of this: a continuous
  per-spot number rather than one global threshold.
* **Seeding the background clip from the previous frame.** 10% faster and a
  different answer - traces shifted by a median of 14%, because the clip only
  ever removes pixels and a seeded level is already clipped.
* **FFT convolution for the background smoothing.** Measured against three
  alternatives; direct convolution on a decimated grid won by a factor of two
  and FFT does not overtake it until about sigma 40.
