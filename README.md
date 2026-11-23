# bridgenorm

Batch correction for cytometry from bridge samples: reference material from the same donor measured in every batch.

The differences between bridge samples measure the technical batch effect. bridgenorm estimates the batch effect from them alone, so study samples never enter a fit and their biology cannot be removed as batch effect. Each batch is corrected toward one anchor batch. A fitted model is frozen with a digest, so it applies to new data without changing earlier corrections.

## Methods

`bn_fit_node_quantile()` is the main method. It works on reference cells only, in two passes:

1. A self-organising map (SOM) on the mapping markers splits cells into nodes (cell states). Each batch gets a paired shift per node, which is the mean bridge-sample difference from the anchor (`bn_fit_nodes()`).
2. A second map is fitted on the shifted cells. In each (batch, node), the batch's marker quantiles are mapped onto the anchor's. This corrects changes in the shape of a distribution, not only its location. The offset curve of a small node is shrunk toward the batch-wide curve with weight n / (n + `shrink_k`).

Only batches whose bridge donors equal the anchor's get quantile maps. A batch that shares no bridge donor with the anchor is left raw and flagged.

Simpler models are available as baselines:

- `bn_fit_paired()`: one paired shift per population and marker. Every other model must beat this baseline.
- `bn_fit()`: a linear mixed model per population and marker. Donor effects are fixed, and study and batch effects are random, with variance components from REML. A batch's correction is its batch-effect BLUP minus the anchor's.
- `bn_fit_shared()`: one mixed model per marker, with a batch effect shared across populations.

Every setting that can change an estimate is a required argument with no default.

## Install

```r
remotes::install_github("mdmanurung/bridgenorm")
```

## Use

```r
samples <- data.frame(
  sample_id = ..., batch = ...,
  reference_id = ..., # bridge donor; NA for study samples
  aliquot = ...       # optional
)
bn_bridges(samples, anchor = "b1")

model <- bn_fit_node_quantile(
  x, samples, population = labels, markers = markers,
  transform = list(name = "asinh", cofactor = 6000), anchor = "b1",
  min_cells = 50, unbridged = "abstain", method = "paired", eb_settings = NULL,
  grid = c(6, 6), mapping_markers = mapping, train_cells = 5000, rlen = 100,
  seed = 1, assignment = "nearest", fallback = "global",
  quantile_probs = seq(0.025, 0.975, length.out = 21), tail = "shift",
  shrink_k = 300
)
res <- bn_integrate(model, x, samples, labels, transform = model$transform, unsupported = "abstain")
bn_save(model, "bn_model.rds")
```

The values above show the shape of a call. They are not recommendations. See the Get started article for a worked example.

## Origin

bridgenorm was called heb (hierarchical empirical Bayes) up to version 0.7.2. The mixed-model fits are ported from the `heb` and `heb_shared` engines of SpectraWeaveR (`R/corrector_heb.R`, `R/corrector_heb_shared.R`, `R/corrector_paired_shift.R`, `R/corrector.R`). On the SpectraWeaveR test simulations, all 32 combinations of `pool`, `anchor_effect`, `variance_prior`, `residual_model` and `scale_model` give the same effects and corrected values as SpectraWeaveR, to the last bit. The internal-anchor option was left out: it failed its real-data check (plan A14) and was dropped on 2026-11-20.

One intended difference: with `unbridged = "abstain"`, a batch whose references were never measured in the anchor batch gets no effect. SpectraWeaveR still fitted an effect for such a batch through its donors in other batches, although its documentation said it would not.
