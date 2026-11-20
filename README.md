# heb

Hierarchical empirical Bayes batch correction for cytometry, estimated from repeated reference samples.

For each population and marker, `heb_fit()` fits a linear mixed model to aliquot means of the reference samples. Donor effects are fixed. Study and batch effects are random. Variance components come from REML. The correction of a batch is the BLUP of its batch effect minus that of the anchor batch. `heb_fit_shared()` fits one model per marker instead, with a batch effect shared across populations plus a per-population deviation.

`heb_integrate()` applies a frozen model. A batch the model has not seen gets a paired-reference shift against the frozen anchor means. Earlier effects never change. A batch with no reference shared with the anchor is left raw and flagged.

Every setting that can change an estimate is a required argument with no default.

## Install

```r
remotes::install_github("mdmanurung/bridgenorm")
```

## Use

```r
samples <- data.frame(
  sample_id = ..., batch = ...,
  reference_id = ..., # NA for samples that are not references
  aliquot = ...       # optional
)
heb_bridges(samples, anchor = "b1")

model <- heb_fit(
  x, samples, population = labels, markers = c("CD3", "CD4"),
  transform = list(name = "asinh", cofactor = 6000), anchor = "b1",
  hierarchy = data.frame(batch = ..., study = ...), min_cells = 50,
  unbridged = "abstain", variance_floor = 1e-8, pool = "population",
  anchor_effect = "random", variance_prior = "boundary_avoiding",
  residual_model = "known_plus_extra", scale_model = "none"
)
res <- heb_integrate(model, x, samples, labels, transform = model$transform, unsupported = "abstain")
heb_save(model, "heb_model.rds")
```

The values above show the shape of a call. They are not recommendations.

## Origin

Ported from the `heb` and `heb_shared` engines of SpectraWeaveR (`R/corrector_heb.R`, `R/corrector_heb_shared.R`, `R/corrector_paired_shift.R`, `R/corrector.R`). The internal-anchor option was left out: it failed its real-data check (plan A14) and was dropped on 2026-11-20. On the SpectraWeaveR test simulations, all 32 combinations of `pool`, `anchor_effect`, `variance_prior`, `residual_model` and `scale_model` give the same effects and corrected values as SpectraWeaveR, to the last bit.

One intended difference: with `unbridged = "abstain"`, a batch whose references were never measured in the anchor batch gets no effect. SpectraWeaveR still fitted an effect for such a batch through its donors in other batches, although its documentation said it would not.
