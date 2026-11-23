# bridgenorm

Batch correction for cytometry from bridge samples: reference material from the same donor measured in every batch.

The differences between bridge samples measure the technical batch effect. bridgenorm estimates the batch effect from them alone, so study samples never enter a fit and their biology cannot be removed as batch effect. Each batch is corrected toward one anchor batch. A fitted model is frozen with a digest, so it applies to new data without changing earlier corrections.

## Methods

`bn_fit()` is the main method. It works on reference cells only, in two passes:

1. A self-organising map (SOM) on the mapping markers splits cells into nodes (cell states). Each batch gets a paired shift per node, which is the mean bridge-sample difference from the anchor.
2. A second map is fitted on the shifted cells. In each (batch, node), the batch's marker quantiles are mapped onto the anchor's. This corrects changes in the shape of a distribution, not only its location. The offset curve of a small node is shrunk toward the batch-wide curve with weight n / (n + `shrink_k`).

Only batches whose bridge donors equal the anchor's get quantile maps. A batch that shares no bridge donor with the anchor is left raw and flagged.

`bn_fit_paired()` fits the baseline: one paired shift per population and marker, with no nodes and no quantile maps. Every other model must beat it.

The anchor, SOM grid, mapping markers, seed and transform depend on the data and have no default. The other settings default to the values benchmarked on two spectral cytometry data sets (data set A and data set B).

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

model <- bn_fit(
  x, samples, population = labels, markers = markers,
  transform = list(name = "asinh", cofactor = 6000), anchor = "b1",
  grid = c(6, 6), mapping_markers = mapping, seed = 1
)
res <- bn_apply(model, x, samples, labels)
bn_save(model, "bn_model.rds")
```

The grid and mapping markers above show the shape of a call; choose them for your panel. See the Get started article for a worked example.

## Origin

bridgenorm was called heb (hierarchical empirical Bayes) up to version 0.7.2. Version 0.9.0 removed the mixed-model fits ported from SpectraWeaveR, because no benchmarked method used them; experiments that need them use heb 0.7.x.
