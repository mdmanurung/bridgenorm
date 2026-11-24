# bridgenorm

Batch correction for cytometry from bridge samples: reference material from the same donor measured in every batch.

Site and articles: <https://mdmanurung.github.io/bridgenorm/>

The differences between bridge samples measure the technical batch effect. bridgenorm estimates the batch effect from them alone, so study samples never enter a fit and their biology cannot be removed as batch effect. Each batch is corrected toward one anchor batch. A fitted model is frozen with a digest, so it applies to new data without changing earlier corrections.

## Methods

`bn_fit()` is the main method. It works on reference cells only, in two passes:

1. A self-organising map (SOM) on the mapping markers splits cells into nodes (cell states). Each batch gets a paired shift per node, which is the mean bridge-sample difference from the anchor.
2. A second map is fitted on the shifted cells. In each (batch, node), the batch's marker quantiles are mapped onto the anchor's. This corrects changes in the shape of a distribution, not only its location. The offset curve of a small node is shrunk toward the batch-wide curve with weight n / (n + `shrink_k`).

Only batches whose bridge donors equal the anchor's get quantile maps. A batch that shares no bridge donor with the anchor is left raw and flagged.

Two options are off by default. Estimates at the defaults are unchanged from 0.9.0.

- `target = "virtual"` corrects every batch, the anchor included, toward a virtual reference. Per SOM node, the reference is the equal-weight mean of the batches' quantile curves. It lowered the held-out score on the full data set A data (0.955 to 0.906), but the gain did not hold on subsets: with 8 of the 11 scorable batches it was worse on every data set (+0.03 to +0.05), and leaving out one study changed its direction. The default stays `target = "anchor"`. It does not combine with `shrink_k = "eb"`.
- `shrink_k = "eb"` estimates one pseudo-count per marker by empirical Bayes, never below 300. On data set A (fresh sim seeds 1008-1009) it lowered node5 RMSE by 4.9-6.5% but made confounded-expression recovery slightly worse. On the data set B tuning round it was worse than the matching scalar `shrink_k` (0.0521 vs 0.0444 at grid 8x8, 0.0439 vs 0.0303 at grid 6x6; |1 - recovery|, lower is better).

`bn_fit_paired()` fits the baseline: one paired shift per population and marker, with no nodes and no quantile maps. Every other model must beat it.

The anchor, SOM grid, mapping markers, seed and transform depend on the data and have no default. The other settings default to the values benchmarked on two spectral cytometry data sets (data set A and data set B).

## Comparison with CytoNorm and cyCombine

Held-out score: the median over markers of the worst batch-pair EMD on held-out bridge wells (3 wells; lower is better).

| Data set | bridgenorm | CytoNorm 2.0.10 | cyCombine 0.3.0 | Paired shift | Uncorrected |
|---|---|---|---|---|---|
| data set A (heb 0.7.0) | 0.954 | 1.592 | 3.533 | 2.376 | 7.307 |
| data set B (0.9.0 defaults; grid 6x6) | 0.797; 0.766 | 0.805 | 1.134 | 1.442 | 2.454 |

bridgenorm scored lowest on both data sets. On data set B its lead over CytoNorm is within the spread across wells (SD about 0.2 for both). cyCombine was fitted on the scored cells, which favours it. CytoNorm ran with its defaults (FlowSOM 15x15, 10 clusters).

On the data set B confounded simulation, bridgenorm also recovered expression differences best (|1 - recovery| 0.030 at grid 6x6, 0.044 at defaults, against 0.050 for CytoNorm and 0.118 for cyCombine). These sim seeds also tuned the grid 6x6 setting, so they favour it; a rerun on unused seeds is planned. Cell RMSE to the simulated truth favours bridgenorm by design, because the truth is the anchor batch and CytoNorm targets the batch mean.

The comparison used `target = "anchor"`. The virtual target has not been compared with CytoNorm or cyCombine.

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
