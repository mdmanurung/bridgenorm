# bridgenorm 0.8.0

* The package is renamed from heb to bridgenorm. Exported functions, classes and schema strings change prefix
  from `heb_` to `bn_` (`heb_fit_node_quantile()` is now `bn_fit_node_quantile()`, and so on). The `heb_settings`
  argument of `bn_fit_nodes()` is now `eb_settings`, and the mixed-model method value `"heb"` is now `"eb"`
  (`"heb_shared"` is now `"eb_shared"`). Models saved by heb do not load; refit them. Estimates are unchanged.

# heb 0.7.2

* `heb_fit_node_quantile()` builds quantile tables only for batches whose reference donors (`reference_id`) are
  exactly the anchor's. Before, every batch got tables, so an unbridged batch that pass 1 left raw was still
  quantile-mapped onto the anchor's donors, and a partly bridged batch pooled donors the anchor lacks. Such a
  batch now takes the pass-2 node shift if bridged and stays raw if not, as the README states. Batches in the
  data set A and data set B benchmarks all share the anchor's donors, so their results are unchanged (decision RV2).

# heb 0.7.1

* `heb_fit_node_quantile()` sorts each column of a shrunk quantile table (rearrangement; Chernozhukov,
  Fernandez-Val and Galichon 2010, Econometrica 78:1093). With `shrink_k = 300` the mix of node, anchor and
  batch-wide curves left about 0.4% of table columns non-monotone (worst step -0.037 asinh), so a few cells
  swapped order. The fix changed no benchmark score by more than 0.0007 at SOM seeds 473-475 (tusoskill run
  heb_round3b_20261122, candidate c016; decision R3b12). `shrink_k = 0` is unchanged.

# heb 0.7.0

* `heb_fit_node_quantile()` gains a required `shrink_k`: each non-anchor (batch, pass-2 node) offset curve
  (batch minus anchor quantiles) is pulled toward the batch-wide offset curve with weight n / (n + shrink_k),
  n = the batch's reference cells in the node (tusoskill run heb_confounded_20261121, candidate c018;
  decision TR3k, shrink_k = 300). `shrink_k = 0` reproduces 0.6.0. **Breaking:** calls without `shrink_k`
  now stop. In the round 3 lockbox (sim seeds 1004-1005, SOM seeds 473-475) shrink_k = 300 lowered recovery
  distance on every seed (-0.009 to -0.013) and RMSE to truth, with no expression harm.

# heb 0.6.0

* New `heb_fit_node_quantile()`: two-pass node correction with per-node quantile maps (tusoskill run
  heb_nodes_20261121, candidate c016; decisions T1-T3). Pass 1 is `heb_fit_nodes()`; pass 2 refits the
  map on the shifted reference cells and, in each node where a batch and the anchor both have at least
  `min_cells` reference cells pooled over wells, maps the batch's marker quantiles onto the anchor's
  (linear interpolation, a constant shift beyond the outer quantiles). Other cells take the pass-2 node
  shift. `quantile_probs` and `tail` are required, like every `heb_fit_nodes()` setting.
  On the data set A lockbox simulation it lowered the worst-pair EMD from 0.257 to 0.198 (node-specific
  shifts) and from 0.763 to 0.585 (smooth shifts) against the 6x6 node shift.
* `heb_integrate()` applies the new `heb_node_quantile_model`; its `actions` carry a `pass` column, and
  pass-2 rows mapped by quantiles read `apply / node_quantile`. `heb_save()` and `heb_load()` accept it.
* `heb_fit_nodes()` and its corrections are unchanged.

# heb 0.5.1

Fixes from an independent review (2026-11-21); corrected values are unchanged for valid inputs.

* `heb_fit_nodes()` stops unless every mapping marker is one of `markers`. A mapping marker outside
  them was mapped without the global pre-shift, so cells of a shifted batch landed in the wrong node.
* `heb_integrate()` on a node model: a (batch, node, marker) row reads `apply / node_fallback_global`
  only when every cell in it took the global effect. When the global model covers some of the node's
  populations but not others, the row is `abstain / node_fallback_partial` (the other cells keep raw
  values), so `unsupported = "error"` stops. It used to read `apply`.
* `heb_integrate()` on a node model returns a copy that keeps the effects of batches estimated at
  integrate time, with `parent_digest`, as for population models. It used to return the input model.
* `heb_fit_nodes()` restores the caller's random number stream on exit, and stores the kohonen version
  and its default settings in `settings$som`. The help page states that the map depends on row order
  within a sample for a given seed.
* Batch shrinkage adds `batch_shrinkage_status`: `"single_batch_group"` marks a group with too few
  batches to moderate, which keeps its paired contrast (weight 1); the value is unchanged.
* `n_aliquots` and `n_aliquots_anchor` count aliquots per donor; labels shared across donors were
  counted once.
* A missing batch or study label in `hierarchy` stops with a clear message, instead of a
  "not positive definite" error.

# heb 0.5.0

* New `heb_fit_nodes()`: node-level correction (roadmap R8). A self-organising map (kohonen, now in
  Suggests) is trained on reference cells after the global paired shift, and paired shift or `heb_fit()`
  is fitted per map node. `heb_integrate()` assigns each cell to its nearest node and subtracts that
  node's effect; a node without an effect takes the global effect (`fallback = "global"`) or keeps its
  raw value (`"abstain"`). `heb_save()` and `heb_load()` accept the new `heb_node_model`. All settings
  (`grid`, `mapping_markers`, `train_cells`, `rlen`, `seed`, `assignment`, `fallback`) are required.

# heb 0.4.0

* `heb_fit()` gains shrinkage across batches (roadmap R2), four required settings:
  `batch_shrinkage` ("blup", "standard_eb", "adaptive"), `batch_shrinkage_centre` ("study", "zero"),
  `batch_shrinkage_group` ("study", "all") and `batch_shrinkage_tau2` ("profile", "reml"). Non-"blup"
  values replace the BLUP by the unshrunk paired contrast moderated across batches, with cell plus
  pooled well variance; the BLUP is kept as `estimate_blup`. `"blup"` gives the 0.3.0 effects.

# heb 0.3.0

* `heb_fit()` gains three required settings, `moderation`, `moderation_centre` and
  `moderation_variance`: each batch's marker effects can be shrunk across markers by standard
  empirical Bayes or the adaptive weight of ComBat-refQL (`fit_adaptive_eb`, MIT). Moderated fits
  keep `estimate_raw` and add `moderation_weight`, `moderation_tau2`, `moderation_centre` and
  `moderation_status`. `moderation = "none"` gives the 0.2.0 effects.
* Review fixes (independent review, 2026-11-20):
  * moderation leaves a marker with zero sampling variance unmoderated, where it used to crash;
  * the prior variance is found by a log-grid scan plus a log-scale refine, so a second basin or
    small-unit data no longer stop it at a local optimum or at zero;
  * a (batch, population) with no reference donor shared with the anchor in that population gets
    no effect (it used to get one through other populations' bridges);
  * the `log_sd` scale centre uses only aliquots that pass `min_cells`;
  * `heb_fit()` needs `min_cells >= 2`;
  * new-batch paired rows of a moderated model record identity moderation.

# heb 0.2.0

* New `heb_fit_paired()`: the paired-shift baseline as a frozen model.
* New `heb_reference_noise()`: reference noise per batch to inform the anchor choice
  (report only; roadmap R7).
* A batch whose references all fall below `min_cells` now counts as unbridged.
* A singular mixed-model system stops with the population and marker named.
* `heb_fit_shared()` works with one donor and one population.
