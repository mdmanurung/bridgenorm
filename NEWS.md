# bridgenorm 1.1.0

## Default changes

* `bn_fit()` now defaults to `shrink_target = "value"` and `shrink_k = 1000` (1.0.0: `"quantile"` and 300).
  Corrected values change for every fit that relied on the defaults; pass `shrink_target = "quantile",
  shrink_k = 300` to reproduce 1.0.0. Models saved by 1.0.0 still load and apply unchanged.
* The new defaults were chosen by tuning on two real multi-batch cytometry studies, scored by the earth
  mover's distance between batches on held-out bridge wells. The comparison used every fold, including
  those used to pick the settings, so the gain is an upper bound. Check the defaults on new data.
* `settings$shrink_target` is now always recorded. A model without it used `"quantile"`.

## New

* `bn_fit(shrink_target = )`: `"value"` shrinks each node's offset curve toward the batch-wide map
  evaluated at the node's own target quantiles, so an intensity-dependent effect shared by all nodes is not
  shrunk toward a wrong curve. `"quantile"` keeps the 1.0.0 rule.
* `bn_fit(shrink_reweight = TRUE)` computes the batch-wide curve with weights that give each batch the
  target's cell-state mix of bridge cells. Default `FALSE`.
* `bn_fit(quantile_maps = FALSE)` keeps no quantile tables, so every cell takes the pass-2 node shift (an
  ablation). Default `TRUE`.

# bridgenorm 1.0.0

The model contract is now stable: exported functions, model classes, schema strings and the columns of
`actions` change only with a major version.

## Breaking changes

* Models saved by 0.x do not load; refit them. Schema strings are now `bn_node_quantile_model/2`,
  `bn_node_model/2` and `bn_paired_model/2`.
* Every model has the parent class `bn_model`. `bn_fit_paired()` returns class
  `c("bn_paired_model", "bn_model")` (was `bn_model`).
* `bn_fit()` drops the single-valued arguments `method`, `eb_settings`, `assignment` and `tail`. Their values
  are fixed and still recorded in `settings`.
* In the `actions` of a `bn_fit()` result, the SOM node column is `node` (was `population`).
* Dead fields are removed: `hierarchy` from paired models, the duplicate `global` model from node models,
  `eb_settings` and `assignment` from node settings.
* `seed` must be a whole number. `kohonen` moves from Suggests to Imports.

## New

* `goal_batches`: with `target = "virtual"`, `bn_fit()` and `bn_fit_paired()` map every batch toward the
  equal-weight mean of a chosen set of batches, for example the batches of one study, as CytoNorm's goal
  distribution. Unbridged or donor-mismatched goal batches are dropped; a goal lacking a reference summary the
  anchor has stops the fit. No benchmark benefit is claimed yet.
* `bn_fit()` warns when a batch gets no quantile maps (donors that pass `min_cells` differ from the anchor's)
  and when bridge donors' cell shares differ from the anchor's by more than 10 percentage points.
* `bn_apply()` sends a message naming batches it left entirely raw.

## Fixes

* `bn_apply()` no longer stops when a marker has one value in a node's quantile table: that marker keeps the
  pass-2 node shift (reason code `node_quantile_tied`).
* An empty or blank `sample_id`, `batch`, `reference_id` or `aliquot`, and carriage returns in ids or
  population labels, now stop with the samples named. Before, `""` silently made study samples a bridge donor.
* Quantile tables are built only for batches whose aliquots pass `min_cells`, the gate the shift passes use.
* The fit no longer depends on the collation locale or on the session's RNG kind. Results on data whose ids
  sort the same in every locale are unchanged.
* The digest leaves out timestamps and hashes a canonical form, so two fits on the same data have the same
  digest.
* `bn_reference_noise()` takes `transform` as optional. Choice errors name the function the user called.
  `bn_load()` checks that the file exists.

# bridgenorm 0.10.0

* No default changes. Estimates at the defaults (`shrink_k = 300`, `target = "anchor"`) are unchanged from 0.9.0.
  Both additions are options.
* `bn_fit(shrink_k = "eb")` estimates one pseudo-count per marker by empirical Bayes, with a floor of 300. The
  values are stored in `settings$shrink_k_marker`. It stays optional; in benchmarks it did not beat the scalar
  default.
* `bn_fit(target = "virtual")` and `bn_fit_paired(target = "virtual")` correct every batch, the anchor included,
  toward a virtual reference. In each SOM node the reference is the equal-weight mean of the batches' quantile
  curves (a 1D Wasserstein barycenter). The anchor still sets the bridges and the map.
  `target = "virtual"` does not take `shrink_k = "eb"`, and no batch may be named like the virtual key. In
  `bn_fit_paired()` each anchor reference summary is replaced by its mean over the bridged batches, so the anchor
  gets an effect too. `bn_apply()` reads the target from the model; models from 0.9.0 map to the anchor.

# bridgenorm 0.9.0

* The mixed-model fits are removed: `bn_fit()` (REML), `bn_fit_shared()`, moderation, batch shrinkage and the
  `"eb"` node method. No benchmarked method used them. Experiments that need them stay on heb 0.7.x.
* `bn_fit_node_quantile()` is now `bn_fit()`, the main method. `bn_integrate()` is now `bn_apply()`.
  `bn_fit_nodes()` is internal; it is the pass-1 and pass-2 model of `bn_fit()`.
* Settings that depend on the data stay required: `anchor`, `grid`, `mapping_markers`, `seed`, `transform`,
  `markers` and `population`. The rest default to the benchmarked values: `min_cells = 50`,
  `unbridged = "abstain"`, `train_cells = 500`, `rlen = 100`, 21 quantiles from 0.025 to 0.975, `tail = "shift"`,
  `shrink_k = 300`. `bn_fit_paired()` and `bn_reference_noise()` default to `min_cells = 50`; `bn_apply()`
  defaults to the model's transform and `unsupported = "abstain"`. Estimates are unchanged.

# bridgenorm 0.8.0

* The package is renamed from heb to bridgenorm. Exported functions, classes and schema strings change prefix
  from `heb_` to `bn_` (`heb_fit_node_quantile()` is now `bn_fit_node_quantile()`, and so on). The `heb_settings`
  argument of `bn_fit_nodes()` is now `eb_settings`, and the mixed-model method value `"heb"` is now `"eb"`
  (`"heb_shared"` is now `"eb_shared"`). Models saved by heb do not load; refit them. Estimates are unchanged.

# heb 0.7.2

* `heb_fit_node_quantile()` builds quantile tables only for batches whose reference donors (`reference_id`) are
  exactly the anchor's. Before, every batch got tables, so an unbridged batch that pass 1 left raw was still
  quantile-mapped onto the anchor's donors, and a partly bridged batch pooled donors the anchor lacks. Such a
  batch now takes the pass-2 node shift if bridged and stays raw if not, as the README states.

# heb 0.7.1

* `heb_fit_node_quantile()` sorts each column of a shrunk quantile table (rearrangement; Chernozhukov,
  Fernandez-Val and Galichon 2010, Econometrica 78:1093). With `shrink_k = 300` the mix of node, anchor and
  batch-wide curves left about 0.4% of table columns non-monotone (worst step -0.037 asinh), so a few cells
  swapped order. The fix changes corrected values only slightly. `shrink_k = 0` is unchanged.

# heb 0.7.0

* `heb_fit_node_quantile()` gains a required `shrink_k`: each non-anchor (batch, pass-2 node) offset curve
  (batch minus anchor quantiles) is pulled toward the batch-wide offset curve with weight n / (n + shrink_k),
  n = the batch's reference cells in the node. `shrink_k = 0` reproduces 0.6.0. **Breaking:** calls without `shrink_k`
  now stop. In simulation benchmarks `shrink_k = 300` lowered recovery error and RMSE to truth.

# heb 0.6.0

* New `heb_fit_node_quantile()`: two-pass node correction with per-node quantile maps. Pass 1 is `heb_fit_nodes()`; pass 2 refits the
  map on the shifted reference cells and, in each node where a batch and the anchor both have at least
  `min_cells` reference cells pooled over wells, maps the batch's marker quantiles onto the anchor's
  (linear interpolation, a constant shift beyond the outer quantiles). Other cells take the pass-2 node
  shift. `quantile_probs` and `tail` are required, like every `heb_fit_nodes()` setting.
  In simulation benchmarks it lowered the worst-pair EMD against the node shift.
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
