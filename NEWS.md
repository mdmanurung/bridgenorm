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
