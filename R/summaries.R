# Reference summaries and paired effects for new batches.
#
# Ported from SpectraWeaveR R/corrector_paired_shift.R (.sw_ps_summaries,
# .sw_ps_donor_level, .sw_ps_effects). Cells are measurements, not replicates:
# the replicate unit is the reference aliquot.

# Columns every effects table must have.
.bn_effect_cols <- c(
  "batch",
  "population",
  "marker",
  "estimate",
  "se_donor",
  "n_donors",
  "n_aliquots",
  "n_aliquots_anchor",
  "support",
  "source"
)

.bn_empty_effects <- function() {
  tibble::tibble(
    batch = character(),
    population = character(),
    marker = character(),
    estimate = double(),
    se_donor = double(),
    n_donors = integer(),
    n_aliquots = integer(),
    n_aliquots_anchor = integer(),
    support = character(),
    source = character()
  )
}

.bn_empty_excluded <- function() {
  tibble::tibble(
    reference_id = character(),
    aliquot = character(),
    batch = character(),
    population = character(),
    n_cells = integer(),
    reason = character()
  )
}

# Aliquot summaries in long form, plus the summaries below min_cells.
# value = mean of the cells; var_mean = var(cells) / n, the known sampling variance.
.bn_summaries <- function(xm, cells, min_cells) {
  markers <- colnames(xm)
  ref <- cells[cells$is_reference, , drop = FALSE]
  empty <- tibble::tibble(
    reference_id = character(),
    aliquot = character(),
    batch = character(),
    population = character(),
    marker = character(),
    value = double(),
    n_cells = integer(),
    var_mean = double()
  )
  if (nrow(ref) == 0) {
    return(list(summ = empty, excluded = .bn_empty_excluded()))
  }
  key <- paste(
    ref$reference_id,
    ref$aliquot,
    ref$batch,
    ref$population,
    sep = "\r"
  )
  groups <- split(ref$row, key)
  info <- ref[
    match(names(groups), key),
    c("reference_id", "aliquot", "batch", "population")
  ]
  n <- lengths(groups)
  keep <- n >= min_cells

  excluded <- tibble::tibble(
    info[!keep, , drop = FALSE],
    n_cells = as.integer(n[!keep]),
    reason = sprintf(
      "%d cell(s), below min_cells = %d; excluded from the effect estimate.",
      n[!keep],
      min_cells
    )
  )
  if (!any(keep)) {
    return(list(summ = empty, excluded = excluded))
  }

  # Mean of sorted values, so a summary does not depend on row order.
  means <- do.call(
    rbind,
    lapply(groups[keep], function(ii) {
      apply(xm[ii, , drop = FALSE], 2, function(v) mean(sort(v)))
    })
  )
  var_means <- do.call(
    rbind,
    lapply(groups[keep], function(ii) {
      apply(xm[ii, , drop = FALSE], 2, stats::var) / length(ii)
    })
  )
  k <- which(keep)
  n_m <- length(markers)
  summ <- tibble::tibble(
    reference_id = rep(info$reference_id[k], each = n_m),
    aliquot = rep(info$aliquot[k], each = n_m),
    batch = rep(info$batch[k], each = n_m),
    population = rep(info$population[k], each = n_m),
    marker = rep(markers, length(k)),
    value = as.vector(t(means)),
    n_cells = rep(as.integer(n[k]), each = n_m),
    var_mean = as.vector(t(var_means))
  )
  list(summ = summ, excluded = excluded)
}

# Donor summaries: equal-weight mean of aliquot summaries per reference and batch.
.bn_donor_level <- function(summ) {
  key <- paste(
    summ$reference_id,
    summ$batch,
    summ$population,
    summ$marker,
    sep = "\r"
  )
  g <- unname(split(seq_len(nrow(summ)), key))
  first <- vapply(g, `[`, integer(1), 1)
  tibble::tibble(
    summ[first, c("reference_id", "batch", "population", "marker")],
    value = vapply(g, function(ii) mean(sort(summ$value[ii])), numeric(1)),
    n_aliquots = as.integer(lengths(g))
  )
}

# Anchor targets: the donor summaries of the anchor batch.
.bn_targets <- function(donor, anchor) {
  anc <- donor[donor$batch == anchor, , drop = FALSE]
  if (nrow(anc) == 0) {
    stop(
      sprintf("Anchor batch '%s' has no reference summary.", anchor),
      call. = FALSE
    )
  }
  tibble::tibble(
    reference_id = anc$reference_id,
    population = anc$population,
    marker = anc$marker,
    target = anc$value,
    n_aliquots_anchor = anc$n_aliquots
  )
}

# Paired effects for batches unseen at fit time: mean over shared references of
# (batch summary - frozen anchor summary), no shrinkage.
.bn_paired_effects <- function(donor, targets, source) {
  paired <- dplyr::inner_join(
    donor,
    targets,
    by = c("reference_id", "population", "marker"),
    relationship = "many-to-one"
  )
  if (nrow(paired) == 0) {
    return(.bn_empty_effects())
  }
  diff <- paired$value - paired$target
  key <- paste(paired$batch, paired$population, paired$marker, sep = "\r")
  g <- unname(split(seq_len(nrow(paired)), key))
  first <- vapply(g, `[`, integer(1), 1)
  tibble::tibble(
    batch = paired$batch[first],
    population = paired$population[first],
    marker = paired$marker[first],
    estimate = vapply(g, function(ii) mean(sort(diff[ii])), numeric(1)),
    se_donor = vapply(
      g,
      function(ii) {
        if (length(ii) < 2) {
          NA_real_
        } else {
          stats::sd(sort(diff[ii])) / sqrt(length(ii))
        }
      },
      numeric(1)
    ),
    df = ifelse(lengths(g) >= 2, lengths(g) - 1, NA_real_),
    n_donors = as.integer(lengths(g)),
    n_aliquots = vapply(g, function(ii) sum(paired$n_aliquots[ii]), integer(1)),
    n_aliquots_anchor = vapply(
      g,
      function(ii) sum(paired$n_aliquots_anchor[ii]),
      integer(1)
    ),
    support = "paired_reference",
    source = source
  )
}
