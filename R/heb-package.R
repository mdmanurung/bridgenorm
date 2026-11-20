#' heb: hierarchical empirical Bayes batch correction from reference samples
#'
#' Fit with \code{\link{heb_fit}} (one model per population and marker) or
#' \code{\link{heb_fit_shared}} (one model per marker, batch effect shared
#' across populations). Apply with \code{\link{heb_integrate}}. Store with
#' \code{\link{heb_save}} and \code{\link{heb_load}}. Check the design with
#' \code{\link{heb_bridges}}. Fit the paired-shift baseline with
#' \code{\link{heb_fit_paired}}, and report reference noise per batch with
#' \code{\link{heb_reference_noise}} before choosing the anchor.
#'
#' Every setting that can change an estimate is a required argument with no
#' default.
#'
#' @keywords internal
#' @importFrom dplyr .data
"_PACKAGE"
