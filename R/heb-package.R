#' heb: hierarchical empirical Bayes batch correction from reference samples
#'
#' Fit with \code{\link{heb_fit}} (one model per population and marker) or
#' \code{\link{heb_fit_shared}} (one model per marker, batch effect shared
#' across populations). Apply with \code{\link{heb_integrate}}. Store with
#' \code{\link{heb_save}} and \code{\link{heb_load}}. Check the design with
#' \code{\link{heb_bridges}}.
#'
#' Every setting that can change an estimate is a required argument with no
#' default.
#'
#' @keywords internal
"_PACKAGE"
