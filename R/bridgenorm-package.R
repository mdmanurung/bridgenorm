#' bridgenorm: cytometry batch correction from bridge samples
#'
#' Fit the main two-pass node and quantile model with
#' \code{\link{bn_fit_node_quantile}}, or a node-shift model with
#' \code{\link{bn_fit_nodes}}. Baselines: \code{\link{bn_fit_paired}} (paired
#' shift), \code{\link{bn_fit}} (mixed model per population and marker) and
#' \code{\link{bn_fit_shared}} (batch effect shared across populations).
#' Apply any model with \code{\link{bn_integrate}}; store it with
#' \code{\link{bn_save}} and \code{\link{bn_load}}. Check the design with
#' \code{\link{bn_bridges}}, and report bridge-sample noise per batch with
#' \code{\link{bn_reference_noise}} before choosing the anchor.
#'
#' Every setting that can change an estimate is a required argument with no
#' default.
#'
#' @keywords internal
#' @importFrom dplyr .data
"_PACKAGE"
