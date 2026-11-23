#' bridgenorm: cytometry batch correction from bridge samples
#'
#' Fit the two-pass node and quantile model with \code{\link{bn_fit}}, or the
#' paired-shift baseline with \code{\link{bn_fit_paired}}. Correct data with
#' \code{\link{bn_apply}}; store a model with \code{\link{bn_save}} and
#' \code{\link{bn_load}}. Check the design with \code{\link{bn_bridges}}, and
#' report bridge-sample noise per batch with \code{\link{bn_reference_noise}}
#' before choosing the anchor.
#'
#' Settings that depend on the data (anchor, SOM grid, mapping markers, seed,
#' transform) have no default. The other settings default to the values
#' benchmarked on the data set A and data set B data.
#'
#' @keywords internal
#' @importFrom dplyr .data
"_PACKAGE"
