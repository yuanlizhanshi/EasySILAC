#' EasySILAC: SILAC protein turnover analysis with R and Rust fitting engines
#'
#' @description
#' Fitting (steady-state and non-steady-state), labeling-efficiency (Hill)
#' correction, QC and visualization for pulsed SILAC proteomics data.
#'
#' @keywords internal
#' @import dplyr
#' @import ggplot2
#' @importFrom tibble tibble
#' @importFrom purrr map map_dfr map_dbl
#' @importFrom glue glue
#' @importFrom stats lm coef nls nls.control predict median cor.test setNames
#' @importFrom cowplot plot_grid
#' @importFrom stringr str_extract
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom grDevices cairo_pdf dev.off png tiff
#' @importFrom grid grid.draw
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom utils globalVariables
"_PACKAGE"

utils::globalVariables(c(
  "Dim.1", "Dim.2", "filename", "gene", "intercept", "no_remove", "ntr",
  "old", "r2", "removed_t", "total", "type", "value", ".",
  "E_fit", "E_x", "Peptide", "c_fit", "c_obs", "conc", "max_ratio",
  "mean_ratio", "n_timepoints_detected", "new_ratio", "old_ratio",
  "point_number_old", "point_number", "point_number_total", "r2_for_rank", "r2_old_decay", "P0", "L0", "ratio_fit", "spearman_rho",
  # Rust entry points loaded at runtime by silacfit_load_rust()
  "fit_silac_rs", "fit_silac_nonsteady_rs",
  "fit_silac_batch_rs", "fit_silac_nonsteady_batch_rs"
))
