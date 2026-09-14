# Rust-accelerated fitting engine (rextendr, compiled on first use).
#
# The Rust implementations in inst/rust/silacfit.rs are numerically identical
# to the R versions (closed-form OLS replaces lm(); the bounded least-squares
# closed form replaces nls(port), which is exact because the model is linear
# in the synthesis rate given k and P0).

#' Load the compiled Rust fitting library
#'
#' @description Compiles (first call in a session, cached afterwards by
#'   rextendr) and loads the Rust fitting functions. Called automatically by
#'   the `*_rust` / `*_batch` functions; requires a Rust toolchain and the
#'   rextendr package.
#'
#' @return Invisibly TRUE.
#' @export
silacfit_load_rust <- function() {
  if (!requireNamespace("rextendr", quietly = TRUE)) {
    stop("Package 'rextendr' is required for the Rust fitting engine. ",
         "Install it with install.packages('rextendr').")
  }
  if (!isTRUE(.GlobalEnv$.EasySILAC_rust_loaded)) {
    rextendr::rust_source(
      system.file("rust", "silacfit.rs", package = "EasySILAC", mustWork = TRUE),
      cache_build = TRUE,
      quiet = TRUE,
      env = .GlobalEnv
    )
    .GlobalEnv$.EasySILAC_rust_loaded <- TRUE
  }
  invisible(TRUE)
}

#' Fit SILAC heavy-ratio kinetics, Rust engine
#'
#' @description Rust implementation of [fit_silac()] with identical numerical
#'   output (verified to ~1e-14 on production data).
#'
#' @inheritParams fit_silac
#' @return A one-row tibble, or NULL if not fittable.
#' @export
fit_silac_rust <- function(ratio, t, Peptide, total = NULL) {
  silacfit_load_rust()
  res <- fit_silac_rs(
    as.numeric(ratio), as.numeric(t),
    if (is.null(total)) NULL else as.numeric(total),
    as.character(Peptide)
  )
  if (is.null(res)) return(NULL)
  tibble::as_tibble(res)
}

#' Fit non-steady-state SILAC kinetics, Rust engine
#'
#' @description Rust implementation of [fit_silac_nonsteady()] with identical
#'   numerical output.
#'
#' @inheritParams fit_silac_nonsteady
#' @return A one-row tibble, or NULL if not fittable.
#' @export
fit_silac_nonsteady_rust <- function(ratio, total, t, Peptide,
                                     min_points = 3, min_k = 0.0001) {
  silacfit_load_rust()
  res <- fit_silac_nonsteady_rs(
    as.numeric(ratio), as.numeric(total), as.numeric(t),
    as.character(Peptide)
  )
  if (is.null(res)) return(NULL)
  tibble::as_tibble(res)
}

#' Batch steady-state fit over matrix rows (multithreaded Rust)
#'
#' @description Fits every row of `mtx_ratio` in one call, multithreaded
#'   inside Rust. Rows that fail the fit are dropped (same as
#'   `purrr::map_dfr` over [fit_silac()] skipping NULLs).
#'
#' @param mtx_ratio Numeric matrix (features x samples) of heavy ratios.
#' @param t Numeric vector of labeling times (length `ncol(mtx_ratio)`).
#' @param mtx_total Numeric matrix of total intensities, same dimensions.
#'
#' @return A data.frame with one row per successfully fitted feature;
#'   columns match [fit_silac()].
#' @export
fit_silac_batch <- function(mtx_ratio, t, mtx_total) {
  silacfit_load_rust()
  fit_silac_batch_rs(
    as.numeric(mtx_ratio), nrow(mtx_ratio), ncol(mtx_ratio),
    as.numeric(t), as.numeric(mtx_total), rownames(mtx_ratio)
  )
}

#' Batch non-steady-state fit over matrix rows (multithreaded Rust)
#'
#' @inheritParams fit_silac_batch
#' @return A data.frame; columns match [fit_silac_nonsteady()].
#' @export
fit_silac_nonsteady_batch <- function(mtx_ratio, mtx_total, t) {
  silacfit_load_rust()
  fit_silac_nonsteady_batch_rs(
    as.numeric(mtx_ratio), nrow(mtx_ratio), ncol(mtx_ratio),
    as.numeric(mtx_total), as.numeric(t), rownames(mtx_ratio)
  )
}

#' Enable or disable the Rust fitting engine globally
#'
#' @description When enabled, [fit_silac()] and [fit_silac_nonsteady()]
#'   dispatch to the Rust implementations ([fit_silac_rust()] /
#'   [fit_silac_nonsteady_rust()]). Disabled by default (pure R).
#'
#' @param enable Logical; default TRUE.
#'
#' @return Invisibly the new setting.
#' @export
silac_use_rust <- function(enable = TRUE) {
  options(EasySILAC.use_rust = isTRUE(enable))
  invisible(isTRUE(enable))
}

#' Fit SILAC heavy-ratio kinetics (steady state)
#'
#' @description
#' Fits `log(1 - ratio) ~ t` by linear regression on all replicate points:
#' the heavy fraction follows `ratio(t) = 1 - C * exp(-k * t)`. Returns the
#' degradation rate `k`, half-life, intercept, standard errors and (when
#' `total` is given) the protein abundance `P0` and synthesis rate `k * P0`.
#'
#' By default the pure-R implementation is used; after `silac_use_rust(TRUE)`
#' the numerically identical Rust engine is used (see [fit_silac_rust()]).
#'
#' @param ratio Numeric vector of observed heavy ratios H/(H+L).
#' @param t Numeric vector of labeling times.
#' @param Peptide Peptide/protein identifier (carried into the output).
#' @param total Optional numeric vector of total intensities (H+L).
#'
#' @return A one-row tibble with fit results, or NULL if the fit is not
#'   possible (too few points or non-finite rate).
#' @export
fit_silac <- function(ratio, t, Peptide, total = NULL) {
  if (isTRUE(getOption("EasySILAC.use_rust", FALSE))) {
    return(fit_silac_rust(ratio, t, Peptide, total))
  }
  .fit_silac_R(ratio = ratio, t = t, Peptide = Peptide, total = total)
}

#' Fit SILAC kinetics without the steady-state assumption
#'
#' @description
#' Non-steady-state fit: the degradation rate `k` comes from
#' `lm(log(old) ~ t)` with `old = total * (1 - ratio)`, and the synthesis rate
#' `s` from `nls(total ~ P0*exp(-k*(t-t0)) + (s/k)*(1-exp(-k*(t-t0))))`
#' (port algorithm, `s >= 0`). Reports `r2_old_decay` and total-fit error
#' metrics used for QC filtering.
#'
#' By default the pure-R implementation is used; after `silac_use_rust(TRUE)`
#' the numerically identical Rust engine is used (see
#' [fit_silac_nonsteady_rust()]).
#'
#' @param ratio Numeric vector of observed heavy ratios.
#' @param total Numeric vector of total intensities.
#' @param t Numeric vector of labeling times.
#' @param Peptide Peptide/protein identifier.
#' @param min_points Minimum number of distinct time points (default 3).
#' @param min_k Lower bound applied to the degradation rate (default 1e-4).
#'
#' @return A one-row tibble with fit results, or NULL if not fittable.
#' @export
fit_silac_nonsteady <- function(ratio, total, t, Peptide,
                                min_points = 3, min_k = 0.0001) {
  if (isTRUE(getOption("EasySILAC.use_rust", FALSE))) {
    return(fit_silac_nonsteady_rust(ratio, total, t, Peptide, min_points, min_k))
  }
  .fit_silac_nonsteady_R(ratio = ratio, total = total, t = t,
                         Peptide = Peptide, min_points = min_points, min_k = min_k)
}

# ---------------------------------------------------------------------------
# Helpers for protein-level fitting from peptide matrices
# ---------------------------------------------------------------------------

#' Build a per-gene ratio matrix from peptide-level ratios
#'
#' @description Aggregates peptides to genes with `matrixStats::colMedians`
#'   per sample (single-peptide genes pass through unchanged), matching the
#'   SILAC pipeline protein-level aggregation.
#'
#' @param mtx_ratio Peptide x sample ratio matrix.
#' @param gene_info Data frame with columns `peptide` and `gene`.
#' @param genes Gene order of the output rows (default: unique genes).
#'
#' @return A gene x sample numeric matrix.
#' @export
build_gene_ratio_matrix <- function(mtx_ratio, gene_info, genes = unique(gene_info$gene)) {
  pep_by_gene <- split(gene_info$peptide, gene_info$gene)[genes]
  t(vapply(
    pep_by_gene,
    function(peps) {
      sub <- mtx_ratio[peps, , drop = FALSE]
      if (nrow(sub) == 1) as.numeric(sub[1, ]) else matrixStats::colMedians(sub, na.rm = TRUE)
    },
    numeric(ncol(mtx_ratio))
  ))
}

#' Build the per-gene fit_info table from a gene ratio matrix
#'
#' @param gene_ratio_mtx Gene x sample ratio matrix.
#' @param genes Gene order (rows).
#' @param sample_names Sample names (columns).
#' @param t_vec Labeling times (columns).
#'
#' @return A long tibble with gene, sample_name, t, ratio, adjusted_ratio.
#' @export
build_fit_info <- function(gene_ratio_mtx, genes, sample_names, t_vec) {
  n_s <- ncol(gene_ratio_mtx)
  tibble::tibble(
    gene = rep(genes, each = n_s),
    sample_name = rep(sample_names, times = length(genes)),
    t = rep(t_vec, times = length(genes)),
    ratio = as.numeric(t(gene_ratio_mtx)),
    adjusted_ratio = as.numeric(t(gene_ratio_mtx))
  )
}

#' Batch steady-state protein fit with per-gene list output
#'
#' @description Reproduces the SILAC pipeline's per-gene steady-state fitting
#'   loop as a single batch Rust call, keeping the output structure (an
#'   unnamed list of per-gene `list(fit_results, fit_info)`).
#'
#' @inheritParams build_gene_ratio_matrix
#' @param t_vec Labeling times.
#' @param mtx_total Gene x sample total-intensity matrix.
#' @param sample_names Sample names for the fit_info table.
#'
#' @return An unnamed list, one element per gene.
#' @export
run_steady_fit_block <- function(mtx_ratio, t_vec, mtx_total, gene_info, sample_names) {
  genes <- unique(gene_info$gene)
  gene_ratio_mtx <- build_gene_ratio_matrix(mtx_ratio, gene_info, genes)
  fit_df <- as.data.frame(fit_silac_batch(gene_ratio_mtx, t_vec, mtx_total[genes, ]))
  fit_split <- split(fit_df, fit_df$Peptide)
  fit_info <- build_fit_info(gene_ratio_mtx, genes, sample_names, t_vec)
  fit_info_split <- split(fit_info, fit_info$gene)
  lapply(genes, function(g) {
    list(
      fit_results = if (g %in% names(fit_split)) tibble::as_tibble(fit_split[[g]]) else NULL,
      fit_info = fit_info_split[[g]]
    )
  })
}
