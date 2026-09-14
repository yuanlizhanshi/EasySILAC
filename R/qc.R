#' Quality control for SILAC peptide data
#'
#' @description
#' Filter a SILAC peptide-level SummarizedExperiment by intensity quantile,
#' number of detected time points, ratio range and replicate reproducibility
#' (Spearman), and normalize the total intensity. Returns a list with the
#' filtered assays and QC statistics.
#'
#' @param se A SummarizedExperiment with raw/log2 intensity assays.
#' @param assay_log2 Name of the log2 intensity assay.
#' @param assay_raw Name of the raw intensity assay.
#' @param type_col Column in `colData` marking heavy/light channel.
#' @param time_col Column in `colData` with the labeling time.
#' @param sample_col Column in `colData` with the sample name.
#' @param peptide_col Column in `rowData` with the peptide id.
#' @param heavy_label Label of the heavy channel in `type_col`.
#' @param light_label Label of the light channel in `type_col`.
#' @param intensity_q Intensity quantile cutoff (default 0.01).
#' @param min_detected_timepoints_intensity Minimum number of time points with
#'   detected intensity (default 3).
#' @param min_detected_timepoints_ratio Minimum number of time points with a
#'   valid H/(H+L) ratio (default 3).
#' @param ratio_threshold Minimum ratio value kept (default 0).
#' @param min_spearman_rho Minimum replicate Spearman correlation (default 0).
#' @param norm_total Normalization target for total intensity (default 1e6).
#' @param verbose Print progress messages.
#'
#' @return A list with filtered assays, QC statistics and kept peptides.
#' @export
quality_control_SILAC <- function(
    se,
    assay_log2 = "log2_intensity",
    assay_raw = "raw_intensity",
    type_col = "type",
    time_col = "t",
    sample_col = "sample_name",
    peptide_col = "Peptide",
    heavy_label = "H",
    light_label = "L",
    intensity_q = 0.01,
    min_detected_timepoints_intensity = 3,
    min_detected_timepoints_ratio = 3,
    ratio_threshold = 0,
    min_spearman_rho = 0,
    norm_total = 1e6,
    verbose = TRUE
) {
  msg <- function(..., verbose = TRUE) {
    if (isTRUE(verbose)) {
      message(glue::glue(...))
    }
  }

  
  pct <- function(x, y) {
    if (y == 0) return(NA_real_)
    round(100 * x / y, 1)
  }
  
  stopifnot(inherits(se, "SummarizedExperiment"))
  
  if (!assay_log2 %in% assayNames(se)) {
    stop(glue("assay_log2 '{assay_log2}' not found in assayNames(se)."))
  }
  
  if (!assay_raw %in% assayNames(se)) {
    stop(glue("assay_raw '{assay_raw}' not found in assayNames(se)."))
  }
  
  meta0 <- as.data.frame(colData(se))
  
  need_cols <- c(type_col, time_col)
  missing_cols <- setdiff(need_cols, colnames(meta0))
  if (length(missing_cols) > 0) {
    stop(glue("Missing colData columns: {paste(missing_cols, collapse = ', ')}"))
  }
  
  if (!sample_col %in% colnames(meta0)) {
    meta0[[sample_col]] <- colnames(se)
  }
  
  if (is.null(rownames(se))) {
    stop("rownames(se) must be peptide IDs.")
  }
  
  msg("======================================")
  msg("Starting SILAC ratio processing")
  msg("======================================")
  msg("Input peptides: {nrow(se)}")
  msg("Input samples : {ncol(se)}")
  
  n0 <- nrow(se)
  
  # ----------------------------
  # 1. intensity-level QC
  # ----------------------------
  mtx_log2 <- assay(se, assay_log2)
  
  intensity_log2_value <- tibble::tibble(
    value = as.vector(mtx_log2)
  ) %>%
    dplyr::filter(!is.na(value), value > 0)
  
  q_cutoff <- as.numeric(quantile(
    intensity_log2_value$value,
    probs = intensity_q,
    na.rm = TRUE
  ))
  
  intensity_qc <- count_detected_timepoints(
    mtx_log2,
    time_vec = colData(se)[[time_col]],
    threshold = q_cutoff
  )
  
  intensity_qc_freq <- table(
    intensity_qc$n_timepoints_detected
  ) %>%
    as.data.frame()
  
  valid_pep <- intensity_qc %>%
    dplyr::filter(n_timepoints_detected >= min_detected_timepoints_intensity)
  
  msg("--------------------------------------")
  msg("[1] Intensity-level QC")
  msg("Intensity threshold: q{intensity_q} = {round(q_cutoff, 3)}")
  msg("Detected >= {min_detected_timepoints_intensity} timepoints: {nrow(valid_pep)}/{n0} peptides retained ({pct(nrow(valid_pep), n0)}%)")
  
  if (nrow(valid_pep) == 0) {
    stop("No peptides passed intensity-level QC.")
  }
  
  se_qc <- se[valid_pep[[peptide_col]], ]
  

  # 2. split H / L channels

  se_H <- se_qc[, colData(se_qc)[[type_col]] == heavy_label]
  se_L <- se_qc[, colData(se_qc)[[type_col]] == light_label]
  
  msg("--------------------------------------")
  msg("[2] Channel split")
  msg("H samples: {ncol(se_H)}")
  msg("L samples: {ncol(se_L)}")
  
  if (ncol(se_H) == 0 || ncol(se_L) == 0) {
    stop("H or L channel has zero samples. Check type_col / heavy_label / light_label.")
  }
  
  if (ncol(se_H) != ncol(se_L)) {
    stop(glue("H/L sample number mismatch: H = {ncol(se_H)}, L = {ncol(se_L)}"))
  }
  
  mtx_H <- assay(se_H, assay_raw)
  mtx_L <- assay(se_L, assay_raw)
  
  mtx_H[is.na(mtx_H)] <- 0
  mtx_L[is.na(mtx_L)] <- 0
  

  # 3. normalize total intensity first,
  # then redistribute to H and L

  mtx_total_raw <- mtx_H + mtx_L
  
  total_colsum <- colSums(mtx_total_raw, na.rm = TRUE)
  total_colsum[total_colsum == 0] <- NA_real_
  
  mtx_total_norm <- sweep(
    mtx_total_raw,
    2,
    total_colsum,
    FUN = "/"
  ) * norm_total
  
  mtx_H_prop <- mtx_H / mtx_total_raw
  mtx_L_prop <- mtx_L / mtx_total_raw
  
  mtx_H_prop[!is.finite(mtx_H_prop)] <- 0
  mtx_L_prop[!is.finite(mtx_L_prop)] <- 0
  
  mtx_H_norm <- mtx_total_norm * mtx_H_prop
  mtx_L_norm <- mtx_total_norm * mtx_L_prop
  
  mtx_H_norm[!is.finite(mtx_H_norm)] <- 0
  mtx_L_norm[!is.finite(mtx_L_norm)] <- 0
  
  mtx_total_norm <- mtx_H_norm + mtx_L_norm
  

  # 4. calculate H / (H + L)

  mtx_ratio <- mtx_H / mtx_total_raw
  mtx_ratio[!is.finite(mtx_ratio)] <- NA
  
  ratio_na_fraction <- mean(is.na(mtx_ratio))
  ratio_median <- median(mtx_ratio, na.rm = TRUE)
  
  msg("--------------------------------------")
  msg("[3] Ratio matrix")
  msg("Ratio definition: H / (H + L)")
  msg("Ratio matrix size: {nrow(mtx_ratio)} peptides x {ncol(mtx_ratio)} samples")
  msg("NA fraction: {round(100 * ratio_na_fraction, 2)}%")
  msg("Median ratio: {round(ratio_median, 3)}")
  msg("Total intensity normalized to column sum = {norm_total}")
  

  # 5. build ratio metadata

  sample_names <- colnames(mtx_ratio)
  
  meta_H <- as.data.frame(colData(se_H))
  
  meta_ratio <- meta_H
  meta_ratio$missing_number <- apply(mtx_ratio, 2, function(x) sum(is.na(x)))
  
  if (!sample_col %in% colnames(meta_ratio)) {
    meta_ratio[[sample_col]] <- sample_names
  }
  
  rownames(meta_ratio) <- sample_names
  
  msg("--------------------------------------")
  msg("[4] Ratio colData")
  msg("Ratio samples: {nrow(meta_ratio)}")
  msg("Timepoint distribution:")
  if (isTRUE(verbose)) {
    print(table(meta_ratio[[time_col]]))
  }
  
  # 6. ratio-level QC

  ratio_qc <- count_detected_timepoints(
    mtx_ratio,
    time_vec = colData(se_H)[[time_col]],
    threshold = ratio_threshold
  )
  
  ratio_qc_freq <- table(
    ratio_qc$n_timepoints_detected
  ) %>%
    as.data.frame()
  
  ratio_qc_filter <- ratio_qc %>%
    dplyr::filter(n_timepoints_detected >= min_detected_timepoints_ratio)
  
  msg("--------------------------------------")
  msg("[5] Ratio-level QC")
  msg("Ratio detection threshold: > {ratio_threshold}")
  msg("Detected >= {min_detected_timepoints_ratio} timepoints: {nrow(ratio_qc_filter)}/{nrow(mtx_ratio)} peptides retained ({pct(nrow(ratio_qc_filter), nrow(mtx_ratio))}%)")
  
  if (nrow(ratio_qc_filter) == 0) {
    stop("No peptides passed ratio-level QC.")
  }
  
  mtx_ratio_filter <- mtx_ratio[ratio_qc_filter[[peptide_col]], , drop = FALSE]

  # 7. mean ratio by timepoint

  time_vec <- colData(se_H)[[time_col]]
  time_levels <- unique(time_vec)
  
  mtx_ratio_mean <- sapply(time_levels, function(tt) {
    idx <- which(time_vec == tt)
    rowMeans(mtx_ratio_filter[, idx, drop = FALSE], na.rm = TRUE)
  })
  
  if (is.null(dim(mtx_ratio_mean))) {
    mtx_ratio_mean <- matrix(
      mtx_ratio_mean,
      ncol = 1,
      dimnames = list(rownames(mtx_ratio_filter), time_levels)
    )
  }
  
  colnames(mtx_ratio_mean) <- paste0("ratio_", time_levels)
  
  mtx_ratio_mean <- as.data.frame(mtx_ratio_mean) %>%
    tibble::rownames_to_column(peptide_col)
  
  msg("--------------------------------------")
  msg("[6] Mean ratio by timepoint")
  msg("Timepoints: {paste(time_levels, collapse = ', ')}")
  msg("Peptides in mean ratio matrix: {nrow(mtx_ratio_mean)}")
  

  # 8. correlation QC

  corr_mat <- t(apply(
    mtx_ratio_filter[mtx_ratio_mean[[peptide_col]], , drop = FALSE],
    1,
    calc_corr,
    t = time_vec
  ))
  
  mtx_ratio_mean$spearman_rho <- corr_mat[, "rho"]
  mtx_ratio_mean$p_value <- corr_mat[, "p"]
  
  mtx_ratio_mean_filter <- mtx_ratio_mean %>%
    dplyr::filter(!is.na(spearman_rho)) %>%
    dplyr::filter(spearman_rho > min_spearman_rho)
  
  msg("--------------------------------------")
  msg("[7] Correlation QC")
  msg("Non-NA correlation: {sum(!is.na(mtx_ratio_mean$spearman_rho))}/{nrow(mtx_ratio_mean)} peptides")
  msg("rho > {min_spearman_rho}: {nrow(mtx_ratio_mean_filter)}/{nrow(mtx_ratio_mean)} peptides retained ({pct(nrow(mtx_ratio_mean_filter), nrow(mtx_ratio_mean))}%)")
  
  if (nrow(mtx_ratio_mean_filter) == 0) {
    stop("No peptides passed correlation QC.")
  }
  

  # 9. final peptide set

  final_peptides <- mtx_ratio_mean_filter[[peptide_col]]
  
  mtx_ratio_final <- mtx_ratio_filter[final_peptides, , drop = FALSE]
  mtx_H_final <- mtx_H_norm[final_peptides, , drop = FALSE]
  mtx_L_final <- mtx_L_norm[final_peptides, , drop = FALSE]
  mtx_total_final <- mtx_total_norm[final_peptides, , drop = FALSE]
  
  colnames(mtx_H_final) <- colnames(mtx_ratio_final)
  colnames(mtx_L_final) <- colnames(mtx_ratio_final)
  colnames(mtx_total_final) <- colnames(mtx_ratio_final)
  

  # 10. sample-level ratio median

  sample_ratio_median <- tibble::tibble(
    sample = colnames(mtx_ratio_final),
    ratio_median = matrixStats::colMedians(mtx_ratio_final, na.rm = TRUE),
    missing_number = apply(mtx_ratio_final, 2, function(x) sum(is.na(x)))
  ) %>%
    dplyr::left_join(
      meta_H %>% tibble::rownames_to_column("sample"),
      by = "sample"
    )
  
  msg("--------------------------------------")
  msg("[8] Sample-level QC")
  msg("Sample ratio median range: {round(min(sample_ratio_median$ratio_median, na.rm = TRUE), 3)} - {round(max(sample_ratio_median$ratio_median, na.rm = TRUE), 3)}")
  

  # 11. construct FINAL ratio SE

  se_ratio <- SummarizedExperiment(
    assays = list(
      ratio = mtx_ratio_final,
      H_intensity = mtx_H_final,
      L_intensity = mtx_L_final,
      total_intensity = mtx_total_final
    ),
    rowData = rowData(se_qc[final_peptides, ]),
    colData = S4Vectors::DataFrame(meta_ratio)
  )
  
  msg("======================================")
  msg("Final Summary")
  msg("======================================")
  msg("Initial peptides: {n0}")
  msg("After intensity QC: {nrow(se_qc)}")
  msg("After ratio QC: {nrow(mtx_ratio_filter)}")
  msg("After correlation QC / final se_ratio: {nrow(se_ratio)}")
  msg("Overall final retention: {pct(nrow(se_ratio), n0)}%")
  msg("======================================")
  
  list(
    se_ratio = se_ratio,
    
    mtx_ratio = mtx_ratio,
    mtx_ratio_filter = mtx_ratio_filter,
    mtx_ratio_final = mtx_ratio_final,
    
    mtx_ratio_mean = mtx_ratio_mean,
    mtx_ratio_mean_filter = mtx_ratio_mean_filter,
    
    mtx_H_norm = mtx_H_norm,
    mtx_L_norm = mtx_L_norm,
    mtx_total_norm = mtx_total_norm,
    
    mtx_H_final = mtx_H_final,
    mtx_L_final = mtx_L_final,
    mtx_total_final = mtx_total_final,
    
    intensity_qc = intensity_qc,
    intensity_qc_freq = intensity_qc_freq,
    ratio_qc = ratio_qc,
    ratio_qc_freq = ratio_qc_freq,
    ratio_qc_filter = ratio_qc_filter,
    
    sample_ratio_median = sample_ratio_median,
    meta_ratio = meta_ratio,
    
    q_cutoff = q_cutoff,
    valid_pep = valid_pep,
    final_peptides = final_peptides,
    norm_total = norm_total
  )
}

#' Count detected time points
#'
#' @param mtx Numeric matrix (features x samples).
#' @param time_vec Time point of each sample.
#' @param threshold Detection threshold (values above count as detected).
#'
#' @return Integer vector of per-feature counts of detected time points.
#' @export
count_detected_timepoints <- function(mtx, time_vec, threshold) {
  stopifnot(
    is.matrix(mtx),
    length(time_vec) == ncol(mtx)
  )

  detected <- mtx > threshold

  time_levels <- unique(time_vec)

  time_detect_mat <- sapply(time_levels, function(tp) {
    cols <- which(time_vec == tp)
    matrixStats::rowAnys(detected[, cols, drop = FALSE], na.rm = TRUE)
  })

  n_timepoints <- rowSums(time_detect_mat)

  tibble::tibble(
    Peptide = rownames(mtx),
    n_timepoints_detected = n_timepoints
  )
}

#' Turnover weighted by a quantitative assay
#'
#' @param protein_fit_se SummarizedExperiment with fit results in `rowData`.
#' @param assay Assay used for weighting (default "ratio").
#'
#' @return A tibble of weighted turnover statistics.
#' @export
get_weighted_turnover <-  function(protein_fit_se,assay = 'ratio'){
  total_intensity <- assay(protein_fit_se,'intensity')
  H_ratio <- assay(protein_fit_se,assay)
  total_sample_ratio <- tibble(
    sample = colnames(protein_fit_se),
    ratio = colSums(total_intensity * H_ratio,na.rm = T)/colSums(total_intensity,na.rm = T),
    t = as.numeric(protein_fit_se$t) 
  )
  fit_res <- fit_silac(ratio = total_sample_ratio$ratio,t = total_sample_ratio$t,Peptide = 'Total (weighted mean)')
  return(fit_res)
}

#' Correlation between a vector and time
#'
#' @param x Numeric vector.
#' @param t Time vector.
#'
#' @return Numeric correlation.
#' @export
calc_corr <- function(x, t) {
  keep <- !is.na(x) & !is.na(t)
  if (sum(keep) < 3) return(c(rho = NA, p = NA))

  ct <- suppressWarnings(cor.test(x[keep], t[keep], method = "spearman"))

  c(rho = unname(ct$estimate), p = ct$p.value)
}

#' Extract the new-to-total ratio of one peptide
#'
#' @param se SummarizedExperiment.
#' @param peptide Peptide row name.
#'
#' @return A tibble with time and ratio.
#' @export
extract_ntr_from_se <- function(se, peptide) {
  tibble(
    t = colData(se)$t,
    ntr = as.numeric(assay(se, "ratio")[peptide, ])
  )
}

#' Precision-weighted degradation rate
#'
#' @param k Numeric vector of per-peptide rates.
#' @param stderr Standard errors of `k`.
#' @param min_peptides Minimum number of peptides required (default 3).
#'
#' @return Weighted mean of `k`, or NA if fewer than `min_peptides`.
#' @export
weighted_k <- function(k, stderr, min_peptides = 3) {
  ok <- !is.na(k) & !is.na(stderr) & stderr > 0

  k <- k[ok]
  stderr <- stderr[ok]

  if (length(k) == 0) return(NA_real_)

  if (length(k) < min_peptides) {
    return(mean(k, na.rm = TRUE))
  }

  w <- 1 / stderr
  sum(w * k) / sum(w)
}
