#' Fit SILAC heavy-ratio kinetics (steady state)
#'
#' @description
#' Fits `log(1 - ratio) ~ t` by linear regression on all replicate points:
#' the heavy fraction follows `ratio(t) = 1 - C * exp(-k * t)`. Returns the
#' degradation rate `k`, half-life, intercept, standard errors and (when
#' `total` is given) the protein abundance `P0` and synthesis rate `k * P0`.
#'
#' A Rust implementation with identical output is available via
#' [fit_silac_rust()] / [fit_silac_batch()], or globally with
#' `silac_use_rust(TRUE)`.
#'
#' @param ratio Numeric vector of observed heavy ratios H/(H+L).
#' @param t Numeric vector of labeling times.
#' @param Peptide Peptide/protein identifier (carried into the output).
#' @param total Optional numeric vector of total intensities (H+L).
#'
#' @return A one-row tibble with fit results, or NULL if the fit is not
#'   possible (too few points or non-finite rate).
#' @keywords internal
#' @noRd
.fit_silac_R <- function(ratio, t, Peptide, total = NULL) {
  stopifnot(length(ratio) == length(t))
  if (!is.null(total)) stopifnot(length(total) == length(t))
  
  keep <- is.finite(ratio) & 
    is.finite(t) &
    ratio > 0 &
    ratio < 1
  
  if (!is.null(total)) {
    keep <- keep & is.finite(total)
  }
  
  ratio <- ratio[keep]
  t <- t[keep]
  if (!is.null(total)) total <- total[keep]
  
  if (length(ratio) < 3) return(NULL)
  if (length(unique(t)) < 3) return(NULL)
  
  
  eps <- 1e-6
  ratio <- pmin(pmax(ratio, eps), 1 - eps)
  
  y <- log(1 - ratio)
  
  fit <- lm(y ~ t)
  b <- coef(fit)[1]
  slope <- coef(fit)[2]
  k <- -slope
  
  s <- summary(fit)
  
  
  P0 <- NA_real_
  synthesis_rate <- NA_real_
  
  if (!is.null(total)) {
    t0 <- min(t, na.rm = TRUE)
    idx0 <- which(t == t0)
    
    if (length(idx0) > 0) {
      P0 <- mean(total[idx0], na.rm = TRUE)
      synthesis_rate <- k * P0
    }
  }
  
  tibble::tibble(
    Peptide = Peptide,
    point_number = length(ratio),
    intercept = b,
    k = k,
    half_life = log(2) / k,
    r2 = s$r.squared,
    se_intercept = s$coefficients["(Intercept)", 2],
    se_k = s$coefficients["t", 2],
    P0 = P0,
    synthesis_rate = synthesis_rate,
    H0 = 1 - exp(intercept)
  )
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
#' A Rust implementation with identical output is available via
#' [fit_silac_nonsteady_rust()] / [fit_silac_nonsteady_batch()], or globally
#' with `silac_use_rust(TRUE)`.
#'
#' @param ratio Numeric vector of observed heavy ratios.
#' @param total Numeric vector of total intensities.
#' @param t Numeric vector of labeling times.
#' @param Peptide Peptide/protein identifier.
#' @param min_points Minimum number of distinct time points (default 3).
#' @param min_k Lower bound applied to the degradation rate (default 1e-4).
#'
#' @return A one-row tibble with fit results, or NULL if not fittable.
#' @keywords internal
#' @noRd
.fit_silac_nonsteady_R <- function(
    ratio,
    total,
    t,
    Peptide,
    min_points = 3,
    min_k = 0.0001
) {
  stopifnot(length(ratio) == length(t))
  stopifnot(length(total) == length(t))
  
  df_total <- tibble::tibble(
    ratio = as.numeric(ratio),
    total = as.numeric(total),
    t = as.numeric(t)
  ) %>%
    dplyr::filter(
      is.finite(ratio),
      is.finite(total),
      is.finite(t),
      total > 0
    ) %>%
    dplyr::mutate(
      old = total * (1 - ratio)
    ) %>%
    dplyr::arrange(t)
  
  df_old <- df_total %>%
    dplyr::filter(
      is.finite(old),
      old > 0
    )
  
  if (length(unique(df_old$t)) < min_points) return(NULL)
  
  fit_k <- tryCatch(
    lm(log(old) ~ t, data = df_old),
    error = function(e) NULL
  )
  
  if (is.null(fit_k)) return(NULL)
  
  k_raw <- -coef(fit_k)[["t"]]
  intercept_k <- coef(fit_k)[["(Intercept)"]]
  L0_fit <- exp(intercept_k)
  
  if (!is.finite(k_raw)) return(NULL)
  
  # 给降解速率设置最低值
  k <- max(k_raw, min_k)
  
  t0 <- min(df_total$t, na.rm = TRUE)
  P0 <- mean(df_total$total[df_total$t == t0], na.rm = TRUE)
  
  if (!is.finite(P0) || P0 <= 0) return(NULL)
  
  df_fit_s <- df_total %>%
    dplyr::filter(t > t0)
  
  if (nrow(df_fit_s) < 1) return(NULL)
  
  fit_s <- tryCatch(
    nls(
      total ~ P0 * exp(-k * (t - t0)) +
        (s / k) * (1 - exp(-k * (t - t0))),
      data = df_fit_s,
      start = list(s = k * P0),
      algorithm = "port",
      lower = c(s = 0),
      control = nls.control(maxiter = 200, warnOnly = TRUE)
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit_s)) {
    s_fit <- NA_real_
    r2_total_fit <- NA_real_
    rmse_total_fit <- NA_real_
    nrmse_total_fit <- NA_real_
  } else {
    s_fit <- coef(fit_s)[["s"]]
    
    pred_total <- predict(fit_s)
    obs_total <- df_fit_s$total
    
    ss_res <- sum((obs_total - pred_total)^2, na.rm = TRUE)
    ss_tot <- sum((obs_total - mean(obs_total, na.rm = TRUE))^2, na.rm = TRUE)
    
    r2_total_fit <- ifelse(ss_tot > 0, 1 - ss_res / ss_tot, NA_real_)
    
    rmse_total_fit <- sqrt(mean((obs_total - pred_total)^2, na.rm = TRUE))
    nrmse_total_fit <- rmse_total_fit / mean(obs_total, na.rm = TRUE)
  }
  
  tibble::tibble(
    Peptide = Peptide,
    point_number_total = nrow(df_total),
    point_number_old = nrow(df_old),
    P0 = P0,
    L0_fit = L0_fit,
    k_raw = k_raw,
    k = k,
    k_is_min_capped = k_raw < min_k,
    half_life = log(2) / k,
    synthesis_fit = s_fit,
    synthesis_norm = s_fit / P0,
    balance = s_fit / (k * P0),
    r2_old_decay = summary(fit_k)$r.squared,
    r2_total_fit = r2_total_fit,
    rmse_total_fit = rmse_total_fit,
    nrmse_total_fit = nrmse_total_fit
  )
}

#' Steady-state SILAC fit on per-time-point geometric means
#'
#' @description
#' Like [fit_silac()], but replicates are first averaged per time point with
#' the geometric mean before the linear fit. Note this inflates R-squared
#' relative to fitting all replicate points.
#'
#' @inheritParams fit_silac
#'
#' @return A one-row tibble with fit results, or NULL.
#' @export
fit_silac_mean <- function(ratio, t, Peptide, total = NULL) {
  stopifnot(length(ratio) == length(t))
  if (!is.null(total)) stopifnot(length(total) == length(t))
  
  df <- tibble::tibble(
    ratio = as.numeric(ratio),
    t = as.numeric(t),
    total = if (is.null(total)) NA_real_ else as.numeric(total)
  ) %>%
    dplyr::filter(
      is.finite(ratio),
      is.finite(t),
      ratio > 0,
      ratio < 1
    )
  
  if (!is.null(total)) {
    df <- df %>%
      dplyr::filter(is.finite(total), total > 0)
  }
  
  if (nrow(df) < 3) return(NULL)
  if (length(unique(df$t)) < 3) return(NULL)
  
  df_mean <- df %>%
    dplyr::group_by(t) %>%
    dplyr::summarise(
      ratio = geom_mean(ratio),
      total = if (all(is.na(total))) NA_real_ else geom_mean(total),
      .groups = "drop"
    ) %>%
    dplyr::filter(
      is.finite(ratio),
      ratio > 0,
      ratio < 1,
      is.finite(t)
    ) %>%
    dplyr::arrange(t)
  
  if (nrow(df_mean) < 3) return(NULL)
  if (length(unique(df_mean$t)) < 3) return(NULL)
  
  eps <- 1e-6
  df_mean <- df_mean %>%
    dplyr::mutate(
      ratio = pmin(pmax(ratio, eps), 1 - eps),
      y = log(1 - ratio)
    )
  
  fit <- lm(y ~ t, data = df_mean)
  
  b <- coef(fit)[["(Intercept)"]]
  slope <- coef(fit)[["t"]]
  k <- -slope
  
  if (!is.finite(k)) return(NULL)
  
  s <- summary(fit)
  
  P0 <- NA_real_
  synthesis_rate <- NA_real_
  
  if (!is.null(total)) {
    t0 <- min(df_mean$t, na.rm = TRUE)
    idx0 <- which(df_mean$t == t0)
    
    if (length(idx0) > 0) {
      P0 <- mean(df_mean$total[idx0], na.rm = TRUE)
      synthesis_rate <- k * P0
    }
  }
  
  tibble::tibble(
    Peptide = Peptide,
    point_number = nrow(df_mean),
    intercept = b,
    k = k,
    half_life = log(2) / k,
    r2 = s$r.squared,
    se_intercept = s$coefficients["(Intercept)", 2],
    se_k = s$coefficients["t", 2],
    P0 = P0,
    synthesis_rate = synthesis_rate,
    H0 = 1 - exp(b)
  )
}

#' Non-steady-state SILAC fit on per-time-point geometric means
#'
#' @description Like [fit_silac_nonsteady()], but replicates are first
#'   averaged per time point with the geometric mean.
#'
#' @inheritParams fit_silac_nonsteady
#'
#' @return A one-row tibble with fit results, or NULL.
#' @export
fit_silac_nonsteady_mean <- function(
    ratio,
    total,
    t,
    Peptide,
    min_points = 3,
    min_k = 0.0001
) {
  stopifnot(length(ratio) == length(t))
  stopifnot(length(total) == length(t))
  
  df_raw <- tibble::tibble(
    ratio = as.numeric(ratio),
    total = as.numeric(total),
    t = as.numeric(t)
  ) %>%
    dplyr::filter(
      is.finite(ratio),
      is.finite(total),
      is.finite(t),
      ratio > 0,
      ratio < 1,
      total > 0
    ) %>%
    dplyr::mutate(
      old = total * (1 - ratio)
    ) %>%
    dplyr::filter(
      is.finite(old),
      old > 0
    )
  
  if (nrow(df_raw) < min_points) return(NULL)
  if (length(unique(df_raw$t)) < min_points) return(NULL)
  
  df_total <- df_raw %>%
    dplyr::group_by(t) %>%
    dplyr::summarise(
      ratio = geom_mean(ratio),
      total = geom_mean(total),
      old = geom_mean(old),
      .groups = "drop"
    ) %>%
    dplyr::filter(
      is.finite(ratio),
      is.finite(total),
      is.finite(old),
      is.finite(t),
      ratio > 0,
      ratio < 1,
      total > 0,
      old > 0
    ) %>%
    dplyr::arrange(t)
  
  df_old <- df_total %>%
    dplyr::filter(
      is.finite(old),
      old > 0
    )
  
  if (length(unique(df_old$t)) < min_points) return(NULL)
  
  fit_k <- tryCatch(
    lm(log(old) ~ t, data = df_old),
    error = function(e) NULL
  )
  
  if (is.null(fit_k)) return(NULL)
  
  k_raw <- -coef(fit_k)[["t"]]
  intercept_k <- coef(fit_k)[["(Intercept)"]]
  L0_fit <- exp(intercept_k)
  
  if (!is.finite(k_raw)) return(NULL)
  
  k <- max(k_raw, min_k)
  
  t0 <- min(df_total$t, na.rm = TRUE)
  P0 <- mean(df_total$total[df_total$t == t0], na.rm = TRUE)
  
  if (!is.finite(P0) || P0 <= 0) return(NULL)
  
  df_fit_s <- df_total %>%
    dplyr::filter(t > t0)
  
  if (nrow(df_fit_s) < 1) return(NULL)
  
  fit_s <- tryCatch(
    nls(
      total ~ P0 * exp(-k * (t - t0)) +
        (s / k) * (1 - exp(-k * (t - t0))),
      data = df_fit_s,
      start = list(s = k * P0),
      algorithm = "port",
      lower = c(s = 0),
      control = nls.control(maxiter = 200, warnOnly = TRUE)
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit_s)) {
    s_fit <- NA_real_
    r2_total_fit <- NA_real_
    rmse_total_fit <- NA_real_
    nrmse_total_fit <- NA_real_
  } else {
    s_fit <- coef(fit_s)[["s"]]
    
    pred_total <- predict(fit_s)
    obs_total <- df_fit_s$total
    
    ss_res <- sum((obs_total - pred_total)^2, na.rm = TRUE)
    ss_tot <- sum((obs_total - mean(obs_total, na.rm = TRUE))^2, na.rm = TRUE)
    
    r2_total_fit <- ifelse(ss_tot > 0, 1 - ss_res / ss_tot, NA_real_)
    
    rmse_total_fit <- sqrt(mean((obs_total - pred_total)^2, na.rm = TRUE))
    nrmse_total_fit <- rmse_total_fit / mean(obs_total, na.rm = TRUE)
  }
  
  tibble::tibble(
    Peptide = Peptide,
    point_number_total = nrow(df_total),
    point_number_old = nrow(df_old),
    P0 = P0,
    L0_fit = L0_fit,
    k_raw = k_raw,
    k = k,
    k_is_min_capped = k_raw < min_k,
    half_life = log(2) / k,
    synthesis_fit = s_fit,
    synthesis_norm = s_fit / P0,
    balance = s_fit / (k * P0),
    r2_old_decay = summary(fit_k)$r.squared,
    r2_total_fit = r2_total_fit,
    rmse_total_fit = rmse_total_fit,
    nrmse_total_fit = nrmse_total_fit
  )
}

#' Steady-state fit choosing the best of raw/mean strategies
#'
#' @inheritParams fit_silac
#' @return A one-row tibble with fit results, or NULL.
#' @export
fit_silac_mean_best <- function(ratio, t, Peptide, total = NULL) {
  stopifnot(length(ratio) == length(t))
  if (!is.null(total)) stopifnot(length(total) == length(t))
  
  t_unique <- sort(unique(t[is.finite(t)]))
  
  candidate_drop_t <- c(NA_real_, t_unique)
  
  res_list <- purrr::map(candidate_drop_t, function(drop_t) {
    
    if (is.na(drop_t)) {
      keep <- rep(TRUE, length(t))
    } else {
      keep <- t != drop_t
    }
    
    res <- tryCatch(
      fit_silac_mean(
        ratio = ratio[keep],
        t = t[keep],
        Peptide = Peptide,
        total = if (is.null(total)) NULL else total[keep]
      ),
      error = function(e) NULL
    )
    
    if (is.null(res)) return(NULL)
    
    res %>%
      dplyr::mutate(
        removed_t = drop_t,
        .before = point_number
      )
  })
  
  res_all <- dplyr::bind_rows(res_list)
  
  if (nrow(res_all) == 0) return(NULL)
  
  res_best <- res_all %>%
    dplyr::mutate(
      r2_for_rank = dplyr::if_else(is.finite(r2), r2, -Inf),
      no_remove = is.na(removed_t)
    ) %>%
    dplyr::arrange(
      dplyr::desc(r2_for_rank),
      dplyr::desc(no_remove),
      dplyr::desc(point_number)
    ) %>%
    dplyr::slice(1) %>%
    dplyr::select(-r2_for_rank, -no_remove)
  
  res_best
}

#' Non-steady-state fit choosing the best of raw/mean strategies
#'
#' @inheritParams fit_silac_nonsteady
#' @return A one-row tibble with fit results, or NULL.
#' @export
fit_silac_nonsteady_mean_best <- function(
    ratio,
    total,
    t,
    Peptide,
    min_points = 3,
    min_k = 0.0001
) {
  stopifnot(length(ratio) == length(t))
  stopifnot(length(total) == length(t))
  
  t_unique <- sort(unique(t[is.finite(t)]))
  
  candidate_drop_t <- c(NA_real_, t_unique)
  
  res_list <- purrr::map(candidate_drop_t, function(drop_t) {
    
    if (is.na(drop_t)) {
      keep <- rep(TRUE, length(t))
    } else {
      keep <- t != drop_t
    }
    
    res <- tryCatch(
      fit_silac_nonsteady_mean(
        ratio = ratio[keep],
        total = total[keep],
        t = t[keep],
        Peptide = Peptide,
        min_points = min_points,
        min_k = min_k
      ),
      error = function(e) NULL
    )
    
    if (is.null(res)) return(NULL)
    
    res %>%
      dplyr::mutate(
        removed_t = drop_t,
        .before = point_number_total
      )
  })
  
  res_all <- dplyr::bind_rows(res_list)
  
  if (nrow(res_all) == 0) return(NULL)
  
  res_best <- res_all %>%
    dplyr::mutate(
      r2_for_rank = dplyr::if_else(is.finite(r2_old_decay), r2_old_decay, -Inf),
      no_remove = is.na(removed_t)
    ) %>%
    dplyr::arrange(
      dplyr::desc(r2_for_rank),
      dplyr::desc(no_remove),
      dplyr::desc(point_number_old)
    ) %>%
    dplyr::slice(1) %>%
    dplyr::select(-r2_for_rank, -no_remove)
  
  res_best
}
