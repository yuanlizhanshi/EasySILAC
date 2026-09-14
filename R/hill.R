#' Fit global labeling efficiency with a Hill function
#'
#' @description
#' Fits the labeling-efficiency curve `E(x) = x^h / (K^h + x^h)` jointly to
#' the heavy-ratio time courses of the H100/H75/H50 samples (minpack.lm),
#' returning K, h, the dilution rate lambda and corrected (collapsed) ratios
#' `c(t) = observed / E(x)`.
#'
#' @param full_ratio A tibble with one row per sample and time point.
#' @param sample_col Column marking the labeling-efficiency sample (H100/H75/H50).
#' @param conc_col Column with the nominal heavy concentration.
#' @param time_col Column with the labeling time.
#' @param ratio_col Column with the (mean) heavy ratio.
#' @param start Named list of start values for K, h, lambda.
#' @param lower Named numeric vector of lower bounds.
#' @param upper Named numeric vector of upper bounds.
#' @param conc_check Named numeric vector of nominal concentrations.
#' @param conc_grid Grid of concentrations for evaluating E(x).
#' @param maxiter Maximum nlsLM iterations.
#'
#' @return A list with fit parameters (K, h, lambda), the efficiency function
#'   `E_fun`, corrected data frames, fit QC metrics and diagnostic plots.
#' @export
fit_hill_labelling_efficiency <- function(
    full_ratio,
    sample_col = "sample",
    conc_col = "conc",
    time_col = "t",
    ratio_col = "mean_ratio",
    start = list(
      K = 0.5,
      h = 3,
      lambda = 0.25
    ),
    lower = c(
      K = 1e-6,
      h = 0.1,
      lambda = 1e-6
    ),
    upper = c(
      K = 10,
      h = 20,
      lambda = 10
    ),
    conc_check = c(
      H100 = 1,
      H67 = 0.67,
      H33 = 0.33
    ),
    conc_grid = seq(0.01, 1.2, length.out = 300),
    maxiter = 1000
) {
  
  requireNamespace("dplyr")
  requireNamespace("tibble")
  requireNamespace("ggplot2")
  requireNamespace("minpack.lm")
  

  # 1. Prepare fitting dataframe

  
  df_fit <- full_ratio %>%
    dplyr::mutate(
      conc = as.numeric(.data[[conc_col]]),
      t = as.numeric(.data[[time_col]]),
      mean_ratio = as.numeric(.data[[ratio_col]]),
      sample = as.character(.data[[sample_col]])
    ) %>%
    dplyr::filter(
      is.finite(conc),
      is.finite(t),
      is.finite(mean_ratio),
      conc > 0,
      t > 0,
      mean_ratio >= 0,
      mean_ratio <= 1
    )
  
  if (nrow(df_fit) < 5) {
    stop("Too few valid data points after filtering.")
  }
  
  if (length(unique(df_fit$conc)) < 2) {
    stop("At least two concentration levels are required to fit concentration-dependent E(x).")
  }
  
  if (length(unique(df_fit$t)) < 3) {
    stop("At least three time points are recommended to fit c(t).")
  }
  

  # 2. Fit joint Hill model

  # mean_ratio = E(x) * c(t)
  #
  # E(x) = x^h / (K^h + x^h)
  # c(t) = 1 - exp(-lambda * t)

  
  fit_joint_hill <- minpack.lm::nlsLM(
    mean_ratio ~ 
      (conc^h / (K^h + conc^h)) *
      (1 - exp(-lambda * t)),
    
    data = df_fit,
    
    start = start,
    lower = lower,
    upper = upper,
    
    control = minpack.lm::nls.lm.control(
      maxiter = maxiter,
      ftol = 1e-10,
      ptol = 1e-10
    )
  )
  
  par_hill <- coef(fit_joint_hill)
  
  K_fit <- unname(par_hill["K"])
  h_fit <- unname(par_hill["h"])
  lambda_fit <- unname(par_hill["lambda"])
  
  E_fun_hill <- function(x) {
    x^h_fit / (K_fit^h_fit + x^h_fit)
  }
  
  c_fun <- function(t) {
    1 - exp(-lambda_fit * t)
  }
  
  ratio_fun <- function(conc, t) {
    E_fun_hill(conc) * c_fun(t)
  }

  
  full_ratio_check_hill <- df_fit %>%
    dplyr::mutate(
      E_x = E_fun_hill(conc),
      c_obs = mean_ratio / E_x,
      c_obs_clamped = pmin(pmax(c_obs, 0), 1),
      c_fit = c_fun(t),
      ratio_fit = ratio_fun(conc, t),
      residual = mean_ratio - ratio_fit
    )
  

  
  rss <- sum(full_ratio_check_hill$residual^2, na.rm = TRUE)
  tss <- sum((full_ratio_check_hill$mean_ratio - mean(full_ratio_check_hill$mean_ratio, na.rm = TRUE))^2, na.rm = TRUE)
  r2 <- 1 - rss / tss
  
  rmse <- sqrt(mean(full_ratio_check_hill$residual^2, na.rm = TRUE))
  nrmse_mean <- rmse / mean(full_ratio_check_hill$mean_ratio, na.rm = TRUE)
  
  fit_qc <- tibble::tibble(
    metric = c(
      "n_points",
      "n_samples",
      "n_concentrations",
      "n_timepoints",
      "RSS",
      "RMSE",
      "nRMSE_mean",
      "R2"
    ),
    value = c(
      nrow(df_fit),
      length(unique(df_fit$sample)),
      length(unique(df_fit$conc)),
      length(unique(df_fit$t)),
      rss,
      rmse,
      nrmse_mean,
      r2
    )
  )
  

  fit_summary <- tibble::tibble(
    parameter = c("K", "h", "lambda"),
    value = c(K_fit, h_fit, lambda_fit),
    meaning = c(
      "Half-maximal relative heavy amino-acid concentration",
      "Hill coefficient controlling concentration-response steepness",
      "Shared protein replacement rate in high-turnover reference proteins"
    )
  )
  

  
  max_ratio_by_sample <- df_fit %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(
      conc = dplyr::first(conc),
      max_ratio = max(mean_ratio, na.rm = TRUE),
      .groups = "drop"
    )
  
  E_df_check <- tibble::tibble(
    sample = names(conc_check),
    conc = as.numeric(conc_check)
  ) %>%
    dplyr::mutate(
      E_fit = E_fun_hill(conc)
    ) %>%
    dplyr::left_join(
      max_ratio_by_sample %>%
        dplyr::select(sample, max_ratio),
      by = "sample"
    ) %>%
    dplyr::mutate(
      max_corrected_ratio = max_ratio / E_fit
    )
  
  # -----------------------------
  # 7. Predicted E(x)
  # -----------------------------
  
  E_pred <- tibble::tibble(
    conc = conc_grid
  ) %>%
    dplyr::mutate(
      E_fit = E_fun_hill(conc)
    )
  
  # -----------------------------
  # 8. Plots
  # -----------------------------
  
  p_E <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = E_df_check,
      ggplot2::aes(x = conc, y = E_fit),
      size = 3
    ) +
    ggplot2::geom_line(
      data = E_pred,
      ggplot2::aes(x = conc, y = E_fit),
      linewidth = 1.2
    ) +
    ggplot2::scale_x_continuous(
      breaks = as.numeric(conc_check),
      labels = names(conc_check)
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Relative heavy amino-acid concentration",
      y = "Fitted labelling efficiency E(x)",
      title = "Concentration-dependent labelling efficiency"
    ) +
    ggplot2::theme_test()
  
  p_collapse <- ggplot2::ggplot(
    full_ratio_check_hill,
    ggplot2::aes(x = t)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(y = c_obs, color = sample),
      size = 3
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = c_fit),
      color = "red",
      linewidth = 1.2
    ) +
    ggplot2::labs(
      x = "Time (t)",
      y = "c(t) = observed ratio / E(x)",
      title = "Collapse after labelling-efficiency correction"
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1.1)) +
    ggplot2::theme_test()
  
  p_obs_fit <- ggplot2::ggplot(
    full_ratio_check_hill,
    ggplot2::aes(x = t)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(y = mean_ratio, color = sample),
      size = 3
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = ratio_fit, color = sample),
      linewidth = 1.2
    ) +
    ggplot2::labs(
      x = "Time (t)",
      y = "Observed heavy ratio",
      title = "Observed ratio vs fitted ratio"
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_test()

  
  list(
    fit = fit_joint_hill,
    parameters = par_hill,
    K = K_fit,
    h = h_fit,
    lambda = lambda_fit,
    E_fun = E_fun_hill,
    c_fun = c_fun,
    ratio_fun = ratio_fun,
    df_fit = df_fit,
    df_corrected = full_ratio_check_hill,
    E_df_check = E_df_check,
    E_pred = E_pred,
    fit_summary = fit_summary,
    fit_qc = fit_qc,
    plots = list(
      p_E = p_E,
      p_collapse = p_collapse,
      p_obs_fit = p_obs_fit
    )
  )
}
