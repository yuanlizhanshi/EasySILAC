#' PCA score plot
#'
#' @description Note: unlike `EasyProtein::plot_pca`, this version supports an
#'   additional `shape_by` aesthetic.
#'
#' @param pca.df Data frame of PCA coordinates with sample metadata.
#' @param pca.res PCA result object (e.g. from FactoMineR::PCA).
#' @param colorby Column of `pca.df` mapped to color.
#' @param shape_by Optional column mapped to point shape.
#' @param label Label mode passed to the text layer.
#'
#' @return A ggplot object.
#' @export
plot_pca <- function(pca.df, pca.res, colorby, shape_by = NULL,label) {
  
  base_plot <- ggplot(pca.df, aes(x = Dim.1, y = Dim.2)) +
    geom_point(aes_string(color = colorby,shape = shape_by), size = rel(2.5)) +
    labs(
      x = paste0("PC1 ", round(pca.res$eig[1, 2], 2), "%"),
      y = paste0("PC2 ", round(pca.res$eig[2, 2], 2), "%"),
      color = NULL
    ) +
    theme_test() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold")
    )
  
  no_label <- is.null(label) || isFALSE(label) || identical(label, "NULL")
  
  if (no_label) {
    return(base_plot)
  } else {
    if (!is.character(label) || length(label) != 1 || !nzchar(label)) {
      stop("'label' must be a column name, FALSE, or 'NULL'.")
    }
    if (!label %in% colnames(pca.df)) {
      stop(paste0("Column '", label, "' not found in pca.df"))
    }
    return(
      base_plot +
        ggrepel::geom_text_repel(aes_string(label = label))
    )
  }
}

#' Density scatter plot on log10 scales
#'
#' @description Density-colored scatter (ggpointdensity + viridis) with
#'   log10 axes, identity line and Pearson correlation annotation.
#'
#' @param x Name of the x column.
#' @param y Name of the y column.
#' @param min_value Axis minimum (both axes).
#' @param max_value Axis maximum (both axes).
#' @param xlabs X axis label.
#' @param ylabs Y axis label.
#' @param title Plot title.
#' @param df Data frame.
#'
#' @return A ggplot object.
#' @export
plot_scatter_log10_v1 <- function(
    x,
    y,
    min_value = NULL,
    max_value = NULL,
    xlabs = NULL,
    ylabs = NULL,
    title = NULL,
    df
) {
  df_tmp <- df %>% dplyr::select(x, y)

  if (is.null(min_value)) {
    min_value <- min(df_tmp[[x]], df_tmp[[y]], na.rm = TRUE) * 0.95
  }
  if (is.null(max_value)) {
    max_value <- max(df_tmp[[x]], df_tmp[[y]], na.rm = TRUE) * 1.05
  }
  if (is.null(xlabs)) {
    xlabs <- glue::glue(x)
  }
  if (is.null(ylabs)) {
    ylabs <- glue::glue(y)
  }
  df_tmp <- df_tmp %>%
    dplyr::filter(.data[[x]] > min_value, .data[[y]] > min_value)
  p <- ggplot(df_tmp, aes(!!sym(x), !!sym(y))) +
    geom_point() +
    coord_fixed(ratio = 1) +
    scale_x_continuous(
      limits = c(min_value, max_value),
      labels = scales::label_number(),
      trans = 'log10'
    ) +
    scale_y_continuous(
      limits = c(min_value, max_value),
      labels = scales::label_number(),
      trans = 'log10'
    ) +
    ggpubr::stat_cor() +
    theme_test() +
    geom_abline(slope = 1, color = "blue", linetype = "dashed") +
    ggpointdensity::geom_pointdensity() +
    viridis::scale_color_viridis() +
    guides(color = 'none') +
    labs(
      x = xlabs,
      y = ylabs,
      title = title,
      subtitle = paste0('N = ', nrow(df_tmp))
    )
  return(p)
}

#' Extract observed ratios of one peptide for plotting
#'
#' @param se SummarizedExperiment.
#' @param peptide Peptide row name.
#' @param time_col `colData` column with labeling time.
#' @param ratio_assay Assay with the heavy ratio.
#' @param sample_col `colData` column with the sample group.
#'
#' @return A tibble with time, ratio and sample.
#' @export
get_data_to_plot <- function(
    se,
    peptide,
    time_col   = "t",
    ratio_assay  = "ratio",
    sample_col = "conc"   # H100 / H067 / H033
) {


  df_obs <- tibble(
    t      = as.numeric(colData(se)[[time_col]]),
    ntr    = as.numeric(assay(se, ratio_assay)[peptide, ]),
    sample = colData(se)[[sample_col]]
  ) %>%
    dplyr::filter(!is.na(t), !is.na(ntr))

  df_obs$sample <- factor(df_obs$sample, levels = unique(df_obs$sample))

  if (nrow(df_obs) < 3)
    stop("Not enough points after removing NA")


  rd <- as.data.frame(rowData(se)[peptide, ])

  # 找 k 列
  k_cols <- grep("_k$", colnames(rd), value = TRUE)

  df_pred <- purrr::map_dfr(k_cols, function(k_col) {

    sample_name <- sub("_k$", "", k_col)

    k <- as.numeric(rd[[k_col]])
    a <- as.numeric(rd[[paste0(sample_name, "_intercept")]])

    if (is.na(k) | is.na(a)) return(NULL)

    t_range <- df_obs %>%
      dplyr::filter(sample == sample_name) %>%
      dplyr::pull(t)

    if (length(t_range) < 2) return(NULL)

    t_pred <- seq(min(t_range), max(t_range), length.out = 200)

    C <- exp(a)
    ntr_pred <- 1 - C * exp(-k * t_pred)

    tibble(
      t      = t_pred,
      ntr    = ntr_pred,
      sample = sample_name,
      H0     = 1 - C
    )
  })

  return(list(
    obs  = df_obs,
    pred = df_pred
  ))
}

#' Density plot of a numeric column
#'
#' @param data Data frame.
#' @param column Column name to plot.
#' @param title Plot title.
#' @param add_2_fold Add dashed lines at two-fold changes.
#'
#' @return A ggplot object.
#' @export
plot_density <- function(data, column = NULL, title = NULL, add_2_fold = FALSE) {
  
  if (is.null(column)) {
    x <- data
    colname <- deparse(substitute(data))
  } else {
    colname <- rlang::as_string(rlang::ensym(column))
    x <- data[[colname]]
  }
  
  x <- x[is.finite(x)]
  
  if (length(x) < 2) {
    stop("Fewer than 2 finite values; cannot draw the density plot. Check for NA/NaN/Inf or the column type.")
  }
  
  value_mean <- mean(x)
  value_med  <- median(x)
  subtitle_text <- glue::glue(
    "N = {length(x)}\nmean = {round(value_mean, 2)} | median = {round(value_med, 2)}"
  )

  if (isTRUE(add_2_fold)) {
    fold_range_lower <- min(value_med * 0.5, value_med * 2)
    fold_range_upper <- max(value_med * 0.5, value_med * 2)
    pct_in_2_fold <- mean(x >= fold_range_lower & x <= fold_range_upper) * 100
    subtitle_text <- glue::glue(
      "{subtitle_text}\nwithin 0.5-2x median = {sprintf('%.1f%%', pct_in_2_fold)}"
    )
  }
  
  if (is.null(title)) {
    title <- colname
  }
  
  ggplot(data.frame(x = x), aes(x = x)) +
    geom_density(fill = "#2c7fb8", alpha = 0.4) +
    labs(
      x = colname,
      subtitle = subtitle_text,
      title = title
    ) +
    theme_test()
}

#' Plot one peptide's ratio time course
#'
#' @param se SummarizedExperiment.
#' @param pep Peptide row name.
#'
#' @return A ggplot object.
#' @export
plot_pep <- function(se,pep){
  res <- get_data_to_plot(se, pep)


  p <- ggplot() +
    geom_point(data = res$obs,
               aes(t, ntr, color = sample)) +
    geom_line(data = res$pred,
              aes(t, ntr, color = sample)) +
    ylim(0, 1) +
    scale_x_continuous(
      breaks = c(1, 2, 4, 8),
      labels = c("1", "2", "4", "8")
    )+
    labs(
      x = "Time",
      y = "H/(H+L)",
      title = pep
    ) +
    theme_test()
  return(p)
}

#' Plot a steady-state SILAC fit
#'
#' @description Observed ratios (black points) and the fitted curve
#'   `1 - C * exp(-k * t)` (red line), annotated with k, H0, half-life and R2.
#'
#' @param se SummarizedExperiment with the ratio assay.
#' @param fit_df Data frame with fit results (columns `gene`, `k`,
#'   `intercept`, `r2`, `half_life`).
#' @param peptide Gene/peptide to plot.
#' @param plot_title Plot title.
#' @param time_col `colData` column with labeling time.
#' @param plot_type "ratio" or "intensity".
#' @param ratio_assay Assay with the (possibly corrected) heavy ratio.
#' @param conc_assay Assay with total intensity (intensity mode).
#'
#' @return A ggplot object.
#' @export
plot_protein_silac_fit <- function(
    se,
    fit_df,
    peptide,
    plot_title = NULL,
    time_col = "t",
    plot_type = c("ratio", "intensity"),
    ratio_assay = "ratio",
    conc_assay = "conc"
) {
  
  plot_type <- match.arg(plot_type)
  
  
  t <- as.numeric(SummarizedExperiment::colData(se)[[time_col]])
  
  row <- fit_df %>%
    dplyr::filter(gene == peptide)
  
  if (nrow(row) == 0) {
    stop("Gene not found in fit_df")
  }
  
  k  <- as.numeric(row$k[1])
  a  <- as.numeric(row$intercept[1])
  r2 <- as.numeric(row$r2[1])
  hl <- as.numeric(row$half_life[1])
  
  C <- exp(a)
  
  # ============================================================
  # ratio mode
  # ============================================================
  if (plot_type == "ratio") {
    
    ratio <- as.numeric(
      SummarizedExperiment::assay(se, ratio_assay)[peptide, ]
    )
    
    df_obs <- data.frame(
      t = t,
      ratio = ratio
    ) %>%
      dplyr::filter(
        is.finite(t),
        is.finite(ratio)
      ) %>%
      dplyr::arrange(t)
    
    if (nrow(df_obs) < 3) {
      stop("Not enough points after removing NA")
    }
    
    t_pred <- seq(
      min(df_obs$t),
      max(df_obs$t),
      length.out = 200
    )
    
    ratio_pred <- 1 - C * exp(-k * t_pred)
    
    df_pred <- data.frame(
      t = t_pred,
      ratio = ratio_pred
    )
    
    p <- ggplot() +
      geom_point(
        data = df_obs,
        aes(x = t, y = ratio),
        color = "black",
        size = 3
      ) +
      geom_line(
        data = df_pred,
        aes(x = t, y = ratio),
        color = "red",
        linewidth = 1.2
      ) +
      scale_x_continuous(
        breaks = sort(unique(df_obs$t))
      ) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(
        x = "Time",
        y = "H / (H + L)",
        subtitle = peptide,
        title = plot_title
      ) +
      annotate(
        "text",
        x = max(df_obs$t) * 0.55,
        y = 0.25,
        hjust = 0,
        label = paste0(
          "k = ", signif(k, 3),
          "\nH0 = ", round(1 - C, 3),
          "\nHalf-life = ", signif(hl, 3),
          "\nR\u00B2 = ", round(r2, 2)
        )
      ) +
      theme_test(base_size = 14)
    
    return(p)
  }
  
  # ============================================================
  # intensity mode
  # only use samples with valid ratio, H, and L
  # ============================================================
  if (plot_type == "intensity") {
    
    conc <- as.numeric(
      SummarizedExperiment::assay(se, conc_assay)[peptide, ]
    )
    
    ratio <- as.numeric(
      SummarizedExperiment::assay(se, ratio_assay)[peptide, ]
    )
    
    df_obs <- data.frame(
      t = t,
      conc = conc,
      ratio = ratio
    ) %>%
      dplyr::filter(
        is.finite(t),
        is.finite(conc),
        is.finite(ratio),
        conc > 0
      ) %>%
      dplyr::mutate(
        ratio = pmin(pmax(ratio, 1e-6), 1 - 1e-6),
        new = conc * ratio,
        old = conc * (1 - ratio),
        total = conc
      ) %>%
      dplyr::arrange(t)
    
    if (nrow(df_obs) < 3) {
      stop("Not enough points after removing NA")
    }
    
    message(
      "Intensity mode uses ",
      nrow(df_obs),
      " samples after matching conc-valid and ratio-valid samples."
    )
    
    df_long <- rbind(
      data.frame(
        t = df_obs$t,
        value = df_obs$new,
        type = "new"
      ),
      data.frame(
        t = df_obs$t,
        value = df_obs$old,
        type = "old"
      ),
      data.frame(
        t = df_obs$t,
        value = df_obs$total,
        type = "total"
      )
    )
    
    df_total_mean <- df_obs %>%
      dplyr::group_by(t) %>%
      dplyr::summarise(
        value = mean(total, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(type = "total")
    
    t_pred <- seq(
      min(df_obs$t, na.rm = TRUE),
      max(df_obs$t, na.rm = TRUE),
      length.out = 200
    )
    
    N0 <- df_obs %>%
      dplyr::filter(t == min(t, na.rm = TRUE)) %>%
      dplyr::summarise(
        N0 = mean(total, na.rm = TRUE)
      ) %>%
      dplyr::pull(N0)
    
    old_frac_pred <- exp(a - k * t_pred)
    old_frac_pred <- pmin(pmax(old_frac_pred, 0), 1)
    new_frac_pred <- 1 - old_frac_pred
    
    old_pred <- N0 * old_frac_pred
    new_pred <- N0 * new_frac_pred
    
    df_pred <- rbind(
      data.frame(
        t = t_pred,
        value = new_pred,
        type = "new"
      ),
      data.frame(
        t = t_pred,
        value = old_pred,
        type = "old"
      )
    ) %>%
      dplyr::filter(
        is.finite(value),
        value >= 0
      )
    
    p <- ggplot() +
      geom_point(
        data = df_long,
        aes(x = t, y = value, color = type),
        size = 2.5,
        alpha = 0.85
      ) +
      geom_line(
        data = df_total_mean,
        aes(x = t, y = value, color = type),
        linewidth = 0.9,
        linetype = 2
      ) +
      geom_line(
        data = df_pred,
        aes(x = t, y = value, color = type),
        linewidth = 1.0,
        linetype = 2
      ) +
      scale_color_manual(
        values = c(
          new = "#db6968",
          old = "#88c4e8",
          total = "grey60"
        )
      ) +
      scale_x_continuous(
        breaks = sort(unique(df_obs$t))
      ) +
      labs(
        x = "Time",
        y = "Concentration",
        title = plot_title,
        subtitle = peptide,
        color = NULL
      ) +
      theme_test(base_size = 14)
    
    return(p)
  }
}

#' Plot a non-steady-state SILAC fit
#'
#' @description Observed ratios and the non-steady-state fitted curve,
#'   annotated with k, synthesis rate, half-life, R2 of the old-channel decay
#'   and nRMSE of the total fit.
#'
#' @param se SummarizedExperiment.
#' @param fit_df Data frame with fit results (columns `gene`/`Peptide`, `k`,
#'   `synthesis_fit`, `half_life`, `r2_old_decay`, `nrmse_total_fit`,
#'   `L0_fit`).
#' @param peptide Gene/peptide to plot.
#' @param plot_title Plot title.
#' @param time_col `colData` column with labeling time.
#' @param plot_type "ratio" or "intensity".
#' @param ratio_assay Assay with the heavy ratio.
#' @param conc_assay Assay with total intensity.
#' @param eps Small constant to keep ratios inside (0, 1).
#'
#' @return A ggplot object.
#' @export
plot_protein_silac_nonsteady_fit <- function(
    se,
    fit_df,
    peptide,
    plot_title = NULL,
    time_col = "t",
    plot_type = c("ratio", "intensity"),
    ratio_assay = "ratio",
    conc_assay = "conc",
    eps = 1e-6
) {
  
  plot_type <- match.arg(plot_type)
  t <- as.numeric(SummarizedExperiment::colData(se)[[time_col]])
  
  # 兼容 gene 或 Peptide 列
  if ("gene" %in% colnames(fit_df)) {
    row <- fit_df %>%
      dplyr::filter(gene == peptide)
  } else if ("Peptide" %in% colnames(fit_df)) {
    row <- fit_df %>%
      dplyr::filter(Peptide == peptide)
  } else {
    stop("fit_df must contain either gene or Peptide column.")
  }
  
  if (nrow(row) == 0) {
    stop("Peptide not found in fit_df")
  }
  
  if (nrow(row) > 1) {
    row <- row[1, ]
  }
  
  k <- as.numeric(row$k[1])
  s_fit <- as.numeric(row$synthesis_fit[1])
  hl <- as.numeric(row$half_life[1])
  r2_old <- as.numeric(row$r2_old_decay[1])
  nrmse_total_fit <- as.numeric(row$nrmse_total_fit[1])
  L0_fit <- as.numeric(row$L0_fit[1])
  
  if (!is.finite(k) || k <= 0) stop("Invalid k.")
  if (!is.finite(s_fit)) stop("Invalid synthesis_fit.")
  
  ratio <- as.numeric(SummarizedExperiment::assay(se, ratio_assay)[peptide, ])
  conc <- as.numeric(SummarizedExperiment::assay(se, conc_assay)[peptide, ])
  
  df_obs0 <- data.frame(
    t = t,
    conc = conc,
    ratio = ratio
  ) %>%
    dplyr::filter(
      is.finite(t),
      is.finite(conc),
      is.finite(ratio),
      conc > 0
    ) %>%
    dplyr::mutate(
      ratio = pmin(pmax(ratio, eps), 1 - eps),
      new_ratio = ratio,
      old_ratio = 1 - ratio,
      new = conc * new_ratio,
      old = conc * old_ratio,
      total = conc
    ) %>%
    dplyr::arrange(t)
  
  if (nrow(df_obs0) < 3) {
    stop("Not enough points after removing NA")
  }
  
  t0 <- min(df_obs0$t, na.rm = TRUE)
  
  P0_raw <- df_obs0 %>%
    dplyr::filter(t == t0) %>%
    dplyr::summarise(P0 = mean(total, na.rm = TRUE)) %>%
    dplyr::pull(P0)
  
  L0_raw <- df_obs0 %>%
    dplyr::filter(t == t0) %>%
    dplyr::summarise(L0 = mean(old, na.rm = TRUE)) %>%
    dplyr::pull(L0)
  
  if (!is.finite(P0_raw) || P0_raw <= 0) stop("Invalid P0.")
  if (!is.finite(L0_raw) || L0_raw <= 0) stop("Invalid L0.")
  
  t_pred <- seq(
    min(df_obs0$t, na.rm = TRUE),
    max(df_obs0$t, na.rm = TRUE),
    length.out = 200
  )
  
  tau_pred <- t_pred - t0
  
  # ------------------------------------------------------------
  # 核心修正：
  # L0_fit 来自 lm(log(old) ~ t)，所以 L0_fit 是 t = 0 的 old
  # 不能写成 L0_fit * exp(-k * (t - t0))
  # ------------------------------------------------------------
  if (is.finite(L0_fit) && L0_fit > 0) {
    old_pred_raw <- L0_fit * exp(-k * t_pred)
  } else {
    old_pred_raw <- L0_raw * exp(-k * tau_pred)
  }
  
  # total fit 仍然从 t0 的 P0 开始
  total_pred_raw <- P0_raw * exp(-k * tau_pred) +
    (s_fit / k) * (1 - exp(-k * tau_pred))
  
  new_pred_raw <- total_pred_raw - old_pred_raw
  
  # 防止噪音导致 new prediction 小于 0
  new_pred_raw[new_pred_raw < 0] <- 0
  
  new_ratio_pred <- pmin(pmax(new_pred_raw / total_pred_raw, 0), 1)
  
  if (plot_type == "ratio") {
    
    df_pred_ratio <- data.frame(
      t = t_pred,
      new_ratio = new_ratio_pred
    ) %>%
      dplyr::filter(is.finite(t), is.finite(new_ratio))
    
    p <- ggplot() +
      geom_point(data = df_obs0, aes(x = t, y = new_ratio), color = "black", size = 3) +
      geom_line(data = df_pred_ratio, aes(x = t, y = new_ratio), color = "red", linewidth = 1.2) +
      scale_x_continuous(breaks = sort(unique(df_obs0$t))) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(
        x = "Time",
        y = "New fraction / ratio",
        title = plot_title,
        subtitle = peptide
      ) +
      annotate(
        "text",
        x = max(df_obs0$t, na.rm = TRUE) * 0.55,
        y = 0.45,
        hjust = 0,
        label = paste0(
          "k = ", signif(k, 3),
          "\ns = ", signif(s_fit, 3),
          "\nHalf-life = ", signif(hl, 3),
          "\nR\u00B2 old = ", round(r2_old, 2),
          "\nnRMSE = ", round(nrmse_total_fit, 2)
        )
      ) +
      theme_test(base_size = 14)
    
    return(p)
  }
  
  if (plot_type == "intensity") {
    
    message(
      "Raw concentration-derived intensity mode uses ",
      nrow(df_obs0),
      " samples."
    )
    
    df_obs <- df_obs0
    
    df_long <- rbind(
      data.frame(t = df_obs$t, value = df_obs$new, type = "new"),
      data.frame(t = df_obs$t, value = df_obs$old, type = "old"),
      data.frame(t = df_obs$t, value = df_obs$total, type = "total")
    )
    
    df_total_mean <- df_obs %>%
      dplyr::group_by(t) %>%
      dplyr::summarise(value = mean(total, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(type = "total")
    
    # 这里只画合成 new 和降解 old 的拟合线
    # 不额外画 total fitted，避免图例乱掉
    df_pred <- rbind(
      data.frame(t = t_pred, value = new_pred_raw, type = "new"),
      data.frame(t = t_pred, value = old_pred_raw, type = "old")
    ) %>%
      dplyr::filter(is.finite(value), value >= 0)
    
    p <- ggplot() +
      geom_point(
        data = df_long,
        aes(x = t, y = value, color = type),
        size = 2.5,
        alpha = 0.85
      ) +
      geom_line(
        data = df_total_mean,
        aes(x = t, y = value, color = type),
        linewidth = 0.9,
        linetype = 2
      ) +
      geom_line(
        data = df_pred,
        aes(x = t, y = value, color = type),
        linewidth = 1.0,
        linetype = 2
      ) +
      scale_color_manual(
        values = c(new = "#db6968", old = "#88c4e8", total = "grey60"),
        breaks = c("new", "old", "total")
      ) +
      scale_x_continuous(breaks = sort(unique(df_obs$t))) +
      labs(
        x = "Time",
        y = "Concentration",
        title = plot_title,
        subtitle = peptide,
        color = NULL
      ) +
      theme_test(base_size = 14)
    
    return(p)
  }
}

#' Plot a labeling-efficiency-corrected non-steady-state fit
#'
#' @description Like [plot_protein_silac_nonsteady_fit()], for ratios already
#'   corrected by the global Hill labeling efficiency.
#'
#' @inheritParams plot_protein_silac_nonsteady_fit
#'
#' @return A ggplot object.
#' @export
plot_protein_silac_nonsteady_adjust_fit <-function(
    se,
    fit_df,
    peptide,
    plot_title = NULL,
    time_col = "t",
    plot_type = c("ratio", "intensity"),
    ratio_assay = "ratio",
    conc_assay = "conc",
    eps = 1e-6
) {
  
  plot_type <- match.arg(plot_type)
  
  # ============================================================
  # time
  # ============================================================
  t <- as.numeric(SummarizedExperiment::colData(se)[[time_col]])
  
  # ============================================================
  # get fit parameters
  # 兼容 Peptide / gene 两种列名
  # ============================================================
  if ("Peptide" %in% colnames(fit_df)) {
    row <- fit_df %>% dplyr::filter(Peptide == peptide)
  } else if ("gene" %in% colnames(fit_df)) {
    row <- fit_df %>% dplyr::filter(gene == peptide)
  } else {
    stop("fit_df must contain either 'Peptide' or 'gene' column.")
  }
  
  if (nrow(row) == 0) {
    stop("Peptide not found in fit_df.")
  }
  
  k <- as.numeric(row$k[1])
  s_fit <- as.numeric(row$synthesis_fit[1])
  hl <- as.numeric(row$half_life[1])
  r2_old <- as.numeric(row$r2_old_decay[1])
  nrmse_total_fit <- as.numeric(row$nrmse_total_fit[1])
  L0_fit <- as.numeric(row$L0_fit[1])
  
  if (!is.finite(k) || k <= 0) {
    stop("Invalid k.")
  }
  
  if (!is.finite(s_fit)) {
    stop("Invalid synthesis_fit.")
  }
  
  # ============================================================
  # observed data
  # 注意：这里默认 ratio_assay 已经是矫正后的 ratio
  # ============================================================
  if (!peptide %in% rownames(se)) {
    stop("peptide not found in rownames(se).")
  }
  
  ratio <- as.numeric(
    SummarizedExperiment::assay(se, ratio_assay)[peptide, ]
  )
  
  conc <- as.numeric(
    SummarizedExperiment::assay(se, conc_assay)[peptide, ]
  )
  
  df_obs0 <- data.frame(
    t = t,
    conc = conc,
    ratio = ratio
  ) %>%
    dplyr::filter(
      is.finite(t),
      is.finite(conc),
      is.finite(ratio),
      conc > 0
    ) %>%
    dplyr::mutate(
      ratio = pmin(pmax(ratio, eps), 1 - eps),
      new_ratio = ratio,
      old_ratio = 1 - ratio,
      new = conc * new_ratio,
      old = conc * old_ratio,
      total = conc
    ) %>%
    dplyr::arrange(t)
  
  if (nrow(df_obs0) < 3) {
    stop("Not enough points after removing NA.")
  }
  
  message(
    "Using ",
    nrow(df_obs0),
    " samples after matching conc-valid and ratio-valid samples."
  )
  
  # ============================================================
  # initial values
  # ============================================================
  t0 <- min(df_obs0$t, na.rm = TRUE)
  
  P0_raw <- df_obs0 %>%
    dplyr::filter(t == t0) %>%
    dplyr::summarise(
      P0 = mean(total, na.rm = TRUE)
    ) %>%
    dplyr::pull(P0)
  
  L0_raw_at_t0 <- df_obs0 %>%
    dplyr::filter(t == t0) %>%
    dplyr::summarise(
      L0 = mean(old, na.rm = TRUE)
    ) %>%
    dplyr::pull(L0)
  
  if (!is.finite(P0_raw) || P0_raw <= 0) {
    stop("Invalid P0.")
  }
  
  # ============================================================
  # prediction
  #
  # 关键修正：
  # fit_silac_nonsteady 里面 L0_fit 来自 lm(log(old) ~ t)
  # 所以 L0_fit 是 old(t = 0)，不是 old(t = t0)
  #
  # 因此 old prediction 应该用：
  # old(t) = L0_fit * exp(-k * t)
  #
  # 或者先转换：
  # L0_at_t0 = L0_fit * exp(-k * t0)
  # old(t) = L0_at_t0 * exp(-k * (t - t0))
  # ============================================================
  t_pred <- seq(
    min(df_obs0$t, na.rm = TRUE),
    max(df_obs0$t, na.rm = TRUE),
    length.out = 200
  )
  
  tau_pred <- t_pred - t0
  
  if (is.finite(L0_fit) && L0_fit > 0) {
    L0_at_t0 <- L0_fit * exp(-k * t0)
  } else {
    L0_at_t0 <- L0_raw_at_t0
  }
  
  if (!is.finite(L0_at_t0) || L0_at_t0 <= 0) {
    stop("Invalid L0_at_t0.")
  }
  
  old_pred_raw <- L0_at_t0 * exp(-k * tau_pred)
  
  total_pred_raw <- P0_raw * exp(-k * tau_pred) +
    (s_fit / k) * (1 - exp(-k * tau_pred))
  
  new_pred_raw <- total_pred_raw - old_pred_raw
  
  new_ratio_pred <- new_pred_raw / total_pred_raw
  new_ratio_pred <- pmin(pmax(new_ratio_pred, 0), 1)
  
  # ============================================================
  # ratio mode
  # ============================================================
  if (plot_type == "ratio") {
    
    df_pred_ratio <- data.frame(
      t = t_pred,
      new_ratio = new_ratio_pred
    ) %>%
      dplyr::filter(
        is.finite(t),
        is.finite(new_ratio)
      )
    
    x_text <- min(df_obs0$t, na.rm = TRUE) +
      0.08 * diff(range(df_obs0$t, na.rm = TRUE))
    
    p <- ggplot() +
      geom_point(
        data = df_obs0,
        aes(x = t, y = new_ratio),
        color = "black",
        size = 3
      ) +
      geom_line(
        data = df_pred_ratio,
        aes(x = t, y = new_ratio),
        color = "red",
        linewidth = 1.2
      ) +
      scale_x_continuous(
        breaks = sort(unique(df_obs0$t))
      ) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(
        x = "Time",
        y = "New fraction / ratio",
        title = plot_title,
        subtitle = peptide
      ) +
      annotate(
        "text",
        x = x_text,
        y = 0.82,
        hjust = 0,
        vjust = 1,
        fontface = "bold",
        label = paste0(
          "k = ", signif(k, 3),
          "\ns = ", signif(s_fit, 3),
          "\nHalf-life = ", signif(hl, 3),
          "\nR\u00B2 old = ", round(r2_old, 2),
          "\nnRMSE = ", round(nrmse_total_fit, 2)
        )
      ) +
      theme_test(base_size = 14)
    
    return(p)
  }
  
  # ============================================================
  # intensity mode
  # ============================================================
  if (plot_type == "intensity") {
    
    message(
      "Raw concentration-derived intensity mode uses ",
      nrow(df_obs0),
      " samples."
    )
    
    df_obs <- df_obs0
    
    df_long <- rbind(
      data.frame(
        t = df_obs$t,
        value = df_obs$new,
        type = "new"
      ),
      data.frame(
        t = df_obs$t,
        value = df_obs$old,
        type = "old"
      ),
      data.frame(
        t = df_obs$t,
        value = df_obs$total,
        type = "total"
      )
    )
    
    df_pred <- rbind(
      data.frame(
        t = t_pred,
        value = new_pred_raw,
        type = "new_fit"
      ),
      data.frame(
        t = t_pred,
        value = old_pred_raw,
        type = "old_fit"
      )
    ) %>%
      dplyr::filter(
        is.finite(t),
        is.finite(value)
      )
    
    df_pred_new <- df_pred %>%
      dplyr::filter(type == "new_fit")
    
    df_pred_old <- df_pred %>%
      dplyr::filter(type == "old_fit")
    
    p <- ggplot() +
      # observed points
      geom_point(
        data = df_long,
        aes(x = t, y = value, color = type),
        size = 2.5,
        alpha = 0.85
      ) +
      # observed lines: new / old / total 都只连接真实点
      geom_line(
        data =       data.frame(
          t = df_obs$t,
          value = df_obs$total,
          type = "total"
        ),
        aes(x = t, y = value, color = type, group = type),
        linewidth = 0.7,
        alpha = 0.75
      ) +
      # new fitted line，不进入图例
      geom_line(
        data = df_pred_new,
        aes(x = t, y = value),
        color = "#db6968",
        linewidth = 1.0,
        linetype = 2,
        show.legend = FALSE
      ) +
      # old fitted line，不进入图例
      geom_line(
        data = df_pred_old,
        aes(x = t, y = value),
        color = "#88c4e8",
        linewidth = 1.0,
        linetype = 2,
        show.legend = FALSE
      ) +
      scale_color_manual(
        values = c(
          new = "#db6968",
          old = "#88c4e8",
          total = "grey60"
        ),
        breaks = c(
          "new",
          "old",
          "total"
        )
      ) +
      scale_x_continuous(
        breaks = sort(unique(df_obs$t))
      ) +
      labs(
        x = "Time",
        y = "Concentration",
        title = plot_title,
        subtitle = peptide,
        color = NULL
      ) +
      theme_test(base_size = 14)
    
    return(p)
  }
}
