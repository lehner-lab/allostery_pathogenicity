# ── Shared Helper Functions ───────────────────────────────────────────────────
# Reusable functions extracted from src/*.Rmd scripts.
# Loaded automatically via load.project().

# ── PDB B-factor Mapping ─────────────────────────────────────────────────────

#' Map per-residue values to PDB B-factor column for ChimeraX visualization.
#'
#' Writes a new PDB with B-factors replaced by the values in `data`.
#' Residues not in `data` are set to PDB_OUTLIER_BFACTOR (999) so they
#' can be filtered in ChimeraX with `select @@bfactor>=999`.
#'
#' @param data        Data frame with position and value columns.
#' @param pdb_path    Path to input PDB file.
#' @param output_path Path to write the modified PDB.
#' @param position_col Column name for residue numbers (default "position").
#' @param value_col   Column name for the value to map (default "value").
#' @return Invisible list with min/max values and output path.
map_residuals_to_pdb <- function(data,
                                 pdb_path,
                                 output_path,
                                 position_col = "position",
                                 value_col = "value") {
  pdb <- bio3d::read.pdb(pdb_path, rm.alt = TRUE)
  new_b <- pdb$atom$b

  positions <- data[[position_col]]
  values <- data[[value_col]]

  for (i in seq_along(positions)) {
    idx <- which(pdb$atom$resno == positions[i])
    new_b[idx] <- values[i]
  }

  # Mark unmapped residues with sentinel value
  matched <- which(pdb$atom$resno %in% positions)
  new_b[setdiff(seq_along(new_b), matched)] <- PDB_OUTLIER_BFACTOR

  pdb$atom$b <- new_b
  bio3d::write.pdb(pdb, file = output_path)

  invisible(list(
    min_value = min(values, na.rm = TRUE),
    max_value = max(values, na.rm = TRUE),
    output_file = output_path
  ))
}


# ── Bootstrap Exponential Decay ──────────────────────────────────────────────

#' Fit exponential decay y = a * exp(-b * x) with bootstrap confidence bands.
#'
#' @param data    Data frame.
#' @param x_col   Column name for distance (x-axis).
#' @param y_col   Column name for |residual| (y-axis).
#' @param n_boot  Number of bootstrap replicates.
#' @param seed    Random seed for reproducibility.
#' @return List with: model (nls fit), fit_df (prediction grid with CI),
#'         params (a, b, half_distance with CIs).
bootstrap_exp_decay <- function(data,
                                x_col,
                                y_col,
                                n_boot = 1000,
                                seed = 11) {
  x <- data[[x_col]]
  y <- abs(data[[y_col]])

  # Fit on full data
  fit_data <- data.frame(x = x, y = y)
  model <- minpack.lm::nlsLM(
    y ~ a * exp(-b * x),
    data = fit_data,
    start = list(a = 1, b = 0.1)
  )

  # Prediction grid
  x_vals <- seq(min(x), max(x), length.out = 200)

  # Bootstrap
  set.seed(seed)
  boot_params <- replicate(n_boot, {
    samp <- fit_data[sample(nrow(fit_data), replace = TRUE), ]
    fit <- try(
      minpack.lm::nlsLM(
        y ~ a * exp(-b * x),
        data = samp,
        start = list(a = 1, b = 0.1)
      ),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) c(NA, NA) else coef(fit)
  })
  boot_params <- t(boot_params)
  boot_params <- boot_params[complete.cases(boot_params), ]

  boot_preds <- apply(boot_params, 1, function(p) {
    p[1] * exp(-p[2] * x_vals)
  })

  fit_df <- data.frame(
    x = x_vals,
    y_pred = predict(model, newdata = data.frame(x = x_vals)),
    lower = apply(boot_preds, 1, quantile, probs = 0.025),
    upper = apply(boot_preds, 1, quantile, probs = 0.975)
  )
  names(fit_df)[1:2] <- c(x_col, y_col)

  # Parameter summary
  coefs <- coef(model)
  se <- summary(model)$coefficients[, "Std. Error"]
  half_d <- log(2) / coefs["b"]
  half_d_ci <- quantile(
    log(2) / boot_params[, "b"],
    probs = c(0.025, 0.975)
  )

  list(
    model = model,
    fit_df = fit_df,
    params = list(
      a = c(
        est = coefs["a"],
        lower = coefs["a"] - 1.96 * se["a"],
        upper = coefs["a"] + 1.96 * se["a"]
      ),
      b = c(
        est = coefs["b"],
        lower = coefs["b"] - 1.96 * se["b"],
        upper = coefs["b"] + 1.96 * se["b"]
      ),
      half_distance = c(est = half_d, half_d_ci)
    )
  )
}


# ── Plotting Helpers ─────────────────────────────────────────────────────────

#' Add marginal density plots to a ggplot using project defaults.
#'
#' Wraps ggExtra::ggMarginal with the standard density styling used
#' across all scatter plots in this project.
#'
#' @param p A ggplot object.
#' @return A ggplot with marginal density panels.
add_marginal_density <- function(p) {
  ggExtra::ggMarginal(
    p,
    type = "density",
    margins = "both",
    groupColour = FALSE,
    groupFill = FALSE,
    size = 10,
    colour = "grey",
    fill = "lightgrey"
  )
}


#' Scatter plot with LOESS fit line and residual color gradient.
#'
#' Standard layout used for ddGf vs pathogenicity score plots.
#'
#' @param df          Data frame.
#' @param x_col       Column for x-axis (e.g. "ddG_pred").
#' @param y_col       Column for y-axis (e.g. "function_score" or "ESM1v").
#' @param residual_col Column for color gradient (residuals).
#' @param fit_line_df Data frame with LOESS fit line (x_col, y_col columns).
#' @param title       Plot title.
#' @param x_lab       X-axis label.
#' @param y_lab       Y-axis label.
#' @param color_lab   Color legend label.
#' @param color_limits 2-element vector for color scale limits (optional).
#' @param xlim_vals   X-axis limits (optional).
#' @param ylim_vals   Y-axis limits (optional).
#' @param add_marginal Whether to add marginal density (default TRUE).
#' @return A ggplot (with marginals if requested).
plot_scatter_with_loess <- function(df,
                                    x_col,
                                    y_col,
                                    residual_col,
                                    fit_line_df,
                                    title = "",
                                    x_lab = "Predicted ddGf",
                                    y_lab = "Experimental score",
                                    color_lab = "LOESS residuals",
                                    color_limits = NULL,
                                    xlim_vals = NULL,
                                    ylim_vals = NULL,
                                    add_marginal = TRUE) {
  if (is.null(color_limits)) {
    lim <- max(abs(df[[residual_col]]), na.rm = TRUE)
    color_limits <- c(-lim, lim)
  }

  p <- ggplot(df, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    color = .data[[residual_col]]
  )) +
    geom_point(size = 2, alpha = 0.35) +
    geom_vline(
      xintercept = 0, linetype = "dashed",
      linewidth = 0.5, color = "grey"
    ) +
    geom_hline(
      yintercept = 0, linetype = "dashed",
      linewidth = 0.5, color = "grey"
    ) +
    geom_line(
      data = fit_line_df,
      aes(x = .data[[x_col]], y = .data[[y_col]]),
      inherit.aes = FALSE, color = "black", linewidth = 1
    ) +
    scale_color_gradient2(
      low = "red", mid = "grey", high = "blue",
      midpoint = 0, limits = color_limits
    ) +
    labs(title = title, x = x_lab, y = y_lab, color = color_lab) +
    theme_classic() +
    theme(legend.position = "none")

  if (!is.null(xlim_vals)) p <- p + xlim(xlim_vals)
  if (!is.null(ylim_vals)) p <- p + ylim(ylim_vals)

  if (add_marginal) add_marginal_density(p) else p
}


#' Violin plot with median summary and significance annotations.
#'
#' Standard layout for comparing residual distributions across groups.
#'
#' @param df          Data frame.
#' @param group_col   Column for x-axis grouping.
#' @param value_col   Column for y-axis values.
#' @param colors      Named vector of fill colors.
#' @param title       Plot title.
#' @param y_lab       Y-axis label.
#' @param comparisons List of 2-element vectors for signif bars (optional).
#' @return A ggplot object.
plot_violin_with_stats <- function(df,
                                   group_col,
                                   value_col,
                                   colors,
                                   title = "",
                                   y_lab = "|Residual|",
                                   comparisons = NULL) {
  label_df <- df %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarise(
      n = dplyr::n(),
      median_val = median(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )

  y_max <- max(abs(df[[value_col]]), na.rm = TRUE)

  p <- ggplot(df, aes(
    x = .data[[group_col]],
    y = .data[[value_col]],
    fill = .data[[group_col]]
  )) +
    geom_violin(
      trim = FALSE, scale = "width",
      alpha = 0.8, color = NA
    ) +
    geom_jitter(
      width = 0.15, size = 2,
      alpha = 0.7, color = "lightgrey"
    ) +
    stat_summary(
      fun = median, geom = "crossbar",
      width = 0.4, color = "black", fatten = 1
    ) +
    stat_summary(
      fun = median, geom = "point",
      shape = 23, size = 2, fill = "black",
      color = "black", stroke = 0.7
    ) +
    geom_text(
      data = label_df,
      aes(
        x = .data[[group_col]],
        y = y_max * 1.1,
        label = paste0("n = ", n)
      ),
      inherit.aes = FALSE, size = 4
    ) +
    geom_text(
      data = label_df,
      aes(
        x = .data[[group_col]],
        y = median_val + 0.5,
        label = sprintf("%.2f", median_val)
      ),
      inherit.aes = FALSE, size = 5
    ) +
    scale_fill_manual(values = colors) +
    labs(title = title, x = "", y = y_lab) +
    theme_classic() +
    theme(legend.position = "none")

  if (!is.null(comparisons)) {
    p <- p + ggsignif::geom_signif(
      comparisons = comparisons,
      map_signif_level = FALSE,
      test = "wilcox.test",
      step_increase = 0.1,
      tip_length = 0.01
    )
  }

  p
}
