#' Plot Hypervolume Overlap Dendrogram (Internal)
#'
#' Generates and saves a dendrogram plot from an hclust object derived from
#' hypervolume overlap.
#'
#' @param hclust_object An object of class `hclust`.
#' @param output_path Character. Full path to save the PNG image.
#' @param title Character. Main title for the plot.
#' @param target_species_highlight Character. Optional. Species ID to highlight in the dendrogram.
#' @param global_species_palette Named character vector. Palette for coloring dendrogram labels.
#' @param verbose Logical.
#'
#' @return The `output_path` if plotting was successful, otherwise NULL.
#' @importFrom graphics par plot mtext
#' @importFrom stats as.dendrogram
#' @importFrom dendextend set
#' @noRd
.plot_hypervolume_dendrogram <- function(hclust_object,
                                         output_path,
                                         title = "Niche Overlap Dendrogram",
                                         target_species_highlight = NULL,
                                         global_species_palette = NULL,
                                         verbose = TRUE) {
  if (verbose) message("      Generating dendrogram plot: ", title)
  dend <- stats::as.dendrogram(hclust_object)
  if (!is.null(target_species_highlight) && !is.null(global_species_palette)) {
    label_colors <- ifelse(labels(dend) == target_species_highlight,
                           global_species_palette[target_species_highlight] %||% "red",
                           "black")
    dend <- dendextend::set(dend, "labels_col", label_colors)
  } else if (!is.null(global_species_palette) && length(global_species_palette) >= length(labels(dend))) {
     label_colors_all <- global_species_palette[labels(dend)]
     if(all(!is.na(label_colors_all))) {
        dend <- dendextend::set(dend, "labels_col", label_colors_all)
     }
  }
  dend <- dendextend::set(dend, "labels_cex", 0.7)
  dend <- dendextend::set(dend, "branches_lwd", 1.2)

  grDevices::png(output_path, width = max(1200, length(labels(dend)) * 40), height = 800, res = 150)
  old_par <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(old_par))
  graphics::par(mar = c(8, 4, 4, 2) + 0.1)
  graphics::plot(dend, main = title, ylab = "Height (Dissimilarity)", horiz = FALSE, hang = -1, cex.axis = 0.8)
  grDevices::dev.off()

  return(output_path)
}

#' Plot Pairwise PERMANOVA R2 Results (Internal)
#'
#' Generates a bar plot of R2 values from pairwise PERMANOVA, colored by significance.
#'
#' @param permanova_results_df Tibble. Output from `.run_permanova_pairwise_module`.
#' @param target_id Character. The ID of the target species.
#' @param output_dir Character. Directory to save the plot.
#' @param plot_suffix Character. Suffix for the plot filename.
#' @param verbose Logical.
#' @noRd
.plot_permanova_pairwise_r2 <- function(permanova_results_df, target_id, output_dir, plot_suffix = "", verbose = TRUE) {
  if (is.null(permanova_results_df) || nrow(permanova_results_df) == 0) {
    if (verbose) message("      PERMANOVA R2 Plot: No PERMANOVA results to plot.")
    return()
  }
  plot_data <- permanova_results_df %>%
    dplyr::filter(specie1 == target_id) %>%
    dplyr::filter(!is.na(R2) & !is.na(p_adj)) %>%
    dplyr::mutate(
      significance_level = dplyr::case_when(
        p_adj < 0.001 ~ "*** (p < 0.001)",
        p_adj < 0.01  ~ "** (p < 0.01)",
        p_adj < 0.05  ~ "* (p < 0.05)",
        TRUE          ~ "ns (p >= 0.05)"
      ),
      significance_level = factor(significance_level, levels = c("*** (p < 0.001)", "** (p < 0.01)", "* (p < 0.05)", "ns (p >= 0.05)"))
    ) %>%
    dplyr::arrange(dplyr::desc(R2))
  if (nrow(plot_data) == 0) {
    if (verbose) message("      PERMANOVA R2 Plot: No valid data after filtering for target '", target_id, "'.")
    return()
  }
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = R2, y = stats::reorder(specie2, R2), fill = significance_level)) +
    ggplot2::geom_col(alpha = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(round(R2, 2))), hjust = -0.2, size = 3) +
    ggplot2::scale_fill_brewer(palette = "RdYlBu", direction = -1, name = "Significance (adj. p-value)") +
    ggplot2::labs(
      title = paste("Pairwise PERMANOVA: R-squared Values vs.", target_id, plot_suffix),
      subtitle = "Higher R2 indicates greater dissimilarity explained by species identity",
      x = bquote(R^2 ~ "Value (Effect Size)"),
      y = "Comparison Species"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = "bottom",
                   axis.text.y = ggplot2::element_text(size = 8)) +
    ggplot2::xlim(0, NA)
  plot_filename <- file.path(output_dir, paste0("permanova_pairwise_R2", plot_suffix, "_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".png"))
  ggplot2::ggsave(filename = plot_filename, plot = p, width = 8, height = max(4, nrow(plot_data) * 0.3 + 2), dpi = 300, limitsize = FALSE)
  if (verbose) message("      Pairwise PERMANOVA R2 plot saved to: ", plot_filename)
}

#' Plot Pairwise PERMANOVA P-value Results (Internal)
#'
#' Generates a bar plot of -log10(adjusted p-values) from pairwise PERMANOVA.
#'
#' @param permanova_results_df Tibble. Output from `.run_permanova_pairwise_module`.
#' @param target_id Character. The ID of the target species.
#' @param output_dir Character. Directory to save the plot.
#' @param plot_suffix Character. Suffix for the plot filename.
#' @param verbose Logical.
#' @noRd
.plot_permanova_pairwise_pvalue <- function(permanova_results_df, target_id, output_dir, plot_suffix = "", verbose = TRUE) {
  if (is.null(permanova_results_df) || nrow(permanova_results_df) == 0) {
    if (verbose) message("      PERMANOVA P-value Plot: No PERMANOVA results to plot.")
    return()
  }
  plot_data <- permanova_results_df %>%
    dplyr::filter(specie1 == target_id) %>%
    dplyr::filter(!is.na(p_adj)) %>%
    dplyr::mutate(
      log10_p_adj = -log10(p_adj),
      log10_p_adj = ifelse(is.infinite(log10_p_adj), -log10(.Machine$double.eps*0.1), log10_p_adj),
      signif_status = ifelse(p_adj < 0.05, "Significant (p.adj < 0.05)", "Not Significant (p.adj >= 0.05)"),
      label_text = paste0("R2=", round(R2, 2), "\np.adj=", format.pval(p_adj, digits = 2, eps = 0.001))
    ) %>%
    dplyr::arrange(p_adj)

  if (nrow(plot_data) == 0) {
    if (verbose) message("      PERMANOVA P-value Plot: No valid data after filtering for target '", target_id, "'.")
    return()
  }
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = stats::reorder(specie2, p_adj), y = log10_p_adj, fill = signif_status)) +
    ggplot2::geom_col(alpha = 0.8, width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = label_text), vjust = -0.3, size = 2.5, color = "black") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    ggplot2::annotate("text", x = Inf, y = -log10(0.05), label = "p.adj = 0.05", hjust = 1.1, vjust = -0.5, size = 3, color = "red") +
    ggplot2::scale_fill_manual(values = c("Not Significant (p.adj >= 0.05)" = "skyblue3", "Significant (p.adj < 0.05)" = "salmon3"), name = "Significance") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
    ggplot2::labs(
      title = paste("Pairwise PERMANOVA: Significance vs.", target_id, plot_suffix),
      subtitle = "Higher bars indicate stronger statistical significance of dissimilarity",
      x = "Comparison Species",
      y = bquote(-log[10] ~ "(adjusted p-value)")
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 8),
                   legend.position = "bottom")
  plot_filename <- file.path(output_dir, paste0("permanova_pairwise_pvalue", plot_suffix, "_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".png"))
  ggplot2::ggsave(filename = plot_filename, plot = p, width = max(7, nrow(plot_data) * 0.4 + 2), height = 7, dpi = 300, limitsize = FALSE)
  if (verbose) message("      Pairwise PERMANOVA P-value plot saved to: ", plot_filename)
}