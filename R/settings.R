#' Set Analysis Settings
#'
#' Configures various parameters for the niche analysis workflow.
#'
#' @param target_species Character vector. Optional. A vector of full scientific names
#'        (e.g., "Genus_species") to be treated as target species.
#'        If NULL, `target_genera` must be provided.
#' @param target_genera Character vector. Optional. A vector of genus names. All species
#'        belonging to these genera found in the prepared data will be treated as
#'        target species. If NULL, `target_species` must be provided.
#' @param exclude_congenerics_if_species_target Logical. If `target_species` is used and
#'        this is `TRUE`, other species from the same genus as a specified target species
#'        will be excluded from the set of "comparison" species for that target.
#'        Default is `TRUE`.
#'
#' @param analyses_to_run Character vector. Specifies which analysis modules to execute.
#'        Possible values include: "overlap_indices", "equivalency_test",
#'        "permanova_pairwise", "lda_multigroup", "modularity_network", "hypervolume_analysis",
#'        "variable_stability", "similarity_validation", "categorical_association", "enfa",
#'        "dendrogram_global_hv". Default runs a comprehensive set.
#'
#' @param variable_selection_method Character. Method for selecting environmental variables.
#'        Options: "Boruta", "RF_Importance", "Manual". Default is "Boruta".
#' @param rf_importance_n_top Integer. If `variable_selection_method` is "RF_Importance",
#'        this specifies the number of top variables to select. Default is 8.
#' @param boruta_max_confirmed_vars Integer. If `variable_selection_method` is "Boruta",
#'        this limits the number of confirmed variables. This limit is applied *after*
#'        the main Boruta run (and after stability analysis, if performed). Default is 6.
#' @param boruta_internal_max_runs Integer. The `maxRuns` parameter for each individual
#'        Boruta algorithm execution (main run). Default is 100.
#' @param boruta_max_runs_for_stability Integer. The `maxRuns` parameter for each Boruta
#'        execution *within* the bootstrap stability analysis. Can be lower than
#'        `boruta_internal_max_runs` for speed. Default is 50.
#' @param manual_selected_vars Named list. Optional. If `variable_selection_method` is "Manual",
#'        this list should contain entries for each target species, where each entry is a
#'        character vector of variable names to use.
#'        Example: `list(Tuber_aestivum = c("bio1", "bio12"), Tuber_borchii = c("aspect", "soil_pH"))`.
#'
#' @param generate_plots Logical. Should plots be generated? Default is `TRUE`.
#' @param plot_type Character. For certain plots (e.g., distributions, PCA, LDA),
#'        choose between "all_comparison_species" or "n_closest_species".
#'        Default is "n_closest_species".
#' @param plot_n_closest Integer. If `plot_type` is `"n_closest_species"`, this specifies
#'        how many of the most similar species to include in these plots. Default is 10.
#' @param plot_n_closest_overlap_metric Character. Metric used to determine the "closest"
#'        species for subset plots. Valid options include "composite_overlap",
#'        "mean_schoener_D", "mean_pianka_O", "mean_czek_psi", "mean_hurlbert_L".
#'        Default is "mean_schoener_D".
#' @param generate_report Logical. Should a summary report (e.g., Markdown) be generated?
#'        Default is `TRUE`.
#' @param continuous_dist_plot_type Character. Type of plot for continuous variable distributions.
#'        Options: "violin" (default), "ridgeline".
#' @param plot_unified_distributions Logical. If TRUE, generates unified distribution plots
#'        (across all species in the dataset) for variable types specified in
#'        `plot_unified_distributions_types`. Default is `FALSE`.
#' @param plot_unified_distributions_types Character vector. If `plot_unified_distributions`
#'        is TRUE, this specifies which types of unified plots to generate.
#'        Currently, "categorical" is supported. Default is `c("categorical")`.
#' @param plot_unified_all_species_distributions Logical. If TRUE, generates a single, unified
#'        distribution plot containing *all* species in the dataset, after the main
#'        analysis loop completes. Default is `TRUE`.
#' @param plot_categorical_distributions Logical. If TRUE, generates plots for categorical
#'        variable distributions for each species pair.
#'        Default is `TRUE`.
#'
#' @param filter_by_shared_categorical_classes Logical. If `TRUE`, comparison species will be
#'        further filtered to include only those that share at least one class with the
#'        target species in at least one of the variables specified in
#'        `categorical_vars_for_class_sharing_filter`. Default is `FALSE`.
#' @param categorical_vars_for_class_sharing_filter Character vector. Names of categorical
#'        variables to use for the shared class filter.
#' @param correlation_threshold Numeric. Absolute correlation threshold for identifying
#'        highly correlated variables during pre-analysis. Default is 0.8.
#' @param vif_threshold Numeric. VIF threshold for identifying multicollinearity. Default is 5.
#' @param min_final_variables_for_multivariate Numeric. Minimum number of variables required
#'        after selection to proceed with multivariate analyses (PCA, LDA). Default is 4.
#' @param boruta_perform_stability_analysis Logical. If TRUE and `variable_selection_method` is "Boruta",
#'        a stability analysis will be performed using bootstrapping, and its results will
#'        inform the final variable selection. Default is `TRUE`.
#' @param boruta_stability_freq_threshold Numeric. If `boruta_perform_stability_analysis` is TRUE,
#'        this is the minimum frequency (0-1) a variable must be confirmed in bootstrap runs to be retained. Default is 0.7.
#' @param equiv_dist_similar_threshold Numeric. Centroid distance threshold for equivalency test
#'        to consider species "potentially similar" even if p-value is significant. Default is 2.0.
#' @param modularity_network_overlap_threshold Numeric. Overlap threshold for building
#'        the modularity network. Default is 0.25.
#' @param n_bootstrap_stability Integer. Number of bootstrap iterations for variable
#'        stability analysis (if "Boruta" is used). Default is 50.
#' @param permanova_use_gower Logical. If TRUE, uses Gower distance (suitable for
#'        mixed continuous and categorical variables) for PERMANOVA. If FALSE,
#'        uses Euclidean distance on standardized numeric variables only. Default is `TRUE`.
#' @param validation_cv_folds Integer. Number of cross-validation folds for similarity
#'        validation. Default is 5.
#' @param similarity_validation_p_threshold Numeric. The p-value threshold (from PERMANOVA
#'        or Equivalency Test) used to identify "statistically similar" species for
#'        the similarity validation module. Default is 0.05.
#' @param similarity_validation_source Character vector. Specifies the source(s) for identifying
#'        "similar" species for the validation module. Can be a combination of:
#'        "permanova", "equivalency", "overlap_indices", "hypervolume".
#'        Default is `c("permanova", "equivalency")`.
#' @param similarity_validation_overlap_threshold Numeric. If `similarity_validation_source`
#'        includes "overlap_indices" or "hypervolume", this is the minimum overlap value
#'        (e.g., Schoener's D or Sorensen) for two species to be considered "similar".
#'        Default is 0.6.
#' @param overlap_n_bins Integer. Number of bins to use for discretizing
#'        continuous variables when calculating overlap indices. Default is 10.
#' @param hypervolume_n_vars Integer. Number of top (e.g., from Boruta/RF) or least
#'        correlated numeric variables to use for hypervolume construction. Default is 3.
#' @param hypervolume_method Character. Method for hypervolume construction.
#'        Options: "gaussian", "box", "svm". Default is "gaussian".
#' @param hypervolume_samples_per_point Integer. Number of Monte Carlo samples per
#'        data point for Gaussian KDE. Default is 100.
#' @param hypervolume_plot_overlap_metric Character. Overlap metric to use for selecting
#'        the closest species for hypervolume plots (`_hv_pca_topN`, `_hv_orig_topN`).
#'        Options: "sorensen", "jaccard". Default is "sorensen".
#' @param hypervolume_n_vars_plot Integer. Number of variables to use for hypervolume
#'        plots (e.g., 2D or 3D projections). Default is 3.
#' @param dendrogram_global_hv_metric Character. Overlap metric to use for constructing the
#'        global hypervolume dendrogram. Options: "jaccard", "sorensen". Default is "jaccard".
#' @param dendrogram_global_hv_cluster_method Character. Clustering method for `hclust`
#'        for the global hypervolume dendrogram. Default is "average".
#' @param dendrogram_global_hv_n_vars Integer. Number of numeric variables to use for
#'        constructing the global hypervolume dendrogram. Default is 5.
#'
#' @param enfa_n_axes Integer. Number of ENFA axes to retain and use for plotting.
#'        If NULL, all significant axes (or a default based on data) will be used. Default is 2.
#' @param enfa_plot_type Character. Type of plot for ENFA results. Options: "biplot", "scatter".
#'        "biplot" shows both species scores and environmental vectors. "scatter" shows only species scores.
#'        Default is "biplot".
#' @param enfa_plot_species_labels Logical. If TRUE, labels species points in ENFA plots. Default is TRUE.
#' @param enfa_plot_env_vectors Logical. If TRUE, plots environmental vectors in ENFA biplots. Default is TRUE.
#' @param enfa_plot_env_labels Logical. If TRUE, labels environmental vectors in ENFA biplots. Default is TRUE.
#' @param enfa_marginality_threshold Numeric. A threshold for marginality to highlight species
#'        with high marginality in plots or summaries. If NULL, no specific threshold is applied.
#'        Default is NULL.
#' @param enfa_specialization_threshold Numeric. A threshold for specialization to highlight species
#'        with high specialization in plots or summaries. If NULL, no specific threshold is applied.
#'        Default is NULL.
#'
#' @return A list object of class `niche_analysis_settings` containing all
#'         configured parameters.
#' @export
#'
#' @examples
#' \dontrun{
#'   # Example: Analyze Tuber aestivum and all species of the Quercus genus
#'   # - Variable selection: Boruta
#'   # - Plots: N closest species, ridgeline for continuous, categorical distributions enabled
#'   # - Various technical parameters customized
#'   analysis_config_boruta <- set_analysis_settings(
#'     target_species = "Tuber_aestivum",
#'     target_genera = "Quercus",
#'     exclude_congenerics_if_species_target = TRUE,
#'     analyses_to_run = c("overlap_indices", "permanova_pairwise", "lda_multigroup", "enfa",
#'                         "modularity_network", "hypervolume_analysis", "dendrogram_global_hv"),
#'     variable_selection_method = "Boruta",
#'     boruta_max_confirmed_vars = 5,
#'     generate_plots = TRUE,
#'     plot_type = "n_closest_species",
#'     plot_n_closest = 5,
#'     plot_n_closest_overlap_metric = "mean_schoener_D",
#'     generate_report = TRUE,
#'     continuous_dist_plot_type = "ridgeline",
#'     plot_categorical_distributions = TRUE,
#'     plot_unified_all_species_distributions = TRUE,
#'     filter_by_shared_categorical_classes = FALSE, # Example: not filtering by shared classes
#'     correlation_threshold = 0.75,                # Custom correlation threshold
#'     vif_threshold = 4,                           # Custom VIF threshold
#'     min_final_variables_for_multivariate = 3,    # Custom minimum variables
#'     permanova_use_gower = TRUE                   # Use Gower for PERMANOVA
#'   )
#'
#'   # Example 2: Manual variable selection for specific targets
#'   # - Focuses on overlap indices.
#'   # - Filters comparison species by shared categorical classes.
#'   # - Uses violin plots for continuous distributions.
#'    manual_vars <- list(
#'     Tuber_melanosporum = c("bio1majority", "bio12majority", "lucas_pH_CaClmajority"),
#'     Boletus_edulis = c("aspectmean", "tree_cover", "Soil")
#'   )
#'   analysis_config_manual <- set_analysis_settings(
#'     target_species = c("Tuber_melanosporum", "Boletus_edulis"),
#'    analyses_to_run = c("overlap_indices", "categorical_association"),
#'     variable_selection_method = "Manual",
#'     manual_selected_vars = manual_vars,
#'     filter_by_shared_categorical_classes = TRUE,
#'     categorical_vars_for_class_sharing_filter = c("Soil", "Vegetation"), # Example filter vars
#'     continuous_dist_plot_type = "violin",
#'     plot_categorical_distributions = TRUE
#'   )
#'
#'   # Example 3: Customizing ENFA parameters
#'   # - Activates ENFA and sets specific plotting options and thresholds.
#'   analysis_config_enfa_custom <- set_analysis_settings(
#'     target_species = "Tuber_aestivum",
#'     analyses_to_run = c("enfa"),
#'     generate_plots = TRUE,
#'     enfa_n_axes = 3, # Use 3 ENFA axes
#'     enfa_plot_type = "biplot",
#'     enfa_plot_species_labels = FALSE, # Do not label individual points
#'     enfa_plot_env_vectors = TRUE,
#'     enfa_plot_env_labels = TRUE,
#'     enfa_marginality_threshold = 0.5, # Highlight if marginality > 0.5
#'     enfa_specialization_threshold = 0.7 # Highlight if specialization > 0.7
#'   )
#' }


set_analysis_settings <- function(
    target_species = NULL,
    target_genera = NULL,
    exclude_congenerics_if_species_target = TRUE,
    analyses_to_run = c("overlap_indices", "equivalency_test", "permanova_pairwise",
                        "lda_multigroup", "modularity_network", "hypervolume_analysis", "dendrogram_global_hv",
                        "variable_stability", "similarity_validation", "categorical_association", "enfa"),
    variable_selection_method = "Boruta",
    rf_importance_n_top = 8,
    boruta_max_confirmed_vars = 6,
    categorical_vars_for_class_sharing_filter = NULL,
    boruta_internal_max_runs = 100,
    boruta_max_runs_for_stability = 50,
    boruta_perform_stability_analysis = TRUE,
    boruta_stability_freq_threshold = 0.7,
    manual_selected_vars = NULL,
    generate_plots = TRUE,
    plot_type = "n_closest_species",
    plot_n_closest = 10,
    plot_n_closest_overlap_metric = "mean_schoener_D",
    generate_report = TRUE,
    continuous_dist_plot_type = "violin",
    plot_categorical_distributions = TRUE,
    plot_unified_distributions = FALSE, 
    plot_unified_distributions_types = c("categorical"),
    plot_unified_all_species_distributions = TRUE,
    filter_by_shared_categorical_classes = FALSE,
    correlation_threshold = 0.8,
    vif_threshold = 5,
    min_final_variables_for_multivariate = 4,
    equiv_dist_similar_threshold = 2.0,
    modularity_network_metric = "mean_schoener_D",
    modularity_network_overlap_threshold = 0.25,
    n_bootstrap_stability = 50,
    permanova_use_gower = TRUE,
    similarity_validation_p_threshold = 0.05,
    similarity_validation_source = c("permanova", "equivalency"),
    similarity_validation_overlap_threshold = 0.6,
    validation_cv_folds = 5,
    overlap_n_bins = 10,
    hypervolume_n_vars = 3,
    hypervolume_method = "gaussian",
    hypervolume_samples_per_point = 100,
    hypervolume_plot_overlap_metric = "sorensen",
    hypervolume_n_vars_plot = 3,
    dendrogram_global_hv_metric = "sorensen",
    dendrogram_target_specific_hv_metric = "sorensen",
    dendrogram_global_hv_cluster_method = "average",
    dendrogram_global_hv_n_vars = 5,
    # ENFA Parameters
    enfa_n_axes = 2,
    enfa_plot_type = "biplot", # Default to "biplot"
    enfa_plot_species_labels = TRUE, # Default to TRUE
    enfa_plot_env_vectors = TRUE, # Default to TRUE
    enfa_plot_env_labels = TRUE, # Default to TRUE
    enfa_marginality_threshold = NULL,
    enfa_specialization_threshold = NULL
) {
  # --- 1. Input Validation ---
  if (is.null(target_species) && is.null(target_genera)) {
    stop("Either 'target_species' or 'target_genera' must be provided.")
  }
  if (!is.null(target_species) && !is.character(target_species)) {
    stop("'target_species' must be a character vector or NULL.")
  }
  if (!is.null(target_genera) && !is.character(target_genera)) {
    stop("'target_genera' must be a character vector or NULL.")
  }
  if (!is.logical(exclude_congenerics_if_species_target) || length(exclude_congenerics_if_species_target) != 1) {
    stop("'exclude_congenerics_if_species_target' must be a single logical value.")
  }
  valid_analyses <- c("overlap_indices", "equivalency_test", "permanova_pairwise",
                      "lda_multigroup", "modularity_network", "hypervolume_analysis",
                      "dendrogram_global_hv", "variable_stability", "similarity_validation",
                      "categorical_association", "enfa")
  if (!all(analyses_to_run %in% valid_analyses)) {
    stop("Invalid 'analyses_to_run'. Possible values are: ", paste(valid_analyses, collapse = ", "))
  }
  valid_var_selection <- c("Boruta", "RF_Importance", "Manual")
  if (!variable_selection_method %in% valid_var_selection) {
    stop("Invalid 'variable_selection_method'. Options are: ", paste(valid_var_selection, collapse = ", "))
  }
  if (variable_selection_method == "Manual" && is.null(manual_selected_vars)) {
    warning("Variable selection method is 'Manual' but 'manual_selected_vars' is NULL. This might lead to errors later if not handled.")
  }
  if (!is.null(manual_selected_vars) && !is.list(manual_selected_vars)) {
    stop("'manual_selected_vars' must be a named list if provided.")
  }
  if (!is.logical(generate_plots) || length(generate_plots) != 1) {
    stop("'generate_plots' must be a single logical value.")
  }
  valid_plot_types <- c("all_comparison_species", "n_closest_species")
  if (!plot_type %in% valid_plot_types) {
    stop("Invalid 'plot_type'. Options are: ", paste(valid_plot_types, collapse = ", "))
  }
  if (plot_type == "n_closest_species" && (is.null(plot_n_closest) || !is.numeric(plot_n_closest) || plot_n_closest < 1)) {
    stop("'plot_n_closest' must be a positive integer when plot_type is 'n_closest_species'.")
  }
  if (!is.numeric(correlation_threshold) || length(correlation_threshold) != 1 || correlation_threshold < 0 || correlation_threshold > 1) {
      stop("'correlation_threshold' must be a single numeric value between 0 and 1.")
  }
  if (!is.numeric(vif_threshold) || length(vif_threshold) != 1 || vif_threshold < 1) {
      stop("'vif_threshold' must be a single numeric value >= 1.")
  }
  if (!is.numeric(min_final_variables_for_multivariate) || length(min_final_variables_for_multivariate) != 1 || min_final_variables_for_multivariate < 0) {
      stop("'min_final_variables_for_multivariate' must be a single non-negative integer.")
  }
  if (!is.numeric(equiv_dist_similar_threshold) || length(equiv_dist_similar_threshold) != 1 || equiv_dist_similar_threshold < 0) {
      stop("'equiv_dist_similar_threshold' must be a single non-negative numeric value.")
  }
  if (!is.numeric(modularity_network_overlap_threshold) || length(modularity_network_overlap_threshold) != 1 || modularity_network_overlap_threshold < 0 || modularity_network_overlap_threshold > 1) {
      stop("'modularity_network_overlap_threshold' must be a single numeric value between 0 and 1.")
  }
  if (!is.numeric(n_bootstrap_stability) || length(n_bootstrap_stability) != 1 || n_bootstrap_stability < 10) {
      warning("'n_bootstrap_stability' should ideally be >= 100 for robust results. Current value: ", n_bootstrap_stability)
  }
  if (!is.numeric(validation_cv_folds) || length(validation_cv_folds) != 1 || validation_cv_folds < 2) {
      stop("'validation_cv_folds' must be a single integer >= 2.")
  }
  if (!is.numeric(similarity_validation_p_threshold) || length(similarity_validation_p_threshold) != 1 || similarity_validation_p_threshold < 0 || similarity_validation_p_threshold > 1) {
    stop("'similarity_validation_p_threshold' must be a single numeric value between 0 and 1.")
  }
  valid_sim_sources <- c("permanova", "equivalency", "overlap_indices", "hypervolume")
  if (!is.character(similarity_validation_source) || !all(similarity_validation_source %in% valid_sim_sources)) {
    stop("Invalid 'similarity_validation_source'. Options are: ", paste(valid_sim_sources, collapse = ", "))
  }
  if (!is.numeric(similarity_validation_overlap_threshold) || length(similarity_validation_overlap_threshold) != 1 || similarity_validation_overlap_threshold < 0 || similarity_validation_overlap_threshold > 1) {
    stop("'similarity_validation_overlap_threshold' must be a single numeric value between 0 and 1.")
  }
  if (!is.numeric(overlap_n_bins) || length(overlap_n_bins) != 1 || overlap_n_bins < 2) {
      stop("'overlap_n_bins' must be a single integer >= 2.")
  }
  if (!is.numeric(hypervolume_n_vars) || length(hypervolume_n_vars) != 1 || hypervolume_n_vars < 1) {
      stop("'hypervolume_n_vars' must be a single integer >= 1.")
  }
  if (!is.numeric(hypervolume_samples_per_point) || length(hypervolume_samples_per_point) != 1 || hypervolume_samples_per_point < 1) {
      stop("'hypervolume_samples_per_point' must be a single integer >= 1.")
  }
  if (!hypervolume_plot_overlap_metric %in% c("jaccard", "sorensen")) {
    stop("'hypervolume_plot_overlap_metric' must be 'jaccard' or 'sorensen'.")
  }
  if (!is.numeric(hypervolume_n_vars_plot) || length(hypervolume_n_vars_plot) != 1 || hypervolume_n_vars_plot < 0) {
      stop("'hypervolume_n_vars_plot' must be a single non-negative integer.")
  }
  if (!dendrogram_global_hv_metric %in% c("jaccard", "sorensen")) { 
    stop("'dendrogram_global_hv_metric' must be 'jaccard' or 'sorensen'.")
  }
  if (!dendrogram_target_specific_hv_metric %in% c("jaccard", "sorensen")) {
    stop("'dendrogram_target_specific_hv_metric' must be 'jaccard' or 'sorensen'.")
  }
  if (!is.numeric(dendrogram_global_hv_n_vars) || length(dendrogram_global_hv_n_vars) != 1 || dendrogram_global_hv_n_vars < 2) {
    stop("'dendrogram_global_hv_n_vars' must be a single integer >= 2.")
  }
  if (!is.null(enfa_n_axes) && (!is.numeric(enfa_n_axes) || length(enfa_n_axes) != 1 || enfa_n_axes < 1)) {
    stop("'enfa_n_axes' must be a single positive integer or NULL.")
  }
  valid_enfa_plot_types <- c("biplot", "scatter")
  if (!enfa_plot_type %in% valid_enfa_plot_types) {
    stop("Invalid 'enfa_plot_type'. Options are: ", paste(valid_enfa_plot_types, collapse = ", "))
  }
  if (!is.logical(enfa_plot_species_labels) || length(enfa_plot_species_labels) != 1) {
    stop("'enfa_plot_species_labels' must be a single logical value.")
  }
  if (!is.logical(enfa_plot_env_vectors) || length(enfa_plot_env_vectors) != 1) {
    stop("'enfa_plot_env_vectors' must be a single logical value.")
  }
  if (!is.logical(enfa_plot_env_labels) || length(enfa_plot_env_labels) != 1) {
    stop("'enfa_plot_env_labels' must be a single logical value.")
  }
  if (!is.null(enfa_marginality_threshold) && (!is.numeric(enfa_marginality_threshold) || length(enfa_marginality_threshold) != 1 || enfa_marginality_threshold < 0)) {
    stop("'enfa_marginality_threshold' must be a single non-negative numeric value or NULL.")
  }
  if (!is.null(enfa_specialization_threshold) && (!is.numeric(enfa_specialization_threshold) || length(enfa_specialization_threshold) != 1 || enfa_specialization_threshold < 0)) {
    stop("'enfa_specialization_threshold' must be a single non-negative numeric value or NULL.")
  }
  # --- 2. Construct Settings List ---
  settings <- list(
    # Target Specification
    target_species_input = target_species,
    target_genera_input = target_genera,
    exclude_congenerics_if_species_target = exclude_congenerics_if_species_target,
    # Analyses Configuration
    analyses_to_run = unique(analyses_to_run),
    # Variable Selection Configuration
    variable_selection_method = variable_selection_method,
    rf_importance_n_top = as.integer(rf_importance_n_top),
    boruta_max_confirmed_vars = as.integer(boruta_max_confirmed_vars),
    boruta_internal_max_runs = as.integer(boruta_internal_max_runs),
    boruta_max_runs_for_stability = as.integer(boruta_max_runs_for_stability),
    boruta_perform_stability_analysis = as.logical(boruta_perform_stability_analysis),
    boruta_stability_freq_threshold = as.numeric(boruta_stability_freq_threshold),
    manual_selected_vars = manual_selected_vars,
    # Plotting Configuration (General)
    generate_plots = generate_plots,
    plot_type = plot_type,
    plot_n_closest = if (plot_type == "n_closest_species") as.integer(plot_n_closest) else NULL,
    plot_n_closest_overlap_metric = plot_n_closest_overlap_metric,
    # Report Configuration (General)
    generate_report = generate_report,
    continuous_dist_plot_type = continuous_dist_plot_type,
    plot_unified_distributions = as.logical(plot_unified_distributions),
    plot_unified_all_species_distributions = as.logical(plot_unified_all_species_distributions),
    plot_unified_distributions_types = if(!is.null(plot_unified_distributions_types)) as.character(plot_unified_distributions_types) else NULL,
    plot_categorical_distributions = plot_categorical_distributions,
    filter_by_shared_categorical_classes = as.logical(filter_by_shared_categorical_classes),
    categorical_vars_for_class_sharing_filter = if(!is.null(categorical_vars_for_class_sharing_filter)) as.character(categorical_vars_for_class_sharing_filter) else NULL,
    # Technical Parameters for Analyses
    correlation_threshold = correlation_threshold,
    vif_threshold = vif_threshold,
    min_final_variables_for_multivariate = as.integer(min_final_variables_for_multivariate),
    equiv_dist_similar_threshold = equiv_dist_similar_threshold,
    modularity_network_metric = modularity_network_metric,
    modularity_network_overlap_threshold = modularity_network_overlap_threshold,
    n_bootstrap_stability = as.integer(n_bootstrap_stability),
    validation_cv_folds = as.integer(validation_cv_folds),
    similarity_validation_p_threshold = as.numeric(similarity_validation_p_threshold),
    similarity_validation_source = similarity_validation_source,
    similarity_validation_overlap_threshold = as.numeric(similarity_validation_overlap_threshold),
    permanova_use_gower = as.logical(permanova_use_gower),
    overlap_n_bins = as.integer(overlap_n_bins),
    hypervolume_n_vars = as.integer(hypervolume_n_vars),
    hypervolume_method = hypervolume_method,
    hypervolume_samples_per_point = as.integer(hypervolume_samples_per_point),
    hypervolume_plot_overlap_metric = hypervolume_plot_overlap_metric,
    hypervolume_n_vars_plot = as.integer(hypervolume_n_vars_plot),
    # Dendrogram Global Hypervolume Settings
    dendrogram_global_hv_metric = dendrogram_global_hv_metric,
    dendrogram_target_specific_hv_metric = dendrogram_target_specific_hv_metric,
    dendrogram_global_hv_cluster_method = dendrogram_global_hv_cluster_method,
    dendrogram_global_hv_n_vars = as.integer(dendrogram_global_hv_n_vars),
    # ENFA Settings
    enfa_n_axes = if(!is.null(enfa_n_axes)) as.integer(enfa_n_axes) else NULL,
    enfa_plot_type = enfa_plot_type,
    enfa_plot_species_labels = as.logical(enfa_plot_species_labels),
    enfa_plot_env_vectors = as.logical(enfa_plot_env_vectors),
    enfa_plot_env_labels = as.logical(enfa_plot_env_labels),
    enfa_marginality_threshold = enfa_marginality_threshold,
    enfa_specialization_threshold = enfa_specialization_threshold
  )
  class(settings) <- "niche_analysis_settings"
  message("Analysis settings configured.")
  return(settings)
}
