#' Resolve Actual Target Species List (Internal)
#'
#' Determines the final list of target species IDs based on user settings
#' (specific species and/or genera) and the available species in the prepared data.
#'
#' @param settings An object of class 'niche_analysis_settings'.
#' @param available_species_ids Character vector. All unique species IDs present
#'        in the `niche_data_object$data`.
#' @param data_frame_for_genus_lookup Tibble/data.frame. The `niche_data_object$data`
#'        containing at least the species ID column and the genus column.
#' @param species_id_col_name Character. Name of the species ID column in `data_frame_for_genus_lookup`.
#' @param genus_col_name_processed Character. Name of the processed genus column
#'        in `data_frame_for_genus_lookup` (e.g., "Genus_processed_internal"
#'        or the name saved in `niche_data_object$metadata$genus_col_name_processed`).
#' @param verbose Logical. If TRUE, prints messages.
#'
#' @return A character vector of unique, valid target species IDs.
#' @noRd
.resolve_target_species_list <- function(settings,
                                         available_species_ids,
                                         data_frame_for_genus_lookup,
                                         species_id_col_name,
                                         genus_col_name_processed,
                                         verbose = TRUE) {
  target_species_from_input <- character(0)
  # 1. Process specific target species names
  if (!is.null(settings$target_species_input) && length(settings$target_species_input) > 0) {
    target_species_from_input <- unique(settings$target_species_input)
    if (verbose) message("  User-specified target species: ", paste(target_species_from_input, collapse = ", "))
  }
  # 2. Process target genera
  species_from_genera <- character(0)
  if (!is.null(settings$target_genera_input) && length(settings$target_genera_input) > 0) {
    if (!genus_col_name_processed %in% names(data_frame_for_genus_lookup)) {
      warning("Genus column '", genus_col_name_processed, "' not found in the provided data for resolving target genera. Skipping genus-based targets.")
    } else {
      if (verbose) message("  Resolving species for target genera: ", paste(settings$target_genera_input, collapse = ", "))
      species_from_genera <- data_frame_for_genus_lookup %>%
        dplyr::filter(!!dplyr::sym(genus_col_name_processed) %in% settings$target_genera_input) %>%
        dplyr::pull(!!dplyr::sym(species_id_col_name)) %>%
        unique()
      if (verbose) message("    Found ", length(species_from_genera), " species from specified genera: ", paste(species_from_genera, collapse = ", "))
    }
  }
  # 3. Combine and filter by availability
  potential_targets <- unique(c(target_species_from_input, species_from_genera))
  if (length(potential_targets) == 0) {
    if (verbose) message("  No potential target species identified from inputs.")
    return(character(0))
  }
  actual_targets <- intersect(potential_targets, available_species_ids)
  if (length(actual_targets) == 0 && verbose) {
    message("  Warning: None of the specified target species or species from target genera were found in the prepared dataset.")
  }
  return(actual_targets)
}
#' Prepare Data for a Single Target Species Analysis (Internal)
#'
#' Subsets the main niche data for a given target species and its comparison species.
#' Optionally excludes congeneric species from the comparison set if specified in settings.
#'
#' @param current_target_id Character. The ID of the current target species.
#' @param niche_data_object An object of class 'niche_data'.
#' @param settings An object of class 'niche_analysis_settings'.
#' @param verbose Logical. If TRUE, prints messages.
#'
#' @return A list containing:
#'   \item{data_subset_for_analysis}{A tibble filtered for the target and its comparison species.}
#'   \item{comparison_species_ids}{A character vector of species IDs used for comparison.}
#'   \item{target_genus}{The genus of the current target species (or NA if not found).}
#'         Returns NULL if data preparation fails (e.g., target not found).
#' @noRd
.prepare_target_comparison_data <- function(current_target_id,
                                            niche_data_object,
                                            settings,
                                            verbose = TRUE) {
  full_data <- niche_data_object$data
  species_col <- niche_data_object$metadata$species_id_col
  genus_col <- niche_data_object$metadata$genus_col_name_processed
  if (!current_target_id %in% full_data[[species_col]]) {
    if (verbose) message("  Error: Current target '", current_target_id, "' not found in the dataset.")
    return(NULL)
  }
  all_other_species <- setdiff(unique(full_data[[species_col]]), current_target_id)
  target_genus_value <- NA_character_
  if (genus_col %in% names(full_data)) {
    target_genus_value <- unique(full_data[[genus_col]][full_data[[species_col]] == current_target_id])
    if(length(target_genus_value) > 1) {
      if(verbose) warning("  Multiple genera found for target '", current_target_id, "'. Using the first one: '", target_genus_value[1], "'.")
      target_genus_value <- target_genus_value[1]
    } else if (length(target_genus_value) == 0) {
      target_genus_value <- NA_character_
    }
  } else {
    if (verbose) warning("  Genus column '", genus_col, "' not found for congeneric exclusion logic.")
  }
  comparison_species_final <- all_other_species
  # Exclude congenerics if:
  # 1. Setting is TRUE
  # 2. The current_target_id was specified directly by the user (not derived from a genus target)
  # 3. The target_genus_value is valid
  if (settings$exclude_congenerics_if_species_target &&
      current_target_id %in% (settings$target_species_input %||% character(0)) &&
      !is.na(target_genus_value) && genus_col %in% names(full_data)) {
    congenerics_to_exclude <- full_data %>%
      dplyr::filter(!!dplyr::sym(genus_col) == target_genus_value & !!dplyr::sym(species_col) != current_target_id) %>%
      dplyr::pull(!!dplyr::sym(species_col)) %>%
      unique()
    if (length(congenerics_to_exclude) > 0 && verbose) {
      message("  Excluding ", length(congenerics_to_exclude), " congeneric species for target '", current_target_id, "': ", paste(congenerics_to_exclude, collapse = ", "))
    }
    comparison_species_final <- setdiff(all_other_species, congenerics_to_exclude)
  }
  if (length(comparison_species_final) == 0) {
    if (verbose) message("  No comparison species remaining for target '", current_target_id, "' after filtering.")
    return(list(
      data_subset_for_analysis = full_data %>% dplyr::filter(!!dplyr::sym(species_col) == current_target_id),
      comparison_species_ids = character(0),
      target_genus = target_genus_value
    ))
  }
  if (settings$filter_by_shared_categorical_classes &&
      !is.null(settings$categorical_vars_for_class_sharing_filter) &&
      length(settings$categorical_vars_for_class_sharing_filter) > 0 &&
      length(comparison_species_final) > 0) {
    vars_for_filter_actual <- intersect(settings$categorical_vars_for_class_sharing_filter,
                                        niche_data_object$metadata$env_vars$categorical)
    if (length(vars_for_filter_actual) == 0) {
      if (verbose) warning("  No valid categorical variables found from 'categorical_vars_for_class_sharing_filter' in the dataset. Skipping shared class filter.")
    } else {
      if (verbose) {
        message("  Applying filter: shared categorical classes using variables: ",
                paste(vars_for_filter_actual, collapse = ", "))
      }
      target_classes_by_var <- list()
      for (cat_var_filter in vars_for_filter_actual) {
        if (cat_var_filter %in% names(full_data)) {
          target_classes_by_var[[cat_var_filter]] <- unique(as.character(full_data[[cat_var_filter]][full_data[[species_col]] == current_target_id]))
          target_classes_by_var[[cat_var_filter]] <- target_classes_by_var[[cat_var_filter]][!is.na(target_classes_by_var[[cat_var_filter]])]
        }
      }
      if (length(unlist(target_classes_by_var)) == 0) {
        if (verbose) warning("  Target species '", current_target_id, "' does not occupy any classes in the specified filter variables (", paste(vars_for_filter_actual, collapse=", "),"). No comparison species will be selected by this filter.")
        comparison_species_final <- character(0)
      } else {
        species_passing_class_filter <- purrr::map_lgl(comparison_species_final, function(comp_sp_id) {
          shares_any_class_with_target <- FALSE
          for (cat_var_filter in names(target_classes_by_var)) { 
            if (cat_var_filter %in% names(full_data) && length(target_classes_by_var[[cat_var_filter]]) > 0) {
              comp_classes_for_var <- unique(as.character(full_data[[cat_var_filter]][full_data[[species_col]] == comp_sp_id]))
              comp_classes_for_var <- comp_classes_for_var[!is.na(comp_classes_for_var)]
              if (length(intersect(target_classes_by_var[[cat_var_filter]], comp_classes_for_var)) > 0) {
                return(TRUE)
              }
            }
          }
          return(FALSE)
        })
        comparison_species_final <- comparison_species_final[species_passing_class_filter]
        if (verbose) message("    ", length(comparison_species_final), " comparison species remaining after shared categorical class filter.")
      }
    }
  }
  data_subset <- full_data %>%
    dplyr::filter(!!dplyr::sym(species_col) %in% c(current_target_id, comparison_species_final))
  return(list(
    data_subset_for_analysis = data_subset,
    comparison_species_ids = comparison_species_final,
    target_genus = target_genus_value
  ))
}
"%||%" <- function(a, b) if (!is.null(a)) a else b

#' Perform Variable Selection (Internal)
#'
#' Orchestrates the variable selection process using the specified method
#' (Boruta, RF Importance, Manual) and applies post-selection filters
#' (correlation, VIF).
#'
#' @param data_for_selection Tibble. The data subset for the current target species
#'        and its comparison species. Must contain the species ID column and
#'        all potential environmental variables.
#' @param target_id Character. The ID of the current target species.
#' @param all_potential_env_vars Character vector. Names of all environmental
#'        variables available in `data_for_selection` before selection.
#' @param categorical_vars Character vector. Names of categorical environmental variables.
#' @param species_col_actual_name Character. The actual name of the species ID column in `data_for_selection`.
#' @param metadata_labels List. Optional. A list containing user-defined labels for variables, typically from `niche_data_object$metadata$labels`.
#' @param settings An object of class 'niche_analysis_settings'.
#' @param output_dir_target_specific Character. Path to the output directory
#'        specific to the current target species.
#' @param verbose Logical. If TRUE, prints progress messages.
#' @noRd
#' @return A list containing:
#'   \item{selected_variables}{Character vector of final selected variable names.}
#'   \item{initial_selection_method}{Character. The method used for initial selection.}
#'   \item{initial_selected_variables}{Character vector. Variables after initial selection.}
#'   \item{details}{A list containing method-specific results (e.g., Boruta stats)
#'                  and filter results (correlation matrix, VIF values).}
.perform_variable_selection <- function(data_for_selection,
                                        target_id,
                                        all_potential_env_vars,
                                        categorical_vars,
                                        species_col_actual_name,
                                        settings,
                                        metadata_labels = NULL, 
                                        output_dir_target_specific,
                                        verbose = TRUE) {
  if (verbose) message("  Starting variable selection for target: ", target_id, " using method: ", settings$variable_selection_method)
  selection_details <- list()
  # --- 1. Initial Variable Selection based on method ---
  # Operate selection ONLY on continuous variables
  continuous_vars_available <- setdiff(all_potential_env_vars, categorical_vars)
  user_specified_categorical_vars_in_data <- intersect(all_potential_env_vars, categorical_vars) 
  if (settings$variable_selection_method == "Boruta") {
    selection_details$stability_analysis_performed = FALSE
    if (verbose) message("    Running Boruta selection module...")
    boruta_output <- .run_boruta_selection_module(
      data_for_boruta = data_for_selection, 
      target_id = target_id,
      potential_vars = continuous_vars_available, 
      species_col_actual_name = species_col_actual_name,
      settings = settings,
      output_dir_boruta = output_dir_target_specific,
      metadata_labels = metadata_labels,
      save_outputs = TRUE, 
      apply_max_confirmed_limit = FALSE,
      verbose = verbose
    )
    initial_selected_continuous_vars <- boruta_output$confirmed_variables
    selection_details$boruta_stats <- boruta_output$boruta_stats_df
    selection_details$boruta_tentative_variables <- boruta_output$tentative_variables
    selection_details$boruta_plot_path <- boruta_output$boruta_plot_path
    if (settings$boruta_perform_stability_analysis && length(initial_selected_continuous_vars) > 0) {
      if (verbose) message("    Performing Boruta stability analysis...")
      selection_details$stability_analysis_performed <- TRUE
      stability_results_df <- .run_variable_stability_module(
        data_for_stability = data_for_selection,
        target_id = target_id,
        potential_vars_for_stability = continuous_vars_available,
        species_col_actual_name = species_col_actual_name, 
        settings = settings, 
        output_dir_module = output_dir_target_specific,
        verbose = verbose
      )
      selection_details$stability_frequencies = stability_results_df
      stable_and_confirmed_vars <- stability_results_df %>%
        dplyr::filter(Selection_Frequency >= settings$boruta_stability_freq_threshold) %>%
        dplyr::pull(Variable)
      initial_selected_continuous_vars <- intersect(initial_selected_continuous_vars, stable_and_confirmed_vars)
      if (verbose) message("      Continuous variables after stability filter (freq >= ", settings$boruta_stability_freq_threshold, "): ", paste(initial_selected_continuous_vars, collapse = ", "))
    }
    if (!is.null(settings$boruta_max_confirmed_vars) && length(initial_selected_continuous_vars) > settings$boruta_max_confirmed_vars) {
      if (!is.null(selection_details$boruta_stats) && "Original_Variable" %in% names(selection_details$boruta_stats) && "meanImp" %in% names(selection_details$boruta_stats)) {
        top_vars_after_stability <- selection_details$boruta_stats %>%
          dplyr::filter(Original_Variable %in% initial_selected_continuous_vars) %>%
          dplyr::slice_max(order_by = meanImp, n = settings$boruta_max_confirmed_vars) %>%
          dplyr::pull(Original_Variable)
        initial_selected_continuous_vars <- top_vars_after_stability
      } else {
        initial_selected_continuous_vars <- head(initial_selected_continuous_vars, settings$boruta_max_confirmed_vars)
      }
      if (verbose) message("      Continuous variables after applying boruta_max_confirmed_vars limit: ", paste(initial_selected_continuous_vars, collapse = ", "))
    }
    initial_selected_vars <- initial_selected_continuous_vars
  } else if (settings$variable_selection_method == "RF_Importance") {
    if (verbose) message("    Running RF Importance selection module...")
    rf_output <- .run_rf_importance_module(
      data_for_rf = data_for_selection,
      target_id = target_id,
      potential_vars = continuous_vars_available,
      species_col_actual_name = species_col_actual_name,
      settings = settings,
      output_dir_rf = output_dir_target_specific,
      verbose = verbose
    )
    initial_selected_vars <- rf_output$selected_variables
    selection_details$rf_importance_scores <- rf_output$importance_scores
    selection_details$rf_model_oob_error <- rf_output$rf_model_oob_error
  } else if (settings$variable_selection_method == "Manual") {
    if (verbose) message("    Using manual variable selection...")
    if (!is.null(settings$manual_selected_vars) && target_id %in% names(settings$manual_selected_vars)) {
      initial_selected_vars <- settings$manual_selected_vars[[target_id]]
      initial_selected_vars <- intersect(initial_selected_vars, all_potential_env_vars)
      if (verbose) message("      Manually selected variables for ", target_id, ": ", paste(initial_selected_vars, collapse = ", "))
    } else {
      warning("Manual selection chosen, but no variables specified for target '", target_id, "' or target not in list. Using all available variables (continuous and categorical).")
      initial_selected_vars <- all_potential_env_vars
    }
    selection_details$manual_selection_list_used <- settings$manual_selected_vars[[target_id]] %||% "N/A"
  } else {
    warning("Unknown variable selection method: ", settings$variable_selection_method, ". Using all available non-categorical variables.")
    initial_selected_vars <- continuous_vars_available
  }
  if (length(initial_selected_vars) == 0) {
    if (verbose) message("    No variables selected by the initial method for ", target_id, ". Cannot proceed with filters.")
    return(list(
      selected_variables = character(0),
      initial_selection_method = settings$variable_selection_method,
      initial_selected_variables = character(0),
      details = selection_details
    ))
  }
  
  # --- 2. Post-selection Filters (Correlation and VIF) ---
  clean_all_potential_env_vars <- all_potential_env_vars[!is.na(all_potential_env_vars)]
  clean_categorical_vars <- categorical_vars[!is.na(categorical_vars)]
  clean_initial_selected_vars <- initial_selected_vars[!is.na(initial_selected_vars)]
  continuous_vars_available_cleaned <- setdiff(clean_all_potential_env_vars, clean_categorical_vars)
  continuous_vars_from_initial_selection <- intersect(clean_initial_selected_vars, continuous_vars_available_cleaned)
  # Correctly handle categorical variables based on selection method
  if (settings$variable_selection_method == "Manual") {
    # If manual, only keep the categorical vars that were manually selected
    categorical_vars_to_keep <- intersect(clean_initial_selected_vars, clean_categorical_vars)
  } else {
    # For Boruta/RF, keep all available categorical vars as they were not part of the selection
    categorical_vars_to_keep <- intersect(clean_all_potential_env_vars, clean_categorical_vars)
  }
  vars_after_correlation <- continuous_vars_from_initial_selection
  if (length(continuous_vars_from_initial_selection) >= 2) {
    cor_filter_output <- .filter_by_correlation(
      data_for_correlation = data_for_selection,
      selected_vars = continuous_vars_from_initial_selection,
      correlation_threshold = settings$correlation_threshold,
      output_dir_corr = output_dir_target_specific,
      target_id = target_id,
      verbose = verbose
    )
    vars_after_correlation <- cor_filter_output$remaining_variables
    selection_details$correlation_matrix <- cor_filter_output$correlation_matrix
    selection_details$removed_by_correlation <- cor_filter_output$removed_variables
  } else {
    if (verbose) message("    Skipping correlation filter: less than 2 numeric variables initially selected.")
    selection_details$correlation_details <- "Skipped ( < 2 continuous vars from initial selection)"
    selection_details$removed_by_correlation <- character(0)
  }
  if (length(vars_after_correlation) == 0 && length(continuous_vars_from_initial_selection) > 0) {
    if (verbose) message("    No numeric variables remaining after correlation filter for ", target_id)
  } else if (length(vars_after_correlation) > 0) {
    if (verbose) message("    Numeric variables after correlation filter: ", paste(vars_after_correlation, collapse = ", "))
  }
  final_numeric_vars <- vars_after_correlation
  if (length(vars_after_correlation) >= 2) {
    vif_filter_output <- .filter_by_vif(
      data_for_vif = data_for_selection,
      selected_vars = vars_after_correlation,
      vif_threshold = settings$vif_threshold,
      output_dir_vif = output_dir_target_specific,
      target_id = target_id,
      verbose = verbose
    )
    final_numeric_vars <- vif_filter_output$remaining_variables
    selection_details$vif_values <- vif_filter_output$vif_values
    selection_details$removed_by_vif <- vif_filter_output$removed_variables
  } else {
    if (verbose) message("    Skipping VIF filter: less than 2 numeric variables after correlation filter.")
    selection_details$vif_details <- "Skipped ( < 2 continuous vars after correlation)"
    selection_details$removed_by_vif <- character(0)
  }
  final_selected_vars <- unique(c(final_numeric_vars, categorical_vars_to_keep))
  return(list(
    selected_variables = final_selected_vars,
    initial_selection_method = settings$variable_selection_method,
    initial_selected_variables = initial_selected_vars,
    details = selection_details
  ))
}

#' Run Target-Specific Hypervolume Dendrogram Module (Internal)
#'
#' Constructs hypervolumes for a target species and its comparison set using
#' the variables selected for that target, calculates pairwise overlap, and
#' generates a dendrogram.
#'
#' @param data_for_module Tibble. Data for the target and its comparison species,
#'        containing the species ID column and the selected variables for this target.
#' @param target_id Character. The ID of the target species.
#' @param comparison_ids Character vector. IDs of comparison species for this target.
#' @param vars_for_hv Character vector. Numeric variables selected for this target's analyses.
#' @param species_id_col Character. Name of the species ID column.
#' @param settings An object of class 'niche_analysis_settings'.
#' @param output_dir_target_specific Character. Output directory for this target.
#' @param global_species_palette Named character vector. Palette for coloring.
#' @param verbose Logical.
#'
#' @return A list containing the plot path, overlap matrix, and hclust object, or NULL.
#' @noRd
.run_dendrogram_target_specific_module <- function(data_for_module,
                                                   target_id,
                                                   comparison_ids,
                                                   vars_for_hv,
                                                   species_id_col,
                                                   settings,
                                                   output_dir_target_specific,
                                                   global_species_palette,
                                                   verbose = TRUE) {
  if (verbose) message("    Running Target-Specific Hypervolume Dendrogram module for: ", target_id)
  
  species_in_this_analysis <- unique(c(target_id, comparison_ids))
  if (length(species_in_this_analysis) < 2) {
    if (verbose) message("      Target-Specific HV Dendro: Less than 2 species for this target's dendrogram. Skipping.")
    return(NULL)
  }
  if (length(vars_for_hv) < 2) {
    if (verbose) message("      Target-Specific HV Dendro: Less than 2 variables for hypervolumes. Skipping.")
    return(NULL)
  }
  
  hv_list_target_specific <- list()
  for (sp_id in species_in_this_analysis) {
    sp_data <- data_for_module %>%
      dplyr::filter(!!dplyr::sym(species_id_col) == sp_id) %>%
      dplyr::select(dplyr::all_of(vars_for_hv)) %>%
      stats::na.omit()
    if (nrow(sp_data) <= length(vars_for_hv)) {
      if (verbose) warning("      Target-Specific HV Dendro: Insufficient data for ", sp_id, ". Skipping this species.")
      next
    }
    sp_data_scaled <- scale(sp_data, center = TRUE, scale = TRUE)
    hv_list_target_specific[[sp_id]] <- tryCatch({
      suppressMessages(hypervolume::hypervolume_gaussian(sp_data_scaled, name = sp_id, samples.per.point = settings$hypervolume_samples_per_point, verbose = FALSE))
    }, error = function(e) { NULL })
  }
  hv_list_target_specific <- hv_list_target_specific[!sapply(hv_list_target_specific, is.null)]
  if (length(hv_list_target_specific) < 2) {
    if (verbose) message("      Target-Specific HV Dendro: Less than 2 valid hypervolumes. Skipping.")
    return(NULL)
  }
  
  metric_for_target_dendro <- settings$dendrogram_target_specific_hv_metric %||% settings$dendrogram_global_hv_metric
  overlap_matrix_target <- .calculate_pairwise_hv_overlap_matrix(hv_list_target_specific, metric_for_target_dendro)
  if(is.null(overlap_matrix_target) || all(is.na(overlap_matrix_target[lower.tri(overlap_matrix_target)]))) {
    if (verbose) message("      Target-Specific HV Dendro: Overlap matrix is NULL or all NA. Skipping.")
    return(NULL)
  }
  dist_matrix_target <- stats::as.dist(1 - overlap_matrix_target)
  hclust_obj_target <- stats::hclust(dist_matrix_target, method = settings$dendrogram_global_hv_cluster_method %||% "average")
  
  plot_path_target <- .plot_hypervolume_dendrogram(
    hclust_object = hclust_obj_target,
    output_path = file.path(output_dir_target_specific, paste0("dendrogram_target_specific_hv_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".png")),
    title = paste("Niche Dendrogram for", target_id, "and Comparisons\n(Hypervolume", metric_for_target_dendro, "Overlap)"),
    target_species_highlight = target_id,
    global_species_palette = global_species_palette,
    verbose = verbose
  )
  if (verbose) message("    Target-Specific Hypervolume Dendrogram module completed. Plot saved to: ", plot_path_target)
  return(list(plot_path = plot_path_target, overlap_matrix = overlap_matrix_target, hclust_object = hclust_obj_target))
}

#' Run Global Hypervolume Dendrogram Module (Internal)
#'
#' Constructs hypervolumes for all species in the dataset using a common set of variables,
#' calculates pairwise overlap, and generates a dendrogram based on hierarchical clustering.
#'
#' @param niche_data_object An object of class 'niche_data'.
#' @param settings An object of class 'niche_analysis_settings'.
#' @param output_dir_base Character. The main output directory for the analysis.
#'        A subdirectory "global_dendrogram_hv" will be created here.
#' @param global_species_palette Named character vector. Palette for coloring dendrogram labels if needed.
#' @param verbose Logical.
#'
#' @return A list containing the path to the saved dendrogram plot, the overlap matrix,
#'         and the hclust object, or NULL if the analysis cannot be performed.
#' @noRd
.run_dendrogram_global_hv_module <- function(niche_data_object,
                                             settings,
                                             output_dir_base,
                                             global_species_palette, 
                                             verbose = TRUE) {
  if (verbose) message("  Running Global Hypervolume Dendrogram module...")
  module_output_dir <- file.path(output_dir_base, "global_dendrogram_hv")
  if (!dir.exists(module_output_dir)) {
    dir.create(module_output_dir, recursive = TRUE)
  }
  all_species_ids <- unique(niche_data_object$data[[niche_data_object$metadata$species_id_col]])
  if (length(all_species_ids) < 2) {
    if (verbose) message("    Global HV Dendrogram: Less than 2 species in the dataset. Skipping.")
    return(NULL)
  }
  continuous_vars_all <- niche_data_object$metadata$env_vars$continuous
  if (length(continuous_vars_all) == 0) {
    if (verbose) message("    Global HV Dendrogram: No continuous variables available. Skipping.")
    return(NULL)
  }
  vars_for_hv_global <- continuous_vars_all
  if (length(continuous_vars_all) > settings$dendrogram_global_hv_n_vars) {
    if (verbose) message("    Global HV Dendrogram: More than ", settings$dendrogram_global_hv_n_vars, " continuous variables available. Selecting top ", settings$dendrogram_global_hv_n_vars, " (currently by order).")
    vars_for_hv_global <- head(continuous_vars_all, settings$dendrogram_global_hv_n_vars)
  }
  if (length(vars_for_hv_global) < 2) {
    if (verbose) message("    Global HV Dendrogram: Less than 2 variables selected/available for hypervolume construction. Skipping.")
    return(NULL)
  }
  if (verbose) message("    Global HV Dendrogram: Using variables: ", paste(vars_for_hv_global, collapse = ", "))
  hv_list_global <- list()
  for (sp_id in all_species_ids) {
    sp_data <- niche_data_object$data %>%
      dplyr::filter(!!dplyr::sym(niche_data_object$metadata$species_id_col) == sp_id) %>%
      dplyr::select(dplyr::all_of(vars_for_hv_global)) %>%
      stats::na.omit()
    if (nrow(sp_data) <= length(vars_for_hv_global)) {
      if (verbose) warning("      Global HV Dendrogram: Insufficient data points for ", sp_id, " (", nrow(sp_data), " points for ", length(vars_for_hv_global), " vars). Skipping this species.")
      next
    }
    sp_data_scaled <- scale(sp_data, center = TRUE, scale = TRUE)
    hv_list_global[[sp_id]] <- tryCatch({
      suppressMessages(hypervolume::hypervolume_gaussian(sp_data_scaled, name = sp_id, samples.per.point = settings$hypervolume_samples_per_point, verbose = FALSE))
    }, error = function(e) {
      if (verbose) warning("      Global HV Dendrogram: Hypervolume construction failed for ", sp_id, ": ", e$message)
      NULL
    })
  }
  hv_list_global <- hv_list_global[!sapply(hv_list_global, is.null)]
  if (length(hv_list_global) < 2) {
    if (verbose) message("    Global HV Dendrogram: Less than 2 valid hypervolumes constructed. Skipping.")
    return(NULL)
  }
  overlap_matrix <- .calculate_pairwise_hv_overlap_matrix(hv_list_global, settings$dendrogram_global_hv_metric)
  if(is.null(overlap_matrix)){
    if(verbose) message("    Global HV Dendrogram: Failed to calculate overlap matrix. Skipping.")
    return(NULL)
  }
  utils::write.csv(overlap_matrix, file.path(module_output_dir, "global_hypervolume_overlap_matrix.csv"), row.names = TRUE)
  if (all(is.na(overlap_matrix[lower.tri(overlap_matrix)]))) {
    if (verbose) message("    Global HV Dendrogram: Overlap matrix consists only of NAs (excluding diagonal). Cannot perform clustering. Skipping.")
    return(NULL)
  }
  dist_matrix <- stats::as.dist(1 - overlap_matrix)
  if (all(is.na(dist_matrix)) || all(!is.finite(dist_matrix))) {
    if (verbose) message("    Global HV Dendrogram: Distance matrix is all NA or non-finite after creation. Cannot perform clustering. Skipping.")
    return(NULL)
  }
  hclust_obj <- stats::hclust(dist_matrix, method = settings$dendrogram_global_hv_cluster_method %||% "average")
  saveRDS(hclust_obj, file.path(module_output_dir, "global_hypervolume_hclust_object.rds"))
  # Plot dendrogram
  plot_path <- .plot_hypervolume_dendrogram(hclust_object = hclust_obj,
                                            output_path = file.path(module_output_dir, "dendrogram_global_hypervolume_overlap.png"),
                                            title = paste("Global Niche Dendrogram (Hypervolume", settings$dendrogram_global_hv_metric, "Overlap)"),
                                            target_species_highlight = NULL,
                                            global_species_palette = global_species_palette,
                                            verbose = verbose)
  if (verbose) message("    Global Hypervolume Dendrogram module completed. Plot saved to: ", plot_path)
  return(list(plot_path = plot_path, overlap_matrix = overlap_matrix, hclust_object = hclust_obj, vars_used = vars_for_hv_global))
}

#' Calculate Pairwise Hypervolume Overlap Matrix (Internal)
#' @noRd
.calculate_pairwise_hv_overlap_matrix <- function(hv_list, overlap_metric_name) {
  species_for_matrix <- names(hv_list)
  overlap_matrix <- matrix(NA, nrow = length(species_for_matrix), ncol = length(species_for_matrix),
                           dimnames = list(species_for_matrix, species_for_matrix))
  for (i in 1:length(species_for_matrix)) {
    for (j in i:length(species_for_matrix)) {
      if (i == j) {
        overlap_matrix[i, j] <- 1.0
      } else {
        hv_set <- tryCatch(hypervolume::hypervolume_set(hv_list[[species_for_matrix[i]]], hv_list[[species_for_matrix[j]]], check.memory = FALSE, verbose = FALSE), error = function(e) NULL)
        if (!is.null(hv_set)) {
          overlap_stats_set <- hypervolume::hypervolume_overlap_statistics(hv_set)
          overlap_value <- overlap_stats_set[[overlap_metric_name]] %||% NA_real_
          overlap_matrix[i, j] <- overlap_matrix[j, i] <- overlap_value
        }
      }
    }
  }
  return(overlap_matrix)
}

# --- Overlap Module Helper Functions ---
#' Normalize Continuous Data into Bins (Internal)
#' @noRd
.normalize_to_bins <- function(x, n_bins = 10) {
  x_clean <- stats::na.omit(x)
  if (length(unique(x_clean)) < 2) {
    if(length(x_clean) > 0) {
      counts <- rep(0, n_bins)
      counts[ceiling(n_bins/2)] <- length(x_clean)
      return(counts / sum(counts))
    } else {
      return(rep(0, n_bins))
    }
  }
  min_x <- min(x_clean)
  max_x <- max(x_clean)
  if (min_x == max_x) {
    counts <- rep(0, n_bins)
    counts[ceiling(n_bins/2)] <- length(x_clean)
    return(counts / sum(counts))
  }
  breaks <- seq(min_x, max_x, length.out = n_bins + 1)
  breaks[n_bins + 1] <- breaks[n_bins + 1] + .Machine$double.eps
  h <- graphics::hist(x_clean, breaks = breaks, plot = FALSE, include.lowest = TRUE, right = FALSE)
  total_counts <- sum(h$counts)
  if (total_counts == 0) return(rep(0, n_bins))
  h$counts / total_counts
}

#' Calculate Schoener's D for Continuous Data (Internal)
#' @noRd
.calculate_schoener_D_continuous <- function(freq1, freq2) {
  if (length(freq1) != length(freq2) || any(is.na(c(freq1, freq2)))) return(NA_real_)
  1 - 0.5 * sum(abs(freq1 - freq2))
}

#' Calculate Pianka's O for Continuous Data (Internal)
#' @noRd
.calculate_pianka_O_continuous <- function(freq1, freq2) {
  if (length(freq1) != length(freq2) || any(is.na(c(freq1, freq2)))) return(NA_real_)
  numerator <- sum(freq1 * freq2)
  denominator <- sqrt(sum(freq1^2) * sum(freq2^2))
  if (denominator == 0) return(NA_real_) else numerator / denominator
}

#' Calculate Czekanowski's PSI for Proportions (Internal)
#' @noRd
.calculate_czekanowski_psi_proportions <- function(prop1, prop2) {
  if (length(prop1) != length(prop2) || any(is.na(c(prop1, prop2)))) return(NA_real_)
  sum(pmin(prop1, prop2))
}

#' Calculate Hurlbert's L for Proportions (Internal)
#' @noRd
.calculate_hurlbert_L_proportions <- function(prop1, prop2) {
  if (length(prop1) != length(prop2) || any(is.na(c(prop1, prop2)))) return(NA_real_)
  sum_p1_sq <- sum(prop1^2)
  sum_p2_sq <- sum(prop2^2)
  if (sum_p1_sq == 0 || sum_p2_sq == 0) return(NA_real_)
  B1 <- 1 / sum_p1_sq 
  B2 <- 1 / sum_p2_sq
  psi <- .calculate_czekanowski_psi_proportions(prop1, prop2)
  denominator <- sqrt(B1 * B2)
  if(is.na(denominator) || denominator == 0) NA_real_ else psi / denominator
}

#' Calculate Overlap Indices for a Single Categorical Variable (Internal)
#' @noRd
.calculate_categorical_var_overlap <- function(data_pair, species_col, var_cat, sp1_id, sp2_id) {
  if (!is.factor(data_pair[[var_cat]])) data_pair[[var_cat]] <- as.factor(data_pair[[var_cat]])
  data_filtered <- data_pair %>%
    dplyr::filter(!!dplyr::sym(species_col) %in% c(sp1_id, sp2_id)) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(var_cat), ~forcats::fct_drop(.x)))
  if (nrow(data_filtered) < 2 || length(levels(data_filtered[[var_cat]])) < 1) {
    return(list(schoener = NA_real_, pianka = NA_real_, czek_psi = NA_real_, hurlbert_l = NA_real_))
  }
  prop_table <- data_filtered %>%
    dplyr::count(!!dplyr::sym(species_col), !!dplyr::sym(var_cat)) %>%
    dplyr::group_by(!!dplyr::sym(species_col)) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup() %>%
    tidyr::complete(!!dplyr::sym(species_col), !!dplyr::sym(var_cat), fill = list(n = 0, prop = 0))
  prop1 <- prop_table %>% dplyr::filter(!!dplyr::sym(species_col) == sp1_id) %>% dplyr::pull(prop)
  prop2 <- prop_table %>% dplyr::filter(!!dplyr::sym(species_col) == sp2_id) %>% dplyr::pull(prop)
  if(length(prop1) == 0 || length(prop2) == 0 || length(prop1) != length(prop2)){ 
    return(list(schoener = NA_real_, pianka = NA_real_, czek_psi = NA_real_, hurlbert_l = NA_real_))
  }
  tryCatch({
    list(
      schoener = 1 - 0.5 * sum(abs(prop1 - prop2)),
      pianka = sum(prop1 * prop2) / sqrt(sum(prop1^2) * sum(prop2^2)),
      czek_psi = .calculate_czekanowski_psi_proportions(prop1, prop2),
      hurlbert_l = .calculate_hurlbert_L_proportions(prop1, prop2)
    )
  }, error = function(e) {
    list(schoener = NA_real_, pianka = NA_real_, czek_psi = NA_real_, hurlbert_l = NA_real_)
  })
}

#' Run Niche Overlap Indices Module (Internal)
#'
#' Calculates various niche overlap indices between a target species and
#' a set of comparison species.
#'
#' @param data_for_module Tibble. Data containing species ID and selected environmental variables.
#' @param target_id Character. ID of the target species.
#' @param comparison_ids Character vector. IDs of comparison species.
#' @param selected_vars Character vector. Names of selected environmental variables.
#' @param species_id_col Character. Name of the species ID column in `data_for_module`.
#' @param categorical_vars_names Character vector. Names of categorical variables among `selected_vars`.
#' @param output_dir_module Character. Path to save the CSV output.
#' @param settings An object of class 'niche_analysis_settings'.
#' @param verbose Logical.
#'
#' @return A tibble with pairwise overlap indices.
#' @noRd
.run_overlap_module <- function(data_for_module, target_id, comparison_ids,
                                selected_vars, species_id_col, categorical_vars_names,
                                output_dir_module, settings, verbose = TRUE) {
  
  if (verbose) message("    Calculating niche overlap indices...")
  # Filter selected_vars by type
  continuous_vars_to_process <- intersect(selected_vars, setdiff(names(data_for_module), categorical_vars_names)) 
  categorical_vars_to_process <- intersect(selected_vars, categorical_vars_names) 
  n_bins <- settings$overlap_n_bins
  if (verbose) message("      Overlap: Processing ", length(continuous_vars_to_process), " continuous and ", length(categorical_vars_to_process), " categorical variables.")
  all_overlap_results <- purrr::map_dfr(comparison_ids, function(comp_id) {
    if (verbose) message("      Overlap: ", target_id, " vs ", comp_id)
    data_pair <- data_for_module %>%
      dplyr::filter(!!dplyr::sym(species_id_col) %in% c(target_id, comp_id))
    per_variable_overlaps <- list()
    # Continuous variables
    if (length(continuous_vars_to_process) > 0) {
      for (var_cont in continuous_vars_to_process) {
        d1 <- data_pair[[var_cont]][data_pair[[species_id_col]] == target_id]
        d2 <- data_pair[[var_cont]][data_pair[[species_id_col]] == comp_id]
        if(length(stats::na.omit(d1)) == 0 || length(stats::na.omit(d2)) == 0) {
          per_variable_overlaps[[var_cont]] <- list(schoener_D = NA_real_, pianka_O = NA_real_, czek_psi = NA_real_, hurlbert_L = NA_real_)
          next
        }
        freq1 <- .normalize_to_bins(d1, n_bins)
        freq2 <- .normalize_to_bins(d2, n_bins)
        per_variable_overlaps[[var_cont]] <- list(
          schoener_D = .calculate_schoener_D_continuous(freq1, freq2),
          pianka_O = .calculate_pianka_O_continuous(freq1, freq2),
          czek_psi = .calculate_czekanowski_psi_proportions(freq1, freq2), 
          hurlbert_L = .calculate_hurlbert_L_proportions(freq1, freq2)    
        )
      }
    }
    # Categorical variables
    if (length(categorical_vars_to_process) > 0) {
      for (var_cat in categorical_vars_to_process) {
        per_variable_overlaps[[var_cat]] <- .calculate_categorical_var_overlap(data_pair, species_id_col, var_cat, target_id, comp_id)
      }
    }
    if (length(per_variable_overlaps) == 0) {
      return(tibble::tibble(target_species = target_id, comparison_species = comp_id,
                            mean_schoener_D = NA_real_, mean_pianka_O = NA_real_,
                            mean_czek_psi = NA_real_, mean_hurlbert_L = NA_real_, composite_overlap = NA_real_))
    }
    # Calculate mean indices
    mean_schoener <- mean(purrr::map_dbl(per_variable_overlaps, "schoener_D", .default = NA_real_), na.rm = TRUE)
    mean_pianka <- mean(purrr::map_dbl(per_variable_overlaps, "pianka_O", .default = NA_real_), na.rm = TRUE)
    mean_czek <- mean(purrr::map_dbl(per_variable_overlaps, "czek_psi", .default = NA_real_), na.rm = TRUE)
    mean_hurlbert <- mean(purrr::map_dbl(per_variable_overlaps, "hurlbert_L", .default = NA_real_), na.rm = TRUE)
    # Composite index (weights from your script)
    composite <- 0.4 * mean_schoener + 0.3 * mean_czek + 0.3 * mean_hurlbert
    tibble::tibble(target_species = target_id, comparison_species = comp_id,
                   mean_schoener_D = mean_schoener, mean_pianka_O = mean_pianka,
                   mean_czek_psi = mean_czek, mean_hurlbert_L = mean_hurlbert,
                   composite_overlap = composite)
  })
  # Save to CSV
  csv_filename <- file.path(output_dir_module, paste0("overlap_indices_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".csv"))
  utils::write.csv(all_overlap_results, csv_filename, row.names = FALSE)
  if (verbose) message("    Overlap indices saved to: ", csv_filename)
  return(all_overlap_results)
}

# # --- PERMANOVA Module Helper Function ---

#' Run Pairwise PERMANOVA Module (Internal)
#'
#' Performs pairwise PERMANOVA tests between a target species and comparison species
#' using either Gower or Euclidean distance.
#'
#' @param data_for_module Tibble. Data containing species ID and selected environmental variables.
#' @param target_id Character. ID of the target species.
#' @param comparison_ids Character vector. IDs of comparison species.
#' @param selected_vars Character vector. Names of selected environmental variables.
#' @param species_id_col Character. Name of the species ID column in `data_for_module`.
#' @param categorical_vars_names Character vector. Names of categorical variables among `selected_vars`.
#' @param output_dir_module Character. Path to save the CSV output.
#' @param settings An object of class 'niche_analysis_settings'.
#' @param verbose Logical.
#'
#' @return A tibble with pairwise PERMANOVA results.
#' @noRd
.run_permanova_pairwise_module <- function(data_for_module, target_id, comparison_ids,
                                           selected_vars, species_id_col, categorical_vars_names,
                                           output_dir_module, settings, verbose = TRUE) {
  
  if (verbose) message("    Calculating pairwise PERMANOVA tests...")
  use_gower <- settings$permanova_use_gower
  n_permutations <- settings$n_bootstrap_stability %||% 999 
  if (use_gower) {
    vars_for_permanova <- selected_vars 
    if (verbose) message("      Using Gower distance (mixed variables).")
  } else {
    vars_for_permanova <- setdiff(selected_vars, categorical_vars_names)
    if (verbose) message("      Using Euclidean distance (standardized numeric variables).")
    if (length(vars_for_permanova) < 1) {
      if (verbose) message("      PERMANOVA: No numeric variables selected for Euclidean distance. Skipping module.")
      return(tibble::tibble(
        specie1 = target_id,
        specie2 = comparison_ids,
        F_value = NA_real_, R2 = NA_real_, p_value = NA_real_, p_betadisper = NA_real_,
        n_tot = NA_integer_, n1 = NA_integer_, n2 = NA_integer_,
        warning = "No numeric variables available for Euclidean", p_adj = NA_real_
      ))
    }
  }
  if (verbose) message("      Variables used for PERMANOVA: ", paste(vars_for_permanova, collapse = ", "))
  all_permanova_results <- purrr::map_dfr(comparison_ids, function(comp_id) {
    if (verbose) message("      PERMANOVA: ", target_id, " vs ", comp_id)
    # Subset data for the pair and selected variables
    sub_data <- data_for_module %>%
      dplyr::filter(!!dplyr::sym(species_id_col) %in% c(target_id, comp_id)) %>%
      dplyr::select(!!dplyr::sym(species_id_col), dplyr::all_of(vars_for_permanova)) %>%
      stats::na.omit()
    # Check for sufficient data points per group
    n1 <- sum(sub_data[[species_id_col]] == target_id)
    n2 <- sum(sub_data[[species_id_col]] == comp_id)
    if (n1 < 2 || n2 < 2 || nrow(sub_data) < 3) {
      return(tibble::tibble(specie1 = target_id, specie2 = comp_id, F_value = NA_real_, R2 = NA_real_, p_value = NA_real_, p_betadisper = NA_real_, n_tot = nrow(sub_data), n1 = n1, n2 = n2, warning = "Insufficient samples (<2 per group or <3 total)"))
    }
    # Prepare input based on distance metric
    permanova_input <- NULL
    error_msg <- NA_character_
    dist_matrix_for_betadisper <- NULL
    if (use_gower) {
      # Gower distance requires factors for categorical variables
      data_for_dist <- sub_data %>%
        dplyr::select(dplyr::all_of(vars_for_permanova)) %>%
        dplyr::mutate(dplyr::across(dplyr::intersect(vars_for_permanova, categorical_vars_names), as.factor))
      # Check for zero variance in numeric columns before distance calculation
      numeric_cols <- data_for_dist %>% dplyr::select(where(is.numeric))
      if (ncol(numeric_cols) > 0) {
        sd_check <- apply(numeric_cols, 2, stats::sd, na.rm = TRUE)
        if (any(sd_check == 0, na.rm = TRUE)) {
          error_msg <- "Zero variance in one or more numeric variables for Gower"
        }
      }
      if (is.na(error_msg)) {
        dist_matrix <- tryCatch({
          cluster::daisy(data_for_dist, metric = "gower", warnBin = FALSE, warnAsym = FALSE)
        }, error = function(e) {
          warning(paste("Error calculating Gower distance between", target_id, "and", comp_id, ":", e$message))
          NULL
        })
        if (is.null(dist_matrix)) {
          error_msg <- "Error calculating Gower distance"
        } else {
          permanova_input <- dist_matrix
          dist_matrix_for_betadisper <- dist_matrix
        }
      }
    } else {
      data_for_dist <- sub_data %>% dplyr::select(dplyr::all_of(vars_for_permanova)) 
      # Check for zero variance
      sd_check <- apply(data_for_dist, 2, stats::sd, na.rm = TRUE)
      if (any(sd_check == 0, na.rm = TRUE)) {
        error_msg <- "Zero variance in one or more numeric variables for Euclidean"
      } else {
        # Standardize data
        permanova_input <- vegan::decostand(data_for_dist, "standardize")
        # Calculate Euclidean distance for betadisper
        dist_matrix_for_betadisper <- tryCatch({
          vegan::vegdist(permanova_input, method = "euclidean")
        }, error = function(e) {
          warning(paste("Error calculating Euclidean distance for Betadisper between", target_id, "and", comp_id, ":", e$message))
          NULL
        })
        if (is.null(dist_matrix_for_betadisper)) error_msg <- "Error calculating Euclidean dist (Betadisper)"
      }
    }
    if (!is.na(error_msg)) {
      return(tibble::tibble(specie1 = target_id, specie2 = comp_id, F_value = NA_real_, R2 = NA_real_, p_value = NA_real_, p_betadisper = NA_real_, n_tot = nrow(sub_data), n1 = n1, n2 = n2, warning = error_msg))
    }
    # Perform adonis2
    perm_res <- tryCatch({
      data_for_formula <- sub_data[, species_id_col, drop = FALSE]
      formula_str <- paste0("permanova_input ~ `", species_id_col, "`")
      current_formula <- stats::as.formula(formula_str)
      vegan::adonis2(current_formula, data = data_for_formula, permutations = n_permutations)
    }, error = function(e) {
      detailed_error_msg <- paste("Error during adonis2 between", target_id, "and", comp_id, "using vars:", paste(vars_for_permanova, collapse=", "), ". Error: ", e$message)
      if(verbose) warning(detailed_error_msg)
      NULL
    })
    betadisper_p <- NA_real_
    if (!is.null(dist_matrix_for_betadisper) && length(unique(sub_data[[species_id_col]])) > 1) {
      betadisp_res <- tryCatch({
        vegan::betadisper(dist_matrix_for_betadisper, sub_data[[species_id_col]])
      }, error = function(e) {
        warning(paste("Error during betadisper between", target_id, "and", comp_id, ":", e$message))
        NULL
      })
      if (!is.null(betadisp_res)) {
        betadisp_anova <- tryCatch({
          stats::anova(betadisp_res)
        }, error = function(e) {
          warning(paste("Error during anova(betadisper) between", target_id, "and", comp_id, ":", e$message))
          NULL
        })
        if (!is.null(betadisp_anova) && "Pr(>F)" %in% names(betadisp_anova)) {
          betadisper_p <- betadisp_anova$`Pr(>F)`[1]
        }
      }
    }
    if (is.null(perm_res) || !is.list(perm_res) || !all(c("F", "R2", "Pr(>F)") %in% names(perm_res)) || is.null(perm_res$F) || is.null(perm_res$`Pr(>F)`)) {
      return(tibble::tibble(specie1 = target_id, specie2 = comp_id, F_value = NA_real_, R2 = NA_real_, p_value = NA_real_, p_betadisper = betadisper_p, n_tot = nrow(sub_data), n1 = n1, n2 = n2, warning = "Error or invalid adonis2 result"))
    } else {
      return(tibble::tibble(
        specie1 = target_id,
        specie2 = comp_id,
        F_value = perm_res$F[1],
        R2 = perm_res$R2[1],
        p_value = perm_res$`Pr(>F)`[1],
        p_betadisper = betadisper_p,
        n_tot = nrow(sub_data),
        n1 = n1,
        n2 = n2,
        warning = NA_character_
      ))
    }
  })
  if (nrow(all_permanova_results) > 0 && any(!is.na(all_permanova_results$p_value))) {
    all_permanova_results <- all_permanova_results %>%
      dplyr::mutate(p_adj = stats::p.adjust(p_value, method = "BH"))
  } else if (nrow(all_permanova_results) > 0) {
    all_permanova_results$p_adj <- NA_real_ # Add p_adj column even if all are NA
  }
  csv_filename <- file.path(output_dir_module, paste0("permanova_pairwise_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".csv"))
  utils::write.csv(all_permanova_results, csv_filename, row.names = FALSE)
  if (verbose) message("    Pairwise PERMANOVA results saved to: ", csv_filename)
  
  return(all_permanova_results)
}

# --- LDA Module Helper Function ---
#' Run Linear Discriminant Analysis (LDA) Module (Internal)
#'
#' Performs LDA to discriminate between the target species and comparison species
#' (or among all species in the provided data subset).
#'
#' @param data_for_module Tibble. Data containing species ID and selected environmental variables.
#'        Should include the target and all relevant comparison species.
#' @param target_id Character. ID of the target species (used for context, not for filtering if
#'        `comparison_ids` includes all species for a multigroup LDA).
#' @param comparison_ids Character vector. IDs of comparison species. For multigroup LDA,
#'        this might be all species present in `data_for_module` excluding the target,
#'        or all species including the target if `data_for_module` is already subsetted.
#' @param selected_vars Character vector. Names of selected environmental variables.
#' @param species_id_col Character. Name of the species ID column in `data_for_module`.
#' @param categorical_vars_names Character vector. Names of categorical variables (to be excluded).
#' @param output_dir_module Character. Path to save the CSV outputs and plots.
#' @param settings An object of class 'niche_analysis_settings'.
#' @param verbose Logical.
#'
#' @return A list containing the LDA model, predictions, confusion matrix, accuracy,
#'         and variable importance, or NULL if LDA fails.
#' @noRd
.run_lda_module <- function(data_for_module, target_id, comparison_ids,
                            selected_vars, species_id_col, categorical_vars_names,
                            output_dir_module, settings, verbose = TRUE) {
  
  if (verbose) message("    Running Linear Discriminant Analysis (LDA) module...")
  
  numeric_vars_lda <- setdiff(selected_vars, categorical_vars_names)
  if (length(numeric_vars_lda) < 1) {
    if (verbose) message("      LDA: No numeric variables selected. Skipping module.")
    return(NULL)
  }
  lda_data_analysis <- data_for_module %>%
    dplyr::select(dplyr::all_of(species_id_col), dplyr::all_of(numeric_vars_lda)) %>%
    stats::na.omit() %>%
    dplyr::mutate(!!dplyr::sym(species_id_col) := factor(!!dplyr::sym(species_id_col)))
  group_counts <- lda_data_analysis %>% dplyr::count(!!dplyr::sym(species_id_col))
  if (any(group_counts$n < 2) || nrow(lda_data_analysis) < (length(numeric_vars_lda) + nrow(group_counts)) || nrow(group_counts) < 2) {
    if (verbose) warning("      LDA: Insufficient data for LDA (groups < 2, total obs < vars + groups, or < 2 groups). Skipping.")
    return(NULL)
  }
  sd_check <- apply(lda_data_analysis[, numeric_vars_lda, drop = FALSE], 2, stats::sd, na.rm = TRUE)
  if (any(sd_check == 0, na.rm = TRUE)) {
    if (verbose) warning("      LDA: Zero variance in one or more numeric variables. Skipping.")
    return(NULL)
  }
  lda_formula_str <- paste0("`", species_id_col, "` ~ ", paste0("`", numeric_vars_lda, "`", collapse = " + "))
  lda_formula <- stats::as.formula(lda_formula_str)
  lda_model <- tryCatch({
    MASS::lda(lda_formula, data = lda_data_analysis, CV = FALSE)
  }, error = function(e) {
    if (verbose) warning("      LDA model fitting failed: ", e$message)
    NULL
  })
  if (is.null(lda_model)) return(NULL)
  lda_cv_predictions <- tryCatch({
    MASS::lda(lda_formula, data = lda_data_analysis, CV = TRUE)$class
  }, error = function(e) {
    if (verbose) warning("      LDA cross-validation failed: ", e$message)
    NULL
  })
  confusion_matrix_cv <- if (!is.null(lda_cv_predictions)) table(Observed = lda_data_analysis[[species_id_col]], Predicted = lda_cv_predictions) else NULL
  accuracy_cv <- if (!is.null(confusion_matrix_cv)) sum(diag(confusion_matrix_cv)) / sum(confusion_matrix_cv) else NA_real_
  predictions_train <- stats::predict(lda_model)
  lda_scores_df <- as.data.frame(predictions_train$x)
  utils::write.csv(lda_scores_df, file.path(output_dir_module, paste0("lda_scores_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".csv")), row.names = FALSE)
  variable_importance_ld1 <- NULL
  if (!is.null(lda_model$scaling) && ncol(lda_model$scaling) >= 1) {
    variable_importance_ld1 <- data.frame(
      Variable = rownames(lda_model$scaling),
      LD1_Loading = lda_model$scaling[, 1]
    ) %>% dplyr::arrange(dplyr::desc(abs(LD1_Loading)))
    utils::write.csv(variable_importance_ld1, file.path(output_dir_module, paste0("lda_variable_importance_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".csv")), row.names = FALSE)
  }
  if (verbose) message("    LDA module completed. Accuracy (CV): ", round(accuracy_cv, 3))
  return(list(
    lda_model = lda_model,
    predictions_on_training_data = predictions_train,
    cv_confusion_matrix = confusion_matrix_cv,
    cv_accuracy = accuracy_cv,
    ld_scores = lda_scores_df,
    variable_importance_ld1 = variable_importance_ld1,
    data_used_for_lda = lda_data_analysis, 
    numeric_vars_in_lda_model = numeric_vars_lda 
  ))
}

# --- Hypervolume Module Helper Function ---
#' Run Niche Hypervolume Analysis Module (Internal)
#'
#' Constructs hypervolumes for target and comparison species and calculates overlap.
#'
#' @param data_for_module Tibble. Data containing species ID and selected environmental variables.
#' @param target_id Character. ID of the target species.
#' @param comparison_ids Character vector. IDs of comparison species.
#' @param selected_vars Character vector. Names of all selected environmental variables.
#' @param species_id_col Character. Name of the species ID column in `data_for_module`.
#' @param categorical_vars_names Character vector. Names of categorical variables (to be excluded).
#' @param variable_selection_details List. Output from `.perform_variable_selection()`,
#'        containing details like Boruta stats or RF importance to select vars for hypervolume.
#' @param output_dir_module Character. Path to save hypervolume objects and CSV outputs.
#' @param settings An object of class 'niche_analysis_settings'.
#' @param verbose Logical.
#'
#' @return A list containing paths to saved hypervolumes, and a tibble of pairwise
#'         hypervolume overlap statistics, or NULL if analysis fails.
#' @noRd
.run_hypervolume_module <- function(data_for_module, target_id, comparison_ids,
                                    selected_vars, species_id_col, categorical_vars_names,
                                    variable_selection_details, 
                                    output_dir_module, settings, verbose = TRUE) {
  if (verbose) message("    Running Niche Hypervolume Analysis module...")
  numeric_vars_all <- setdiff(selected_vars, categorical_vars_names)
  vars_for_hv <- character(0)
  if (settings$variable_selection_method == "Boruta" && !is.null(variable_selection_details$boruta_stats)) {
    confirmed_boruta <- variable_selection_details$initial_selected_variables
    vars_for_hv <- intersect(confirmed_boruta, numeric_vars_all)
    if(length(vars_for_hv) > settings$hypervolume_n_vars) {
      top_boruta_for_hv <- variable_selection_details$boruta_stats %>%
        dplyr::filter(Original_Variable %in% vars_for_hv) %>%
        dplyr::slice_max(order_by = meanImp, n = settings$hypervolume_n_vars) %>%
        dplyr::pull(Original_Variable)
      vars_for_hv <- top_boruta_for_hv
    }
  } else if (settings$variable_selection_method == "RF_Importance" && !is.null(variable_selection_details$rf_importance_scores)) {
    vars_for_hv <- variable_selection_details$initial_selected_variables 
    vars_for_hv <- intersect(vars_for_hv, numeric_vars_all) 
  }
  if (length(vars_for_hv) == 0 || length(vars_for_hv) < 2) {
    vars_for_hv <- head(numeric_vars_all, settings$hypervolume_n_vars)
  }
  vars_for_hv <- head(vars_for_hv, settings$hypervolume_n_vars) 
  if (length(vars_for_hv) < 2) {
    if (verbose) message("      Hypervolume: Less than 2 numeric variables selected for hypervolume construction. Skipping module.")
    return(NULL)
  }
  if (verbose) message("      Hypervolume: Using variables: ", paste(vars_for_hv, collapse = ", "))
  hv_objects_dir <- file.path(output_dir_module, "hypervolume_objects")
  if (!dir.exists(hv_objects_dir)) dir.create(hv_objects_dir, recursive = TRUE)
  hypervolume_list <- list()
  build_save_hv <- function(sp_id, data, vars, hv_method, samples_per_point, path_dir, sp_verbose) {
    sp_data <- data %>%
      dplyr::filter(!!dplyr::sym(species_id_col) == sp_id) %>%
      dplyr::select(dplyr::all_of(vars)) %>%
      stats::na.omit()
    if (nrow(sp_data) <= length(vars)) {
      if (sp_verbose) warning("      Hypervolume: Insufficient data points for ", sp_id, " (", nrow(sp_data), " points for ", length(vars), " vars).")
      return(NULL)
    }
    sp_data_scaled <- scale(sp_data, center = TRUE, scale = TRUE)
    hv <- tryCatch({
      suppressMessages({
        hypervolume::hypervolume_gaussian(sp_data_scaled, name = sp_id,
                                          samples.per.point = samples_per_point,
                                          verbose = FALSE)      })
    }, error = function(e) {
      if (sp_verbose) warning("      Hypervolume construction failed for ", sp_id, ": ", e$message)
      NULL
    })
    
    if (!is.null(hv)) {
      saveRDS(hv, file.path(path_dir, paste0("hv_", gsub("[^A-Za-z0-9_.-]", "_", sp_id), ".rds")))
    }
    return(hv)
  }
  hypervolume_list[[target_id]] <- build_save_hv(target_id, data_for_module, vars_for_hv, settings$hypervolume_method, settings$hypervolume_samples_per_point, hv_objects_dir, verbose)
  for (comp_id in comparison_ids) {
    hypervolume_list[[comp_id]] <- build_save_hv(comp_id, data_for_module, vars_for_hv, settings$hypervolume_method, settings$hypervolume_samples_per_point, hv_objects_dir, verbose)
  }
  
  valid_hvs <- hypervolume_list[!sapply(hypervolume_list, is.null)]
  
  if (length(valid_hvs) < 2 || is.null(valid_hvs[[target_id]])) {
    if (verbose) message("      Hypervolume: Not enough valid hypervolumes constructed for overlap analysis.")
    return(list(hypervolume_objects_paths = list.files(hv_objects_dir, pattern = "\\.rds$", full.names = TRUE), overlap_stats = NULL, vars_used = vars_for_hv, hypervolume_objects = valid_hvs))
  }
  hv_target <- valid_hvs[[target_id]]
  overlap_stats_list <- purrr::map_dfr(setdiff(names(valid_hvs), target_id), function(comp_id_hv) {
    hv_comp <- valid_hvs[[comp_id_hv]]
    if (is.null(hv_comp)) return(tibble::tibble(target_species = target_id, comparison_species = comp_id_hv, jaccard=NA_real_, sorensen=NA_real_, frac_unique_target=NA_real_, frac_unique_comp=NA_real_))
    hv_set <- tryCatch(hypervolume::hypervolume_set(hv_target, hv_comp, check.memory = FALSE, verbose = FALSE), error = function(e) NULL)
    if (is.null(hv_set)) return(tibble::tibble(target_species = target_id, comparison_species = comp_id_hv, jaccard=NA_real_, sorensen=NA_real_, frac_unique_target=NA_real_, frac_unique_comp=NA_real_))
    stats <- hypervolume::hypervolume_overlap_statistics(hv_set)
    jaccard_val <- if("jaccard" %in% names(stats)) stats[['jaccard']] else NA_real_
    sorensen_val <- if("sorensen" %in% names(stats)) stats[['sorensen']] else NA_real_
    unique_a_val <- if("unique_a" %in% names(stats)) stats[['unique_a']] else NA_real_
    unique_b_val <- if("unique_b" %in% names(stats)) stats[['unique_b']] else NA_real_
    tibble::tibble(target_species = target_id, comparison_species = comp_id_hv,
                   jaccard = jaccard_val, sorensen = sorensen_val,
                   frac_unique_target = unique_a_val, frac_unique_comp = unique_b_val)
  })
  utils::write.csv(overlap_stats_list, file.path(output_dir_module, paste0("hypervolume_overlap_stats_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".csv")), row.names = FALSE)
  if (verbose) message("    Hypervolume overlap statistics saved.")
  return(list(
    hypervolume_objects_paths = list.files(hv_objects_dir, pattern = "\\.rds$", full.names = TRUE),
    overlap_stats = overlap_stats_list,
    vars_used = vars_for_hv,
    hypervolume_objects = valid_hvs
  ))
} 

# --- Variable Stability Module Helper Function ---

#' Run Variable Stability Analysis using Bootstrapped Boruta (Internal)
#'
#' Assesses the stability of variable selection by running Boruta on multiple
#' bootstrap samples of the data.
#'
#' @param data_for_stability Tibble. Data containing the 'species' column and
#'        the numeric environmental variables to be tested. This should be the
#'        same data used for the main Boruta run (`data_for_selection` in `.perform_variable_selection`).
#' @param target_id Character. ID of the target species.
#' @param potential_vars_for_stability Character vector. Names of numeric environmental
#'        variables to test for stability (typically those fed into the main Boruta).
#' @param settings An object of class 'niche_analysis_settings'.
#' @param output_dir_module Character. Path to save the stability frequency table.
#' @param verbose Logical.
#'
#' @return A tibble with variables and their selection frequency.
#' @noRd
.run_variable_stability_module <- function(data_for_stability,
                                           target_id,
                                           potential_vars_for_stability,
                                           species_col_actual_name,
                                           settings,
                                           output_dir_module,
                                           verbose = TRUE) {
  
  if (verbose) message("    Running Variable Stability Analysis (Bootstrapped Boruta)...")
  n_boot <- settings$n_bootstrap_stability %||% 50
  if (n_boot < 10) {
    if (verbose) warning("      Number of bootstrap iterations for stability (", n_boot, ") is low. Consider increasing for robust results.")
  }
  all_bootstrap_selections <- list()
  if (verbose) pb_boot <- utils::txtProgressBar(min = 0, max = n_boot, style = 3, width = 50, char = ".")
  for (i in 1:n_boot) {
    if (verbose) utils::setTxtProgressBar(pb_boot, i)
    boot_indices <- sample(nrow(data_for_stability), size = nrow(data_for_stability), replace = TRUE)
    boot_data <- data_for_stability[boot_indices, ]
    boruta_boot_output <- .run_boruta_selection_module(
      data_for_boruta = boot_data,
      target_id = target_id,
      potential_vars = potential_vars_for_stability,
      species_col_actual_name = species_col_actual_name, 
      settings = settings,
      output_dir_boruta = tempdir(), 
      save_outputs = FALSE,
      apply_max_confirmed_limit = FALSE, 
      verbose = FALSE, 
      max_runs_override = settings$boruta_max_runs_for_stability %||% 50 
    )
    all_bootstrap_selections[[i]] <- boruta_boot_output$confirmed_variables
  }
  if (verbose) close(pb_boot)
  stability_freq_df <- tibble::tibble(Variable = unlist(all_bootstrap_selections)) %>%
    dplyr::count(Variable, name = "Selection_Count") %>%
    dplyr::mutate(Selection_Frequency = Selection_Count / n_boot) %>%
    dplyr::arrange(dplyr::desc(Selection_Frequency), Variable)
  csv_filename <- file.path(output_dir_module, paste0("variable_stability_freq_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".csv"))
  utils::write.csv(stability_freq_df, csv_filename, row.names = FALSE)
  if (verbose) message("    Variable stability frequencies saved to: ", csv_filename)
  return(stability_freq_df)
}

#' Run Chi-squared/Fisher Tests for Categorical Environmental Variables (Internal)
#' @noRd
.run_env_categorical_tests <- function(data_for_tests, target_id, species_id_col, categorical_vars_to_test, verbose = TRUE) {
  if (verbose) message("        Running association tests for categorical env. variables...")
  if (length(categorical_vars_to_test) == 0) return(NULL)
  # Create a binary group: Target vs. All Others in this subset
  data_tests_prep <- data_for_tests %>%
    dplyr::mutate(TestGroup = factor(ifelse(!!dplyr::sym(species_id_col) == target_id, "Target", "Comparison")))
  results_list <- purrr::map_dfr(categorical_vars_to_test, function(var_cat) {
    if (!var_cat %in% names(data_tests_prep)) return(NULL)
    contingency_table_data <- data_tests_prep %>%
      dplyr::select(TestGroup, !!dplyr::sym(var_cat)) %>%
      stats::na.omit() %>%
      dplyr::mutate(!!dplyr::sym(var_cat) := factor(!!dplyr::sym(var_cat)))
    # Create the contingency table *before* trying to use it.
    contingency_table <- table(contingency_table_data$TestGroup, contingency_table_data[[var_cat]])
    # Create a readable summary string of the contingency table for CSV output
    format_row_to_string <- function(row_name, table_row) {
      paste0(row_name, ": [", paste(names(table_row), table_row, sep = "=", collapse = ", "), "]")
    }
    counts_summary_string <- tryCatch({
      target_str <- if ("Target" %in% rownames(contingency_table)) format_row_to_string("Target", contingency_table["Target", ]) else "Target: []"
      comparison_str <- if ("Comparison" %in% rownames(contingency_table)) format_row_to_string("Comparison", contingency_table["Comparison", ]) else "Comparison: []"
      paste(target_str, comparison_str, sep = "; ")
    }, error = function(e) { "Error generating counts summary" })
    if (nrow(contingency_table_data) < 5 || length(unique(contingency_table_data$TestGroup)) < 2 || length(levels(contingency_table_data[[var_cat]])) < 2) {
      return(tibble::tibble(Variable = var_cat, Test_Type = "Skipped", Chi_Sq_Statistic = NA_real_, P_Value = NA_real_,
                            Counts_Summary = counts_summary_string,
                            Notes = "Insufficient data or groups/levels",
                            Contingency_Table = list(contingency_table)))
    }
    test_result <- tryCatch({
      if (any(contingency_table < 5)) {
        fisher_test <- stats::fisher.test(contingency_table, simulate.p.value = TRUE, B=1999) # B for more stable p-value
        list(test_type = "Fisher's Exact Test", statistic = NA_real_, p_value = fisher_test$p.value)
      } else {
        chi_sq_test <- stats::chisq.test(contingency_table)
        list(test_type = "Chi-squared Test", statistic = chi_sq_test$statistic, p_value = chi_sq_test$p.value)
      }
    }, error = function(e) {
      list(test_type = "Error", statistic = NA_real_, p_value = NA_real_, notes = as.character(e$message))
    })
    tibble::tibble(Variable = var_cat, Test_Type = test_result$test_type,
                   Chi_Sq_Statistic = test_result$statistic, P_Value = test_result$p_value,
                   Counts_Summary = counts_summary_string,
                   Notes = test_result$notes %||% NA_character_,
                   Contingency_Table = list(contingency_table))
  })
  return(results_list)
}


#' Run Categorical Association Module (Internal)
#'
#' Orchestrates Apriori for indicator species and association tests for
#' categorical environmental variables.
#'
#' @param data_for_env_tests Tibble. Data for environmental tests (target + comparisons, selected vars).
#' @param target_id Character. ID of the target species.
#' @param species_id_col Character. Name of the species ID column.
#' @param selected_categorical_vars Character vector. Categorical env. vars selected for analysis.
#' @param output_dir_module Character. Path to save CSV outputs.
#' @param settings An object of class 'niche_analysis_settings'.
#' @param verbose Logical.
#'
#' @return A list containing `env_categorical_test_results`.
#' @noRd
.run_categorical_association_module <- function(data_for_env_tests,
                                                target_id,
                                                species_id_col,
                                                selected_categorical_vars,
                                                output_dir_module, settings, verbose = TRUE) {
  if (verbose) message("    Running Categorical Association module...")
  # Environmental Categorical Variable Association Tests
  env_cat_tests_df <- .run_env_categorical_tests(data_for_env_tests, target_id, species_id_col, selected_categorical_vars, verbose)
  if (!is.null(env_cat_tests_df) && nrow(env_cat_tests_df) > 0) {
    df_for_csv <- env_cat_tests_df %>% dplyr::select(-dplyr::any_of("Contingency_Table"))
    utils::write.csv(df_for_csv, file.path(output_dir_module, paste0("env_categorical_tests_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".csv")), row.names = FALSE)
    if (verbose) message("        Environmental categorical variable association test results saved.")
  } else if (verbose && length(selected_categorical_vars) > 0) {
    message("        No results from environmental categorical variable association tests for ", target_id)
  }
  return(list(
    env_categorical_test_results = env_cat_tests_df
  ))
}

# --- Plotting Orchestrator and Helper Functions ---

#' Get N Closest Species (Internal)
#'
#' Identifies the N closest species to the target based on a specified overlap metric.
#'
#' @param target_id Character. The ID of the target species.
#' @param overlap_results Tibble. Output from `.run_overlap_module()`.
#' @param n_closest Integer. Number of closest species to retrieve.
#' @param overlap_metric Character. Name of the overlap metric column in `overlap_results` to use for ranking.
#' @param verbose Logical.
#' @return Character vector of the N closest species IDs (excluding the target itself).
#' @noRd
.get_n_closest_species <- function(target_id, overlap_results, n_closest, overlap_metric, verbose = TRUE) {
  if (is.null(overlap_results) || nrow(overlap_results) == 0 || !overlap_metric %in% names(overlap_results)) {
    if (verbose) message("      Cannot determine N closest species: overlap results missing or metric '", overlap_metric, "' not found.")
    return(character(0))
  }
  
  closest_species <- overlap_results %>%
    dplyr::filter(!!dplyr::sym("target_species") == target_id) %>% 
    dplyr::arrange(dplyr::desc(!!dplyr::sym(overlap_metric))) %>%
    dplyr::slice_head(n = n_closest) %>%
    dplyr::pull(!!dplyr::sym("comparison_species")) 
  return(closest_species)
}

#' Plot Variable Distributions (Continuous and Categorical) (Internal)
#'
#' Generates violin/boxplot for continuous variables and bar plots for categorical variables.
#'
#' @param data_for_plot Tibble. Data containing species ID and selected environmental variables.
#' @param target_id Character. ID of the target species (for context in titles).
#' @param species_to_include Character vector. Species to include in the plot.
#' @param selected_vars Character vector. Variables to plot.
#' @param species_id_col Character. Name of the species ID column.
#' @param niche_data_metadata List. The `metadata` element from `niche_data_object`.
#' @param global_palette Named character vector. Palette mapping species names to colors.
#' @param output_dir Character. Directory to save the plot.
#' @param plot_suffix Character. Suffix for the plot filename (e.g., "_all", "_subsetN").
#' @param settings_for_plot List. The `analysis_settings` object, or relevant parts for plotting.
#' @param verbose Logical.
#' @noRd
.plot_variable_distributions <- function(data_for_plot, target_id, species_to_include,
                                         selected_vars, species_id_col,
                                         niche_data_metadata, global_palette,
                                         output_dir, plot_suffix = "",
                                         settings_for_plot, verbose = TRUE) {
  if (verbose) message("      Generating distribution plots", plot_suffix, " for target: ", target_id)
  
  data_subset <- data_for_plot %>%
    dplyr::filter(!!dplyr::sym(species_id_col) %in% species_to_include) %>%
    dplyr::select(dplyr::all_of(species_id_col), dplyr::all_of(selected_vars))
  if (verbose && settings_for_plot$plot_categorical_distributions && length(intersect(selected_vars, niche_data_metadata$env_vars$categorical)) == 0) {
    message("        NOTE: Plotting of categorical distributions is enabled, but no categorical variables were found within the 'selected_vars' for this plot group (", plot_suffix, ").")
    message("              Selected vars for this plot: ", paste(selected_vars, collapse=", "))
    message("              Available categoricals in metadata: ", paste(niche_data_metadata$env_vars$categorical, collapse=", "))
  }
  if (nrow(data_subset) == 0) {
    if (verbose) message("        No data available for distribution plots after filtering for species: ", paste(species_to_include, collapse=", "))
    return()
  }
  ordered_species_levels <- intersect(c(target_id, sort(setdiff(species_to_include, target_id))), names(global_palette))
  data_subset[[species_id_col]] <- factor(data_subset[[species_id_col]], levels = ordered_species_levels)
  palette_for_plot <- global_palette[names(global_palette) %in% levels(data_subset[[species_id_col]])]
  continuous_vars_plot <- intersect(selected_vars, niche_data_metadata$env_vars$continuous)
  categorical_vars_plot <- intersect(selected_vars, niche_data_metadata$env_vars$categorical)
  # Continuous variables plot
  if (length(continuous_vars_plot) > 0) {
    df_plot_cont <- data_subset %>%
      dplyr::select(dplyr::all_of(species_id_col), dplyr::all_of(continuous_vars_plot)) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(continuous_vars_plot), names_to = "env_var_name", values_to = "value") %>% 
      dplyr::mutate(value = as.numeric(value)) %>%
      stats::na.omit()
    if (nrow(df_plot_cont) > 0) {
      p_cont <- NULL      
      if (settings_for_plot$continuous_dist_plot_type == "ridgeline" && requireNamespace("ggridges", quietly = TRUE)) {
        current_plot_width_cont <- if (plot_suffix == "_all_comparison") 4 else 5 
        p_cont <- ggplot2::ggplot(df_plot_cont, ggplot2::aes(x = value, y = !!dplyr::sym(species_id_col), fill = !!dplyr::sym(species_id_col), color = !!dplyr::sym(species_id_col))) +
          ggridges::geom_density_ridges(alpha = 0.7, scale = 1.2, rel_min_height = 0.01,
                                        jittered_points = settings_for_plot$ridgeline_show_points %||% TRUE, 
                                        point_shape = '|', point_size = 1.5,                                        
                                        position = ggridges::position_points_jitter(height = 0)) +
          ggplot2::facet_wrap(~ env_var_name, scales = "free_x", ncol = 1, 
                              labeller = ggplot2::as_labeller(niche_data_metadata$labels$continuous_vars)) +
          ggplot2::scale_fill_manual(values = palette_for_plot, guide = "none") +
          ggplot2::scale_color_manual(values = palette_for_plot, guide = "none") + 
          ggplot2::labs(title = paste("Continuous Variable Distributions (Ridgeline)", plot_suffix),
                        subtitle = paste("Target:", target_id), x = "Value", y = "Species") +
          ggplot2::theme_minimal(base_size = 10) + 
          ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey92"),
                         axis.text.y = ggplot2::element_text(size = settings_for_plot$ridgeline_y_axis_text_size %||% 7))
        num_vars_cont_plot <- dplyr::n_distinct(df_plot_cont$env_var_name) 
        num_species_this_plot <- dplyr::n_distinct(df_plot_cont[[species_id_col]])
        
        height_per_species_factor <- settings_for_plot$ridgeline_height_per_species %||% 0.15 
        base_height_per_facet <- settings_for_plot$ridgeline_base_height_per_facet %||% 1.2 
        min_total_height_ridge <- settings_for_plot$ridgeline_min_total_height %||% 7
        max_total_height_ridge <- settings_for_plot$ridgeline_max_total_height %||% 40
        
        calculated_height <- num_vars_cont_plot * (base_height_per_facet + num_species_this_plot * height_per_species_factor)
        plot_height_cont <- min(max(min_total_height_ridge, calculated_height), max_total_height_ridge)
        
      } else { 
        if (settings_for_plot$continuous_dist_plot_type == "ridgeline") {
          if(verbose) warning("Package 'ggridges' not installed. Defaulting to violin plots for continuous variables.")
        }
        current_plot_width_cont <- 15 
        p_cont <- ggplot2::ggplot(df_plot_cont, ggplot2::aes(x = !!dplyr::sym(species_id_col), y = value, fill = !!dplyr::sym(species_id_col))) +
          ggplot2::geom_violin(trim = FALSE, alpha = 0.7, scale = "width", show.legend = FALSE) +
          ggplot2::geom_boxplot(width = 0.15, fill = "white", alpha = 0.6, outlier.shape = NA, coef = 0, show.legend = FALSE) +
          ggplot2::stat_summary(fun = stats::median, geom = "point", shape = 18, size = 2.5, color = "black", show.legend = FALSE) +
          ggplot2::facet_wrap(~ env_var_name, scales = "free_y", ncol = min(3, length(continuous_vars_plot)),
                              labeller = ggplot2::as_labeller(niche_data_metadata$labels$continuous_vars)) +
          ggplot2::scale_fill_manual(values = palette_for_plot) +
          ggplot2::labs(title = paste("Continuous Variable Distributions (Violin)", plot_suffix),
                        subtitle = paste("Target:", target_id), x = "Species", y = "Value") +
          ggplot2::theme_bw(base_size = 10) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
                         legend.position = "none",
                         strip.background = ggplot2::element_rect(fill = "grey92"),
                         panel.spacing = ggplot2::unit(1, "lines"))
        num_facet_rows_cont_violin <- ceiling(length(continuous_vars_plot) / min(3, length(continuous_vars_plot)))
        plot_height_cont <- max(7, num_facet_rows_cont_violin * 4) 
      }
      ggplot2::ggsave(filename = file.path(output_dir, paste0("distributions_continuous", plot_suffix, "_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".png")),
                      plot = p_cont, width = current_plot_width_cont, height = plot_height_cont, dpi = 300, limitsize = FALSE)
    }
  }
  if (settings_for_plot$plot_categorical_distributions && length(categorical_vars_plot) > 0) {
    if (verbose) {
      message("        Attempting to plot categorical distributions for: ", paste(categorical_vars_plot, collapse = ", "))
      message("          Number of rows in data_subset for categorical plotting: ", nrow(data_subset))
    }
    df_plot_cat <- data_subset %>%
      dplyr::select(dplyr::all_of(species_id_col), dplyr::all_of(categorical_vars_plot)) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(categorical_vars_plot), names_to = "env_var_name", values_to = "value") %>%
      dplyr::mutate(value = factor(value)) %>% 
      stats::na.omit()
    df_plot_cat_labeled <- df_plot_cat %>%
      dplyr::group_by(env_var_name) %>%
      dplyr::mutate(value_labeled = factor(value,
                                           levels = levels(value), 
                                           labels = niche_data_metadata$labels$categorical_vars[[unique(env_var_name)]][levels(value)] %||% levels(value)
      )) %>%
      dplyr::ungroup()
    
    if (nrow(df_plot_cat) == 0 && verbose) {
      message("        No data for categorical distributions after pivoting/NA removal for target: ", target_id, plot_suffix)
      return()
    }
    max_levels_cat <- df_plot_cat_labeled %>%
      dplyr::group_by(env_var_name) %>%
      dplyr::summarise(n_levels = dplyr::n_distinct(value_labeled), .groups = 'drop') %>%
      dplyr::pull(n_levels) %>%
      max(na.rm = TRUE)
    num_colors_needed_cat <- max(3, max_levels_cat, na.rm = TRUE)
    palette_name_cat_fill <- if (num_colors_needed_cat > 8 && num_colors_needed_cat <= 12) "Set3" else if (num_colors_needed_cat > 12) "Paired" else "Set2"
    category_colors_fill <- tryCatch({
      RColorBrewer::brewer.pal(num_colors_needed_cat, palette_name_cat_fill)[1:max_levels_cat]
    }, warning = function(w){ randomcoloR::distinctColorPalette(max_levels_cat) }, error = function(e){ randomcoloR::distinctColorPalette(max_levels_cat) })
    if (nrow(df_plot_cat_labeled) > 0) {
      p_cat <- ggplot2::ggplot(df_plot_cat_labeled,
                               ggplot2::aes(x = value_labeled, 
                                            fill = !!dplyr::sym(species_id_col), 
                                            group = !!dplyr::sym(species_id_col))) + 
        ggplot2::geom_bar(position = "stack", alpha = 0.8, show.legend = TRUE) + 
        ggplot2::facet_wrap(~ env_var_name, scales = "free_x", ncol = min(3, length(categorical_vars_plot)),
                            labeller = ggplot2::as_labeller(niche_data_metadata$labels$categorical_vars %||% niche_data_metadata$labels$continuous_vars)) + 
        ggplot2::scale_fill_manual(values = palette_for_plot, name = "Species") + 
        ggplot2::labs(title = paste("Categorical Variable Distributions", plot_suffix),
                      subtitle = paste("Target:", target_id),
                      x = "Category Level", 
                      y = "Number of Occurrences",
                      fill = "Species") + 
        ggplot2::theme_bw(base_size = 10) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
                       legend.position = "bottom",
                       strip.background = ggplot2::element_rect(fill = "grey92"),
                       panel.spacing = ggplot2::unit(1, "lines"))
      plot_height_cat <- if (plot_suffix == "_all_comparison") max(14, ceiling(length(categorical_vars_plot) / min(3, length(categorical_vars_plot))) * 5.5) else max(7, ceiling(length(categorical_vars_plot) / 3) * 4)
      ggplot2::ggsave(filename = file.path(output_dir, paste0("distributions_categorical", plot_suffix, "_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".png")),
                      plot = p_cat, width = 14, height = plot_height_cat, dpi = 300, limitsize = FALSE)
    }
    if (verbose) message("      Distribution plots", plot_suffix, " saved for target: ", target_id)
  }
}

#' Plot Unified Categorical Distributions for All Species in the Dataset (Internal)
#'
#' Generates a single bar plot showing categorical variable distributions for all
#' species present in the dataset. This function is called once after the main analysis loop.
#'
#' @param data_for_plot Tibble. The full `niche_data_object$data`.
#' @param species_to_include Character vector. All species IDs to include in the plot.
#' @param categorical_vars_to_plot Character vector. All categorical variables to plot.
#' @param species_id_col Character. Name of the species ID column.
#' @param niche_data_metadata List. The `metadata` element from `niche_data_object`.
#' @param global_palette Named character vector. Palette mapping all species names to colors.
#' @param output_dir Character. The main output directory (not target-specific).
#' @param verbose Logical.
#' @noRd
.plot_unified_all_species_distributions <- function(data_for_plot,
                                                    species_to_include,
                                                    categorical_vars_to_plot,
                                                    species_id_col,
                                                    niche_data_metadata,
                                                    global_palette,
                                                    output_dir,
                                                    verbose = TRUE) {
  if (verbose) message("\nGenerating unified distribution plot for all species...")

  if (is.null(categorical_vars_to_plot) || length(categorical_vars_to_plot) == 0) {
    if (verbose) message("  Skipping unified plot: No categorical variables found in the dataset.")
    return()
  }

  data_subset <- data_for_plot %>%
    dplyr::filter(!!dplyr::sym(species_id_col) %in% species_to_include) %>%
    dplyr::select(dplyr::all_of(species_id_col), dplyr::all_of(categorical_vars_to_plot))

  if (nrow(data_subset) == 0) {
    if (verbose) message("  Skipping unified plot: No data available for the specified species.")
    return()
  }

  # Prepare data for plotting
  df_plot_cat <- data_subset %>%
    tidyr::pivot_longer(cols = dplyr::all_of(categorical_vars_to_plot), names_to = "env_var_name", values_to = "value") %>%
    dplyr::mutate(value = factor(value)) %>%
    stats::na.omit()

  # Apply pretty labels
  df_plot_cat_labeled <- df_plot_cat %>%
    dplyr::group_by(env_var_name) %>%
    dplyr::mutate(value_labeled = factor(value,
                                         levels = levels(value),
                                         labels = niche_data_metadata$labels$categorical_vars[[unique(env_var_name)]][levels(value)] %||% levels(value)
    )) %>%
    dplyr::ungroup()

  if (nrow(df_plot_cat_labeled) == 0) {
    if (verbose) message("  Skipping unified plot: No data remaining after preparing for plotting.")
    return()
  }

  palette_for_plot <- global_palette[names(global_palette) %in% species_to_include]

  p_cat_unified <- ggplot2::ggplot(df_plot_cat_labeled,
                                   ggplot2::aes(x = value_labeled, fill = !!dplyr::sym(species_id_col))) +
    ggplot2::geom_bar(position = "stack", alpha = 0.8, show.legend = TRUE) +
    ggplot2::facet_wrap(~ env_var_name, scales = "free_x", ncol = min(3, length(categorical_vars_to_plot)),
                        labeller = ggplot2::as_labeller(niche_data_metadata$labels$categorical_vars)) +
    ggplot2::scale_fill_manual(values = palette_for_plot, name = "Species") +
    ggplot2::labs(title = "Unified Categorical Variable Distributions for All Species",
                  x = "Category Level",
                  y = "Number of Occurrences") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
                   legend.position = "bottom",
                   strip.background = ggplot2::element_rect(fill = "grey92"),
                   panel.spacing = ggplot2::unit(1, "lines"))

  plot_height_cat <- max(14, ceiling(length(categorical_vars_to_plot) / min(3, length(categorical_vars_to_plot))) * 5.5)
  plot_filename <- file.path(output_dir, "distributions_categorical_unified_all_species.png")
  ggplot2::ggsave(filename = plot_filename, plot = p_cat_unified, width = 14, height = plot_height_cat, dpi = 300, limitsize = FALSE)
  if (verbose) message("  Unified all-species distribution plot saved to: ", plot_filename)
}

#' Generate All Plots for a Target Species (Internal Orchestrator)
#'
#' Calls specific plotting functions based on available results and settings.
#'
#' @param results_for_target List. Aggregated results for the current target species.
#' @param data_for_plots Tibble. Data subset for plotting (target + all comparisons, selected vars).
#' @param niche_data_object_metadata List. The `metadata` element from `niche_data_object`.
#' @param global_species_palette Named character vector. Palette mapping all species names to colors.
#' @param settings An object of class 'niche_analysis_settings'.
#' @param output_dir_target_specific Character. Path to the target's output directory.
#' @param verbose Logical.
#' @noRd
.generate_all_plots_for_target <- function(results_for_target, data_for_plots,
                                           niche_data_object_metadata, global_species_palette,
                                           settings, output_dir_target_specific, verbose = TRUE) {
  
  if (!settings$generate_plots) return()
  if (verbose) message("    Orchestrating plot generation for target: ", results_for_target$metadata$target_species_id)
  
  target_id <- results_for_target$metadata$target_species_id
  all_comparison_species_ids <- results_for_target$metadata$comparison_species_used
  selected_vars_final <- results_for_target$metadata$final_variables_used
  species_id_col <- niche_data_object_metadata$species_id_col
  # --- Plot distributions for ALL comparison species ---
  tryCatch({ 
    .plot_variable_distributions(
      data_for_plot = data_for_plots, target_id = target_id,
      species_to_include = c(target_id, all_comparison_species_ids),
      selected_vars = selected_vars_final, species_id_col = species_id_col, 
      global_palette = global_species_palette,      
      niche_data_metadata = niche_data_object_metadata, 
      output_dir = output_dir_target_specific, plot_suffix = "_all_comparison",
      settings_for_plot = settings, verbose = verbose 
    )
  }, error = function(e) {
    if(verbose) warning("    WARNING: Distribution plot for _all_comparison failed for target '", target_id, "'. Error: ", e$message, ". Continuing analysis.")
  })
  
  # --- Target-Specific Dendrogram ---
  if (!is.null(results_for_target$hypervolume_analysis) && !is.null(results_for_target$hypervolume_analysis$vars_used) && length(results_for_target$hypervolume_analysis$vars_used) >= 2) {
    tryCatch({
      .run_dendrogram_target_specific_module(
        data_for_module = data_for_plots, 
        target_id = target_id,
        comparison_ids = all_comparison_species_ids,
        vars_for_hv = results_for_target$hypervolume_analysis$vars_used,
        species_id_col = species_id_col,
        settings = settings,
        output_dir_target_specific = output_dir_target_specific,
        global_species_palette = global_species_palette,
        verbose = verbose
      )
    }, error = function(e) {
      if(verbose) warning("    WARNING: Target-specific dendrogram failed for target '", target_id, "'. Error: ", e$message, ". Continuing analysis.")
    })
  }
  # --- Handle N Closest Species Plots ---
  if (settings$plot_type == "n_closest_species" && !is.null(settings$plot_n_closest) && settings$plot_n_closest > 0) {
    if (verbose) {
      message("    Attempting to generate plots for N closest species.")
      message("      plot_n_closest: ", settings$plot_n_closest)
      message("      plot_n_closest_overlap_metric: ", settings$plot_n_closest_overlap_metric)
    }
    n_closest_ids <- .get_n_closest_species(
      target_id = target_id,
      overlap_results = results_for_target$overlap_indices,
      n_closest = settings$plot_n_closest,
      overlap_metric = settings$plot_n_closest_overlap_metric,
      verbose = verbose
    )
    if (verbose) message("      Species identified by .get_n_closest_species: ", paste(n_closest_ids, collapse=", "))
    if (length(n_closest_ids) > 0) {
      species_for_subset_plot <- unique(c(target_id, n_closest_ids))
      if (verbose) message("      Final species for subset plots (target + N closest): ", paste(species_for_subset_plot, collapse=", "))
      suffix_subset <- paste0("_subset", settings$plot_n_closest)
      tryCatch({ 
        .plot_variable_distributions(
          data_for_plot = data_for_plots, target_id = target_id,
          species_to_include = species_for_subset_plot,
          selected_vars = selected_vars_final, species_id_col = species_id_col,
          niche_data_metadata = niche_data_object_metadata, global_palette = global_species_palette,
          output_dir = output_dir_target_specific, plot_suffix = suffix_subset,
          settings_for_plot = settings, verbose = verbose
        )
      }, error = function(e) {
        if(verbose) warning("    WARNING: Distribution plot for ", suffix_subset, " failed for target '", target_id, "'. Error: ", e$message, ". Continuing analysis.")
      })
      # PCA Biplot for SUBSET
      data_pca_subset <- data_for_plots %>%
        dplyr::filter(!!dplyr::sym(species_id_col) %in% species_for_subset_plot) %>%
        dplyr::select(dplyr::all_of(species_id_col), dplyr::all_of(intersect(selected_vars_final, niche_data_object_metadata$env_vars$continuous))) %>%
        stats::na.omit()
      if(nrow(data_pca_subset) > ncol(data_pca_subset) -1 && length(unique(data_pca_subset[[species_id_col]])) > 1){
        tryCatch({ 
          pca_res_subset <- stats::prcomp(data_pca_subset[, -1], scale. = TRUE, center = TRUE)
          .plot_pca_biplot(pca_results = pca_res_subset, data_for_pca_plot = data_pca_subset,
                           target_id = target_id, species_id_col = species_id_col,
                           niche_data_metadata = niche_data_object_metadata, global_palette = global_species_palette,
                           output_dir = output_dir_target_specific, plot_suffix = suffix_subset, verbose = verbose)
        }, error = function(e) {
          if(verbose) warning("    WARNING: PCA biplot for ", suffix_subset, " failed for target '", target_id, "'. Error: ", e$message, ". Continuing analysis.")
        })
      } else { 
        if(verbose) message("      Skipping PCA plot for ", suffix_subset, ": insufficient data or groups after subsetting.")
      }
      # LDA Plot for SUBSET
      if (!is.null(results_for_target$lda_multigroup)) {
        data_lda_subset <- data_for_plots %>% dplyr::filter(!!dplyr::sym(species_id_col) %in% species_for_subset_plot)
        if(!is.null(results_for_target$lda_multigroup$ld_scores)){
          tryCatch({ 
            .plot_lda_results(lda_results_module = results_for_target$lda_multigroup,
                              data_for_lda_plot = data_lda_subset,
                              target_id = target_id, species_id_col = species_id_col,
                              global_palette = global_species_palette, 
                              output_dir = output_dir_target_specific,
                              plot_suffix = suffix_subset, verbose = verbose)
          }, error = function(e) {
            if(verbose) warning("    WARNING: LDA plot for ", suffix_subset, " failed for target '", target_id, "'. Error: ", e$message, ". Continuing analysis.")
          })
        }
      } else {
        if(verbose) message("      Skipping LDA plot for ", suffix_subset, ": LDA results not available.")
      }
      # --- Pairwise PERMANOVA Plots for SUBSET ---
      if (!is.null(results_for_target$permanova_pairwise)) {
        permanova_data_subset <- results_for_target$permanova_pairwise %>%
          dplyr::filter(specie1 == target_id & specie2 %in% n_closest_ids)
        
        if (nrow(permanova_data_subset) > 0) {
          tryCatch({
            .plot_permanova_pairwise_r2(
              permanova_results_df = permanova_data_subset,
              target_id = target_id, output_dir = output_dir_target_specific,
              plot_suffix = suffix_subset, verbose = verbose
            )
            .plot_permanova_pairwise_pvalue(
              permanova_results_df = permanova_data_subset,
              target_id = target_id, output_dir = output_dir_target_specific,
              plot_suffix = suffix_subset, verbose = verbose
            )
          }, error = function(e) {
            if(verbose) warning("    WARNING: Pairwise PERMANOVA plots for ", suffix_subset, " failed for target '", target_id, "'. Error: ", e$message, ". Continuing analysis.")
          })
        } else {
          if(verbose) message("      Skipping Pairwise PERMANOVA plots for ", suffix_subset, ": no PERMANOVA results for the N closest species.")
        }
      } else {
        if(verbose) message("      Skipping Pairwise PERMANOVA plots for ", suffix_subset, ": permanova_pairwise results not available.")
      }
      if (!is.null(results_for_target$modularity_network) && !is.null(results_for_target$modularity_network$network)) {
        if (verbose) message("      Attempting to generate Modularity network plot for ", suffix_subset)
        tryCatch({
          net_global <- results_for_target$modularity_network$network
          groups_global <- results_for_target$modularity_network$groups
          mod_score_global <- results_for_target$modularity_network$modularity_score
          valid_species_for_subset_plot <- intersect(species_for_subset_plot, igraph::V(net_global)$name)
          if (length(valid_species_for_subset_plot) >= 2) { 
            net_subset <- igraph::induced_subgraph(net_global, V(net_global)$name %in% valid_species_for_subset_plot)
            
            groups_subset <- NULL
            if (!is.null(groups_global) && length(groups_global) > 0 && !is.null(names(groups_global))) {
              groups_subset <- groups_global[names(groups_global) %in% igraph::V(net_subset)$name]
              if(length(groups_subset) == 0 || is.null(names(groups_subset)) || !all(names(groups_subset) %in% igraph::V(net_subset)$name) ) {
                groups_subset <- NULL 
              }
            }
            
            modularity_results_subset <- list(
              network = net_subset,
              modularity_score = mod_score_global, 
              groups = groups_subset,
              n_groups = if(!is.null(groups_subset) && length(groups_subset) > 0) length(unique(groups_subset)) else NA_integer_
            )
            
            .plot_modularity_network(
              modularity_results = modularity_results_subset,
              target_id = target_id,
              overlap_threshold = settings$modularity_network_overlap_threshold %||% 0.25, # Use fallback
              overlap_metric_name = settings$modularity_network_metric %||% "mean_schoener_D", # Use the correct metric from settings
              global_palette = global_species_palette,
              output_dir = output_dir_target_specific,
              plot_suffix = suffix_subset, 
              file_path_override = NULL,
              verbose = verbose
            )
          } else {
            if (verbose) message("        Skipping Modularity network plot for ", suffix_subset, ": less than 2 valid species from subset found in the global network or global network is too small.")
          }
        }, error = function(e) {
          if(verbose) warning("    WARNING: Modularity network plot for ", suffix_subset, " failed for target '", target_id, "'. Error: ", e$message, ". Details: ", paste(capture.output(print(e)), collapse = "\n"))
        })
      } else {
        if(verbose) message("      Skipping Modularity network plot for ", suffix_subset, ": Modularity results or global network not available.")
      }
      if (!is.null(results_for_target$hypervolume_analysis) && settings$hypervolume_n_vars_plot > 0) {
        metric_for_hv_ranking <- settings$hypervolume_plot_overlap_metric %||% "sorensen"
        
        # Use HV overlap stats if available and the metric exists, otherwise fallback to n_closest_ids
        if (!is.null(results_for_target$hypervolume_analysis$overlap_stats) &&
            metric_for_hv_ranking %in% names(results_for_target$hypervolume_analysis$overlap_stats)) {
          
          if(verbose) message("      Selecting species for HV plot based on '", metric_for_hv_ranking, "' hypervolume overlap.")
          
          closest_species_by_hv <- results_for_target$hypervolume_analysis$overlap_stats %>%
            dplyr::filter(target_species == target_id) %>%
            dplyr::arrange(dplyr::desc(!!dplyr::sym(metric_for_hv_ranking))) %>%
            dplyr::slice_head(n = (settings$hypervolume_n_vars_plot %||% 3) - 1) %>%
            dplyr::pull(comparison_species)
          
          species_for_hv_plot_subset <- unique(c(target_id, closest_species_by_hv))
          
        } else {
          if(verbose) warning("      Could not use HV overlap for plot species selection (data or metric '", metric_for_hv_ranking, "' missing). Falling back to 'n_closest_ids' (based on Schoener's D).")
          species_for_hv_plot_subset <- head(n_closest_ids, (settings$hypervolume_n_vars_plot %||% 3) - 1)
          species_for_hv_plot_subset <- unique(c(target_id, species_for_hv_plot_subset))
        }
        
        
        
        tryCatch({ 
          .plot_hypervolumes_on_pca(
            data_for_module = data_for_plots,
            target_id = target_id,
            species_for_hv_plot = species_for_hv_plot_subset,
            species_id_col = species_id_col,
            vars_for_hv_pca = results_for_target$hypervolume_analysis$vars_used, 
            settings = settings, 
            global_palette = global_species_palette,
            output_dir = output_dir_target_specific,
            plot_suffix = paste0("_hv_pca_top", length(species_for_hv_plot_subset)),
            verbose = verbose
          )
        }, error = function(e) {
          if(verbose) warning("    WARNING: Hypervolume PCA plot for subset failed for target '", target_id, "'. Error: ", e$message, ". Continuing analysis.")
        })
      } else { 
        if(verbose) message("      Skipping Hypervolume PCA plot for subset: Hypervolume results not available or hypervolume_n_vars_plot <= 0.")
      }
      # Hypervolume on ORIGINAL VARS for SUBSET
      if (!is.null(results_for_target$hypervolume_analysis) &&
          !is.null(results_for_target$hypervolume_analysis$hypervolume_objects) && 
          settings$hypervolume_n_vars_plot > 0) { 
        metric_for_hv_ranking <- settings$hypervolume_plot_overlap_metric %||% "sorensen"
        
        # Use HV overlap stats if available, otherwise fallback to n_closest_ids
        if (!is.null(results_for_target$hypervolume_analysis$overlap_stats) &&
            metric_for_hv_ranking %in% names(results_for_target$hypervolume_analysis$overlap_stats)) {
          
          if(verbose) message("      Selecting species for Original HV plot based on '", metric_for_hv_ranking, "' hypervolume overlap.")
          
          closest_species_by_hv_orig <- results_for_target$hypervolume_analysis$overlap_stats %>%
            dplyr::filter(target_species == target_id) %>%
            dplyr::arrange(dplyr::desc(!!dplyr::sym(metric_for_hv_ranking))) %>%
            dplyr::slice_head(n = (settings$hypervolume_n_vars_plot %||% 3) - 1) %>%
            dplyr::pull(comparison_species)
          
          species_for_hv_orig_subset_plot <- unique(c(target_id, closest_species_by_hv_orig))
          
        } else {
          if(verbose) warning("      Could not use HV overlap for original var plot species selection (data or metric '", metric_for_hv_ranking, "' missing). Falling back to 'n_closest_ids' (based on Schoener's D).")
          species_for_hv_orig_subset_plot <- head(n_closest_ids, settings$hypervolume_n_vars_plot - 1)
          species_for_hv_orig_subset_plot <- unique(c(target_id, species_for_hv_orig_subset_plot))
        }
        
        tryCatch({
          .plot_hypervolumes_original_vars(
            hv_objects_list = results_for_target$hypervolume_analysis$hypervolume_objects,
            target_id = target_id,
            species_to_plot = species_for_hv_orig_subset_plot,
            vars_used_for_hv = results_for_target$hypervolume_analysis$vars_used,
            global_palette = global_species_palette,
            output_dir = output_dir_target_specific,
            plot_suffix = paste0("_hv_orig_top", length(species_for_hv_orig_subset_plot)),
            settings = settings,
            verbose = verbose
          )
        }, error = function(e) {
          if(verbose) warning("    WARNING: Hypervolume plot on original vars for subset failed for target '", target_id, "'. Error: ", e$message, ". Continuing analysis.")
        })
      }
      
    } else { 
      if(verbose) message("      Skipping subset plots: No closest species identified or n_closest_ids is empty.")
    }
  } else { 
    if(verbose) message("      Skipping N-closest-species plots: plot_type is not 'n_closest_species' or plot_n_closest is not set appropriately.")
  }
  # --- Pairwise PERMANOVA Plots (R2 and P-value) for ALL comparison species ---
  # Only generate these if we are NOT in "n_closest_species" plot mode with a valid n.
  if (!(settings$plot_type == "n_closest_species" && 
        !is.null(settings$plot_n_closest) && 
        settings$plot_n_closest > 0)) {
    
    if (!is.null(results_for_target$permanova_pairwise)) {
      tryCatch({
        .plot_permanova_pairwise_r2(
          permanova_results_df = results_for_target$permanova_pairwise,
          target_id = target_id, output_dir = output_dir_target_specific,
          plot_suffix = "_all_comparison", verbose = verbose
        )
        .plot_permanova_pairwise_pvalue(
          permanova_results_df = results_for_target$permanova_pairwise,
          target_id = target_id, output_dir = output_dir_target_specific,
          plot_suffix = "_all_comparison", verbose = verbose
        )
      }, error = function(e) {
        if(verbose) warning("    WARNING: Pairwise PERMANOVA plots for _all_comparison failed for target '", target_id, "'. Error: ", e$message, ". Continuing analysis.")
      })
    }
  } else {
    if (verbose) message("      Skipping _all_comparison PERMANOVA plots because N-closest species plots are active and configured.")
  }
  data_pca_all <- data_for_plots %>%
    dplyr::select(dplyr::all_of(species_id_col), dplyr::all_of(intersect(selected_vars_final, niche_data_object_metadata$env_vars$continuous))) %>%
    stats::na.omit()
  if(nrow(data_pca_all) > ncol(data_pca_all) -1 && length(unique(data_pca_all[[species_id_col]])) > 1){
    tryCatch({
      pca_res_all <- stats::prcomp(data_pca_all[, -1], scale. = TRUE, center = TRUE)
      .plot_pca_biplot(pca_results = pca_res_all, data_for_pca_plot = data_pca_all,
                       target_id = target_id, species_id_col = species_id_col,
                       niche_data_metadata = niche_data_object_metadata, global_palette = global_species_palette,
                       output_dir = output_dir_target_specific, plot_suffix = "_all_comparison", verbose = verbose)
    }, error = function(e) {
      if(verbose) warning("    WARNING: PCA biplot for _all_comparison failed for target '", target_id, "'. Error: ", e$message, ". Continuing analysis.")
    })
  } else {
    if(verbose) message("      Skipping PCA plot for all comparison species: insufficient data or groups.")
    
  }
  # --- LDA Plot for ALL comparison species ---
  if (!is.null(results_for_target$lda_multigroup)) { 
    tryCatch({
      .plot_lda_results(lda_results_module = results_for_target$lda_multigroup,
                        data_for_lda_plot = data_for_plots,
                        target_id = target_id, species_id_col = species_id_col,
                        global_palette = global_species_palette, output_dir = output_dir_target_specific,
                        plot_suffix = "_all_comparison", verbose = verbose)
    }, error = function(e) {
      if(verbose) warning("    WARNING: LDA plot for _all_comparison failed for target '", target_id, "'. Error: ", e$message, ". Continuing analysis.")
    })
  } else { 
    if(verbose) message("      Skipping LDA plot for _all_comparison: LDA results not available.")
  }
  # --- Modularity Network for ALL comparison species ---
  if (!is.null(results_for_target$modularity_network)) {
    tryCatch({ 
      .plot_modularity_network(
        modularity_results = results_for_target$modularity_network,
        target_id = target_id,
        overlap_threshold = settings$modularity_network_overlap_threshold %||% 0.25, # Use fallback
        overlap_metric_name = settings$modularity_network_metric %||% "mean_schoener_D", # Use the correct metric from settings
        global_palette = global_species_palette, 
        output_dir = output_dir_target_specific, 
        plot_suffix = "_all_comparison", 
        file_path_override = NULL, 
        verbose = verbose
      )
    }, error = function(e) {
      if(verbose) warning("    WARNING: Modularity network plot for _all_comparison failed for target '", target_id, "'. Error: ", e$message, ". Continuing analysis. Details: ", paste(capture.output(print(e)), collapse = "\n"))
    })
  } else { 
    if(verbose) message("      Skipping Modularity network plot for _all_comparison: Modularity results not available.")
  }
  if (!is.null(results_for_target$hypervolume_analysis) &&
      !is.null(results_for_target$hypervolume_analysis$hypervolume_objects)) {
    species_for_hv_orig_all_plot <- names(results_for_target$hypervolume_analysis$hypervolume_objects)
    if (length(species_for_hv_orig_all_plot) > 0) {
      tryCatch({
        .plot_hypervolumes_original_vars(
          hv_objects_list = results_for_target$hypervolume_analysis$hypervolume_objects,
          target_id = target_id,
          species_to_plot = species_for_hv_orig_all_plot,
          vars_used_for_hv = results_for_target$hypervolume_analysis$vars_used,
          global_palette = global_species_palette,
          output_dir = output_dir_target_specific,
          plot_suffix = "_hv_orig_all_comparison",
          settings = settings,
          verbose = verbose
        )
      }, error = function(e) {
        if(verbose) warning("    WARNING: Hypervolume plot on original vars for _all_comparison failed for target '", target_id, "'. Error: ", e$message, ". Continuing analysis.")
      })
    }
  }
  if (verbose) message("    Plot generation orchestration complete for: ", target_id)
}

#' Plot Hypervolumes on Original Variables (Internal)
#'
#' Plots hypervolumes (previously constructed on original environmental variables)
#' for a specified set of species.
#'
#' @param hv_objects_list Named list. A list of `Hypervolume` objects, where names are species IDs.
#'        (e.g., from `results_for_target$hypervolume_analysis$hypervolume_objects`).
#' @param target_id Character. ID of the target species.
#' @param species_to_plot Character vector. Species to include in this plot.
#' @param vars_used_for_hv Character vector. Names of original environmental variables
#'        used to construct these hypervolumes (for axis labels).
#' @param global_palette Named character vector. Palette mapping species names to colors.
#' @param output_dir Character. Directory to save the plot.
#' @param plot_suffix Character. Suffix for the plot filename.
#' @param settings An object of class 'niche_analysis_settings'.
#' @param verbose Logical.
#' @noRd
.plot_hypervolumes_original_vars <- function(hv_objects_list,
                                             target_id,
                                             species_to_plot,
                                             vars_used_for_hv,
                                             global_palette,
                                             output_dir,
                                             plot_suffix = "",
                                             settings,
                                             verbose = TRUE) {
  
  if (verbose) message("      Generating hypervolume plot on original variables", plot_suffix, " for target: ", target_id)
  
  if (is.null(hv_objects_list) || length(hv_objects_list) == 0) {
    if (verbose) message("        Original HV Plot: Hypervolume objects list is NULL or empty. Skipping.")
    return()
  }
  if (length(species_to_plot) < 1) { 
    if (verbose) message("        Original HV Plot: No species specified to plot. Skipping.")
    return()
  }
  if (length(vars_used_for_hv) < 2) {
    if (verbose) message("        Original HV Plot: Less than 2 variables were used for HV construction. Cannot create 2D+ plot. Skipping.")
    return()
  }
  # Filter hypervolumes for the species to plot
  hvs_to_plot_actual <- hv_objects_list[names(hv_objects_list) %in% species_to_plot]
  hvs_to_plot_actual <- hvs_to_plot_actual[!sapply(hvs_to_plot_actual, is.null)] 
  
  if (length(hvs_to_plot_actual) == 0) {
    if (verbose) message("        Original HV Plot: No valid hypervolumes found for the selected species: ", paste(species_to_plot, collapse=", "), ". Skipping.")
    return()
  }
  # Determine number of dimensions to display
  dims_to_show_in_plot <- min(settings$hypervolume_n_vars_plot %||% 2, length(vars_used_for_hv), hypervolume::hypervolume_n_dimensions(hvs_to_plot_actual[[1]]))
  dims_to_show_in_plot <- min(settings$hypervolume_n_vars_plot %||% 2, length(vars_used_for_hv), hypervolume::hypervolume_get_dimensionality(hvs_to_plot_actual[[1]]))
  dims_to_show_in_plot <- min(dims_to_show_in_plot, length(vars_used_for_hv))
  if (dims_to_show_in_plot < 2 && length(vars_used_for_hv) >=2 && ncol(hvs_to_plot_actual[[1]]@RandomPoints) >=2 ) dims_to_show_in_plot <- 2
  if (dims_to_show_in_plot < 2) {
    if (verbose) message("        Original HV Plot: Cannot plot with less than 2 dimensions. Skipping.")
    return()
  }
  # Join hypervolumes into a HypervolumeList for plotting
  hv_list_obj_for_plot <- tryCatch(
    hypervolume::hypervolume_join(hvs_to_plot_actual),
    error = function(e) {
      if (verbose) warning("        Original HV Plot: Failed to join hypervolumes for plotting: ", e$message)
      NULL
    }
  )
  if (is.null(hv_list_obj_for_plot)) return()
  plot_filename <- file.path(output_dir, paste0("hypervolumes_original_vars", plot_suffix, "_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".png"))
  palette_for_hv_plot <- global_palette[names(global_palette) %in% names(hvs_to_plot_actual)]
  if (length(palette_for_hv_plot) != length(hvs_to_plot_actual)) {
    if (verbose) warning("        Original HV Plot: Palette issue for hypervolume plot. Using default colors.")
    num_hvs_to_plot <- length(hvs_to_plot_actual)
    palette_for_hv_plot <- stats::setNames(
      RColorBrewer::brewer.pal(max(3, num_hvs_to_plot), "Set2")[1:num_hvs_to_plot],
      names(hvs_to_plot_actual)
    )
    if(length(palette_for_hv_plot) != num_hvs_to_plot) palette_for_hv_plot <- NULL
  }
  axis_labels_plot <- vars_used_for_hv[1:dims_to_show_in_plot]
  if(!is.null(niche_data_object_metadata$labels$continuous_vars)){
    labeled_axes <- niche_data_object_metadata$labels$continuous_vars[axis_labels_plot]
    axis_labels_plot <- ifelse(is.na(labeled_axes), axis_labels_plot, labeled_axes)
  }
  grDevices::png(plot_filename, width = 1200, height = 1000, res = 150, pointsize = 20)
  plot_success_hv_orig <- FALSE
  tryCatch({
    old_par <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mar = c(5, 5, 5, 2) + 0.1)
    hypervolume::plot.HypervolumeList(
      hv_list_obj_for_plot,
      main = paste("Hypervolumes on Original Variables (", dims_to_show_in_plot, "D shown)", plot_suffix, "\nTarget:", target_id),
      colors = palette_for_hv_plot,
      show.legend = TRUE,
      show.random = TRUE, random.points.cex = 0.5,
      show.centroid = TRUE, centroid.pch = 21, centroid.cex = 1.2,
      contour.lwd = 1.5, contour.type = "kde",
      kde.plot.npoints = 50,
      plot.dimensionality = dims_to_show_in_plot 
    )
    plot_success_hv_orig <- TRUE
  }, error = function(e) {
    if (verbose) warning("        Original HV Plot: Error during hypervolume::plot.HypervolumeList: ", e$message)
    graphics::plot(1, type="n", xlab="", ylab="", main="Error in Original Hypervolume Plot")
    graphics::text(1,1, "Could not generate original hypervolume plot.")
  }, finally = {
    grDevices::dev.off()
  })
  
  if (plot_success_hv_orig && verbose) {
    message("      Hypervolume plot on original variables", plot_suffix, " saved to: ", plot_filename)
  } else if (verbose) {
    message("      Hypervolume plot on original variables", plot_suffix, " failed to save or was not generated.")
  }
}

#' Plot PCA Biplot (Internal)
#'
#' Generates a PCA biplot using factoextra.
#'
#' @param pca_results The output from `stats::prcomp`.
#' @param data_for_pca_plot Tibble. The data used to fit the PCA, must contain the `species_id_col`.
#' @param target_id Character. ID of the target species.
#' @param species_id_col Character. Name of the species ID column.
#' @param niche_data_metadata List. The `metadata` element from `niche_data_object`.
#' @param global_palette Named character vector. Palette mapping species names to colors.
#' @param output_dir Character. Directory to save the plot.
#' @param plot_suffix Character. Suffix for the plot filename.
#' @param verbose Logical.
#' @noRd
.plot_pca_biplot <- function(pca_results, data_for_pca_plot, target_id, species_id_col,
                             niche_data_metadata, global_palette,
                             output_dir, plot_suffix = "", verbose = TRUE) {
  
  if (is.null(pca_results)) {
    if (verbose) message("      PCA results are NULL, skipping PCA plot", plot_suffix)
    return()
  }
  if (is.null(pca_results$x) || nrow(pca_results$x) == 0 || ncol(pca_results$x) < 2) {
    if (verbose) message("      PCA results ($x) are insufficient (no data or <2 PCs), skipping PCA plot", plot_suffix)
    return()
  }
  if (is.null(pca_results$rotation) || nrow(pca_results$rotation) == 0 || ncol(pca_results$rotation) < 1) {
    if (verbose) message("      PCA results ($rotation) are insufficient (no variables or PCs), skipping PCA plot", plot_suffix)
    return()
  }
  if (verbose) message("      Generating PCA biplot", plot_suffix, " for target: ", target_id)
  species_in_pca_data <- unique(as.character(data_for_pca_plot[[species_id_col]]))
  if (length(species_in_pca_data) == 0) {
    if(verbose) message("        PCA plot", plot_suffix, ": No species found in data_for_pca_plot. Skipping.")
    return()
  }
  ordered_species_levels_pca <- intersect(c(target_id, sort(setdiff(species_in_pca_data, target_id))), names(global_palette))
  
  if (length(ordered_species_levels_pca) == 0) {
    if(verbose) message("        PCA plot", plot_suffix, ": No species levels for PCA plot after intersecting with palette. Skipping.")
    return()
  }
  palette_for_pca_plot <- global_palette[names(global_palette) %in% ordered_species_levels_pca]
  
  if (length(palette_for_pca_plot) == 0) { 
    if(verbose) message("        PCA plot", plot_suffix, ": Palette for PCA plot is empty. Skipping.")
    return()
  }
  data_for_pca_plot_factored <- data_for_pca_plot
  data_for_pca_plot_factored[[species_id_col]] <- factor(data_for_pca_plot_factored[[species_id_col]], levels = ordered_species_levels_pca)
  data_for_pca_plot_factored <- data_for_pca_plot_factored %>%
    dplyr::filter(!is.na(!!dplyr::sym(species_id_col))) 
  if(nrow(data_for_pca_plot_factored) == 0) {
    if(verbose) message("        PCA plot", plot_suffix, ": No data remaining after factoring species and removing NAs. Skipping.")
    return()
  }
  species_in_factored_data <- unique(as.character(data_for_pca_plot_factored[[species_id_col]]))
  ordered_species_levels_pca_final <- intersect(ordered_species_levels_pca, species_in_factored_data)
  if (length(ordered_species_levels_pca_final) == 0) {
    if(verbose) message("        PCA plot", plot_suffix, ": No species levels remaining after final filtering for PCA. Skipping.")
    return()
  }
  data_for_pca_plot_factored[[species_id_col]] <- factor(data_for_pca_plot_factored[[species_id_col]], levels = ordered_species_levels_pca_final)
  palette_for_pca_plot_final <- palette_for_pca_plot[names(palette_for_pca_plot) %in% ordered_species_levels_pca_final]
  species_counts_for_ellipse <- table(data_for_pca_plot_factored[[species_id_col]])
  can_draw_ellipses <- all(species_counts_for_ellipse[species_counts_for_ellipse > 0] >= 3) && length(species_counts_for_ellipse[species_counts_for_ellipse > 0]) >= 1
  if (!can_draw_ellipses && verbose) {
    message("        PCA plot", plot_suffix, ": Ellipses will not be drawn (some species < 3 obs or only one species group with enough obs).")
  }
  pca_results_labeled_vars <- pca_results
  current_var_names <- rownames(pca_results$rotation)
  labels_for_current_vars <- niche_data_metadata$labels$continuous_vars[current_var_names]
  final_labels <- ifelse(is.na(labels_for_current_vars) | labels_for_current_vars == "", current_var_names, labels_for_current_vars)
  names(final_labels) <- current_var_names
  rownames(pca_results_labeled_vars$rotation) <- final_labels[rownames(pca_results$rotation)]
  p_pca <- tryCatch({
    factoextra::fviz_pca_biplot(pca_results_labeled_vars,
                                geom.ind = "point",
                                col.ind = data_for_pca_plot_factored[[species_id_col]],
                                palette = palette_for_pca_plot_final, 
                                pointshape = 16,
                                addEllipses = can_draw_ellipses,
                                ellipse.type = "norm", ellipse.level = 0.90, ellipse.alpha = 0.1,
                                pointsize = 2.5,
                                geom.var = c("arrow", "text"), col.var = "grey20",
                                repel = TRUE,
                                title = paste("PCA Biplot", plot_suffix, "-", target_id),
                                legend.title = "Species") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "right")
  }, error = function(e) {
    if(verbose) warning("      Error generating PCA biplot", plot_suffix, ": ", e$message)
    NULL
  })
  if (!is.null(p_pca)) {
    tryCatch({
      ggplot2::ggsave(filename = file.path(output_dir, paste0("pca_biplot", plot_suffix, "_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".png")),
                      plot = p_pca, width = 12, height = 8, dpi = 300)
      if (verbose) message("      PCA biplot", plot_suffix, " saved.")
    }, error = function(e_save) {
      if(verbose) warning("      Error SAVING PCA biplot", plot_suffix, " (viewport error likely here if p_pca was bad): ", e_save$message)
      dummy_plot <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle(paste("PCA Biplot Failed for", target_id, plot_suffix))
      ggplot2::ggsave(filename = file.path(output_dir, paste0("pca_biplot_FAILED_", plot_suffix, "_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".png")),
                      plot = dummy_plot, width = 12, height = 8, dpi = 72)
    })
  }
}

#' Plot LDA Results (Internal)
#'
#' Generates scatter plot of LDA results (LD1 vs LD2 or LD1 density).
#'
#' @param lda_results_module List. The output from `.run_lda_module()`.
#' @param data_for_lda_plot Tibble. Data used for LDA, containing species ID and LD scores.
#' @param target_id Character. ID of the target species.
#' @param species_id_col Character. Name of the species ID column.
#' @param global_palette Named character vector. Palette mapping species names to colors.
#' @param output_dir Character. Directory to save the plot.
#' @param plot_suffix Character. Suffix for the plot filename.
#' @param verbose Logical.
#' @noRd
.plot_lda_results <- function(lda_results_module, data_for_lda_plot, target_id, species_id_col,
                              global_palette, output_dir, plot_suffix = "", verbose = TRUE) {
  
  if (is.null(lda_results_module) || is.null(lda_results_module$ld_scores)) {
    if (verbose) message("      LDA results or LD scores are NULL, skipping LDA plot", plot_suffix)
    return()
  }
  if (verbose) message("      Generating LDA plot", plot_suffix, " for target: ", target_id)
  
  lda_scores_df <- lda_results_module$ld_scores
  data_from_lda_model_run <- lda_results_module$data_used_for_lda
  
  if (is.null(data_from_lda_model_run) || nrow(data_from_lda_model_run) == 0) {
    if (verbose) message("      LDA plot", plot_suffix, ": Data used for LDA model is missing or empty. Skipping plot.")
    return()
  }
  if (nrow(data_from_lda_model_run) != nrow(lda_scores_df)) {
    if (verbose) warning("      LDA plot", plot_suffix, ": Mismatch in rows between LDA data (", nrow(data_from_lda_model_run) ,") and LD scores (", nrow(lda_scores_df),"). Skipping plot.")
    return()
  }
  full_lda_data_with_scores <- dplyr::bind_cols(data_from_lda_model_run, lda_scores_df)
  species_to_display_in_current_plot <- unique(as.character(data_for_lda_plot[[species_id_col]]))
  plot_data_final_filtered <- full_lda_data_with_scores %>%
    dplyr::filter(!!dplyr::sym(species_id_col) %in% species_to_display_in_current_plot)
  if (nrow(plot_data_final_filtered) == 0) {
    if (verbose) message("      LDA plot", plot_suffix, ": No data remaining after filtering for species: ", paste(species_to_display_in_current_plot, collapse=", "), ". Skipping plot.")
    return()
  }
  num_unique_species_in_plot <- length(unique(plot_data_final_filtered[[species_id_col]]))
  if (num_unique_species_in_plot == 0) { 
    if (verbose) message("      LDA plot", plot_suffix, ": No species to plot. Skipping.")
    return()
  }
  species_in_plot_data <- unique(as.character(plot_data_final_filtered[[species_id_col]]))
  ordered_species_levels_lda <- intersect(c(target_id, sort(setdiff(species_in_plot_data, target_id))), names(global_palette))
  plot_data_final_filtered[[species_id_col]] <- factor(plot_data_final_filtered[[species_id_col]], levels = ordered_species_levels_lda)
  palette_for_lda_plot <- global_palette[names(global_palette) %in% ordered_species_levels_lda]
  p_lda <- NULL
  if ("LD2" %in% names(lda_scores_df)) {
    if (num_unique_species_in_plot < 2 && verbose) {
      message("      LDA plot", plot_suffix, ": Need at least 2 species for a 2D LDA plot. Skipping.")
    } else {
      p_lda <- ggplot2::ggplot(plot_data_final_filtered, ggplot2::aes(x = LD1, y = LD2, color = !!dplyr::sym(species_id_col))) +
        ggplot2::geom_point(alpha = 0.7, size = 2.5) +
        ggplot2::stat_ellipse(ggplot2::aes(fill = !!dplyr::sym(species_id_col)), geom = "polygon", alpha = 0.1, level = 0.90) +
        ggplot2::scale_color_manual(values = palette_for_lda_plot) +
        ggplot2::scale_fill_manual(values = palette_for_lda_plot) +
        ggplot2::labs(title = paste("LDA Plot (LD1 vs LD2)", plot_suffix, "-", target_id),
                      subtitle = paste("CV Accuracy:", round(lda_results_module$cv_accuracy, 3)),
                      color = "Species", fill = "Species") +
        ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "right")
      current_lda_plot_width <- if (plot_suffix == "_all_comparison") 14 else 10
      
      ggplot2::ggsave(filename = file.path(output_dir, paste0("lda_plot_2D", plot_suffix, "_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".png")),
                      plot = p_lda, width = current_lda_plot_width, height = 7, dpi = 300)
    }
  } else if ("LD1" %in% names(lda_scores_df)) {
    if (num_unique_species_in_plot < 1 && verbose) {
      message("      LDA plot", plot_suffix, ": Need at least 1 species for a 1D LDA plot. Skipping.")
    } else {
      p_lda <- ggplot2::ggplot(plot_data_final_filtered, ggplot2::aes(x = LD1, y = !!dplyr::sym(species_id_col), color = !!dplyr::sym(species_id_col))) +
        ggridges::geom_density_ridges(ggplot2::aes(fill = !!dplyr::sym(species_id_col)), alpha = 0.7, scale = 0.9) +
        ggplot2::scale_color_manual(values = palette_for_lda_plot) +
        ggplot2::scale_fill_manual(values = palette_for_lda_plot) +
        ggplot2::labs(title = paste("LDA Plot (LD1 Density)", plot_suffix, "-", target_id),
                      subtitle = paste("CV Accuracy:", round(lda_results_module$cv_accuracy, 3)),
                      x = "LD1 Score", y = "Species") +
        ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "none")
      ggplot2::ggsave(filename = file.path(output_dir, paste0("lda_plot_1D", plot_suffix, "_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".png")),
                      plot = p_lda, width = 8, height = max(6, num_unique_species_in_plot*0.8), dpi = 300, limitsize = FALSE)
    }
  }
  if (!is.null(p_lda) && verbose) message("      LDA plot", plot_suffix, " saved.")
}

#' Plot Modularity Network (Internal)
#'
#' Visualizes the niche overlap network with detected communities.
#'
#' @param modularity_results List. Output from `.run_modularity_module()`.
#' @param target_id Character. ID of the target species.
#' @param overlap_threshold Numeric. The overlap threshold used (for title).
#' @param overlap_metric_name Character. Name of the overlap metric used (for title).
#' @param output_dir Character. Directory to save the plot.
#' @param plot_suffix Character. Suffix for the plot filename.
#' @param file_path_override Character. Optional. Full path to save the plot, overriding default naming.
#' @param verbose Logical.
#' @noRd
.plot_modularity_network <- function(modularity_results, target_id,
                                     overlap_threshold, overlap_metric_name,
                                     global_palette,
                                     output_dir, 
                                     plot_suffix = "", file_path_override = NULL, verbose = TRUE) {
  if (is.null(modularity_results) || is.null(modularity_results$network)) {
    if (verbose) message("      Modularity network plot", plot_suffix, ": Modularity results or network object are NULL. Skipping plot.")
    return()
  }
  if (verbose) message("      Generating modularity network plot", plot_suffix, " for target: ", target_id)
  net <- modularity_results$network
  groups <- modularity_results$groups
  mod_score <- modularity_results$modularity_score
  if (igraph::gsize(net) == 0) {
    if(verbose) message("      Modularity network plot", plot_suffix, ": Network is empty (no nodes). Skipping plot.")
    return()
  }
  igraph::V(net)$size <- 6 
  igraph::V(net)$frame.color <- "grey30" 
  igraph::V(net)$label.cex <- 0.7 
  igraph::V(net)$label.color <- "black"
  igraph::V(net)$shape <- "circle" 
  if (!is.null(groups) && length(groups) > 0) {
    if(is.null(names(groups)) || !all(names(groups) %in% igraph::V(net)$name)) {
      if(verbose) warning("      Modularity network plot", plot_suffix, ": Group names do not match network vertex names. Using default colors.")
      groups <- NULL 
    } else {
      unique_group_ids <- sort(unique(as.character(groups)))
      num_groups_plot <- length(unique_group_ids)
      if (num_groups_plot > 0) {
        palette_groups <- character(0)
        if (num_groups_plot <= 8) {
          palette_groups <- RColorBrewer::brewer.pal(max(3, num_groups_plot), "Set2")
        } else {
          palette_groups <- randomcoloR::distinctColorPalette(num_groups_plot)
        }
        group_color_map <- stats::setNames(palette_groups[1:num_groups_plot], as.character(unique_group_ids))
        igraph::V(net)$color <- group_color_map[as.character(groups[igraph::V(net)$name])]
        igraph::V(net)$color[is.na(igraph::V(net)$color)] <- "grey" 
      } else { 
        if(verbose) warning("      Modularity network plot", plot_suffix, ": No unique group IDs found despite groups object being present. Using default node color.")
        igraph::V(net)$color <- "skyblue"
      }    }
  } else {
    igraph::V(net)$color <- "skyblue" 
  }
  if (target_id %in% igraph::V(net)$name) {
    igraph::V(net)$shape[igraph::V(net)$name == target_id] <- "square" 
    igraph::V(net)$size[igraph::V(net)$name == target_id] <- 9 
    igraph::V(net)$frame.color[igraph::V(net)$name == target_id] <- "black"
    igraph::V(net)$frame.width <- ifelse(igraph::V(net)$name == target_id, 2, 1)
  } else { 
    igraph::V(net)$frame.width <- 1
  }
  if (igraph::gsize(net) > 0 && "weight" %in% igraph::edge_attr_names(net)) {
    valid_weights <- igraph::E(net)$weight[!is.na(igraph::E(net)$weight)]
    if (length(valid_weights) > 0) {
      min_w <- min(valid_weights, na.rm = TRUE)
      max_w <- max(valid_weights, na.rm = TRUE)
      if (max_w > min_w) {
        igraph::E(net)$width <- 1 + 4 * (igraph::E(net)$weight - min_w) / (max_w - min_w)
      } else {
        igraph::E(net)$width <- 2.5 
      }
      igraph::E(net)$width[is.na(igraph::E(net)$width)] <- 1
    } else { igraph::E(net)$width <- 1.5 }
  } else {
    igraph::E(net)$width <- 1.5
  }
  igraph::E(net)$color <- grDevices::adjustcolor("gray50", alpha.f = 0.6)
  plot_layout <- igraph::layout_with_fr(net)
  actual_plot_filename <- if (!is.null(file_path_override)) {
    file_path_override
  } else {
    file.path(output_dir, paste0("modularity_network", gsub("[^A-Za-z0-9_]", "", plot_suffix), "_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".png"))
  }
  plot_title <- paste("Niche Modularity Network (", gsub("_", " ", plot_suffix), ")\nRef. Target: ", target_id,
                      " | Metric: ", gsub("mean_", "", overlap_metric_name), " > ", round(overlap_threshold, 2),
                      " | Modularity: ", round(mod_score, 3))
  grDevices::png(actual_plot_filename, width = 1200, height = 1000, res = 120)
  plot_success <- FALSE
  tryCatch({
    old_par <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(old_par), add = TRUE)     
    graphics::par(mar = c(1, 1, 3, 1))
    igraph::plot.igraph(net, layout = plot_layout,
                        vertex.label = igraph::V(net)$name,
                        main = plot_title,
                        cex.main = 0.9)
    leg_labels_edge_width <- NULL
    widths_for_legend <- NULL
    if (igraph::gsize(net) > 0 && "weight" %in% igraph::edge_attr_names(net) && any(!is.na(igraph::E(net)$weight))) {
      valid_overlaps_legend <- igraph::E(net)$weight[!is.na(igraph::E(net)$weight)]
      if(length(valid_overlaps_legend) > 0) {
        min_w_val <- min(valid_overlaps_legend, na.rm = TRUE)
        med_w_val <- stats::median(valid_overlaps_legend, na.rm = TRUE)
        max_w_val <- max(valid_overlaps_legend, na.rm = TRUE)
        all_edge_widths <- igraph::E(net)$width[!is.na(igraph::E(net)$width)]
        if(length(all_edge_widths) > 0) {
          min_lwd <- min(all_edge_widths, na.rm = TRUE)
          med_lwd <- stats::median(all_edge_widths, na.rm = TRUE)
          max_lwd <- max(all_edge_widths, na.rm = TRUE)
          candidate_widths <- unique(c(min_lwd, med_lwd, max_lwd))
          candidate_labels <- round(unique(c(min_w_val, med_w_val, max_w_val)), 2)
          widths_for_legend <- candidate_widths[is.finite(candidate_widths)]
          leg_labels_edge_width <- candidate_labels[is.finite(candidate_labels)]
          if(length(widths_for_legend) > 0 && length(leg_labels_edge_width) > 0 && length(widths_for_legend) != length(leg_labels_edge_width)) {
            widths_for_legend <- widths_for_legend[1:min(length(widths_for_legend), length(leg_labels_edge_width))]
            leg_labels_edge_width <- leg_labels_edge_width[1:min(length(widths_for_legend), length(leg_labels_edge_width))]
          }
        } else {
          widths_for_legend <- c(1.5)
          leg_labels_edge_width <- round(med_w_val, 2)
          if(is.na(leg_labels_edge_width) || !is.finite(leg_labels_edge_width)) leg_labels_edge_width <- NULL 
        }
      }
    }
    if (!is.null(groups) && length(unique(groups)) > 1 && exists("group_color_map") && exists("unique_group_ids")) {
      graphics::legend("bottomleft", legend = paste("Module", unique_group_ids), fill = group_color_map[as.character(unique_group_ids)],
                       bty = "n", title = "Modules", cex = 0.8)
    }
    if (!is.null(leg_labels_edge_width) && length(leg_labels_edge_width) > 0 &&
        !is.null(widths_for_legend) && length(widths_for_legend) > 0 &&
        length(leg_labels_edge_width) == length(widths_for_legend)) {
      graphics::legend("bottomright", title = paste("Overlap (", overlap_metric_name, ")", sep=""), legend = leg_labels_edge_width, lwd = widths_for_legend, col = grDevices::adjustcolor("gray50", alpha.f = 0.6), bty = "n", cex = 0.7, seg.len=1.5)
    }
    plot_success <- TRUE 
  }, error = function(e) {
    if(verbose) {
      warning("      ERROR during igraph plot/legend for modularity network '", plot_suffix, "': ", e$message)
      warning("      Full error details for modularity plot: ", paste(capture.output(print(e)), collapse="\n"))
    } 
    graphics::plot(1, type="n", main="Modularity Plot Failed", xlab="", ylab="")
    graphics::text(1,1, "Error generating network plot.")
  }, finally = {
    grDevices::dev.off()
    if (!plot_success && verbose) {
      message("      Modularity network plot '", plot_suffix, "' for target '", target_id, "' FAILED to generate due to an error during plotting operations.")
    }
  })
  if (plot_success && verbose) {
    message("      Modularity network plot", plot_suffix, " saved to: ", actual_plot_filename)
  }}

#' Generate Summary Tables for Top N Similar Species (Internal)
#'
#' Creates CSV summary tables for different analyses, focusing on the top N
#' most similar species to the target.
#'
#' @param target_id Character. ID of the target species.
#' @param all_results_for_target List. Aggregated results for the current target.
#'        This list is expected to contain elements like `overlap_indices`,
#'        `permanova_pairwise`, `equivalency_test`, `lda_multigroup`,
#'        `hypervolume_analysis$overlap_stats`.
#' @param settings An object of class 'niche_analysis_settings'.
#' @param output_dir Character. Path to the target-specific output directory.
#' @param niche_data_object_metadata List. Metadata from `niche_data_object` for labels.
#' @param verbose Logical.
#' @noRd

.generate_top_n_summary_tables_per_analysis <- function(target_id,
                                                        all_results_for_target,
                                                        settings,
                                                        output_dir,
                                                        niche_data_object_metadata,
                                                        verbose = TRUE) {
  if (verbose) message("    Generating summary tables for top N similar species...")
  
  if (is.null(all_results_for_target$overlap_indices) || nrow(all_results_for_target$overlap_indices) == 0) {
    if (verbose) message("      No overlap indices found for target '", target_id, "'. Skipping summary tables.")
    return()
  }
  n_closest <- settings$plot_n_closest %||% 10
  overlap_metric_for_ranking <- settings$plot_n_closest_overlap_metric %||% "composite_overlap"
  if (!overlap_metric_for_ranking %in% names(all_results_for_target$overlap_indices)) {
    if (verbose) warning("      Overlap metric '", overlap_metric_for_ranking, "' not found in overlap_indices. Using 'composite_overlap' or first available.")
    overlap_metric_for_ranking <- intersect(c("composite_overlap", "mean_schoener_D", "mean_pianka_O", "mean_czek_psi", "mean_hurlbert_L"), names(all_results_for_target$overlap_indices))[1]
    if(is.na(overlap_metric_for_ranking)){
      if(verbose) message("      No suitable overlap metric found for ranking. Skipping summary tables.")
      return()
    }
  }
  top_n_species_df <- all_results_for_target$overlap_indices %>%
    dplyr::filter(target_species == target_id) %>%
    dplyr::arrange(dplyr::desc(!!dplyr::sym(overlap_metric_for_ranking))) %>%
    dplyr::slice_head(n = n_closest)
  top_n_comparison_ids <- top_n_species_df$comparison_species
  if (length(top_n_comparison_ids) == 0) {
    if (verbose) message("      No similar species identified for target '", target_id, "' based on overlap ranking. Skipping summary tables.")
    return()
  }
  # --- 1. Overlap Indices Summary Table ---
  summary_overlap <- top_n_species_df %>%
    dplyr::select(target_species, comparison_species, dplyr::starts_with("mean_"), dplyr::starts_with("composite_"))
  utils::write.csv(summary_overlap, file.path(output_dir, paste0("summary_top", n_closest, "_overlap_indices_", target_id, ".csv")), row.names = FALSE)
  # --- 2. PERMANOVA Pairwise Summary Table ---
  if (!is.null(all_results_for_target$permanova_pairwise)) {
    summary_permanova <- all_results_for_target$permanova_pairwise %>%
      dplyr::filter(specie1 == target_id & specie2 %in% top_n_comparison_ids) %>%
      dplyr::left_join(top_n_species_df %>% dplyr::select(comparison_species, !!dplyr::sym(overlap_metric_for_ranking)), by = c("specie2" = "comparison_species")) %>%
      dplyr::arrange(dplyr::desc(!!dplyr::sym(overlap_metric_for_ranking)))
    utils::write.csv(summary_permanova, file.path(output_dir, paste0("summary_top", n_closest, "_permanova_", target_id, ".csv")), row.names = FALSE)
  }
  # --- 3. Equivalency Test Summary Table ---
  if (!is.null(all_results_for_target$equivalency_test)) {
    summary_equivalency <- all_results_for_target$equivalency_test %>%
      dplyr::filter(target_species == target_id & comparison_species %in% top_n_comparison_ids) %>%
      dplyr::left_join(top_n_species_df %>% dplyr::select(comparison_species, !!dplyr::sym(overlap_metric_for_ranking)), by = "comparison_species") %>%
      dplyr::arrange(dplyr::desc(!!dplyr::sym(overlap_metric_for_ranking)))
    utils::write.csv(summary_equivalency, file.path(output_dir, paste0("summary_top", n_closest, "_equivalency_", target_id, ".csv")), row.names = FALSE)
  }
  # --- 4. LDA Summary (e.g., CV accuracy, top discriminating vars if applicable) ---
  if (!is.null(all_results_for_target$lda_multigroup)) {
    lda_summary_info <- list(
      target_species = target_id,
      overall_cv_accuracy = all_results_for_target$lda_multigroup$cv_accuracy,
      top_discriminating_vars_ld1 = if(!is.null(all_results_for_target$lda_multigroup$variable_importance_ld1)) {
        all_results_for_target$lda_multigroup$variable_importance_ld1 %>%
          dplyr::slice_head(n = 5) %>% # Top 5
          dplyr::mutate(Variable_Label = niche_data_object_metadata$labels$continuous_vars[Variable] %||% Variable) %>%
          dplyr::pull(Variable_Label) %>% paste(collapse = "; ")
      } else { "N/A" }
    )
    utils::write.csv(as.data.frame(lda_summary_info), file.path(output_dir, paste0("summary_lda_overall_for_top", n_closest, "_context_", target_id, ".csv")), row.names = FALSE)
  }
  # --- 5. Hypervolume Overlap Summary Table ---
  if (!is.null(all_results_for_target$hypervolume_analysis$overlap_stats)) {
    summary_hypervolume <- all_results_for_target$hypervolume_analysis$overlap_stats %>%
      dplyr::filter(target_species == target_id & comparison_species %in% top_n_comparison_ids) %>%
      dplyr::left_join(top_n_species_df %>% dplyr::select(comparison_species, !!dplyr::sym(overlap_metric_for_ranking)), by = "comparison_species") %>%
      dplyr::arrange(dplyr::desc(!!dplyr::sym(overlap_metric_for_ranking)))
    utils::write.csv(summary_hypervolume, file.path(output_dir, paste0("summary_top", n_closest, "_hypervolume_overlap_", target_id, ".csv")), row.names = FALSE)
  }
  # --- 6. Categorical Association Summary (if applicable) ---
  if (!is.null(all_results_for_target$categorical_association$env_categorical_test_results)) {
    sig_cat_associations_target <- all_results_for_target$categorical_association$env_categorical_test_results %>%
      dplyr::filter(P_Value < 0.05) %>%
      dplyr::mutate(Variable_Label = niche_data_object_metadata$labels$categorical_vars[[unique(Variable)]][Variable] %||% Variable)
    if(nrow(sig_cat_associations_target) > 0){
      utils::write.csv(sig_cat_associations_target, file.path(output_dir, paste0("summary_significant_categorical_associations_for_", target_id, ".csv")), row.names = FALSE)
    }
  }
  if (verbose) message("    Summary tables for top N species generated for target '", target_id, "'.")
}

#' Plot Hypervolumes on PCA (Internal)
#'
#' Performs PCA on selected variables for a subset of species,
#' then constructs and plots hypervolumes on the first few PCA dimensions.
#'
#' @param data_for_module Tibble. Data for the current target and all its comparison species,
#'        containing the species ID column and all *original selected environmental variables*.
#' @param target_id Character. ID of the target species.
#' @param species_for_hv_plot Character vector. Species to include in this plot (target + N closest).
#' @param species_id_col Character. Name of the species ID column.
#' @param vars_for_hv_pca Character vector. Names of *original environmental variables*
#'        to be used for the PCA (typically top N numeric vars from selection).
#' @param settings An object of class 'niche_analysis_settings'.
#' @param global_palette Named character vector. Palette mapping species names to colors.
#' @param output_dir Character. Directory to save the plot.
#' @param plot_suffix Character. Suffix for the plot filename.
#' @param verbose Logical.
#' @noRd
.plot_hypervolumes_on_pca <- function(data_for_module, target_id, species_for_hv_plot,
                                      species_id_col, vars_for_hv_pca,
                                      settings, global_palette,
                                      output_dir, plot_suffix = "", verbose = TRUE) {
  if (verbose) message("      Generating hypervolume plot on PCA dimensions", plot_suffix, " for target: ", target_id)
  if (is.null(data_for_module) || nrow(data_for_module) == 0) {
    if(verbose) message("        HV PCA Plot: Input data_for_module is NULL or empty. Skipping.")
    return()
  }
  if (length(species_for_hv_plot) < 2) {
    if(verbose) message("        HV PCA Plot: Less than 2 species to plot. Skipping.")
    return()
  }
  if (length(vars_for_hv_pca) < 2) {
    if(verbose) message("        HV PCA Plot: Less than 2 variables for PCA. Skipping.")
    return()
  }
  # 1. Filter data for selected species and variables for PCA
  data_pca_hv <- data_for_module %>%
    dplyr::filter(!!dplyr::sym(species_id_col) %in% species_for_hv_plot) %>%
    dplyr::select(dplyr::all_of(c(species_id_col, vars_for_hv_pca))) %>%
    stats::na.omit()

  if (nrow(data_pca_hv) == 0) {
    if(verbose) message("        HV PCA Plot: No data remaining after filtering for species/vars and NA removal. Skipping.")
    return()
  }
  if (nrow(data_pca_hv) <= length(vars_for_hv_pca) || dplyr::n_distinct(data_pca_hv[[species_id_col]]) < 2) {
    if(verbose) message("        HV PCA Plot: Insufficient data after filtering for PCA (obs <= vars, or < 2 groups). Skipping.")
    return()
  }
  numeric_data_for_pca <- data_pca_hv[, vars_for_hv_pca, drop = FALSE]
  species_ids_for_scores <- data_pca_hv[[species_id_col]] # Extract species IDs from the cleaned data
  sd_check <- apply(numeric_data_for_pca, 2, stats::sd, na.rm = TRUE)
  if (any(sd_check == 0, na.rm = TRUE) && ncol(numeric_data_for_pca) > 0) {
    if(verbose) warning("        HV PCA Plot: Zero variance in one or more variables for PCA. Skipping. Vars: ", paste(names(sd_check[sd_check == 0]), collapse=", "))
    return()
  }
  # 2. Perform PCA
  pca_res_hv <- tryCatch({
    stats::prcomp(numeric_data_for_pca, scale. = TRUE, center = TRUE)
  }, error = function(e) {
    if(verbose) warning("        HV PCA Plot: PCA failed: ", e$message)
    NULL
  })
  if (is.null(pca_res_hv)) return()
  if(verbose){
    message("        HV PCA Plot: PCA successful. Summary of PCA results:")
    print(summary(pca_res_hv))
    message("        HV PCA Plot: Dimensions of PCA scores (pca_res_hv$x): ", paste(dim(pca_res_hv$x), collapse="x"))
  }
  # 3. Select top N PCA dimensions for hypervolume construction
  n_dims_for_hv_on_pca <- min(settings$hypervolume_n_vars_plot %||% 2, ncol(pca_res_hv$x))
  if(verbose) message("        HV PCA Plot: Number of PCA dimensions to use for HV: ", n_dims_for_hv_on_pca)
  if (n_dims_for_hv_on_pca < 2) {
    if(verbose) message("        HV PCA Plot: Less than 2 PCA dimensions available/selected for hypervolume. Skipping.")
    return()
  }
  pca_scores_for_hv_build <- as.data.frame(pca_res_hv$x[, 1:n_dims_for_hv_on_pca, drop = FALSE])
  pca_scores_for_hv_build[[species_id_col]] <- species_ids_for_scores # Use the correctly ordered species IDs
  # 4. Calculate hypervolumes on these PCA scores
  hv_list_on_pca <- list()
  if(verbose) message("        HV PCA Plot: Building hypervolumes on ", n_dims_for_hv_on_pca, " PCA dimensions for ", length(species_for_hv_plot), " species.")
  if(verbose) message("        HV PCA Plot: Species for HV plot: ", paste(species_for_hv_plot, collapse=", "))
  for (sp_name in species_for_hv_plot) {
    sp_pca_scores <- pca_scores_for_hv_build %>%
      dplyr::filter(!!dplyr::sym(species_id_col) == sp_name) %>%
      dplyr::select(-dplyr::all_of(species_id_col))
    if(verbose) message("          HV PCA Plot: Processing species '", sp_name, "'. Number of PCA score rows: ", nrow(sp_pca_scores))
    if (nrow(sp_pca_scores) > n_dims_for_hv_on_pca) {
      if(verbose) message("            HV PCA Plot: Attempting hypervolume_gaussian for '", sp_name, "'.")
      hv_list_on_pca[[sp_name]] <- tryCatch({
        suppressMessages({
          hypervolume::hypervolume_gaussian(sp_pca_scores, 
                                            name = sp_name,
                                            samples.per.point = settings$hypervolume_samples_per_point %||% 50, 
                                            verbose = FALSE)        
        })
      }, error = function(e_hv) {
        if(verbose) warning("        HV PCA Plot: Error calculating hypervolume for ", sp_name, " on PCA scores: ", e_hv$message)
        NULL
      })
    } else {
      if(verbose) warning("        HV PCA Plot: Insufficient data points for ", sp_name, " (", nrow(sp_pca_scores), ") to build hypervolume on ", n_dims_for_hv_on_pca, " PCA dimensions.")
    }
  }
  hv_list_on_pca <- hv_list_on_pca[!sapply(hv_list_on_pca, is.null)]
  if(verbose){
    message("        HV PCA Plot: Number of valid hypervolumes constructed: ", length(hv_list_on_pca))
    if(length(hv_list_on_pca) > 0) message("        HV PCA Plot: Names of valid hypervolumes: ", paste(names(hv_list_on_pca), collapse=", "))
  }
  if (length(hv_list_on_pca) < 1) {
    if(verbose) message("        HV PCA Plot: Less than 1 valid hypervolume constructed on PCA scores. Skipping plot.")
    return()
  }
  # 5. Plot the HypervolumeList
  hv_list_obj_for_plot <- tryCatch(hypervolume::hypervolume_join(hv_list_on_pca), error = function(e) NULL)
  if(verbose && !is.null(hv_list_obj_for_plot)) message("        HV PCA Plot: hypervolume_join successful.")
  if(is.null(hv_list_obj_for_plot)){
    if(verbose) warning("        HV PCA Plot: Failed to join hypervolumes for plotting.")
    return()
  }
  plot_filename <- file.path(output_dir, paste0("hypervolumes_on_pca", plot_suffix, "_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".png"))
  palette_for_hv_plot <- global_palette[names(global_palette) %in% names(hv_list_on_pca)]
  if(length(palette_for_hv_plot) == 0 || length(palette_for_hv_plot) != length(hv_list_on_pca)) {
    if(verbose) warning("        HV PCA Plot: Palette issue. Using default colors.")
    num_hvs_to_plot <- length(hv_list_on_pca)
    palette_for_hv_plot <- stats::setNames(RColorBrewer::brewer.pal(max(3, num_hvs_to_plot), "Set2")[1:num_hvs_to_plot], names(hv_list_on_pca))
    if(length(palette_for_hv_plot) != num_hvs_to_plot) palette_for_hv_plot <- NULL
  }
  grDevices::png(plot_filename, width = 1200, height = 1000, res = 150, pointsize = 20)
  plot_success_hv_pca <- FALSE
  tryCatch({
    if(verbose) message("        HV PCA Plot: Attempting hypervolume::plot.HypervolumeList...")
    old_par <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mar = c(5, 5, 5, 2) + 0.1)
    hypervolume::plot.HypervolumeList(
      hv_list_obj_for_plot,
      main = paste("Hypervolumes on PCA (", n_dims_for_hv_on_pca, "D)", plot_suffix, "\nTarget:", target_id),
      xlab = paste0("PC1 (", round(summary(pca_res_hv)$importance[2,1]*100,1),"%)"),
      ylab = if(n_dims_for_hv_on_pca >= 2) paste0("PC2 (", round(summary(pca_res_hv)$importance[2,2]*100,1),"%)") else "",
      colors = palette_for_hv_plot,
      show.legend = TRUE,
      show.random = TRUE, random.points.cex = 0.5,
      show.centroid = TRUE, centroid.pch = 21, centroid.cex = 1.2,
      contour.lwd = 1.5, contour.type = "kde",
      kde.plot.npoints = 50
    )
    plot_success_hv_pca <- TRUE
  }, error = function(e) {
    if(verbose) warning("        HV PCA Plot: Error during hypervolume::plot.HypervolumeList: ", e$message)
    graphics::plot(1, type="n", xlab="", ylab="", main="Error in Hypervolume PCA Plot")
    graphics::text(1,1, "Could not generate hypervolume plot on PCA.")
  }, finally = {
    grDevices::dev.off()
  })
  if (plot_success_hv_pca && verbose) {
    message("        HV PCA Plot: Hypervolume plot on PCA dimensions", plot_suffix, " saved to: ", plot_filename)
  } else if (verbose) {
    message("        HV PCA Plot: Hypervolume plot on PCA dimensions", plot_suffix, " failed to save or was not generated.")
  }
}

# --- Report Generation Helper Functions ---

#' Generate Markdown Report for a Target Species (Internal)
#'
#' Creates a Markdown report summarizing the analysis results for a single target species.
#'
#' @param results_for_target List. Aggregated results for the current target species.
#' @param niche_data_object_metadata List. The `metadata` element from `niche_data_object`.
#' @param settings An object of class 'niche_analysis_settings'.
#' @param output_dir_target_specific Character. Path to the target's output directory.
#' @param verbose Logical.
#' @noRd
.generate_markdown_report_for_target <- function(results_for_target,
                                                 niche_data_object_metadata,
                                                 settings,
                                                 output_dir_target_specific,
                                                 verbose = TRUE) {
  
  if (verbose) message("    Generating Markdown report for: ", results_for_target$metadata$target_species_id)
  
  target_id <- results_for_target$metadata$target_species_id
  report_content <- c()
  add_md_section <- function(title, level = 1) paste0("\n", paste(rep("#", level), collapse = ""), " ", title, "\n")
  add_md_text <- function(...) paste0(..., "\n")
  add_md_table <- function(df, caption = "") {
    if (!is.null(df) && inherits(df, "data.frame") && nrow(df) > 0) {
      df_labeled <- df
      if ("Variable" %in% names(df_labeled) && !is.null(niche_data_object_metadata$labels$continuous_vars)) {
        df_labeled$Variable <- niche_data_object_metadata$labels$continuous_vars[df_labeled$Variable] %||% df_labeled$Variable
      }
      if ("Original_Variable" %in% names(df_labeled) && !is.null(niche_data_object_metadata$labels$continuous_vars)) {
        df_labeled$Original_Variable <- niche_data_object_metadata$labels$continuous_vars[df_labeled$Original_Variable] %||% df_labeled$Original_Variable
      }
      table_md <- knitr::kable(df_labeled, format = "pipe", caption = caption, digits = 3)
      return(paste(c(table_md, "\n"), collapse = "\n"))
    }
    return(paste0(caption, "\n\nNo data available for this table.\n"))
  }
  add_md_plot <- function(plot_filename_base, caption = "", plot_suffix = "") {
    full_plot_filename <- paste0(plot_filename_base, plot_suffix, "_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".png")
    if (file.exists(file.path(output_dir_target_specific, full_plot_filename))) {
      return(paste0("!", caption, "\n\n*Figure: ", caption, "*\n"))
    }
    return(paste0("\n*Plot '", full_plot_filename, "' not found for ", caption, ".*\n"))
  }
  # --- Report Content ---
  report_content <- c(report_content, add_md_section(paste("Niche Analysis Report for:", target_id), level = 1))
  report_content <- c(report_content, add_md_text(paste("Report generated on:", Sys.Date())))
  report_content <- c(report_content, add_md_section("Analysis Metadata", level = 2))
  meta <- results_for_target$metadata
  report_content <- c(report_content, add_md_text(paste("- **Target Species:**", meta$target_species_id)))
  report_content <- c(report_content, add_md_text(paste("- **Variable Selection Method:**", meta$analysis_settings_used$variable_selection_method)))
  report_content <- c(report_content, add_md_text(paste("- **Final Variables Used:**", paste(meta$final_variables_used, collapse = ", "))))
  report_content <- c(report_content, add_md_text(paste("- **Comparison Species Count:**", length(meta$comparison_species_used))))
  if (!is.null(meta$variable_selection_details)) {
    report_content <- c(report_content, add_md_section("Variable Selection Details", level = 2))
    if (meta$analysis_settings_used$variable_selection_method == "Boruta" && !is.null(meta$variable_selection_details$details$boruta_stats)) {
      report_content <- c(report_content, add_md_text("Boruta selection results:"))
      report_content <- c(report_content, add_md_plot("boruta_plot", caption = "Boruta Variable Importance Plot"))
      report_content <- c(report_content, add_md_table(meta$variable_selection_details$details$boruta_stats, caption = "Boruta Statistics"))
    }
    if (!is.null(meta$variable_selection_details$details$correlation_matrix)) {
      report_content <- c(report_content, add_md_plot("correlation_matrix", caption = "Correlation Matrix of Selected Numeric Variables"))
    }
    if (!is.null(meta$variable_selection_details$details$vif_values)) {
      report_content <- c(report_content, add_md_table(meta$variable_selection_details$details$vif_values, caption = "VIF Values for Final Numeric Variables"))
    }
  }
  if ("overlap_indices" %in% names(results_for_target)) {
    report_content <- c(report_content, add_md_section("Niche Overlap Indices", level = 2))
    report_content <- c(report_content, add_md_table(results_for_target$overlap_indices, caption = "Pairwise Niche Overlap Indices"))
    report_content <- c(report_content, add_md_plot("distributions_continuous", caption = "Continuous Variable Distributions (All Comparison Species)", plot_suffix = "_all_comparison"))
    report_content <- c(report_content, add_md_plot("distributions_categorical", caption = "Categorical Variable Distributions (All Comparison Species)", plot_suffix = "_all_comparison"))
    if (settings$plot_type == "n_closest_species") {
      report_content <- c(report_content, add_md_plot("distributions_continuous", caption = paste("Continuous Variable Distributions (Subset N=",settings$plot_n_closest,")", sep=""), plot_suffix = paste0("_subset", settings$plot_n_closest)))
      report_content <- c(report_content, add_md_plot("distributions_categorical", caption = paste("Categorical Variable Distributions (Subset N=",settings$plot_n_closest,")", sep=""), plot_suffix = paste0("_subset", settings$plot_n_closest)))
    }
  }
  if ("equivalency_test" %in% names(results_for_target)) {
    report_content <- c(report_content, add_md_section("Niche Equivalency Test", level = 2))
    report_content <- c(report_content, add_md_table(results_for_target$equivalency_test, caption = "Pairwise Niche Equivalency Test Results"))
  }
  if ("permanova_pairwise" %in% names(results_for_target)) {
    report_content <- c(report_content, add_md_section("Pairwise PERMANOVA", level = 2))
    report_content <- c(report_content, add_md_table(results_for_target$permanova_pairwise, caption = "Pairwise PERMANOVA Results"))
    report_content <- c(report_content, add_md_plot("permanova_pairwise_R2", caption = "Pairwise PERMANOVA R-squared Values", plot_suffix = "_all_comparison"))
    report_content <- c(report_content, add_md_plot("permanova_pairwise_pvalue", caption = "Pairwise PERMANOVA Significance (-log10 p-value)", plot_suffix = "_all_comparison"))
  }
  if ("lda_multigroup" %in% names(results_for_target) && !is.null(results_for_target$lda_multigroup)) {
    report_content <- c(report_content, add_md_section("Linear Discriminant Analysis (LDA)", level = 2))
    report_content <- c(report_content, add_md_text(paste("Cross-validated Accuracy:", round(results_for_target$lda_multigroup$cv_accuracy, 3))))
    report_content <- c(report_content, add_md_table(as.data.frame.matrix(results_for_target$lda_multigroup$cv_confusion_matrix), caption = "LDA Confusion Matrix (Cross-Validated)"))
    report_content <- c(report_content, add_md_table(results_for_target$lda_multigroup$variable_importance_ld1, caption = "LDA Variable Importance (Loadings on LD1)"))
    report_content <- c(report_content, add_md_plot("lda_plot_2D", caption = "LDA Plot (LD1 vs LD2 - All Comparison Species)", plot_suffix = "_all_comparison"))
    report_content <- c(report_content, add_md_plot("lda_plot_1D", caption = "LDA Plot (LD1 Density - All Comparison Species)", plot_suffix = "_all_comparison"))
    if (settings$plot_type == "n_closest_species") {
      report_content <- c(report_content, add_md_plot("lda_plot_2D", caption = paste("LDA Plot (LD1 vs LD2 - Subset N=",settings$plot_n_closest,")", sep=""), plot_suffix = paste0("_subset", settings$plot_n_closest)))
      report_content <- c(report_content, add_md_plot("lda_plot_1D", caption = paste("LDA Plot (LD1 Density - Subset N=",settings$plot_n_closest,")", sep=""), plot_suffix = paste0("_subset", settings$plot_n_closest)))
    }
  }
  if ("modularity_network" %in% names(results_for_target) && !is.null(results_for_target$modularity_network)) {
    report_content <- c(report_content, add_md_section("Niche Modularity Network", level = 2))
    report_content <- c(report_content, add_md_text(paste("Modularity Score:", round(results_for_target$modularity_network$modularity_score, 3))))
    report_content <- c(report_content, add_md_text(paste("Number of Detected Modules:", results_for_target$modularity_network$n_groups)))
    report_content <- c(report_content, add_md_plot("modularity_network", caption = "Niche Overlap Network and Detected Modules", plot_suffix = "_all_comparison"))
  }
  if ("hypervolume_analysis" %in% names(results_for_target) && !is.null(results_for_target$hypervolume_analysis)) {
    report_content <- c(report_content, add_md_section("Hypervolume Analysis", level = 2))
    report_content <- c(report_content, add_md_text(paste("Variables used for hypervolumes:", paste(results_for_target$hypervolume_analysis$vars_used, collapse = ", "))))
    report_content <- c(report_content, add_md_table(results_for_target$hypervolume_analysis$overlap_stats, caption = "Pairwise Hypervolume Overlap Statistics"))
    report_content <- c(report_content, add_md_plot("hypervolumes_on_pca", caption = paste("Hypervolume Comparison (Top ", settings$hypervolume_n_vars_plot, " species)", sep=""), plot_suffix = paste0("_hv_pca_top", settings$hypervolume_n_vars_plot)))
    report_content <- c(report_content, add_md_plot("dendrogram_target_specific_hv", caption = "Target-Specific Niche Dendrogram (Hypervolume Overlap)"))
  }
  if ("variable_stability" %in% names(results_for_target) && !is.null(results_for_target$variable_stability)) {
    report_content <- c(report_content, add_md_section("Variable Selection Stability (Boruta Bootstrap)", level = 2))
    report_content <- c(report_content, add_md_table(results_for_target$variable_stability, caption = "Variable Selection Frequencies from Bootstrap"))
  }
  if ("similarity_validation" %in% names(results_for_target) && !is.null(results_for_target$similarity_validation)) {
    if (!is.null(results_for_target$similarity_validation$mean_cv_accuracy)) {
      report_content <- c(report_content, add_md_section("Similarity Validation (Target vs. Statistically Similar)", level = 2))
      report_content <- c(report_content, add_md_text(paste("Mean Cross-Validated Accuracy (Random Forest):", round(results_for_target$similarity_validation$mean_cv_accuracy, 3))))
      report_content <- c(report_content, add_md_text(paste("Species validated against:", paste(results_for_target$similarity_validation$similar_species_validated_against, collapse = ", "))))
      report_content <- c(report_content, add_md_table(results_for_target$similarity_validation$cv_accuracy_results, caption = "Accuracy per CV Fold"))
    } else if (!is.null(results_for_target$similarity_validation$status_message)) {
      report_content <- c(report_content, add_md_section("Similarity Validation (Target vs. Statistically Similar)", level = 2))
      report_content <- c(report_content, add_md_text(paste0("*", results_for_target$similarity_validation$status_message, "*")))
    }
  }
  if ("categorical_association" %in% names(results_for_target) && !is.null(results_for_target$categorical_association)) {
    report_content <- c(report_content, add_md_section("Categorical Associations", level = 2))
    if (!is.null(results_for_target$categorical_association$indicator_species_rules)) {
      report_content <- c(report_content, add_md_text("Indicator Species Association Rules (Apriori):"))
      report_content <- c(report_content, add_md_table(results_for_target$categorical_association$indicator_species_rules, caption = paste("Top association rules for target", target_id)))
    }
    if (!is.null(results_for_target$categorical_association$env_categorical_test_results)) {
      report_content <- c(report_content, add_md_text("Association with Categorical Environmental Variables:"))
      report_content <- c(report_content, add_md_table(results_for_target$categorical_association$env_categorical_test_results, caption = "Chi-sq/Fisher Tests for Categorical Env. Variables"))
    }
  }
  # --- ENFA Section ---
  if (!is.null(results_for_target$enfa)) {
    report_content <- c(report_content, add_md_section("Ecological Niche Factor Analysis (ENFA)", level = 2))
    report_content <- c(report_content, add_md_text("ENFA identifies the main ecological gradients (factors) that explain the species' distribution. It separates these factors into **Marginality** (how much the species' niche deviates from the average conditions) and **Specialization** (how narrow the species' niche is)."))
    
    # Summary Table
    summary_table_path <- fs::path_join(c("ENFA_Results", paste0("ENFA_comparison_summary_", target_id, ".csv")))
    if (fs::file_exists(fs::path_join(c(output_dir_target_specific, summary_table_path)))) {
      summary_table <- utils::read.csv(fs::path_join(c(output_dir_target_specific, summary_table_path)))
      report_content <- c(report_content, add_md_section("Summary of ENFA Metrics", level = 3))
      report_content <- c(report_content, add_md_table(summary_table, caption = "Comparison of Marginality and Specialization values."))
    }
    
    # Target Species Plot
    target_plot_path <- fs::path_join(c("ENFA_Results", paste0("ENFA_", target_id, ".png")))
    if (fs::file_exists(fs::path_join(c(output_dir_target_specific, target_plot_path)))) {
      report_content <- c(report_content, add_md_section("ENFA Biplot for Target Species", level = 3))
      report_content <- c(report_content, paste0("![ENFA Biplot for ", target_id, "](", target_plot_path, ")\n\n*Figure: ENFA Biplot for ", target_id, "*\n"))
      report_content <- c(report_content, add_md_text("The plot shows the species' presence points (colored cloud) and its average position (centroid) in the ecological space defined by marginality and specialization. Blue arrows represent the influence of environmental variables."))
    }
    
    # Comparison Plot
    comparison_plot_path <- fs::path_join(c("ENFA_Results", paste0("ENFA_comparison_", target_id, ".png")))
    if (fs::file_exists(fs::path_join(c(output_dir_target_specific, comparison_plot_path)))) {
      report_content <- c(report_content, add_md_section("ENFA Comparison Biplot", level = 3))
      report_content <- c(report_content, add_md_plot(fs::path_join("ENFA_Results", paste0("ENFA_comparison_", target_id)), caption = "ENFA Comparison Biplot"))
      report_content <- c(report_content, add_md_text("This plot projects the target species and its closest neighbors onto the same ecological space. It allows for a direct visual comparison of their niche positions and breadths relative to the main environmental gradients."))
    }
    
    if (!is.null(results_for_target$enfa_closest_species) && length(results_for_target$enfa_closest_species) > 0) {
      report_content <- c(report_content, add_md_text("Individual ENFA plots for each neighbor species are also available in the `ENFA_Results` directory."))
    }
  }
  report_filepath <- file.path(output_dir_target_specific, paste0("analysis_report_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".md"))
  writeLines(unlist(report_content), report_filepath)
  if (verbose) message("    Markdown report saved to: ", report_filepath)
}

# --- Similarity Validation Module Helper Function ---

#' Run Similarity Validation using Cross-Validated Random Forest (Internal)
#'
#' Assesses how well the target species can be distinguished from a set of
#' specified "similar" species using k-fold cross-validation with Random Forest.
#'
#' @param data_for_validation Tibble. Data containing species ID and selected environmental variables.
#'        Should include the target and the similar species.
#' @param target_id Character. ID of the target species.
#' @param similar_species_ids Character vector. IDs of species considered highly similar
#'        to the target for this validation.
#' @param selected_vars Character vector. Names of selected environmental variables to use in RF.
#' @param species_id_col Character. Name of the species ID column.
#' @param similar_species_details Tibble. A detailed dataframe explaining why each
#'        species was considered similar.
#' @param categorical_vars_names Character vector. Names of categorical variables (RF handles factors).
#' @param output_dir_module Character. Path to save the CV accuracy results.
#' @param settings An object of class 'niche_analysis_settings'.
#' @param verbose Logical.
#'
#' @return A list containing:
#'   \item{cv_accuracy_results}{Tibble with accuracy per fold.}
#'   \item{mean_cv_accuracy}{Overall mean accuracy from cross-validation.}
#'   \item{n_observations_used}{Number of observations used in the validation.}
#'   \item{similar_species_validated_against}{Vector of species IDs used as "similar".}
#'         Returns NULL if validation cannot be performed.
#' @noRd
.run_similarity_validation_module <- function(data_for_validation,
                                              target_id,
                                              similar_species_ids,
                                              similar_species_details,
                                              selected_vars,
                                              species_id_col,
                                              categorical_vars_names, 
                                              output_dir_module,
                                              settings,
                                              verbose = TRUE) {
  if (verbose) message("    Running Similarity Validation module...")
  
  # Define the output path once at the top
  output_csv_path <- file.path(output_dir_module, paste0("similarity_validation_accuracy_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".csv"))
  details_csv_path <- file.path(output_dir_module, paste0("similarity_validation_details_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".csv"))

  # Save the details file if details are provided
  if (!is.null(similar_species_details) && nrow(similar_species_details) > 0) {
    # Consolidate reasons for each species
    consolidated_details <- similar_species_details %>%
      dplyr::mutate(reason = paste0(source, " (", metric, "=", round(value, 3), ")")) %>%
      dplyr::group_by(comparison_species) %>%
      dplyr::summarise(reasons_for_inclusion = paste(reason, collapse = "; "), .groups = 'drop')
    utils::write.csv(consolidated_details, details_csv_path, row.names = FALSE)
  }

  # Helper function to write a placeholder CSV and return a status message for the report
  write_placeholder_and_return <- function(reason, verbose_msg) {
    if (verbose) message(verbose_msg)
    placeholder_df <- tibble::tibble(
        Status = "Skipped",
        Reason = reason,
        Timestamp = as.character(Sys.time())
    )
    utils::write.csv(placeholder_df, output_csv_path, row.names = FALSE)
    # Return a list that the report generator can understand
    return(list(status_message = paste("Module skipped:", reason)))
  }

  if (is.null(similar_species_ids) || length(similar_species_ids) == 0) {
    return(write_placeholder_and_return("No statistically similar species were identified for validation.", "      Similarity Validation: No statistically similar species found. Creating placeholder file and skipping."))
  }

  data_val_prep <- data_for_validation %>%
    dplyr::filter(!!dplyr::sym(species_id_col) %in% c(target_id, similar_species_ids)) %>%
    dplyr::select(dplyr::all_of(species_id_col), dplyr::all_of(selected_vars)) %>%
    dplyr::mutate(ValidationGroup = factor(ifelse(!!dplyr::sym(species_id_col) == target_id, "Target", "Similar"))) %>%
    dplyr::select(-dplyr::all_of(species_id_col)) %>% 
    stats::na.omit()
  vars_to_factor_val <- intersect(selected_vars, categorical_vars_names)
  if(length(vars_to_factor_val) > 0) {
    data_val_prep <- data_val_prep %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(vars_to_factor_val), as.factor))
  }
  if (nrow(data_val_prep) < 20 || length(unique(data_val_prep$ValidationGroup)) < 2 || any(table(data_val_prep$ValidationGroup) < settings$validation_cv_folds)) {
    return(write_placeholder_and_return("Insufficient data for cross-validation after filtering.", "      Similarity Validation: Insufficient data or groups for cross-validation. Creating placeholder file and skipping."))
  }
  if (!requireNamespace("caret", quietly = TRUE)) {
    return(write_placeholder_and_return("Required package 'caret' is not installed.", "      Similarity Validation: Package 'caret' not found. Creating placeholder file and skipping."))
  }
  set.seed(123)
  train_control <- caret::trainControl(method = "cv", number = settings$validation_cv_folds)
  rf_formula <- stats::as.formula("ValidationGroup ~ .")
  cv_model <- tryCatch({
    caret::train(rf_formula, data = data_val_prep, method = "rf", 
                 trControl = train_control,
                 metric = "Accuracy", 
                 ntree = 999) 
  }, error = function(e) {
    # Use the placeholder function on error
    write_placeholder_and_return(paste("caret::train failed with error:", e$message), paste0("      Similarity Validation: caret::train failed: ", e$message, ". Creating placeholder file."))
    NULL
  })
  if (is.null(cv_model)) {
    # The error case in tryCatch already handles file writing, just return a status for the report
    return(list(status_message = "Module failed: caret::train returned NULL."))
  }
  mean_accuracy <- cv_model$results$Accuracy[which.max(cv_model$results$Accuracy)]
  accuracy_per_fold_df <- cv_model$resample %>% dplyr::select(Fold = Resample, Accuracy)
  if (verbose) message("    Similarity Validation mean CV Accuracy: ", round(mean_accuracy, 3))
  utils::write.csv(accuracy_per_fold_df, output_csv_path, row.names = FALSE)
  return(list(
    cv_accuracy_results = accuracy_per_fold_df,
    mean_cv_accuracy = mean_accuracy,
    n_observations_used = nrow(data_val_prep),
    similar_species_validated_against = similar_species_ids,
    rf_model_object = cv_model
  ))
}

# --- Modularity Network Module Helper Function ---

#' Run Niche Modularity Analysis Module (Internal)
#'
#' Constructs a niche overlap network and performs modularity analysis.
#'
#' @param data_for_module Tibble. Data containing species ID and selected environmental variables.
#'        Should include the target and all relevant comparison species.
#' @param target_id Character. ID of the target species.
#' @param selected_vars Character vector. Names of selected environmental variables.
#' @param species_id_col Character. Name of the species ID column in `data_for_module`.
#' @param categorical_vars_names Character vector. Names of categorical variables.
#' @param output_dir_module Character. Path to save CSV outputs and plots.
#' @param settings An object of class 'niche_analysis_settings'.
#' @param verbose Logical.
#'
#' @return A list containing the network object, modularity score, group memberships,
#'         and number of groups, or NULL if analysis fails.
#' @noRd
.run_modularity_module <- function(data_for_module, target_id,
                                   selected_vars, species_id_col, categorical_vars_names,
                                   output_dir_module, settings, verbose = TRUE) {
  if (verbose) message("    Running Niche Modularity Analysis module...")
  # Use the metric from settings, with a fallback to the original default
  modularity_overlap_metric <- settings$modularity_network_metric %||% "mean_schoener_D"
  overlap_threshold <- settings$modularity_network_overlap_threshold %||% 0.25
  if (verbose) message("      Modularity: Using '", modularity_overlap_metric, "' with threshold > ", overlap_threshold)
  all_species_in_module <- unique(data_for_module[[species_id_col]])
  if (length(all_species_in_module) < 2) {
    if (verbose) message("      Modularity: Less than 2 species in the data for module. Skipping network construction.")
    return(list(network = igraph::make_empty_graph(directed = FALSE), modularity_score = NA_real_, groups = NULL, n_groups = NA_integer_))
  }
  # --- RECALCULATE ALL-VS-ALL OVERLAPS ---
  if (verbose) message("      Modularity: Recalculating all-vs-all overlaps for network construction.") 
  all_pairs <- utils::combn(all_species_in_module, 2, simplify = FALSE) 
  n_bins_mod <- settings$overlap_n_bins
  all_pairwise_overlaps_for_modularity <- purrr::map_dfr(all_pairs, function(pair) {
    sp1_id <- pair[1]
    sp2_id <- pair[2]
    data_pair_mod <- data_for_module %>% dplyr::filter(!!dplyr::sym(species_id_col) %in% c(sp1_id, sp2_id))
    per_variable_overlaps_mod <- list()
    continuous_vars_mod <- intersect(selected_vars, setdiff(names(data_for_module), categorical_vars_names))
    categorical_vars_mod_actual <- intersect(selected_vars, categorical_vars_names)
    if (length(continuous_vars_mod) > 0) {
      for (var_cont in continuous_vars_mod) {
        d1 <- data_pair_mod[[var_cont]][data_pair_mod[[species_id_col]] == sp1_id]
        d2 <- data_pair_mod[[var_cont]][data_pair_mod[[species_id_col]] == sp2_id]
        if(length(stats::na.omit(d1))==0 || length(stats::na.omit(d2))==0) {
          per_variable_overlaps_mod[[var_cont]] <- list(schoener_D=NA_real_, pianka_O=NA_real_, czek_psi=NA_real_, hurlbert_L=NA_real_)
          next
        }
        freq1 <- .normalize_to_bins(d1, n_bins_mod); freq2 <- .normalize_to_bins(d2, n_bins_mod)
        per_variable_overlaps_mod[[var_cont]] <- list(
          schoener_D=.calculate_schoener_D_continuous(freq1,freq2),
          pianka_O=.calculate_pianka_O_continuous(freq1,freq2),
          czek_psi=.calculate_czekanowski_psi_proportions(freq1,freq2),
          hurlbert_L=.calculate_hurlbert_L_proportions(freq1,freq2)
        )
      }
    }
    if (length(categorical_vars_mod_actual) > 0) {
      for (var_cat in categorical_vars_mod_actual) {
        per_variable_overlaps_mod[[var_cat]] <- .calculate_categorical_var_overlap(data_pair_mod, species_id_col, var_cat, sp1_id, sp2_id)
      }
    }
    
    if (length(per_variable_overlaps_mod) == 0) {
      return(tibble::tibble(from=sp1_id, to=sp2_id, weight=NA_real_, mean_schoener_D=NA_real_, mean_pianka_O=NA_real_, mean_czek_psi=NA_real_, mean_hurlbert_L=NA_real_, composite_overlap=NA_real_))
    }
    mean_schoener_mod <- mean(purrr::map_dbl(per_variable_overlaps_mod, "schoener_D",.default=NA_real_),na.rm=TRUE)
    mean_pianka_mod <- mean(purrr::map_dbl(per_variable_overlaps_mod, "pianka_O",.default=NA_real_),na.rm=TRUE)
    mean_czek_mod <- mean(purrr::map_dbl(per_variable_overlaps_mod, "czek_psi",.default=NA_real_),na.rm=TRUE)
    mean_hurlbert_mod <- mean(purrr::map_dbl(per_variable_overlaps_mod, "hurlbert_L",.default=NA_real_),na.rm=TRUE)
    composite_mod <- 0.4 * mean_schoener_mod + 0.3 * mean_czek_mod + 0.3 * mean_hurlbert_mod
    tibble::tibble(from = sp1_id, to = sp2_id,
                   mean_schoener_D=mean_schoener_mod, mean_pianka_O=mean_pianka_mod,
                   mean_czek_psi=mean_czek_mod, mean_hurlbert_L=mean_hurlbert_mod,
                   composite_overlap=composite_mod)
  })
  edges_df <- all_pairwise_overlaps_for_modularity %>%
    # Dynamically rename the selected metric column to 'weight' for igraph
    dplyr::rename(weight = !!dplyr::sym(modularity_overlap_metric)) %>%
    dplyr::filter(!is.na(weight) & weight > overlap_threshold)
  if (nrow(edges_df) == 0) {
    if (verbose) message("      Modularity: No species pairs found above the overlap threshold. Network will be empty.")
  }
  nodes_df <- data.frame(name = all_species_in_module)
  niche_net <- igraph::graph_from_data_frame(d = edges_df, vertices = nodes_df, directed = FALSE)
  if (igraph::gsize(niche_net) == 0) {
    if (verbose) message("      Modularity: Network has no edges. Modularity analysis not meaningful.")
    return(list(network = niche_net, modularity_score = NA_real_, groups = stats::setNames(rep(1, length(all_species_in_module)), all_species_in_module), n_groups = 1, plot_path = NULL))
  }
  communities <- tryCatch({
    igraph::cluster_walktrap(niche_net, weights = if("weight" %in% igraph::edge_attr_names(niche_net)) igraph::E(niche_net)$weight else NULL)
  }, error = function(e) {
    if (verbose) warning("      Modularity: Community detection (cluster_walktrap) failed: ", e$message)
    NULL
  })
  if (is.null(communities)) return(list(network = niche_net, modularity_score = NA_real_, groups = NULL, n_groups = NA_integer_))
  modularity_score_val <- igraph::modularity(communities)
  group_membership <- igraph::membership(communities)
  if (verbose) message("      Modularity score: ", round(modularity_score_val, 3), ", Number of groups: ", length(unique(group_membership)))
  groups_df <- tibble::tibble(species = names(group_membership), module_id = as.integer(group_membership)) %>%
    dplyr::arrange(module_id, species)
  utils::write.csv(groups_df, file.path(output_dir_module, paste0("modularity_groups_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".csv")), row.names = FALSE)
  return(list(
    network = niche_net,
    modularity_score = modularity_score_val, 
    groups = group_membership, 
    n_groups = length(unique(group_membership))
  ))
}

#' Run Niche Equivalency Test Module (Internal)
#'
#' Performs a niche equivalency test between a target species and comparison species
#' based on centroid distances in standardized environmental space.
#'
#' @param data_for_module Tibble. Data containing species ID and selected environmental variables.
#' @param target_id Character. ID of the target species.
#' @param comparison_ids Character vector. IDs of comparison species.
#' @param selected_vars Character vector. Names of selected environmental variables.
#' @param species_id_col Character. Name of the species ID column in `data_for_module`.
#' @param categorical_vars_names Character vector. Names of categorical variables (to be excluded).
#' @param output_dir_module Character. Path to save the CSV output.
#' @param settings An object of class 'niche_analysis_settings'.
#' @param verbose Logical.
#'
#' @return A tibble with pairwise equivalency test results.
#' @noRd
.run_equivalency_module <- function(data_for_module, target_id, comparison_ids,
                                    selected_vars, species_id_col, categorical_vars_names,
                                    output_dir_module, settings, verbose = TRUE) {
  
  if (verbose) message("    Calculating niche equivalency tests...")
  numeric_vars <- setdiff(selected_vars, categorical_vars_names)
  if (length(numeric_vars) < 1) {
    if (verbose) message("      Niche Equivalency: No numeric variables selected. Skipping module.")
    return(tibble::tibble(
      target_species = target_id,
      comparison_species = comparison_ids,
      centroid_distance = NA_real_,
      p_value = NA_real_,
      n_permutations = NA_integer_,
      variables_used = paste(numeric_vars, collapse = "; "),
      notes = "No numeric variables available"
    ))
  }
  n_permutations <- settings$n_bootstrap_stability %||% 999
  all_equivalency_results <- purrr::map_dfr(comparison_ids, function(comp_id) {
    if (verbose) message("      Equivalency: ", target_id, " vs ", comp_id)
    sp1_data_raw <- data_for_module %>%
      dplyr::filter(!!dplyr::sym(species_id_col) == target_id) %>%
      dplyr::select(dplyr::all_of(numeric_vars)) %>%
      stats::na.omit()
    sp2_data_raw <- data_for_module %>%
      dplyr::filter(!!dplyr::sym(species_id_col) == comp_id) %>%
      dplyr::select(dplyr::all_of(numeric_vars)) %>%
      stats::na.omit()
    if (nrow(sp1_data_raw) < 2 || nrow(sp2_data_raw) < 2) {
      return(tibble::tibble(target_species = target_id, comparison_species = comp_id, centroid_distance = NA_real_, p_value = NA_real_, n_permutations = n_permutations, notes = "Insufficient data for one or both species"))
    }
    combined_data_for_scaling <- rbind(sp1_data_raw, sp2_data_raw)
    sd_check <- apply(combined_data_for_scaling, 2, stats::sd, na.rm = TRUE)
    if (any(sd_check == 0, na.rm = TRUE)) {
      return(tibble::tibble(target_species = target_id, comparison_species = comp_id, centroid_distance = NA_real_, p_value = NA_real_, n_permutations = n_permutations, notes = "Zero variance in one or more variables"))
    }
    combined_std <- scale(combined_data_for_scaling, center = TRUE, scale = TRUE)
    sp1_std <- combined_std[1:nrow(sp1_data_raw), , drop = FALSE]
    sp2_std <- combined_std[(nrow(sp1_data_raw) + 1):nrow(combined_std), , drop = FALSE]
    centroid1 <- colMeans(sp1_std, na.rm = TRUE)
    centroid2 <- colMeans(sp2_std, na.rm = TRUE)
    observed_distance <- stats::dist(rbind(centroid1, centroid2))[1]
    perm_distances <- purrr::map_dbl(1:n_permutations, function(i) {
      permuted_indices <- sample(nrow(combined_std))
      p1 <- combined_std[permuted_indices[1:nrow(sp1_std)], , drop = FALSE]
      p2 <- combined_std[permuted_indices[(nrow(sp1_std) + 1):nrow(combined_std)], , drop = FALSE]
      stats::dist(rbind(colMeans(p1, na.rm = TRUE), colMeans(p2, na.rm = TRUE)))[1]
    })
    p_value <- (sum(perm_distances >= observed_distance, na.rm = TRUE) + 1) / (n_permutations + 1)
    
    tibble::tibble(target_species = target_id, comparison_species = comp_id,
                   centroid_distance = observed_distance, p_value = p_value, n_permutations = n_permutations, notes = NA_character_)
  })
  all_equivalency_results <- all_equivalency_results %>%
    dplyr::mutate(variables_used = paste(numeric_vars, collapse = "; "))
  csv_filename <- file.path(output_dir_module, paste0("equivalency_test_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".csv"))
  utils::write.csv(all_equivalency_results, csv_filename, row.names = FALSE)
  if (verbose) message("    Equivalency test results saved to: ", csv_filename)
  return(all_equivalency_results)
}


#' Run Boruta Variable Selection Module (Internal)
#'
#' Performs variable selection using the Boruta algorithm.
#'
#' @param data_for_boruta Tibble. Data containing the 'species' column (factor)
#'        and the numeric environmental variables to be tested.
#' @param target_id Character. The ID of the current target species.
#' @param potential_vars Character vector. Names of numeric environmental variables to test.
#' @param species_col_actual_name Character. The actual name of the species ID column in `data_for_boruta`.
#' @param settings An object of class 'niche_analysis_settings'.
#' @param output_dir_boruta Character. Path to save Boruta plots and stats.
#' @param save_outputs Logical. If TRUE, saves plot and stats CSV.
#' @param apply_max_confirmed_limit Logical. If TRUE, applies `settings$boruta_max_confirmed_vars`.
#' @param metadata_labels List. Optional. A list containing user-defined labels for variables, typically from `niche_data_object$metadata$labels`.
#' @param verbose Logical. If TRUE, prints progress messages.
#'
#' @return A list containing:
#'   \item{confirmed_variables}{Character vector of confirmed variable names.}
#'   \item{tentative_variables}{Character vector of tentative variable names.}
#'   \item{boruta_stats_df}{Dataframe of Boruta statistics.}
#'   \item{boruta_plot_path}{Path to the saved Boruta plot.}
#' @noRd
.run_boruta_selection_module <- function(data_for_boruta,
                                         target_id,
                                         potential_vars,
                                         species_col_actual_name,
                                         settings,
                                         output_dir_boruta,
                                         save_outputs = TRUE,
                                         apply_max_confirmed_limit = TRUE,
                                         metadata_labels = NULL,
                                         verbose = TRUE,
                                         max_runs_override = NULL) { 
  if (verbose && save_outputs) message("      Executing Boruta for target: ", target_id)
  
  if(!species_col_actual_name %in% names(data_for_boruta)){
    warning("      Boruta: Species ID column not found in data_for_boruta. Skipping.")
    return(list(confirmed_variables = potential_vars, tentative_variables = character(0), boruta_stats_df = NULL, boruta_plot_path = NULL))
  }
  boruta_data_prep <- data_for_boruta %>%
    dplyr::select(dplyr::all_of(species_col_actual_name), dplyr::all_of(potential_vars)) %>%
    stats::na.omit()
  if (nrow(boruta_data_prep) < 20) {
    if (verbose) warning("      Boruta: Less than 20 complete rows after NA omit for target '", target_id, "'. Returning all potential vars.")
    return(list(confirmed_variables = potential_vars, tentative_variables = character(0), boruta_stats_df = NULL, boruta_plot_path = NULL))
  }
  boruta_data_analysis <- boruta_data_prep %>%
    dplyr::mutate(TargetGroup = factor(ifelse(!!dplyr::sym(species_col_actual_name) == target_id, "Target", "Other"))) %>%
    dplyr::select(-dplyr::all_of(species_col_actual_name)) 
  group_counts <- table(boruta_data_analysis$TargetGroup)
  if (length(group_counts) < 2 || any(group_counts < 5)) {
    if (verbose) warning("      Boruta: Insufficient or imbalanced groups for target '", target_id, "' (counts: ", paste(names(group_counts), group_counts, collapse=", "), "). Returning all potential vars.")
    return(list(confirmed_variables = potential_vars, tentative_variables = character(0), boruta_stats_df = NULL, boruta_plot_path = NULL))
  }
  original_names_map <- stats::setNames(potential_vars, make.names(potential_vars, unique = TRUE))
  
  data_boruta_renamed <- boruta_data_analysis
  # Ensure correct sanitized names are used for Boruta formula
  sanitized_predictor_names <- names(original_names_map)[match(potential_vars, original_names_map)]
  colnames(data_boruta_renamed) <- c(sanitized_predictor_names, "TargetGroup")
  
  boruta_formula <- stats::as.formula("TargetGroup ~ .")
  actual_max_runs <- max_runs_override %||% settings$boruta_internal_max_runs %||% 100
  boruta_run <- tryCatch({
    Boruta::Boruta(boruta_formula, data = data_boruta_renamed, doTrace = 0,
                   maxRuns = actual_max_runs)
  }, error = function(e) {
    if (verbose) warning("      Boruta execution failed for target '", target_id, "': ", e$message)
    NULL
  })
  if (is.null(boruta_run)) {
    return(list(confirmed_variables = potential_vars, tentative_variables = character(0), boruta_stats_df = NULL, boruta_plot_path = NULL))
  }
  # Attempt to fix missing rownames in ImpHistory, which can occur in some cases
  if (!is.null(boruta_run$ImpHistory)) {
    current_rownames <- rownames(boruta_run$ImpHistory)
    if (is.null(current_rownames) || all(current_rownames == "")) {
      expected_rows <- length(sanitized_predictor_names) + 3 # For shadow vars
      actual_rows <- nrow(boruta_run$ImpHistory)
      if (actual_rows == expected_rows) {
        rownames(boruta_run$ImpHistory) <- c(sanitized_predictor_names, "shadowMin", "shadowMean", "shadowMax")
      } else if (actual_rows == length(sanitized_predictor_names)) {
        rownames(boruta_run$ImpHistory) <- sanitized_predictor_names
      }
    }
  }
  
  
  if (verbose && save_outputs) message("      Main Boruta run completed for target: ", target_id)
  confirmed_vars_sanitized <- Boruta::getSelectedAttributes(boruta_run, withTentative = FALSE)
  tentative_vars_sanitized <- Boruta::getSelectedAttributes(boruta_run, withTentative = TRUE)
  confirmed_variables <- original_names_map[confirmed_vars_sanitized]
  tentative_variables <- original_names_map[tentative_vars_sanitized]
  confirmed_variables <- confirmed_variables[!is.na(confirmed_variables)] 
  tentative_variables <- tentative_variables[!is.na(tentative_variables)]
  boruta_stats_raw <- Boruta::attStats(boruta_run)
  boruta_stats_df <- as.data.frame(boruta_stats_raw) %>%
    tibble::rownames_to_column("Sanitized_Variable") %>%
    dplyr::mutate(Original_Variable = original_names_map[Sanitized_Variable]) %>%
    dplyr::select(Original_Variable, Sanitized_Variable, dplyr::everything()) %>%
    dplyr::arrange(dplyr::desc(meanImp))
  boruta_plot_path <- NULL
  if (save_outputs) {
    
    plot_boruta_obj_for_plot <- boruta_run
    
    if (!is.null(plot_boruta_obj_for_plot$ImpHistory) && ncol(plot_boruta_obj_for_plot$ImpHistory) > 0) {
      
      sanitized_var_names_from_boruta <- rownames(plot_boruta_obj_for_plot$ImpHistory)
      # Only proceed with label mapping if sanitized_var_names_from_boruta is not empty and seems valid
      if (length(sanitized_var_names_from_boruta) > 0 && !all(sanitized_var_names_from_boruta == "")) {
        final_labels_for_plot <- character(length(sanitized_var_names_from_boruta))
        user_labels_applied_count <- 0
                
        for(i in seq_along(sanitized_var_names_from_boruta)){
          s_name <- sanitized_var_names_from_boruta[i]
          original_name <- original_names_map[s_name] # Mappa sanitizzato -> originale
          # Default iniziale all'originale (o al sanitizzato se la mappa originale fallisce)        
          current_label <- if (!is.na(original_name) && nzchar(original_name)) original_name else s_name
          
          # Tenta di sovrascrivere con l'etichetta utente
          if (!is.null(metadata_labels) && !is.null(metadata_labels$continuous_vars) &&
              !is.na(original_name) && nzchar(original_name) && # Assicurati che original_name sia valido per il lookup
              original_name %in% names(metadata_labels$continuous_vars)) {
            user_label_candidate <- metadata_labels$continuous_vars[[original_name]]
            if (!is.na(user_label_candidate) && nzchar(user_label_candidate)) {
              current_label <- user_label_candidate
              if (current_label != original_name && (!is.na(original_name) && nzchar(original_name))) { # Conta solo se l'etichetta utente è diversa dall'originale VALIDO
                user_labels_applied_count <- user_labels_applied_count + 1
              }
            }
          }
          final_labels_for_plot[i] <- current_label
        }
        
        # Controllo di unicità e fallback
        if (length(final_labels_for_plot) == length(unique(final_labels_for_plot))) {
          rownames(plot_boruta_obj_for_plot$ImpHistory) <- final_labels_for_plot
          if (verbose) {
            if (user_labels_applied_count > 0) {
              message("      Boruta plot: Successfully applied user-defined labels for target '", target_id, "' to Boruta object.")
            } else {
              message("      Boruta plot: Using original variable names for target '", target_id, "' in Boruta object (user-defined labels not found or not different).")
            }
          }
        } else if (length(final_labels_for_plot) > 0) { 
          original_names_from_map_direct <- as.vector(original_names_map[sanitized_var_names_from_boruta])
          original_names_from_map_direct[is.na(original_names_from_map_direct) | !nzchar(original_names_from_map_direct)] <- 
            sanitized_var_names_from_boruta[is.na(original_names_from_map_direct) | !nzchar(original_names_from_map_direct)] 
          
          if (length(original_names_from_map_direct) == length(unique(original_names_from_map_direct))) {
            rownames(plot_boruta_obj_for_plot$ImpHistory) <- original_names_from_map_direct
            if (verbose) warning(paste("      Boruta plot: User-defined labels for target '", target_id, "' resulted in non-unique names.",
                                       "Using original variable names for Boruta object plot labels instead."))
          } else {
            rownames(plot_boruta_obj_for_plot$ImpHistory) <- sanitized_var_names_from_boruta
            if (verbose) warning(paste("      Boruta plot: Neither user-defined nor original variable names were unique for target '", target_id, "'.",
                                       "Using sanitized (formula-compatible) names in the Boruta object plot."))
          }       
        }
      } else { # End of: if (length(sanitized_var_names_from_boruta) > 0 && !all(sanitized_var_names_from_boruta == ""))
        if (verbose) warning("      Boruta Plot: rownames(ImpHistory) are missing or invalid. Standard Boruta plot might not use custom labels. Attempting ggplot fallback if standard plot fails.")
      }
    }
    plot_filename <- paste0("boruta_plot_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".png")
    boruta_plot_path <- file.path(output_dir_boruta, plot_filename)
    plot_generation_successful <- FALSE
    
    # Attempt standard Boruta plot first
    if (length(rownames(plot_boruta_obj_for_plot$ImpHistory)) > 0 && !all(rownames(plot_boruta_obj_for_plot$ImpHistory) == "")) {
      grDevices::png(boruta_plot_path, width = 1000, height = 700, res = 100) 
      old_par_settings <- NULL
      tryCatch({ 
        old_par_settings <- graphics::par(mar = c(10, 12, 4, 2) + 0.1) 
        plot(plot_boruta_obj_for_plot, xlab = "", 
             main = paste("Boruta Variable Importance for", target_id),
             las = 2, cex.axis = 0.8)        
        plot_generation_successful <- TRUE
      }, error = function(e_plot) {
        if (verbose) warning("      Error plotting Boruta (standard method) for target '", target_id, "': ", e_plot$message)
        # Ensure device is off if plot fails
        # tryCatch({ if(!is.null(old_par_settings)) graphics::par(old_par_settings) }, error = function(e_par_reset){})
        # grDevices::dev.off() # dev.off() will be called in finally
      }, finally = {
        if (!is.null(old_par_settings)) {
          try(graphics::par(old_par_settings), silent = TRUE) 
        }
        grDevices::dev.off() 
      })
    } else {
      if (verbose) warning("      Skipping standard Boruta plot due to missing/invalid ImpHistory rownames.")
    }
    
    # If standard plot failed or was skipped, and we have boruta_stats_df, try ggplot
    if (!plot_generation_successful && !is.null(boruta_stats_df) && nrow(boruta_stats_df) > 0 && requireNamespace("ggplot2", quietly = TRUE)) {
      if (verbose) message("      Attempting ggplot2 fallback for Boruta plot for target '", target_id, "'.")
      
      plot_data_ggplot <- boruta_stats_df %>%
        dplyr::mutate(
          VariableLabel = dplyr::coalesce(
            if (!is.null(metadata_labels) && !is.null(metadata_labels$continuous_vars)) metadata_labels$continuous_vars[Original_Variable] else NULL,
            Original_Variable,
            Sanitized_Variable
          ),
          Decision = decision # Assuming 'decision' column exists from attStats
        ) %>%
        dplyr::filter(!is.na(meanImp)) %>% # Ensure meanImp is not NA
        dplyr::arrange(meanImp) %>%
        dplyr::mutate(VariableLabel = factor(VariableLabel, levels = VariableLabel)) # Order by importance
      
      if (nrow(plot_data_ggplot) > 0) {
        p_gg <- ggplot2::ggplot(plot_data_ggplot, ggplot2::aes(x = VariableLabel, y = meanImp, fill = Decision)) +
          ggplot2::geom_bar(stat = "identity") +
          ggplot2::coord_flip() +
          ggplot2::scale_fill_manual(values = c("Confirmed" = "darkgreen", "Tentative" = "orange", "Rejected" = "darkred", "ShadowMax" = "grey50", "ShadowMean" = "grey70", "ShadowMin" = "grey90"),
                                     limits = c("Confirmed", "Tentative", "Rejected", "ShadowMax", "ShadowMean", "ShadowMin"), drop = FALSE) +
          ggplot2::labs(title = paste("Boruta Variable Importance (ggplot2) for", target_id),
                        x = "Variable", y = "Mean Importance (MeanDecreaseAccuracy)") +
          ggplot2::theme_minimal(base_size = 10) +
          ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8), legend.position = "bottom")
        
        tryCatch({
          ggplot2::ggsave(boruta_plot_path, plot = p_gg, width = 10, height = max(7, nrow(plot_data_ggplot) * 0.3), dpi = 100, limitsize = FALSE)
          plot_generation_successful <- TRUE
          if (verbose) message("      Boruta plot (ggplot2) saved to: ", boruta_plot_path)
        }, error = function(e_ggsave) {
          if (verbose) warning("      Error saving ggplot2 Boruta plot for target '", target_id, "': ", e_ggsave$message)
        })
      } else {
        if (verbose) warning("      No data for ggplot2 Boruta plot after filtering for target '", target_id, "'.")
      }
    }
    # Final fallback if all plotting fails
    if (!plot_generation_successful) {
      tryCatch({
        grDevices::png(boruta_plot_path, width = 1000, height = 700, res = 100)
        graphics::plot(1, type = "n", xlab = "", ylab = "", axes = FALSE, main = "Boruta Plot Error")
        graphics::text(1, 1, paste("Failed to generate Boruta plot for", target_id), cex = 0.9)
      }, error = function(e_final_fallback) {}, finally = { grDevices::dev.off() })
    }
    
    if (plot_generation_successful && verbose) {
      message("      Boruta plot saved to: ", boruta_plot_path)
    } else if (verbose) {
      message("      Boruta plot for '", target_id, "' could not be generated or an error occurred. Check: ", boruta_plot_path)
    }
    stats_filename <- paste0("boruta_stats_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".csv")
    utils::write.csv(boruta_stats_df, file.path(output_dir_boruta, stats_filename), row.names = FALSE)
    if (verbose) message("      Boruta stats saved to: ", file.path(output_dir_boruta, stats_filename))
  }
  
  if (apply_max_confirmed_limit) {
    if (!is.null(settings$boruta_max_confirmed_vars) && length(confirmed_variables) > settings$boruta_max_confirmed_vars) {
      if (verbose && save_outputs) message("      Limiting confirmed variables to top ", settings$boruta_max_confirmed_vars, " based on mean importance.")
      
      if (!is.null(boruta_stats_df) && nrow(boruta_stats_df) > 0) {
        confirmed_stats_for_ranking <- boruta_stats_df %>%
          dplyr::filter(Original_Variable %in% confirmed_variables)
        if (nrow(confirmed_stats_for_ranking) > 0) {
          vars_with_imp <- confirmed_stats_for_ranking %>%
            dplyr::filter(!is.na(meanImp)) %>%
            dplyr::arrange(dplyr::desc(meanImp))
          vars_without_imp <- confirmed_stats_for_ranking %>%
            dplyr::filter(is.na(meanImp))
          selected_vars_list <- head(vars_with_imp$Original_Variable, settings$boruta_max_confirmed_vars)
          num_still_needed <- settings$boruta_max_confirmed_vars - length(selected_vars_list)
          if (num_still_needed > 0 && nrow(vars_without_imp) > 0) {
            additional_vars_to_consider <- intersect(confirmed_variables, vars_without_imp$Original_Variable)
            selected_vars_list <- c(selected_vars_list, head(additional_vars_to_consider, num_still_needed))
          }
          confirmed_variables <- unique(selected_vars_list)
        } else {
          if(verbose) warning("      Could not find stats for initially confirmed Boruta variables. Resulting selection might be empty.")
          confirmed_variables <- character(0)
        }
      } else {
        if(verbose) warning("      Cannot apply boruta_max_confirmed_vars: boruta_stats_df is NULL or empty. Taking first N confirmed arbitrarily.")
        current_confirmed_count <- length(confirmed_variables)
        confirmed_variables <- head(confirmed_variables, settings$boruta_max_confirmed_vars)
      }
    }
  }
  return(list(
    confirmed_variables = confirmed_variables,
    tentative_variables = tentative_variables,
    boruta_stats_df = boruta_stats_df,
    boruta_plot_path = boruta_plot_path
  ))
}


#' Run Random Forest Importance Selection Module (Internal)
#'
#' Performs variable selection based on Random Forest variable importance.
#'
#' @param data_for_rf Tibble. Data containing the 'species' column (factor)
#'        and the numeric environmental variables to be tested.
#' @param target_id Character. The ID of the current target species.
#' @param potential_vars Character vector. Names of numeric environmental variables to test.
#' @param species_col_actual_name Character. The actual name of the species ID column in `data_for_rf`.
#' @param settings An object of class 'niche_analysis_settings'.
#' @param output_dir_rf Character. Path to save RF importance scores.
#' @param verbose Logical. If TRUE, prints progress messages.
#'
#' @return A list containing:
#'   \item{selected_variables}{Character vector of top N important variable names.}
#'   \item{importance_scores}{Dataframe of RF importance scores for all tested variables.}
#'   \item{rf_model_oob_error}{Out-of-Bag error rate of the RF model.}
#' @noRd
.run_rf_importance_module <- function(data_for_rf,
                                      target_id,
                                      potential_vars,
                                      species_col_actual_name,
                                      settings,
                                      output_dir_rf,
                                      verbose = TRUE) {
  if (verbose) message("      Executing Random Forest Importance for target: ", target_id)
  if(!species_col_actual_name %in% names(data_for_rf)){
    warning("      RF Importance: Species ID column not found in data_for_rf. Skipping.")
    return(list(selected_variables = potential_vars, importance_scores = NULL, rf_model_oob_error = NA_real_))
  }
  rf_data_prep <- data_for_rf %>%
    dplyr::select(dplyr::all_of(species_col_actual_name), dplyr::all_of(potential_vars)) %>%
    stats::na.omit()
  if (nrow(rf_data_prep) < 20) {
    if (verbose) warning("      RF Importance: Less than 20 complete rows after NA omit for target '", target_id, "'. Returning all potential vars.")
    return(list(selected_variables = potential_vars, importance_scores = NULL, rf_model_oob_error = NA_real_))
  }
  rf_data_analysis <- rf_data_prep %>%
    dplyr::mutate(TargetGroup = factor(ifelse(!!dplyr::sym(species_col_actual_name) == target_id, "Target", "Other"))) %>%
    dplyr::select(-dplyr::all_of(species_col_actual_name)) 
  group_counts <- table(rf_data_analysis$TargetGroup)
  if (length(group_counts) < 2 || any(group_counts < 5)) {
    if (verbose) warning("      RF Importance: Insufficient or imbalanced groups for target '", target_id, "' (counts: ", paste(names(group_counts), group_counts, collapse=", "), "). Returning all potential vars.")
    return(list(selected_variables = potential_vars, importance_scores = NULL, rf_model_oob_error = NA_real_))
  }
  rf_formula <- stats::as.formula("TargetGroup ~ .")
  rf_model <- tryCatch({
    randomForest::randomForest(rf_formula, data = rf_data_analysis, ntree = 500, importance = TRUE, proximity = FALSE)
  }, error = function(e) {
    if (verbose) warning("      Random Forest execution failed for target '", target_id, "': ", e$message)
    NULL
  })
  if (is.null(rf_model)) {
    return(list(selected_variables = potential_vars, importance_scores = NULL, rf_model_oob_error = NA_real_))
  }
  if (verbose) message("      Random Forest completed for target: ", target_id)
  # Get importance scores (MeanDecreaseAccuracy or MeanDecreaseGini)
  importance_raw <- randomForest::importance(rf_model, type = 1) 
  importance_scores_df <- as.data.frame(importance_raw) %>%
    tibble::rownames_to_column("Variable") %>%
    dplyr::rename(MeanDecreaseAccuracy = MeanDecreaseAccuracy) %>% 
    dplyr::arrange(dplyr::desc(MeanDecreaseAccuracy))
  stats_filename <- paste0("rf_importance_scores_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".csv")
  utils::write.csv(importance_scores_df, file.path(output_dir_rf, stats_filename), row.names = FALSE)
  if (verbose) message("      RF importance scores saved to: ", file.path(output_dir_rf, stats_filename))
  n_top <- settings$rf_importance_n_top %||% length(potential_vars)
  selected_variables <- head(importance_scores_df$Variable, n_top)
  return(list(
    selected_variables = selected_variables,
    importance_scores = importance_scores_df,
    rf_model_oob_error = rf_model$err.rate[nrow(rf_model$err.rate), "OOB"]
  ))
}

#' Filter Variables by Correlation (Internal)
#'
#' Removes highly correlated variables from a given set.
#'
#' @param data_for_correlation Tibble. Data containing the numeric variables.
#' @param selected_vars Character vector. Names of numeric variables to check for correlation.
#' @param correlation_threshold Numeric. Absolute correlation threshold above which
#'        one variable from a pair will be removed.
#' @param output_dir_corr Character. Path to save the correlation matrix.
#' @param target_id Character. ID of the target species, for naming output files.
#' @param verbose Logical. If TRUE, prints messages.
#'
#' @return A list containing:
#'   \item{remaining_variables}{Character vector of variable names after filtering.}
#'   \item{correlation_matrix}{The calculated correlation matrix.}
#'   \item{removed_variables}{Character vector of variables removed due to high correlation.}
#' @noRd
.filter_by_correlation <- function(data_for_correlation,
                                   selected_vars,
                                   correlation_threshold,
                                   output_dir_corr,
                                   target_id,
                                   verbose = TRUE) {
  if (verbose) message("      Applying correlation filter (threshold: ", correlation_threshold, ")")
  if (length(selected_vars) < 2) {
    if (verbose) message("      Less than 2 variables, skipping correlation filter.")
    return(list(remaining_variables = selected_vars, correlation_matrix = NULL, removed_variables = character(0)))
  }
  data_numeric <- data_for_correlation %>%
    dplyr::select(dplyr::all_of(selected_vars)) %>%
    stats::na.omit()
  if (nrow(data_numeric) < 2) {
    if (verbose) warning("      Correlation filter: Not enough complete observations after NA removal. Skipping.")
    return(list(remaining_variables = selected_vars, correlation_matrix = NULL, removed_variables = character(0)))
  }
  cor_matrix <- stats::cor(data_numeric, use = "pairwise.complete.obs")
  cor_matrix_path <- file.path(output_dir_corr, paste0("correlation_matrix_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".csv"))
  utils::write.csv(cor_matrix, cor_matrix_path, row.names = TRUE)
  if (verbose) message("      Correlation matrix saved to: ", cor_matrix_path)
  vars_to_keep <- selected_vars
  vars_removed_by_cor <- character(0)
  # Iteratively remove variables
  # Find the pair with the highest correlation, remove one, recalculate, repeat
  # A simpler approach: find all pairs above threshold, then decide which one of each pair to remove
  # based on mean absolute correlation with other variables.
  if(requireNamespace("caret", quietly = TRUE)){
    highly_correlated_indices <- caret::findCorrelation(cor_matrix, cutoff = correlation_threshold, verbose = FALSE, names = FALSE)
    if(length(highly_correlated_indices) > 0){
      vars_removed_by_cor <- colnames(cor_matrix)[highly_correlated_indices]
      vars_to_keep <- setdiff(selected_vars, vars_removed_by_cor)
      if (verbose) message("      Variables removed by caret::findCorrelation: ", paste(vars_removed_by_cor, collapse = ", "))
    } else {
      if (verbose) message("      No variables removed by caret::findCorrelation.")
    }
  } else {
    if (verbose) warning("      Package 'caret' not found. Using a simpler correlation filter. Results might differ.")
    temp_cor_matrix <- cor_matrix
    diag(temp_cor_matrix) <- 0
    while(max(abs(temp_cor_matrix), na.rm = TRUE) > correlation_threshold && length(vars_to_keep) > 1) {
      max_cor_val <- max(abs(temp_cor_matrix), na.rm = TRUE)
      if(!is.finite(max_cor_val) || max_cor_val <= correlation_threshold) break
      pair_indices <- which(abs(temp_cor_matrix) == max_cor_val, arr.ind = TRUE)[1,]
      var1_name <- rownames(temp_cor_matrix)[pair_indices[1]]
      var2_name <- colnames(temp_cor_matrix)[pair_indices[2]]
      mean_abs_cor_var1 <- mean(abs(temp_cor_matrix[var1_name, vars_to_keep[vars_to_keep != var1_name]]), na.rm = TRUE)
      mean_abs_cor_var2 <- mean(abs(temp_cor_matrix[var2_name, vars_to_keep[vars_to_keep != var2_name]]), na.rm = TRUE)
      var_to_remove <- if (mean_abs_cor_var1 >= mean_abs_cor_var2) var1_name else var2_name
      if (verbose) message("      High correlation (", round(max_cor_val,3), ") between ", var1_name, " and ", var2_name, ". Removing: ", var_to_remove)
      vars_to_keep <- setdiff(vars_to_keep, var_to_remove)
      vars_removed_by_cor <- c(vars_removed_by_cor, var_to_remove)
      if(length(vars_to_keep) < 2) break
      temp_cor_matrix <- stats::cor(data_numeric[, vars_to_keep, drop=FALSE], use = "pairwise.complete.obs")
      diag(temp_cor_matrix) <- 0
    }
  }
  return(list(
    remaining_variables = vars_to_keep,
    correlation_matrix = cor_matrix,
    removed_variables = unique(vars_removed_by_cor)
  ))
}

#' Filter Variables by VIF (Variance Inflation Factor) (Internal)
#'
#' Iteratively removes variables with VIF above a specified threshold.
#'
#' @param data_for_vif Tibble. Data containing the numeric variables.
#' @param selected_vars Character vector. Names of numeric variables to check for VIF.
#' @param vif_threshold Numeric. VIF threshold above which variables will be removed.
#' @param output_dir_vif Character. Path to save VIF values.
#' @param target_id Character. ID of the target species, for naming output files.
#' @param verbose Logical. If TRUE, prints messages.
#'
#' @return A list containing:
#'   \item{remaining_variables}{Character vector of variable names after VIF filtering.}
#'   \item{vif_values}{Dataframe of final VIF values for remaining variables.}
#'   \item{removed_variables}{Character vector of variables removed due to high VIF.}
#' @noRd
.filter_by_vif <- function(data_for_vif,
                           selected_vars,
                           vif_threshold,
                           output_dir_vif,
                           target_id,
                           verbose = TRUE) {
  if (verbose) message("      Applying VIF filter (threshold: ", vif_threshold, ")")
  if (length(selected_vars) < 2) {
    if (verbose) message("      Less than 2 variables, skipping VIF filter.")
    return(list(remaining_variables = selected_vars, vif_values = NULL, removed_variables = character(0)))
  }


  
  data_numeric <- data_for_vif %>%
    dplyr::select(dplyr::all_of(selected_vars)) %>%
    stats::na.omit()
  if (nrow(data_numeric) <= length(selected_vars)) {
    if (verbose) warning("      VIF filter: Not enough observations (", nrow(data_numeric), ") for the number of variables (", length(selected_vars), ") after NA removal. Skipping VIF.")
    return(list(remaining_variables = selected_vars, vif_values = NULL, removed_variables = character(0)))
  }
  if(requireNamespace("usdm", quietly = TRUE)){
    vif_result <- tryCatch({
      usdm::vifstep(data_numeric, th = vif_threshold)
    }, error = function(e){
      if(verbose) warning("      Error during usdm::vifstep: ", e$message, ". Skipping VIF filter.")
      NULL
    })
    if(!is.null(vif_result) && inherits(vif_result, "VIF")){
      remaining_vars_vif <- vif_result@variables
      removed_vars_vif <- setdiff(selected_vars, remaining_vars_vif)
      final_vif_values_df <- vif_result@results 
      if (verbose && length(removed_vars_vif) > 0) message("      Variables removed by usdm::vifstep: ", paste(removed_vars_vif, collapse = ", "))
      if (verbose && length(removed_vars_vif) == 0) message("      No variables removed by usdm::vifstep.")
      vif_path <- file.path(output_dir_vif, paste0("vif_values_", gsub("[^A-Za-z0-9_.-]", "_", target_id), ".csv"))
      utils::write.csv(final_vif_values_df, vif_path, row.names = FALSE)
      if (verbose) message("      Final VIF values saved to: ", vif_path)
      return(list(
        remaining_variables = remaining_vars_vif,
        vif_values = final_vif_values_df,
        removed_variables = removed_vars_vif
      ))
    } else {
      if (verbose) warning("      usdm::vifstep did not complete successfully. Skipping VIF filter.")
      return(list(remaining_variables = selected_vars, vif_values = NULL, removed_variables = character(0)))
    }
  } else {
    if (verbose) warning("      Package 'usdm' not available for iterative VIF. Skipping VIF filter. Please install 'usdm' for this feature.")
    return(list(remaining_variables = selected_vars, vif_values = NULL, removed_variables = character(0)))
  }
}

#' Extract Genus and Species from Scientific Names (Internal)
#'
#' Splits scientific names into genus and species epithet.
#' Assumes the first word is the genus and the second is the species epithet.
#' Handles cases with only one word (genus only, species becomes NA).
#' Also handles cases where scientific_names might contain NAs or empty strings.
#'
#' @param scientific_names A character vector of scientific names.
#' @return A list with two character vectors: `genus` and `species`.
#' @noRd
.extract_genus_species <- function(scientific_names) {
  # Initialize vectors for genus and species with NAs
  genus_col <- rep(NA_character_, length(scientific_names))
  species_col <- rep(NA_character_, length(scientific_names))

  # Process only non-NA and non-empty names
  valid_indices <- !is.na(scientific_names) & nzchar(trimws(scientific_names))
  valid_names <- scientific_names[valid_indices]

  if (length(valid_names) > 0) {
    # Split the scientific name into parts based on the first space
    name_parts_list <- strsplit(valid_names, " ", fixed = TRUE)
    genus_col[valid_indices] <- sapply(name_parts_list, function(parts) if(length(parts) >= 1) parts[1] else NA_character_)
    species_col[valid_indices] <- sapply(name_parts_list, function(parts) if(length(parts) >= 2) parts[2] else NA_character_)
  }
  return(list(genus = genus_col, species = species_col))
}