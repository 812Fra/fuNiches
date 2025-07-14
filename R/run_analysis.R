#' Run Niche Analysis Workflow
#'
#' Orchestrates the complete niche analysis workflow based on prepared data
#' and user-defined settings.
#'
#' @param niche_data_object An object of class 'niche_data_object' (output from \code{\link{prepare_niche_data}}).
#' @param analysis_settings An object of class 'niche_analysis_settings' (output from \code{\link{set_analysis_settings}}).
#' @param output_dir Character. Path to the main directory where all outputs (CSVs, plots, reports) will be saved.
#'        Subdirectories for each target species will be created here.
#' @param verbose Logical. If TRUE, prints progress messages to the console. Default is TRUE.
#'
#' @return A list (invisibly) containing detailed results for each target species analyzed.
#'         The primary output is written to files in the specified `output_dir`.
#' @export
#' @importFrom dplyr %>% filter mutate select group_by summarise n distinct across all_of left_join count sym pull arrange desc slice_head bind_rows case_when where
#' @examples
#' \dontrun{
#' # --- 0. Define Paths and Basic Parameters ---
#' my_gpkg_path <- "path/to/your/occurrences.gpkg" # MODIFY
#' my_raster_dir <- "path/to/your/raster_layers/"  # MODIFY
#' my_output_dir <- "path/to/your/analysis_output/" # MODIFY
#'
#' # Ensure output directory exists
#' if (!dir.exists(my_output_dir)) {
#'   dir.create(my_output_dir, recurive = TRUE)
#' }
#'
#' # --- 1. Prepare Niche Data ---
#' # Assuming 'SoilType' and 'LandCover' are categorical raster layers in my_raster_dir
#' # and the GPKG has 'scientificName', 'Genus', 'Species' columns.
#' prepared_niche_data <- prepare_niche_data(
#'   occurrences_gpkg_path = my_gpkg_path,
#'   raster_dir_path = my_raster_dir,
#'   categorical_var_names = c("SoilType", "LandCover"), # Names of categorical rasters
#'   scientific_name_col = "scientificName",
#'   genus_col_name = "Genus",
#'   species_col_name = "Species",
#'   min_occurrences_per_species = 10,
#'   continuous_var_labels = list(bio1 = "Annual Mean Temperature (°C)"),
#'   categorical_var_labels = list(SoilType = c("1" = "Sandy", "2" = "Clay"))
#' )
#'
#' # Check if data preparation was successful
#' if (!is.null(prepared_niche_data) && nrow(prepared_niche_data$data) > 0) {
#'
#'   # --- 2. Set Analysis Settings ---
#'   # Example: Analyze Tuber aestivum and all species of Quercus genus
#'   # using Boruta, generating plots for the top 5 closest species.
#'   analysis_configuration <- set_analysis_settings(
#'     target_species = "Tuber_aestivum", # Example target species
#'     target_genera = "Quercus",         # Example target genus
#'     exclude_congenerics_if_species_target = TRUE,
#'     analyses_to_run = c("overlap_indices", "permanova_pairwise", "lda_multigroup", "dendrogram_global_hv", "enfa"),
#'     variable_selection_method = "Boruta",
#'     boruta_max_confirmed_vars = 5,
#'     generate_plots = TRUE,
#'     plot_type = "n_closest_species",
#'     plot_n_closest = 5,
#'     plot_n_closest_overlap_metric = "mean_schoener_D",
#'     enfa_run_for_closest_species = TRUE, # Set to TRUE to run ENFA on closest species
#'     enfa_n_closest = 3,                  # Number of closest species for ENFA
#'     enfa_plot_presence_points = TRUE,    # Whether to plot the cloud of presence points
#'     generate_report = TRUE,
#'     continuous_dist_plot_type = "ridgeline",
#'     plot_categorical_distributions = TRUE
#'   )
#'
#'   # --- 3. Run the Niche Analysis ---
#'   analysis_results <- run_niche_analysis(
#'     niche_data_object = prepared_niche_data,
#'     analysis_settings = analysis_configuration,
#'     output_dir = my_output_dir,
#'     verbose = TRUE
#'   )
#'
#'   # The 'analysis_results' object (returned invisibly) contains detailed outputs.
#'   # All primary outputs (CSVs, plots, reports) are saved in 'my_output_dir'.
#'
#' } else {
#'   message("Niche data preparation failed or resulted in no data. Analysis cannot proceed.")
#' }
#' }
#' @importFrom utils txtProgressBar setTxtProgressBar write.csv head read.csv
#' @importFrom knitr kable
#' @importFrom stats setNames na.omit cor dist prcomp lm median quantile fisher.test chisq.test p.adjust as.formula
#' @importFrom methods as
#' @importFrom igraph graph_from_data_frame E V gorder gsize edge_attr_names cluster_walktrap modularity membership layout_with_fr plot.igraph
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices png dev.off adjustcolor
#' @importFrom ade4 dudi.pca
#' @importFrom adehabitatHS enfa
#' @importFrom fs path_join dir_exists dir_create
#' @importFrom ggplot2 ggsave guides guide_legend
#' @importFrom rmarkdown render
#' @importFrom graphics par legend

run_niche_analysis <- function(niche_data_object,
                               analysis_settings,
                               output_dir,
                               verbose = TRUE) {

  # --- 1. Validazione Input Essenziali ---
  if (!inherits(niche_data_object, "niche_data")) { 
    stop("'niche_data_object' must be an object prepared by prepare_niche_data().")
  }
  if (!inherits(analysis_settings, "niche_analysis_settings")) {
    stop("'analysis_settings' must be an object configured by set_analysis_settings().")
  }
  if (!is.character(output_dir) || length(output_dir) != 1 || !nzchar(output_dir)) {
    stop("'output_dir' must be a valid character string path.")
  }
  if (!dir.exists(output_dir)) {
    if (verbose) message("Output directory '", output_dir, "' does not exist. Creating it.")
    tryCatch(dir.create(output_dir, recursive = TRUE),
             error = function(e) stop("Failed to create output directory: ", e$message))
  }

  # --- 2. Risoluzione Specie Target Effettive ---
actual_target_species_ids <- .resolve_target_species_list(
  settings = analysis_settings,
  available_species_ids = unique(niche_data_object$data[[niche_data_object$metadata$species_id_col]]),
  data_frame_for_genus_lookup = niche_data_object$data, 
  species_id_col_name = niche_data_object$metadata$species_id_col, # This is correct
  genus_col_name_processed = niche_data_object$metadata$genus_col_name_processed, 
  verbose = verbose
)
  if (is.null(actual_target_species_ids) || length(actual_target_species_ids) == 0) {
    stop("No valid target species could be identified from the settings and the prepared data.")
  }
  if (verbose) message("Targets for analysis: ", paste(actual_target_species_ids, collapse = ", "))
  # --- Create Global Species Color Palette ---
  all_species_in_data <- unique(niche_data_object$data[[niche_data_object$metadata$species_id_col]])
  global_species_color_palette <- NULL
  if (length(all_species_in_data) > 0 && requireNamespace("randomcoloR", quietly = TRUE)) {
    set.seed(42) 
    global_species_color_palette <- randomcoloR::distinctColorPalette(length(all_species_in_data), runTsne = FALSE) 
    names(global_species_color_palette) <- all_species_in_data
    if (verbose) message("Global species color palette created for ", length(all_species_in_data), " species.")
  }
  # --- 3. Initialization Overall Results List ---
    all_run_results <- list()
  # --- 4. Main Loop for Target Species ---
  if (verbose) message("\nStarting analysis loop for ", length(actual_target_species_ids), " target(s)...")

  for (i in seq_along(actual_target_species_ids)) {
    current_target_id <- actual_target_species_ids[i]
    if (verbose) message(paste0("\n\n===== PROCESSING TARGET (", i, "/", length(actual_target_species_ids), "): ", current_target_id, " ====="))
    target_specific_output_dir <- file.path(output_dir, gsub("[^A-Za-z0-9_.-]", "_", current_target_id))
    if (!dir.exists(target_specific_output_dir)) {
      dir.create(target_specific_output_dir, recursive = TRUE)
    }

    # 4.1 Data compilation for current target (target vs. comparison)
    target_comparison_data_list <- .prepare_target_comparison_data(
      current_target_id = current_target_id,
      niche_data_object = niche_data_object,
      settings = analysis_settings,
      verbose = verbose
    )
    if (is.null(target_comparison_data_list) || nrow(target_comparison_data_list$data_subset_for_analysis) == 0 || length(target_comparison_data_list$comparison_species_ids) == 0) {
      if (verbose) message("  Skipping ", current_target_id, ": insufficient data or no comparison species after filtering.")
      all_run_results[[current_target_id]] <- list(
        status = "skipped_insufficient_data_for_comparison",
        target_id = current_target_id,
                species_col_actual_name = niche_data_object$metadata$species_id_col,
        timestamp = Sys.time()
      )
      next
    }

    # 4.2 Variable Selection for the current target
    variable_selection_output <- .perform_variable_selection(
      data_for_selection = target_comparison_data_list$data_subset_for_analysis,
      target_id = current_target_id,
      all_potential_env_vars = niche_data_object$metadata$env_vars$all,
      categorical_vars = niche_data_object$metadata$env_vars$categorical,
      species_col_actual_name = niche_data_object$metadata$species_id_col, 
      settings = analysis_settings,
      metadata_labels = niche_data_object$metadata$labels,
      output_dir_target_specific = target_specific_output_dir,
      verbose = verbose
    )

    final_selected_variables <- variable_selection_output$selected_variables
    if (length(final_selected_variables) < analysis_settings$min_final_variables_for_multivariate) {
      if (verbose) message("  Skipping ", current_target_id, ": less than ",
                           analysis_settings$min_final_variables_for_multivariate,
                           " variables remaining after selection (", length(final_selected_variables), " found: ",
                           paste(final_selected_variables, collapse=", "),").")
      all_run_results[[current_target_id]] <- list(
        status = "skipped_insufficient_variables_after_selection",
        target_id = current_target_id,
        variables_initial_selection = variable_selection_output$initial_selection,
        variables_after_filters = final_selected_variables,
        timestamp = Sys.time()
      )
      next
    }
    if (verbose) message("  Final selected variables for ", current_target_id, ": ", paste(final_selected_variables, collapse = ", "))
    # Initializes the result list for the current target
    current_target_analysis_results <- list(
      metadata = list(
        target_species_id = current_target_id,
        run_timestamp = Sys.time(),
        analysis_settings_used = analysis_settings,
        variable_selection_details = variable_selection_output,
        final_variables_used = final_selected_variables,
        comparison_species_used = target_comparison_data_list$comparison_species_ids,
        niche_data_object_source = niche_data_object$metadata$parameters_used$occurrences_gpkg_path
      )
    )
    data_for_final_analyses <- target_comparison_data_list$data_subset_for_analysis %>%
      dplyr::select(
        dplyr::all_of(niche_data_object$metadata$species_id_col),
        dplyr::all_of(final_selected_variables)
      )
    # 4.3 Running the selected analysis modules
    if ("overlap_indices" %in% analysis_settings$analyses_to_run) {
      if (verbose) message("  Running: Overlap Indices module...")
      current_target_analysis_results$overlap_indices <- .run_overlap_module(
        data_for_module = data_for_final_analyses,
        target_id = current_target_id,
        comparison_ids = target_comparison_data_list$comparison_species_ids,
        selected_vars = final_selected_variables,
        species_id_col = niche_data_object$metadata$species_id_col, 
        categorical_vars_names = niche_data_object$metadata$env_vars$categorical, 
        output_dir_module = target_specific_output_dir,
        settings = analysis_settings,
        verbose = verbose
      )
    }
    if ("equivalency_test" %in% analysis_settings$analyses_to_run) {
      if (verbose) message("  Running: Niche Equivalency Test module...")
      current_target_analysis_results$equivalency_test <- .run_equivalency_module(
        data_for_module = data_for_final_analyses,
        target_id = current_target_id,
        comparison_ids = target_comparison_data_list$comparison_species_ids,
        selected_vars = final_selected_variables,
        species_id_col = niche_data_object$metadata$species_id_col, 
        categorical_vars_names = niche_data_object$metadata$env_vars$categorical, 
        output_dir_module = target_specific_output_dir,
        settings = analysis_settings,
        verbose = verbose
      )
    }
    if (!is.null(current_target_analysis_results$overlap_indices)) {
      .generate_top_n_summary_tables_per_analysis(
        target_id = current_target_id,
        all_results_for_target = current_target_analysis_results,
        settings = analysis_settings,
        output_dir = target_specific_output_dir,
        niche_data_object_metadata = niche_data_object$metadata,
        verbose = verbose
      )
    }
        # PERMANOVA Pairwise Module
    if ("permanova_pairwise" %in% analysis_settings$analyses_to_run) {
      if (verbose) message("  Running: Pairwise PERMANOVA module...")
      current_target_analysis_results$permanova_pairwise <- .run_permanova_pairwise_module(
        data_for_module = data_for_final_analyses,
        target_id = current_target_id,
        comparison_ids = target_comparison_data_list$comparison_species_ids,
        selected_vars = final_selected_variables,
        species_id_col = niche_data_object$metadata$species_id_col,
        categorical_vars_names = niche_data_object$metadata$env_vars$categorical,
        output_dir_module = target_specific_output_dir,
        settings = analysis_settings,
        verbose = verbose
      )
    }
    # LDA Module (Multigroup)
    if ("lda_multigroup" %in% analysis_settings$analyses_to_run) {
      if (verbose) message("  Running: Linear Discriminant Analysis (LDA) module...")
      current_target_analysis_results$lda_multigroup <- .run_lda_module(
        data_for_module = data_for_final_analyses,
        target_id = current_target_id, 
        comparison_ids = target_comparison_data_list$comparison_species_ids, 
        selected_vars = final_selected_variables,
        species_id_col = niche_data_object$metadata$species_id_col,
        categorical_vars_names = niche_data_object$metadata$env_vars$categorical,
        output_dir_module = target_specific_output_dir,
        settings = analysis_settings,
        verbose = verbose
      )
    }
    # Modularity Network Module
    if ("modularity_network" %in% analysis_settings$analyses_to_run) {
      if (verbose) message("  Running: Niche Modularity Network module...")
      if (!is.null(current_target_analysis_results$overlap_indices)) {
        current_target_analysis_results$modularity_network <- .run_modularity_module(
          data_for_module = data_for_final_analyses, 
          target_id = current_target_id,
          comparison_ids = target_comparison_data_list$comparison_species_ids,
          selected_vars = final_selected_variables, 
          species_id_col = niche_data_object$metadata$species_id_col,
          categorical_vars_names = niche_data_object$metadata$env_vars$categorical, 
          overlap_results = current_target_analysis_results$overlap_indices, 
          output_dir_module = target_specific_output_dir,
          settings = analysis_settings,
          verbose = verbose
        )
      } else {
        if (verbose) message("    Skipping Modularity Network: Overlap indices results not available.")
      }
    }
    # Hypervolume Module
    if ("hypervolume_analysis" %in% analysis_settings$analyses_to_run) {
      if (verbose) message("  Running: Niche Hypervolume Analysis module...")
      current_target_analysis_results$hypervolume_analysis <- .run_hypervolume_module(
        data_for_module = data_for_final_analyses,
        target_id = current_target_id,
        comparison_ids = target_comparison_data_list$comparison_species_ids,
        selected_vars = final_selected_variables,
        species_id_col = niche_data_object$metadata$species_id_col,
        categorical_vars_names = niche_data_object$metadata$env_vars$categorical,
        variable_selection_details = current_target_analysis_results$metadata$variable_selection_details, 
        output_dir_module = target_specific_output_dir,
        settings = analysis_settings,
        verbose = verbose
      )
    }
    # Variable Stability Module (if Boruta was used)
    if ("variable_stability" %in% analysis_settings$analyses_to_run &&
        analysis_settings$variable_selection_method == "Boruta") {
      if (verbose) message("  Running: Variable Stability Analysis module...")
      current_target_analysis_results$variable_stability <- .run_variable_stability_module(
        data_for_stability = target_comparison_data_list$data_subset_for_analysis, 
        target_id = current_target_id,
        potential_vars_for_stability = setdiff(niche_data_object$metadata$env_vars$all, niche_data_object$metadata$env_vars$categorical),
        species_col_actual_name = niche_data_object$metadata$species_id_col,
        settings = analysis_settings,
        output_dir_module = target_specific_output_dir,
        verbose = verbose
      )
    } else if ("variable_stability" %in% analysis_settings$analyses_to_run && analysis_settings$variable_selection_method != "Boruta") {
        if(verbose) message("    Skipping Variable Stability: Module selected, but variable selection method was not Boruta.")
    }
    # Similarity Validation Module
    if ("similarity_validation" %in% analysis_settings$analyses_to_run) {
      if (verbose) message("  Running: Similarity Validation module...")

      similar_species_details_list <- list()
      validation_sources <- analysis_settings$similarity_validation_source %||% c("permanova", "equivalency")

      if ("permanova" %in% validation_sources) {
        if (!is.null(current_target_analysis_results$permanova_pairwise) && nrow(current_target_analysis_results$permanova_pairwise) > 0 && "p_adj" %in% names(current_target_analysis_results$permanova_pairwise)) {
            p_thresh <- analysis_settings$similarity_validation_p_threshold %||% 0.05
            perm_similar_df <- current_target_analysis_results$permanova_pairwise %>%
              dplyr::filter(p_adj >= p_thresh) %>%
              dplyr::transmute(comparison_species = specie2, source = "PERMANOVA", metric = "p_adj", value = p_adj, threshold_used = paste(">=", p_thresh))
            if(nrow(perm_similar_df) > 0) similar_species_details_list[['permanova']] <- perm_similar_df
            if (verbose) message("    -> Found ", nrow(perm_similar_df), " similar species from PERMANOVA (p.adj >= ", p_thresh, ").")
        }
      }
      if ("equivalency" %in% validation_sources) {
        if (!is.null(current_target_analysis_results$equivalency_test) && nrow(current_target_analysis_results$equivalency_test) > 0 && "p_value" %in% names(current_target_analysis_results$equivalency_test)) {
            p_thresh <- analysis_settings$similarity_validation_p_threshold %||% 0.05
            equiv_similar_df <- current_target_analysis_results$equivalency_test %>%
              dplyr::filter(p_value >= p_thresh) %>%
              dplyr::transmute(comparison_species = comparison_species, source = "Equivalency Test", metric = "p_value", value = p_value, threshold_used = paste(">=", p_thresh))
            if(nrow(equiv_similar_df) > 0) similar_species_details_list[['equivalency']] <- equiv_similar_df
            if (verbose) message("    -> Found ", nrow(equiv_similar_df), " similar species from Equivalency Test (p >= ", p_thresh, ").")
        }
      }
      if ("overlap_indices" %in% validation_sources) {
        overlap_thresh <- analysis_settings$similarity_validation_overlap_threshold %||% 0.6
        overlap_metric <- analysis_settings$plot_n_closest_overlap_metric %||% "mean_schoener_D"
        if (!is.null(current_target_analysis_results$overlap_indices) && nrow(current_target_analysis_results$overlap_indices) > 0 && overlap_metric %in% names(current_target_analysis_results$overlap_indices)) {
            overlap_similar_df <- current_target_analysis_results$overlap_indices %>%
              dplyr::filter(!!dplyr::sym(overlap_metric) >= overlap_thresh) %>%
              dplyr::transmute(comparison_species = comparison_species, source = "Overlap Index", metric = overlap_metric, value = !!dplyr::sym(overlap_metric), threshold_used = paste(">=", overlap_thresh))
            if(nrow(overlap_similar_df) > 0) similar_species_details_list[['overlap_indices']] <- overlap_similar_df
            if (verbose) message("    -> Found ", nrow(overlap_similar_df), " similar species from Overlap Indices (", overlap_metric, " >= ", overlap_thresh, ").")
        } else if (verbose) {
            warning("    -> Could not use 'overlap_indices' source: results are missing or metric '", overlap_metric, "' not found.")
        }
      }
      if ("hypervolume" %in% validation_sources) {
        overlap_thresh <- analysis_settings$similarity_validation_overlap_threshold %||% 0.6
        hv_metric <- analysis_settings$hypervolume_plot_overlap_metric %||% "sorensen"
        if (!is.null(current_target_analysis_results$hypervolume_analysis$overlap_stats) && nrow(current_target_analysis_results$hypervolume_analysis$overlap_stats) > 0 && hv_metric %in% names(current_target_analysis_results$hypervolume_analysis$overlap_stats)) {
            hv_similar_df <- current_target_analysis_results$hypervolume_analysis$overlap_stats %>%
              dplyr::filter(!!dplyr::sym(hv_metric) >= overlap_thresh) %>%
              dplyr::transmute(comparison_species = comparison_species, source = "Hypervolume", metric = hv_metric, value = !!dplyr::sym(hv_metric), threshold_used = paste(">=", overlap_thresh))
            if(nrow(hv_similar_df) > 0) similar_species_details_list[['hypervolume']] <- hv_similar_df
            if (verbose) message("    -> Found ", nrow(hv_similar_df), " similar species from Hypervolume overlap (", hv_metric, " >= ", overlap_thresh, ").")
        } else if (verbose) {
            warning("    -> Could not use 'hypervolume' source: results are missing or metric '", hv_metric, "' not found.")
        }
      }

      similar_species_details_df <- dplyr::bind_rows(similar_species_details_list)
      similar_species_for_validation <- unique(similar_species_details_df$comparison_species)

      if (length(similar_species_for_validation) == 0 && verbose) {
        message("    -> No species met the specified similarity criteria. Similarity Validation module will be skipped.")
      }

      # The helper function will now handle all cases, including when no similar species are found.
      current_target_analysis_results$similarity_validation <- .run_similarity_validation_module(
        data_for_validation = data_for_final_analyses,
        target_id = current_target_id,
        similar_species_ids = similar_species_for_validation,
        similar_species_details = similar_species_details_df,
        selected_vars = final_selected_variables,
        species_id_col = niche_data_object$metadata$species_id_col,
        categorical_vars_names = niche_data_object$metadata$env_vars$categorical,
        output_dir_module = target_specific_output_dir,
        settings = analysis_settings,
        verbose = verbose
      )
    }
    # Categorical Association Module (if selected and categorical vars exist)
    if ("categorical_association" %in% analysis_settings$analyses_to_run) {
      if (verbose) message("  Running: Categorical Association module...")      
      final_categorical_vars_for_tests <- intersect(final_selected_variables, niche_data_object$metadata$env_vars$categorical)
            if (length(final_categorical_vars_for_tests) == 0) {
        if (verbose) message("    Skipping Categorical Association: No categorical variables selected or available in final data.")
      } else {
      current_target_analysis_results$categorical_association <- .run_categorical_association_module(
        data_for_env_tests = data_for_final_analyses,  
        target_id = current_target_id,
        species_id_col = niche_data_object$metadata$species_id_col,
        selected_categorical_vars = final_categorical_vars_for_tests,
        output_dir_module = target_specific_output_dir,
        settings = analysis_settings,
        verbose = verbose
      )
    }}
    # ENFA Module
    if ("enfa" %in% analysis_settings$analyses_to_run) {
      if (verbose) message("  Running: Ecological Niche Factor Analysis (ENFA) module...")
      enfa_results_for_target <- .run_enfa_module(
        data_for_module = data_for_final_analyses,
        target_id = current_target_id,
        selected_vars = final_selected_variables,
        categorical_vars_names = niche_data_object$metadata$env_vars$categorical,
        species_id_col = niche_data_object$metadata$species_id_col,
        metadata_labels = niche_data_object$metadata$labels,
        settings = analysis_settings,
        output_dir_module = target_specific_output_dir, # Pass for potential internal use, though plots saved externally
        verbose = verbose
      )
      if (!is.null(enfa_results_for_target)) {
        current_target_analysis_results$enfa <- enfa_results_for_target

        # Save ENFA plot directly here to ensure it's always handled,
        # regardless of the structure of .generate_all_plots_for_target.
        if (analysis_settings$generate_plots && !is.null(enfa_results_for_target$plot)) {
          enfa_plot_dir <- fs::path_join(c(target_specific_output_dir, "ENFA_Results"))
          if (!fs::dir_exists(enfa_plot_dir)) {
            fs::dir_create(enfa_plot_dir, recurse = TRUE)
          }
          plot_filename <- fs::path_join(c(enfa_plot_dir, paste0("ENFA_", current_target_id, ".png")))
          if (verbose) message("    Saving ENFA plot for ", current_target_id, " to: ", plot_filename)
          ggplot2::ggsave(plot_filename, plot = enfa_results_for_target$plot, width = 8, height = 6, units = "in", dpi = 300)
        }
      }
      # --- Run ENFA for N closest species (if requested) ---
      n_closest_for_enfa <- analysis_settings$enfa_n_closest %||% 3 # Default to 3 if not specified
      if ((analysis_settings$enfa_run_for_closest_species %||% TRUE) && n_closest_for_enfa > 0) {
        if (verbose) message("  -> Running ENFA for the ", n_closest_for_enfa, " closest species to ", current_target_id, "...")
        if (is.null(current_target_analysis_results$overlap_indices)) {
          if (verbose) warning("    Cannot run ENFA for closest species: overlap indices not calculated or available.")
        } else {
          closest_species_ids <- .get_n_closest_species(
            target_id = current_target_id,
            overlap_results = current_target_analysis_results$overlap_indices,
            n_closest = n_closest_for_enfa,
            overlap_metric = analysis_settings$plot_n_closest_overlap_metric %||% "composite_overlap",
            verbose = verbose
          )
          if (length(closest_species_ids) > 0) {
            current_target_analysis_results$enfa_closest_species <- list()
            for (closest_sp_id in closest_species_ids) {
              if (verbose) message("    --> Running ENFA for closest species: ", closest_sp_id)
              enfa_results_for_closest <- .run_enfa_module(
                data_for_module = data_for_final_analyses,
                target_id = closest_sp_id, # The "target" for this run is the closest species
                selected_vars = final_selected_variables, # Use the same variables as the main target
                categorical_vars_names = niche_data_object$metadata$env_vars$categorical,
                species_id_col = niche_data_object$metadata$species_id_col,
                metadata_labels = niche_data_object$metadata$labels,
                settings = analysis_settings,
                output_dir_module = target_specific_output_dir,
                verbose = verbose
              )
              if (!is.null(enfa_results_for_closest)) {
                current_target_analysis_results$enfa_closest_species[[closest_sp_id]] <- enfa_results_for_closest
                if (analysis_settings$generate_plots && !is.null(enfa_results_for_closest$plot)) {
                  enfa_plot_dir <- fs::path_join(c(target_specific_output_dir, "ENFA_Results"))
                  if (!fs::dir_exists(enfa_plot_dir)) {
                    fs::dir_create(enfa_plot_dir, recurse = TRUE)
                  }
                  # Create a descriptive filename
                  plot_filename_closest <- fs::path_join(c(enfa_plot_dir, paste0("ENFA_closest_to_", current_target_id, "_FOR_", closest_sp_id, ".png")))
                  if (verbose) message("      Saving ENFA plot for ", closest_sp_id, " to: ", plot_filename_closest)
                  ggplot2::ggsave(plot_filename_closest, plot = enfa_results_for_closest$plot, width = 8, height = 6, units = "in", dpi = 300)
                }
              }
            }
          } else {
            if (verbose) message("    No closest species found to run additional ENFA.")
          }
        }
      }
      # --- Generate Comparison Plot and Summary Table ---
      # Check if there are results for both the target and at least one closest species
      if (!is.null(current_target_analysis_results$enfa) && 
          !is.null(current_target_analysis_results$enfa_closest_species) && 
          length(current_target_analysis_results$enfa_closest_species) > 0) {
        
        if(verbose) message("  -> Generating ENFA comparison plot and summary table...")
        
        # Generate and save the combined plot
        .generate_enfa_comparison_plot(current_target_analysis_results, target_specific_output_dir, global_species_color_palette, verbose)
        
        # Generate and save the summary table
        .generate_enfa_summary_table(current_target_analysis_results, target_specific_output_dir, verbose)
        
      }
    }

    # --- 4.4 Plot and Report Generation ---
    # Orchestrate all plot generation for the current target
    if (analysis_settings$generate_plots) {
      # This function is assumed to exist elsewhere and handles all plotting logic
      .generate_all_plots_for_target(
        results_for_target = current_target_analysis_results,
        data_for_plots = data_for_final_analyses,
        niche_data_object_metadata = niche_data_object$metadata,
        global_species_palette = global_species_color_palette,
        settings = analysis_settings,
        output_dir = target_specific_output_dir,
        verbose = verbose
      )
    }

    # Generate the final report for the current target
    if (analysis_settings$generate_report) {
      if (verbose) message("  Generating final report for ", current_target_id, "...")
      # This function is assumed to exist elsewhere and handles report compilation
      .generate_report_for_target(
        results_for_target = current_target_analysis_results,
        settings = analysis_settings,
        output_dir = target_specific_output_dir,
        verbose = verbose
      )
    }

    # --- 4.5 Save and Store Results ---
    # Save the complete results object for the current target as an .rds file
    results_rds_path <- file.path(target_specific_output_dir, paste0("analysis_results_", current_target_id, ".rds"))
    if (verbose) message("  Saving detailed results object to: ", results_rds_path)
    saveRDS(current_target_analysis_results, file = results_rds_path)

    # Add the results for the current target to the main list that will be returned
    all_run_results[[current_target_id]] <- current_target_analysis_results
  }

  # --- Generate Unified Plot for All Species in the Dataset (if requested) ---
  if (analysis_settings$generate_plots && (analysis_settings$plot_unified_all_species_distributions %||% FALSE)) {
    # Start with all species
    species_to_include_in_unified_plot <- all_species_in_data

    # Apply congeneric exclusion if requested, respecting the 'if_species_target' condition
    if (analysis_settings$exclude_congenerics_if_species_target) {
      if (verbose) message("\nApplying congeneric exclusion for the unified plot...")

      # Identify which of the actual targets were specified directly as species
      targets_specified_as_species <- intersect(actual_target_species_ids, analysis_settings$target_species_input %||% character(0))

      if (length(targets_specified_as_species) > 0) {
        # Find the genera of these specific targets
        target_genera_for_exclusion <- niche_data_object$data %>%
          dplyr::filter(!!dplyr::sym(niche_data_object$metadata$species_id_col) %in% targets_specified_as_species) %>%
          dplyr::pull(!!dplyr::sym(niche_data_object$metadata$genus_col_name_processed)) %>%
          unique() %>%
          stats::na.omit()

        if (length(target_genera_for_exclusion) > 0) {
          if (verbose) message("  Genera from user-specified species targets identified for exclusion: ", paste(target_genera_for_exclusion, collapse = ", "))

          # Identify all comparison species that are congenerics of these specific targets
          congenerics_to_exclude <- niche_data_object$data %>%
            dplyr::filter(
              !!dplyr::sym(niche_data_object$metadata$genus_col_name_processed) %in% target_genera_for_exclusion &
                !(!!dplyr::sym(niche_data_object$metadata$species_id_col) %in% actual_target_species_ids) # Still exclude from all comparison species
            ) %>%
            dplyr::pull(!!dplyr::sym(niche_data_object$metadata$species_id_col)) %>%
            unique()

          if (length(congenerics_to_exclude) > 0) {
            if (verbose) message("  Excluding ", length(congenerics_to_exclude), " congeneric comparison species from unified plot: ", paste(congenerics_to_exclude, collapse = ", "))
            species_to_include_in_unified_plot <- setdiff(species_to_include_in_unified_plot, congenerics_to_exclude)
          } else {
            if (verbose) message("  No congeneric comparison species found to exclude.")
          }
        } else {
          if (verbose) message("  Could not determine genera for the user-specified species targets. Skipping congeneric exclusion.")
        }
      } else {
        if (verbose) message("  No targets were specified directly as species (all came from genera). Skipping congeneric exclusion for unified plot.")
      }
    }
    tryCatch({
      .plot_unified_all_species_distributions(
        data_for_plot = niche_data_object$data,
        species_to_include = species_to_include_in_unified_plot, # Use the potentially filtered list
        categorical_vars_to_plot = niche_data_object$metadata$env_vars$categorical,
        species_id_col = niche_data_object$metadata$species_id_col,
        niche_data_metadata = niche_data_object$metadata,
        global_palette = global_species_color_palette,
        output_dir = output_dir, # Save in the main output directory
        verbose = verbose
      )
    }, error = function(e) {
      if (verbose) warning("Unified all-species distribution plot failed. Error: ", e$message)
    })
  }

  # --- Generate Global Hypervolume Dendrogram (if requested) ---
  if ("dendrogram_global_hv" %in% analysis_settings$analyses_to_run) {
    if (verbose) message("\nGenerating Global Hypervolume Dendrogram for all species...")
    if (is.null(all_run_results[["global_analysis_summaries"]])) {
      all_run_results$global_analysis_summaries <- list()
    }
    all_run_results$global_analysis_summaries$dendrogram_global_hv <- tryCatch({
      .run_dendrogram_global_hv_module(
        niche_data_object = niche_data_object,
        settings = analysis_settings,
        output_dir_base = output_dir,
        global_species_palette = global_species_color_palette,
        verbose = verbose
      )
    }, error = function(e) {
      if (verbose) warning("  Global Hypervolume Dendrogram module failed: ", e$message)
      list(status = "error", message = e$message, timestamp = Sys.time())
    })
    saveRDS(all_run_results, file = file.path(output_dir, "all_targets_analysis_results_with_global_dendro.rds"))
  }

  # --- Finalization ---
  if (verbose) message("\n\n===== Analysis Workflow Complete =====")

  # Save the entire results object to a single file for later use.
  # This is a good practice and provides a complete record before the function returns.
  final_results_path <- file.path(output_dir, "all_targets_analysis_results.rds")
  if (verbose) message("  Saving complete results object for all targets to: ", final_results_path)
  tryCatch({
    saveRDS(all_run_results, file = final_results_path)
    if (verbose) message("  Successfully saved complete results object.")
  }, error = function(e) {
    if (verbose) warning("  Failed to save the complete results object. Error: ", e$message)
  })

  if (verbose) message("  Function finished. Returning results invisibly.")

  return(invisible(all_run_results))
}

#' @title .run_enfa_module
#' @description Helper function to perform Ecological Niche Factor Analysis (ENFA) for a single species.
#' @param data_for_module A data frame containing species occurrence and environmental data for the current target and comparison species.
#' @param target_id Character. The ID of the current target species.
#' @param selected_vars Character vector. Names of the selected environmental variables to use for ENFA.
#' @param species_id_col Character. The name of the column containing species IDs in `data_for_module`.
#' @param settings A list of analysis settings, typically generated by `set_analysis_settings`.
#' @param output_dir_module Character. The target-specific output directory (used for messages, not direct saving here).
#' @param verbose Logical. If TRUE, prints progress messages.
#' @return A list containing ENFA results for the target species, including the `enfa_object`,
#'         marginality, specialization, eigenvalues, species scores, environmental loadings, and a ggplot object.
#'         Returns NULL if ENFA cannot be performed or fails.
#' @importFrom ade4 dudi.pca
#' @importFrom adehabitatHS enfa
#' @importFrom ggplot2 ggplot aes geom_point geom_segment geom_text labs theme_minimal arrow unit
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr %>% select all_of mutate
#' @keywords internal
.run_enfa_module <- function(data_for_module, target_id, selected_vars, categorical_vars_names, species_id_col, metadata_labels, settings, output_dir_module, verbose) {
  if (verbose) message("    Performing ENFA for: ", target_id)

  # 1. Prepare data: ENFA works only on numeric variables
  numeric_vars_for_enfa <- setdiff(selected_vars, categorical_vars_names)
  if (verbose) message("      ENFA DEBUG: Found ", length(numeric_vars_for_enfa), " numeric variables: ", paste(numeric_vars_for_enfa, collapse=", "))
  if (length(numeric_vars_for_enfa) < 2) {
    if (verbose) message("      ENFA SKIPPED: Not enough numeric variables (need at least 2).")
    return(NULL)
  }

  # 2. Prepare a clean data frame with only the necessary columns and handle NAs
  data_for_enfa_clean <- data_for_module %>%
    dplyr::select(dplyr::all_of(c(species_id_col, numeric_vars_for_enfa))) %>%
    stats::na.omit() %>%
    as.data.frame() # Coerce to base data.frame for maximum compatibility with older packages

  if (verbose) message("      ENFA DEBUG: Rows in data_for_module: ", nrow(data_for_module), ". Rows after na.omit: ", nrow(data_for_enfa_clean), ".")
  if (nrow(data_for_enfa_clean) == 0) {
    if (verbose) message("      ENFA SKIPPED: No complete observations found for the selected numeric variables.")
    return(NULL)
  }

  # 3. Create presence vector from the CLEANED data
  # Coerce to numeric (0/1) for maximum compatibility with older packages like adehabitatHS
  presence_vector <- as.numeric(data_for_enfa_clean[[species_id_col]] == target_id)
  if (verbose) message("      ENFA DEBUG: Number of presence points for '", target_id, "' in clean data: ", sum(presence_vector))

  # Add a more robust check: we need both presences and absences for the analysis to be meaningful.
  if (length(unique(presence_vector)) < 2) {
    if (verbose) message("      ENFA SKIPPED: The presence vector for '", target_id, "' does not contain both presences and absences in the clean data subset. It's either all presences or all absences.")
    return(NULL)
  }
  if (sum(presence_vector) < 5) { # ENFA can be unstable with very few points
    if (verbose) message("      ENFA SKIPPED: Too few presence points (< 5) for a stable analysis for '", target_id, "'.")
    return(NULL)
  }

  # 4. Get the final environmental data (now guaranteed to be clean) and check for constants
  env_data_filtered <- data_for_enfa_clean %>%
    dplyr::select(dplyr::all_of(numeric_vars_for_enfa))

  constant_cols <- sapply(env_data_filtered, function(x) length(unique(x)) < 2)
  if (any(constant_cols)) {
    if (verbose) message("      ENFA DEBUG: Constant columns found and removed: ",
                                paste(names(env_data_filtered)[constant_cols], collapse = ", "))
    env_data_filtered <- env_data_filtered[, !constant_cols, drop = FALSE]
  }
  if (verbose) message("      ENFA DEBUG: Number of varying numeric variables after constant check: ", ncol(env_data_filtered))
  if (ncol(env_data_filtered) < 2) {
    if (verbose) message("      ENFA SKIPPED: Not enough varying environmental variables remaining.")
    return(NULL)
  }

  # 5. Perform ENFA using the standard and most stable two-step method: PCA then ENFA.
  # The error "object of class dudi expected" is a strong indication that this is the required workflow.
  pca_env <- tryCatch({
    ade4::dudi.pca(env_data_filtered, scannf = FALSE, nf = min(5, ncol(env_data_filtered)))
  }, error = function(e) {
    # Use message() for more immediate and clear feedback in the console
    message(paste0("\n\n*** ENFA ERROR in PCA step for '", target_id, "' ***"))
    message("    Error message: ", e$message)
    message("    This is often due to issues in the input data 'env_data_filtered'. Skipping ENFA module.\n")
    return(NULL)
  })
  if (is.null(pca_env)) {
    return(NULL)
  }
  
  # --- CRITICAL SANITY CHECK for the PCA object ---
  # This is the most likely point of failure. Check for NaNs in the PCA output.
  if (verbose) message("      ENFA DEBUG: Performing sanity check on PCA output...")
  # Use unlist() to safely check for NaNs inside lists, data.frames, or vectors, which avoids the "method not implemented for type 'list'" error.
  if (any(is.nan(unlist(pca_env$li))) || any(is.nan(unlist(pca_env$co))) || any(is.nan(unlist(pca_env$eig)))) {
    message(paste0("\n\n*** ENFA CRITICAL FAILURE for '", target_id, "' ***"))
    message("    The PCA produced NaN (Not a Number) values. This is almost always due to")
    message("    high multicollinearity (perfectly correlated variables) in the input data.")
    message("    Check the correlation matrix for the following variables:")
    message("    ", paste(colnames(env_data_filtered), collapse = ", "))
    message("    Is NaN in site scores ($li)? ", any(is.nan(unlist(pca_env$li))))
    message("    Is NaN in var loadings ($co)? ", any(is.nan(unlist(pca_env$co))))
    message("    Is NaN in eigenvalues ($eig)? ", any(is.nan(unlist(pca_env$eig))))
    message("    Skipping ENFA module.\n")
    return(NULL)
  }

  enfa_obj <- tryCatch({
    # Pass the dudi object from the PCA to enfa.
    # Determine the number of axes to request, ensuring it's not more than available.
    nf_requested <- if (!is.null(settings$enfa_n_axes)) settings$enfa_n_axes else 2
    nf_available <- pca_env$nf
    nf_to_use <- min(nf_requested, nf_available)
    if (nf_requested > nf_available && verbose) {
      warning(paste0("Requested ", nf_requested, " ENFA axes, but only ", nf_available, " are available from PCA. Using ", nf_to_use, "."))
    }
    
    adehabitatHS::enfa(pca_env, presence_vector, scannf = FALSE, nf = nf_to_use)
  }, error = function(e) {
    # Use message() for more immediate and clear feedback in the console
    message(paste0("\n\n*** ENFA ERROR in ENFA analysis step for '", target_id, "' ***"))
    message("    Error message: ", e$message)
    message("    This can happen if the PCA object or presence vector is not as expected. Skipping ENFA module.\n")
    return(NULL)
  })

  if (is.null(enfa_obj)) {
    return(NULL)
  }

  # 6. DEBUG and Validate the structure of the enfa_obj before proceeding
  if (verbose) {
    message("      ENFA DEBUG: Structure of the returned enfa object:")
    utils::str(enfa_obj, max.level = 1)
    if (!is.null(enfa_obj$mar)) message("        - Marginality per variable (mar): ", paste(round(enfa_obj$mar, 2), collapse = ", "))
    if (!is.null(enfa_obj$s)) message("        - Eigenvalues (s): ", paste(round(enfa_obj$s, 2), collapse = ", "))
    if (!is.null(enfa_obj$li)) message("        - Site scores (li) dimensions: ", paste(dim(enfa_obj$li), collapse = "x"))
    if (!is.null(enfa_obj$co)) message("        - Env loadings (co) dimensions: ", paste(dim(enfa_obj$co), collapse = "x"))
  }

  # More robust validation checks
  # In modern adehabitatHS, $m is the global marginality (a single value), while $mar is a vector of marginalities per variable.
  # We check $m for the summary, and $mar for its structure.
  if (is.null(enfa_obj$m) || !is.numeric(enfa_obj$m) || length(enfa_obj$m) != 1 || !is.finite(enfa_obj$m)) {
    if (verbose) message("      ENFA SKIPPED: Produced a non-numeric, non-finite, or invalid global marginality value ($m). This often indicates an issue with the presence/absence data distribution relative to the environmental variables.")
    return(NULL)
  }
  # The eigenvalues are in the '$s' component.
  if (is.null(enfa_obj$s) || !is.numeric(enfa_obj$s) || !is.vector(enfa_obj$s) || any(!is.finite(enfa_obj$s))) {
    if (verbose) message("      ENFA SKIPPED: Produced non-numeric or non-finite eigenvalues ($s). This can happen if specialization axes cannot be computed.")
    return(NULL)
  }
  # Check that scores/loadings are either a data.frame or matrix with numeric content
  if (is.null(enfa_obj$li) || !(is.data.frame(enfa_obj$li) || is.matrix(enfa_obj$li)) || !all(sapply(enfa_obj$li, is.numeric))) {
    if (verbose) message("      ENFA SKIPPED: Produced invalid or non-numeric site scores (li). Expected a data.frame or matrix with numeric columns.")
    return(NULL)
  }
  if (is.null(enfa_obj$co) || !(is.data.frame(enfa_obj$co) || is.matrix(enfa_obj$co)) || !all(sapply(enfa_obj$co, is.numeric))) {
    if (verbose) message("      ENFA SKIPPED: Produced invalid or non-numeric environmental loadings (co). Expected a data.frame or matrix with numeric columns.")
    return(NULL)
  }

  # 7. Extract results (now safe to do so)
  marginality <- enfa_obj$m # Use $m for global marginality
  # The eigenvalues are in $s. The first corresponds to the marginality axis.
  # Specialization is calculated on the remaining axes.
  eigenvalues <- enfa_obj$s
  specialization_eigenvalues <- if (length(eigenvalues) > 1) eigenvalues[-1] else numeric(0)
  specialization <- if (length(specialization_eigenvalues) > 0) sqrt(mean(specialization_eigenvalues)) else 0

  # Get scores only for the presence points for plotting
  presence_scores <- enfa_obj$li[presence_vector == 1, , drop = FALSE]
  env_loadings <- enfa_obj$co # Loadings of environmental variables on ENFA axes

  # 8. Generate Plot
  n_axes_to_plot <- min(ncol(presence_scores), ncol(env_loadings))
  enfa_plot <- NULL
  plot_data_species <- NULL # Initialize to ensure it's always returned
  if (n_axes_to_plot < 2) {
    if (verbose) warning(paste0("Not enough ENFA axes (", n_axes_to_plot, ") to create a 2D plot for species ", target_id, ". Skipping plot generation."))
  } else {
    plot_data_species <- as.data.frame(presence_scores[, 1:2, drop = FALSE])
    colnames(plot_data_species) <- paste0("Axis", 1:2)
    plot_data_env <- as.data.frame(env_loadings[, 1:2])
    colnames(plot_data_env) <- paste0("Axis", 1:2)
    
    # Use pretty labels for variables if available
    var_names <- rownames(env_loadings)
    # The `rlang::`%||%` operator is a safe way to provide a default value if the label is NULL.
    # We need to ensure rlang is available or use a base R alternative. Let's use base R for robustness.
    get_label <- function(var_name) {
      label <- if (!is.null(metadata_labels$continuous_vars)) metadata_labels$continuous_vars[[var_name]] else NULL
      if (is.null(label) || !is.character(label) || length(label) != 1) var_name else label
    }
    if (!is.null(metadata_labels)) {
      pretty_var_names <- sapply(var_names, get_label, USE.NAMES = FALSE)
      plot_data_env$Variable <- pretty_var_names
    } else {
      plot_data_env$Variable <- var_names
    }

    explained_var <- round(eigenvalues / sum(eigenvalues) * 100, 2)
    # The first axis in ENFA is always the marginality axis.
    # The second is the first specialization axis.
    xlab_label <- paste0("Marginality Axis (", explained_var[1], "%)")
    ylab_label <- paste0("First Specialization Axis (", explained_var[2], "%)")

    subtitle_text <- paste0("Marginality: ", round(marginality, 2), ", Specialization: ", round(specialization, 2))

    enfa_plot <- ggplot2::ggplot() +
      ggplot2::labs(title = paste0("ENFA for ", target_id),
                    subtitle = subtitle_text,
                    x = xlab_label,
                    y = ylab_label) +
      ggplot2::theme_minimal()

    if (settings$enfa_plot_type == "biplot" && settings$enfa_plot_env_vectors) {
      # Scale vectors for better visualization
      scaling_factor <- max(abs(range(plot_data_species[, 1:2], na.rm = TRUE))) / max(abs(range(plot_data_env[, 1:2], na.rm = TRUE))) * 0.8
      plot_data_env <- plot_data_env %>% dplyr::mutate(Axis1_scaled = Axis1 * scaling_factor, Axis2_scaled = Axis2 * scaling_factor)
      enfa_plot <- enfa_plot +
        ggplot2::geom_segment(data = plot_data_env,
                              ggplot2::aes(x = 0, y = 0, xend = Axis1_scaled, yend = Axis2_scaled),
                              arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")),
                              color = "darkblue", alpha = 0.7)
      if (settings$enfa_plot_env_labels) {
        enfa_plot <- enfa_plot +
          ggplot2::geom_text(data = plot_data_env,
                             ggplot2::aes(x = Axis1_scaled, y = Axis2_scaled, label = Variable),
                             vjust = -0.5, hjust = 0.5, color = "darkblue", size = 3)
      }
    }
    # Optionally plot the cloud of individual presence points
    if (settings$enfa_plot_presence_points %||% TRUE) {
      enfa_plot <- enfa_plot +
        ggplot2::geom_point(data = plot_data_species,
                            ggplot2::aes(x = Axis1, y = Axis2),
                            color = "darkseagreen", alpha = 0.6, size = 1.5)
    }
    
    # Plot the species label and a point at the centroid of the presence points
    if (settings$enfa_plot_species_labels) {
      # Label only the centroid of the presence points for clarity
      if (nrow(plot_data_species) > 0) {
        centroid_data <- data.frame(Axis1 = mean(plot_data_species$Axis1, na.rm = TRUE), Axis2 = mean(plot_data_species$Axis2, na.rm = TRUE), label = target_id)
        
        # Add a distinct point for the centroid
        enfa_plot <- enfa_plot +
          ggplot2::geom_point(data = centroid_data,
                              ggplot2::aes(x = Axis1, y = Axis2),
                              color = "darkgreen", size = 4, shape = 18) # shape=18 is a diamond
        
        # Add the text label for the centroid, ensuring it doesn't overlap the point
        enfa_plot <- enfa_plot +
          ggrepel::geom_text_repel(data = centroid_data,
                                   ggplot2::aes(x = Axis1, y = Axis2, label = label), fontface = "bold", color = "darkgreen", size = 4,
                                   box.padding = 0.5, point.padding = 0.5, min.segment.length = 0)
      }
    }
  }

  # 9. Save numerical results
  enfa_results_dir <- fs::path_join(c(output_dir_module, "ENFA_Results"))
  if (!fs::dir_exists(enfa_results_dir)) {
    fs::dir_create(enfa_results_dir, recurse = TRUE)
  }
  # Robustly create the summary data frame to avoid length mismatch errors
  summary_list <- list(
    list(Metric = "Marginality", Value = marginality),
    list(Metric = "Specialization", Value = specialization)
  )
  if (length(eigenvalues) > 0) {
    eigen_list <- lapply(seq_along(eigenvalues), function(i) {
      list(Metric = paste0("Eigenvalue_Axis_", i), Value = eigenvalues[i])
    })
    summary_list <- c(summary_list, eigen_list)
  }
  summary_df <- dplyr::bind_rows(summary_list)

  summary_filename <- fs::path_join(c(enfa_results_dir, paste0("ENFA_summary_", target_id, ".csv")))
  utils::write.csv(summary_df, summary_filename, row.names = FALSE)
  utils::write.csv(as.data.frame(enfa_obj$li), fs::path_join(c(enfa_results_dir, paste0("ENFA_site_scores_", target_id, ".csv"))), row.names = TRUE)
  utils::write.csv(as.data.frame(env_loadings), fs::path_join(c(enfa_results_dir, paste0("ENFA_env_loadings_", target_id, ".csv"))), row.names = TRUE)
  if (verbose) message("      ENFA numerical results saved to: ", enfa_results_dir)

  # Remove the 'call' element, which can cause issues with saving/storing the object
  enfa_obj$call <- NULL

  return(list(
    enfa_object = enfa_obj,
    marginality = marginality,
    specialization = specialization,
    eigenvalues = eigenvalues,
    site_scores = enfa_obj$li, # Note: these are scores for ALL sites, not just presence
    env_loadings = env_loadings,
    plot = enfa_plot,
    presence_scores_for_plot = plot_data_species
  ))
}

#' @title .generate_enfa_comparison_plot
#' @description Generates a single ENFA biplot comparing the target species with its closest neighbors.
#' @param results A list of ENFA results for the target and closest species.
#' @param output_dir Character. The directory to save the plot.
#' @param palette A named color palette for species.
#' @param verbose Logical. If TRUE, prints messages.
#' @keywords internal
.generate_enfa_comparison_plot <- function(results, output_dir, palette, verbose) {
  target_enfa <- results$enfa
  target_id <- results$metadata$target_species_id
  
  # Start with the base plot of the target species (which includes the environmental vectors)
  comparison_plot <- target_enfa$plot
  if (is.null(comparison_plot)) {
    if (verbose) warning("Cannot generate ENFA comparison plot: base plot for target is missing.")
    return(NULL)
  }
  
  # --- 1. Clean the base plot: remove all species-specific points and labels ---
  # We keep only the axes and the environmental vectors.
  layers_to_remove_indices <- c()
  for (i in seq_along(comparison_plot$layers)) {
    layer <- comparison_plot$layers[[i]]
    # Remove any layer that is a point, or a text label.
    if (inherits(layer$geom, "GeomPoint") || inherits(layer$geom, "GeomTextRepel")) {
      layers_to_remove_indices <- c(layers_to_remove_indices, i)
    }
  }
  if (length(layers_to_remove_indices) > 0) {
    comparison_plot$layers[layers_to_remove_indices] <- NULL
  }

  comparison_plot <- comparison_plot + ggplot2::labs(title = paste("ENFA Comparison for", target_id, "and Neighbors"))

  # --- 2. Prepare data for all species ---
  all_points_data <- list()
  all_centroids_data <- list()
  
  # Add target species data first
  target_presence_data <- target_enfa$presence_scores_for_plot
  if (!is.null(target_presence_data) && nrow(target_presence_data) > 0) {
    target_presence_data$Species <- target_id
    all_points_data[[target_id]] <- target_presence_data
    
    all_centroids_data[[target_id]] <- data.frame(
      Axis1 = mean(target_presence_data$Axis1, na.rm = TRUE),
      Axis2 = mean(target_presence_data$Axis2, na.rm = TRUE),
      Species = target_id
    )
  }

  # Add neighbor species data
  for (sp_id in names(results$enfa_closest_species)) {
    sp_results <- results$enfa_closest_species[[sp_id]]
    
    # Project neighbor presence onto target's ENFA space
    presence_vector_neighbor <- as.numeric(sp_results$enfa_object$pr)
    presence_scores_neighbor <- target_enfa$enfa_object$li[presence_vector_neighbor == 1, , drop = FALSE]
    
    if (nrow(presence_scores_neighbor) > 0) {
      plot_data_neighbor <- as.data.frame(presence_scores_neighbor[, 1:2, drop = FALSE])
      colnames(plot_data_neighbor) <- c("Axis1", "Axis2")
      plot_data_neighbor$Species <- sp_id
      all_points_data[[sp_id]] <- plot_data_neighbor
      
      all_centroids_data[[sp_id]] <- data.frame(
        Axis1 = mean(plot_data_neighbor$Axis1, na.rm = TRUE),
        Axis2 = mean(plot_data_neighbor$Axis2, na.rm = TRUE),
        Species = sp_id
      )
    }
  }
  
  # Combine into single data frames
  final_points_df <- dplyr::bind_rows(all_points_data)
  final_centroids_df <- dplyr::bind_rows(all_centroids_data)

  # --- 3. Draw everything from scratch ---
  if (nrow(final_points_df) > 0) {
    comparison_plot <- comparison_plot +
      ggplot2::geom_point(data = final_points_df, ggplot2::aes(x = Axis1, y = Axis2, color = Species), alpha = 0.5, size = 1.5)
  }
  
  if (nrow(final_centroids_df) > 0) {
    comparison_plot <- comparison_plot +
      ggplot2::geom_point(data = final_centroids_df, ggplot2::aes(x = Axis1, y = Axis2, color = Species), size = 4, shape = 18) +
      ggrepel::geom_text_repel(data = final_centroids_df,
                               ggplot2::aes(x = Axis1, y = Axis2, label = Species, color = Species),
                               fontface = "bold", box.padding = 0.5, point.padding = 0.5,
                               min.segment.length = 0, show.legend = FALSE)
  }
  
  # Add final touches
  comparison_plot <- comparison_plot +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::guides(color = ggplot2::guide_legend(title = "Species", override.aes = list(size = 3, alpha = 1)))
  
  # Save the plot
  plot_filename <- fs::path_join(c(output_dir, "ENFA_Results", paste0("ENFA_comparison_", target_id, ".png")))
  ggplot2::ggsave(plot_filename, plot = comparison_plot, width = 10, height = 8, units = "in", dpi = 300)
  if (verbose) message("    ENFA comparison plot saved to: ", plot_filename)
}

#' @title .generate_enfa_summary_table
#' @description Creates a CSV file comparing ENFA metrics for a target and its neighbors.
#' @param results A list of ENFA results for the target and closest species.
#' @param output_dir Character. The directory to save the CSV.
#' @param verbose Logical. If TRUE, prints messages.
#' @keywords internal
.generate_enfa_summary_table <- function(results, output_dir, verbose) {
  summary_list <- list()
  
  # Add target species data
  summary_list[[results$metadata$target_species_id]] <- data.frame(
    Species = results$metadata$target_species_id,
    Role = "Target",
    Marginality = results$enfa$marginality,
    Specialization = results$enfa$specialization
  )
  
  # Add neighbor species data
  for (sp_id in names(results$enfa_closest_species)) {
    summary_list[[sp_id]] <- data.frame(
      Species = sp_id,
      Role = "Neighbor",
      Marginality = results$enfa_closest_species[[sp_id]]$marginality,
      Specialization = results$enfa_closest_species[[sp_id]]$specialization
    )
  }
  
  summary_df <- dplyr::bind_rows(summary_list)
  table_filename <- fs::path_join(c(output_dir, "ENFA_Results", paste0("ENFA_comparison_summary_", results$metadata$target_species_id, ".csv")))
  utils::write.csv(summary_df, table_filename, row.names = FALSE)
  if (verbose) message("    ENFA comparison summary table saved to: ", table_filename)
}

#' @title .generate_report_for_target
#' @description Renders an R Markdown report for a single target species analysis.
#' @param results_for_target A list containing all analysis results for the target.
#' @param settings The analysis settings object.
#' @param output_dir The directory where the report will be saved.
#' @param verbose Logical. If TRUE, prints messages.
#' @keywords internal
.generate_report_for_target <- function(results_for_target, settings, output_dir, verbose) {

  # Define the name for the output report file
  report_filename <- paste0("Niche_Analysis_Report_", results_for_target$metadata$target_species_id, ".md")
  
  # Locate the R Markdown template within the package
  # This assumes a standard package structure where the template is in inst/rmarkdown/...
  template_path <- system.file(
    "rmarkdown", "templates", "niche_analysis_report", "skeleton", "skeleton.Rmd",
    package = "fuNiches" 
    )

  if (!nzchar(template_path)) {
    if (verbose) warning("  Report generation skipped: R Markdown template not found in package 'fuNiches'. Expected at: 'inst/rmarkdown/templates/niche_analysis_report/skeleton/skeleton.Rmd'")
    return(NULL)
  }

  # Render the report, passing the results and settings as parameters
  tryCatch({
    rmarkdown::render(
      input = template_path,
      output_file = report_filename,
      output_dir = output_dir,
      params = list(
        analysis_results = results_for_target,
        analysis_settings = settings
      ),
      quiet = TRUE # Keep console clean
    )
    if (verbose) message("  Report successfully generated: ", file.path(output_dir, report_filename))
  }, error = function(e) {
    if (verbose) {
      warning("  Failed to generate report for ", results_for_target$metadata$target_species_id, ". R Markdown error: ", e$message)
    }
  })
}
