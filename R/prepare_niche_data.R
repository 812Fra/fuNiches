#' Prepare Niche Data for Analysis
#'
#' This function loads species occurrence data (GPKG) and environmental raster layers,
#' processes them, and returns a structured object ready for niche analysis.
#'
#' @param occurrences_gpkg_path Character. Path to the GPKG file containing species occurrences.
#'        The GPKG must contain point geometries.
#' @param scientific_name_col Character. Name of the column in the GPKG's attribute table
#'        containing the full scientific name (e.g., "Genus species"). Used if `genus_col_name`
#'        and `species_col_name` are not found or not provided. Default is "scientificName".
#' @param genus_col_name Character. Optional. Name of the column for Genus.
#'        If provided and found, this column will be used. Otherwise, the function
#'        attempts to find a column named "Genus", or derive it from `scientific_name_col`.
#'        Default is "Genus".
#' @param species_col_name Character. Optional. Name of the column for the specific epithet.
#'        If provided and found, this column will be used. Otherwise, the function
#'        attempts to find a column named "Species", or derive it from `scientific_name_col`.
#'        Default is "Species".
#' @param species_id_col Character. Name of the column that will be created (or used if existing and correctly formatted) by combining
#'        Genus and Species (e.g., "Genus_species"). Default is "species_id".
#' @param raster_dir_path Character. Path to the directory containing environmental raster files (TIF format).
#' @param categorical_var_names Character vector. Names of the environmental variables
#'        (matching raster layer names without .tif extension) that should be treated as categorical.
#'
#' @param min_occurrences_per_species Integer. Minimum number of valid occurrences for a
#'        species to be retained in the final dataset. Default is 10.
#' @param continuous_var_labels Named list. Optional. Custom display labels for continuous variables.
#'        Example: `list(bio1 = "Annual Mean Temp", aspect = "Aspect (degrees)")`.
#'        Variable names in the list must match raster layer names.
#' @param categorical_var_labels Named list of named vectors/lists. Optional. Custom display labels for
#'        categorical variable levels. Example:
#'        `list(SoilType = c("1" = "Sandy", "2" = "Clay"), LandCover = c("1" = "Quercus ilex Forest", "2" = "Pinus radiata Forest))`.
#'        Variable names in the list must match categorical raster layer names.
#' @param trophic_mode_col Character. Optional. Name of the column in the GPKG containing
#'        trophic mode information.
#' @param filtered_trophic_modes Character vector. Optional. If `trophic_mode_col` is provided,
#'        only occurrences with these trophic modes will be kept.
#' @param genera_to_exclude Character vector. Optional. List of genera to exclude from the analysis.
#'        Default is NULL.
#'
#' @return A list object of class `niche_data` containing:
#'   \item{data}{A tibble with processed occurrence data, extracted environmental values,
#'               the site ID column, and the species ID column.}
#'   \item{metadata}{A list containing information about the processed data:
#'     \itemize{
#'       \item `site_id_col`: Name of the site ID column.
#'       \item `genus_col_name_processed`: Name of the genus column used/created.
#'       \item `species_epithet_col_name_processed`: Name of the species epithet column used/created.
#'       \item `species_id_col`: Name of the species ID column.
#'       \item `occurrence_crs`: Original CRS of the occurrence data.
#'       \item `env_vars`: A list with names of all, continuous, and categorical environmental variables used.
#'       \item `labels`: A list containing `continuous_var_labels` and `categorical_var_labels` (merged user input with defaults).
#'       \item `categorical_vars_info`: A list where each element is a categorical variable name,
#'              and its value is a named vector of its levels and their labels.
#'       \item `parameters_used`: A list of key parameters used during data preparation.
#'     }}
#'   \item{raster_stack_info}{Information about the loaded raster stack (e.g., names, resolution, extent),
#'                            or NULL if no raster layers were processed.}
#'
#' @export
#' @export
#' @importFrom sf st_read st_coordinates st_drop_geometry st_crs
#' @importFrom raster stack subset extract res extent crs nlayers
#' @importFrom dplyr %>% filter mutate select group_by summarise n distinct across all_of left_join count sym
#' @importFrom tidyr separate unite drop_na pivot_longer
#' @importFrom tibble as_tibble tibble rownames_to_column
#' @importFrom purrr map map_chr map_dbl map_dfr
#' @importFrom stats setNames na.omit cor dist prcomp lm median quantile fisher.test chisq.test p.adjust as.formula
#' @importFrom utils modifyList txtProgressBar setTxtProgressBar write.csv head read.csv
#' @importFrom methods as
#' @importFrom grDevices png dev.off adjustcolor
#' @importFrom graphics par hist legend mtext
#'
#' @examples
#' \dontrun{
#'   # Assuming you have:
#'   # occurrences.gpkg in "path/to/data"
#'   # raster_layers/ (containing bio1.tif, soil_type.tif) in "path/to/data"
#'
#'   niche_object <- prepare_niche_data(
#'     occurrences_gpkg_path = "path/to/data/occurrences.gpkg",
#'     raster_dir_path = "path/to/data/raster_layers/",
#'     categorical_var_names = c("soil_type"),
#'     scientific_name_col = "nome_scientifico", # if different from default
#'     genus_col_name = "genus",             # if different from default
#'     species_col_name = "species",           # if different from default
#'     min_occurrences_per_species = 5,
#'     continuous_var_labels = list(bio1 = "Annual Mean Temperature (°C)"),
#'     categorical_var_labels = list(soil_type = c("1" = "Sandy Loam", "2" = "Clay Loam"))
#'   )
#'
#'   # Access the processed data
#'   # head(niche_object$data)
#'   # print(niche_object$metadata$env_vars)
#' }
prepare_niche_data <- function(occurrences_gpkg_path,
                               raster_dir_path = NULL, 
                               categorical_var_names = NULL,
                               scientific_name_col = "scientificName",
                               genus_col_name = "Genus",
                               species_col_name = "Species",
                               species_id_col = "species_id",
                               min_occurrences_per_species = 10,
                               continuous_var_labels = list(),
                               categorical_var_labels = list(),
                               trophic_mode_col = NULL,
                               filtered_trophic_modes = NULL,
                               genera_to_exclude  = NULL) {
  # --- 1. Argument Validation ---
  if (!file.exists(occurrences_gpkg_path)) {
    stop("Occurrences GPKG file not found: ", occurrences_gpkg_path)
  }
  # Check raster_dir_path only if it's provided and not NULL
  if (!is.null(raster_dir_path) && !dir.exists(raster_dir_path)) {
    # It's not a fatal error if raster_dir_path is invalid but no rasters are needed/found
    # However, if categorical_var_names are provided, rasters are expected.
    if (!is.null(categorical_var_names) && length(categorical_var_names) > 0) {
      warning("Raster directory not found: ", raster_dir_path,
              ". Categorical variables might not be processed correctly if they come from rasters.")
    } else {
      message("Raster directory not found: ", raster_dir_path, ". No raster-based environmental variables will be processed.")
    }
  }
  if (length(categorical_var_names) == 0 && !is.null(categorical_var_names)) {
    message("`categorical_var_names` is an empty vector. No variables will be treated as categorical from rasters unless already factors in GPKG.")
  }
  # --- 2. Load Raw Spatial and Raster Data ---
  raw_data_list <- .load_raw_spatial_and_raster(gpkg_path = occurrences_gpkg_path, raster_dir = raster_dir_path)
  env_raster_stack <- raw_data_list$env_stack
  loaded_layer_names <- raw_data_list$layer_names
  original_crs <- raw_data_list$original_crs
  occurrences_sf <- raw_data_list$occurrences
  if (nrow(occurrences_sf) == 0) {
    stop("No occurrence data loaded or all occurrences were invalid.")
  }
  # --- 3. Process Occurrence Data (Names, Trophic Mode) ---
  processed_occurrences_list <- .process_occurrence_attributes(
    occurrences_sf = occurrences_sf,
    scientific_name_col = scientific_name_col,
    genus_col_name = genus_col_name,
    species_col_name = species_col_name,
    final_species_id_col = species_id_col,
    trophic_mode_col = trophic_mode_col,
    filtered_trophic_modes = filtered_trophic_modes,
    genera_to_exclude = genera_to_exclude
  )
  processed_occurrences <- processed_occurrences_list$data
  # These are the names of the columns that .process_occurrence_attributes decided to use
  # for genus and species epithet, or created if extracted.
  g_col_source <- processed_occurrences_list$internal_genus_col
  s_col_source <- processed_occurrences_list$internal_species_epithet_col # Assume .process_occurrence_attributes returns this

  # .process_occurrence_attributes is now responsible for creating the final species_id_col
  # (named according to the 'species_id_col' parameter of prepare_niche_data)
  # and ensuring it's correctly formatted (e.g., "Genus_Epithet").
  # We just check if the necessary info was returned for metadata and if species_id_col exists.
  if (is.null(g_col_source) || !nzchar(g_col_source) || 
      is.null(s_col_source) || !nzchar(s_col_source)) {
    stop(paste0("Internal error: .process_occurrence_attributes did not return valid source column names ",
                "for genus ('", g_col_source, "') or species epithet ('", s_col_source, "'). This is needed for metadata."))
  }
  
  if (!species_id_col %in% names(processed_occurrences)) {
    stop(paste0("Internal error: The designated species ID column '", species_id_col, 
                "' was not found in the data returned by .process_occurrence_attributes. ",
                "This column should have been created by the helper function."))
  }  
  
  if (nrow(processed_occurrences) == 0) {
    stop("No occurrence data remaining after attribute processing by .process_occurrence_attributes (check names, trophic mode, exclusions, and species ID creation).")
  }
  message(nrow(processed_occurrences), " records remaining after attribute processing.")
  # --- 4. Extract Environmental Data and Combine ---
  vars_for_extraction <- loaded_layer_names
  combined_data_tibble <- .extract_env_and_combine(
    occurrences_with_ids = processed_occurrences,
    env_raster_stack = env_raster_stack,
    vars_to_extract_from_raster = vars_for_extraction,
    categorical_var_names = categorical_var_names, species_id_col = species_id_col
  )
  if (nrow(combined_data_tibble) == 0) {
    stop("No data remaining after environmental value extraction and NA removal.")
  }
  message(nrow(combined_data_tibble), " records remaining after environmental extraction and NA handling.")
  # --- 5. Filter by Minimum Occurrences per Species ---
  species_counts <- combined_data_tibble %>%
    dplyr::count(!!sym(species_id_col))
  species_to_keep <- species_counts %>%
    dplyr::filter(n >= min_occurrences_per_species) %>%
    dplyr::pull(!!sym(species_id_col))
  final_data_tibble <- combined_data_tibble %>%
    dplyr::filter(!!sym(species_id_col) %in% species_to_keep)
  if (nrow(final_data_tibble) == 0) {
    stop("No species meet the minimum occurrence threshold of ", min_occurrences_per_species)
  }
  message(
    length(unique(final_data_tibble[[species_id_col]])), " species remaining after filtering by min occurrences (>= ",
    min_occurrences_per_species, "). Total records: ", nrow(final_data_tibble)
  )
  # --- 6. Finalize Variable Information and Labels ---
  all_env_vars_in_final_data <- setdiff(names(final_data_tibble),
                                        c(names(attributes(occurrences_sf)$sf_column),
                                          species_id_col,
                                          scientific_name_col, genus_col_name, species_col_name,
                                          "Genus_derived", "Species_epithet_derived",
                                          trophic_mode_col))
  all_env_vars_in_final_data <- all_env_vars_in_final_data[all_env_vars_in_final_data %in% loaded_layer_names | all_env_vars_in_final_data %in% categorical_var_names]
  final_categorical_vars <- intersect(all_env_vars_in_final_data, categorical_var_names)
  final_continuous_vars <- sort(setdiff(all_env_vars_in_final_data, final_categorical_vars))
  default_cont_labels_vec <- stats::setNames(final_continuous_vars, final_continuous_vars)
  if (length(continuous_var_labels) > 0) { 
    name_mapping_from_original_to_cleaned <- raw_data_list$name_map 
    for (user_provided_original_name in names(continuous_var_labels)) {
      user_provided_label_value <- continuous_var_labels[[user_provided_original_name]]
      cleaned_name_for_this_label <- NA_character_
      if (!is.null(name_mapping_from_original_to_cleaned) && user_provided_original_name %in% names(name_mapping_from_original_to_cleaned)) {
        cleaned_name_for_this_label <- name_mapping_from_original_to_cleaned[user_provided_original_name]
      } else if (user_provided_original_name %in% final_continuous_vars) {
        cleaned_name_for_this_label <- user_provided_original_name 
      }
      if (!is.na(cleaned_name_for_this_label) && cleaned_name_for_this_label %in% names(default_cont_labels_vec)) {
        default_cont_labels_vec[cleaned_name_for_this_label] <- user_provided_label_value
      } else {
        warning("Label provided for continuous variable '", user_provided_original_name, 
                "' could not be matched to a final continuous variable or its cleaned name. Label ignored.")
      }
    }
  }
  final_cont_labels_vector_for_output <- default_cont_labels_vec
  final_cat_labels_processed <- list()
  final_cat_vars_info <- list()
  for (cat_var in final_categorical_vars) {
    if (cat_var %in% names(final_data_tibble) && is.factor(final_data_tibble[[cat_var]])) {
      actual_levels <- levels(final_data_tibble[[cat_var]])
      default_level_labels <- stats::setNames(actual_levels, actual_levels)
      user_provided_level_labels <- list()
      if (cat_var %in% names(categorical_var_labels)) {
        user_provided_level_labels <- categorical_var_labels[[cat_var]]
        if((!is.list(user_provided_level_labels) && !is.vector(user_provided_level_labels)) || is.null(names(user_provided_level_labels))){
          warning("Labels for categorical variable '", cat_var, "' are not a named list/vector. Using default level names.")
          user_provided_level_labels <- list()
        }
      }
  merged_level_labels <- utils::modifyList(as.list(default_level_labels), as.list(user_provided_level_labels))
      final_level_labels_for_var <- stats::setNames(
        purrr::map_chr(actual_levels, ~ ifelse(.x %in% names(merged_level_labels), merged_level_labels[[.x]], .x)),
        actual_levels
      )
      final_cat_labels_processed[[cat_var]] <- final_level_labels_for_var
      final_cat_vars_info[[cat_var]] <- final_level_labels_for_var
    }
  }
  # --- 7. Construct the Output Object ---
  niche_data_object <- list(
    data = final_data_tibble,
    metadata = list(
      species_id_col = species_id_col,
      occurrence_crs = original_crs,
      genus_col_name_processed = g_col_source, # The source column for genus
      species_epithet_col_name_processed = s_col_source, # The source column for species epithet
      env_vars = list(
        all = all_env_vars_in_final_data,
        continuous = final_continuous_vars,
        categorical = final_categorical_vars
      ),
      labels = list(
        continuous_vars = final_cont_labels_vector_for_output, # Use the processed vector
        categorical_vars = final_cat_labels_processed
      ),
      categorical_vars_info = final_cat_vars_info,
      parameters_used = list(
        occurrences_gpkg_path = occurrences_gpkg_path,
        raster_dir_path = raster_dir_path,
        min_occurrences_per_species = min_occurrences_per_species,
        genera_to_exclude = genera_to_exclude
      )
    ),
    raster_stack_info = if (!is.null(env_raster_stack)) {
      list(names = names(env_raster_stack),
           resolution = raster::res(env_raster_stack),
           extent = raster::extent(env_raster_stack),
           crs = raster::crs(env_raster_stack))
    } else { NULL }
  )
  class(niche_data_object) <- c("niche_data", "list")
  message("Niche data preparation complete.")
  return(niche_data_object)
}
