#' Load Raw Spatial Occurrences and Raster Stack (Internal)
#'
#' Loads occurrence data from a GPKG file, generates a unique site ID
#' based on coordinates, and loads environmental raster layers.
#'
#' @param gpkg_path Character. Path to the GPKG file.
#' @param raster_dir Character. Path to the directory containing raster files.
#' @return A list containing:
#'   \item{occurrences}{An sf object with species occurrences and the new site ID column.}
#'   \item{env_stack}{A RasterStack object of the environmental layers, or NULL if none found/loaded.}
#'   \item{layer_names}{A character vector of the names of the loaded raster layers.}
#'   \item{original_crs}{The CRS of the original occurrence data.}
#'   \item{name_map}{A named character vector mapping original raster layer names to cleaned names, or NULL.}
#' @noRd
.load_raw_spatial_and_raster <- function(gpkg_path, raster_dir) {
  occurrences_raw <- tryCatch({
    sf::st_read(gpkg_path, quiet = TRUE)
  }, error = function(e) {
    stop("Failed to load GPKG file '", gpkg_path, "': ", e$message)
  })
  if (nrow(occurrences_raw) == 0) {
    stop("No features found in GPKG file: ", gpkg_path)
  }
  original_crs <- sf::st_crs(occurrences_raw)
  env_stack <- NULL
  layer_names <- character(0)
  name_map <- NULL
  if (dir.exists(raster_dir)) {
    env_files <- list.files(raster_dir, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)
    if (length(env_files) > 0) {
      if(length(env_files) > 0) message("Found ", length(env_files), " .tif files in raster directory: ", raster_dir)
      env_stack <- tryCatch({
        raster::stack(env_files)
      }, error = function(e_stack) {
        warning("Error loading rasters into a stack from '", raster_dir, "': ", e_stack$message,
                "\nCheck that all TIF files have the same extent, resolution, and projection. Proceeding without rasters.")
        NULL
      })
      if (!is.null(env_stack)) {
        names_before_cleaning <- names(env_stack)
        cleaned_raster_names <- make.names(names_before_cleaning, unique = TRUE, allow_ = TRUE)
        if (!identical(names_before_cleaning, cleaned_raster_names)) {
          renamed_info <- paste0("'", names_before_cleaning[names_before_cleaning != cleaned_raster_names],
                                 "' -> '", cleaned_raster_names[names_before_cleaning != cleaned_raster_names], "'",
                                 collapse = "; ")
          warning("Some raster layer names were sanitized for R compatibility: ", renamed_info)
          names(env_stack) <- cleaned_raster_names
        }
        layer_names <- names(env_stack)
        name_map <- stats::setNames(cleaned_raster_names, names_before_cleaning)
        message("Raster layers loaded: ", paste(layer_names, collapse = ", "))
              } else {
        warning("No raster layers were successfully loaded from '", raster_dir, "'. Check file integrity and format.")
      }
    } else {
      message("No .tif files found in raster directory: ", raster_dir)
    }
  } else {
    message("Raster directory '", raster_dir, "' not found or not specified. No rasters will be loaded by this helper.")
  }
  return(list(
    occurrences = occurrences_raw,
    env_stack = env_stack,
    layer_names = layer_names,
    original_crs = original_crs,
    name_map = name_map
  ))
}

#' Process Occurrence Attributes (Internal)
#'
#' Handles species name standardization (Genus, Species, species_id),
#' filters by trophic mode, and excludes specified genera.
#'
#' @param occurrences_sf sf object. Raw occurrences from .load_raw_spatial_and_raster.
#' @param scientific_name_col Character. Name of the column containing full scientific names.
#' @param genus_col_name Character. Name of the column for Genus.
#' @param species_col_name Character. Name of the column for the specific epithet.
#' @param final_species_id_col Character. Name for the final species ID column (e.g., "Genus_species").
#' @param trophic_mode_col Character. Optional. Name of the trophic mode column.
#' @param filtered_trophic_modes Character vector. Optional. Vector of trophic modes to keep.
#' @param genera_to_exclude Character vector. Optional. List of genera to exclude.
#' @return A list containing the processed sf object (`data`), 
#'         the name of the genus column used/created (`internal_genus_col`),
#'         and the name of the species epithet column used/created (`internal_species_epithet_col`).
#' @return An sf object with processed attributes, including the new `final_species_id_col`.
#' @noRd
.process_occurrence_attributes <- function(occurrences_sf,
                                           scientific_name_col,
                                           genus_col_name,
                                           species_col_name,
                                           final_species_id_col,
                                           trophic_mode_col,
                                           filtered_trophic_modes,
                                           genera_to_exclude) {

  data_intermediate <- occurrences_sf
  # These will store the names of the columns ultimately used for genus and species epithet
  actual_genus_col_source_name <- NULL
  actual_species_epithet_col_source_name <- NULL

  # --- Priority Logic for Genus and Species ---
  # 1. User-specified genus_col_name and species_col_name
  if (!is.null(genus_col_name) && genus_col_name %in% names(data_intermediate) &&
      !is.null(species_col_name) && species_col_name %in% names(data_intermediate) &&
      !all(is.na(data_intermediate[[genus_col_name]])) && !all(data_intermediate[[genus_col_name]] == "") &&
      !all(is.na(data_intermediate[[species_col_name]])) && !all(data_intermediate[[species_col_name]] == "")) {
    message("Using existing columns '", genus_col_name, "' as Genus and '", species_col_name, "' as Species.")
    actual_genus_col_source_name <- genus_col_name
    actual_species_epithet_col_source_name <- species_col_name

  # 2. Extraction from scientific_name_col
  } else if (!is.null(scientific_name_col) && scientific_name_col %in% names(data_intermediate)) {
    message("Attempting to extract Genus and Species from column '", scientific_name_col, "'.")
    data_intermediate[[scientific_name_col]] <- as.character(data_intermediate[[scientific_name_col]])
    extracted_parts <- .extract_genus_species(data_intermediate[[scientific_name_col]])
    
    # Define names for the new columns that will hold extracted parts
    extracted_genus_col_name <- "funiches_genus_extracted"
    extracted_species_epithet_col_name <- "funiches_species_epithet_extracted"
    
    data_intermediate[[extracted_genus_col_name]] <- extracted_parts$genus
    data_intermediate[[extracted_species_epithet_col_name]] <- extracted_parts$species
    
    # Check if extraction was successful for at least some rows
    if (any(!is.na(data_intermediate[[extracted_genus_col_name]])) && any(!is.na(data_intermediate[[extracted_species_epithet_col_name]]))) {
      actual_genus_col_source_name <- extracted_genus_col_name
      actual_species_epithet_col_source_name <- extracted_species_epithet_col_name
      message("  Successfully extracted Genus and Species to new internal columns: '", actual_genus_col_source_name, "' and '", actual_species_epithet_col_source_name, "'.")
    } else {
      message("  Extraction from '", scientific_name_col, "' did not yield any valid genus/species. Will check for default 'Genus'/'Species' columns if present.")
      # Fallback to default "Genus" and "Species" if extraction fails
      if ("Genus" %in% names(data_intermediate) && "Species" %in% names(data_intermediate) &&
          !all(is.na(data_intermediate[["Genus"]])) && !all(data_intermediate[["Genus"]] == "") &&
          !all(is.na(data_intermediate[["Species"]])) && !all(data_intermediate[["Species"]] == "")) {
        message("Using existing columns 'Genus' as Genus and 'Species' as Species.")
        actual_genus_col_source_name <- "Genus"
        actual_species_epithet_col_source_name <- "Species"
      }
    }
  # 3. Default "Genus" and "Species" columns (if not already handled by fallback from extraction)
  } else if ("Genus" %in% names(data_intermediate) && "Species" %in% names(data_intermediate) &&
             !all(is.na(data_intermediate[["Genus"]])) && !all(data_intermediate[["Genus"]] == "") &&
             !all(is.na(data_intermediate[["Species"]])) && !all(data_intermediate[["Species"]] == "")) {
    message("Using existing columns 'Genus' as Genus and 'Species' as Species.") # This was the message you saw
    actual_genus_col_source_name <- "Genus"
    actual_species_epithet_col_source_name <- "Species"  } else {
    stop(paste0("Cannot determine Genus and Species. Please provide valid '",genus_col_name,"'/'",species_col_name,
                "' columns, or a valid '",scientific_name_col,"' column, or ensure default 'Genus' and 'Species' columns exist and are populated."))
  }

  # Validate that source columns were successfully identified
  if (is.null(actual_genus_col_source_name) || !nzchar(actual_genus_col_source_name) ||
      is.null(actual_species_epithet_col_source_name) || !nzchar(actual_species_epithet_col_source_name)) {
    stop("Failed to identify or create valid source columns for Genus and Species epithet.")
  }

  # Ensure source columns are character and filter NAs/empty
  data_intermediate[[actual_genus_col_source_name]] <- as.character(data_intermediate[[actual_genus_col_source_name]])
  data_intermediate[[actual_species_epithet_col_source_name]] <- as.character(data_intermediate[[actual_species_epithet_col_source_name]])

  data_intermediate <- data_intermediate %>%
    dplyr::filter(
      !is.na(!!sym(actual_genus_col_source_name)) & !!sym(actual_genus_col_source_name) != "",
      !is.na(!!sym(actual_species_epithet_col_source_name)) & !!sym(actual_species_epithet_col_source_name) != ""
    )
  if(nrow(data_intermediate) == 0) {
    warning("No valid records after attempting to derive/process Genus and Species names.")
    return(list(data = data_intermediate, internal_genus_col = actual_genus_col_source_name, internal_species_epithet_col = actual_species_epithet_col_source_name))
  }
  # --- Filter by trophic mode ---
  if (!is.null(trophic_mode_col) && trophic_mode_col %in% names(data_intermediate)) {
    if (!is.null(filtered_trophic_modes) && length(filtered_trophic_modes) > 0) {
      message("Filtering by trophic modes: ", paste(filtered_trophic_modes, collapse=", "))
      data_intermediate <- data_intermediate %>%
        dplyr::mutate(..TrophicModeCleaned = toupper(trimws(as.character(!!sym(trophic_mode_col))))) %>%
        dplyr::filter(..TrophicModeCleaned %in% toupper(trimws(filtered_trophic_modes))) %>%
        dplyr::select(-..TrophicModeCleaned)
      message(nrow(data_intermediate), " records remaining after trophic mode filtering.")
    }
  } else if (!is.null(trophic_mode_col)) {
    warning("Trophic mode column '", trophic_mode_col, "' not found. Skipping trophic mode filtering.")
  }
  # --- Exclude specified genera ---
  if (!is.null(genera_to_exclude) && length(genera_to_exclude) > 0 && any(nzchar(genera_to_exclude))) {
    message("Excluding genera: ", paste(genera_to_exclude, collapse=", "))
    data_intermediate <- data_intermediate %>%
      dplyr::filter(!(!!sym(actual_genus_col_source_name) %in% genera_to_exclude))
    message(nrow(data_intermediate), " records remaining after genus exclusion.")
  }
  # --- Create final species_id ---
  data_intermediate <- data_intermediate %>%
    dplyr::mutate(!!sym(final_species_id_col) := paste(trimws(!!sym(actual_genus_col_source_name)), 
                                                      trimws(!!sym(actual_species_epithet_col_source_name)), 
                                                      sep = "_"))
  # Clean the final_species_id_col
  data_intermediate <- data_intermediate %>%
    dplyr::filter(!is.na(!!sym(final_species_id_col)) & 
                  (!!sym(final_species_id_col) != "NA_NA") &
                  !startsWith(!!sym(final_species_id_col), "_") & 
                  !endsWith(!!sym(final_species_id_col), "_"))
  return(list(    data = data_intermediate,
    internal_genus_col = actual_genus_col_source_name,
    internal_species_epithet_col = actual_species_epithet_col_source_name    
    ))
}

#' Extract Environmental Data, Combine, and Clean (Internal)
#'
#' Extracts environmental data from a raster stack for occurrence locations,
#' combines it with occurrence data, handles NA values, and converts
#' specified categorical variables to factors.
#'
#' @param occurrences_with_ids sf object. Processed occurrences from
#'        `.process_occurrence_attributes`, must include `site_id_col` and `species_id_col`.
#' @param env_raster_stack RasterStack object. Environmental layers.
#' @param vars_to_extract_from_raster Character vector. Names of raster layers to extract.
#'        These should match names in `env_raster_stack`.
#' @param categorical_var_names Character vector. Names of variables to be treated as categorical.
#' @param species_id_col Character. Name of the species ID column
#'        (from rasters or already in `occurrences_with_ids`).
#'
#' @return A tibble with combined occurrence and environmental data, with NAs removed
#'         from environmental variables and categorical variables converted to factors.
#' @noRd
.extract_env_and_combine <- function(occurrences_with_ids,
                                     env_raster_stack,
                                     vars_to_extract_from_raster, 
                                     categorical_var_names,
                                     species_id_col) {
  cols_to_keep_initial <- c(species_id_col, names(attributes(occurrences_with_ids)$sf_column))
  data_for_extraction <- occurrences_with_ids 
  env_values_extracted_df <- NULL
  if (!is.null(env_raster_stack) && raster::nlayers(env_raster_stack) > 0 &&
      !is.null(vars_to_extract_from_raster) && length(vars_to_extract_from_raster) > 0) {
    actual_vars_in_stack <- intersect(vars_to_extract_from_raster, names(env_raster_stack))
    if (length(actual_vars_in_stack) > 0) {
      message("Extracting values for ", length(actual_vars_in_stack), " raster layers: ", paste(actual_vars_in_stack, collapse=", "))
      env_raster_stack_subset <- raster::subset(env_raster_stack, actual_vars_in_stack)
      extracted_values_matrix <- tryCatch({
        raster::extract(env_raster_stack_subset, data_for_extraction)
      }, error = function(e) {
        warning("Error during raster extraction: ", e$message, ". Environmental variables from rasters will be NA.")
        NULL
      })
      if (!is.null(extracted_values_matrix) && nrow(extracted_values_matrix) == nrow(data_for_extraction)) {
        env_values_extracted_df <- tibble::as_tibble(extracted_values_matrix)
        if (ncol(env_values_extracted_df) == length(actual_vars_in_stack)) {
          names(env_values_extracted_df) <- actual_vars_in_stack
        } else {
          warning("Mismatch in number of extracted columns and expected raster variables. Naming might be incorrect.")
        }
      } else {
        warning("Raster extraction failed or returned inconsistent number of rows.")
        if (length(actual_vars_in_stack) > 0) {
          env_values_extracted_df <- tibble::as_tibble(
            matrix(NA_real_, nrow = nrow(data_for_extraction), ncol = length(actual_vars_in_stack))
          )
          names(env_values_extracted_df) <- actual_vars_in_stack
        }
      }
    } else {
      message("No matching raster layers found in stack for extraction based on 'vars_to_extract_from_raster'.")
    }
  } else {
    message("No raster stack provided, or no variables specified for extraction from rasters, or raster stack is empty.")
  }
  combined_data_no_geom <- sf::st_drop_geometry(data_for_extraction)
  if (!is.null(env_values_extracted_df) && nrow(env_values_extracted_df) == nrow(combined_data_no_geom)) {
    new_env_cols <- names(env_values_extracted_df)[!names(env_values_extracted_df) %in% names(combined_data_no_geom)]
    if(length(new_env_cols) > 0) {
      combined_data_tibble <- dplyr::bind_cols(combined_data_no_geom, env_values_extracted_df[, new_env_cols, drop = FALSE])
    } else if (ncol(env_values_extracted_df) > 0) {
      warning("Column name clash or no new columns to add from raster extraction. Check input variable names.")
      combined_data_tibble <- tibble::as_tibble(combined_data_no_geom)
    } else {
      combined_data_tibble <- tibble::as_tibble(combined_data_no_geom)
    }
  } else {
    combined_data_tibble <- tibble::as_tibble(combined_data_no_geom)
  }
  env_vars_for_na_check <- character(0)
  if (!is.null(vars_to_extract_from_raster)) {
    env_vars_for_na_check <- c(env_vars_for_na_check,
                               intersect(vars_to_extract_from_raster, names(combined_data_tibble)))
  }
  if (!is.null(categorical_var_names)) {
    env_vars_for_na_check <- c(env_vars_for_na_check,
                               intersect(categorical_var_names, names(combined_data_tibble)))
  }
  env_vars_for_na_check <- unique(env_vars_for_na_check)
  if (length(env_vars_for_na_check) > 0) {
    rows_before_na_omit <- nrow(combined_data_tibble)
    combined_data_tibble <- combined_data_tibble %>%
      tidyr::drop_na(dplyr::all_of(env_vars_for_na_check))
    rows_after_na_omit <- nrow(combined_data_tibble)
    if (rows_before_na_omit > rows_after_na_omit) {
      message(rows_before_na_omit - rows_after_na_omit,
              " rows removed due to NAs in one or more of the following environmental variables: ",
              paste(env_vars_for_na_check, collapse=", "))
    }
  }
  if (!is.null(categorical_var_names) && length(categorical_var_names) > 0) {
    vars_to_factor <- intersect(categorical_var_names, names(combined_data_tibble))
    if (length(vars_to_factor) > 0) {
      message("Converting to factor: ", paste(vars_to_factor, collapse = ", "))
      existing_vars_to_factor <- vars_to_factor[vars_to_factor %in% names(combined_data_tibble)]
      if(length(existing_vars_to_factor) > 0){
        combined_data_tibble <- combined_data_tibble %>%
          dplyr::mutate(dplyr::across(dplyr::all_of(existing_vars_to_factor), as.factor))
      }
    }
  }
  return(combined_data_tibble)
}
