#!/usr/bin/env Rscript

# R Script: Record Lookup and ANC Data Extraction
# Purpose: Extract absolute neutrophil count (ANC) data for BNHL records at specific time points
# Author: Generated for AG Rejeski analysis
# Date: 2025-12-13

# Load required libraries
library(tidyverse)
library(readr)

# Define file paths
master_file <- "/Users/linus/Library/Mobile Documents/com~apple~CloudDocs/10 AG Rejeski/Paper/Eosinophile_ SS/master_lmu_06_11_2025.csv"
lab_file <- "/Users/linus/Library/Mobile Documents/com~apple~CloudDocs/10 AG Rejeski/Paper/Eosinophile_ SS/lab_export.csv"
eosinophils_file <- "/Users/linus/Library/Mobile Documents/com~apple~CloudDocs/10 AG Rejeski/Paper/Eosinophile_ SS/eosinophils abs_continous.csv"
output_dir <- "/Users/linus/Library/Mobile Documents/com~apple~CloudDocs/10 AG Rejeski/Paper/Eosinophile_ SS"
output_file <- file.path(output_dir, "absolute_neutrophil_count.csv")

# Define time points
time_labels <- c("-90", "-60", "-30", "-21", "-14", "-5", "0", "3", "7", "14", "21",
                 "M1", "M2", "M3", "M6", "M12", "M18", "M24", "M30", "M36", "M42", "M48", "M54")
days_numeric <- c(-90, -60, -30, -21, -14, -5, 0, 3, 7, 14, 21,
                  30, 60, 90, 180, 360, 540, 720, 900, 1080, 1260, 1440, 1620)

# Tolerance window for matching (±3 days)
tolerance <- 3

# Step 1: Read eosinophils file to get entity_ids from row 1 (column headers)
cat("Reading eosinophils file to extract entity_ids...\n")
eosinophils_headers <- names(read_csv(eosinophils_file, n_max = 0, show_col_types = FALSE))

# Extract entity_ids (bnhl_x identifiers) - exclude Time and days columns
entity_ids <- eosinophils_headers[!eosinophils_headers %in% c("Time", "days")]
cat(sprintf("Found %d entity_ids in eosinophils file: %s\n",
            length(entity_ids),
            paste(head(entity_ids, 10), collapse = ", ")))
if (length(entity_ids) > 10) {
  cat(sprintf("... and %d more\n", length(entity_ids) - 10))
}

# Step 2: Read master file and match entity_id to record_id
cat("\nReading master file...\n")
master_data <- read_csv(master_file, show_col_types = FALSE)

# Match entity_ids to record_ids
cat("Matching entity_ids to record_ids...\n")
record_to_bnhl <- master_data %>%
  filter(entity_id %in% entity_ids) %>%
  select(record_id, entity_id) %>%
  distinct() %>%
  filter(!is.na(record_id), !is.na(entity_id)) %>%
  rename(bnhl_id = entity_id)

cat(sprintf("Found %d unique record IDs matching entity_ids\n", nrow(record_to_bnhl)))

# Check for entity_ids that were not found
missing_entities <- setdiff(entity_ids, record_to_bnhl$bnhl_id)
if (length(missing_entities) > 0) {
  cat(sprintf("WARNING: %d entity_ids from eosinophils file not found in master file:\n",
              length(missing_entities)))
  cat(paste(head(missing_entities, 10), collapse = ", "), "\n")
  if (length(missing_entities) > 10) {
    cat(sprintf("... and %d more\n", length(missing_entities) - 10))
  }
}

# Step 3: Read lab export file and filter for relevant record_ids
cat("\nReading lab export file...\n")
lab_data <- read_csv(lab_file, show_col_types = FALSE)

# Filter for relevant record_ids
lab_filtered <- lab_data %>%
  filter(record_id %in% record_to_bnhl$record_id) %>%
  select(record_id, day_therapy, anc_masch_abs) %>%
  filter(!is.na(day_therapy), !is.na(anc_masch_abs))

cat(sprintf("Found %d lab measurements for %d unique records\n",
            nrow(lab_filtered),
            n_distinct(lab_filtered$record_id)))

# Step 4: Match time points with tolerance window
cat("\nMatching time points (tolerance = ±", tolerance, "days)...\n")

# Initialize tracking for imputations
imputation_log <- list()
total_matches <- 0
exact_matches <- 0
tolerance_matches <- 0

# Create output data frame
output_data <- data.frame(
  Time = time_labels,
  days = days_numeric,
  stringsAsFactors = FALSE
)

# Process each record
for (i in 1:nrow(record_to_bnhl)) {
  rec_id <- record_to_bnhl$record_id[i]
  bnhl_id <- record_to_bnhl$bnhl_id[i]

  # Get data for this record
  rec_data <- lab_filtered %>%
    filter(record_id == rec_id)

  # Initialize column for this bnhl_id
  anc_values <- rep(NA_real_, length(days_numeric))

  # Match each time point
  for (j in 1:length(days_numeric)) {
    target_day <- days_numeric[j]

    # Find measurements within tolerance window
    matches <- rec_data %>%
      filter(abs(day_therapy - target_day) <= tolerance)

    if (nrow(matches) > 0) {
      total_matches <- total_matches + 1

      # Check if any are exact matches
      exact <- matches %>% filter(day_therapy == target_day)

      if (nrow(exact) > 0) {
        exact_matches <- exact_matches + 1
        # Take lowest value if multiple exact matches
        anc_values[j] <- min(exact$anc_masch_abs)
        if (nrow(exact) > 1) {
          imputation_log[[length(imputation_log) + 1]] <- sprintf(
            "%s at day %d: %d exact matches found, took minimum (%.2f)",
            bnhl_id, target_day, nrow(exact), anc_values[j]
          )
        }
      } else {
        tolerance_matches <- tolerance_matches + 1
        # Take lowest value from tolerance matches
        anc_values[j] <- min(matches$anc_masch_abs)

        # Log this imputation
        closest_day <- matches$day_therapy[which.min(abs(matches$day_therapy - target_day))]
        imputation_log[[length(imputation_log) + 1]] <- sprintf(
          "%s at day %d: No exact match, used day %d (within tolerance, value=%.2f)",
          bnhl_id, target_day, closest_day, anc_values[j]
        )
      }
    }
  }

  # Add column to output
  output_data[[bnhl_id]] <- anc_values
}

# Step 5: Save output file
cat("\nSaving output to:", output_file, "\n")
write_csv(output_data, output_file)

# Step 6: Print summary statistics
cat("\n=== SUMMARY ===\n")
cat(sprintf("Total records processed: %d\n", nrow(record_to_bnhl)))
cat(sprintf("Total time points: %d\n", length(days_numeric)))
cat(sprintf("Possible data points: %d\n", nrow(record_to_bnhl) * length(days_numeric)))
cat(sprintf("Total matches found: %d\n", total_matches))
cat(sprintf("  - Exact matches: %d\n", exact_matches))
cat(sprintf("  - Tolerance matches (±%d days): %d\n", tolerance, tolerance_matches))
cat(sprintf("Missing data points: %d\n",
            nrow(record_to_bnhl) * length(days_numeric) - total_matches))

# Print imputation log
if (length(imputation_log) > 0) {
  cat("\n=== IMPUTATION LOG ===\n")
  cat(sprintf("Total imputations/warnings: %d\n\n", length(imputation_log)))
  for (log_entry in imputation_log) {
    cat(log_entry, "\n")
  }
}

cat("\n=== COMPLETE ===\n")
cat("Output file created successfully!\n")
