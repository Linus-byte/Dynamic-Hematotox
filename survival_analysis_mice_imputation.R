# ============================================
# SURVIVAL ANALYSIS BY CAR-HEMATOTOX TRAJECTORY
# Baseline ht_rg + Day 14 Score Stratification
# 4 Groups: high/high, high/low, low/high, low/low
# WITH MICE MULTIPLE IMPUTATION
# ============================================

library(tidyverse)
library(survival)
library(survminer)
library(cowplot)
library(broom)
library(mice)  # For multiple imputation

# ============================================
# CONFIGURATION
# ============================================

main_input_path <- '/Users/linus/Library/Mobile Documents/com~apple~CloudDocs/10 AG Rejeski/Dynamic_hemato/R/master_lmu_17_11_2025_lk.csv'
lab_input_path <- "/Users/linus/Library/Mobile Documents/com~apple~CloudDocs/10 AG Rejeski/dynamic_hemato/R/lab_export.csv"
output_path <- "/Users/linus/Library/Mobile Documents/com~apple~CloudDocs/10 AG Rejeski/dynamic_hemato/R"

# Imputation settings
imp_nr <- 10  # Number of imputed datasets (as in supplemental)

# ============================================
# COLOR SCHEME (4 groups)
# ============================================

colors_trajectory <- c(
  "high/high" = "#9E1972",
  "low/high"  = "#E67BB0",
  "high/low"  = "#0084BC",
  "low/low"   = "#00588B"
)

# ============================================
# HELPER FUNCTIONS
# ============================================

format_pval_short <- function(p) {
  suppressWarnings(ifelse(is.na(p), NA, ifelse(p < 0.001, "<0.001", format.pval(p, digits = 3))))
}

calculate_hematotox <- function(plt, anc, hb, crp, ferritin) {

  score <- 0

  if (!is.na(plt)) {
    if (plt > 175) score <- score + 0
    else if (plt >= 75) score <- score + 1
    else score <- score + 2
  } else {
    return(NA_real_)
  }

  if (!is.na(anc)) {
    if (anc > 1.2) score <- score + 0
    else score <- score + 1
  } else {
    return(NA_real_)
  }

  if (!is.na(hb)) {
    if (hb > 9.0) score <- score + 0
    else score <- score + 1
  } else {
    return(NA_real_)
  }

  if (!is.na(crp)) {
    if (crp < 3.0) score <- score + 0
    else score <- score + 1
  } else {
    return(NA_real_)
  }

  if (!is.na(ferritin)) {
    if (ferritin < 650) score <- score + 0
    else if (ferritin <= 2000) score <- score + 1
    else score <- score + 2
  } else {
    return(NA_real_)
  }

  return(score)
}

# ============================================
# STEP 1: LOAD MAIN DATASET
# ============================================

cat("=== STEP 1: LOADING MAIN DATASET ===\n\n")

df_main <- read_csv2(main_input_path)
cat("Main dataset loaded: n =", nrow(df_main), "\n")

df_main <- df_main %>%
  mutate(
    record_id = as.character(record_id),
    ht_rg = as.numeric(ht_rg),
    pfs_ev = as.numeric(pfs_ev),
    os_ev = as.numeric(os_ev),
    pfs_d = as.numeric(pfs_d),
    os_d = as.numeric(os_d),
    baseline_group = case_when(
      ht_rg == 1 ~ "high",
      ht_rg == 0 ~ "low",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(baseline_group)) %>%
  select(record_id, baseline_group, ht_rg, pfs_ev, os_ev, pfs_d, os_d)

cat("Patients with valid ht_rg: n =", nrow(df_main), "\n")
cat("\nBaseline ht_rg distribution:\n")
print(table(df_main$baseline_group, useNA = "ifany"))
cat("\n")

# ============================================
# STEP 2: LOAD AND WRANGLE LAB DATA
# ============================================

cat("=== STEP 2: LOADING LAB DATA ===\n\n")

df_lab <- read_csv2(lab_input_path)
cat("Lab dataset loaded: n =", nrow(df_lab), "\n")

df_lab <- df_lab %>%
  mutate(
    record_id = as.character(record_id),
    day_therapy = as.numeric(day_therapy),
    plt = case_when(
      str_trim(as.character(plt)) == "s.Bem." ~ NA_character_,
      TRUE ~ as.character(plt)
    ),
    plt = as.numeric(str_replace(plt, ",", ".")),
    anc_masch_abs = case_when(
      str_trim(as.character(anc_masch_abs)) == "<0,10" ~ "0.1",
      TRUE ~ as.character(anc_masch_abs)
    ),
    anc_masch_abs = as.numeric(str_replace(anc_masch_abs, ",", ".")),
    crp = as.numeric(str_replace(as.character(crp), ",", ".")),
    ferritin = as.numeric(str_replace(as.character(ferritin), ",", ".")),
    hb = as.numeric(str_replace(as.character(hb), ",", "."))
  )

cat("Lab values cleaned\n\n")

# ============================================
# STEP 3: DIAGNOSE DATA AVAILABILITY
# ============================================

cat("=== STEP 3: DIAGNOSING DATA AVAILABILITY ===\n\n")

patients_with_labs <- df_lab %>%
  filter(record_id %in% df_main$record_id) %>%
  pull(record_id) %>%
  unique()

cat("Patients in main dataset:", nrow(df_main), "\n")
cat("Patients with ANY lab data:", length(patients_with_labs), "\n")
cat("Patients without lab data:", nrow(df_main) - length(patients_with_labs), "\n\n")

df_lab_d14_window <- df_lab %>%
  filter(day_therapy >= 10 & day_therapy <= 18) %>%
  filter(record_id %in% df_main$record_id)

patients_in_window <- unique(df_lab_d14_window$record_id)
cat("Patients with lab data in day 10-18 window:", length(patients_in_window), "\n\n")

cat("Variable availability across ALL timepoints for patients in main:\n")
all_patient_labs <- df_lab %>% filter(record_id %in% df_main$record_id)

for (var in c("plt", "anc_masch_abs", "hb", "crp", "ferritin")) {
  pts_with_var <- all_patient_labs %>%
    filter(!is.na(!!sym(var))) %>%
    pull(record_id) %>%
    unique() %>%
    length()
  cat("  ", var, ": ", pts_with_var, "/", nrow(df_main), " patients have at least 1 measurement\n")
}
cat("\n")

# ============================================
# STEP 4: EXTRACT DAY 14 VALUES
# ============================================

cat("=== STEP 4: EXTRACTING DAY 14 VALUES ===\n\n")

get_d14_value <- function(days, values) {
  df <- tibble(day = days, val = values) %>%
    filter(!is.na(val))

  if (nrow(df) == 0) return(NA_real_)

  d14_exact <- df %>% filter(day == 14)
  if (nrow(d14_exact) > 0) {
    return(d14_exact$val[1])
  }

  df <- df %>%
    mutate(dist_from_14 = abs(day - 14)) %>%
    arrange(dist_from_14)

  return(df$val[1])
}

d14_values <- df_lab_d14_window %>%
  group_by(record_id) %>%
  summarise(
    plt_d14 = get_d14_value(day_therapy, plt),
    anc_d14 = get_d14_value(day_therapy, anc_masch_abs),
    hb_d14 = get_d14_value(day_therapy, hb),
    crp_d14 = get_d14_value(day_therapy, crp),
    ferritin_d14 = get_d14_value(day_therapy, ferritin),
    .groups = 'drop'
  )

cat("Day 14 values extracted for", nrow(d14_values), "patients\n\n")

cat("Day 14 value availability (from window 10-18):\n")
cat("  plt:      ", sum(!is.na(d14_values$plt_d14)), "/", nrow(d14_values), "\n")
cat("  anc:      ", sum(!is.na(d14_values$anc_d14)), "/", nrow(d14_values), "\n")
cat("  hb:       ", sum(!is.na(d14_values$hb_d14)), "/", nrow(d14_values), "\n")
cat("  crp:      ", sum(!is.na(d14_values$crp_d14)), "/", nrow(d14_values), "\n")
cat("  ferritin: ", sum(!is.na(d14_values$ferritin_d14)), "/", nrow(d14_values), "\n\n")

d14_values <- d14_values %>%
  mutate(
    n_missing_initial = is.na(plt_d14) + is.na(anc_d14) + is.na(hb_d14) +
                        is.na(crp_d14) + is.na(ferritin_d14)
  )

cat("BEFORE imputation:\n")
cat("  Complete (0 missing):  n =", sum(d14_values$n_missing_initial == 0), "\n")
cat("  1 missing value:       n =", sum(d14_values$n_missing_initial == 1), "\n")
cat("  2 missing values:      n =", sum(d14_values$n_missing_initial == 2), "\n")
cat("  3+ missing values:     n =", sum(d14_values$n_missing_initial >= 3), "\n\n")

# ============================================
# STEP 5: MICE MULTIPLE IMPUTATION
# ============================================

cat("=== STEP 5: MICE MULTIPLE IMPUTATION ===\n\n")

# Merge day 14 values with main data for imputation
# Including baseline_group as it may inform imputation
df_for_imputation <- df_main %>%
  inner_join(d14_values, by = "record_id") %>%
  mutate(
    # Convert baseline_group to numeric for imputation
    baseline_numeric = ifelse(baseline_group == "high", 1, 0)
  )

cat("Patients available for imputation: n =", nrow(df_for_imputation), "\n\n")

# Define variables for MICE imputation
# imputerVars: Variables that inform the imputation (predictors)
# imputedVars: Variables to be imputed

imputerVars <- c("plt_d14", "anc_d14", "hb_d14", "crp_d14", "ferritin_d14",
                 "baseline_numeric", "ht_rg")
imputedVars <- c("plt_d14", "anc_d14", "hb_d14", "crp_d14", "ferritin_d14")

# Prepare data for mice - select only relevant columns
mice_data <- df_for_imputation %>%
  select(record_id, all_of(imputerVars), pfs_ev, os_ev, pfs_d, os_d)

# Identify which variables actually have missing values
missVars <- names(mice_data)[colSums(is.na(mice_data)) > 0]
cat("Variables with missing values:", paste(missVars, collapse = ", "), "\n\n")

# Update imputedVars to only include those with actual missingness
imputedVars <- intersect(imputedVars, missVars)
cat("Variables to be imputed:", paste(imputedVars, collapse = ", "), "\n\n")

# Build predictor matrix (following supplemental file approach)
all_columns <- colnames(mice_data)
predictorMatrix <- matrix(0, ncol = length(all_columns), nrow = length(all_columns))
rownames(predictorMatrix) <- all_columns
colnames(predictorMatrix) <- all_columns

# Create imputer matrix (which variables can be used as predictors)
# Lab values + baseline info can predict each other
imputerMatrix <- predictorMatrix
for (col in intersect(imputerVars, all_columns)) {
  imputerMatrix[, col] <- 1
}

# Create imputed matrix (which variables will be imputed)
imputedMatrix <- predictorMatrix
for (row in intersect(imputedVars, all_columns)) {
  imputedMatrix[row, ] <- 1
}

# Construct full predictor matrix
predictorMatrix <- imputerMatrix * imputedMatrix

# Diagonals must be zeros (a variable cannot impute itself)
diag(predictorMatrix) <- 0

# record_id should never be used as predictor or be imputed
predictorMatrix["record_id", ] <- 0
predictorMatrix[, "record_id"] <- 0

cat("Predictor matrix constructed\n")
cat("Variables used as predictors for imputation:\n")
for (var in imputedVars) {
  predictors <- names(which(predictorMatrix[var, ] == 1))
  cat("  ", var, "predicted by:", paste(predictors, collapse = ", "), "\n")
}
cat("\n")

# Set up imputation methods
# Use "pmm" (predictive mean matching) for continuous variables
# Leave blank for variables that don't need imputation
methods <- rep("", length(all_columns))
names(methods) <- all_columns
methods[imputedVars] <- "pmm"

cat("Performing MICE imputation with", imp_nr, "datasets...\n")
cat("Method: Predictive Mean Matching (pmm)\n\n")

# Perform the actual imputation
set.seed(42)  # For reproducibility
imputedDataSets <- mice(
  data = mice_data,
  m = imp_nr,                      # Number of imputed datasets
  predictorMatrix = predictorMatrix,
  method = methods,
  maxit = 50,                      # Maximum iterations
  print = FALSE
)

cat("MICE imputation completed!\n\n")

# Check convergence
cat("Imputation convergence check:\n")
print(imputedDataSets$loggedEvents)
cat("\n")

# ============================================
# STEP 6: ANALYZE IMPUTATION RESULTS
# ============================================

cat("=== STEP 6: ANALYZING IMPUTATION RESULTS ===\n\n")

# Get summary of imputations
cat("Summary of imputed values:\n")
for (var in imputedVars) {
  n_imputed <- sum(is.na(mice_data[[var]]))
  if (n_imputed > 0) {
    cat("\n", var, "(", n_imputed, "values imputed):\n")

    # Get imputed values across all datasets
    imputed_vals <- sapply(1:imp_nr, function(i) {
      complete(imputedDataSets, i)[[var]][is.na(mice_data[[var]])]
    })

    if (is.matrix(imputed_vals)) {
      cat("  Mean across imputations:", round(mean(imputed_vals), 2), "\n")
      cat("  SD across imputations:", round(sd(imputed_vals), 2), "\n")
      cat("  Range:", round(min(imputed_vals), 2), "-", round(max(imputed_vals), 2), "\n")
    }
  }
}
cat("\n")

# ============================================
# STEP 7: CALCULATE HT SCORES FOR EACH IMPUTED DATASET
# ============================================

cat("=== STEP 7: CALCULATING HT SCORES ===\n\n")

# Function to process a single imputed dataset
process_imputed_dataset <- function(imp_data, df_main) {

  # Merge with main data
  df_analysis <- df_main %>%
    inner_join(
      imp_data %>% select(record_id, plt_d14, anc_d14, hb_d14, crp_d14, ferritin_d14),
      by = "record_id"
    ) %>%
    mutate(
      ht_score_d14 = mapply(calculate_hematotox, plt_d14, anc_d14, hb_d14, crp_d14, ferritin_d14),
      d14_group = case_when(
        ht_score_d14 >= 3 ~ "high",
        ht_score_d14 <= 2 ~ "low",
        TRUE ~ NA_character_
      ),
      trajectory_group = paste0(baseline_group, "/", d14_group)
    ) %>%
    filter(!is.na(ht_score_d14)) %>%
    mutate(
      pfs_d = pmin(pfs_d, 1095, na.rm = TRUE),
      os_d = pmin(os_d, 1095, na.rm = TRUE),
      trajectory_group = factor(trajectory_group,
                                levels = c("high/high", "low/high", "high/low", "low/low"))
    )

  return(df_analysis)
}

# Process all imputed datasets
imputed_analyses <- lapply(1:imp_nr, function(i) {
  imp_data <- complete(imputedDataSets, i)
  process_imputed_dataset(imp_data, df_main)
})

cat("HT scores calculated for all", imp_nr, "imputed datasets\n\n")

# Check trajectory group distributions across imputations
cat("Trajectory group distributions across imputed datasets:\n")
for (i in 1:imp_nr) {
  cat("Dataset", i, ":", table(imputed_analyses[[i]]$trajectory_group), "\n")
}
cat("\n")

# ============================================
# STEP 8: POOLED SURVIVAL ANALYSIS
# ============================================

cat("=== STEP 8: POOLED SURVIVAL ANALYSIS ===\n\n")

# Function to run Cox regression on one imputed dataset
run_cox_analysis <- function(df, outcome = "pfs") {
  if (outcome == "pfs") {
    cox_model <- coxph(Surv(pfs_d, pfs_ev) ~ trajectory_group, data = df)
  } else {
    cox_model <- coxph(Surv(os_d, os_ev) ~ trajectory_group, data = df)
  }
  return(cox_model)
}

# Run Cox regression on all imputed datasets
cat("Running Cox regression on all imputed datasets...\n\n")

# PFS analysis
pfs_models <- lapply(imputed_analyses, function(df) run_cox_analysis(df, "pfs"))

# OS analysis
os_models <- lapply(imputed_analyses, function(df) run_cox_analysis(df, "os"))

# Pool results using Rubin's rules
# Extract coefficients and variances from each model
pool_cox_results <- function(models) {

  # Get coefficient names from first model
  coef_names <- names(coef(models[[1]]))
  m <- length(models)

  # Extract coefficients and variances
  coefs <- sapply(models, coef)
  vars <- sapply(models, function(m) diag(vcov(m)))

  if (is.vector(coefs)) {
    coefs <- matrix(coefs, nrow = 1)
    vars <- matrix(vars, nrow = 1)
  }

  # Rubin's rules
  # Q_bar = mean of estimates
  Q_bar <- rowMeans(coefs)

  # U_bar = mean of within-imputation variances
  U_bar <- rowMeans(vars)

  # B = between-imputation variance
  B <- apply(coefs, 1, var)

  # Total variance
  T_var <- U_bar + (1 + 1/m) * B

  # Standard errors
  se <- sqrt(T_var)

  # Hazard ratios and CIs
  hr <- exp(Q_bar)
  hr_lower <- exp(Q_bar - 1.96 * se)
  hr_upper <- exp(Q_bar + 1.96 * se)

  # P-values (Wald test)
  z <- Q_bar / se
  p_value <- 2 * pnorm(-abs(z))

  results <- data.frame(
    term = coef_names,
    estimate = hr,
    conf.low = hr_lower,
    conf.high = hr_upper,
    std.error = se,
    p.value = p_value,
    row.names = NULL
  )

  return(results)
}

# Pool PFS results
pfs_pooled <- pool_cox_results(pfs_models)
cat("POOLED PFS RESULTS (Hazard Ratios, ref: high/high):\n")
print(pfs_pooled %>%
        mutate(across(c(estimate, conf.low, conf.high), ~round(., 3)),
               p.value = format_pval_short(p.value)))
cat("\n")

# Pool OS results
os_pooled <- pool_cox_results(os_models)
cat("POOLED OS RESULTS (Hazard Ratios, ref: high/high):\n")
print(os_pooled %>%
        mutate(across(c(estimate, conf.low, conf.high), ~round(., 3)),
               p.value = format_pval_short(p.value)))
cat("\n")

# ============================================
# STEP 9: SURVIVAL PLOTTING (using first imputed dataset as representative)
# ============================================

cat("=== STEP 9: CREATING SURVIVAL PLOTS ===\n\n")

# Use first imputed dataset for visualization
# (Alternative: could use majority vote or average survival curves)
df_analysis_plot <- imputed_analyses[[1]]

cat("Using first imputed dataset for visualization\n")
cat("Final cohort for survival analysis: n =", nrow(df_analysis_plot), "\n\n")

cat("Trajectory group distribution:\n")
print(table(df_analysis_plot$trajectory_group))
cat("\n")

# Plotting function (modified to use pooled results)
plot_survival_trajectory <- function(fit, df, outcome_name, time_var, event_var,
                                     file_prefix, pooled_results) {

  horizon_days <- 1095

  survdiff_res <- survdiff(Surv(df[[time_var]], df[[event_var]]) ~ trajectory_group, data = df)
  logrank_p <- 1 - pchisq(survdiff_res$chisq, df = length(survdiff_res$n) - 1)

  surv_tbl <- as.data.frame(summary(fit)$table)
  if (!("strata" %in% names(surv_tbl))) {
    surv_tbl <- surv_tbl %>% rownames_to_column("strata")
  }

  summary_stats <- surv_tbl %>%
    mutate(
      group = str_remove(strata, "^trajectory_group="),
      n = .data[[intersect(c("n", "records"), names(surv_tbl))[1]]],
      events = .data[[intersect(c("events", "event"), names(surv_tbl))[1]]],
      median_days = .data[[intersect(c("median", "Median"), names(surv_tbl))[1]]],
      lower = .data[[intersect(c("0.95LCL", "X0.95LCL", "lower"), names(surv_tbl))[1]]],
      upper = .data[[intersect(c("0.95UCL", "X0.95UCL", "upper"), names(surv_tbl))[1]]]
    ) %>%
    transmute(
      group,
      n = as.numeric(n),
      events = as.numeric(events),
      median_days = round(as.numeric(median_days), 1),
      lower = round(as.numeric(lower), 1),
      upper = round(as.numeric(upper), 1)
    )

  p <- ggsurvplot(
    fit,
    data = df,
    risk.table = TRUE,
    conf.int = FALSE,
    palette = colors_trajectory,
    legend.title = "Trajectory\n(Baseline/Day14)",
    legend.labs = levels(df$trajectory_group),
    title = paste0(outcome_name, " by CAR-HEMATOTOX Trajectory (MICE Imputed)"),
    subtitle = paste0("Baseline ht_rg + Day 14 HT Score | Multiple Imputation (m=", imp_nr, ")"),
    xlab = "Time (months)",
    ylab = paste0(outcome_name, " probability"),
    xlim = c(0, horizon_days),
    break.time.by = 180,
    surv.median.line = "hv",
    pval = FALSE,
    ggtheme = theme_minimal(base_size = 14),
    risk.table.height = 0.28,
    tables.theme = theme_cleantable()
  )

  median_text <- paste0(
    "Median ", outcome_name, ":\n",
    paste0(
      summary_stats$group, ": ",
      ifelse(is.na(summary_stats$median_days), "NR",
             paste0(round(summary_stats$median_days / 30.4, 1), " mo")),
      collapse = "\n"
    )
  )

  # Use pooled results for HR annotation
  hr_lines <- c("Pooled HRs (ref: high/high):")
  if (nrow(pooled_results) > 0) {
    for (i in seq_len(nrow(pooled_results))) {
      hr_lines <- c(hr_lines,
                    paste0(
                      "  ", str_remove(pooled_results$term[i], "^trajectory_group"),
                      ": ", round(pooled_results$estimate[i], 2),
                      " (", round(pooled_results$conf.low[i], 2), "-",
                      round(pooled_results$conf.high[i], 2), ")"
                    ))
    }
  }
  hr_text <- paste(hr_lines, collapse = "\n")

  annot_text <- paste0(
    "Log-rank p = ", format_pval_short(logrank_p), "\n\n",
    hr_text, "\n\n",
    median_text
  )

  p$plot <- p$plot +
    scale_x_continuous(
      breaks = seq(0, horizon_days, by = 180),
      labels = c("0", "6", "12", "18", "24", "30", "36")
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "gray40") +
    annotate(
      "text", x = horizon_days * 0.02, y = 0.08,
      label = annot_text, hjust = 0, vjust = 0, size = 2.6,
      color = "black", family = "mono"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = element_blank(),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    )

  p$table <- p$table + theme(axis.title.y = element_blank())

  combined_plot <- plot_grid(p$plot, p$table, ncol = 1, rel_heights = c(3, 1))

  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    file.path(output_path, paste0(file_prefix, "_km_plot.png")),
    plot = combined_plot,
    width = 10, height = 9, dpi = 300
  )

  cat("Plot saved:", file_prefix, "_km_plot.png\n")

  return(list(logrank_p = logrank_p, summary_stats = summary_stats))
}

# ============================================
# STEP 10: RUN SURVIVAL ANALYSES
# ============================================

cat("=== STEP 10: RUNNING SURVIVAL ANALYSES ===\n\n")

cat("Running PFS analysis...\n")
pfs_fit <- survfit(Surv(pfs_d, pfs_ev) ~ trajectory_group, data = df_analysis_plot)
pfs_results <- plot_survival_trajectory(pfs_fit, df_analysis_plot, "Progression-Free Survival",
                                         "pfs_d", "pfs_ev", "km_pfs_mice_imputed",
                                         pfs_pooled)

cat("\nRunning OS analysis...\n")
os_fit <- survfit(Surv(os_d, os_ev) ~ trajectory_group, data = df_analysis_plot)
os_results <- plot_survival_trajectory(os_fit, df_analysis_plot, "Overall Survival",
                                        "os_d", "os_ev", "km_os_mice_imputed",
                                        os_pooled)

# ============================================
# STEP 11: EXPORT FILES
# ============================================

cat("\n=== STEP 11: EXPORTING FILES ===\n\n")

# Export pooled results
pooled_export <- bind_rows(
  pfs_pooled %>% mutate(outcome = "PFS"),
  os_pooled %>% mutate(outcome = "OS")
)
write_csv2(pooled_export, file.path(output_path, "mice_pooled_cox_results.csv"))
cat("Exported: mice_pooled_cox_results.csv\n")

# Export imputation diagnostics
mice_summary <- data.frame(
  variable = imputedVars,
  n_imputed = sapply(imputedVars, function(v) sum(is.na(mice_data[[v]]))),
  method = "pmm",
  n_datasets = imp_nr
)
write_csv2(mice_summary, file.path(output_path, "mice_imputation_summary.csv"))
cat("Exported: mice_imputation_summary.csv\n")

# Export patient lists by group (from first imputation)
for (grp in levels(df_analysis_plot$trajectory_group)) {
  grp_data <- df_analysis_plot %>%
    filter(trajectory_group == grp) %>%
    select(record_id, baseline_group, ht_score_d14, d14_group, trajectory_group)

  grp_filename <- paste0("trajectory_mice_", gsub("/", "_", grp), "_patients.csv")
  write_csv2(grp_data, file.path(output_path, grp_filename))
  cat("Exported:", grp_filename, "- n =", nrow(grp_data), "\n")
}

# ============================================
# SUMMARY
# ============================================

cat("\n", strrep("=", 70), "\n")
cat("MICE IMPUTATION ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n\n")

cat("COHORT SUMMARY:\n")
cat("  Patients in main dataset:        ", nrow(df_main), "\n")
cat("  Patients with day 14 window data:", nrow(d14_values), "\n")
cat("  Patients in survival analysis:   ", nrow(df_analysis_plot), "\n\n")

cat("IMPUTATION SUMMARY:\n")
cat("  Method: MICE with Predictive Mean Matching (pmm)\n")
cat("  Number of imputed datasets: ", imp_nr, "\n")
cat("  Variables imputed: ", paste(imputedVars, collapse = ", "), "\n\n")

cat("STATISTICAL APPROACH:\n")
cat("  Cox regression run on each imputed dataset\n")
cat("  Results pooled using Rubin's rules\n")
cat("  Survival curves from representative dataset\n\n")

cat("OUTPUT FILES:\n")
cat("  - km_pfs_mice_imputed_km_plot.png\n")
cat("  - km_os_mice_imputed_km_plot.png\n")
cat("  - mice_pooled_cox_results.csv\n")
cat("  - mice_imputation_summary.csv\n")
cat("  - trajectory_mice_*_patients.csv (4 files)\n")
cat("\nOutput directory:", output_path, "\n")
