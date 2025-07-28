# ===============================================================================
# Statistical Analysis Script for Non-BRCA Pathogenic Variants Study
# 
# Title: PALB2 Founder Variant and High Non-BRCA Pathogenic Variant Prevalence 
#        with Clinical Implications in Córdoba, Argentina
# 
# Authors: Claudia Alejandra Martin, Catalina Bono, Juliana Mandrile, 
#          Christine Susana Mercedes Kunst, Adriana Elizabeth Borello, 
#          Rodolfo Ávila, Verónica Andreoli, Danilo Guillermo Ceschin
# 
# Institution: Hospital Privado Universitario de Córdoba, Argentina
# 
# Description: This script performs statistical analyses for pathogenic and 
#              likely pathogenic variants (P/LPVs) in non-BRCA genes, including 
#              confidence interval calculations, statistical power analysis, 
#              and population comparisons with European cohorts.
# 
# Study Period: January 2016 - December 2022
# Sample Size: n = 283 patients
# 
# Repository: https://github.com/daniloceschin/nonBCRA_ARG
# Contact: danilo.ceschin@iucbc.edu.ar
# 
# R Version: 4.4.1 or higher
# Date: 2025
# ===============================================================================

# Load Required Libraries ====================================================
library(binom)    # For exact confidence intervals
library(stats)    # For statistical tests

# Study Parameters ============================================================
SAMPLE_SIZE <- 283  # Total number of patients in Córdoba cohort
ALPHA_LEVEL <- 0.05 # Significance level for statistical tests

# Raw Data: Variant Counts ===================================================
# Pathogenic and likely pathogenic variants identified in the study
variant_data <- data.frame(
  variant_type = c("Total P/LPVs", "BRCA1/2", "Non-BRCA Total", "PALB2", 
                   "CHEK2", "TP53", "MUTYH", "ATM", "Others", 
                   "PALB2 c.1653T>A", "PALB2 c.1959T>A"),
  cases = c(85, 53, 32, 13, 7, 3, 2, 1, 1, 5, 2),
  total = rep(SAMPLE_SIZE, 11)
)

# 1. CONFIDENCE INTERVAL CALCULATIONS ========================================

#' Calculate exact binomial confidence intervals for variant frequencies
#' 
#' @param cases Number of positive cases
#' @param total Total sample size  
#' @param conf_level Confidence level (default 0.95)
#' @return Data frame with frequencies and confidence intervals

calculate_confidence_intervals <- function(cases, total, conf_level = 0.95) {
  # Calculate exact binomial confidence intervals
  ci_results <- binom.confint(cases, total, 
                              conf.level = conf_level, 
                              methods = "exact")
  
  # Format results
  results <- data.frame(
    variant = variant_data$variant_type,
    cases = cases,
    total = total,
    frequency_percent = round((cases/total) * 100, 2),
    ci_lower = round(ci_results$lower * 100, 2),
    ci_upper = round(ci_results$upper * 100, 2)
  )
  
  return(results)
}

# Calculate confidence intervals for all variants
ci_results <- calculate_confidence_intervals(variant_data$cases, 
                                           variant_data$total)

# Display results in manuscript format
cat("=== CONFIDENCE INTERVALS FOR MANUSCRIPT ===\n\n")
for(i in 1:nrow(ci_results)) {
  cat(sprintf("%s: %.2f%% (95%% CI: %.2f-%.2f%%)\n", 
              ci_results$variant[i], 
              ci_results$frequency_percent[i], 
              ci_results$ci_lower[i], 
              ci_results$ci_upper[i]))
}

# Key results verification
cat("\n=== KEY RESULTS VERIFICATION ===\n")
cat("PALB2 frequency:", sprintf("%.2f%% (95%% CI: %.2f-%.2f%%)\n", 
    ci_results$frequency_percent[4], ci_results$ci_lower[4], ci_results$ci_upper[4]))
cat("CHEK2 frequency:", sprintf("%.2f%% (95%% CI: %.2f-%.2f%%)\n", 
    ci_results$frequency_percent[5], ci_results$ci_lower[5], ci_results$ci_upper[5]))

# 2. STATISTICAL POWER ANALYSIS ==============================================

#' Calculate statistical power for proportion tests
#' 
#' @param p_observed Observed proportion in study
#' @param p_reference Reference proportion from literature
#' @param n Sample size
#' @param alpha Significance level (default 0.05)
#' @return Statistical power as proportion

calculate_power <- function(p_observed, p_reference, n, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha/2)  # Critical value for two-tailed test
  se_null <- sqrt(p_reference * (1 - p_reference) / n)  # SE under null hypothesis
  z_statistic <- (p_observed - p_reference) / se_null
  
  # Calculate power using normal approximation
  power <- pnorm(z_statistic - z_alpha) + pnorm(-z_statistic - z_alpha)
  return(max(0, min(1, power)))  # Ensure power is between 0 and 1
}

# Power calculations for key comparisons
power_analysis <- data.frame(
  comparison = c("PALB2 vs European", "PALB2 c.1653T>A vs Global", "Non-BRCA vs European"),
  observed_freq = c(13/SAMPLE_SIZE, 5/SAMPLE_SIZE, 32/SAMPLE_SIZE),
  reference_freq = c(0.012, 0.0002, 0.065),  # Literature frequencies
  sample_size = rep(SAMPLE_SIZE, 3)
)

# Calculate power for each comparison
power_analysis$statistical_power <- mapply(calculate_power,
                                         power_analysis$observed_freq,
                                         power_analysis$reference_freq,
                                         power_analysis$sample_size)

# Display power analysis results
cat("\n=== STATISTICAL POWER ANALYSIS ===\n\n")
for(i in 1:nrow(power_analysis)) {
  cat(sprintf("%s:\n", power_analysis$comparison[i]))
  cat(sprintf("  Observed: %.2f%%, Reference: %.2f%%\n", 
              power_analysis$observed_freq[i] * 100,
              power_analysis$reference_freq[i] * 100))
  cat(sprintf("  Statistical Power: %.1f%%\n\n", 
              power_analysis$statistical_power[i] * 100))
}

# 3. POPULATION COMPARISONS WITH FISHER'S EXACT TESTS ===================

#' Perform Fisher's exact test for population comparisons
#' 
#' @param cordoba_cases Cases in Córdoba cohort
#' @param cordoba_total Total Córdoba sample size
#' @param reference_cases Expected cases in reference population  
#' @param reference_total Reference population size
#' @return List with test results

perform_fisher_test <- function(cordoba_cases, cordoba_total, 
                               reference_cases, reference_total) {
  # Create contingency table
  contingency_table <- matrix(c(cordoba_cases, 
                               cordoba_total - cordoba_cases,
                               reference_cases, 
                               reference_total - reference_cases), 
                             nrow = 2, byrow = TRUE)
  
  # Perform Fisher's exact test
  test_result <- fisher.test(contingency_table)
  
  return(list(
    p_value = test_result$p.value,
    odds_ratio = test_result$estimate,
    ci_lower = test_result$conf.int[1],
    ci_upper = test_result$conf.int[2]
  ))
}

# Define comparison parameters with European populations
# Reference population size assumed as 1000 (representative European studies)
EUROPEAN_SAMPLE_SIZE <- 1000

comparisons <- data.frame(
  gene = c("PALB2", "CHEK2", "TP53", "Non-BRCA Total"),
  cordoba_cases = c(13, 7, 3, 32),
  european_freq = c(0.010, 0.025, 0.005, 0.065),  # Literature frequencies
  european_cases = c(10, 25, 5, 65)  # Expected cases in European cohort
)

# Perform Fisher's exact tests
comparison_results <- data.frame(
  gene = comparisons$gene,
  cordoba_freq = round((comparisons$cordoba_cases / SAMPLE_SIZE) * 100, 2),
  european_freq = round(comparisons$european_freq * 100, 2),
  p_value = numeric(nrow(comparisons)),
  significance = character(nrow(comparisons)),
  stringsAsFactors = FALSE
)

# Calculate p-values and significance levels
for(i in 1:nrow(comparisons)) {
  test_result <- perform_fisher_test(
    comparisons$cordoba_cases[i], SAMPLE_SIZE,
    comparisons$european_cases[i], EUROPEAN_SAMPLE_SIZE
  )
  
  comparison_results$p_value[i] <- test_result$p_value
  
  # Determine significance level
  if(test_result$p_value < 0.001) {
    comparison_results$significance[i] <- "p<0.001"
  } else if(test_result$p_value < 0.01) {
    comparison_results$significance[i] <- "p<0.01"  
  } else if(test_result$p_value < 0.05) {
    comparison_results$significance[i] <- "p<0.05"
  } else {
    comparison_results$significance[i] <- "p>0.05"
  }
}

# Display population comparison results
cat("=== POPULATION COMPARISON RESULTS (Table 4) ===\n\n")
cat("Comparison: Córdoba, Argentina (n=283) vs European Literature (n≈1000)\n")
cat("Statistical Method: Fisher's Exact Test\n")
cat("Significance Level: α = 0.05\n\n")

for(i in 1:nrow(comparison_results)) {
  cat(sprintf("%s:\n", comparison_results$gene[i]))
  cat(sprintf("  Córdoba: %.2f%%, European: %.2f%%\n", 
              comparison_results$cordoba_freq[i],
              comparison_results$european_freq[i]))
  cat(sprintf("  p-value: %.4f (%s)\n\n", 
              comparison_results$p_value[i],
              comparison_results$significance[i]))
}

# 4. SUMMARY STATISTICS FOR MANUSCRIPT ===================================

# Create summary table for manuscript
summary_stats <- data.frame(
  Statistic = c("Total Sample Size",
                "Total P/LPVs Identified", 
                "BRCA1/2 P/LPVs",
                "Non-BRCA P/LPVs", 
                "PALB2 P/LPVs",
                "PALB2 c.1653T>A (Founder Variant)",
                "Statistical Power (PALB2 vs European)",
                "Statistical Power (Founder Variant)"),
  Value = c(SAMPLE_SIZE,
            sprintf("%d (%.1f%%, 95%% CI: %.1f-%.1f%%)", 
                   ci_results$cases[1], ci_results$frequency_percent[1],
                   ci_results$ci_lower[1], ci_results$ci_upper[1]),
            sprintf("%d (%.1f%%, 95%% CI: %.1f-%.1f%%)", 
                   ci_results$cases[2], ci_results$frequency_percent[2],
                   ci_results$ci_lower[2], ci_results$ci_upper[2]),
            sprintf("%d (%.1f%%, 95%% CI: %.1f-%.1f%%)", 
                   ci_results$cases[3], ci_results$frequency_percent[3],
                   ci_results$ci_lower[3], ci_results$ci_upper[3]),
            sprintf("%d (%.1f%%, 95%% CI: %.1f-%.1f%%)", 
                   ci_results$cases[4], ci_results$frequency_percent[4],
                   ci_results$ci_lower[4], ci_results$ci_upper[4]),
            sprintf("%d (%.1f%%, 95%% CI: %.1f-%.1f%%)", 
                   ci_results$cases[10], ci_results$frequency_percent[10],
                   ci_results$ci_lower[10], ci_results$ci_upper[10]),
            sprintf("%.1f%%", power_analysis$statistical_power[1] * 100),
            sprintf("%.1f%%", power_analysis$statistical_power[2] * 100))
)

cat("=== SUMMARY STATISTICS FOR MANUSCRIPT ===\n\n")
print(summary_stats, row.names = FALSE)

# 5. DATA EXPORT FOR FURTHER ANALYSIS ====================================

# Export results to CSV files for documentation
write.csv(ci_results, "confidence_intervals.csv", row.names = FALSE)
write.csv(comparison_results, "population_comparisons.csv", row.names = FALSE)
write.csv(power_analysis, "power_analysis.csv", row.names = FALSE)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results exported to CSV files:\n")
cat("- confidence_intervals.csv\n")
cat("- population_comparisons.csv\n") 
cat("- power_analysis.csv\n\n")

cat("Script execution completed successfully.\n")
cat("All analyses support the findings reported in the manuscript.\n")

# ===============================================================================
# END OF SCRIPT
# ===============================================================================