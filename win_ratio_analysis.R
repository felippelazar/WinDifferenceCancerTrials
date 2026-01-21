# ============================================================================ #
# Win-Ratio as an Alternative Outcome
# Author: Felippe Lazar, Universidade de Sao Paulo, 2024 #
# ============================================================================ #

# ============================================================================ #
# Loading Libraries

library(IPDfromKM)
library(gtsummary)
library(tidyverse)
library(readxl)
library(broom)
library(ggsurvfit)
library(survRM2)
library(rio)
library(survival)
library(BuyseTest)
library(patchwork)
source('tryRetry.R')
source('win_ratio_functions.R')

# Loading the Metadata
metadata <- rio::import("data/harvard_misc_dataset/PhaseIII_ClinicalTrials_metadata_F.xlsx") %>%
  janitor::clean_names()

dim(metadata) # 161  18

# Checking if all have at least one figure
sum(!is.na(metadata$os_figure) | !is.na(metadata$surrogate_event_free_survival_figure)) # 161 

metadata_long <- metadata %>%
  tidyr::pivot_longer(cols = matches("figure"), names_to = "figure_outcome_ipd", values_to = "figure_ipd")

surrogate_trials <- metadata_long %>%
  dplyr::filter(figure_outcome_ipd == "surrogate_event_free_survival_figure") %>%
  dplyr::filter(!is.na(figure_ipd)) %>%
  dplyr::select(trial_name, pub_med_id, trial_endpoint = surrogate_event_free_survival_term, 
                trial_metastatic = includes_cancers_with_distant_metastases_no_0_yes_1,
                trial_scale = surrogate_event_free_survival_scale_days_d_weeks_w_months_m_years_y,
                trial_endpoint_superior = superior_surrogate_event_free_survival_no_0_yes_1, 
                figure_outcome_ipd, figure_ipd) %>%
  dplyr::mutate(trial_filename = paste0(trial_name, "_", figure_ipd))

nrow(surrogate_trials) # 146

survival_trials <- metadata_long %>%
  dplyr::filter(figure_outcome_ipd == "os_figure") %>%
  dplyr::filter(!is.na(figure_ipd)) %>%
  dplyr::mutate(trial_endpoint = "OS") %>%
  dplyr::select(trial_name, pub_med_id, trial_endpoint,
                trial_metastatic = includes_cancers_with_distant_metastases_no_0_yes_1,
                trial_scale = os_scale_days_d_weeks_w_months_m_years_y, 
                trial_endpoint_superior = superior_os_no_0_yes_1,
                figure_outcome_ipd, figure_ipd) %>%
  dplyr::mutate(trial_filename = paste0(trial_name, "_", figure_ipd))

nrow(survival_trials) # 116

metadata_tidied <- bind_rows(surrogate_trials, survival_trials)

# Loading the IPD Data from Harvard Dataset
ipd_files <- list.files("data/harvard_misc_dataset/", pattern = ".*csv", full.names = TRUE)

# Joining with the Metadata and Transforming the Outcome Time in Weeks
ipd_data_list <- lapply(ipd_files, function(x){
  
  data <- rio::import(x)
  data <- data %>%
    dplyr::mutate(trial_filename = stringr::str_extract(x, "\\/\\/(.*?)[.]csv$", group = 1)) %>%
    dplyr::left_join(metadata_tidied) %>%
    dplyr::mutate(
      TimeWeeks = case_when(
        trial_scale == "W" ~ Time,
        trial_scale == "M" ~ Time * 4.345,   # average weeks per month
        trial_scale == "Y" ~ Time * 52.1775, # average weeks per year
        TRUE ~ NA_real_
      )
    )
  
  if(!is.na(first(data$trial_name))) return(data)
  
})


# Getting Only Survival IPD
ipd_data_list <- purrr::compact(ipd_data_list)
names_files <- sapply(ipd_data_list, function(x) first(x$trial_filename))
names(ipd_data_list) <- names_files

ipd_survival_estimates <- pbapply::pblapply(ipd_data_list, getSurvivalEstimates)
names(ipd_survival_estimates) <- names_files
saveRDS(ipd_survival_estimates, "data/ipd_survival_estimates.rds")
ipd_survival_estimates <- readRDS("data/ipd_survival_estimates.rds")

ipd_winratio_estimates <- pbapply::pblapply(ipd_data_list, function(x) getWinRatioBatchManual(x) %try% NULL)
names(ipd_winratio_estimates) <- names_files
saveRDS(ipd_winratio_estimates, "data/ipd_winratio_estimates.rds")
ipd_winratio_estimates <- readRDS("data/ipd_winratio_estimates.rds")

ipd_peron_estimates <- pbapply::pblapply(ipd_data_list, function(x) getBuyseTestPeron(x) %try% NULL)
names(ipd_peron_estimates) <- names_files
saveRDS(ipd_peron_estimates, "data/ipd_peron_estimates.rds")
ipd_peron_estimates <- readRDS("data/ipd_peron_estimates.rds")

ipd_gehan_estimates <- pbapply::pblapply(ipd_data_list, function(x) getBuyseTestGehan(x) %try% NULL)
names(ipd_gehan_estimates) <- names_files
saveRDS(ipd_gehan_estimates, "data/ipd_gehan_estimates.rds")
ipd_gehan_estimates <- readRDS("data/ipd_gehan_estimates.rds")

# Combining all Results
survival_estimates_df <- do.call(bind_rows, ipd_survival_estimates)
ipd_winratio_estimates_df <- do.call(bind_rows, ipd_winratio_estimates)
# ipd_gehan_estimates_df <- do.call(bind_rows, ipd_gehan_estimates)
ipd_peron_estimates_df <- do.call(bind_rows, ipd_peron_estimates)

# # Getting Only Part of the IPD Gehan and Peron Estimates
# ipd_gehan_estimates_wide_df <- ipd_gehan_estimates_df %>%
#   dplyr::filter(gehan_min_time %in% c(0, 4, 12, 26, 52, 104)) %>%
#   tidyr::pivot_wider(
#     id_cols = c(trial_name, trial_filename),
#     names_from = gehan_min_time,
#     values_from = matches('gehan_'),
#     names_glue = "{.value}_{gehan_min_time}wk"
#   )

ipd_peron_estimates_wide_df <- ipd_peron_estimates_df %>%
  dplyr::filter(peron_min_time %in% c(0, 4, 12, 26, 52, 104)) %>%
  tidyr::pivot_wider(
    id_cols = c(trial_name, trial_filename),
    names_from = peron_min_time,
    values_from = matches('peron_'),
    names_glue = "{.value}_{peron_min_time}wk"
  )

# Now Combining All of Them
final_results_df <- survival_estimates_df %>%
  # dplyr::left_join(ipd_gehan_estimates_wide_df, by = c("trial_name", "trial_filename")) %>%
  dplyr::left_join(ipd_peron_estimates_wide_df, by = c("trial_name", "trial_filename")) %>%
  dplyr::left_join(metadata_tidied,
                   by = c("trial_name", "trial_filename"))


# Define the week thresholds
week_thresholds <- c(0, 4, 12, 26, 52, 104)

# Create the wide dataset
ipd_peron_estimates_wide_df <- ipd_peron_estimates_df %>%
      dplyr::filter(peron_min_time %in% week_thresholds) %>%
      tidyr::pivot_wider(
            id_cols = c(trial_name, trial_filename),
            names_from = peron_min_time,
            values_from = matches('peron_'),
            names_glue = "{.value}_{peron_min_time}wk"
      )

# Combine all datasets
final_results_df <- survival_estimates_df %>%
      dplyr::left_join(ipd_peron_estimates_wide_df, by = c("trial_name", "trial_filename")) %>%
      dplyr::left_join(metadata_tidied, by = c("trial_name", "trial_filename"))

# Set Peron variables to NA if median survival < weeks threshold
final_results_df <- final_results_df %>%
      dplyr::mutate(
            across(
                  matches("^peron_.*_4wk$"),
                  ~dplyr::if_else(median_fup_estimate < 4, NA_real_, .)
            ),
            across(
                  matches("^peron_.*_12wk$"),
                  ~dplyr::if_else(median_fup_estimate < 12, NA_real_, .)
            ),
            across(
                  matches("^peron_.*_26wk$"),
                  ~dplyr::if_else(median_fup_estimate < 26, NA_real_, .)
            ),
            across(
                  matches("^peron_.*_52wk$"),
                  ~dplyr::if_else(median_fup_estimate < 52, NA_real_, .)
            ),
            across(
                  matches("^peron_.*_104wk$"),
                  ~dplyr::if_else(median_fup_estimate < 104, NA_real_, .)
            )
      )

colnames(final_results_df)

# Getting Only those with OS Results
final_results_os_df <- final_results_df %>%
  dplyr::filter(trial_endpoint == 'OS') %>%
  dplyr::filter(trial_metastatic == 1)

# ========================================================================================================
# Creating Visualizations and Validation of Values

order_trials <- final_results_os_df %>%
  arrange(desc(cox_estimate)) %>%
  pull(trial_filename) %>%
  unique()

metrics_long <- final_results_os_df %>%
  dplyr::mutate(trial_filename = factor(trial_filename, levels = order_trials)) %>%
  #dplyr::mutate(across(matches('gehan_nb_estimate'), ~ .*100)) %>%
  dplyr::mutate(across(matches('peron_nb_estimate'), ~ .*100)) %>%
  dplyr::mutate(across(matches('landmark_diff'), ~ .*100)) %>%
  tidyr::pivot_longer(
    cols = c(matches('peron_nb_estimate'), # matches('gehan_nb_estimate'), 
             matches('landmark_diff'), 'median_diff', 'max_rmean_diff', 'cox_estimate'),
    names_to = 'metric_name',
    values_to = 'metric_value'
  ) %>%
  dplyr::select(trial_name, trial_filename, metric_name, metric_value) %>%
  dplyr::mutate(metric_name = factor(metric_name, levels = c("cox_estimate", "median_diff", "max_rmean_diff", 
                                                             "peron_nb_estimate_0wk", "peron_nb_estimate_4wk", "peron_nb_estimate_12wk", "peron_nb_estimate_26wk", "peron_nb_estimate_52wk", "peron_nb_estimate_104wk",
                                                             "gehan_nb_estimate_0wk", "gehan_nb_estimate_4wk", "gehan_nb_estimate_12wk", "gehan_nb_estimate_26wk", "gehan_nb_estimate_52wk", "gehan_nb_estimate_104wk",
                                                             "landmark_diff_6mo", "landmark_diff_12mo", "landmark_diff_18mo", "landmark_diff_24mo", "landmark_diff_30mo", "landmark_diff_36mo"))) %>%
  group_by(metric_name) %>%
  mutate(metric_rank = ifelse(metric_name == 'cox_estimate', min_rank(metric_value), min_rank(desc(metric_value)))) %>% 
  ungroup()

metrics_long %>%
  dplyr::distinct(trial_filename, metric_name, .keep_all = TRUE) %>%
  ggplot(aes(x = metric_name, y = trial_filename, label = round(metric_value, 2))) +
  geom_tile(aes(fill = metric_rank), color = 'gray50', alpha = 0.5) +
  scale_fill_distiller(palette = "Spectral", direction = 1, name = "Rank") +
  geom_text(fontface = "bold", size = 7) + 
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0))

ggsave(plot = last_plot(), 'figures/os_metrics_heatmap_value.pdf', width = 52, height = 72, unit = 'cm')

metrics_long %>%
  dplyr::distinct(trial_filename, metric_name, .keep_all = TRUE) %>%
  ggplot(aes(x = metric_name, y = trial_filename, label = round(metric_rank, 2))) +
  geom_tile(aes(fill = metric_rank), color = 'gray50', alpha = 0.5) +
  scale_fill_distiller(palette = "Spectral", direction = 1, name = "Rank") +
  geom_text(fontface = "bold", size = 7) + 
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0))

ggsave(plot = last_plot(), 'figures/os_metrics_heatmap_rank.pdf', width = 52, height = 72, unit = 'cm')

# ========================================================================================================
# Creating Scatter Plot of Values - NB Benefit Against the Others

scatter_data_value <- metrics_long %>%
  rename(metric_x = metric_name, value_x = metric_value) %>%
  inner_join(
    metrics_long %>% rename(metric_y = metric_name, value_y = metric_value),
    by = c("trial_filename", "trial_name"), 
    relationship = "many-to-many" # Required in newer dplyr versions for this expansion
  ) %>%
  filter(
    grepl("nb_estimate", metric_y),
    !grepl("nb_estimate", metric_x)
  )

cor_data <- scatter_data_value %>%
  group_by(metric_x, metric_y) %>%
  summarize(
    # Calculate R, handling potential NAs
    pearson_r = cor(value_x, value_y, method = "pearson", use = "complete.obs"),
    # Create the label string
    label = paste0("R = ", round(pearson_r, 2)),
    .groups = "drop"
  )

scatter_data_value %>%
  ggplot(aes(x = value_x, y = value_y)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.5) +
  geom_text(
    data = cor_data,
    aes(label = label),
    x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, 
    size = 3, fontface = "bold", inherit.aes = FALSE
  ) +
  
  facet_grid(metric_y ~ metric_x, scales = "free") +
  theme_bw() +
  labs(x = "Comparison Metrics", y = "NB Statistics") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave('figures/scatter_metrics_nb_value.pdf', width = 30, height = 42, unit = 'cm')

cor_data %>%
  ggplot(aes(x = metric_x, y = metric_y, fill = abs(pearson_r), label = round(pearson_r, 2))) +
  geom_tile() + 
  geom_text() + 
  scale_fill_distiller(palette = "Spectral", direction = -1) + 
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0))

ggsave('figures/corplot_metrics_nb_value.pdf', width = 21, height = 18, unit = 'cm')

# ========================================================================================================
# Creating Scatter Plot of Ranks - NB Benefit Against the Others

scatter_data_rank <- metrics_long %>%
  rename(metric_x = metric_name, value_x = metric_rank) %>%
  inner_join(
    metrics_long %>% rename(metric_y = metric_name, value_y = metric_rank),
    by = c("trial_filename", "trial_name"), 
    relationship = "many-to-many" # Required in newer dplyr versions for this expansion
  ) %>%
  filter(
    grepl("nb_estimate", metric_y),
    !grepl("nb_estimate", metric_x)
  )

cor_data <- scatter_data_value %>%
  group_by(metric_x, metric_y) %>%
  summarize(
    # Calculate R, handling potential NAs
    spearman_r = cor(value_x, value_y, method = "spearman", use = "complete.obs"),
    # Create the label string
    label = paste0("R = ", round(spearman_r, 2)),
    .groups = "drop"
  )

scatter_data_rank %>%
  ggplot(aes(x = value_x, y = value_y)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.5) +
  geom_text(
    data = cor_data,
    aes(label = label),
    x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, 
    size = 3, fontface = "bold", inherit.aes = FALSE
  ) +
  
  facet_grid(metric_y ~ metric_x, scales = "free") +
  theme_bw() +
  labs(x = "Comparison Metrics", y = "NB Statistics") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave('figures/scatter_metrics_nb_rank.pdf', width = 30, height = 42, unit = 'cm')

cor_data %>%
  ggplot(aes(x = metric_x, y = metric_y, fill = abs(spearman_r), label = round(spearman_r, 2))) +
  geom_tile() + 
  geom_text() + 
  scale_fill_distiller(palette = "Spectral", direction = -1) + 
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0))

ggsave('figures/corplot_metrics_nb_rank.pdf', width = 21, height = 18, unit = 'cm')

# ========================================================================================================
# Creating Scatter Plot of Values for the Traditional Methods

scatter_data_value <- metrics_long %>%
  rename(metric_x = metric_name, value_x = metric_value) %>%
  inner_join(
    metrics_long %>% rename(metric_y = metric_name, value_y = metric_value),
    by = c("trial_filename", "trial_name"), 
    relationship = "many-to-many" # Required in newer dplyr versions for this expansion
  ) %>%
  filter(
    !grepl("nb_estimate", metric_y),
    !grepl("nb_estimate", metric_x)
  )

cor_data <- scatter_data_value %>%
  group_by(metric_x, metric_y) %>%
  summarize(
    # Calculate R, handling potential NAs
    pearson_r = cor(value_x, value_y, method = "pearson", use = "complete.obs"),
    # Create the label string
    label = paste0("R = ", round(pearson_r, 2)),
    .groups = "drop"
  )

scatter_data_value %>%
  ggplot(aes(x = value_x, y = value_y)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.5) +
  geom_text(
    data = cor_data,
    aes(label = label),
    x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, 
    size = 3, fontface = "bold", inherit.aes = FALSE
  ) +
  
  facet_grid(metric_y ~ metric_x, scales = "free") +
  theme_bw() +
  labs(x = "Comparison Metrics", y = "NB Statistics") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave('figures/scatter_metrics_standard_value.pdf', width = 30, height = 42, unit = 'cm')

cor_data %>%
  ggplot(aes(x = metric_x, y = metric_y, fill = abs(pearson_r), label = round(pearson_r, 2))) +
  geom_tile() + 
  geom_text() + 
  scale_fill_distiller(palette = "Spectral", direction = -1) + 
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0))

ggsave('figures/corplot_metrics_standard_value.pdf', width = 21, height = 18, unit = 'cm')

# ========================================================================================================
# Creating Scatter Plot of Ranks for the Traditional Methods

scatter_data_rank <- metrics_long %>%
  rename(metric_x = metric_name, value_x = metric_rank) %>%
  inner_join(
    metrics_long %>% rename(metric_y = metric_name, value_y = metric_rank),
    by = c("trial_filename", "trial_name"), 
    relationship = "many-to-many" # Required in newer dplyr versions for this expansion
  ) %>%
  filter(
    !grepl("nb_estimate", metric_y),
    !grepl("nb_estimate", metric_x)
  )

cor_data <- scatter_data_value %>%
  group_by(metric_x, metric_y) %>%
  summarize(
    # Calculate R, handling potential NAs
    spearman_r = cor(value_x, value_y, method = "spearman", use = "complete.obs"),
    # Create the label string
    label = paste0("R = ", round(spearman_r, 2)),
    .groups = "drop"
  )

scatter_data_rank %>%
  ggplot(aes(x = value_x, y = value_y)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.5) +
  geom_text(
    data = cor_data,
    aes(label = label),
    x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, 
    size = 3, fontface = "bold", inherit.aes = FALSE
  ) +
  
  facet_grid(metric_y ~ metric_x, scales = "free") +
  theme_bw() +
  labs(x = "Comparison Metrics", y = "NB Statistics") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave('figures/scatter_metrics_standard_rank.pdf', width = 30, height = 42, unit = 'cm')

cor_data %>%
  ggplot(aes(x = metric_x, y = metric_y, fill = abs(spearman_r), label = round(spearman_r, 2))) +
  geom_tile() + 
  geom_text() + 
  scale_fill_distiller(palette = "Spectral", direction = -1) + 
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0))

ggsave('figures/corplot_metrics_standard_rank.pdf', width = 21, height = 18, unit = 'cm')

# ========================================================================================================
# Creating Individual Figures for Each Trial
names(ipd_data_list) <- sapply(ipd_data_list, function(x) first(x$trial_filename))

# Selecting Trials to Save
trial_filenames <- final_results_df %>%
  dplyr::arrange(cox_estimate) %>%
  dplyr::filter(trial_endpoint == 'OS') %>%
  pull(trial_filename)

v_thresholds <- c(26, 52, 104)
milestones <- c(26, 52, 104)

semantic_colors <- c(
  "percWin"              = "#00E5FF", 
  "percLoss"             = "#3F51B5", 
  "percTrueTieEvent"     = "#7F8C8D", 
  "percTrueTieCensored"  = "#BDC3C7", 
  "percFalseTieEvent"    = "#F39C12", 
  "percFalseTieCensored" = "#F1C40F"  
)

i = 0

for(trial_filename_it in trial_filenames){
  
  {
    
    i = i + 1
    pmid <- metadata_tidied$pub_med_id[metadata_tidied$trial_filename == trial_filename_it]
    
    cox_graph <- survfit2(Surv(TimeWeeks, Event) ~ Arm, data = ipd_data_list[[trial_filename_it]]) %>%
      ggsurvfit() +
      add_risktable_strata_symbol(symbol = "\U25CF", size = 10) +
      add_risktable(
        risktable_height = 0.20,
        size = 4, # increase font size of risk table statistics
        risktable_stats = 'n.risk',
        theme = list(
          theme_risktable_default()
        )
      ) + 
      add_censor_mark() +
      add_quantile() + 
      labs(x = 'Weeks') + 
      add_pvalue(location  = "caption") +
      labs(title = glue::glue('Overall Survival - {trial_filename_it}')) +
      scale_color_manual(values = c("#3F51B5", "#00E5FF"))
    
    cox_graph <- ggsurvfit_build(cox_graph)
    
    censoring_graph <- survfit2(Surv(TimeWeeks, Event == 0) ~ Arm, data = ipd_data_list[[trial_filename_it]]) %>%
      ggsurvfit() +
      add_risktable_strata_symbol(symbol = "\U25CF", size = 10) +
      add_risktable(
        risktable_height = 0.20,
        size = 4, # increase font size of risk table statistics
        risktable_stats = 'n.risk',
        theme = list(
          theme_risktable_default()
        )
      ) + 
      add_censor_mark() +
      add_quantile() + 
      labs(x = 'Weeks') + 
      add_pvalue(location = "caption") +
      labs(title = glue::glue('Inverse-KM - PMID: {pmid}')) + 
      scale_color_manual(values = c("#3F51B5", "#00E5FF"))
    
    censoring_graph <- ggsurvfit_build(censoring_graph)
    
    # tile_statistics <- metrics_long %>%
    #   dplyr::filter(trial_filename == trial_filename_it) %>%
    #   dplyr::mutate(fake_y = 'metric_value') %>%
    #   ggplot(aes(x = metric_name, y = fake_y, label = round(metric_value, 2))) +
    #   geom_tile(fill = 'white', color = 'gray50') +
    #   geom_text(size = 5) +
    #   theme(axis.text.y = element_blank(),
    #         plot.background = element_blank(),
    #         panel.background = element_blank()) + 
    #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    #   scale_x_discrete(position = "top")
    
    gehan_winratio <- ipd_gehan_estimates[[trial_filename_it]] %>%
      dplyr::mutate(
        NetBenefit = gehan_nb_estimate,
        NegLoss = -(gehan_total_loss)/(gehan_total_wins + gehan_total_loss + gehan_total_ties),
        percWin = gehan_total_wins/(gehan_total_wins + gehan_total_loss + gehan_total_ties),
        statsSign = ifelse(gehan_wr_p.value < 0.05, "Significant", "Not Significant")
      )
    
    gehan_net_benefit_graph <- gehan_winratio %>%
      ggplot(aes(x = gehan_min_time)) +
      geom_col(aes(y = percWin, fill = "Favors Experimental Arm"), alpha = 0.3) +
      geom_col(aes(y = NegLoss, fill = "Favors Standard of Care"), alpha = 0.3) +
      geom_hline(yintercept = 0, color = "gray40", size = 0.5) +
      geom_col(aes(y = NetBenefit, fill = statsSign), alpha = 0.9) +
      geom_point(aes(y = NetBenefit, fill = "Net Benefit (Surplus)"), alpha = 0.9, size = 0.7) +
      geom_line(aes(y = NetBenefit, fill = "Net Benefit (Surplus)"), alpha = 0.9) +
      geom_vline(xintercept = c(26, 52, 104), color = "gray60", size = 0.5, linetype = 2) +
      ggrepel::geom_label_repel(
        data = gehan_winratio %>% filter(gehan_min_time %in% milestones) %>%
          mutate(nudge_amount = ifelse(NetBenefit >= 0, 0.03, -0.03)),
        aes(
          x = gehan_min_time, 
          y = NetBenefit, 
          label = scales::percent(NetBenefit, accuracy = 0.1)
        ),
        fontface = "bold", 
        size = 3.5, 
        color = "#006064",
        direction = "y",
        nudge_y = gehan_winratio %>% filter(gehan_min_time %in% milestones) %>%
          mutate(nudge_amount = ifelse(NetBenefit >= 0, 0.01, -0.01)) %>%
          pull(nudge_amount),
        segment.color = "#006064",
        segment.size = 0.4,
        min.segment.length = 0
      ) +
      # SCALES
      coord_cartesian(ylim = c(-0.3, 0.3)) +
      scale_y_continuous(labels = scales::percent, breaks = c(-0.4, -0.3, -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.3, 0.4)) +
      scale_fill_manual(
        values = c(
          "Similar Outcome (Tie)"       = "#CFD8DC", # Light Blue-Grey
          "Favors Experimental Arm"     = "#00E5FF", # Cyan (Your preference)
          "Favors Standard of Care"     = "#5C6BC0", # Soft Indigo (Instead of Scary Red)
          "Net Benefit (Surplus)"       = "#006064",  # Deep Teal/Cyan (Solid)
          "Significant"       = "#006064",
          "Not Significant" = "gray75"
        ),
        breaks = c("Favors Experimental Arm", "Net Benefit (Surplus)", 
                   "Favors Standard of Care", "Similar Outcome (Tie)")
      ) +
      # SCALES
      theme_minimal(base_size = 10) +
      labs(
        x = "Survival Time Threshold (Weeks)",
        y = "NTB Gehan (%)"
      ) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = "none"
      ) + 
      # We place them at the top of the plot (e.g., 0.45 or 45%)
      annotate(
        "text", 
        x = v_thresholds + 5, # nudge slightly right
        y = 0.45, 
        label = paste(v_thresholds, "wks"),
        angle = 90, size = 3, color = "gray60"
      ) +
      # Ensure the labels aren't cut off
      coord_cartesian(clip = "off") +
      theme(plot.margin = margin(5, 40, 5, 5))
    
    peron_winratio <- ipd_peron_estimates[[trial_filename_it]] %>%
      dplyr::mutate(
        NetBenefit = peron_nb_estimate,
        NegLoss = -(peron_total_loss)/(peron_total_wins + peron_total_loss + peron_total_ties),
        percWin = peron_total_wins/(peron_total_wins + peron_total_loss + peron_total_ties),
        statsSign = ifelse(peron_wr_p.value < 0.05, "Significant", "Not Significant")
      )
    
    peron_net_benefit_graph <- peron_winratio %>%
      ggplot(aes(x = peron_min_time)) +
      geom_col(aes(y = percWin, fill = "Favors Experimental Arm"), alpha = 0.3) +
      geom_col(aes(y = NegLoss, fill = "Favors Standard of Care"), alpha = 0.3) +
      geom_hline(yintercept = 0, color = "gray40", size = 0.5) +
      geom_col(aes(y = NetBenefit, fill = statsSign), alpha = 0.9) +
      geom_point(aes(y = NetBenefit, fill = "Net Benefit (Surplus)"), alpha = 0.9, size = 0.7) +
      geom_line(aes(y = NetBenefit, fill = "Net Benefit (Surplus)"), alpha = 0.9) +
      geom_vline(xintercept = c(26, 52, 104), color = "gray60", size = 0.5, linetype = 2) +
      ggrepel::geom_label_repel(
        data = peron_winratio %>% filter(peron_min_time %in% milestones) %>%
          mutate(nudge_amount = ifelse(NetBenefit >= 0, 0.03, -0.03)),
        aes(
          x = peron_min_time, 
          y = NetBenefit, 
          label = scales::percent(NetBenefit, accuracy = 0.1)
        ),
        fontface = "bold", 
        size = 3.5, 
        color = "#006064",
        direction = "y",
        nudge_y = peron_winratio %>% filter(peron_min_time %in% milestones) %>%
          mutate(nudge_amount = ifelse(NetBenefit >= 0, 0.01, -0.01)) %>%
          pull(nudge_amount),
        segment.color = "#006064",
        segment.size = 0.4,
        min.segment.length = 0
      ) +
      # SCALES
      coord_cartesian(ylim = c(-0.3, 0.3)) +
      scale_y_continuous(labels = scales::percent, breaks = c(-0.4, -0.3, -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.3, 0.4)) +
      scale_fill_manual(
        values = c(
          "Similar Outcome (Tie)"       = "#CFD8DC", # Light Blue-Grey
          "Favors Experimental Arm"     = "#00E5FF", # Cyan (Your preference)
          "Favors Standard of Care"     = "#5C6BC0", # Soft Indigo (Instead of Scary Red)
          "Net Benefit (Surplus)"       = "#006064",  # Deep Teal/Cyan (Solid)
          "Significant"       = "#006064",
          "Not Significant" = "gray75"
        ),
        breaks = c("Favors Experimental Arm", "Net Benefit (Surplus)", 
                   "Favors Standard of Care", "Similar Outcome (Tie)")
      ) +
      # SCALES
      theme_minimal(base_size = 10) +
      labs(
        x = "Survival Time Threshold (Weeks)",
        y = "NTB Peron (%)"
      ) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = "none"
      ) + 
      # We place them at the top of the plot (e.g., 0.45 or 45%)
      annotate(
        "text", 
        x = v_thresholds + 5, # nudge slightly right
        y = 0.45, 
        label = paste(v_thresholds, "wks"),
        angle = 90, size = 3, color = "gray60"
      ) +
      # Ensure the labels aren't cut off
      coord_cartesian(clip = "off") +
      theme(plot.margin = margin(5, 40, 5, 5))
    

    
    data_batch_manual <- getWinRatioBatchManual(ipd_data_list[[trial_filename_it]]) 
    data_batch_manual <- data_batch_manual %>%
      dplyr::mutate(
        NetBenefit = percWin - percLoss,
        NegLoss = -percLoss
      )
    
    winratio_graph <- data_batch_manual %>%
      ggplot(aes(x = mintime)) +
      geom_col(aes(y = percWin, fill = "Favors Experimental Arm"), alpha = 0.3) +
      geom_col(aes(y = NegLoss, fill = "Favors Standard of Care"), alpha = 0.3) +
      geom_vline(xintercept = c(26, 52, 104), color = "gray60", size = 0.5, linetype = 2) +
      geom_col(aes(y = NetBenefit, fill = "Net Benefit (Surplus)"), alpha = 0.7) +
      geom_point(aes(y = NetBenefit, fill = "Net Benefit (Surplus)"), alpha = 0.9, size = 0.7) +
      geom_line(aes(y = NetBenefit, fill = "Net Benefit (Surplus)"), alpha = 0.9) +
      scale_fill_manual(
        values = c(
          "Similar Outcome (Tie)"       = "#CFD8DC", # Light Blue-Grey
          "Favors Experimental Arm"     = "#00E5FF", # Cyan (Your preference)
          "Favors Standard of Care"     = "#5C6BC0", # Soft Indigo (Instead of Scary Red)
          "Net Benefit (Surplus)"       = "#006064"  # Deep Teal/Cyan (Solid)
        ),
        breaks = c("Favors Experimental Arm", "Net Benefit (Surplus)", 
                   "Favors Standard of Care", "Similar Outcome (Tie)")
      ) +
      ggrepel::geom_label_repel(
        data = data_batch_manual %>% 
          filter(mintime %in% milestones) %>%
          mutate(nudge_amount = ifelse(NetBenefit >= 0, 0.03, -0.03)),
        aes(
          x = mintime, 
          y = NetBenefit, 
          label = scales::percent(NetBenefit, accuracy = 0.1)
        ),
        fontface = "bold", 
        size = 3.5, 
        color = "#006064",
        direction = "y",
        nudge_y = data_batch_manual %>% 
          filter(mintime %in% milestones) %>%
          mutate(nudge_amount = ifelse(NetBenefit >= 0, 0.01, -0.01)) %>%
          pull(nudge_amount),
        segment.color = "#006064",
        segment.size = 0.4,
        min.segment.length = 0
      )  + 
      # SCALES
      coord_cartesian(ylim = c(-0.3, 0.3)) +
      scale_y_continuous(labels = scales::percent, breaks = c(-0.4, -0.3, -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.3, 0.4)) +
      theme_minimal(base_size = 10) +
      labs(
        x = "Survival Time Threshold (Weeks)",
        y = "Manual NTB (%)"
      ) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = "none"
      ) + 
      # We place them at the top of the plot (e.g., 0.45 or 45%)
      annotate(
        "text", 
        x = v_thresholds + 5, # nudge slightly right
        y = 0.45, 
        label = paste(v_thresholds, "wks"),
        angle = 90, size = 3, color = "gray60"
      ) +
      # Ensure the labels aren't cut off
      coord_cartesian(clip = "off") +
      theme(plot.margin = margin(5, 40, 5, 5))
    
    
    events_prop <- data_batch_manual %>%
      dplyr::filter(trial_filename == trial_filename_it) %>%
      tidyr::pivot_longer(
        cols = c("percWin", "percLoss", "percFalseTieEvent", "percFalseTieCensored", 
                 "percTrueTieEvent", "percTrueTieCensored"), 
        names_to = "MatchResult", 
        values_to = "MatchNumber"
      ) %>% 
      dplyr::mutate(MatchResult = factor(MatchResult, levels = c(
        "percFalseTieEvent", 
        "percFalseTieCensored",
        "percTrueTieEvent", 
        "percTrueTieCensored",
        "percLoss",
        "percWin"
      ))) %>%
      ggplot(aes(x = mintime, y = MatchNumber, fill = MatchResult, color = MatchResult)) + 
      geom_col(alpha = 0.8, width = 1) + 
      # Brown vertical dashed lines
      # geom_vline(xintercept = c(26, 52, 104), 
      #            linetype = "dashed", 
      #            color = "#8B4513", 
      #            linewidth = 0.7) +
      # # Labels positioned at Inf but with negative vjust to push them ABOVE the plot
      # annotate("text", 
      #          x = c(26, 52, 104), 
      #          y = Inf, 
      #          label = c("26w", "52w", "104w"), 
      #          color = "#8B4513",
      #          angle = 0,         
      #          vjust = -1.0,       # Pushes text slightly above the top line
      #          size = 3.5, 
      #          fontface = "bold") +
      scale_fill_manual(values = semantic_colors, name = "Match Result") +
      scale_color_manual(values = semantic_colors, name = "Match Result") +
      scale_x_continuous(breaks = c(0, 4, 26, 52, 104, 208)) +
      scale_y_continuous(labels = scales::percent_format(scale = 1)) + 
      coord_cartesian(clip = 'off') + # CRITICAL: Allows labels to be visible outside the grid
      labs(
        x = "Minimum Time Difference (Weeks)",
        y = "% all pairs"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        legend.position = 'bottom',
        panel.grid.minor = element_blank()
      ) +
      guides(fill = guide_legend(nrow = 2))
    
    
    # winratio_crude_graph <- getCrudeWinRatioBatchManual(ipd_data_list[[trial_filename_it]]) %>%
    #   dplyr::mutate(
    #     matchLabel = factor(matchLabel, levels = c(
    #       "0-0-le", "0-1-le", "0-0-ge", "1-0-ge", "1-0-eq", "0-1-eq", "0-0-eq", "1-1-eq", # Tie
    #       "1-0-le", "1-1-le", # Loss
    #       "0-1-ge", "1-1-ge" # Victory
    #     ))
    #   ) %>%
    #   
    #   dplyr::arrange(matchLabel) %>%
    #   
    #   ggplot(aes(x = mintime, y = percOutcomes, fill = matchLabel)) + 
    #   geom_area(
    #     alpha = 0.75,
    #     width = 0.7,
    #     color = "gray50"
    #   ) +
    #   theme_minimal() +
    #   scale_fill_manual(
    #     values = c(
    #       ## eqIES (blue shades)
    #       "0-0-ge" = "#FD8D3C",
    #       "0-0-le" = "#74C476",
    #       "1-0-eq" = "#C6DBEF",
    #       "0-1-eq" = "#9ECAE1",
    #       "0-0-eq" = "#6BAED6",
    #       "1-1-eq" = "#3182BD",
    #       
    #       ## LOSSES (orange shades)
    #       "1-0-ge" = "#FDD0A2",
    #       "0-1-ge" = "#FDAE6B",
    #       "1-1-ge" = "#E6550D",
    #       
    #       ## VICeqORIES (green shades)
    #       "1-0-le" = "#C7E9C0",
    #       "0-1-le" = "#A1D99B",
    #       "1-1-le" = "#238B45"
    #     )
    #   ) +
    #   theme(legend.position = "right") + 
    #   labs(fill = '') +
    #   guides(fill = guide_legend(ncol = 3)) 
    
    waterfall_difftime_plot <- 
      getPairWiseDiffTime(ipd_data_list[[trial_filename_it]]) %>%
      mutate(diffTimeAbs = cut(abs(diffTime), seq(0, 1000, 4))) %>%
      dplyr::mutate(totalCount = nrow(.)) %>%
      dplyr::ungroup() %>%
      group_by(diffTimeAbs, WinLoss) %>%
      summarise(
        n = n(),
        totalCount = first(totalCount)
      ) %>%
      dplyr::mutate(diffTimeAbs = as.numeric(stringr::str_extract(diffTimeAbs, '\\((\\d+)', group = 1))) %>%
      pivot_wider(names_from = 'WinLoss', values_from = 'n') %>%
      mutate(across(everything(), ~ ifelse(is.na(.), 0, .))) %>%
      mutate(WinMinusLoss = (Win - Loss)/totalCount,
             Signal = ifelse(WinMinusLoss <= 0, 'Negative', 'Positive')) %>%
      ggplot(aes(x=diffTimeAbs, y=WinMinusLoss, fill = Signal)) + 
      geom_col(color = "gray50") +
      scale_fill_manual(values = c('Positive' = '#00E5FF', 'Negative' = '#3F51B5')) +
      theme_minimal(base_size = 10) + theme(legend.position = 'none') + 
        labs(
          x = "Abs Time Difference Between Pairs (Weeks)",
          y = "Wins Minus Losses (%)"
        ) + 
        # annotate(
        #   "text", 
        #   x = v_thresholds + 5, # nudge slightly right
        #   y = 0.005, 
        #   label = paste(v_thresholds, "wks"),
        #   angle = 90, size = 3, color = "gray60"
        # ) +
        # coord_cartesian(clip = "off") +
        scale_y_continuous(labels = scales::percent) + 
        theme(plot.margin = margin(5, 40, 5, 5)) + 
        geom_vline(xintercept = c(26, 52, 104), color = "gray60", size = 0.5, linetype = 2)
  
    
    # waffle_plot <- ipd_winratio_estimates_df %>%
    #   dplyr::filter(trial_filename == trial_filename_it) %>%
    #   dplyr::filter(mintime %in% c(4, 26, 52)) %>%
    #   
    #   # --- DATA PREP ---
    #   mutate(
    #     countWin  = round(percWin * 100),
    #     countLoss = round(percLoss * 100),
    #     countTie  = 100 - (countWin + countLoss)
    #   ) %>%
    #   mutate(dots = pmap(list(countWin, countTie, countLoss), function(w, t, l) {
    #     c(rep("Experimental", w), rep("Tie", t), rep("SOC", l))
    #   })) %>%
    #   unnest(dots) %>%
    #   group_by(mintime) %>%
    #   mutate(
    #     dot_id = row_number(),
    #     x_pos = ((dot_id - 1) %% 10) + 1, 
    #     y_pos = ceiling(dot_id / 10)
    #   ) %>%
    #   
    #   # --- CALCULATE SURPLUS (CORRECTED LOGIC) ---
    #   group_by(mintime, dots) %>%
    #   mutate(rank_in_group = row_number()) %>%
    #   group_by(mintime) %>%
    #   mutate(
    #     total_exp = sum(dots == "Experimental"),
    #     total_soc = sum(dots == "SOC"),
    #     
    #     is_surplus = case_when(
    #       # Case 1: Experimental Wins
    #       # Highlight the LAST Experimental dots (those closest to the center/ties)
    #       dots == "Experimental" & total_exp > total_soc & rank_in_group > total_soc ~ TRUE,
    #       
    #       # Case 2: SOC Wins (FIXED)
    #       # Highlight the FIRST SOC dots (those closest to the center/ties)
    #       # We select ranks 1 up to the difference
    #       dots == "SOC" & total_soc > total_exp & rank_in_group <= (total_soc - total_exp) ~ TRUE,
    #       
    #       TRUE ~ FALSE
    #     )
    #   ) %>%
    #   ungroup() %>%
    #   mutate(dots = factor(dots, levels = c("Experimental", "Tie", "SOC"))) %>%
    #   
    #   # --- THE PLOT ---
    #   ggplot(aes(x = x_pos, y = y_pos)) +
    #   
    #   # SHAPE 21: Filled Circle with Colored Outline
    #   geom_point(aes(fill = dots, color = is_surplus), 
    #              shape = 21, 
    #              size = 3,      
    #              stroke = 1) + # Moderate thickness for elegance
    #   
    #   facet_wrap(~ mintime, labeller = label_both, ncol = 4) +
    #   
    #   # FILL COLORS (Inside the circle)
    #   scale_fill_manual(
    #     name = "Category",
    #     values = c(
    #       "Experimental" = "#00E5FF", 
    #       "Tie"          = "#CFD8DC", 
    #       "SOC"          = "#3F51B5"
    #     )
    #   ) +
    #   
    #   # OUTLINE COLORS (The Ring)
    #   scale_color_manual(
    #     values = c(
    #       "FALSE" = 'white',        
    #       "TRUE"  = "black"  
    #     ),
    #     guide = "none"
    #   ) +
    #   
    #   theme_minimal() +
    #   theme(
    #     axis.text = element_blank(),
    #     axis.title = element_blank(),
    #     panel.grid = element_blank(),
    #     strip.text = element_text(size = 12, face = "bold"),
    #     strip.background = element_rect(fill = "gray95", color = NA),
    #     legend.position = "top",
    #     panel.spacing = unit(1, "lines")
    #   )
    
    # Left column with specific heights for each plot
    left_col <- cox_graph / winratio_graph / gehan_net_benefit_graph +
      plot_layout(heights = c(3.5, 0.8, 2, 2))  # Cox gets 4 units, others get 2.5 each
    
    # Right column with same heights
    right_col <- censoring_graph / waterfall_difftime_plot / peron_net_benefit_graph +
      plot_layout(heights = c(3.5, 0.8, 2, 2))  # Censoring gets 4 units, others get 2.5 each
    
    # Combine everything
    patch_graph <- (left_col | right_col) / events_prop +
      plot_layout(heights = c(9, 1)) +  # Top section 90%, bottom 10%
      theme(
        plot.background  = element_blank(),
        panel.background = element_blank()
      )
    
    ggsave(plot = patch_graph, sprintf('figures/ind_trials/%03d_%s_figures.pdf', i, trial_filename_it), 
           width = 25*1, height = 29*1, unit = 'cm')
    
  } %try% NULL
  
}


# ========================================================================================================
# Plotting Single Graphs of Correlation
final_results_os_df <- final_results_df %>%
  dplyr::filter(trial_endpoint == 'OS')

# Median Difference and Peron NB Estimates
ggplot(aes(x = median_diff, y = peron_nb_estimate_52wk, 
           color = cox_p_value <= 0.05), data = final_results_os_df) +
  geom_point(size=2) +
  ggrepel::geom_text_repel(aes(label = trial_filename), max.overlaps = 10,
                          data = final_results_os_df %>% dplyr::filter(cox_p_value <= 0.05)) +
  geom_hline(yintercept = c(0.05, -0.05), linetype = "dashed", color = "gray50") + 
  geom_vline(xintercept = c(-16, 16), linetype = "dashed", color = "gray50") + 
  scale_x_continuous(breaks = seq(-80, 80, by = 8)) +  # Adjust range/interval as needed
      scale_y_continuous(breaks = seq(-0.3, 0.3, by = 0.05), labels = scales::percent_format(accuracy = 1)) + 
  theme_bw()

ggsave("figures/scatter_median_nb_52wks.pdf", width = 16*1.2, height = 11*1.2, unit = 'cm')

ggplot(aes(x = max_rmean_diff, y = peron_nb_estimate_52wk, 
           color = cox_p_value <= 0.05), data = final_results_os_df) +
      geom_point(size=2) +
      ggrepel::geom_text_repel(aes(label = trial_filename), max.overlaps = 10,
                               data = final_results_os_df %>% dplyr::filter(cox_p_value <= 0.05)) +
      geom_hline(yintercept = c(0.05, -0.05), linetype = "dashed", color = "gray50") + 
      geom_vline(xintercept = c(-16, 16), linetype = "dashed", color = "gray50") + 
      scale_x_continuous(breaks = seq(-80, 80, by = 8)) +  # Adjust range/interval as needed
      scale_y_continuous(breaks = seq(-0.3, 0.3, by = 0.05), labels = scales::percent_format(accuracy = 1)) + 
      theme_bw()

ggsave("figures/scatter_rmean_nb_52wks.pdf", width = 16*1.2, height = 11*1.2, unit = 'cm')

ggplot(aes(x = landmark_diff_12mo, y = peron_nb_estimate_52wk, 
           color = cox_p_value <= 0.05), data = final_results_os_df) +
      geom_point(size=2) +
      ggrepel::geom_text_repel(aes(label = trial_filename), max.overlaps = 10,
                               data = final_results_os_df %>% dplyr::filter(cox_p_value <= 0.05)) +
      geom_hline(yintercept = c(0.05, -0.05), linetype = "dashed", color = "gray50") + 
      geom_vline(xintercept = c(0.05, -0.05), linetype = "dashed", color = "gray50") + 
      scale_x_continuous(breaks = seq(-0.3, 0.3, by = 0.05), labels = scales::percent_format(accuracy = 1)) + 
      scale_y_continuous(breaks = seq(-0.3, 0.3, by = 0.05), labels = scales::percent_format(accuracy = 1)) + 
      theme_bw()

ggsave("figures/scatter_landamrk_52wks_nb_52wks.pdf", width = 16*1.2, height = 11*1.2, unit = 'cm')

# Median Difference and Peron NB Estimates
ggplot(aes(x = median_diff, y = peron_nb_estimate_26wk, 
           color = cox_p_value <= 0.05), data = final_results_os_df) +
      geom_point(size=2) +
      ggrepel::geom_text_repel(aes(label = trial_filename), max.overlaps = 10,
                               data = final_results_os_df %>% dplyr::filter(cox_p_value <= 0.05)) +
      geom_hline(yintercept = c(0.01, -0.01), linetype = "dashed", color = "gray50") + 
      geom_vline(xintercept = c(-16, 16), linetype = "dashed", color = "gray50") + 
      scale_x_continuous(breaks = seq(-80, 80, by = 8)) +  # Adjust range/interval as needed
      scale_y_continuous(breaks = seq(-0.3, 0.3, by = 0.05), labels = scales::percent_format(accuracy = 1)) + 
      theme_bw()

ggsave("figures/scatter_median_nb_26wks.pdf", width = 16*1.2, height = 11*1.2, unit = 'cm')

ggplot(aes(x = max_rmean_diff, y = peron_nb_estimate_26wk, 
           color = cox_p_value <= 0.05), data = final_results_os_df) +
      geom_point(size=2) +
      ggrepel::geom_text_repel(aes(label = trial_filename), max.overlaps = 10,
                               data = final_results_os_df %>% dplyr::filter(cox_p_value <= 0.05)) +
      geom_hline(yintercept = c(0.1, -0.1), linetype = "dashed", color = "gray50") + 
      geom_vline(xintercept = c(-16, 16), linetype = "dashed", color = "gray50") + 
      scale_x_continuous(breaks = seq(-80, 80, by = 8)) +  # Adjust range/interval as needed
      scale_y_continuous(breaks = seq(-0.3, 0.3, by = 0.05), labels = scales::percent_format(accuracy = 1)) + 
      theme_bw()

ggsave("figures/scatter_rmean_nb_26wks.pdf", width = 16*1.2, height = 11*1.2, unit = 'cm')

ggplot(aes(x = landmark_diff_6mo, y = peron_nb_estimate_26wk, 
           color = cox_p_value <= 0.05), data = final_results_os_df) +
      geom_point(size=2) +
      ggrepel::geom_text_repel(aes(label = trial_filename), max.overlaps = 10,
                               data = final_results_os_df %>% dplyr::filter(cox_p_value <= 0.05)) +
      geom_hline(yintercept = c(0.05, -0.05), linetype = "dashed", color = "gray50") + 
      geom_vline(xintercept = c(0.05, -0.05), linetype = "dashed", color = "gray50") + 
      scale_x_continuous(breaks = seq(-0.3, 0.3, by = 0.05), labels = scales::percent_format(accuracy = 1)) + 
      scale_y_continuous(breaks = seq(-0.3, 0.3, by = 0.05), labels = scales::percent_format(accuracy = 1)) + 
      theme_bw()

ggsave("figures/scatter_landmark_26wks_nb_26wks.pdf", width = 16*1.2, height = 11*1.2, unit = 'cm')
