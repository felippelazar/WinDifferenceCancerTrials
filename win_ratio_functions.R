# ============================================================================ #
# Win-Ratio Functions
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
library(survRM2)
source('R/tryRetry.R')

# ============================================================================ #
# 1. Getting Strata Names from a survfit Object

tidy_mstrata_survfit <- function(survfit_obj, stat_name = 'median', times_interval = c(6, 12, 18, 24, 30, 36)){
  
  if(stat_name == 'median'){
    table <- as.data.frame(summary(survfit_obj)$table)
    table <- table %>% mutate(strata = rownames(table)) %>% 
      relocate(strata, .before = records) %>%
      rename(conf.low = `0.95LCL`,
             conf.high = `0.95UCL`) %>%
      mutate(med_est_interval = sprintf('%.2f (%.2f - %.2f)', median, conf.low, conf.high))
    
    rownames(table) <- NULL
    return(table)
    
  }else if(stat_name == 'times'){
    summary_survfit_obj <- summary(survfit_obj, times = times_interval)
    if(!'strata' %in% names(summary_survfit_obj)){summary_survfit_obj$strata <- rep('all', length(summary_survfit_obj$time))}
    table <- data.frame(strata = summary_survfit_obj$strata,
                        time = summary_survfit_obj$time, 
                        n.risk = summary_survfit_obj$n.risk,
                        n.event = summary_survfit_obj$n.event,
                        n.censor = summary_survfit_obj$n.censor,
                        surv_estimate = summary_survfit_obj$surv,
                        surv_lower = summary_survfit_obj$lower,
                        surv_higher = summary_survfit_obj$upper)
    
    table <- table %>%
      mutate(surv_est_interval = sprintf('%.2f (%.2f - %.2f)', surv_estimate, surv_lower, surv_higher)) %>%
      relocate(surv_est_interval, .before = surv_estimate)
    return(table)
  }
  
}

# ============================================================================ #
# 2. Defining Win/Loss/Tie Functions

defWinLossComplete <- function(time1, outcome1, time2, outcome2, mintime = 0, tie_threshold = 4) {
  # difftime > 0 means time1 is later than time2
  difftime <- time1 - time2 
  abs_difftime <- abs(difftime)
  
  comparisonResult <- case_when(
    
    # 1. WIN CASES: Prioritize all Win conditions
    # 1a. S1 Censored vs S2 Died, but S1 observed significantly longer
    (outcome1 == 0) & (outcome2 == 1) & (difftime > mintime) ~ 'Win',
    # 1b. Both Died, but S1 survived significantly longer
    (outcome1 == 1) & (outcome2 == 1) & (difftime > mintime) ~ 'Win',
    
    # 2. LOSS CASES: Prioritize all Loss conditions
    # 2a. S1 Died vs S2 Censored, but S1 died significantly earlier
    (outcome1 == 1) & (outcome2 == 0) & (difftime < -mintime) ~ 'Loss',
    # 2b. Both Died, but S1 died significantly earlier
    (outcome1 == 1) & (outcome2 == 1) & (difftime < -mintime) ~ 'Loss',
    
    # --- TIE CLASSIFICATIONS ---
    
    # 3. TRUE TIE 1: Both Censored, and follow-up times are very similar
    # The absolute difference is less than or equal to the tie_threshold (e.g., 30 days)
    (outcome1 == 0) & (outcome2 == 0) & (abs_difftime <= tie_threshold) ~ 'True Tie (Censored)',
    
    # 4. TRUE TIE 2: Both Died, but the difference does not meet the mintime threshold
    # Note: This is now explicitly pulled out, though the current logic handles it with the next TRUE condition.
    # We only need this if you want to explicitly name the tie. Otherwise, it's covered by the final TRUE.
    (outcome1 == 1) & (outcome2 == 1) & (abs_difftime <= mintime) ~ 'True Tie (Event)',
    
    # 5. FALSE TIE (Ambiguous/Default Tie)
    # This catches:
    # - Censored/Death where time difference was not significant enough for Win/Loss.
    # - Censored/Censored where follow-up difference was large (abs_difftime > tie_threshold).
    # - Death/Censored where time difference was not significant enough for Win/Loss.
    (outcome1 == 1) & (outcome2 == 0) ~ 'False Tie (Any Event)',
    (outcome1 == 0) & (outcome2 == 1) ~ 'False Tie (Any Event)',
    (outcome1 == 0) & (outcome2 == 0) ~ 'False Tie (Censored)',
    TRUE ~ NA_character_
  )
  
  return(comparisonResult)
  
}

defWinLoss <- function(time1, outcome1, time2, outcome2, mintime = 0) {
  # difftime > 0 means time1 is later than time2
  difftime <- time1 - time2 
  
  comparisonResult <- case_when(
    
    # 1. Subject 1 Censored (0) vs. Subject 2 Died (1) -> Win only if observed time is significantly longer
    (outcome1 == 0) & (outcome2 == 1) & (difftime > mintime) ~ 'Win',
    
    # 2. Subject 1 Died (1) vs. Subject 2 Censored (0) -> Loss only if death time is significantly earlier
    (outcome1 == 1) & (outcome2 == 0) & (difftime < -mintime) ~ 'Loss',
    
    # 3. Both Died (1 vs 1) -> Win if Subject 1 died LATER
    (outcome1 == 1) & (outcome2 == 1) & (difftime > mintime) ~ 'Win',
    
    # 4. Both Died (1 vs 1) -> Loss if Subject 1 died EARLIER
    (outcome1 == 1) & (outcome2 == 1) & (difftime < -mintime) ~ 'Loss',
    
    # 5. Catch-all for Tie
    # This includes:
    # - Both Censored (0 vs 0).
    # - Censored/Death comparisons where the time difference was small (i.e., |difftime| <= mintime).
    # - Death/Death comparison where the time difference was small (i.e., |difftime| <= mintime).
    TRUE ~ 'Tie'
  )
  
  return(comparisonResult)
}

# ============================================================================ #
# 3. Getting Survival Estimates

getSurvivalEstimates <- function(data){
  
  print(first(data$trial_name))
  
  # Refactoring the Treatment Arm
  data <- data %>%
    dplyr::mutate(Arm = factor(Arm, levels = unique(.$Arm)))
  
  # Median Survival Estimates
  median_data <- survfit(Surv(TimeWeeks, Event) ~ Arm, data = data) %>%
    tidy_mstrata_survfit(stat_name = "median") %>%
    summarise(
      median_trt_estimate = median[[1]],
      median_ctl_estimate = median[[2]],
      median_diff = median[[1]] - median[[2]],
    ) %>%
    dplyr::mutate(log_rank_p_value = survdiff(Surv(TimeWeeks, Event) ~ Arm, data = data)$p) %>%
    dplyr::mutate(inf_censoring_p_value = survdiff(Surv(TimeWeeks, Event==0) ~ Arm, data = data)$p %try% NULL)
  
  fit_fup  <- survfit(Surv(TimeWeeks, Event == 0) ~ 1, data = data)
  fit_arm  <- survfit(Surv(TimeWeeks, Event) ~ Arm, data = data)
  
  n_events_trt <- summary(fit_arm)$table[1, "events"]
  n_events_ctl <- summary(fit_arm)$table[2, "events"]
  n_events_total <- n_events_trt + n_events_ctl
  
  n_patients_trt <- summary(fit_arm)$table[1, "records"]
  n_patients_ctl <- summary(fit_arm)$table[2, "records"]
  n_patients_total <- nrow(data)
  
  # Median Survival Estimates
  fup_data <- data.frame(
    median_fup_estimate   = summary(fit_fup)$table["median"],
    median_fup_conf_low   = summary(fit_fup)$table["0.95LCL"],
    median_fup_conf_high  = summary(fit_fup)$table["0.95UCL"],
    n_events_total        = n_events_total,
    n_events_trt          = n_events_trt,
    n_events_ctl          = n_events_ctl,
    n_patients_total      = n_patients_total,
    n_patients_trt        = n_patients_trt,
    n_patients_ctl        = n_patients_ctl,
    rate_events_total     = n_events_total / n_patients_total,
    rate_events_trt       = n_events_trt / n_patients_trt,
    rate_events_ctl       = n_events_ctl / n_patients_ctl
  )
  
  tau_points <- c(52/2, 52, 52*2, 52*3, 52*4, 52*5, min(tapply(data$TimeWeeks, data$Arm, max))) # 6mo, 12mo, 24mo, 36mo, 48mo
  
  results <- lapply(tau_points, function(t) {fit <- rmst2(time = data$TimeWeeks, status = data$Event, arm = as.numeric(as.factor(data$Arm)) - 1, tau = t) %try% NULL})
  
  names(results) <- c("26wk", "52wk", "104wk", "156wk", "208wk", "260wk", "max")
  
  rmean_landmark_table <- imap_dfr(results, function(fit, label) {
    if (is.null(fit)) return(NULL) # Handles the %try% NULL cases
    
    # Extract the difference row from the unadjusted results table
    # Row 2 is typically the "RMST (arm=1)-(arm=0)" row
    diff_row <- as.data.frame(fit$unadj) %>% 
      slice(2) %>% 
      select(Est., p) %>%
      rename(rmean_diff = Est., p_value = p)
    
    mutate(diff_row, landmark = label)
  }) %>%
    select(landmark, rmean_diff, p_value) %>%
    mutate(landmark = factor(landmark, levels = c("26wk", "52wk", "104wk", "156wk", "208wk", "260wk", "max"))) %>%
    pivot_wider(
      names_from = landmark,
      values_from = c(rmean_diff, p_value),
      names_glue = "{landmark}_{.value}"
    )
  
  # Landmark Survival Estimates
  landmark_data <- lapply(list(seq(6, 36, 6), seq(6, 30, 6), seq(6, 24, 6), seq(6, 18, 6), seq(6, 12, 6), seq(6, 6, 6)), function(times) {
    
    {survfit(Surv(TimeWeeks*0.2301496, Event) ~ Arm, data = data) %>%
        tidy_mstrata_survfit(stat_name = "times", times_interval = times) %>%
        group_by(time) %>%
        summarise(
          landmark_trt = surv_estimate[[1]],
          landmark_ctl = surv_estimate[[2]],
          landmark_diff = surv_estimate[[1]] - surv_estimate[[2]]
        ) %>%
        pivot_wider(names_from = "time", values_from = c("landmark_trt", "landmark_ctl", "landmark_diff")) %>%
        setNames(., paste0(colnames(.), "mo"))} %try% NULL
    
  }) 
  
  landmark_data <- purrr::compact(landmark_data)[[1]]
  
  # Getting HR Cox Proportional Data
  cox_data <- coxph(Surv(TimeWeeks, Event) ~ Arm, data = data %>% dplyr::mutate(Arm = factor(Arm, levels = rev(unique(.$Arm))))) %>%
    broom.helpers::tidy_plus_plus(exponentiate = T) %>%
    mutate(est.conf.interval = sprintf('%.2f (%.2f - %.2f)', estimate, conf.low, conf.high)) %>%
    dplyr::filter(!is.na(statistic)) %>%
    dplyr::select(
      cox_term = term,
      cox_estimate = estimate,
      cox_conf_low = conf.low,
      cox_conf_high = conf.high, 
      cox_p_value = p.value,
      cox_est_conf_interval = est.conf.interval
    ) %>%
    dplyr::mutate(cox_zph = cox.zph(coxph(Surv(TimeWeeks, Event) ~ Arm, data = data %>% dplyr::mutate(Arm = factor(Arm, levels = rev(unique(.$Arm))))))$table[1,3])
  
  trial_metadata <- data %>%
    dplyr::distinct(trial_name, trial_filename, pub_med_id)
  
  return(cbind(trial_metadata, median_data, rmean_landmark_table, landmark_data, cox_data, fup_data))
}

# ============================================================================ #
# 4. Get Win Ratio Batch on Files

getWinRatioBatchManual <- function(data){
  
  max_fup_time <- max(data$TimeWeeks, na.rm=TRUE)
  max_time <- ceiling(max_fup_time/2)*2
  
  # Getting Complete Win Ratio Data
  groupped_data <- data %>%
    dplyr::mutate(ID = row_number()) %>%
    dplyr::mutate(Arm = factor(Arm, levels = unique(.$Arm))) %>%
    dplyr::select(ID, TimeWeeks, Event, Arm) %>%
    group_by(Arm) %>%
    group_split()
  
  groupped_data <- lapply(1:length(groupped_data), function(x){
    setNames(groupped_data[[x]],  paste0(colnames(groupped_data[[x]]), x))
  })
  
  expanded_grid <- expand.grid(groupped_data[[1]]$ID1, groupped_data[[2]]$ID2) %>%
    dplyr::rename(ID1 = Var1, ID2 = Var2) %>%
    dplyr::left_join(groupped_data[[1]]) %>%
    dplyr::left_join(groupped_data[[2]]) %>%
    tidyr::crossing(mintime = seq(0, max_time, 2)) %>%
    dplyr::mutate(
      WinLoss = defWinLossComplete(
        TimeWeeks1, Event1,
        TimeWeeks2, Event2,
        mintime = mintime
      )
    )
  
  win_ratio_data <-  expanded_grid %>%
    group_by(mintime) %>%
    summarise(
      totalMatches = n(),
      totalWin = sum(WinLoss == "Win"),
      totalLoss = sum(WinLoss == "Loss"),
      totalFalseTie = sum(WinLoss == "False Tie (Any Event)" | WinLoss == "False Tie (Censored)"),
      totalFalseTieEvent = sum(WinLoss == "False Tie (Any Event)"),
      totalFalseTieCensored = sum(WinLoss == "False Tie (Censored)"),
      totalTrueTieCensored = sum(WinLoss == "True Tie (Censored)"),
      totalTrueTieEvent = sum(WinLoss == "True Tie (Event)"),
      totalTrueTie = sum(WinLoss == "True Tie (Event)" | WinLoss == "True Tie (Censored)"),
      totalTie = sum(WinLoss == "True Tie (Event)" | WinLoss == "True Tie (Censored)" | WinLoss == "False Tie (Any Event)" | WinLoss == "False Tie (Censored)"),
      
      percWin = totalWin/totalMatches,
      percLoss = totalLoss/totalMatches,
      percFalseTie = totalFalseTie/totalMatches,
      percFalseTieCensored = totalFalseTieCensored/totalMatches,
      percFalseTieEvent = totalFalseTieEvent/totalMatches,
      percTrueTieCensored = totalTrueTieCensored/totalMatches,
      percTrueTieEvent = totalTrueTieEvent/totalMatches,
      percTrueTie = totalTrueTie/totalMatches,
      percTie = totalTie/totalMatches,
      
      winRatioNoTies = totalWin/totalLoss,
      successRatio = (totalWin/totalMatches + (totalLoss/totalMatches*0.5)) / (totalLoss/totalMatches + (totalLoss/totalMatches*0.5)),
      winDifference = (totalWin - totalLoss)/totalMatches
    )
  
  trial_metadata <- data %>%
    dplyr::distinct(trial_name, trial_filename, pub_med_id)
  
  return(cbind(win_ratio_data, trial_metadata))
}


# ============================================================================ #
# 5. Get Crude Win Ratio Estimates on Batch

getCrudeWinRatioBatchManual <- function(data){
  
  max_fup_time <- max(data$TimeWeeks, na.rm=TRUE)
  max_time <- ceiling(max_fup_time/2)*2
  
  # Getting Complete Win Ratio Data
  groupped_data <- data %>%
    dplyr::mutate(ID = row_number()) %>%
    dplyr::mutate(Arm = factor(Arm, levels = unique(.$Arm))) %>%
    dplyr::select(ID, TimeWeeks, Event, Arm) %>%
    group_by(Arm) %>%
    group_split()
  
  groupped_data <- lapply(1:length(groupped_data), function(x){
    setNames(groupped_data[[x]],  paste0(colnames(groupped_data[[x]]), x))
  })
  
  expanded_grid <- expand.grid(groupped_data[[1]]$ID1, groupped_data[[2]]$ID2) %>%
    dplyr::rename(ID1 = Var1, ID2 = Var2) %>%
    dplyr::left_join(groupped_data[[1]]) %>%
    dplyr::left_join(groupped_data[[2]]) %>%
    tidyr::crossing(mintime = seq(0, max_time, 2)) %>%
    dplyr::mutate(
      diffTime = TimeWeeks1 - TimeWeeks2,
      outcomeMatch = case_when(
        diffTime >= mintime ~ 'ge',
        diffTime <= -mintime ~ 'le',
        TRUE ~ 'eq'
      ),
      matchLabel = paste0(Event1, "-", Event2, "-", outcomeMatch),
      WinLoss = defWinLossComplete(
        TimeWeeks1, Event1,
        TimeWeeks2, Event2,
        mintime = mintime
      )
    )
  
  crude_win_ratio_data <-  expanded_grid %>%
    group_by(mintime) %>%
    dplyr::mutate(totalMatches = n()) %>%
    ungroup() %>%
    group_by(mintime, matchLabel) %>%
    dplyr::summarise(
      totalOutcomes = n(),
      percOutcomes = n()/first(totalMatches)
    )
  
  trial_metadata <- data %>%
    dplyr::distinct(trial_name, trial_filename, pub_med_id)
  
  return(cbind(crude_win_ratio_data, trial_metadata))
}

# ============================================================================ #
# 6. Get Pairwise Diff Time

getPairWiseDiffTime <- function(data){
  
  # Getting Complete Win Ratio Data
  groupped_data <- data %>%
    dplyr::mutate(ID = row_number()) %>%
    dplyr::mutate(Arm = factor(Arm, levels = unique(.$Arm))) %>%
    dplyr::select(ID, TimeWeeks, Event, Arm) %>%
    group_by(Arm) %>%
    group_split()
  
  groupped_data <- lapply(1:length(groupped_data), function(x){
    setNames(groupped_data[[x]],  paste0(colnames(groupped_data[[x]]), x))
  })
  
  expanded_grid <- expand.grid(groupped_data[[1]]$ID1, groupped_data[[2]]$ID2) %>%
    dplyr::rename(ID1 = Var1, ID2 = Var2) %>%
    dplyr::left_join(groupped_data[[1]]) %>%
    dplyr::left_join(groupped_data[[2]]) %>%
    tidyr::crossing(mintime = 1) %>%
    dplyr::mutate(
      diffTime = TimeWeeks1 - TimeWeeks2,
      outcomeMatch = case_when(
        diffTime >= mintime ~ '1V',
        diffTime <= -mintime ~ '1L',
        TRUE ~ 'T'
      ),
      matchLabel = paste0(Event1, "-", Event2, "-", outcomeMatch),
      WinLoss = defWinLossComplete(
        TimeWeeks1, Event1,
        TimeWeeks2, Event2,
        mintime = mintime
      )
    )
  
  trial_metadata <- data %>%
    dplyr::distinct(trial_name, trial_filename, pub_med_id)
  
  return(cbind(expanded_grid, trial_metadata))
  
}

# ============================================================================ #
# 7. Get BuyseTest Peron Method

getBuyseTestPeron <- function(data) {
  
  max_fup_time <- max(data$TimeWeeks, na.rm=TRUE)
  max_time <- ceiling(max_fup_time/2)*2
  
  data <- data %>%
    dplyr::mutate(Arm = factor(Arm, levels = rev(unique(.$Arm))))
  
  models_peron <- lapply(seq(0, max_time, 2), function(x) {
    
    model <- BuyseTest(
      formula          = Arm ~ TTE(TimeWeeks, status = Event, threshold = x), 
      data             = data, 
      scoring.rule     = "Peron", 
      method.inference = "u-statistic"
    )})
  
  confint_netbenefit <- lapply(models_peron, function(model_peron) {
    confint(model_peron, statistic = "netBenefit") %>%
      setNames(., paste0("peron_nb_", names(.))) %>%
      dplyr::mutate(peron_min_time = as.numeric(model_peron@threshold))
  })
  
  confint_winratio <- lapply(models_peron, function(model_peron) {
    confint(model_peron, statistic = "winRatio") %>%
      setNames(., paste0("peron_wr_", names(.))) %>%
      dplyr::mutate(peron_min_time = as.numeric(model_peron@threshold))
  })
  
  counts_matches <- lapply(models_peron, function(model_peron) {
    
    data.frame(
      peron_min_time = model_peron@threshold,
      peron_total_wins = sum(model_peron@count.favorable),
      peron_total_loss = sum(model_peron@count.unfavorable),
      peron_total_ties = sum(model_peron@count.neutral),
      peron_total_uninformative = sum(model_peron@count.uninf)
    )
    
  })
  
  final_table <- do.call(bind_rows, counts_matches) %>%
    dplyr::left_join(do.call(bind_rows, confint_netbenefit)) %>%
    dplyr::left_join(do.call(bind_rows, confint_winratio)) %>%
    dplyr::mutate(
      trial_name = first(data$trial_name),
      trial_filename = first(data$trial_filename)
    )
  
  return(final_table)
  
}

# ============================================================================ #
# 8. Get BuyseTest Gehan Method

getBuyseTestGehan <- function(data){
  
  max_fup_time <- max(data$TimeWeeks, na.rm=TRUE)
  max_time <- ceiling(max_fup_time/2)*2
  
  data <- data %>%
    dplyr::mutate(Arm = factor(Arm, levels = rev(unique(.$Arm))))
  
  models_gehan <- lapply(seq(0, max_time, 2), function(x) {
    
    model <- BuyseTest(
      formula          = Arm ~ TTE(TimeWeeks, status = Event, threshold = x), 
      data             = data, 
      scoring.rule     = "Gehan", 
      method.inference = "u-statistic"
    )})
  
  confint_netbenefit <- lapply(models_gehan, function(model_gehan) {
    confint(model_gehan, statistic = "netBenefit") %>%
      setNames(., paste0("gehan_nb_", names(.))) %>%
      dplyr::mutate(gehan_min_time = as.numeric(model_gehan@threshold))
  })
  
  confint_winratio <- lapply(models_gehan, function(model_gehan) {
    confint(model_gehan, statistic = "winRatio") %>%
      setNames(., paste0("gehan_wr_", names(.))) %>%
      dplyr::mutate(gehan_min_time = as.numeric(model_gehan@threshold))
  })
  
  counts_matches <- lapply(models_gehan, function(model_gehan) {
    
    data.frame(
      gehan_min_time = model_gehan@threshold,
      gehan_total_wins = sum(model_gehan@count.favorable),
      gehan_total_loss = sum(model_gehan@count.unfavorable),
      gehan_total_ties = sum(model_gehan@count.neutral),
      gehan_total_uninformative = sum(model_gehan@count.uninf)
    )
    
  })
  
  final_table <- do.call(bind_rows, counts_matches) %>%
    dplyr::left_join(do.call(bind_rows, confint_netbenefit)) %>%
    dplyr::left_join(do.call(bind_rows, confint_winratio)) %>%
    dplyr::mutate(
      trial_name = first(data$trial_name),
      trial_filename = first(data$trial_filename)
    )
  
  return(final_table)
  
  
}


