# Three-model Cox analysis of prepared two-year landmark data; covariates are already imputed.
# Set BMI_T2DM_ROOT to the absolute project directory and prepare the input files.

library(dplyr)
library(survival)
library(broom)
library(data.table)
library(future.apply)
library(survminer)

setwd(file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/observational"))

landmark_data <- fread(file.path(Sys.getenv("BMI_T2DM_ROOT"), "data/observational/proteins_E11_landmark2y.csv"))

protein_names <- c(
  "ADH1B",
  "AP3B1",
  "APOE",
  "ENTPD6",
  "GRP",
  "RBP1"
)

# Input is already prepared: exclude baseline diabetes and early events/censoring;
# followup_landmark_years starts at the two-year landmark.
# Protein standardization and covariate imputation occur upstream.
data0 <- landmark_data %>%
  select(
    eid,
    event,
    followup_landmark_years,
    all_of(protein_names))

cat("2-year landmark participants: ",nrow(data0),"\n")
cat("E11 events after landmark: ",sum(data0$event == 1, na.rm = TRUE),"\n")

cov <- fread(
  file.path(Sys.getenv("BMI_T2DM_ROOT"), "data/observational/micecov.csv")) %>%
  select(-18) %>%
  dplyr::rename(
    TG_lip = TG)
head(cov)

data <- merge(data0,cov,by = "eid")
cat("Participants after covariate merge: ",nrow(data),"\n")

cat_vars <- c("sex","activity","Ethnicity","smoking","education","alcohol","healthy_diet","income")
var1 <- c("age","sex","Ethnicity")
var2 <- c("TDindex","education")
var3 <- c("smoking","alcohol","healthy_diet")

analysis_vars <- unique(c("eid","event","followup_landmark_years",protein_names,var1,var2,var3))

missing_vars <- setdiff(analysis_vars,names(data))
if (length(missing_vars) > 0) {
  stop(
    paste("Missing variables: ",paste(missing_vars, collapse = ", ")))
}

data_analysis <- data %>%select(all_of(analysis_vars))

na_counts <- colSums(is.na(data_analysis))
na_cols <- na_counts[na_counts > 0]
cat("\nVariables with missing values: \n")
print(na_cols)

n_before_complete <- nrow(data_analysis)
data1 <- data_analysis %>%filter(complete.cases(.))
n_after_complete <- nrow(data1)
n_missing_excluded <-n_before_complete -n_after_complete
cat( "\nParticipants before complete-case filtering: ",n_before_complete,"\n")
cat("Participants excluded for missing values: ",n_missing_excluded,"\n")
cat("Complete-case participants: ",n_after_complete,"\n")

cat_vars_present <- intersect(cat_vars,names(data1))
data2 <- data1 %>%
  mutate(across(all_of(cat_vars_present),as.factor))

cat("\nFinal Cox analysis participants: ",nrow(data2), "\n")
cat("Final E11 events: ",sum(data2$event == 1),"\n")
cat("Final non-events/censored participants: ",sum(data2$event == 0),"\n")

stopifnot(all( data2$followup_landmark_years > 0))

model_formulas <- list(
  model1 = paste(
    "Surv(followup_landmark_years, event) ~",
    paste(
      var1,
      collapse = "+")),
  model2 = paste(
    "Surv(followup_landmark_years, event) ~",
    paste(c(var1,var2),collapse = "+")),
  model3 = paste(
    "Surv(followup_landmark_years, event) ~",
    paste(c(var1,var2,var3),collapse = "+")))

print(model_formulas)

all_model_res <- lapply(
  names(model_formulas),
  function(m_name) {
    message(paste(">>> Starting model:",m_name))
    
    res_list <- lapply(
      protein_names,
      function(protein) {
        formula_str <- paste0(model_formulas[[m_name]]," + ",protein)
        tryCatch(
          {
            fit <- coxph(as.formula(formula_str),data = data2)
            broom::tidy(fit,exponentiate = TRUE,conf.int = TRUE) %>%
              filter(term == protein) %>%
              select(
                term,
                estimate,
                std.error,
                p.value,
                conf.low,
                conf.high) %>%
              mutate(
                model_name = m_name,
                protein_name = protein,
                n = fit$n,
                events = fit$nevent)
            
          },
          error = function(e) {
            data.frame(
              model_name = m_name,
              protein_name = protein,
              term = protein,
              estimate = NA,
              std.error = NA,
              p.value = NA,
              conf.low = NA,
              conf.high = NA,
              n = NA,
              events = NA,
              error = conditionMessage(e))
            
          }
        )
      }
    )
    model_res <- bind_rows(res_list) %>%
      arrange(protein_name)
    message(paste(">>> Completed model:",m_name))
    return(model_res)
  }
)

final_res <- bind_rows(all_model_res)

final_res <- final_res %>%
  select(
    model_name,
    protein_name,
    term,
    estimate,
    conf.low,
    conf.high,
    std.error,
    p.value,
    n,
    events,
    everything())
print(final_res)

out_path <- paste0(file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/observational/"),"e11_landmark2y_all_models_results.csv")

cat("\nCox output path (writing disabled in source): ",out_path,"\n")

analysis_data_path <- paste0(file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/observational/"),"e11_landmark2y_cox_population.csv")

sample_summary <- data.table(
  item = c(
    "2-year landmark population",
    "After merging covariates",
    "Excluded because of missing data",
    "Final Cox analysis population",
    "E11 events in final Cox analysis",
    "Non-events or censored participants in final Cox analysis"),
  n = c(
    nrow(data0),
    nrow(data),
    n_missing_excluded,
    nrow(data2),
    sum(data2$event == 1),
    sum(data2$event == 0)))
print(sample_summary)
summary_path <- paste0(file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/observational/"),"e11_landmark2y_sample_summary.csv")

cat("\n========================================\n")
cat("2-year landmark Cox analysis completed\n")
cat("========================================\n")
cat("Original landmark participants: ",nrow(data0),"\n")
cat("Participants after covariate merge: ",nrow(data),"\n")
cat("Participants excluded for missing values: ",n_missing_excluded,"\n")
cat("Final Cox participants: ",nrow(data2),"\n")
cat("E11 events: ",sum(data2$event == 1),"\n")
cat("Non-events/censored participants: ",sum(data2$event == 0),"\n")
cat("========================================\n")
