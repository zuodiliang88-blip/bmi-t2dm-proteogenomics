# Four-model linear regression; Model4 is additional to the three manuscript models.
# Set BMI_T2DM_ROOT to the absolute project directory and prepare the input files.

library(dplyr)
library(data.table)
library(broom)
library(purrr)

proteins <- fread(file.path(Sys.getenv("BMI_T2DM_ROOT"), "data/observational/2920proteins.csv"))%>%
  select(eid,ADH1B,AP3B1,APOE,ENTPD6,GRP,RBP1)
cov<-fread(file.path(Sys.getenv("BMI_T2DM_ROOT"), "data/observational/micecov.csv"))%>%
  select(-18)%>%
  dplyr::rename(TG_lip=TG)
head(cov)
data<-merge(proteins,cov,by="eid")
na_counts <- colSums(is.na(data))
na_cols <- na_counts[na_counts > 0]
data1<-data%>%
  filter(complete.cases(.))
cat_vars <- c("sex", "activity","Ethnicity", "smoking", "education","alcohol","healthy_diet","income")
data2 <- data1%>% mutate(across(all_of(cat_vars), as.factor))

protein_names <- colnames(data2)[2:7]
var1 <- c("age", "sex","Ethnicity")
var2 <- c("TDindex", "education")
var3 <- c( "smoking", "alcohol","healthy_diet")
var4 <- c("TG_lip", "HDL", "LDL", "SBP", "DBP")
# Models 1-3 correspond to the manuscript; Model4 is an additional adjustment model.
# Protein standardization and covariate imputation occur upstream.
cov_list <- list(
  Model1 = var1,
  Model2 = c(var1, var2),
  Model3 = c(var1, var2, var3),
  Model4 = c(var1, var2, var3, var4)
)
results_list <- list()
run_lm <- function(protein, covars, dat) {
  fm <- as.formula(
    paste("BMI ~", paste(c(protein, covars), collapse = " + "))
  )
  
  fit <- lm(fm, data = dat)
  
  broom::tidy(fit, conf.int = TRUE) %>%
    filter(term == protein) %>%
    mutate(
      protein = protein,
      n = nobs(fit),
      r.squared = summary(fit)$r.squared,
      adj.r.squared = summary(fit)$adj.r.squared
    ) %>%
    select(protein, term, estimate, conf.low, conf.high,
           std.error, statistic, p.value, n, r.squared, adj.r.squared)
}
result_list <- imap(cov_list, ~{
  map_dfr(protein_names, run_lm, covars = .x, dat = data2) %>%
    mutate(model = .y, .before = 1) %>%
    arrange(protein)
})
final_result <-bind_rows(result_list)
fwrite(final_result, file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/observational/protein_bmi_linear_regression_results.csv"))
