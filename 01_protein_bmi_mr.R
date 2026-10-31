options(timeout = 15)
library(TwoSampleMR)
library(MRInstruments)
library(ieugwasr)
library(plinkbinr)
library(dplyr)
library(data.table)
library(RadialMR)
rm(list = ls())
gc()
set.seed(20261019)
dir.create(file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/mr"), showWarnings = FALSE, recursive = TRUE)

columns = c("id.exposure","id.outcome","outcome","exposure","method","nsnp","b","se","pval","lo_ci","up_ci","or","or_lci95","or_uci95")
MRResults = data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(MRResults) = columns
MRResults$id <- character(0)

hete <- c("id.exposure","id.outcome","outcome","exposure","method","Q","Q_df","Q_pval")
heterogeneity <- data.frame(matrix(nrow = 0, ncol = length(hete)))
colnames(heterogeneity) <- hete
heterogeneity$id <- character(0)

ply <- c("id.exposure","id.outcome","outcome","exposure","egger_intercept","se","pval")
pleiotropy <- data.frame(matrix(nrow = 0, ncol = length(ply)))
colnames(pleiotropy) <- ply
pleiotropy$id <- character(0)

trait1.trait2.dat_mr <- data.frame()

sig_files <- list.files(file.path(Sys.getenv("BMI_T2DM_ROOT"), "data/pqtl/protein_file_list/"), full.names = FALSE)
sig_files <- paste0(file.path(Sys.getenv("BMI_T2DM_ROOT"), "data/pqtl/ukb_cis_1mb/"), sub("\\.gz$", "", sig_files))

outcomedata <- read_outcome_data(
  filename = file.path(Sys.getenv("BMI_T2DM_ROOT"), "data/gwas/BMI_GCST90475156.tsv.gz"),
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "freq",
  pval_col = "p",
  samplesize_col = "n"
)

for (file_name in sig_files) {

  tryCatch({

    if (!file.exists(file_name)) {
      cat("File does not exist:", file_name, "\n")
      next
    }
    id <- paste0(basename(file_name), ".gz")
    trait1.exposure_raw <- fread(file_name)

    print(file_name)
    trait1.exposure_data <- read_exposure_data(file_name,
                                               sep = "\t",
                                               snp_col = "SNP",
                                               beta_col = "Beta",
                                               se_col = "Se",
                                               eaf_col = "Eaf",
                                               effect_allele_col = "Effect_allele",
                                               other_allele_col = "Other_allele",
                                               pval_col = "P",
                                               samplesize_col = "N",
                                               gene_col = "Protein"
    )

    trait1.exposure_data <- trait1.exposure_data[trait1.exposure_data$pval.exposure < 5e-8, ]

    if (nrow(trait1.exposure_data) == 0) next

    dat <- trait1.exposure_data[, c("SNP", "pval.exposure")]
    colnames(dat)[1:2] <- c("rsid", "pval")
    pval_raw <- as.numeric(trait1.exposure_raw$P[match(tolower(dat$rsid), tolower(trait1.exposure_raw$SNP))])
    abs_beta_se <- abs(trait1.exposure_data$beta.exposure / trait1.exposure_data$se.exposure)
    stopifnot(all(is.finite(pval_raw)), all(pval_raw >= 0), all(is.finite(abs_beta_se[pval_raw == 0])))
    clump_order <- order(pval_raw, ifelse(pval_raw == 0, -abs_beta_se, 0), seq_along(pval_raw))

    dat$pval[clump_order] <- seq_len(nrow(dat)) / (nrow(dat) + 1) * 5e-8
    expo_clump <- ld_clump(dat = dat,
                           plink_bin = "plink",
                           bfile = file.path(Sys.getenv("BMI_T2DM_ROOT"), "data/ld/EUR"))
    trait1.exposure_data <- trait1.exposure_data[trait1.exposure_data$SNP %in% expo_clump$rsid, ]
    if (nrow(trait1.exposure_data) == 0) next

    trait1.trait2.outcome_data <- outcomedata[outcomedata$SNP %in% trait1.exposure_data$SNP, ]

    trait1.trait2.dat <- harmonise_data(exposure_dat = trait1.exposure_data, outcome_dat = trait1.trait2.outcome_data)
    trait1.trait2.dat$id.outcome <- "BMI"
    trait1.trait2.dat<-steiger_filtering(trait1.trait2.dat)
    trait1.trait2.dat <- subset(trait1.trait2.dat, mr_keep == TRUE)

    if (nrow(trait1.trait2.dat) == 0) next
    if (nrow(trait1.trait2.dat) >= 2) {
      Q <- mr_heterogeneity(trait1.trait2.dat)
      Q_pval <- Q$Q_pval[Q$method == "Inverse variance weighted"][1]
      if (is.finite(Q_pval) && Q_pval < 0.05) {
        radial_ivw <- tryCatch({
          dat_mr_input <- tsmr_to_rmr_format(trait1.trait2.dat)
          ivw_radial(dat_mr_input, alpha = 0.05)
        }, error = function(e) {
          cat("Radial未完成，保留原SNP：", conditionMessage(e), "\n")
          NULL
        })
        if (!is.null(radial_ivw) && is.data.frame(radial_ivw$outliers)) {
          heterogeneous_snps <- radial_ivw$outliers$SNP[!is.na(radial_ivw$outliers$p.value) & radial_ivw$outliers$p.value < 0.05]
          heterogeneous_snps <- intersect(heterogeneous_snps, trait1.trait2.dat$SNP)
          if (length(heterogeneous_snps) > 0 && length(heterogeneous_snps) < nrow(trait1.trait2.dat)) {
            trait1.trait2.dat <- trait1.trait2.dat[!trait1.trait2.dat$SNP %in% heterogeneous_snps, ]
          }
        }
      }
    }

    if (nrow(trait1.trait2.dat) > 0) {
      trait1.trait2.dat$id <- id
      trait1.trait2.dat_mr <- rbind(trait1.trait2.dat_mr, trait1.trait2.dat)
    }

    nsnp <- nrow(trait1.trait2.dat)

    temp_mr <- NULL
    temp_hete <- NULL
    temp_pleio <- NULL

    if (nsnp >= 3) {
      temp_hete <- mr_heterogeneity(trait1.trait2.dat); temp_hete$id <- id
      temp_pleio <- mr_pleiotropy_test(trait1.trait2.dat); temp_pleio$id <- id

      temp_mr <- mr(trait1.trait2.dat)
      temp_mr <- TwoSampleMR::generate_odds_ratios(temp_mr)
      temp_mr$id <- id
    }

    if (nsnp == 2) {
      temp_mr <- mr(trait1.trait2.dat)
      temp_mr <- TwoSampleMR::generate_odds_ratios(temp_mr)
      temp_mr$id <- id

      temp_hete <- mr_heterogeneity(trait1.trait2.dat)
      temp_hete$id <- id
    }

    if (nsnp == 1) {
      temp_mr <- mr(trait1.trait2.dat)
      temp_mr <- TwoSampleMR::generate_odds_ratios(temp_mr)
      temp_mr$id <- id
    }

    if (!is.null(temp_mr)) MRResults <- rbind(MRResults, temp_mr)
    if (!is.null(temp_hete)) heterogeneity <- rbind(heterogeneity, temp_hete)
    if (!is.null(temp_pleio)) pleiotropy <- rbind(pleiotropy, temp_pleio)

  }, error = function(e) {
    cat("Error occurred while processing file:", file_name, "\n")
    cat("Error message:", conditionMessage(e), "\n")
  })
}

write.table(MRResults,file=file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/mr/bmi_result"),col.names = TRUE,row.names = FALSE,sep = "\t",quote = FALSE)
write.table(heterogeneity,file=file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/mr/bmi_hete"),col.names = TRUE,row.names = FALSE,sep = "\t",quote = FALSE)
write.table(pleiotropy,file=file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/mr/bmi_ple"),col.names = TRUE,row.names = FALSE,sep = "\t",quote = FALSE)
write.table(trait1.trait2.dat_mr,file=file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/mr/bmi_trait1_trait2_dat_mr"),col.names = TRUE,row.names = FALSE,sep = "\t",quote = FALSE)
