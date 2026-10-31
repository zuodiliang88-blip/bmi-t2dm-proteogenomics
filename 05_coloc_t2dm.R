# Protein/T2DM colocalization; T2DM is a case-control trait.
# Set BMI_T2DM_ROOT to the absolute project directory and prepare the input files.

library(dplyr)
library(data.table)
library(coloc)

rm(list = ls())
gc()

# Candidate proteins must be selected from the primary MR results before this step.
# Retain Wald ratio/IVW results with BH FDR < 0.05 and Q/Egger-intercept P >= 0.05
# where available, then retain proteins with concordant BMI and T2DM effect directions.
# Supply their original protein-file identifiers in the id column of final_55_proteins.csv.
A <- fread(file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/mr/final_55_proteins.csv"))
stopifnot(nrow(A) == 55L, !anyDuplicated(A$id))
base_path <- file.path(Sys.getenv("BMI_T2DM_ROOT"), "data/pqtl/ukb_cis_1mb/")

file_names <- A$id
file_names <- sub("\\.gz$", "", file_names)

sig_files <- paste0(base_path, file_names)
sig_files <- unique(sig_files)

gwas <- fread(file.path(Sys.getenv("BMI_T2DM_ROOT"), "data/gwas/t2dm_EUR.txt.gz"))
gwas$varbeta <- gwas$se^2
gwas <- subset(gwas, !duplicated(SNP))

gwas$z <- gwas$beta / gwas$se
gwas$p <- as.numeric(gwas$p)

gwas <- gwas %>% na.omit() %>% dplyr::filter(p >= 0)

all_results <- list()

for (file_path in sig_files) {
  tryCatch({
  cat("COLOC START", basename(file_path), as.character(Sys.time()), "\n")
  pqtl <- fread(file_path)
  pqtl <- pqtl[, c("Chr", "Pos", "SNP", "Effect_allele", "Other_allele", "Eaf", "Beta", "Se", "P", "N")]
  colnames(pqtl) <- c("CHR", "BP", "SNP", "A1", "A2", "freq", "beta", "se", "p", "n")
  pqtl$varbeta <- pqtl$se^2
  pqtl$z <- pqtl$beta / pqtl$se

  pqtl <- pqtl[is.finite(p) & p >= 0 & is.finite(z) & is.finite(BP)]
  if (nrow(pqtl) == 0) stop("No valid pQTL variants")
  lead <- pqtl[order(p, ifelse(p == 0, -abs(z), 0), seq_len(nrow(pqtl)))]
  leadchr <- lead$CHR[1]
  leadbp <- as.numeric(lead$BP[1])
  cat("LEAD", lead$SNP[1], "CHR", leadchr, "BP", leadbp, "P", lead$p[1],
      "absZ", abs(lead$z[1]), "\n")
  QTLdata <- pqtl[pqtl$CHR == leadchr & pqtl$BP > leadbp - 5e5 & pqtl$BP < leadbp + 5e5, ]
  QTLdata <- subset(QTLdata, !duplicated(SNP))
  QTLdata <- na.omit(QTLdata)

  sameSNP <- intersect(QTLdata$SNP, gwas$SNP)
  QTLdata <- QTLdata[QTLdata$SNP %in% sameSNP, ]
  GWASdata <- gwas[gwas$SNP %in% sameSNP, ]

  dat_merge <- merge(QTLdata, GWASdata, by = "SNP", suffixes = c("_pqtl", "_gwas"))
  dat_merge$label_a1a2 <- ifelse(dat_merge$A1_pqtl == dat_merge$A2_gwas, "yes", "no")
  dat_merge$label_a2a1 <- ifelse(dat_merge$A2_pqtl == dat_merge$A1_gwas, "yes", "no")
  dat_merge$A1_gwas2 <- ifelse(dat_merge$label_a1a2 == "yes" & dat_merge$label_a2a1 == "yes",
                               dat_merge$A1_pqtl, dat_merge$A1_gwas)
  dat_merge$A2_gwas2 <- ifelse(dat_merge$label_a1a2 == "yes" & dat_merge$label_a2a1 == "yes",
                               dat_merge$A2_pqtl, dat_merge$A2_gwas)
  dat_merge$beta_gwas2 <- ifelse(dat_merge$label_a1a2 == "yes" & dat_merge$label_a2a1 == "yes",
                                 -dat_merge$beta_gwas, dat_merge$beta_gwas)
  dat_harmonise <- dat_merge[!duplicated(dat_merge$SNP), ]

  dat_harmonise$MAF_pqtl <- ifelse(dat_harmonise$freq_pqtl < 0.5, dat_harmonise$freq_pqtl, 1 - dat_harmonise$freq_pqtl)
  dat_harmonise$MAF_gwas <- ifelse(dat_harmonise$freq_gwas < 0.5, dat_harmonise$freq_gwas, 1 - dat_harmonise$freq_gwas)

  dat_harmonise <- dat_harmonise[
    is.finite(beta_pqtl) & is.finite(varbeta_pqtl) & is.finite(MAF_pqtl) &
      is.finite(beta_gwas) & is.finite(varbeta_gwas) & is.finite(MAF_gwas) &
      varbeta_pqtl > 0 & varbeta_gwas > 0 &
      MAF_pqtl > 0 & MAF_pqtl < 1 &
      MAF_gwas > 0 & MAF_gwas < 1,
  ]
  if (nrow(dat_harmonise) == 0) stop("No valid shared variants")
  cat("SHARED_SNPS", nrow(dat_harmonise), "\n")

  pqtl_form<- list(snp =dat_harmonise$SNP,

                   beta = dat_harmonise$beta_pqtl,
                   varbeta = dat_harmonise$varbeta_pqtl,
                   MAF = dat_harmonise$MAF_pqtl,
                   N = dat_harmonise$n_pqtl,
                   type = "quant")
  gwas_form<- list(snp = dat_harmonise$SNP,

                   beta = dat_harmonise$beta_gwas2,
                   varbeta = dat_harmonise$varbeta_gwas,
                   type = "cc")

  coloc_result <- tryCatch({
    coloc.abf(dataset1 = pqtl_form, dataset2 = gwas_form)
  }, error = function(e) {
    message(paste("Error in file:", file_path, conditionMessage(e)))
    return(NULL)
  })

  if (!is.null(coloc_result)) {
    res_df <- as.data.frame(t(coloc_result$summary))
    res_df$file <- basename(file_path)
    all_results[[length(all_results) + 1]] <- res_df
    cat("COLOC FINISH", basename(file_path), "PP.H4", coloc_result$summary["PP.H4.abf"], "\n")
  }
  }, error = function(e) {
    message("Error in file: ", file_path, " : ", conditionMessage(e))
  })
}

final_results <- bind_rows(all_results)

write.table(final_results, file = file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/mr/t2dm_ukb_all_coloc_5e5.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
stopifnot(nrow(final_results) == length(sig_files), !anyDuplicated(final_results$file))
stopifnot(all(abs(rowSums(final_results[, c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")]) - 1) < 1e-8))
cat("ALL COLOC COMPLETED T2DM", nrow(final_results), as.character(Sys.time()), "\n")
