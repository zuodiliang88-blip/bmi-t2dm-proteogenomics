# Three-trait colocalization with the original trait-type helper included.
# Set BMI_T2DM_ROOT to the absolute project directory and prepare the input files.

library(dplyr)
library(data.table)
library(moloc)
rm(list = ls())
gc()
# Required study helper: moloc_typed.R
moloc_test_typed <- moloc::moloc_test
old_assignment <- quote(listData[[i]]$type <- ifelse("Ncases" %in% names(listData[[i]]), "cc", "quant"))
new_assignment <- quote(stopifnot(length(unique(listData[[i]]$type)) == 1L, unique(listData[[i]]$type) %in% c("quant", "cc")))
changed <- 0L
patch_type <- function(expr) {
  if (identical(expr, old_assignment)) { changed <<- changed + 1L; return(new_assignment) }
  if (is.call(expr)) for (j in seq_along(expr)) {
    if (is.call(expr[[j]])) expr[[j]] <- patch_type(expr[[j]])
  }
  expr
}
body(moloc_test_typed) <- patch_type(body(moloc_test_typed))
stopifnot(changed == 1L)
rm(old_assignment, new_assignment, changed, patch_type)
workspace_dir <- file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/mr/moloc_results")

bmi_path  <- file.path(Sys.getenv("BMI_T2DM_ROOT"), "data/gwas/BMI_GCST90475156.tsv.gz")
t2dm_path <- file.path(Sys.getenv("BMI_T2DM_ROOT"), "data/gwas/t2dm_EUR.txt.gz")

dir.create(file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/mr/moloc_results"), showWarnings = FALSE)

# Candidate proteins must be selected from the primary MR results before this step.
# Retain Wald ratio/IVW results with BH FDR < 0.05 and Q/Egger-intercept P >= 0.05
# where available, then retain proteins with concordant BMI and T2DM effect directions.
# Supply their original protein-file identifiers in the id column of final_55_proteins.csv.
A <- fread(file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/mr/final_55_proteins.csv"))
base_path <- file.path(Sys.getenv("BMI_T2DM_ROOT"), "data/pqtl/ukb_cis_1mb/")

file_names <- A$id
file_names <- sub("\\.gz$", "", file_names)

sig_files <- paste0(base_path, file_names)
sig_files <- unique(sig_files)
stopifnot(length(sig_files) == 55L, all(file.exists(sig_files)))

cat("Loading BMI and T2DM GWAS data...\n")
bmi_full  <- fread(bmi_path)
setnames(bmi_full, old=c("chr","bp","p","beta","se","A1","A2","freq","n"),
         new=c("CHR","BP","P","BETA","SE","A1","A2","EAF","N"), skip_absent=TRUE)

t2dm_full <- fread(t2dm_path)
setnames(t2dm_full, old=c("chr","bp","p","beta","se","A1","A2","freq","n"),
         new=c("CHR","BP","P","BETA","SE","A1","A2","EAF","N"), skip_absent=TRUE)

moloc_summary_list <- list()

for(  p_file  in sig_files) {

  p_name <- basename(p_file)
  cat("MOLOC START", p_name, as.character(Sys.time()), "\n")
  tryCatch({

    p_dat <- fread(p_file)
    setnames(p_dat, old=c("Chr","Pos","P","Beta","Se","Effect_allele","Other_allele","Eaf","N"),
             new=c("CHR","BP","P","BETA","SE","A1","A2","EAF","N"), skip_absent=TRUE)

    p_dat <- p_dat[is.finite(P) & P >= 0 & is.finite(BETA/SE) & is.finite(BP)]
    if (nrow(p_dat) == 0L) stop("No valid pQTL variants")
    top_snp_row <- p_dat[order(P, ifelse(P == 0, -abs(BETA/SE), 0), seq_len(nrow(p_dat)))][1]
    t_chr <- top_snp_row$CHR
    t_pos <- top_snp_row$BP
    top_snp_id <- top_snp_row$SNP
    cat("LEAD", top_snp_id, "CHR", t_chr, "BP", t_pos, "\n")

    p_reg   <- p_dat[CHR == t_chr & BP > (t_pos - 5e5) & BP < (t_pos + 5e5)]
    sameSNP <- Reduce(intersect, list(p_reg$SNP, bmi_full$SNP, t2dm_full$SNP))

    p_reg     <-  p_reg  [SNP %in% sameSNP]
    bmi_reg    <-  bmi_full[SNP %in% sameSNP]
    t2d_reg    <- t2dm_full[SNP %in% sameSNP]

    m1 <- merge(p_reg, bmi_reg, by="SNP", suffixes=c("_p", "_bmi"))
    m2 <- merge(m1, t2d_reg, by="SNP")

    setnames(m2, old=c("BETA","SE","A1","A2","P","N","EAF"),
             new=c("BETA_t2d","SE_t2d","A1_t2d","A2_t2d","P_t2d","N_t2d","EAF_t2d"), skip_absent=TRUE)

    if(nrow(m2) >= 50) {

      m2$flip_bmi <- 0
      m2$flip_bmi[m2$A1_p == m2$A1_bmi & m2$A2_p == m2$A2_bmi] <- 1
      m2$flip_bmi[m2$A1_p == m2$A2_bmi & m2$A2_p == m2$A1_bmi] <- -1

      m2$flip_t2d <- 0
      m2$flip_t2d[m2$A1_p == m2$A1_t2d & m2$A2_p == m2$A2_t2d] <- 1
      m2$flip_t2d[m2$A1_p == m2$A2_t2d & m2$A2_p == m2$A1_t2d] <- -1

      m_final <- m2[flip_bmi != 0 & flip_t2d != 0]

      if(nrow(m_final) >= 50) {

        maf_vec <- ifelse(m_final$EAF_p < 0.5, m_final$EAF_p, 1 - m_final$EAF_p)

        df_p <- data.frame(SNP=as.character(m_final$SNP), BETA=as.numeric(m_final$BETA_p),
                           SE=as.numeric(m_final$SE_p), N=as.numeric(m_final$N_p),
                           MAF=maf_vec, type="quant", stringsAsFactors=F)

        df_bmi <- data.frame(SNP=as.character(m_final$SNP), BETA=as.numeric(m_final$BETA_bmi * m_final$flip_bmi),
                             SE=as.numeric(m_final$SE_bmi), N=as.numeric(m_final$N_bmi),
                             MAF=maf_vec, type="quant", stringsAsFactors=F)

        df_t2d <- data.frame(SNP=as.character(m_final$SNP), BETA=as.numeric(m_final$BETA_t2d * m_final$flip_t2d),
                             SE=as.numeric(m_final$SE_t2d), N=as.numeric(m_final$N_t2d),
                             MAF=maf_vec, type="cc", stringsAsFactors=F)

        data_single <- list(pQTL=df_p, BMI=df_bmi, T2DM=df_t2d)

        stopifnot(all(df_p$type == "quant"), all(df_bmi$type == "quant"), all(df_t2d$type == "cc"))
        res_moloc <- moloc_test_typed(data_single, priors = c(1e-4, 1e-5, 1e-6))

        ppa_table <- res_moloc$priors_lkl_ppa

        if (!is.null(ppa_table)) {

          res_summary <- as.data.frame(t(ppa_table$PPA))
          colnames(res_summary) <- paste0("PPA_", rownames(ppa_table))

          pp_abc_val <- ppa_table["abc", "PPA"]
          res_summary$PP_abc_shared <- pp_abc_val

          res_summary$Protein <- p_name
          res_summary$n_snps   <- res_moloc$nsnps

          if(!is.null(res_moloc$best_snp)) {
            res_summary$Best_SNP <- tail(res_moloc$best_snp$best.snp.coloc, 1)
          }

          moloc_summary_list[[p_name]] <- res_summary
          cat(sprintf("   - moloc completed. Three-trait posterior (PPA_abc): %.4f | SNP: %s\n",
                      pp_abc_val, res_summary$Best_SNP))
        } else {
          cat("   - Error: ppa_table was unavailable\n")
        }
      } else {
        cat("   - Skipped: fewer than 50 aligned shared SNPs\n")
      }
    } else {
      cat("   - Skipped: fewer than 50 shared SNPs after merging\n")
    }

  }, error = function(e) {
    cat(sprintf("MOLOC ERROR %s: %s\n", p_name, e$message))
  })
}

if (length(moloc_summary_list) > 0) {
  final_moloc_df <- rbindlist(moloc_summary_list, fill=TRUE)
  stopifnot(!anyDuplicated(final_moloc_df$Protein),
            all(is.finite(final_moloc_df$PP_abc_shared)),
            all(abs(rowSums(as.matrix(final_moloc_df[, grep("^PPA_", names(final_moloc_df)), with=FALSE])) - 1) < 1e-6))

  dir.create(file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/mr/moloc_results"), showWarnings = FALSE)

  stopifnot(nrow(final_moloc_df) == length(sig_files))

  fwrite(final_moloc_df, file = file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/mr/moloc_results/bfd__ALL_MOLOC_SUMMARY.csv"), sep = ",")

  significant_moloc <- final_moloc_df[PP_abc_shared > 0.7]
  fwrite(significant_moloc, file = file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/mr/moloc_results/bfd__SIGNIFICANT_MOLOC_H4_07.csv"), sep = ",")

  cat("\n=======================================================\n")
  cat(sprintf("Analysis completed. Number of significant proteins: %d\n", nrow(significant_moloc)))
} else {
  cat("\nNo valid results were generated.\n")
}

stopifnot(length(moloc_summary_list) == length(sig_files))
cat("ALL MOLOC COMPLETED", length(moloc_summary_list), as.character(Sys.time()), "\n")
