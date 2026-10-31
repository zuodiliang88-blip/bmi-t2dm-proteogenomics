library(dplyr)
library(data.table)
library(CARMA)

rm(list = ls())
setwd(file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/finemapping/ADH1B/numcausal5"))
set.seed(20261018)
source(file.path(Sys.getenv("BMI_T2DM_ROOT"), "CARMA_strict.R"))
source(file.path(Sys.getenv("BMI_T2DM_ROOT"), "finemap_helpers.R"))

bmi_path  <- file.path(Sys.getenv("BMI_T2DM_ROOT"), "data/gwas/BMI_GCST90475156.tsv.gz")
t2dm_path <- file.path(Sys.getenv("BMI_T2DM_ROOT"), "data/gwas/t2dm_EUR.txt.gz")
p_file    <- file.path(Sys.getenv("BMI_T2DM_ROOT"), "data/pqtl/ukb_cis_1mb/adh1b_p00325_oid30612_v1_1MB.txt")

plink_bin <- "plink"
ref_bfile <- file.path(Sys.getenv("BMI_T2DM_ROOT"), "data/ld/EUR")

outdir <- file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/finemapping/ADH1B/numcausal5")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

bmi_full <- fread(bmi_path)
setnames(
  bmi_full,
  old = c("chr","bp","p","beta","se","A1","A2","freq","n"),
  new = c("CHR","BP","P","BETA","SE","A1","A2","EAF","N"),
  skip_absent = TRUE
)

t2dm_full <- fread(t2dm_path)
setnames(
  t2dm_full,
  old = c("chr","bp","p","beta","se","A1","A2","freq","n"),
  new = c("CHR","BP","P","BETA","SE","A1","A2","EAF","N"),
  skip_absent = TRUE
)

p_dat <- fread(p_file)
setnames(
  p_dat,
  old = c("Chr","Pos","P","Beta","Se","Effect_allele","Other_allele","Eaf","N"),
  new = c("CHR","BP","P","BETA","SE","A1","A2","EAF","N"),
  skip_absent = TRUE
)

for (dt in list(bmi_full, t2dm_full, p_dat)) {
  if ("BP" %in% names(dt)) dt[, BP := as.numeric(BP)]
  if ("P" %in% names(dt)) dt[, P := as.numeric(P)]
  if ("BETA" %in% names(dt)) dt[, BETA := as.numeric(BETA)]
  if ("SE" %in% names(dt)) dt[, SE := as.numeric(SE)]
  if ("EAF" %in% names(dt)) dt[, EAF := as.numeric(EAF)]
  if ("N" %in% names(dt)) dt[, N := as.numeric(N)]
}

bmi_full  <- unique(bmi_full, by = "SNP")
t2dm_full <- unique(t2dm_full, by = "SNP")
p_dat     <- unique(p_dat, by = "SNP")

target_snp <- "rs1229984"
window_size <- 250000

protein_anchor <- p_dat[SNP == target_snp]
bmi_anchor     <- bmi_full[SNP == target_snp]
t2dm_anchor    <- t2dm_full[SNP == target_snp]

cat("Protein target SNP:", target_snp, " CHR:", protein_anchor$CHR[1], " BP:", protein_anchor$BP[1], "\n")
cat("BMI target SNP:", target_snp, " CHR:", bmi_anchor$CHR[1], " BP:", bmi_anchor$BP[1], "\n")
cat("T2DM target SNP:", target_snp, " CHR:", t2dm_anchor$CHR[1], " BP:", t2dm_anchor$BP[1], "\n")

p_chr <- protein_anchor$CHR[1]
p_pos <- protein_anchor$BP[1]

b_chr <- bmi_anchor$CHR[1]
b_pos <- bmi_anchor$BP[1]

t_chr <- t2dm_anchor$CHR[1]
t_pos <- t2dm_anchor$BP[1]

p_reg <- p_dat[
  CHR == p_chr &
    BP > (p_pos - window_size) &
    BP < (p_pos + window_size)
]

bmi_reg <- bmi_full[
  CHR == b_chr &
    BP > (b_pos - window_size) &
    BP < (b_pos + window_size)
]

t2d_reg <- t2dm_full[
  CHR == t_chr &
    BP > (t_pos - window_size) &
    BP < (t_pos + window_size)
]

p_reg   <- unique(p_reg, by = "SNP")
bmi_reg <- unique(bmi_reg, by = "SNP")
t2d_reg <- unique(t2d_reg, by = "SNP")

cat("Protein SNP count:", nrow(p_reg), "\n")
cat("BMI SNP count:", nrow(bmi_reg), "\n")
cat("T2DM SNP count:", nrow(t2d_reg), "\n")

run_carma_one <- function(z_vec, ld_matrix, snp_vec, num_causal = 5) {

  n_snps <- length(z_vec)
  stopifnot(n_snps == nrow(ld_matrix), n_snps == ncol(ld_matrix),
            n_snps == length(snp_vec), identical(names(z_vec), snp_vec))

  z.list <- list(as.numeric(z_vec[1:n_snps]))
  ld.list <- list(ld_matrix[1:n_snps, 1:n_snps, drop = FALSE])
  lambda.list <- list(1 / sqrt(n_snps))

  fit <- CARMA_strict(
    z.list = z.list,
    ld.list = ld.list,
    lambda.list = lambda.list,
    effect.size.prior = "Spike-slab",
    num.causal = num_causal,
    outlier.switch = TRUE,
    printing.log = FALSE,
    output.labels = file.path(outdir, paste0("CARMA_native_", carma_trait))
  )

  pips <- fit[[1]]$PIPs
  verified <- check_model_space(file.path(outdir, paste0("CARMA_native_", carma_trait)), pips, num_causal)
  fwrite(verified, file.path(outdir, paste0(carma_trait, "_model_check.tsv")), sep="\t")

  res <- data.table(
    SNP = snp_vec[1:n_snps],
    PIP = pips
  )

  return(list(
    fit = fit,
    pip = res
  ))
}

protein_snp_file  <- file.path(outdir, "protein_snps.txt")
protein_ld_prefix <- file.path(outdir, "protein_ld")

write.table(
  p_reg$SNP,
  protein_snp_file,
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)

system2(
  plink_bin,
  args = c(
    "--bfile", ref_bfile,
    "--keep-allele-order",
    "--extract", protein_snp_file,
    "--r", "square",
    "--freq",
    "--write-snplist",
    "--out", protein_ld_prefix
  )
)

protein_ld_matrix <- as.matrix(fread(paste0(protein_ld_prefix, ".ld"), header = FALSE))
protein_ld_snps   <- fread(paste0(protein_ld_prefix, ".snplist"), header = FALSE)[[1]]

p_reg <- p_reg[SNP %in% protein_ld_snps]
p_reg <- p_reg[match(protein_ld_snps, SNP)]
p_reg <- p_reg[!is.na(SNP)]

protein_keep_idx <- match(p_reg$SNP, protein_ld_snps)
protein_ld_matrix <- protein_ld_matrix[protein_keep_idx, protein_keep_idx, drop = FALSE]

protein_ld_matrix <- (protein_ld_matrix + t(protein_ld_matrix)) / 2
diag(protein_ld_matrix) <- 1

aligned <- align_to_ld(p_reg, protein_ld_matrix, p_reg$SNP, protein_ld_prefix)
p_reg <- aligned$dat
protein_ld_matrix <- aligned$ld
protein_z <- aligned$z
stopifnot(target_snp %in% p_reg$SNP)

carma_trait <- "protein"
cat("CARMA START protein", as.character(Sys.time()), "\n")
carma_protein <- run_carma_one(
  z_vec = protein_z,
  ld_matrix = protein_ld_matrix,
  snp_vec = p_reg$SNP,
  num_causal = 5
)
saveRDS(list(fit = carma_protein$fit, snp_order = p_reg$SNP, z = protein_z, num.causal = 5, seed = 20261018), file.path(outdir, "protein_CARMA_fit.rds"))
cat("CARMA FINISH protein", as.character(Sys.time()), "\n")

protein_pip <- copy(carma_protein$pip)
setnames(protein_pip, "PIP", "PIP_protein")

protein_result <- data.table(
  SNP  = p_reg$SNP,
  CHR  = p_reg$CHR,
  BP   = p_reg$BP,
  A1   = p_reg$A1,
  A2   = p_reg$A2,
  P    = p_reg$P,
  BETA = p_reg$BETA,
  Z_LD = p_reg$Z_LD,
  LD_A1 = p_reg$LD_A1, LD_A2 = p_reg$LD_A2,
  SE   = p_reg$SE
)

protein_result <- merge(protein_result, protein_pip, by = "SNP", all.x = TRUE)
setorder(protein_result, -PIP_protein)
protein_result <- protein_result[!is.na(PIP_protein)]
cs_table <- native_cs(carma_protein$fit, p_reg$SNP,
                      file.path(outdir, "protein_CARMA_native_CS99.tsv"))
protein_result[, CS_ID := vapply(SNP, function(snp) {
  ids <- cs_table[SNP == snp, CS_ID]
  if (length(ids)) paste(ids, collapse=";") else NA_character_
}, character(1))]

fwrite(
  protein_result,
  file.path(outdir, "protein_CARMA_result.txt"),
  sep = "\t"
)

bmi_snp_file  <- file.path(outdir, "BMI_snps.txt")
bmi_ld_prefix <- file.path(outdir, "BMI_ld")

write.table(
  bmi_reg$SNP,
  bmi_snp_file,
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)

system2(
  plink_bin,
  args = c(
    "--bfile", ref_bfile,
    "--keep-allele-order",
    "--extract", bmi_snp_file,
    "--r", "square",
    "--freq",
    "--write-snplist",
    "--out", bmi_ld_prefix
  )
)

bmi_ld_matrix <- as.matrix(fread(paste0(bmi_ld_prefix, ".ld"), header = FALSE))
bmi_ld_snps   <- fread(paste0(bmi_ld_prefix, ".snplist"), header = FALSE)[[1]]

bmi_reg <- bmi_reg[SNP %in% bmi_ld_snps]
bmi_reg <- bmi_reg[match(bmi_ld_snps, SNP)]
bmi_reg <- bmi_reg[!is.na(SNP)]

bmi_keep_idx <- match(bmi_reg$SNP, bmi_ld_snps)
bmi_ld_matrix <- bmi_ld_matrix[bmi_keep_idx, bmi_keep_idx, drop = FALSE]

bmi_ld_matrix <- (bmi_ld_matrix + t(bmi_ld_matrix)) / 2
diag(bmi_ld_matrix) <- 1

aligned <- align_to_ld(bmi_reg, bmi_ld_matrix, bmi_reg$SNP, bmi_ld_prefix)
bmi_reg <- aligned$dat
bmi_ld_matrix <- aligned$ld
bmi_z <- aligned$z
stopifnot(target_snp %in% bmi_reg$SNP)

carma_trait <- "bmi"
cat("CARMA START bmi", as.character(Sys.time()), "\n")
carma_bmi <- run_carma_one(
  z_vec = bmi_z,
  ld_matrix = bmi_ld_matrix,
  snp_vec = bmi_reg$SNP,
  num_causal = 5
)
saveRDS(list(fit = carma_bmi$fit, snp_order = bmi_reg$SNP, z = bmi_z, num.causal = 5, seed = 20261018), file.path(outdir, "bmi_CARMA_fit.rds"))
cat("CARMA FINISH bmi", as.character(Sys.time()), "\n")

bmi_pip <- copy(carma_bmi$pip)
setnames(bmi_pip, "PIP", "PIP_BMI")

bmi_result <- data.table(
  SNP  = bmi_reg$SNP,
  CHR  = bmi_reg$CHR,
  BP   = bmi_reg$BP,
  A1   = bmi_reg$A1,
  A2   = bmi_reg$A2,
  P    = bmi_reg$P,
  BETA = bmi_reg$BETA,
  Z_LD = bmi_reg$Z_LD,
  LD_A1 = bmi_reg$LD_A1, LD_A2 = bmi_reg$LD_A2,
  SE   = bmi_reg$SE
)

bmi_result <- merge(bmi_result, bmi_pip, by = "SNP", all.x = TRUE)
setorder(bmi_result, -PIP_BMI)
bmi_result <- bmi_result[!is.na(PIP_BMI)]
cs_table <- native_cs(carma_bmi$fit, bmi_reg$SNP,
                      file.path(outdir, "BMI_CARMA_native_CS99.tsv"))
bmi_result[, CS_ID := vapply(SNP, function(snp) {
  ids <- cs_table[SNP == snp, CS_ID]
  if (length(ids)) paste(ids, collapse=";") else NA_character_
}, character(1))]

fwrite(
  bmi_result,
  file.path(outdir, "BMI_CARMA_result.txt"),
  sep = "\t"
)

t2dm_snp_file  <- file.path(outdir, "T2DM_snps.txt")
t2dm_ld_prefix <- file.path(outdir, "T2DM_ld")

write.table(
  t2d_reg$SNP,
  t2dm_snp_file,
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)

system2(
  plink_bin,
  args = c(
    "--bfile", ref_bfile,
    "--keep-allele-order",
    "--extract", t2dm_snp_file,
    "--r", "square",
    "--freq",
    "--write-snplist",
    "--out", t2dm_ld_prefix
  )
)

t2dm_ld_matrix <- as.matrix(fread(paste0(t2dm_ld_prefix, ".ld"), header = FALSE))
t2dm_ld_snps   <- fread(paste0(t2dm_ld_prefix, ".snplist"), header = FALSE)[[1]]

t2d_reg <- t2d_reg[SNP %in% t2dm_ld_snps]
t2d_reg <- t2d_reg[match(t2dm_ld_snps, SNP)]
t2d_reg <- t2d_reg[!is.na(SNP)]

t2dm_keep_idx <- match(t2d_reg$SNP, t2dm_ld_snps)
t2dm_ld_matrix <- t2dm_ld_matrix[t2dm_keep_idx, t2dm_keep_idx, drop = FALSE]

t2dm_ld_matrix <- (t2dm_ld_matrix + t(t2dm_ld_matrix)) / 2
diag(t2dm_ld_matrix) <- 1

aligned <- align_to_ld(t2d_reg, t2dm_ld_matrix, t2d_reg$SNP, t2dm_ld_prefix)
t2d_reg <- aligned$dat
t2dm_ld_matrix <- aligned$ld
t2dm_z <- aligned$z
stopifnot(target_snp %in% t2d_reg$SNP)

carma_trait <- "t2dm"
cat("CARMA START t2dm", as.character(Sys.time()), "\n")
carma_t2dm <- run_carma_one(
  z_vec = t2dm_z,
  ld_matrix = t2dm_ld_matrix,
  snp_vec = t2d_reg$SNP,
  num_causal = 5
)
saveRDS(list(fit = carma_t2dm$fit, snp_order = t2d_reg$SNP, z = t2dm_z, num.causal = 5, seed = 20261018), file.path(outdir, "t2dm_CARMA_fit.rds"))
cat("CARMA FINISH t2dm", as.character(Sys.time()), "\n")

t2dm_pip <- copy(carma_t2dm$pip)
setnames(t2dm_pip, "PIP", "PIP_T2DM")

t2dm_result <- data.table(
  SNP  = t2d_reg$SNP,
  CHR  = t2d_reg$CHR,
  BP   = t2d_reg$BP,
  A1   = t2d_reg$A1,
  A2   = t2d_reg$A2,
  P    = t2d_reg$P,
  BETA = t2d_reg$BETA,
  Z_LD = t2d_reg$Z_LD,
  LD_A1 = t2d_reg$LD_A1, LD_A2 = t2d_reg$LD_A2,
  SE   = t2d_reg$SE
)

t2dm_result <- merge(t2dm_result, t2dm_pip, by = "SNP", all.x = TRUE)
setorder(t2dm_result, -PIP_T2DM)
t2dm_result <- t2dm_result[!is.na(PIP_T2DM)]
cs_table <- native_cs(carma_t2dm$fit, t2d_reg$SNP,
                      file.path(outdir, "T2DM_CARMA_native_CS99.tsv"))
t2dm_result[, CS_ID := vapply(SNP, function(snp) {
  ids <- cs_table[SNP == snp, CS_ID]
  if (length(ids)) paste(ids, collapse=";") else NA_character_
}, character(1))]

fwrite(
  t2dm_result,
  file.path(outdir, "T2DM_CARMA_result.txt"),
  sep = "\t"
)

cat("All done.\n")

outdir <- file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/finemapping/ADH1B/numcausal5")

protein_result <- fread(file.path(outdir, "protein_CARMA_result.txt"))
bmi_result     <- fread(file.path(outdir, "BMI_CARMA_result.txt"))
t2dm_result    <- fread(file.path(outdir, "T2DM_CARMA_result.txt"))

merged_result <- merge(
  protein_result[, .(SNP, PIP_protein, CS_protein=CS_ID)],
  bmi_result[, .(SNP, PIP_BMI, CS_BMI=CS_ID)], by="SNP", all=TRUE
)
merged_result <- merge(merged_result,
  t2dm_result[, .(SNP, PIP_T2DM, CS_T2DM=CS_ID)], by="SNP", all=TRUE)
setorder(merged_result, -PIP_protein)
fwrite(merged_result, file.path(outdir,"merged_CARMA_result.txt"), sep="\t", na="NA")
writeLines(capture.output(sessionInfo()), file.path(outdir,"sessionInfo.txt"))
cat("RUN COMPLETED", as.character(Sys.time()), "\n")
