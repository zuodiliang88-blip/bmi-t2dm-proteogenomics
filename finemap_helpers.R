library(Matrix)
ref_bim <- fread(file.path(Sys.getenv("BMI_T2DM_ROOT"), "data/ld/EUR.bim"), header = FALSE)
stopifnot(!anyDuplicated(ref_bim$V2))

align_to_ld <- function(dat, ld_matrix, snp_order, prefix) {
  stopifnot(identical(as.character(dat$SNP), as.character(snp_order)),
            nrow(dat) == nrow(ld_matrix), nrow(ld_matrix) == ncol(ld_matrix),
            !anyDuplicated(dat$SNP), all(is.finite(ld_matrix)))
  ref <- ref_bim[match(dat$SNP, V2)]
  freq <- fread(paste0(prefix, ".frq"))
  freq <- freq[match(dat$SNP, SNP)]
  a1 <- toupper(dat$A1); a2 <- toupper(dat$A2)
  r1 <- toupper(ref$V5); r2 <- toupper(ref$V6)
  stopifnot(!anyNA(ref$V2), identical(as.character(freq$SNP), as.character(dat$SNP)))
  rf <- as.numeric(freq$MAF)

  f_rev <- toupper(freq$A1) == r2 & toupper(freq$A2) == r1
  f_same <- toupper(freq$A1) == r1 & toupper(freq$A2) == r2
  stopifnot(all(f_same | f_rev))
  rf[f_rev] <- 1 - rf[f_rev]
  same <- a1 == r1 & a2 == r2
  rev <- a1 == r2 & a2 == r1
  ca1 <- chartr("ACGT", "TGCA", a1); ca2 <- chartr("ACGT", "TGCA", a2)
  snv <- nchar(a1) == 1L & nchar(a2) == 1L & a1 %in% c("A","C","G","T") & a2 %in% c("A","C","G","T")
  csame <- snv & ca1 == r1 & ca2 == r2
  crev <- snv & ca1 == r2 & ca2 == r1
  pal <- snv & ((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
                (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C"))
  sign <- rep(NA_real_, nrow(dat))
  sign[!pal & (same | csame)] <- 1
  sign[!pal & (rev | crev)] <- -1
  eaf <- if ("EAF" %in% names(dat)) as.numeric(dat$EAF) else rep(NA_real_, nrow(dat))
  dplus <- abs(eaf - rf); dminus <- abs((1 - eaf) - rf)

  inform <- pal & (same | rev) & is.finite(eaf) & eaf > 0 & eaf < 1 &
            is.finite(rf) & pmin(eaf, 1-eaf) < 0.42 & pmin(rf, 1-rf) < 0.42
  sign[inform & dplus <= 0.10 & dminus - dplus >= 0.10] <- 1
  sign[inform & dminus <= 0.10 & dplus - dminus >= 0.10] <- -1
  valid <- is.finite(dat$BETA) & is.finite(dat$SE) & dat$SE > 0
  keep <- !is.na(sign) & valid
  reason <- ifelse(!valid, "invalid_beta_se", ifelse(is.na(sign),
           ifelse(pal, "palindromic_unresolved", "allele_mismatch"),
           ifelse(sign == -1, "kept_Z_flipped", "kept_Z_same")))
  fwrite(data.table(SNP=dat$SNP, A1=dat$A1, A2=dat$A2, LD_A1=r1, LD_A2=r2,
                    EAF=eaf, LD_A1_FREQ=rf, sign=sign, keep=keep, reason=reason),
         paste0(prefix, "_allele_alignment.tsv"), sep="\t")
  cat("ALIGNMENT", basename(prefix), "input", nrow(dat), "retained", sum(keep),
      "flipped", sum(sign[keep] == -1), "removed", sum(!keep), "\n")
  print(table(reason))
  stopifnot(sum(keep) > 1L)
  dat <- copy(dat[keep])
  dat[, Z_LD := (BETA / SE) * sign[keep]]
  dat[, `:=`(LD_A1=r1[keep], LD_A2=r2[keep])]
  ld_matrix <- ld_matrix[keep, keep, drop=FALSE]
  stopifnot(all(is.finite(dat$Z_LD)), max(abs(ld_matrix-t(ld_matrix))) < 1e-8,
            max(abs(diag(ld_matrix)-1)) < 1e-8)
  z <- dat$Z_LD; names(z) <- dat$SNP
  fwrite(dat, paste0(prefix, "_aligned_input.tsv"), sep="\t")
  list(dat=dat, ld=ld_matrix, z=z)
}

native_cs <- function(fit, snps, filename) {
  sets <- fit[[1]][["Credible set"]][[2]]
  if (length(sets) == 0L) {
    ans <- data.table(CS_ID=integer(), SNP=character(), PIP=numeric())
  } else {
    ans <- rbindlist(lapply(seq_along(sets), function(i) {
      idx <- as.integer(sets[[i]])
      stopifnot(all(idx >= 1L & idx <= length(snps)))
      data.table(CS_ID=i, SNP=snps[idx], PIP=fit[[1]]$PIPs[idx])
    }))
  }
  fwrite(ans, filename, sep="\t")
  ans
}

check_model_space <- function(directory, pip, cap) {
  models <- readMM(file.path(directory, "post_locus_1_poi_gamma.mtx"))
  score <- scan(file.path(directory, "post_locus_1_poi_likeli.txt"), quiet=TRUE)
  stopifnot(nrow(models)==length(score), ncol(models)==length(pip), all(is.finite(score)))
  sizes <- Matrix::rowSums(models)
  w <- exp(score-max(score)); w <- w/sum(w)
  rebuilt <- as.numeric(Matrix::crossprod(models,w))
  stopifnot(max(sizes)<=cap, sum(pip)<=cap+1e-8, max(abs(rebuilt-pip))<1e-8)
  cat("VERIFIED cap", cap, "models", nrow(models), "max_model", max(sizes),
      "sumPIP", format(sum(pip),digits=15), "PIP_rebuild_error", max(abs(rebuilt-pip)), "\n")
  data.table(num_causal=cap, n_snps=length(pip), n_models=nrow(models),
             max_model_size=max(sizes), sum_PIP=sum(pip), pip_rebuild_error=max(abs(rebuilt-pip)))
}
