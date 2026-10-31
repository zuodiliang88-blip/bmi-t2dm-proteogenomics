library(data.table)

select_proteins <- function(result, hete, ple) {
  x <- fread(result)[method %in% c("Wald ratio", "Inverse variance weighted")]
  x[, fdr_p := p.adjust(pval, "fdr")]
  h <- fread(hete)[method == "Inverse variance weighted", .(id, Q_pval)]
  e <- fread(ple)[, .(id, ple_pval = pval)]
  x <- merge(merge(x, h, by="id", all.x=TRUE), e, by="id", all.x=TRUE)
  x[fdr_p < .05 & (is.na(Q_pval) | Q_pval >= .05) & (is.na(ple_pval) | ple_pval >= .05)]
}
a <- select_proteins(file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/mr/bmi_result"), file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/mr/bmi_hete"), file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/mr/bmi_ple"))
b <- select_proteins(file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/mr/t2dm_step2_result"), file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/mr/t2dm_step2_hete"), file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/mr/t2dm_step2_ple"))
A <- merge(a[, .(id, b_bmi=b)], b[, .(id, b_t2dm=b)], by="id")
A <- A[b_bmi * b_t2dm > 0]
stopifnot(nrow(A)==55L, !anyDuplicated(A$id))
fwrite(A, file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/mr/final_55_proteins.csv"))
cat("BMI selected", nrow(a), "T2DM selected", nrow(b), "overlap same direction", nrow(A), "\n")

expected <- fread(file.path(Sys.getenv("BMI_T2DM_ROOT"), "expected_dual_role_candidates.csv"))
stopifnot(setequal(A$id, expected$id))
