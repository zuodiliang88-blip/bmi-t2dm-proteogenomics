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
