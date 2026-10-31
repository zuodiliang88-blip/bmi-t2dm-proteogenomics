library(CARMA)
library(data.table)
args <- commandArgs(TRUE)
root <- file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/finemapping")

defs <- Filter(function(x) {
  is.call(x) && identical(x[[1]],as.name('<-')) &&
    identical(x[[2]],as.name('credible.set.fun.improved'))
}, as.list(body(CARMA::CARMA)))
stopifnot(length(defs)==1L)
fun_expr <- defs[[1]][[3]]
cs_fun <- eval(fun_expr)
all_tables <- list(); set_tables <- list()
for (protein in c('ENTPD6','ADH1B')) {
  results <- file.path(root,protein,'numcausal5')
  tables <- list()
  for (trait in c('protein','BMI','T2DM')) {
    tag <- tolower(trait)
    obj <- readRDS(file.path(results,paste0(tag,'_CARMA_fit.rds')))
    prefix <- switch(trait,protein='protein_ld',BMI='BMI_ld',T2DM='T2DM_ld')
    dat <- fread(file.path(results,paste0(prefix,'_aligned_input.tsv')))
    snps <- fread(file.path(results,paste0(prefix,'.snplist')),header=FALSE)[[1]]
    ld <- as.matrix(fread(file.path(results,paste0(prefix,'.ld')),header=FALSE))
    idx <- match(obj$snp_order,snps)
    stopifnot(!anyNA(idx),identical(dat$SNP,obj$snp_order),obj$num.causal==5)
    ld <- ld[idx,idx,drop=FALSE]
    pip <- obj$fit[[1]]$PIPs
    cs99 <- cs_fun(pip,ld,rho=.99)
    stopifnot(identical(cs99[[2]],obj$fit[[1]][['Credible set']][[2]]))
    cs95 <- cs_fun(pip,ld,rho=.95)
    member <- integer(length(pip))
    if (length(cs95[[2]])) member[unique(unlist(cs95[[2]]))] <- 1L
    out <- dat[,.(SNP,CHR,BP,A1,A2)]
    out[, PIP:=pip];out[, CS95:=member]
    out[, `:=`(protein=protein,trait=trait)]
    tables[[trait]] <- out
    cat('EXPORTED',protein,trait,'N',nrow(out),'CS95_GROUPS',length(cs95[[2]]),
        'CS95_MEMBERS',sum(member),'\n',file=stderr())
    if (length(cs95[[2]])) {
      for (j in seq_along(cs95[[2]])) {
        ids <- cs95[[2]][[j]]
        set_tables[[length(set_tables)+1L]] <- data.table(protein=protein,trait=trait,
          CS_ID=j,SNP=obj$snp_order[ids],PIP=pip[ids])
      }
    }
  }
  all_tables[[protein]] <- rbindlist(tables)
}

cat(jsonlite::toJSON(list(data=rbindlist(all_tables),sets=rbindlist(set_tables)),
                    dataframe='rows',auto_unbox=TRUE,na='null',digits=NA))
