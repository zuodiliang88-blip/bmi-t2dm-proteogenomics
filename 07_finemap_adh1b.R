# One CARMA example at ADH1B; original strict implementation and LD helpers included.
# Set BMI_T2DM_ROOT to the absolute project directory and prepare the input files.

library(dplyr)
library(data.table)
library(CARMA)

rm(list = ls())
setwd(file.path(Sys.getenv("BMI_T2DM_ROOT"), "results/finemapping/ADH1B/numcausal5"))
set.seed(20261018)
# Required study helper: CARMA_strict.R
# Modified CARMA implementation by the study authors; GPL-3.0-or-later.
# Original CARMA: Zikun Yang, Chen Wang and Iuliana Ionita-Laza.
# Original authors: Zikun Yang, Chen Wang and Iuliana Ionita-Laza.
# See LICENSE-CARMA for the original licence.
CARMA_strict<-function(z.list,ld.list,w.list=NULL,lambda.list=NULL,output.labels='.',label.list=NULL,
                                 effect.size.prior='Spike-slab',rho.index=0.99,BF.index=10,EM.dist='Logistic',
                                Max.Model.Dim=2e+5,all.iter=3,all.inner.iter=10,input.alpha=0,epsilon.threshold=1e-5,printing.log=F,
                                 num.causal=10,y.var=1,tau=0.04,outlier.switch=T,outlier.BF.index=1/3.2,prior.prob.computation='Logistic'){

  Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

  {
    log.2pi<-log(2*pi)
    L<-length(z.list)
    p.list<-list()
    for(i in 1:L){
       z.list[[i]]<-as.matrix(z.list[[i]])
      p.list[[i]]<-nrow(z.list[[i]])
    }
    B<-Max.Model.Dim
    all.B.list<-list()
    for(i in 1:L){
    all.B.list[[i]]<-list()
    all.B.list[[i]][[1]]<-integer(0)
    all.B.list[[i]][[2]]<-Matrix(nrow = 0,ncol=p.list[[i]],data=0,sparse = T)
    }
    q.list<-list()
    if(!is.null(w.list)){
      for(i in 1:L){

        q.list[[i]]<-ncol(w.list[[i]])
        invariant.var.index<-which((apply(as.matrix(w.list[[i]][,-1]),2,sd))==0)
        if(length(invariant.var.index)!=0){
          invariant.var<-w.list[[i]][,invariant.var.index+1]
          w.list[[i]]<-as.matrix(cbind(1,scale(w.list[[i]][,-1])))
          w.list[[i]][,invariant.var.index+1]<-invariant.var
        }else{
          w.list[[i]]<-as.matrix(cbind(1,scale(w.list[[i]][,-1])))
        }

      }
    }
    if(is.null(label.list)){
      for(i in 1:L){
        label.list[[i]]=paste0('locus_',i)
      }
    }
    Sigma.list<-list()
    for(i in 1:L){
    Sigma.list[[i]]<-as.matrix(ld.list[[i]])
    }
    S.list<-list()
    for(i in 1:L){
    S.list[[i]]<-integer(0)
    }

    all.C.list<-list()
    for(i in 1:L){
      all.C.list[[i]]<-list()
      all.C.list[[i]][[1]]<-integer(0)
      all.C.list[[i]][[2]]<-Matrix(nrow = 0,ncol=p.list[[i]],data=0,sparse = T)
    }

    all.epsilon.threshold<-0
    epsilon.list<-list()
    for(i in 1:L){
      epsilon.list[[i]]<-epsilon.threshold*p.list[[i]]
      all.epsilon.threshold<-all.epsilon.threshold+  epsilon.threshold*p.list[[i]]
    }
    model.prior='Poisson'
    standardize.model.space=T

  }

  EM.M.step.func<-function(input.response=NULL,w=w,input.alpha=0.5,EM.dist='Logistic'){
    if(EM.dist=='Poisson'){
      count.index<-input.response
      cv.poisson<-cv.glmnet(w,count.index,family = 'poisson',alpha=input.alpha,type.measure='deviance' )
      cv.index<-which(cv.poisson$lambda==cv.poisson$lambda.min)
      glm.beta<-as.matrix(c(cv.poisson$glmnet.fit$a0[cv.index],cv.poisson$glmnet.fit$beta[-1,cv.index]))
    }
    if(EM.dist=='Logistic'){
      response.matrix<-matrix(c(1-input.response, input.response),length(input.response),2)
      cv.logistic<-cv.glmnet(x=w,y=response.matrix, family='binomial',alpha=input.alpha,type.measure = 'deviance')
      cv.index<-which(cv.logistic$lambda==cv.logistic$lambda.min)
      glm.beta<-as.matrix(c(cv.logistic$glmnet.fit$a0[cv.index],cv.logistic$glmnet.fit$beta[-1,cv.index]))

    }

    return(glm.beta=glm.beta)
  }

  credible.set.fun.improved<-function(pip,ld,true.beta=NULL,rho=0.99){
    candidate.r<-c((seq(from=.5,0.95,by=0.05)),seq(0.96,0.99,0.01))

    snp.list<-list()
    colnames(ld)<-rownames(ld)<-1:nrow(ld)
    for(r in 1:length(candidate.r)){
      working.ld<-ld
      cor.threshold<-candidate.r[r]
      pip.order<-order(pip,decreasing = T)
      snp.list[[r]]<-list()
      group.index<-1
      s<-1
      while(sum(pip[pip.order[s:length(pip.order)]])>rho){

        cor.group<-as.numeric(names(which(abs(working.ld[which(pip.order[s]==as.numeric(colnames(working.ld))),])>cor.threshold)))
        if(sum(pip[cor.group])>rho){
          group.pip<- pip[cor.group]
          snp.index<- cor.group[order(group.pip,decreasing = T)][1: min(which(cumsum( sort(group.pip,decreasing = T))>rho))]
          snp.list[[r]][[group.index]]<-snp.index
          group.index<-group.index+1
          pip.order<-pip.order[-match(snp.index,pip.order)]
          working.ld<-working.ld[-match(snp.index,as.numeric(colnames(working.ld))),-match(snp.index,as.numeric(colnames(working.ld)))]
        }else{
          s=s+1
        }
      }

    }
    if(sum(sapply(snp.list,length))!=0){
      group.index<-max(which(sapply(snp.list,length)==max(sapply(snp.list,length))))
      credible.set.list<-snp.list[[group.index]]
      if(!is.null(true.beta)){
        purity<-c()
        for(s in 1:length(credible.set.list)){
          purity<-c(purity, (mean(ld[ credible.set.list[[s]], credible.set.list[[s]]]^2)))
        }
        causal.snp<-ifelse(length(na.omit(match(unlist(credible.set.list),true.beta)))!=0,
                           length(na.omit(match(unlist(credible.set.list),true.beta))),0)
        length.credible<-length(credible.set.list)
        return(list(c(causal.snp,
                      ifelse(length.credible==0,NA,length.credible),
                      mean(sapply(credible.set.list,length)),
                      mean(purity)),credible.set.list))

      }else{
        purity<-c()
        for(s in 1:length(credible.set.list)){
          purity<-c(purity, (mean(ld[ credible.set.list[[s]], credible.set.list[[s]]]^2)))
        }
        length.credible<-length(credible.set.list)
        return(list(c(ifelse(length.credible==0,NA,length.credible),
                      mean(sapply(credible.set.list,length)),
                      mean(purity)),
                    credible.set.list))
      }
    }else{
      return(list(rep(0,4),list()))
    }
  }

  credible.model.fun<-function(likelihood,model.space,bayes.threshold=10){
    na.index<-which(is.na(likelihood))
    if(length(na.index)!=0){
      likelihood<-likelihood[-na.index]
      model.space<-model.space[-na.index,]
    }
    post.like.temp<-likelihood-likelihood[1]
    post.prob<-exp(post.like.temp)/(sum(exp(post.like.temp)))
    bayes.f<-post.prob[1]/post.prob
    candidate.model<-1

    credible.model.list<-list()
    credible.model.list[[1]]<-list()
    input.rs<-c()
    for(ss in 1:length(which(bayes.f<bayes.threshold))){
      credible.model.list[[1]][[ss]]<-which(model.space[ss,]==1)
      input.rs<-c(input.rs,which(model.space[ss,]==1))
    }
    credible.model.list[[2]]<-data.frame(Posterior.Prob=post.prob[which(bayes.f<bayes.threshold)])
    credible.model.list[[3]]<-unique(input.rs)
    return(credible.model.list )

  }

  Module.Cauchy.Shotgun<-function(z,ld.matrix,Max.Model.Dim=1e+4,input.S=NULL,lambda,label,printing.log=F,
                                              num.causal=10,output.labels,y.var=1,effect.size.prior=effect.size.prior,model.prior=model.prior,
                                              outlier.switch,input.conditional.S.list=NULL,tau=1/0.05^2,
                                              C.list=NULL,prior.prob=NULL,epsilon=1e-3,inner.all.iter=10){
       {

         prob.list<-list()
         p<-nrow(z);
         if(model.prior=='input.prob'){
            posi.log.pro<-log(prior.prob)
            nega.log.pro<-log(1-prior.prob)
            input.prior.dist<-function(x){
              variable.index<-which(x==1)
              if(any(variable.index)){
                return( sum(posi.log.pro[variable.index])+sum(nega.log.pro[(1:p)[-variable.index]])-sum(nega.log.pro))
              }else{
                return(sum(nega.log.pro))
              }
            }
            prior.dist<-input.prior.dist
          }

       if(model.prior=='Poisson'){
        Poisson.prior.dist<-function(t){
          dim.model<-sum(t)
          result<-dim.model*log(lambda)+lfactorial(p-dim.model)-lfactorial(p)
          return(result)
        }
        prior.dist<-Poisson.prior.dist
       }
      if(model.prior=='beta-binomial'){
        beta.binomial.dist<-function(t){
          dim.model<-sum(t)
          result<- lbeta(dim.model+1,p-dim.model+9)-lbeta(1,p+9)
          return(result)
        }
        prior.dist<-beta.binomial.dist
      }

         if(effect.size.prior=='Cauchy'){
           marginal_likelihood=Cauchy_fixed_sigma_marginal
           tau.sample<-rgamma(1e+5,0.5,rate=0.5)
           if(outlier.switch){
             outlier_likelihood=outlier_Cauchy_fixed_sigma_marginal
             outlier.tau=tau.sample
           }
         }
         if(effect.size.prior=='Hyper-g'){
           marginal_likelihood=hyper_g_fixed_sigma_marginal
           tau.sample<-rgamma(1e+5,0.5,rate=0.5)
         }
         if(effect.size.prior=='Normal'){
           marginal_likelihood=Normal_fixed_sigma_marginal
           tau.sample<-tau
           if(outlier.switch){
             outlier_likelihood=outlier_Normal_fixed_sigma_marginal
             outlier.tau=tau.sample
           }
         }
         if(effect.size.prior=='Spike-slab'){
           marginal_likelihood=ind_Normal_fixed_sigma_marginal
           tau.sample<-tau
           if(outlier.switch){

               outlier_likelihood=outlier_ind_Normal_marginal
               outlier.tau=tau.sample

           }
         }

         p<-nrow(z);
      log.2pi<-log(2*pi)

      B<-Max.Model.Dim
      stored.result.prob<-rep(0,p)
      stored.bf<-0
      Sigma<-as.matrix(ld.matrix)

      if(!is.null(input.S)){
        S<-input.S
      }else{
        S<-integer(0)
      }
      conditional.S=NULL
      null.model<-Matrix(nrow = 1,ncol=p,data=0,sparse = T)
      null.margin<-prior.dist(null.model)
      if(is.null(C.list)){
        C.list<-list()
        C.list[[2]]<-list()
        C.list[[1]]<-list()
      }
           B.list<-list()
           B.list[[1]]<-prior.dist(null.model)
           B.list[[2]]<-Matrix(nrow = 1,ncol=p,data=0,sparse = T)
         if(length(input.conditional.S.list)==0){
           conditional.S.list<-list()
           conditional.S=NULL
         }else{
           conditional.S.list<-input.conditional.S.list
           conditional.S<-input.conditional.S.list$Index
           conditional.S<-unique(conditional.S)
           S<-conditional.S
         }
    }

    set.gamma.func<-function(input.S,condition.index=NULL){
       set.gamma.func.base<-function(S){
         add.function<-function(y){results<-(apply(as.matrix(S_sub),1,function(x){return(sort(c(x,y)))}))
         return(t(results))
         }

         set.gamma<-list()
         for(i in 1:3){
           set.gamma[[i]]<-c()
         }

         if(length(S)==0){
           S_sub<-(1:p)
           set.gamma[[1]]<-c()
           set.gamma[[2]]<-apply(as.matrix(S_sub),1,function(x){return(c(x,S))})
           set.gamma[[2]]<-as.matrix(set.gamma[[2]])
           set.gamma[[3]]<-c()
         }
         if(length(S)>1){
           S_sub<-(1:p)[-S]
           set.gamma[[1]]<-t(combn(S,length(S)-1))
           if(length(S)>2){
             set.gamma[[1]]<-t(apply(as.matrix(set.gamma[[1]]),1,sort))
           }

           set.gamma[[2]]<-apply(as.matrix(S_sub),1,function(x){return(sort(c(x,S)))})
           set.gamma[[2]]<-t(set.gamma[[2]])

           set.gamma[[3]]<-add.function(set.gamma[[1]][1,])
           for(i in 2:nrow(set.gamma[[1]])){
             set.gamma[[3]]<-rbind(set.gamma[[3]],add.function(set.gamma[[1]][i,]))
           }
         }
         if(length(S)==1){
           S_sub<-(1:p)[-S]
           set.gamma[[1]]<-t(combn(S,length(S)-1))
           set.gamma[[2]]<-apply(as.matrix(S_sub),1,function(x){return(sort(c(x,S)))})
           set.gamma[[2]]<-t(set.gamma[[2]])
           set.gamma[[3]]<-t(add.function(set.gamma[[1]][1,]))
         }
         return(set.gamma)
       }
       set.gamma.func.conditional<-function(input.S,condition.index){
         add.function<-function(y){results<-(apply(as.matrix(S_sub),1,function(x){return(sort(c(x,y)))}))
         return(t(results))
         }

         set.gamma<-list()
         for(i in 1:3){
           set.gamma[[i]]<-c()
         }
         S=input.S[-match(condition.index,input.S)]

         if(length(S)==0){
           S=integer(0)
           S_sub<-(1:p)[-condition.index]
           set.gamma[[1]]<-c()
           set.gamma[[2]]<-apply(as.matrix(S_sub),1,function(x){return(c(x,S))})
           set.gamma[[2]]<-as.matrix(set.gamma[[2]])
           set.gamma[[3]]<-c()
         }
         if(length(S)==1){
           S_sub<-(1:p)[-input.S]
           set.gamma[[1]]<-t(combn(S,length(S)-1))
           set.gamma[[2]]<-apply(as.matrix(S_sub),1,function(x){return(sort(c(x,S)))})
           set.gamma[[2]]<-as.matrix((t(set.gamma[[2]])))
           set.gamma[[3]]<-t(add.function(set.gamma[[1]][1,]))
         }
         if(length(S)>1){
           S_sub<-(1:p)[-input.S]
           if(length(S)>2){
             set.gamma[[1]]<-t(combn(S,length(S)-1))
             set.gamma[[1]]<-t(apply(as.matrix(set.gamma[[1]]),1,sort))
           }else{
             set.gamma[[1]]<-t(combn(S,length(S)-1))
           }
           set.gamma[[2]]<-apply(as.matrix(S_sub),1,function(x){return(sort(c(x,S)))})
           set.gamma[[2]]<-as.matrix(t(set.gamma[[2]]))
           set.gamma[[3]]<-add.function(set.gamma[[1]][1,])
           factorial(3)
           for(i in 2:nrow(set.gamma[[1]])){
             set.gamma[[3]]<-rbind(set.gamma[[3]],add.function(set.gamma[[1]][i,]))
           }
         }

         return(set.gamma)
       }

      if(is.null(condition.index)){
        results<-set.gamma.func.base(input.S)
      }else{
        results<-set.gamma.func.conditional(input.S,condition.index)
      }
      return(results)
    }
    duplicated.dgCMatrix <- function (dgCMat, MARGIN) {
      MARGIN <- as.integer(MARGIN)
      n <- nrow(dgCMat)
      p <- ncol(dgCMat)
      J <- rep(1:p, diff(dgCMat@p))
      I <- dgCMat@i + 1
      x <- dgCMat@x
      if (MARGIN == 1L) {

        names(x) <- J
        RowLst <- split(x, I)
        is_empty <- setdiff(1:n, I)
        result <- duplicated.default(RowLst)
      } else if (MARGIN == 2L) {

        names(x) <- I
        ColLst <- split(x, J)
        is_empty <- setdiff(1:p, J)
        result <- duplicated.default(ColLst)
      } else {
        warning("invalid MARGIN; return NULL")
        result <- NULL
      }

      if(any(is_empty)){
        out <- logical(if(MARGIN == 1L) n else p)
        out[-is_empty] <- result
        if(length(is_empty) > 1)
          out[is_empty[-1]] <- TRUE
        result <- out
      }

      result
    }
    match.dgCMatrix <- function (dgCMat1,dgCMat2) {

      n1 <- nrow(dgCMat1)
      p1 <- ncol(dgCMat1)
      J1 <- rep(1:p1, diff(dgCMat1@p))
      I1 <- dgCMat1@i + 1
      x1 <- dgCMat1@x
      n2 <- nrow(dgCMat2)
      p2 <- ncol(dgCMat2)
      J2 <- rep(1:p2, diff(dgCMat2@p))
      I2 <- dgCMat2@i + 1
      x2 <- dgCMat2@x

      names(x1) <- J1
      RowLst1 <- split(J1, I1)
      is_empty1<- setdiff(1:n1, I1)

      names(x2) <- J2
      RowLst2 <- split(J2, I2)
      is_empty2 <- setdiff(1:n2, I2)
      result<-(match(RowLst2,RowLst1))
      if(any(which(result>is_empty1))){
        result[which(result>=is_empty1)]=result[which(result>=is_empty1)]+1
      }
      if(any(is_empty1)){
        if(any(is_empty2)){
          result<-c(is_empty1,result)
        }
      }else{
        if(any(is_empty2)){
          result<-c(NA,result)
        }
      }
      return(result)

    }

    PIP.func<-function(likeli,model.space){
      infi.index<-which(is.infinite(likeli))
      if(length(infi.index)!=0){
        likeli<-likeli[-infi.index]
        model.space<-model.space[-infi.index,]
      }
      na.index<-which(is.na(likeli))
      if(length(na.index)!=0){
        likeli<-likeli[-na.index]
        model.space<-model.space[-na.index,]
      }
      aa<-likeli-max(likeli,na.rm=T)
      prob.sum<-sum(exp(aa))
      result.prob<-rep(NA,p)
      for(i in 1:p){
        result.prob[i]<-sum(exp(aa[ which(model.space[,i]==1)]))/prob.sum
      }
      return(result.prob)
    }
    index.fun.inner<-function(x){
      m<-as(as(as(matrix(0,nrow=nrow(x),ncol=p), "dMatrix"), "generalMatrix"), "TsparseMatrix")
      m@i<-as.integer(rep(1:nrow(x)-1,each=ncol(x)))
      m@j<-as.integer(c(t(x))-1)
      m@x=rep(1,nrow(x)*ncol(x))
      m<-as(m,"CsparseMatrix")
      return(m)
    }
    index.fun<-function(outer.x,Max.Model.Dimins=10){
     if(nrow(outer.x)>1000){
    index.bins<-which((1:nrow(outer.x))%%floor(nrow(outer.x)/Max.Model.Dimins)==0)
    result.m<-index.fun.inner(outer.x[1:index.bins[1],,drop=F])
    for(b in 1:(length(index.bins)-1)){
      result.m<-rbind(result.m,index.fun.inner(outer.x[(index.bins[b]+1):index.bins[b+1],,drop=F]))
    }
    if(index.bins[length(index.bins)]!=nrow(outer.x)){
      result.m<-rbind(result.m,index.fun.inner(outer.x[(index.bins[length(index.bins)]+1):nrow(outer.x),,drop=F]))
    }
     }else{
         result.m<-index.fun.inner(outer.x)
     }
    return(result.m)
    }
    ridge.fun<-function(x){
      temp.ld.S<-x*modi.ld.S+(1-x)*diag(nrow(modi.ld.S))
      temp.Sigma[test.S,test.S]<-temp.ld.S
      return( outlier_likelihood(test.S,temp.Sigma,z,outlier.tau,length(test.S),1) )
    }

    for(l in 1:inner.all.iter){
      for(h in 1:10){

        {
          set.gamma<-set.gamma.func(S,conditional.S)
          if(is.null(conditional.S)){
            working.S=S
            base.model<-null.model
            base.model.margin<-null.margin
          }else{
            working.S=S[-match(conditional.S,S)]
            if(length(working.S)!=0){
            base.model<-as(Matrix(nrow=1,ncol=p,sparse = T,data=0),'CsparseMatrix')
            base.model[,working.S]<-1
            p_S=length(working.S);
            base.model.margin<-marginal_likelihood(working.S,Sigma,z,tau=tau.sample,p_S=p_S,y.var)+prior.dist(base.model)
            }else{
              base.model<-null.model
              base.model.margin<-null.margin
            }
            }
          set.gamma.margin<-list()
          set.gamma.prior<-list()
          matrix.gamma<-list()
          if(length(working.S)!=0){
            S.model<-as(Matrix(nrow=1,ncol=p,sparse = T,data=0),'CsparseMatrix')
            S.model[,working.S]<-1
            p_S=length(working.S);
            current.log.margin<-marginal_likelihood(working.S,Sigma,z,tau=tau.sample,p_S=p_S,y.var)+prior.dist(S.model)
          }else{
            current.log.margin<-prior.dist(null.model)
          }

          if(length(working.S)>1){
            for(i in  1:length(set.gamma)){
              t0=Sys.time()
              matrix.gamma[[i]]<-index.fun(set.gamma[[i]])

              if(length(C.list[[2]])<ncol(set.gamma[[i]])){
                C.list[[2]][[ncol(set.gamma[[i]])]]<-as(Matrix(nrow=0,ncol=p,sparse = T,data=0),'CsparseMatrix')
                C.list[[1]][[ncol(set.gamma[[i]])]]<-integer(0)
                computed.index<-integer(0)
              }else{
                computed.index<-match.dgCMatrix(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]])
              }

              p_S=dim(set.gamma[[i]])[2]
              if(length(na.omit(computed.index))==0){
                set.gamma.margin[[i]]<-apply(set.gamma[[i]],1,marginal_likelihood,Sigma=Sigma,z=z,tau=tau.sample,p_S=p_S,y_sigma=y.var)
                C.list[[1]][[ncol(set.gamma[[i]])]]<-c(C.list[[1]][[ncol(set.gamma[[i]])]],set.gamma.margin[[i]])
                C.list[[2]][[ncol(set.gamma[[i]])]]<-rbind(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]])
                set.gamma.prior[[i]]<-apply(matrix.gamma[[i]],1,prior.dist)
                set.gamma.margin[[i]]=set.gamma.prior[[i]]+set.gamma.margin[[i]]
              }else{
                set.gamma.margin[[i]]<-rep(NA,nrow(matrix.gamma[[i]]))
                set.gamma.margin[[i]][!is.na(computed.index)]<-C.list[[1]][[ncol(set.gamma[[i]])]][na.omit(computed.index)]
                if(sum(is.na(computed.index))!=0){
                  set.gamma.margin[[i]][is.na(computed.index)]<-apply(set.gamma[[i]][is.na(computed.index),,drop=F],1,marginal_likelihood,Sigma=Sigma,z=z,tau=tau.sample,p_S=p_S,y_sigma=y.var)
                }
                C.list[[1]][[ncol(set.gamma[[i]])]]<-c(C.list[[1]][[ncol(set.gamma[[i]])]],set.gamma.margin[[i]][is.na(computed.index)])
                C.list[[2]][[ncol(set.gamma[[i]])]]<-rbind(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]][is.na(computed.index),,drop=F])
                set.gamma.margin[[i]]<- set.gamma.margin[[i]]+apply(matrix.gamma[[i]],1,prior.dist)
              }
              t1=Sys.time()-t0
            }

            add.B<-list()
            add.B[[1]]<-c(set.gamma.margin[[1]],
                          set.gamma.margin[[2]],
                          set.gamma.margin[[3]])
            add.B[[2]]<-as(Matrix(nrow=0,ncol=p,sparse = T,data=0),'CsparseMatrix')
            for(i in 1:3){
              add.B[[2]]<-rbind(add.B[[2]],matrix.gamma[[i]])
            }

          }

          if(length(working.S)==1){
            set.gamma.margin[[1]]<-null.margin
            matrix.gamma[[1]]<-null.model
            for(i in 2:3){

              matrix.gamma[[i]]<-index.fun(set.gamma[[i]])

              if(length(C.list[[2]])<ncol(set.gamma[[i]])){
                C.list[[2]][[ncol(set.gamma[[i]])]]<-as(Matrix(nrow=0,ncol=p,sparse = T,data=0),'CsparseMatrix')
                C.list[[1]][[ncol(set.gamma[[i]])]]<-integer(0)
                computed.index<-integer(0)
              }else{
                computed.index<-match.dgCMatrix(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]])
              }

              p_S=dim(set.gamma[[i]])[2]
              if(length(na.omit(computed.index))==0){
                set.gamma.margin[[i]]<-apply(set.gamma[[i]],1,marginal_likelihood,Sigma=Sigma,z=z,tau=tau.sample,p_S=p_S,y_sigma=y.var)
                C.list[[1]][[ncol(set.gamma[[i]])]]<-c(C.list[[1]][[ncol(set.gamma[[i]])]],set.gamma.margin[[i]])
                C.list[[2]][[ncol(set.gamma[[i]])]]<-rbind(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]])
                set.gamma.prior[[i]]<-apply(matrix.gamma[[i]],1,prior.dist)
                set.gamma.margin[[i]]=set.gamma.prior[[i]]+set.gamma.margin[[i]]
              }else{
                set.gamma.margin[[i]]<-rep(NA,nrow(matrix.gamma[[i]]))
                set.gamma.margin[[i]][!is.na(computed.index)]<-C.list[[1]][[ncol(set.gamma[[i]])]][na.omit(computed.index)]
                if(sum(is.na(computed.index))!=0){
                  set.gamma.margin[[i]][is.na(computed.index)]<-apply(set.gamma[[i]][is.na(computed.index),,drop=F],1,marginal_likelihood,Sigma=Sigma,z=z,tau=tau.sample,p_S=p_S,y_sigma=y.var)
                }
                C.list[[1]][[ncol(set.gamma[[i]])]]<-c(C.list[[1]][[ncol(set.gamma[[i]])]],set.gamma.margin[[i]][is.na(computed.index)])
                C.list[[2]][[ncol(set.gamma[[i]])]]<-rbind(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]][is.na(computed.index),,drop=F])
                set.gamma.margin[[i]]<- set.gamma.margin[[i]]+apply(matrix.gamma[[i]],1,prior.dist)
              }
            }

            add.B<-list()
            add.B[[1]]<-c(set.gamma.margin[[1]],
                          set.gamma.margin[[2]],
                          set.gamma.margin[[3]])
            add.B[[2]]<-as(Matrix(nrow=0,ncol=p,sparse = T,data=0),'CsparseMatrix')
            for(i in 1:3){
              add.B[[2]]<-rbind(add.B[[2]],matrix.gamma[[i]])
            }
           }
          if(length(working.S)==0){

            for(i in 2){
              matrix.gamma[[i]]<-index.fun(set.gamma[[i]])

              if(length(C.list[[2]])<ncol(set.gamma[[i]])){
                C.list[[2]][[ncol(set.gamma[[i]])]]<-as(Matrix(nrow=0,ncol=p,sparse = T,data=0),'CsparseMatrix')
                C.list[[1]][[ncol(set.gamma[[i]])]]<-integer(0)
                computed.index<-integer(0)
              }else{
                computed.index<-match.dgCMatrix(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]])
              }

              p_S=dim(set.gamma[[i]])[2]
              if(length(na.omit(computed.index))==0){
                set.gamma.margin[[i]]<-apply(set.gamma[[i]],1,marginal_likelihood,Sigma=Sigma,z=z,tau=tau.sample,p_S=p_S,y_sigma=y.var)
                C.list[[1]][[ncol(set.gamma[[i]])]]<-c(C.list[[1]][[ncol(set.gamma[[i]])]],set.gamma.margin[[i]])
                C.list[[2]][[ncol(set.gamma[[i]])]]<-rbind(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]])
                set.gamma.prior[[i]]<-apply(matrix.gamma[[i]],1,prior.dist)
                set.gamma.margin[[i]]=set.gamma.prior[[i]]+set.gamma.margin[[i]]
              }else{
                set.gamma.margin[[i]]<-rep(NA,nrow(matrix.gamma[[i]]))
                set.gamma.margin[[i]][!is.na(computed.index)]<-C.list[[1]][[ncol(set.gamma[[i]])]][na.omit(computed.index)]
                if(sum(is.na(computed.index))!=0){
                  set.gamma.margin[[i]][is.na(computed.index)]<-apply(set.gamma[[i]][is.na(computed.index),,drop=F],1,marginal_likelihood,Sigma=Sigma,z=z,tau=tau.sample,p_S=p_S,y_sigma=y.var)
                }
                C.list[[1]][[ncol(set.gamma[[i]])]]<-c(C.list[[1]][[ncol(set.gamma[[i]])]],set.gamma.margin[[i]][is.na(computed.index)])
                C.list[[2]][[ncol(set.gamma[[i]])]]<-rbind(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]][is.na(computed.index),,drop=F])
                set.gamma.margin[[i]]<- set.gamma.margin[[i]]+ apply(matrix.gamma[[i]],1,prior.dist)
              }
            }

            add.B<-list()
            add.B[[1]]<-c(set.gamma.margin[[2]])
            add.B[[2]]<-matrix.gamma[[2]]
          }

          keep.model <- Matrix::rowSums(add.B[[2]]) <= num.causal
          add.B[[1]] <- add.B[[1]][keep.model]
          add.B[[2]] <- add.B[[2]][keep.model, , drop = FALSE]
          stopifnot(nrow(add.B[[2]]) > 0L)
          add.index<-match.dgCMatrix(B.list[[2]],add.B[[2]])
          if(length(which(!is.na(add.index)))>10){
            check.index<-sample(which(!is.na(add.index)),10)

          }
          if(length(na.omit(add.index))!=0){
            B.list[[1]]<-c((B.list[[1]]),(add.B[[1]][is.na(add.index)]))
            B.list[[2]]<-rbind(B.list[[2]],add.B[[2]][is.na(add.index),,drop=F])
          }else{
            B.list[[1]]<-c((B.list[[1]]),(add.B[[1]]))
            B.list[[2]]<-rbind(B.list[[2]],add.B[[2]])
          }
          B.list[[2]]<-B.list[[2]][order(B.list[[1]],decreasing = T),]
          B.list[[1]]<-B.list[[1]][order(B.list[[1]],decreasing = T)]

          if(length(working.S)!=0){
            set.star<-data.frame(set.index=1:3,gamma.set.index=rep(NA,3),margin=rep(NA,3))
            for(i in 1){
              aa<-set.gamma.margin[[i]]-current.log.margin
              aa<-aa-aa[which.max(aa)]
              if(length(which(is.nan(aa)))!=0){
              aa[which(is.nan(aa))]<-min(aa)
              }
              set.star$gamma.set.index[i] <-c(sample(1:length(set.gamma.margin[[i]]),1,prob=exp(aa)))
              set.star$margin[i]<-set.gamma.margin[[i]][  set.star$gamma.set.index[i]]
              rm(aa)
            }

            if(outlier.switch){
            for(i in 2:length(set.gamma)){
              repeat{

              aa<-set.gamma.margin[[i]]-current.log.margin
              aa<-aa-aa[which.max(aa)]
              if(length(which(is.nan(aa)))!=0){
                aa[which(is.nan(aa))]<-min(aa[!is.na(aa)])
              }

              set.star$gamma.set.index[i]<-c(sample((1:length(set.gamma.margin[[i]])),
                                                    1,prob=exp(aa)))
              set.star$margin[i]<-set.gamma.margin[[i]][  set.star$gamma.set.index[i]]

              test.S<-set.gamma[[i]][set.star$gamma.set.index[i],]

              modi.Sigma<-Sigma
              temp.Sigma<-Sigma
              if(length(test.S)>1){

              modi.ld.S<- modi.Sigma[test.S,test.S]

                opizer<-optimize(ridge.fun,interval=c(0,1),maximum = T)
                modi.ld.S<-opizer$maximum*modi.ld.S+(1-opizer$maximum)*diag(nrow(modi.ld.S))

              modi.Sigma[test.S,test.S]<-modi.ld.S

              test.log.BF<-outlier_likelihood(test.S,Sigma,z,outlier.tau,length(test.S),1)-outlier_likelihood(test.S,modi.Sigma,z,outlier.tau,length(test.S),1)
              test.log.BF<--abs(test.log.BF)
              if(printing.log==T){
              print(paste0('Outlier BF: ', test.log.BF))
              print(test.S)
              print(paste0('This is xi hat: ', opizer))
              }
              }

              if(exp(test.log.BF)<outlier.BF.index){
                set.gamma[[i]]<-set.gamma[[i]][-set.star$gamma.set.index[i],]
                set.gamma.margin[[i]]<-set.gamma.margin[[i]][-set.star$gamma.set.index[i]]
                conditional.S<-c(conditional.S,test.S[is.na(match(test.S,working.S))])
                conditional.S<-unique(conditional.S)
              }else{
                break
              }
              }
              rm(aa)
            }
            }else{
              for(i in 2:length(set.gamma)){

                  aa<-set.gamma.margin[[i]]-current.log.margin
                  aa<-aa-aa[which.max(aa)]
                  if(length(which(is.nan(aa)))!=0){
                    aa[which(is.nan(aa))]<-min(aa[!is.na(aa)])
                  }

                  set.star$gamma.set.index[i]<-c(sample((1:length(set.gamma.margin[[i]])),
                                                        1,prob=exp(aa)))
                  set.star$margin[i]<-set.gamma.margin[[i]][  set.star$gamma.set.index[i]]
                rm(aa)
              }
            }
            if(printing.log==T){
           print(set.star)
          }
            if(length(working.S)>=num.causal){
              set.star<-set.star[-2,]
              aa<-set.star$margin-current.log.margin-max(set.star$margin-current.log.margin)
              sec.sample<-sample(c(1,3),1,prob=exp(aa))
              S<-set.gamma[[sec.sample]][set.star$gamma.set.index[[which(sec.sample==set.star$set.index)]] ,]
            }else{
              aa<-set.star$margin-current.log.margin-max(set.star$margin-current.log.margin)
            sec.sample<-sample(1:3,1,prob=exp(aa) )
            S<-set.gamma[[sec.sample]][set.star$gamma.set.index[[sec.sample]] ,]
            }

            }else{
                set.star<-data.frame(set.index=rep(1,3),gamma.set.index=rep(NA,3),margin=rep(NA,3))
                         aa<-set.gamma.margin[[2]]-current.log.margin
                         aa<-aa-aa[which.max(aa)]
                         if(length(which(is.nan(aa)))!=0){
                           aa[which(is.nan(aa))]<-min(aa)
                         }

                         set.star$gamma.set.index[2] <-c(sample((1:length(set.gamma.margin[[2]]))[order(exp(aa),decreasing = T)[1:(min(length(aa),floor(p/2)))]],
                                                             1,prob=exp(aa)[order(exp(aa),decreasing = T)[1:(min(length(aa),floor(p/2)))]]))
                         set.star$margin[2]<-set.gamma.margin[[2]][  set.star$gamma.set.index[2]]

                         S<-set.gamma[[2]][set.star$gamma.set.index[2],]
                         if(printing.log==T){
                           print(set.star)
                         }
            }
          if(printing.log==T){
          print(paste0('this is running S: ',paste0(S,collapse = ',')))
          }
          S<-unique(c(S,conditional.S))
        }

        }

      result.B.list<-list()
      if(!is.null(conditional.S)){
      all.c.index<-c()

        for(tt in conditional.S){
          c.index<-(B.list[[2]]@i[min(length(B.list[[2]]@i),(B.list[[2]]@p[tt]+1)):B.list[[2]]@p[tt+1]])+1
          all.c.index<-c(all.c.index,c.index)
          }

      all.c.index<-unique(all.c.index)
      temp.B.list<-list()
      temp.B.list[[1]]<-B.list[[1]][-all.c.index]
      temp.B.list[[2]]<-B.list[[2]][-all.c.index,]
      }else{
        temp.B.list<-list()
        temp.B.list[[1]]<-B.list[[1]]
        temp.B.list[[2]]<-B.list[[2]]
      }
      result.B.list<-list()
      result.B.list[[1]]<-temp.B.list[[1]][(1:min(B,nrow(temp.B.list[[2]])))]
      result.B.list[[2]]<-temp.B.list[[2]][(1:min(B,nrow(temp.B.list[[2]]))),]

      stopifnot(all(Matrix::rowSums(result.B.list[[2]]) <= num.causal))
      if(num.causal==1){
        single.set<-matrix(1:p,p,1)
        single.marginal<-apply(single.set,1,marginal_likelihood,Sigma=Sigma,z=z,tau=tau.sample,p_S=p_S,y_sigma=y.var)
        aa<-single.marginal-max(single.marginal,na.rm=T)
        prob.sum<-sum(exp(aa))
        result.prob<-(exp(aa))/prob.sum
      }else{
        result.prob<-PIP.func(result.B.list[[1]],result.B.list[[2]])
      }
      conditional.S.list<-data.frame(Index=conditional.S,Z=z[conditional.S,])
      if(!is.null(output.labels)){
        if(dir.exists(output.labels)==F ){
          dir.create(output.labels,recursive = T)
        }
      write.table(result.B.list[[1]],file=paste0(output.labels,'/post_',label,'_poi_likeli','.txt'),row.names = F,col.names = F)
      writeMM(result.B.list[[2]],file=paste0(output.labels,'/post_',label,'_poi_gamma','.mtx'))
      write.table((result.prob),file=paste0(output.labels,'/post_', label,'.txt'),row.names = F,append = F,col.names = F)
      if(outlier.switch){
        saveRDS(conditional.S.list,file=paste0(output.labels,'/post_', label,'_','outliers.rds'))
      }
      }

      difference<-abs(mean(result.B.list[[1]][1:round(quantile(1:length(result.B.list[[1]]),probs = 0.25))])-stored.bf)

       if(difference<epsilon){
         break
       }else{
         stored.bf<-mean(result.B.list[[1]][1:round(quantile(1:length(result.B.list[[1]]),probs = 0.25))])
       }
    }
    return(list(result.B.list,C.list,result.prob,conditional.S.list,prob.list))
  }

  previous.result<-list()

  for(i in 1:L){
    t0=Sys.time()
   all.C.list[[i]]<-Module.Cauchy.Shotgun(z.list[[i]],ld.list[[i]],epsilon=epsilon.list[[i]],
                                           Max.Model.Dim=Max.Model.Dim,lambda = lambda.list[[i]],
                                          outlier.switch=outlier.switch,tau=tau,
                                           num.causal = num.causal,y.var=y.var,printing.log=printing.log,
                                           label = label.list[[i]],output.labels = output.labels,
                                           effect.size.prior=effect.size.prior,model.prior=model.prior,inner.all.iter = all.inner.iter)
   t1=Sys.time()-t0
   print(paste0('This is locus ',i,' burning time'))
   print((t1))

  }

  for(g in 1:all.iter){
    if(outlier.switch){
      delete.list<-list()
      for(i in 1:L){
        delete.list[[i]]<-integer(0)
        if(nrow(all.C.list[[i]][[4]])!=0){
        temp.delete.list<-c(all.C.list[[i]][[4]]$Index)
        delete.list[[i]]<-temp.delete.list
        }
      }
    }else{
      delete.list<-list()
      for(i in 1:L){
      delete.list[[i]]<-integer(0)
      }
      }

    if(!is.null(w.list)){
      w<-matrix(NA,nrow=0,ncol=ncol(w.list[[1]]))
      colnames(w)<-colnames(w.list[[1]])
      for(i in 1:L){
        if(length(delete.list[[i]])!=0){
          w<-rbind(w,w.list[[i]][-delete.list[[i]],])
        }else{
          w<-rbind(w,w.list[[i]])
        }
      }
    }

    for(i in 1:L){
      previous.result[[i]]<-mean(all.C.list[[i]][[1]][[1]][1:round(quantile(1:length(all.C.list[[i]][[1]][[1]]),probs = 0.25))])
    }
    if(!is.null(w.list)){

    if(EM.dist=='Poisson'){
    if(!standardize.model.space){
      model.space.count<-c()
      for(i in 1:L){
        if(length(delete.list[[i]])!=0){
        model.space.count<-c(model.space.count,colSums(all.C.list[[i]][[1]][[2]][,-delete.list[[i]]]))
        }else{
          model.space.count<-c(model.space.count,colSums(all.C.list[[i]][[1]][[2]]))
        }

      }
    }else{
      model.space.count<-c()
      for(i in 1:L){
        if(length(delete.list[[i]])!=0){
        indi.count<-colSums(all.C.list[[i]][[1]][[2]][,-delete.list[[i]]])
        indi.count<-floor(colSums(all.C.list[[i]][[1]][[2]][,-delete.list[[i]]])/nrow(all.C.list[[i]][[1]][[2]][,-delete.list[[i]]])*Max.Model.Dim)
        }else{
          indi.count<-colSums(all.C.list[[i]][[1]][[2]])
          indi.count<-floor(colSums(all.C.list[[i]][[1]][[2]])/nrow(all.C.list[[i]][[1]][[2]])*Max.Model.Dim)
        }
        model.space.count<-c(model.space.count,indi.count)
      }
    }
    M.step.response=model.space.count

    }

    if(EM.dist=='Logistic'){
      M.step.response<-c()
      for(i in 1:L){
        if(length(delete.list[[i]])!=0){
          indi.pip<-(all.C.list[[i]][[3]][-delete.list[[i]]])
         }else{
          indi.pip<-(all.C.list[[i]][[3]])
        }
        M.step.response<-c(M.step.response,indi.pip)
      }
    }

    try.index<-try(glm.beta<-EM.M.step.func(input.response  = M.step.response ,w=w,input.alpha=input.alpha,EM.dist=EM.dist))
    prior.prob.list<-list()

    if(class(try.index)[1]!='try-error'){
      for(i in 1:L){
        if(prior.prob.computation=='Intercept.approx'){
          glm.beta[1]=log((min(Max.Model.Dim,nrow(all.C.list[[i]][[1]][[2]])))*lambda.list[[i]]/(lambda.list[[i]]+p.list[[i]]))
          prior.prob.list[[i]]<-(exp(w.list[[i]]%*%glm.beta)/(max(1+max(exp(w.list[[i]]%*%glm.beta)),min(Max.Model.Dim,nrow(all.C.list[[i]][[1]][[2]])))))

        }

        if(prior.prob.computation=='Logistic'){
          prior.prob.list[[i]]<-plogis(w.list[[i]]%*%glm.beta)

        }
        if(!is.null(output.labels)){
        write.table((glm.beta),file=paste0(output.labels,'/post_', label.list[[i]],'_theta.txt'),row.names = F,append = F,col.names = F)
        }
      }
     model.prior='input.prob'

    }else{
      model.prior='Poisson'
      prior.prob.list<-list()
      for(i in 1:L){
        prior.prob.list[[i]]<-list(NULL)
      }

    }

    }else{
      prior.prob.list<-list()
      for(i in 1:L){
        prior.prob.list[[i]]<-list(NULL)
      }
    }

    for(i in 1:L){
      t0=Sys.time()
      all.C.list[[i]]<-Module.Cauchy.Shotgun(z=z.list[[i]],ld.list[[i]],input.conditional.S.list = all.C.list[[i]][[4]],
                                             Max.Model.Dim=Max.Model.Dim,y.var=y.var,num.causal = num.causal,epsilon=epsilon.list[[i]],
                                             C.list = all.C.list[[i]][[2]],prior.prob = prior.prob.list[[i]],
                                             outlier.switch=outlier.switch,tau=tau,printing.log=printing.log,
                                             lambda = lambda.list[[i]], label = label.list[[i]],output.labels = output.labels,
                                             effect.size.prior=effect.size.prior,model.prior=model.prior,inner.all.iter = all.inner.iter)
      t1=Sys.time()-t0
      print(paste0('This is locus ',i,' computing time'))
      print((t1))
    }

    difference<-0
     for(i in 1:L){
       difference<-difference+abs(previous.result[[i]]-mean(all.C.list[[i]][[1]][[1]][1:round(quantile(1:length(all.C.list[[i]][[1]][[1]]),probs = 0.25))]))
     }
    if(printing.log==T){
      print(paste0('This is difference; ',difference))
    }

     if(difference<all.epsilon.threshold){
       break
     }
  }

  results.list<-list()
  for(i in 1:L){
    results.list[[i]]<-list()
    pip=all.C.list[[i]][[3]]
    credible.set<-credible.set.fun.improved(pip,ld.list[[i]],rho=rho.index)
    credible.model<-credible.model.fun(all.C.list[[i]][[1]][[1]],all.C.list[[i]][[1]][[2]],bayes.threshold = BF.index)
    results.list[[i]][[1]]<-pip
    results.list[[i]][[2]]<-credible.set
    results.list[[i]][[3]]<-credible.model
    results.list[[i]][[4]]<-all.C.list[[i]][[4]]
    names(results.list[[i]])<-c('PIPs','Credible set','Credible model','Outliers')

    }
  return(results.list)
}

environment(CARMA_strict) <- asNamespace("CARMA")
# Required study helper: finemap_helpers.R
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
