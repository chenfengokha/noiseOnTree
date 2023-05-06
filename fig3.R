library(plyr)
library(dplyr)
library(Seurat)
library(tidyverse)
library(parallel)
library(segmented)
#library(ggtree)
library(ape)
library(lmtest)
options(dplyr.summarise.inform = FALSE)
## ========= define and source functions =============
check.mean.dist.nor.perGene <- function(raw.dist.df, ran.dist){
  mean.ran.dist <- mean(ran.dist)
  median.ran.dist <- median(ran.dist)
  re.depth.mean.expdist <- 
    raw.dist.df %>% dplyr::group_by(mrca, treeDist.asDepth) %>% 
    dplyr::summarise(mrca.mean.dist = mean(exp.dist),
                     mrca.median.dist = median(exp.dist)) %>%
    dplyr::group_by(treeDist.asDepth) %>% 
    dplyr::summarise(mean.median.dist = median(mrca.mean.dist)/mean.ran.dist)
  return(re.depth.mean.expdist)
  
}



test.fit.gene <- function(raw.depth.exp.dist, depth.cut = NULL, test.psi = -3){
  # test fit segment
  test.df <- raw.depth.exp.dist
  ran.dist <- sample(test.df$exp.dist %>% unlist(), 
                     size = nrow(test.df), replace = T)
  if (!is.null(depth.cut)){
    test.df$treeDist.asDepth[which(test.df$treeDist.asDepth > depth.cut)] <- depth.cut
  }
  #test.sampleB <- test.sampleB %>% filter(., mrca != 1155)
  test.depth.mean.expdist <- 
    check.mean.dist.nor.perGene(raw.dist.df = test.df,
                                ran.dist = ran.dist)
  ## fit segmented
  x <- test.depth.mean.expdist$treeDist.asDepth
  y <- test.depth.mean.expdist$mean.median.dist
  ## first fit a no change model
  o <- lm(y ~ 1)
  seg.fit <- 
    tryCatch({
      xx <- -x
      o.seg <- segmented(o, seg.Z = ~xx, psi=list(xx=test.psi))
      if (!is.null(o.seg)) {
        o.seg
      } else {
        NULL
      }
    }, error = function(e){
      NULL
    })
  if (is.null(seg.fit) | !exists("seg.fit")) seg.fit <- NULL
  ## fit linear 
  liner.lm <- lm(formula = mean.median.dist ~ treeDist.asDepth, data = test.depth.mean.expdist)
  ## compare three models
  # calculate all AIC values
  # 1. linear model
  aic.lm <- AIC(liner.lm,k=log(length(x)))
  loglik.linear <- logLik(liner.lm) %>% as.data.frame() %>% unlist()
  # 2. no change model
  aic.lm0 <- AIC(o, k=log(length(x)))
  loglik.0 <- logLik(o) %>% as.data.frame() %>% unlist()
  # 3. segment
  if (is.null(seg.fit)){
    aic.seg <- 9999
    loglik.segfit <- 0
    ## compare seg and linear
    segVSliner.lr.chisq.value = 0
    segVSliner.lr.chi.P.value = 1
    segVSzero.lr.chisq.value = 0
    segVSzero.lr.chi.P.value = 1
    ##lm(y~1) and lm(y~x)
    zeroVSliner.lr <- lrtest(o, liner.lm)
    zeroVSliner.lr.chisq.value <- zeroVSliner.lr$Chisq[!is.na(zeroVSliner.lr$Chisq)]
    zeroVSliner.lr.chi.P.value <- zeroVSliner.lr$`Pr(>Chisq)`[!is.na(zeroVSliner.lr$`Pr(>Chisq)`)]
    if (length(zeroVSliner.lr.chisq.value) == 0) {
      zeroVSliner.lr.chisq.value <- 0
      zeroVSliner.lr.chi.P.value <- 1
    }
    #### == extract slope of gene fit ====
    re.df <- 
      data.frame(aic.seg = aic.seg,
                 aic.lm = aic.lm,
                 aic.lm0 = aic.lm0,
                 loglik.linear = loglik.linear,
                 loglik.0 = loglik.0,
                 loglik.segfit = loglik.segfit,
                 segVSliner.lr.chisq.value = segVSliner.lr.chisq.value,
                 segVSliner.lr.chi.P.value = segVSliner.lr.chi.P.value,
                 segVSzero.lr.chisq.value = segVSzero.lr.chisq.value,
                 segVSzero.lr.chi.P.value = segVSzero.lr.chi.P.value,
                 zeroVSliner.lr.chisq.value = zeroVSliner.lr.chisq.value,
                 zeroVSliner.lr.chi.P.value = zeroVSliner.lr.chi.P.value,
                 slope=999,
                 slopse=999,
                 slopP=999,
                 inter=999,
                 interse=999,
                 interP=999,
                 mybreak=999,
                 mybreakse=999)
  } else {
    aic.seg <- AIC(seg.fit,k=log(length(x)))
    loglik.segfit <- logLik(seg.fit)
    ## compare seg and linear
    segVSliner.lr <- lrtest(seg.fit, liner.lm)
    segVSliner.lr.chisq.value <- segVSliner.lr$Chisq[!is.na(segVSliner.lr$Chisq)]
    segVSliner.lr.chi.P.value <- segVSliner.lr$`Pr(>Chisq)`[!is.na(segVSliner.lr$`Pr(>Chisq)`)]
    ## compare seg and nochange
    segVSzero.lr <- lrtest(seg.fit, o)
    segVSzero.lr.chisq.value <- segVSzero.lr$Chisq[!is.na(segVSzero.lr$Chisq)]
    segVSzero.lr.chi.P.value <- segVSzero.lr$`Pr(>Chisq)`[!is.na(segVSzero.lr$`Pr(>Chisq)`)]
    ##lm(y~1) and lm(y~x)
    zeroVSliner.lr <- lrtest(o, liner.lm)
    zeroVSliner.lr.chisq.value <- zeroVSliner.lr$Chisq[!is.na(zeroVSliner.lr$Chisq)]
    zeroVSliner.lr.chi.P.value <- zeroVSliner.lr$`Pr(>Chisq)`[!is.na(zeroVSliner.lr$`Pr(>Chisq)`)]
    #### == extract slope of gene fit ====
    seg.sumy <- summary(seg.fit)
    coef.df <- seg.sumy$coefficients %>% as.data.frame()
    psi.df <- seg.sumy$psi %>% as.data.frame()
    if(nrow(psi.df)==1){
      re.df <- 
        data.frame(aic.seg = aic.seg,
                   aic.lm = aic.lm,
                   aic.lm0 = aic.lm0,
                   loglik.linear = loglik.linear,
                   loglik.0 = loglik.0,
                   loglik.segfit = loglik.segfit,
                   segVSliner.lr.chisq.value = segVSliner.lr.chisq.value,
                   segVSliner.lr.chi.P.value = segVSliner.lr.chi.P.value,
                   segVSzero.lr.chisq.value = segVSzero.lr.chisq.value,
                   segVSzero.lr.chi.P.value = segVSzero.lr.chi.P.value,
                   zeroVSliner.lr.chisq.value = zeroVSliner.lr.chisq.value,
                   zeroVSliner.lr.chi.P.value = zeroVSliner.lr.chi.P.value,
                   slope=-coef.df[2,1],
                   slopse=coef.df[2,2],
                   slopP=pt(coef.df[2,3],6),
                   inter=coef.df[1,1],
                   interse=coef.df[1,2],
                   interP=coef.df[1,4],
                   mybreak=-psi.df[1,2],
                   mybreakse=psi.df)
    } else{
      re.df <- 
        data.frame(aic.seg = aic.seg,
                   aic.lm = aic.lm,
                   aic.lm0 = aic.lm0,
                   loglik.linear = loglik.linear,
                   loglik.0 = loglik.0,
                   loglik.segfit = loglik.segfit,
                   segVSliner.lr.chisq.value = segVSliner.lr.chisq.value,
                   segVSliner.lr.chi.P.value = segVSliner.lr.chi.P.value,
                   segVSzero.lr.chisq.value = segVSzero.lr.chisq.value,
                   segVSzero.lr.chi.P.value = segVSzero.lr.chi.P.value,
                   zeroVSliner.lr.chisq.value = zeroVSliner.lr.chisq.value,
                   zeroVSliner.lr.chi.P.value = zeroVSliner.lr.chi.P.value,
                   slope=999,
                   slopse=999,
                   slopP=999,
                   inter=999,
                   interse=999,
                   interP=999,
                   mybreak=999,
                   mybreakse=999)
    }
    
  }
  return(re.df)
}


## ============= save and load vars =================
#pair.dist.path1 <- "/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/new_samples_include10x_adjust_results/rawBC/CF2/CF2_perGene_pair_dist.Rds"
pair.dist.path1 <- "/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/sampleABC_re_saturation_results/sampleABC_exp_saturation/sampleA/sampleA_perGene_pair_dist.Rds" 


pergene.dist.df <- readRDS(pair.dist.path1)
#pergene.dist.df <- pergene.dist.df %>% dplyr::select(-exp.dist)

# tmp <- data.frame(stringsAsFactors = F,psi=rep(-(3:6),rep(5,4)),depth.cut=rep(c(7:10,15),4))
# y=7
tmp1 <- data.frame(stringsAsFactors = F,psi=-4,depth.cut=8)


tmp2 <- mclapply(mc.cores = 10,1:20,function(i){
  pergene.fit.df <- 
    mclapply(seq(6, ncol(pergene.dist.df)), function(g) {
      sub.dist.df <- pergene.dist.df[,c(1:5, g)]
      gene.name <- colnames(sub.dist.df)[ncol(sub.dist.df)]
      #message(gene.name)
      # 1. change last colnames as exp.dist
      colnames(sub.dist.df)[ncol(sub.dist.df)] <- "exp.dist"
      gene.fit.re <- 
        test.fit.gene(raw.depth.exp.dist = sub.dist.df,
                      depth.cut = tmp1$depth.cut, 
                      test.psi = tmp1$psi)
      gene.fit.re$gene <- gene.name
      gene.fit.re$i <- i
      return(gene.fit.re)
    }, mc.cores = 6) %>% bind_rows(.)
  
}) %>% rbind.fill() %>% cbind(tmp1)

#saveRDS(tmp2,file="~/project/293Tcelllineagetree/rev/res.C.saturatedGene.all.Rds")
saveRDS(tmp2,file="~/project/293Tcelllineagetree/rev/res.A.saturatedGene.all.Rds")