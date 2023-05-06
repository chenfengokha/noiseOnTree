library(phytools)
library(plyr)
library(dplyr)
library(Seurat)
library(tidyverse)
library(parallel)
library(segmented)
library(ggtree)
library(ape)
library(lmtest)
options(dplyr.summarise.inform = FALSE)
source("~/Rfunction/style.print.R")
###
distance <- (readRDS("~/project/293Tcelllineagetree/plot/sampleA_raw_exp_depth_euclidean_dist_nocc.Rds")$real.dist)[,c(1:3,5)]
mytree <- read.tree("~/project/293Tcelllineagetree/plot/test/adjust_filterMt10.nwk")
#OU_bw_tmp <- read.csv("~/project/293Tcelllineagetree/rev/100.saturatedgene.csv",header = T)
aatmp <- readRDS("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/new_CF2_sampleA_OU/sampleA_all_gene_OUvsBW.Rds")

tmpalthe <- lapply(1:15,function(x){
  data.frame(stringsAsFactors = F,alpha=mean(aatmp$alpha),theta=0.005*1.64^(x-1))
}) %>% rbind.fill()

myall <- 1000
simulateOU_BW_Sat <- mclapply(mc.cores = 15,1:nrow(tmpalthe),function(i){
  OU_BW1 <- tmpalthe[i,]
  simuexp<- mclapply(mc.cores = 5,1:myall,function(n){
    #expression
    #sim OU
    a <- fastBM(mytree,a=OU_BW1$theta,mu=2,sig2=aatmp$sigma %>% mean(),model="OU",theta=OU_BW1$theta,alpha=OU_BW1$alpha) %>% as.data.frame()
    names(a)[1] <-"exp"
    a$type <- "OU"
    a$node <- rownames(a)
    #sim BM
    b <- fastBM(mytree,a=OU_BW1$theta,mu=2,sig2=aatmp$sigma %>% mean(),bounds=c(-Inf,Inf),theta=OU_BW1$theta) %>% as.data.frame()
    names(b)[1] <-"exp"
    b$type <- "Brownian"
    b$node <- rownames(b)
    rbind(a,b) %>% cbind(data.frame(stringsAsFactors = F,gene=n))
    
  })%>% rbind.fill()
  
  #best fit saturation
  tree <- as_tibble(mytree)
  sat <- mclapply(mc.cores = 5,1:myall,function(y){
    simuexp0 <- simuexp %>% dplyr::filter(gene==y)
    
    lapply(1:2, function(z){
      exptmp <- simuexp0 %>% dplyr::filter(type==c("OU","Brownian")[z]) 
      distance$exp1 <- exptmp$exp[match(distance$node1,exptmp$node)]
      distance$exp2 <- exptmp$exp[match(distance$node2,exptmp$node)]
      distance$dis <- abs(distance$exp1-distance$exp2)
      tt <- distance %>%
        group_by(mrca,treeDist.asDepth) %>%
        dplyr::summarize(me=mean(dis))
      tt$treeDist.asDepth[which(tt$treeDist.asDepth>8)] <- 8
      
      aa <- tt %>% group_by(treeDist.asDepth,mrca) %>%
        dplyr::summarize(expdis1=mean(me)) %>%
        group_by(treeDist.asDepth) %>%
        dplyr::summarize(expdis2=median(expdis1)) %>%
        as.data.frame()
      
      treedepth <- aa$treeDist.asDepth
      me <- aa$expdis2
      o <-lm(me~1)
      
      seg.fit <- 
        tryCatch({
          xx <- -treedepth
          o2 <- segmented(o,seg.Z=~xx,psi=list(xx=-3))
          if (!is.null(o2)) {
            o2
          } else {
            NULL
          }
        }, error = function(e){
          NULL
        })
      
      ## fit linear 
      liner.lm <- lm(formula = me ~ treedepth, data = aa)
      ## compare three models
      # calculate all AIC values
      # 1. linear model
      aic.lm <- AIC(liner.lm,k=log(length(treedepth)))
      loglik.linear <- logLik(liner.lm) %>% as.data.frame() %>% unlist()
      # 2. no change model
      aic.lm0 <- AIC(o, k=log(length(treedepth)))
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
        aic.seg <- AIC(seg.fit,k=log(length(treedepth)))
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
                       mybreakse=psi.df[1,3])
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
      re.df$type <- unique(exptmp$type)
      return(re.df)
      
    }) %>% rbind.fill() 
  }) %>% rbind.fill()  
  
  sat <- sat %>% cbind(OU_BW1)
  sat$alpha[which(sat$type=="Brownian")] <- 0.5*1.31^(-1)
  sat
  
}) %>% rbind.fill() %>% as.data.frame()

save(simulateOU_BW_Sat,file="~/project/293Tcelllineagetree/rev/res.simulateOU_BW_Sat.V2.1000.Rdata")

# 
load("~/project/293Tcelllineagetree/rev/res.simulateOU_BW_Sat.V2.1000.Rdata")

simulateOU_BW_Sat$type <- factor(simulateOU_BW_Sat$type,levels = c("OU","Brownian"))

simulateOU_BW_Sat %>% 
  group_by(theta,alpha,type) %>% 
  dplyr::summarize(nn=length(which((aic.seg<aic.lm) & (aic.seg<aic.lm0) & segVSliner.lr.chi.P.value<0.05 & segVSzero.lr.chi.P.value <0.05))/length(aic.seg)) %>%
  as.data.frame() %>%
  ggplot(aes(x=theta,y=nn,fill=type,color=type))+
  geom_bar(stat = "identity",position = "dodge")+
  scale_x_continuous(trans = "log2",breaks = c(0.01,0.1,1))+
  #geom_errorbar(aes(ymax=m+se,ymin=m-se),position = "dodge")+
  labs(x="θ",y="Frequency of genes best fitted to saturation model")+
  style.print()