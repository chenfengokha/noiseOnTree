library(plyr)
library(dplyr)
library(parallel)
library(ape)
library(tidytree)
library(ggtree)
library(phyloTop)
library(rpart)
library(stringr)
library(readr)
library(scales)
library(data.table)
library(tidyr)
#1)mother-border-distance ~ son-mean of pairwise distance
#sample1
exp <- readRDS("~/project/293Tcelllineagetree/plot/test/adjust_filterMt10_comSCOREandUMI.oneCell.exp.Rds") %>% as.data.frame()
mytree <- read.tree("~/project/293Tcelllineagetree/plot/test/adjust_filterMt10.nwk")
cfraw <- readRDS("~/project/293Tcelllineagetree/plot/sampleA_raw_exp_depth_euclidean_dist_nocc.Rds")$real.dist
#sample2
exp <- readRDS("~/project/293Tcelllineagetree/plot/test/CF2_oneCellexp.Rds") %>% as.data.frame()
mytree <- read.tree("~/project/293Tcelllineagetree/plot/test/CF2.union.adjust.nwk")
cfraw <- readRDS("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/new_samples_include10x_adjust_results/rawBC/CF2/CF2_depth_euclidean_exp_dist_nocc.Rds")$real.dist
#saturated gene
olgene <- read.csv("~/saturatedgeneAC.csv",stringsAsFactors = F,sep = " ")

exp <- exp[which(rownames(exp) %in% olgene$gene),]

tree <- as_tibble(mytree)
tree2 <- tree %>% as.data.frame()

#all mother cell expression were estimated by all of her daughter cell.
inter1 <- cfraw[,c(3,5)] %>% unique() %>% dplyr::filter(treeDist.asDepth==13)
n <- nrow(inter1)
nodeexp <- mclapply(mc.cores = 60,1:n, function(x){
  int <- inter1$mrca[x]
  son <- tree2 %>% dplyr::filter(parent==int)
  myexp <- exp
  myexpress1 <- myexp[,which(c(colnames(myexp) %in% son$label))] %>% as.data.frame()
  myexpress2 <- (internodeexp %>% dplyr::filter((int %in% (son %>% dplyr::filter(label=="1"))$node)))[,c(-1)] %>% t() %>% as.data.frame()
  myexpress <- cbind(myexpress1,myexpress2)
  #myexpress <- myexpress1
  expmean <- mclapply(mc.cores = 1,1:nrow(myexpress), function(z){
    aa <- myexpress[z,] %>% as.numeric() %>% c() %>% mean()
    data.frame(aa,stringsAsFactors = F)
  }) %>% rbind.fill()
  cbind(data.frame(int,stringsAsFactors = F),expmean %>% t() %>% as.data.frame())
  
})

aa <- lapply(1:length(nodeexp), function(x){
  nodeexp[[x]]
}) %>% rbind.fill()
names(aa)[c(-1)] <- rownames(exp)
internodeexp <- internodeexp %>% rbind(aa)
#names(internodeexp)[c(-1)] <- rownames(exp)
#saveRDS(internodeexp,file = "~/project/293Tcelllineagetree/plot/test/06.real.motherexp.Acpm.boddisandexpdis.Rds")
#saveRDS(internodeexp,file = "~/project/293Tcelllineagetree/plot/test/06.real.motherexp.cf2cpm.boddisandexpdis.Rds")
#saveRDS(internodeexp,file = "~/project/293Tcelllineagetree/plot/test/06.real.motherexp.A.boddisandexpdis.Rds")
saveRDS(internodeexp,file = "~/project/293Tcelllineagetree/plot/test/06.real.motherexp.cf2.boddisandexpdis.Rds")

##1.1)border distance using saturated genes
exp$mean <- apply(exp[,1:ncol(exp)],1,mean) %>% as.numeric()
exp$ma <- apply(exp[,1:(ncol(exp)-1)], 1, max)

exp$gene <- rownames(exp)
row.names(exp) <- NULL
exp11 <- exp %>% arrange(mean)
seqexpmean <- seq(1,nrow(olgene),by=floor(nrow(olgene)/10))[1:10]

lowbottom <- min(exp11$mean[which(exp11$mean!=0)])

seqbor <- mclapply(mc.cores = 1,1:10,function(i){
  
  mygene <- exp11[seqexpmean[i]:nrow(exp),]$gene
  mclapply(mc.cores = 60,1:1000,function(x){
    mygene <- sample(mygene,1000,replace = T)
    
    internodeexp1 <- internodeexp 
    internodeexp1 <- internodeexp1[,c(1,which(colnames(internodeexp1) %in% mygene))]
    mydata <- mclapply(mc.cores = 1,1:nrow(internodeexp1),function(y){
      mydf <- internodeexp1[y,] %>% t() %>% as.data.frame()
      names(mydf) <- "V1"
      mydf$gene <- rownames(mydf)
      mydf1 <- mydf[c(-1),]
      mydf1$ma <- exp11$ma[match(mydf1$gene,exp11$gene)]
      
      taa <- abs(mydf1$V1-mydf1$ma)/mydf1$ma
      
      bordis <- 10^(sum(log(taa+lowbottom,10))/length(taa))
      data.frame(int=internodeexp1[y,1],bordis,stringsAsFactors = F)
      
    }) %>% rbind.fill()
    
    
    myexp <- exp[which(exp$gene %in% mygene),]
    myexpdis <- mclapply(mc.cores =1,1:nrow(cfraw),function(z){
      bb <- cfraw[z,]
      
      dd1 <- myexp[,which(names(myexp)==bb$node1)] %>% as.numeric()
      dd2 <- myexp[,which(names(myexp)==bb$node2)] %>% as.numeric()
      data.frame(tmp3=sum(abs(dd1-dd2)^2)^0.5,stringsAsFactors = F,int=bb$mrca)
    }) %>% rbind.fill() %>% group_by(int) %>%
      dplyr::summarize(expdis=mean(tmp3))
    mydata$expdis <- myexpdis$expdis[match(mydata$int,myexpdis$int)]
    t1 <- cor.test(mydata$bordis,mydata$expdis,method = "s")
    
    data.frame(rho=as.numeric(t1$estimate),p=t1$p.value,i,x,stringsAsFactors = F)
  }) %>% rbind.fill()
  
}) %>% rbind.fill()
#saveRDS(mydata,file = "~/project/293Tcelllineagetree/plot/test/06.real.all.cor.A.boddisandexpdis.Rds") 
#saveRDS(mydata,file = "~/project/293Tcelllineagetree/plot/test/06.real.all.cor.cf2.boddisandexpdis.Rds")

#saveRDS(seqbor,file = "~/project/293Tcelllineagetree/plot/test/06.real.cor.A.boddisandexpdis.Rds") 
saveRDS(seqbor,file = "~/project/293Tcelllineagetree/plot/test/06.real.cor.cf2.boddisandexpdis.Rds") 
seqbor$i <- factor(seqbor$i,levels = c(1:10))
source("~/Rfunction/style.print.R")
seqbor %>% ggplot(aes(x=i,y=rho,group=i))+geom_violin()+
  geom_boxplot(width = 0.1, outlier.colour = NA)+
  style.print()


#4)mother-border-distance ~ son-mean of pairwise distance, simulation
##mother express random 1000
#sample 1
exp <- readRDS("~/project/293Tcelllineagetree/plot/test/adjust_filterMt10_comSCOREandUMI.oneCell.exp.Rds") %>% as.data.frame()
mytree <- read.tree("~/project/293Tcelllineagetree/plot/test/adjust_filterMt10.nwk")

cfraw <- readRDS("~/project/293Tcelllineagetree/plot/sampleA_raw_exp_depth_euclidean_dist_nocc.Rds")$real.dist
#sample 2
exp <- readRDS("~/project/293Tcelllineagetree/plot/test/CF2_oneCellexp.Rds") %>% as.data.frame()
mytree <- read.tree("~/project/293Tcelllineagetree/plot/test/CF2.union.adjust.nwk")
cfraw <- readRDS("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/new_samples_include10x_adjust_results/rawBC/CF2/CF2_depth_euclidean_exp_dist_nocc.Rds")$real.dist
#overlapgene
olgene <- read.csv("~/saturatedgeneAC.csv",stringsAsFactors = F,sep = " ")

exp <- exp[which(rownames(exp) %in% olgene$gene),]
tmp <- cfraw[,c(3,5)] %>% unique()

tree <- as_tibble(mytree)
tree2 <- tree %>% as.data.frame()
#all mother cell expression were estimated by all of her daughter cell.
tmp0 <- data.frame(raw=colnames(exp)[1:873],stringsAsFactors = F)
tmp1 <- lapply(1:1000, function(x){
  data.frame(v=sample(tmp0[,1],873,replace = F),stringsAsFactors = F) %>% t() %>% as.data.frame()
}) %>% rbind.fill()
cbind(tmp0,tmp1 %>% t() %>% as.data.frame()) -> tmp2
#saveRDS(tmp2,file = "~/project/293Tcelllineagetree/plot/test/06.random.rank.A.boddisandexpdis.Rds")
saveRDS(tmp2,file = "~/project/293Tcelllineagetree/plot/test/06.random.rank.cf2.boddisandexpdis.Rds")

tmp$treeDist.asDepth %>% max()
inter1 <- tmp %>% dplyr::filter(treeDist.asDepth==13)
n <- nrow(inter1)

nodeexp <- mclapply(mc.cores = 1,1:n, function(x){
  int <- inter1$mrca[x]
  son <- tree2 %>% dplyr::filter(parent==int)
  mclapply(mc.cores = 60,1:1000,function(y){
    myexp <- exp
    tmp3 <- tmp2[,c(1,y+1)]
    names(tmp3)[2] <- "V"
    aa <- as.vector(tmp3$V)[match(colnames(myexp),as.vector(tmp3$raw))]
    colnames(myexp) <- aa
    myexpress1 <- myexp[,which(c(colnames(myexp) %in% son$label))] %>% as.data.frame()
    myexpress2 <- (internodeexp %>% dplyr::filter((int %in% (son %>% dplyr::filter(label=="1"))$node) & (type==y)))[,c(-1,-2)] %>% t() %>% as.data.frame()
    myexpress <- cbind(myexpress1,myexpress2)
    #myexpress <- myexpress1
    expmean <- mclapply(mc.cores = 1,1:nrow(myexpress), function(z){
      aa <- myexpress[z,] %>% as.numeric() %>% c() %>% mean()
      data.frame(aa,stringsAsFactors = F)
    }) %>% rbind.fill()
    cbind(data.frame(int,type=y,stringsAsFactors = F),expmean %>% t() %>% as.data.frame())
    
  }) %>% rbind.fill()
  
})

aa <- lapply(1:length(nodeexp), function(x){
  nodeexp[[x]]
}) %>% rbind.fill()
names(aa)[c(-1,-2)] <- rownames(exp)
internodeexp <- internodeexp %>% rbind(aa)

#saveRDS(internodeexp,file = "~/project/293Tcelllineagetree/plot/test/06.random.motherexp.A.boddisandexpdis.Rds")
saveRDS(internodeexp,file = "~/project/293Tcelllineagetree/plot/test/06.random.motherexp.cf2.boddisandexpdis.Rds")

##4.1)border distance

exp$mean <- apply(exp[,1:ncol(exp)],1,mean) %>% as.numeric()

exp$ma <- apply(exp[,1:(ncol(exp)-1)], 1, max)

exp$gene <- rownames(exp)
row.names(exp) <- NULL
exp11 <- exp %>% arrange(mean)
seqexpmean <- seq(1,nrow(olgene),by=floor(nrow(olgene)/10))[1:10]
# internodeexp$mean <- exp11$mean[match(internodeexp$gene,exp11$gene)]
lowbottom <- min(exp11$mean[which(exp11$mean!=0)])
#load("~/lineagepaper/06.fig6.borderdistanc.tmp2.Rdata")
seqbor1 <- mclapply(mc.cores = 1,1:10,function(i){
  
  mygene <- exp11[seqexpmean[i]:nrow(exp),]$gene
  
  aa <- mclapply(mc.cores = 60,1:1000,function(x){
    internodeexp1 <- internodeexp %>% dplyr::filter(type==x)
    internodeexp1 <- internodeexp1[,c(1:2,which(colnames(internodeexp1) %in% mygene))]
    mydata <- mclapply(mc.cores = 1,1:nrow(internodeexp1),function(y){
      mydf <- internodeexp1[y,] %>% t() %>% as.data.frame()
      names(mydf) <- "V1"
      mydf$gene <- rownames(mydf)
      mydf1 <- mydf[c(-1,-2),]
      mydf1$ma <- exp11$ma[match(mydf1$gene,exp11$gene)]
      
      taa <- abs(mydf1$V1-mydf1$ma)/mydf1$ma
      
      bordis <- 10^(sum(log(taa+lowbottom,10))/length(taa))
      internodeexp1[y,1:2] %>% cbind(bordis)
    }) %>% rbind.fill()
    myexp <- exp[which(exp$gene %in% mygene),]
    myexpdis <- mclapply(mc.cores = 1,1:nrow(cfraw),function(z){
      bb <- cfraw[z,]
      cc <- tmp2[,c(1,x+1)]
      names(cc)[2] <- "A"
      dd1 <- myexp[,which(names(myexp)==(as.vector(cc$A)[which(as.vector(cc$raw)==bb$node1)]))] %>% as.numeric()
      dd2 <- myexp[,which(names(myexp)==(as.vector(cc$A)[which(as.vector(cc$raw)==bb$node2)]))] %>% as.numeric()
      data.frame(tmp3=sum(abs(dd1-dd2)^2)^0.5,stringsAsFactors = F) %>% cbind(int=bb$mrca)
    }) %>% rbind.fill() %>% group_by(int) %>%
      dplyr::summarize(expdis=mean(tmp3))
    mydata$expdis <- myexpdis$expdis[match(mydata$int,myexpdis$int)]
    t1 <- cor.test(mydata$bordis,mydata$expdis,method = "s")
    
    data.frame(rho=as.numeric(t1$estimate),p=t1$p.value,type=x,stringsAsFactors = F)
  }) %>% rbind.fill() %>% cbind(i)
  
}) %>% rbind.fill()
#saveRDS(seqbor1,file = "~/project/293Tcelllineagetree/plot/test/06.random.cor.A.boddisandexpdis.Rds") 
saveRDS(seqbor1,file = "~/project/293Tcelllineagetree/plot/test/06.random.cor.cf2.boddisandexpdis.Rds") 

###########################################
##fig4.B
#sample1
sampA <- readRDS("~/project/293Tcelllineagetree/plot/test/06.real.all.cor.A.boddisandexpdis.Rds") %>%
  group_by(int,bordis) %>%
  dplyr::summarize(m=mean(expdis)) %>% as.data.frame()
#sample2
samp2 <- readRDS("~/project/293Tcelllineagetree/plot/test/06.real.all.cor.cf2.boddisandexpdis.Rds") %>%
  group_by(int,bordis) %>%
  dplyr::summarize(m=mean(expdis)) %>% as.data.frame()
bordA <- readRDS("~/project/293Tcelllineagetree/plot/test/06.real.motherexp.Anew.bordis.Rds")
bord2 <-  readRDS("~/project/293Tcelllineagetree/plot/test/06.real.motherexp.cf2new.bordis.Rds")

sampA$newdis <- bordA$bordis[match(sampA$int,bordA$int)]
samp2$newdis <- bord2$bordis[match(samp2$int,bord2$int)]
cor.test(sampA$newdis,sampA$m,method = "s")$p.value
cor.test(samp2$newdis,samp2$m,method = "s")$p.value

ggplot()+
  geom_point(data=sampA,aes(newdis,m),alpha=0.5,color="#F8766D") +style.print()
ggplot()+  
  geom_point(data=samp2,aes(newdis,m),alpha=0.5,color="#00BFC4") +
  style.print()
##fig4.c
#sample1
realA <- readRDS("~/project/293Tcelllineagetree/plot/test/06.real.cor.A.boddisandexpdis.Rds") %>% group_by(i) %>% dplyr::summarize(mean=mean(rho),sd=sd(rho))
ranA <- readRDS("~/project/293Tcelllineagetree/plot/test/06.random.cor.A.boddisandexpdis.Rds")
#sample2
realcf2 <- readRDS("~/project/293Tcelllineagetree/plot/test/06.real.cor.cf2.boddisandexpdis.Rds") %>% group_by(i) %>% dplyr::summarize(mean=mean(rho),sd=sd(rho))
rancf2 <- readRDS("~/project/293Tcelllineagetree/plot/test/06.random.cor.cf2.boddisandexpdis.Rds")

real <- rbind(cbind(realA,type2="A"),cbind(realcf2,type2="cf2")) %>% dplyr::filter(i!=10)
real$i <- factor(real$i,levels = as.character(1:9))
ran <- rbind(cbind(ranA,type2="A"),cbind(rancf2,type2="cf2")) %>% dplyr::filter(i!=10)

ran$i <- factor(ran$i,levels = as.character(1:9))

ggplot()+
  geom_point(data = real,aes(x=i,y=mean,color=type2,group=type2))+
  geom_errorbar(data=real,aes(x=i,y=mean,ymin=mean-sd,ymax=mean+sd,color=type2,group=type2), width = 0.25)+
  geom_violin(data = ran,aes(x=i,y=rho,color=type2),width=0.8)+
  geom_boxplot(data = ran,aes(x=i,y=rho,color=type2),width=0.4, outlier.colour = NA)+
  scale_y_continuous(limits = c(-0.2,0.95),breaks = c(0,0.4,0.8))+
  style.print()



