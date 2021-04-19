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
exp <- readRDS("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/adjust_filterMt10_comSCOREandUMI.oneCell.exp.Rds") %>% as.data.frame()
cyclegene <- read.table("/mnt/data/home/lzz/project/2020-6-18-IVDD_scRNA/material/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")
exp <- exp[which(!(rownames(exp) %in% as.vector(cyclegene$V1))),]
exp$gene <- row.names(exp)
row.names(exp) <- NULL
load("~/lineagepaper/06.distanceexpress.Rdata")
tmp <- distanceexpress[,c(4:5)] %>% unique()
mytree <- read.tree("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/adjust_filterMt10.nwk")
tree <- as_tibble(mytree)
tree2 <- tree %>% as.data.frame()
#all mother cell expression were estimated by all of her daughter cell.
inter1 <- tmp %>% dplyr::filter(treedepth==18)
n <- nrow(inter1)
# tt <- c(1,6,7)
# aa <- 1:n
# # tt <- c(6,9,36,42,43)
# aa2 <- aa[which(!(aa %in% tt))]
nodeexp <- mclapply(mc.cores = 20,1:n, function(x){
  int <- inter1$internode[x]
  son <- tree2 %>% dplyr::filter(parent==int)
  sonexp1 <- exp[,c(ncol(exp),which(colnames(exp) %in% son$label))]
  son2 <- son %>% dplyr::filter(label==1)
  sonexp2 <- internodeexp[,-2] %>% dplyr::filter(intnode %in% son$node) %>%
    reshape2::dcast(gene~intnode,value.var="exp1")
  sonexpraw <- merge(sonexp2, sonexp1, by = "gene",all = T)
  sonexpzhishu <- merge((internodeexp[,-1] %>% dplyr::filter(intnode %in% son$node) %>%
                           reshape2::dcast(gene~intnode,value.var="exp2")), sonexp1, by = "gene",all = T)
  lapply(1:nrow(sonexpraw), function(y){
    gene <- sonexpraw[y,1]
    a <- sonexpraw[y,-1] %>% as.numeric()
    b <- sonexpzhishu[y,-1] %>% as.numeric()
    exp1 <- mean(c(a))
    if(length(which(b==a))==0){
      exp2 <- mean(c(b))
    } else {
      b[which(b==a)] <- 2^(b[which(b==a)])-1
      exp2 <- mean(c(b))
    }
    c <- data.frame(exp1,exp2,gene)
    c$gene <- as.vector(c$gene)
    c
  }) %>% rbind.fill() %>% cbind(intnode=int,depth=unique(inter1$treedepth))
}) %>% rbind.fill()
internodeexp <- internodeexp %>% rbind(nodeexp)
save(internodeexp,file = "~/lineagepaper/06.fig6.allmotherexpress.Rdata")
##
# lapply(1:length(nodeexp), function(x){
#   if(is.null(nrow(nodeexp[[x]]))){
#     a <- data.frame(x,type=0)
#   } else{
#     a <- data.frame(x,type=1)
#   }
#   a
# }) %>% rbind.fill() %>% dplyr::filter(type==0)

###
#inter1 <- tmp %>% dplyr::filter(treedepth==12)
# n <- nrow(inter1)
int <- inter1$internode[n]
son <- tree2 %>% dplyr::filter(parent==int)
sonexp1 <- exp[,c(ncol(exp),which(colnames(exp) %in% son$label))]
sonexp1
son2 <- son %>% dplyr::filter(label==1)
sonexp2 <- internodeexp[,-2] %>% dplyr::filter(intnode %in% son$node) %>%
  reshape2::dcast(gene~intnode,value.var="exp1")

sonexp2$exp1 <- apply(sonexp2[,2:ncol(sonexp2)], 1, mean)

sonexp3 <- internodeexp[,-1] %>% dplyr::filter(intnode %in% son$node) %>%
  reshape2::dcast(gene~intnode,value.var="exp2")
sonexp3$exp2 <- apply(sonexp3[,2:ncol(sonexp3)], 1, mean)
merge(sonexp2[,c(1,ncol(sonexp2))], sonexp3[,c(1,ncol(sonexp3))], by = "gene",all = T) %>% cbind(intnode=int,depth=unique(inter1$treedepth)) -> nodeexp
head(nodeexp)
internodeexp <- internodeexp %>% rbind(nodeexp)
save(internodeexp,file = "~/lineagepaper/06.fig6.allmotherexpress.Rdata")

# sonexpraw <- merge(sonexp2, sonexp1, by = "gene",all = T)
# sonexpzhishu <- merge((internodeexp[,-1] %>% dplyr::filter(intnode %in% son$node) %>%
#                          reshape2::dcast(gene~intnode,value.var="exp2")), sonexp1, by = "gene",all = T)



sonexp1$exp2 <- 2^sonexp1[,2]-1
names(sonexp1)[2] <- "exp1"
sonexp1 %>% cbind(intnode=int,depth=unique(inter1$treedepth)) -> nodeexp
internodeexp <- internodeexp %>% rbind(nodeexp)
save(internodeexp,file = "~/lineagepaper/06.fig6.allmotherexpress.Rdata")


# sonexp1$exp2 <- 2^sonexp1[,2]-1
# names(sonexp1)[2] <- "exp1"

lapply(1:nrow(sonexp1), function(y){
  gene <- sonexp1[y,1]
  a <- sonexp1[y,-1] %>% as.numeric()

  c <- data.frame(exp1=mean(c(a)),exp2=mean(c(2^a-1)),gene)
  c$gene <- as.vector(c$gene)
  c
}) %>% rbind.fill() %>% cbind(intnode=int,depth=unique(inter1$treedepth)) -> nodeexp


##2)border distance with mean using 492 saturated genes
load("~/lineagepaper/06.fig6.allmotherexpress.Rdata")
load("~/lineagepaper/06.gene_slope.all.mean0.1.Rdata")
gene_slope %>% dplyr::filter(seg_lm<0.05) -> gene_slope
internodeexp <- internodeexp %>% dplyr::filter(gene %in% gene_slope$gene)
load("~/lineagepaper/06.distanceexpress3.Rdata")
tmp <- distanceexpress %>% separate(3,c("a","b"))
mytree <- read.tree("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/adjust_filterMt10.nwk")
tree <- as_tibble(mytree)
tree2 <- tree %>% as.data.frame()
tmp$node1 <- tree2$label[match(tmp$a,tree2$node)]
tmp$node2 <- tree2$label[match(tmp$b,tree2$node)]


exp <- readRDS("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/adjust_filterMt10_comSCOREandUMI.oneCell.exp.Rds") %>% as.data.frame()

exp <- exp[(rownames(exp) %in% gene_slope$gene),]

exp$mean <- apply(exp[,1:ncol(exp)],1,mean) %>% as.numeric()

exp$ma <- apply(exp[,1:(ncol(exp)-2)], 1, max)
lowbottom <- min(exp$mean)

exp11 <- exp[,c(874:875)] %>% arrange(mean)
exp11$gene <- rownames(exp11)
seqexpmean <- seq(1,492,by=25)[1:10]
internodeexp$mean <- exp11$mean[match(internodeexp$gene,rownames(exp11))]


seqbor <- mclapply(mc.cores = 3,1:length(seqexpmean),function(y){
  interexp1 <- internodeexp %>% dplyr::filter(mean >= exp11$mean[seqexpmean[y]])
  mclapply(mc.cores = 20,1:1000, function(i){
    gene1 <- sample(unique(interexp1$gene),length(unique(interexp1$gene)),replace = T)
    interexp <- interexp1 %>% dplyr::filter(gene %in% gene1)
    exptmp <- exp %>% dplyr::filter(rownames(exp) %in% gene1)
    exptmp1 <- exptmp[,c(1:873)]
    mydistanc <- mclapply(mc.cores = 1,1:nrow(tmp),function(z){
      tmp1 <- tmp[z,]
      expressa <- exptmp1[,which(colnames(exptmp1)==tmp1$node1)] %>% as.numeric()
      expressb <- exptmp1[,which(colnames(exptmp1)==tmp1$node2)] %>% as.numeric()
      tmp3 <- sum(abs(expressa-expressb)^2)^0.5

      data.frame(expdis=tmp3) %>% cbind(tmp1[,c(6:7)])
    }) %>% rbind.fill() %>% group_by(treedepth,internode) %>%
      dplyr::summarize(me=mean(expdis),.groups = 'drop')

    bordis <- mclapply(mc.cores = 1,1:nrow(mydistanc),function(x){
      a <- mydistanc[x,]
      moexp <- interexp %>% dplyr::filter(intnode==a$internode)
      allexp <- merge(moexp, exp11, by = "gene",all = T) %>%
        dplyr::filter(!is.na(exp1)) %>%
        group_by(gene) %>%
        dplyr::mutate(bordistoma=abs(exp1-ma))

      ma <- allexp$bordistoma
      bordistancetoma <- 10^(sum(log(ma+lowbottom,10))/length(ma))

      data.frame(bordistancetoma) %>% cbind(a %>% as.data.frame())
    }) %>% rbind.fill()
    # ##2) border distance ~ expression distance, 492 genes
    #
    # save(bordis,file = "~/lineagepaper/06.bordis.expdis.492genes.Rdata")
    # load("~/lineagepaper/06.bordis.expdis.492genes.Rdata")
    # source("~/Rfunction/style.print.R")
    # t2 <- cor.test(bordis$bordistancetoma,bordis$me,method = "s")
    # p2 <- bordis %>%
    #   ggplot(aes(bordistancetoma,me))+
    #   geom_point(alpha=0.3)+
    #   scale_x_continuous(limits = c(1.8,2.127))+
    #   #labs(x="to max border",y="pairwise expression distance",title = paste("rho = ",round(as.numeric(t2$estimate),3),"P = ",round(as.numeric(t2$p.value),3)))+
    #   style.print()

    t1 <- cor.test(bordis$bordistancetoma,bordis$me,method = "s")
    data.frame(arho=as.numeric(t1$estimate),ap=t1$p.value,type=y,stringsAsFactors = F)
  }) %>% rbind.fill()

}) %>% rbind.fill()
save(seqbor,file = "~/lineagepaper/06.seqbor.492genes.10part.1000.Rdata")

#3)mother-border-distance ~ son-mean of pairwise distance, simulation
##mother express random 1000
exp <- readRDS("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/adjust_filterMt10_comSCOREandUMI.oneCell.exp.Rds") %>% as.data.frame()
cyclegene <- read.table("/mnt/data/home/lzz/project/2020-6-18-IVDD_scRNA/material/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")
exp <- exp[which(!(rownames(exp) %in% as.vector(cyclegene$V1))),]
exp$gene <- row.names(exp)
row.names(exp) <- NULL
load("~/lineagepaper/06.distanceexpress.Rdata")
tmp <- distanceexpress[,c(4:5)] %>% unique()
mytree <- read.tree("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/adjust_filterMt10.nwk")
tree <- as_tibble(mytree)
tree2 <- tree %>% as.data.frame()
#all mother cell expression were estimated by all of her daughter cell.
tmp0 <- data.frame(raw=colnames(exp)[1:873])
tmp1 <- lapply(1:1000, function(x){
  data.frame(v=sample(tmp0[1,],873,replace = F)) %>% t() %>% as.data.frame()
}) %>% rbind.fill()
cbind(tmp0,tmp1 %>% t() %>% as.data.frame()) -> tmp2
save(tmp2,file = "~/lineagepaper/06.fig6.borderdistanc.tmp2.Rdata")
myexp <- exp[,1:873]

inter1 <- tmp %>% dplyr::filter(treedepth==18)
n <- nrow(inter1)

nodeexp <- mclapply(mc.cores = 1,1:n, function(x){
  int <- inter1$internode[x]
  son <- tree2 %>% dplyr::filter(parent==int)
  mclapply(mc.cores = 7,1:1000,function(y){
    tmp3 <- tmp2[,c(1,y+1)]
    names(tmp3)[2] <- "V"
    aa <- as.vector(tmp3$V)[match(colnames(myexp),as.vector(tmp3$raw))]
    colnames(myexp) <- aa
    myexpress1 <- myexp[,which(c(colnames(myexp) %in% son$label))] %>% as.data.frame()
    myexpress2 <- (internodeexp %>% dplyr::filter((int %in% (son %>% dplyr::filter(label=="1"))$node) & (type==y)))[,c(-1,-2)] %>% t() %>% as.data.frame()
    myexpress <- cbind(myexpress1,myexpress2)
    expmean <- mclapply(mc.cores = 7,1:nrow(myexpress), function(z){
      aa <- myexpress[z,] %>% as.numeric() %>% c() %>% mean()
      data.frame(aa)
    }) %>% rbind.fill()
    cbind(data.frame(int,type=y),expmean %>% t() %>% as.data.frame())

  }) %>% rbind.fill()

})
aa <- lapply(1:length(nodeexp), function(x){
  nodeexp[[x]]
}) %>% rbind.fill()

internodeexp <- internodeexp %>% rbind(aa)
save(internodeexp,file = "~/lineagepaper/06.fig6.allmotherexpress.1000.Rdata")

##4.1)border distance with expdis
load("~/lineagepaper/06.fig6.allmotherexpress.1000.Rdata")
load("~/lineagepaper/06.gene_slope.all.mean0.1.Rdata")
gene_slope %>% dplyr::filter(seg_lm<0.05) -> gene_slope
exp <- readRDS("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/adjust_filterMt10_comSCOREandUMI.oneCell.exp.Rds") %>% as.data.frame()
cyclegene <- read.table("/mnt/data/home/lzz/project/2020-6-18-IVDD_scRNA/material/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")
exp <- exp[which(!(rownames(exp) %in% as.vector(cyclegene$V1))),] 
colnames(internodeexp)[c(-1,-2)] <- rownames(exp)
internodeexp <- internodeexp[,c(1:2,which(colnames(internodeexp) %in% gene_slope$gene))]

load("~/lineagepaper/06.distanceexpress3.Rdata")
tmp <- distanceexpress %>% separate(3,c("a","b"))
mytree <- read.tree("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/adjust_filterMt10.nwk")
tree <- as_tibble(mytree)
tree2 <- tree %>% as.data.frame()
tmp$node1 <- tree2$label[match(tmp$a,tree2$node)]
tmp$node2 <- tree2$label[match(tmp$b,tree2$node)]  

exp <- exp[which(rownames(exp) %in% gene_slope$gene),]
exp$mean <- apply(exp[,1:ncol(exp)],1,mean) %>% as.numeric()

exp$ma <- apply(exp[,1:(ncol(exp)-2)], 1, max)

exp$gene <- rownames(exp)
row.names(exp) <- NULL
exp11 <- exp[,c(874:876)] %>% arrange(mean)
seqexpmean <- seq(1,492,by=25)[1:10]
# internodeexp$mean <- exp11$mean[match(internodeexp$gene,exp11$gene)]
lowbottom <- min(exp11$mean[which(exp11$mean!=0)])
load("~/lineagepaper/06.fig6.borderdistanc.tmp2.Rdata")
seqbor1 <- mclapply(mc.cores = 1,1:10,function(i){
  
  mygene <- exp11[seqexpmean[i]:nrow(exp11),]$gene
  
  aa <- mclapply(mc.cores = 60,1:1000,function(x){
    internodeexp1 <- internodeexp %>% dplyr::filter(type==x)
    internodeexp1 <- internodeexp1[,c(1:2,which(colnames(internodeexp1) %in% mygene))]
    mydata <- mclapply(mc.cores = 1,1:nrow(internodeexp1),function(y){
      mydf <- internodeexp1[y,] %>% t() %>% as.data.frame()
      names(mydf) <- "V1"
      mydf$gene <- rownames(mydf)
      mydf1 <- mydf[c(-1,-2),]
      mydf1$ma <- exp11$ma[match(mydf1$gene,exp11$gene)]
      
      taa <- abs(mydf1$V1-mydf1$ma)
      
      bordis <- 10^(sum(log(taa+lowbottom,10))/length(taa))
      internodeexp1[y,1:2] %>% cbind(bordis)
    }) %>% rbind.fill()
    myexp <- exp[which(exp$gene %in% mygene),]
    myexpdis <- mclapply(mc.cores = 1,1:nrow(tmp),function(z){
      bb <- tmp[z,]
      cc <- tmp2[,c(1,x+1)]
      names(cc)[2] <- "A"
      dd1 <- myexp[,which(names(myexp)==(as.vector(cc$A)[which(as.vector(cc$raw)==bb$node1)]))] %>% as.numeric()
      dd2 <- myexp[,which(names(myexp)==(as.vector(cc$A)[which(as.vector(cc$raw)==bb$node2)]))] %>% as.numeric()
      data.frame(tmp3=sum(abs(dd1-dd2)^2)^0.5) %>% cbind(int=bb$internode)
    }) %>% rbind.fill() %>% group_by(int) %>%
      dplyr::summarize(expdis=mean(tmp3))
    mydata$expdis <- myexpdis$expdis[match(mydata$int,myexpdis$int)]
    t1 <- cor.test(mydata$bordis,mydata$expdis,method = "s")
    
    data.frame(rho=as.numeric(t1$estimate),p=t1$p.value,type=x,stringsAsFactors = F)
  }) %>% rbind.fill() %>% cbind(i)
  
}) %>% rbind.fill()
save(seqbor1,file = "~/lineagepaper/06.fig6.seqbor.492.sample1000.Rdata") 

##draw picture
load("~/lineagepaper/06.seqbor.492genes.10part.1000.Rdata")
load("~/lineagepaper/06.fig6.seqbor.492.sample1000.Rdata")
seqbor %>% group_by(type) %>% dplyr::summarize(mean=mean(arho),sd=sd(arho)) -> seqbor0
ggplot()+
  geom_point(data = seqbor0,aes(x=type,y=mean),color="#F8766D")+
  geom_errorbar(data=seqbor0,aes(x=type,y=mean,ymin=mean-sd,ymax=mean+sd), width = 0.25,color="#F8766D")+
  #geom_violin(data = seqbor,aes(x=type,y=arho,group=type))+
  geom_violin(data = seqbor1,aes(x=i,y=rho,group=i),color="#00BFC4")+
  geom_boxplot(data = seqbor1,aes(x=i,y=rho,group=i),width = 0.1,color="#00BFC4", outlier.colour = NA)+
  scale_x_continuous(breaks = c(1,3,5,7,9))+
  scale_y_continuous(limits = c(-0.15,0.95),breaks = c(0,0.3,0.6,0.9))+
  style.print()








