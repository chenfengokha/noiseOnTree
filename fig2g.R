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
library(segmented)



##########3)add leaves
exp <- readRDS("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/adjust_filterMt10_comSCOREandUMI.oneCell.exp.Rds") %>% as.data.frame()
cyclegene <- read.table("/mnt/data/home/lzz/project/2020-6-18-IVDD_scRNA/material/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")
exp <- exp[which(!(rownames(exp) %in% as.vector(cyclegene$V1))),] 
mytree <- read.tree("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/modify_iqtreeInput_rerunIQtree/ran_add_tree1_new.treefile")
tree22 <- as_tibble(mytree) %>% as.data.frame()

#tree depth of each internode
internode <- unique((tree22 %>% dplyr::filter(is.na(label)))$node)
internodedepth <- mclapply(mc.cores = 10,1:length(internode), function(j){
  
  leavesdataall <- offspring(mytree,internode[j])
  depth <- lapply(1:length(leavesdataall),function(z){
    tmpdd <- leavesdataall[z]
    data.frame(z,dep=which(ancestor(mytree, tmpdd)==internode[j]))
  }) %>% rbind.fill() %>%
    dplyr::filter(dep==max(dep)) %>% unique()
  data.frame(treedepth=unique(depth$dep),internode=internode[j],stringsAsFactors = F)
  
}) %>% rbind.fill()

#cusin of new leaf  
newleaf <- unique((tree22 %>% dplyr::filter(substr(label,1,1)=="r"))$label)
oldleaf <- unique((tree22 %>% dplyr::filter(substr(label,(nchar(label)-1),nchar(label))=="_1") %>% dplyr::filter(substr(label,1,1)=="N"))$label)

cusinORsisternode <- mclapply(mc.cores = 10,1:length(newleaf),function(x){
  mynewleaf <- newleaf[x]
  depth.2 <- mclapply(mc.cores = 5,1:length(oldleaf),function(y){
    myoldleaf <- oldleaf[y]
    mynewnode <- tree22$node[which(tree22$label==mynewleaf)]
    myoldnode <- tree22$node[which(tree22$label==myoldleaf)]
    a <- ancestor(mytree, mynewnode)
    b <- ancestor(mytree, myoldnode)
    interab <- a[which(a %in% b)[1]]
    data.frame(stringsAsFactors = F,mynewleaf,myoldleaf,mynewnode,myoldnode,depth=internodedepth$treedepth[internodedepth$internode==interab])
  }) %>% rbind.fill()
  
  # #depth =1 or depth =2 is the new expression
  # if(2 %in% (depth.2$depth)){
  #   tmpnodedata <- (depth.2 %>% dplyr::filter(depth==2))
  #   tmpnode <- tmpnodedata[sample(1:nrow(tmpnodedata),1000,replace = T),] %>% cbind(type=1:1000)
  #   #newexp <- exp[,which(colnames(exp)==tmpnode$myoldleaf)]
  # }else{
  #   tmpnodedata <- (depth.2 %>% dplyr::filter(depth==1))
  #   tmpnode <- tmpnodedata[sample(1:nrow(tmpnodedata),1000,replace = T),] %>% cbind(type=1:1000)
  #   #newexp <- exp[,which(colnames(exp)==tmpnode$myoldleaf)]
  # }
  mclapply(mc.cores = 5,1:1000,function(i){
    if(sample(1:20,1) %in% c(1,4,6,10,15,16)){
      tmpdepth <- sort(unique(depth.2$depth),decreasing = F)[2]
      tmpdepth1 <- (depth.2 %>% dplyr::filter(depth==tmpdepth))
      tmpnodedata <- tmpdepth1[sample(1:nrow(tmpdepth1),1),]
    } else {
      tmpdepth <- sort(unique(depth.2$depth),decreasing = F)[1]
      tmpnodedata <- (depth.2 %>% dplyr::filter(depth==tmpdepth))[1,]
    }
    
    tmpnodedata %>% cbind(data.frame(stringsAsFactors = F,type=i))
  }) %>% rbind.fill()
  
  
}) %>% rbind.fill()

save(cusinORsisternode,file = "~/project/293Tcelllineagetree/data.01.cusinORsisternode.Rdata")

##new expression and treedepth
allleaf <- c(newleaf,oldleaf)
tmprank <- combn(allleaf,2) %>% as.data.frame() %>% t()
testres <- lapply(1:1000, function(i){
  tmpcusin <- cusinORsisternode %>% dplyr::filter(type==i)
  
  newdepth.exp <- mclapply(mc.cores = 50,1:nrow(tmprank),function(x){
    leaf1 <- as.character(tmprank[x,1])
    leaf2 <- as.character(tmprank[x,2])
    a <- ancestor(mytree, tree22$node[which(tree22$label==leaf1)])
    b <- ancestor(mytree, tree22$node[which(tree22$label==leaf2)])
    interab <- a[which(a %in% b)[1]]
    depthab <- internodedepth$treedepth[internodedepth$internode==interab]
    if(substr(leaf1,1,1)=="r"){
      exp1 <- exp[,which(colnames(exp)==tmpcusin$myoldleaf[which(tmpcusin$mynewleaf==leaf1)])]
    } else{
      exp1 <- exp[,which(colnames(exp)==leaf1)]
    }
    if(substr(leaf2,1,1)=="r"){
      exp2 <- exp[,which(colnames(exp)==tmpcusin$myoldleaf[which(tmpcusin$mynewleaf==leaf2)])]
    } else{
      exp2 <- exp[,which(colnames(exp)==leaf2)]
    }
    data.frame(stringsAsFactors = F,expdis=sum(abs(exp1-exp2)^2)^0.5,expR=1-as.numeric(as.vector(cor.test(exp1,exp2,method = "p")$estimate)),treedepth=depthab,internode=interab)
    
  }) %>% rbind.fill()
  save(newdepth.exp,file = paste("~/project/293Tcelllineagetree/data.01.newdepth.exp",i,".Rdata",sep = ""))
  aa <- newdepth.exp
  aa$treedepth[which(aa$treedepth>10)] <- 10
  aa <- aa %>% group_by(treedepth,internode) %>%
    dplyr::summarize(dis1=mean(expdis),R1=mean(expR)) %>%
    group_by(treedepth) %>%
    dplyr::summarize(dis2=median(dis1)/mean(aa$expdis),R2=median(R1)/mean(aa$expR)) %>%
    as.data.frame()
  
  #1)depth-expdis(Euclidean distance)
  treedepth <- aa$treedepth
  med <- aa$dis2
  od <- lm(med~1)
  xxd <- -treedepth
  
  t1 <- summary(segmented(od,seg.Z=~xxd,psi=list(xxd=-3)))
  #1)depth-expdis(Euclidean distance)
  
  meR <- aa$R2
  oR <-lm(meR~1)
  xxR <- -treedepth
  
  tR <- summary(segmented(oR,seg.Z=~xxR,psi=list(xxR=-3)))
  
  data.frame(mybreakD=abs(t1$psi[1,2]),sebreakD=t1$psi[1,3],i,mybreakR=abs(tR$psi[1,2]),sebreakR=tR$psi[1,3],stringsAsFactors = F)
  
}) %>% rbind.fill()


save(testres,file = "~/project/293Tcelllineagetree/data.01.testres.Rdata")

