suppressMessages({
  library(velocyto.R);
  library(tidyverse);
  library(pagoda2);
  library(parallel);
  library(matlib);
  library(ggtree);
  library(ape)
})
## ========  some functions =======
get.inTree.mx.oneCell <- function(info.df, raw.mx){
  NodeName <- info.df$NodeName[match(colnames(raw.mx), info.df$BC)]
  colnames(raw.mx) <- NodeName
  re.mx <- raw.mx[, grepl("_1$", NodeName)]
  return(re.mx)
}

##  calculate every cell distance to boundary
cell.distTObundary <- function(exp.matrix){
  max.boundary <- apply(exp.matrix, 1, max) 
  pseudo.count <- min(rowMeans(exp.matrix))
  node.dist <- apply(exp.matrix, 2, function(x) 10^(sum(log(abs(x-max.boundary)+pseudo.count,10))/length(max.boundary)))
  node.dist.df <- node.dist %>% as.data.frame()
  colnames(node.dist.df) <- "dist"
  return(node.dist.df)
}

## from velocity result calculate angle
raw.arrow.angle <- function(rvl.re, pair.df){
  current.mx <- rvl.re$current
  project.mx <- rvl.re$projected
  diff.mx <- current.mx - project.mx
  all.info <- mclapply(seq(nrow(pair.df)), function(i){
    n1 <- pair.df$node1[i]
    n2 <- pair.df$node2[i]
    ## 
    n1.project <- project.mx[,n1] %>% as.vector()
    n2.project <- project.mx[,n2] %>% as.vector()
    n1.diff <- diff.mx[,n1] %>% as.vector()
    n2.diff <- diff.mx[,n2] %>% as.vector()
    pair.angle <- angle(n1.diff, n2.diff)
    arrow.dist <- 10^(sum(log(abs(n1.project-n2.project) + 0.001, 10))/100)
    return(data.frame(node1 = n1,
                      node2 = n2,
                      angle = pair.angle,
                      project.dist = arrow.dist,
                      stringsAsFactors = F))
  }, mc.cores=60) %>% bind_rows()
  return(all.info)
}

raw.angle.dist <- function(pair.angle, cell.dist.df, mother.dist.df){
  pair.angle$node1.dist <- cell.dist.df[pair.angle$node1, "dist"]
  pair.angle$node2.dist <- cell.dist.df[pair.angle$node2, "dist"]
  pair.angle$mother.dist <- mother.dist.df$mother.dist[match(pair.angle$node1, mother.dist.df$label)]
  return(pair.angle)
}


## compare border vs inside angles                    
bVSi.angle.cut <- function(pair.angle, all.dist, quantile.cut = seq(0.05,1,0.05), border.cut = seq(0.05,0.5,0.05)){
  # 1. get all distance to border
  quantile.dist <- quantile(all.dist, probs = quantile.cut)
  #print(quantile.dist)
  re.df <- lapply(seq(length(border.cut)), function(i){
    toBorder <- paste0(border.cut[i]*100, "%")
    toInside <- paste0((1-border.cut[i])*100, "%")
    border.info <- pair.angle %>% filter(node1.dist <= quantile.dist[toBorder] & node2.dist <= quantile.dist[toBorder])
    inside.info <- pair.angle %>% filter(node1.dist >= quantile.dist[toInside] & node2.dist >= quantile.dist[toInside])
    border.angle.mean <- mean(border.info$angle)
    inside.angle.mean <- mean(inside.info$angle)
    # ks test and wilcox test
    ks.p <- ks.test(border.info$angle, inside.info$angle, alternative = "greater")$p.value
    wilcox.p <- wilcox.test(border.info$angle, inside.info$angle, alternative = "less")$p.value
    #compare.ratio <- (border.len.mean - inside.len.mean) / inside.len.mean
    return(data.frame(cut.to.border = toBorder,
                      cut.to.inside = toInside,
                      border.num = nrow(border.info),
                      inside.num = nrow(inside.info),
                      border.angle.mean = border.angle.mean,
                      inside.angle.mean = inside.angle.mean,
                      ks.p = ks.p,
                      wilcox.p = wilcox.p))
  }) %>% bind_rows()
  return(re.df)
}


mother.bVSi.angle.cut <- function(pair.angle, all.dist, quantile.cut = seq(0.05,1,0.05), border.cut = seq(0.05,0.5,0.05)){
  # 1. get all distance to border
  quantile.dist <- quantile(all.dist, probs = quantile.cut)
  #print(quantile.dist)
  re.df <- lapply(seq(length(border.cut)), function(i){
    toBorder <- paste0(border.cut[i]*100, "%")
    toInside <- paste0((1-border.cut[i])*100, "%")
    border.info <- pair.angle %>% filter(mother.dist <= quantile.dist[toBorder])
    inside.info <- pair.angle %>% filter(mother.dist >= quantile.dist[toInside])
    border.angle.mean <- mean(border.info$angle)
    inside.angle.mean <- mean(inside.info$angle)
    # ks test and wilcox test
    ks.p <- ks.test(border.info$angle, inside.info$angle, alternative = "greater")$p.value
    wilcox.p <- wilcox.test(border.info$angle, inside.info$angle, alternative = "less")$p.value
    #compare.ratio <- (border.len.mean - inside.len.mean) / inside.len.mean
    return(data.frame(cut.to.border = toBorder,
                      cut.to.inside = toInside,
                      border.num = nrow(border.info),
                      inside.num = nrow(inside.info),
                      border.angle.mean = border.angle.mean,
                      inside.angle.mean = inside.angle.mean,
                      ks.p = ks.p,
                      wilcox.p = wilcox.p))
  }) %>% bind_rows()
  return(re.df)
}




#  ======== load global vars ========
sampleA.loom <- read.loom.matrices("/mnt/data/home/lzz/project/2019-7-31-10xSampleA/results/293T_20190731sampleA/velocyto/293T_20190731sampleA.loom")
alleleInfo <- read.csv("~/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/adjust_filterMt10.AllelesInfo.txt", 
                       sep = '\t', stringsAsFactors = F)
oneCellNode.pair <- read.csv("~/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/adjust_filterMt10_comSCOREandUMI.oneCellNode.pair_Dist.txt",
                            sep = '\t', stringsAsFactors = F)

sisterPair <- filter(oneCellNode.pair, treeDist == 1)



#####
emat <- sampleA.loom$spliced
rownames(emat) <- make.unique(rownames(emat))
nmat <- sampleA.loom$unspliced
rownames(nmat) <- make.unique(rownames(nmat))
# 1. change emat and nmat colname
colnames(emat) <- gsub("293T_20190731sampleA:|x", "", colnames(emat))
colnames(nmat) <- gsub("293T_20190731sampleA:|x", "", colnames(nmat)) 
# 2. change to node name and get results
emat.onecell <- get.inTree.mx.oneCell(alleleInfo, emat)
nmat.onecell <- get.inTree.mx.oneCell(alleleInfo, nmat)



# =============  run ====================
# load 492 seuration genes
s492.gene <- readLines("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/material/11.txt")
length(s492.gene)

emat.492 <- emat.onecell[rownames(emat.onecell) %in% s492.gene,]
rownames(emat.492) <- make.unique(rownames(emat.492))
nmat.492 <- nmat.onecell[rownames(nmat.onecell) %in% s492.gene,]
rownames(nmat.492) <- make.unique(rownames(nmat.492))

noPCA492ve <- gene.relative.velocity.estimates(emat=emat.492,
                                                nmat=nmat.492,
                                               deltaT=1,
                                               #fit.quantile = 1,
                                                kCells=1,
                                                #cell.dist=cell.dist.710
                                               min.nmat.emat.slope = -999,
                                              min.nmat.emat.correlation = -1
                                                )


### sister pair border vs inside 
# load boundary
load("/mnt/data/home/chenfeng/lineagepaper/06.bordis.expdis.492genes.Rdata")
# bordis
tree <- read.tree("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/adjust_filterMt10.nwk")
tree_tibble <- as_tibble(tree)
node.mother <- filter(tree_tibble, grepl("_1$", label)) %>% select(label, parent)
## keep only have onecell node pair mother
node.mother <- node.mother %>% filter(label %in% c(sisterPair$node1, sisterPair$node2))
node.mother$mother.dist <- bordis$bordistancetoma[match(node.mother$parent, bordis$internode)]

sisterPair.dist <- mclapply(seq(nrow(sisterPair)), function(i){
  sub.df <- sisterPair[i,]
  n1 <- sub.df$node1[1]
  n2 <- sub.df$node2[1]
  # get mother dist info
  n1.m.info <- filter(node.mother, label == n1)
  n2.m.info <- filter(node.mother, label == n2)
  # return result
  return(data.frame(node1 = n1,
                    node2 = n2,
                    n1.mother = n1.m.info$parent[1],
                    n2.mother = n2.m.info$parent[1],
                    n1.mother.dist = n1.m.info$mother.dist[1],
                    n2.mother.dist = n2.m.info$mother.dist[1],
                    stringsAsFactors = F))
}, mc.cores=60) %>% bind_rows()

noPCA492ve.angle <- raw.arrow.angle(noPCA492ve, sisterPair.dist)

head(noPCA492ve.angle)

noPCA492ve.angle.info <- left_join(noPCA492ve.angle, sisterPair.dist)

noPCA492ve.angle.info <- noPCA492ve.angle.info %>% group_by(n1.mother) %>% do({
  mydf <- .
  mydf$mother.mean.angle <- mean(mydf$angle)
  mydf
})

noPCA492ve.angle.gy.mother <- noPCA492ve.angle.info %>% group_by(n1.mother) %>% summarise(mother.dist = mean(n1.mother.dist), angle.mean = mean(angle))

cor.test(noPCA492ve.angle.gy.mother$mother.dist, noPCA492ve.angle.gy.mother$angle.mean, method = "spearman")
cor.test(noPCA492ve.angle.gy.mother$mother.dist, noPCA492ve.angle.gy.mother$angle.mean, method = "pearson")

pdf("~/project/2019-8-22-PacBio.SampleA/paper_writing/velocity_results/noPCA492ve_492boundary_dist_angle.pdf")
options(repr.plot.width=14, repr.plot.height=10)
noPCA492ve.angle.gy.mother %>% ggplot() +
  geom_point(aes(x=mother.dist, y=angle.mean), alpha = 0.5) +
  annotate("text", x=1.85, y=90, label=paste0("ρ = 0.37\nP = 2.382e-09"), size=5) +
  labs(x="Distance to expression boundary", y= "Angular distance of RNA velocity with noise saturation between sister cells") +
  scale_x_continuous(limits = c(1.8,2.127)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 15)
        #             legend.text = element_text(size=15),
        #             legend.title = element_text(size=15),
        #             legend.position = "right"
  ) 
dev.off()



















