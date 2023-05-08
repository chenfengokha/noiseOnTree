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
#####################1)down sample
suppressMessages({
  library(ape)
  library(ggtree)
  library(tidyverse)
})

suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(parallel)
  library(segmented)
  library(ggtree)
  library(ape)
})

check.mean.dist <- function(raw.dist.df){
  re.depth.mean.expdist <- 
    raw.dist.df %>% group_by(mrca, treeDist.asDepth) %>% 
    summarise(mrca.mean.dist = mean(exp.dist),
              mrca.median.dist = median(exp.dist),
              mrca.sd.dist = sd(exp.dist)) %>%
    group_by(treeDist.asDepth) %>% 
    summarise(mean.mean.dist = mean(mrca.mean.dist),
              mean.median.dist = median(mrca.mean.dist),
              median.mean.dist = mean(mrca.median.dist),
              median.median.dist = median(mrca.median.dist))
  return(re.depth.mean.expdist)
  
}

check.mean.dist.nor <- function(raw.dist.df, ran.dist.df){
  mean.ran.dist <- mean(ran.dist.df$ran.exp.dist)
  median.ran.dist <- median(ran.dist.df$ran.exp.dist)
  re.depth.mean.expdist <- 
    raw.dist.df %>% dplyr::group_by(mrca, treeDist.asDepth) %>% 
    dplyr::summarise(mrca.mean.dist = mean(exp.dist),
              mrca.median.dist = median(exp.dist),
              mrca.sd.dist = sd(exp.dist)) %>%
    dplyr::group_by(treeDist.asDepth) %>% 
    dplyr::summarise(mean.mean.dist = mean(mrca.mean.dist)/mean.ran.dist,
              mean.median.dist = median(mrca.mean.dist)/median.ran.dist,
              median.mean.dist = mean(mrca.median.dist)/mean.ran.dist,
              median.median.dist = median(mrca.median.dist)/median.ran.dist)
  return(re.depth.mean.expdist)
  
}

view.tree.nodes <- function(phylo.obj, select.nodes, view.as = "internal"){
  if (view.as == "internal"){
    phylo.obj2 <- groupClade(phylo.obj, .node=select.nodes, group_name = "diffs")
  } else {
    phylo.obj2 <- groupOTU(phylo.obj, select.nodes, group_name = "diffs")
  }
  p <- 
    ggtree(phylo.obj2, layout = "circular", aes(color = diffs), branch.length = "none") + xlim(-10, NA) +
    scale_color_manual(values=c("black", "red"),
                       label = c("others", "selected")) +
    theme(legend.position = c(0.5, 0.5))
  return(p)
}


## get all subtree depth, and BI score to know how multifurcation it is
# BIscore = (Interal node num + 1) / leaves
# if tree is bifurcate, BIscore is 1
subtree.BIscore <- function(tree_path){
  ## run
  tree <- read.tree(tree_path)
  # get oneCell node pair
  oneCellNode <- tree$tip.label[grepl("_1$", tree$tip.label)]
  oneCellNode.pair <- combn(oneCellNode, 2)
  oneCellNode.pair <- data.frame(node1=oneCellNode.pair[1,], node2=oneCellNode.pair[2,], stringsAsFactors = F)
  
  # get tree tibble
  inSubtrees <- subtrees(tree)
  distMatrix <- dist.nodes(tree)
  treeTibble <- as_tibble(tree)
  
  getDist <- function(nodeName, node){
    n2 <- dplyr::filter(treeTibble, label == nodeName) %>% `[`("node") %>% unlist()
    dis <- distMatrix[node, n2]
    return(dis)
  }
  ## get all subtree depth
  all.subtree.info <- lapply(seq(length(inSubtrees)), function(s){
    subtree <- inSubtrees[[s]]
    s.name <- subtree$name
    # depth
    treeDeep <- Map(getDist, subtree$tip.label, s.name) %>% unlist() %>% as.vector() %>% max()
    # bifurcate score
    bi.score <- (length(subtree$node.label) + 1) / length(subtree$tip.label)
    return(data.frame(subtree = s.name,
                      subtree.depth = treeDeep,
                      bi.score = bi.score,
                      leaves.num = length(subtree$tip.label),
                      stringsAsFactors = F))
  }) %>% bind_rows()
  return(all.subtree.info)
}

treeDepth.dist.no.cc <- function(tree_path, exp_path, 
                                 cc.gene_path = "/mnt/data/home/lzz/project/2020-6-18-IVDD_scRNA/material/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt", 
                                 use.method = "Euclidean",
                                 use.cores=40){
  ## run
  tree <- read.tree(tree_path)
  exp.mx <- readRDS(exp_path)
  cyclegene <- read.table(cc.gene_path)
  exp.mx <- exp.mx[which(!(rownames(exp.mx) %in% as.vector(cyclegene$V1))),] 
  # get oneCell node pair
  oneCellNode <- tree$tip.label[grepl("_1$", tree$tip.label)]
  oneCellNode.pair <- combn(oneCellNode, 2)
  oneCellNode.pair <- data.frame(node1=oneCellNode.pair[1,], node2=oneCellNode.pair[2,], stringsAsFactors = F)
  
  # get tree tibble
  inSubtrees <- subtrees(tree)
  distMatrix <- dist.nodes(tree)
  treeTibble <- as_tibble(tree)
  
  getDist <- function(nodeName, node){
    n2 <- dplyr::filter(treeTibble, label == nodeName) %>% `[`("node") %>% unlist()
    dis <- distMatrix[node, n2]
    return(dis)
  }
  ## get all subtree depth
  all.subtree.depth <- lapply(seq(length(inSubtrees)), function(s){
    subtree <- inSubtrees[[s]]
    s.name <- subtree$name
    treeDeep <- Map(getDist, subtree$tip.label, s.name) %>% unlist() %>% as.vector() %>% max()
    return(data.frame(subtree = s.name,
                      subtree.depth = treeDeep,
                      stringsAsFactors = F))
  }) %>% bind_rows()
  # get node pair dist dataframe
  # choice distance method
  if (use.method == "Euclidean") {
    cal.dist <- function(e1, e2) return(sum(abs(e1-e2)^2)^0.5)
  } else if (use.method == "Pearson") {
    cal.dist <- function(e1, e2) return(1-cor.test(e1, e2, method = "pearson")$estimate)
  } else {
    stop("error. use.method must one of Euclidean and Pearson")
  }
  ##
  oneCellNode.pair.dist <- mclapply(seq(nrow(oneCellNode.pair)), function(i){
    node1 <- oneCellNode.pair$node1[i]
    node2 <- oneCellNode.pair$node2[i]
    n1 <- filter(treeTibble, label == node1) %>% `[`("node") %>% unlist()
    n2 <- filter(treeTibble, label == node2) %>% `[`("node") %>% unlist()
    # dist1
    mrca <- getMRCA(tree, c(node1, node2))
    treedist <- max(distMatrix[n1, mrca], distMatrix[n2, mrca])
    # dist normalize as subtree depth
    treedist2 <- filter(all.subtree.depth, subtree == mrca) %>% `[`("subtree.depth") %>% unlist()
    # calculate exp distance
    n1.exp <- exp.mx[,node1] %>% unlist() %>% as.numeric()
    n2.exp <- exp.mx[,node2] %>% unlist() %>% as.numeric()
    exp.dist <- cal.dist(n1.exp, n2.exp)
    # return result
    return(data.frame(node1=node1,
                      node2=node2,
                      mrca=mrca,
                      treeDist=treedist,
                      treeDist.asDepth = treedist2,
                      exp.dist = exp.dist,
                      stringsAsFactors = F))
  }, mc.cores = use.cores) %>% bind_rows(.)
  rownames(oneCellNode.pair.dist) <- seq(1, nrow(oneCellNode.pair.dist))
  ### run_random exp
  ran <- lapply(1:1000,function(r){
    ran1 <- exp.mx[,sample(1:ncol(exp.mx),1,replace=T)] %>% as.numeric()
    ran2 <- exp.mx[,sample(1:ncol(exp.mx),1,replace=T)] %>% as.numeric()
    ran.dist <- cal.dist(ran1, ran2)
    return(data.frame(ran.exp.dist=ran.dist,
                      ran.time=r))
  }) %>% bind_rows()
  # ran.mean <- ran %>% dplyr::summarize(ran=mean(ran.exp.dist))
  # ranexp.dist <- ran$ran %>% unlist()
  #oneCellNode.pair.dist$ran.exp.dist <- ranexp.dist
  return(list(real.dist = oneCellNode.pair.dist,
              ran.dist = ran))
}


treeDepth.dist.no.cc.besOnFilter <- function(oneCellNode_pair_dist_path, exp_path, 
                                             cc.gene_path = "/mnt/data/home/lzz/project/2020-6-18-IVDD_scRNA/material/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt", 
                                             use.method = "Euclidean",
                                             use.cores=40){
  ## run
  exp.mx <- readRDS(exp_path)
  cyclegene <- read.table(cc.gene_path)
  exp.mx <- exp.mx[which(!(rownames(exp.mx) %in% as.vector(cyclegene$V1))),] 
  # get oneCell node pair
  oneCellNode <- tree$tip.label[grepl("_1$", tree$tip.label)]
  oneCellNode.pair <- combn(oneCellNode, 2)
  oneCellNode.pair <- data.frame(node1=oneCellNode.pair[1,], node2=oneCellNode.pair[2,], stringsAsFactors = F)
  
  # get tree tibble
  inSubtrees <- subtrees(tree)
  distMatrix <- dist.nodes(tree)
  treeTibble <- as_tibble(tree)
  
  getDist <- function(nodeName, node){
    n2 <- dplyr::filter(treeTibble, label == nodeName) %>% `[`("node") %>% unlist()
    dis <- distMatrix[node, n2]
    return(dis)
  }
  ## get all subtree depth
  all.subtree.depth <- lapply(seq(length(inSubtrees)), function(s){
    subtree <- inSubtrees[[s]]
    s.name <- subtree$name
    treeDeep <- Map(getDist, subtree$tip.label, s.name) %>% unlist() %>% as.vector() %>% max()
    return(data.frame(subtree = s.name,
                      subtree.depth = treeDeep,
                      stringsAsFactors = F))
  }) %>% bind_rows()
  # get node pair dist dataframe
  # choice distance method
  if (use.method == "Euclidean") {
    cal.dist <- function(e1, e2) return(sum(abs(e1-e2)^2)^0.5)
  } else if (use.method == "Pearson") {
    cal.dist <- function(e1, e2) return(1-cor.test(e1, e2, method = "pearson")$estimate)
  } else {
    stop("error. use.method must one of Euclidean and Pearson")
  }
  ##
  oneCellNode.pair.dist <- mclapply(seq(nrow(oneCellNode.pair)), function(i){
    node1 <- oneCellNode.pair$node1[i]
    node2 <- oneCellNode.pair$node2[i]
    n1 <- filter(treeTibble, label == node1) %>% `[`("node") %>% unlist()
    n2 <- filter(treeTibble, label == node2) %>% `[`("node") %>% unlist()
    # dist1
    mrca <- getMRCA(tree, c(node1, node2))
    treedist <- max(distMatrix[n1, mrca], distMatrix[n2, mrca])
    # dist normalize as subtree depth
    treedist2 <- filter(all.subtree.depth, subtree == mrca) %>% `[`("subtree.depth") %>% unlist()
    # calculate exp distance
    n1.exp <- exp.mx[,node1] %>% unlist() %>% as.numeric()
    n2.exp <- exp.mx[,node2] %>% unlist() %>% as.numeric()
    exp.dist <- cal.dist(n1.exp, n2.exp)
    # return result
    return(data.frame(node1=node1,
                      node2=node2,
                      mrca=mrca,
                      treeDist=treedist,
                      treeDist.asDepth = treedist2,
                      exp.dist = exp.dist,
                      stringsAsFactors = F))
  }, mc.cores = use.cores) %>% bind_rows(.)
  rownames(oneCellNode.pair.dist) <- seq(1, nrow(oneCellNode.pair.dist))
  ### run_random exp
  ran <- lapply(1:1000,function(r){
    ran1 <- exp.mx[,sample(1:ncol(exp.mx),1,replace=T)] %>% as.numeric()
    ran2 <- exp.mx[,sample(1:ncol(exp.mx),1,replace=T)] %>% as.numeric()
    ran.dist <- cal.dist(ran1, ran2)
    return(data.frame(ran.exp.dist=ran.dist,
                      ran.time=r))
  }) %>% bind_rows()
  # ran.mean <- ran %>% dplyr::summarize(ran=mean(ran.exp.dist))
  # ranexp.dist <- ran$ran %>% unlist()
  #oneCellNode.pair.dist$ran.exp.dist <- ranexp.dist
  return(list(real.dist = oneCellNode.pair.dist,
              ran.dist = ran))
}

treeDepth.dist.no.cc.besOnFilter2 <- function(tree_path, exp_mx, 
                                              cc.gene_path = "/mnt/data/home/lzz/project/2020-6-18-IVDD_scRNA/material/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt", 
                                              use.method = "Euclidean",
                                              use.cores=40){
  ## run
  tree <- read.tree(tree_path)
  exp.mx <- exp_mx
  cyclegene <- read.table(cc.gene_path)
  exp.mx <- exp.mx[which(!(rownames(exp.mx) %in% as.vector(cyclegene$V1))),] 
  # get oneCell node pair
  oneCellNode <- tree$tip.label[grepl("_1$", tree$tip.label)]
  oneCellNode.pair <- combn(oneCellNode, 2)
  oneCellNode.pair <- data.frame(node1=oneCellNode.pair[1,], node2=oneCellNode.pair[2,], stringsAsFactors = F)
  oneCellNode.pair <- oneCellNode.pair %>% filter(., node1 %in% colnames(exp.mx) & node2 %in% colnames(exp.mx))
  # get tree tibble
  inSubtrees <- subtrees(tree)
  distMatrix <- dist.nodes(tree)
  treeTibble <- as_tibble(tree)
  
  getDist <- function(nodeName, node){
    n2 <- dplyr::filter(treeTibble, label == nodeName) %>% `[`("node") %>% unlist()
    dis <- distMatrix[node, n2]
    return(dis)
  }
  ## get all subtree depth
  all.subtree.depth <- lapply(seq(length(inSubtrees)), function(s){
    subtree <- inSubtrees[[s]]
    s.name <- subtree$name
    treeDeep <- Map(getDist, subtree$tip.label, s.name) %>% unlist() %>% as.vector() %>% max()
    return(data.frame(subtree = s.name,
                      subtree.depth = treeDeep,
                      stringsAsFactors = F))
  }) %>% bind_rows()
  # get node pair dist dataframe
  # choice distance method
  if (use.method == "Euclidean") {
    cal.dist <- function(e1, e2) return(sum(abs(e1-e2)^2)^0.5)
  } else if (use.method == "Pearson") {
    cal.dist <- function(e1, e2) return(1-cor.test(e1, e2, method = "pearson")$estimate)
  } else {
    stop("error. use.method must one of Euclidean and Pearson")
  }
  ##
  oneCellNode.pair.dist <- mclapply(seq(nrow(oneCellNode.pair)), function(i){
    node1 <- oneCellNode.pair$node1[i]
    node2 <- oneCellNode.pair$node2[i]
    n1 <- filter(treeTibble, label == node1) %>% `[`("node") %>% unlist()
    n2 <- filter(treeTibble, label == node2) %>% `[`("node") %>% unlist()
    # dist1
    mrca <- getMRCA(tree, c(node1, node2))
    treedist <- max(distMatrix[n1, mrca], distMatrix[n2, mrca])
    # dist normalize as subtree depth
    treedist2 <- filter(all.subtree.depth, subtree == mrca) %>% `[`("subtree.depth") %>% unlist()
    # calculate exp distance
    n1.exp <- exp.mx[,node1] %>% unlist() %>% as.numeric()
    n2.exp <- exp.mx[,node2] %>% unlist() %>% as.numeric()
    exp.dist <- cal.dist(n1.exp, n2.exp)
    # return result
    return(data.frame(node1=node1,
                      node2=node2,
                      mrca=mrca,
                      treeDist=treedist,
                      treeDist.asDepth = treedist2,
                      exp.dist = exp.dist,
                      stringsAsFactors = F))
  }, mc.cores = use.cores) %>% bind_rows(.)
  rownames(oneCellNode.pair.dist) <- seq(1, nrow(oneCellNode.pair.dist))
  ### run_random exp
  ran <- lapply(1:1000,function(r){
    ran1 <- exp.mx[,sample(1:ncol(exp.mx),1,replace=T)] %>% as.numeric()
    ran2 <- exp.mx[,sample(1:ncol(exp.mx),1,replace=T)] %>% as.numeric()
    ran.dist <- cal.dist(ran1, ran2)
    return(data.frame(ran.exp.dist=ran.dist,
                      ran.time=r))
  }) %>% bind_rows()
  # ran.mean <- ran %>% dplyr::summarize(ran=mean(ran.exp.dist))
  # ranexp.dist <- ran$ran %>% unlist()
  #oneCellNode.pair.dist$ran.exp.dist <- ranexp.dist
  return(list(real.dist = oneCellNode.pair.dist,
              ran.dist = ran))
}


test.fit.segment <- function(raw.depth.exp.dist.list, depth.cut = NULL, test.psi = -3){
  # test fit segment
  test.df <- raw.depth.exp.dist.list$real.dist
  ran.df <- raw.depth.exp.dist.list$ran.dist
  if (!is.null(depth.cut)){
    test.df$treeDist.asDepth[which(test.df$treeDist.asDepth > depth.cut)] <- depth.cut
  }
  test.df$nor.exp <- test.df$exp.dist / median(ran.df$ran.exp.dist)
  #test.sampleB <- test.sampleB %>% filter(., mrca != 1155)
  test.depth.mean.expdist <- 
    check.mean.dist.nor(raw.dist.df = test.df,
                        ran.dist.df = ran.df)
  ## fit segmented
  x <- test.depth.mean.expdist$treeDist.asDepth
  y <- test.depth.mean.expdist$mean.median.dist
  o <- lm(y ~ 1)
  xx <- -x
  o2 <- segmented(o, seg.Z = ~xx, psi=list(xx=test.psi))
  my.fitted <- fitted(o2)
  message("======== check pvalue of segmented =============")
  print(summary(o2))
  message("======== check fitted values =============")
  print(my.fitted)
  my.fitted.plot <- my.fitted %>% as.data.frame()
  my.fitted.plot$tree.depth <- rownames(my.fitted.plot) %>% as.numeric()
  colnames(my.fitted.plot)[1] <- "fit.val"
  test.plot <- 
    test.df %>%
    ggplot(aes(x=treeDist.asDepth, y=nor.exp,group=treeDist.asDepth))+
    geom_boxplot(color="grey",scale="width",adjust=0.5,trim=F) +
    # stat_summary(fun.data = "mean", geom = "crossbar",
    #              colour = "red", width = 0.2) +
    geom_line(data = my.fitted.plot, aes(x= tree.depth, y=fit.val,group="1")) + 
    #geom_line(aes(treeDist.asDepth,mediansim,color="sim",group="1"))+
    #scale_y_continuous(limits = c(20,50),breaks = c(20,30,40,50))+
    labs(x="Depth of each tree",y="Mean value of pairwise expression distance in each tree")
  ## return results
  return(list(segmented.fit.model = o2,
              fit.df = my.fitted.plot,
              use.fit.raw.df = test.df,
              use.fit.mean.df = test.depth.mean.expdist,
              test.plot = test.plot))
}



### ========== set some condition to filter ==================
## 1. filter and re normalized
condition.exp.mx.fromSeurat <- function(gene.ratio.cut = NULL, cells.features.cut, raw.seurat, all.gene.ratio){
  use.seurat <- subset(raw.seurat, subset = nFeature_RNA > cells.features.cut)
  # use gene set
  cut.genes <- all.gene.ratio %>% filter(., exp.ratio > gene.ratio.cut) %>% rownames()
  message("kee gene number: ", length(cut.genes))
  use.seurat <- use.seurat[cut.genes,]
  use.seurat <- NormalizeData(use.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  exp.mx <- as.matrix(use.seurat@assays$RNA@data)
  exp.mx <- exp.mx[Matrix::rowSums(exp.mx) > 0, Matrix::colSums(exp.mx) > 0]
  ##
  return(exp.mx)
}

## 2. normalized in all cells seurat, just extract values from exp matrix

condition.exp.mx.fromMatrix <- function(gene.ratio.cut = NULL, cells.features.cut, raw.seurat, raw.oneCell.exp, allele.infos, all.gene.ratio){
  use.seurat <- subset(raw.seurat, subset = nFeature_RNA > cells.features.cut)
  use.genes <- rownames(use.seurat)
  use.cells <- gsub("*_|-1", "", colnames(use.seurat))
  # use gene set
  cut.genes <- all.gene.ratio %>% filter(., exp.ratio > gene.ratio.cut) %>% rownames()
  message("kee gene number: ", length(cut.genes))
  use.genes <- use.genes[use.genes %in% cut.genes]
  ## filter
  allele.infos <- filter(allele.infos, cellNum == 1)
  use.cells.nodename <- allele.infos$NodeName[match(use.cells, allele.infos$BC)]
  use.cells.nodename <- na.omit(use.cells.nodename)
  exp.mx <- raw.oneCell.exp[use.genes, use.cells.nodename]
  exp.mx <- exp.mx[Matrix::rowSums(exp.mx) > 0, Matrix::colSums(exp.mx) > 0]
  ##
  return(exp.mx)
}


## filter scale matrix
condition.exp.scale.mx <- function(gene.ratio.cut = NULL, cells.features.cut, raw.seurat,allele.infos, all.gene.ratio){
  use.seurat <- subset(raw.seurat, subset = nFeature_RNA > cells.features.cut)
  # use gene set
  cut.genes <- all.gene.ratio %>% filter(., exp.ratio > gene.ratio.cut) %>% rownames()
  message("kee gene number: ", length(cut.genes))
  use.seurat <- use.seurat[cut.genes,]
  use.seurat <- use.seurat[Matrix::rowSums(use.seurat) != 0, Matrix::colSums(use.seurat) > 0]
  use.seurat <- NormalizeData(use.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  use.seurat <- ScaleData(use.seurat, features = rownames(use.seurat))
  scale.mx <- as.matrix(use.seurat[["RNA"]]@scale.data)
  colnames(scale.mx) <- gsub("*_|-1", "", colnames(scale.mx))
  allele.infos <- filter(allele.infos, cellNum == 1)
  colnames(scale.mx) <- allele.infos$NodeName[match(colnames(scale.mx), allele.infos$BC)]
  return(scale.mx)
}


## filter PCA matrix
condition.exp.PCA.mx <- function(gene.ratio.cut = NULL, cells.features.cut, raw.seurat, allele.infos, all.gene.ratio){
  use.seurat <- subset(raw.seurat, subset = nFeature_RNA > cells.features.cut)
  # use gene set
  cut.genes <- all.gene.ratio %>% filter(., exp.ratio > gene.ratio.cut) %>% rownames()
  message("kee gene number: ", length(cut.genes))
  use.seurat <- use.seurat[cut.genes,]
  use.seurat <- use.seurat[Matrix::rowSums(use.seurat) != 0, Matrix::colSums(use.seurat) > 0]
  use.seurat <- NormalizeData(use.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  use.seurat <- ScaleData(use.seurat, features = rownames(use.seurat))
  #     scale.mx <- as.matrix(use.seurat[["RNA"]]@scale.data)
  use.seurat <-
    RunPCA(use.seurat,
           features = rownames(use.seurat),
           verbose = F,
           npcs = 100)
  PCA.mx <- use.seurat@reductions$pca@cell.embeddings %>% t() %>% as.matrix()
  colnames(PCA.mx) <- colnames(use.seurat)
  colnames(PCA.mx) <- gsub("*_|-1", "", colnames(PCA.mx))
  allele.infos <- filter(allele.infos, cellNum == 1)
  colnames(PCA.mx) <- allele.infos$NodeName[match(colnames(PCA.mx), allele.infos$BC)]
  return(PCA.mx)
}

###########################
### =============================== define and source functions ======================
source("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/scripts/R_scripts/expression_saturation/1a.saturation_analysis_utils.R")

test.fit.pd.ran.segment <- function(raw.depth.exp.dist.list, depth.cut = NULL, test.psi = -3){
  # test fit segment
  test.df <- raw.depth.exp.dist.list$real.dist
  ran.df <- raw.depth.exp.dist.list$ran.dist
  if (!is.null(depth.cut)){
    test.df$treeDist.asDepth[which(test.df$treeDist.asDepth > depth.cut)] <- depth.cut
  }
  ## for Euclidean distance 
  median.ran.eu.dist <- median(ran.df$ran.eu.dist)
  median.ran.pd.dist <- median(ran.df$ran.pd.dist)
  re.depth.mean.dist <- 
    test.df %>% dplyr::group_by(mrca, treeDist.asDepth) %>% 
    dplyr::summarise(mrca.mean.eu.dist = mean(eu.dist),
                     mrca.mean.pd.dist = mean(pd.dist)) %>%
    dplyr::group_by(treeDist.asDepth) %>% 
    dplyr::summarise(mean.median.eu.dist = median(mrca.mean.eu.dist)/median.ran.eu.dist,
                     mean.median.pd.dist = median(mrca.mean.pd.dist)/median.ran.pd.dist)
  ## ==== fit segmented ======
  ## fit 
  # X
  x <- re.depth.mean.dist$treeDist.asDepth
  xx <- -x
  # Y
  y.eu <- re.depth.mean.dist$mean.median.eu.dist
  y.pd <- re.depth.mean.dist$mean.median.pd.dist
  ## ========== for EU ==========
  eu.fit <- tryCatch({
    o.eu <- lm(y.eu ~ 1)
    o2.eu <- segmented(o.eu, seg.Z = ~xx, psi=list(xx=test.psi))
    eu.breakPoint <- o2.eu$psi %>% as.data.frame()
    if (!is.null(eu.breakPoint)) {
      eu.breakPoint
    } else {
      NULL
    }
  }, error = function(e){
    NULL
  })
  if (is.null(eu.fit) | !exists("eu.fit") | nrow(eu.fit) == 0) eu.fit <-  data.frame(Initial = -3, Est. = -1, St.Err = 999)
  colnames(eu.fit) <- paste0("Eu.", colnames(eu.fit))
  ## for PD
  pd.fit <- tryCatch({
    o.pd <- lm(y.pd ~ 1)
    o2.pd <- segmented(o.pd, seg.Z = ~xx, psi=list(xx=test.psi))
    pd.breakPoint <- o2.pd$psi %>% as.data.frame()
    if (!is.null(pd.breakPoint)) {
      pd.breakPoint
    } else {
      NULL
    }
  }, error = function(e){
    NULL
  })
  if (is.null(pd.fit) | !exists("pd.fit") | nrow(pd.fit) == 0) pd.fit <-  data.frame(Initial = -3, Est. = -1, St.Err = 999)
  colnames(pd.fit) <- paste0("Pd.", colnames(pd.fit))
  ## return results
  fit.df <- bind_cols(eu.fit, pd.fit)
  return(fit.df)
}



## ======= run ===============
## set vars
tree.path <- "/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/adjust_filterMt10.nwk"
oneCell.exp.path <- "/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/adjust_filterMt10_comSCOREandUMI.oneCell.exp.Rds"

down_sample_dir <- "/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/sampleABC_re_saturation_results/sampleABC_exp_saturation/sampleA/down_sample_test/trim_nwk_files/trim0.3/"

raw.exp.euclidean.dist_path <- 
  "/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/sampleABC_re_saturation_results/sampleABC_exp_saturation/sampleA/sampleA_raw_exp_depth_euclidean_dist_nocc.Rds"

raw.exp.pearson.dist_path <- 
  "/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/sampleABC_re_saturation_results/sampleABC_exp_saturation/sampleA/sampleA_raw_exp_depth_pearson_dist_nocc.Rds"

out_rds_dir <- 
  "/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/sampleABC_re_saturation_results/sampleABC_exp_saturation/sampleA/down_sample_test/Rds_data/"



###### run
raw.eu.list <- readRDS(raw.exp.euclidean.dist_path)
raw.eu.df <- raw.eu.list$real.dist
rm(raw.eu.list)

raw.pd.list <- readRDS(raw.exp.pearson.dist_path)
raw.pd.df <- raw.pd.list$real.dist
rm(raw.pd.list)

## =========== test =======================
whole.tree.path = tree.path

whole.tree.obj <- read.tree(whole.tree.path)

fit.depth.cut = 6
sim.times = 1000

trim.ratio = 0.1
down_sample_dir <- "/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/sampleABC_re_saturation_results/sampleABC_exp_saturation/sampleA/down_sample_test/trim_nwk_files/"
tmp.dir = paste0(down_sample_dir, "trim", trim.ratio, "/")

## set output results path
out_dist_df_path <- paste0(out_rds_dir, trim.ratio, "_101_1000_euANDpd_dist_df.Rds")

out_fit_df_path <- paste0(out_rds_dir, trim.ratio, "_101_1000_fit_df.Rds")


# ### create down sample trees
# system.time(
#   mclapply(seq(sim.times), function(s){
#     # 1. get keep leaves 
#     keep.leaves <- sample(x = whole.tree.obj$tip.label,
#                           size = round(length(whole.tree.obj$tip.label)*(1-trim.ratio)))
#     # 2. build
#     keep.tree.path <- paste0(tmp.dir, "trim", trim.ratio, "_", s, "_tree.nwk")
#     
#     trim.cmd <- paste0("/mnt/data/home/lzz/anaconda3/bin/python ~/project/2019-8-22-PacBio.SampleA/scripts/PacBio_Analysis_scripts/trim_nodes_and_their_internal_nodes.py -t ",
#                        whole.tree.path,
#                        " -n ", paste(keep.leaves, collapse = ","),
#                        " -o ", keep.tree.path)
#     
#     system(trim.cmd)
#   }, mc.cores = 80)
# )
# 
## ======================== get keep node compare dist ===================================
message("=============== get simulation dist tables =============")
sim.dist.list <- lapply(seq(101,sim.times), function(s){
  # 1. get keep leaves 
  # keep.leaves <- sample(x = whole.tree.obj$tip.label,
  #                       size = round(length(whole.tree.obj$tip.label)*(1-trim.ratio)))
  # 2. build
  keep.tree.path <- paste0(tmp.dir, "trim", trim.ratio, "_", s, "_tree.nwk")
  
  # trim.cmd <- paste0("/mnt/data/home/lzz/anaconda3/bin/python ~/project/2019-8-22-PacBio.SampleA/scripts/PacBio_Analysis_scripts/trim_nodes_and_their_internal_nodes.py -t ",
  #                    whole.tree.path,
  #                    " -n ", paste(keep.leaves, collapse = ","),
  #                    " -o ", keep.tree.path)
  # 
  # system(trim.cmd)
  
  # 3.  calculate dist fit
  keep.tree <- read.tree(keep.tree.path)
  
  keep.eu.df <- raw.eu.df %>% filter((node1 %in% keep.tree$tip.label) & (node2 %in% keep.tree$tip.label))
  keep.pd.df <- raw.pd.df %>% filter((node1 %in% keep.tree$tip.label) & (node2 %in% keep.tree$tip.label))
  
  inSubtrees <- subtrees(keep.tree)
  distMatrix <- dist.nodes(keep.tree)
  treeTibble <- as_tibble(keep.tree)
  
  getDist <- function(nodeName, node){
    n2 <- dplyr::filter(treeTibble, label == nodeName) %>% `[`("node") %>% unlist()
    dis <- distMatrix[node, n2]
    return(dis)
  }
  
  ## get all subtree depth
  
  all.subtree.depth <- mclapply(seq(length(inSubtrees)), function(t){
    subtree <- inSubtrees[[t]]
    s.name <- subtree$name
    treeDeep <- Map(getDist, subtree$tip.label, s.name) %>% unlist() %>% as.vector() %>% max()
    return(data.frame(subtree = s.name,
                      subtree.depth = treeDeep,
                      stringsAsFactors = F))
  }, mc.cores = 80) %>% bind_rows()
  
  keep.dist.df <- 
    mclapply(seq(nrow(keep.eu.df)), function(i){
      #i = 1
      node1 <- keep.eu.df$node1[i]
      node2 <- keep.eu.df$node2[i]
      # n1 <- filter(treeTibble, label == node1) %>% `[`("node") %>% unlist()
      # n2 <- filter(treeTibble, label == node2) %>% `[`("node") %>% unlist()
      # dist1
      mrca <- getMRCA(keep.tree, c(node1, node2))
      # treedist <- max(distMatrix[n1, mrca], distMatrix[n2, mrca])
      # dist normalize as subtree depth
      treedist2 <- filter(all.subtree.depth, subtree == mrca) %>% `[`("subtree.depth") %>% unlist()
      # calculate exp distance
      # n1.exp <- exp.mx[,node1] %>% unlist() %>% as.numeric()
      # n2.exp <- exp.mx[,node2] %>% unlist() %>% as.numeric()
      eu.dist <- keep.eu.df$exp.dist[i]
      pd.dist <- keep.pd.df$exp.dist[i]
      # return result
      return(data.frame(node1=node1,
                        node2=node2,
                        mrca=mrca,
                        #treeDist=treedist,
                        treeDist.asDepth = treedist2,
                        eu.dist = eu.dist,
                        pd.dist = pd.dist,
                        stringsAsFactors = F))
    }, mc.cores = 80) %>% bind_rows(.)
  
  rownames(keep.dist.df) <- seq(1, nrow(keep.dist.df))
  return(keep.dist.df)
})

message("=============== done and save simulation tables=============")
saveRDS(sim.dist.list,
        file = out_dist_df_path)

###
message("================ get fit breakpoints results =============")
sim.re <- mclapply(seq(101,sim.times), function(sim.tmp){
  sim <- sim.tmp - 100
  keep.dist.df <- sim.dist.list[[sim]]
  ### run_random exp
  # keep.nodes <- unique(c(keep.dist.df$node1, keep.dist.df$node2))
  ran <- lapply(1:1000,function(r){
    ran.eu.dist <- sample(keep.dist.df$eu.dist, 1, replace = TRUE)
    ran.pd.dist <- sample(keep.dist.df$pd.dist, 1, replace = TRUE)
    return(data.frame(ran.eu.dist=ran.eu.dist,
                      ran.pd.dist=ran.pd.dist,
                      ran.time=r))
  }) %>% bind_rows()
  
  trim.list <- 
    list(real.dist = keep.dist.df,
         ran.dist = ran)
  
  fit.breakPoint <- 
    test.fit.pd.ran.segment(raw.depth.exp.dist.list = trim.list,
                            depth.cut = fit.depth.cut, 
                            test.psi = -3)
  
  fit.breakPoint <- fit.breakPoint %>%
    mutate(sim.time = sim.tmp,
           trim.ratio = trim.ratio)
  ## === return results ====
  return(fit.breakPoint)
}, mc.cores = 80) %>% bind_rows(.) 

message("================ done and save breakpoints results =============")
saveRDS(sim.re,
        file = out_fit_df_path)

##########2)up sample
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
    dplyr::summarize(dis1=mean(expdis)) %>%
    group_by(treedepth) %>%
    dplyr::summarize(dis2=median(dis1)/mean(aa$expdis)) %>%
    as.data.frame()
  
  #1)depth-expdis(Euclidean distance)
  treedepth <- aa$treedepth
  med <- aa$dis2
  od <- lm(med~1)
  xxd <- -treedepth
  
  t1 <- summary(segmented(od,seg.Z=~xxd,psi=list(xxd=-3)))
  
  
  data.frame(mybreakD=abs(t1$psi[1,2]),sebreakD=t1$psi[1,3],i,stringsAsFactors = F)
  
}) %>% rbind.fill()

save(testres,file = "~/project/293Tcelllineagetree/data.01.testres1.Rdata")
