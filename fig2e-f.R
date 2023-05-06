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


tree.path <- "tree.nwk"
oneCell.exp.path <- "oneCell.exp.Rds"


tree <- read.tree(tree.path)
oneCell.exp.mx <- readRDS(oneCell.exp.path)

raw.depth.exp.dist <- 
  readRDS("depth_dist_nocc.Rds")

all.gene.ratio <- apply(oneCell.exp.mx, 1, function(x) sum(x > 0) / length(x)) %>% as.data.frame()


raw.depth.euclidean.exp.dist.nocc <- 
  treeDepth.dist.no.cc(tree_path = tree.path,
                       exp_path = oneCell.exp.path,
                       use.method = "Euclidean",
                       use.cores = 80)


raw.exp.fit <- 
    test.fit.segment(raw.depth.exp.dist.list = raw.depth.euclidean.exp.dist.nocc,
                     depth.cut = 6,
                     test.psi = -3)