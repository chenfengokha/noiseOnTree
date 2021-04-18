# load pkgs
library(argparse)
library(dplyr)
library(parallel)

# some functions
get.exp.dist <- function(treeDist.path, exp.matrix){
  treeDist.df <- read.csv(treeDist.path, sep = '\t', stringsAsFactors = F)
  treeDistAndExpDist <- mclapply(seq(nrow(treeDist.df)), function(i){
    node1 <- treeDist.df$node1[i]
    node2 <- treeDist.df$node2[i]
    treeDist <- treeDist.df$treeDist[i]
    # get exp 
    node1.exp <- exp.matrix[, node1] %>% as.numeric()
    node2.exp <- exp.matrix[, node2] %>% as.numeric()
    expDist <- dist(rbind(node1.exp, node2.exp), method = "euclidean")[1]
    return(data.frame(node1=node1, node2=node2, treeDist = treeDist, expDist = expDist))
  }, mc.cores = 18) %>% bind_rows(.)
  return(treeDistAndExpDist)
}

get.treeDistVSexpDist <- function(treeANDexpDis.df){
  treeDistVSexpDist <- treeANDexpDis.df %>% group_by(treeDist) %>% summarise(exp.dis.mean=mean(expDist), exp.dis.sd=sd(expDist))
  return(treeDistVSexpDist)
}

## accept argument from command line
parser <- ArgumentParser(description = "calculate the correlation of tree distance and expression distance.")
parser$add_argument("--inTreeDist", type="character", help="PATH of input tree distance file")
parser$add_argument("--inExpMx", type="character", help="PATH of input exp matrix")
parser$add_argument("-o", "--out", type="character", help="PATH of output results RData")

args <- parser$parse_args()
#

## run
#test 
exp.matrix <- readRDS(args$inExpMx)
exp.tree.dist <- get.exp.dist(args$inTreeDist, exp.matrix)
treeDistVSexpDist <- get.treeDistVSexpDist(exp.tree.dist)
out_results <- list(raw=exp.tree.dist,
                    treeVSexpMean=treeDistVSexpDist)
saveRDS(out_results, args$out)





