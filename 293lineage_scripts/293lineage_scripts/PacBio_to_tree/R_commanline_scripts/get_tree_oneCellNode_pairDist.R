library(dplyr)
library(parallel)
library(ape)
library(ggtree)
library(argparse)


getTreeDist <- function(pair.df, intree){
  distMatrix <- dist.nodes(intree)
  df.tree <- as_tibble(intree)
  oneCellNode.pair.dist <- mclapply(seq(nrow(pair.df)), function(i){
    node1 <- pair.df$node1[i]
    node2 <- pair.df$node2[i]
    n1 <- filter(df.tree, label == node1) %>% `[`("node") %>% unlist()
    n2 <- filter(df.tree, label == node2) %>% `[`("node") %>% unlist()
    mrca <- getMRCA(intree, c(node1, node2))
    treedist <- max(distMatrix[n1, mrca], distMatrix[n2, mrca])
    return(data.frame(node1=node1,
                      node2=node2,
                      treeDist=treedist,
                      stringsAsFactors = F))
  }, mc.cores = 20) %>% bind_rows(.)
  return(oneCellNode.pair.dist)
}


## accept argument from command line
parser <- ArgumentParser(description = "get oneCelll node pair distance in tree")
parser$add_argument("-i", "--intree", type="character", help="PATH of tree, nwk")
parser$add_argument("-o", "--outDist", type="character", help="PATH of output tree dist")

args <- parser$parse_args()

## run
tree <- read.tree(args$intree)
# get oneCell node pair
oneCellNode <- tree$tip.label[grepl("_1$", tree$tip.label)]
oneCellNode.pair <- combn(oneCellNode, 2)
oneCellNode.pair <- data.frame(node1=oneCellNode.pair[1,], node2=oneCellNode.pair[2,], stringsAsFactors = F)
# get tree dist
oneCellNode_dist <- getTreeDist(oneCellNode.pair, tree)
# write out results
write.table(oneCellNode_dist, args$outDist, sep = '\t', quote = F, row.names = F)




