library(argparse)
library(dplyr)
library(parallel)

######
################## some functions####
checkCon2 <- function(Node1, Node2, alleleInfo, rawEvents){
  con1 <- filter(alleleInfo, NodeName == Node1)
  con2 <- filter(alleleInfo, NodeName == Node2)
  # get cell BC
  BC1 <- con1$BC[1]
  BC2 <- con2$BC[1]
  # get raw events
  raw1 <- filter(rawEvents, BC == BC1)[,c(paste0("target", seq(1,13)))]
  raw2 <- filter(rawEvents, BC == BC2)[,c(paste0("target", seq(1,13)))]
  # overlap edits
  overlap <- dplyr::intersect(raw1, raw2)
  return(data.frame(node1 = Node1,
                    node2 = Node2,
                    overlapE = nrow(overlap), stringsAsFactors = F))
  
  
}

getOverlapInfo.dist <- function(Pair.df, alleleInfo, rawEvents){
  OverlapInfo <- mclapply(seq(nrow(Pair.df)), function(i){
    n1 <- Pair.df$node1[i]
    n2 <- Pair.df$node2[i]
    treeD <- Pair.df$treeDist[i]
    subCheck <- checkCon2(n1, n2, alleleInfo, rawEvents)
    subCheck$treeDist <- treeD
    return(subCheck)
  }, mc.cores = 50) %>% bind_rows(.)
  return(OverlapInfo)
}

## accept argument from command line
parser <- ArgumentParser(description = "Analysis data, get the overlap infos and other results.")
parser$add_argument("--allele", type="character", help="PATH of alleleinfo file")
parser$add_argument("--moreInfo", type="character", help="PATH of moreInfo file")
parser$add_argument("--pairDist", type="character", help="PATH of oneCell Node pair tree dist file")
parser$add_argument("-o", "--out", type="character", help="PATH of output results Rds out")

useArgs <- parser$parse_args()

###### run function
run_get_results <- function(alleleInfos_path, moreInfos, oneCellNode.pair.dist_path, checkString){
  alleleInfos <- read.csv(alleleInfos_path, sep = '\t', stringsAsFactors = F)
  oneCellNode.pair.dist <- read.csv(oneCellNode.pair.dist_path, sep = '\t', stringsAsFactors = F)
  # check have anc cells and oneCell node
  haveANC.cells <- moreInfos %>% filter(., (!!sym(checkString)) == "anc") %>% `[`("BC") %>% unlist() %>% unique() # 
  haveANC.oneCellNode <- alleleInfos %>% filter(., BC %in% haveANC.cells & cellNum == 1) %>% `[`("NodeName") %>% unlist() # 
  # check overlap infos
  moreInfos.no.dis <- moreInfos %>% filter(., (!!sym(checkString)) != "dis")
  overlapInfo.dist <- getOverlapInfo.dist(oneCellNode.pair.dist, alleleInfos, moreInfos.no.dis)
  # check share_ratio
  if (length(haveANC.oneCellNode) == 1){
    haveANC.share_ratio <- NULL
  } else {
    haveANC.pair.dist <- oneCellNode.pair.dist %>% filter(., node1 %in% haveANC.oneCellNode & node2 %in% haveANC.oneCellNode)
    haveANC.overlapInfo <- getOverlapInfo.dist(haveANC.pair.dist, alleleInfos, moreInfos.no.dis)
    haveANC.share_ratio <- haveANC.overlapInfo %>% dplyr::group_by(treeDist) %>% dplyr::summarise(share_ratio=sum(overlapE > 0) / n())
  }
  return(list(haveANC_cells=haveANC.cells,
         haveANC_oneCellNode=haveANC.oneCellNode,
         Cells=nrow(alleleInfos),
         Nodes=length(unique(alleleInfos$NodeName)),
         oneCellNodes=nrow(filter(alleleInfos, cellNum == 1)),
         overlapInfo_dist=overlapInfo.dist,
         haveANC_share_ratio=haveANC.share_ratio))
}

### run
moreInfo <- read.csv(useArgs$moreInfo, sep = '\t', stringsAsFactors = F)

result_list <- run_get_results(useArgs$allele, moreInfo, useArgs$pairDist, "checkInfo")


saveRDS(result_list, useArgs$out)








