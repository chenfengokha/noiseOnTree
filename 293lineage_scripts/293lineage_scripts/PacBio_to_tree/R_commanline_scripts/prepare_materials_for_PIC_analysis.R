library(ape)
library(ggtree)
library(castor)
library(dplyr)
library(argparse)




#========== get pic mx ===========
getPICmx <- function(exp.mx.path, intree){
  exp.mx <- readRDS(exp.mx.path)
  t.exp.mx <- t(exp.mx)
  # order by tree
  t.exp.mx <- t.exp.mx[intree$tip.label, ]
  pic.result <- get_independent_contrasts(intree, t.exp.mx)
  pic.mx <- pic.result$PICs
  return(pic.mx)
}



## accept argument from command line
parser <- ArgumentParser(description = "get oneCelll node pair distance in tree")
parser$add_argument("-i", "--intree", type="character", help="PATH of tree, nwk")
parser$add_argument("--oneCellExp", type="character", help="PATH of one cell node exp matrix")
parser$add_argument("--alltreeExp", type="character", help="PATH of all tree node exp matrix")
parser$add_argument("-o", "--outDir", type="character", help="PATH of output directory")

args <- parser$parse_args()



### test sampleA ####
alltree <- read.tree("~/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/adjust_filterMt10.nwk")
oneCellNode <- alltree$tip.label[grepl("_1$", alltree$tip.label)]
multiCellNode <- alltree$tip.label[!grepl("_1$", alltree$tip.label)]
oneCellTree <- drop.tip(alltree, multiCellNode, trim.internal = FALSE)








# one cell Node
oneCell.exp.mx <- readRDS("~/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/adjust_filterMt10_comSCOREandUMI.oneCell.exp.Rds")
t.oneCell.exp.mx <- t(oneCell.exp.mx)
t.oneCell.exp.mx <- t.oneCell.exp.mx[oneCellTree$tip.label, ]
oneCellNode.pic.result <- get_independent_contrasts(oneCellTree,t.oneCell.exp.mx)
oneCellNode.pic.mx <- oneCellNode.pic.result$PICs
# all cell Node
alltree.exp.mx <- readRDS("~/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/adjust_filterMt10_comSCOREandUMI.allTree.Exp.Rds")
t.alltree.exp.mx <- t(alltree.exp.mx)
t.alltree.exp.mx <- t.alltree.exp.mx[alltree$tip.label, ]
alltree.pic.result <- get_independent_contrasts(alltree, t.alltree.exp.mx)
alltree.pix.mx <- alltree.pic.result$PICs

#========  get gene pair ===========
get.genePair <- function(t.exp.mx, stringDB_info, CORUM_info){
  genePair <- combn(colnames(t.exp.mx), 2)
  genePair <- data.frame(gene1=genePair[1,], gene2=genePair[2,], stringsAsFactors = F)
  ### join more infos
  # join stringDB info
  new.genePair <- left_join(genePair, stringDB_info)
  new.genePair <- new.genePair[!duplicated(new.genePair[,1:2]),]
  # join CORUM info
  new.genePair <- left_join(new.genePair, CORUM_info)
  genePair.MoreInfo <- new.genePair[!duplicated(new.genePair[,1:2]),]
  genePair.MoreInfo[is.na(genePair.MoreInfo)] <- "N"
  print(nrow(genePair) == nrow(genePair.MoreInfo))
  return(genePair.MoreInfo)
}




oneCell.gene.pair <- combn(colnames(t.oneCell.exp.mx), 2)
oneCell.gene.pair <- data.frame(gene1 = oneCell.gene.pair[1,], gene2 = oneCell.gene.pair[2,], stringsAsFactors = F)

new.oneCell.gene.pair <- left_join(oneCell.gene.pair, stringDB_double)
new.oneCell.gene.pair <- new.oneCell.gene.pair[!duplicated(new.oneCell.gene.pair),]
new.oneCell.gene.pair <- left_join(new.oneCell.gene.pair, CORUM_data)
oneCell.gene.pair.MoreInfo <- new.oneCell.gene.pair[!duplicated(new.oneCell.gene.pair[,1:2]),]
oneCell.gene.pair.MoreInfo[is.na(oneCell.gene.pair.MoreInfo)] <- "N"
##

alltree.gene.pair <- combn(colnames(t.alltree.exp.mx),2)
alltree.gene.pair <- data.frame(gene1 = alltree.gene.pair[1,], gene2 = alltree.gene.pair[2,], stringsAsFactors = F)
new.alltree.gene.pair <- left_join(alltree.gene.pair, stringDB_double)
new.alltree.gene.pair <- new.alltree.gene.pair[!duplicated(new.alltree.gene.pair[,1:2]), ]
new.alltree.gene.pair <- left_join(new.alltree.gene.pair, CORUM_data)
alltree.gene.pair.MoreInfo <- new.alltree.gene.pair[!duplicated(new.alltree.gene.pair[,1:2]),]
alltree.gene.pair.MoreInfo[is.na(alltree.gene.pair.MoreInfo)] <- "N"
saveRDS(alltree.gene.pair.MoreInfo, "~/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/pic_analysis/sampleA_filterMt10_alltree_genePair.Rds")
























# load stringDB data
library(dplyr)
stringDB <- read.csv("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/material/stringDB_knownGene_corSig.txt", sep = '\t', stringsAsFactors = F)
stringDB2 <- stringDB[,c(2,1,3)]
colnames(stringDB2) <- c("gene1", "gene2", "coexpression")
stringDB_double <- bind_rows(stringDB, stringDB2)
rm(list = c("stringDB", "stringDB2"))
stringDB_double <- stringDB_double[,1:2]
stringDB_double$inString <- "Y"

# load CORUM data
CORUM_data <- read.csv("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/material/human_complex_CORUM_df.txt", sep = '\t', stringsAsFactors = F)



tmp_result <- left_join(oneCell.Exp.corSig_infos, CORUM_data)
tmp_result <- tmp_result[!duplicated(tmp_result),]










