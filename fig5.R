##pic cor

library(dplyr)
library(parallel)


get_pic_cor <- function(gene.pair, pic.mx, cors=20){
  results_df <- mclapply(seq(nrow(gene.pair)), function(i){
    sub.df <- gene.pair[i,]
    g1 <- sub.df$gene1[1]
    g2 <- sub.df$gene2[1]
    stringDB_info <- sub.df$inString[1]
    CORUM_info <- sub.df$complexPro[1]
    # get pic
    g1.pic <- pic.mx[, g1]
    g2.pic <- pic.mx[, g2]
    fit <- lm(g1.pic ~ g2.pic -1)
    c <- summary(fit)
    p.value <- as.data.frame(c$coefficients)$`Pr(>|t|)`
    slope <- as.data.frame(c$coefficients)$Estimate
    r.squared <- c$r.squared
    adj.r.squared <- c$adj.r.squared
    sub.df <- data.frame(gene1 = g1, 
                         gene2 = g2,
                         slope = slope,
                         p_value = p.value,
                         r.squared = r.squared,
                         adj.r.squared = adj.r.squared,
                         inString = stringDB_info,
                         complexPro = CORUM_info,
                         stringsAsFactors = F)
    return(sub.df)
  }, mc.cores = cors) %>% bind_rows()
  results_df$P_adjust <- p.adjust(results_df$p_value, method = "BH")
  return(results_df)
}



# calculate whole tree node gene pair pic correlation
# allNode.pair <- readRDS("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/pic_analysis/sampleA_filterMt10_alltree_genePair.Rds")
# wholeNode.picMx <- readRDS("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/pic_analysis/sampleA_filterMt10_alltree_PICmx.Rds")
# allNode.result <- get_pic_cor(allNode.pair, wholeNode.picMx)
# saveRDS(allNode.result, "/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/pic_analysis/sampleA_filterMt10_alltree_pic_cor.Rds")
# rm(list = c("allNode.pair", "wholeNode.picMx", "allNode.result"))
# gc()
# message("finish alltree pic cor analysis")

# calculate one cell node gene pair pic correlation
result_dir = 
    "/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/new_samples_include10x_adjust_results/rawBC/sample1/"

oneCell.genePair <- readRDS(paste0(result_dir, "PIC/sample1_oneCell_genePair_moreInfo.Rds"))
oneCellNode.picMx <- readRDS(paste0(result_dir, "PIC/sample1_oneCell_picMx.Rds"))
# oneCellNode.result <- get_pic_cor(oneCellNode.pair, oneCellNode.picMx)

# gene.pic.info <- get_pic_cor(gene.pair =  oneCellNode.pair,
#                                   pic.mx = oneCellNode.picMx,
#                                   cors = 20)
# #saveRDS(select.gene.pic.info, file = "~/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/pic_analysis/20210329_select_gene_pic_info.Rds")
# saveRDS(gene.pic.info, file = paste0(result_dir, "PIC/CF2_oneCell_picCor_info.Rds"))



## for part 1
message("run subpair part1")
sub.pair <- oneCell.genePair[1:30113424,]
rownames(sub.pair) <- seq(nrow(sub.pair))
pic.cor <- 
  get_pic_cor(sub.pair, oneCellNode.picMx, cors=40)

saveRDS(part1.cor,
        file = paste0(result_dir, "PIC/newsample1_oneCell_picCor_info.Rds"))


##original cor
library(dplyr)
library(parallel)

# calculate wholeNodeExp
get_exp_cor <- function(gene.pair, exp.mx, cors=40){
  results_df <- mclapply(seq(nrow(gene.pair)), function(i){
    sub.df <- gene.pair[i,]
    g1 <- sub.df$gene1[1]
    g2 <- sub.df$gene2[1]
    stringDB_info <- sub.df$inString[1]
    CORUM_info <- sub.df$complexPro[1]
    # get pic
    g1.exp <- exp.mx[g1, ] %>% unlist() %>% as.numeric()
    g2.exp <- exp.mx[g2, ] %>% unlist() %>% as.numeric()
    # fit pic cor
    ## fit -1
    fit1 <- lm(g1.exp ~ g2.exp -1)
    c1 <- as.data.frame(summary(fit1)$coefficients)
    ## fit 
    fit2 <- lm(g1.exp ~ g2.exp)
    c2 <- as.data.frame(summary(fit2)$coefficients)
    ## spearman
    spearman.cor <- cor.test(g1.exp, g2.exp, method = "spearman")
    ## pearson
    pearson.cor <- cor.test(g1.exp, g2.exp, method = "pearson")
    # result
    sub.df <- data.frame(gene1 = g1, 
                         gene2 = g2,
                         exp.cor1 = c1$Estimate,
                         p_value1 = c1$`Pr(>|t|)`,
                         exp.cor2 = c2$Estimate,
                         p_value2 = c2$`Pr(>|t|)`,
                         spe.cor = spearman.cor$estimate,
                         spe.pvalue = spearman.cor$p.value,
                         pea.cor = pearson.cor$estimate,
                         pea.pvalue = pearson.cor$p.value,
                         inString = stringDB_info,
                         complexPro = CORUM_info,
                         stringsAsFactors = F)
    return(sub.df)
    gc()
  }, mc.cores = cors) %>% bind_rows()
  results_df$P_adjust1 <- p.adjust(results_df$p_value1, method = "BH")
  results_df$P_adjust2 <- p.adjust(results_df$p_value2, method = "BH")
  results_df$spe.P_adjust <- p.adjust(results_df$spe.pvalue, method = "BH")
  results_df$pea.P_adjust <- p.adjust(results_df$pea.pvalue, method = "BH")
  return(results_df)
} 


# one cell node
result_dir = 
    "/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/new_samples_include10x_adjust_results/"

oneCell.genePair <- readRDS("gene_pair.Rds")
oneCellNode.picMx <- readRDS("oneCell_picMx.Rds")
oneCell.exp <- readRDS("oneCellexp.Rds"))

rest.cor <- 
  get_exp_cor(oneCell.genePair, oneCell.exp, cors=40)

saveRDS(rest.cor,
        file = paste0(result_dir, "oneCell_expCor_info_rest.Rds"))

message("finish one cell tree exp cor")



