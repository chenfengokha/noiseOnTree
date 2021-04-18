library(dplyr)
library(parallel)



# calculate one cell node gene pair pic correlation
oneCellNode.pair <- readRDS("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/pic_analysis/sampleA_filterMt10_oneCellNode_genePair.Rds")
oneCellNode.picMx <- readRDS("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/pic_analysis/sampleA_filterMt10_oneCell_PICmx.Rds")

allComplexes.df <- read.csv("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/material/allComplexes.txt", sep = '\t', stringsAsFactors = F)
allComplexes.df <- filter(allComplexes.df, Organism == "Human")
###

get_pic_cor.info <- function(gene.pair, pic.mx, complex.df, cores=80){
  results_df <- mclapply(seq(nrow(gene.pair)), function(i){
    g1 <- gene.pair$gene1[i]
    g2 <- gene.pair$gene2[i]
    # get pic
    g1.pic <- pic.mx[, g1]
    g2.pic <- pic.mx[, g2]
    # get complex info
    g1_g2.comp <- complex.df$ComplexName[mapply(grepl, g1, complex.df$subunits.Gene.name.) & mapply(grepl, g2, complex.df$subunits.Gene.name.)]
    if (length(g1_g2.comp) == 0) g1_g2.comp <- "N"
    # fit pic cor
    # fit pic cor
    fit <- lm(g1.pic ~ g2.pic -1)
    c <- summary(fit)
    p.value <- as.data.frame(c$coefficients)$`Pr(>|t|)`
    slope <- as.data.frame(c$coefficients)$Estimate
    r.squared <- c$r.squared
    adj.r.squared <- c$adj.r.squared
    # return sub df
    sub.df <- data.frame(gene1 = g1, 
                         gene2 = g2,
                         slope = slope,
                         p_value = p.value,
                         r.squared = r.squared,
                         adj.r.squared = adj.r.squared,
                         complexPro = g1_g2.comp,
                         stringsAsFactors = F)
    return(sub.df)
  }, mc.cores = cores) %>% bind_rows()
  results_df$P_adjust <- p.adjust(results_df$p_value, method = "BH")
  return(results_df)
}


gene.pic.info <- get_pic_cor.info(gene.pair =  oneCellNode.pair,
                                  pic.mx = oneCellNode.picMx,
                                  complex.df = allComplexes.df,
                                  cores = 60)
#saveRDS(select.gene.pic.info, file = "~/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/pic_analysis/20210329_select_gene_pic_info.Rds")
saveRDS(gene.pic.info, file = "~/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/pic_analysis/20210330_sampleA_filterMt10_pic_exp_info.Rds")



message("finish one cell node tree pic cor analysis")










