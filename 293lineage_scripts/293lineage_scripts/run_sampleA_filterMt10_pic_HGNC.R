library(tidyverse)
library(parallel)

oneCell.exp <- readRDS("~/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/adjust_filterMt10_comSCOREandUMI.oneCell.exp.Rds")
oneCellNode.picMx <- readRDS("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/pic_analysis/sampleA_filterMt10_oneCell_PICmx.Rds")

allComplexes.df <- read.csv("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/material/allComplexes.txt", sep = '\t', stringsAsFactors = F)
allComplexes.df <- filter(allComplexes.df, Organism == "Human")

###  ======================   HGNC group info calculate
HGNC.group.df <- read.csv("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/material/HGUC_group.csv", stringsAsFactors = F)

use.group.df <- filter(HGNC.group.df, parent.group.ID != "#N/A")

use.gene.name <- lapply(seq(nrow(use.group.df)), function(i){
  As <- use.group.df$Approved.symbol[i]
  # previous symbols
  ps <- use.group.df$Previous.symbols[i] %>% toupper()
  # Alias symbols
  als <- use.group.df$Alias.symbols[i] %>% toupper()
  # 
  if (As %in% rownames(oneCell.exp)) {
    re.symbol <- As
  } else if (ps != "" & any(unlist(strsplit(ps, ",")) %in% rownames(oneCell.exp))) {
    ps_vector <- unlist(strsplit(ps, ","))
    re.symbol <- ps_vector[ps_vector %in% rownames(oneCell.exp)][1]
  } else if (als != "" & any(unlist(strsplit(als, ",")) %in% rownames(oneCell.exp))) {
    als_vector <- unlist(strsplit(als, ","))
    re.symbol <- als_vector[als_vector %in% rownames(oneCell.exp)][1]
  } else {
    return("N")
  }
}) %>% unlist()

use.group.df <- cbind(use.gene.name, use.group.df)


gene.combn.df <- use.group.df %>% group_by(parent.group.ID) %>% do({
  myDF <- .
  use.gene <- myDF$use.gene.name[myDF$use.gene.name != 'N']
  if (length(use.gene) > 1) {
    use.combn <- combn(use.gene, 2)
    use.combn <- data.frame(gene1=use.combn[1,], gene2=use.combn[2,], stringsAsFactors = F)
    use.combn$gene1.Group.ID <- myDF$Group.ID[match(use.combn$gene1, myDF$use.gene.name)]
    use.combn$gene2.Group.ID <- myDF$Group.ID[match(use.combn$gene2, myDF$use.gene.name)]
    use.combn$gene1.Group.name <- myDF$Group.name[match(use.combn$gene1, myDF$use.gene.name)]
    use.combn$gene2.Group.name <- myDF$Group.name[match(use.combn$gene2, myDF$use.gene.name)]
    use.combn$parent.group.ID <- myDF$parent.group.ID[1]
    use.combn
  } else {data.frame(NA)}
})

gene.combn.df <- gene.combn.df %>% select(-ncol(gene.combn.df))
gene.combn.df <- gene.combn.df[complete.cases(gene.combn.df),]

get_cor.info.HGNC <- function(gene.pair, use.mx, complex.df, type="exp", cores=80){
  #result_list <- list()
  results_df <-  mclapply(seq(nrow(gene.pair)), function(i){
    nvector <- gene.pair[i,]
    g1 <- nvector$gene1
    g2 <- nvector$gene2
    # add some info
    g1.group.id <- nvector$gene1.Group.ID
    g2.group.id <- nvector$gene2.Group.ID
    g1.group.name <- nvector$gene1.Group.name
    g2.group.name <- nvector$gene2.Group.name
    parent.group.id <- nvector$parent.group.ID
    # get exp
    if (type == "exp") {
      g1.vector <- use.mx[g1, ] %>% unlist()
      g2.vector <- use.mx[g2, ] %>% unlist()
    } else {
      g1.vector <- use.mx[, g1] %>% unlist()
      g2.vector <- use.mx[, g2] %>% unlist()
    }
    
    # get complex info
    g1_g2.comp <- complex.df$ComplexName[mapply(grepl, g1, complex.df$subunits.Gene.name.) & mapply(grepl, g2, complex.df$subunits.Gene.name.)]
    if (length(g1_g2.comp) == 0) g1_g2.comp <- "N"
    # fit exp cor
    fit <- lm(g1.vector ~ g2.vector -1)
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
                         gene1.Group.ID = g1.group.id,
                         gene2.Group.ID = g2.group.id,
                         gene1.Group.name = g1.group.name,
                         gene2.Group.name = g2.group.name,
                         parent.group.ID = parent.group.id,
                         stringsAsFactors = F)
    #result_list[[i]] <<- sub.df
  }, mc.cores = cores) %>% bind_rows()
  results_df$P_adjust <- p.adjust(results_df$p_value, method = "BH")
  #return(result_list)
  return(results_df)
}


# ## exp cor 
# select.exp.cor.info.HGNC <- get_exp_cor.info.HGNC(gene.pair = gene.combn.df,
#                                                   exp.mx = oneCell.exp,
#                                                   complex.df = allComplexes.df,
#                                                   cores = 20)



load("/mnt/data/home/chenfeng/lineagepaper/complextest.Rdata")
small.combn.df <- filter(gene.combn.df, parent.group.ID %in% aa$mother)
## pic cor
select.pic.cor.info.HGNC <- get_cor.info.HGNC(gene.pair = small.combn.df,
                                              pic.mx = oneCellNode.picMx,
                                              complex.df = allComplexes.df,
                                              type = "pic"
                                              cores = 28)




saveRDS(select.pic.cor.info.HGNC, file = "/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/pic_analysis/select_gene_picANDexp/HGNC_pic.Rds")
