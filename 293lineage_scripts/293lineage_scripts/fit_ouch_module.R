library(ouch)
library(ape)
library(lmtest)
library(dplyr)
library(ggplot2)
library(BiocGenerics)
library(reshape2)
library(tidyr)
# ======= save and load image =========
save.image("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/OU_model/fit_ou.RData")
load("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/OU_model/fit_ou.RData")

### =========================== our data process ===================
# read exp mx
oneCell.expMx <- readRDS("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/adjust_filterMt10_comSCOREandUMI.oneCell.exp.Rds")
oneCell.expMx <- t(oneCell.expMx)
oneCell.expMx <- as.data.frame(oneCell.expMx)
# get oneCell node tree
alltree <- read.tree("~/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/adjust_filterMt10.nwk")
oneCellNode <- alltree$tip.label[grepl("_1$", alltree$tip.label)]
multiCellNode <- alltree$tip.label[!grepl("_1$", alltree$tip.label)]
oneCellTree <- drop.tip(alltree, multiCellNode, trim.internal = FALSE)


# 1.trans ape phylo to ouchtree obj
oneCellTree <- ape2ouch(oneCellTree, scale = FALSE)
# 2.merge datafram and tree
oneCellTree.df <- as(oneCellTree, "data.frame")
oneCell.expMx$labels <- rownames(oneCell.expMx)
oneCellTree.df <- merge(oneCellTree.df, oneCell.expMx, by="labels", all=TRUE)
rownames(oneCellTree.df) <- oneCellTree.df$nodes
# 3.remake the ouch tree
oneCellTree <- with(oneCellTree.df, ouchtree(nodes=nodes,ancestors=ancestors,times=times,labels=labels))

# 4.test to fit module 
genes <- colnames(oneCell.expMx)[-ncol(oneCell.expMx)]
oneCellTree.df$regimes <- as.factor("global")
## run on all genes
library(parallel)
all.gene.OUvsBW <- mclapply(seq(length(genes)), function(i){
  select_gene <- genes[i]
  # fit bw module
  bm <- brown(tree=oneCellTree,data=oneCellTree.df[select_gene])
  # fit ou module
  ou <- hansen(
    tree=oneCellTree,
    data=oneCellTree.df[select_gene],
    regimes=oneCellTree.df["regimes"],
    sqrt.alpha=1,
    sigma=1
  )
  # 5.calculate the likelihood ratio test pvalue
  LR.value <- 2*(ou@loglik - bm@loglik)
  LR.p.value <- 1 - pchisq(2*(ou@loglik - bm@loglik), 1)
  # 6.extract coefficient
  theta <- coef(ou)$theta %>% unlist() %>% as.vector()
  alpha <- coef(ou)$alpha.matrix %>% unlist() %>% as.vector()
  sigma <- coef(ou)$sigma %>% unlist() %>% as.vector()
  # 7.return values
  return(data.frame(gene = select_gene,
                    LR.value = LR.value,
                    LR.p.value = LR.p.value,
                    theta = theta,
                    alpha = alpha,
                    sigma = sigma, stringsAsFactors = F))
}, mc.cores = 80) %>% bind_rows()

LR.p.adjust <- p.adjust(all.gene.OUvsBW$LR.p.value, method = "BH")
all.gene.OUvsBW$LR.p.adjust <- LR.p.adjust
filter(all.gene.OUvsBW, LR.p.adjust <= 0.05) %>% nrow()


## load  slope genes
slope.492.genes <- readLines("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/material/11.txt")
slope.492.genes.OUvsBW <- all.gene.OUvsBW %>% filter(gene %in% slope.492.genes)

### ====== use alpha cutoff to get.FC enrich =====

test.enrich <- function(our.all.genes, our.genes, other.genes){
  all.gene.num <- c(our.all.genes, other.genes) %>% unique() %>% length()
  TP <- sum(other.genes %in% our.genes)
  FP <- sum(!(other.genes %in% our.genes))
  FN <- length(our.genes) - TP
  TN <- all.gene.num - TP - FN - FP
  phyper.re <- phyper(TP-1, TP+FN, FP+TN, TP+FP, lower.tail = F, log.p = F)
  fc <- (TP*all.gene.num) / ((TP+FP)*(TP+FN))
  chisq.re <- chisq.test(matrix(c(TP, FN, FP, TN), nrow = 2, ncol = 2))
  return(data.frame(phyper.p = phyper.re,
                    fc = fc,
                    chisq.p = chisq.re$p.value))
}

### ============== select gene exp in cells cutoff =================
enrich.result <- function(test.gene, all.df, exp.cutoff=NULL, alpha.cutoffs=NULL, cores=20){
  if (exp.cutoff) {
    all.sig.df <- all.df %>% filter(LR.p.adjust < 0.05 & exp.percent >= exp.cutoff)
  } else {
    all.sig.df <- all.df %>% filter(LR.p.adjust < 0.05)
  }
  if (is.null(alpha.cutoffs)) alpha.cutoffs <- unique(all.sig.df$alpha)[unique(all.sig.df$alpha) %>% order()]
  #alpha.cutoffs <- seq(0,30,0.1)
  testVSall <- mclapply(seq(length(alpha.cutoffs)), function(i){
    alpha.cut <- alpha.cutoffs[i]
    sig.gene <- filter(all.sig.df, alpha >= alpha.cut) %>% `[`("gene") %>% unlist()
    # all 
    enrich.result <- test.enrich(our.all.genes = unlist(all.df$gene),
                                 our.genes = sig.gene,
                                 other.genes = test.gene)
    #
    return(data.frame(alpha.cut = alpha.cut,
                      fc.enrich = enrich.result$fc,
                      phyper.p = enrich.result$phyper.p,
                      chisq.p = enrich.result$chisq.p,
                      all.sig.num = length(sig.gene),
                      stringsAsFactors = F))
  }, mc.cores = cores) %>% bind_rows()
  return(testVSall)
}


alphaVSexpratio.enrich <- function(test.gene, all.df, exp.cut.vector=NULL, alpha.cut.vector=NULL, cores=20){
  alphaVSexp.cut.df <- lapply(seq(length(exp.cut.vector)), function(x){
    exp.ratio <- exp.cut.vector[x]
    sub.df <- enrich.result(test.gene = test.gene,
                            all.df = all.df,
                            exp.cutoff = exp.ratio,
                            alpha.cutoffs = alpha.cut.vector,
                            cores = cores)
    sub.df$exp.ratio <- exp.ratio
    return(sub.df)
  }) %>% bind_rows()
  return(alphaVSexp.cut.df)
}

# load exp matrix use to exp ratio cutoff
oneCell.exp <- readRDS("~/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/comprehen_SCOREandUMI_results_with_filter/adjust_filterMt10_comSCOREandUMI.oneCell.exp.Rds")
slope.492.exp <- oneCell.exp[slope.492.genes, ]
slope.492.expPercent <- apply(slope.492.exp, 1, function(x) return(sum(x != 0)/length(x))) %>% as.data.frame()
colnames(slope.492.expPercent) <- "ratio"
slope.492.expPercent %>% ggplot() + geom_histogram(aes(x=ratio, y=..count../sum(..count..)))


###############
##
## ==== heatmap plot for alpha cut and exp ratio cut =============
##
#################
use.alpht.cut <- seq(2,10,1)
use.exp.cut <- seq(0.1,0.9,0.1)

s492VSall.aplhaVSexpratio.df <- alphaVSexpratio.enrich(test.gene = slope.492.genes,
                                                       all.df = all.gene.OUvsBW,
                                                       exp.cut.vector = use.exp.cut,
                                                       alpha.cut.vector = use.alpht.cut,
                                                       cores = 20)
filter(s492VSall.aplhaVSexpratio.df, phyper.p > 0.05)
filter(s492VSall.aplhaVSexpratio.df, fc.enrich < 1)

rownames(s492VSall.aplhaVSexpratio.df) <- seq(1, nrow(s492VSall.aplhaVSexpratio.df))

max(plot.fd.p.df)
source("/mnt/data/home/chenfeng/Rfunction/style.print.R")
pdf("~/project/2019-8-22-PacBio.SampleA/paper_writing/plots_And_materials/OU_model_heatmap.pdf")
s492VSall.aplhaVSexpratio.df %>% 
  ggplot(aes(x=alpha.cut,fill= fc.enrich,y= exp.ratio,color=phyper.p < 0.05)) +
  geom_tile(size=0) +
  scale_fill_distiller(type="seq", direction = 1, name="Fold Enrichment") +
  geom_rect(xmin = 9.525, xmax = 10.475, ymin=0.0525, ymax = 0.1475, fill=NA, color="gray", size=1) +
  scale_color_manual(values=c("TRUE" = NA, "FALSE"="gray"), guide=NULL) +
  scale_x_continuous(breaks = use.alpht.cut, labels = paste0(">", use.alpht.cut), expand = c(0,0)) +
  scale_y_continuous(breaks = use.exp.cut, labels = paste0(use.exp.cut * 100, "%"), expand = c(0,0)) +
  labs(x = "α value in fitted Ornstein Uhlenbeck models",
       y = "Minimum fraction of expressed cells") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        legend.position = "top"
        ) 

dev.off()


