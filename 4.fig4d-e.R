#RNA velocity
suppressMessages({
    library(velocyto.R);
    library(tidyverse);
    #library(pagoda2);
    library(parallel);
    library(matlib)
    library(ggtree)
    library(ape)
    library(plyr)
})

raw.arrow.len <- function(rvl.re){
    current.mx <- rvl.re$current
    project.mx <- rvl.re$projected
    arrow.len.powers <- (current.mx - project.mx) ^ 2
    arrow.len <- sqrt(colSums(arrow.len.powers))
    arrow.len.df <- arrow.len %>% as.data.frame()
    colnames(arrow.len.df) <- "raw.len"
    return(arrow.len.df)
}

angle2 <- function(u, v){
    return(acos(sum(u*v)/(sqrt(sum(u^2)) * sqrt(sum(v^2)))) * 180 / pi)
}


raw.arrow.angle <- function(rvl.re, pair.df){
    current.mx <- rvl.re$current
    project.mx <- rvl.re$projected
    diff.mx <- current.mx - project.mx
    all.info <- mclapply(seq(nrow(pair.df)), function(i){
        n1 <- pair.df$node1[i]
        n2 <- pair.df$node2[i]
        ## 
        n1.project <- project.mx[,n1] %>% as.vector()
        n2.project <- project.mx[,n2] %>% as.vector()
        n1.diff <- diff.mx[,n1] %>% as.vector()
        n2.diff <- diff.mx[,n2] %>% as.vector()
        pair.angle <- angle(n1.diff, n2.diff)
        pair.angle2 <- angle2(n1.diff, n2.diff)
        arrow.dist <- 10^(sum(log(abs(n1.project-n2.project) + 0.001, 10))/100)
        return(data.frame(node1 = n1,
                          node2 = n2,
                          angle = pair.angle,
                          angle2 = pair.angle2,
                          project.dist = arrow.dist,
                          stringsAsFactors = F))
    }, mc.cores=60) %>% bind_rows()
    return(all.info)
}


get_mother_exp <- function(node_pair, exp_mx, tree_df, use_genes = NULL){
    # set vars 
    if (is.null(use_genes)){
        exp <- exp_mx
    } else {
        exp <- exp_mx[which(rownames(exp_mx) %in% use_genes),]
    }
    all_depths <- node_pair$treeDist.asDepth%>% unique() %>% unlist() %>% as.numeric() %>% sort()
    # init_exp
    use_depth = 1
    internode_exp <- NULL
    current_depth_mrca <- filter(node_pair, treeDist.asDepth==use_depth) %>% `[`("mrca") %>% unlist() %>% unique()
    init_exp <-
        mclapply(mc.cores=80,seq(length(current_depth_mrca)), function(m){
            #print(current_depth_mrca[m])
            sons <- tree_df %>% filter(.,parent==current_depth_mrca[m])
            #print("done")
            # 
            sub_exp <- exp
            leaves_exp <- sub_exp[,which(c(colnames(sub_exp) %in% sons$label))]
            current_mrca_exp <- leaves_exp %>% rowMeans() %>% as.matrix()
            colnames(current_mrca_exp) <- current_depth_mrca[m]
            return(current_mrca_exp)
        }) %>% bind_cols() %>% as.matrix()
    rownames(init_exp) <- rownames(exp)
    internode_exp <- init_exp
    # other exp
    for (use_depth in all_depths[2:length(all_depths)]){
        current_depth_mrca <- filter(node_pair, treeDist.asDepth==use_depth) %>% `[`("mrca") %>% unlist() %>% unique()
        current_depth_mrca_exp <-
            mclapply(seq(length(current_depth_mrca)), function(m){
                #print(m)
                ##print(current_depth_mrca[m])
                sons <- tree_df %>% filter(.,parent==current_depth_mrca[m])
                sub_exp <- exp
                leaves_exp <- sub_exp[,which(c(colnames(sub_exp) %in% sons$label))]
                #print("done")
                #
                inter_sons <- (sons %>% dplyr::filter(label=="1"))$node %>% unlist()
                if (is_empty(inter_sons)) {
                    current_mrca_exp <- leaves_exp %>% rowMeans() %>% as.matrix()
                    colnames(current_mrca_exp) <- current_depth_mrca[m]
                } else {
                    sub_inter_exp <-internode_exp[,which(c(colnames(internode_exp) %in% inter_sons))]
                    all_sons_exp <- bind_cols(leaves_exp, sub_inter_exp) %>% as.matrix()
                    rownames(all_sons_exp) <- rownames(leaves_exp)
                    current_mrca_exp <- all_sons_exp %>% rowMeans() %>% as.matrix()
                    colnames(current_mrca_exp) <- current_depth_mrca[m]
                }
                return(current_mrca_exp)
            }, mc.cores=80) %>% bind_cols() %>% as.matrix()
        internode_exp <- bind_cols(current_depth_mrca_exp, internode_exp) %>% as.matrix()
    }
    rownames(internode_exp) <- rownames(exp)
    return(internode_exp)
}

##############

raw.arrow.len <- function(rvl.re){
    current.mx <- rvl.re$current
    project.mx <- rvl.re$projected
    arrow.len.powers <- (current.mx - project.mx) ^ 2
    arrow.len <- sqrt(colSums(arrow.len.powers))
    arrow.len.df <- arrow.len %>% as.data.frame()
    colnames(arrow.len.df) <- "raw.len"
    return(arrow.len.df)
}

angle2 <- function(u, v){
    return(acos(sum(u*v)/(sqrt(sum(u^2)) * sqrt(sum(v^2)))) * 180 / pi)
}


raw.arrow.angle <- function(rvl.re, pair.df){
    current.mx <- rvl.re$current
    project.mx <- rvl.re$projected
    diff.mx <- current.mx - project.mx
    all.info <- mclapply(seq(nrow(pair.df)), function(i){
        n1 <- pair.df$node1[i]
        n2 <- pair.df$node2[i]
        ## 
        n1.project <- project.mx[,n1] %>% as.vector()
        n2.project <- project.mx[,n2] %>% as.vector()
        n1.diff <- diff.mx[,n1] %>% as.vector()
        n2.diff <- diff.mx[,n2] %>% as.vector()
        pair.angle <- angle(n1.diff, n2.diff)
        pair.angle2 <- angle2(n1.diff, n2.diff)
        arrow.dist <- 10^(sum(log(abs(n1.project-n2.project) + 0.001, 10))/100)
        return(data.frame(node1 = n1,
                          node2 = n2,
                          angle = pair.angle,
                          angle2 = pair.angle2,
                          project.dist = arrow.dist,
                          stringsAsFactors = F))
    }, mc.cores=60) %>% bind_rows()
    return(all.info)
}

get_mother_DistToBoundary <- 
    function(inter_exp, raw_exp) {
        select.exp.mx <- exp.mx[rownames(inter_exp),]
        max.boundary <- apply(select.exp.mx, 1, max)
        
        pseudo.count <- 0.001

        mrca.dist2boundary <- 
            lapply(seq(ncol(inter_exp)), function(m){
                sub.df <- inter_exp[,m, drop=F]
                mrca <- colnames(sub.df)
                exp.vector <- sub.df[names(max.boundary),] %>% unlist()
                exp.dist <- 10^(sum(log(abs(exp.vector-max.boundary)+pseudo.count,10))/length(max.boundary))
                return(data.frame(mrca=mrca,
                                  boundary.dist = exp.dist,
                                  stringsAsFactors = F))
            }) %>% bind_rows()
        return(mrca.dist2boundary)
    }

node_pair <- readRDS("/data02/backup/02.result/lizz/project/sampleA_results/exp_and_other_materials/sample1_depth_euclidean_exp_dist_nocc.Rds")$real.dist
exp.mx <- 
 readRDS("/data02/backup/02.result/lizz/project/sampleA_results/exp_and_other_materials/sample1_oneCellexp.Rds")

selectgene <- 
    read.csv("/public/group/lizz/project/sampleA_results/1aaaa_final_fig4DE/use_genes.csv", sep=" ", stringsAsFactors = F)

selectgene <- selectgene$gene %>% unlist()

tree <- read.tree("/public/group/lizz/project/sampleA_results/sample1.union.adjust.nwk")
# exp$gene <- rownames(exp)
tree <- as_tibble(tree)
tree2 <- tree %>% as.data.frame()

### calculate velocity
get.inTree.mx.oneCell <- function(info.df, raw.mx){
  NodeName <- info.df$NodeName[match(colnames(raw.mx), info.df$BC)]
  colnames(raw.mx) <- NodeName
  re.mx <- raw.mx[, grepl("_1$", NodeName)]
  return(re.mx)
}

sample1.loom <- read.loom.matrices("/public/group/lizz/project/sampleA_results/1a-velocity/sample1.loom")
alleleInfo <- read.csv("/data02/backup/02.result/lizz/trans_tmp/sample1.AllelesInfo.txt", 
                       sep = '\t', stringsAsFactors = F)


emat <- sample1.loom$spliced
rownames(emat) <- make.unique(rownames(emat))
nmat <- sample1.loom$unspliced
rownames(nmat) <- make.unique(rownames(nmat))

# 2. change to node name and get results
emat.onecell <- get.inTree.mx.oneCell(alleleInfo, emat)
nmat.onecell <- get.inTree.mx.oneCell(alleleInfo, nmat)

emat.select1 <- emat.onecell[rownames(emat.onecell) %in% selectgene,]
rownames(emat.select1) <- make.unique(rownames(emat.select1))
nmat.select1 <- nmat.onecell[rownames(nmat.onecell) %in% selectgene,]
rownames(nmat.select1) <- make.unique(rownames(nmat.select1))

sample1.velo <- 
    gene.relative.velocity.estimates(emat=emat.select1,
                                     nmat=nmat.select1,
                                     deltaT=1,
                                     #fit.quantile = 1,
                                     kCells=1,
                                     #cell.dist=cell.dist.710
                                     min.nmat.emat.slope = -999,
                                     min.nmat.emat.correlation = -1
                                                )


select.angle <- 
    raw.arrow.angle(rvl.re = sample1.velo,
                    pair.df = node_pair)

select.angle <- left_join(select1.angle, node_pair)


select_gene_inter_exp <-
    get_mother_exp(node_pair = node_pair,
                   exp_mx = exp.mx,
                   tree_df = tree2,
                   use_genes = selectgene)

select_mrca.dist2boundary <- 
    get_mother_DistToBoundary(inter_exp = select_gene_inter_exp,
                              raw_exp = exp.mx)


select.angle$mrca <- as.character(select.angle$mrca)

select.angle$new.mrcaDist2boundary <- 
    select_mrca.dist2boundary$boundary.dist[match(select.angle$mrca, select_mrca.dist2boundary$mrca)]

select_gene.mdistVSangle <- 
    select.angle %>% dplyr::group_by(mrca, new.mrcaDist2boundary) %>% 
    dplyr::summarise(angle.mean = mean(angle),
              angle.sd = sd(angle))

a <- cor.test(select_gene.mdistVSangle$new.mrcaDist2boundary, select_gene.mdistVSangle$angle.mean, method = 'spearman')

#OU
suppressMessages({
    library(parallel)
    library(ouch)
    library(ape)
    library(lmtest)
    library(dplyr)
    library(ggplot2)
    library(BiocGenerics)
    library(reshape2)
    library(tidyr)
    library(DT)
})

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


enrich.expANDvar.rank.cut2 <- function(test.gene, all.df, exp.cutoff=NULL, 
                               cutoffs=seq(0.1,1,0.1), cutoff.column="theta",cores=20){
   all.df <- all.df %>% arrange(desc(exp.mean))
  #all.df <- all.df[order(all.df$exp.mean, decreasing = T),]
  if (!is.null(exp.cutoff)) {
    use.all.df <- all.df[1:round(nrow(all.df) * (exp.cutoff)),]
    test.gene <- test.gene[test.gene %in% use.all.df$gene]
    all.sig.df <- use.all.df %>% filter(LR.p.adjust < 0.05)
    
  } else {
    use.all.df <- all.df
    all.sig.df <- all.df %>% filter(LR.p.adjust < 0.05)
  }
  all.sig.df <- all.sig.df %>% arrange(desc(!!sym(cutoff.column)))
  #alpha.cutoffs <- seq(0,30,0.1)
  testVSall <- mclapply(seq(length(cutoffs)), function(i){
    use.cut <- cutoffs[i]
    cut.ou.sig <- all.sig.df[1:round(nrow(all.sig.df)*use.cut),]
    sig.gene <- cut.ou.sig$gene %>% unlist()
    # all 
    enrich.re <- test.enrich(our.all.genes = unlist(all.df$gene),
                                 our.genes = sig.gene,
                                 other.genes = test.gene)
    #
    return(data.frame(cutoff = use.cut,
                      cut.var = cutoff.column,
                      fc.enrich = enrich.re$fc,
                      phyper.p = enrich.re$phyper.p,
                      chisq.p = enrich.re$chisq.p,
                      all.sig.num = length(sig.gene),
                      stringsAsFactors = F))
  }, mc.cores = cores) %>% bind_rows()
  return(testVSall)
}

varsVSexpMean.enrich2 <- function(test.gene, all.df, 
                                  exp.cut.vector=seq(0.1,1,0.1), 
                                  cut.vector=seq(0.1,1,0.1), 
                                  cut.vars = "theta",
                                  cores=20){
  alphaVSexp.cut.df <- lapply(seq(length(exp.cut.vector)), function(x){
    exp.ratio <- exp.cut.vector[x]
    sub.df <- enrich.expANDvar.rank.cut2(test.gene = test.gene,
                            all.df = all.df,
                            exp.cutoff = exp.ratio,
                            cutoffs = cut.vector,
                            cutoff.column = cut.vars,
                            cores = cores)
    sub.df$exp.ratio <- exp.ratio
    return(sub.df)
  }) %>% bind_rows()
  return(alphaVSexp.cut.df)
}

all.gene.OUvsBW <-  readRDS("sample1_all_gene_OUvsBW.Rds")

selectgene <- 
    read.csv("/public/group/lizz/project/sampleA_results/1aaaa_final_fig4DE/use_genes.csv", sep=" ", stringsAsFactors = F)

selectgene <- selectgene$gene %>% unlist()


slelect1.alphaANDexp.cut.enrich2 <- 
    varsVSexpMean.enrich2(test.gene = selectgene,
                         all.df = all.gene.OUvsBW,
                         exp.cut.vector = seq(0.1,1,0.05),
                         cut.vector = seq(0.1,1,0.05),
                         cut.vars = "alpha",
                         cores = 40)

varsVSexpMean.enrich2.sim <- 
    function(test.gene, 
             all.df, 
             exp.cut.vector=seq(0.1,1,0.1), 
             cut.vector=seq(0.1,1,0.1), 
             cut.vars = "theta",
             cores=20,
             per.time = 1000
             ){
    alphaVSexp.cut.df <- lapply(seq(length(exp.cut.vector)), function(x){
        exp.ratio <- exp.cut.vector[x]
        per.df <- 
            mclapply(mc.cores=80,seq(per.time), function(per){
                per_all.df <- sample_frac(all.df,0.8, replace=F)
                per_test.gene <- test.gene[test.gene %in% per_all.df$gene]
                sub.df <- 
                    enrich.expANDvar.rank.cut2(test.gene = per_test.gene,
                                                all.df = per_all.df,
                                                exp.cutoff = exp.ratio,
                                                cutoffs = cut.vector,
                                                cutoff.column = cut.vars,
                                                cores = cores)
                sub.df$exp.ratio <- exp.ratio
                sub.df$times <- per
                return(sub.df)
            }) %>% bind_rows()
        return(per.df)
  }) %>% bind_rows()
  return(alphaVSexp.cut.df)
}

slelect1.alphaANDexp.cut.enrich2.sim <- 
    varsVSexpMean.enrich2.sim(test.gene = selectgene,
                             all.df = all.gene.OUvsBW,
                             exp.cut.vector = seq(0.1,1,0.05),
                             cut.vector = seq(0.1,1,0.05),
                             cut.vars = "alpha",
                             cores = 2,
                             per.time=1000)

##########################
##load data and plot
#fig4d
#sample1
sampA <- readRDS("~/project/293Tcelllineagetree/rev/fig4d-e/sampleA_999_angle.Rds") %>% 
  group_by(mrca,new.mrcaDist2boundary) %>% dplyr::summarize(m=mean(angle))
#sample2
samp2 <- readRDS("~/project/293Tcelllineagetree/rev/fig4d-e/CF2_999_angle.Rds") %>% 
  group_by(mrca,new.mrcaDist2boundary) %>% dplyr::summarize(m=mean(angle))

bordA <- readRDS("~/project/293Tcelllineagetree/plot/test/06.real.motherexp.Anew.bordis.Rds")
bord2 <-  readRDS("~/project/293Tcelllineagetree/plot/test/06.real.motherexp.cf2new.bordis.Rds")

sampA$newdis <- bordA$bordis[match(sampA$mrca,bordA$int)]
samp2$newdis <- bord2$bordis[match(samp2$mrca,bord2$int)]

sampA %>% as.data.frame() %>% cbind(type="Sample 1") %>%
  #%>% rbind(samp2 %>% as.data.frame() %>% cbind(type="Sample 2")) %>%
  ggplot(aes(newdis,m))+
  geom_point(alpha=0.5,color="#F8766D")+
  scale_y_continuous(limits = c(10,120),breaks = c(20,60,100))+
  facet_grid(~type)+
  style.print()
samp2 %>% as.data.frame() %>% cbind(type="Sample 2") %>%
  ggplot(aes(newdis,m))+
  geom_point(alpha=0.5,color="#00BFC4")+
  scale_y_continuous(limits = c(10,120),breaks = c(20,60,100))+
  facet_grid(~type)+
  style.print()

#fig4e
testdata1 <- readRDS("~/project/293Tcelllineagetree/rev/fig4d-e/sampleA_OUenrich_cut_by_expANDalpha_selectGene1_sim1000_raw.Rds")
testdata2 <- readRDS("~/project/293Tcelllineagetree/rev/fig4d-e/CF2_OUenrich_cut_by_expANDalpha_selectGene1_sim1000_raw.Rds")

testdata <- rbind(testdata1,testdata2) %>% dplyr::filter(cut.var=="alpha") %>% dplyr::filter(cutoff %in% c("0.5","0.2")) %>%
  group_by(cutoff,exp.ratio,sample.name,cut.var) %>% dplyr::summarize(m=mean(fc.enrich),sd=1.96*sd(fc.enrich)) %>% 
  dplyr::filter(exp.ratio %in% c((1:10)/10,"0.3","0.7"))

testdata$cutoff <- as.character(testdata$cutoff)
testdata$cutoff <- factor(testdata$cutoff,levels = c("0.5","0.2"))
testdata$exp.ratio <- factor(testdata$exp.ratio,levels = as.character((20:1)/20))
testdata$sample.name <- factor(testdata$sample.name,levels = c("sampleA","CF2"))

ggplot(data=testdata,aes(x=exp.ratio,y=m,linetype=cutoff,colour=sample.name,group=paste(cutoff,sample.name)))+
  geom_line(position=position_dodge(width=0.5))+
  scale_y_continuous(limits = c(0.9,17),breaks = c(1,4,13),trans = "log2")+
  geom_point(position=position_dodge(width=0.5),aes(fill="ptype"))+
  geom_errorbar(position=position_dodge(width=0.5),aes(ymax=m+sd,ymin=m-sd))+
  scale_shape_manual(values = c(1,2))+
  style.print()

