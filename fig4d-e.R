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

