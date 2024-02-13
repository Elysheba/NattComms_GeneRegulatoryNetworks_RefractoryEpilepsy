library(beeswarm)
library(BiocParallel)
library(DT)
library(edgeR)
library(ggplot2)
library(here)
library(limma)
library(MASS)
library(openxlsx)
library(plotly)
library(RColorBrewer)
library(readxl)
library(scales)
library(tibble)
library(tximport)
library(reshape2)
library(WGCNA)
library(CoReMo)
library(lme4)
library(ggpubr)
library(patchwork)
library(ggrepel)
library(tidyverse)
library(CoReMo)
library(limma)
library(BiocParallel)
library(biomaRt)
library(GraphAT)
library(DT)
library(RColorBrewer)
library(UpSetR)
library(cluster)
library(gplots)
library(dplyr)
library(BED)
library(RTBKM)

load(here("Data/04_QC/Preprocessing.rda"))
load(here("Data/04_QC/HS_expression_matrix_batch_corrected.rda"))
load(here("Data/06_Modules/HS/Module-Data/spearman-KeptProbes-HS.rda"))

sampleInfo <- sampleInfo %>%
  filter(PA_Diagnosis %in% c("HS")) %>%
  mutate(log10_seizure_freq = log10(Seizure_Frequency_months))
d <- bc_vg[colnames(corMat),sampleInfo$GS_RNASEQ_ID]

##
wgcna.clMeth <- "ward.D"
wgcna.corMeth <- "spearman"
wgcna.softPower <- 6
wgcna.minClusSize <- 1
wgcna.k <- 22

if(file.exists(here("Data/06_Modules/HS/Module-Data/LOO_HS_Cluster.rda"))){
  load(here("Data/06_Modules/HS/Module-Data/LOO_HS_Cluster.rda"))
}else{
  set.seed(73)
  LOO.clusters <- tibble(ensgene = rownames(corMat))
  for(i in 0:ncol(d)){
    LOO.name <- paste0("BS_",i)
    print(LOO.name)
    if(i == 0) {
      d.in <- d[rownames(corMat),]
      wgcna.cormat.tmp <- corMat
    } else {
      d.in <- d[rownames(corMat),-i]
      wgcna.cormat.tmp <- d.in %>% t %>% cor(method = wgcna.corMeth)
    }
    LOO.tmp <- {1-abs(wgcna.cormat.tmp)^wgcna.softPower} %>%
      as.dist %>% stats::hclust(d = ., method = wgcna.clMeth)
    LOO.tmp <- LOO.tmp %>%
      msCutree(tree    = .,
               k       = wgcna.k,
               minsize = wgcna.minClusSize,
               d       = d,
               corMeth = wgcna.corMeth) %>%
      tibble(A = names(.),
             B = .) %>%
      `colnames<-`(c("ensgene", LOO.name))
    LOO.clusters <- LOO.clusters %>%
      left_join(LOO.tmp)
    rm(LOO.tmp, LOO.name, d.in, wgcna.cormat.tmp); gc()
  }
  LOO.clusters %>% head
  save(LOO.clusters, file = here("Data/06_Modules/HS/Module-Data/LOO_HS_Cluster.rda"))
}


if(file.exists(here("Data/06_Modules/HS/Module-Data/AdjacencyMatrix_BS_HS.rda"))){
  load(here("Data/06_Modules/HS/Module-Data/AdjacencyMatrix_BS_HS.rda"))
}else{
  # Take Clusters and create adjacency matrix
  # Sum adjacency matricies to look at all permutations
  adj.matrix <- GraphAT::clust2Mat({LOO.clusters[[2]] %>%
      `names<-`(LOO.clusters$ensgene)}) %>%
    `colnames<-`(LOO.clusters$ensgene) %>%
    `rownames<-`(LOO.clusters$ensgene)
  for(i in 3:ncol(LOO.clusters)) {
    print(paste0("Computing ", i))
    adj.matrix <- adj.matrix + {clust2Mat({LOO.clusters[[i]] %>%
        `names<-`(LOO.clusters$ensgene)}) %>%
        `colnames<-`(LOO.clusters$ensgene) %>%
        `rownames<-`(LOO.clusters$ensgene)}
  }
  save(adj.matrix, file = here("Data/06_Modules/HS/Module-Data/AdjacencyMatrix_BS_HS.rda"))
}
#
#
# Look for dense modules >900 members and consider
# them junk modules. Consider genes "junk", if the
# R2 of the module is below 0.05 and they are
# assigned to that module in >50% of permutations.
junk.genes      <- c()
mod.decision    <- c()
for(i in 2:ncol(LOO.clusters)) {
  perm.in         <- colnames(LOO.clusters)[i]
  message(paste0("Bootstrap: ", perm.in))

  for(j in 1:wgcna.k) {
    message(paste0("Module: ", j, " of ",wgcna.k))

    # get genes in module
    subset.genes  <- LOO.clusters[,c(1,i)] %>%
      filter_at(2, all_vars(. == j)) %>%
      pull(ensgene)

    # subset correlation matrix
    cormat.subset <- corMat[subset.genes,subset.genes] %>%
      .[upper.tri(.)] %>% abs %>%
      {. ^ 2}
    # get mean
    # get median
    # get size
    # append to mod.decision
    mod.decision  <- bind_rows(mod.decision,
                               tibble(Bootstrap = perm.in,
                                      `Module ID` = j,
                                      `Mean R2` = mean(cormat.subset),
                                      `Median R2` = median(cormat.subset),
                                      `Module Size` = length(subset.genes)))
  }

  junk.modules    <- mod.decision %>%
    filter(`Median R2` <= 0.05,
           Bootstrap   == perm.in) %>%
    pull(`Module ID`)

  junk.genes.tmp  <- LOO.clusters %>% filter((!!as.name(perm.in)) %in% junk.modules) %>% pull(ensgene)
  junk.genes      <- c(junk.genes, junk.genes.tmp)
}
junk.genes.thres  <- (ncol(LOO.clusters) - 1) * .5
junk.genes.out    <- junk.genes %>% table %>% {.[. >= junk.genes.thres]} %>% names
#
#
# Remove Junk genes from summarised adjancency
# Matrix

similarity.mat  <- adj.matrix[-match(junk.genes.out, rownames(adj.matrix)),
                              -match(junk.genes.out, colnames(adj.matrix))]
#
save(similarity.mat, mod.decision, junk.genes.out, junk.genes,junk.genes.thres,
     file = here("Data/06_Modules/HS/Module-Data/Bootstrapping_Similarity_HS.rda"))

