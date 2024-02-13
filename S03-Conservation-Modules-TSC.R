########################################@
## Assess conservation of modules ----

rm(list=ls())
gc()

library(CoReMo)
library(randomcoloR)
library(TBTools)
library(RTBKM)
library(here)
library(parallel)
library(tidyverse)
library(DT)

#*******************************************************************************@
## Reference = TSC modules ----
modAname <- "Correlated genes in TSC"
rpname <- "AMC"
refscope <- list(be = "Gene",
                 source = "Ens_gene",
                 organism = "human")
refList <- do.call(focusOnScope, c(
  list(
    x=getBekBeids(
      analysisNames=modAname,
      analysisProjects=rpname
    ),
    force=FALSE
  ),
  refscope
))
names(refList) <- metadata(refList)$name
varGenes <- focusOnScope(
                    getBekBeids(analysisNames = "Variable genes in TSC"), force = FALSE,
                    scope = refscope)

#*******************************************************************************@
## AMC conservation  ----
print("AMC other conditions")
an <- listAnalyses(projects = c("AMC"), types = "Quantification")
an <- getAnalyses(projects = c("AMC"), types = "Quantification", names = an$name)
an <- inner_join(an$analyses %>% dplyr::select(analysisName = name,
                                               project,
                                               uuid),
                 an$inputs$Conditions %>% dplyr::select(condition = name,
                                                        # project,
                                                        description,
                                                        uuid),
                 by = "uuid") %>%
  inner_join(an$inputs$Datasets %>% dplyr::select(dataset = name,
                                                  # project,
                                                  uuid),
             by = "uuid") %>%
  inner_join(an$features$`BE scope`,
             by = "uuid") %>%
  distinct() %>%
  filter(!condition %in% c("TSC", "HS", "HS_disease", "FCD2b_disease", "TSC_disease", "mTOR_disease", "GG",
                           "GNT_disease", "GNT_control", "FCD 2b", "Control", "Dravet", "DNT"))
## getAnalyses doesn't list the correct dataset to use for mTOR, TSC and HS
an <- an %>%
  mutate(dataset = case_when(condition == "TSC_control" ~ "Batch corrected RNASEQ in TSC",
                             condition == "mTOR_control" ~ "Batch corrected RNASEQ in mTORopathies",
                             condition == "HS_control" ~ "Batch corrected RNASEQ in hippocampal sclerosis",
                             TRUE ~ dataset))
toRet <- lapply(1:nrow(an),
                function(tc){
                  ##
                  san <- an %>% slice(tc)
                  print(paste("Project", san$project, "- Condition:", san$condition))

                  ## correlation matrix
                  samples <- listSamples(projects = san$project, condition = san$condition)

                  ## Expression matrix
                  varGenes <- focusOnScope(varGenes,
                                           be = san$be,
                                           source = san$source,
                                           organism = san$organism)

                  expr <- getDataset(project = san$project,
                                     name = san$dataset,
                                     samples = samples$name,
                                     feat = unlist(varGenes))

                  ## reference
                  # crefList <- focusOnScope(refList,
                  #                          be = san$be,
                  #                          source = san$source,
                  #                          organism = san$organism)
                  crefList <- refList
                  ## remove duplicated entries from the converted reference module list
                  dup <- unlist(crefList)[duplicated(unlist(crefList))]
                  length(dup)
                  crefList <- filterByBEID(crefList,
                                           intersect(setdiff(unlist(crefList), dup), rownames(expr)))
                  ## Only keep module genes
                  expr <- expr[unlist(crefList),]

                  ## Conservation strategy
                  cons <- lmcTest(l = crefList, d = expr, f = corSummary,
                                  method="spearman", sf = median, square = TRUE, byrow = TRUE,
                                  alternative = "greater",
                                  perm = 1000, mc.cores = 20)

                  toRet <- list(
                    conservation = cons,
                    info = san,
                    bescope = list(be = san$be,
                                   source = san$source,
                                   organism = san$organism),
                    refList = crefList)
                  rm(expr, feat, samples, cons, crossTab, crefList)
                  return(toRet)
                })

names(toRet) <- paste(an %>% pull(project),
                      an %>% pull(condition),
                      sep = "_")
print("Saving conservation AMC")

save(toRet, file = here("Data/06_Modules/TSC/Conservation/TSC_AMC_Control_Conservation.rda"))
