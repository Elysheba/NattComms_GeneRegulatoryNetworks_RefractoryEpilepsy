## Libraries
library(here)
library(dplyr)
library(parallel)
library(here)
library(tidyverse)
library(DT)

pname <- "AMC"
cname <- "TSC"

## get data
samples <- listSamples(projects = pname , condition = c("TSC_control")) %>%
  mutate(cohort = "Control") %>%
  bind_rows(listSamples(projects = pname , condition = c("TSC_disease")) %>%
              mutate(cohort = "Disease")) %>%
  select(sample_id = name, cohort)
contrasts <- "Disease-Control"
##
modAname <- "Correlated genes in TSC"
rpname <- "AMC"
refscope <- list(be = "Gene",
                 source = "Ens_gene",
                 organism = "human")
moduleList <- do.call(focusOnScope, c(
  list(
    x=getBekBeids(
      analysisNames=modAname,
      analysisProjects=rpname
    ),
    force=FALSE
  ),
  refscope
))
names(moduleList) <- metadata(moduleList)$name
##
varGenes <- getBekBeids(analysisNames = "Variable genes in TSC", )
varGenes <- focusOnScope(varGenes, scope = refscope, force = FALSE)
##
expression <- getDataset(project = pname,
                   name = "Batch corrected RNASEQ in TSC",
                   features = unlist(varGenes))
samples <- samples %>%
  filter(sample_id %in% colnames(expression))

diffco <- dcTest(sampleTable = samples, contrasts = contrasts, moduleList = moduleList,
                 expression = expression, corMeth = "spearman", perm = 1000, mc.cores = 5)
save(diffco, file = here("Data/06_Modules/TSC/Conservation/DifferentialCoexpression_TSC.rda"))

