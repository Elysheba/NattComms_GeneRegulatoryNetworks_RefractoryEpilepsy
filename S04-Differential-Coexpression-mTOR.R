## Libraries
library(here)
library(dplyr)
library(parallel)
library(here)
library(tidyverse)
library(DT)

pname <- "AMC"
cname <- "mTOR"

#################################@
## Compare cohorts
## FCD2a-FCD2b, FCD2a-TSC, FCD2b-TSC

## get data
samples <- listSamples(projects = pname , condition = c("mTOR_control")) %>%
  mutate(cohort = "Control") %>%
  bind_rows(listSamples(projects = pname , condition = c("mTOR_disease")) %>%
              mutate(cohort = "Disease")) %>%
  select(sample_id = name, cohort)
contrasts <- c("FCD2a-FCD2b", "FCD2a-TSC", "FCD2b-TSC")

##
load(here("Data/04_QC/Preprocessing.rda"))
sampleInfo <- sampleInfo %>%
  filter(PA_Diagnosis %in% c("Control", "FCD 2a", "FCD 2b", "TSC") &
           !Area_of_Resection %in% c("Hippocampus")) %>%
  filter(!GS_RNASEQ_ID %in% c("103866-001-192")) %>%
  mutate(log10_seizure_freq = log10(Seizure_Frequency_months)) %>%
  filter(Age_At_Time_of_Operation <= 46) %>%
  select(GS_RNASEQ_ID, PA_Diagnosis)
samples <- inner_join(samples, sampleInfo,
                      by = c("sample_id" = "GS_RNASEQ_ID")) %>%
  select(sample_id, cohort = PA_Diagnosis) %>%
  filter(cohort != "Control") %>%
  mutate(cohort = gsub(" ", "", cohort))

##
modAname <- "Correlated genes in mTOR"
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
varGenes <- getBekBeids(analysisNames = "Variable genes in mTORopathies", )
varGenes <- focusOnScope(varGenes, scope = refscope, force = FALSE)
##
expression <- getDataset(project = pname,
                         name = "Batch corrected RNASEQ in mTORopathies",
                         features = unlist(varGenes))
samples <- samples %>%
  filter(sample_id %in% colnames(expression))

diffco <- lapply(contrasts, function(x){
  toRet <- dcTest(sampleTable = samples, contrasts = x, moduleList = moduleList,
                  expression = expression, corMeth = "spearman", perm = 1000, mc.cores = 5)
  return(toRet)
})

save(diffco, file = here("Data/06_Modules/mTOR/Module-Data/DifferentialCoexpression_mTOR.rda"))











