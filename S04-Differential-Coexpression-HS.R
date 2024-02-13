## Libraries
library(here)
library(dplyr)
library(parallel)

## Load files
load(here("Data/04_QC/Preprocessing.rda"))
sampleInfo <- sampleInfo %>%
  filter(PA_Diagnosis %in% c("Control", "HS") &
           Area_of_Resection %in% c("Hippocampus", "Temporal")) %>%
  filter(Library_Prep_batch != "A") %>%
  mutate(log10_seizure_freq = log10(Seizure_Frequency_months))
all_samples <- sampleInfo
load(here("Data/06_Modules/HS/Module-Data/HS_CorrMembership.rda"))
load(here("Data/04_QC/HS_expression_matrix_batch_corrected.rda"))
load(here("Data","06_Modules","HS/Module-Data/Modules-HS.rda"))

## Parameters
corMeth <- "spearman"
indTable <- all_samples %>% select(GS_RNASEQ_ID, PA_Diagnosis)
genesKept <- rownames(data.cormat)
pm <- 10000


## Load clusters
scope <- c("indTable", "bc_vg", "genesKept", "moduleList")
cl <- makeCluster(5)

## Run analyses
clusterExport(cl, scope)
rDC <- do.call(cbind, parLapply(cl, 1:pm, function(i) {
  randInd <- indTable
  randInd$PA_Diagnosis <- sample(randInd$PA_Diagnosis)
  rdiv0 <- cor(t(bc_vg[genesKept,
                       randInd$GS_RNASEQ_ID[randInd$PA_Diagnosis == "Control"]]))
  rdiv21 <- cor(t(bc_vg[genesKept,
                        randInd$GS_RNASEQ_ID[randInd$PA_Diagnosis == "HS"]]))
  toAdd <- unlist(lapply(
    moduleList,
    function(x){median(as.dist(rdiv21[unlist(x),unlist(x)])^2) -
        median(as.dist(rdiv0[unlist(x),unlist(x)])^2)
    }))
  rm(rdiv21, rdiv0)
  return(toAdd)
}))
stopCluster(cl)
save(rDC, file = here("Data/06_Modules/HS/Conservation/DifferentialCoexpression_HS.rda"))
