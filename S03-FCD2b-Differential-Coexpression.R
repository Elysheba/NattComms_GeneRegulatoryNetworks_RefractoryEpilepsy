## Libraries
library(here)
library(dplyr)
library(parallel)

## Load files
load(here("Data/04_QC/Preprocessing.rda"))
sampleInfo <- sampleInfo %>%
  filter(PA_Diagnosis %in% c("Control", "FCD 2b") &
           Area_of_Resection %in% c("Frontal", "Occipital", "Parietal", "Temporal")) %>%
  mutate(log10_seizure_freq = log10(Seizure_Frequency_months))
all_samples <- sampleInfo
load(here("Data/06_Modules/FCD2b/Module-Data/FCD2b_CorrMembership.rda"))
load(here("Data","06_Modules","FCD2b/Module-Data/Modules-FCD2b.rda"))

## Parameters
corMeth <- "spearman"
indTable <- all_samples %>% select(GS_RNASEQ_ID, PA_Diagnosis)
genesKept <- rownames(data.cormat)
pm <- 10000

## Load clusters
scope <- c("indTable", "vg", "genesKept", "moduleList")
cl <- makeCluster(5)

## Run analyses
clusterExport(cl, scope)
rDC <- do.call(cbind, parLapply(cl, 1:pm, function(i) {
  randInd <- indTable
  randInd$PA_Diagnosis <- sample(randInd$PA_Diagnosis)
  rdiv0 <- cor(t(vg$E[genesKept,
                       randInd$GS_RNASEQ_ID[randInd$PA_Diagnosis == "Control"]]))
  rdiv21 <- cor(t(vg$E[genesKept,
                        randInd$GS_RNASEQ_ID[randInd$PA_Diagnosis == "FCD 2b"]]))
  toAdd <- unlist(lapply(
    moduleList,
    function(x){median(as.dist(rdiv21[unlist(x),unlist(x)])^2) -
        median(as.dist(rdiv0[unlist(x),unlist(x)])^2)
    }))
  rm(rdiv21, rdiv0)
  return(toAdd)
}))
stopCluster(cl)
save(rDC, file = here("Data/06_Modules/FCD2b/Conservation/DifferentialCoexpression_FCD2b.rda"))
