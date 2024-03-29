# MOSD(Multi-omics integration with weighted affinity and self-diffusion  applied  for cancer subtyping)
## MOSD:a R package for multi-omics integration 

Depends: R(>3.6)

Import packages: survival, survminer, ggplot2

### Overview
Characterizing cancer molecular subtyping  is crucical for improving  prognosis and individualized treatment.Integrative analysis of genomic data has become an important approach for cancer subtyping.Here,we present multi-omics integration via weighted affinity and self-diffusion(MOSD) to dissect cancer heterogeneity.MOSD first construct local scaling affinity on every data type,and then integrate the affinities by weighted linear combination followed by the self-diffusion  to further improve the samples similarities for the downstream clustering analysis.
## Installation

You can install MOSD from GitHub directly using devtools.
- install.packages("devtools")
- devtools::install_github("DXCODEE/MOSD")

## Example
  ```    
  Since mRNA is too large, we used BRCA DNA methylation and miRNA to demonstrate the MOSD 
  data(BRCA)
  d<-list(x[[1]],x[[2]])
  test<-MOSD(d,6)
 ```
## Compare subtypes outcome

  ```
  library(survival)
  library(survminer)
  library(ggplot2)
  brca.clin<-x[[3]]
  clin <- Surv(as.numeric(brca.clin$futime), as.numeric(brca.clin$fustat))
  plot.KM(clin,as.integer(test$clu))

 ```

