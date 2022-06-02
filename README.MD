# MOSD

## Multi-omics integration with weighted affinity and self-diffusion  applied  for cancer subtyping
### Overview
Cancer subtyping is critical for precise patients prognosis and treatment.Integrative analysis of genomic data has become an important approach for cancer subtyping.Here,we present multi-omics integration via weighted affinity and self-diffusion(MOSD) to dissect cancer heterogeneity.MOSD first construct local scaling affinity on every datatype and integrate the affinities by weighted linear combination.followed by the self-diffusion  to further improve the samples similarities for the downstream clustering analysis.
## Installation

You can install MOSD from GitHub directly using devtools.
- install.packages("devtools")
- devtools::install_github("DXCODEE/MOSD")

## Example
  ```
	 x<-readRDS("data/COADREAD.rds")
	 tes<-MOSD(x)
 ```
## compare subtypes outcome

  ```
	  library(survcomp)
          library(survival)
          library(survminer)
	  library(ggplots)
	  clin <- Surv(as.numeric(survival$Survival), as.numeric(survival$Death))
	  plot.KM(clin,as.integer(tes$clu))

 ```

