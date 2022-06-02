# LSSD

## Self-diffusion on local scaling affinity for single cell clustering
### Overview
Single-cell RNA sequencing(scRNA-seq) has enabled gene transcriptomic profiling to be studied at individual cell level, advancing our knowledge of the cellular heterogeneity and underlying mechanisms. A key step of this methodology is to perform unsupervised clustering, which, however, often suffers challenges of zero count observations, high level of noise, as well as high dimensionality. To overcome the challenges, we proposed self-diffusion on local scaling affinity (LSSD) to enhance cells similarities metric learning for dissecting cellular heterogeneity. Local scaling infers the self-tuning of cell-to-cell distances that are used to construct cells affinity. Our approach implements self-diffusion process by propagating the affinity matrices to further improve the cell similarities for the downstream clustering analysis.
## Installation

You can install LSSDsc from GitHub directly using devtools.
- install.packages("devtools")
- devtools::install_github("DuanX8/LSSD")

## Example
  ```
	 x<-readRDS("data/pollen.rds")
	 data<-logcounts(x)
	 A<-affs(t(data))
	 Sda<-self.diffusion(A,4)
	 label<-spec.clu(Sda,11)

 ```
## Visualization
We used the diffused graph incorporating UMAP dimensionality reduction method to visualize the clusters.

  ```
	  library(umap)
	  um1<-umap(Sda)
	  lab_true<-data$cell_type1
	  ft<-data.frame(um1$layout[,1],um1$layout[,2],lab_true)
	  colnames(ft)<-c("u1","u2","lab")
	  ggplot(ft, aes(ft$u1,ft$u2))+geom_point(aes(colour = factor(ft$lab)))+scale_color_manual(values=c("#1F77B4FF","#FF7F0EFF","#2CA02CFF","#D62728FF","#9467BDFF","#8C564BFF","#E377C2FF","#7F7F7FFF"))+labs(x="LSSD_UMAP1",y="LSSD_UMAP2")+theme(legend.position = "none")

 ```
 
 ![visualization](img/LSSD_umap.png)
