
#' @name MOSD
#' @title Multi-omics integration with weighted local scaling affinity and self-diffusion  applied  for cancer subtyping
#' @description Cancer subtypes identifition by multi-omics integration
#' @param x A list  of matrices contains multi-omics datasets where row is samples and column is features
#' @param k Clustering number;if k==NULL,the clustering number is estimated by seperarion cost method.

#' @return a list;
#' (1)clu:the clustering result;
#' (2)S:the final integrated diffusion matrix;
#' (3)es.num:the estimated clustering number;
#' @export
#' @examples
#' data(COADREAD)
#' d<-list(mydatGE,mydatME,mydatMI)
#' test <- MOSD(d)

MOSD<-function(x,k=NULL){
    s<-sapply(1:length(x), function(i){
      dim(x[[i]])[2]
    })
    w<-perc(s)
    we<-exp(w)
    p<-perc(we)
    A<-lapply(1:length(x),function(i){
      affs(x[[i]])
    })
    AL<-lapply(1:length(x), function(i){
      A[[i]]*p[i]
    })
    AA<-Reduce('+',AL)

    sd_matrix<-self.diffusion(AA,3)

    if(is.null(k)){
      es<-es.num.graph(sd_matrix,2:10)
      k<-es[[1]]
    }

    lab<-spec.clu(sd_matrix,k)
    results <- list()

    results[["clu"]] <- lab
    results[["S"]] <- sd_matrix

    if(!is.null(es[[1]])){

    results[["es"]]<-es

    }

    return(results)
  }

