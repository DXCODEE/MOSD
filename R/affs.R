  #' affs Function
  #'
  #' constructing  local scaling affinity
  #' @param x a matrix where the rows represent samples, and the columns represents the features
  #' @return affinity graph
  #' @export
 

  affs<-function(x,na.action = na.omit)
  {
    x <- na.action(x)
    rown <- rownames(x)
    x <- as.matrix(x)
    m <- nrow(x)

    s <- rep(0,m)
    dota <- rowSums(x*x)/2
    dis <- crossprod(t(x))
    for (i in 1:m)
      dis[i,]<- 2*(-dis[i,] + dota + rep(dota[i],m))
    ## fix numerical prob.
    dis[dis < 0] <- 0
    for (i in 1:m)
      s[i] <- median(sort(sqrt(dis[i,]))[1:5])
    ## Compute Affinity Matrix
    km <- exp(-dis / s%*%t(s))
    km
    # return(list(km=km))

  }
