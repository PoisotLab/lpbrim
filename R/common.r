#' @export
spread = function (v, m = 0, M = 1)
{
   v <- v - min(v)
   v <- v/max(v)
   v <- v * (M - m)
   v <- v + m
   return(v)
}

#' @export
mostFrequent = function(vec, w=NA)
{
   if(is.na(w)) w <- rep(1, length(vec)) # TODO use weights at some point?
   tvec <- table(vec)
   nvec <- as.vector(tvec[as.factor(vec)])
   return(which.max(nvec))
}

#' @export
netToAdjacency = function(net){
    x <- as.factor(net[,1])
    y <- as.factor(net[,2])
    z <- net[,3]
    adj <- matrix(0,
                  nrow=nlevels(x),
                  ncol=nlevels(y),
                  dimnames=list(levels(x), levels(y)))
    adj[cbind(x, y)] <- z
    return(adj)
}
