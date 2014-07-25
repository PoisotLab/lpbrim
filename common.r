spread = function (v, m = 0, M = 1) 
{
   v <- v - min(v)
   v <- v/max(v)
   v <- v * (M - m)
   v <- v + m
   return(v)
}

mostFrequent = function(vec, w=NA)
{
   if(is.na(w)) w <- rep(1, length(vec)) # what is this supposed to do? w isn't used anywhere else in the function!
   tvec <- table(vec)
   nvec <- as.vector(tvec[as.factor(vec)])
   return(which.is.max(nvec))
}

netToAdjacency = function(net){
    x <- as.factor(net[,1])
    y <- as.factor(net[,2])
    z <- net[,3]
    adj <- matrix(0,
                  nrow=nlevels(x),
                  ncol=nlevels(y),
                  dimnames=list(levels(x), levels(y)))
    adj[cbind(x, y)] <- z
    adj
}
