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
   if(is.na(w)) w <- rep(1, length(vec))
   nvec <- vec
   names(nvec) <- names(vec)
   for(i in 1:length(vec)) nvec[i] <- sum(vec==vec[i])
   return(which.is.max(nvec))
}
