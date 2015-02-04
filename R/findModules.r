#' @title Find modules
#' @description This function takes a matrix, and performs a search for the best
#' partition in modules.
#' @export
findModules = function(M,iter=50,sparse=TRUE, ...)
{
   if(is.null(rownames(M))) rownames(M) <- paste('r',c(1:NROW(M)),sep='')
   if(is.null(colnames(M))) colnames(M) <- paste('c',c(1:NCOL(M)),sep='')
   if(sparse) M <- Matrix(M, sparse=TRUE)
   ModulOutput <- alply(c(1:iter),1, function(x) bBRIM(M), ...)
   Qs <- unlist(lapply(ModulOutput,function(x)x$Q))
   maxQs <- which.is.max(Qs)
   return(ModulOutput[[maxQs]])
}
