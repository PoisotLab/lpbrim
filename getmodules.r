## Extract the within-modules networks
getmodules = function(mod)
{
	Lev1 <- mod$S[c(1:nrow(mod$M)),]
	Lev2 <- mod$S[c((nrow(mod$M)+1):nrow(mod$S)),]
	LoW = NULL
	for(i in 1:ncol(mod$S))
	{
		LoW[[i]] <- mod$M[(Lev1[,i]==1),(Lev2[,i]==1)]
	}
	return(LoW)
}