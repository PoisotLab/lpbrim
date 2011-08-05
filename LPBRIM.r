library(nnet)
library(latticeExtra)
library(snow)
library(snowfall)

source('common.r')
source('plots.r')

bLP = function (x,as.adjacency=TRUE) {
	cat('Label propagation starting\n')
	if(as.adjacency) x[x>0] <- 1
    OrderVec <- c(rownames(x),colnames(x))
  	x <- x[sample(c(1:nrow(x))),sample(c(1:ncol(x)))]
   	# HTL labels
   	lT <- c(1:nrow(x))
   	names(lT) <- rownames(x)
   	# LTL labels
   	## Set to NA for initialization of the graph
   	lB <- rep(0,ncol(x))
   	names(lB) <- colnames(x)
 	## Seeding the initial modularity
 	oldQ <- 0
 	newQ <- 1e-10
 	Nsteps <- 0
   	# Labels propagation loop
   	while(oldQ<newQ)
   	{
   		Nsteps <- Nsteps + 1
   		oldQ <- newQ
   		## Step 1 : update lB
   		for(lsp in 1:ncol(x))
   		{
   			Nei <- rownames(x)[x[,lsp]>0]
	  		NeiLab <- lT[Nei]
   			if(as.adjacency)
   			{
	   			lB[lsp] <- NeiLab[mostFrequent(NeiLab)]
   			} else {
   				LS <- x[Nei,lsp]
				lB[lsp] <- NeiLab[mostFrequentW(NeiLab,LS)]
   			}   			
   		}
   		names(lB) <- colnames(x)
   		## Step 2 : update lT
   		for(tsp in 1:nrow(x))
   		{
   			Nei <- colnames(x)[x[tsp,]>0]
     		NeiLab <- lB[Nei]
			if(as.adjacency)
   			{
	   			lT[tsp] <- NeiLab[mostFrequent(NeiLab)]
   			} else {
   				LS <- x[tsp,Nei]
				lT[tsp] <- NeiLab[mostFrequentW(NeiLab,LS)]
   			}
   		}
   		names(lT) <- rownames(x)
	   ## Shaping the vectors
   		Modules <- c(lT,lB)
	   	Comms <- unique(Modules)
	 	NComm <- length(Comms)
	 	Smat = matrix(0,ncol=NComm,nrow=sum(dim(x)))
	 	colnames(Smat) <- Comms
	 	rownames(Smat) <- c(rownames(x),colnames(x))
		for(i in 1:length(Modules)) Smat[names(Modules)[i],as.character(Modules[i])]<-1
	   	newQ <- Qbip(x,Smat)
	   	cat('\rCurrent Q\t',newQ,'\t')
   	}
	Modules <- as.numeric(as.factor(Modules))
	names(Modules) <- c(names(lT),names(lB))
	cat('\n')
	cat('Label propagation ended after',Nsteps,'step(s)\n')
	return(Modules[OrderVec])
}



Qbip = function(x,s)
{
	Q <- 0
	x[x>0] <- 1
 	p <- nrow(x)
 	h <- ncol(x)
 	m <- sum(x)
 	nc <- ncol(s)
 	A <- x
 	P <- x
 	for(i in 1:nrow(x)) for(j in 1:ncol(x)) P[i,j] <- (sum(x[i,])*sum(x[,j]))/m
 	B <- A-P	
 	Rm <- s[c(1:p),]
 	Tm <- s[c(p+(1:h)),]
 	## Induce Qr from T
 	BT <- B%*%Tm
 	Isum <- NULL
 	for(i in 1:p)
 	{
 		Ksum <- NULL
 		for(k in 1:nc)
 		{
 			Ksum[k] <- Rm[i,k]*(BT[i,k])
 		}
   		Isum[i] <- sum(Ksum)
 	}
 	Q = (1/m)*sum(Isum)	 
	return(Q)
}

bBRIM = function(x)
{
 	Nsteps <- 0
	CommDiv <- bLP(x)
	x[x>0] <- 1
 	Comms <- unique(CommDiv)
 	NComm <- length(Comms)
 	Smat = matrix(0,ncol=NComm,nrow=sum(dim(x)))
 	colnames(Smat) <- Comms
 	rownames(Smat) <- c(rownames(x),colnames(x))
	## Fill the S matrix
 	for(i in 1:length(CommDiv)) Smat[names(CommDiv)[i],as.character(CommDiv[i])]<-1
	## Initial modularity 
	FromR <- TRUE
	cBM <- -10
	preBM <- 10
 	## Some important values
 	p <- nrow(x)
 	h <- ncol(x)
 	m <- sum(x)
 	nc <- ncol(Smat)
 	A <- x
 	P <- x
 	for(i in 1:nrow(x)) for(j in 1:ncol(x)) P[i,j] <- (sum(x[i,])*sum(x[,j]))/m
 	B <- A-P
 	## Optimization loop
  	cat('\nBRIM optimization starting\n')
	while((Nsteps<3)|(preBM!=cBM))
	{
 		Nsteps <- Nsteps + 1
  		preBM <- Qbip(x,Smat)
  		cat('\rCurrent Q\t',preBM,'\t')
		## Modularity matrix
		Rm <- as.matrix(Smat[c(1:p),])
 		Tm <- as.matrix(Smat[c(p+(1:h)),])
 		## Product matrix for T & R
 		rBT <- B%*%Tm
		tBT <- t(B)%*%Rm
	
 		## Optimization
		if(FromR)
		{
			for(i in 1:nrow(x))
			{
				Rm[i,] <- rep(0,ncol(Smat))
				Rm[i,which.max(rBT[i,])] <- 1
			}
		} else {
			for(i in 1:ncol(x))
			{
				Tm[i,] <- rep(0,ncol(Smat))
				Tm[i,which.max(tBT[i,])] <- 1
			}	
		}
   		Smat[c(1:nrow(x)),] <- Rm
  		Smat[(nrow(x)+c(1:ncol(x))),] <- Tm
   		## New bipartition
  		cBM <- Qbip(x,Smat)
   		FromR <- !FromR
	}
	cat('\n')
	cat('BRIM convergence reached in',Nsteps,'step(s)\n')
	return(list(S=Smat,M=x,Q=cBM,c=ncol(Smat)))
}

findModules = function(M,iter=50,cpu=1)
{
	usePar <- ifelse(cpu>1,TRUE,FALSE)
	sfInit(parallel=usePar,cpus=cpu,type='SOCK')
 	if(usePar)
 	{
 		sfLibrary(nnet)
   		sfExportAll()
 	}
 	ModulOutput <- sfLapply(c(1:iter),function(x) bBRIM(M))
  	sfStop()	
 	Qs <- unlist(lapply(ModulOutput,function(x)x$Q))
	#cs <- unlist(lapply(ModulOutput,function(x)x$c))
   	#trellis.device('pdf')
   	#print(xyplot(Qs~cs,type=c('g','p','a'),jitter.x=TRUE))
   	#dev.off()
   	maxQs <- which.is.max(Qs)
 	return(ModulOutput[[maxQs]])
}