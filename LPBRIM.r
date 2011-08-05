library(nnet)
library(latticeExtra)
library(snow)
library(snowfall)

spread = function (v, m = 0, M = 1) 
{
    v <- v - min(v)
    v <- v/max(v)
    v <- v * (M - m)
    v <- v + m
    return(v)
}

mostFrequent = function(vec)
{
	nvec <- vec
 	names(nvec) <- names(vec)
	for(i in 1:length(vec)) nvec[i] <- sum(vec==vec[i])
 	return(which.is.max(nvec))
}

mostFrequentW = function(v,w)
{
	nvec <- v
 	names(nvec) <- names(v)
	for(i in 1:length(v)) nvec[i] <- sum(w[v==v[i]])
 	return(which.is.max(nvec))
}

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

plotModules = function(mod)
{
	opar <- par(no.readonly=TRUE)
	##
	#TX <- spread(x[x>0],0.5,3)
 	x <- mod$M
	x[x>0] <- 1
	##
 	par(mar=c(0,0,0,0),xaxt='n',yaxt='n')
	##
 	S <- mod$S
 	Modules <- numeric(nrow(S))
 	names(Modules) <- rownames(S)
 	for(i in 1:nrow(S))
 	{
 		Modules[i] <- as.numeric(colnames(S)[S[i,]==1])
 	}
 	ModNum <- as.numeric(as.factor(Modules))
 	names(ModNum) <- names(Modules)
 	Tm <- ModNum[rownames(x)]
 	Bm <- ModNum[colnames(x)]
 	## Initate Plot
 	CommColor <- colorRampPalette(brewer.pal(11,'Spectral'))(ncol(S))
 	plot(0,xlim=c(-1,2),ylim=c(0,1),pch=NA,bty='n')
 	##
 	SeqT <- rank(Tm,ties.method='random')
 	SeqB <- rank(Bm,ties.method='random')
 	## Within each module …
 	for(cmod in unique(ModNum))
 	{
 		## … degree-ranking of the HTL species
 		Tcurr <- SeqT[Tm==cmod]
   		if(length(Tcurr)>1)
   		{
   			Gen <- apply(x[names(Tcurr),],1,function(x) sum(x>0))
	   		WCrank <- rank(Gen,ties.method='random')
	   		Tcurr <- (min(Tcurr)-1)+WCrank
	   		SeqT[Tm==cmod] <- Tcurr
   		}
   		
   		## … degree-ranking of the HTL species
 		Bcurr <- SeqB[Bm==cmod]
   		if(length(Bcurr)>1)
   		{
	 		Gen <- apply(x[,names(Bcurr)],2,function(x) sum(x>0))
	   		WCrank <- rank(Gen,ties.method='random')
	   		Bcurr <- (min(Bcurr)-1)+WCrank
	   		SeqB[Bm==cmod] <- Bcurr
   		}
 	}
 
 	SeqT <- spread(SeqT)
 	SeqB <- spread(SeqB)
 
 	for(ts in 1:nrow(x)) for(bs in 1:ncol(x)) if(x[ts,bs]>0)
 	{
 		ccol <- ifelse(Tm[ts]==Bm[bs],CommColor[Tm[ts]],'grey')
 		clwd <- ifelse(Tm[ts]==Bm[bs],2,1)
 		segments(0,SeqB[bs],1,SeqT[ts],lwd=clwd,col=ccol)
 	}
 
 	## Plot Top Species
 	points(rep(1,nrow(x)),SeqT,col=CommColor[Tm],pch=19,cex=1.4)
 	text(rep(1.1,nrow(x)),SeqT,rownames(x),cex=0.8,adj=c(0,0.5))
 	## Plot Bottom Species
 	points(rep(0,ncol(x)),SeqB,col=CommColor[Bm],pch=19,cex=1.4)
	text(rep(-0.1,ncol(x)),SeqB,colnames(x),cex=0.8,adj=c(1,0.5))
 	## 
 	par(opar)
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

plotMatrixModules = function (mod,mode='blocks') {	
	x <- mod$M
 	S <- mod$S
 	Modules <- numeric(nrow(S))
 	names(Modules) <- rownames(S)
 	for(i in 1:nrow(S))
 	{
 		Modules[i] <- as.numeric(colnames(S)[S[i,]==1])
 	}
 	ModNum <- as.numeric(as.factor(Modules))
 	names(ModNum) <- names(Modules)
 	Tm <- ModNum[rownames(x)]
 	Bm <- ModNum[colnames(x)]
	
	opar <- par(no.readonly=TRUE)
 	par(xaxt='n',yaxt='n',mar=c(0,0,0,0),mai=c(0.2,0.2,0.2,0.2))
 
 	CommColor <- colorRampPalette(brewer.pal(11,'Paired'))(ncol(S))
 
 	SeqT <- rank(Tm,ties.method='random')
 	SeqB <- rank(Bm,ties.method='random')
 
 	## Within each module …
 	for(cmod in unique(ModNum))
 	{
 		## … degree-ranking of the HTL species
 		Tcurr <- SeqT[Tm==cmod]
   		if(length(Tcurr)>1)
   		{
   			Gen <- apply(x[names(Tcurr),],1,function(x) sum(x>0))
	   		WCrank <- rank(Gen,ties.method='random')
	   		Tcurr <- (min(Tcurr)-1)+WCrank
	   		SeqT[Tm==cmod] <- Tcurr
   		}
   		
   		## … degree-ranking of the HTL species
 		Bcurr <- SeqB[Bm==cmod]
   		if(length(Bcurr)>1)
   		{
	 		Gen <- apply(x[,names(Bcurr)],2,function(x) sum(x>0))
	   		WCrank <- rank(Gen,ties.method='random')
	   		Bcurr <- (min(Bcurr)-1)+WCrank
	   		SeqB[Bm==cmod] <- Bcurr
   		}
 	}
 
 	Tm <- Tm[names(SeqT)]
 	Bm <- Bm[names(SeqB)]
 
plot(0,pch=NA,xlim=c(0.9,ncol(x)+0.1),ylim=c(0.9,nrow(x)+0.1),asp=1,xlab='',ylab='',bty='n')
	rect(0.2,0.2,ncol(x)+0.8,nrow(x)+0.8)
   	## FRAMES
   	if(mode%in%c('frames','both'))
   	{
   		Ucom <- unique(Tm)
   		for(cc in 1:length(Ucom))
   		{
   			Ccom <- Ucom[cc]
   			Tpos <- range(SeqT[Tm==Ccom])
	  		Bpos <- range(SeqB[Bm==Ccom])
	  		Tpos <- Tpos + c(-0.45,0.45)
	  		Bpos <- Bpos + c(-0.45,0.45)
	  		## Module
	  		if(mode=='frames') rect(Bpos[1],Tpos[1],Bpos[2],Tpos[2],lwd=2)
	  		if(mode=='both') rect(Bpos[1],Tpos[1],Bpos[2],Tpos[2],border=CommColor[Ccom],lwd=2)
   		}
   	}
   	## Squares
	for(i in 1:length(SeqT))
	{
	  for(j in 1:length(SeqB))
	  {
	  	Tcoord <- SeqT[i]
		Bcoord <- SeqB[j]
	  	ccol <- 'grey'
	  	if(mode!='frames')
	  	{
	  		ccol <- ifelse(Tm[i]==Bm[j],CommColor[Tm[i]],'grey')
	  	} else {
	  		ccol <- ifelse(Tm[i]==Bm[j],'black','grey')
	  	}
	  	if(x[i,j]>0) symbols(Bcoord,Tcoord,squares=0.7,add=TRUE,inches=FALSE,bg=ccol,fg=NA)
	  }
	} 
	par(opar)
}