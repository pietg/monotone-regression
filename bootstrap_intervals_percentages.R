########################################################
####   MONOTONE REGRESSION        ######
########################################################

	library(Rcpp)
	sourceCpp("CI_SLSE.cpp")
	sourceCpp("CI_NW.cpp")
  	NumIt = 1000
  	ngrid = 100
  	n = 100
  	m = 2
  	sigma = 0.1

  	f <- function(x) {x^2+x/5}
  	#f <- function(x) {exp(4 * (x - 0.5))/(1+exp(4 * (x - 0.5)))}
  	x0 <-seq(0,1,by=0.01)
   	y0<-f(x0)
   	
   	percentages1 <- integer(ngrid-1)
   	percentages2 <- integer(ngrid-1)
   	
   	for (k in 1:(ngrid-1))
   	{
   		percentages1[k]<-0
   		percentages2[k]<-0	
   	}
   	
   	grid <- numeric(ngrid)
   	
   	for (i in (0:ngrid))
   	{
   		grid[i] = i/ngrid
   	}
   	
   	B <- matrix(0, nrow = ngrid, ncol = 4)
   	C <- matrix(0, nrow = ngrid, ncol = 4)
  	
 for (j in 1: NumIt)
 {
  	sim = 1000+j
  	
  	print(j)
   
  	set.seed(sim)
  	
  	a = runif(n,0,1) 
  	b = f(a) + rnorm(n,0,sigma)
  
  	X = matrix(cbind(a,b),n,m, byrow = FALSE)
	
  	output1 <- CI_SLSE(X,n,sim)
  	output2 <- CI_NW(X,n,sim)
  	
  	B <- output1$CI_SLSE
  	C <- output2$CI_NW
  	
  	for (k in 1:(ngrid-1))
  	{
  		if (B[k,3]>f(grid[k]) || B[k,4]<f(grid[k]))
  		{
  			percentages1[k] <- percentages1[k]+1
  		}
  	}
  	
  	for (k in 1:(ngrid-1))
  	{
  		if (C[k,3]>f(grid[k]) || C[k,4]<f(grid[k]))
  		{
  			percentages2[k] <- percentages2[k]+1
  		}
  	}
  }
   	
   	z1<-1-percentages1/NumIt
   	z2<-1-percentages2/NumIt
    	x1 <- seq(0.01,0.99,by=0.01) 
   	  	
   	plot(c(-10000,-10000),xlim=c(0.0,1), ylim=c(0.85,1), main= "", ylab="",xlab="",bty="n",las=1)
	lines(x1,z1,lty=1,lwd=2,col="blue")
	lines(x1,z2,lty=1,lwd=2,col="red")
	lines(c(0.0,1),c(0.95,0.95),lwd=2)
	



