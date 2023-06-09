########################################################
####   MONOTONE REGRESSION        ######
########################################################

	library(Rcpp)
	sourceCpp("SLSE.cpp")
	sourceCpp("Nadaraya_Watson.cpp")
  	NumIt = 1000
  	ngrid = 100
  	n = 500
  	m = 2
  	sigma = 0.1

  	f <- function(x) {x^2+x/5}
  	#f <- function(x) {exp(4 * (x - 0.5))/(1+exp(4 * (x - 0.5)))}
  	x0 <-seq(0,1,by=0.01)
   	y0<-f(x0)
   	
   	percentages <- integer(ngrid-1)
   	grid <- numeric(ngrid)
   	
   	for (i in (0:ngrid))
   	{
   		grid[i] = i/ngrid
   	}
   	B <- matrix(0, nrow = ngrid, ncol = 2)
   	
   	outputMat <- matrix(0, nrow= NumIt, ncol= 2, byrow = FALSE)
   	colnames(outputMat) <- c("SLSE", "NW")
  	
 for (j in 1: NumIt)
 {
  	sim = j
  	
  	print(j)
   
  	set.seed(sim)
  	
  	a = runif(n,0,1) 
  	b = f(a) + rnorm(n,0,sigma)
  
  	X = matrix(cbind(a,b),n,m, byrow = FALSE)
	
  	output <- SLSE(X,n)
  	output2 <- NW(X,n)
  	
  	outputMat[j,] <- c(output$SLSE[60,2],output2$SLSE[60,2])
  	
  }
  
  outputMat

	pdf("BoxPlot_value.pdf")
	boxplot(outputMat,las=1)
	abline(h=f(grid[60]),lwd=2,lty=1,col = "red")
	dev.off()


