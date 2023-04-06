########################################################
####   MONOTONE REGRESSION        ######
########################################################

	library(Rcpp)
	sourceCpp("bandwidth_choice.cpp")

  	n = 1000
  	m = 2
  	sigma = 0.1

  	f <- function(x) {x^2+x/5}
  	x0 <-seq(0,1,by=0.01)
   	y0<-f(x0)
   	
  	
  	sim = 1000
   
  	set.seed(sim)
  	
  	a = runif(n,0,1) 
  	b = f(a) + rnorm(n,0,sigma)
  
  	X = matrix(cbind(a,b),n,m, byrow = FALSE)
	
  	output <- bandwidth(X,n,sim)
  	
  	B <- output$MSE
  	     	
   	x<-B[,1]
   	y<-B[,2]
   	  	
   	plot(c(-10000,-10000),xlim=c(0,1), ylim=c(min(y),max(y)), main= "", ylab="",xlab="",bty="n",las=1)
	lines(x,y,lty=1,lwd=2,col="blue")
	
	
	