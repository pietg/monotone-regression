########################################################
####   MONOTONE REGRESSION        ######
########################################################

	library(Rcpp)
	sourceCpp("bandwidth_choice.cpp")

  	n = 500
  	m = 2
  	sigma = 0.1
  	B=10000

  	f <- function(x) {x^2+x/5}
  	#f <- function(x) {exp(4 * (x - 0.5))/(1+exp(4 * (x - 0.5)))}
  	x0 <-seq(0,1,by=0.01)
   	y0<-f(x0)
   	
  	
  	sim = 1000
   
  	set.seed(sim)
  	
  	a = runif(n,0,1) 
  	b = f(a) + rnorm(n,0,sigma)
  
  	X = matrix(cbind(a,b),n,m, byrow = FALSE)
	
  	output <- bandwidth(X,n,sim,B)
  	
  	C <- output$MSE
  	     	
   	x<-C[,1]
   	y<-C[,2]
   	  	
   	plot(c(-10000,-10000),xlim=c(0,1), ylim=c(min(y),max(y)), main= "", ylab="",xlab="",bty="n",las=1)
	lines(x,y,lty=1,lwd=2,col="blue")
	
	
	