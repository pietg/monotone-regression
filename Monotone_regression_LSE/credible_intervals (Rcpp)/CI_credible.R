    library(Rcpp)
	sourceCpp("credible_intervals.cpp")
	
	f <- function(x) {x^2+x/5}
	#f <- function(x) {exp(4*(x-0.5))/(1+exp(4* (x-0.5)))}
   	x0 <-seq(0,1,by=0.01)
   	y0<-f(x0)

	
	n=1000
	J=floor(n^(1/3)*log(n))
	#J=floor(n^(1/2))
	zeta <- rep(0,J)
	#zeta <- f(seq(0,1,length=J)) 
	lambda_sq <- rep(100,J)
	
  	output <- credible(n,J,zeta,lambda_sq)
   
	#B <- output$posterior_mean
	#C <- output$CI_LSE
	#D <- output$percentages
	
	B<-read.table("posterior_mean.txt")
  	C<-read.table("CI_LSE.txt")
  	D<-read.table("percentages.txt")
  
   x1<-B[,1]
   z1<-B[,2]
   x2<-C[,1]
   y1<-C[,2]
   y2<-C[,3]
   y3<-C[,4]
   
     
   plot(c(-100,-100),xlim=c(0,1), ylim=c(min(z1),max(y3)), main= "",ylab="",xlab="",bty="n",las=1)
   lines(x1,z1,lwd=2,type ="s",col="blue")
   lines(x2,y1,lwd=2,col="blue",type='s')
   lines(x0, y0,lwd=2,lty=2,col="red")
   lines(x2,y2,lwd=1,lty=1,type='s')
   lines(x2,y3,lwd=1,lty=1,type='s')
   segments(x2,y2,x2,y3)

   
    
    u<-D[,1]
   	v<-1-D[,2]
   	  	
   	plot(c(-10000,-10000),xlim=c(0.0,1.0), ylim=c(0.75,1), main= "", ylab="",xlab="",bty="n",las=1)
	lines(u,v,lty=1)
	#lines(c(0,1),c(0.95,0.95),col="red")
	lines(c(0,1),c(0.96324,0.96324),lwd=2,col="red")
