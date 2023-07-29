	library(Rcpp)
	sourceCpp("correlation_percentile.cpp")
	
	n <- 1000
	
 	output <- correlation_percentile(n)
   
	B <- output$LSE
	C <- output$CI_LSE
	D <- output$percentages
	E <- output$data
	
	#B<-read.table("SLSE.txt")
  	#C<-read.table("CI_LSE.txt")
  	#D<-read.table("percentages.txt")
  	#E<-read.table("data.txt")
  
   x1<-B[,1]
   z1<-B[,2]
   x2<-C[,1]
   y1<-C[,2]
   y2<-C[,3]
   y3<-C[,4]
   u1<-E[,1]
   v1<-E[,2]
   
   f <- function(x) {x^2+x/5}
   x0 <-seq(0,1,by=0.01)
   y0<-f(x0)
   
   plot(c(-100,-100),xlim=c(0,1), ylim=c(min(y2),max(y3)), main= "",ylab="",xlab="",bty="n",las=1)
   #lines(x1,z1,lwd=2,type='s')
   lines(x2,y1,lwd=2,lty=1,col="blue")
   lines(x0, y0,lwd=2,lty=2,col="red")
   lines(x2,y2,lwd=1,lty=1,type='s')
   lines(x2,y3,lwd=1,lty=1,type='s')
   segments(x2,y2,x2,y3)
   #points(u1,v1,pch = 20) 

   
   
    u<-D[,1]
   	v<-1-D[,2]
   	  	
   	plot(c(-10000,-10000),xlim=c(0.0,1.0), ylim=c(0.5,1), main= "", ylab="",xlab="",bty="n",las=1)
	lines(u,v,lty=1)
	lines(c(0,1),c(0.95,0.95),col="red")
	#lines(c(0,1),c(0.96324,0.96324),col="red")
