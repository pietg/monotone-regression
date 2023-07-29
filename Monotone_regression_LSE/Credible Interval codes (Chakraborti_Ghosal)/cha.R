library("msm")
library("fdrtool")
library("dplyr")
library("isotone")

load("Ainv.rda")

# isotone2 function: to evaluate the isotonic regression estimator
# the vector y should be sorted wrt x's
#return values: r gives the isotonized step-function, t is its value at t1
isotone2 <- function(y, w, t1){
  J <- length(w)
  xval <- c(1:J) / J
  
  index.w <- which(w==0)
  w0 <- w
  
  y <- y[w>0]
  xval <- xval[w>0]
  w <- w[w>0]
  
  fit <- gpava(xval, y, weights = w, solver = weighted.mean,
                                      ties = "primary", p = NA)
  
  r <- stepfun(xval, c(fit[[1]][1], fit[[1]]))
  t <- r(t1) #value of the isotonized step function at given t1
  return(list(r, t))
}
####################################################################

# true regression functions 
#define f1
f1 <- function(x){
  return (x^2 + x/5)
}

#define f2
f2 <- function(x){
  xx <- exp(4 * (x - 0.5))
  return (xx / (1 + xx))
}
##########################

# obtain Nj, yBar and sigmahat
extract.nj.sigmahat <- function(x, y, J, zi, tausq){
  n <- length(x)
  ybar <- c(1:J); xbar <- c(1:J) / J; nj <- c(1:J);
  sigmaHatSq <- NULL;
  
  for(j in 1:J){
    xInIj <- x > (j-1)/J & x <= j/J #indices present in j-th interval.
    nj[j] <- sum(xInIj) #number of Xj's in j-th interval.
    
    if(nj[j] == 0){
      ybar[j] <- 0
      next # go to next j if ther's no data in this interval
    }
    else{
      yInIj <-  y[xInIj]
      ybar[j] <- mean(yInIj) # define ybar
      sigmaHatSqIj <- t(yInIj-zi[j])%*%solve(diag(nj[j])+
                                               matrix(tausq[j],nj[j],nj[j]))%*%(yInIj-zi[j])
      sigmaHatSq <- c(sigmaHatSq, sigmaHatSqIj)
    }
  }
  sigmaHatSq <- sum(sigmaHatSq)/n
  sigmaHat <- sqrt(sigmaHatSq)
  return(list(nj, ybar, sigmaHat))
}
##########################################################################################

#main function to obtain a projection-posterior credible interval for a given x_0
projection.credible.interval <- function(x, y, x0 = 0.5, credibility.level=0.95,
                                         plot.function=TRUE, nPostSamp=1000){
  n <- length(x)
  
  #sort y according to x
  xy.matrix.sorted <- arrange(as.data.frame(cbind(x, y)), x)
  x <- xy.matrix.sorted[, 1]
  y <- xy.matrix.sorted[, 2]
  
  # set J 
  #J <- floor((n ^ (1/3))*log(n))
  #J <- floor((n ^ (2/3))*log(n))
  J<-n
  
  #prior hyperparameters
  zi <- rep(0,J); tausq <- rep(100, J)
  
  # calculate Nj, YjBar and sigmahat
  obj1 <- extract.nj.sigmahat(x, y, J, zi, tausq)
  nj <- obj1[[1]]
  ybar <- obj1[[2]]
  sigmaHat <- obj1[[3]]
  
  alpha.level <-  (1 - credibility.level) / 2
  corrected.cred <- A.inv(credibility.level)
  correct.level <- (1-corrected.cred) / 2
  
  # posterior mean and variance
  postMean <- (nj*ybar + zi/tausq)/(nj + 1/tausq)
  postSd <- sigmaHat/sqrt(nj + 1/tausq)
  
  if(plot.function==TRUE){
    # plot the data 
    plot(x, y, xlab="x", ylab="y",type='n')
    xx <- seq(0.001, 0.999, 0.001)
    function.sample <- NULL
    sample.fx0 <- NULL
    for(i in 1: nPostSamp){
      theta <- rnorm(J, postMean, postSd)
      theta1<- isotone2(theta, nj, x0) 
      sample.fx0 <- c(sample.fx0, theta1[[2]])
      
      yy1 <- theta1[[1]](xx)
      function.sample <- cbind(function.sample, yy1)
      #lines(xx, yy1, xlab="x", ylab="fstar(x)", col="green")
    }
    #lines(xas,f1(xas))
    lines(xx, apply(function.sample,1, FUN=function(x){quantile(x, 1-alpha.level)}))
    lines(xx, apply(function.sample,1, FUN=function(x){quantile(x, alpha.level)}))
    lines(xx, rowMeans(function.sample), lty="dashed")
    legend(0.02, max(y), legend=c(paste(100*credibility.level, "% projection-posterior credible limits"),"Projection-posterior mean"), lty=1:2, cex=0.7)
  }
  
  if(plot.function != TRUE){
    sample.fx0 <- NULL
    # draw samples from projection-posterior
    for(i in 1: nPostSamp){
      theta <- rnorm(J, postMean, postSd)
      theta1<- isotone2(theta, nj, x0) 
      sample.fx0 <- c(sample.fx0, theta1[[2]])
    }
  }
    
  lower.limit <- quantile(sample.fx0, alpha.level)
  upper.limit <- quantile(sample.fx0, 1 - alpha.level)
  interval.obs <- cbind(lower.limit, upper.limit)
  
  lower.limit <- quantile(sample.fx0, correct.level)
  upper.limit <- quantile(sample.fx0, 1 - correct.level)
  interval.correct <- cbind(lower.limit, upper.limit)
  
  return(list(interval.obs, interval.correct))
    
}

set.seed(1010)
n <- 1000
sigma <- 0.1
x0 <- 0.5
f <- f2
x <- runif(n)
y <- f(x) + rnorm(n, 0, sigma)
cred0 <- projection.credible.interval(x, y, x0 = x0, credibility.level=0.95,
                           plot.function=TRUE, nPostSamp=1000)

#The $100(1-\alpha)\%$ projection-posterior credible interval is
 cred0[[1]] # unadjusted credible interval
                           
#The corrected projection-posterior credible interval with asymptotic coverage  $100(1-\alpha)\%$ is
cred0[[2]]    # adjusted credible interval

#The true value of $f(x_0)$ is
f(x0)

                     
