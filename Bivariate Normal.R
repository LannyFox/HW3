bivariate.normal <- function(M1, M2, S1, S2, p1){
  if (runif(1)<= p1){
    return(rnorm(1,mean=M1,sd=sqrt(S1)))
  } else {
    return(rnorm(1,mean=M2,sd=sqrt(S2)))
  }}
X <- NULL
set.seed(40)
for (i in 1:100){
  X[i] <- bivariate.normal(M1=0, M2=3,S1=1,S2=4,p1=.7)
}
plot(density(X), xlab='Kernel Density plot of X',main='')