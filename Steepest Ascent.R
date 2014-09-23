steepest.ascent <- function(X, start, step=.0005){ 
  
  #X is the vector with all sample values xi
  #start should be given as a 5-element vector c(p, mu1, mu2, sigma1, sigma2)
  #I only had one p value because p2 = 1-p1
  #M[1]=p, M[2]=mu1, M[3]=mu2, M[4]=sigma1, M[5]=sigma2
  
  M <- start #M matrix will store each step in a separate row
  
  # Q will be the denominator in a lot of the gradient's elements, 
  # so calculating it separately will simplify the code overall
  Q <- M[1]*((2*pi*M[4])^-.5)*exp(-(X-M[2])^2/(2*M[4])) + (1-M[1])*((2*pi*M[5])^-.5)*exp(-(X-M[3])^2/(2*M[5]))
  
  
  grad <- c(sum((((2*pi*M[4])^-.5)*exp(-(X-M[2])^2/(2*M[4])) - ((2*pi*M[5])^-.5)*exp(-(X-M[3])^2/(2*M[5])))/Q), #p portion of gradient
            sum((M[1]*((2*pi*M[4])^-.5)*exp(-(X-M[2])^2/(2*M[4]))*2*(X-M[2]))/Q), #mu1 portion of gradient
            sum((1-M[1])*((2*pi*M[5])^-.5)*exp(-(X-M[3])^2/(2*M[5]))*2*(X-M[3])/Q), #mu2 portion of gradient
            sum((M[1]/(sqrt(2*pi)*Q*M[4]^4))*exp(-(X-M[2])^2/(2*M[4]))*(M[4]^2.5*(-X+M[2])-.5*M[4]^2.5 +.5*M[4]^1.5*(X-M[2])^2)),#sigma1
            sum(((1-M[1])/(sqrt(2*pi)*Q*M[5]^4))*exp(-(X-M[3])^2/(2*M[5]))*(M[5]^2.5*(-X+M[3])-.5*M[5]^2.5 +.5*M[5]^1.5*(X-M[3])^2))#sigma2
  )
  
  norm.grad <- sqrt(sum(grad^2))
  slope <- grad/norm.grad
  
  ## The count variable will be used to make the program stop in the event that it gets stuck between two values
  count <- 1 
  while (norm.grad > .01 & count<10000){ #while loop will run until normalized gradient is sufficiently small
    M <- M+slope*step
    Q <- M[1]*((2*pi*M[4])^-.5)*exp(-(X-M[2])^2/(2*M[4])) + (1-M[1])*((2*pi*M[5])^-.5)*exp(-(X-M[3])^2/(2*M[5]))
    grad <- c(sum((((2*pi*M[4])^-.5)*exp(-(X-M[2])^2/(2*M[4])) - ((2*pi*M[5])^-.5)*exp(-(X-M[3])^2/(2*M[5])))/Q), #p
              sum((M[1]*((2*pi*M[4])^-.5)*exp(-(X-M[2])^2/(2*M[4]))*2*(X-M[2]))/Q), #mu1
              sum((1-M[1])*((2*pi*M[5])^-.5)*exp(-(X-M[3])^2/(2*M[5]))*2*(X-M[3])/Q), #mu2
              sum((M[1]/(sqrt(2*pi)*Q*M[4]^4))*exp(-(X-M[2])^2/(2*M[4]))*(M[4]^2.5*(-X+M[2])-.5*M[4]^2.5 +.5*M[4]^1.5*(X-M[2])^2)),#sigma1
              sum(((1-M[1])/(sqrt(2*pi)*Q*M[5]^4))*exp(-(X-M[3])^2/(2*M[5]))*(M[5]^2.5*(-X+M[3])-.5*M[5]^2.5 +.5*M[5]^1.5*(X-M[3])^2))#sigma2
    )
    norm.grad <- sqrt(sum(grad^2))
    slope <- grad/norm.grad
   # print(norm.grad) # this print statement I used to diagnose problems with the model by printing different varaibles as it ran.
    count <- count+1
  }
  return(data.frame(p=M[1], Mu1=M[2],Mu2=M[3],sigma1=M[4],sigma2=M[5]))  
}