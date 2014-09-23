steepest.ascent <- function(X, start, step=.0001){ 
  
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
            sum((M[1]/Q)*(exp(-(X-M[2])/(2*M[4]))*(-pi)*(2*pi*M[4])^(-1.5)+(2*pi*M[4])^-.5*exp(-(X-M[2])^2/(2*M[4]))*(X-M[2])^2/(2*M[4]^2))),#sigma1 portion
            sum(((1-M[1])/Q)*(exp(-(X-M[3])/(2*M[5]))*(-pi)*(2*pi*M[5])^(-1.5)+(2*pi*M[5])^-.5*exp(-(X-M[3])^2/(2*M[5]))*(X-M[3])^2/(2*M[5]^2)))#sigma2 portion
  )
  
  norm.grad <- sqrt(sum(grad^2))
  slope <- grad/norm.grad
  
  while (norm.grad> .001){ #while loop will run until normalized gradient is sufficiently small
    M <- M+slope*step
    Q <- M[1]*((2*pi*M[4])^-.5)*exp(-(X-M[2])^2/(2*M[4])) + (1-M[1])*((2*pi*M[5])^-.5)*exp(-(X-M[3])^2/(2*M[5]))
    grad <- c(sum((((2*pi*M[4])^-.5)*exp(-(X-M[2])^2/(2*M[4])) - ((2*pi*M[5])^-.5)*exp(-(X-M[3])^2/(2*M[5])))/Q), #p
              sum((M[1]*((2*pi*M[4])^-.5)*exp(-(X-M[2])^2/(2*M[4]))*2*(X-M[2]))/Q), #mu1
              sum((1-M[1])*((2*pi*M[5])^-.5)*exp(-(X-M[3])^2/(2*M[5]))*2*(X-M[3])/Q), #mu2
              sum((M[1]/Q)*(exp(-(X-M[2])/(2*M[4]))*(-pi)*(2*pi*M[4])^(-1.5)+(2*pi*M[4])^-.5*(X-M[2])^2*exp(-(X-M[2])^2/(2*M[4]))/(2*M[4]^2))),#sigma1
              sum(((1-M[1])/Q)*(exp(-(X-M[3])/(2*M[5]))*(-pi)*(2*pi*M[5])^(-1.5)+(2*pi*M[5])^-.5*exp(-(X-M[3])^2/(2*M[5]))*(X-M[3])^2/(2*M[5]^2))) #sigma2
    )
    norm.grad <- sqrt(sum(grad^2))
    slope <- grad/norm.grad
    print(M)
  }
  return(M)  
}