# Parameters ----
# k is gamma shape parameter
# beta is gamma scale parameter
# lambda is claim frequency per exposure
# fixed cost is operation cost
# r should be 0

# Fix some surplus values surpluses used 
x_surplus <- c(100, 5000, 15000, 20000) #, round(beta*k*lambda*N*demand(1), digits = -1))

# Define first claim

N <- c(10000, 10000)
k <- c(2, 2)
beta <- c(500, 500)
lambda <- c(0.08, 0.08)
fixed_cost <-   0.2*k*beta*N*lambda*0.4
r <- c(0,0) 

nu <- 1

S_ <- function(x, a, b){
  
  return(c(1-pgamma(q = x , shape = a[1], scale = b[1]), 1-pgamma(q = x , shape = a[2], scale = b[2])))
}


#' @param b1 - logit parameter intercept
#' @param b2 - logit parameter loading slope
demand_1 <- function(theta){
  b1 <- -0.6
  b2 <- 4
  
  1/(1+exp(b1+b2*theta))
}

#' @param b1 - logit parameter intercept
#' @param b2 - logit parameter loading slope
demand_2 <- function(theta){
  b1 <- -0.6
  b2 <- 4.5
  
  1/(1+exp(b1+b2*theta))
}

demand <- function(theta){
  b1 <- c(-0.6, 4)
  b2 <- c(-0.6, 4.5)
  return(c(1/(1+exp(b1[1]+b1[2]*theta)), 1/(1+exp(b2[1]+b2[2]*theta))))
  #return(c(1,1))
}


# Probability Inside integral
S_marginal <-  function(x, shape, scale) 1-pgamma(q = x,shape = shape, scale = scale)
S_ <- function(x, a, b){
  
  return(c(1-pgamma(q = x , shape = a[1], scale = b[1]), 1-pgamma(q = x , shape = a[2], scale = b[2])))
}


# Color codes ------
my_blue <- colorRampPalette(brewer.pal(n = 9, "Blues")[3:7])(8)[1:8]
my_red <- colorRampPalette(brewer.pal(n = 9, "Reds")[4:8])(8)[1:8]
my_grey <- brewer.pal(n = 9, "Greys")[9]