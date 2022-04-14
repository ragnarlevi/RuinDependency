# Parameters ----
# k is gamma shape parameter
# beta is gamma scale parameter
# lambda is claim frequency per exposure
# fixed cost is operation cost
# r should be 0

# Fix  surplus values
x_surplus <- c(2, 5, 10)

kendell <- c(0.05, 0.15,  0.4,  0.8)
theta_grid <- seq(from = 0.15, to = 0.4, by = 0.01)
theta_grid_ruin <- c(seq(from = 0.06, to = 0.28, by = 0.01), seq(from = 0.281, to = 0.31, by = 0.001), seq(from = 0.32, to = 0.44, by = 0.01)) 


# Define claims
N <- c(1, 1)
k <- c(1, 2)
beta <- c(1, 3/2)
lambda <- c(10000, 10000/3)
fixed_cost <-   c(1000/3, 1000/3)
r <- c(0,0) 

nu <- 1



#' @param b1 - logit parameter intercept
#' @param b2 - logit parameter loading slope
demand_1 <- function(theta){
  p <- 1/3
  tt <- 1/5
  b1 <- log((1-p)/p) - 1/(1-p)
  b2 <- 1/(tt*(1-p))
  
  1/(1+exp(b1+b2*theta))
}

#' @param b1 - logit parameter intercept
#' @param b2 - logit parameter loading slope
demand_2 <- function(theta){
  p <- 1/3
  tt <- 1/5
  b1 <- log((1-p)/p) - 1/(1-p)
  b2 <- 1/(tt*(1-p))
  
  1/(1+exp(b1+b2*theta))
}

demand <- function(theta){
  return(c(demand_1(theta), demand_2(theta)))
}

demand2var <- function(theta1, theta2){
  return(c(demand_1(theta1), demand_2(theta2)))
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