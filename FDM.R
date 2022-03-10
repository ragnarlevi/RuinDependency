#' library(tidyverse)
#' library(plotly)
#' 
#' ## Exponential claim correct
#' 
#' #'@param mu Rate paramter
#' ruin_exp <- function(theta, a, b, r, mu, surplus){
#'   
#' 
#'   lambda <- a- b*theta
#'   p <- (1+theta)*lambda/mu-r
#'   const <- lambda/(p*mu)
#'   R <- (mu-lambda/p)
#'   
#'   ruin <- const*exp(-R*surplus)
#'   
#'   
#'   return(ruin)
#'   
#' }
#' 
#' ruin_exp_logit <- function(theta, demand, D, l, mu, r, surplus){
#'   
#'   
#'   lambda <- D*demand(theta)*l
#'   p <- (1+theta)*lambda*mu-r
#'   
#'   const <- lambda/(p*mu)
#'   R <- (mu-lambda/p)
#'   
#'   ruin <- const*exp(-R*surplus)
#'   
#'   
#'   return(ruin)
#'   
#' }
#' 
#' 
#' ##### infinite time Ruin ####
#' 
#' 
#' 
#' #### Polynomial approximation #####
#' 
#' infinite_ruin_linear_interpolation <- function(to, h_x, a,b,r, mu, theta, f_, F_, F_bar, F_2bar, E_exp_trunc,...){
#'   
#'   lambda <- a-b*theta
#'   p <- lambda*(1+theta)*mu-r
#'   stopifnot(p>lambda*mu)
#'   
#'   
#'   # we expect v(t,x) -> 0 as x->infinity
#'   x_steps <- seq(from = 0,  to = to, by = h_x)
#'   N_x <- length(x_steps)
#'   
#'   
#'   # We subract a triangle with sides equal to N_t
#'   # as we expect the 
#'   V <- matrix(NA, nrow = N_x , ncol =1)
#'   
#'   V[1] <- lambda*mu/p
#'   
#'   ratio <- lambda/p
#'   
#'   for(i in 1:(N_x-1)){
#'     
#'     if(i == 1){
#'       
#'       # sum independent of V_i
#'       s_indp <- 0
#'       # constants infront of V_i
#'       V_i <- 0
#'       # Single integral
#'       s_indp <- s_indp + V[1]*h_x - (V[1]/h_x)*0.5*(h_x)^2
#'       # a_r <- + (V[1]*h_x - (V[1]/h_x)*0.5*(h_x)^2)*ratio
#'       V_i <- V_i - (1/h_x)*0.5*(h_x)^2
#'       # a_l <- (1/h_x)*0.5*(h_x)^2*ratio
#'       
#'       # add truncation
#'       s_indp <- s_indp - (E_exp_trunc(x = x_steps[2], ...)- E_exp_trunc(x = x_steps[1], ...))
#'       
#'       # double integration
#'       tmp1 <- V[2-1]*(F_bar(x_steps[2], ...) - F_bar(x_steps[2-1], ...)) - F_(x_steps[2-1], ...)*V[2-1]*h_x
#'       tmp2 <- -x_steps[1]*(F_bar(x_steps[2], ...)-F_bar(x_steps[1], ...)) -
#'         0.5*F_(x_steps[1], ...)*((x_steps[2]-2*x_steps[1])^2-(x_steps[1]-2*x_steps[1])^2) +
#'         F_2bar(x_steps[2],...) - F_2bar(x_steps[1],...) - F_bar(x_steps[1], ...)*h_x
#'       
#'       #d_r <- (tmp1 -(V[1]/h_x)*tmp2)*ratio
#'       #d_l <- (1/h_x)*tmp2*ratio
#'       
#'       #  (V[i-1] + a_r - b_r - d_r)/(1-a_l+d_l)
#'       
#'       
#'       # correct above  tmp1 + tmp2*(0.998-1)/h_x
#'       
#'       s_indp <- s_indp - (tmp1 - (V[1]/h_x)*tmp2)
#'       
#'       V_i <- V_i + tmp2/h_x
#'       
#'       # multiply with ratio
#'       
#'       s_indp <- s_indp*ratio + V[i-1+1]
#'       V_i <- V_i*ratio + 1 # plus one because the equation contaions V_i
#'       
#'       V[i+1] <- s_indp/V_i
#'       
#'       
#'       
#'     }else if(i == 2){
#'       
#'       
#'       # sum independent of V_i
#'       s_indp <- 0
#'       # constants infront of V_i
#'       V_i <- 0
#'       
#'       # Single integral
#'       s_indp <- s_indp + V[i-1+1]*h_x - (V[i-1+1]/h_x)*0.5*(h_x)^2
#'       V_i <- V_i - (1/h_x)*0.5*(h_x)^2
#'       
#'       # add truncation
#'       s_indp <- s_indp -(E_exp_trunc(x = x_steps[i+1], ...)- E_exp_trunc(x = x_steps[i-1+1], ...))
#'       
#'       # double integration x0 t0 x1
#'       tmp1 <- V[i-1+1]*(F_(x_steps[2], ...) - F_(x_steps[1], ...))*h_x
#'       tmp2 <- 0.5*F_(x_steps[2], ...)*((x_steps[3]-x_steps[2]-x_steps[2])^2-(x_steps[2]-x_steps[2]-x_steps[2])^2) - 
#'         0.5*F_(x_steps[1], ...)*((x_steps[3]-x_steps[1]-x_steps[2])^2-(x_steps[2]-x_steps[1]-x_steps[2])^2) + 
#'         (F_bar(x_steps[2], ...)-F_bar(x_steps[1], ...))*h_x
#'       
#'       s_indp <- s_indp - tmp1 + tmp2*V[2]/h_x
#'       V_i <- V_i + tmp2/h_x
#'       # correct above
#'       
#'       
#'       # double integration x1 t0 z
#'       tmp1 <- V[1]*(F_bar(x_steps[3], ...) - F_bar(x_steps[2], ...)) - V[1]*F_(x_steps[2], ...)*h_x + 
#'         ((V[2]-V[1])/h_x)*(-x_steps[1]*(F_bar(x_steps[3], ...) - F_bar(x_steps[2], ...)) -
#'                              0.5*F_(x_steps[2], ...)*((x_steps[3] - x_steps[2] -x_steps[1])^2-(x_steps[2] - x_steps[2] -x_steps[1])^2) +
#'                              F_2bar(x_steps[3], ...) - F_2bar(x_steps[2], ...) -F_bar(x_steps[2], ...)*h_x
#'                            
#'         ) 
#'       
#'       # correct above 
#'       
#'       s_indp <- s_indp - tmp1
#'       
#'       
#'       # multiply with ratio
#'       
#'       s_indp <- s_indp*ratio + V[i-1+1]
#'       V_i <- V_i*ratio + 1 # plus one because the equation contaions V_i
#'       
#'       V[i+1] <- s_indp/V_i
#'       
#'     }else{
#'       
#'       # sum independent of V_i
#'       s_indp <- 0
#'       # constants infront of V_i
#'       V_i <- 0
#'       
#'       # Single integral
#'       s_indp <- s_indp + V[i-1+1]*h_x - (V[i-1+1]/h_x)*0.5*(h_x)^2
#'       V_i <- V_i - (1/h_x)*0.5*(h_x)^2
#'       
#'       # add truncation
#'       s_indp <- s_indp -(E_exp_trunc(x = x_steps[i+1], ...)- E_exp_trunc(x = x_steps[i-1+1], ...))
#'       
#'       # double integration x0 t0 x1
#'       tmp1 <- V[i-1+1]*(F_(x_steps[2], ...) - F_(x_steps[1], ...))*h_x
#'       tmp2 <- 0.5*F_(x_steps[2], ...)*((x_steps[i+1]-x_steps[2]-x_steps[i-1+1])^2-(x_steps[i-1+1]-x_steps[2]-x_steps[i-1+1])^2) - 
#'         0.5*F_(x_steps[1], ...)*((x_steps[i+1]-x_steps[1]-x_steps[i-1+1])^2-(x_steps[i-1+1]-x_steps[1]-x_steps[i-1+1])^2) + 
#'         (F_bar(x_steps[2], ...)-F_bar(x_steps[1], ...))*h_x
#'       
#'       s_indp <- s_indp - tmp1 + tmp2*(V[i-1+1]/h_x)
#'       V_i <- V_i +tmp2/h_x
#'       # correct above
#'       
#'       #double integration from xj-1 to xj
#'       
#'       for(j in 2:(i-1)){
#'         
#'         tmp1 <- V[i-j+1]*(F_(x_steps[j+1], ...)-F_(x_steps[j], ...))*h_x + ((V[i-j+1+1]-V[i-j+1])/h_x)*(
#'           0.5*F_(x_steps[j+1], ...)*((x_steps[i+1]-x_steps[j+1]-x_steps[i-j+1])^2-(x_steps[i-1+1]-x_steps[j+1]-x_steps[i-j+1])^2)-
#'             0.5*F_(x_steps[j], ...)*((x_steps[i+1]-x_steps[j]-x_steps[i-j+1])^2-(x_steps[i-1+1]-x_steps[j]-x_steps[i-j+1])^2) +
#'             (F_bar(x_steps[j+1], ...)-F_bar(x_steps[j], ...))*h_x
#'         )
#'         
#'         s_indp <- s_indp - tmp1 
#'       }
#'       
#'       
#'       # double integration x1 t0 z
#'       tmp1 <- V[1]*(F_bar(x_steps[i+1], ...) - F_bar(x_steps[i-1+1], ...)) - V[1]*F_(x_steps[i-1+1], ...)*h_x + 
#'         ((V[2]-V[1])/h_x)*(-x_steps[1]*(F_bar(x_steps[i+1], ...) - F_bar(x_steps[i-1+1], ...)) -
#'                              0.5*F_(x_steps[i-1+1], ...)*((x_steps[i+1] - x_steps[i-1+1] -x_steps[1])^2-(x_steps[i-1+1] - x_steps[i-1+1] -x_steps[1])^2) +
#'                              F_2bar(x_steps[i+1], ...) - F_2bar(x_steps[i-1+1], ...) - F_bar(x_steps[i-1+1], ...)*h_x
#'                            
#'         ) 
#'       
#'       # correct above 
#'       
#'       s_indp <- s_indp - tmp1
#'       
#'       
#'       # multiply with ratio
#'       
#'       s_indp <- s_indp*ratio + V[i-1+1]
#'       V_i <- V_i*ratio + 1 # plus one because the equation contaions V_i
#'       
#'       V[i+1] <- s_indp/V_i
#'       
#'     }
#'     
#'     
#'   }
#'   
#'   
#'   ret <- list()
#'   ret$V <- V
#'   ret$x <- x_steps
#'   return(ret)
#'   
#' }
#' 
#' infinite_ruin_linear_interpolation2 <- function(to, h_x, a,b,r, mu, theta, f_, ...){
#'   
#'   lambda <- a-b*theta
#'   p <- lambda*(1+theta)*mu-r
#'   stopifnot(p>lambda*mu)
#'   
#'   F_ <- function(x, ...){integrate(f_, lower = 0, upper = x, ...)$value}
#'   F_bar <- function(x,...){
#'     fun <- Vectorize(F_, vectorize.args = "x")
#'     integrate(fun, lower = 0, upper = x, ...)$value
#'   }
#'   F_2bar <- function(x,...){
#'     fun <- Vectorize(F_bar, vectorize.args = "x")
#'     integrate(fun, lower = 0, upper = x, ...)$value
#'   }
#'   E <- function(x1, x2, ...){
#'     fun <- function(x,...){1-F_(x,...)}
#'     fun <- Vectorize(fun, vectorize.args = "x")
#'     integrate(fun, lower = x1, upper = x2, ...)$value
#'   }
#' 
#'   # F_t(2, rate = 1)
#'   # F_(2, rate = 1)
#'   # 
#'   # 
#'   # F_bart(2, rate = 1)
#'   # F_bar(2, rate = 1)
#'   # 
#'   # F_2bart(2, rate = 1)
#'   # F_2bar(2, rate = 1)
#'   # 
#'   # Et(0.1, 0.2, rate = 1)
#'   # E_exp_trunc(0.2, rate = 1) - E_exp_trunc(0.1, rate = 1)
#'   
#'   
#'   
#'   # we expect v(t,x) -> 0 as x->infinity
#'   x_steps <- seq(from = 0,  to = to, by = h_x)
#'   N_x <- length(x_steps)
#'   
#'   
#'   # We subract a triangle with sides equal to N_t
#'   # as we expect the 
#'   V <- matrix(NA, nrow = N_x , ncol =1)
#'   
#'   V[1] <- lambda*mu/p
#'   
#'   ratio <- lambda/p
#'   
#'   for(i in 1:(N_x-1)){
#'     
#'     if(i == 1){
#'       
#'       # sum independent of V_i
#'       s_indp <- 0
#'       # constants infront of V_i
#'       V_i <- 0
#'       # Single integral
#'       s_indp <- s_indp + V[1]*h_x - (V[1]/h_x)*0.5*(h_x)^2
#'       # a_r <- + (V[1]*h_x - (V[1]/h_x)*0.5*(h_x)^2)*ratio
#'       V_i <- V_i - (1/h_x)*0.5*(h_x)^2
#'       # a_l <- (1/h_x)*0.5*(h_x)^2*ratio
#'       
#'       # add truncation
#'       s_indp <- s_indp - E(x1 = x_steps[1], x2 = x_steps[2],...)
#'       
#'       # double integration
#'       tmp1 <- V[2-1]*(F_bar(x_steps[2], ...) - F_bar(x_steps[2-1], ...)) - F_(x_steps[2-1], ...)*V[2-1]*h_x
#'       tmp2 <- -x_steps[1]*(F_bar(x_steps[2], ...)-F_bar(x_steps[1], ...)) -
#'         0.5*F_(x_steps[1], ...)*((x_steps[2]-2*x_steps[1])^2-(x_steps[1]-2*x_steps[1])^2) +
#'         F_2bar(x_steps[2],...) - F_2bar(x_steps[1],...) - F_bar(x_steps[1], ...)*h_x
#'       
#'       # correct above  tmp1 + tmp2*(0.998-1)/h_x
#'       
#'       s_indp <- s_indp - (tmp1 - (V[1]/h_x)*tmp2)
#'       
#'       V_i <- V_i + tmp2/h_x
#'       
#'       # multiply with ratio
#'       
#'       s_indp <- s_indp*ratio + V[i-1+1]
#'       V_i <- V_i*ratio + 1 # plus one because the equation contaions V_i
#'       
#'       V[i+1] <- s_indp/V_i
#'       
#'       
#'       
#'     }else if(i == 2){
#'       
#'       
#'       # sum independent of V_i
#'       s_indp <- 0
#'       # constants infront of V_i
#'       V_i <- 0
#'       
#'       # Single integral
#'       s_indp <- s_indp + V[i-1+1]*h_x - (V[i-1+1]/h_x)*0.5*(h_x)^2
#'       V_i <- V_i - (1/h_x)*0.5*(h_x)^2
#'       
#'       # add truncation
#'       s_indp <- s_indp - E(x1 = x_steps[i-1+1], x2 = x_steps[i+1],...)
#'       
#'       # double integration x0 t0 x1
#'       tmp1 <- V[i-1+1]*(F_(x_steps[2], ...) - F_(x_steps[1], ...))*h_x
#'       tmp2 <- 0.5*F_(x_steps[2], ...)*((x_steps[3]-x_steps[2]-x_steps[2])^2-(x_steps[2]-x_steps[2]-x_steps[2])^2) - 
#'         0.5*F_(x_steps[1], ...)*((x_steps[3]-x_steps[1]-x_steps[2])^2-(x_steps[2]-x_steps[1]-x_steps[2])^2) + 
#'         (F_bar(x_steps[2], ...)-F_bar(x_steps[1], ...))*h_x
#'       
#'       s_indp <- s_indp - tmp1 + tmp2*V[2]/h_x
#'       V_i <- V_i + tmp2/h_x
#'       # correct above
#'       
#'       
#'       # double integration x1 t0 z
#'       tmp1 <- V[1]*(F_bar(x_steps[3], ...) - F_bar(x_steps[2], ...)) - V[1]*F_(x_steps[2], ...)*h_x + 
#'         ((V[2]-V[1])/h_x)*(-x_steps[1]*(F_bar(x_steps[3], ...) - F_bar(x_steps[2], ...)) -
#'                              0.5*F_(x_steps[2], ...)*((x_steps[3] - x_steps[2] -x_steps[1])^2-(x_steps[2] - x_steps[2] -x_steps[1])^2) +
#'                              F_2bar(x_steps[3], ...) - F_2bar(x_steps[2], ...) -F_bar(x_steps[2], ...)*h_x
#'                            
#'         ) 
#'       
#'       # correct above 
#'       
#'       s_indp <- s_indp - tmp1
#'       
#'       
#'       # multiply with ratio
#'       
#'       s_indp <- s_indp*ratio + V[i-1+1]
#'       V_i <- V_i*ratio + 1 # plus one because the equation contaions V_i
#'       
#'       V[i+1] <- s_indp/V_i
#'       
#'     }else{
#'       
#'       # sum independent of V_i
#'       s_indp <- 0
#'       # constants infront of V_i
#'       V_i <- 0
#'       
#'       # Single integral
#'       s_indp <- s_indp + V[i-1+1]*h_x - (V[i-1+1]/h_x)*0.5*(h_x)^2
#'       V_i <- V_i - (1/h_x)*0.5*(h_x)^2
#'       
#'       # add truncation
#'       s_indp <- s_indp - E(x1 = x_steps[i-1+1], x2 = x_steps[i+1],...)
#'       
#'       # double integration x0 t0 x1
#'       tmp1 <- V[i-1+1]*(F_(x_steps[2], ...) - F_(x_steps[1], ...))*h_x
#'       tmp2 <- 0.5*F_(x_steps[2], ...)*((x_steps[i+1]-x_steps[2]-x_steps[i-1+1])^2-(x_steps[i-1+1]-x_steps[2]-x_steps[i-1+1])^2) - 
#'         0.5*F_(x_steps[1], ...)*((x_steps[i+1]-x_steps[1]-x_steps[i-1+1])^2-(x_steps[i-1+1]-x_steps[1]-x_steps[i-1+1])^2) + 
#'         (F_bar(x_steps[2], ...)-F_bar(x_steps[1], ...))*h_x
#'       
#'       s_indp <- s_indp - tmp1 + tmp2*(V[i-1+1]/h_x)
#'       V_i <- V_i +tmp2/h_x
#'       # correct above
#'       
#'       #double integration from xj-1 to xj
#'       
#'       for(j in 2:(i-1)){
#'         
#'         tmp1 <- V[i-j+1]*(F_(x_steps[j+1], ...)-F_(x_steps[j], ...))*h_x + ((V[i-j+1+1]-V[i-j+1])/h_x)*(
#'           0.5*F_(x_steps[j+1], ...)*((x_steps[i+1]-x_steps[j+1]-x_steps[i-j+1])^2-(x_steps[i-1+1]-x_steps[j+1]-x_steps[i-j+1])^2)-
#'             0.5*F_(x_steps[j], ...)*((x_steps[i+1]-x_steps[j]-x_steps[i-j+1])^2-(x_steps[i-1+1]-x_steps[j]-x_steps[i-j+1])^2) +
#'             (F_bar(x_steps[j+1], ...)-F_bar(x_steps[j], ...))*h_x
#'         )
#'         
#'         s_indp <- s_indp - tmp1 
#'       }
#'       
#'       
#'       # double integration x1 t0 z
#'       tmp1 <- V[1]*(F_bar(x_steps[i+1], ...) - F_bar(x_steps[i-1+1], ...)) - V[1]*F_(x_steps[i-1+1], ...)*h_x + 
#'         ((V[2]-V[1])/h_x)*(-x_steps[1]*(F_bar(x_steps[i+1], ...) - F_bar(x_steps[i-1+1], ...)) -
#'                              0.5*F_(x_steps[i-1+1], ...)*((x_steps[i+1] - x_steps[i-1+1] -x_steps[1])^2-(x_steps[i-1+1] - x_steps[i-1+1] -x_steps[1])^2) +
#'                              F_2bar(x_steps[i+1], ...) - F_2bar(x_steps[i-1+1], ...) - F_bar(x_steps[i-1+1], ...)*h_x
#'                            
#'         ) 
#'       
#'       # correct above 
#'       
#'       s_indp <- s_indp - tmp1
#'       
#'       
#'       # multiply with ratio
#'       
#'       s_indp <- s_indp*ratio + V[i-1+1]
#'       V_i <- V_i*ratio + 1 # plus one because the equation contaions V_i
#'       
#'       V[i+1] <- s_indp/V_i
#'       
#'     }
#'     
#'     
#'   }
#'   
#'   
#'   ret <- list()
#'   ret$V <- V
#'   ret$x <- x_steps
#'   return(ret)
#'   
#' }
#' 
#' 
#' 
#' 
#' # Infnite Survival ----
#' 
#' infnite_suvival <- function(to, h_x, theta, a, b, r, mu, S_bar, S_2bar, ...){
#'   
#'   lambda <- a-b*theta
#'   p <- (1+theta)*lambda*mu-r
#'   stopifnot(p>lambda*mu)
#'   # we expect v(t,x) -> 0 as x->infinity
#'   x_steps <- seq(from = 0,to = to, by = h_x)
#'   N_x <- length(x_steps)
#'   
#'   
#'   V_bar <- matrix(NA, nrow = N_x , ncol =1)
#'   
#'   # Do survival
#'   V_bar[1] <- 1-lambda*mu/p
#'  # V_bar[2] <- 1-1/(1+theta)
#'   
#'   ratio <- lambda/p
#'  # V_bar[3] <- V_bar[2]+ ratio*(V_bar[2]*S_(x_steps[3], ...)*h_x/2+V_bar[2]*S_(x_steps[2], ...)*h_x/2)
#'   for(i in 1:(N_x-1)){
#'     if(i==1){
#'       
#'       tmp1 <- V_bar[i-1+1]*(S_bar(x_steps[i+1], ...)-S_bar(x_steps[i+1-1], ...))
#'       tmp2 <- (x_steps[i+1] - x_steps[i+1] -x_steps[i+1-1])*S_bar(x_steps[i+1], ...) - (x_steps[i+1]-2*x_steps[i+1-1])*S_bar(x_steps[i+1-1],...) +
#'         S_2bar(x_steps[i+1], ...)-S_2bar(x_steps[i+1-1], ...)
#'       
#'       V_bar[1+1] <- (V_bar[0+1]+ratio*(tmp1-tmp2*V_bar[1-1+1]/h_x))/(1-ratio*tmp2/h_x)
#'       
#'     }else{
#'       
#'       s <- 0
#'       for(j in 2:i){
#'         s <- s + V_bar[i-j+1]*(S_bar(x_steps[j+1], ...)-S_bar(x_steps[j-1+1], ...)) + 
#'           ((V_bar[i-j+1+1]-V_bar[i-j+1])/h_x)*((x_steps[i+1] - x_steps[j+1] - x_steps[i-j+1])*S_bar(x_steps[j+1], ...)-
#'                                                (x_steps[i+1] - x_steps[j-1+1] - x_steps[i-j+1])*S_bar(x_steps[j-1+1],...) + 
#'                                                 S_2bar(x_steps[j+1], ...)-S_2bar(x_steps[j-1+1], ...))
#'       }
#'       
#'       tmp1 <- V_bar[i-1+1]*(S_bar(x_steps[1+1], ...)-S_bar(x_steps[0+1], ...))
#'       tmp2 <- (x_steps[i+1] - x_steps[1+1] -x_steps[i-1+1])*S_bar(x_steps[1+1], ...) - (x_steps[i+1]-x_steps[0+1]-x_steps[i-1+1])*S_bar(x_steps[0+1],...) + 
#'         S_2bar(x_steps[1+1], ...)-S_2bar(x_steps[0+1], ...)
#'       
#'       V_bar[i+1] <- (V_bar[1]+ratio*s+ratio*(tmp1-tmp2*V_bar[i-1+1]/h_x))/(1-ratio*tmp2/h_x)
#'       
#'       
#'     }
#'     
#'   }
#' 
#'   ret <- list()
#'   ret$V_bar <- V_bar
#'   ret$V <- 1-V_bar
#'   ret$x <- x_steps
#'   return(ret)
#' }
#' 
#' #' @param lambda_true is the actual intensity of the compound Poisson
#' #' @param mu mean of the processes
#' infnite_suvival2 <- function(to, h_x, theta, a, b, r, mu, f_, lower = 0, lambda_true = NA, mean_true = NA, ...){
#' 
#'   lambda <- sum(a-b*theta)
#'   p <- sum((1+theta)*lambda*mu-r)
#'   # We a surplus horzion as computers are finite, we expect v(t,x) -> 0 as x->infinity
#'   x_steps <- seq(from = 0,to = to, by = h_x)
#'   N_x <- length(x_steps)
#' 
#'   
#'   S_ <- function(x,...){
#'     1-integrate(f_, lower = lower, upper = x, ...)$value
#'   }
#'   S_ <- Vectorize(S_, vectorize.args = "x")
#'   
#'   S_bar <- function(x,...){
#'     integrate(fun, lower = lower, upper = x, ...)$value
#'   }
#'   S_bar <- Vectorize(S_bar, vectorize.args = "x")
#'   
#'   S_2bar <- function(x,...){
#'     integrate(fun, lower = lower, upper = x, ...)$value
#'   }
#'   S_2bar <- Vectorize(S_bar, vectorize.args = "x")
#'   
#'   # S_t(2, 1)
#'   # S_(2,1)
#'   # 
#'   # S_bar(2, rate = 1)
#'   # S_bart(2, rate = 1)
#'   # 
#'   # S_2bar(2,1)
#'   # S_2bart(2,1)
#'   
#'   
#'   V_bar <- matrix(NA, nrow = N_x , ncol =1)
#'   
#'   if(is.na(mean_true)){
#'     fun <- Vectorize(S_, vectorize.args = "x")
#'     mean_true <- integrate(fun, lower = 0, upper = 30, ...)$value
#'     
#'     print(mean_true)
#'   }
#'   
#'   # Do survival
#'   
#'   if(is.na(lambda_true)){
#'     lambda_true <- lambda
#'   }
#'   V_bar[1] <- 1-lambda_true*mean_true/p
#'   # V_bar[2] <- 1-1/(1+theta)
#'   
#'   ratio <- lambda_true/p
#'   # V_bar[3] <- V_bar[2]+ ratio*(V_bar[2]*S_(x_steps[3], ...)*h_x/2+V_bar[2]*S_(x_steps[2], ...)*h_x/2)
#'   for(i in 1:(N_x-1)){
#'     print(paste0(i, " of ", N_x, sep = ""))
#'     if(i==1){
#'       
#'       tmp1 <- V_bar[i-1+1]*(S_bar(x_steps[i+1], ...)-S_bar(x_steps[i+1-1], ...))
#'       tmp2 <- (x_steps[i+1] - x_steps[i+1] -x_steps[i+1-1])*S_bar(x_steps[i+1], ...) - (x_steps[i+1]-2*x_steps[i+1-1])*S_bar(x_steps[i+1-1],...) + 
#'         S_2bar(x_steps[i+1], ...)-S_2bar(x_steps[i+1-1], ...)
#'       
#'       V_bar[1+1] <- (V_bar[0+1]+ratio*(tmp1-tmp2*V_bar[1-1+1]/h_x))/(1-ratio*tmp2/h_x)
#'       
#'     }else{
#'       
#'       s <- 0
#'       for(j in 2:i){
#'         s <- s + V_bar[i-j+1]*(S_bar(x_steps[j+1], ...)-S_bar(x_steps[j-1+1], ...)) + 
#'           ((V_bar[i-j+1+1]-V_bar[i-j+1])/h_x)*((x_steps[i+1] - x_steps[j+1] - x_steps[i-j+1])*S_bar(x_steps[j+1], ...)-
#'                                                  (x_steps[i+1] - x_steps[j-1+1] - x_steps[i-j+1])*S_bar(x_steps[j-1+1],...) + 
#'                                                  S_2bar(x_steps[j+1], ...)-S_2bar(x_steps[j-1+1], ...))
#'       }
#'       
#'       tmp1 <- V_bar[i-1+1]*(S_bar(x_steps[1+1], ...)-S_bar(x_steps[0+1], ...))
#'       tmp2 <- (x_steps[i+1] - x_steps[1+1] -x_steps[i-1+1])*S_bar(x_steps[1+1], ...) - (x_steps[i+1]-x_steps[0+1]-x_steps[i-1+1])*S_bar(x_steps[0+1],...) + 
#'         S_2bar(x_steps[1+1], ...)-S_2bar(x_steps[0+1], ...)
#'       
#'       V_bar[i+1] <- (V_bar[1]+ratio*s+ratio*(tmp1-tmp2*V_bar[i-1+1]/h_x))/(1-ratio*tmp2/h_x)
#'       
#'       
#'     }
#'     
#'   }
#'   
#'   ret <- list()
#'   ret$V_bar <- V_bar
#'   ret$V <- 1-V_bar
#'   ret$x <- x_steps
#'   return(ret)
#' }
#' 
#' # Try to make 2 a lot more efficient
infnite_survival3 <- function(to, h_x, p, lambda_true, mean_true, f_ = "", S_ = "", print_time = FALSE, lower = 0, ...){


  # We a surplus horzion as computers are finite, we expect v(t,x) -> 0 as x->infinity
  x_steps <- seq(from = 0,to = to, by = h_x)
  N_x <- length(x_steps)
  V_bar <- matrix(NA, nrow = N_x , ncol =1)


  # Do survival


  V_bar[1] <- 1-lambda_true*mean_true/p

  if(V_bar[1]>=1 | V_bar[1]<= 0){
    warning(paste0("V_bar is: ", V_bar[1]))
  }

  if(is.character(f_) & is.character(S_)){
    stop("Either f_ or S_has to be specified")
  }

  if(is.character(S_)){
    S_ <- function(x,...){
      1-integrate(f_, lower = lower, upper = x, ...)$value
    }
    S_ <- Vectorize(S_, vectorize.args = "x")
  }


  S_bar <- function(x,...){
    integrate(S_, lower = lower, upper = x, ...)$value
  }
  S_bar <- Vectorize(S_bar, vectorize.args = "x")

  S_2bar <- function(x,...){
    integrate(S_bar, lower = lower, upper = x, ...)$value
  }
  S_2bar <- Vectorize(S_2bar, vectorize.args = "x")

  ttt <- Sys.time()
 # print(ttt)
  # Calculate all values before looping
  S_bar_v <- S_bar(x_steps, ...)
  S_2bar_v <- S_2bar(x_steps, ...)

  if(print_time)print(Sys.time() - ttt)




  # V_bar[2] <- 1-1/(1+theta)

  ratio <- lambda_true/p
  # V_bar[3] <- V_bar[2]+ ratio*(V_bar[2]*S_(x_steps[3], ...)*h_x/2+V_bar[2]*S_(x_steps[2], ...)*h_x/2)
  for(i in 1:(N_x-1)){
    # print(paste0(i, " of ", N_x, sep = ""))

    if(i==1){

      tmp1 <- V_bar[i-1+1]*(S_bar_v[i+1]-S_bar_v[i+1-1])
      tmp2 <- (x_steps[i+1] - x_steps[i+1] -x_steps[i+1-1])*S_bar_v[i+1] - (x_steps[i+1]-2*x_steps[i+1-1])*S_bar_v[i+1-1] +
        S_2bar_v[i+1] - S_2bar_v[i+1-1]

      V_bar[1+1] <- (V_bar[0+1]+ratio*(tmp1-tmp2*V_bar[1-1+1]/h_x))/(1-ratio*tmp2/h_x)

    }else{

      s <- 0
      for(j in 2:i){
        s <- s + V_bar[i-j+1]*(S_bar_v[j+1]-S_bar_v[j-1+1]) +
          ((V_bar[i-j+1+1]-V_bar[i-j+1])/h_x)*((x_steps[i+1] - x_steps[j+1] - x_steps[i-j+1])*S_bar_v[j+1]-
                                                 (x_steps[i+1] - x_steps[j-1+1] - x_steps[i-j+1])*S_bar_v[j-1+1] +
                                                 S_2bar_v[j+1] - S_2bar_v[j-1+1])
      }

      tmp1 <- V_bar[i-1+1]*(S_bar_v[1+1] - S_bar_v[0+1])
      tmp2 <- (x_steps[i+1] - x_steps[1+1] -x_steps[i-1+1])*S_bar_v[1+1] - (x_steps[i+1]-x_steps[0+1]-x_steps[i-1+1])*S_bar_v[0+1] +
        S_2bar_v[1+1] - S_2bar_v[0+1]

      V_bar[i+1] <- (V_bar[1]+ratio*s+ratio*(tmp1-tmp2*V_bar[i-1+1]/h_x))/(1-ratio*tmp2/h_x)


    }

  }

  ret <- list()
  ret$V_bar <- V_bar
  ret$V <- 1-V_bar
  ret$x <- x_steps
  return(ret)
}
#' 
#' 
#' # Instead of taking in a, b,r etc.., Take p and lambda
#' infnite_suvival4 <- function(to, h_x, p, lambda, mu, S_bar, S_2bar, ...){
#'   
#'   stopifnot(p>lambda*mu)
#'   # we expect v(t,x) -> 0 as x->infinity
#'   x_steps <- seq(from = 0,to = to, by = h_x)
#'   N_x <- length(x_steps)
#'   
#'   
#'   V_bar <- matrix(NA, nrow = N_x , ncol =1)
#'   
#'   # Do survival
#'   V_bar[1] <- 1-lambda*mu/p
#'   # V_bar[2] <- 1-1/(1+theta)
#'   
#'   ratio <- lambda/p
#'   # V_bar[3] <- V_bar[2]+ ratio*(V_bar[2]*S_(x_steps[3], ...)*h_x/2+V_bar[2]*S_(x_steps[2], ...)*h_x/2)
#'   for(i in 1:(N_x-1)){
#'     if(i==1){
#'       
#'       tmp1 <- V_bar[i-1+1]*(S_bar(x_steps[i+1], ...)-S_bar(x_steps[i+1-1], ...))
#'       tmp2 <- (x_steps[i+1] - x_steps[i+1] -x_steps[i+1-1])*S_bar(x_steps[i+1], ...) - (x_steps[i+1]-2*x_steps[i+1-1])*S_bar(x_steps[i+1-1],...) +
#'         S_2bar(x_steps[i+1], ...)-S_2bar(x_steps[i+1-1], ...)
#'       
#'       V_bar[1+1] <- (V_bar[0+1]+ratio*(tmp1-tmp2*V_bar[1-1+1]/h_x))/(1-ratio*tmp2/h_x)
#'       
#'     }else{
#'       
#'       s <- 0
#'       for(j in 2:i){
#'         s <- s + V_bar[i-j+1]*(S_bar(x_steps[j+1], ...)-S_bar(x_steps[j-1+1], ...)) + 
#'           ((V_bar[i-j+1+1]-V_bar[i-j+1])/h_x)*((x_steps[i+1] - x_steps[j+1] - x_steps[i-j+1])*S_bar(x_steps[j+1], ...)-
#'                                                  (x_steps[i+1] - x_steps[j-1+1] - x_steps[i-j+1])*S_bar(x_steps[j-1+1],...) + 
#'                                                  S_2bar(x_steps[j+1], ...)-S_2bar(x_steps[j-1+1], ...))
#'       }
#'       
#'       tmp1 <- V_bar[i-1+1]*(S_bar(x_steps[1+1], ...)-S_bar(x_steps[0+1], ...))
#'       tmp2 <- (x_steps[i+1] - x_steps[1+1] -x_steps[i-1+1])*S_bar(x_steps[1+1], ...) - (x_steps[i+1]-x_steps[0+1]-x_steps[i-1+1])*S_bar(x_steps[0+1],...) + 
#'         S_2bar(x_steps[1+1], ...)-S_2bar(x_steps[0+1], ...)
#'       
#'       V_bar[i+1] <- (V_bar[1]+ratio*s+ratio*(tmp1-tmp2*V_bar[i-1+1]/h_x))/(1-ratio*tmp2/h_x)
#'       
#'       
#'     }
#'     
#'   }
#'   
#'   ret <- list()
#'   ret$V_bar <- V_bar
#'   ret$V <- 1-V_bar
#'   ret$x <- x_steps
#'   return(ret)
#' }
#'  
#' infnite_suvival_copula <- function(to, h_x, theta, a, b, r, mu, f_, lower = 0, ...){
#'   
#'   lambda <- sum(a-b*theta)
#'   p <- sum((1+theta)*lambda*mu-r)
#'   # We a surplus horzion as computers are finite, we expect v(t,x) -> 0 as x->infinity
#'   x_steps <- seq(from = 0,to = to, by = h_x)
#'   N_x <- length(x_steps)
#'   
#'   
#'   S_ <- function(x,...){
#'     1-integrate(f_, lower = lower, upper = x, ...)$value
#'   }
#'   S_bar <- function(x,...){
#'     fun <- Vectorize(S_, vectorize.args = "x")
#'     integrate(fun, lower = lower, upper = x, ...)$value
#'   }
#'   S_2bar <- function(x,...){
#'     fun <- Vectorize(S_bar, vectorize.args = "x")
#'     integrate(fun, lower = lower, upper = x, ...)$value
#'   }
#'   
#'   # S_t(2, 1)
#'   # S_(2,1)
#'   # 
#'   # S_bar(2, rate = 1)
#'   # S_bart(2, rate = 1)
#'   # 
#'   # S_2bar(2,1)
#'   # S_2bart(2,1)
#'   
#'   
#'   V_bar <- matrix(NA, nrow = N_x , ncol =1)
#'   
#'   # Do survival
#'   V_bar[1] <- 1-lambda*mu/p
#'   # V_bar[2] <- 1-1/(1+theta)
#'   
#'   ratio <- lambda/p
#'   # V_bar[3] <- V_bar[2]+ ratio*(V_bar[2]*S_(x_steps[3], ...)*h_x/2+V_bar[2]*S_(x_steps[2], ...)*h_x/2)
#'   for(i in 1:(N_x-1)){
#'     if(i==1){
#'       
#'       tmp1 <- V_bar[i-1+1]*(S_bar(x_steps[i+1], ...)-S_bar(x_steps[i+1-1], ...))
#'       tmp2 <- (x_steps[i+1] - x_steps[i+1] -x_steps[i+1-1])*S_bar(x_steps[i+1], ...) - (x_steps[i+1]-2*x_steps[i+1-1])*S_bar(x_steps[i+1-1],...) + 
#'         S_2bar(x_steps[i+1], ...)-S_2bar(x_steps[i+1-1], ...)
#'       
#'       V_bar[1+1] <- (V_bar[0+1]+ratio*(tmp1-tmp2*V_bar[1-1+1]/h_x))/(1-ratio*tmp2/h_x)
#'       
#'     }else{
#'       
#'       s <- 0
#'       for(j in 2:i){
#'         s <- s + V_bar[i-j+1]*(S_bar(x_steps[j+1], ...)-S_bar(x_steps[j-1+1], ...)) + 
#'           ((V_bar[i-j+1+1]-V_bar[i-j+1])/h_x)*((x_steps[i+1] - x_steps[j+1] - x_steps[i-j+1])*S_bar(x_steps[j+1], ...)-
#'                                                  (x_steps[i+1] - x_steps[j-1+1] - x_steps[i-j+1])*S_bar(x_steps[j-1+1],...) + 
#'                                                  S_2bar(x_steps[j+1], ...)-S_2bar(x_steps[j-1+1], ...))
#'       }
#'       
#'       tmp1 <- V_bar[i-1+1]*(S_bar(x_steps[1+1], ...)-S_bar(x_steps[0+1], ...))
#'       tmp2 <- (x_steps[i+1] - x_steps[1+1] -x_steps[i-1+1])*S_bar(x_steps[1+1], ...) - (x_steps[i+1]-x_steps[0+1]-x_steps[i-1+1])*S_bar(x_steps[0+1],...) + 
#'         S_2bar(x_steps[1+1], ...)-S_2bar(x_steps[0+1], ...)
#'       
#'       V_bar[i+1] <- (V_bar[1]+ratio*s+ratio*(tmp1-tmp2*V_bar[i-1+1]/h_x))/(1-ratio*tmp2/h_x)
#'       
#'       
#'     }
#'     
#'   }
#'   
#'   ret <- list()
#'   ret$V_bar <- V_bar
#'   ret$V <- 1-V_bar
#'   ret$x <- x_steps
#'   return(ret)
#' }
#' 
#' ##### Probability method ----
#' 
#' # @h - steps of x considered
#' # @t - starting time
#' # @Time - final time
#' # @x - final surplus considered
#' # @theta - loading
#' # @dist - distribution function of claims
#' # @mu - mean of the claim frequency
#' # @a - claim frequency if theta = 0
#' # @b - for every unit of theta lambda decreases by b
#' # @r - cost op operation
#' # @... - additional argument into dist
#'   Finite_Ruin_Probability <- function(h, t, Time, x, theta, dist, mu,a,b,r, V_T = 0, ... ){
#'   
#'   
#'   
#'   # h = T/K  
#'   t_steps <- seq(from = t, to = Time, by = h)
#'   
#'   lambda <- max(0,a-b*theta)
#'   p <- (1+theta)*lambda*mu-r
#'   
#'   x_steps <- seq(from = 0, to = x, by = p*h)
#'   
#'   K <- length(t_steps) - 1
#'   M <- length(x_steps) - 1
#'   
#'   V <- matrix(NA, nrow = K+1, ncol = M+1)
#'   V[K+1,] <- V_T
#' 
#'   
#'   # Note R starts indexing at 1
#'   for(i in seq(from = K-1, to = 0, by = -1)){
#'     for(j in 0:M){
#'       if(j< M && i <= K-1){
#'         V[i+1,j+1] <- (1-pexp(h, rate = lambda))*V[i+1+1, j+1+1] + 
#'           pexp(h, rate = lambda)*(sum((dist(p*h*(1:(j+1)), ...)-dist(p*h*(0:j), ...))*V[i+2, rev(0:j)+1]) 
#'                                   + (1-dist((j+1)*p*h, ...)))
#'       }else if(j == M && i <= K-1){
#'         # Assume V[i+1+1, j+1+1] = 0
#'         V[i+1,j+1] <- (1-pexp(h, rate = lambda))*0 + 
#'           pexp(h, rate = lambda)*(sum((dist(p*h*(1:(j+1)), ...)-dist(p*h*(0:j), ...))*V[i+2, rev(0:j)+1]) 
#'                                   + (1-dist((j+1)*p*h, ...)))
#'       }
#'       
#'     }
#'   }
#'   
#'   
#'   
#'   ret <- list()
#'   ret$V <- V
#'   ret$V_bar <- 1-V
#'   ret$t <- t_steps
#'   ret$x <- x_steps
#'   return(ret)
#'   
#'   
#' }
#' 
#' 
#' Finite_Ruin_Probability_2 <- function(h, t, Time, x, theta, dist, mu, r, demand, D,l, V_T = 0, func = FALSE, ... ){
#'   
#'   #k_1 = 2
#'   #beta_1 = 500
#'   
#'   
#'   # h = Time/K  
#'   t_steps <- seq(from = t, to = Time, by = h)
#' 
#'   lambda <- D*l*demand(theta)
#'   p <- (1+theta)*lambda*mu-r
#'   #print(p)
#'   if(p >0){
#'     x_steps <- seq(from = 0, to = x, by = p*h)
#'     
#'     K <- length(t_steps) - 1
#'     M <- length(x_steps) - 1
#'     
#'     V <- matrix(NA, nrow = K+1, ncol = M+1)
#'     if(func){
#'       V[K+1,] <- V_T(x_steps)
#'     }else{
#'       V[K+1,] <- V_T
#'     }
#'     
#'     
#'     # Note R starts indexing at 1
#'     for(i in seq(from = K-1, to = 0, by = -1)){
#'       # print(i)
#'       for(j in 0:M){
#'         if(j< M && i <= K-1){
#'           V[i+1,j+1] <- (1-pexp(h, rate = lambda))*V[i+1+1, j+1+1] +
#'             pexp(h, rate = lambda)*(sum((dist(p*h*(1:(j+1)), ...)-dist(p*h*(0:j), ...))*V[i+2, rev(0:j)+1])
#'                                     + (1-dist((j+1)*p*h, ...)))
#'         }else if(j == M && i <= K-1){
#'           # Assume V[i+1+1, j+1+1] = 0
#'           V[i+1,j+1] <- (1-pexp(h, rate = lambda))*0 +
#'             pexp(h, rate = lambda)*(sum((dist(p*h*(1:(j+1)), ...)-dist(p*h*(0:j), ...))*V[i+2, rev(0:j)+1])
#'                                     + (1-dist((j+1)*p*h, ...)))
#'         }
#'         
#'       }
#'     }
#'   }else if(p<0){
#'     #print(p)
#'     x_steps <- rev(seq(from = x, to = 0, by = p*h))
#' 
#'     K <- length(t_steps) - 1
#'     M <- length(x_steps) - 1
#' 
#'     V <- matrix(NA, nrow = K+1, ncol = M+1)
#'     if(func){
#'       V[K+1,] <- V_T(x_steps)
#'     }else{
#'       V[K+1,] <- V_T
#'     }
#' 
#' 
#' 
#'     # Note R starts indexing at 1
#'     for(i in seq(from = K-1, to = 0, by = -1)){
#'       # print(i)
#'       for(j in 0:M){
#'         if(j> 1 && i <= K-1){
#'           V[i+1,j+1] <- (1-pexp(h, rate = lambda))*V[i+1+1, j+1-1] +
#'             pexp(h, rate = lambda)*(sum((dist(-p*h*(1:(j-1)), ...)-dist(-p*h*(0:(j-2)), ...))*V[i+2, rev(0:(j-2))+1])
#'                                     + (1-dist((j-1)*(-p)*h, ...)))
#'         }else if(j == 0 && i <= K-1){
#'           # Assume V[i+1+1, j+1+1] = 0
#'           V[i+1,j+1] <- 1
#'         }else if(j == 1){
#'           V[i+1,j+1] <- 1
#'         }
#' 
#'       }
#'     }
#'   }else{
#'     # print(p)
#'     x_steps <- seq(from = 0, to = x, by = p_0)
#' 
#'     K <- length(t_steps) - 1
#'     M <- length(x_steps) - 1
#' 
#'     V <- matrix(NA, nrow = K+1, ncol = M+1)
#'     if(func){
#'       V[K+1,] <- V_T(x_steps)
#'     }else{
#'       V[K+1,] <- V_T
#'     }
#' 
#' 
#'     # Note R starts indexing at 1
#'     for(i in seq(from = K-1, to = 0, by = -1)){
#'       # print(i)
#'       for(j in 0:M){
#' 
#'         V[i+1,j+1] <- (1-pexp(h, rate = lambda))*V[i+1+1, j+1] +
#'             pexp(h, rate = lambda)*(sum((dist(p_0*(1:(j+1)), ...)-dist(p_0*(0:j), ...))*V[i+2, rev(0:j)+1])
#'                                     + (1-dist((j+1)*p_0, ...)))
#' 
#'       }
#'     }
#' 
#'   }
#' 
#'   
#'   
#'   
#'   ret <- list()
#'   ret$V <- V
#'   ret$V_bar <- 1-V
#'   ret$t <- t_steps
#'   ret$x <- x_steps
#'   return(ret)
#'   
#'   
#' }
#' 
#' 
#' 
#' Finite_Ruin_Probability_3 <- function(h, t, Time, x, theta, dist, mu, r, a, b, V_T = 0, func = FALSE, p_0 =0.005, ... ){
#'   
#'   #k_1 = 2
#'   #beta_1 = 500
#'   
#'   
#'   # h = Time/K  
#'   t_steps <- seq(from = t, to = Time, by = h)
#'   
#'   lambda <- max(0, a-b*theta)
#'   p <- (1+theta)*lambda*mu-r
#'   
#'   if(p>0){
#'     #print(p) 
#'     x_steps <- seq(from = 0, to = x, by = p*h)
#' 
#'     K <- length(t_steps) - 1
#'     M <- length(x_steps) - 1
#'     
#'     V <- matrix(NA, nrow = K+1, ncol = M+1)
#'     
#'     if(func){
#'       V[K+1,] <- V_T(x_steps)
#'     }else{
#'       V[K+1,] <- V_T
#'     }
#'     
#'   
#'     
#'     # Note R starts indexing at 1
#'     for(i in seq(from = K-1, to = 0, by = -1)){
#'       # print(i)
#'       for(j in 0:M){
#'         if(j< M && i <= K-1){
#'           V[i+1,j+1] <- (1-pexp(h, rate = lambda))*V[i+1+1, j+1+1] +
#'             pexp(h, rate = lambda)*(sum((dist(p*h*(1:(j+1)), ...)-dist(p*h*(0:j), ...))*V[i+2, rev(0:j)+1])
#'                                     + (1-dist((j+1)*p*h, ...)))
#'         }else if(j == M && i <= K-1){
#'           # Assume V[i+1+1, j+1+1] = 0
#'           V[i+1,j+1] <- (1-pexp(h, rate = lambda))*0 +
#'             pexp(h, rate = lambda)*(sum((dist(p*h*(1:(j+1)), ...)-dist(p*h*(0:j), ...))*V[i+2, rev(0:j)+1])
#'                                     + (1-dist((j+1)*p*h, ...)))
#'         }
#'         
#'       }
#'     }
#'     
#'   }
#'   
#'   # else if(p<0){
#'   #   #print(p) 
#'   #   x_steps <- rev(seq(from = x, to = 0, by = p*h))
#'   #   
#'   #   K <- length(t_steps) - 1
#'   #   M <- length(x_steps) - 1
#'   #   
#'   #   V <- matrix(NA, nrow = K+1, ncol = M+1)
#'   #   if(func){
#'   #     V[K+1,] <- V_T(x_steps)
#'   #   }else{
#'   #     V[K+1,] <- V_T
#'   #   }
#'   #   
#'   #   
#'   #   
#'   #   # Note R starts indexing at 1
#'   #   for(i in seq(from = K-1, to = 0, by = -1)){
#'   #     # print(i)
#'   #     for(j in 0:M){
#'   #       if(j> 0 && i <= K-1){
#'   #         V[i+1,j+1] <- (1-pexp(h, rate = lambda))*V[i+1+1, j+1-1] +
#'   #           pexp(h, rate = lambda)*(sum((dist(p*h*(1:(j)), ...)-dist(p*h*(0:(j-1)), ...))*V[i+2, rev(0:(j-1))+1])
#'   #                                   + (1-dist((j)*p*h, ...)))
#'   #       }else if(j == 0 && i <= K-1){
#'   #         # Assume V[i+1+1, j+1+1] = 0
#'   #         V[i+1,j+1] <- 1
#'   #       }
#'   #       
#'   #     }
#'   #   }
#'   # }else{
#'   #   # print(p) 
#'   #   x_steps <- seq(from = 0, to = x, by = p_0)
#'   #   
#'   #   K <- length(t_steps) - 1
#'   #   M <- length(x_steps) - 1
#'   #   
#'   #   V <- matrix(NA, nrow = K+1, ncol = M+1)
#'   #   if(func){
#'   #     V[K+1,] <- V_T(x_steps)
#'   #   }else{
#'   #     V[K+1,] <- V_T
#'   #   }
#'   #   
#'   #   
#'   #   # Note R starts indexing at 1
#'   #   for(i in seq(from = K-1, to = 0, by = -1)){
#'   #     # print(i)
#'   #     for(j in 0:M){
#'   #       
#'   #       V[i+1,j+1] <- (1-pexp(h, rate = lambda))*V[i+1+1, j+1] +
#'   #           pexp(h, rate = lambda)*(sum((dist(p_0*(1:(j+1)), ...)-dist(p_0*(0:j), ...))*V[i+2, rev(0:j)+1])
#'   #                                   + (1-dist((j+1)*p_0, ...)))
#'   #       
#'   #     }
#'   #   }
#'   #   
#'   # }
#'   
#'   
#'   
#'   ret <- list()
#'   ret$V <- V
#'   ret$V_bar <- 1-V
#'   ret$t <- t_steps
#'   ret$x <- x_steps
#'   return(ret)
#'   
#'   
#' }
#' 
#' 
#' 
#' ### Derivative probability method ----
#' Finite_Ruin_theta <- function(theta, h, t, Time, x, mu,a,b,r, rate = 1/mu){
#'   
#'   # theta = theta_min
#'   # h = 0.1
#'   # t = 0
#'   # Time = 2
#'   # x = 0.5
#'   # mu = 1
#'   # a = 1
#'   # b = 0.6
#'   # r = 0.25
#'   # rate = 1
#'   # x_calc = 0.5
#'   # t_calc = 0
#'   
#'   # h = T/K  
#'   t_steps <- seq(from = t, to = Time, by = h)
#'   
#'   lambda <- a-b*theta
#'   p <- (1+theta)*lambda*mu-r
#'   
#'   x_steps <- seq(from = 0, to = x, by = p*h)
#'   
#'   K <- length(t_steps) - 1
#'   M <- length(x_steps) - 1
#'   
#'   
#'   V_d <- matrix(NA, nrow = K+1, ncol = M+1)
#'   V_d[K+1,] <- 0
#'   
#'   V <- matrix(NA, nrow = K+1, ncol = M+1)
#'   V[K+1,] <- 0
#'   
#'   h_d <- 0.01
#'   
#'   # Note R starts indexing at 1
#'   for(i in seq(from = K-1, to = 0, by = -1)){
#'     #(M+i-K+1)
#'     for(j in 0:(M+i-K+1)){
#'       if(i == K-1){
#'         p <- (1+theta)*(mu)*(a-b*(theta))-r
#'         p_r <- (1+theta+h_d)*(mu)*(a-b*(theta+h_d))-r
#'         p_l <- (1+theta-h_d)*(mu)*(a-b*(theta-h_d))-r
#'         
#'         
#'         V[i+1,j+1] <- pexp(h, rate = lambda)*(1-pexp((j+1)*p*h, rate = 1/mu))
#'         V_d[i+1,j+1] <-   (pexp(h, rate = a-b*(theta+h_d))-pexp(h, rate = a-b*(theta-h_d)))*(1-pexp((j+1)*p*h, rate =1/mu))/(2*h_d) +
#'           pexp(h, rate = a-b*theta)*((1-pexp((j+1)*p_r*h, rate = 1/mu))-(1-pexp((j+1)*p_l*h, rate = 1/mu)))/(2*h_d)
#'         
#'         
#'       }else if(j< M && i < K-1){
#'         
#'         p <- (1+theta)*(mu)*(a-b*(theta))-r
#'         p_r <- (1+theta+h_d)*(mu)*(a-b*(theta+h_d))-r
#'         p_l <- (1+theta-h_d)*(mu)*(a-b*(theta-h_d))-r
#'         
#'         V[i+1,j+1] <- (1-pexp(h, rate = lambda))*V[i+1+1, j+1+1] + 
#'           pexp(h, rate = lambda)*(sum((pexp(p*h*(1:(j+1)), rate = 1/mu)-pexp(p*h*(0:j), rate= 1/mu))*V[i+2, rev(0:j)+1]) 
#'                                   + (1-pexp((j+1)*p*h, rate = 1/mu)))
#'         
#'         
#'         tmp_1 <- ((1-pexp(h, rate = a-b*(theta+h_d)))-(1-pexp(h, rate = a-b*(theta-h_d))))*V[i+1+1, j+1+1]/(2*h_d)   + 
#'           (1-pexp(h, rate = a-b*theta))*V_d[i+1+1,j+1+1] 
#'         
#'         
#'         tmp_2 <- (pexp(h, rate = a-b*(theta+h_d))-pexp(h, rate = a-b*(theta-h_d)))*
#'           (sum((pexp(p*h*(1:(j+1)), 1/mu)-pexp(p*h*(0:j), 1/mu))*V[i+2, rev(0:j)+1]) + (1-pexp((j+1)*p*h, 1/mu)))/(2*h_d) +
#'           pexp(h, rate = theta)*sum(((pexp(p_r*h*(1:(j+1)), rate = 1/mu)-pexp(p_r*h*(0:j), rate = 1/mu))-
#'                                        (pexp(p_l*h*(1:(j+1)), 1/mu)-pexp(p_l*h*(0:j), 1/mu)))*V[i+2, rev(0:j)+1]/(2*h_d) +
#'                                       ((1-pexp((j+1)*p_r*h, 1/mu))-(1-pexp((j+1)*p_l*h, 1/mu)))/(2*h_d) +
#'                                       sum((pexp(p*h*(1:(j+1)), 1/mu)-pexp(p*h*(0:j), 1/mu))*V_d[i+2, rev(0:j)+1]))
#'         
#'         V_d[i+1, j+1] <-   tmp_1 + tmp_2
#'       }
#'       # else if(j == M && i < K-1){
#'       # 
#'       # p <- (1+theta)*(mu)*(a-b*(theta))-r
#'       # p_r <- (1+theta+h_d)*(mu)*(a-b*(theta+h_d))-r
#'       # p_l <- (1+theta-h_d)*(mu)*(a-b*(theta-h_d))-r
#'       # 
#'       # 
#'       # 
#'       #   V[i+1,j+1] <- (1-pexp(h, rate = lambda))*0 +
#'       #     pexp(h, rate = lambda)*(sum((pexp(p*h*(1:(j+1)), 1/mu)-pexp(p*h*(0:j), 1/mu))*V[i+2, rev(0:j)+1])
#'       #                             + (1-pexp((j+1)*p*h, 1/mu)))
#'       # 
#'       # 
#'       #   V_d[i+1, j+1] <- #((1-pexp(h, rate = a-b*(theta+0.01)))-(1-pexp(h, rate = a-b*(theta-0.01))))*V[i+1, j+1]/0.01   +
#'       #     #(1-pexp(h, rate = a-b*theta))*V_d[i+1+1,j+1+1] +
#'       #     (pexp(h, rate = a-b*(theta+0.01))-pexp(h, rate = a-b*(theta-0.01)))*(sum((pexp(p*h*(1:(j+1)), rate = 1/mu)-
#'       #                                                                                 pexp(p*h*(0:j), rate = 1/mu))*V[i+2, rev(0:j)+1])
#'       #                                                                          + (1-pexp((j+1)*p*h, rate = 1/mu)))/0.01 +
#'       #     pexp(h, rate = theta)*sum(((pexp(p_r*h*(1:(j+1)), rate = 1/mu)-pexp(p_r*h*(0:j), rate = 1/mu))-
#'       #                                  (pexp(p_l*h*(1:(j+1)), 1/mu)-pexp(p_l*h*(0:j), 1/mu)))*V[i+2, rev(0:j)+1]/0.01 +
#'       #                                 ((1-pexp((j+1)*p_r*h, 1/mu))-(1-pexp((j+1)*p_l*h, 1/mu)))/0.01 +
#'       #                                 sum((pexp(p*h*(1:(j+1)), 1/mu)-pexp(p*h*(0:j), 1/mu))*V_d[i+2, rev(0:j)+1])
#'       #     )
#'       # 
#'       # }
#'       
#'       
#'     }
#'   }
#'   
#'   
#'   
#'   ret <- list()
#'   ret$V <- V
#'   ret$V_bar <- 1-V
#'   ret$V_d <- V_d
#'   ret$t <- t_steps
#'   ret$x <- x_steps
#'   return(ret)
#'   
#' }
#' 
#' 
#' Finite_Ruin_Probability_2_old <- function(h, t, Time, x, theta, dist, mu, r, demand, D,l, V_T = 0, func = FALSE, ... ){
#'   
#'   #k_1 = 2
#'   #beta_1 = 500
#'   
#'   
#'   # h = Time/K  
#'   t_steps <- seq(from = t, to = Time, by = h)
#'   
#'   lambda <- D*l*max(demand(theta)-0.05,0)
#'   p <- (1+theta)*lambda*mu-r
#'   #print(p)
#'   if(p <=0){
#'     ret <- list()
#'     ret$V <- matrix(1, nrow = 100, ncol = 100)
#'     ret$V_bar <- 1-V
#'     ret$t <- t_steps
#'     ret$x <- seq(from = 0, by = 1, length.out = 100)
#'     return(ret)
#'   }
#'   
#'   x_steps <- seq(from = 0, to = x, by = p*h)
#'   
#'   
#'   
#'   K <- length(t_steps) - 1
#'   M <- length(x_steps) - 1
#'   
#'   V <- matrix(NA, nrow = K+1, ncol = M+1)
#'   if(func){
#'     V[K+1,] <- V_T(x_steps)
#'   }else{
#'     V[K+1,] <- V_T
#'   }
#'   
#'   
#'   
#'   # Note R starts indexing at 1
#'   # for(i in seq(from = K-1, to = 0, by = -1)){
#'   #   print(paste0(i, " of ", length(seq(from = K-1, to = 0, by = -1))))
#'   #   for(j in 0:M){
#'   #     print(paste0(j, " of ", M))
#'   #     if(j< M && i <= K-1){
#'   #       V[i+1,j+1] <- (1-pexp(h, rate = lambda))*V[i+1+1, j+1+1] + 
#'   #         pexp(h, rate = lambda)*(sum((dist(p*h*(1:(j+1)), shape = k_1, scale = beta_1)-dist(p*h*(0:j), shape = k_1, scale = beta_1))*V[i+2, rev(0:j)+1]) 
#'   #                                 + (1-dist((j+1)*p*h, shape = k_1, scale = beta_1)))
#'   #     }else if(j == M && i <= K-1){
#'   #       # Assume V[i+1+1, j+1+1] = 0
#'   #       V[i+1,j+1] <- (1-pexp(h, rate = lambda))*0 + 
#'   #         pexp(h, rate = lambda)*(sum((dist(p*h*(1:(j+1)), shape = k_1, scale = beta_1)-dist(p*h*(0:j), shape = k_1, scale = beta_1))*V[i+2, rev(0:j)+1]) 
#'   #                                 + (1-dist((j+1)*p*h, shape = k_1, scale = beta_1)))
#'   #     }
#'   #     
#'   #   }
#'   # }
#'   
#'   # fig <- plot_ly(showscale = TRUE,
#'   #                type = "surface",
#'   #                x = ~ x_steps,
#'   #                y = ~ t_steps,
#'   #                z = ~ V, 
#'   #                colorscale = "Viridis", colorbar=list(title="") ) #list(c(0, 0.3, 1),c("blue","purple", "red"))
#'   # fig <- fig %>% layout(
#'   #   scene = list(
#'   #     xaxis = list(nticks = 10, title = "Surplus", range = c(0,x),tickfont = list(size = 15), titlefont = list(size = 20)),
#'   #     zaxis = list(tickvals = c(0.2, 0.4, 0.6, 0.8, 1), title= "Probability of Ruin", range = c(0,1),tickfont = list(size = 15), titlefont = list(size = 20)),
#'   #     yaxis = list(title = "Time",tickfont = list(size = 15), titlefont = list(size = 20)),
#'   #     camera = list(eye = list(x = 0, y = -1, z = 0.5)),
#'   #     aspectratio = list(x = 1, y = 1, z = 0.7)))
#'   # 
#'   # fig
#'   
#'   # Note R starts indexing at 1
#'   for(i in seq(from = K-1, to = 0, by = -1)){
#'     # print(i)
#'     for(j in 0:M){
#'       if(j< M && i <= K-1){
#'         V[i+1,j+1] <- (1-pexp(h, rate = lambda))*V[i+1+1, j+1+1] +
#'           pexp(h, rate = lambda)*(sum((dist(p*h*(1:(j+1)), ...)-dist(p*h*(0:j), ...))*V[i+2, rev(0:j)+1])
#'                                   + (1-dist((j+1)*p*h, ...)))
#'       }else if(j == M && i <= K-1){
#'         # Assume V[i+1+1, j+1+1] = 0
#'         V[i+1,j+1] <- (1-pexp(h, rate = lambda))*0 +
#'           pexp(h, rate = lambda)*(sum((dist(p*h*(1:(j+1)), ...)-dist(p*h*(0:j), ...))*V[i+2, rev(0:j)+1])
#'                                   + (1-dist((j+1)*p*h, ...)))
#'       }
#'       
#'     }
#'   }
#'   
#'   
#'   
#'   ret <- list()
#'   ret$V <- V
#'   ret$V_bar <- 1-V
#'   ret$t <- t_steps
#'   ret$x <- x_steps
#'   return(ret)
#'   
#'   
#' }
#' 
#' 
#' 
#' 
#' # Ordinay copula functions ----
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
