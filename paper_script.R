library(tidyverse)
library(ggnewscale) # add new color scale -> new_scale_color
# Set color brewer
library(RColorBrewer)
library(plotly)

source("FDM.R")
source("paper_function_per_levy_copula.R")


# Color codes ------
my_blue <- colorRampPalette(brewer.pal(n = 9, "Blues")[3:7])(8)[1:8]
my_red <- colorRampPalette(brewer.pal(n = 9, "Reds")[4:8])(8)[1:8]
my_grey <- brewer.pal(n = 9, "Greys")[9] 

plot_all <- function(ruin_lines, profit_lines, theta_lim = c(0,1), optimal_expected_profit, optimal_theta = NA, opt_ruin_theta = NA, text_offset = 0,
                     breaks = c(0, 10000, 20000, 30000, 40000), size = 2, y_lim = c(0,1.05) ){
  #my_blue = brewer.pal(n = 9, "Blues")[2:9] 
  my_blue <- colorRampPalette(brewer.pal(n = 9, "Blues")[3:7])(8)[1:8]
  my_grey <- brewer.pal(n = 9, "Greys")[9] 
  gg <- ggplot() +
    
    geom_line(aes(x = x, y = value, color = surplus), ruin_lines, size = 1) +
    scale_color_manual(values = my_blue) +
    
    # geom_point(mapping = aes(x = x, y = value, color = surplus), ruin_point, size = 2) +
    labs( color = c("Surplus")) +
    new_scale_color() + 
    
    geom_line(aes(x = x, y = y, color = id), profit_lines, size = 1, alpha = 0.6) +
    # geom_point(mapping = aes(x = x, y = y, color = id), data = profit_point, size = 2, shape = 15, alpha = 0.6) +
    #new_scale_color() +
    scale_color_manual(values = "#000000") +
    
    geom_hline(yintercept = 0, color = "black", size = 1) +
    geom_vline(xintercept = 0, color = "black", size = 1) +
    # geom_vline(xintercept = optimal_theta, color = "black", size = 0.5) +
    annotate("text", x = optimal_theta-text_offset , y = 1.03, label = bquote(theta[profit]^"*" == .(optimal_theta)), size = size ) +
    annotate("text", x = opt_ruin_theta + text_offset , y = 1.03, label = bquote(theta[ruin]^"*" == .(opt_ruin_theta)), size = size) +
    geom_segment(aes(x = optimal_theta, xend = optimal_theta, y = 0, yend = 1), linetype = 2) +
    geom_segment(aes(x = opt_ruin_theta, xend = opt_ruin_theta, y = 0, yend = 1), linetype = 2) +
    labs(x = expression(theta), y = "Probability of Ruin \n", color = "") +
    #ggtitle("") +
    theme_bw(base_size = 20) +
    scale_x_continuous(expand = c(0.0015, 0),limits = theta_lim) +
    scale_y_continuous(expand = c(0.002, 0), 
                       limits = y_lim,
                       sec.axis = sec_axis(~.*optimal_expected_profit, 
                                           name = "Expected Profit \n", breaks = breaks), label = scales::comma) +
    theme(legend.justification = "top")
  return(gg)
}



# Declare Constants ----
N <- c(1, 1)
k <- c(2, 2)
beta <- c(500, 500)
lambda <- c(800, 800)
fixed_cost <-  0.2*0.4*k*beta*lambda
r <- c(0,0) 

kendell <- c(0.05, 0.1, 0.25, 0.4, 0.65, 0.8)

# Todo, take in copula calculate the probabilities
demand <- function(theta){
  b1 <- c(-0.6, 4)
  b2 <- c(-0.6, 4.5)
  return(c(1/(1+exp(b1[1]+b1[2]*theta)), 1/(1+exp(b2[1]+b2[2]*theta))))
  #return(c(1,1))
}

#' @param b1 - logit parameter intercept
#' @param b2 - logit parameter loading slope
demand_1 <- function(theta){
  b1 <- -0.6
  b2 <- 4
  
  1/(1+exp(b1+b2*theta))
}




# Independent case ----

S_ <- function(x, a, b){
  
  return(c(1-pgamma(q = x , shape = a[1], scale = b[1]), 1-pgamma(q = x , shape = a[2], scale = b[2])))
}

# 

out_one_indp <- one_loading_inference_indp(N = N, 
                                           r = r, 
                                           fixed_cost = fixed_cost,
                                           lambda = lambda, 
                                           claim_mean = k*beta,
                                           demand = demand, 
                                           h_x = 100,
                                           x_surplus = c(100),# 5000, 15000, 20000, 40000), 
                                           theta_finess = 0.01,
                                           S_ = S_, 
                                           a = k, 
                                           b = beta)



# Dependent Case Clayton -----

theta_grid <- c(0.25, 0.3, 0.32, 0.34,0.36,0.38,0.39, 0.392, 0.394, 0.396, 0.398, 0.4, 0.402, 0.404, 0.406, 0.408, 0.41,0.42, 0.44, 0.46, 0.5)

x_surplus <- c(100)#, 2000, 5000, 10000)

nu <- 1

# independent copula ----

ind_copula <- function(a,b, cop_par) return(a*b)

clayton_indp <- list()
clayton_indp$indp <- list()
for(x in x_surplus){
  print(paste0("Surplus: ", x))
  clayton_indp$indp[[paste0(x)]] <- one_loading_inference_clayton(N = N, 
                                                                  r = r, 
                                                                  fixed_cost = fixed_cost, 
                                                                  lambda = lambda, 
                                                                  nu = nu, 
                                                                  k = k, 
                                                                  beta = beta,
                                                                  x_surplus = x,
                                                                  demand = demand, 
                                                                  h_x = 100, 
                                                                  f_z_max = 20000, 
                                                                  theta_grid = c(0.4),
                                                                  f_z_limit = 20000, 
                                                                  ord_copula = ind_copula,
                                                                  cop_par = 1) 
  
}


# Gumbel copula nu = 1 -----

gumbel <- function(a,b,cop_par){
  return(exp(-((-log(a))^cop_par + (-log(b))^cop_par )^(1/cop_par)))
}


cop_params <- 1/(1-kendell)
clayton_gumbel <- list()


for(i in cop_params){
  print(paste0("copula parameter: ", i))
  clayton_gumbel[[paste0(i)]] <- list()
  for(x in c(5000)){
    print(paste0("Surplus: ", x))
    clayton_gumbel[[paste0(i)]][[paste0(x)]] <- one_loading_inference_clayton(N = N, 
                                                                        r = r, 
                                                                        fixed_cost = fixed_cost,
                                                                        lambda = lambda, 
                                                                        nu = 1, 
                                                                        k = k, 
                                                                        beta = beta,
                                                                        x_surplus = x,
                                                                        demand = demand, 
                                                                        h_x = 100, 
                                                                        f_z_max = 20000, 
                                                                        theta_grid = theta_grid,
                                                                        f_z_limit = 20000,
                                                                        ord_copula = gumbel,
                                                                        cop_par = i)
    
    
  }
  
}


ruin_values <- list()


for( i in 1:length(cop_params)){

  ruin_values[[i]] <-clayton_gumbel[[i]]$`100`$df
  ruin_values[[i]]$kendell <- as.character(kendell[i])
  ruin_values[[i]]$cop_param <- cop_params[i]

}

ruin_values <- as.data.frame(do.call(rbind,ruin_values))

ggplot() + geom_line(aes(x = x, y = value, colour =  kendell), data = ruin_values, size = 1 )+theme_bw(24)+
  xlab("Loading")+
  xlab("Loading")+
  ylab("Probability of Ruin") +
  ggtitle("Gumbel Acquisition Structure")+
  labs(colour = "Kendell's tau")

#### Clayton nu = 1 ----


clayton <- function(a,b,cop_par){
  return((a^(-cop_par)+b^(-cop_par)-1)^(-1/cop_par))
}


cop_params <- 2*kendell/(1-kendell)
clayton_clayton <- list()


for(i in cop_params[2]){
  print(paste0("copula parameter: ", i))
  clayton_clayton[[paste0(i)]] <- list()
  for(x in x_surplus[3]){
    print(paste0("Surplus: ", x))
    clayton_clayton[[paste0(i)]][[paste0(x)]] <- one_loading_inference_clayton(N = N, 
                                                                              r = r, 
                                                                              fixed_cost = fixed_cost,
                                                                              lambda = lambda, 
                                                                              nu = 1, 
                                                                              k = k, 
                                                                              beta = beta,
                                                                              x_surplus = x,
                                                                              demand = demand, 
                                                                              h_x = 100, 
                                                                              f_z_max = 20000, 
                                                                              theta_grid = theta_grid,
                                                                              f_z_limit = 20000,
                                                                              ord_copula = clayton,
                                                                              cop_par = i)
    
    
  }
  
}

ruin_values <- list()


for( i in 1:length(cop_params)){

  
  ruin_values[[i]] <-clayton_clayton[[i]]$`5000`$df
  ruin_values[[i]]$kendell <- as.character(kendell[i])
  ruin_values[[i]]$cop_param <- cop_params[i]
  
}

ruin_values <- as.data.frame(do.call(rbind,ruin_values))

ggplot() + geom_line(aes(x = x, y = value, colour =  kendell), data = ruin_values, size = 1 )+theme_bw(24)+
  xlab("Loading")+
  xlab("Loading")+
  ylab("Probability of Ruin") +
  ggtitle("Clayton Acquisition Structure")+
  labs(colour = "Kendell's tau")



# Gumbel copula nu = 0.5  -----



cop_params <- 1/(1-kendell)
clayton_gumbel <- list()


for(i in cop_params){
  print(paste0("copula parameter: ", i))
  clayton_gumbel[[paste0(i)]] <- list()
  for(x in x_surplus){
    print(paste0("Surplus: ", x))
    clayton_gumbel[[paste0(i)]][[paste0(x)]] <- one_loading_inference_clayton(N = N, 
                                                                              r = r, 
                                                                              fixed_cost = fixed_cost,
                                                                              lambda = lambda, 
                                                                              nu = 0.5, 
                                                                              k = k, 
                                                                              beta = beta,
                                                                              x_surplus = x,
                                                                              demand = demand, 
                                                                              h_x = 100, 
                                                                              f_z_max = 20000, 
                                                                              theta_grid = theta_grid,
                                                                              f_z_limit = 20000,
                                                                              ord_copula = gumbel,
                                                                              cop_par = i)
    
    
  }
  
}




#### Clayton nu = 0.5 ----



cop_params <- 2*kendell/(1-kendell)
clayton_clayton <- list()


for(i in cop_params[1]){
  print(paste0("copula parameter: ", i))
  clayton_clayton[[paste0(i)]] <- list()
  for(x in x_surplus[3]){
    print(paste0("Surplus: ", x))
    clayton_clayton[[paste0(i)]][[paste0(x)]] <- one_loading_inference_clayton(N = N, 
                                                                               r = r, 
                                                                               fixed_cost = fixed_cost,
                                                                               lambda = lambda, 
                                                                               nu = 0.5, 
                                                                               k = k, 
                                                                               beta = beta,
                                                                               x_surplus = x,
                                                                               demand = demand, 
                                                                               h_x = 100, 
                                                                               f_z_max = 20000, 
                                                                               theta_grid = theta_grid,
                                                                               f_z_limit = 20000,
                                                                               ord_copula = clayton,
                                                                               cop_par = i)
    
    
  }
  
}


