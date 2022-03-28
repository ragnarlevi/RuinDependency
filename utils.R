#' @param N - number of customers
#' @param r - prop op cost
#' @param fixed_cost - fixed op cost
#' @param lambda - inesity ber customer
#' @param claim_mean - Severity mean
#' @param x_surplus - surplusus to calculate ruin probabiity as function of theta
#' @param demand function
#' @param S_ severity survival
#' @param h_x pde step size 
#' 
#' 

library(progress)
library(tidyverse)

plot_ruin <- function(V, t1, t2, type = "surface", t1_opt = NA, t2_opt = NA, v_min = NA ){
  fig_ruin <- plot_ly(showscale = TRUE,
                      type = type,
                      y = ~ t1,
                      x = ~ t2,
                      z = ~ V,
                      colorbar=list(title="Ruin Prob.", 
                                    tickfont = list(size = 20),
                                    titlefont = list(size = 24)))
  
  fig_ruin <- fig_ruin %>% layout(yaxis = list(title = "&#952;<sub>1</sub>", 
                                               tickfont = list(size = 22), 
                                               titlefont = list(size = 24), 
                                               range = c(0, 0.7)),
                                  xaxis = list(title = "&#952;<sub>2</sub>", 
                                               tickfont = list(size = 22), 
                                               titlefont = list(size = 24),
                                               range = c(0, 0.7)))
  
  
  text <- paste0("(",t2_opt, ",", t1_opt, ")-> ",v_min, sep ="" )
  
  fig_ruin <- fig_ruin %>% colorbar(limits = c(0,1) ) %>%
    add_annotations(x = t2_opt, y = t1_opt, text = text, font = list(color = "black", 
                                                                     size = 20), arrowcolor="black")
  
  
  
  fig_ruin
}

plot_profit<- function(V, t1, t2, type = "surface", t1_opt = NA, t2_opt = NA, prof_max = NA ){
  fig_ruin <- plot_ly(showscale = TRUE,
                      type = type,
                      y = ~ t1,
                      x = ~ t2,
                      z = ~ V,
                      colorbar=list(title="Expected Profit", 
                                    tickfont = list(size = 20), 
                                    titlefont = list(size = 24)  ) )
  fig_ruin <- fig_ruin %>% layout(yaxis = list(title = "&#952;<sub>1</sub>", 
                                               tickfont = list(size = 22),
                                               titlefont = list(size = 24)),
                                  xaxis = list(title = "&#952;<sub>2</sub>", 
                                               tickfont = list(size = 22), 
                                               titlefont = list(size = 24)))
  
  text <- paste0("(",t2_opt, ",", t1_opt, ")-> ", prof_max, sep ="" )
  
  fig_ruin <- fig_ruin  %>%
    add_annotations(x = t2_opt, y = t1_opt, text = text, font = list(color = "black", size = 20),
                    arrowcolor="black")
  
  
  
  fig_ruin
}



infinite_survival3 <- function(to, h_x, p, lambda_true, mean_true, f_ = "", S_ = "", print_time = FALSE, lower = 0, ...){
  
  
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


plot_all <- function(ruin_lines, profit_lines, theta_lim = c(0,1), optimal_expected_profit, optimal_theta = NA, opt_ruin_theta = NA, text_offset = 0,
                     breaks = c(0, 10000, 20000, 30000, 40000), size = 2, y_lim = c(0,1.05) ){
  #my_blue = brewer.pal(n = 9, "Blues")[2:9] 
  my_blue <- colorRampPalette(brewer.pal(n = 9, "Blues")[3:7])(8)[1:8]
  my_grey <- brewer.pal(n = 9, "Greys")[9] 
  gg <- ggplot() +
    
    geom_line(aes(x = x, y = value, color = factor(surplus, levels = x_surplus, ordered = TRUE)), ruin_lines, size = 1) +
    scale_color_manual(values = my_blue[c(1, 3,6,8)]) +
    
    # geom_point(mapping = aes(x = x, y = value, color = surplus), ruin_point, size = 2) +
    labs( color = c("Initial Surplus")) +
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





one_loading_inference_indp <- function(N, r, fixed_cost, lambda = rep(1, length(N)), claim_mean, 
                                       x_surplus, demand, h_x = 10, theta_finess = 0.05, S_, ...){
  
  
  
  # Find minimum theta
  thetas <- seq(from = -1, to = 5, by = 0.001)
  # if demand returns a vector then doing vector thingy can lear to some errors
  expected_income <- list()
  for( i in 1:length(thetas)){
    
    l <- N*demand(thetas[i])*lambda/sum(N*demand(thetas[i])*lambda)
    claim_mean_tmp <- sum(claim_mean*l)
    expected_income[[i]] <- (1+thetas[i])*sum(N*demand(thetas[i])*(lambda*claim_mean-r)) - sum(fixed_cost) - claim_mean_tmp*sum(lambda*N*demand(thetas[i]))
    }
  expected_income <- unlist(expected_income)
  
  
  if(all(expected_income<=0) ){
    stop("Always ruin")
  }
  
  
  theta_optimal <- thetas[expected_income == max(expected_income)]
  theta_min <- thetas[min(expected_income[expected_income>0 & thetas < theta_optimal]) == expected_income]
  # print(theta_min)
  theta_max <- min(thetas[min(expected_income[expected_income>0 & thetas > theta_optimal]) == expected_income],5)
  
  
  if(theta_max == thetas[length(thetas)]){
    warning("No theta above optimal. theta that gives negative expected value")
  }
  
  
  
  # Ruin probability as function of theta for given surplus
  print(theta_min)
  df_dep_one <- list()
  thetas_ok <- seq(from = theta_min, to = theta_max, by = theta_finess)
  for(i in 1:length(x_surplus)){
    print(i)
    thetas_to_df <- list()
    value_to_df <- list()
    for(j in 1:length(thetas_ok)){

      
      #vector lambda 
      l <- N*demand(thetas_ok[j])*lambda/sum(N*demand(thetas_ok[j])*lambda)

      
      
      # The survival distribution is a sum
      S_tmp <- function(x, ...){
        
        return(sum(l*S_(x, ...)))
        
      }
      
      S_tmp <- Vectorize(S_tmp, vectorize.args = "x")
      
      
      
      mean_true <-  sum(claim_mean*l)
      lambda_true <-  sum(lambda*N*demand(thetas_ok[j])) 
      p <-  (1+thetas_ok[j])*sum(N*demand(thetas_ok[j])*(lambda*claim_mean)) - sum(fixed_cost)
      #    print(paste("p ", p))
      #     print(paste("mean true ", mean_true))
      #     print(paste("lambda true ", lambda_true))
      if(p<=mean_true*lambda_true){
        warning(paste0("Ruin will happen ", thetas_ok[j]))
        
        thetas_to_df[[j]] <- NA
        value_to_df[[j]] <- NA
      }else {
        out <- infinite_survival3(to = x_surplus[i], 
                                h_x = h_x, 
                                p = p,
                                lambda_true = lambda_true,
                                mean_true = mean_true,
                                S_ = S_tmp,
                                ...)
        
        length_out <- length(out$V)
        
        thetas_to_df[[j]] <- thetas_ok[j]
        value_to_df[[j]] <- out$V[length_out]
      }
      
      
      
    }
    
    df_dep_one[[i]] <- data.frame(x = unlist(thetas_to_df), value = unlist(value_to_df), surplus = paste(x_surplus[i]))
  }
  
  df_dep_one <- do.call(what = rbind, df_dep_one)
  
  
  
  # plot
  # gg <- plot_all(ruin_lines = df_dep_one, 
  #                profit_lines = data.frame(x = thetas, y = expected_income/max(expected_income), id = "Expected Profit"), 
  #                theta_lim = c(0,theta_max), 
  #                optimal_expected_profit = max(expected_income), optimal_theta = theta_optimal)
  # 
  # 
  surpluses <- unique(df_dep_one$surplus)
  optimal <- list()
  for( i in 1:length(surpluses)){
    
    optimal[[i]] <- thetas_ok[df_dep_one$value[df_dep_one$surplus == surpluses[i]] 
                              == min(df_dep_one$value[df_dep_one$surplus == surpluses[i]])]
    
  }
  optimal <- unlist(optimal)
  if(length(unique(optimal)) != 1) warning("more than one optimal value for ruin")
  
  ret <- list()
  #ret$plot <- gg
  ret$df <- df_dep_one
  ret$expected_proift <- expected_income
  ret$thetas <- thetas
  ret$thetas_ok <- thetas_ok
  ret$optimal_expected_profit <- max(expected_income)
  ret$optimal_theta <- theta_optimal
  ret$opt_ruin_theta <- optimal[1]
  return(ret)
  
  
}



# Assuming only two
one_loading_inference_dep <- function(N, r, fixed_cost, lambda = rep(1, length(N)), nu, k, beta,  x_surplus, demand, 
                                      h_x = 10, f_z_max = 10000, theta_finess = 0.1, f_z_limit = f_z_max, verbose = TRUE){
  
  
  
  
  
  
  S_common <- function(x1,x2,a,b, l1, l2, l_common, nu){
    
    tmp <- (((l1*(1-pgamma(q = x1, shape = a[1], scale = b[1])))^(-nu) + 
               (l2*(1-pgamma(q = x2, shape = a[2], scale = b[2])))^(-nu))^(-1/nu))/l_common
    return(tmp)
  }
  
  
  # We use the clayton copula
  clayton_copula <- function(x1,x2, nu){
    tmp <- ((x1)^(-nu)+(x2)^(-nu))^(-1/nu)
    return(tmp)
  }
  
  
  # assume exponential margin, 
  f_bivariate <- function(x1, x2, nu, a, b, l1, l2, l_common){
    
    copula_density <- function(u1, u2, nu){
      
      return((1+nu)*(u1^(-nu)+u2^(-nu))^(-1/nu - 2)*(u1^(-nu-1)*u2^(-nu-1)))
    }
    
    return(l1*l2*dgamma(x = x1, shape = a[1], scale = b[1])*dgamma(x = x2, shape = a[2], scale = b[2])*
             copula_density(u1 = l1*(1-pgamma(q = x1, shape = a[1], scale = b[1])), 
                            u2 = l2*(1-pgamma(q = x2, shape = a[2], scale = b[2])), nu)/l_common)
    
  }
  
  
  # Density of independent jumps
  f_indp_perp <- function(x, a,b, nu, l1, l2, l_common){
    l1_perp <- l1-l_common
    
    tmp1 <- l1*dgamma(x = x, shape = a, scale = b)/l1_perp
    tmp2 <- ((l1*(1-pgamma(q = x, shape = a, scale = b)))^(-nu) + l2^(-nu))^(-1/nu-1)*
      (l1*(1-pgamma(q = x, shape = a, scale = b)))^(-nu-1)*l1*dgamma(x = x, shape = a, scale = b)/l1_perp
    
    if(x< f_z_limit){
      return(tmp1 - tmp2)
    }else{
      return(0)
    }
    
  }
  
  f_indp_perp <- Vectorize(f_indp_perp, vectorize.args = "x")
  # density of the sum of common jumps is
  f_z <- function(z, a,b, nu, l1, l2, l_common){
    # print(z)  
    int_function <- function(x) f_bivariate(x, z-x, nu = nu, a = a, b = b, l1, l2, l_common)
    
    if(f_z_limit<z){
      return(0)
    }else{
      return(integrate(int_function, lower = 0, upper = z)$value) 
    }
  }
  f_z <- Vectorize(f_z, vectorize.args = "z")
  
  
  # Density of the marginal common jump
  f_1perp <- function(x, l1, l2, a,b, nu){
    
    l_common <- clayton_copula(l1,l2,nu)
    
    ((l1*(1-pgamma(q = x, shape = a, scale = b)))^(-nu) + l2^(-nu))^(-1/nu-1)*
      (l1*(1-pgamma(q = x, shape = a, scale = b)))^(-nu-1)*l1*dgamma(x = x, shape = a, scale = b)/l_common
    
  }
  f_1perp <- Vectorize(f_1perp, vectorize.args = "x")
  
  
  
  
  f_combined <- function(x, l1, l2, nu,a, b, p_1, p_1_0, p_2, p_0_2){
    
    l_common <- clayton_copula(l1,l2,nu)
    l1_perp <- l1 - l_common
    l2_perp <- l2 - l_common
    
    l_tilde <- p_1*l1_perp + p_2*l2_perp + (p_1_0 + p_0_2 + p_1*p_2)*l_common
    
    
    p_1*l1_perp*f_indp_perp(x = x, a = a[1], b = b[1], nu = nu, l1 = l1, l2 = l2, l_common = l_common)/l_tilde + 
      p_2*l2_perp*f_indp_perp(x = x, a = a[2], b = b[2], nu = nu, l1 = l2, l2 = l1, l_common = l_common)/l_tilde +
      p_1_0*l_common*f_1perp(x = x, l1 = l1, l2 = l2, a = a[1], b = b[1], nu = nu)/l_tilde +
      p_0_2*l_common*f_1perp(x = x, l1 = l2, l2 = l1, a = a[2], b = b[2], nu = nu)/l_tilde +
      p_1*p_2*l_common*f_z(z = x, a = a, b = b, nu = nu, l1 = l1, l2 = l2, l_common = l_common )/l_tilde
    
    
  }
  f_combined <- Vectorize(f_combined, vectorize.args = "x")
  
  
  
  mean_func <- function(x ,l1, l2, nu,a, b, p_1, p_1_0, p_2, p_0_2 ){
    x*f_combined(x, l1, l2, nu,a, b, p_1, p_1_0, p_2, p_0_2 )
  }
  
  
  
  # Find minimum theta
  thetas <- seq(from = 0, to = 2, by = 0.01)
  # if demand returns a vector then doing vector thingy can lear to some errors
  expected_income <- list()
  for( i in 1:length(thetas)){
    demands <- demand(thetas[i])
    #l1 <- N[1]*demands[1]*lambda[1]
    #l2 <- N[2]*demands[2]*lambda[2]
    l1 <- N[1]*lambda[1]
    l2 <- N[2]*lambda[2]
    
    l_common <- clayton_copula(l1,l2,nu)
    l1_perp <- l1 - l_common
    l2_perp <- l2 - l_common
    
    # Define probabilities
    p_1 <- demands[1]
    p_2 <- demands[2]
    p_1_0 <- demands[1]*(1-demands[2])
    p_0_2 <- demands[2]*(1-demands[1])
    p_1_2 <- demands[1]*demands[2]
    
    
    sum_to_one <- integrate(f_combined,lower = 0, upper = f_z_max, 
                            l1 = l1, l2 = l2, nu = nu,a = k, b = beta, 
                            p_1 = p_1, p_1_0 = p_1_0, p_2 = p_2, p_0_2 = p_0_2)$value
    
    
    
    
    if(sum_to_one <0.9999){
      stop(paste0("f_z sum is ", sum_to_one))
    }
    
    l_all <- p_1*l1_perp + p_2*l2_perp + (p_1_0 + p_0_2 + p_1*p_2)*l_common
    
    claim_mean_all <- integrate(mean_func, lower = 0, upper = f_z_max, 
                                l1 = l1, l2 = l2, nu = nu,a = k, b = beta, 
                                p_1 = p_1, p_1_0 = p_1_0, p_2 = p_2, p_0_2 = p_0_2)$value
    
    expected_income[[i]] <- (1+thetas[i])*sum(c(l1*demands[1], l2*demands[2])*k*beta ) - sum(fixed_cost) - claim_mean_all*l_all
    
  }
  
  print("Expected profit done.")
  expected_income <- unlist(expected_income)
  plot(thetas, expected_income)
  
  
  if(all(expected_income<=0) ){
    stop("Always ruin")
  }
  
  
  theta_optimal <- thetas[expected_income == max(expected_income)]
  theta_min <- thetas[min(expected_income[expected_income>0 & thetas < theta_optimal]) == expected_income]
  #  print(theta_min)
  theta_max <- thetas[min(expected_income[expected_income>0 & thetas > theta_optimal]) == expected_income]
  
  if(theta_max == thetas[length(thetas)]){
    warning("No theta above optimal theta that gives negative expected value")
  }
  
  
  # Ruin probability as function of theta for given surplus
  
  df_dep_one <- list()
  thetas_ok <- seq(from = theta_min, to = theta_max, by = theta_finess)
  if(verbose){
    pb <- progress_bar$new(format = "  downloading [:bar] :percent eta: :eta",
                           total = length(x_surplus)*length(thetas_ok))
  }
  for(i in 1:length(x_surplus)){
    print(i)
    thetas_to_df <- list()
    value_to_df <- list()
    for(j in 1:length(thetas_ok)){
      
      
      demands <- demand(thetas_ok[j])
      l1 <- N[1]*lambda[1]
      l2 <- N[2]*lambda[2]
      
      l_common <- clayton_copula(l1,l2,nu)
      l1_perp <- l1 - l_common
      l2_perp <- l2 - l_common
      
      
      # Define probabilities
      p_1 <- demands[1]
      p_2 <- demands[2]
      p_1_0 <- demands[1]*(1-demands[2])
      p_0_2 <- demands[2]*(1-demands[1])
      p_1_2 <- demands[1]*demands[2]
      
      
      
      
      l_all <- p_1*l1_perp + p_2*l2_perp + (p_1_0 + p_0_2 + p_1*p_2)*l_common
      
      
      mean_true <-  integrate(mean_func, lower = 0, upper = f_z_max, 
                              l1 = l1, l2 = l2, nu = nu,a = k, b = beta, 
                              p_1 = p_1, p_1_0 = p_1_0, p_2 = p_2, p_0_2 = p_0_2)$value
      # print("mean_true calculated")
      lambda_true <-  l_all
      p <-  (1+thetas_ok[j])*sum(c(l1*demands[1], l2*demands[2])*k*beta ) - sum(fixed_cost)
      
      if(p<=mean_true*lambda_true){
        warning(paste0("Ruin will happen ", thetas_ok[j]))
        thetas_to_df[[j]] <- NA
        value_to_df[[j]] <- NA
      }else{
        out <- infinite_survival3(to = x_surplus[i], 
                                 h_x = h_x, 
                                 p = p,
                                 lambda_true = lambda_true,
                                 mean_true = mean_true,
                                 f_ = f_combined,
                                 l1 = l1, l2 = l2, nu = nu,a = k, b = beta, p_1 = p_1, p_1_0 = p_1_0, p_2 = p_2, p_0_2 = p_0_2)
        
        length_out <- length(out$V)
        #print("out$V calculated")
        
        thetas_to_df[[j]] <- thetas_ok[j]
        value_to_df[[j]] <- out$V[length_out] 
      }
      
      if(verbose){
        pb$tick()
      }
      
    }
    
    df_dep_one[[i]] <- data.frame(x = unlist(thetas_to_df), value = unlist(value_to_df), surplus = paste(x_surplus[i]))
    
  }
  
  df_dep_one <- do.call(what = rbind, df_dep_one)
  
  
  
  # plot
  # gg <- plot_all(ruin_lines = df_dep_one, 
  #                profit_lines = data.frame(x = thetas, y = expected_income/max(expected_income), id = "Expected Profit"), 
  #                theta_lim = c(0,theta_max), 
  #                optimal_expected_profit = max(expected_income))
  
  surpluses <- unique(df_dep_one$surplus)
  optimal <- list()
  for( i in 1:length(surpluses)){
    
    optimal[[i]] <- thetas_ok[df_dep_one$value[df_dep_one$surplus == surpluses[i]] 
                              == min(df_dep_one$value[df_dep_one$surplus == surpluses[i]])]
    
  }
  optimal <- unlist(optimal)
  if(length(unique(optimal)) != 1) warning("more than one optimal value for ruin")
  
  
  ret <- list()
  # ret$plot <- gg
  ret$df <- df_dep_one
  ret$expected_proift <- expected_income
  ret$thetas <- thetas
  ret$thetas_ok <- thetas_ok
  ret$optimal_expected_profit <- max(expected_income)
  ret$optimal_theta <- theta_optimal
  ret$opt_ruin_theta <- optimal[1]
  return(ret)
  
  
}



# Assuming only two
one_loading_inference_clayton <- function(N, r, fixed_cost, lambda = rep(1, length(N)), nu, k, beta,  x_surplus, demand, 
                                      h_x = 10, f_z_max = 10000, f_z_limit = f_z_max, ord_copula, theta_grid, cop_par, verbose = TRUE){
  
  "Similar to one_loading_inference_dep but can take copula for acquisition"
  
  
  
  
  S_common <- function(x1,x2,a,b, l1, l2, l_common, nu){
    
    tmp <- (((l1*(1-pgamma(q = x1, shape = a[1], scale = b[1])))^(-nu) + 
               (l2*(1-pgamma(q = x2, shape = a[2], scale = b[2])))^(-nu))^(-1/nu))/l_common
    return(tmp)
  }
  

  
  # We use the clayton Levy copula
  clayton_copula <- function(x1,x2, nu){
    tmp <- ((x1)^(-nu)+(x2)^(-nu))^(-1/nu)
    return(tmp)
  }
  
  
  # assume gamma margin, straight forward to adjust for more general case
  f_bivariate <- function(x1, x2, nu, a, b, l1, l2, l_common){
    
    copula_density <- function(u1, u2, nu){
      
      return((1+nu)*(u1^(-nu)+u2^(-nu))^(-1/nu - 2)*(u1^(-nu-1)*u2^(-nu-1)))
    }
    
    val <- l1*l2*dgamma(x = x1, shape = a[1], scale = b[1])*dgamma(x = x2, shape = a[2], scale = b[2])*
      copula_density(u1 = l1*(1-pgamma(q = x1, shape = a[1], scale = b[1])), 
                     u2 = l2*(1-pgamma(q = x2, shape = a[2], scale = b[2])), nu)/l_common
    
    val[is.nan(val)] <- 0

    
    return(val)
    
  }
  
  
  # Density of independent jumps
  f_indp_perp <- function(x, a,b, nu, l1, l2, l_common){
    l1_perp <- l1-l_common
    
    tmp1 <- l1*dgamma(x = x, shape = a, scale = b)/l1_perp
    tmp2 <- ((l1*(1-pgamma(q = x, shape = a, scale = b)))^(-nu) + l2^(-nu))^(-1/nu-1)*
      (l1*(1-pgamma(q = x, shape = a, scale = b)))^(-nu-1)*l1*dgamma(x = x, shape = a, scale = b)/l1_perp
    
    if(x < f_z_limit){
      return(tmp1 - tmp2)
    }else{
      return(0)
    }
    
  }
  
  f_indp_perp <- Vectorize(f_indp_perp, vectorize.args = "x")
  # density of the sum of common jumps is
  f_z <- function(z, a,b, nu, l1, l2, l_common){
    int_function <- function(x){
      val <- f_bivariate(x, z-x, nu = nu, a = a, b = b, l1, l2, l_common)
      
      }
    
    if(f_z_limit<z){
      return(0)
    }else{
      return(integrate(int_function, lower = 0, upper = z)$value) 
    }
  }
  
  f_z <- Vectorize(f_z, vectorize.args = "z")
  
  
  # Density of the marginal common jump
  
  f_1perp <- function(x, l1, l2, a,b, nu){
    
    l_common <- clayton_copula(l1,l2,nu)
    
    ((l1*(1-pgamma(q = x, shape = a, scale = b)))^(-nu) + l2^(-nu))^(-1/nu-1)*
      (l1*(1-pgamma(q = x, shape = a, scale = b)))^(-nu-1)*l1*dgamma(x = x, shape = a, scale = b)/l_common
    
  }
  
  f_1perp <- Vectorize(f_1perp, vectorize.args = "x")
  
  #integrate(f_1perp, lower = 0, upper = 100, l1 = 3, l2 = 4, a = 2, b = 5, nu = 0.2)
  
  
  f_combined <- function(x, l1, l2, nu, a, b, p_1, p_1_0, p_2, p_0_2, p_1_2 = p_1*p_2){
    
    l_common <- clayton_copula(l1,l2,nu)
    l1_perp <- l1 - l_common
    l2_perp <- l2 - l_common
    
    l_tilde <- p_1*l1_perp + p_2*l2_perp + (p_1_0 + p_0_2 + p_1_2)*l_common
    
    
    val <- p_1*l1_perp*f_indp_perp(x = x, a = a[1], b = b[1], nu = nu, l1 = l1, l2 = l2, l_common = l_common)/l_tilde + 
      p_2*l2_perp*f_indp_perp(x = x, a = a[2], b = b[2], nu = nu, l1 = l2, l2 = l1, l_common = l_common)/l_tilde +
      p_1_0*l_common*f_1perp(x = x, l1 = l1, l2 = l2, a = a[1], b = b[1], nu = nu)/l_tilde +
      p_0_2*l_common*f_1perp(x = x, l1 = l2, l2 = l1, a = a[2], b = b[2], nu = nu)/l_tilde +
      p_1_2*l_common*f_z(z = x, a = a, b = b, nu = nu, l1 = l1, l2 = l2, l_common = l_common )/l_tilde
    # print(val)
    
    val[is.nan(val)] <- 0
    val
    
    
    
    
    
  }
  f_combined <- Vectorize(f_combined, vectorize.args = "x")
  
  
  
  mean_func <- function(x ,l1, l2, nu,a, b, p_1, p_1_0, p_2, p_0_2, p_1_2 = p_1*p_2 ){
    x*f_combined(x, l1, l2, nu,a, b, p_1, p_1_0, p_2, p_0_2,p_1_2 )
  }
  
  
  
  
  
  
  # Find minimum theta
  thetas <- seq(from = 0, to = 2, by = 0.01)
  # if demand returns a vector then doing vector thingy can lear to some errors
  expected_income <- list()
  for( i in 1:length(thetas)){
  #  print(i)
    demands <- demand(thetas[i])

    
    l1 <- N[1]*lambda[1]
    l2 <- N[2]*lambda[2]
    
    l_common <- clayton_copula(l1,l2,nu)
    l1_perp <- l1 - l_common
    l2_perp <- l2 - l_common

    # Define probabilities
    p_1 <- demands[1]
    p_2 <- demands[2]
    
    
    p_1_0 <- ord_copula(1,1-p_2, cop_par)  - ord_copula(1-p_1, 1-p_2, cop_par)
    p_0_2 <- ord_copula(1-p_1,1, cop_par) - ord_copula(1-p_1, 1-p_2, cop_par)
    p_1_2 <- ord_copula(1,1, cop_par) - ord_copula(1-p_1,1, cop_par) - ord_copula(1,1-p_2, cop_par) + ord_copula(1-p_1,1-p_2, cop_par)
    
    p_1_0
    p_0_2
    p_1_2
    
    # p_1_0 <- demands[1]*(1-demands[2])
    # p_0_2 <- demands[2]*(1-demands[1])
    # p_1_2 <- demands[1]*demands[2]
    # 
    # p_1_0
    # p_0_2
    # p_1_2
    
    
    if(any(c(p_1_0>1,p_0_2>1, p_1_2>1))) stop("ordinary copula probability bigger than 1")
    if(any(c(p_1_0<0,p_0_2<0, p_1_2<0))) stop("ordinary copula probability less than 0")
    

    # if(p_1_0 != demands[1]*(1-demands[2])) print(thetas[i])
    # if(p_0_2 != demands[2]*(1-demands[1])) warning("p_0_2 not same as indp")
    # if(p_1_2 != demands[1]*demands[2]) warning("p_1_2 not same as indp")
    sum_to_one <- integrate(f_combined,lower = 0, upper = f_z_max, 
                            l1 = l1, l2 = l2, nu = nu,a = k, b = beta, 
                            p_1 = p_1, p_1_0 = p_1_0, p_2 = p_2, p_0_2 = p_0_2, p_1_2 = p_1_2,subdivisions = 2000)$value
    
    
    
    
    if(sum_to_one <0.9999){
      stop(paste0("f_z sum is ", sum_to_one))
    }
    
    l_all <- p_1*l1_perp + p_2*l2_perp + (p_1_0 + p_0_2 + p_1_2)*l_common
    
    claim_mean_all <- integrate(mean_func, lower = 0, upper = f_z_max, 
                                l1 = l1, l2 = l2, nu = nu,a = k, b = beta, 
                                p_1 = p_1, p_1_0 = p_1_0, p_2 = p_2, p_0_2 = p_0_2, p_1_2 = p_1_2)$value
    
    expected_income[[i]] <- (1+thetas[i])*sum(c(l1*demands[1], l2*demands[2])*k*beta ) - sum(fixed_cost) - claim_mean_all*l_all
   
  }
  
  print("Expected profit done.")
  expected_income <- unlist(expected_income)
  plot(thetas, expected_income)
  

  
  
  if(all(expected_income<=0) ){
    stop("Always ruin")
  }
  
  
  theta_optimal <- thetas[expected_income == max(expected_income)]
  theta_min <- thetas[min(expected_income[expected_income>0 & thetas < theta_optimal]) == expected_income]
  #  print(theta_min)
  theta_max <- thetas[min(expected_income[expected_income>0 & thetas > theta_optimal]) == expected_income]
  
  if(theta_max == thetas[length(thetas)]){
    warning("No theta above optimal theta that gives negative expected value")
  }
  
  
  
  thetas_ok <- theta_grid
  # Ruin probability as function of theta for given surplus
  if(verbose){
    pb <- progress_bar$new(format = "  downloading [:bar] :percent eta: :eta",
                           total = length(x_surplus)*length(thetas_ok))
  }
  
  df_dep_one <- list()
  for(i in 1:length(x_surplus)){
    print(i)
    thetas_to_df <- list()
    value_to_df <- list()
    for(j in 1:length(thetas_ok)){
      
      
      demands <- demand(thetas_ok[j])
      l1 <- N[1]*lambda[1]
      l2 <- N[2]*lambda[2]
      
      l_common <- clayton_copula(l1,l2,nu)
      l1_perp <- l1 - l_common
      l2_perp <- l2 - l_common
      
      
      # Define probabilities
      p_1 <- demands[1]
      p_2 <- demands[2]
      zero_bid <- 1-demand(0)
      
      
      p_1_0 <- ord_copula(1,1-p_2, cop_par) - ord_copula(1-p_1, 1-p_2, cop_par)
      p_0_2 <- ord_copula(1-p_1,1, cop_par) - ord_copula(1-p_1, 1-p_2, cop_par)
      p_1_2 <- ord_copula(1,1, cop_par) - ord_copula(1-p_1,1, cop_par) - ord_copula(1,1-p_2, cop_par) + ord_copula(1-p_1,1-p_2, cop_par)
      
      
      if(any(c(p_1_0>1,p_0_2>1, p_1_2>1))){
        warning("ordinary copula probability bigger than 1") 
        next
      }
      if(any(c(p_1_0<0,p_0_2<0, p_1_2<0))){
        warning("ordinary copula probability less than 0")
        next
        
      } 
  
      
      l_all <- p_1*l1_perp + p_2*l2_perp + (p_1_0 + p_0_2 + p_1_2)*l_common
      
      
      
      mean_true <-  integrate(mean_func, lower = 0, upper = f_z_max, 
                              l1 = l1, l2 = l2, nu = nu,a = k, b = beta, 
                              p_1 = p_1, p_1_0 = p_1_0, p_2 = p_2, p_0_2 = p_0_2, p_1_2 = p_1_2)$value
     # print("mean_true calculated")
      lambda_true <-  l_all
      p <-  (1+thetas_ok[j])*sum(c(l1*demands[1], l2*demands[2])*k*beta ) - sum(fixed_cost)
      
      if(p<=mean_true*lambda_true){
        warning(paste0("Ruin will happen ", thetas_ok[j]))
        thetas_to_df[[j]] <- NA
        value_to_df[[j]] <- 1
      }else{
        out <- infinite_survival3(to = x_surplus[i], 
                                h_x = h_x, 
                                p = p,
                                lambda_true = lambda_true,
                                mean_true = mean_true,
                                f_ = f_combined,
                                l1 = l1, l2 = l2, nu = nu,a = k, b = beta, p_1 = p_1, p_1_0 = p_1_0, p_2 = p_2,
                                p_0_2 = p_0_2, p_1_2 = p_1_2)
        
        length_out <- length(out$V)
        print(thetas_ok[j])
        
        thetas_to_df[[j]] <- thetas_ok[j]
        value_to_df[[j]] <- out$V[length_out] 
      }
      
      if(verbose){
        pb$tick()
      }
      
    }
    
    df_dep_one[[i]] <- data.frame(x = unlist(thetas_to_df), value = unlist(value_to_df), surplus = paste(x_surplus[i]))
  }
  
  df_dep_one <- do.call(what = rbind, df_dep_one)
  
  
  
  # plot
  # gg <- plot_all(ruin_lines = df_dep_one, 
  #                profit_lines = data.frame(x = thetas, y = expected_income/max(expected_income), id = "Expected Profit"), 
  #                theta_lim = c(0,theta_max), 
  #                optimal_expected_profit = max(expected_income))
  
  surpluses <- unique(df_dep_one$surplus)
  optimal <- list()
  for( i in 1:length(surpluses)){
    
    optimal[[i]] <- thetas_ok[df_dep_one$value[df_dep_one$surplus == surpluses[i]] 
                              == min(df_dep_one$value[df_dep_one$surplus == surpluses[i]])]
    
  }
  optimal <- unlist(optimal)
  if(length(unique(optimal)) != 1) warning("more than one optimal value for ruin")
  
  
  ret <- list()
  # ret$plot <- gg
  ret$df <- df_dep_one
  ret$expected_proift <- expected_income
  ret$thetas <- thetas
  ret$thetas_ok <- thetas_ok
  ret$optimal_expected_profit <- max(expected_income)
  ret$optimal_theta <- theta_optimal
  ret$opt_ruin_theta <- optimal[1]
  return(ret)
  
  
}


#' @param claim_mean - claim mean of the individual claims
two_loadings_indep <- function(N, r, fixed_cost, lambda, k, beta, x_surplus, claim_mean, demand,  h_x = 10, S_, theta_finess = 0.01, verbose = TRUE, ...){
  
  
  
  # 
  #   N <- c(10000, 10000)
  #   k <- c(2,2)
  #   beta <- c(500,500)
  #   fixed_cost <-  c(100,100)
  #   lambda <- c(0.08, 0.08)
  #   r <- 0.05*k*beta*lambda # 5% of mean claims, minus a fixed cost
  #   nu <- 0.5
  #   demand <- function(theta1, theta2){
  #     b1 <- c(-0.5, 3)
  #     b2 <- c(-0.5, 4)
  #     return(c(1/(1+exp(b1[1]+b1[2]*theta1)), 1/(1+exp(b2[1]+b2[2]*theta2))))
  #   }
  
  
  thetas_1 <- seq(from = 0, to = 0.7, by = theta_finess )
  thetas_2 <- seq(from = 0, to = 0.7, by = theta_finess )
  expected_income <- list()
  for(i in 1:length(thetas_1)){
    expected_income[[i]] <- list()
  }
  
  # intesity for claim mean
  # l <- lambda*N*demand(thetas[i])
  #  claim_mean <- l*k*beta/sum(l)
  
  for(i in 1:length(thetas_1)){
    for(j in 1:length(thetas_2)){
      
      l <- N*demand(thetas_1[i], thetas_2[j])*lambda/sum(N*demand(thetas_1[i], thetas_2[j])*lambda)
      claim_mean_tmp <- sum(l*claim_mean)
      
      expected_income[[i]][[j]] <- sum((1+c(thetas_1[i], thetas_1[j]))*N*demand(thetas_1[i], thetas_2[j])*(lambda*k*beta-r)) - sum(fixed_cost) - claim_mean_tmp*sum(lambda*N*demand(thetas_1[i], thetas_2[j]))
      expected_income[[i]] <- unlist(expected_income[[i]])
      
    }
  }
  
  
  
  expected_income <- do.call(rbind, expected_income)
  
  if(verbose){
    pb <- progress_bar$new(format = "  downloading [:bar] :percent eta: :eta",
                           total = length(thetas_1)*length(thetas_2))
  }
  
  
  V <- matrix(NA, nrow = length(thetas_1), ncol = length(thetas_2))
  for(i in 1:length(thetas_1)){
    for(j in 1:length(thetas_2)){
      # print(paste0( j + (i-1)*length(thetas_2), " of ", length(thetas_2)*length(thetas_1)))
      
      l <- N*demand(thetas_1[i], thetas_2[j])*lambda/sum(N*demand(thetas_1[i], thetas_2[j])*lambda)
      
      S_tmp<- function(x, ...){
        
        return(sum(l*S_(x, ...)))
        
      }
      S_tmp <- Vectorize(S_tmp, vectorize.args = "x")
      
      lambda_true <- sum(lambda*N*demand(thetas_1[i], thetas_2[j]))
      mean_true <- sum(l*claim_mean)
      p <-  sum((1+c(thetas_1[i], thetas_2[j]))*N*demand(thetas_1[i], thetas_2[j])*(lambda*k*beta-r)) - sum(fixed_cost)
      
      if(p<mean_true*lambda_true){
        V[i,j] <- 1
      }else{
        
        
        out <- infinite_survival3(to = x_surplus, 
                                h_x = h_x, 
                                p = p,
                                lambda_true = lambda_true,
                                mean_true = mean_true,
                                S_ = S_tmp,
                                ...)
        V[i,j] <- out$V[length(out$V)]  
      }
      
      if(verbose){
        pb$tick()
      }
      
    }
  }
  
  
  V[is.na(V)] <- 99
  v_min <- min(V)
  t_opt <- which(V == v_min, arr.ind = TRUE)
  t1_opt <- thetas_1[t_opt[1]]
  t2_opt <- thetas_2[t_opt[2]]
  V[V == 99] <- NA
  
  
  expected_income[is.na(expected_income)] <- -99999
  profit_max <- max(expected_income)
  t_opt_prof <- which(expected_income == profit_max, arr.ind = TRUE)
  t1_opt_prof <- thetas_1[t_opt_prof[1]]
  t2_opt_prof <- thetas_2[t_opt_prof[2]]
  expected_income[expected_income == -99999] <- NA
  
  ret <- list()
  ret$V <- V
  ret$theta_1 <- thetas_1
  ret$theta_2 <- thetas_2
  ret$expected_income <- expected_income
  ret$t1_opt <- t1_opt
  ret$t2_opt <- t2_opt
  ret$t1_opt_prof <- t1_opt_prof
  ret$t2_opt_prof <- t2_opt_prof
  
  return(ret)
  
}


two_loadings_dep <- function(N, r, fixed_cost, lambda, k, beta, x_surplus, demand,  h_x = 10, 
                             theta_finess = 0.05, f_z_max = 10000, nu, f_z_limit = f_z_max, verbose = TRUE){
  
  
  
  # 
  #   N <- c(10000, 10000)
  #   k <- c(2,2)
  #   beta <- c(500,500)
  #   fixed_cost <-  c(100,100)
  #   lambda <- c(0.08, 0.08)
  #   r <- 0.05*k*beta*lambda # 5% of mean claims, minus a fixed cost
  #   nu <- 0.5
  #   demand <- function(theta1, theta2){
  #     b1 <- c(-0.5, 3)
  #     b2 <- c(-0.5, 4)
  #     return(c(1/(1+exp(b1[1]+b1[2]*theta1)), 1/(1+exp(b2[1]+b2[2]*theta2))))
  #   }
  
  
  
  S_common <- function(x1,x2,a,b, l1, l2, l_common, nu){
    
    tmp <- (((l1*(1-pgamma(q = x1, shape = a[1], scale = b[1])))^(-nu) + 
               (l2*(1-pgamma(q = x2, shape = a[2], scale = b[2])))^(-nu))^(-1/nu))/l_common
    return(tmp)
  }
  
  # F_common <- function(x1,x2, a,b,l1,l2, nu){
  #   l_common <- clayton_copula(l1,l2,nu)
  #   print(l_common)
  #   1-clayton_copula(l1*(1-pgamma(q = x1, shape = a[1], scale = b[1])), l2, nu)/l_common -
  #     clayton_copula(, l1, nu)/l_common +
  #     clayton_copula(l1*(1-pgamma(q = x1, shape = a[1], scale = b[1])), 0, nu)/l_common 
  #   
  # }
  # 
  # F_common(10,1000, c(2,2),c(5,5), 3,4, 0.5)
  # S_common(0,0,c(2,2), c(5,5), 2,3, clayton_copula(2,3,0.2), 0.2)
  
  # We use the clayton copula
  clayton_copula <- function(x1,x2, nu){
    tmp <- ((x1)^(-nu)+(x2)^(-nu))^(-1/nu)
    return(tmp)
  }
  
  
  # assume exponential margin, straight forward to adjust for more general case
  f_bivariate <- function(x1, x2, nu, a, b, l1, l2, l_common){
    
    copula_density <- function(u1, u2, nu){
      
      return((1+nu)*(u1^(-nu)+u2^(-nu))^(-1/nu - 2)*(u1^(-nu-1)*u2^(-nu-1)))
    }
    
    return(l1*l2*dgamma(x = x1, shape = a[1], scale = b[1])*dgamma(x = x2, shape = a[2], scale = b[2])*
             copula_density(u1 = l1*(1-pgamma(q = x1, shape = a[1], scale = b[1])), 
                            u2 = l2*(1-pgamma(q = x2, shape = a[2], scale = b[2])), nu)/l_common)
    
  }
  
  
  # Density of independent jumps
  f_indp_perp <- function(x, a,b, nu, l1, l2, l_common){
    l1_perp <- l1-l_common
    
    tmp1 <- l1*dgamma(x = x, shape = a, scale = b)/l1_perp
    tmp2 <- ((l1*(1-pgamma(q = x, shape = a, scale = b)))^(-nu) + l2^(-nu))^(-1/nu-1)*
      (l1*(1-pgamma(q = x, shape = a, scale = b)))^(-nu-1)*l1*dgamma(x = x, shape = a, scale = b)/l1_perp
    
    if(x< f_z_limit){
      return(tmp1 - tmp2)
    }else{
      return(0)
    }
    
  }
  
  f_indp_perp <- Vectorize(f_indp_perp, vectorize.args = "x")
  # density of the sum of common jumps is
  f_z <- function(z, a,b, nu, l1, l2, l_common){
    # print(z)  
    int_function <- function(x) f_bivariate(x, z-x, nu = nu, a = a, b = b, l1, l2, l_common)
    
    if(f_z_limit<z){
      return(0)
    }else{
      return(integrate(int_function, lower = 0, upper = z)$value) 
    }
  }
  
  f_z <- Vectorize(f_z, vectorize.args = "z")
  
  
  # Density of the marginal common jump
  
  f_1perp <- function(x, l1, l2, a,b, nu){
    
    l_common <- clayton_copula(l1,l2,nu)
    
    ((l1*(1-pgamma(q = x, shape = a, scale = b)))^(-nu) + l2^(-nu))^(-1/nu-1)*
      (l1*(1-pgamma(q = x, shape = a, scale = b)))^(-nu-1)*l1*dgamma(x = x, shape = a, scale = b)/l_common
    
  }
  
  f_1perp <- Vectorize(f_1perp, vectorize.args = "x")
  
  integrate(f_1perp, lower = 0, upper = 100, l1 = 3, l2 = 4, a = 2, b = 5, nu = 0.2)
  
  # 
  # f_indp_perp_tilde_1 <- function(x, p_1, p_1_0, a,b, nu, l1, l2, l_common){
  #   l1_perp <- l1-l_common
  #   w_indp <- (p_1*l1_perp)/(p_1*l1_perp + p_1_0*l_common)
  #   w_com <- (p_1_0*l_common)/(p_1*l1_perp + p_1_0*l_common)
  #   
  #   
  #   return(w_indp*f_indp_perp(x, a[1], b[1], nu, l1, l2, l_common) + w_com*f_z(x, a, b, nu, l1,l2, l_common))
  #   
  # }
  # f_indp_perp_tilde_1 <- Vectorize(f_indp_perp_tilde_1, vectorize.args = 'x' )
  # 
  # f_indp_perp_tilde_2 <- function(x, p_2, p_0_2, a,b, nu, l1, l2, l_common){
  #   l2_perp <- l2-l_common
  #   w_indp <- (p_2*l2_perp)/(p_2*l2_perp + p_0_2*l_common)
  #   w_com <- (p_0_2*l_common)/(p_2*l2_perp + p_0_2*l_common)
  #   
  #   
  #   return(w_indp*f_indp_perp(x, a[2], b[2], nu, l2, l1, l_common) + w_com*f_z(x, a, b, nu, l1,l2, l_common))
  #   
  # }
  # f_indp_perp_tilde_2 <- Vectorize(f_indp_perp_tilde_2, vectorize.args = 'x' )
  
  
  
  
  # f_combined <- function(x,l1_perp_tilde,l2_perp_tilde, l_common_tilde,  l1, l2  l_common, a,b,  nu){
  #   lambda <- l1_perp_tilde + l2_perp_tilde + l_common_tilde
  #   tmp1 <- l1_perp*f_indp_perp_tilde(x, p_1)/lambda
  #   tmp2 <- l2_perp*f_indp_perp(x = x, a = a[2], b = b[2], nu = nu, l1 = l2, l2 = l1, l_common = l_common)/lambda
  #   tmp3 <- l_common*f_z(z = x, a = a, b = b, nu = nu, l1 = l1, l2 = l2, l_common = l_common )/lambda
  #   return(tmp1 + tmp2 + tmp3)
  #   
  # }
  # f_combined <- Vectorize(f_combined, vectorize.args = "x")
  
  
  
  # f_combined <- function(x, l1_perp_tilde, l2_perp_tilde, l_common_tilde,  l1, l2,  l_common, a,b,  nu, p_1, p_1_0, p_2, p_0_2){
  #   
  #   lambda <- l1_perp_tilde + l2_perp_tilde + l_common_tilde
  #   
  #   tmp1 <- l1_perp_tilde*f_indp_perp_tilde_1(x = x, p_1 = p_1, p_1_0 = p_1_0, a = a, b = b, 
  #                                     nu = nu, l1 = l1, l2 = l2, l_common = l_common)/lambda
  #   
  #   tmp2 <- l2_perp_tilde*f_indp_perp_tilde_2(x = x, p_2 = p_2, p_0_2 = p_0_2, a = a, b = b, 
  #                                     nu = nu, l1 = l2, l2 = l1, l_common = l_common)/lambda
  #   
  #   tmp3 <- l_common_tilde*f_z(z = x, a = a, b = b, nu = nu, l1 = l1, l2 = l2, l_common = l_common )/lambda
  #   
  #   
  #   return(tmp1 + tmp2 + tmp3)
  # 
  # }
  
  f_combined <- function(x, l1, l2, nu,a, b, p_1, p_1_0, p_2, p_0_2){
    
    l_common <- clayton_copula(l1,l2,nu)
    l1_perp <- l1 - l_common
    l2_perp <- l2 - l_common
    
    l_tilde <- p_1*l1_perp + p_2*l2_perp + (p_1_0 + p_0_2 + p_1*p_2)*l_common
    
    
    p_1*l1_perp*f_indp_perp(x = x, a = a[1], b = b[1], nu = nu, l1 = l1, l2 = l2, l_common = l_common)/l_tilde + 
      p_2*l2_perp*f_indp_perp(x = x, a = a[2], b = b[2], nu = nu, l1 = l2, l2 = l1, l_common = l_common)/l_tilde +
      p_1_0*l_common*f_1perp(x = x, l1 = l1, l2 = l2, a = a[1], b = b[1], nu = nu)/l_tilde +
      p_0_2*l_common*f_1perp(x = x, l1 = l2, l2 = l1, a = a[2], b = b[2], nu = nu)/l_tilde +
      p_1*p_2*l_common*f_z(z = x, a = a, b = b, nu = nu, l1 = l1, l2 = l2, l_common = l_common )/l_tilde
    
    # p_1*l1_perp/l_tilde + 
    #   p_2*l2_perp/l_tilde +
    #   p_1_0*l_common/l_tilde +
    #   p_0_2*l_common/l_tilde +
    #   p_1*p_2*l_common/l_tilde
    
    
    
    
  }
  f_combined <- Vectorize(f_combined, vectorize.args = "x")
  
  
  
  mean_func <- function(x ,l1, l2, nu,a, b, p_1, p_1_0, p_2, p_0_2 ){
    x*f_combined(x, l1, l2, nu,a, b, p_1, p_1_0, p_2, p_0_2 )
  }
  
  
  
  
  thetas_1 <- seq(from = 0, to = 0.5, by = 0.01 )
  thetas_2 <- seq(from = 0, to = 0.5, by = 0.01 )
  expected_income <- list()
  for(i in 1:length(thetas_1)){
    expected_income[[i]] <- list()
  }
  
  # intesity for claim mean
  # l <- lambda*N*demand(thetas[i])
  #  claim_mean <- l*k*beta/sum(l)
  
  for(i in 1:length(thetas_1)){
    for(j in 1:length(thetas_2)){
      print(paste(i, j))
      
      demands <- demand(thetas_1[i], thetas_2[j])
      
      l1 <- N[1]*lambda[1]
      l2 <- N[2]*lambda[2]
      
      l_common <- clayton_copula(l1,l2,nu)
      l1_perp <- l1 - l_common
      l2_perp <- l2 - l_common
      
      # Define probabilities
      p_1 <- demands[1]
      p_2 <- demands[2]
      p_1_0 <- demands[1]*(1-demands[2])
      p_0_2 <- demands[2]*(1-demands[1])
      p_1_2 <- demands[1]*demands[2]
      
      
      # 
      # # Adjust 
      # l1_perp_tilde <- demands[1]*l1_perp + p_1_0*l_common
      # 
      # l2_perp_tilde <- demands[2]*l2_perp + p_0_2*l_common
      # l_common_tilde <- p_1_2*l_common 
      
      
      # f_combined(5000, l1_perp_tilde = l1_perp_tilde, l_common_tilde = l_common_tilde, l2_perp_tilde = l2_perp_tilde, l1 = l1, l2 = l2,
      #            l_common = l_common, a = a, b = b, nu = nu, p_1 = p_1, p_1_0 = p_1_0, p_2 = p_2, p_0_2 = p_0_2
      # )
      
      
      l_all <- p_1*l1_perp + p_2*l2_perp + (p_1_0 + p_0_2 + p_1*p_2)*l_common
      
      claim_mean_all <- integrate(mean_func, lower = 0, upper = f_z_max, 
                                  l1 = l1, l2 = l2, nu = nu,a = k, b = beta, 
                                  p_1 = p_1, p_1_0 = p_1_0, p_2 = p_2, p_0_2 = p_0_2)$value
      
      
      
      expected_income[[i]][[j]] <- sum((1+c(thetas_1[i], thetas_1[j]))*N*demand(thetas_1[i], thetas_2[j])*(lambda*k*beta)) - sum(fixed_cost) - claim_mean_all*l_all
      expected_income[[i]] <- unlist(expected_income[[i]])
      
    }
  }
  expected_income <- do.call(rbind, expected_income)
  theta_1_ei <- thetas_1
  theta_2_ei <- thetas_2
  
  thetas_1 <- c(seq(from = 0, to = 029, by = 0.05 ), 
                seq(from = 0.295, to = 0.35, by = 0.01), 
                seq(from = 0.355, to = 0.7, by = 0.05 ))
  thetas_2 <- c(seq(from = 0, to = 029, by = 0.05 ), 
                seq(from = 0.295, to = 0.35, by = 0.01), 
                seq(from = 0.355, to = 0.7, by = 0.05 ))
  
  
  
  if(verbose){
    pb <- progress_bar$new(format = "  downloading [:bar] :percent eta: :eta",
                           total = length(thetas_1)*length(thetas_2))
  }
  
  V <- matrix(NA, nrow = length(thetas_1), ncol = length(thetas_2))
  for(i in 1:length(thetas_1)){
    for(j in 1:length(thetas_2)){
      # print(paste0( j + (i-1)*length(thetas_2), " of ", length(thetas_2)*length(thetas_1)))
      
      demands <- demand(thetas_1[i], thetas_2[j])
      
      l1 <- N[1]*lambda[1]
      l2 <- N[2]*lambda[2]
      
      l_common <- clayton_copula(l1,l2,nu)
      l1_perp <- l1 - l_common
      l2_perp <- l2 - l_common
      
      # Define probabilities
      p_1 <- demands[1]
      p_2 <- demands[2]
      p_1_0 <- demands[1]*(1-demands[2])
      p_0_2 <- demands[2]*(1-demands[1])
      p_1_2 <- demands[1]*demands[2]
      
      l_all <- p_1*l1_perp + p_2*l2_perp + (p_1_0 + p_0_2 + p_1*p_2)*l_common
      
      claim_mean_all <- integrate(mean_func, lower = 0, upper = f_z_max, 
                                  l1 = l1, l2 = l2, nu = nu,a = k, b = beta, 
                                  p_1 = p_1, p_1_0 = p_1_0, p_2 = p_2, p_0_2 = p_0_2)$value
      
      
      lambda_true <-  l_all
      p <-  sum((1+c(thetas_1[i], thetas_2[j]))*N*demand(thetas_1[i], thetas_2[j])*(lambda*k*beta)) - sum(fixed_cost)
      
      
      if(p<claim_mean_all*lambda_true){
        V[i,j] <- 1
      }else{
        
        out <- infinite_survival3(to = 5000, 
                                h_x = 500, 
                                p = p,
                                lambda_true = lambda_true,
                                mean_true = claim_mean_all,
                                f_ = f_combined,
                                l1 = l1, l2 = l2, nu = nu,a = k, b = beta, p_1 = p_1, p_1_0 = p_1_0, p_2 = p_2, p_0_2 = p_0_2)
        
        V[i,j] <- out$V[length(out$V)]  
      }
      
      if(verbose){
        pb$tick()
      }
      
    }
  }
  
  
  V[is.na(V)] <- 99
  v_min <- min(V)
  t_opt <- which(V == v_min, arr.ind = TRUE)
  t1_opt <- thetas_1[t_opt[1]]
  t2_opt <- thetas_2[t_opt[2]]
  V[V == 99] <- NA
  
  
  expected_income[is.na(expected_income)] <- -99999
  profit_max <- max(expected_income)
  t_opt_prof <- which(expected_income == profit_max, arr.ind = TRUE)
  t1_opt_prof <- theta_1_ei[t_opt_prof[1]]
  t2_opt_prof <- theta_1_ei[t_opt_prof[2]]
  expected_income[expected_income == -99999] <- NA
  
  
  
  ret <- list()
  ret$V <- V
  ret$theta_1 <- thetas_1
  ret$theta_2 <- thetas_2
  ret$t1_opt <- t1_opt
  ret$t2_opt <- t2_opt
  ret$expected_income <- expected_income
  ret$theta_1_ei <- theta_1_ei
  ret$theta_2_ei <- theta_2_ei
  ret$t1_opt_prof <- t1_opt_prof
  ret$t2_opt_prof <- t2_opt_prof
  
  return(ret)
  
}
