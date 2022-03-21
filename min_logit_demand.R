library(tidyverse)
library(ggnewscale) # add new color scale -> new_scale_color
# Set color brewer
library(RColorBrewer)
library(plotly)


source("utils.R")  # functions
source("Parameters.R")  # Parameters used in numerical analysis


# First Marginal ----

out_1 <- one_loading_inference_indp(N = N[1], 
                                    r = r[1], 
                                    fixed_cost = fixed_cost[1],
                                    h_x = 100,
                                    lambda = lambda[1], 
                                    claim_mean = k[1]*beta[1],
                                    x_surplus = x_surplus, 
                                    demand = demand_1, 
                                    theta_finess = 0.01,
                                    S_ = S_marginal, 
                                    shape = k[1], 
                                    scale = beta[1])

# Fin theta which minimizes the the first marginal
v_min_1 <- out_1$df %>% group_by(surplus) %>% summarise(V_min = min(value), theta_min = x[which(value == min(value))])
#out_1$plot



plot_1 <- plot_all(ruin_lines = out_1$df, 
         profit_lines = data.frame(x = out_1$thetas, y = out_1$expected_proift/out_1$optimal_expected_profit, id = "Expected Profit"), 
         optimal_expected_profit = out_1$optimal_expected_profit, 
         optimal_theta = out_1$optimal_theta, theta_lim = c(0, 0.75), 
         y_lim = c(0,1.07),
         opt_ruin_theta = out_1$opt_ruin_theta, text_offset = 0.05,
         breaks = c(0, 5000, 10000, 15000, 20000, 30000, round(out_1$optimal_expected_profit)),
         size = 5.5)

my_blue <- colorRampPalette(brewer.pal(n = 9, "Blues")[3:7])(8)[1:8]
v_min_1$id <- as.character(round(v_min_1$V_min , digits = 3))
v_min_1$id <- factor(v_min_1$id, levels = v_min_1$id)
plot_1 + new_scale_color() +
  geom_point(aes(x = theta_min, y = V_min, color = id, digits = 3), 
             v_min_1, size = 3) +
  scale_color_manual(values = my_blue, name = "Min. Ruin Prob.")

# Second Marginal ----


out_2 <- one_loading_inference_indp(N = N[2], 
                                    r = r[2], 
                                    fixed_cost = fixed_cost[2],
                                    lambda = lambda[2], 
                                    h_x = 100,
                                    claim_mean = k[2]*beta[2], 
                                    theta_finess = 0.005,
                                    x_surplus = x_surplus, 
                                    demand = demand_2, 
                                    S_ = S_marginal, 
                                    shape = k[2], 
                                    scale = beta[2])

v_min_2 <- out_2$df %>% group_by(surplus) %>% 
  summarise(V_min = min(value), theta_min = x[which(value == min(value))])


plot_2 <- plot_all(ruin_lines = out_2$df, 
         profit_lines = data.frame(x = out_2$thetas, y = out_2$expected_proift/out_2$optimal_expected_profit, id = "Expected Profit"), 
         theta_lim = c(0, 0.55), 
         optimal_expected_profit = out_2$optimal_expected_profit, 
         optimal_theta = out_2$optimal_theta,
         opt_ruin_theta = out_2$opt_ruin_theta, text_offset = 0.04, 
         breaks = c(0, 5000, 10000, 20000, 30000, round(out_2$optimal_expected_profit)),
         y_lim = c(0,1.07),size = 5.5)

v_min_2$id <- as.character(round(v_min_2$V_min , digits = 3))
v_min_2$id <- factor(v_min_2$id, levels = v_min_2$id)
plot_2 + new_scale_color() +
  geom_point(aes(x = theta_min, y = V_min, color = id, digits = 3), 
             v_min_2, size = 3) +
  scale_color_manual(values = my_blue, name = "Min. Ruin Prob.")




# Joint but indpendent ----


out_one_indp <- one_loading_inference_indp(N = N, 
                                           r = r, 
                                           fixed_cost = fixed_cost,
                                           lambda = lambda, 
                                           claim_mean = k*beta,
                                           demand = demand, 
                                           h_x = 100,
                                           x_surplus = x_surplus,
                                           theta_finess = 0.01,
                                           S_ = S_, 
                                           a = k, 
                                           b = beta)

v_min_indp <- out_one_indp$df %>% group_by(surplus) %>% 
  summarise(V_min = min(value), theta_min = x[which(value == min(value))])

  plot_indp <- plot_all(ruin_lines = out_one_indp$df, 
                   profit_lines = data.frame(x = out_one_indp$thetas, y = out_one_indp$expected_proift/out_one_indp$optimal_expected_profit, id = "Expected Profit"), 
                   theta_lim = c(0, 0.6), 
                   optimal_expected_profit = out_one_indp$optimal_expected_profit, 
                   optimal_theta = out_one_indp$optimal_theta,
                   opt_ruin_theta = out_one_indp$opt_ruin_theta,
                   text_offset = 0.04, 
                   breaks = c(0, 10000, 20000, 30000, 40000, 50000, out_one_indp$optimal_expected_profit), 
                   size = 5.5, y_lim = c(0, 1.08))
plot_indp
v_min_indp$id <- as.character(round(v_min_indp$V_min , digits = 3))
v_min_indp$id <- factor(v_min_indp$id, levels = v_min_indp$id)
plot_indp + new_scale_color() +
  geom_point(aes(x = theta_min, y = V_min, color = id, digits = 3), 
             v_min_indp, size = 3) +
  scale_color_manual(values = my_blue, name = "Min. Ruin Prob.")



# compare ruin thetas ----
out_1$opt_ruin_theta
out_2$opt_ruin_theta
out_one_indp$opt_ruin_theta


# Joint with Levy copula dependence ----
# These take some time to calculate, The higher the surplus, the higher the complexity.


  
out_both_one_dep_1 <- one_loading_inference_dep(N = N, 
                                              r = r, 
                                              fixed_cost = fixed_cost, 
                                              lambda = lambda, 
                                              nu = nu, 
                                              k = k, 
                                              beta = beta,
                                              x_surplus = x_surplus[1],
                                              demand = demand, 
                                              h_x = 100, 
                                              f_z_max = 20000, 
                                              theta_finess = 0.01,
                                              f_z_limit = 20000)


 out_both_one_dep_2 <- one_loading_inference_dep(N = N, 
                                                r = r, 
                                                fixed_cost = fixed_cost, 
                                                lambda = lambda, 
                                                nu = nu, 
                                                k = k, 
                                                beta = beta,
                                                x_surplus = x_surplus[2],
                                                demand = demand, 
                                                h_x = 100, 
                                                f_z_max = 20000, 
                                                theta_finess = 0.01)
 
save.image("./Data_clayton.RData")

out_both_one_dep_3 <- one_loading_inference_dep(N = N, 
                                                r = r, 
                                                fixed_cost = fixed_cost, 
                                                lambda = lambda, 
                                                nu = nu, 
                                                k = k, 
                                                beta = beta,
                                                x_surplus = x_surplus[3],
                                                demand = demand, 
                                                h_x = 100, 
                                                f_z_max = 20000, 
                                                theta_finess = 0.02)

save.image("./Data_clayton.RData")

out_both_one_dep_4 <- one_loading_inference_dep(N = N, 
                                                r = r, 
                                                fixed_cost = fixed_cost, 
                                                lambda = lambda, 
                                                nu = nu, 
                                                k = k, 
                                                beta = beta,
                                                x_surplus = 20000,
                                                demand = demand, 
                                                h_x = 500, 
                                                f_z_max = 10000, 
                                                theta_finess = 0.02)

save.image("./Data_clayton.RData")

#save(out_both_one_dep_3, file = "dep_3.RData")

df_dep <- rbind(out_both_one_dep_1$df, out_both_one_dep_2$df, out_both_one_dep_3$df, out_both_one_dep_4$df)
profit_lines_dep <-  data.frame(x = out_both_one_dep_1$thetas, y = out_both_one_dep_1$expected_proift/out_one_indp$optimal_expected_profit, id = "Expected Profit Dependence")
profit_lines_indp <- data.frame(x = out_one_indp$thetas, y = out_one_indp$expected_proift/out_one_indp$optimal_expected_profit, id = "Expected Profit Independence")
v_min_dep <- df_dep %>% group_by(surplus) %>% 
  summarise(V_min = min(value), theta_min = x[which(value == min(value))])

plot_dep <- plot_all(ruin_lines = df_dep, 
                     profit_lines = profit_lines_dep,
                     theta_lim = c(0,0.6),
                     optimal_expected_profit = out_both_one_dep_1$optimal_expected_profit, 
                      optimal_theta = out_both_one_dep_1$optimal_theta,
                     opt_ruin_theta = out_both_one_dep_1$opt_ruin_theta, 
                     text_offset = 0.07,
                     breaks = c(0, 10000, 20000, 30000, 40000, 50000, out_both_one_dep_1$optimal_expected_profit),
                    size = 5.5, y_lim = c(0, 1.08))

my_blue <- colorRampPalette(brewer.pal(n = 9, "Blues")[3:7])(8)[1:8]
v_min_dep$id <- as.character(round(v_min_dep$V_min , digits = 3))
v_min_dep$id <- factor(v_min_dep$id, levels = v_min_dep$id)
plot_dep + new_scale_color() +
  geom_point(aes(x = theta_min, y = V_min, color = id, digits = 3), 
             v_min_dep, size = 3) +
  scale_color_manual(values = my_blue, name = "Min. Ruin Prob.")
  



df_dep_tmp <- df_dep
df_dep_tmp$surplus <- as.character(df_dep_tmp$surplus)
surplus_tmp <- unique(df_dep_tmp$surplus)
for( j in 1:length(surplus_tmp)){
  tmp_min <- round(v_min_dep$V_min[v_min_dep$surplus == surplus_tmp[j]], digits = 2)
  df_dep_tmp$surplus[df_dep_tmp$surplus == surplus_tmp[j]] <-  paste0(surplus_tmp[j], " -> ", tmp_min, sep = "")
}
df_dep_tmp$surplus <- factor(df_dep_tmp$surplus, 
                              levels =c("100 -> 0.92",
                                        "5000 -> 0.61", 
                                        "15000 -> 0.26",
                                        "20000 -> 0.15") )

df_indp_tmp <- out_one_indp$df
df_indp_tmp$surplus <- as.character(df_indp_tmp$surplus)
surplus_tmp <- unique(df_indp_tmp$surplus)
for( j in 1:length(surplus_tmp)){
  tmp_min <- round(v_min_indp$V_min[v_min_indp$surplus == surplus_tmp[j]], digits = 2)
  df_indp_tmp$surplus[df_indp_tmp$surplus == surplus_tmp[j]] <-  paste0(surplus_tmp[j], " -> ", tmp_min, sep = "")
}
df_indp_tmp$surplus <- factor(df_indp_tmp$surplus, 
                             levels =c("100 -> 0.92",
                                       "5000 -> 0.57", 
                                       "15000 -> 0.21",
                                       "20000 -> 0.13",
                                       "30000 -> 0.05",
                                       "40000 -> 0.02") )

# my_blue = brewer.pal(n = 9, "Blues")[2:9] 
my_blue <- colorRampPalette(brewer.pal(n = 9, "Blues")[3:7])(8)[1:8]
my_red <- colorRampPalette(brewer.pal(n = 9, "Reds")[4:8])(8)[1:8]
my_grey <- brewer.pal(n = 9, "Greys")[9] 

gg <- ggplot() +
  
  geom_line(aes(x = x, y = y, color = id), profit_lines_dep, size = 1, linetype = "dashed") +
  scale_color_manual(values = "#882020", name = "") +
  new_scale_color() + 
  geom_line(aes(x = x, y = y, color = id), profit_lines_indp, size = 0.5, linetype = "dashed") +
  # geom_point(mapping = aes(x = x, y = y, color = id), data = profit_point, size = 2, shape = 15, alpha = 0.6) +
  #new_scale_color() +
  scale_color_manual(values = "#36454f", name = "") +
  new_scale_color()+
  
  geom_line(aes(x = x, y = value, color = surplus), df_indp_tmp, size = 1) +
  scale_color_manual(values = my_blue[c(1,3,4,5,6,8)], name = "Surplus, Ind. Claims") +
  new_scale_color() +
  
  
  # geom_point(mapping = aes(x = x, y = value, color = surplus), ruin_point, size = 2) +
  # new_scale_color() + 
  # geom_line(aes(x = x, y = y, color = id), out_one_indp$expected_proift, size = 1, alpha = 0.6)+
  # # geom_point(mapping = aes(x = x, y = y, color = id), data = profit_point, size = 2, shape = 15, alpha = 0.6) +
  # scale_color_manual(values = my_grey) +
  
  geom_line(aes(x = x, y = value, color = surplus), df_dep_tmp, size = 1) +
  scale_color_manual(values = my_red[c(1,3,5)], name = "Surplus, Dep. Claims" )  +
  
  # geom_point(mapping = aes(x = x, y = value, color = surplus), ruin_point, size = 2) +
  # new_scale_color() + 
  # geom_line(aes(x = x, y = y, color = id), out_both_one_dep$expected_proift, size = 1, alpha = 0.6)+
  # # geom_point(mapping = aes(x = x, y = y, color = id), data = profit_point, size = 2, shape = 15, alpha = 0.6) +
  # #new_scale_color() +
  # scale_color_manual(values = my_green) +
  
  
  # geom_hline(yintercept = 0, color = "black", size = 1) +
  # geom_vline(xintercept = 0, color = "black", size = 1) +
  # # geom_vline(xintercept = optimal_theta, color = "black", size = 0.5) +
  # annotate("text", x = out_both_one_dep_1$optimal_theta-0.03 , y = 1.03,
  #          label = bquote(theta[profit]^"*" == .(out_both_one_dep_1$optimal_theta)), size = 5.5)  +
  # annotate("text", x = out_both_one_dep_1$opt_ruin_theta + 0.03 , y = 1.03, 
  #          label = bquote(theta[ruin_dep]^"*" == .(out_both_one_dep_1$opt_ruin_theta)), size = 5.5) +
  # geom_segment(aes(x = out_both_one_dep_1$optimal_theta, xend = out_both_one_dep_1$optimal_theta, y = 0, yend = 1), linetype = 2) +
  # geom_segment(aes(x = out_both_one_dep_1$opt_ruin_theta, xend = out_both_one_dep_1$opt_ruin_theta, y = 0, yend = 1), linetype = 2) +
  labs(x = expression(theta), y = "Probability of Ruin \n", color = "") +
  #ggtitle("") +
  theme_bw(base_size = 20) +
  scale_x_continuous(expand = c(0.0015, 0),limits = c(0.15, 0.6)) +
   scale_y_continuous(expand = c(0.002, 0), 
                      limits = c(0,1.07),
                      sec.axis = sec_axis(~.*out_both_one_dep_1$optimal_expected_profit, 
                                          name = "Expected Profit \n", breaks = c(0, 10000, out_both_one_dep_1$optimal_expected_profit)), label = scales::comma) +
  theme(legend.justification = "top") 
gg




# Both, two, indp, ----  
  


S_ <- function(x, a, b){
  
  return(c(1-pgamma(q = x , shape = a[1], scale = b[1]), 1-pgamma(q = x , shape = a[2], scale = b[2])))
}
  
out_two_loading_indp_15000 <-  two_loadings_indep(N = N, r = r, fixed_cost = fixed_cost, lambda = lambda, k = k, beta = beta, claim_mean = k*beta, demand = demand,
                                               x_surplus = 15000, h_x = 100, S_ = S_, a = k, b = beta)
out_two_loading_indp_2000 <-  two_loadings_indep(N = N, r = r, fixed_cost = fixed_cost, lambda = lambda, k = k, 
                                                 beta = beta, claim_mean = k*beta, demand = demand,
                                                 x_surplus = 5000, h_x = 100, S_ = S_, a = k, b = beta)
out_two_loading_indp_30000 <-  two_loadings_indep(N = N, r = r, fixed_cost = fixed_cost, lambda = lambda, k = k, 
                                                  beta = beta, claim_mean = k*beta, demand = demand,
                                                  x_surplus = 30000, h_x = 100, S_ = S_, a = k, b = beta)
  
  
  
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
                                                 range = c(0, 0.8)),
                                    xaxis = list(title = "&#952;<sub>2</sub>", 
                                                 tickfont = list(size = 22), 
                                                 titlefont = list(size = 24),
                                                 range = c(0, 0.8)))
    
    
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

  
  

# 2000
plot_ruin(out_two_loading_indp_2000$V, out_two_loading_indp_2000$theta_1, out_two_loading_indp_2000$theta_2,
          type = "contour", t1_opt = out_two_loading_indp_2000$t1_opt, t2_opt = out_two_loading_indp_2000$t2_opt, 
          v_min = round(min(out_two_loading_indp_2000$V[!is.na(out_two_loading_indp_2000$V)]), digits = 2))

plot_profit(out_two_loading_indp_2000$expected_income, 
            out_two_loading_indp_2000$theta_1,
            out_two_loading_indp_2000$theta_2,
            type = "contour", 
            t1_opt = out_two_loading_indp_2000$t1_opt_prof, 
            t2_opt = out_two_loading_indp_2000$t2_opt_prof,
            prof_max =  round(max(out_two_loading_indp_2000$expected_income[!is.na(out_two_loading_indp_2000$expected_income)]), digits = 0))
# 15000
plot_ruin(out_two_loading_indp_15000$V, out_two_loading_indp_15000$theta_1, out_two_loading_indp_15000$theta_2,
          type = "contour", t1_opt = out_two_loading_indp_15000$t1_opt, t2_opt = out_two_loading_indp_15000$t2_opt, 
          v_min = round(min(out_two_loading_indp_15000$V[!is.na(out_two_loading_indp_15000$V)]), digits = 2))

plot_profit(out_two_loading_indp_15000$expected_income, 
            out_two_loading_indp_15000$theta_1,
            out_two_loading_indp_15000$theta_2,
            type = "contour", 
            t1_opt = out_two_loading_indp_15000$t1_opt_prof, 
            t2_opt = out_two_loading_indp_15000$t2_opt_prof,
            prof_max =  round(max(out_two_loading_indp_15000$expected_income[!is.na(out_two_loading_indp_15000$expected_income)]), digits = 0))
# 30000
plot_ruin(out_two_loading_indp_30000$V, out_two_loading_indp_30000$theta_1, out_two_loading_indp_30000$theta_2,
          type = "contour", t1_opt = out_two_loading_indp_30000$t1_opt, t2_opt = out_two_loading_indp_30000$t2_opt, 
          v_min = round(min(out_two_loading_indp_30000$V[!is.na(out_two_loading_indp_30000$V)]), digits = 2))

plot_profit(out_two_loading_indp_30000$expected_income, 
            out_two_loading_indp_30000$theta_1,
            out_two_loading_indp_30000$theta_2,
            type = "contour", 
            t1_opt = out_two_loading_indp_30000$t1_opt_prof, 
            t2_opt = out_two_loading_indp_30000$t2_opt_prof,
            prof_max =  round(max(out_two_loading_indp_30000$expected_income[!is.na(out_two_loading_indp_30000$expected_income)]), digits = 0))

  
# Both two  dependent -----
  
out_two_loading_dep_15000 <-  two_loadings_dep(N = N, r = r, fixed_cost = fixed_cost, lambda = lambda, k = k, beta = beta, 
                                             x_surplus = 5000, demand = demand, h_x = 100, theta_finess = 0.05,
                                             f_z_max = 20000, nu = nu)

plot_ruin(out_two_loading_dep_15000$V, out_two_loading_dep_15000$theta_1, out_two_loading_dep_15000$theta_2,
          type = "contour", t1_opt = out_two_loading_dep_15000$t1_opt, t2_opt = out_two_loading_dep_15000$t2_opt, 
          v_min = round(min(out_two_loading_dep_15000$V[!is.na(out_two_loading_dep_15000$V)]), digits = 2))

plot_profit(out_two_loading_dep_15000$expected_income, 
            out_two_loading_dep_15000$theta_1_ei,
            out_two_loading_dep_15000$theta_2_ei,
            type = "contour", 
            t1_opt = out_two_loading_dep_15000$t1_opt_prof, 
            t2_opt = out_two_loading_dep_15000$t2_opt_prof,
            prof_max =  round(max(out_two_loading_dep_15000$expected_income[!is.na(out_two_loading_dep_15000$expected_income)]), digits = 0))



#keep_logit <- out_1
#save(keep_logit, file = "keep_logit.RData")




# Test ----

ruin_exp_multi <- function(theta, a, b, r, claim_mean, surplus){
  
  
  lambda <- sum(a)- sum(b*theta)
  mu <- sum((a-b*theta)*claim_mean)/lambda
  p <- sum((1+theta)*(a-b*theta)*claim_mean)-sum(r)
  const <- lambda*mu/(p)
  R <- (1/mu-lambda/p)
  
  ruin <- const*exp(-R*surplus)
  
  return(ruin)
}




surplus <- seq(ffrom = 0, to = 20000, by = 100)
theta <- 0.4


  #l <- N_1*demand_1(thetas_ok[i])*lambda_1/sum(N*demand(thetas_ok[i])*lambda)
  mean_true <-  1500
  lambda_true <-  sum(lambda_1*N_1*demand_1(theta))
  p <-  (1+theta)*sum(N_1*demand_1(theta)*(lambda_1*mean_true-r_1)) - sum(fixed_cost_1)
  
  
  
  const <- lambda_true*mean_true/(p)
  R <- (1/mean_true-lambda_true/p)
  ruin_1 <- const*exp(-R*surplus)
  
  
  
  l <-  N*demand(theta)*lambda/sum(N*demand(theta)*lambda)
  lambda_true <-  sum(lambda*N*demand(theta))
  p <-  (1+theta)*sum(N*demand(theta)*(lambda*mean_true-r)) - sum(fixed_cost)
  
  
  
  const <- lambda_true*mean_true/(p)
  R <- (1/mean_true-lambda_true/p)
  ruin_2 <- const*exp(-R*surplus)


  plot(surplus, ruin_1, type = "l", ylim = c(0,1))
  lines(surplus, ruin_2)
  
  

# Save ----
  
  
   load("min_logit.RData")
  
  
  #save.image("min_logit.RData")






