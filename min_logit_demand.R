library(tidyverse)
library(ggnewscale) # add new color scale -> new_scale_color
# Set color brewer
library(RColorBrewer)
library(plotly)
library(progress)


source("utils.R")  # functions
source("parameters_diff_demand.R")  # Parameters used in numerical analysis


file = "./Data_clayton_1_diff_loadings_new.RData"
print(nu)
h_x <- 0.1
# First Marginal ----


out_1 <- one_loading_inference_indp(N = N[1], 
                                    r = r[1], 
                                    fixed_cost = fixed_cost[1],
                                    h_x = h_x,
                                    lambda = lambda[1], 
                                    claim_mean = k[1]*beta[1],
                                    x_surplus = x_surplus, 
                                    demand = demand_1, 
                                    theta_finess = 0.001,
                                    infinite_survival = infinite_survival1,
                                    theta_grid_ruin = theta_grid_ruin,
                                    S_ = S_marginal, 
                                    shape = k[1], 
                                    scale = beta[1])


# Fin theta which minimizes the the first marginal
v_min_1 <- out_1$df %>% group_by(surplus) %>% summarise(V_min = min(value), theta_min = x[which(value == min(value))])
#out_1$plot

v_min_1 <- v_min_1[order(v_min_1$V_min, decreasing = TRUE),]



plot_1 <- plot_all(ruin_lines = out_1$df, 
         profit_lines = data.frame(x = out_1$thetas, y = out_1$expected_proift/out_1$optimal_expected_profit, id = "Expected Profit"), 
         optimal_expected_profit = out_1$optimal_expected_profit, 
         optimal_theta = out_1$optimal_theta, 
         theta_lim = c(0, 0.5), 
         y_lim = c(0,1.07),
         opt_ruin_theta = round(out_1$opt_ruin_theta,3), 
         text_offset = 0.05,
         breaks = c(0, 50, 100, 150, 200, 250, 300, round(out_1$optimal_expected_profit), 350),
         size = 5.5)

my_blue <- colorRampPalette(brewer.pal(n = 9, "Blues")[3:7])(8)[c(1, 4, 8)]
v_min_1$id <- as.character(round( v_min_1$V_min , digits = 3))
v_min_1$id <- factor(v_min_1$id, levels = v_min_1$id)
plot_1 + new_scale_color() +
  geom_point(aes(x = theta_min, y = V_min, color = id, digits = 3), 
             v_min_1, size = 3) +
  scale_color_manual(values = my_blue, name = "Min. Ruin Prob.") +theme_bw(20)

# Second Marginal ----


out_2 <- one_loading_inference_indp(N = N[2], 
                                    r = r[2], 
                                    fixed_cost = fixed_cost[2],
                                    lambda = lambda[2], 
                                    h_x = h_x,
                                    claim_mean = k[2]*beta[2], 
                                    theta_finess = 0.001,
                                    x_surplus = x_surplus, 
                                    demand = demand_2, 
                                    S_ = S_marginal,
                                    infinite_survival = infinite_survival3,
                                    theta_grid_ruin = theta_grid_ruin,
                                    shape = k[2], 
                                    scale = beta[2])

v_min_2 <- out_2$df %>% group_by(surplus) %>% 
  summarise(V_min = min(value), theta_min = x[which(value == min(value))])
v_min_2 <- v_min_2[order(v_min_2$V_min, decreasing = TRUE),]

plot_2 <- plot_all(ruin_lines = out_2$df, 
         profit_lines = data.frame(x = out_2$thetas, y = out_2$expected_proift/out_2$optimal_expected_profit, id = "Expected Profit"), 
         theta_lim = c(0, 0.55), 
         optimal_expected_profit = out_2$optimal_expected_profit, 
         optimal_theta = out_2$optimal_theta,
         y_lim = c(0,1.07),
         opt_ruin_theta = round(out_2$opt_ruin_theta,3), 
         text_offset = 0.04, 
         breaks = c(0, 50, 100, 150, 200, 250, 300, round(out_2$optimal_expected_profit), 350),
         size = 5.5)

my_blue <- colorRampPalette(brewer.pal(n = 9, "Blues")[3:7])(8)[c(1, 4, 8)]
v_min_2$id <- as.character(round( v_min_2$V_min , digits = 3))
v_min_2$id <- factor(v_min_2$id, levels = v_min_2$id)
plot_2 + new_scale_color() +
  geom_point(aes(x = theta_min, y = V_min, color = id, digits = 3), 
             v_min_2, size = 3) +
  scale_color_manual(values = my_blue, name = "Min. Ruin Prob.") +theme_bw(20)




# Joint but indpendent ----


out_one_indp <- one_loading_inference_indp(N = N, 
                                           r = r, 
                                           fixed_cost = fixed_cost,
                                           lambda = lambda, 
                                           claim_mean = k*beta,
                                           demand = demand, 
                                           h_x = h_x,
                                           x_surplus = x_surplus,
                                           theta_finess = 0.001,
                                           theta_grid_ruin = theta_grid_ruin,
                                           infinite_survival = infinite_survival3,
                                           S_ = S_, 
                                           a = k, 
                                           b = beta)

v_min_indp <- out_one_indp$df %>% group_by(surplus) %>% 
  summarise(V_min = min(value), theta_min = x[which(value == min(value))])

  plot_indp <- plot_all(ruin_lines = out_one_indp$df, 
                   profit_lines = data.frame(x = out_one_indp$thetas, y = out_one_indp$expected_proift/out_one_indp$optimal_expected_profit, id = "Expected Profit"), 
                   theta_lim = c(0, 0.6), 
                   optimal_expected_profit = out_one_indp$optimal_expected_profit, 
                   optimal_theta = round(out_one_indp$optimal_theta, 2),
                   opt_ruin_theta = round(out_one_indp$opt_ruin_theta,2),
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

  
out_both_one_dep_1 <- one_loading_inference_dep2(N = N, 
                                              r = r, 
                                              fixed_cost = fixed_cost, 
                                              lambda = lambda, 
                                              nu = nu, 
                                              k = k, 
                                              beta = beta,
                                              x_surplus = x_surplus[1],
                                              demand = demand,
                                              theta_grid_ruin = theta_grid_ruin,
                                              infinite_survival = infinite_survival3,
                                              h_x = h_x, 
                                              f_z_max = 30, 
                                              theta_finess = 0.01,
                                              f_z_limit = 30)

save.image(file)
 out_both_one_dep_2 <- one_loading_inference_dep2(N = N, 
                                                r = r, 
                                                fixed_cost = fixed_cost, 
                                                lambda = lambda, 
                                                nu = nu, 
                                                k = k, 
                                                beta = beta,
                                                x_surplus = x_surplus[2],
                                                theta_grid_ruin = theta_grid_ruin,
                                                infinite_survival = infinite_survival3,
                                                demand = demand, 
                                                h_x = h_x, 
                                                f_z_max = 30, 
                                                theta_finess = 0.01)
 
save.image(file)

out_both_one_dep_3 <- one_loading_inference_dep2(N = N, 
                                                r = r, 
                                                fixed_cost = fixed_cost, 
                                                lambda = lambda, 
                                                nu = nu, 
                                                k = k, 
                                                beta = beta,
                                                x_surplus = x_surplus[3],
                                                theta_grid_ruin = theta_grid_ruin,
                                                infinite_survival = infinite_survival3,
                                                demand = demand, 
                                                h_x = h_x, 
                                                f_z_max = 30, 
                                                theta_finess = 0.01)

save.image(file)



#save(out_both_one_dep_3, file = "dep_3.RData")

df_dep <- rbind(out_both_one_dep_1$df, out_both_one_dep_2$df, out_both_one_dep_3$df)
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
  tmp_min_theta <- round(v_min_dep$theta_min[v_min_dep$surplus == surplus_tmp[j]], digits = 3)
  df_dep_tmp$surplus[df_dep_tmp$surplus == surplus_tmp[j]] <-  paste0(surplus_tmp[j], " -> ", tmp_min, ", \U03B8 = ", tmp_min_theta,  sep = "")
}
df_dep_tmp$surplus <- factor(df_dep_tmp$surplus, 
                              levels =c("2 -> 0.78, \U03B8 = 0.297",
                                        "5 -> 0.65, \U03B8 = 0.3", 
                                        "10 -> 0.48, \U03B8 = 0.3") )


df_indp_tmp <- out_one_indp$df
df_indp_tmp$surplus <- as.character(df_indp_tmp$surplus)
surplus_tmp <- unique(df_indp_tmp$surplus)
for( j in 1:length(surplus_tmp)){
  tmp_min <- round(v_min_indp$V_min[v_min_indp$surplus == surplus_tmp[j]], digits = 2)
  tmp_min_theta <- round(v_min_indp$theta_min[v_min_indp$surplus == surplus_tmp[j]], digits = 3)
  df_indp_tmp$surplus[df_indp_tmp$surplus == surplus_tmp[j]] <-  paste0(surplus_tmp[j], " -> ", tmp_min, ", \U03B8 = ", tmp_min_theta,  sep = "")
}
df_indp_tmp$surplus <- factor(df_indp_tmp$surplus, 
                             levels =c("2 -> 0.77, \U03B8 = 0.292",
                                       "5 -> 0.63, \U03B8 = 0.292", 
                                       "10 -> 0.45, \U03B8 = 0.292") )

# my_blue = brewer.pal(n = 9, "Blues")[2:9] 
my_blue <- colorRampPalette(brewer.pal(n = 9, "Blues")[3:7])(8)[c(1,4,8)]
my_red <- colorRampPalette(brewer.pal(n = 9, "Reds")[4:8])(8)[c(1,4,8)]
my_grey <- brewer.pal(n = 9, "Greys")[9] 

gg <- ggplot() +
  
  geom_line(aes(x = x, y = y), profit_lines_dep, size = 1, linetype = "dashed") +
  scale_color_manual(values = "#882020", name = "") +
  new_scale_color() + 
  geom_line(aes(x = x, y = y), profit_lines_indp, size = 0.5, linetype = "dashed") +
  # geom_point(mapping = aes(x = x, y = y, color = id), data = profit_point, size = 2, shape = 15, alpha = 0.6) +
  #new_scale_color() +
  scale_color_manual(values = "#36454f", name = "") +
  new_scale_color()+
  
  geom_line(aes(x = x, y = value, color = surplus), df_indp_tmp, size = 1) +
  scale_color_manual(values = my_blue, name = "Surplus, Ind. Claims") +
  new_scale_color() +
  
  
  # geom_point(mapping = aes(x = x, y = value, color = surplus), ruin_point, size = 2) +
  # new_scale_color() + 
  # geom_line(aes(x = x, y = y, color = id), out_one_indp$expected_proift, size = 1, alpha = 0.6)+
  # # geom_point(mapping = aes(x = x, y = y, color = id), data = profit_point, size = 2, shape = 15, alpha = 0.6) +
  # scale_color_manual(values = my_grey) +
  
  geom_line(aes(x = x, y = value, color = surplus), df_dep_tmp, size = 1) +
  scale_color_manual(values = my_red, name = "Surplus, Dep. Claims" )  +
  
  # geom_point(mapping = aes(x = x, y = value, color = surplus), ruin_point, size = 2) +
  # new_scale_color() + 
  # geom_line(aes(x = x, y = y, color = id), out_both_one_dep$expected_proift, size = 1, alpha = 0.6)+
  # # geom_point(mapping = aes(x = x, y = y, color = id), data = profit_point, size = 2, shape = 15, alpha = 0.6) +
  # #new_scale_color() +
  # scale_color_manual(values = my_green) +
  
  
  geom_hline(yintercept = 0, color = "black", size = 1) +
  geom_vline(xintercept = 0, color = "black", size = 1) +
  # geom_vline(xintercept = optimal_theta, color = "black", size = 0.5) +
  annotate("text", x = out_both_one_dep_1$optimal_theta-0.03 , y = 1.03,
           label = bquote(theta[profit]^"*" == .(out_both_one_dep_1$optimal_theta)), size = 5.5)  +
  annotate("text", x = out_both_one_dep_1$opt_ruin_theta + 0.03 , y = 1.03,
           label = bquote(theta[ruin_dep]^"*" == .(out_both_one_dep_1$opt_ruin_theta)), size = 5.5) +
  geom_segment(aes(x = out_both_one_dep_1$optimal_theta, xend = out_both_one_dep_1$optimal_theta, y = 0, yend = 1), linetype = 2) +
  geom_segment(aes(x = out_both_one_dep_1$opt_ruin_theta, xend = out_both_one_dep_1$opt_ruin_theta, y = 0, yend = 1), linetype = 2) +
  labs(x = expression(theta), y = "Probability of Ruin \n", color = "") +
  #ggtitle("") +
  theme_bw(base_size = 20) +
  scale_x_continuous(expand = c(0.0015, 0),limits = c(0.05, 0.45)) +
   scale_y_continuous(expand = c(0.002, 0), 
                      limits = c(0,1.07),
                      sec.axis = sec_axis(~.*out_both_one_dep_1$optimal_expected_profit, 
                                          name = "Expected Profit \n", breaks = c(0, 200, 400, 600, out_both_one_dep_1$optimal_expected_profit, 800), label = scales::comma)) +
  theme(legend.justification = "top") 
gg




# Both, two, indp, ----  
  

  
out_two_loading_indp_1 <-  two_loadings_indep(N = N, r = r, fixed_cost = fixed_cost, lambda = lambda, k = k, beta = beta, 
                                              claim_mean = k*beta, demand = demand2var,
                                              x_surplus = x_surplus[2], h_x = h_x, S_ = S_, a = k, b = beta)
  
  
  
save.image(file)
  
plot_ruin(out_two_loading_indp_1$V, out_two_loading_indp_1$theta_1, out_two_loading_indp_1$theta_2,
          type = "contour", t1_opt = out_two_loading_indp_1$t1_opt, t2_opt = out_two_loading_indp_1$t2_opt, 
          v_min = round(min(out_two_loading_indp_1$V[!is.na(out_two_loading_indp_1$V)]), digits = 2))

plot_profit(out_two_loading_indp_1$expected_income, 
            out_two_loading_indp_1$theta_1,
            out_two_loading_indp_1$theta_2,
            type = "contour", 
            t1_opt = out_two_loading_indp_1$t1_opt_prof, 
            t2_opt = out_two_loading_indp_1$t2_opt_prof,
            prof_max =  round(max(out_two_loading_indp_1$expected_income[!is.na(out_two_loading_indp_1$expected_income)]), digits = 0))


  
# Both two  dependent -----
  
# out_two_loading_dep_1 <-  two_loadings_dep(N = N, r = r, fixed_cost = fixed_cost, lambda = lambda, k = k, beta = beta, 
#                                              x_surplus = x_surplus[2], demand = demand2var, h_x = 100, theta_finess = 0.05,
#                                              f_z_max = 20000, nu = nu)
# 
# plot_ruin(out_two_loading_dep_1$V, out_two_loading_dep_1$theta_1, out_two_loading_dep_1$theta_2,
#           type = "contour", t1_opt = out_two_loading_dep_1$t1_opt, t2_opt = out_two_loading_dep_1$t2_opt, 
#           v_min = round(min(out_two_loading_dep_1$V[!is.na(out_two_loading_dep_1$V)]), digits = 2))
# 
# plot_profit(out_two_loading_dep_1$expected_income, 
#             out_two_loading_dep_1$theta_1_ei,
#             out_two_loading_dep_1$theta_2_ei,
#             type = "contour", 
#             t1_opt = out_two_loading_dep_1$t1_opt_prof, 
#             t2_opt = out_two_loading_dep_1$t2_opt_prof,
#             prof_max =  round(max(out_two_loading_dep_1$expected_income[!is.na(out_two_loading_dep_1$expected_income)]), digits = 0))
# 
# save.image(file)

#keep_logit <- out_1
#save(keep_logit, file = "keep_logit.RData")


# Depednent acquisition ----


# Dependent Case Clayton -----




# Gumbel copula -----

gumbel <- function(a,b,cop_par){
  return(exp(-((-log(a))^cop_par + (-log(b))^cop_par )^(1/cop_par)))
}


cop_params <- 1/(1-kendell)
clayton_gumbel <- list()

source("utils.R")
for(i in cop_params){
  print(paste0("copula parameter: ", i))
  clayton_gumbel[[paste0(i)]] <- list()
  for(x in c(x_surplus[2])){
    print(paste0("Surplus: ", x))
    clayton_gumbel[[paste0(i)]][[paste0(x)]] <- one_loading_inference_clayton2(N = N, 
                                                                              r = r, 
                                                                              fixed_cost = fixed_cost,
                                                                              lambda = lambda, 
                                                                              nu = nu, 
                                                                              k = k, 
                                                                              beta = beta,
                                                                              x_surplus = x,
                                                                              demand = demand, 
                                                                              h_x = h_x, 
                                                                              f_z_max = 30, 
                                                                              theta_grid = theta_grid_ruin,
                                                                              f_z_limit = 30,
                                                                              ord_copula = gumbel,
                                                                              cop_par = i)
    save.image(file)
    
    
  }
  
}


save.image(file)



ruin_values <- list()
for( i in 1:length(cop_params)){
  
  ruin_values[[i]] <-clayton_gumbel[[i]]$`5`$df
  ruin_values[[i]]$kendell <- as.character(kendell[i])
  ruin_values[[i]]$cop_param <- cop_params[i]
  
}

ruin_values <- as.data.frame(do.call(rbind,ruin_values))

ggplot() + geom_line(aes(x = x, y = value, colour =  kendell), data = ruin_values, size = 1 )+theme_bw(24)+
  xlab("\U03B8")+
  ylab("Probability of Ruin") +
  ggtitle("Gumbel Acquisition Structure")+
  labs(colour = "Kendell's tau")

#### Clayton ----


clayton <- function(a,b,cop_par){
  return((a^(-cop_par)+b^(-cop_par)-1)^(-1/cop_par))
}


cop_params <- 2*kendell/(1-kendell)
clayton_clayton <- list()


for(i in cop_params){
  print(paste0("copula parameter: ", i))
  clayton_clayton[[paste0(i)]] <- list()
  for(x in x_surplus[2]){
    print(paste0("Surplus: ", x))
    clayton_clayton[[paste0(i)]][[paste0(x)]] <- one_loading_inference_clayton2(N = N, 
                                                                               r = r, 
                                                                               fixed_cost = fixed_cost,
                                                                               lambda = lambda, 
                                                                               nu = nu, 
                                                                               k = k, 
                                                                               beta = beta,
                                                                               x_surplus = x,
                                                                               demand = demand, 
                                                                               h_x = h_x, 
                                                                               f_z_max = 30, 
                                                                               theta_grid = theta_grid,
                                                                               f_z_limit = 30,
                                                                               ord_copula = clayton,
                                                                               cop_par = i)
    save.image(file)
    
    
  }
  
}

ruin_values <- list()


for( i in 1:length(cop_params)){
  
  
  ruin_values[[i]] <-clayton_clayton[[i]]$`5`$df
  ruin_values[[i]]$kendell <- as.character(kendell[i])
  ruin_values[[i]]$cop_param <- cop_params[i]
  
}

ruin_values <- as.data.frame(do.call(rbind,ruin_values))

ggplot() + geom_line(aes(x = x, y = value, colour =  kendell), data = ruin_values, size = 1 )+theme_bw(24)+
  xlab("\U03B8")+
  ylab("Probability of Ruin") +
  ggtitle("Clayton Acquisition Structure")+
  labs(colour = "Kendell's tau")






# Save ----
  

save.image(file)

  
  