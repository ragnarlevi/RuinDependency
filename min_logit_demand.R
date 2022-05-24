library(tidyverse)
library(ggnewscale) # add new color scale -> new_scale_color
# Set color brewer
library(RColorBrewer)
library(plotly)
library(progress)


source("utils.R")  # functions
source("parameters_diff_demand.R")  # Parameters used in numerical analysis


file = "./Data_clayton_1_diff_sens_may.RData"
print(nu)
h_x <- 0.3

# Explore Sensitivity ----

out_sens <- list()
for(x in c(x_surplus,20,30)){
  out_sens[[paste(x)]] <- sensitivity_test(fixed_cost = fixed_cost[1],
                               h_x = h_x,
                               lambda = lambda[1],
                               claim_mean = k[1]*beta[1],
                               x_surplus = x,
                               infinite_survival = infinite_survival3,
                               beta_grid = seq(from = 0.05, to = 0.7, by = 0.01),
                               S_ = S_marginal,
                               shape = k[1],
                               scale = beta[1])

}

out_sens <- bind_rows(out_sens)
out_sens$opt_ruin_val[out_sens$opt_ruin_val>1] <- 1
out_sens$x <- factor(out_sens$x)
out_sens$x <- as.character(out_sens$x)
out_sens$x[out_sens$x == "2"] = "initial surplus=2" 
out_sens$x[out_sens$x == "5"] = "initial surplus=5" 
out_sens$x[out_sens$x == "10"] = "initial surplus=10" 
out_sens$x[out_sens$x == "20"] = "initial surplus=20" 

out_sens$x <- factor(out_sens$x, levels = c("initial surplus=2","initial surplus=5","initial surplus=10","initial surplus=20"))
out_sens2 <-  out_sens %>%  group_by(x) %>% mutate(profit_scaled=opt_profit_val*1.3/max(opt_profit_val) )
out_sens2 <- out_sens2[out_sens2$x == "initial surplus=10",]
out_sens2$x <- "Maximal profit"


ggplot() + geom_line(aes(x = beta, y = opt_ruin_val, color = x), data = out_sens) +
  scale_color_discrete(name = "Minimal ruin probability:")+
  new_scale_color() + 
  geom_line(aes(x = beta, y = profit_scaled, color = x), data = out_sens2, linetype = "dashed") +
  scale_x_continuous(limits = c(2,10), name = bquote(beta[1]) ) +
  scale_color_discrete(name = "")+
  new_scale_color() +
  geom_line(aes(x = beta, y = opt_ruin, color = "Optimal ruin loading"), data = out_sens[out_sens$x =="initial surplus=5",], linetype = "longdash", size = 1) +
  scale_color_manual(name = "", values = c("Black"), labels = c(expression(paste(theta['ruin']))))+
  new_scale_color() +
  geom_line(aes(x = beta, y = opt_profit, color = "Optimal profit loading"), data = out_sens[out_sens$x =="initial surplus=5",], linetype = "dotdash", size = 1 ) +
  scale_color_manual(name = "", values = c( "orange"), labels = c(expression(paste(theta['profit']))))+
  scale_y_continuous(name = "Ruin probability / Loading", expand = c(0.002, 0), 
                     limits = c(0,1.6),
                     sec.axis = sec_axis(~.*max(out_sens$opt_profit_val)/1.3, 
                                         name = "Profit \n", breaks = c(0, 500, 1000, 1500,  2000, 2500, 3000), label = scales::comma)) +
  theme(legend.justification = "top") +
  theme_bw(20)
# 
# 
# 
# # Explore shape scale ----
# 
shape_test <- list()
for(x in c(x_surplus, 10, 20)){
  print(x)

  for(shape in 1:50){
    b1 <- log((1-p_demand)/p_demand) - 1/(1-p_demand)
    b2 <- 1/(tt*(p_demand))

    theta <- (1/b2)*(log(lambda[1]/(fixed_cost[1]*b2)) - b1)

    mean_true <-  1
    lambda_true <- lambda[1]*demand_1(theta)
    income <-  sum((1+theta)*demand_1(theta)*(lambda[1])) - fixed_cost[1]

    if(x >= 10){
      h_x1 <- 0.1
    }else{
      h_x1 <- 0.01
    }
    out_v <- infinite_survival3(to = x,
                               h_x = h_x1,
                               p = income,
                               lambda_true = lambda_true,
                               mean_true = mean_true,
                               S_ = S_marginal,
                               shape = shape,
                               scale = 1/shape)

    shape_test[[paste(x,shape)]] <- data.frame(x = x, shape = shape, scale = 1/shape, V = out_v$V[length(out_v$V)])

  }
}

shape_test <- bind_rows(shape_test)
shape_test$x <- factor(shape_test$x)

ggplot(shape_test) + geom_line(aes(x = scale, y = V, color = x)) +
  scale_color_discrete(name = "Initial surplus") +
  scale_x_continuous(name = "Scale k") +
  scale_y_continuous(name = "Probability of ruin") +
  theme_bw(20)

# 
# 

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
         profit_lines = data.frame(x = out_1$thetas, y = out_1$expected_proift/500, id = "Expected Profit"), 
         optimal_expected_profit = out_1$optimal_expected_profit, 
         optimal_theta = out_1$optimal_theta, 
         theta_lim = c(0, 0.5), 
         y_lim = c(0,1.07),
         opt_ruin_theta = round(out_1$opt_ruin_theta,3), 
         text_offset = 0.05,
         breaks = c(0, 100, 200, 300, 400, 500),
         size = 5.5)

my_blue <- colorRampPalette(brewer.pal(n = 9, "Blues")[3:7])(8)[c(1, 4, 8)]
v_min_1$id <- as.character(round( v_min_1$V_min , digits = 3))
v_min_1$id <- factor(v_min_1$id, levels = v_min_1$id)
plot_1 + new_scale_color() +
  geom_point(aes(x = theta_min, y = V_min, color = id, digits = 3), 
             v_min_1, size = 3) +
  scale_color_manual(values = my_blue, name = "Min. ruin prob.") +theme_bw(20)

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

out_2$df$value[out_2$df$value == 1] <- NA

plot_2 <- plot_all(ruin_lines = out_2$df, 
         profit_lines = data.frame(x = out_2$thetas, y = out_2$expected_proift/500, id = "Expected Profit"), 
         theta_lim = c(0, 0.5), 
         optimal_expected_profit = out_2$optimal_expected_profit, 
         optimal_theta = out_2$optimal_theta,
         y_lim = c(0,1.07),
         opt_ruin_theta = round(out_2$opt_ruin_theta,3), 
         text_offset = 0.04, 
         breaks = c(0, 100, 200, 300, 400, 500),
         size = 5.5)

my_blue <- colorRampPalette(brewer.pal(n = 9, "Blues")[3:7])(8)[c(1, 4, 8)]
v_min_2$id <- as.character(round( v_min_2$V_min , digits = 3))
v_min_2$id <- factor(v_min_2$id, levels = v_min_2$id)
plot_2 + new_scale_color() +
  geom_point(aes(x = theta_min, y = V_min, color = id, digits = 3), 
             v_min_2, size = 3) +
  scale_color_manual(values = my_blue, name = "Min. ruin prob.") +theme_bw(20)




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





# Joint with Levy copula dependence ----
# These take some time to calculate, The higher the surplus, the higher the complexity.

#   
# out_both_one_dep_1 <- one_loading_inference_dep2(N = N, 
#                                               r = r, 
#                                               fixed_cost = fixed_cost, 
#                                               lambda = lambda, 
#                                               nu = nu, 
#                                               k = k, 
#                                               beta = beta,
#                                               x_surplus = x_surplus[1],
#                                               demand = demand,
#                                               theta_grid_ruin = theta_grid_ruin,
#                                               infinite_survival = infinite_survival3,
#                                               h_x = h_x, 
#                                               f_z_max = 30, 
#                                               theta_finess = 0.01,
#                                               f_z_limit = 30)
# 
# save.image(file)
#  out_both_one_dep_2 <- one_loading_inference_dep2(N = N, 
#                                                 r = r, 
#                                                 fixed_cost = fixed_cost, 
#                                                 lambda = lambda, 
#                                                 nu = nu, 
#                                                 k = k, 
#                                                 beta = beta,
#                                                 x_surplus = x_surplus[2],
#                                                 theta_grid_ruin = theta_grid_ruin,
#                                                 infinite_survival = infinite_survival3,
#                                                 demand = demand, 
#                                                 h_x = h_x, 
#                                                 f_z_max = 30, 
#                                                 theta_finess = 0.01)
#  
# save.image(file)
# 
# out_both_one_dep_3 <- one_loading_inference_dep2(N = N, 
#                                                 r = r, 
#                                                 fixed_cost = fixed_cost, 
#                                                 lambda = lambda, 
#                                                 nu = nu, 
#                                                 k = k, 
#                                                 beta = beta,
#                                                 x_surplus = x_surplus[3],
#                                                 theta_grid_ruin = theta_grid_ruin,
#                                                 infinite_survival = infinite_survival3,
#                                                 demand = demand, 
#                                                 h_x = h_x, 
#                                                 f_z_max = 30, 
#                                                 theta_finess = 0.01)
# 
# save.image(file)
# 
# 
# 
# #save(out_both_one_dep_3, file = "dep_3.RData")
# 
# df_dep <- rbind(out_both_one_dep_1$df, out_both_one_dep_2$df, out_both_one_dep_3$df)
# profit_lines_dep <-  data.frame(x = out_both_one_dep_1$thetas, y = out_both_one_dep_1$expected_proift/out_one_indp$optimal_expected_profit, id = "Expected Profit Dependence")
# profit_lines_indp <- data.frame(x = out_one_indp$thetas, y = out_one_indp$expected_proift/out_one_indp$optimal_expected_profit, id = "Expected Profit Independence")
# v_min_dep <- df_dep %>% group_by(surplus) %>% 
#   summarise(V_min = min(value), theta_min = x[which(value == min(value))])
# 
# plot_dep <- plot_all(ruin_lines = df_dep, 
#                      profit_lines = profit_lines_dep,
#                      theta_lim = c(0,0.6),
#                      optimal_expected_profit = out_both_one_dep_1$optimal_expected_profit, 
#                       optimal_theta = out_both_one_dep_1$optimal_theta,
#                      opt_ruin_theta = out_both_one_dep_1$opt_ruin_theta, 
#                      text_offset = 0.07,
#                      breaks = c(0, 10000, 20000, 30000, 40000, 50000, out_both_one_dep_1$optimal_expected_profit),
#                     size = 5.5, y_lim = c(0, 1.08))
# 
# my_blue <- colorRampPalette(brewer.pal(n = 9, "Blues")[3:7])(8)[1:8]
# v_min_dep$id <- as.character(round(v_min_dep$V_min , digits = 3))
# v_min_dep$id <- factor(v_min_dep$id, levels = v_min_dep$id)
# plot_dep + new_scale_color() +
#   geom_point(aes(x = theta_min, y = V_min, color = id, digits = 3), 
#              v_min_dep, size = 3) +
#   scale_color_manual(values = my_blue, name = "Min. Ruin Prob.")
#   
# 
# 
# 
# df_dep_tmp <- df_dep
# df_dep_tmp$surplus <- as.character(df_dep_tmp$surplus)
# surplus_tmp <- unique(df_dep_tmp$surplus)
# for( j in 1:length(surplus_tmp)){
#   tmp_min <- round(v_min_dep$V_min[v_min_dep$surplus == surplus_tmp[j]], digits = 2)
#   tmp_min_theta <- round(v_min_dep$theta_min[v_min_dep$surplus == surplus_tmp[j]], digits = 3)
#   df_dep_tmp$surplus[df_dep_tmp$surplus == surplus_tmp[j]] <-  paste0(surplus_tmp[j], " -> ", tmp_min, ", \U03B8 = ", tmp_min_theta,  sep = "")
# }
# df_dep_tmp$surplus <- factor(df_dep_tmp$surplus, 
#                               levels =c("2 -> 0.71, \U03B8 = 0.275",
#                                         "5 -> 0.49, \U03B8 = 0.278", 
#                                         "10 -> 0.27, \U03B8 = 0.279") )
# 
# 
# df_indp_tmp <- out_one_indp$df
# df_indp_tmp$surplus <- as.character(df_indp_tmp$surplus)
# surplus_tmp <- unique(df_indp_tmp$surplus)
# for( j in 1:length(surplus_tmp)){
#   tmp_min <- round(v_min_indp$V_min[v_min_indp$surplus == surplus_tmp[j]], digits = 2)
#   tmp_min_theta <- round(v_min_indp$theta_min[v_min_indp$surplus == surplus_tmp[j]], digits = 3)
#   df_indp_tmp$surplus[df_indp_tmp$surplus == surplus_tmp[j]] <-  paste0(surplus_tmp[j], " -> ", tmp_min, ", \U03B8 = ", tmp_min_theta,  sep = "")
# }
# df_indp_tmp$surplus <- factor(df_indp_tmp$surplus, 
#                              levels =c("2 -> 0.69, \U03B8 = 0.268",
#                                        "5 -> 0.45, \U03B8 = 0.268", 
#                                        "10 -> 0.23, \U03B8 = 0.268") )
# 
# # my_blue = brewer.pal(n = 9, "Blues")[2:9] 
# my_blue <- colorRampPalette(brewer.pal(n = 9, "Blues")[3:7])(8)[c(1,4,8)]
# my_red <- colorRampPalette(brewer.pal(n = 9, "Reds")[4:8])(8)[c(1,4,8)]
# my_grey <- brewer.pal(n = 9, "Greys")[9] 
# 
# gg <- ggplot() +
#   
#   geom_line(aes(x = x, y = y), profit_lines_dep, size = 1, linetype = "dashed") +
#   scale_color_manual(values = "#882020", name = "") +
#   new_scale_color() + 
#   geom_line(aes(x = x, y = y), profit_lines_indp, size = 0.5, linetype = "dashed") +
#   # geom_point(mapping = aes(x = x, y = y, color = id), data = profit_point, size = 2, shape = 15, alpha = 0.6) +
#   #new_scale_color() +
#   scale_color_manual(values = "#36454f", name = "") +
#   new_scale_color()+
#   
#   geom_line(aes(x = x, y = value, color = surplus), df_indp_tmp, size = 1) +
#   scale_color_manual(values = my_blue, name = "Surplus, Ind. Claims") +
#   new_scale_color() +
#   
#   
#   # geom_point(mapping = aes(x = x, y = value, color = surplus), ruin_point, size = 2) +
#   # new_scale_color() + 
#   # geom_line(aes(x = x, y = y, color = id), out_one_indp$expected_proift, size = 1, alpha = 0.6)+
#   # # geom_point(mapping = aes(x = x, y = y, color = id), data = profit_point, size = 2, shape = 15, alpha = 0.6) +
#   # scale_color_manual(values = my_grey) +
#   
#   geom_line(aes(x = x, y = value, color = surplus), df_dep_tmp, size = 1) +
#   scale_color_manual(values = my_red, name = "Surplus, Dep. Claims" )  +
#   
#   # geom_point(mapping = aes(x = x, y = value, color = surplus), ruin_point, size = 2) +
#   # new_scale_color() + 
#   # geom_line(aes(x = x, y = y, color = id), out_both_one_dep$expected_proift, size = 1, alpha = 0.6)+
#   # # geom_point(mapping = aes(x = x, y = y, color = id), data = profit_point, size = 2, shape = 15, alpha = 0.6) +
#   # #new_scale_color() +
#   # scale_color_manual(values = my_green) +
#   
#   
#   geom_hline(yintercept = 0, color = "black", size = 1) +
#   geom_vline(xintercept = 0, color = "black", size = 1) +
#   # geom_vline(xintercept = optimal_theta, color = "black", size = 0.5) +
#   annotate("text", x = out_both_one_dep_1$optimal_theta-0.03 , y = 1.03,
#            label = bquote(theta[profit]^"*" == .(out_both_one_dep_1$optimal_theta)), size = 5.5)  +
#   annotate("text", x = out_both_one_dep_1$opt_ruin_theta + 0.03 , y = 1.03,
#            label = bquote(theta[ruin_dep]^"*" == .(out_both_one_dep_1$opt_ruin_theta)), size = 5.5) +
#   geom_segment(aes(x = out_both_one_dep_1$optimal_theta, xend = out_both_one_dep_1$optimal_theta, y = 0, yend = 1), linetype = 2) +
#   geom_segment(aes(x = out_both_one_dep_1$opt_ruin_theta, xend = out_both_one_dep_1$opt_ruin_theta, y = 0, yend = 1), linetype = 2) +
#   labs(x = expression(theta), y = "Probability of Ruin \n", color = "") +
#   #ggtitle("") +
#   theme_bw(base_size = 20) +
#   scale_x_continuous(expand = c(0.0015, 0),limits = c(0.05, 0.45)) +
#    scale_y_continuous(expand = c(0.002, 0), 
#                       limits = c(0,1.07),
#                       sec.axis = sec_axis(~.*out_both_one_dep_1$optimal_expected_profit, 
#                                           name = "Expected Profit \n", breaks = c(0, 200, 400, 600, out_both_one_dep_1$optimal_expected_profit, 800), label = scales::comma)) +
#   theme(legend.justification = "top") 
# gg




# Both, two, indp, ----  
#   
# 
# out_two_loading_indp_1 <-  two_loadings_indep(N = N, r = r, fixed_cost = fixed_cost, lambda = lambda, 
#                                               k = k, beta = beta, 
#                                               claim_mean = k*beta, demand = demand2var, 
#                                               theta_grid = theta_grid_two, theta_grid_profit = theta_grid_profit_two,
#                                               x_surplus = x_surplus[1], h_x = h_x, S_ = S_, a = k, b = beta)
#   
#   
#   
# save.image(file)
#   
# plot_ruin(out_two_loading_indp_1$V, out_two_loading_indp_1$theta_1, out_two_loading_indp_1$theta_2,
#           type = "contour", t1_opt = out_two_loading_indp_1$t1_opt, t2_opt = out_two_loading_indp_1$t2_opt, 
#           v_min = round(min(out_two_loading_indp_1$V[!is.na(out_two_loading_indp_1$V)]), digits = 2))
# 
# plot_profit(out_two_loading_indp_1$expected_income, 
#             out_two_loading_indp_1$theta_1,
#             out_two_loading_indp_1$theta_2,
#             type = "contour", 
#             t1_opt = out_two_loading_indp_1$t1_opt_prof, 
#             t2_opt = out_two_loading_indp_1$t2_opt_prof,
#             prof_max =  round(max(out_two_loading_indp_1$expected_income[!is.na(out_two_loading_indp_1$expected_income)]), digits = 0))
# 
# 
#   
# # Both two  dependent -----
#   
# out_two_loading_dep_1 <-  two_loadings_dep(N = N, r = r, fixed_cost = fixed_cost, lambda = lambda, k = k, beta = beta,
#                                              x_surplus = x_surplus[1], demand = demand2var, h_x = h_x, 
#                                            theta_grid = theta_grid_two, 
#                                            theta_grid_profit = theta_grid_profit_two,
#                                              f_z_max = 30, nu = nu)
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


# Optimize Find Best Independent ----


df_indp_two_opt <- list()
for(x in c(30)){
  print(x)
  
  if(x >= 30){
    h_x <- 1
  }else if(x >= 21){
    h_x <- 0.6
    
  }else if( x >= 9){
    
    h_x <- 0.3
  }else{
    h_x <- 0.1
  }
  out_independent <- optim(par = c(0.3, 0.28), fn = two_loadings_indep_opt_ruin,
                           fixed_cost = fixed_cost,
                           lambda = lambda, 
                           k = k, 
                           beta = beta,
                           x_surplus = x, 
                           claim_mean = beta*k,
                           demand = demand2var,
                           h_x = h_x,
                           S_ = S_)
  
  profit_ruin <- -two_loadings_indep_opt_profit(out_independent$par,
                                               fixed_cost = fixed_cost,
                                               lambda = lambda, 
                                               k = k, 
                                               beta = beta,
                                               x_surplus = x, 
                                               claim_mean = beta*k,
                                               demand = demand2var,
                                               h_x = h_x,
                                               S_ = S_)
  
  out_independent_profit <- optim(par = c(0.3, 0.28), fn = two_loadings_indep_opt_profit,
                                  fixed_cost = fixed_cost,
                                  lambda = lambda, 
                                  k = k, 
                                  beta = beta,
                                  x_surplus = x, 
                                  claim_mean = beta*k,
                                  demand = demand2var,
                                  h_x = h_x,
                                  S_ = S_)
  
  df_indp_two_opt[[paste(x)]] <- data.frame(theta1_ruin = out_independent$par[1],
                                            theta2_ruin = out_independent$par[2],
                                            ruin_val = out_independent$value,
                                            profit_ruin = profit_ruin,
                                            theta1_profit = out_independent_profit$par[1],
                                            theta2_profit = out_independent_profit$par[2],
                                            profit = -out_independent_profit$value,
                                            surplus = x)
  
  # save.image(file)
}




df_indp_two_opt <- bind_rows(df_indp_two_opt)
df_indp_two_opt$Case <- "Indep. claim"

save.image(file)




# Optimize Find Best Dependent ----
df_dep_two_opt <- list()
for(x in c(30)){
  print(x)
  
  if (x >= 30){
    h_x <- 1
    
  }else if (x >= 21){
    h_x <- 0.6
    
  }else if( x >= 9){
    
    h_x <- 0.3
  }else{
    h_x <- 0.1
  }
  
  init_par <- c(df_indp_two_opt$theta1_ruin[df_indp_two_opt$surplus == x], df_indp_two_opt$theta2_ruin[df_indp_two_opt$surplus == x])
  out_dependent <- optim(par = init_par, fn = two_loadings_dep_opt_ruin,
                         fixed_cost = fixed_cost,
                         lambda = lambda, 
                         k = k, 
                         beta = beta,
                         x_surplus = x, 
                         demand = demand2var,
                         h_x = h_x,
                         f_z_max = 30,
                         nu = nu,
                         control = list(trace = 1))
  
  
  
  profit_ruin <- -two_loadings_dep_opt_profit2(out_dependent$par,
                                              fixed_cost = fixed_cost,
                                              lambda = lambda, 
                                              k = k, 
                                              beta = beta,
                                              x_surplus = x, 
                                              demand = demand2var,
                                              h_x = h_x,
                                              f_z_max = 30,
                                              nu = nu)
  
  
  
  init_par <- c(df_indp_two_opt$theta1_profit[df_indp_two_opt$surplus == x], df_indp_two_opt$theta2_profit[df_indp_two_opt$surplus == x])
  out_dependent_profit <- optim(par = init_par, fn = two_loadings_dep_opt_profit,
                                fixed_cost = fixed_cost,
                                lambda = lambda, 
                                k = k, 
                                beta = beta,
                                x_surplus = x, 
                                demand = demand2var,
                                h_x = h_x,
                                f_z_max = 30,
                                nu = nu,
                                control = list(trace = 1))
  
  
  df_dep_two_opt[[paste(x)]] <- data.frame(theta1_ruin = out_dependent$par[1],
                                            theta2_ruin = out_dependent$par[2],
                                            ruin_val = out_dependent$value,
                                           profit_ruin = profit_ruin,
                                            theta1_profit = out_dependent_profit$par[1],
                                            theta2_profit = out_dependent_profit$par[2],
                                            profit = -out_dependent_profit$value,
                                            surplus = x)

  save.image(file)
}



df_dep_two_opt <- bind_rows(df_dep_two_opt)
df_dep_two_opt$Case <- "Dep. claim"
save.image(file)

# Optimize Find Best Dependent and Dependent Acquisition ----

cop_params <- 1/(1-kendell)

df_dep_acq_two_opt <- list()
for(x in c(30)){
  for(idx_cop in c(3)){
    
    if (x >= 30){
      h_x <- 1
      
    }else if (x >= 21){
      h_x <- 0.6
      
    }else if( x >= 9){
      
      h_x <- 0.3
    }else{
      h_x <- 0.1
    }
    
    print(paste(x, " ", idx_cop))
    init_par <- c(df_indp_two_opt$theta1_ruin[df_indp_two_opt$surplus == x], df_indp_two_opt$theta2_ruin[df_indp_two_opt$surplus == x])
    out_dependent_acq <- optim(par = init_par, fn = two_loadings_dep_acq_opt_ruin,
                           fixed_cost = fixed_cost,
                           lambda = lambda, 
                           k = k, 
                           beta = beta,
                           x_surplus = x, 
                           demand = demand2var,
                           h_x = h_x,
                           f_z_max = 30,
                           nu = nu,
                           ord_copula = gumbel, cop_par = cop_params[idx_cop],
                           control = list(trace = 1))
    
    profit_ruin <- -two_loadings_dep_acq_opt_profit(out_dependent_acq$par,
                                                    fixed_cost = fixed_cost,
                                                    lambda = lambda, 
                                                    k = k, 
                                                    beta = beta,
                                                    x_surplus = x, 
                                                    demand = demand2var,
                                                    h_x = h_x,
                                                    f_z_max = 30,
                                                    nu = nu,
                                                    ord_copula = gumbel, cop_par = cop_params[idx_cop])
    
    init_par <- c(df_indp_two_opt$theta1_profit[df_indp_two_opt$surplus == x], df_indp_two_opt$theta2_profit[df_indp_two_opt$surplus == x])
    out_dependent_acq_profit <- optim(par = init_par, fn = two_loadings_dep_acq_opt_profit,
                                  fixed_cost = fixed_cost,
                                  lambda = lambda, 
                                  k = k, 
                                  beta = beta,
                                  x_surplus = x, 
                                  demand = demand2var,
                                  h_x = h_x,
                                  f_z_max = 30,
                                  nu = nu,
                                  ord_copula = gumbel, cop_par = cop_params[idx_cop],
                                  control = list(trace = 1))
    
    
    df_dep_acq_two_opt[[paste(x,idx_cop)]] <- data.frame(theta1_ruin = out_dependent_acq$par[1],
                                                 theta2_ruin = out_dependent_acq$par[2],
                                                 ruin_val = out_dependent_acq$value,
                                                 profit_ruin = profit_ruin,
                                                 theta1_profit = out_dependent_acq_profit$par[1],
                                                 theta2_profit = out_dependent_acq_profit$par[2],
                                                 profit = -out_dependent_acq_profit$value,
                                                 surplus = x,
                                                 cop_param = cop_params[idx_cop],
                                                 kendall= kendell[idx_cop])
    
    save.image(file)
    
  }
  
  
}

df_dep_acq_two_opt <- bind_rows(df_dep_acq_two_opt)
df_dep_acq_two_opt$Case <- "Dep. claim, Dep. acquisition"
save.image(file)


all_opt <- bind_rows(df_indp_two_opt, 
                     df_dep_two_opt,
                    df_dep_acq_two_opt[df_dep_acq_two_opt$kendall == 0.15,]
                     )


# PLot opt ----

ggplot(df_indp_two_opt) + geom_line(aes(x = surplus, y = theta1_ruin, color = "1"), size = 1.02)+
  geom_line(aes(x = surplus, y = theta2_ruin, color = "2"),size = 1.02)+
  scale_color_manual(name = "", labels = c(expression(paste(theta['1,ruin'])), 
                                           expression(paste(theta['2,ruin']))),
                     values = c("lightblue", "Darkblue")) +
  new_scale_color()+
  geom_line(aes(x = surplus, y = theta1_profit, color = "Theta 1 w.r.t. profit"), linetype = "dotted", size = 1.01)+
  geom_line(aes(x = surplus, y = theta2_profit, color = "Theta 2 w.r.t. profit"), linetype = "dotted", size = 1.01)+
  # scale_color_manual(name = "", labels = c(expression(paste("   ",theta['1,profit'])==paste(theta['2,profit'])) ),
  #                    values = c("DarkRed")) +
  scale_color_manual(name = "", labels = c(expression(paste(theta['1,profit'])), expression(paste(theta['2,profit'])) ),
                                         values = c("#FF7F7F", "DarkRed")) +
  new_scale_color()+
  geom_line(aes(x = surplus, y = profit_ruin/500*0.325, color = "Profit w.r.t. ruin loadings"), linetype = "dashed")+
  scale_color_manual(name = "", labels = expression(paste("Expected profit at ", theta['ruin']^'*')),  
                                                    values = c("orange")) +
  scale_y_continuous(lim = c(0.15, 0.37), sec.axis = sec_axis(~.*500/0.325, 
                                         name = "Expected profit \n", breaks = c(0, 250, 300,350, 400, 450, 500, 550, 600), 
                                         label = scales::comma),
                     name =   bquote(theta))  +
  scale_x_continuous(name  = "Initial surplus")+

  theme_bw(20)


    ggplot(df_dep_two_opt) + geom_line(aes(x = surplus, y = theta1_ruin, color = "1"), size = 1.02)+
    geom_line(aes(x = surplus, y = theta2_ruin, color = "2"),size = 1.02)+
    scale_color_manual(name = "", labels = c(expression(paste(theta['1,ruin'])), 
                                             expression(paste(theta['2,ruin']))),
                       values = c("lightblue", "Darkblue")) +
    new_scale_color()+
    geom_line(aes(x = surplus, y = theta1_profit, color = "Theta 1 w.r.t. profit"), linetype = "dotted", size = 1.01)+
    geom_line(aes(x = surplus, y = theta2_profit, color = "Theta 2 w.r.t. profit"), linetype = "dotted", size = 1.01)+
    # scale_color_manual(name = "", labels = c(expression(paste("   ",theta['1,profit'])==paste(theta['2,profit'])) ),
    #                    values = c("DarkRed")) +
    scale_color_manual(name = "", labels = c(expression(paste(theta['1,profit'])), expression(paste(theta['2,profit'])) ),
                       values = c("#FF7F7F", "DarkRed")) +
    new_scale_color()+
    geom_line(aes(x = surplus, y = profit_ruin/500*0.325, color = "Profit w.r.t. ruin loadings"), linetype = "dashed")+
    scale_color_manual(name = "", labels = expression(paste("Expected profit at ", theta['ruin']^'*')),  
                       values = c("orange")) +
    scale_y_continuous(lim = c(0.15, 0.37), sec.axis = sec_axis(~.*500/0.325, 
                                                                name = "Expected profit \n", breaks = c(0, 250, 300,350, 400, 450, 500, 550, 600), 
                                                                label = scales::comma),
                       name =   bquote(theta))  +
    scale_x_continuous(name  = "Initial surplus")+
    
    theme_bw(20)


df_dep_acq_two_opt_plot <- df_dep_acq_two_opt[df_dep_acq_two_opt$kendall == 0.4,]
ggplot(df_dep_acq_two_opt_plot) + geom_line(aes(x = surplus, y = theta1_ruin, color = "1"), size = 1.02)+
  geom_line(aes(x = surplus, y = theta2_ruin, color = "2"),size = 1.02)+
  scale_color_manual(name = "", labels = c(expression(paste(theta['1,ruin'])), 
                                           expression(paste(theta['2,ruin']))),
                     values = c("lightblue", "Darkblue")) +
  new_scale_color()+
  geom_line(aes(x = surplus, y = theta1_profit, color = "Theta 1 w.r.t. profit"), linetype = "dotted", size = 1.01)+
  geom_line(aes(x = surplus, y = theta2_profit, color = "Theta 2 w.r.t. profit"), linetype = "dotted", size = 1.01)+
  # scale_color_manual(name = "", labels = c(expression(paste("   ",theta['1,profit'])==paste(theta['2,profit'])) ),
  #                    values = c("DarkRed")) +
  scale_color_manual(name = "", labels = c(expression(paste(theta['1,profit'])), expression(paste(theta['2,profit'])) ),
                     values = c("#FF7F7F", "DarkRed")) +
  new_scale_color()+
  geom_line(aes(x = surplus, y = profit_ruin/500*0.325, color = "Profit w.r.t. ruin loadings"), linetype = "dashed")+
  scale_color_manual(name = "", labels = expression(paste("Expected profit at ", theta['ruin']^'*')),  
                     values = c("orange")) +
  scale_y_continuous(lim = c(0.15, 0.37), sec.axis = sec_axis(~.*500/0.325, 
                                                              name = "Expected profit \n", breaks = c(0, 250, 300,350, 400, 450, 500, 550, 600), 
                                                              label = scales::comma),
                     name =   bquote(theta))  +
  scale_x_continuous(name  = "Initial surplus")+
  
  theme_bw(20)
  
all_opt$Case[all_opt$Case == "Indep. claim"] <- "Indep. risks"
all_opt$Case[all_opt$Case == "Dep. claim"] <- "Dep. risks"
all_opt$Case[all_opt$Case == "Dep. claim, Dep. acquisition"] <- "Dep. risks, dep. acquisitions"
all_opt$Case <- factor(all_opt$Case, levels = c("Indep. risks", "Dep. risks", "Dep. risks, dep. acquisitions"))

  ggplot(all_opt) + geom_line(aes(x = surplus, y = ruin_val , color = Case))+
    scale_color_discrete(name = "Cases") +
    scale_x_continuous(name = "Initial surplus") +
    scale_y_continuous(name = "Probability of ruin", lim = c(0,1))+
    theme_bw(20)

  
  
  ggplot(all_opt) + geom_point(aes(x = theta1_ruin, y = theta2_ruin , color = factor(surplus), shape = Case), size = 3)+
    geom_point(aes(x = x, y = y, shape = "Marginal"),
               data.frame(x = out_1$opt_ruin_theta,y = out_2$opt_ruin_theta), 
               color = 'black', 
               shape = "cross", size = 3)+
    scale_color_discrete(name = "Surplus")+
    scale_shape_discrete(name = "Cases")+
    scale_x_continuous(name = bquote(theta[1]), lim = c(0.27, 0.3) ,breaks = seq(0.27, 0.3, 0.01)) +
    scale_y_continuous(name =  bquote(theta[2]), lim = c(0.24, 0.27) ,breaks = seq(0.24, 0.27, 0.01)) +
    
    theme_bw(20)
  


  
    dep_acq_data_plot <- df_dep_acq_two_opt[df_dep_acq_two_opt$surplus == 9,]
  ggplot(dep_acq_data_plot) + geom_line(aes(x = kendall, y = theta1_ruin, color = "1"), size = 1.05) +
    geom_line(aes(x = kendall, y = theta2_ruin, color = "2"), size = 1.05)+
    scale_color_manual(name = "", labels = c(expression(paste(theta['1,ruin'])), 
                                             expression(paste(theta['2,ruin']))),
                       values = c("lightblue", "Darkblue")) +
    new_scale_color()+
    geom_line(aes(x = kendall, y = ruin_val/2.5, color = "12"),linetype = "dashed", size = 1.05) +
    scale_color_manual(name = "", labels = expression(paste("Optimal ruin probability ", theta['ruin'])),  
                       values = c("orange")) +
    scale_y_continuous(lim = c(0.12, 0.4), sec.axis = sec_axis(~.*(1*2.5),
                                           name = "Ruin probability \n",
                                           label = scales::comma),
                       name =   bquote(theta),
                       breaks = c(seq(0.1, 0.4, 0.05)))  +
    scale_x_continuous(name  = expression(paste("Kendall's ",tau)), breaks = c(seq(0, 0.8, 0.2)), lim = c(0, 0.8))+
    
    theme_bw(20)
    
  
  
# 
# 
# # Gumbel copula -----
# 
# gumbel <- function(a,b,cop_par){
#   return(exp(-((-log(a))^cop_par + (-log(b))^cop_par )^(1/cop_par)))
# }
# 
# 
# cop_params <- 1/(1-kendell)
# clayton_gumbel <- list()
# 
# source("utils.R")
# for(i in cop_params){
#   print(paste0("copula parameter: ", i))
#   clayton_gumbel[[paste0(i)]] <- list()
#   for(x in c(x_surplus[2])){
#     print(paste0("Surplus: ", x))
#     clayton_gumbel[[paste0(i)]][[paste0(x)]] <- one_loading_inference_clayton(N = N, 
#                                                                               r = r, 
#                                                                               fixed_cost = fixed_cost,
#                                                                               lambda = lambda, 
#                                                                               nu = nu, 
#                                                                               k = k, 
#                                                                               beta = beta,
#                                                                               x_surplus = x,
#                                                                               demand = demand, 
#                                                                               h_x = h_x, 
#                                                                               f_z_max = 30, 
#                                                                               theta_grid = theta_grid_ruin,
#                                                                               f_z_limit = 30,
#                                                                               ord_copula = gumbel,
#                                                                               cop_par = i)
#     save.image(file)
#     
#     
#   }
#   
# }
# 
# 
# save.image(file)
# 
# 
# 
# ruin_values <- list()
# for( i in 1:length(cop_params)){
#   
#   ruin_values[[i]] <-clayton_gumbel[[i]]$`5`$df
#   ruin_values[[i]]$kendell <- as.character(kendell[i])
#   ruin_values[[i]]$cop_param <- cop_params[i]
#   
# }
# 
# ruin_values <- as.data.frame(do.call(rbind,ruin_values))
# 
# ggplot() + geom_line(aes(x = x, y = value, colour =  kendell), data = ruin_values, size = 1 )+theme_bw(24)+
#   xlab("\U03B8")+
#   ylab("Probability of Ruin") +
#   ggtitle("Gumbel Acquisition Structure")+
#   labs(colour = "Kendell's tau")
# 
# #### Clayton ----
# 
# 
# clayton <- function(a,b,cop_par){
#   return((a^(-cop_par)+b^(-cop_par)-1)^(-1/cop_par))
# }
# 
# 
# cop_params <- 2*kendell/(1-kendell)
# clayton_clayton <- list()
# 
# 
# for(i in cop_params){
#   print(paste0("copula parameter: ", i))
#   clayton_clayton[[paste0(i)]] <- list()
#   for(x in x_surplus[2]){
#     print(paste0("Surplus: ", x))
#     clayton_clayton[[paste0(i)]][[paste0(x)]] <- one_loading_inference_clayton(N = N, 
#                                                                                r = r, 
#                                                                                fixed_cost = fixed_cost,
#                                                                                lambda = lambda, 
#                                                                                nu = nu, 
#                                                                                k = k, 
#                                                                                beta = beta,
#                                                                                x_surplus = x,
#                                                                                demand = demand, 
#                                                                                h_x = h_x, 
#                                                                                f_z_max = 30, 
#                                                                                theta_grid = theta_grid_ruin,
#                                                                                f_z_limit = 30,
#                                                                                ord_copula = clayton,
#                                                                                cop_par = i)
#     save.image(file)
#     
#     
#   }
#   
# }
# 
# ruin_values <- list()
# 
# 
# for( i in 1:length(cop_params)){
#   
#   
#   ruin_values[[i]] <-clayton_clayton[[i]]$`5`$df
#   ruin_values[[i]]$kendell <- as.character(kendell[i])
#   ruin_values[[i]]$cop_param <- cop_params[i]
#   
# }
# 
# ruin_values <- as.data.frame(do.call(rbind,ruin_values))
# 
# ggplot() + geom_line(aes(x = x, y = value, colour =  kendell), data = ruin_values, size = 1 )+theme_bw(24)+
#   xlab("\U03B8")+
#   ylab("Probability of Ruin") +
#   ggtitle("Clayton Acquisition Structure")+
#   labs(colour = "Kendell's tau")
# 
# 
# 
# 
# 
# 
# # Save ----
#   
# 
# save.image(file)
# 
# 
# ## Test optimality
# 
# alpha_fun2var(out$par)
# 
# alpha_fun2var <- function(param){
#   
#   denom <- sum(lambda*demand2var(param[1],param[2]))
#   num <- sum((1+param)*lambda*demand2var(param[1],param[2])*k*beta -fixed_cost)
#   
#   return(denom/num)
#   
#   
# }
# 
# out <- optim(c(0.3, 0.3), alpha_fun2var, method = "BFGS")
# out
# 
# 
# test_v <- infinite_survival3(to = 5,h_x = h_x, p = )
# 
# 
# #vector lambda 
# 
# 
# l <- N*demand2var(out$par[1], out$par[2])*lambda/sum(N*demand2var(out$par[1], out$par[2])*lambda)
# 
# 
# # The survival distribution is a sum
# S_tmp <- function(x, ...){
#   
#   return(sum(l*S_(x, ...)))
#   
# }
# 
# S_tmp <- Vectorize(S_tmp, vectorize.args = "x")
# 
# 
# 
# mean_true <-  sum(k*beta*l)
# lambda_true <-  sum(lambda*N*demand2var(out$par[1], out$par[2])) 
# p <-  sum((1+out$par)*N*demand2var(out$par[1], out$par[2])*(lambda*k*beta)) - sum(fixed_cost)
# 
# out_v <- infinite_survival3(to = 5, 
#                          h_x = 0.01, 
#                          p = p,
#                          lambda_true = lambda_true,
#                          mean_true = mean_true,
#                          S_ = S_tmp,
#                          a = k,
#                          b = beta)
# 
# out_v$V
# 
# 
# alpha_fun1 <- function(param){
#   
#   denom <- lambda[1]*demand_1(param)
#   num <- (1+param[1])*lambda[1]*demand_1(param)*k[1]*beta[1] -fixed_cost[1]
#   
#   return(denom/num)
#   
#   
# }
# 
# out1 <- optim(c(0.2), alpha_fun1)
# out1
# 
# 
# alpha_fun2 <- function(param){
#   
#   denom <- lambda[2]*demand_2(param)
#   num <- (1+param)*lambda[2]*demand_2(param)*k[2]*beta[2] -fixed_cost[2]
#   
#   return(denom/num)
#   
#   
# }
# 
# out2 <- optim(c(0.25), alpha_fun2)
# out2
# 
# 
# 


  
  