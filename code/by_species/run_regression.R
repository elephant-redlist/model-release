
library(plyr)
library(rstan)
library(ggplot2)
library(gridExtra)
library(reshape2)

###########
# PRELIMS #
###########

load("inputs.rda")

#reg_input <- reg_dat
#reg_data  <- with(reg_dat, data.frame(partition = reg_dat$partition[XP], site = sites[XS], year = years[XY], y = y))
#save(reg_data, reg_input, file = "inputs.rda")

dir_res <- file.path("..", "..", "results", "by_species")

dir.create(dir_res, recursive = TRUE)

stamps  <- c("constant", "linear")

for (stamp in stamps) {
  
  # get and compile model code
  reg_code <- readLines(file.path("../stan", paste0("regression_", stamp,".stan")))
  reg      <- stan_model(model_code = reg_code, model_name = stamp)
  
  #######
  # RUN #
  #######
  
  # get initial values
  load(paste0("inits_", stamp, ".rda"))
  
  # MCMC sampling
  reg_out <- sampling(reg, data = reg_input, init = function() reg_ini, iter = 2000, thin = 2, chains = 4, pars = c("eta", "scalar_log", "xt_log", "alpha_re"), include = FALSE) 
  
  save(reg_out, file = file.path(dir_res, paste0("fit_", stamp, ".rda")))
  
  # SUMMARY TRACE PLOTS
  x <- extract(reg_out, pars = "traces", permuted = FALSE)
  dimnames(x) <- list(iter = 1:dim(x)[1], chain = 1:4, parameter = c("initial density", "slope (constant)", "slope (linear)", "random site (constant)", "random site (linear)", "random error", "dispersion", "trend"))
  x <- reshape2::melt(x)
  
  gg <- ggplot(x) + 
    geom_line(aes(iter, value, col = as.factor(chain))) + 
    facet_wrap(~parameter, scales = "free_y") + theme_bw(base_size = 16) + 
    guides(col = "none") + theme(axis.text = element_blank()) +
    labs(x = "Iteration", y = "")
  ggsave(gg, file = file.path(dir_res, paste0("traces_", stamp, ".pdf")), width = 10)
  
  ############
  # FIGURE 1 #
  ############
  
  # fit to data
  # posterior prediction of data
  dfr <- extract(reg_out, "y_sim")[[1]]
  dfr <- apply(dfr, 2, function(x) c(value = mean(x), hat = median(x), low = as.numeric(quantile(x, 0.025)), upp = as.numeric(quantile(x, 0.975))))
  dfr <- cbind(reg_data, t(dfr))
  
  dfr$within_ci <- ifelse(dfr$upp > dfr$y & dfr$y > dfr$low, TRUE, FALSE) 
  
  gg <- ggplot(subset(dfr, y > 0)) + 
    geom_point(aes(x = y / 1e3, y = hat / 1e3, col = within_ci)) + 
    geom_errorbar(aes(x = y / 1e3, ymax = upp / 1e3, ymin = low / 1e3, col = within_ci), alpha = 0.2) +
    facet_grid(toupper(partition)~.) +
    labs(x = "Observed", y = "Predicted") + 
    geom_abline(intercept = 0, slope = 1, col = "black") + guides(col = "none") +
    theme_bw(base_size = 22) + theme(strip.background = element_blank(), strip.text = element_blank())
  gg <- gg + scale_y_log10(breaks = c(0.01, 1, 100), labels = c("0.01", "1", "100")) + scale_x_log10(breaks = c(0.01, 1, 100), labels = c("0.01", "1", "100"))

  gg_fit <- gg
  
  # trend
  dfr <- extract(reg_out, "trend_per_partition")[[1]]
  dimnames(dfr) <- list(iter = 1:dim(dfr)[1], partition = toupper(reg_input$partition))
  dfr <- reshape2::melt(dfr)
  
  rects     <- data.frame(xstart = c(Inf, 0.8, 0.5, 0.2), xend = c(0.8, 0.5, 0.2, -Inf), col = c("Least concern", "Vulnerable", "Endangered", "Critical"))
  rects$col <- factor(rects$col, levels = c("Least concern", "Vulnerable", "Endangered", "Critical"))
  
  dfr_med  <- ddply(dfr, .(partition), summarise, hat = median(value))
  dfr_mean <- ddply(dfr, .(partition), summarise, hat = mean(value))
  
  gg <- ggplot(dfr) + 
    geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = col), alpha = 0.3) +
    scale_fill_manual(values = c("green", "yellow", "orange", "red")) +
    geom_density(aes(x = value, y = after_stat(density)), fill = "grey") +
    geom_vline(data = dfr_med,  aes(xintercept = hat), linewidth = 1) +
    geom_vline(data = dfr_mean, aes(xintercept = hat), linewidth = 1, linetype = "dashed") +
    xlim(0.0, 1.0) + facet_grid(partition ~ ., scales = "free_y") + guides(fill = "none") + 
    labs(x = "Trend", y = "", fill = "") + theme_bw(base_size = 20) + theme(axis.title.y = element_blank(), axis.text.y = element_blank())
  
  gg_post <- gg
  
  gg <- grid.arrange(gg_fit, gg_post, ncol = 2, widths = c(7, 7)); dev.off()
  ggsave(gg, file = file.path(dir_res, paste0("trend_fit_", stamp, ".pdf")), width = 12, height = 5)
  
  ###############################
  # TABULATE DECLINE (TABLE 2)  #
  ###############################
  
  ff1 <- function(x) round(mean(x, na.rm = TRUE), 2)
  ff2 <- function(x) round(median(x, na.rm = TRUE), 2)
  ff3 <- function(x) paste0("(", round(quantile(x, 0.05, na.rm = TRUE), 2), ",", round(quantile(x, 0.95, na.rm = TRUE), 2), ")")
  
  tab_trend  <- data.frame(model = stamp, partition = reg_input$partition, trend_mean = NA, trend_median = NA, trend_ci = NA)
  tab_trend$partition <- factor(tab_trend$partition, levels = reg_input$partition)
  
  tmp <- extract(reg_out,  pars = "trend_per_partition")[[1]]
  
  tab_trend <- split(tab_trend, ~ partition)
  
  for (i in 1:length(tab_trend)) {
    
    # summary values
    tab_trend[[i]][1, "trend_mean"]   <- ff1(tmp[,i])
    tab_trend[[i]][1, "trend_median"] <- ff2(tmp[,i])
    tab_trend[[i]][1, "trend_ci"]     <- ff3(tmp[,i])
    
  }
  
  tab_trend <- ldply(tab_trend)[,-1]
  
  # print and save
  print(tab_trend)
    
  save(tab_trend, file = file.path(dir_res, paste0("tab_trend_", stamp, ".rda")))
  
  # tidy up
  rm(reg_ini, reg_out, gg, gg_fit, gg_post, tab_trend); gc()
}

