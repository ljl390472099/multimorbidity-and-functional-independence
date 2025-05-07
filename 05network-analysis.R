# Network analysis

# Load libraries
library(networktools)
library(bootnet)
library(mgm)
library(dplyr)
library(magrittr)
library(NetworkComparisonTest)
library(psych)
library(qgraph)
library(IsingSampler)
library(IsingFit)
library(tidyverse)

# Prepare datasets (residualizing covariates)
prepare_data <- function(data, disease_vars) {
  data |>
    select(all_of(c(disease_vars, "adl", "age", "sex", "marr", "edu", "exercise", "smoke", "drink", "income.c", "bmi.c"))) |>
    mutate(across(everything(), as.numeric)) |>
    na.omit()
}

df_net        <- prepare_data(elsa.long,   c("hb", "diabetes", "cancer", "lung", "heart_new", "stroke", "psy.dep", "arth", "mmse.c"))
df_net.charls <- prepare_data(charls.long, c("hb", "diabetes", "cancer", "lung", "heart", "stroke", "psy.dep", "arth", "memary.mmse"))
df_net.hrs    <- prepare_data(hrs.long,    c("hp", "diabetes", "cancer", "lung", "heart_new", "stroke", "psy.dep", "arth", "memary.mmse"))
df_net.mhas   <- prepare_data(mhas.long,   c("hb", "diabetes", "cancer", "lung", "hrtatt", "stroke", "dep", "arth", "mmse.c"))
df_net.share  <- prepare_data(share.long,  c("hb", "diabetes", "cancer", "lung", "heart", "stroke", "mental", "arth", "cognitive"))

# Generate network function
generate_network <- function(data, set.layout = FALSE) {
  set.seed(123)
  names(data)[1:10] <- paste0("D", 1:10)
  
  # Residualize each disease variable
  residual_mat <- matrix(nrow = nrow(data), ncol = 10) |> as.data.frame()
  colnames(residual_mat) <- paste0("D", 1:10)
  for (i in 1:10) {
    formula <- as.formula(paste0("D", i, " ~ age + sex + marr + edu + exercise + smoke + drink + income.c + bmi.c"))
    residual_mat[[i]] <- glm(formula, data = data, family = "binomial")$residuals
  }
  
  network_data <- as.matrix(residual_mat)
  p <- ncol(network_data)
  fit_obj <- mgm(data = network_data,
                 type = rep('g', p),
                 level = rep(1, p),
                 lambdaSel = 'CV',
                 ruleReg = 'OR',
                 pbar = TRUE)
  
  pred_obj <- predict(fit_obj, data = network_data, errorCon = 'R2')
  
  network_est <- estimateNetwork(data = network_data,
                                 default = "EBICglasso",
                                 tuning = 0.5,
                                 corMethod = "cor",
                                 corArgs = list(method = "spearman", use = "pairwise.complete.obs"))
  
  layout_to_use <- if (set.layout) p.network.elsa[[2]]$layout else "spring"
  
  p_net <- plot(network_est,
                layout = layout_to_use,
                groups = groups,
                label.cex = 1,
                label.color = 'black',
                negDashed = TRUE,
                legend = TRUE,
                nodeNames = items,
                legend.cex = 0.5,
                legend.mode = 'style1',
                pie = pred_obj$error[,2])
  
  return(list(network_est, p_net))
}

# Node labels and groupings
items <- list("Hypertension", "Diabetes", "Cancer", "Lung disease", "Heart disease",
              "Stroke", "Psychological disorder", "Arthritis", "Cognitive disorder", "Disability")
groups <- list("Chronic Diseases" = 1:9, "Outcome" = 10)

# Run network estimation
p.network.elsa   <- generate_network(df_net)
p.network.charls <- generate_network(df_net.charls, set.layout = TRUE)
p.network.hrs    <- generate_network(df_net.hrs, set.layout = TRUE)
p.network.mhas   <- generate_network(df_net.mhas, set.layout = TRUE)
p.network.share  <- generate_network(df_net.share, set.layout = TRUE)

# Save network plots
save_network_plot <- function(name, plot_object) {
  pdf(paste0("./results/network/", name, ".pdf"), width = 8, height = 5)
  plot(plot_object)
  dev.off()
}

save_network_plot("1charls", p.network.charls[[2]])
save_network_plot("2elsa",   p.network.elsa[[2]])
save_network_plot("3hrs",    p.network.hrs[[2]])
save_network_plot("4mhas",   p.network.mhas[[2]])
save_network_plot("5share",  p.network.share[[2]])

# Save centrality plots
save_centrality_plot <- function(name, network_object) {
  pdf(paste0("./results/network/", name, ".pdf"), width = 8, height = 5)
  centralityPlot(network_object, labels = items, include = c("Strength", "Closeness", "Betweenness"))
  dev.off()
}

save_centrality_plot("1_1charls", p.network.charls[[1]])
save_centrality_plot("2_1elsa",   p.network.elsa[[1]])
save_centrality_plot("3_1hrs",    p.network.hrs[[1]])
save_centrality_plot("4_1mhas",   p.network.mhas[[1]])
save_centrality_plot("5_1share",  p.network.share[[1]])

summary(p.network.share[[1]])

# Sensitivity analysis
network_sensitive <- function(network_data, nBoots = 100, nCores = 6) {
  set.seed(123)
  boot_result <- bootnet(network_data,
                         statistics = c("Strength", "Closeness", "Betweenness", "edge"),
                         nBoots = nBoots,
                         nCores = nCores)
  
  # Plot node strength and edge stability
  plot1 <- plot(boot_result, order = "sample") +
    theme_minimal(base_size = 16) +
    theme(axis.text = element_text(colour = "black"),
          legend.text = element_text(size = 14),
          strip.text = element_text(size = 16, face = "bold"))
  
  # Plot difference test for node strength
  plot2 <- plot(boot_result, "strength", plot = "difference") +
    theme_minimal(base_size = 16)
  
  # Edge stability plot
  plot3 <- plot(boot_result, "edge", plot = "difference") +
    theme_minimal(base_size = 16)
  
  # Correlation stability coefficient
  cs <- corStability(boot_result)
  cat("Correlation stability coefficient:", cs, "\n")
  
  return(list(boot = boot_result, plot1 = plot1, plot2 = plot2, plot3 = plot3, cs = cs))
}

# Example usage
 result_sens <- network_sensitive(p.network.elsa[[1]])
