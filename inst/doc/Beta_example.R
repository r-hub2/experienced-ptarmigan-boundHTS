## ----setup, message = FALSE, warning = FALSE----------------------------------
set.seed(123)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(rjags)
library(coda)
library(splines)


n_obs   <- 400
burn_in <- 150
time    <- seq_len(n_obs)

## ----bottom-level-ar1---------------------------------------------------------
# Bottom-level nodes: AA and AB
# Modelled as correlated AR(1) processes on the logit scale

beta0 <- c(0.2, 0.1)
beta1 <- c(0.65, 0.55)

# Stationary mean of each AR(1) process
mu <- beta0 / (1 - beta1)

# Innovation covariance
sd_AA <- sd_AB <- 0.15
rho   <- -0.6

Sigma <- matrix(
  c(sd_AA^2, rho * sd_AA * sd_AB,
    rho * sd_AA * sd_AB, sd_AB^2),
  nrow = 2
)

logit_mat <- matrix(NA_real_, n_obs, 2)
logit_mat[1, ] <- rnorm(2)

# Simulate correlated AR(1) dynamics on the logit scale
for (t in 2:n_obs) {
  logit_mat[t, ] <-
    mu +
    beta1 * (logit_mat[t - 1, ] - mu) +
    MASS::mvrnorm(1, mu = c(0, 0), Sigma = Sigma)
}

# Remove burn-in
logit_mat <- logit_mat[-seq_len(burn_in), ]

# Transform to proportions
AA <- plogis(logit_mat[, 1])
AB <- plogis(logit_mat[, 2])

## ----dependence-check---------------------------------------------------------
# Dependence between bottom-level nodes
stats::cor.test(logit_mat[, 1], logit_mat[, 2], method = "kendall")

## ----aggregation-A------------------------------------------------------------
# Aggregation at node A with additional variability
A_mean <- 0.5 * AA + 0.5 * AB

kappa_A <- 200
A <- stats::rbeta(
  n      = length(A_mean),
  shape1 = kappa_A * A_mean,
  shape2 = kappa_A * (1 - A_mean)
)

## ---- aggregation-diagnostics, fig.height=3 ----------------------------------
plot(
  A - A_mean,
  type = "l",
  ylab = expression(A[t] - E(A[t])),
  xlab = "Time"
)
abline(h = 0, lty = 2)

# Top-level split
B <- 1 - A


## ----top-level----------------------------------------------------------------

# Final data set
sim_data <- tibble(
  Time = 1:c(n_obs-burn_in),
  A  = A,
  B  = B,
  AA = AA, 
  AB = AB
)


## ----jags-settings------------------------------------------------------------
n_chains <- 3

# Initial long run (performed once)
init_burnin_iter <- 4000
init_sample_iter <- 4000
init_thin <- 2

# Subsequent rolling updates
update_short <- 1000
sample_short <- 2000
thin_short <- 1

# helper function
samples_to_matrix <- function(samples_mcmc_list, varname) {
  mat_list <- lapply(samples_mcmc_list, as.matrix)
  fullmat  <- do.call(rbind, mat_list)
  cols <- grep(paste0("^", varname, "(\\b|\\[)"),
               colnames(fullmat), value = TRUE)
  fullmat[, cols, drop = FALSE]
}


## ----eval=FALSE---------------------------------------------------------------
# preds_post_all <- preds_mean_all <- data.frame()
# n_testyears <- c(max(sim_data$Time)-50):max(sim_data$Time)
# 
# # -------------------------
# # Prepare initial training window (first t) to compile models once
# # -------------------------
# first_t <- min(n_testyears)
# 
# train_idx0 <- which(sim_data$Time < first_t)
# training_data0 <- sim_data[train_idx0, , drop = FALSE]
# yA0 <- training_data0[, "A", drop = FALSE]
# 
# # Top model initial data (for model compilation) ---
# N_pred0_top <- nrow(yA0) + 1
# X0 <- splines::bs(1:N_pred0_top, df = 10)
# K0 <- ncol(X0)
# 
# jags_data_top_env <- list2env(list(
#   n_bottom = 1,
#   X = X0,
#   K = K0,
#   n_obs = nrow(yA0),
#   N_pred = N_pred0_top,
#   y = as.matrix(yA0)
# ), parent = emptyenv())
# 
# 
# # Compile top model once
# cat("Compiling TOP model (this happens once)...\n")
# mod_top <- rjags::jags.model(
#   file = system.file("model", "beta_sim_top_series_model.txt", package = "boundHTS"),
#   data = jags_data_top_env,
#   n.chains = n_chains,
#   n.adapt = 500  # small adapt at compile; we'll burn-in explicitly below
# )
# 
# # Bottom model initial data
# # Use the first training dataset bottom columns
# yB0 <- training_data0[, c("AA", "AB"), drop = FALSE]
# N_pred0_bottom <- nrow(yB0) + 1
# 
# jags_data_bottom_env <- list2env(list(
#   n_bottom = ncol(yB0),
#   priormu_b0 = c(-0.85, -0.5), # default priors for compilation
#   priormu_b1 = c(0.5, 0.4),
#   priorsigma_b0 = c(0.5, 0.5),
#   priorsigma_b1 = c(0.3, 0.3),
#   prior_phi_alpha = 2,
#   prior_phi_beta = 2,
#   n_obs = nrow(yB0),
#   N_pred = N_pred0_bottom,
#   y = as.matrix(yB0)
# ), parent = emptyenv())
# 
# cat("Compiling BOTTOM model (this happens once)...\n")
# mod_bottom <- rjags::jags.model(
#   file = system.file("model", "beta_sim_bottom_series_model.txt", package = "boundHTS"),
#   data = jags_data_bottom_env,
#   n.chains = n_chains,
#   n.adapt = 500
# )

## ----eval=FALSE---------------------------------------------------------------
# # -------------------------
# # First (initial) warm long run for both models (long burnin + sampling)
# # -------------------------
# cat("Initial burn-in & sampling for TOP model...\n")
# 
# # burn-in
# update(mod_top, init_burnin_iter)
# 
# # sample
# params_top <- c("phi", "mu", "beta0", "beta_coef", "y_pred")
# samps_top_init <- coda.samples(mod_top, variable.names = params_top,
#                                n.iter = init_sample_iter, thin = init_thin)
# 
# # get mean of y_pred and phi and mu
# y_pred_top_mat <- samples_to_matrix(samps_top_init, "y_pred")
# phi_top_mat <- samples_to_matrix(samps_top_init, "phi")
# mu_top_mat <- samples_to_matrix(samps_top_init, "mu")
# 
# y_pred_top_mean <- colMeans(y_pred_top_mat)
# phi_top_mean <- colMeans(phi_top_mat)
# mu_top_mean <- colMeans(mu_top_mat)
# 
# cat("Initial burn-in & sampling for BOTTOM model...\n")
# params_bottom <- c("phi", "mu", "beta0", "beta_coef", "y_pred")
# update(mod_bottom, init_burnin_iter)
# samps_bottom_init <- coda.samples(mod_bottom, variable.names = params_bottom,
#                                   n.iter = init_sample_iter, thin = init_thin)
# 
# y_pred_bottom_mat <- samples_to_matrix(samps_bottom_init, "y_pred")
# phi_bottom_mat <- samples_to_matrix(samps_bottom_init, "phi")
# mu_bottom_mat <- samples_to_matrix(samps_bottom_init, "mu")
# 
# y_pred_bottom_mean <- colMeans(y_pred_bottom_mat)
# phi_bottom_mean <- colMeans(phi_bottom_mat)
# mu_bottom_mean <- colMeans(mu_bottom_mat)
# 
# # store initial mu estimates if first_t matches your save criteria
# if (first_t == 200) {
#   mu_A <- mu_top_mean # mean
#   mu_A_post <- mu_top_mat # posterior samples
#   mu_bottom <- mu_bottom_mean
#   mu_bottom_post <- mu_bottom_mat
# }
# 
# # Append results for first_t
# preds_mean_all <- bind_rows(preds_mean_all, as.data.frame(
#   cbind(Time = first_t,
#         A = y_pred_top_mean[ncol(y_pred_top_mat) ],      # last entry corresponds to N_pred
#         pred_AA = y_pred_bottom_mean[ ncol(y_pred_bottom_mat) ],
#         pred_AB = y_pred_bottom_mean[ ncol(y_pred_bottom_mat) - 1 ], # adjust if ordering differs
#         phi_A = mean(phi_top_mean),
#         phi_AA = phi_bottom_mean[1],
#         phi_AB = phi_bottom_mean[2])
# ))
# 
# # save posterior draws (for first_t)
# # For consistent column ordering create a data.frame with Time and draws
# n_draws_init <- nrow(y_pred_top_mat)
# Npred_index_top <- ncol(y_pred_top_mat)  # last column is the N_pred
# Npred_index_bottom <- ncol(y_pred_bottom_mat)
# 
# pred_post_init <- data.frame(
#   Time = rep(first_t, n_draws_init),
#   A = y_pred_top_mat[, Npred_index_top],
#   AA = y_pred_bottom_mat[, Npred_index_bottom - 1],  # adjust if necessary
#   AB = y_pred_bottom_mat[, Npred_index_bottom]
# )
# preds_post_all <- bind_rows(preds_post_all, pred_post_init)

## ----eval=FALSE---------------------------------------------------------------
# remaining_t <- setdiff(n_testyears, first_t)
# 
# for (t in remaining_t) {
#   cat("Processing t =", t, "...\n")
# 
#   train_idx <- which(sim_data$Time < t)
#   training_data <- sim_data[train_idx, , drop = FALSE]
# 
#   # TOP model: update environment variables for new training data -----
#   yA <- training_data[, "A", drop = FALSE]
#   N_pred_top <- nrow(yA) + 1
#   X_top <- splines::bs(1:N_pred_top, df = 10)
#   K_top <- ncol(X_top)
# 
#   # replace data in the top model's environment
#   assign("X", X_top, envir = jags_data_top_env)
#   assign("K", K_top, envir = jags_data_top_env)
#   assign("n_obs", nrow(yA), envir = jags_data_top_env)
#   assign("N_pred", N_pred_top, envir = jags_data_top_env)
#   assign("y", as.matrix(yA), envir = jags_data_top_env)
# 
#   # short warm-up and sample
#   update(mod_top, update_short)
#   samps_top <- coda.samples(mod_top, variable.names = params_top,
#                             n.iter = sample_short, thin = thin_short)
# 
#   # extract means and posterior predictive draws
#   y_pred_top_mat <- samples_to_matrix(samps_top, "y_pred")
#   phi_top_mat <- samples_to_matrix(samps_top, "phi")
#   mu_top_mat <- samples_to_matrix(samps_top, "mu")
# 
#   y_pred_top_mean <- colMeans(y_pred_top_mat)
#   phi_top_mean <- colMeans(phi_top_mat)
#   mu_top_mean <- colMeans(mu_top_mat)
# 
#   # BOTTOM model: update data env -----
#   yB <- training_data[, c("AA", "AB"), drop = FALSE]
#   N_pred_bottom <- nrow(yB) + 1
#   assign("n_bottom", ncol(yB), envir = jags_data_bottom_env)
#   assign("n_obs", nrow(yB), envir = jags_data_bottom_env)
#   assign("N_pred", N_pred_bottom, envir = jags_data_bottom_env)
#   assign("y", as.matrix(yB), envir = jags_data_bottom_env)
# 
#   update(mod_bottom, update_short)
#   samps_bottom <- coda.samples(mod_bottom, variable.names = params_bottom,
#                                n.iter = sample_short, thin = thin_short)
# 
#   y_pred_bottom_mat <- samples_to_matrix(samps_bottom, "y_pred")
#   phi_bottom_mat <- samples_to_matrix(samps_bottom, "phi")
#   mu_bottom_mat <- samples_to_matrix(samps_bottom, "mu")
# 
#   y_pred_bottom_mean <- colMeans(y_pred_bottom_mat)
#   phi_bottom_mean <- colMeans(phi_bottom_mat)
#   mu_bottom_mean <- colMeans(mu_bottom_mat)
# 
#   # create summary row
#   Npred_idx_top <- ncol(y_pred_top_mat)
#   Npred_idx_bottom <- ncol(y_pred_bottom_mat)
# 
#   preds_mean <- cbind(
#     Time = t,
#     A      = y_pred_top_mean[Npred_idx_top],
#     pred_AA = y_pred_bottom_mean[Npred_idx_bottom - 1],
#     pred_AB = y_pred_bottom_mean[Npred_idx_bottom],
#     phi_A  = mean(phi_top_mean),
#     phi_AA = phi_bottom_mean[1],
#     phi_AB = phi_bottom_mean[2]
#   )
# 
#   # posterior draws for this t
#   n_draws <- nrow(y_pred_top_mat)
#   pred_post <- data.frame(
#     Time = rep(t, n_draws),
#     A = y_pred_top_mat[, Npred_idx_top],
#     AA = y_pred_bottom_mat[, Npred_idx_bottom - 1],
#     AB = y_pred_bottom_mat[, Npred_idx_bottom]
#   )
# 
#   preds_mean_all <- bind_rows(preds_mean_all, as.data.frame(preds_mean))
#   preds_post_all <- bind_rows(preds_post_all, pred_post)
# }
# 

## ----eval=FALSE---------------------------------------------------------------
# 
# # Save estimates
# mu_estimates <- list(mean = list(mu_A, mu_bottom),
#                      post = list(mu_A_post, mu_bottom_post))
# 
# 

## ----eval=FALSE---------------------------------------------------------------
# 
# mu_mean <- mu_estimates[[1]] # means
# 
# # ---- Top-level (A) ------------------------------------------------
# mu_top <- tibble(
#   param = names(mu_mean[[1]]),
#   estimate = mu_mean[[1]]
# ) %>%
#   separate_wider_delim(
#     col = param,
#     names = c("Time", "index_top"),
#     delim = ","
#   ) %>%
#   mutate(
#     Time = gsub("mu\\[", "", Time),
#     index_top = gsub("\\]", "", index_top),
#     Node = ifelse(index_top == 1, "A", NA)
#   ) %>%
#   select(Time, Node, estimate) %>%
#   pivot_wider(names_from = Node, values_from = estimate)
# 
# # ---- Bottom-level (AA, AB) ----------------------------------------
# mu_bottom <- tibble(
#   param = names(mu_mean[[2]]),
#   estimate = mu_mean[[2]]
# ) %>%
#   separate_wider_delim(
#     col = param,
#     names = c("Time", "index_bottom"),
#     delim = ","
#   ) %>%
#   mutate(
#     Time = gsub("mu\\[", "", Time),
#     index_bottom = gsub("\\]", "", index_bottom),
#     Node = ifelse(index_bottom == 1, "AA", "AB")
#   ) %>%
#   select(Time, Node, estimate) %>%
#   pivot_wider(names_from = Node, values_from = estimate)
# 
# # Combined latent mean estimates
# mu_mean_df <- left_join(mu_top, mu_bottom)

## ----eval=FALSE---------------------------------------------------------------
# 
# mu_post <- mu_estimates[[2]]
# 
# # ---- Top-level posterior samples ----------------------------------
# mu_top_post <- as_tibble(mu_post[[1]]) %>%
#   mutate(sample_id = row_number()) %>%
#   pivot_longer(
#     cols = everything(),
#     names_to = "param",
#     values_to = "estimate"
#   ) %>%
#   separate_wider_delim(
#     col = param,
#     names = c("Time", "index_top"),
#     delim = ","
#   ) %>%
#   mutate(
#     Time = as.numeric(gsub("mu\\[", "", Time)),
#     index_top = gsub("\\]", "", index_top),
#     Node = ifelse(index_top == 1, "A", NA)
#   ) %>%
#   select(Time, sample_id, Node, estimate) %>%
#   pivot_wider(names_from = Node, values_from = estimate) %>%
#   filter(Time < 200)
# 
# # ---- Bottom-level posterior samples -------------------------------
# mu_bottom_post <- as_tibble(mu_post[[2]]) %>%
#   mutate(sample_id = row_number()) %>%
#   pivot_longer(
#     cols = everything(),
#     names_to = "param",
#     values_to = "estimate"
#   ) %>%
#   separate_wider_delim(
#     col = param,
#     names = c("Time", "index_bottom"),
#     delim = ","
#   ) %>%
#   mutate(
#     Time = as.numeric(gsub("mu\\[", "", Time)),
#     index_bottom = gsub("\\]", "", index_bottom),
#     Node = ifelse(index_bottom == 1, "AA", "AB")
#   ) %>%
#   select(Time, sample_id, Node, estimate) %>%
#   pivot_wider(names_from = Node, values_from = estimate) %>%
#   filter(Time < 200)
# 
# # Combined posterior samples
# mu_post_df <- left_join(mu_top_post, mu_bottom_post)
# 

## ----eval=FALSE---------------------------------------------------------------
# 
# pred_ysim <- preds_mean_all %>%
#   as_tibble() %>%
#   rename(
#     AA = pred_AA,
#     AB = pred_AB
#   ) %>%
#   select(Time, A, AA, AB, phi_A, phi_AA, phi_AB)
# 

## ----eval=FALSE---------------------------------------------------------------
# post_ysim <- preds_post_all %>%
#   select(Time, A, AA, AB) %>%
#   group_by(Time) %>%
#   mutate(sample_id = row_number())
# 

## ----eval=FALSE---------------------------------------------------------------
# z_values <- seq(0, 1, length.out = 1000) # density grid
# top_node <- "A"
# bottom_nodes <- c("AA", "AB")
# phi_bottom_nodes <- c("phi_AA", "phi_AB")   # bottom-level variance parameters
# weights_bottom <- c(0.5, 0.5)
# 
# 
# # Using the test set of years
# test_years <- sim_data %>%
#   select(Time) %>%
#   filter(Time >= 200) %>% # validation test set
#   unlist() %>%
#   as.vector()
# 
# p <- ncol(sim_data)-1 # N tilted_nodes in series
# groups <- list(2, c(2, 2))  # Hierarchical structure
# 
# # Create results dataframe to store results in
# rec_df <- tibble()
# 
# # Set up tilting inputs
# f_tilde_exp <- list()
# nu_exp <- list()
# exp_samps <- list()
# nu_star_mat <- matrix(NA, nrow=2, ncol=3) # n_bottom x n_series
# 
# # Define node set for A-group
# all_nodes <- c("A", "AA", "AB")
# tilted_nodes <- c("A")
# 

## ----eval=FALSE---------------------------------------------------------------
# # Loop over prediction time points
# for (i in seq_along(test_years)) {
# 
#   t_i <- test_years[i]
# 
#   # Extract posterior samples
#   P_top     <- post_ysim %>% filter(Time==t_i) %>% ungroup() %>% select(all_of(top_node))
#   P_bottom  <- post_ysim %>% filter(Time==t_i)  %>% ungroup() %>% select(all_of(bottom_nodes))
# 
#   # Extract mean phi for top + bottom
#   phi_top    <- pred_ysim[i, "phi_A"]
#   phi_bottom <- as.numeric(pred_ysim[i, phi_bottom_nodes])
# 
#   # Beta parameters for bottom distributions
#   alpha_bottom <- P_bottom * phi_bottom
#   beta_bottom  <- (1 - P_bottom) * phi_bottom
# 
#   # Sample bottom
#   weighted_samps <- array(NA, dim = c(nrow(alpha_bottom), 1, 2))
# 
#   for (j in 1:2) {
#     for(s in 1:nrow(alpha_bottom)) {
#       weighted_samps[s, 1, j] <- ExtDist::rBeta_ab(
#         1,
#         alpha_bottom[s,j],
#         beta_bottom[s,j],
#         0,
#         weights_bottom[j]
#       )
#     }
#   }
# 
#   # Build marginal densities for unweighted bottom nodes
#   bottom_dens <- apply(P_bottom, 2, density)
# 
#   for (b in 1:2) {
#     dens_i <- bottom_dens[[b]]
# 
#     tmp <- tibble(
#       Node    = bottom_nodes[b],
#       Time    = t_i,
#       Z       = z_values,
#       Density = approx(dens_i$x, dens_i$y,
#                        xout = z_values, rule = 2)$y
#     )
# 
#     rec_df <- bind_rows(rec_df, tmp)
#   }
# 
#   # Convolution to get A distribution
#   Density_top <- Beta_convolution_density_parallel(
#     z_values,
#     alpha_bottom,
#     beta_bottom,
#     weighted_samps,
#     weights_bottom
#   )
# 
#   tmp_top <- tibble(
#     Node    = top_node,
#     Time    = t_i,
#     Z       = z_values,
#     Density = Density_top
#   )
# 
#   rec_df <- bind_rows(rec_df, tmp_top)
# }
# 
# 

## ----eval=FALSE---------------------------------------------------------------
# for (t in seq_along(test_years)) {
#   # Extract convolution densities for these nodes
#   rec_dens_df <- rec_df %>%
#     filter(Time == test_years[t]) %>%
#     pivot_wider(names_from = Node, values_from = Density) %>%
#     select(Z, all_of(tilted_nodes))
# 
#   # Predicted means
#   mu_theory <- as.numeric(pred_ysim[t, tilted_nodes])
# 
#   # y grid
#   y_vals <- sort(unique(rec_dens_df$Z))
# 
#   # Initialize results
#   nu_star_vec <- numeric(length(tilted_nodes))
#   tilted_samps <- matrix(NA, nrow = 5000, ncol = length(tilted_nodes))
#   f_tilted <- matrix(NA, nrow = length(z_values), ncol = length(tilted_nodes))
#   colnames(tilted_samps) <- colnames(f_tilted) <- tilted_nodes
# 
#   for (i in seq_along(tilted_nodes)) {
#     name <- tilted_nodes[i]
# 
#     # Extract and normalise convolution base density
#     f_y <- as.numeric(rec_dens_df[[name]])
#     f_y <- f_y / pracma::trapz(y_vals, f_y)
# 
#     # Base mean
#     mu_base <- pracma::trapz(y_vals, y_vals * f_y)
# 
#     cat("Node:", name, "\n")
#     cat("mu base:", mu_base, "\n")
#     cat("mu theory:", mu_theory[i], "\n")
# 
#     # Root finding bracket search
#     nu_grid <- seq(-2000, 2000, length.out = 4000)
#     vals <- sapply(nu_grid, moment_condition_tilting,
#                    f_y = f_y,
#                    y_vals = y_vals,
#                    mu_theory = mu_theory[i])
# 
#     idx <- which(diff(sign(vals)) != 0)
# 
#     if (length(idx) == 0) {
#       warning("No sign change for ", name,
#               " — using no tilting (nu=1).")
#       nu_star <- 0
#     } else {
#       lower <- nu_grid[idx[1]]
#       upper <- nu_grid[idx[1] + 1]
#       nu_star <- uniroot(moment_condition_tilting,
#                          lower = lower,
#                          upper = upper,
#                          f_y = f_y,
#                          y_vals = y_vals,
#                          mu_theory = mu_theory[i])$root
#     }
# 
#     nu_star_vec[i] <- nu_star
# 
#     # Tilted density + samples
#     f_tilt_i <- tilted_density_cont(nu_star, f_y, y_vals)
#     f_tilt_i <- f_tilt_i / pracma::trapz(y_vals, f_tilt_i)
#     f_tilted[, i] <- approx(y_vals, f_tilt_i, xout = z_values, rule = 2)$y
# 
#     tilted_samps[, i] <- sample(y_vals, 5000, replace = TRUE, prob = f_tilt_i)
# 
#     cat("Tilted mean:", pracma::trapz(y_vals, y_vals * f_tilt_i), "\n")
#   }
# 
#   # Save per-t results
#   nu_exp[[t]]      <- nu_star_vec
#   exp_samps[[t]]   <- tilted_samps
#   f_tilde_exp[[t]] <- f_tilted
# }
# 
# 

## ----echo=FALSE, out.width = "100%",  fig.align = "center"--------------------
knitr::include_graphics("visualisations/example_beta_tilted_density_testset_A.png")

## ----echo=FALSE, out.width = "100%",  fig.align = "center"--------------------
knitr::include_graphics("visualisations/example_beta_tilted_density_testset_A_ridges.png")


