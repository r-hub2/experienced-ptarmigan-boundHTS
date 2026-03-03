## ----setup, message = FALSE, warning = FALSE----------------------------------
set.seed(123)
n_obs   <- 400
burn_in <- 150
time    <- seq_len(n_obs)
nodes <- c("A", "AA", "AB")

## ----bottom-level-ar1---------------------------------------------------------
# Bottom-level nodes: AA and AB
target_mean <- c(0.10, 0.05) # we want lower proportions for this simulation
mu <- qlogis(target_mean) # work backwards to estimate coefficients

beta1 <- c(0.05, 0.05)
beta0 <- mu * (1 - beta1) # ensures mean is stationary

# AR noise
sd_AA <- sd_AB <- 0.15
rho <- -0.6

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

# Add in zero inflation
zi_AA <- 0.05  # P(x=0) = 5%
zi_AB <- 0.05  # P(x=0) = 5%

inflate_with_zero <- function(mu, zi) {
  u <- runif(length(mu))
  y <- numeric(length(mu))

  # Structural zero
  y[u < zi] <- 0

  # Latent AR part
  idx <- u >= zi
  y[idx] <- mu[idx] 
  return(y)
}

AA_zero <- inflate_with_zero(
  mu  = AA,
  zi = zi_AA
)

AB_zero <-inflate_with_zero(
  mu  = AB,
  zi = zi_AB
)

plot(1:250, AA_zero, type='l')
lines(1:250, AA, col='red')
plot(1:250, AB_zero, type='l')


## ----dependence-check---------------------------------------------------------
# Dependence between bottom-level nodes
stats::cor.test(AA_zero, AB_zero, method = "kendall")

## ----aggregation-A------------------------------------------------------------
# Aggregation at node A with additional variability
A_mean <- 0.5 * AA_zero + 0.5 * AB_zero

kappa_A <- 200 # very tight to reduce noise
A <- stats::rbeta(
  n      = length(A_mean),
  shape1 = kappa_A * A_mean,
  shape2 = kappa_A * (1 - A_mean)
)

## Check white noise added to aggregation (+/- 5%)
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
sim_data <- dplyr::tibble(
  Time = 1:c(n_obs-burn_in),
  A  = A,
  B  = B,
  AA = AA_zero, 
  AB = AB_zero
)


## ----brms-settings, eval=FALSE------------------------------------------------
# training_data <- sim_data %>% filter(Time < 249)
# test_times <- sort(unique(sim_data$Time[sim_data$Time >= 249]))
# 
# nodes <- c("A", "AA", "AB")
# 
# form_base <- bf(value ~ s(Time, k = 4),
#                 phi ~ 1,
#                 zi ~ 1)
# 
# prior <- c(
#   brms::prior(normal(0, 1), class = "Intercept"),
#   brms::prior(normal(log(30), 0.5), class = "Intercept", dpar = "phi"),
#   brms::prior(normal(-4, 1), class = "Intercept", dpar = "zi"),
#   brms::prior(exponential(2), class = "sds")
# )
# 
# 
# fits <- list(
#   A  = brms::brm(brms::update(form_base, value ~ .),
#                  family = zero_inflated_beta(),
#                  data = training_data %>% rename(value = A),
#                  prior = prior,
#                  chains = 4, cores = 4, iter = 4000, warmup = 1000,
#                  backend = "cmdstanr",
#                  control = list(adapt_delta = 0.995, max_treedepth = 12)),
# 
#   AA = brms::brm(brms::update(form_base, value ~ .),
#                  family = zero_inflated_beta(),
#                  data = training_data %>% rename(value = AA),
#                  prior = prior,
#                  chains = 4, cores = 4, iter = 4000, warmup = 1000,
#                  backend = "cmdstanr",
#                  control = list(adapt_delta = 0.995, max_treedepth = 12)),
# 
#   AB = brms::brm(brms::update(form_base, value ~ .),
#                  family = zero_inflated_beta(),
#                  data = training_data %>% rename(value = AB),
#                  prior = prior,
#                  chains = 4, cores = 4, iter = 4000, warmup = 1000,
#                  backend = "cmdstanr",
#                  control = list(adapt_delta = 0.995, max_treedepth = 12)),
# )
# 

## ----eval=FALSE---------------------------------------------------------------
# cmdstanr::summary(fits$A)
# plot(fits$A)
# brms::pp_check(fits$A)
# 
# cmdstanr::summary(fits$AA)
# plot(fits$AA)
# brms::pp_check(fits$AA)
# 
# cmdstanr::summary(fits$AB)
# plot(fits$AB)
# brms::pp_check(fits$AB)
# 

## ----eval=FALSE---------------------------------------------------------------
# n_draws <- 1000
# n_nodes <- length(nodes)
# n_times <- length(test_times)
# 
# forecast_results <- vector("list", n_times)
# names(forecast_results) <- test_times
# 
# for (i in seq_along(test_times)) {
# 
#   t <- test_times[i]
#   message("Forecasting time ", t)
# 
#   data_up_to_t <- sim_data %>% dplyr::filter(Time <= t)
# 
#   years <- sort(unique(data_up_to_t$Time))
#   n_years <- length(years)
# 
#   pred_array <- array(NA, c(n_draws, n_nodes, n_years),
#                       dimnames = list(NULL, nodes, years))
#   phi_array  <- pred_array
#   zi_array  <- pred_array
# 
#   for (s in seq_along(nodes)) {
# 
#     node <- nodes[s]
#     fit  <- fits[[node]]
# 
#     node_data <- data_up_to_t %>%
#       dplyr::select(Time, value = all_of(node))
# 
#     pred_array[, s, ] <-
#       brms::posterior_predict(fit, newdata = node_data)[1:n_draws, ]
# 
#     phi_array[, s, ] <-
#       brms::posterior_epred(fit, dpar = "phi", newdata = node_data)[1:n_draws, ]
# 
#     zi_array[, s, ] <-
#       brms::posterior_epred(fit, dpar = "zi", newdata = node_data)[1:n_draws, ]
#   }
# 
#   forecast_results[[i]] <- list(
#     Time       = t,
#     post_draws = pred_array,
#     phi_draws  = phi_array,
#     zi_draws  = zi_array)
# }
# 
# 
# 

## ----eval=FALSE---------------------------------------------------------------
# mean_draws_1 <- apply(forecast_results[[1]]$post_draws, c(2,3), mean)
# mean_draws_2 <- apply(forecast_results[[2]]$post_draws, c(2,3), mean)
# 
# mu_mean_df <- dplyr::as_tibble(rbind(t(mean_draws_1), t(mean_draws_2)[250,])) %>%
#   dplyr::mutate(Time = 1:250)
# 

## ----eval=FALSE---------------------------------------------------------------
# # Combined posterior samples
# mu_post_df <- abind::abind(forecast_results[[1]]$post_draws, forecast_results[[2]]$post_draws[,,250])
# 

## ----eval=FALSE---------------------------------------------------------------
# phi_post_df <- abind::abind(forecast_results[[1]]$phi_draws, forecast_results[[2]]$phi_draws[,,250])
# zi_post_df <- abind::abind(forecast_results[[1]]$zi_draws, forecast_results[[2]]$zi_draws[,,250])
# 

## ----eval=FALSE---------------------------------------------------------------
# z_values <- seq(0, 1, length.out = 1000) # density grid
# top_node <- "A"
# bottom_nodes <- c("AA", "AB")
# phi_bottom_nodes <- c("phi_AA", "phi_AB")   # bottom-level precision parameters
# weights_bottom <- c(0.5, 0.5)
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
# 
# rec_df <- data.frame()
# 
# # Loop over prediction time points
# for (i in seq_along(test_times)) {
# 
#   t_i <- test_times[i]
# 
#   bottom_data <- sim_data %>% dplyr::filter(Time <= t_i) %>% dplyr::select('AA', 'AB') %>% as.matrix()
# 
#   # Sample bottom
#   weighted_samps <- rZIB_4p(n_mc = 100,
#                              sub_obs_data = bottom_data,
#                              phi_array = phi_post_df[,c(2,3),1:t_i],
#                              zi_array = zi_post_df[,c(2,3),1:t_i],
#                              weights = weights_bottom)
# 
# 
#   # Build marginal densities for unweighted bottom nodes
#   bottom_dens <- apply(mu_post_df[,c(2,3),], 2, density)
# 
#   for (b in 1:2) {
#     dens_i <- bottom_dens[[b]]
# 
#     tmp <- dplyr::tibble(
#       Node    = bottom_nodes[b],
#       Time    = t_i,
#       Z       = z_values,
#       Density = stats::approx(dens_i$x, dens_i$y,
#                        xout = z_values, rule = 2)$y
#     )
# 
#     rec_df <- dplyr::bind_rows(rec_df, tmp) # add density of bottom series to results dataframe
#   }
# 
# 
#   # Convolution to get A distribution
#   Density_top <- ZIB_convolution_density(
#     Y_mc = weighted_samps[,,t_i],
#     phi_array = phi_post_df[, c(2,3), t_i],
#     zi_array = zi_post_df[, c(2,3), t_i],
#     weights = weights_bottom,
#     z_values = z_values,
#     n_mc = 100
#   )
# 
#   tmp_top <- dplyr::tibble(
#     Node    = top_node,
#     Time    = t_i,
#     Z       = z_values,
#     Density = Density_top
#   )
# 
#   rec_df <- dplyr::bind_rows(rec_df, tmp_top) # add density of convoluted top series to results dataframe
# }
# 

## ----eval=FALSE---------------------------------------------------------------
# for (t in seq_along(test_times)) {
#   # Extract convolution densities for these nodes
#   rec_dens_df <- rec_df %>%
#     dplyr::filter(Time == test_times[t]) %>%
#     dplyr::pivot_wider(names_from = Node, values_from = Density) %>%
#     dplyr::select(Z, all_of(tilted_nodes))
# 
#   # Predicted means
#   mu_theory <- as.numeric(mu_mean_df[test_times[t], tilted_nodes])
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
#       nu_star <- stats::uniroot(moment_condition_tilting,
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
knitr::include_graphics("visualisations/example_ZIB_tilted_density_testset_A_ridges.png")


