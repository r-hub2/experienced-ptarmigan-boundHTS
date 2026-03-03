## ----setup, message = FALSE, warning = FALSE----------------------------------
set.seed(123)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(rjags)
library(coda)
library(splines)

## -----------------------------------------------------------------------------
N <- 2000 # Number of observations
M <- 4 # Number of bottom level series


## -----------------------------------------------------------------------------
b0 <- runif(M, min=0.4, max = 0.5) # M covariates
b1 <- runif(M, min=0.4, max = 0.5)
b2 <- runif(M, min=0.4, max = 0.5)
x <- runif(N, min = 0, max = 1) # N observations

## ----echo=FALSE---------------------------------------------------------------
cat("b0 coefficients: ", b0)
cat("\n b1 coefficients: ", b1)
cat("\n b2 coefficients: ", b2)
cat("\n x covariate: ", head(x, n=10))


## -----------------------------------------------------------------------------
lambda <- matrix(NA, nrow = N, ncol = M)
y <- matrix(NA, nrow = N, ncol = M)

for(m in 1:M) {
  lambda[, m] <- exp(b0[m] + b1[m] * x + b2[m] * x^2)
  y[, m] <- rpois(N, lambda = lambda[, m])
}

colnames(y) <- c("AA", "AB", "BA", "BB")
head(y)


## -----------------------------------------------------------------------------
eps <- matrix(NA, nrow = N, ncol = M)

eps[,1] <- sample(c(-1,0,1), N, replace = TRUE, prob = c(0.3, 0.4, 0.3))
eps[,2] <- sample(c(-1,0,1), N, replace = TRUE, prob = c(0.5, 0.4, 0.1))
eps[,3] <- sample(c(-1,0,1), N, replace = TRUE, prob = c(0.33, 0.33, 0.34))
eps[,4] <- sample(c(-1,0,1), N, replace = TRUE, prob = c(0.2, 0.3, 0.5))


## -----------------------------------------------------------------------------
y_star <- matrix(NA, nrow = N, ncol = M)

for(m in 1:M) {
  y_star[, m] <- pmax(0, y[, m] + eps[, m])  # ensure non-negative counts
}

colnames(y_star) <- c("AA", "AB", "BA", "BB")
head(y_star)


## -----------------------------------------------------------------------------
hist(y_star[,4], main = "Histogram of series BB", xlab = "Counts")


## ----top-level----------------------------------------------------------------

Tot = apply(y, 1, sum) # Top-level sum of the undisturbed series
Tot = ifelse(Tot < 0, 0, Tot) # Ensure positivity

# Final data set
sim_data <- tibble(
  X=x,
  Tot  = Tot,
  AA = y_star[,1], 
  AB = y_star[,2],
  BA = y_star[,3],
  BB = y_star[,4]
)


## -----------------------------------------------------------------------------
hts_data_long <- tidyr::pivot_longer(sim_data, cols = -c(X), names_to = "Level", values_to = "Value")

ggplot(hts_data_long, aes(x = X, y = Value, color = Level)) +
  geom_line() +
  labs(title = "Simulated hierarchical count", x = "X", y = "Count") +
  theme_minimal() + 
  theme(legend.position = "none") +
  facet_wrap(~Level, scales='free')


## -----------------------------------------------------------------------------
N <- nrow(sim_data) # 2000 rows
m <- 4 # bottom series
n_series <- ncol(sim_data[,-1])
n_train <- c(1:c(N-50)) # withhold the last 50 observations for validation
test_indices <- c(length(n_train)+1):N
n_samples <- length(test_indices)
sum_bottom <- c("AA", "AB", "BA", "BB")
top_y_vals <- seq(from = 0, to = c(max(sim_data$Tot)+10))
ally <- sim_data[,-c(1)] # remove covariates


## ----poisson-regression-------------------------------------------------------
# Containers for results
poiss_reg <- list()   # fitted GLMs at each step
poiss_GLM_reg <- list()
lambda_list <- list()   # estimated Poisson rates (lambda)
fitted_list <- list()   # predictive Poisson samples

# Bottom-level series
sum_bottom <- c("AA", "AB", "BA", "BB")

for (t in seq_along(test_indices)) {

  # Training indices (up to time t-1)
  train_idx <- seq_len(test_indices[t] - 1)

  # Store estimated lambdas and predictive samples
  lambda_est <- fits <- matrix(
    NA,
    nrow = length(train_idx) + 1,
    ncol = n_series
  )

  colnames(lambda_est) <- colnames(fits) <- colnames(ally)

  # Fit Poisson GLM separately for each bottom-level series
  for (series in colnames(ally)) {

    # Quadratic Poisson regression
    formula_str <- paste(series, "~ X + I(X^2)")
    fit <- glm(
      formula = as.formula(formula_str),
      data    = sim_data[train_idx, c("X", series)],
      family = "poisson"
    )

    # Store fitted model
    poiss_reg[[series]] <- fit

    # In-sample fitted Poisson means
    lambda_est[train_idx, series] <- fit$fitted.values

    # One-step-ahead Poisson mean
    lambda_est[test_indices[t], series] <-
      predict(fit, newdata = sim_data[test_indices[t], c("X", series)])

    # Draw predictive samples from Poisson distribution
    fits[, series] <- rpois(
      n = test_indices[t],
      lambda = lambda_est[, series]
    )
  }

  # Save results for this forecast origin
  lambda_list[[t]]     <- lambda_est
  fitted_list[[t]]     <- fits
  poiss_GLM_reg[[t]]   <- poiss_reg
}


## ----eval=FALSE---------------------------------------------------------------
# f_tilde_exp <- list()
# nu_exp <- list()
# f_y <- vector()
# bottom_sum <- c("AA", "AB", "BA", "BB")
# pmf_values <- vector()

## ----eval=FALSE---------------------------------------------------------------
# # Apply to predictions
# for(t in 1:length(test_indices)) {
#   # lambda values
#   lambda_vals <- as.data.frame(lambda_list[[t]])
#   mu_theory <- as.vector(unlist(lambda_vals[test_indices[t], ])) # predictive mean
# 
#   lambda_bseries <- lambda_vals[test_indices[t],bottom_sum]
#   lambda_conv <- sum(lambda_bseries) # sum poissons
#   lambda_vec <- c(lambda_conv, as.numeric(lambda_bseries)) # convolution and bottom series lambda
# 
#   # fitted values
#   fitted_vals <- as_tibble(fitted_list[[t]])
#   fitted_bseries <- fitted_vals[test_indices[t],bottom_sum]
# 
#   # Construct tilted pmf for top and bottom series
#   f_tilt <- matrix(NA, nrow = length(top_y_vals), ncol = length(lambda_vec))
#   for(k in 1:length(lambda_vec)) {
#     # Convolution step
#     f_y <- dpois(top_y_vals, lambda_vec[k]) # density of convolution
# 
#     # Solve for nu
#     nu_star <- uniroot(moment_condition_tilting, interval = c(-1, 1),
#                        f_y = f_y, y_vals = top_y_vals,
#                        mu_theory = mu_theory[k])$root
# 
#     f_tilt[,k] <- tilted_density_discrete(nu_star, f_y, top_y_vals)
#   }
#   colnames(f_tilt) <- colnames(fitted_vals)
# 
#   nu_exp[[t]] <- nu_star
#   f_tilde_exp[[t]] <- f_tilt # density of top level
# }
# 

## ----echo=FALSE, out.width = "100%",  fig.align = "center"--------------------
knitr::include_graphics("visualisations/example_pois_tilted_density_testset_Tot.png")


