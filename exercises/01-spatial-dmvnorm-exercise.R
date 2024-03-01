# Goals:
# - Code an exponential correlation function
# - Learn to fit a Gaussian random field directly with RTMB
# - Observe computational performance of this approach

library(RTMB)
library(ggplot2)
options(ggplot2.continuous.colour = "viridis")
options(ggplot2.continuous.fill = "viridis")
theme_set(theme_light())

# Simulate some data --------------------------------------------------------

N_POINTS <- 80

set.seed(123)
loc <- data.frame(X = runif(N_POINTS), Y = runif(N_POINTS))
prediction_grid <- expand.grid(
  X = seq(0, 1, length.out = 10),
  Y = seq(0, 1, length.out = 10)
)
d <- rbind(loc, prediction_grid)
.mesh <- sdmTMB::make_mesh(d, xy_cols = c("X", "Y"), cutoff = 0.05)
dat <- sdmTMB::sdmTMB_simulate(
  formula = ~1,
  data = d,
  mesh = .mesh,
  range = 0.8,
  phi = 0.25,
  sigma_O = 0.4,
  seed = 1,
  B = 0.2 # intercept
)
dat$true <- dat$mu
dat$mu <- NULL
dat$observed[(N_POINTS + 1):nrow(dat)] <- NA
observed_df <- dat[1:N_POINTS, ]
prediction_grid <- dat[(N_POINTS + 1):nrow(dat), ]

# Look at the dat ----------------------------------------------------------

# our observations:
ggplot(observed_df, aes(X, Y, size = observed, colour = observed)) +
  geom_point()

# a grid we will predict on for visualization:
ggplot(prediction_grid, aes(X, Y)) +
  geom_point(pch = 21)

# true process:
ggplot(prediction_grid, aes(X, Y, fill = true)) +
  geom_raster()

nrow(dat)

# RTMB model ----------------------------------------------------------------

nll <- function(pars) {
  getAll(pars, dat)
  sigma_rf <- exp(log_sigma_rf)
  ADREPORT(sigma_rf)
  rho <- exp(log_rho)
  ADREPORT(rho)

  # calculate correlation_matrix for exponential correlation decay:
  N <- length(rf)
  correlation_matrix <- matrix(0, nrow = N, ncol = N)
  for (i in 1:N) {
    for (j in 1:N) {
      s1 <- dat[i, c("X", "Y")] # grab spatial coordinates
      s2 <- dat[j, c("X", "Y")] # grab spatial coordinates
      d <- sqrt(sum((s1 - s2)^2)) # Euclidean distance
      correlation_matrix[i, j] <-  ### Exercise: fill in the exponential correlation function (use d and rho)
    }
  }
  # convert from correlation to covariance matrix:
  covariance_matrix <- sigma_rf^2 * correlation_matrix

  # random field likelihood:
  rf %~% dmvnorm(mu, covariance_matrix)

  # data likelihood:
  not_NA <- !is.na(observed) # the NA observed are the grid
  observed[not_NA] %~% dnorm(rf[not_NA], exp(log_sigma))
}

pars <- list(
  rf = numeric(nrow(dat)), # our spatial random field values
  log_sigma_rf = 0, # controls random field wiggle magnitude
  log_rho = -1, # controls correlation with distance
  mu = 0, # mean
  log_sigma = 0 # log observation SD
)

obj <- MakeADFun(nll, pars, random = "rf")
opt <- nlminb(obj$par, obj$fn, obj$gr)

sdr <- sdreport(obj, opt$par)
sdr_est <- as.list(sdr, "Estimate")
sdr_se <- as.list(sdr, "Std. Error")

# Inspecting the output -----------------------------------------------------

# Question: how sparse or dense is the Hessian?
Matrix::image(obj$env$spHess(random = TRUE))

# Our fixed effect estimates:
sdr

# grab estimates and standard errors at prediction locations:
prediction_grid$est <- sdr_est$rf[is.na(dat$observed)]
prediction_grid$se <- sdr_se$rf[is.na(dat$observed)]

# underlying truth:
g1 <- prediction_grid |>
  ggplot(aes(X, Y, fill = true)) +
  geom_raster() +
  scale_fill_viridis_c(limits = range(prediction_grid$true))

# prediction with matching colour scale:
g2 <- prediction_grid |>
  ggplot(aes(X, Y, fill = est)) +
  geom_raster() +
  geom_point(aes(X, Y), data = observed_df, inherit.aes = FALSE) +
  scale_fill_viridis_c(limits = range(prediction_grid$true))

cowplot::plot_grid(g1, g2, nrow = 1)

# standard error:
prediction_grid |>
  ggplot(aes(X, Y, fill = se)) +
  geom_raster() +
  geom_point(aes(X, Y), data = observed_df, inherit.aes = FALSE) +
  scale_fill_viridis_c(option = "A")

prediction_grid |>
  ggplot(aes(true, est)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

# Exercises -----------------------------------------------------------------

# 1. Fill in the exponential correlation and fit the model.
# 2. Try simulating and fitting different sample sizes. How well does this scale?
# 3. How sparse or dense is the Hessian?
# 4. Discuss why the prediction does or does not match the truth well.
# 5. Discuss why the standard error is larger in some locations than others.
