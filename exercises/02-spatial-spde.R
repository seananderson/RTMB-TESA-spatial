# Goals:
# - Work through a complete example of spatial modelling with the SPDE approach
#   to GMRFs inspecting all the steps along the way.

library(RTMB)
library(ggplot2)
options(ggplot2.continuous.colour = "viridis")
options(ggplot2.continuous.fill = "viridis")
theme_set(theme_light())

# Simulate some fake data ---------------------------------------------------

N_POINTS <- 1200
set.seed(123)
predictor_dat <- data.frame(X = runif(N_POINTS), Y = runif(N_POINTS))
prediction_grid <- expand.grid(
  X = seq(0, 1, length.out = 60),
  Y = seq(0, 1, length.out = 60)
)
d <- rbind(predictor_dat, prediction_grid)
.mesh <- sdmTMB::make_mesh(d, xy_cols = c("X", "Y"), cutoff = 0.05)
sim_dat <- sdmTMB::sdmTMB_simulate(
  formula = ~1,
  data = d,
  mesh = .mesh,
  range = 0.3,
  phi = 0.25,
  sigma_O = 0.4,
  seed = 123,
  B = 0.2 # intercept
)
sim_dat$true <- sim_dat$mu
sim_dat$mu <- NULL
sim <- sim_dat[1:N_POINTS, ]
grid_dat <- sim_dat[(N_POINTS + 1):nrow(sim_dat), ]

# truth
ggplot(grid_dat, aes(X, Y, fill = true)) +
  geom_raster()

# observations
ggplot(sim, aes(X, Y, colour = observed, size = observed)) +
  geom_point()

# Mesh and SPDE matrix construction -----------------------------------------

# make a 2D 'mesh' object:
mesh <- fmesher::fm_mesh_2d(
  loc = as.matrix(sim[, c("X", "Y")]),
  max.edge = c(0.1, 0.3), # maximum triangle edge length; inside + border
  cutoff = 0.06 # minimum triangle edge length
)

plot(mesh)
points(sim[, c("X", "Y")], col = "red", pch = 21, cex = 0.5)
mesh$n
mesh$loc

# compute 'finite element' matrices for SPDE approach:
spde <- fmesher::fm_fem(mesh)
# see equation between eq. (3) and eq. (4) in Lindgren and Rue (2015)
# https://doi.org/10.18637/jss.v063.i19
# Based on Lindgren et al. 2011
# https://doi.org/10.1111/j.1467-9868.2011.00777.x
spde$c0
spde$g1
spde$g2

# compute bilinear interpolation matrix from mesh to data:
interpolator_data <- fmesher::fm_basis(
  mesh,
  loc = as.matrix(sim[, c("X", "Y")])
)
nrow(sim)
nrow(interpolator_data)

mesh$n
ncol(interpolator_data)

interpolator_data[519, ]
# so, for row of data 519, bilinear interpolation is combination of these vertices:
vert <- which(interpolator_data[519, ] != 0)
vert
interpolator_data[519, vert]
plot(mesh)
points(mesh$loc[vert,], col = "blue", pch = 20)
points(sim[519, c("X", "Y")], col = "red", pch = 20, cex = 0.8)

# compute bilinear interpolation matrix from mesh to a prediction grid:
interpolator_prediction <- fmesher::fm_basis(
  mesh,
  loc = as.matrix(prediction_grid[, c("X", "Y")])
)

# RTMB setup  ---------------------------------------------------------------

nll <- function(par) {
  getAll(par, tmb_data)
  tau <- exp(log_tau)
  kappa <- exp(log_kappa)

  # compute precision matrix Q, e.g. Lindgren + Rue 2015 JSS p4:
  Q <- tau^2 * (kappa^4 * spde$c0 + 2 * kappa^2 * spde$g1 + spde$g2)
  rf %~% dgmrf(mu, Q)

  # project random effects from vertices to data locations:
  rf_at_observed <- interpolator_data %*% rf

  # data likelihood:
  observed %~% dnorm(as.vector(rf_at_observed), exp(log_obs_sd))

  # project random effects from vertices to prediction locations:
  predictions <- interpolator_prediction %*% rf

  # ADREPORT(predictions) # slower, but gets you SEs on the prediction grid
  REPORT(predictions) # fast!

  # derived values
  # distance correlation ~0.1 (p4 in Lindgren et al. 2011):
  range <- sqrt(8) / kappa
  ADREPORT(range)

  # marginal SD (p5 in Lindgren et al. 2011):
  sigma <- 1 / sqrt(4 * pi * exp(2 * log_tau + 2 * log_kappa))
  ADREPORT(sigma)
}

tmb_data <- list(observed = sim$observed)
tmb_data$mesh <- mesh
tmb_data$spde <- spde
tmb_data$interpolator_data <- interpolator_data
tmb_data$interpolator_prediction <- interpolator_prediction

par <- list(
  log_tau = 0,
  log_kappa = 0,
  rf = numeric(mesh$n),
  mu = 0,
  log_obs_sd = 0
)

obj <- MakeADFun(nll, par, random = "rf")

# Fit model -----------------------------------------------------------------

opt <- nlminb(obj$par, obj$fn, obj$gr)
# most time is spent on ADREPORT(predictions):
sdrep <- sdreport(obj)

Matrix::image(obj$env$spHess(random = TRUE))

# Look at fit ---------------------------------------------------------------

pl <- as.list(sdrep, "Estimate", report = TRUE)
plsd <- as.list(sdrep, "Std. Error", report = TRUE)

pl$range
plsd$range
pl$sigma
plsd$sigma

r <- obj$report(obj$env$last.par.best)

grid_dat$estimate <- r$predictions[,1]

g1 <- ggplot(grid_dat, aes(X, Y, fill = true)) +
  geom_raster() +
  ggtitle("True")
g2 <- ggplot(grid_dat, aes(X, Y, fill = estimate)) +
  geom_raster() +
  ggtitle("Estimated")
cowplot::plot_grid(g1, g2, nrow = 1)

ggplot(grid_dat, aes(true, estimate)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, colour = "red")

if (FALSE) { # only available if ADREPORT(predictions) uncommented
  grid_dat$se <- plsd$predictions
  ggplot(grid_dat, aes(X, Y, fill = se)) +
    geom_raster() +
    scale_fill_viridis_c(option = "C") +
    scale_size_area() +
    geom_point(
      data = sim, mapping = aes(X, Y, size = abs(observed)),
      inherit.aes = FALSE, pch = 21
    ) +
    ggtitle("Standard error")
}

est <- as.list(sdrep, "Estimate")
# Discussion:
# What is est$rf?
# What are the locations of those values?
est$rf
mesh$loc
mesh$n
dat <- data.frame(mesh$loc, rf = est$rf)
ggplot(dat, aes(X1, X2, colour = rf)) + geom_point()

# Fit same model with sdmTMB ------------------------------------------------

mesh0 <- sdmTMB::make_mesh(sim, c("X", "Y"), mesh = mesh)
fit <- sdmTMB::sdmTMB(observed ~ 1, data = sim, mesh = mesh0)
fit

p <- predict(fit, newdata = grid_dat)
ggplot(p, aes(X, Y, fill = est)) + geom_raster()

# same:
logLik(fit)
opt$objective

AIC(fit)
2 * length(opt$par) - 2 * -opt$objective

fit$sd_report
sdrep

sdmTMB::tidy(fit, "ran_pars", conf.int = TRUE)
pl$range
pl$sigma

# Exercises -----------------------------------------------------------------

# 1. Experiment with different numbers of data observations and
#    range sizes.
# 2. Experiment with different mesh resolutions.
# 3. Try switching `ADREPORT(predictions)` on to plot SEs in space.

