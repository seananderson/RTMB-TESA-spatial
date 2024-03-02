# Goals:
# - Practice filling some key parts of an SPDE-based spatial GMRF model
# - Gain experience with a custom non-linear spatial model that's not easily
#   fit with any package such as sdmTMB.
# - Gain additional familiarity interpreting output from these models
# - Investigate the consequences of allowing early life growth rate to vary
#   spatially.

library(RTMB)
library(ggplot2)
library(fmesher)
options(ggplot2.continuous.colour = "viridis")
theme_set(theme_light())

source(here::here("exercises/03-growth.R"))
pop <- readRDS(here::here("data/pop.rds"))

# Plot raw data -------------------------------------------------------------

ggplot(pop, aes(X, Y, colour = length)) + geom_point()
ggplot(pop, aes(X, Y, colour = age)) + geom_point()

grid_dat <- sdmTMB::qcs_grid
ggplot(grid_dat, aes(X, Y, fill = depth)) + geom_raster() +
  scale_fill_viridis_c(option = "G", direction = -1, trans = "log10") +
  coord_fixed()

# Mesh and SPDE setup -------------------------------------------------------

# make a 2D 'mesh' object:
mesh <- fmesher::fm_mesh_2d(
  loc = pop[, c("X", "Y")], # UTM coordinates in km
  max.edge = c(30, 50),     # maximum triangle edge length; inside + border
  offset = c(15, 40),       # boundary around data and in border
  cutoff = 3                # minimum triangle edge length
)
plot(mesh)
points(pop[, c("X", "Y")])

# compute 'finite element' matrices for SPDE approach:
spde <- fmesher::fm_fem(mesh)

# compute bilinear interpolation matrix from mesh to data:
interpolator_data <- fmesher::fm_basis(mesh, loc = as.matrix(pop[, c("X", "Y")]))

# compute bilinear interpolation matrix from mesh to a prediction grid:
interpolator_prediction <- fmesher::fm_basis(mesh, as.matrix(grid_dat))

# RTMB model ----------------------------------------------------------------

tmb_data <- list(
  age = pop$age,
  length = pop$length,
  X = pop$X,
  Y = pop$Y,
  spde = spde,
  interpolator_data = interpolator_data,
  interpolator_prediction = interpolator_prediction
)

par <- list(
  log_L_inf = log(max(pop$length)),
  log_tau = 0,
  log_kappa = 0,
  log_omega = numeric(mesh$n),
  mu_log_omega = 0,
  t0 = 0,
  log_obs_sd = 0
)

nll <- function(par) {
  getAll(par, tmb_data)

  L_inf <- exp(log_L_inf)
  tau <- exp(log_tau)
  kappa <- exp(log_kappa)

  # GMRF for k:
  # SPDE-based precision matrix, Q:
  Q <- tau^2 * (kappa^4 * spde$c0 + 2 * kappa^2 * spde$g1 + spde$g2)

  # omega: growth rate of length measurement units per year in early life
  # omega = k * L_inf

  # random effects likelihood:
  log_omega %~% dgmrf(mu_log_omega, Q)

  # project from knots to data locations and exp():
  omega <- exp(interpolator_data %*% log_omega)

  # growth curve
  # Gallucci and Quinn (1979)
  # Cahill et al. 2020 https://doi.org/10.1139/cjfas-2019-0434
  pred <- L_inf * (1 - exp(-(omega / L_inf) * (age - t0)))

  # project omega to prediction grid for visualization:
  omega_pred <- exp(interpolator_prediction %*% log_omega)
  ADREPORT(omega_pred)
  # REPORT(omega_pred) # faster, but no SEs on predictions

  # data likelihood:
  log(length) %~% dnorm(as.vector(log(pred)), exp(log_obs_sd))
  # derived values:

  # distance correlation ~0.1 (p4 in Lindgren et al. 2011):
  range <- sqrt(8) / kappa
  ADREPORT(range)
  # marginal SD (p5 in Lindgren et al. 2011):
  sigma <- 1 / sqrt(4 * pi * exp(2 * log_tau + 2 * log_kappa))
  ADREPORT(sigma)
}

obj <- MakeADFun(nll, par, random = "log_omega")

spatial_fit <- nlminb(obj$par, obj$fn, obj$gr)

sdrep <- sdreport(obj)
sdrep

# Inspect parameter estimates -----------------------------------------------

est <- as.list(sdrep, "Estimate", report = TRUE)
se <- as.list(sdrep, "Std. Error", report = TRUE)
est
se

# Question: What does the range mean here?
#           What are the units of the range? Note we're working with UTMs in km.

est2 <- as.list(sdrep, "Estimate")
est2$mu_log_omega
est2$log_omega

# Question: What is mu_log_omega?
# Question: What is log_omega?
# Question: How many log_omega's are there and why?

mesh$n
mesh$loc
dat <- data.frame(mesh$loc, log_omega = est2$log_omega)
ggplot(dat, aes(X1, X2, colour = log_omega)) + geom_point()

# Question: What does the above represent?

pred <- grid_dat
pred$omega_pred <- est$omega_pred
pred$omega_se <- se$omega_pred

ggplot(pred, aes(X, Y, fill = omega_pred)) +
  geom_raster() +
  geom_point(data = pop, mapping = aes(X, Y), inherit.aes = FALSE) +
  scale_fill_viridis_c(option = "A") +
  coord_fixed()

ggplot(pred, aes(X, Y, fill = omega_se)) +
  geom_raster() +
  geom_point(data = pop, mapping = aes(X, Y), inherit.aes = FALSE) +
  scale_fill_viridis_c(option = "B") +
  coord_fixed()

# Compare to the non-spatial model ------------------------------------------

# Compare the marginal AIC values:

AIC_non_spatial <- 2 * length(non_spatial_fit$par) - 2 * -non_spatial_fit$objective
AIC_spatial <- 2 * length(spatial_fit$par) - 2 * -spatial_fit$objective

AIC_non_spatial
AIC_spatial

# Question: Which model does marginal AIC favour as more parsimonious?

# Question: How similar are the parameter estimates?
sdrep
sdrep_non_spatial

# Plot curves ---------------------------------------------------------------

p <- as.list(sdrep, "Estimate")
L_inf <- exp(p$log_L_inf)
t0 <- p$t0
omega <- exp(p$mu_log_omega)

ages <- seq(min(pop$age), max(pop$age), length.out = 100)
lengths <- L_inf * (1 - exp(-(omega / L_inf) * (ages - t0)))

plot(pop[, c("age", "length")], col = "#00000010", pch = 19)
lines(ages, lengths, col = "red", lwd = 2)
o <- order(pop$age)
lines(pop$age[o], pred_non_spatial[o], col = "blue", lwd = 2)

# You made it to the end. Great! Try these extra-credit exercises:

# Extra exercise 1: Check the sensitivity to mesh resolution.

# Extra exercise 2: How different does the growth curve look with the maximum
# vs. minimum 'omega' spatially varying values?

# Extra exercise 3: Try estimating the growth curve and its standard error
# within the RTMB model. Estimate the curve for the average log_omega
# (without the random fields included in the prediction). Is this more
# uncertain than the model that doesn't estimate random fields?
