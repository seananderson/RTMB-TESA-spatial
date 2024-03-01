# Goals:
# - Gain exposure to how these models can be made spatiotemporal.

library(RTMB)
library(ggplot2)
options(ggplot2.continuous.colour = "viridis")
options(ggplot2.continuous.fill = "viridis")
theme_set(theme_light())

# Simulate some fake data ---------------------------------------------------

dat <- sdmTMB::pcod
prediction_grid <- sdmTMB::qcs_grid

# observations
ggplot(dat, aes(X, Y, colour = present)) +
  geom_point() +
  facet_wrap(~year)

# Mesh and SPDE matrix construction -----------------------------------------

# make a 2D 'mesh' object:
mesh <- fmesher::fm_mesh_2d(
  loc = as.matrix(dat[, c("X", "Y")]),
  cutoff = 10 # minimum triangle edge length
)
plot(mesh)
# compute 'finite element' matrices for SPDE approach:
spde <- fmesher::fm_fem(mesh)
# compute bilinear interpolation matrix from mesh to data:
interpolator_data <- fmesher::fm_basis(
  mesh,
  loc = as.matrix(dat[, c("X", "Y")])
)

# RTMB setup  ---------------------------------------------------------------

nll <- function(par) {
  getAll(par, tmb_data)
  tau <- exp(log_tau)
  kappa <- exp(log_kappa)

  # compute precision matrix Q, e.g. Lindgren + Rue 2015 JSS p4:
  Q <- tau^2 * (kappa^4 * spde$c0 + 2 * kappa^2 * spde$g1 + spde$g2)
  rf[,1] %~% dgmrf(0, Q)
  for (i in 2:n_t) {
    rf[,i] %~% dgmrf(rf[,i-1], Q)
  }

  # project random effects from vertices to data locations:
  rf_at_observations <- matrix(0, nrow = length(observed), ncol = ncol(rf))
  for (i in 1:n_t) {
    rf_at_observations[,i] <- as.vector(interpolator_data %*% rf[,i])
  }

  # pick out the appropriate time slice for each row of data:
  eta <- numeric(length(observed))
  for (i in 1:length(observed)) {
    eta[i] <- mu + rf_at_observations[i,t_i[i]] # linear predictor in link space
  }
  REPORT(eta)

  # data likelihood:
  observed %~% dbinom_robust(size = 1, logit_p = as.vector(eta), log = TRUE)
  # which is just a more robust version of:
  # observed %~% dbinom(size = 1, prob = plogis(eta), log = TRUE) # plogis = inverse logit

  range <- sqrt(8) / kappa
  ADREPORT(range)
  sigma <- 1 / sqrt(4 * pi * exp(2 * log_tau + 2 * log_kappa))
  ADREPORT(sigma)
}

tmb_data <- list(observed = dat$present)
tmb_data$mesh <- mesh
tmb_data$spde <- spde
tmb_data$interpolator_data <- interpolator_data
tmb_data$n_t <- length(unique(dat$year))
tmb_data$t_i <- as.numeric(as.factor(dat$year))

par <- list(
  log_tau = 0,
  log_kappa = 0,
  rf = matrix(0, ncol = tmb_data$n_t, nrow = mesh$n),
  mu = 0
)

obj <- MakeADFun(nll, par, random = "rf")

# Fit model -----------------------------------------------------------------

opt <- nlminb(obj$par, obj$fn, obj$gr)
sdrep <- sdreport(obj)
sdrep

r <- obj$report()
plr <- as.list(sdrep, "Estimate", report = TRUE)
names(plr)
plr$range
plr$sigma

dat$encounter_prob <- plogis(r$eta)
ggplot(dat, aes(X, Y, colour = present)) +
  geom_point() +
  facet_wrap(~year)

ggplot(dat, aes(X, Y, colour = encounter_prob)) +
  geom_point() +
  facet_wrap(~year) +
  scale_colour_viridis_c(option = "C")

# Same with sdmTMB ----------------------------------------------------------

library(sdmTMB)
.mesh = make_mesh(dat, c("X", "Y"), mesh = mesh)
m <- sdmTMB(
  present ~ 1,
  family = binomial(),
  data = dat,
  time = "year",
  spatial = "off",
  spatiotemporal = "rw",
  mesh = .mesh
)
m
m$sd_report
tidy(m, "ran_pars", conf.int = TRUE)

# Question:
# What is wrong with our model still?
plot(r$eta)
