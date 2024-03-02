# Goals:
# - See how we can add a spatial random field to our spatiotemporal model

library(RTMB)
library(ggplot2)
options(ggplot2.continuous.colour = "viridis")
options(ggplot2.continuous.fill = "viridis")
theme_set(theme_light())

# Load data -----------------------------------------------------------------

dat <- sdmTMB::pcod

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
  tau0 <- exp(log_tau0)
  tau <- exp(log_tau)
  kappa <- exp(log_kappa)

  # compute precision matrix Q, e.g. Lindgren + Rue 2015 JSS p4:
  Q0 <- tau0^2 * (kappa^4 * spde$c0 + 2 * kappa^2 * spde$g1 + spde$g2) #< NEW
  Q <- tau^2 * (kappa^4 * spde$c0 + 2 * kappa^2 * spde$g1 + spde$g2)

  rf0 %~% dgmrf(0, Q0) #< NEW
  rf[,1] %~% dgmrf(0, Q)
  for (i in seq(2, n_t)) {
    rf[,i] %~% dgmrf(rf[,i-1], Q)
  }

  # project random effects from vertices to data locations:
  rf_at_observations0 <- interpolator_data %*% rf0  #< NEW
  rf_at_observations <- matrix(0, nrow = length(observed), ncol = n_t)
  for (i in seq_len(n_t)) {
    rf_at_observations[,i] <- as.vector(interpolator_data %*% rf[,i])
  }

  # pick out the appropriate time slice for each row of data
  # and add on any other components:
  eta <- numeric(length(observed))
  for (i in seq_along(observed)) {
    eta[i] <- mu + rf_at_observations0[i] + #< NEW
      rf_at_observations[i,t_i[i]]
  }
  REPORT(rf_at_observations0) #< NEW
  REPORT(rf_at_observations)
  REPORT(eta)

  # data likelihood:
  observed %~% dbinom_robust(size = 1, logit_p = eta, log = TRUE)
  # which is just a more robust version of:
  # observed %~% dbinom(size = 1, prob = plogis(eta), log = TRUE)

  range <- sqrt(8) / kappa
  ADREPORT(range)
  sigma0 <- 1 / sqrt(4 * pi * exp(2 * log_tau0 + 2 * log_kappa))
  sigma <- 1 / sqrt(4 * pi * exp(2 * log_tau + 2 * log_kappa))
  ADREPORT(sigma0) #< NEW
  ADREPORT(sigma)
}

tmb_data <- list(
  observed = dat$present,
  spde = spde,
  interpolator_data = interpolator_data,
  n_t = length(unique(dat$year)),
  t_i = as.numeric(as.factor(dat$year))
)

par <- list(
  log_tau0 = 0, #< NEW
  log_tau = 0,
  log_kappa = 0,
  rf0 = rep(0, length = mesh$n), #< NEW
  rf = matrix(0, ncol = tmb_data$n_t, nrow = mesh$n),
  mu = 0
)

obj <- MakeADFun(nll, par, random = c("rf0", "rf"))  #< NEW: rf0 added

# Fit model -----------------------------------------------------------------

opt <- nlminb(obj$par, obj$fn, obj$gr)
sdrep <- sdreport(obj)
sdrep

plr <- as.list(sdrep, "Estimate", report = TRUE)
names(plr)
plr$range
plr$sigma0
plr$sigma
r <- obj$report()

dat$encounter_prob <- plogis(r$eta)
dat$rf_spatial <- r$rf_at_observations0[,1]

ggplot(dat, aes(X, Y, colour = present)) +
  geom_point() +
  facet_wrap(~year)

ggplot(dat, aes(X, Y, colour = encounter_prob)) +
  geom_point() +
  facet_wrap(~year)

ggplot(dat, aes(X, Y, colour = rf_spatial)) +
  geom_point() +
  facet_wrap(~year)

ggplot(dat, aes(X, Y, colour = rf_spatial)) +
  geom_point() +
  facet_wrap(~year)

Matrix::image(obj$env$spHess(random = TRUE))

# Compare same model in sdmTMB ----------------------------------------------

library(sdmTMB)
mesh0 <- make_mesh(dat, c("X", "Y"), mesh = mesh)
m <- sdmTMB(
  present ~ 1,
  family = binomial(),
  data = dat,
  time = "year",
  spatial = "on", #< NEW
  spatiotemporal = "rw",
  mesh = mesh0
)
m
m$sd_report
tidy(m, "ran_pars", conf.int = TRUE)

plot(r$eta)
