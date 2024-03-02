# Goals:
# - Learn to fit spatially varying coefficient (SVC) models with an SPDE
#   approach in RTMB.
#
# See:
# Barnett, L.A.K., Ward, E.J., and Anderson, S.C. 2021. Improving estimates of
# species distribution change by incorporating local trends. Ecography 44(3):
# 427–439. https://doi.org/10.1111/ecog.05176.
# for the type of local spatially varying local trend model we are fitting.
#
# Also see:
#
# Gelfand, A.E., Kim, H.-J., Sirmans, C.F., and Banerjee, S. 2003. Spatial
# modeling with spatially varying coefficient processes. Journal of the American
# Statistical Association 98(462): 387–396.
# doi:10.1198/016214503000170.

# Thorson, J.T., Barnes, C.L., Friedman, S.T., Morano, J.L., and Siple, M.C. 2023.
# Spatially varying coefficients can improve parsimony and descriptive power for
# species distribution models. Ecography: e06510. doi:10.1111/ecog.06510.

library(RTMB)
library(ggplot2)
options(ggplot2.continuous.colour = "viridis")
theme_set(theme_light())

# Look at raw data ----------------------------------------------------------

dat <- readRDS(here::here("data/dogfish-pdo.rds"))
dat$year_centered <- dat$year - 2013
dat$pdo_centered <- dat$pdo - mean(dat$pdo)

# observations
ggplot(dat, aes(X, Y, colour = log(catch_weight))) +
  geom_point() +
  facet_wrap(~year)

ggplot(dat, aes(X, Y, colour = pdo_centered)) +
  geom_point() +
  facet_wrap(~year)

ggplot(dat, aes(X, Y, colour = year_centered)) +
  geom_point() +
  facet_wrap(~year)

# Mesh and SPDE matrix construction -----------------------------------------

# make a 2D 'mesh' object:
mesh <- fmesher::fm_mesh_2d(
  loc = as.matrix(dat[, c("X", "Y")]),
  cutoff = 5 # minimum triangle edge length
)
plot(mesh)
points(dat[, c("X", "Y")], col = "red", pch = 21, cex = 0.5)

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
  Q1 <- tau[1]^2 * (kappa^4 * spde$c0 + 2 * kappa^2 * spde$g1 + spde$g2)
  Q2 <- tau[2]^2 * (kappa^4 * spde$c0 + 2 * kappa^2 * spde$g1 + spde$g2)

  rf1 %~% dgmrf(0, Q1)
  rf2 %~% dgmrf(0, Q2)

  # project random effects from vertices to data locations:
  rf1_at_observations <- interpolator_data %*% rf1
  rf2_at_observations <- interpolator_data %*% rf2

  # create linear predictor:
  eta <- # linear predictor in link (log) space
    log(area_swept) + # 'offset' such that we're modelling catch per 1 unit area swept
    b[1] + # intercept
    rf1_at_observations + # "intercept" spatial random field
    b[2] * covariate + # overall trend
    rf2_at_observations * covariate # spatially varying coefficient trend

  # data likelihood:
  observed %~% dtweedie(
    exp(as.vector(eta)), # our expected mean (as.vector() is removing sparse matrix class)
    phi = exp(log_phi), # Tweedie dispersion parameter
    p = plogis(tweedie_p) + 1, # Tweedie power parameter; must be 1 < p < 2
    log = TRUE
  )

  # reporting:
  spatial_intercept <- b[1] + rf1_at_observations
  spatial_slope <- b[2] + rf2_at_observations
  REPORT(spatial_intercept)
  ADREPORT(spatial_slope)
}

tmb_data <- list(
  observed = dat$catch_weight,
  area_swept = dat$area_swept,
  spde = spde,
  interpolator_data = interpolator_data,
  covariate = dat$year_centered
)

par <- list(
  log_tau = rep(0, 2),
  log_kappa = 0,
  rf1 = rep(0, length = mesh$n),
  rf2 = rep(0, length = mesh$n),
  b = rep(0, 2),
  log_phi = 0,
  tweedie_p = 0
)

obj <- MakeADFun(nll, par, random = c("rf1", "rf2"))

# Fit model -----------------------------------------------------------------

opt <- nlminb(obj$par, obj$fn, obj$gr)
sdrep <- sdreport(obj)
sdrep

# Question: What does each parameter mean?

# Extract estimates ---------------------------------------------------------

pl <- as.list(sdrep, "Estimate")
plr <- as.list(sdrep, "Estimate", report = TRUE)
plrsd <- as.list(sdrep, "Std. Error", report = TRUE)
r <- obj$report(obj$env$last.par.best)

dat$intercept <- r$spatial_intercept[, 1]
dat$svc <- plr$spatial_slope[,1]
dat$svc_se <- plrsd$spatial_slope[,1]
dat$svc_lwr <- exp(dat$svc - 1.96 * dat$svc_se)
dat$svc_upr <- exp(dat$svc + 1.96 * dat$svc_se)

# Visualize the estimates ---------------------------------------------------

lims <- range(c(dat$svc_lwr, dat$svc_upr)) # for plotting
lims[2] <- min(lims[2], 5) # in case there are giant outlying values

ggplot(dat, aes(X, Y, colour = intercept)) +
  geom_point()

# Question: What does this intercept random field represent?
# At what value of our predictor is this based?

g0 <- ggplot(dat, aes(X, Y, colour = exp(svc))) +
  geom_point() +
  scale_color_gradient2(limits = lims, midpoint = 1) +
  coord_fixed()+
  ggtitle("SVC estimate")

g1 <- ggplot(dat, aes(X, Y, colour = svc_se)) +
  geom_point() +
  coord_fixed()+
  scale_fill_viridis_c(option = "C") +
  ggtitle("SVC standard error")

g2 <- ggplot(dat, aes(X, Y, colour = svc_lwr)) +
  geom_point() +
  scale_color_gradient2(limits = lims, midpoint = 1)+
  coord_fixed()+
  ggtitle("SVC lower CI")

g3 <- ggplot(dat, aes(X, Y, colour = svc_upr)) +
  geom_point() +
  scale_color_gradient2(limits = lims, midpoint = 1)+
  coord_fixed()+
  ggtitle("SVC upper CI")

cowplot::plot_grid(g0, g1, g2, g3, nrow = 2, align = "v")

# Question: What do the above plots tell us?

ggplot(dat, aes(X, Y, colour = svc_upr < 1)) +
  geom_point() +
  coord_fixed() +
  scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "grey50")) +
  ggtitle("SVC upper CI below 1")

# Question:
# If the predictor (year or PDO) increased by 1 unit, what would you expect at
# a given location in space?

# Compare same model in sdmTMB ----------------------------------------------

library(sdmTMB)
mesh0 <- make_mesh(dat, c("X", "Y"), mesh = mesh)
m <- sdmTMB(
  catch_weight ~ year_centered,
  # catch_weight ~ pdo_centered,
  spatial_varying = ~year_centered,
  # spatial_varying = ~pdo_centered,
  offset = log(dat$area_swept),
  family = tweedie(),
  data = dat,
  spatial = "on",
  mesh = mesh0
)
sanity(m)
m
m$sd_report
tidy(m, conf.int = TRUE)
tidy(m, "ran_pars", conf.int = TRUE)

logLik(m)
opt$objective

# Exercise ------------------------------------------------------------------

# Swap out `year_centered` for `pdo_centered` above and re-run the above code.
# Answer the above questions with your partner as you work through the code.
