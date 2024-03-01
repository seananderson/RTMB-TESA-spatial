library(RTMB)
pop <- readRDS(here::here("data/pop.rds"))
pop$log_length <- log(pop$length)

par <- list(
  log_L_inf = log(max(pop$length)),
  log_omega = 0,
  t0 = 0,
  log_obs_sd = 0
)

nll <- function(par) {
  getAll(par, pop)
  L_inf <- exp(log_L_inf)
  # omega: growth rate of length measurement units per year in early life
  # omega = k * L_inf
  omega <- exp(log_omega)
  # Gallucci and Quinn (1979)
  # Cahill et al. 2020 https://doi.org/10.1139/cjfas-2019-0434
  pred <- L_inf * (1 - exp(-(omega / L_inf) * (age - t0)))
  REPORT(pred)

  obs_sd <- exp(log_obs_sd)
  log_length %~% dnorm(log(pred), obs_sd)
}

obj <- MakeADFun(nll, par)

non_spatial_fit <- nlminb(obj$par, obj$fn, obj$gr)

sdrep_non_spatial <- sdreport(obj)
sdrep_non_spatial

r <- obj$report(obj$env$last.par.best)
pred_non_spatial <- r$pred

plot(pop[, c("age", "length")], col = "#00000010", pch = 19)
o <- order(pop$age)
lines(pop$age[o], pred_non_spatial[o], col = "red", lwd = 2)
