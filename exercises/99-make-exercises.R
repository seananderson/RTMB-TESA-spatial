setwd(here::here())
setwd("exercises")
x <- readLines("03-growth-spde.R")

x[grep("Q <- ", x)] <- "  Q <-  ### Exercise: fill in precision matrix calculation"

x[grep("log_omega %~% ", x)] <- "  log_omega %~%  ### Exercise: fill in the GRMF density with mu_log_omega and Q"

x[grep("  omega <- exp", x)] <- "  omega <- exp() ### Exercise: fill in the bilinear interpolation, to project `log_omega` using `interpolator_data`; wrap it in exp()"

x[grep("Cahill et al. 2020 ", x) + 1] <- "  pred <-  ### Exercise: fill in the growth curve, see the non-spatial growth .R file"

x[grep("range <-", x)] <- "  range <-  ### Exercise: fill in range calculation"

x[grep("sigma <-", x)] <- "  sigma <-  ### Exercise: fill in SD calculation"

# i <- grep("AIC_non_spatial <- ", x)
# j <- grep("AIC_spatial$", x)

# x <- x[c(1:(i-1), (j-1):length(x))]
# x[i] <- "AIC_non_spatial <-  ### Exercise fill in AIC calculation"
# x[i+1] <- "AIC_spatial <-  ### Exercise fill in AIC calculation"

# x[grep("mesh\\$n$", x)] <- ""
# x[grep("mesh\\$loc$", x)] <-  ""
# x[grep("ggplot\\(dat, aes\\(X1, X2", x)] <-  ""
# x[grep("dat <- data.frame\\(mesh\\$loc", x)] <-  ""

# x <- c(x, "\n# You made it to the end. Great! Try these extra-credit exercises:\n")
# x <- c(x, "# Extra exercise 1: Check the sensitivity to mesh resolution.\n")
# x <- c(x, "# Extra exercise 2: How different does the growth curve look with the maximum vs. minimum 'omega' spatially varying values?\n")
# x <- c(x, "# Extra exercise 2: Try estimating the growth curve and its standard error within the RTMB model. Estimate the curve for the average log_omega (without the random fields included in the prediction). Is this more uncertain than the model that doesn't estimate random fields?\n")

writeLines(x, "03-growth-spde-exercise.R")


x <- readLines("01-spatial-dmvnorm.R")
x[grep("correlation_matrix\\[i, j\\]", x)] <- "      correlation_matrix[i, j] <-  ### Exercise: fill in the exponential correlation function (use d and rho)"
writeLines(x, "01-spatial-dmvnorm-exercise.R")

setwd(here::here())
