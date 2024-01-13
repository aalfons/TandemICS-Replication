# Tolerance ellipses of different scatter and shape matrices -------------------

# Packages and functions -------------------------------------------------------
library("ICSClust")
library("ggthemes")  # colorblind palette

# from function tolEllipsePlot() from rrcov
# Arguments:
# - loc: location estimate
# - cov: scatter estimate
ellipse <- function(loc, cov) {
  dist <- sqrt(qchisq(0.975, 2))
  A <- solve(cov)
  eA <- eigen(A)
  ev <- eA$values
  lambda1 <- max(ev)
  lambda2 <- min(ev)
  eigvect <- eA$vectors[, order(ev)[2]]
  z <- seq(0, 2 * pi, by = 0.01)
  z1 <- dist/sqrt(lambda1) * cos(z)
  z2 <- dist/sqrt(lambda2) * sin(z)
  alfa <- atan(eigvect[2]/eigvect[1])
  r <- matrix(c(cos(alfa), -sin(alfa), sin(alfa), cos(alfa)), ncol = 2)
  t(loc + t(cbind(z1, z2) %*% r))
}

# function for parsing axis labels
parse_labels <- function(labels, ...) parse(text = labels, ...)


# Example ellipses -------------------------------------------------------------

# set seed of random number generator
set.seed(20220907)

# simulate data
pct_clusters <- c(1/3, 1/3, 1/3)
df <- mixture_sim(pct_clusters = pct_clusters, n = 1000, p = 2, delta = 10)

# graphical parameters
colors <- scales::hue_pal()(5)[c(2, 4, 5)]
names(colors) <- c("Group1", "Group2", "Group3")
col_vect <- colors[df$cluster]

# compute the ellipse for the sample mean and the sample covariance matrix
X <- df[, -1]
m <- colMeans(X, na.rm = TRUE)
S <- var(X, na.rm = TRUE)
ellipse_cov <- ellipse(loc = m, cov = S)


# Local scatters -----

# open png file
png("figures/examples/ellipses_local.png", width = 5, height = 5,
    units = "in", res = 250)
par(mar = c(4, 5, 0.5, 0.1))

# plot the simulated data
plot(X, col = col_vect, cex.axis = 1.25, cex.lab = 1.5,
     xlab = parse_labels("X[1]"), ylab = parse_labels("X[2]"),
     xlim = range(ellipse_cov[,1])*2, ylim = range(ellipse_cov[,2])*1.5)

# add tolerance ellipse for the sample covariance matrix
points(ellipse_cov, type = "l", col = "black")

# compute scatters
scatters <-  c("ICS_mlc", "ICS_lcov",  "ICS_tcov")
scatters_labels <- c( "MLC", "LCOV", "TCOV")
scatters_args_list <- list(`ICS_mlc` = list(location = TRUE, maxiter = 1000),
                           `ICS_lcov` = list(),
                           `ICS_tcov` = list())
n_scatters <- length(scatters)
col_scatters <- ggthemes::colorblind_pal()(n_scatters+1)[-1]

# loop over scatters
for (i in seq_along(scatters) ){
  # compute scatter estimate
  estimate <- do.call(scatters[i],
                      append(list(x = X), scatters_args_list[[i]]))
  print(scatters[i])
  print(estimate$scatter)
  # make sure we have a location estimate
  if (is.null(estimate$location)) {
    # use sample mean for LCOV and TCOV
    estimate$location <- m
  }
  # compute ellipse and add it to plot
  current_ellipse <- ellipse(loc = estimate$location, cov = estimate$scatter)
  points(current_ellipse, type = "l", col = col_scatters[i])
}

# add legend
legend("topright", legend = c("COV", parse_labels(scatters_labels)),
       col = c("black", col_scatters), lwd = 2, cex = 1.25)

# close file
dev.off()


# Global scatters -----

# open png file
png("figures/examples/ellipses_global.png", width = 5, height = 5,
    units = "in", res = 250)
par(mar = c(4, 5, 0.5, 0.1))

# plot the simulated data
plot(X, col = col_vect, cex.axis = 1.25, cex.lab = 1.5,
     xlab = parse_labels("X[1]"), ylab = parse_labels("X[2]"),
     xlim = range(ellipse_cov[,1])*2, ylim = range(ellipse_cov[,2])*1.5)

# add tolerance ellipse for the sample covariance matrix
points(ellipse_cov, type = "l", col = "black")

# compute scatters
scatters <-  c("ICS_cov4", "ICS_ucov")
scatters_labels <- c("COV[4]", "UCOV")
n_scatters <- length(scatters)
col_scatters <- ggthemes::colorblind_pal()(n_scatters+1)[-1]

# loop over scatters
for (i in seq_along(scatters) ){
  # compute scatter estimate
  estimate <- do.call(scatters[i], list(x = X))
  print(scatters[i])
  print(estimate$scatter)
  # compute ellipse and add it to plot
  current_ellipse <- ellipse(loc = estimate$location, cov = estimate$scatter)
  points(current_ellipse, type = "l", col = col_scatters[i])
}

# add legend
legend("topright", legend = c("COV", parse_labels(scatters_labels)),
       col = c("black", col_scatters), lwd = 2, cex = 1.25)

# close file
dev.off()


# MCD -----

# open png file
png("figures/examples/ellipses_MCD.png", width = 5, height = 5,
    units = "in", res = 250)
par(mar = c(4, 5, 0.5, 0.1))

# plot the simulated data
plot(X, col = col_vect, cex.axis = 1.25, cex.lab = 1.5,
     xlab = parse_labels("X[1]"), ylab = parse_labels("X[2]"),
     xlim = range(ellipse_cov[,1])*2, ylim = range(ellipse_cov[,2])*1.5)

# add tolerance ellipse for the sample covariance matrix
points(ellipse_cov, type = "l", col = "black")

# subset sizes
alphas <- c(0.75, 0.50, 0.25, 0.10)
n_scatters <- length(alphas)
col_scatters <- ggthemes::colorblind_pal()(n_scatters+2)[-c(1, 5)]

# loop over subset sizes
for (i in seq_along(alphas) ){
  # compute scatter estimate
  estimate <- ICS_mcd_raw(X, location = TRUE, nsamp = 500, alpha = alphas[i])
  print(paste("MCD", alphas[i], sep = "_"))
  print(estimate$scatter)
  # compute ellipse and add it to plot
  current_ellipse <- ellipse(loc = estimate$location, cov = estimate$scatter)
  points(current_ellipse, type = "l", col = col_scatters[i])
}

# add legend
legend("topright", legend = parse_labels(c("COV", paste0("MCD[", alphas, "]"))),
       col = c("black", col_scatters), lwd = 2, cex = 1.25)

# close file
dev.off()


# RMCD -----

# open png file
png("figures/examples/ellipses_RMCD.png", width = 5, height = 5,
    units = "in", res = 250)
par(mar = c(4, 5, 0.5, 0.1))

# plot the simulated data
plot(X, col = col_vect, cex.axis = 1.25, cex.lab = 1.5,
     xlab = parse_labels("X[1]"), ylab = parse_labels("X[2]"),
     xlim = range(ellipse_cov[,1])*2, ylim = range(ellipse_cov[,2])*1.5)

# add tolerance ellipse for the sample covariance matrix
points(ellipse_cov, type = "l", col = "black")

# subset sizes
alphas <- c(0.75, 0.50, 0.25, 0.10)
n_scatters <- length(alphas)
col_scatters <- ggthemes::colorblind_pal()(n_scatters+2)[-c(1, 5)]

# loop over subset sizes
for (i in seq_along(alphas) ){
  # compute scatter estimate
  estimate <- ICS_mcd_rwt(X, location = TRUE, nsamp = 500, alpha = alphas[i])
  print(paste("RMCD", alphas[i], sep = "_"))
  print(estimate$scatter)
  # compute ellipse and add it to plot
  current_ellipse <- ellipse(loc = estimate$location, cov = estimate$scatter)
  points(current_ellipse, type = "l", col = col_scatters[i])
}

# add legend
legend("topright", legend = parse_labels(c("COV", paste0("RMCD[", alphas, "]"))),
       col = c("black", col_scatters), lwd = 2, cex = 1.25)

# close file
dev.off()
