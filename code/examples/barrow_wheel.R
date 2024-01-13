# Barrow wheel distribution ----------------------------------------------------

# Packages ---------------------------------------------------------------------
library("ICSClust")
library("robustX")      # rbwheel
library("ggplot2")
library("EMMIXmfa")     # MFA
library("clustvarsel")  # clustvarsel


# Data  ------------------------------------------------------------------------

n <- 1000   # number of observations
p <- 3      # number of variables
pct <- 0.2  # percentage of outliers

set.seed(072023)
X <-  rbwheel(n = n, p = p, frac = pct,
              sig1 = 0.1, sig2 = 0.2,
              fullResult = TRUE)$X
clusters <- rep("Group1", n)
clusters[((n*(1-pct))+1):n] <-
  ifelse(X[((n*(1-pct))+1):n,1] < 0, "Group2", "Group3")

data <- cbind(cluster = clusters, data.frame(X))
nb_clusters <- length(unique(data$cluster))
table(clusters)

# Plot simulated data ----------------------------------------------------------

# prepare data frame with nicer columns names
df_plot <- data[, -1]
names(df_plot) <- gsub("([0-9])", "[\\1]", names(df_plot))

# create plot of simulated data
colors <- scales::hue_pal()(5)[c(5, 2, 4)]
plot_data <- component_plot(df_plot, clusters = as.factor(data$cluster)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors)

# file name for plot
file_plot <- "figures/examples/barrow_wheel_data.%s"
# save plot to pdf
pdf(sprintf(file_plot, "pdf"), width = 8, height = 5)
print(plot_data)
dev.off()
# save plot to png
png(sprintf(file_plot, "png"), width = 8, height = 5,
    unit = "in", res = 250)
print(plot_data)
dev.off()


# Clustering on observed data --------------------------------------------------

# kmeans
set.seed(072023)
res <- kmeans_clust(data[, -1], k = 3)
table(res$clusters, data[, 1])
print("kmeans")
print(mclust::adjustedRandIndex(data[, 1], res$clusters))

# mclust
set.seed(072023)
res <- mclust_clust(data[, -1], k = 3)
table(res$clusters, data[, 1])
print("mclust")
print(mclust::adjustedRandIndex(data[, 1], res$clusters))

# rmclust
set.seed(072023)
res <- rmclust_clust(data[, -1], k = 3)
table(res$clusters, data[, 1])
print("rmclust")
print(mclust::adjustedRandIndex(data[, 1], res$clusters))


# Clustering after ICS with COV-COV4 -------------------------------------------

# kmeans
set.seed(072023)
out <- ICSClust(data[, -1], nb_select = 1, nb_clusters = 3,
                method = "kmeans_clust")
plot(out)
table(out$clusters, data[, 1])
mclust::adjustedRandIndex(data[, 1], out$clusters)

# mclust
set.seed(072023)
out <- ICSClust(data[, -1], nb_select = 1, nb_clusters = 3,
                method = "mclust_clust")
plot(out)
table(out$clusters, data[, 1])
mclust::adjustedRandIndex(data[, 1], out$clusters)

# rmclust
set.seed(072023)
out <- ICSClust(data[, -1], nb_select = 1, nb_clusters = 3,
                method = "rmclust_clust")
plot(out)
table(out$clusters, data[, 1])
mclust::adjustedRandIndex(data[, 1], out$clusters)


# Clustering after ICS with TCOV-COV -------------------------------------------

# kmeans
set.seed(072023)
out <- ICSClust(data[, -1], nb_select = 1, nb_clusters = 3,
                ICS_args = list(S1 = ICS_tcov, S2 = ICS_cov),
                method = "kmeans_clust")
plot(out)
table(out$clusters, data[, 1])
print("ICS + kmeans")
print(mclust::adjustedRandIndex(data[, 1], out$clusters))

# mclust
set.seed(072023)
out <- ICSClust(data[, -1], nb_select = 1, nb_clusters = 3,
                ICS_args = list(S1 = ICS_tcov, S2 = ICS_cov),
                method = "mclust_clust")
plot(out)
table(out$clusters, data[, 1])
print("ICS + mclust")
print(mclust::adjustedRandIndex(data[, 1], out$clusters))

# create plot of IC's
plot_IC <- component_plot(out$ICS_out, clusters = data[, 1]) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors)

# file name for plot
file_plot <- "figures/examples/barrow_wheel_IC.%s"
# save plot to pdf
pdf(sprintf(file_plot, "pdf"), width = 8, height = 5)
print(plot_IC)
dev.off()
# save plot to png
png(sprintf(file_plot, "png"), width = 8, height = 5,
    unit = "in", res = 250)
print(plot_IC)
dev.off()

# rmclust
set.seed(072023)
out <- ICSClust(data[, -1], nb_select = 1, nb_clusters = 3,
                ICS_args = list(S1 = ICS_tcov, S2 = ICS_cov),
                method = "rmclust_clust")
plot(out)
table(out$clusters, data[, 1])
print("ICS + rmclust")
print(mclust::adjustedRandIndex(data[, 1], out$clusters))


# MFA --------------------------------------------------------------------------
# McLachlan et al. (2003): implemented in function mfa() in package EMMIXmfa
set.seed(072023)
out <- mfa(Y = data[, -1], g = 3, q = 1)
table(out$clus, data[, 1])
print("MFA")
print(mclust::adjustedRandIndex(data[, 1], out$clus))


# clustvarsel ------------------------------------------------------------------
# Raftery & Dean (2006): implemented in function clustvarsel() in package clustvarsel
set.seed(072023)
res <- clustvarsel(data = data[, -1], G = 3)
table(res$model$classification, data[, 1])
print("clustvarsel")
print(mclust::adjustedRandIndex(data[, 1], res$model$classification))
