# *********************************************
# Authors: Aurore Archimbaud and Andreas Alfons
#          Erasmus University Rotterdam
# *********************************************


# Pairwise distances when structure is in a low-dimensional subspace -----------

# Packages and functions -------------------------------------------------------
library("ICSClust")
library("biotools")
library("ggplot2")
library("fpc")
library("magrittr")

# function for parsing axis labels
parse_labels <- function(labels, ...) parse(text = labels, ...)


# Data -------------------------------------------------------------------------

# set seed of random number generator
set.seed(20220907)

# simulate data from a mixture of 2 Gaussian distributions
pct_clusters <- c(0.85, 0.15)
X_gen <- mixture_sim(pct_clusters = pct_clusters, n = 1000, p = 2, delta = 10)
X <- X_gen[, 2:3]
true_groups <- X_gen$cluster


# standardization
Xstd <- scale(X, center = TRUE, scale = TRUE)
colnames(Xstd) <- paste0("Std.X", 1:ncol(Xstd), "")

# create plot of standardized data
text_size_factor <- 20/6.5
plot_Xstd <- as.data.frame(Xstd) %>%
  ggplot(aes(x = Std.X1, y = Std.X2)) +
  geom_point() +
  theme_bw() +
  theme(axis.title = element_text(size = 11 * text_size_factor),
        axis.text = element_text(size = 9 * text_size_factor),
        aspect.ratio = 1) +
  labs(x = expression(Std~X[1]), y = expression(Std~X[2]))

# file name for plot
file_plot <- "figures/examples/Xstd.%s"
# # save plot to pdf
# pdf(sprintf(file_plot, "pdf"), width = 6, height = 6)
# print(plot_Xstd)
# dev.off()
# save plot to png
png(sprintf(file_plot, "png"), width = 6, height = 6,
    unit = "in", res = 250)
print(plot_Xstd)
dev.off()


# Clustering with kmeans -------------------------------------------------------

## On standardized data ----

# apply kmeans
res_kmeans <- kmeansCBI(data = Xstd, krange = 2, runs = 1)
res_kmeans$partition

# create plot of standardized data with cluster information
plot_kmeans_Xstd <- data.frame(Xstd,
                               kmeans_Xstd = factor(res_kmeans$partition) ) %>%
  ggplot(aes(x = Std.X1, y = Std.X2, color = kmeans_Xstd)) +
  geom_point(show.legend = FALSE) +
  theme_bw() +
  theme(axis.title = element_text(size = 11 * text_size_factor),
        axis.text = element_text(size = 9 * text_size_factor),
        aspect.ratio = 1) +
  labs(x = expression(Std~X[1]), y = expression(Std~X[2]))

# file name for plot
file_plot <- "figures/examples/kmeans_Xstd.%s"
# # save plot to pdf
# pdf(sprintf(file_plot, "pdf"), width = 6, height = 6)
# print(plot_kmeans_Xstd)
# dev.off()
# save plot to png
png(sprintf(file_plot, "png"), width = 6, height = 6,
    unit = "in", res = 250)
print(plot_kmeans_Xstd)
dev.off()

# compute ARI
table(true_groups, res_kmeans$partition)
mclust::adjustedRandIndex(true_groups, res_kmeans$partition)


## On reduced data: PCA ----

# apply PCA
res_PCA <- princomp(Xstd)
colnames(res_PCA$scores) <- gsub("Comp", "PC", colnames(res_PCA$scores))

# apply kmeans
res_kmeans_PCA <- kmeansCBI(data = res_PCA$scores[, 1], krange = 2)
res_kmeans_PCA$partition

# create plot of principal components
plot_PCA <- as.data.frame(res_PCA$scores) %>%
  ggplot(aes(x = PC.1, y = PC.2)) +
  geom_point() +
  theme_bw() +
  theme(axis.title = element_text(size = 11 * text_size_factor),
        axis.text = element_text(size = 9 * text_size_factor),
        aspect.ratio = 1) +
  labs(x = parse_labels("PC[1]"), y = parse_labels("PC[2]"))

# file name for plot
file_plot <- "figures/examples/PCA.%s"
# # save plot to pdf
# pdf(sprintf(file_plot, "pdf"), width = 6, height = 6)
# print(plot_PCA)
# dev.off()
# save plot to png
png(sprintf(file_plot, "png"), width = 6, height = 6,
    unit = "in", res = 250)
print(plot_PCA)
dev.off()

# create plot of standardized data with cluster information
plot_kmeans_PCA <- data.frame(Xstd,
                              kmeans_PCA = factor(res_kmeans_PCA$partition)) %>%
  ggplot(aes(x = Std.X1, y = Std.X2, color = kmeans_PCA)) +
  geom_point(show.legend = FALSE) +
  theme_bw() +
  theme(axis.title = element_text(size = 11 * text_size_factor),
        axis.text = element_text(size = 9 * text_size_factor),
        aspect.ratio = 1) +
  labs(x = expression(Std~X[1]), y = expression(Std~X[2]))

# file name for plot
file_plot <- "figures/examples/kmeans_PCA.%s"
# # save plot to pdf
# pdf(sprintf(file_plot, "pdf"), width = 6, height = 6)
# print(plot_kmeans_PCA)
# dev.off()
# save plot to png
png(sprintf(file_plot, "png"), width = 6, height = 6,
    unit = "in", res = 250)
print(plot_kmeans_PCA)
dev.off()


## On reduced data: ICS ----

# apply ICS
res_ICS <- ICS(Xstd)

# apply kmeans
res_kmeans_ICS <- kmeansCBI(data = components(res_ICS, select = 1), krange = 2)
res_kmeans_ICS$partition

# create plot of invariant coordinates
plot_ICS <- as.data.frame(components(res_ICS)) %>%
  ggplot(aes(x = IC.1, y = IC.2)) +
  geom_point() +
  theme_bw() +
  theme(axis.title = element_text(size = 11 * text_size_factor),
        axis.text = element_text(size = 9 * text_size_factor),
        aspect.ratio = 1) +
  labs(x = parse_labels("IC[1]"), y = parse_labels("IC[2]"))

# file name for plot
file_plot <- "figures/examples/ICS.%s"
# # save plot to pdf
# pdf(sprintf(file_plot, "pdf"), width = 6, height = 6)
# print(plot_ICS)
# dev.off()
# save plot to png
png(sprintf(file_plot, "png"), width = 6, height = 6,
    unit = "in", res = 250)
print(plot_ICS)
dev.off()

# create plot of standardized data with cluster information
plot_kmeans_ICS <- data.frame(Xstd,
                              kmeans_ICS = factor(res_kmeans_ICS$partition)) %>%
  ggplot(aes(x = Std.X1, y = Std.X2, color = kmeans_ICS)) +
  geom_point(show.legend = FALSE) +
  theme_bw() +
  theme(axis.title = element_text(size = 11 * text_size_factor),
        axis.text = element_text(size = 9 * text_size_factor),
        aspect.ratio = 1) +
  labs(x = expression(Std~X[1]), y = expression(Std~X[2]))

# file name for plot
file_plot <- "figures/examples/kmeans_ICS.%s"
# # save plot to pdf
# pdf(sprintf(file_plot, "pdf"), width = 6, height = 6)
# print(plot_kmeans_ICS)
# dev.off()
# save plot to png
png(sprintf(file_plot, "png"), width = 6, height = 6,
    unit = "in", res = 250)
print(plot_kmeans_ICS)
dev.off()


# Distances ----

# On standardized data ----

# compute Euclidean distances
distances_Xstd <- as.vector(dist(Xstd, method = "euclidean"))

# create plot of distances
plot_distances_Xstd <- data.frame(index = seq_along(distances_Xstd),
                                  distance = distances_Xstd) %>%
  ggplot(aes(x = index, y = distance)) +
  geom_point(alpha = 0.2, shape = 16) +
  theme_bw() +
  theme(axis.title = element_text(size = 11 * text_size_factor),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9 * text_size_factor),
        axis.ticks.x = element_blank(),
        aspect.ratio = 1) +
  labs(y = "Distances", x = "Index")

# file name for plot
file_plot <- "figures/examples/distances_Xstd.%s"
# # save plot to pdf
# pdf(sprintf(file_plot, "pdf"), width = 6, height = 6)
# print(plot_distances_Xstd)
# dev.off()
# save plot to png
png(sprintf(file_plot, "png"), width = 6, height = 6,
    unit = "in", res = 250)
print(plot_distances_Xstd)
dev.off()


# On reduced data: PCA ----

# compute Euclidean distances
distances_PCA <- as.vector(dist(res_PCA$scores[, 1], method = "euclidean"))

# create plot of distances
plot_distances_PCA <- data.frame(index = seq_along(distances_PCA),
                                  distance = distances_PCA) %>%
  ggplot(aes(x = index, y = distance)) +
  geom_point(alpha = 0.2, shape = 16) +
  theme_bw() +
  theme(axis.title = element_text(size = 11 * text_size_factor),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9 * text_size_factor),
        axis.ticks.x = element_blank(),
        aspect.ratio = 1) +
  labs(y = "Distances after PCA", x = "Index")

# file name for plot
file_plot <- "figures/examples/distances_PCA.%s"
# # save plot to pdf
# pdf(sprintf(file_plot, "pdf"), width = 6, height = 6)
# print(plot_distances_PCA)
# dev.off()
# save plot to png
png(sprintf(file_plot, "png"), width = 6, height = 6,
    unit = "in", res = 250)
print(plot_distances_PCA)
dev.off()


# On reduced data: ICS ----

# compute Euclidean distances
distances_ICS <- as.vector(dist(components(res_ICS, select = 1),
                                method = "euclidean"))

# create plot of distances
plot_distances_ICS <- data.frame(index = seq_along(distances_ICS),
                                 distance = distances_ICS) %>%
  ggplot(aes(x = index, y = distance)) +
  geom_point(alpha = 0.2, shape = 16) +
  theme_bw() +
  theme(axis.title = element_text(size = 11 * text_size_factor),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9 * text_size_factor),
        axis.ticks.x = element_blank(),
        aspect.ratio = 1) +
  labs(y = "Distances after ICS", x = "Index")

# file name for plot
file_plot <- "figures/examples/distances_ICS.%s"
# # save plot to pdf
# pdf(sprintf(file_plot, "pdf"), width = 6, height = 6)
# print(plot_distances_ICS)
# dev.off()
# save plot to png
png(sprintf(file_plot, "png"), width = 6, height = 6,
    unit = "in", res = 250)
print(plot_distances_ICS)
dev.off()


# ICS with other scatter pairs -------------------------------------------------

# check that results for ICS do not depend on the choice of scatter pair

# control parameters for ICS scatters
ICS_scatters_list <- list(
  `COV-COV4` = list(S1 = ICS_cov, S2 = ICS_cov4),
  `MLC-COV` = list(S1 = ICS_mlc, S2 = ICS_cov),
  `LCOV-COV` = list(S1 = ICS_lcov, S2 = ICS_cov,
                    S1_args = list(mscatter = "cov", proportion = 0.1)),
  `TCOV-COV` = list(S1 = ICS_tcov, S2 = ICS_cov,
                    S1_args = list(beta = 2)),
  `TCOV-UCOV` = list(S1 = ICS_tcov, S2 = ICS_ucov,
                     S1_args = list(beta = 2), S2_args = list(beta = 0.2)),
  `MCD25-MCD95` = list(S1 = ICS_mcd_raw, S2 = ICS_mcd_raw,
                       S1_args = list(alpha = 0.25, nsamp = 500),
                       S2_args = list(alpha = 0.95, nsamp = 500)),
  `MCD10-COV` = list(S1 = ICS_mcd_raw, S2 = ICS_cov,
                     S1_args = list(alpha = 0.10, nsamp = 500)),
  `MCD15-COV` = list(S1 = ICS_mcd_raw, S2 = ICS_cov,
                     S1_args = list(alpha = 0.15, nsamp = 500)),
  `MCD20-COV` = list(S1 = ICS_mcd_raw, S2 = ICS_cov,
                     S1_args = list(alpha = 0.20, nsamp = 500)),
  `MCD25-COV` = list(S1 = ICS_mcd_raw, S2 = ICS_cov,
                     S1_args = list(alpha = 0.25, nsamp = 500)),
  `MCD50-COV` = list(S1 = ICS_mcd_raw, S2 = ICS_cov,
                     S1_args = list(alpha = 0.50, nsamp = 500)),
  `MCD75-COV` = list(S1 = ICS_mcd_raw, S2 = ICS_cov,
                     S1_args = list(alpha = 0.75, nsamp = 500)),
  `RMCD25-RMCD95` = list(S1 = ICS_mcd_rwt, S2 = ICS_mcd_rwt,
                         S1_args = list(alpha = 0.25, nsamp = 500),
                         S2_args = list(alpha = 0.95, nsamp = 500)),
  `RMCD10-COV` = list(S1 = ICS_mcd_rwt, S2 = ICS_cov,
                      S1_args = list(alpha = 0.10, nsamp = 500)),
  `RMCD20-COV` = list(S1 = ICS_mcd_rwt, S2 = ICS_cov,
                      S1_args = list(alpha = 0.20, nsamp = 500)),
  `RMCD25-COV` = list(S1 = ICS_mcd_rwt, S2 = ICS_cov,
                      S1_args = list(alpha = 0.25, nsamp = 500)),
  `RMCD50-COV` = list(S1 = ICS_mcd_rwt, S2 = ICS_cov,
                      S1_args = list(alpha = 0.50, nsamp = 500)),
  `RMCD75-COV` =  list(S1 = ICS_mcd_rwt, S2 = ICS_cov,
                       S1_args = list(alpha = 0.75, nsamp = 500))
)


# loop over scatter pairs and apply kmeans after ICS
res_list <- mapply(function(scatter, ICS_scatters) {

  # apply ICS
  res_ICS <- do.call(ICS, append(list(X = Xstd), ICS_scatters))

  # apply kmeans
  res_kmeans_ICS <- kmeansCBI(data = components(res_ICS, select = 1),
                              krange = 2)

  # # create plot of invariant coordinates
  # plot_ICS <- as.data.frame(components(res_ICS)) %>%
  #   ggplot(aes(x = IC.1, y = IC.2)) +
  #   geom_point() +
  #   theme_bw() +
  #   theme(axis.title = element_text(size = 11 * text_size_factor),
  #         axis.text = element_text(size = 9 * text_size_factor),
  #         aspect.ratio = 1) +
  #   labs(x = parse_labels("IC[1]"), y = parse_labels("IC[2]"))
  #
  # # file name for plot
  # file_plot <- "figures/examples/ICS_%s.%s"
  # # # save plot to pdf
  # # pdf(sprintf(file_plot, scatter, "pdf"), width = 6, height = 6)
  # # print(plot_ICS)
  # # dev.off()
  # # save plot to png
  # png(sprintf(file_plot, scatter, "png"), width = 6, height = 6,
  #     unit = "in", res = 250)
  # print(plot_ICS)
  # dev.off()
  #
  # # create plot of standardized data with cluster information
  # plot_kmeans_ICS <- data.frame(Xstd,
  #                               kmeans_ICS = factor(res_kmeans_ICS$partition)) %>%
  #   ggplot(aes(x = Std.X1, y = Std.X2, color = kmeans_ICS)) +
  #   geom_point(show.legend = FALSE) +
  #   theme_bw() +
  #   theme(axis.title = element_text(size = 11 * text_size_factor),
  #         axis.text = element_text(size = 9 * text_size_factor),
  #         aspect.ratio = 1) +
  #   labs(x = expression(Std~X[1]), y = expression(Std~X[2]))
  #
  # # file name for plot
  # file_plot <- "figures/examples/kmeans_ICS_%s.%s"
  # # # save plot to pdf
  # # pdf(sprintf(file_plot, scatter, "pdf"), width = 6, height = 6)
  # # print(plot_kmeans_ICS)
  # # dev.off()
  # # save plot to png
  # png(sprintf(file_plot, scatter, "png"), width = 6, height = 6,
  #     unit = "in", res = 250)
  # print(plot_kmeans_ICS)
  # dev.off()

  # compute adjusted Rand index
  ARI <- mclust::adjustedRandIndex(true_groups, res_kmeans_ICS$partition)
  data.frame(scatter = scatter, ARI = ARI)

}, scatter = names(ICS_scatters_list), ICS_scatters = ICS_scatters_list,
SIMPLIFY = FALSE, USE.NAMES = FALSE)

# combine results into one data frame
res_df <-  do.call(rbind, res_list)
xtable::xtable(res_df)
