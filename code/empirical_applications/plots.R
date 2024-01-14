
# Packages and functions -------------------------------------------------------

# preprocessing
library("magrittr")
library("dplyr")
library("stringr")

# visualization
library("ICSClust")
library("ggplot2")
library("ggthemes")  # scale_color_colorblind ()
library("scales")
#library("ggh4x")     #  facet_nested()

# latex tables
library("xtable")

# data
library("MASS")      # crabs data
library("cellWise")  # philips data

# function to parse axis labels
parse_labels <- function(text, ...) parse(text = text, ...)


# Import results ---------------------------------------------------------------

# load results from files
file_names <- c("results/empirical_applications/results_crabs_log.RData",
                "results/empirical_applications/results_iris.RData",
                "results/empirical_applications/results_philips.RData")
res_df <- do.call(rbind, lapply(file_names, function(file_name) {
  load(file_name)
  results
})) %>%
  mutate(
    # recode names of selection criteria
    criterion = dplyr::recode(criterion,
                              "normal" = "Normal",
                              "med" = "Med",
                              "var" = "Var",
                              "discriminatory" = "Oracle"),
    # recode names of method labels
    method = dplyr::recode(method,
                           "pam" = "PAM",
                           "reduced_kmeans" = "reduced~kmeans",
                           "factorial_kmeans" = "factorial~kmeans",
                           "mfa" = "MFA"),
    # recode names of data sets
    name = dplyr::recode(name,
                         "crabs_log" = "Crabs",
                         "iris" = "Iris",
                         "philips" = "Philips"),
    # add variable indicating ICS or PCA
    type = ifelse(scatter %in% c("COV", "RMCD[0.75]"), "PCA",
                  ifelse(scatter == "Observed~data", "", "ICS"))
  )

# select results on initial data
res_df_initial <- res_df %>%
  filter(scatter == "Observed~data") %>%
  dplyr::select(name, method, ARI, eta2) %>%
  rename(ARI_initial = ARI, eta2_initial = eta2)

# evaluation measure to be plotted
measure <- "ARI"

# labels to be used in plots
clusters_label <- "Cluster sizes"
criterion_label <- "Component selection"
measure_label <- "Adjusted Rand Index"

# colors to be used in plots (following previous plots)
colors_ICS <- ggthemes::colorblind_pal()(4)[c(3, 1, 2, 4)]
names(colors_ICS) <- c("Oracle", "Var", "Med", "Normal")
colors_PCA <- scales::hue_pal()(4)[c(4, 3)]
names(colors_PCA) <- c("80%", "k-1")
colors <- c(colors_ICS, colors_PCA)


# Main results -----------------------------------------------------------------

# filter on data sets, component selection criteria and scatter matrices
# (based on simulations results), and clustering methods
keep_crit <- c("Med", "Var", "Normal", "80%")
keep_scatter <- c("LCOV-COV", "TCOV-COV",  "TCOV-UCOV",
                  "MCD[0.10]-COV", "MCD[0.25]-COV", "MCD[0.50]-COV",
                  "RMCD[0.75]")
keep_method <- c("PAM", "kmeans", "tkmeans", "mclust", "rmclust")

# select subset of results and merge with results on initial data
res_df_selected <- res_df %>%
  filter(criterion %in% keep_crit,
         scatter %in% keep_scatter,
         method %in% keep_method) %>%
  merge(res_df_initial,
        by = c("name", "method"))

# create plot
plot_selected <-
  res_df_selected %>%
  mutate(criterion = factor(criterion, levels = keep_crit),
         scatter = factor(scatter, levels = keep_scatter),
         method = factor(method, levels = keep_method)) %>%
  ggplot(mapping = aes_string(x = "scatter", y = measure, color = "criterion")) +
  geom_hline(aes(yintercept = ARI_initial))+
  geom_point(position = position_jitterdodge(jitter.width = 0,
                                             jitter.height = 0,
                                             dodge.width = 2/3),
             size = 2) +
  coord_flip() +
  scale_x_discrete(limits = rev, labels = parse_labels) +
  scale_color_manual(name = criterion_label,  values = colors[keep_crit]) +
  theme_bw() +
  theme(axis.title = element_text(size = 11),
        axis.text = element_text(size = 9),
        legend.position = "top",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 9),
        strip.text = element_text(size = 10)) +
  labs(x = NULL, y = measure_label) +
  facet_grid(method ~ name)

# file name for plot
file_plot = "figures/empirical_applications/%s_selected.%s"
# save plot to pdf
pdf(sprintf(file_plot, measure, "pdf"), width = 6.5, height = 5.75)
print(plot_selected)
dev.off()
# save plot to png
png(sprintf(file_plot, measure, "png"), width = 6.5, height = 5.75,
    unit = "in", res = 250)
print(plot_selected)
dev.off()


# Comparison with methods that integrate clustering and dimension reduction ----

# filter on data sets, component selection criteria and scatter matrices
# (based on simulations results), and clustering methods
keep_crit <- "Med"
keep_scatter <- c("Observed~data", "TCOV-COV", "MCD[0.50]-COV")
keep_method <- c("PAM", "kmeans", "tkmeans", "mclust", "rmclust",
                 "reduced~kmeans", "factorial~kmeans", "MFA", "clustvarsel")

# select subset of results
res_df_selected <- res_df %>%
  filter(is.na(criterion) | criterion %in% keep_crit,
         scatter %in% keep_scatter,
         method %in% keep_method)

# additional color to be used in plot
colors <- c("black", colors[keep_crit], colors[keep_crit])
names(colors) <- keep_scatter

# different plot shapes to be used in plot
shapes <- c(15, 19, 17)
names(shapes) <- keep_scatter

# create plot
plot_comparison <- res_df_selected %>%
  mutate(criterion = factor(criterion, levels = keep_crit),
         scatter = factor(scatter, levels = keep_scatter),
         method = factor(method, levels = keep_method)) %>%
  ggplot(mapping = aes_string(x = "method", y = measure, color = "scatter",
                              shape = "scatter")) +
  geom_point(position = position_jitterdodge(jitter.width = 0,
                                             jitter.height = 0,
                                             dodge.width = 2/3),
             size = 2) +
  coord_flip() +
  scale_x_discrete(limits = rev, labels = parse_labels) +
  scale_color_manual(name = "", values = colors, labels = parse_format()) +
  scale_shape_manual(name = "", values = shapes, labels = parse_format()) +
  theme_bw() +
  theme(axis.title = element_text(size = 11),
        axis.text = element_text(size = 9),
        legend.position = "top",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 9),
        strip.text = element_text(size = 10)) +
  labs(x = NULL, y = measure_label) +
  facet_grid(. ~ name)

# file name for plot
file_plot = "figures/empirical_applications/%s_comparison.%s"
# save plot to pdf
pdf(sprintf(file_plot, measure, "pdf"), width = 6, height = 3)
print(plot_comparison)
dev.off()
# save plot to png
png(sprintf(file_plot, measure, "png"), width = 6, height = 3,
    unit = "in", res = 250)
print(plot_comparison)
dev.off()


# Data sets --------------------------------------------------------------------

# colors for different clusters
colors <- scales::hue_pal()(5)[c(2, 4, 5, 1)]

## Crabs data -----

# Load data
data("crabs", package = "MASS")

# Log transformation, no standardization
X <- data.frame(cluster = factor(paste(crabs$sp, crabs$sex, sep = "_")),
                apply(crabs[, 4:8], 2, log))

# ICS with TCOV-COV
out <- ICS(X[,-1], S1 = ICS_tcov, S2 = ICS_cov, S1_args = list(beta = 2))

# create plot
plot_crabs <- component_plot(out, select = c(1:2, 4:5), clusters = X$cluster) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors)

# file name for plot
file_plot = "figures/empirical_applications/crabs_%s-%s.%s"
# save plot to pdf
pdf(sprintf(file_plot, out$S1_label, out$S2_label, "pdf"),
    width = 8, height = 5)
print(plot_crabs)
dev.off()
# save plot to png
png(sprintf(file_plot, out$S1_label, out$S2_label, "png"),
    width = 8, height = 5, unit = "in", res = 250)
print(plot_crabs)
dev.off()


## Iris data -----

# Load data
data("iris", package = "datasets")
X <- data.frame(cluster = iris$Species, iris[,1:4])

# ICS with TCOV-COV
out <- ICS(X[,-1], S1 = ICS_tcov, S2 = ICS_cov, S1_args = list(beta = 2))

# create plot
plot_iris <- component_plot(out, clusters = X$cluster) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors)

# file name for plot
file_plot = "figures/empirical_applications/iris_%s-%s.%s"
# save plot to pdf
pdf(sprintf(file_plot, out$S1_label, out$S2_label, "pdf"),
    width = 8, height = 5)
print(plot_iris)
dev.off()
# save plot to png
png(sprintf(file_plot, out$S1_label, out$S2_label, "png"),
    width = 8, height = 5, unit = "in", res = 250)
print(plot_iris)
dev.off()


## Philips data -----

# Load data
data("data_philips", package = "cellWise")

# Define clusters: based on paper: that observations 491 to 565 form a group while observations 1 to 100 differ from the rest and that these 2 phenomena were explained by the experts
clusters <- rep("Group1", nrow(data_philips))
clusters[1:100] <- "Group2"
clusters[491:565] <- "Group3"
clusters <- factor(clusters)
X <- data.frame(cluster = clusters, data_philips)

# ICS with MCD_0.5-COV
out <- ICS(X[,-1], S1 = ICS_mcd_raw, S2 = ICS_cov,
           S1_args = list(alpha = 0.50, nsamp = 500))

# create plot
plot_philips <- component_plot(out, select = c(1:2, 8:9), clusters = X$cluster) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors)

# file name for plot
file_plot = "figures/empirical_applications/philips_MCD50-COV.%s"
# save plot to pdf
pdf(sprintf(file_plot, "pdf"), width = 8, height = 5)
print(plot_philips)
dev.off()
# save plot to png
png(sprintf(file_plot, "png"), width = 8, height = 5,
    unit = "in", res = 250)
print(plot_philips)
dev.off()
