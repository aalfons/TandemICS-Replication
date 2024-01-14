
# load packages
library("dplyr")
library("ggplot2")
library("ggthemes")
library("magrittr")



# Utility functions ------------------------------------------------------------
# function for parsing axis labels
parse_labels <- function(labels, ...) parse(text = labels, ...)


# Colors and labels for plots --------------------------------------------------

# colors to be used in plots
colors_ICS <- ggthemes::colorblind_pal()(4)[c(3, 1, 2, 4)]
names(colors_ICS) <- c("Oracle", "Var", "Med", "Normal")
colors_PCA <- scales::hue_pal()(4)[c(4, 3)]
names(colors_PCA) <- c("80%", "k-1")
colors <- c(colors_ICS, colors_PCA)

# additional color to be used in plots
color_NA <- "black"

# labels to be used in plots
clusters_label <- "Cluster sizes"
criterion_label <- "Component selection"


# Import results for clustering ------------------------------------------------

# import results and renaming certain variables, dimension reduction methods,
# component selection criteria, etc.

# load results from file
p = 10
load( sprintf("results/simulations/results_clustering_p=%d_r=100.RData", p))

# rename some stuff and filter relevant results
res_df <- results %>%
  # rename some variables for consistency
  rename(outliers = epsilon) %>%
  # recode names of cluster settings, scatter matrices, and selection criteria
  mutate(outliers = ifelse(outliers == 0, "No outliers",
                           sprintf("%d%% outliers", 100 * outliers)),
         criterion = recode(criterion,
                            "discriminatory" = "Oracle",
                            "var" = "Var",
                            "med" = "Med",
                            "normal" = "Normal"),
         method = recode(method,
                         "pam" = "PAM",
                         "mfa" = "MFA")
  ) %>%
  # for clustering on observed data, keep only certain standardization
  filter(
    # for kmeans on observed data, keep only results on standardized data
    !(method == "kmeans" &
        scatter %in% c("Observed~data", "Observed~data~robstd")),
    # for tkmeans and PAM on observed data, keep only results on robustly
    # standardized data
    !(method %in% c("tkmeans", "PAM") &
        scatter %in% c("Observed~data", "Observed~data~std")),
    # for mclust and rmclust on observed data, keep only results without
    # standardization
    !(method %in% c("mclust", "rmclust") &
        scatter %in% c("Observed~data~std", "Observed~data~robstd")),
  ) %>%
  # clean up label for observed data and add variable indicating ICS or PCA
  mutate(scatter = recode(scatter,
                          "Observed~data~std" = "Observed~data",
                          "Observed~data~robstd" = "Observed~data"),
         type = ifelse(scatter %in% c("COV", "RMCD[0.75]"), "PCA",
                       ifelse(scatter %in%
                                c("Observed~data", "Observed~data~std",
                                  "Observed~data~robstd"), "",
                              "ICS")))


# ARI of best performing dimension reduction methods ---------------------------

# evaluation measure to be plotted
measure <- "ARI"
measure_label <- "Adjusted Rand index"

# which methods to select
keep_outliers <- c("No outliers", "2% outliers", "5% outliers")
keep_crit <- c("Med", "80%")
keep_scatter <- c("Observed~data", "LCOV-COV", "TCOV-COV", "TCOV-UCOV",
                  "RMCD[0.75]")
keep_method <- c("PAM", "kmeans", "tkmeans", "mclust", "rmclust")


## overall results across cluster settings

# select subset of results
res_df_selected <- res_df %>%
  filter(outliers %in% keep_outliers,
         is.na(criterion) | criterion %in% keep_crit,
         scatter %in% keep_scatter,
         method %in% keep_method)

# create plot
text_size_factor <- 8/6.5
plot_best <- res_df_selected %>%
  mutate(outliers = factor(outliers, levels = keep_outliers),
         criterion = factor(criterion, levels = keep_crit),
         scatter = factor(scatter, levels = keep_scatter),
         method = factor(method, levels = keep_method)) %>%
  ggplot(mapping = aes_string(x = "scatter", y = measure, color = "criterion",
                              fill = "criterion")) +
  geom_boxplot(alpha = 0.4) +
  scale_x_discrete(labels = parse_labels) +
  scale_color_manual(criterion_label, values = colors[keep_crit],
                     na.value = color_NA) +
  scale_fill_manual(criterion_label, values = colors[keep_crit],
                    na.value = color_NA) +
  theme_bw() +
  theme(axis.title = element_text(size = 11 * text_size_factor),
        axis.text = element_text(size = 9 * text_size_factor),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.title = element_text(size = 11 * text_size_factor),
        legend.text = element_text(size = 9 * text_size_factor),
        strip.text.x = element_text(size = 10 * text_size_factor),
        strip.text.y = element_text(size = 9 * text_size_factor),
        legend.position = "top") +
  labs(x = NULL, y = measure_label) +
  facet_grid(outliers ~ method)

# file name for plot
file_plot = "figures/simulations/clustering/%s_best.%s"
# save plot to pdf
pdf(sprintf(file_plot, measure, "pdf"), width = 8, height = 4.8)
print(plot_best)
dev.off()
# save plot to png
png(sprintf(file_plot, measure, "png"), width = 8, height = 4.8,
    unit = "in", res = 250)
print(plot_best)
dev.off()


## results for selected cluster settings

# to make sure cluster settings are in a certain order
cluster_levels <- c("50-50", "70-30", "80-20", "90-10", "95-5",
                    "33-33-33", "20-50-30", "10-80-10",
                    "20-20-20-20-20", "10-10-20-20-40")

# loop over contamination level and create plot
outlier_suffix <- c("no_outliers", "0.02_outliers", "0.05_outliers")
for (i in seq_along(keep_outliers)) {

  # select subset of cluster settings
  res_df_clusters <- res_df_selected %>%
    filter(outliers %in% keep_outliers[i])

  # create plot
  text_size_factor <- 8/6.5
  current_plot <- res_df_clusters %>%
    mutate(clusters = factor(clusters, levels = cluster_levels),
           criterion = factor(criterion, levels = keep_crit),
           scatter = factor(scatter, levels = keep_scatter),
           method = factor(method, levels = keep_method)) %>%
    ggplot(mapping = aes_string(x = "clusters", y = measure,
                                color = "criterion", fill = "criterion")) +
    geom_boxplot(alpha = 0.4) +
    scale_color_manual(criterion_label, values = colors[keep_crit],
                       na.value = color_NA) +
    scale_fill_manual(criterion_label, values = colors[keep_crit],
                      na.value = color_NA) +
    theme_bw() +
    theme(axis.title = element_text(size = 11 * text_size_factor),
          axis.text = element_text(size = 9 * text_size_factor),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "top",
          legend.title = element_text(size = 11 * text_size_factor),
          legend.text = element_text(size = 9 * text_size_factor),
          strip.text = element_text(size = 10 * text_size_factor)) +
    labs(x = clusters_label, y = measure_label) +
    facet_grid(method ~ scatter, labeller = label_parsed)

  # file name for plot
  file_plot = "figures/simulations/clustering/%s_best_clusters_%s.%s"
  # save plot to pdf
  pdf(sprintf(file_plot, measure, outlier_suffix[i], "pdf"),
      width = 8, height = 7)
  print(current_plot)
  dev.off()
  # save plot to png
  png(sprintf(file_plot, measure, outlier_suffix[i], "png"),
      width = 8, height = 7, unit = "in", res = 250)
  print(current_plot)
  dev.off()

}


# ARI and computation time for comparison with integrated GMM approaches -------

# evaluation measures to be plotted
measures <- c("ARI", "time")
measure_labels <- c("Adjusted Rand index", "Computation time (seconds)")

# which methods to select
keep_outliers <- c("No outliers", "2% outliers", "5% outliers")
keep_crit <- "Med"
keep_scatter <- "TCOV-COV"
keep_method_ICS <- c("mclust", "rmclust")
keep_method_GMM <- c("MFA", "clustvarsel")
keep_method <- c(keep_method_ICS, keep_method_GMM)

# different labels for methods
method_labels <- c(mclust = "ICS + mclust", rmclust = "ICS + rmclust",
                   MFA = "MFA", clustvarsel = "clustvarsel")

# select subset of results
res_df_ICS <- res_df %>%
  filter(outliers %in% keep_outliers,
         is.na(criterion) | criterion %in% keep_crit,
         scatter %in% keep_scatter,
         method %in% keep_method_ICS)
res_df_GMM <- res_df %>%
  filter(outliers %in% keep_outliers,
         method %in% keep_method_GMM)
res_df_selected <- bind_rows(res_df_ICS, res_df_GMM) %>%
  select(Run:method, type, all_of(measures)) %>%
  tidyr::pivot_longer(cols = all_of(measures),
                      names_to = "measure",
                      values_to = "value")

# create plot
text_size_factor <- 8/6.5
plot_GMM <-
  res_df_selected %>%
  mutate(outliers = factor(outliers, levels = keep_outliers),
         criterion = factor(criterion, levels = keep_crit),
         scatter = factor(scatter, levels = keep_scatter),
         method = factor(method, levels = rev(keep_method),
                         labels = rev(method_labels)),
         measure = factor(measure, levels = measures,
                          labels = measure_labels)) %>%
  ggplot(mapping = aes_string(x = "method", y = "value", color = "criterion",
                              fill = "criterion")) +
  geom_boxplot(alpha = 0.4, show.legend = FALSE) +
  scale_x_discrete(labels = parse_labels) +
  scale_color_manual(criterion_label, values = colors[keep_crit],
                     na.value = color_NA) +
  scale_fill_manual(criterion_label, values = colors[keep_crit],
                    na.value = color_NA) +
  theme_bw() +
  theme(axis.title = element_text(size = 11 * text_size_factor),
        axis.text = element_text(size = 9 * text_size_factor),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.title = element_text(size = 11 * text_size_factor),
        legend.text = element_text(size = 9 * text_size_factor),
        strip.text.x = element_text(size = 10 * text_size_factor),
        strip.text.y = element_text(size = 9 * text_size_factor),
        legend.position = "top") +
  labs(x = NULL, y = measure_label) +
  facet_grid(outliers ~ measure, scales = "free_x") +
  coord_flip()

# file name for plot
file_plot = "figures/simulations/clustering/%s_GMM.%s"
# save plot to pdf
pdf(sprintf(file_plot, paste(measures, collapse = "_"), "pdf"),
    width = 8, height = 4.8)
print(plot_GMM)
dev.off()
# save plot to png
png(sprintf(file_plot, paste(measures, collapse = "_"), "png"),
    width = 8, height = 4.8, unit = "in", res = 250)
print(plot_GMM)
dev.off()
