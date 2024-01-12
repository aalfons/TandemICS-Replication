
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

# labels to be used in plots
clusters_label <- "Cluster sizes"
criterion_label <- "Component selection"


# Import results for dimension reduction ---------------------------------------

# import results and renaming certain variables, dimension reduction methods,
# component selection criteria, etc.

# load file
p = 10
load(sprintf("results/simulations/results_reduc_dim_p=%d_r=100.RData", p))

# rename some stuff and filter relevant results
res_df <- results %>%
  # rename some variables for consistency
  rename(outliers = epsilon) %>%
  # recode names of cluster settings, scatter matrices, and selection criteria
  mutate(outliers = ifelse(outliers == 0, "No outliers",
                           sprintf("%d%% outliers", 100 * outliers)),
         criterion = recode(criterion,
                            "discriminatory_crit" = "Oracle",
                            "var_crit" = "Var",
                            "med_crit" = "Med",
                            "normal_crit" = "Normal")
  ) %>%
  # add variable indicating ICS or PCA
  mutate(type = ifelse(scatter %in% c("COV", "RMCD[0.75]"), "PCA",
                       ifelse(scatter %in%
                                c("Observed~data", "Observed~data~std",
                                  "Observed~data~robstd"), "",
                              "ICS"))) %>%
  mutate(criterion = ifelse(is.na(criterion), scatter, criterion))


# evaluation measure to be plotted
measure <- "eta2"
measure_label <- "Discriminatory power"


# Determine for which methods to take a closer look ----------------------------

# filter on component selection criteria
keep_outliers <- "No outliers"
keep_crit <- c("Oracle", "Med", "Var", "Normal", "80%", "k-1")
res_df_overview <- res_df %>%
  filter(outliers %in% keep_outliers,
         criterion %in% keep_crit)

# order methods by average variance based on different selection criteria
# (this ordering is not used directly, but serves as a basis for a manual
# ordering in the plot with an overview of all dimension reduction methods
# and all component selection criteria)
order_df <- res_df_overview %>%
  filter(criterion != "Oracle") %>%
  group_by(scatter, criterion) %>%
  summarize(variance = var(.data[[measure]], na.rm = TRUE),
            .groups = "drop") %>%
  group_by(scatter) %>%
  summarize(variance = min(variance),
            .groups = "drop") %>%
  arrange(variance)
order_df$scatter

# order of dimension reduction methods in plot
scatter_ordered <- c(
  # no dimension reduction
  "Observed~data",
  # best performing scatter pairs for ICS
  "LCOV-COV", "TCOV-COV", "TCOV-UCOV",
  # MCD-COV scatter pairs
  sprintf("MCD[0.%d]-COV", c(10, 20, 25, 50, 75)),
  sprintf("RMCD[0.%d]-COV", c(10, 20, 25, 50, 75)),
  # MCD-MCD scatter pairs
  "MCD[0.25]-MCD[0.95]", "RMCD[0.25]-RMCD[0.95]",
  # other scatter pairs for ICS
  "MLC-COV", "COV-COV[4]",
  # scatter matrices for PCA
  "RMCD[0.75]", "COV"
)

# create plot
text_size_factor <- 10/6.5
plot_all <- res_df_overview %>%
  mutate(scatter = factor(scatter, levels = scatter_ordered),
         criterion = factor(criterion, levels = keep_crit)) %>%
  ggplot(mapping = aes_string(x = "scatter", y = measure, color = "criterion",
                              fill = "criterion")) +
  geom_boxplot(alpha = 0.4, position = position_dodge2(reverse = TRUE)) +
  coord_flip() +
  scale_x_discrete(limits = rev, labels = parse_labels) +
  scale_color_manual(criterion_label,  values = colors[keep_crit]) +
  scale_fill_manual(criterion_label, values = colors[keep_crit]) +
  theme_bw() +
  theme(axis.title = element_text(size = 11 * text_size_factor),
        axis.text = element_text(size = 9 * text_size_factor),
        legend.title = element_text(size = 11 * text_size_factor),
        legend.text = element_text(size = 9 * text_size_factor),
        strip.text = element_text(size = 10 * text_size_factor)) +
  labs(x = NULL, y = measure_label) +
  facet_grid(type ~ ., scales = "free_y", space = "free_y")

# file name for plot
file_plot <- "figures/simulations/reduc_dim/%s_all.%s"
# save plot to pdf
pdf(sprintf(file_plot, measure, "pdf"), width = 10, height = 9.5)
print(plot_all)
dev.off()
# save plot to png
png(sprintf(file_plot, measure, "png"), width = 10, height = 9.5,
    unit = "in", res = 250)
print(plot_all)
dev.off()


# Discriminatory power of best performing dimension reduction methods ----------

# which methods to select
keep_outliers <- c("No outliers", "2% outliers", "5% outliers")
keep_crit <- c("Oracle", "Med", "80%")
keep_scatter <- c("LCOV-COV", "TCOV-COV", "TCOV-UCOV", "RMCD[0.75]")
keep_clusters <- c("50-50", "70-30", "80-20", "90-10", "95-5",
                   "33-33-33", "20-50-30", "10-80-10",
                   "20-20-20-20-20", "10-10-20-20-40")

## overall results across cluster settings

# select subset of methods
res_df_selected <- res_df %>%
  filter(outliers %in% keep_outliers,
         criterion %in% keep_crit,
         scatter %in% keep_scatter)

# create plot
plot_best <- res_df_selected %>%
  mutate(outliers = factor(outliers, levels = keep_outliers),
         criterion = factor(criterion, levels = keep_crit),
         scatter = factor(scatter, levels = keep_scatter)) %>%
  ggplot(mapping = aes_string(x = "scatter", y = measure, color = "criterion",
                              fill = "criterion")) +
  geom_boxplot(alpha = 0.4) +
  scale_x_discrete(labels = parse_labels) +
  scale_color_manual(criterion_label, values = colors[keep_crit]) +
  scale_fill_manual(criterion_label, values = colors[keep_crit]) +
  theme_bw() +
  theme(axis.title = element_text(size = 11),
        axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 9)) +
  labs(x = NULL, y = measure_label) +
  facet_grid(. ~ outliers)

# file name for plot
file_plot <- "figures/simulations/reduc_dim/%s_best.%s"
# save plot to pdf
pdf(sprintf(file_plot, measure, "pdf"), width = 6.5, height = 3)
print(plot_best)
dev.off()
# save plot to png
png(sprintf(file_plot, measure, "png"), width = 6.5, height = 3,
    unit = "in", res = 250)
print(plot_best)
dev.off()


## results for selected cluster settings

# select subset of cluster settings
res_df_clusters <- res_df_selected %>%
  filter(clusters %in% keep_clusters)

# create plot
text_size_factor <- 8/6.5
plot_best_clusters <- res_df_clusters %>%
  mutate(outliers = factor(outliers, levels = keep_outliers),
         clusters = factor(clusters, levels = keep_clusters),
         criterion = factor(criterion, levels = keep_crit),
         scatter = factor(scatter, levels = keep_scatter)) %>%
  ggplot(mapping = aes_string(x = "clusters", y = measure,
                              color = "criterion", fill = "criterion")) +
  geom_boxplot(alpha = 0.4) +
  scale_color_manual(criterion_label, values = colors[keep_crit]) +
  scale_fill_manual(criterion_label, values = colors[keep_crit]) +
  theme_bw() +
  theme(axis.title = element_text(size = 11 * text_size_factor),
        axis.text = element_text(size = 9 * text_size_factor),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "top",
        legend.title = element_text(size = 11 * text_size_factor),
        legend.text = element_text(size = 9 * text_size_factor),
        panel.spacing.y = unit(0.5, "lines"),
        strip.text = element_text(size = 10 * text_size_factor)) +
  labs(x = clusters_label, y = measure_label) +
  facet_grid(outliers ~ scatter, labeller = labeller(scatter = label_parsed))

# file name for plot
file_plot <- "figures/simulations/reduc_dim/%s_best_clusters.%s"
# save plot to pdf
pdf(sprintf(file_plot, measure, "pdf"), width = 8, height = 8)
print(plot_best_clusters)
dev.off()
# save plot to png
png(sprintf(file_plot, measure, "png"), width = 8, height = 8,
    unit = "in", res = 250)
print(plot_best_clusters)
dev.off()


# Discriminatory power of MCD-COV scatter pairs --------------------------------

# which methods to select
keep_outliers <- c("No outliers", "2% outliers", "5% outliers")
keep_crit <- c("Oracle", "Med")
keep_scatter <- sprintf("MCD[0.%d]-COV", c(10, 25, 50))
keep_clusters <- c("50-50", "55-45", "60-40", "70-30", "80-20", "90-10", "95-5",
                   "33-33-33", "10-80-10", "20-20-20-20-20", "10-10-20-20-40")

## overall results across cluster settings

# select subset of methods
res_df_selected <- res_df %>%
  filter(outliers %in% keep_outliers,
         criterion %in% keep_crit,
         scatter %in% keep_scatter)

# # create plot
# plot_MCD <- res_df_selected %>%
#   mutate(outliers = factor(outliers, levels = keep_outliers),
#          criterion = factor(criterion, levels = keep_crit),
#          scatter = factor(scatter, levels = keep_scatter)) %>%
#   ggplot(mapping = aes_string(x = "scatter", y = measure, color = "criterion",
#                               fill = "criterion")) +
#   # geom_boxplot(alpha = 0.4, position = position_dodge2(reverse = TRUE)) +
#   # coord_flip() +
#   # scale_x_discrete(limits = rev, labels = parse_labels) +
#   geom_boxplot(alpha = 0.4) +
#   scale_x_discrete(labels = parse_labels) +
#   scale_color_manual(criterion_label, values = colors[keep_crit]) +
#   scale_fill_manual(criterion_label, values = colors[keep_crit]) +
#   theme_bw() +
#   theme(axis.title = element_text(size = 11),
#         axis.text = element_text(size = 9),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#         legend.title = element_text(size = 11),
#         legend.text = element_text(size = 9)) +
#   labs(x = NULL, y = measure_label) +
#   facet_grid(. ~ outliers)
#
# # file name for plot
# file_plot <- "figures/simulations/reduc_dim/%s_MCD.%s"
# # save plot to pdf
# pdf(sprintf(file_plot, measure, "pdf"), width = 6.5, height = 3)
# print(plot_MCD)
# dev.off()
# # save plot to png
# png(sprintf(file_plot, measure, "png"), width = 6.5, height = 3,
#     unit = "in", res = 250)
# print(plot_MCD)
# dev.off()


## results for selected cluster settings

# select subset of cluster settings
res_df_clusters <- res_df_selected %>%
  filter(clusters %in% keep_clusters)

# create plot
text_size_factor <- 8/6.5
plot_MCD_clusters <- res_df_clusters %>%
  mutate(outliers = factor(outliers, levels = keep_outliers),
         clusters = factor(clusters, levels = keep_clusters),
         criterion = factor(criterion, levels = keep_crit),
         scatter = factor(scatter, levels = keep_scatter)) %>%
  ggplot(mapping = aes_string(x = "clusters", y = measure,
                              color = "criterion", fill = "criterion")) +
  geom_boxplot(alpha = 0.4) +
  scale_color_manual(criterion_label, values = colors[keep_crit]) +
  scale_fill_manual(criterion_label, values = colors[keep_crit]) +
  theme_bw() +
  theme(axis.title = element_text(size = 11 * text_size_factor),
        axis.text = element_text(size = 9 * text_size_factor),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "top",
        legend.title = element_text(size = 11 * text_size_factor),
        legend.text = element_text(size = 9 * text_size_factor),
        strip.text = element_text(size = 10 * text_size_factor)) +
  labs(x = clusters_label, y = measure_label) +
  facet_grid(outliers ~ scatter, labeller = labeller(scatter = label_parsed))

# file name for plot
file_plot <- "figures/simulations/reduc_dim/%s_MCD_clusters.%s"
# save plot to pdf
pdf(sprintf(file_plot, measure, "pdf"), width = 8, height = 8)
print(plot_MCD_clusters)
dev.off()
# save plot to png
png(sprintf(file_plot, measure, "png"), width = 8, height = 8,
    unit = "in", res = 250)
print(plot_MCD_clusters)
dev.off()


# Discriminatory power of COV-COV4  scatter pair -------------------------------

# which methods to select
keep_outliers <- c("No outliers", "2% outliers", "5% outliers")
keep_crit <- c("Oracle", "Med")
keep_scatter <- "COV-COV[4]"

# select subset of cluster settings
res_df_clusters <- res_df %>%
  filter(outliers %in% keep_outliers,
         criterion %in% keep_crit,
         scatter %in% keep_scatter)

# order cluster settings according to number of clusters and variation in
# cluster sizes
cluster_df <- unique(res_df_clusters[, c("q", "clusters")])
cluster_list <- split(cluster_df, cluster_df$q)
order_list <- lapply(cluster_list, function(df) {
  mat <- stringr::str_split(df$clusters, pattern = "-", simplify = TRUE)
  mode(mat) <- "numeric"  # keeps information on dimensions
  variation <- apply(mat, 1, function(x) 1 - sum((x/100)^2))
  cbind(df, variation)
})
order_df <- do.call(rbind, order_list)
cluster_order <- order(order_df$q, -order_df$variation)
cluster_levels <- order_df$clusters[cluster_order]

# create plot
text_size_factor <- 8/6.5
plot_COVCOV4_clusters <- res_df_clusters %>%
  mutate(outliers = factor(outliers, levels = keep_outliers),
         clusters = factor(clusters, levels = cluster_levels),
         criterion = factor(criterion, levels = keep_crit),
         scatter = factor(scatter, levels = keep_scatter)) %>%
  ggplot(mapping = aes_string(x = "clusters", y = measure,
                              color = "criterion", fill = "criterion")) +
  geom_boxplot(alpha = 0.4) +
  scale_color_manual(criterion_label, values = colors[keep_crit]) +
  scale_fill_manual(criterion_label, values = colors[keep_crit]) +
  theme_bw() +
  theme(axis.title = element_text(size = 11 * text_size_factor),
        axis.text = element_text(size = 9 * text_size_factor),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "top",
        legend.title = element_text(size = 11 * text_size_factor),
        legend.text = element_text(size = 9 * text_size_factor),
        strip.text = element_text(size = 10 * text_size_factor)) +
  labs(x = clusters_label, y = measure_label) +
  facet_grid(outliers ~ scatter, labeller = labeller(scatter = label_parsed))

# file name for plot
file_plot <- "figures/simulations/reduc_dim/%s_COV-COV4_clusters.%s"
# save plot to pdf
pdf(sprintf(file_plot, measure, "pdf"), width = 8, height = 8)
print(plot_COVCOV4_clusters)
dev.off()
# save plot to png
png(sprintf(file_plot, measure, "png"), width = 8, height = 8,
    unit = "in", res = 250)
print(plot_COVCOV4_clusters)
dev.off()


# check selected components
df_selection <- res_df_clusters %>%
  mutate(outliers = factor(outliers, levels = keep_outliers),
         clusters = factor(clusters, levels = cluster_levels),
         criterion = factor(criterion, levels = keep_crit),
         scatter = factor(scatter, levels = keep_scatter),
         selected = as.factor(selected)) %>%
  group_by(outliers, clusters, criterion, scatter, selected) %>%
  summarize(Frequency = n(),
            .groups = "drop") %>%
  tidyr::pivot_wider(names_from = criterion, values_from = Frequency)
