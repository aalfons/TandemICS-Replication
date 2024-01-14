
# Packages ----------------------------------------------------------------

library("clustrd")
library("parallel")
library("ICS")
library("rrcov")
library("ICSClust")
library("cellWise")  # philips data
library("EMMIXmfa")
library("clustvarsel")


# Parameters --------------------------------------------------------------
# control parameters for data generation
seed <- 20230508  # seed of the random number generator

# control parameters for ICS scatters
ICS_scatters_list <- list(
  `COV-COV[4]` = list(S1 = ICS_cov, S2 = ICS_cov4),
  `MLC-COV` = list(S1 = ICS_mlc, S2 = ICS_cov),
  `LCOV-COV` = list(S1 = ICS_lcov, S2 = ICS_cov,
                    S1_args = list(mscatter = "cov", proportion = 0.1)),
  `TCOV-COV` = list(S1 = ICS_tcov, S2 = ICS_cov,
                    S1_args = list(beta = 2)),
  `TCOV-UCOV` = list(S1 = ICS_tcov, S2 = ICS_ucov,
                     S1_args = list(beta = 2), S2_args = list(beta = 0.2)),
  `MCD[0.25]-MCD[0.95]` = list(S1 = ICS_mcd_raw, S2 = ICS_mcd_raw,
                               S1_args = list(alpha = 0.25,
                                              nsamp = 500),
                               S2_args = list(alpha = 0.95,
                                              nsamp = 500)),
  `MCD[0.10]-COV` = list(S1 = ICS_mcd_raw, S2 = ICS_cov,
                         S1_args = list(alpha = 0.10,
                                        nsamp = 500)),
  `MCD[0.20]-COV` = list(S1 = ICS_mcd_raw, S2 = ICS_cov,
                         S1_args = list(alpha = 0.20,
                                        nsamp = 500)),
  `MCD[0.25]-COV` = list(S1 = ICS_mcd_raw, S2 = ICS_cov,
                         S1_args = list(alpha = 0.25,
                                        nsamp = 500)),
  `MCD[0.50]-COV` = list(S1 = ICS_mcd_raw, S2 = ICS_cov,
                         S1_args = list(alpha = 0.50,
                                        nsamp = 500)),
  `MCD[0.75]-COV` = list(S1 = ICS_mcd_raw, S2 = ICS_cov,
                         S1_args = list(alpha = 0.75,
                                        nsamp = 500)),
  `RMCD[0.25]-RMCD[0.95]` = list(S1 = ICS_mcd_rwt, S2 = ICS_mcd_rwt,
                                 S1_args = list(alpha = 0.25,
                                                nsamp = 500),
                                 S2_args = list(alpha = 0.95,
                                                nsamp = 500)),
  `RMCD[0.10]-COV` = list(S1 = ICS_mcd_rwt, S2 = ICS_cov,
                          S1_args = list(alpha = 0.10,
                                         nsamp = 500)),
  `RMCD[0.20]-COV` = list(S1 = ICS_mcd_rwt, S2 = ICS_cov,
                          S1_args = list(alpha = 0.20,
                                         nsamp = 500)),
  `RMCD[0.25]-COV` = list(S1 = ICS_mcd_rwt, S2 = ICS_cov,
                          S1_args = list(alpha = 0.25,
                                         nsamp = 500)),
  `RMCD[0.50]-COV` = list(S1 = ICS_mcd_rwt, S2 = ICS_cov,
                          S1_args = list(alpha = 0.50,
                                         nsamp = 500)),
  `RMCD[0.75]-COV` =  list(S1 = ICS_mcd_rwt, S2 = ICS_cov,
                           S1_args = list(alpha = 0.75,
                                          nsamp = 500))
)


# control parameters for ICS criteria
ICS_criteria <- c("normal_crit", "med_crit", "var_crit", "discriminatory_crit")
ICS_criteria_args <- list(
  normal_crit = list(level = 0.05,  test = "agostino.test"),
  med_crit = list(nb_select = c()),
  var_crit = list(nb_select = c()),
  discriminatory_crit = list(clusters = vector(), nb_select = c())
)


# control parameters for PCA criteria
PCA_criteria <- c("pct_crit", "min_crit")
PCA_criteria_args <- list(
  `80%` = c(pct = 0.8),
  `k-1` = c(nb_select = c())
)

min_crit <- function(object, nb_select = NULL) {
  colnames(object$scores)[0:nb_select]
}

pct_crit <- function(object, pct = 0.8) {
  # select the names of the components which explains xxx pct of the total inertia
  nc <- which(cumsum(object$eigenvalues) > pct*sum(object$eigenvalues))[1]
  colnames(object$scores)[0:nc]
}

# control parameters for clustering methods
clustering <- c("kmeans_clust", "tkmeans_clust", "pam_clust", "rimle_clust",
                "mclust_clust", "rmclust_clust"
)
clustering_all <- c("kmeans_clust", "tkmeans_clust", "pam_clust", "rimle_clust",
                    "mclust_clust", "rmclust_clust", "reduced_kmeans_clust",
                    "factorial_kmeans_clust", "mfa_clust", "clustvarsel_clust"
)

reduced_kmeans_clust <- function (X, k, clusters_only = FALSE, ...){
  clust <- clustrd::cluspca(data = X, nclus = k, ndim = k-1,
                            method = "RKM", ...)
  out <- clust$cluster
  if (!clusters_only)
    out <- append(list(clust_method = "reduced_kmeans", clusters = out),
                  clust)
  out
}

factorial_kmeans_clust <- function (X, k, clusters_only = FALSE, ...){
  clust <- clustrd::cluspca(data = X, nclus = k, ndim = k-1,
                            method = "FKM", ...)
  out <- clust$cluster
  if (!clusters_only)
    out <- append(list(clust_method = "factorial_kmeans", clusters = out),
                  clust)
  out
}

clustvarsel_clust <- function (X, k, clusters_only = FALSE, ...){
  clust <- clustvarsel::clustvarsel(data = X,
                                    G = k,
                                    emModels1 = "V",
                                    emModels2 = "VVV",
                                    verbose = FALSE, ...)

  out <- clust$model$classification
  if (!clusters_only)
    out <- append(list(clust_method = "clustvarsel", clusters = out),
                  clust)
  out
}

mfa_clust <- function (X, k, clusters_only = FALSE, ...){
  clust <- EMMIXmfa::mfa(Y = X, g = k, q = k-1,
                         nkmeans = 20, nrandom = 10,
                         sigma_type = 'unique', D_type = 'unique', ...)
  out <- clust$clus
  if (!clusters_only)
    out <- append(list(clust_method = "mfa", clusters = out),
                  clust)
  out
}

clustering_args <- list(
  kmeans_clust = c(iter.max = 100, nstart = 20),
  tkmeans_clust = c(alpha = 0.05),
  pam_clust = c(),
  rimle_clust = c(npr.max = 0.05),
  mclust_clust = c(),
  rmclust_clust = c(),
  reduced_kmeans_clust = c(),
  factorial_kmeans_clust = c(),
  mfa_clust = c(),
  clustvarsel_clust = c()
)


cat(paste(Sys.time(), ": starting ...\n"))
set.seed(seed)


# Data ----

# Load data
data("data_philips", package = "cellWise")

# Define clusters: based on paper: that observations 491 to 565 form a group while observations 1 to 100 differ from the rest and that these 2 phenomena were explained by the experts
clusters <- rep("Group1", nrow(data_philips))
clusters[1:100] <- "Group2"
clusters[491:565] <- "Group3"
clusters <- factor(clusters)
data <- data.frame(cluster = clusters, data_philips)
data_name <- "philips"
nb_clusters <- length(unique(data$cluster))
# Define nb_select by default equals to the number of clusters -1
nb_select <- nb_clusters-1
true_clusters <- data[,1]

# Update some parameters for criteria
PCA_criteria_args$`k-1` <- c(nb_select = nb_select)
ICS_criteria_args$med_crit <- c(nb_select = nb_select)
ICS_criteria_args$var_crit <- c(nb_select = nb_select)
ICS_criteria_args$discriminatory_crit$nb_select <- nb_select
ICS_criteria_args$discriminatory_crit$clusters <- true_clusters

# Details
n = nrow(data)
p = ncol(data)-1
info <- data.frame(name = data_name, n = n, p = p, q = nb_clusters)

# No Dimension reduction ----
# Additional info
criterion <- NA_character_
scatter <- "Observed~data"
selected <- paste(colnames(data[,-1]), collapse = ",")


results_ARI_no_reduction <- lapply(clustering_all, function(method) {
  time_reduction <- system.time({
    if (method == "kmeans_clust"){
      reduced_df <- scale(data[,-1], center = TRUE, scale = TRUE)
    }else if (method == "reduced_kmeans_clust"){

      reduced_df <- data[,-1]
    }else{ # robustly center
      reduced_df <- scale(data[,-1],
                          center = apply(data[,-1], 2, median),
                          scale = apply(data[,-1], 2, mad))
    }
  })[["elapsed"]]


  time_selection <- system.time({
    nb_select <- ncol(reduced_df)
  })[["elapsed"]]

  # Compute discriminatory power
  eta2 <- ICSClust:::eta2_power(reduced_df, clusters = true_clusters,
                                select = 1:nb_select)
  time_clustering <- system.time({
    ARI <- tryCatch({
      clusters <- do.call(method,
                          append(list(X = reduced_df, k = nb_clusters,
                                      clusters_only = TRUE),
                                 clustering_args[[method]]))
      mclust::adjustedRandIndex(true_clusters, clusters)
    }, error = function(e) 0, warning = function(w) 0)

  })[["elapsed"]]
  cbind(info, criterion = criterion, scatter = scatter,
        method = gsub("_clust", "", method), ARI = ARI, eta2 = eta2,
        nb_select = nb_select, selected = selected,
        time_reduction = time_reduction,
        time_selection = time_selection,
        time_clustering = time_clustering)
})

# combine results
df_ARI_no_reduction <- do.call(rbind, results_ARI_no_reduction)

# PCA ----
## classical ----
time_reduction <- system.time({
  PCA_out <- rrcov::PcaClassic(data[,-1], k = p, kmax = p, scale = TRUE)
})[["elapsed"]]
scatter <- "COV"
### criteria ----
results_ARI_PCA_crit <- lapply(1:length(PCA_criteria), function(i) {
  # Select the components
  criterion <- names(PCA_criteria_args)[i]
  time_selection <- system.time({
    select <- do.call(PCA_criteria[i],
                      append(list(object = PCA_out),
                             PCA_criteria_args[[i]]))
  })[["elapsed"]]
  nb_select <- length(select)


  # Compute discriminatory power
  eta2 <-  tryCatch({ICSClust:::eta2_power(PCA_out@scores,
                                           clusters = true_clusters, select = select)
  },error = function(e) 0, warning = function(w) 0)
  reduced_df <- PCA_out@scores[,select, drop = FALSE]
  selected <- paste(colnames(reduced_df), collapse = ",")

  # clustering ----
  results_ARI_PCA <- lapply(clustering, function(method) {
    time_clustering <- system.time({
      ARI <- tryCatch({
        clusters <- do.call(method,
                            append(list(X = reduced_df, k = nb_clusters,
                                        clusters_only = TRUE),
                                   clustering_args[[method]]))
        mclust::adjustedRandIndex(true_clusters, clusters)
      }, error = function(e) 0, warning = function(w) 0)
    })[["elapsed"]]
    cbind(info, criterion = criterion, scatter = scatter,
          method = gsub("_clust", "", method), ARI = ARI, eta2 = eta2,
          nb_select = nb_select, selected = selected,
          time_reduction = time_reduction,
          time_selection = time_selection,
          time_clustering = time_clustering)
  })
  do.call(rbind, results_ARI_PCA)
})

df_results_ARI_PCA <- do.call(rbind, results_ARI_PCA_crit)

## robust -----
time_reduction <- system.time({
  rob_PCA_out <- rrcov::PcaCov(data[,-1], k = p, kmax = p, scale = TRUE,
                               cov.control =
                                 rrcov::CovControlMcd(alpha = 0.75,
                                                      nsamp = "deterministic"))
})[["elapsed"]]
scatter <- "RMCD[0.75]"
### criteria ----
results_ARI_rob_PCA_crit <- lapply(1:length(PCA_criteria), function(i) {
  # Select the components
  criterion <- names(PCA_criteria_args)[i]
  time_selection <- system.time({
    select <- do.call(PCA_criteria[i],
                      append(list(object = rob_PCA_out),
                             PCA_criteria_args[[i]]))
  })[["elapsed"]]
  nb_select <- length(select)


  # Compute discriminatory power
  eta2 <-  tryCatch({ICSClust:::eta2_power(PCA_out@scores,
                                           clusters = true_clusters, select = select)
  },error = function(e) 0, warning = function(w) 0)
  reduced_df <- PCA_out@scores[, select, drop = FALSE]
  selected <- paste(colnames(reduced_df), collapse = ",")

  # clustering ----
  results_ARI_PCA <- lapply(clustering, function(method) {
    time_clustering <- system.time({
      ARI <- tryCatch({
        clusters <- do.call(method,
                            append(list(X = reduced_df, k = nb_clusters,
                                        clusters_only = TRUE),
                                   clustering_args[[method]]))
        mclust::adjustedRandIndex(true_clusters, clusters)
      }, error = function(e) 0, warning = function(w) 0)
    })[["elapsed"]]
    cbind(info, criterion = criterion, scatter = scatter,
          method = gsub("_clust", "", method), ARI = ARI, eta2 = eta2,
          nb_select = nb_select, selected = selected,
          time_reduction = time_reduction,
          time_selection = time_selection,
          time_clustering = time_clustering)
  })
  do.call(rbind, results_ARI_PCA)
})

df_results_ARI_rob_PCA <- do.call(rbind, results_ARI_rob_PCA_crit)

# ICS ----
## scatters ------
results_ARI_ICS_scatters <- lapply(1:length(ICS_scatters_list),
                                   function(i) {
                                     scatter = names(ICS_scatters_list)[i]
                                     time_reduction <- system.time({
                                       ICS_out <- tryCatch({
                                         do.call(ICS::ICS,
                                                 append(list(X = data[,-1]),
                                                        ICS_scatters_list[[i]]))
                                       },error = function(e) NULL)
                                     })[["elapsed"]]

                                     ## criteria ----
                                     results_ARI_ICS_crit <- lapply(ICS_criteria, function(criterion) {
                                       # Select the components
                                       time_selection <- system.time({
                                         select <-  tryCatch({do.call(criterion, append(list(object = ICS_out,
                                                                                             select_only = TRUE),
                                                                                        ICS_criteria_args[[criterion]]))
                                         },error = function(e) NULL, warning = function(w) NULL)
                                       })[["elapsed"]]
                                       nb_select <- length(select)

                                       # Compute discriminatory power
                                       eta2 <-  tryCatch({ICSClust:::eta2_power(ICS::components(ICS_out),
                                                                                clusters = true_clusters,
                                                                                select = select)
                                       },error = function(e) 0, warning = function(w) 0)

                                       reduced_df <- tryCatch({ICS::components(ICS_out, select = select)
                                       },error = function(e) data.frame(), warning = function(w) data.frame())
                                       selected <- paste(colnames(reduced_df), collapse = ",")

                                       # clustering ----
                                       results_ARI_ICS_clust <- lapply(clustering, function(method) {
                                         time_clustering <- system.time({
                                           if(ncol(reduced_df)>0){
                                             ARI <- tryCatch({
                                               clusters <- do.call(method,
                                                                   append(list(X = reduced_df, k = nb_clusters,
                                                                               clusters_only = TRUE),
                                                                          clustering_args[[method]]))
                                               mclust::adjustedRandIndex(true_clusters, clusters)
                                             }, error = function(e) 0, warning = function(w) 0)
                                           }else{
                                             ARI <- 0
                                           }

                                         })[["elapsed"]]

                                         cbind(info, criterion =  gsub("_crit", "", criterion),
                                               scatter = scatter, method = gsub("_clust", "", method),
                                               ARI = ARI, eta2 = eta2,
                                               nb_select = nb_select, selected = selected,
                                               time_reduction = time_reduction,
                                               time_selection = time_selection,
                                               time_clustering = time_clustering)
                                       })

                                       do.call(rbind, results_ARI_ICS_clust)
                                     })

                                     # combine results from current simulation run into data frame
                                     do.call(rbind, results_ARI_ICS_crit)
                                   })

df_ARI_ICS <- do.call(rbind, results_ARI_ICS_scatters)

# combine results into data frame
results <- rbind(df_ARI_no_reduction, df_results_ARI_PCA,  df_results_ARI_rob_PCA,
                 df_ARI_ICS)

# compute the global time
results$time <- rowSums(results[, c("time_reduction", "time_selection", "time_clustering")])

# save results to file
file_results <- "results/empirical_applications/results_%s.RData"
save(results, seed, file = sprintf(file_results, data_name))

# print message that simulation is done
cat(paste(Sys.time(), ": finished.\n"))
