# Tandem clustering with invariant coordinate selection

## About tandem clustering with ICS
Tandem clustering is a well-known technique aiming to improve cluster
identification through initial dimension reduction. However, the usual approach 
using principal component analysis (PCA) has been criticized for focusing 
solely on inertia so that the first components do not necessarily retain the 
structure of interest for clustering. To address this limitation, we propose a 
new tandem clustering approach based on invariant coordinate selection (ICS).

More information can be found in our article:

A. Alfons, A. Archimbaud, K. Nordhausen and A. Ruiz-Gazen (2024).
Tandem clustering with invariant coordinate selection. 
arXiv preprint [arXiv:2212.06108](https://arxiv.org/abs/2212.06108).

## Reproduce results
This repository provides a collection of [R](https://CRAN.R-project.org/) 
scripts to reproduce all examples, simulations and figures in our article.

The easiest way to reproduce the results is to clone this repository with 
[RStudio](https://rstudio.com/products/rstudio/download/).  Running the 
scripts within the resulting RStudio project ensures that there are no issues 
with file paths for storing or reading results, or for producing files 
containing plots.  In addition, the RStudio project uses 
[renv](https://rstudio.github.io/renv/) to make sure that the correct 
versions of all required R packages are used.  After opening the RStudio 
project for the fist time, please type `renv::restore()` on the R command 
line to retrieve the correct versions of all required packages.

Please note that this repository is rather large because it also contains R 
data files with all simulation results.  This way, if you only want to quickly
reproduce the figures with simulation results, you do not actually need to run 
the simulations first.
