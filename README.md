# tb_nanostring_shiny_h.s

This app explores gene expression in 179 participants with active
tuberculosis (TB) or latent TB infection (LTBI), using 107 genes.

Users can choose PCA, UMAP, or NMF and select the number of k-means
clusters. The app displays a scatterplot colored by cluster and shaped
by TB status, cluster centers, a cluster-versus-status table, and the
adjusted Rand index (ARI).

## Files

- app.R — Shiny application
- TBnanostring.rds — course-provided dataset

## Install the required packages

Run the following in the R Console:

```r
install.packages(c(
  "shiny", "ggplot2", "uwot", "mclust", "BiocManager"
))

BiocManager::install("Biobase", update = FALSE, ask = FALSE)

install.packages("NMF")
```

## Run the app

Download and unzip the repository. Keep app.R and TBnanostring.rds
in the same folder.

Open app.R in RStudio and click Run App.

Alternatively, set the R working directory to that folder and run:

```r
shiny::runApp()
```

Startup may take a little time while UMAP and NMF are calculated.
Changing the cluster count reuses the calculated coordinates.

## Methods

PCA uses centered, unscaled expression values. UMAP uses 15 neighbors.
NMF uses rank 2 and 10 runs after shifting each gene to a minimum of zero.
K-means uses the two displayed dimensions and 25 random starts.
TB status is used for comparison after clustering.
