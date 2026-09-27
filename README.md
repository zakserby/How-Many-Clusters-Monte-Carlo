# How Many Clusters Monte Carlo

R script for AEM 6850 (Empirical Methods, Cornell, Fall 2025) that uses a Monte Carlo simulation to ask how many clusters clustered standard errors need.

script_files:
- cluster_se_sim.R : simulates data with a regressor and error that are correlated within clusters, for 100 to 10,000 observations and 2 to 100 clusters, and plots the average clustered standard error divided by the true standard deviation of the estimate
- How-Many-Clusters-Monte-Carlo.Rproj : RStudio project, open this first so the script finds output_figure/

output_figure:
- How_Many_Clusters.png : ratio against number of clusters for each sample size

Other files:
- readme.rtf : full readme with general, methodological and data-specific information
