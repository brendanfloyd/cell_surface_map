# cell_surface_map
This repo contains all the code for analysis of the cell surface proximity labeling data from the Bertozzi lab published in Floyd et al. Mapping the nanoscale organization of the human cell surface proteome reveals new functional associations and surface antigen clusters. BioRxiv (2025). Code contains a focus on critical analyses included the unsupervised hierarchical clustering for functional domain analysis, the random forest classifier for predicting cell surface protein associations, and the strategy for making a 3D network heatmap. Link to the current version of the manuscript can be found here: pending link

# Directory guide
# Preprocessing
Preprocessing contains scripts to take you from the DIA-NN output (used to process all DIA data in this work) to a functioning Perseus software input. An additional script "cleans" the Perseus output and adds in imputated intensity values. The final script uses multiple results files to make a fold change matrix that can be used for all subsequent analyses.

# Prey-centric analysis scripts
These scripts use the fold change matrix to estimate a similarity matrix, compare subunit distances to a known protein complex (a T-cell receptor example is provided), and hierarchically cluster the data with cutoffs at each hierarchy.

# Random forest surface association prediction
These scripts will first generate a combined all-by-all fold change delta matrix and similarity matrix. This will then be used to predict cell surface protein associations using the random forest classifier script. A compiled gold-standard PPI dataset is included in this directory.

# 3D network heatmap
This script will generate the 3D network heatmaps shown in Figure 6 of this work. In addition a graph connectivity score will be estimated for the sampled protein.

We hope that you find these scripts and files useful and do not hesitate to reach out to us with any questions or comments!
