# MRPATH


## Setup

To install the package, run the following command in R.

```
library(devtools)
install_github("danieliong/MRPATH")
```
## How to use MRPATH

The main functions in **MRPATH** are `MR_PATH`, `MRPATH_optimizeInitVals`, and `MRPATH_selectModel`. Each of these functions provide increasing levels of automation in the model-fitting process. In `MR_PATH`, the user must input data, the desired number of clusters, and initial values. In `MRPATH_optimizeInitVals`, the user has to input data and the desired number of clusters. In `MRPATH_selectModel`, the user just has to input data.

Auxiliary functions available in **MRPATH**: `sampleBetas` and `computeClusterMembProb`. These functions are used after obtaining the MC-EM model fit to investigate SNP-specific latent variables in the MR-PATH model.

Plotting functions available in **MRPATH**: `MRPATH_scatterplot`, `MRPATH_barplot`. 
