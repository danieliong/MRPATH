# MRPATH

A preprint for the MR-PATH paper is available at https://arxiv.org/abs/2007.06476.


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

## Examples

You can find examples using the HDl-CHD dataset described in the [preprint](https://arxiv.org/abs/2007.06476) by running any of the following R commands.

```
library(MRPATH)

example("MR_PATH")

example("MRPATH_optimizeInitVals")

example("MRPATH_selectModel")

example("getImportanceSamples")

example("sampleBetas")

example("computeClusterMembProb")

example("MRPATH_scatterplot")

example("MRPATH_barplot")
```

## Troubleshooting

Please report issues and suggestions for the **MRPATH** R package using the [Github issues tracker](https://github.com/danieliong/mr.path/issues).
