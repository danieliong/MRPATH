# MRPATH


## Setup

To install the package, run the following command in R.

```
library(devtools)
install_github("danieliong/MRPATH")
```
## How to use MRPATH

The main functions in **MRPATH** are `MR_PATH`, `MRPATH_optimizeInitVals`, `MRPATH_selectModel`. Each of these functions provide increasing levels of automation in the model-fitting process.

In `MR_PATH`, the user must input data, the desired number of clusters, and initial values. In `MRPATH_optimizeInitVals`, the user has to input data and the desired number of clusters. In `MRPATH_selectModel`, the user just has to input data.  
