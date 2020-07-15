\name{computeClusterMembProb}
\alias{computeClusterMembProb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Cluster membership probabilities in MR-PATH
}
\description{
  Computes SNP-specific cluster membership probabilities from the MR-PATH model.
}
\usage{
    computeClusterMembProb(data, Nsamples = 50000,
        impt_samples = NULL, MCEM_fit = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{data}{A data frame. (see Details in \code{\link{MR_PATH}}).}
    \item{Nsamples}{Number of desired samples.}
    \item{impt_samples}{Importance samples from \code{\link{getImportanceSamples}}. If \code{NULL}, obtain importance samples within function.}
    \item{MCEM_fit}{MC-EM fit from \code{MR_PATH}. If \code{impt_samples = NULL}, use it for obtaining importance samples. Default is \code{NULL}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
If \code{data} contains a column named \code{SNP}, the rows of the output matrix will be labeled by SNP name.
}
\value{
Returns a \code{p} by \code{K} matrix containing cluster membership probabilities.
}
\references{
%% ~put references to the literature/web site here ~
\url{https://arxiv.org/abs/2007.06476}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{MR_PATH}}, \code{\link{sampleBetas}}
}
\examples{
### Load data
data(hdl_chd)

### Filter weak instruments
hdl_chd = hdl_chd[hdl_chd$pval.selection < 5e-8,]

### Set your own K.
# For data-driven model selection, use MRPATH_selectModel
K = 2

### Set your own initial values.
# For automatic initial value selection, use MRPATH_optimizeInitVals
initVals = list("m_X" = mean(hdl_chd$beta.exposure),
                "lambdaX" = sd(hdl_chd$beta.exposure),
                "pis" = rep(1/K, K),
                "mus" = c(-.9, .4),
                "sds" = c(.9, .4))

### Run MC-EM algorithm
MCEM_fit = MR_PATH(K, hdl_chd, initVals)

### Compute SNP-specific cluster membership probabilities 
clustermemb_prob = computeClusterMembProb(hdl_chd, MCEM_fit = MCEM_fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
