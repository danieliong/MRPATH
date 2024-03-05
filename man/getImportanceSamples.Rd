\name{getImportanceSamples}
\alias{getImportanceSamples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Importance Sampling of SNP-specific latent variables in MR-PATH
}
\description{
  Obtain importance samples of SNP-specific latent variables in the MR-PATH model
}
\usage{
getImportanceSamples(data, MCEM_fit, Nsamples = 50000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Nsamples}{Number of desired samples}
  \item{data}{A data frame. (see Details in \code{\link{MR_PATH}}).}
  \item{MCEM_fit}{Output from \code{\link{MR_PATH}}.}
}
% \details{
%%  ~~ If necessary, more details than the description above ~~
% }
\value{
  \item{thetaX}{\code{p} by \code{Nsamples} matrix of importance samples of true marginal SNP-exposure effects}
  \item{beta}{\code{p} by \code{Nsamples} matrix of importance samples of SNP-specific causal effects}
  \item{alpha}{\code{p} by \code{Nsamples} by \code{K} array of importance samples of cluster membership probabilities.}
  \item{W}{\code{p} by \code{Nsamples} matrix of importance weights}
 }
% \references{
%% ~put references to the literature/web site here ~
% }
% \author{
%%  ~~who you are~~
% }
% \note{
%%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{sampleBetas}}, \code{\link{computeClusterMembProb}}
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

### Obtain importance samples
impt_samples = getImportanceSamples(hdl_chd, MCEM_fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
