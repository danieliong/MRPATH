\name{MR_PATH}
\alias{MR_PATH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  MR-PATH Model Fitting
}
\description{
  Fits MR-PATH model with MC-EM algorithm.
}
\usage{
MR_PATH(K, data, initVals, overDispersedY = FALSE,
    equalSds = FALSE, computeSE = TRUE, Nstart_MC = 500L,
    M = 4L, max_Nsamples = 500000L, min_iters = 2L,
    max_iters = 100L, alpha = 0.05, gamma = 0.05, eps = 0.005,
    verbose = FALSE, saveTraj = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{K}{Number of desired clusters (see Details).}
  \item{data}{A data frame (see Details)}
  \item{initVals}{List of initial values (see Details).}
  \item{overDispersedY}{If \code{TRUE}, estimates multiplicative over-dispersion parameter for SNP-outcome effects. Default is \code{FALSE}. (This is still being tested.)}
  \item{equalSds}{If \code{TRUE}, assume mixture components have the same standard deviation. Default is \code{FALSE}.}
  \item{computeSE}{If \code{TRUE}, compute standard errors of parameter estimates. Default is \code{TRUE}.}
  \item{Nstart_MC}{Initial Monte-Carlo sample size for E-step MC approximation. Default is \code{500}.}
  \item{M}{Geometric rate at which MC sample size is increased. Default is \code{4}.}
  \item{max_Nsamples}{Max MC sample size for E-step MC approximation. Default is \code{500000}.}
  \item{min_iters}{Min. number of iterations. Default is \code{2}.}
  \item{max_iters}{Max number of iterations. Default is \code{100}.}
  \item{alpha}{Threshold for Type 1 error rate in E-step approximation. Default is \code{0.05}.}
  \item{gamma}{Threshold for Type 1 error rate in convergence criterion. Default is \code{0.05}.}
  \item{eps}{Threshold for convergence criterion. Default is \code{0.005}.}
  \item{verbose}{If \code{TRUE}, prints output at each iteration. Default is \code{FALSE}.}
  \item{saveTraj}{If \code{TRUE}, save MC-EM trajectory for each parameter. Default is \code{FALSE}.}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
\code{K} can be chosen by the user or by the modified BIC criterion with \code{\link{MRPATH_selectModel}}.


The input data frame must contain the following variables:
\enumerate{
    \item beta.exposure
    \item beta.outcome
    \item se.exposure
    \item se.outcome
}

\code{initVals} must be a list with the following elements:
\enumerate{
    \item \code{m_X}
    \item \code{lambdaX}
    \item \code{pi}
    \item \code{mus}
    \item \code{sds}
}

Since the MC-EM algorithm reaches local optima, it is recommended to run \code{MR_PATH} with a couple different initial values and pick the fit with the highest \code{completeDataLogLik}. One way to pick initial values is to first visualize the data with \code{\link{MRPATH_scatterplot}}. For automatic initial value selection, see \code{\link{MRPATH_optimizeInitVals}}.


\code{MR_PATH} returns the following information about convergence of the MC-EM algorithm:
\enumerate{
    \item \code{Niters}: Number of iterations it takes to converge.
    \item \code{N_MC_end}: Number of MC samples at convergence.
    \item \code{completeDataLogLik}: complete-data log-likelihood value at convergence.
}
}
\value{
Returns a list containing
    \item{paramEst}{List of parameter estimates.}
    \item{standardErrors}{Vector of standard errors.}
    \item{convergenceInfo}{List of information about convergence of MC-EM (see Details).}
    \item{paramTraj}{If \code{saveTraj = TRUE}, returns matrix containing MC-EM trajectory for each parameter.}

    \insertNoCite{iong2020latent}{MRPATH}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\donttest{
require(MRPATH)
### Load HDL-CHD data
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
}
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{MRPATH_optimizeInitVals}}, \code{\link{MRPATH_selectModel}},
\code{\link{MRPATH_scatterplot}}
}
\references{
%% ~put references to the literature/web site here ~
    \insertAllCited{}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
