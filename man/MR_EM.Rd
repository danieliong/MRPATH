\name{MR_EM}
\alias{MR_EM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  MC-EM algorithm for Bayesian Mendelian Randomization
}
\description{
  MC-EM algorithm for Bayesian model that accounts for effect heterogeneity in two-sample summary-data Mendelian Randomization analysis
}
\usage{
MR_EM(K, initVals, Nstart_MC, X, Y, seX, seY, overDispersedY = FALSE, equalSds = TRUE, invGammaPrior = FALSE, prior_shape = 0.01, prior_scale = 0.01, M = 2, max_Nsamples = 1e6, max_iters = 100, alpha = 0.05, gamma = 0.05, eps = 0.01, verbose = FALSE, ncores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{K}{Number of mixture components}
  \item{initVals}{List of initial values with elements named m_X, lambdaX, pis, mus, sds, tau}
  \item{Nstart_MC}{Initial Monte-Carlo sample size for E-step MC approximation}
  \item{X}{Vector containing SNP-exposure effects}
  \item{Y}{Vector containing SNP-outcome effects}
  \item{seX}{Vector containing standard errors of SNP-exposure effects}
  \item{seY}{Vector containing standard errors of SNP-outcome effects}
  \item{overDispersedY}{If TRUE, estimates multiplicative over-dispersion parameter for SNP-outcome effects (Not recommended)}
  \item{equalSds}{If TRUE, assume mixture components have the same standard deviation}
  \item{invGammaPrior}{If TRUE, specify Inverse-Gamma prior on mixture component variances with shape parameter prior_shape and scale parameter prior_scale. (Not recommended)}
  \item{prior_shape}{Shape parameter for Inverse-Gamma prior}
  \item{prior_scale}{Scale parameter for Inverse-Gamma prior}
  \item{M}{Geometric rate at which MC sample size is increased}
  \item{max_Nsamples}{Max MC sample size for E-step MC approximation}
  \item{max_iters}{Max number of iterations}
  \item{alpha}{Threshold for Type 1 error rate in E-step approximation}
  \item{gamma}{Threshold for Type 1 error rate in convergence criterion}
  \item{eps}{Threshold for convergence criterion}
  \item{verbose}{If TRUE, prints output at each iteration}
  \item{ncores}{Number of cores used in E-step sampling. (Parallelization not currently implemented.)}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{m_X}{Empirical Bayes estimate of true instrument strength mean}
  \item{lambdaX}{Empirical Bayes estimate of true instrument strength SD}
  \item{pis}{Empirical Bayes estimate of mixture component proportions}
  \item{mus}{Empirical Bayes estimate of mixture component means}
  \item{sds}{Empirical Bayes estimate of mixture component SDs}
  \item{tau}{Empirical Bayes estimate of over-dispersion parameter}
  \item{m_X_updates}{Trace for m_X}
  \item{lambdaX_updates}{Trace for lambdaX}
  \item{pis_updates}{Trace for pis}
  \item{mus_updates}{Trace for mus}
  \item{sds_updates}{Trace for sds}
  \item{N_iters}{Number of iterations until convergence}
  \item{N_end_MC}{MC sample size at convergence}
  \item{equalSds}{}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
