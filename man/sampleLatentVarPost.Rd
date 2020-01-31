\name{sampleLatentVarPost}
\alias{sampleLatentVarPost}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Importance Sampling of Latent Variable Posterior 
}
\description{
  Samples from the posterior distribution of the latent variables in Bayesian MR model
}
\usage{
sampleLatentVarPost(N_samples, X, Y, seX, seY, params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N_samples}{Number of desired samples}
  \item{X}{Vector containing SNP-exposure effects}
  \item{Y}{Vector containing SNP-outcome effects}
  \item{seX}{Vector containing standard errors of SNP-exposure effects}
  \item{seY}{Vector containing standard errors of SNP-outcome effects}
  \item{params}{List containing prior parameters. Recommended: MR_EM output list}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{muX_samps}{p by N_samples matrix of importance samples from posterior of true instrument strengths}
  \item{beta_samps}{p by N_samples matrix of importance samples from posterior of instrument-specific causal effects}
  \item{alpha_samps}{p by N_samples by K array of importance samples of posterior probabilities}
  \item{W}{Importance weights}
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
