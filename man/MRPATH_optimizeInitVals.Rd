\name{MRPATH_optimizeInitVals}
\alias{MRPATH_optimizeInitVals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  MR-PATH Initial Value Optimizer
}
\description{
  Runs \code{\link{MR_PATH}} multiple times and picks fit with highest complete-data log-likelihood.
}
\usage{
    MRPATH_optimizeInitVals(K, data, Nreps = 10,
        verbose = FALSE, altModel = FALSE, init_seed = 8686, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{K}{Number of desired clusters (see Details).}
    \item{data}{A data frame. (see Details in \code{\link{MR_PATH}}).}
    \item{Nreps}{Number of repetitions to run \code{MR_PATH}.}
    \item{verbose}{If \code{TRUE}, returns complete-data log-likelihood at each repetition.}
    \item{altModel}{If \code{TRUE}, fits alternative model \emph{(This is still work in progress)}. Default is \code{FALSE}.}
    \item{init_seed}{Initial seed (for reproducibility).}
    \item{...}{Additional parameters to be passed to \code{MR_PATH}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
\code{K} can be chosen by the user or by the modified BIC criterion with \code{\link{MRPATH_selectModel}}.
}
\value{
Returns a list with
\item{fit}{MC-EM fit with optimal initial values.}
\item{initVals}{Optimal initial values.}

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

Nreps = 10 # Testing
# Nreps = 30 # Takes longer, but more stable results
verbose = TRUE # Print output
initValsObj = MRPATH_optimizeInitVals(K = K, data = hdl_chd, Nreps = Nreps, verbose=verbose)

# See optimized initial values
initValsObj$initVals

# See MC-EM fit with optimized initial values
initValsObj$fit
}
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{MR_PATH}}, \code{\link{MRPATH_selectModel}}
}

\references{
%% ~put references to the literature/web site here ~
\insertAllCited{}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
