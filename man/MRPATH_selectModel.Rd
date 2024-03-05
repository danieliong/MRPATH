\name{MRPATH_selectModel}
\alias{MRPATH_selectModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Model selection for MR-PATH model using a modified BIC criterion.
}
\description{
  Runs \code{\link{MR_PATH}} multiple times and picks fit with highest complete-data log-likelihood.
}
\usage{
    MRPATH_selectModel(data, K_range = 1:3, Nreps = 20,
        altModel = FALSE, verbose=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{data}{A data frame. (see Details in \code{\link{MR_PATH}}).}
    \item{K_range}{Range of \code{K} values to select from.}
    \item{Nreps}{Number of repetitions for \code{\link{MRPATH_optimizeInitVals}}.}
    \item{altModel}{If \code{TRUE}, fits alternative model \emph{(This is still work in progress)}. Default is \code{FALSE}.}
    \item{verbose}{If \code{TRUE}, prints BIC for each value in \code{K_range}.}
    \item{...}{Additional parameters to be passed to \code{\link{MRPATH_optimizeInitVals}}.}
}
% \details{
%%  ~~ If necessary, more details than the description above ~~
% }
\value{
Returns a list with
\item{bestK}{\code{K} with highest \code{BIC}.}
\item{bestFit}{MC-EM fit with highest \code{BIC}.}
\item{Q}{Vector of complete-data log-likelihoods for each \code{K} in \code{K_range}.}
\item{BIC}{Vector of BIC values for each \code{K} in \code{K_range}.}
}
\references{
%% ~put references to the literature/web site here ~
\url{https://arxiv.org/abs/2007.06476}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{MR_PATH}}, \code{\link{MRPATH_optimizeInitVals}}
}
\examples{
\donttest{
require(MRPATH)
### Load data
data(hdl_chd)

### Filter weak instruments
hdl_chd = hdl_chd[hdl_chd$pval.selection < 5e-8,]

Nreps = 10 # Testing
# Nrep = 30 # Takes longer, but more stable results
verbose = TRUE # Print output
K_range = 1:3 # Set range of K
modSelectionObj = MRPATH_selectModel(hdl_chd, K_range = K_range,
                                     Nreps = Nreps, verbose = verbose)

# See optimal K
modSelectionObj$bestK

# See optimal fit
modSelectionObj$bestFit

# See vector of complete-data log-likelihood values
modSelectionObj$Q

# See vector of BIC values
modSelectionObj$BIC
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
