\name{MRPATH_barplot}
\alias{MRPATH_barplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  MR-PATH Bar Plot
}
\description{
  Plots a barplot for credible intervals and cluster membership probabilities with SNPs ordered by median of SNP-specific causal effects.
}
\usage{
    MRPATH_barplot(data, MCEM_fit, ret.snps = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{data}{A data frame. (see Details in \code{\link{MR_PATH}}).}
    \item{MCEM_fit}{Output from \code{\link{MR_PATH}}.}
    \item{ret.snps}{If \code{TRUE}, returns both plot and a list of ordered SNPs. Default is \code{FALSE}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a 95\% credible interval plot (top) and cluster membership probabilities barplot (bottom).
}
\references{
%% ~put references to the literature/web site here ~
\url{https://arxiv.org/abs/2007.06476}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{MRPATH_scatterplot}}
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

### Plot barplots
MRPATH_barplot(hdl_chd, MCEM_fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
