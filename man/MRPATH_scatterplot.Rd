\name{MRPATH_scatterplot}
\alias{MRPATH_scatterplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  MR-PATH Scatter Plot
}
\description{
  Scatterplot with error bars and results from \code{\link{MR_PATH}}.
}
\usage{
    MRPATH_scatterplot(data, MCEM_fit = NULL,
    exposure_name = "exposure", outcome_name = "outcome", overDispersedY = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{data}{A data frame. (see Details in \code{\link{MR_PATH}}).}
    \item{MCEM_fit}{Output from \code{\link{MR_PATH}}. If \code{NULL}, just plot the data and error bars.}
    \item{exposure_name}{Name of risk exposure variable.}
    \item{outcome_name}{Name of disease outcome variable.}
    \item{overDispersedY}{If \code{TRUE}, replaces \code{se.exposure} with \code{tau * se.exposure}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Each point is colored according to the cluster with the highest cluster membership probability. Solid lines represent \code{mus} and shaded regions represent 68\% interval for each cluster.
}
\value{
Returns a \code{ggplot} scatterplot with error bars and results from \code{\link{MR_PATH}}.
}
\references{
%% ~put references to the literature/web site here ~
\url{https://arxiv.org/abs/2007.06476}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{MRPATH_barplot}}
}
\examples{
### Load data
data(hdl_chd)

### Filter weak instruments
hdl_chd = hdl_chd[hdl_chd$pval.selection < 5e-8,]

### Plot data without MC-EM fit results
MRPATH_scatterplot(hdl_chd)

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

### Plot scatterplot with MC-EM fit results
MRPATH_scatterplot(hdl_chd, MCEM_fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
