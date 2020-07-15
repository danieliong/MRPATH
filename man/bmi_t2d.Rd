\docType{data}
\name{bmi_t2d}
\alias{bmi_t2d}
\title{Effect of Body Mass Index (BMI) on Type-2 Diabetes (T2D)}
\format{A \code{data.frame} with 60 rows and 6 variables.}
\usage{
data(bmi_t2d)
}
\description{
This dataset is created from three genome-wide association studies using the three-sample summary-data MR design \insertCite{zhao2019powerful}{MRPATH}:
\enumerate{
\item \strong{Selection}: \insertCite{akiyama2017genome;textual}{MRPATH}
\item \strong{Exposure}: \insertCite{locke2015genetic;textual}{MRPATH}
\item \strong{Outcome}: \insertCite{mahajan2018fine;textual}{MRPATH}
}
The 60 SNPs selected are independent (distance \eqn{\ge 10} mega base pairs, \eqn{R^2 \le 0.001} in a reference panel) and are associated with T2D (p-value less than \eqn{5*10^{-8}}).
}
\references{
    \insertAllCited{}
}
\keyword{datasets}
