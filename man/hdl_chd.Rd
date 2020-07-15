\docType{data}
\name{hdl_chd}
\alias{hdl_chd}
\title{Effect of HDL Cholesterol (HDL-C) on Coronary Heart Disease (CHD)}
\format{A \code{data.frame} with 151 rows and 6 variables.}
\usage{
data(hdl_chd)
}
\description{
This dataset is created from three genome-wide association studies using the three-sample summary-data MR design \insertCite{zhao2019powerful}{MRPATH}:
\enumerate{
\item \strong{Selection}: GWAS of HDL-C by \insertCite{teslovich2010biological;textual}{MRPATH}
\item \strong{Exposure}: GWAS of lipoprotein subfractions by \insertCite{kettunen2016genome;textual}{MRPATH}
\item \strong{Outcome}: The CARDIoGRAMplusC4D with 1000 Genome Project imputation GWAS of CAD \insertCite{nikpay2015comprehensive}{MRPATH}.
}
The 151 SNPs selected are independent (distance \eqn{\ge 10} mega base pairs, \eqn{R^2 \le 0.001} in a reference panel) and are associated with at least one plasma lipid trait (the minimum p-value with HDl-C, LDL-C, and triglycerides is less than \eqn{10^{-4}}). 
}
\references{
    \insertAllCited{}
}
\keyword{datasets}
