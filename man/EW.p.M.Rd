\name{EW.p.M}
\alias{EW.p.M}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calibration Weights
}
\description{
Calculate the calibration weights given the selected covariates and the sample moments in the target population.
}
\usage{
EW.p.M(df,varX,gf,M)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{df}{a data.frame including "A": treatment with values 1,2,...,K; selected variables.}
\item{varX}{a vector including names of the calibrated covariates, should be a subset of the selected variables included in df.}
\item{gf}{a vector function of calibrated covariates (e.g., moment function).}
\item{M}{sample moments according to gf in the target population.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a vector of the calibration weights, with order matching the individuals in df.
}
\examples{

}

