\name{W.GF.Est}
\alias{W.GF.Est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weighted G-Formula RMST Estimator
}
\description{
Estimate the weight-adjusted G-Formula (GF) RMSTs, RMST difference between two groups and their variances.
}
\usage{
W.GF.Est(df, p, varX, gf, tau)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{df}{a data.frame including "A": treatment with two groups A=1 and A=2; selected variables.}
\item{p}{a vector of weights (e.g., calibration weights, inverse propensity score weights).}
\item{varX}{a vector including names of the calibrated covariates, should be a subset of the selected variables included in df.}
\item{gf}{a vector function of calibrated covariates (e.g., moment function).}
\item{tau}{a scalar of truncation time for RMST.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with
\item{mu}{weighted GF RMSTs for A=1, A=2, and their difference.}
\item{sd}{standard deviations of mu.}
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\examples{

}

