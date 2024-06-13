\name{W.HJ.Est}
\alias{W.HJ.Est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weighted Hajek RMST Estimator
}
\description{
Estimate the weight-adjusted Hajek (HJ) RMSTs, RMST difference between two groups and their variances.
}
\usage{
W.HJ.Est(df, p, tau)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{df}{a data.frame including "A": treatment with two groups A=1 and A=2; selected variables.}
\item{p}{a vector of weights (e.g., calibration weights, inverse propensity score weights).}
\item{tau}{a scalar of truncation time for RMST.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with
\item{mu}{weighted HJ RMSTs for A=1, A=2, and their difference.}
\item{sd}{standard deviations of mu.}
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\examples{

}

