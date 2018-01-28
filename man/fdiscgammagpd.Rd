\name{fdiscgammagpd}
\alias{fdiscgammagpd}
\title{Fit the discrete gamma-GPD spliced threshold model}
\usage{
fdiscgammagpd(x, useq, shift = NULL, pvector=NULL,
std.err = TRUE, method = "Nelder-Mead", ...)
}
\description{
This function takes count data and fits the gamma-GPD spliced threshold model to it. The model consists of a discrete truncated gamma as the bulk distribution, up to the threshold, and a discrete GPD at and above the threshold. The 'shift' is ideally the minimum count in the sample.
}
\arguments{
    \item{x}{A vector of count data.}
    \item{useq}{A vector of possible thresholds to search over. These should be discrete numbers.}
    \item{shift}{The amount the distribution is shifted. It is recommended to use the minimum number in the count data when modeling the clone size distribution of the TCR repertoire.}
    \item{pvector}{A vector of 5 elements corresponding to the initial parameter estimates. These 5 initial values are for the gamma shape and rate, the threshold, and the GPD sigma and xi. If they are not prespecified, the function computes pvector automatically.}
    \item{std.err}{Logical. Should the standard errors on the estimates be computed from the Hessian matrix?}
    \item{method}{Character string listing optimization method fed to \link[stats]{optim}. Defaults to Nelder-Mead.}
}
\value{
    \item{x}{Numerical vector of the original data input}
    \item{init}{Numerical vector of the initial values of the parameter estimates. This is the same as pvector.}
    \item{useq}{Numerical vector containing the thresholds the grid search was performed over.}
    \item{nllhuseq}{Numerical vector of negative log likelihoods computed at each threshold in useq.}
    \item{optim}{Output from optim for the bulk and tail distributions.}
    \item{nllh}{The negative log likelihood corresponding to the maximum likelihood fitted distribution.}
    \item{mle}{A numerical vector containing the estimates for phi, shape, rate, threshold, sigma, and xi.}
    \item{fisherInformation}{The Fisher information matrix computed from the Hessian output from optim.}
}
\author{\email{hbk5086@psu.edu}}
\examples{
data("repertoires")
thresholds1 <- unique(quantile(repertoires[[1]], c(.75,.8,.85,.9,.95)))
thresholds2 <- unique(quantile(repertoires[[2]], c(.75,.8,.85,.9,.95)))

fit1 <- fdiscgammagpd(repertoires[[1]], useq = thresholds1, shift = min(repertoires[[1]]))
fit2 <- fdiscgammagpd(repertoires[[2]], useq = thresholds1, shift = min(repertoires[[2]]))

fit1
fit2
}
