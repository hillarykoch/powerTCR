\name{discgammagpd}
\alias{ddiscgammagpd}
\alias{pdiscgammagpd}
\alias{qdiscgammagpd}
\alias{rdiscgammagpd}


\title{The discrete gamma-GPD spliced threshold distribution}
\description{
Density, distribution function, quantile function and random generation for the discrete gamma-GPD spliced threshold distribution. The distribution has gamma bulk with shape equal to \code{shape} and rate equal to \code{rate}. It is spliced at a threshold equal to \code{u} and has a GPD tail with sigma equal to \code{sigma} and xi equal to \code{xi}. The proportion of data above the threshold phi is equal to \code{phiu} and the data are shifted according to \code{shift}.
}
\usage{
ddiscgammagpd(x, shape, rate, u, sigma, xi, phiu = NULL, shift = 0, log = FALSE)
pdiscgammagpd(q, shape, rate, u, sigma, xi, phiu = NULL, shift = 0)
qdiscgammagpd(p, shape, rate, u, sigma, xi, phiu = NULL, shift = 0)
rdiscgammagpd(n, shape, rate, u, sigma, xi, phiu = NULL, shift = 0)
}
\arguments{
    \item{x, q}{vector of quantiles.}
    \item{p}{vector of probabilities.}
    \item{n}{number of observations.}
    \item{shape}{shape parameter alpha of the truncated gamma distribution.}
    \item{rate}{rate parameter beta of the gamma distribution.}
    \item{u}{threshold.}
    \item{sigma}{scale parameter sigma of the GPD.}
    \item{xi}{shape parameter xi of the GPD}
    \item{phiu}{Propotion of data greater than or equal to threshold u.}
    \item{shift}{value the complete distribution is shifted by. Ideally, this is the smallest value of the count data from one sample.}
    \item{log}{Logical; if TRUE, probabilities p are given as log(p).}
}
\details{
The shape, rate, u, sigma, and xi parameters must be specified by the user. If \code{phiu} is left unspecified, it defaults to 1 minus the distribution function of a discrete gamma distribution (not a discrete truncated gamma) evaluated at \code{u-1}.
}

\value{
\code{ddiscgammagpd} gives the density, \code{pdiscgammagpd} gives the distribution function, \code{qdiscgammagpd} gives the quantile function, and \code{rdiscgammagpd} generates random variables from the described distribution.
}
\author{\email{hbk5086@psu.edu}}
\examples{
# Generate and sort a random sample for a log-log plot
d <- rdiscgammagpd(100, shape = 5, rate = .25, u = 25, sigma = 15, xi = .5, shift = 1)
d <- sort(d, decreasing = TRUE)
plot(log(d), log(1:100))

# When phiu is specified to .2, exactly 80\% of the data are below the threshold u
pdiscgammagpd(24, shape = 5, rate = .25, u = 25, sigma = 15, xi = .5, phiu = .2, shift = 1)

# Plot simulated data versus theoretical quantiles
quantiles <- qdiscgammagpd((100:1)/101, shape = 5, rate = .25, u = 25, sigma = 15, xi = .5, shift = 1)
plot(d, quantiles)

# Density below shift value is 0
ddiscgammagpd(0, shape = 5, rate = .25, u = 25, sigma = 15, xi = .5, shift = 1)
}
