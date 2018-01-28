\name{JS_spliced}
\alias{JS_spliced}
\title{
Compute the Jensen-Shannon distance between two fitted discrete gamma-GPD spliced threshold distributions
}
\description{
After models have been fit to your samples, the pairwise JS distance can be computed between them. This function takes the fitted model parameters from 2 distributions and computes the JS distance between them. When all pairwise distances have been computed, they can be used to do hierarchical clustering. This function assumes you denote one distribution as P and one as Q.
}
\usage{
JS_spliced(grid, shiftp, shiftq, phip, phiq, shapep, shapeq,
           ratep, rateq, threshp, threshq, sigmap, sigmaq, xip, xiq)
}
\arguments{
  \item{grid}{Vector of integers over which to compute the JS distance. The minimum of the grid is ideally the minimum count of all samples being compared. The maximum is ideally something very large (e.g. 100,000) in order to all or nearly all of the model density. The grid should include every integer in its range. See Examples.}
  \item{shiftp}{The shift for distribution P.}
  \item{shiftq}{The shift for distribution Q.}
  \item{phip}{The estimated phi for distribution P.}
  \item{phiq}{The estimated phi for distribution Q.}
  \item{shapep}{The estimated gamma shape parameter for distribution P.}
  \item{shapeq}{The estimated gamma shape parameter for distribution Q.}
  \item{ratep}{The estimated gamma rate parameter for distribution P.}
  \item{rateq}{The estimated gamma rate parameter for distribution Q.}
  \item{threshp}{The estimated threshold for distribution P.}
  \item{threshq}{The estimated threshold for distribution Q.}
  \item{sigmap}{The estimated parameter sigma for distribution P.}
  \item{sigmaq}{The estimated parameter sigma for distribution Q.}
  \item{xip}{The estimated parameter xi for distribution P.}
  \item{xiq}{The estimated parameter xi for distribution Q.}
}
\details{
For 2 discrete distributions P and Q, the Jensen-Shannon distance between them is

\eqn{JSD(P,Q) = \sqrt.5 * [\sum(P_i log P_i/M_i)] + \sum(Q_i log Q_i/M_i)}

where

\eqn{M_i= .5 * (P_i + Q_i)}.
}
\value{
The function directly returns the Jensen-Shannon distance between two fitted distributions P and Q.
}
\author{\email{hbk5086@psu.edu}}
\seealso{
\code{\link{JS_spliced}}
}
\examples{
data("repertoires")

# Fit the discrete gamma-gpd spliced model at some selected threshold on 2 samples
fit1 <- fdiscgammagpd(repertoires[[1]],
                        useq = quantile(repertoires[[1]], .8),
                        shift = min(repertoires[[1]]))
fit2 <- fdiscgammagpd(repertoires[[2]],
                        useq = quantile(repertoires[[2]], .8),
                        shift = min(repertoires[[2]]))

# Create a grid of every integer from the minimum count to a large value
# The chosen "large value" here is only 1,000, for the sake of quick computation.
# Ideally, the large value will be at least 100,000
grid <- min(c(repertoires[[1]], repertoires[[2]])):1000

# Compute the Jensen-Shannon distance between fit1 and fit2
dist <- JS_spliced(grid,
                    shiftp = min(repertoires[[1]]), shiftq = min(repertoires[[2]]),
                    phip = fit1$mle['phi'], phiq = fit2$mle['phi'],
                    shapep = fit1$mle['shape'], shapeq = fit2$mle['shape'],
                    ratep = fit1$mle['rate'], rateq = fit2$mle['rate'],
                    threshp = fit1$mle['thresh'], threshq = fit2$mle['thresh'],
                    sigmap = fit1$mle['sigma'], sigmaq = fit2$mle['sigma'],
                    xip = fit1$mle['xi'], xiq = fit2$mle['xi'])

dist
}
