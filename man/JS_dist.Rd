\name{JS_dist}
\alias{JS_dist}
\title{
Compute the Jensen-Shannon distance between two model fits
}
\description{
This function is a convenient wrapper for JS_spliced and JS_desponds. After models have been fit to your samples, the pairwise JS distance can be computed between them. This function takes two model fits and outputs the JS distance between them. The model fits must be of the same type. That is, they are both fits from the spliced threshold model, or they are both fits from the Desponds et al. model. When all pairwise distances have been computed, they can be used to do hierarchical clustering.
}
\usage{
JS_dist(fit1, fit2, grid, modelType = "Spliced")
}
\arguments{
  \item{fit1}{A fit from the specified modelType.}
  \item{fit2}{A fit from the specified modelType.}
  \item{grid}{Vector of integers over which to compute the JS distance. The minimum of the grid is ideally the minimum count of all samples being compared. The maximum is ideally something very large (e.g. 100,000) in order to all or nearly all of the model density. The grid should include every integer in its range. See Examples.}
  \item{modelType}{The type of model fit1 and fit2 are from. If they were generated using fdiscgammagpd, the type of model is "Spliced". If they were generated using fdesponds, the type of model is "Desponds". Defaults to "Spliced".}
}
\details{
For 2 discrete distributions P and Q, the Jensen-Shannon distance between them is

\eqn{JSD(P,Q) = \sqrt.5 * [\sum(P_i log P_i/M_i)] + \sum(Q_i log Q_i/M_i)}

where

\eqn{M_i= .5 * (P_i + Q_i)}.
}
\value{
The function directly returns the Jensen-Shannon distance between two fitted distributions.
}
\author{\email{hbk5086@psu.edu}}
\seealso{
\code{\link{JS_spliced}, \link{JS_desponds}}
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
dist <- JS_dist(fit1, fit2, grid, "Spliced")
dist
}
