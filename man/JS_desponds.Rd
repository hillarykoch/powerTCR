\name{JS_desponds}
\alias{JS_desponds}
\title{
Compute the Jensen-Shannon distance between two fitted distributions using the Desponds model
}
\description{
After the Desponds et al. (2016) model havs been fit to your samples, the pairwise JS distance can be computed between them. This function takes the fitted model parameters from 2 distributions and computes the JS distance between them. When all pairwise distances have been computed, they can be used to do hierarchical clustering. This function assumes you denote one distribution as P and one as Q.
}
\usage{
JS_desponds(grid, Cminp, Cminq, alphap, alphaq)
}
\arguments{
\item{grid}{Vector of integers over which to compute the JS distance. The minimum of the grid is ideally the minimum count of all samples being compared. The maximum is ideally something very large (e.g. 100,000) in order to all or nearly all of the model density. The grid should include every integer in its range. See Examples.}
\item{Cminp}{The estimated threshold for distribution P.}
\item{Cminq}{The estimated threshold for distribution Q.}
\item{alphap}{The estimated parameter alpha for distribution P.}
\item{alphaq}{The estimated parameter alpha for distribution Q.}
}
\details{
For 2 discrete distributions P and Q, the Jensen-Shannon distance between them is

\eqn{JSD(P,Q) =\sqrt .5 * \int [P(t) log P(t)/M(t)] + \int [Q(t) log Q(t)/M(t)] dt}

where

\eqn{M(t) = .5 * (P(t) + Q(t))}.

}
\value{
The function directly returns the Jensen-Shannon distance between two fitted distributions P and Q.
}
\references{
Desponds, Jonathan, Thierry Mora, and Aleksandra M. Walczak. "Fluctuating fitness shapes the clone-size distribution of immune repertoires." Proceedings of the National Academy of Sciences 113.2 (2016): 274-279.
APA
}
\author{\email{hbk5086@psu.edu}}
\seealso{
\code{\link{JS_spliced}, \link{JS_dist}}
}
\examples{
data("repertoires")

# Fit the discrete gamma-gpd spliced model at some selected threshold on 2 samples
fit1 <- fdesponds(repertoires[[1]])
fit2 <- fdesponds(repertoires[[2]])

# Create a grid of every integer from the minimum threshold to a large value
# When comparing many distributions in advance of clustering,
# the same grid should be used across every comparison
# The chosen "large value" here is only 1,000, for the sake of quick computation.
# Ideally, the large value will be at least 100,000
grid <- min(c(fit1['Cmin'], fit2['Cmin'])):1000

# Compute the Jensen-Shannon distance between fit1 and fit2
dist <- JS_desponds(grid,
                    Cminp = fit1['Cmin'],
                    Cminq = fit2['Cmin'],
                    alphap = fit1['pareto.alpha'],
                    alphaq = fit2['pareto.alpha'])

dist
}
