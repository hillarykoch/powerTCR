\name{get_distances}
\alias{get_distances}
\title{
Compute matrix of pairwise Jensen-Shannon Distances
}
\description{
For a list of model fits (either the spliced model or the Desponds et al. model), compute the matrix of Jensen-Shannon distances. This can then be used for clustering or multi-dimensional scaling.
}
\usage{
get_distances(fits, grid, modelType = "Spliced")
}
\arguments{
  \item{fits}{A list of fits output from either fdiscgammagpd or fdesponds.}
  \item{grid}{Vector of integers over which to compute the JS distance. The minimum of the grid is ideally the minimum count of all samples being compared. The maximum is ideally something very large (e.g. 100,000) in order to all or nearly all of the model density. The grid should include every integer in its range. See \code{\link{JS_dist}}.}
  \item{modelType}{Either "Spliced" or "Desponds", depending on what sort of fits you are supplying. Defaults to "Spliced".}
}
\value{
A symmetric matrix of pairwise Jensen-Shannon distances, with 0 on the diagonal.
}

\author{\email{hbk5086@psu.edu}}
\seealso{
\code{\link{JS_dist}, \link{fdiscgammagpd}, \link{fdesponds}}
}
\examples{
# Simulate 3 datasets
set.seed(123)
s1 <- rdiscgammagpd(1000, shape = 3, rate = .15, u = 25, sigma = 15,
                    xi = .5, shift = 1)
s2 <- rdiscgammagpd(1000, shape = 3.1, rate = .14, u = 26, sigma = 15,
                    xi = .6, shift = 1)
s3 <- rdiscgammagpd(1000, shape = 10, rate = .3, u = 45, sigma = 20,
                    xi = .7, shift = 1)
                    
# Fit the spliced model to each
# Here, we use true thresholds for fast computation for this example
# In practice, you need to select a whole sequence of potential thresholds
sim_fits <- list("s1" = fdiscgammagpd(s1, useq = 25),
                 "s2" = fdiscgammagpd(s2, useq = 26),
                 "s3" = fdiscgammagpd(s3, useq = 45))
                 
# Compute the pairwise JS distance between 3 fitted models
grid <- min(c(s1,s2,s3)):10000
distances <- get_distances(sim_fits, grid, modelType="Spliced")
distances
}
