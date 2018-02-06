\name{clusterPlot}
\alias{clusterPlot}
\title{
Visualize hierarchical clustering of samples
}
\description{
This function is just a simple wrapper for the \link{hclust} function. It takes
a symmetrix matrix displaying pairwise distances between samples and outputs
a plot of the hierarchical clustering using specified linkage. Linkage
method defaults to \code{complete}. Note that the distances must be given as
a matrix object, not a distance object.
}
\usage{
clusterPlot(distances, method = c("complete", "ward.D", "ward.D2", "single",
            "average", "mcquitty", "median", "centroid"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distances}{A symmetric matrix containined the Jensen-Shannon distance between pairs of distributions.}
  \item{method}{Linkage method, as in \link{hclust}}
}
\value{A basic plot of the induced hierarchical clustering.}
\author{\email{hbk5086@psu.edu}}
\note{
The distances must be given as a matrix object, not a distance object. The
distance between a distribution and itself is 0. This corresponds to a
matrix diagonal of 0.
}
\seealso{
\link{JS_spliced}, \link{JS_desponds}
}
\examples{

# Simulate 3 sampled individuals
set.seed(123)
s1 <- rdiscgammagpd(1000, shape = 3, rate = .15, u = 25, sigma = 15, xi = .5, shift = 1)
s2 <- rdiscgammagpd(1000, shape = 3.1, rate = .14, u = 26, sigma = 15, xi = .6, shift = 1)
s3 <- rdiscgammagpd(1000, shape = 10, rate = .3, u = 45, sigma = 20, xi = .7, shift = 1)

# Fit model to the data at the true thresholds
fits <- list("fit1" = fdiscgammagpd(s1, useq = 25),
             "fit2" = fdiscgammagpd(s2, useq = 26),
             "fit3" = fdiscgammagpd(s3, useq = 45))

# Compute the pairwise JS distance between 3 fitted models
distances <- matrix(rep(0, 9), nrow = 3)
colnames(distances) <- rownames(distances) <- c("s1", "s2","s3")
grid <- min(c(s1,s2,s3)):10000
for(i in 1:2){
    for(j in (i+1):3){
        distances[i,j] <- JS_spliced(grid,
                                     shiftp = min(fits[[i]]$x),
                                     shiftq = min(fits[[j]]$x),
                                     phip = fits[[i]]$mle['phi'],
                                     phiq = fits[[j]]$mle['phi'],
                                     shapep = fits[[i]]$mle['shape'],
                                     shapeq = fits[[j]]$mle['shape'],
                                     ratep = fits[[i]]$mle['rate'],
                                     rateq = fits[[j]]$mle['rate'],
                                     threshp = fits[[i]]$mle['thresh'],
                                     threshq = fits[[j]]$mle['thresh'],
                                     sigmap = fits[[i]]$mle['sigma'],
                                     sigmaq = fits[[j]]$mle['sigma'],
                                     xip = fits[[i]]$mle['xi'], 
                                     xiq = fits[[j]]$mle['xi'])
    }
}

# Distances are symmetric
distances <- distances + t(distances)

# Perform clustering. Note that s1 and s2 were generated using similar
# parameters, so we might expect them to be clustered together
\dontrun{clusterPlot(distances, method = c("complete"))}
}
