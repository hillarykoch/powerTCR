\name{get_bootstraps}
\alias{get_bootstraps}
\title{
Get bootstrapped fits for a list of fitted models
}
\description{
In order to get confidence bands on parameter estimates, a parametric bootstrap is recommended. This bootstrapping procedure takes bootstraps above and below the threshold separately, retaining the correct proportion of data that are above or below the threshold.
}
\usage{
get_bootstraps(fits, resamples = 1000, cores = 1, gridStyle = "copy")
}
\arguments{
  \item{fits}{A list of fits output from fdiscgammagpd.}
  \item{resamples}{Number of bootstrap replicates to execute for each model fit. Defaults to 1000.}
  \item{cores}{Number of cores to use, if running in parallel. Defaults to 1.}
  \item{gridStyle}{Defines how the sequence of thresholds is selected in the bootstrap fits. If the default "copy", each bootstrapped fit will be computed using the same grid of thresholds from the original fit. Otherwise, an integer can be supplied. If an integer is supplied, the bootstraps will be fit using a grid of thresholds defined by the originally estimated threshold plus or minus the supplied integer.}
}

\value{
If only one fit is passed, get_bootstraps returns a list of length resamples, where each element is a bootstrapped fit output from fdiscgammagpd. If a list of fits is passed, then the output is a list of lists. Each element of that list is a list of length resamples, where each element is a bootstrapped fit output from fdiscgammagpd.
}
\author{\email{hbk5086@psu.edu}}
\seealso{
\code{\link{fdiscgammagpd}}
}
\examples{
data(repertoires)
fits <- lapply(repertoires,
                function(X) fdiscgammagpd(X, useq = unique(round(quantile(X, c(.75,.8,.85,.9))))))
names(fits) <- names(repertoires)

# You should in practice use a large number of resamples, say, 1000
boot <- get_bootstraps(fits, resamples = 10)

mles <- get_mle(boot[[1]])
xi_CI <- quantile(unlist(purrr::map(mles, 'xi')), c(.025, .5, .975))
xi_CI
}

