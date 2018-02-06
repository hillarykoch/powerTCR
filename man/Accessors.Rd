\name{Accessors}
\alias{get_mle}
\alias{get_nllh}
\title{
Grab important output from a list of fits
}
\description{
These are convenient accessors that will grab important output from a list of fits from fdiscgammagpd. They will grab the maximum likelihood estimates and/or the negative log likelihood for the maximum likelihood estimates.
}
\usage{
get_mle(fits)
get_nllh(fits)
}
\arguments{
  \item{fits}{A list of fits output from fdiscgammagpd.}
}
\value{
A list of out either maximum likelihood estimates (get_mle) or negative log likelihoods (get nllh) corresponding to the list of fits.
}
\author{\email{hbk5086@psu.edu}}
\seealso{
\code{\link{fdiscgammagpd}}
}
\examples{
# Here is a good workflow using fdiscgammagpd:
# Choose quantiles for every sample repertoire in the same manner.
# Then fit the model in the same manner as well.

data("repertoires")

thresholds <- list()
fits <- list()
for(i in 1:2){
    thresholds[[i]] <- unique(round(quantile(repertoires[[i]], c(.8,.85,.9,.95))))
    fits[[i]] <- fdiscgammagpd(repertoires[[i]], useq = thresholds[[i]],
                                shift = min(repertoires[[i]]))
}
names(fits) <- c("fit1", "fit2")

mles <- get_mle(fits)
nllhs <- get_nllh(fits)

mles
nllhs
}