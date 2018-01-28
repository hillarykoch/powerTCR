\name{repertoires}
\alias{repertoires}
\docType{data}
\title{
Two toy examples of sample TCR repertoires.
}
\description{
This data set gives to toy examples of TCR repertoires. Sample "samp1" contains 1,000 clones with a total of 26,288 sequenced T cells. Sample "samp2" contains 800 clones with a total of 24,267 sequenced T cells. These samples have been sorted from largest to smallest clone size.
}
\usage{data("repertoires")}
\format{
The format is:

List of 2
    $ samp1: num [1:1000] 1445 451 309 ...
    $ samp2: num [1:800] 2781 450 447 ...
}
\examples{
data(repertoires)
n1 <- length(repertoires$samp1)
n2 <- length(repertoires$samp2)

# Generates plot on log-log scale
par(mfrow = c(2,1))
plot(log(repertoires$samp1), log(1:n1))
plot(log(repertoires$samp2), log(1:n2))
}
\keyword{datasets}
