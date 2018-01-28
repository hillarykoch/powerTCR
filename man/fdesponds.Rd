\name{fdesponds}
\alias{fdesponds}
\title{
Fit the type-I Pareto distribution as according to Desponds et al. (2016)
}
\description{
This function fits a continuous type-I pareto distribution to a vector of count data. Given data x, a threshold Cmin, and letting n be the number of clones greater than u, the shape parameter alpha is computed as

\eqn{n * 1/(\sum log(x/Cmin)) + 1}.

The method considers every possible threshold (that is, every element of the vector unique(x)). The threshold and alpha which minimize the Kolmogorov-Smirnov statistic are selected.
}
\usage{
fdesponds(x)
}
\arguments{
\item{x}{vector of counts.}
}
\value{
\item{min.KS}{The value of the KS statistic for the fitted Pareto distribution.}
\item{Cmin}{The inferred threshold.}
\item{powerlaw.exponent}{The powerlaw exponent. This is equal to \code{pareto.alpha + 1}}
\item{pareto.alpha}{The inferred shape parameter alpha of the fitted Pareto distribution.}
}
\references{
Desponds, Jonathan, Thierry Mora, and Aleksandra M. Walczak. "Fluctuating fitness shapes the clone-size distribution of immune repertoires." Proceedings of the National Academy of Sciences 113.2 (2016): 274-279.
APA
}
\author{\email{hbk5086@psu.edu}}
\examples{
# Fit the model to sample data

data("repertoires")

fit1 <- fdesponds(repertoires[[1]])
fit2 <- fdesponds(repertoires[[2]])

fit1
fit2
}
