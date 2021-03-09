\name{parseFile}
\alias{parseFile}
\alias{parseFolder}

\title{
Load in and parse TCR data files for use by powerTCR
}
\description{
This function leverages  \link[immunarch]{repLoad} functions from the immunarch package. They are wrappers that output data in the format taken by powerTCR.
}
\usage{
parseFile(path, inframe = TRUE)
}

\arguments{
  \item{path}{Path to input file with TCR repertoire sample data.}
  \item{inframe}{Logical. Should counts only from in-frame sequences be returned? Defaults to TRUE.}
}
\value{
\code{parseFile} returns a vector of counts corresponding to the sample repertoire.
}
\details{
See the immunarch package on GitHub: https://github.com/immunomind/immunarch
}

\references{
Nazarov, Vadim I., et al. "tcR: an R package for T cell receptor repertoire advanced data analysis." BMC bioinformatics 16.1 (2015): 175.
}
