\name{parseFile}
\alias{parseFile}
\alias{parseFolder}

\title{
Load in and parse TCR data files for use by powerTCR
}
\description{
These functions leverage the parse.file and \link[tcR]{parse.folder} functions from the tcR package. They are wrappers that output data in the format taken by powerTCR.
}
\usage{
parseFile(file, format = c('mitcr', 'mitcrbc', 'migec', 'vdjtools',
                           'immunoseq', 'mixcr', 'imseq'),
                            inframe = TRUE)
parseFolder(folder, format = c('mitcr', 'mitcrbc', 'migec', 'vdjtools',
                               'immunoseq', 'mixcr','imseq'),
                                inframe = TRUE)
}

\arguments{
  \item{file}{Path to input file with TCR repertoire sample data.}
  \item{folder}{Path to input folder with one or more input files of TCR repertoire sample data.}
  \item{format}{String specifying the input format. The formats MiTCR ("mitcr"), MiTCR with UMIs ("mitcrbc"), MiGEC ("migec"), VDJtools ("vdjtools"), ImmunoSEQ (formats 1, 2, and 3), MiXCR ("mixcr"), and IMSEQ ("imseq") are currently supported. For more information on formatting, see \link[tcR]{parse.folder}.}
  \item{inframe}{Logical. Should counts only from in-frame sequences be returned? Defaults to TRUE.}
}
\value{
\item{\code{parseFolder} returns a list. Every element of the list is a vector of counts corresponding to a sample repertoire.}

\item{\code{parseFile} returns a vector of counts corresponding to the sample repertoire.}
}
\references{
Nazarov, Vadim I., et al. "tcR: an R package for T cell receptor repertoire advanced data analysis." BMC bioinformatics 16.1 (2015): 175.
APA
}
