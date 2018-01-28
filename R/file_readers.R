##--------------------------------------------------------------------
## Borrow functionality from the tcR package to help folks read files
##--------------------------------------------------------------------

parseFile <- function(file, format = c('mitcr', 'mitcrbc', 'migec', 'vdjtools',
                                       'immunoseq', 'mixcr', 'imseq'),
                      inframe = TRUE){
    data <- tcR::parse.file(.filename = file, .format = format)

    if(inframe){
        data <- get.inframes(data)
    }

    counts <- data$Read.count
    counts
}

parseFolder <- function(folder, format = c('mitcr', 'mitcrbc', 'migec',
                                           'vdjtools', 'immunoseq', 'mixcr',
                                           'imseq'),
                        inframe = TRUE){
    dataList <- parse.folder(.folderpath = folder, .format = format)

    if(inframe){
        dataList <- map(dataList, get.inframes)
    }

    counts <- map(dataList, "Read.count")
    counts
}
