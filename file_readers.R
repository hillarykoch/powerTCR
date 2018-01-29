##--------------------------------------------------------------------
## Borrow functionality from the tcR package to help folks read files
##--------------------------------------------------------------------

parseFile <- function(file, format = c('mitcr', 'mitcrbc', 'migec', 'vdjtools',
                                       'immunoseq', 'mixcr', 'imseq'),
                      inframe = TRUE){
    if(!(is.character(file))){
        stop("file must be a character string.")
    }
    
    if(!(format %in%  c('mitcr', 'mitcrbc', 'migec', 'vdjtools', 'immunoseq',
                        'mixcr', 'imseq'))){
        stop("invalid format specified.")
    }
    
    if(!(is.logical(inframe))){
        stop("inframe must be TRUE or FALSE.")
    }
    
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
    if(!(is.character(folder))){
        stop("folder must be a character string.")
    }
    
    if(!(format %in%  c('mitcr', 'mitcrbc', 'migec', 'vdjtools', 'immunoseq',
                        'mixcr', 'imseq'))){
        stop("invalid format specified.")
    }
    
    if(!(is.logical(inframe))){
        stop("inframe must be TRUE or FALSE.")
    }
    
    dataList <- parse.folder(.folderpath = folder, .format = format)

    if(inframe){
        dataList <- map(dataList, get.inframes)
    }

    counts <- map(dataList, "Read.count")
    counts
}
