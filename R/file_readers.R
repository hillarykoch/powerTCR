##--------------------------------------------------------------------
## Borrow functionality from the immunarch package to read in files
##--------------------------------------------------------------------

parseFile <- function(path, inframe = TRUE){
    if(!(is.character(path))){
        stop("file must be a character string.")
    }
    
    if(!(is.logical(inframe))){
        stop("inframe must be TRUE or FALSE.")
    }
    
    data <- immunarch::repLoad(.path = path)$data[[1]]
    counts <- data$Clones

    if(inframe){
        data <- immunarch::coding(data)
    }

    counts <- data$Clones
    counts
}