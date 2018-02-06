get_nllh <- function(fits){
    if(all(c("x", "shift", "init", "useq", "nllhuseq", "nllh",
         "optim", "mle") %in% names(fits))){
        nllhs <- fits$nllh
    } else{
        nameList <- map(fits, names)
        strCheck <- vapply(fits, function(X) all(c("x", "init", "useq",
                            "nllhuseq", "nllh", "optim", "mle") %in% names(X)),
                                logical(1))
        if(any(!strCheck)){
            stop("One or more elements of \"fits\" is not a model fit
                 from fdiscgammagpd.")
        }
        nllhs <- map(fits, 'nllh')
    }
    nllhs
}

get_mle <- function(fits){
    if(all(c("x", "shift", "init", "useq", "nllhuseq", "nllh",
             "optim", "mle") %in% names(fits))){
        mles <- fits$mle
    } else{
        nameList <- map(fits, names)
        strCheck <- vapply(fits, function(X) all(c("x", "init", "useq",
                            "nllhuseq", "nllh", "optim", "mle") %in% names(X)),
                               logical(1))
        if(any(!strCheck)){
            stop("One or more elements of \"fits\" is not a model fit
                 from fdiscgammagpd.")
        }
        mles <- map(fits, 'mle')
    }
    mles
}
