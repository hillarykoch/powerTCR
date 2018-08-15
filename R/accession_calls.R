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

get_bootstraps <- function(fits,
                           resamples = 1000,
                           cores = 1,
                           gridStyle = "copy"){
    if(!all(c(length(resamples), length(cores)) == 1)){
        stop("cores and resamples must be of length 1.")
    }
    if(!all(is.numeric(c(cores, resamples)))){
        stop("cores and resamples must be integers.")
    }


    if(all(c("x", "shift", "init", "useq", "nllhuseq", "nllh",
             "optim", "mle") %in% names(fits))){

        t1 <- sum(fits$x >= fits$mle['thresh'])
        b1 <- sum(fits$x < fits$mle['thresh'])

        resamp.t1 <- sample(fits$x[fits$x >= fits$mle['thresh']],
                            t1*resamples, replace = T)
        resamp.b1 <- sample(fits$x[fits$x < fits$mle['thresh']],
                            b1*resamples, replace = T)

        resample1 <- matrix(rep(NA, length(fits$x)*resamples), ncol = resamples)
        resample1[1:t1,] <- resamp.t1
        resample1[(t1+1):nrow(resample1),] <- resamp.b1
        resample1 <- data.frame(sapply(seq(resamples),
                                       function(X) sort(resample1[,X],
                                                        decreasing = TRUE)))

        if(gridStyle=="copy"){
            thresh.seq <- fits$useq
        } else if(is.numeric(gridStyle)){
            thresh.seq <- (fits$mle['thresh']-round(gridStyle)):
                (fits$mle['thresh']+round(gridStyle))
            thresh.seq <- thresh.seq[thresh.seq %in% fits$useq]
        } else{
            stop("gridStyle must be either \"copy\" or a positive integer")
        }

        bootstraps <- list()
        if(cores > 1){
            cluster <- makeCluster(cores)
            registerDoParallel(cluster)
            clusterEvalQ(cluster, library(powerTCR))

            bootstraps <- foreach(j = 1:resamples,
                                  .packages = "foreach") %dopar% {
                out <- fdiscgammagpd(resample1[,j],
                                     useq = thresh.seq,
                                     shift = min(resample1[,j]))
                out
            }
            stopCluster(cluster)
        } else if(cores == 1){
            for(j in seq(resamples)){
                bootstraps[[j]] <- fdiscgammagpd(resample1[,j],
                                                 useq = thresh.seq,
                                                 shift = min(resample1[,j]))
            }
        } else {
            stop("cores must be a positive integer")
        }
    } else{
        nameList <- map(fits, names)
        strCheck <- vapply(fits, function(X) {
            all(c("x", "init", "useq","nllhuseq",
                  "nllh", "optim", "mle") %in% names(X))
            },
                           logical(1))
        if(any(!strCheck)){
            stop("One or more elements of \"fits\" is not a model fit
                 from fdiscgammagpd.")
        }

        bootstraps <- list()
        for(i in seq_along(fits)){
            t1 <- sum(fits[[i]]$x >= fits[[i]]$mle['thresh'])
            b1 <- sum(fits[[i]]$x < fits[[i]]$mle['thresh'])

            resamp.t1 <- sample(fits[[i]]$x[fits[[i]]$x >=
                                                fits[[i]]$mle['thresh']],
                                t1*resamples, replace = T)
            resamp.b1 <- sample(fits[[i]]$x[fits[[i]]$x <
                                                fits[[i]]$mle['thresh']],
                                b1*resamples, replace = T)

            resample1 <- matrix(rep(NA, length(fits[[i]]$x)*resamples),
                                ncol = resamples)
            resample1[1:t1,] <- resamp.t1
            resample1[(t1+1):nrow(resample1),] <- resamp.b1
            resample1 <- data.frame(sapply(seq(resamples),
                                           function(X) sort(resample1[,X],
                                                            decreasing = TRUE)))

            if(gridStyle=="copy"){
                thresh.seq <- fits[[i]]$useq
            } else if(is.numeric(gridStyle)){
                thresh.seq <- (fits[[i]]$mle['thresh']-round(gridStyle)):
                    (fits[[i]]$mle['thresh']+round(gridStyle))
                thresh.seq <- thresh.seq[thresh.seq %in% fits[[i]]$useq]
            } else{
                stop("gridStyle must be either \"copy\" or a positive integer")
            }

            if(cores > 1){
                cluster <- makeCluster(cores)
                registerDoParallel(cluster)
                clusterEvalQ(cluster, library(powerTCR))

                bootstraps[[i]] <- foreach(j = 1:resamples,
                                      .packages = "foreach") %dopar% {
                                          out <-
                                              fdiscgammagpd(resample1[,j],
                                                            useq = thresh.seq,
                                                            shift = min(
                                                                resample1[,j])
                                                            )
                                          out
                                      }
                stopCluster(cluster)
                names(bootstraps) <0
            } else if(cores == 1){
                bootstraps[[i]] <- list()
                for(j in seq(resamples)){
                    bootstraps[[i]][[j]] <- fdiscgammagpd(resample1[,j],
                                                     useq = thresh.seq,
                                                     shift = min(resample1[,j]))
                }
            } else {
                stop("cores must be a positive integer")
            }
        }
        names(bootstraps) <- names(fits)
    }
    bootstraps
}

get_distances <- function(fits, grid, modelType = "Spliced"){
    if(!is(grid, "numeric")){
        stop("grid must be numeric.")
    }
    if(any(grid != round(grid))){
        stop("all elements in grid must be integers.")
    }
    if(!(modelType %in% c("Spliced", "Desponds"))){
        stop("modelType must be either \"Spliced\" or \"Desponds\".")
    }

    distances <- matrix(rep(0, length(fits)^2), nrow = length(fits))
    if(!is.null(names(fits))){
        rownames(distances) <- colnames(distances) <- names(fits)
    }

    for(i in seq_len((length(fits)-1))){
        for(j in (i+1):length(fits)){
            distances[i,j] <- JS_dist(fits[[i]],
                                      fits[[j]],
                                      grid,
                                      modelType = modelType)
        }
    }
    distances <- distances + t(distances)
    distances
}

get_diversity <- function(fits){
    rich <- sapply(fits, function(X) specnumber(X$x))
    shannon <- sapply(fits, function(X) diversity(X$x, index = "shannon"))
    clon <- sapply(fits, function(X) 1-vegan::diversity(X$x)/log(vegan::specnumber(X$x)))
    tailprop <- sapply(fits, function(X) sum(X$x[X$x >= X$mle['thresh']])/sum(X$x))

    data.frame("richness" = rich, "shannon" = shannon, "clonality" = clon, "prop_stim" = tailprop)
}

