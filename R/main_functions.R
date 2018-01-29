##-----------------------------------------------------------------------------------------------------
## Fitting function, random number generation, quantile function, density, and distribution functions
##-----------------------------------------------------------------------------------------------------

fdiscgammagpd <- function(x, useq, shift = NULL, pvector=NULL,
                          std.err = TRUE, method = "Nelder-Mead", ...){
    if(!is(x, "numeric")){
        stop("x must be numeric.")
    }
    
    if(!is(shift, "numeric")){
        stop("shift must be numeric.")
    }
    
    if(!is(useq, "numeric")){
        stop("useq must be numeric.")
    }
    
    if(any(x != round(x))){
        stop("all elements in x must be integers.")
    }
    
    if(any(shift != round(shift))){
        stop("shift must be an integer.")
    }
    
    if(any(useq != round(shift))){
        stop("all elements in useq must be integers.")
    }
    
    if(!(length(pvector) == 5)){
        stop("pvector must contain 5 elements.")
    }
    
    if(!(is.logical(std.err))){
        stop("std.err must be TRUE or FALSE.")
    }
    
    if(!(method %in% c("Nelder-Mead","BFGS",
                       "CG", "L-BFGS-B", "SANN", "Brent"))){
        stop("invalid method supplied.")
    }
    
    if(is.null(shift)){
        shift <- min(x)
    }

    if (is.null(pvector)) {
        pvector <- rep(NA,5)
        s <- log(mean(x+0.5))-mean(log(x+0.5))
        k <- (3-s + sqrt((s-3)^2 + 24*s))/12/s
        pvector[1] <- k
        pvector[2] <- k/mean(x)
        pvector[3] <- as.vector(quantile(x, 0.9))

        xu <- x[x>=pvector[3]]
        initfgpd <- evmix::fgpd(xu, min(xu)-10^(-5))
        pvector[4] <- initfgpd$mle[1]
        pvector[5] <- initfgpd$mle[2]
    }

    bulk <- lapply(1:length(useq),
                   function(idx,x,useq) x < useq[idx], x=x, useq=useq)
    tail <- lapply(1:length(useq),
                   function(idx,x,useq) x >= useq[idx], x=x, useq=useq)
    phiu <- lapply(1:length(useq),
                   function(idx,tail) mean(tail[[idx]]), tail=tail)

    gammfit <- list()
    gpdfit <- list()
    nllhu <- rep(NA, length(useq))
    for(i in 1:length(useq)){
        gammfit[[i]] <- tryCatch(expr = fdiscgamma(pvector[1:2],x[bulk[[i]]],
                                                   useq[i],
                                                   phiu[[i]],
                                                   shift,
                                                   method=method),
                                 error = function(err)
                                     message("gamma part could not be fit at a
                                             specified threshold.")
                                 NA)
        gpdfit[[i]] <- tryCatch(expr = fdiscgpd(pvector[4:5],
                                                x[tail[[i]]],
                                                useq[i],
                                                phiu[[i]],
                                                method=method),
                                error = function(err) {
                                    pvec3 <- as.vector(quantile(x,1-phiu[[i]]))
                                    xu <- x[x>=pvec3]
                                    initfgpd.adj <-
                                        evmix::fgpd(x, min(xu)-10^(-5))
                                    pvec4 <- initfgpd.adj$mle[1]
                                    pvec5 <- initfgpd.adj$mle[2]
                                    tryCatch(expr = fdiscgpd(c(pvec4,pvec5),
                                                             x[tail[[i]]],
                                                             useq[i],
                                                             phiu[[i]],
                                                             method=method),
                                             error = function(err2) 
                                             message("GPD part could not be fit
                                                     at a specified
                                                     threshold.")
                                             NA)
                                })
        nllhu[i] <- tryCatch(expr = gammfit[[i]]$value + gpdfit[[i]]$value,
                             error = function(err) NA)
    }

    bestfit <- which.min(nllhu)
    fit.out <- list(gammfit[[bestfit]], gpdfit[[bestfit]])
    names(fit.out) <- c("bulk","tail")
    mle <- c(mean(x >= useq[bestfit]),
             exp(fit.out$bulk$par),
             useq[bestfit],
             exp(fit.out$tail$par[1]),
             fit.out$tail$par[2])
    names(mle) <- c("phi","shape","rate","thresh","sigma","xi")
    if(std.err){
        H <- fit.out$bulk$hessian %>% rbind(matrix(rep(0,4),nrow = 2)) %>%
            cbind(rbind(matrix(rep(0,4),nrow = 2),fit.out$tail$hessian))
        fisherInf <- tryCatch(expr = solve(H), error = function(err) NA)
        out <- list(x = as.vector(x), init = as.vector(pvector),
                    useq = useq, nllhuseq = nllhu,
                    optim = fit.out, nllh = nllhu[bestfit],
                    mle=mle, fisherInformation = fisherInf)
    } else{
        out <- list(x = as.vector(x), init = as.vector(pvector),
                    useq = useq, nllhuseq = nllhu,
                    optim = fit.out, nllh = nllhu[bestfit],
                    mle=mle)
    }
    out
}

# Need this for rgammagpd
qdiscgammagpd <- function(p, shape, rate, u, sigma, xi, phiu=NULL, shift = 0){
    if(!is(x, "numeric")){
        stop("x must be numeric.")
    }
    
    if(!is(shift, "numeric")){
        stop("shift must be numeric.")
    }
    
    if(any(x != floor(x))){
        stop("x must be an integer")
    }
    
    if(any(shift != round(shift))){
        stop("shift must be an integer.")
    }
    
    if(any(c(shape, rate, sigma) <= 0)){
        stop("shape, rate, and sigma must all be positive.")
    }
    
    if(phiu < 0 | phiu > 1){
        stop("phiu must be in [0,1].")
    }
    
    if(p < 0 | p > 1){
        stop("p must be in [0,1].")
    }
    
    if(is.null(phiu)){
        phiu <- 1-pdiscgamma(u-1, shape=shape, rate=rate,
                             thresh=Inf, phiu = 0, shift=shift)
    }
    phib <- (1-phiu)
    nb <- sum(p < phib)
    nu <- sum(p >= phib)

    q <- p
    if(nb > 0){
        q[p < phib] <- qdiscgamma(p[which(p < phib)], shape=shape,
                                  rate=rate, thresh=u, phiu=phiu, shift=shift)
    }
    if(nu > 0){
        q[p >= phib] <- qdiscgpd((p[p >= phib]-phib), thresh=u,
                                 sigma=sigma, xi=xi, phiu=phiu)
    }
    q
}

# Generate data from the model!!
rdiscgammagpd <- function(n, shape, rate, u, sigma, xi, phiu=NULL, shift = 0){
    if(!is(x, "numeric")){
        stop("x must be numeric.")
    }
    
    if(!is(shift, "numeric")){
        stop("shift must be numeric.")
    }
    
    if(any(x != floor(x))){
        stop("x must be an integer")
    }
    
    if(any(shift != round(shift))){
        stop("shift must be an integer.")
    }
    
    if(any(c(shape, rate, sigma) <= 0)){
        stop("shape, rate, and sigma must all be positive.")
    }
    
    if(phiu < 0 | phiu > 1){
        stop("phiu must be in [0,1].")
    }
    
    if(is.null(phiu)){
        phiu <- 1-pdiscgamma(u-1, shape=shape, rate=rate,
                             thresh=Inf, phiu = 0, shift=shift)
    }
    p <- runif(n)
    qdiscgammagpd(p, shape=shape, rate=rate, u=u,
                  sigma=sigma, xi=xi, phiu=phiu, shift=shift)
}

ddiscgammagpd <- function(x, shape, rate, u, sigma, xi,
                          phiu=NULL, shift = 0, log = FALSE){
    if(!is(x, "numeric")){
        stop("x must be numeric.")
    }
    
    if(!is(shift, "numeric")){
        stop("shift must be numeric.")
    }
    
    if(any(x != floor(x))){
        stop("x must be an integer")
    }
    
    if(any(shift != round(shift))){
        stop("shift must be an integer.")
    }
    
    if(any(c(shape, rate, sigma) <= 0)){
        stop("shape, rate, and sigma must all be positive.")
    }
    
    if(phiu < 0 | phiu > 1){
        stop("phiu must be in [0,1].")
    }
    
    if(is.null(phiu)){
        phiu <- 1-pdiscgamma(u-1, shape=shape, rate=rate,
                             thresh=Inf, phiu = 0, shift=shift)
    }

    out <- rep(NA, length(x))

    if(sum(x>=u) != 0){
        out[x>=u] <- ddiscgpd(x[x>=u], u, sigma, xi, phiu, log=log)
    }

    if(sum(x<u) != 0){
        out[x<u] <- ddiscgamma(x[x<u], shape, rate, u, phiu, shift, log=log)
    }
    out
}

pdiscgammagpd <- function(q, shape, rate, u, sigma, xi, phiu=NULL, shift=0){
    if(!is(x, "numeric")){
        stop("x must be numeric.")
    }
    
    if(!is(shift, "numeric")){
        stop("shift must be numeric.")
    }
    
    if(any(x != floor(x))){
        stop("x must be an integer")
    }
    
    if(any(shift != round(shift))){
        stop("shift must be an integer.")
    }
    
    if(any(c(shape, rate, sigma) <= 0)){
        stop("shape, rate, and sigma must all be positive.")
    }
    
    if(phiu < 0 | phiu > 1){
        stop("phiu must be in [0,1].")
    }
    
    probs <- sapply(q, function(q)
        ddiscgammagpd(0:q, shape, rate, u, sigma, xi, phiu, shift) %>%
            sum)
    probs
}
