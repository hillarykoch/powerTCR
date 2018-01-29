##----------------------------------------------------------------------------------
## Density, distribution, and quantile functions, random number generation
## for discrete truncated gamma and discrete gpd
##----------------------------------------------------------------------------------

ddiscgamma <- function(x, shape, rate, thresh, phiu, shift = 0, log = FALSE){
    if(any(x != floor(x))){
        stop("x must be an integer")
    }

    out <- rep(0, length(x))

    up <- pgamma(x+1-shift, shape=shape, rate=rate)
    down <- pgamma(x-shift, shape=shape, rate=rate)

    if(!log){
        b <- pgamma(thresh-shift, shape=shape, rate=rate)
        out[x < thresh] <- ((1-phiu)*(up-down)/b)[x < thresh]
    } else{
        b <- pgamma(thresh-shift, shape=shape, rate=rate, log = TRUE)
        out[x < thresh] <- (log(1-phiu)+log(up-down) - b)[x < thresh]
    }
    out
}

pdiscgamma <- function(q, shape, rate, thresh, phiu, shift = 0){
    probs <- ddiscgamma(0:q, shape, rate, thresh, phiu, shift)
    sum(probs)
}

qdiscgamma <- function(p, shape, rate, thresh, phiu, shift = 0){
    qtrunc(p/(1-phiu), spec = "gamma", a=0,
           b=thresh-shift, shape=shape,rate=rate) %>% floor + shift
}

rdiscgamma <- function(n, shape, rate, thresh, shift = 0){
    rtrunc(n, spec = "gamma", a=0,
           b=thresh-shift, shape=shape, rate=rate) %>% floor + shift
}

ddiscgpd <- function(x, thresh, sigma, xi, phiu, log = FALSE){
    up <- evmix::pgpd(x+1, u=thresh, sigmau=sigma, xi=xi)
    down <- evmix::pgpd(x, u=thresh, sigmau=sigma, xi=xi)

    if(!log){
        phiu*(up-down)
    } else{
        log(phiu) + log(up-down)
    }
}

pdiscgpd <- function(q, thresh, sigma, xi, phiu){
    probs <- ddiscgpd(thresh:q, thresh, sigma, xi, phiu)
    sum(probs)
}

qdiscgpd <- function(p, thresh, sigma, xi, phiu){
    evmix::qgpd(p/phiu, u=thresh, sigmau = sigma, xi=xi) %>% floor
}

rdiscgpd <- function(n, thresh, sigma, xi){
    evmix::rgpd(n, u=thresh, sigmau=sigma, xi=xi) %>% floor
}


##----------------------------------------------------------------------------------
## negative log likelihood and parameter estimation functions
## for discrete truncated gamma and discrete gpd
##----------------------------------------------------------------------------------

discgammanll <- function(param, dat, thresh, phiu, shift=0){
    shape <- exp(param[1])
    rate <- exp(param[2])

    if(any(dat > thresh-1)){ warning("data must be less than the threshold") }

    ll <- log(ddiscgamma(dat, shape, rate, thresh, phiu, shift))

    sum(-ll)
}

discgpdnll <- function(param, dat, thresh, phiu){
    sigma <- exp(param[1])
    xi <- param[2]
    ll <- log(ddiscgpd(dat, thresh, sigma, xi, phiu))

    sum(-ll)
}

fdiscgamma <- function(param, dat, thresh, phiu, shift = 0, method, ...){
    opt <- optim(log(param), discgammanll, dat=dat, thresh=thresh,
                 phiu=phiu, shift=shift, method=method, hessian = TRUE, ...)
    opt
}

fdiscgpd <- function(param, dat, thresh, phiu, method, ...){
    opt <- optim(c(log(param[1]),param[2]), discgpdnll, dat=dat, thresh=thresh,
                 phiu=phiu, method=method, hessian = TRUE, ...)
    opt
}

