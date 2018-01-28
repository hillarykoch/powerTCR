# This function does model fitting and threshold selection as described in Desponds et al (2016)
# Details on page 5 of supplement (right column)
# x is a vector of clone sizes

fdesponds <- function(x){
    x <- sort(x)
    Cmins <- unique(x)
    alpha <- rep(NA, length(Cmins))
    KS <- rep(NA, length(Cmins))

    n <- rep(NA, length(Cmins))
    for(i in 1:length(Cmins)){
        d <- x[x > Cmins[i]]

        n[i] <- length(d)
        alpha[i] <- n[i]*(sum(log(d/Cmins[i])))^(-1) + 1

        empir.cdf <- rep(NA, length(d))
        for(j in 1:length(d)){
            empir.cdf[j] <- mean(d <= d[j])
        }

        est.cdf <- VGAM::ppareto(d, scale=Cmins[i], shape = alpha[i]-1)

        KS[i] <- max(abs(est.cdf - empir.cdf))
    }

    min.KS <- min(KS, na.rm = TRUE)
    Cmin <- Cmins[KS == min(KS, na.rm = TRUE)][1]
    powerlaw.exponent <- alpha[KS == min(KS, na.rm = TRUE)]
    powerlaw.exponent <- powerlaw.exponent[!is.na(powerlaw.exponent)]
    out <- c(min.KS, Cmin, powerlaw.exponent, powerlaw.exponent-1)
    names(out) <- c("min.KS", "Cmin", "powerlaw.exponent", "pareto.alpha")
    out
}
