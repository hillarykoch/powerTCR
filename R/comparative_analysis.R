# Computes Jensen-Shannon divergence
# of 2 CONTINUOUS distributions at a given point
eval_desponds <- function(t, Cminp, Cminq, alphap, alphaq){
    M <- 0.5*(VGAM::dpareto(t, scale=Cminp, shape=alphap) +
                  VGAM::dpareto(t, scale=Cminq, shape=alphaq))

    one <- VGAM::dpareto(t, scale=Cminp, shape=alphap)
    two <- VGAM::dpareto(t, scale=Cminq, shape=alphaq)

    if(one == 0 & two == 0){
        out <- 0
    } else if(one == 0 & two != 0){
        out <- VGAM::dpareto(t, scale=Cminq, shape=alphaq) *
            (log(VGAM::dpareto(t, scale=Cminq, shape=alphaq))-log(M))
    } else if(one != 0 & two == 0){
        out <- VGAM::dpareto(t, scale=Cminp, shape=alphap) *
            (log(VGAM::dpareto(t, scale=Cminp, shape=alphap))-log(M))
    } else{
        out <- VGAM::dpareto(t, scale=Cminp, shape=alphap) *
            (log(VGAM::dpareto(t, scale=Cminp, shape=alphap))-log(M)) +
            VGAM::dpareto(t, scale=Cminq, shape=alphaq) *
            (log(VGAM::dpareto(t, scale=Cminq, shape=alphaq))-log(M))
    }
    out
}

JS_desponds <- function(grid, Cminp, Cminq, alphap, alphaq){
    lower = min(grid)
    upper = max(grid)

    out <- adaptIntegrate(eval_desponds,
                          lowerLimit = lower, upperLimit = upper,
                          Cminp = Cminp, Cminq = Cminq,
                          alphap = alphap, alphaq = alphaq)$integral %>% sqrt
    out
}

# Computes Jensen-Shannon divergence of 2 distributions
JS_spliced <- function(grid, shiftp, shiftq, phip, phiq, shapep, shapeq, ratep,
                       rateq, threshp, threshq, sigmap, sigmaq, xip, xiq){
    K <- max(grid)

    P <- ddiscgammagpd(min(grid):K, shape = shapep, rate = ratep,
                       u=threshp, sigma = sigmap,
                       xi = xip, phiu = phip, shift=shiftp,
                       log = FALSE)
    adjp <- which(P == 0)
    if(length(adjp) != 0){
        P[adjp] <- evmix::dgpd(adjp+0.5, u=threshp,
                              sigmau = sigmap, xi = xip, phiu = phip)
    }

    Q <- ddiscgammagpd(min(grid):K, shape = shapeq, rate = rateq,
                       u=threshq, sigma = sigmaq,
                       xi = xiq, phiu = phiq, shift=shiftq,
                       log = FALSE)
    adjq <- which(Q == 0)
    if(length(adjq) != 0){
        Q[adjq] <- evmix::dgpd(adjq+0.5, u=threshq,
                              sigmau = sigmaq, xi = xiq, phiu = phiq)
    }

    M <- 0.5*(P+Q)
    pzero <- which(P == 0)
    qzero <- which(Q == 0)

    sum1 <- sum2 <- rep(NA, length(grid))
    sum1 <- P*(log(P) - log(M))
    sum2 <- Q*(log(Q) - log(M))

    if(length(intersect(pzero, qzero)) != 0){
        sum1[intersect(pzero, qzero)] <- 0
        sum2[intersect(pzero, qzero)] <- 0
    }
    if(length(setdiff(pzero, qzero)) != 0){
        sum1[setdiff(pzero, qzero)] <- 0
    }
    if(length(setdiff(qzero, pzero)) != 0){
        sum2[setdiff(qzero, pzero)] <- 0
    }

    out <- sqrt(0.5*(sum(sum1) + sum(sum2)))
    out
}

clusterPlot <- function(distances, method = c("complete", "ward.D", "ward.D2",
                                              "single", "average", "mcquitty",
                                              "median", "centroid")){
    if(length(method) > 1){
        warning("More than one clustering method given. Only the first listed
                will be executed.")
        method <- method[1]
    }
    plot(hclust(as.dist(distances), method = method), sub = "", xlab = "")
}

