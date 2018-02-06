test_density <- function(){
    checkEqualsNumeric(ddiscgammagpd(-1, shape=1, rate=1, u=5, sigma=10, xi=.5), 0)
    checkEqualsNumeric(ddiscgammagpd(2, shape=1, rate=1, u=5,
                              sigma=10, xi=.5, shift = 3), 0)
}

test_distribution <- function(){
    checkEqualsNumeric(pdiscgammagpd(24, shape=5, rate=.25, u=25,
                        sigma=15, xi=.5, phiu=.2, shift=1), .8)
    checkEqualsNumeric(pdiscgammagpd(1000000, shape=5, rate=.25, u=25,
                                     sigma=15, xi=.5, phiu=.2, shift=1), 1, tol = .001)
}

test_quantiles <- function(){
    q <- qdiscgammagpd(.1, shape = 5, rate = .25, u = 25,
                       sigma = 15, xi = .5, shift = 1)
    p <- pdiscgammagpd(q, shape=5, rate=.25, u=25, sigma=15, xi=.5, shift=1)
    checkEqualsNumeric(p,.1, tol = .1)
}

test_fitpass <- function(){
    wrong <- list("a" = 1, "b" = 2)
    right1 <- fdiscgammagpd(repertoires[[1]], useq = 15, shift = 1)
    right2 <- fdiscgammagpd(repertoires[[2]], useq = 15, shift = 1)
    checkException(rdiscgammagpd(10, fit=wrong),
                    msg = "\"fit\" is not of the correct structure. It must be
                        one model fit from fdiscgammagpd.")
    checkTrue(is.numeric(JS_dist(right1, right2, 1:1000, "Spliced")))
    checkException(JS_dist(right1, right2, 1:1000, "Desponds"),
                   msg = "\"fit1\" is not of the correct structure. It must be a model
                         fit from fdesponds.")
}

