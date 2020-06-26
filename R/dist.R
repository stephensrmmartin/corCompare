##' Density, distribution function, quantile function, and random generation for the transformed beta distribution with mean equal to mu and shape parameter psi.
##'
##' The transformed beta distribution is the beta distribution, shifted and scaled to be defined between -1 and 1.
##' The tbeta distribution is therefore useful for distributions on individual correlations.
##' Importantly, it is not useful as a distribution for elements within a correlation matrix beyond 2x2, as it has no guarantee that the matrix is positive definite.
##'
##' \deqn{E(x) = mu}
##' \deqn{V(x) = \frac{(1 + \mu)(1 - \mu)}{\psi + 1}}
##'
##' The tbeta distribution is a rescaled, shifted beta distribution under the mean-shape parameterization.
##' \deqn{p(x|mu, psi) = \frac{(.5(x + 1))^{.5(\mu + 1)\psi - 1}(1 - .5(x + 1))^{\psi(1 - .5(\mu + 1)) - 1}}{2\Beta(.5(\mu + 1)\psi, \psi(1 - .5(\mu + 1)))}}
##'
##' Or more simply: It is a beta distribution, shifted by -.5, and multiplied by 2.
##' The density is therefore equivalent to inverse-transforming the variate and mu, and solving for the Beta alpha and beta parameters.
##' The jacobian of this adjustment is .5.
##'
##' The tbeta_diff distribution is the distribution of differences in tbeta-distributed variables.
##' That is, if X ~ tbeta(mu_x, psi_x) and Y ~ tbeta(mu_y, psi_y), then (X - Y) ~ tbeta_diff(mu_x, mu_y, psi_x, psi_y).
##' E(X - Y) = mu_x - mu_y
##' V(X - Y) = V(X) + V(Y), where V() is defined above.
##' 
##' @title Transformed beta distribution
##' @param x Numeric vector. (-1, 1) for tbeta, and (-2, 2) for tbeta_diff.
##' @param mu Numeric vector (-1, 1). Mean parameter.
##' @param psi Numeric vector (0, Inf). Shape parameter.
##' @param log Logical (Default: FALSE). Whether to take the log of the probability.
##' @param p Numeric vector. Probabilities.
##' @param q Numeric vector. Quantiles.
##' @param n Integer. Number of random samples.
##' @param mus Numeric vector (Length 2). Defined in (-1, 1). Mu parameters for the differenced beta variates.
##' @param psis Numeric vector (length 2). Defined in (0, Inf). Psi parameters for the differenced beta variates.
##' @author Stephen R. Martin
##' @name tbeta
NULL

##' @export
##' @rdname tbeta
dtbeta <- function(x, mu, psi, log = FALSE) {
    xt <- (x + 1) / 2
    mut <- (mu + 1) / 2
    a <- mut * psi
    b <- psi - a
    d <- dbeta(xt, a, b, log = TRUE)
    d <- d + log(.5)
    if(log) {
        return(d)
    } else {
        return(exp(d))
    }
}
##' @rdname tbeta
##' @export
qtbeta <- function(p, mu, psi) {
    mut <- (mu + 1) / 2
    a <- mut * psi
    b <- psi - a
    q.raw <- qbeta(p, a, b)
    qt <- q.raw * 2 - 1
    return(qt)
}
##' @rdname tbeta
##' @export
ptbeta <- function(q, mu, psi) {
    qt <- (q + 1) / 2
    mut <- (mu + 1) / 2
    a <- mut * psi
    b <- psi - a
    p <- pbeta(qt, a, b)
    return(p)
}

##' @export
##' @rdname tbeta
rtbeta <- function(n, mu, psi) {
    mut <- (mu + 1) / 2
    a <- mut * psi
    b <- psi - a
    samps <- rbeta(n, a, b)
    samps <- (samps - .5) * 2
    samps
}

##' @rdname tbeta
##' @export
dtbeta_diff <- function(x, mus, psis, log = FALSE) {
    out <- .tbeta_prior_diff_exact(x, mus, psis)
    if(log) {
        out <- log(out)
    }
    return(out)
}

##' @rdname tbeta
##' @export
ptbeta_diff <- function(q, mus, psis, lower.tail = TRUE) {
    sapply(q, function(x) {
        if(lower.tail) {
            integrate(dtbeta_diff, -2, x, mus = mus, psis = psis)$value
        } else {
            integrate(dtbeta_diff, x, 2, mus = mus, psis = psis)$value
        }
    })
    
}

##' @rdname tbeta
##' @export
qtbeta_diff <- function(p, mus, psis) {
    f <- function(q, mus, psis, p) {
        abs(ptbeta_diff(q, mus, psis) - p)
    }
    opt <- sapply(p, function(x){
        optimize(f, c(-2, 2), mus = mus, psis = psis, p = x)$minimum
    })
    return(opt)

}
