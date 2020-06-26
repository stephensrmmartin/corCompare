##' Density and sampling function for transformed beta distribution.
##'
##' The transformed beta distribution is a mean-scale shifted beta distribution ranging from -1 to 1.
##' The tbeta distribution is useful as a prior for a correlation.
##' @title Transformed Beta distribution.
##' @param x Numeric vector.
##' @param mu Numeric vector [-1, 1]. Mean parameter.
##' @param psi Numeric vector [0, Inf]. Shape parameter.
##' @param log Logical. If TRUE, return log density.
##' @return Numeric vector. Densities
##' @author Stephen R. Martin
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
##' @param p Numeric vector. Probabilities.
##' @inheritParams dtbeta mu psi
##' @return 
##' @author Stephen Martin
##' @rdname tbeta
qtbeta <- function(p, mu, psi) {
    mut <- (mu + 1) / 2
    a <- mut * psi
    b <- psi - a
    q.raw <- qbeta(p, a, b)
    qt <- q.raw * 2 - 1
    return(qt)
}
##' @title CDF for tbeta distribution.
##' @param q Numeric vector. Quantiles.
##' @inheritParams dtbeta mu psi
##' @return Numeric vector.
##' @author Stephen R. Martin
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

##' @title Transformed Beta distribution.
##' @param n Number of samples.
##' @inheritParams dtbeta mu psi
##' @return Numeric vector.
##' @author Stephen R. Martin
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
##' Implied PDF of difference in two tbeta random variates.
##'
##' @title Difference-in-two-tbeta-variates distribution.
##' @param x Numeric vector [-2, 2].
##' @param mus Numeric vector (length 2). Mu-parameters to two tbeta distributions.
##' @param psis Numeric vector (length 2). Psi-parameters to two tbeta distributions.
##' @param log Logical (Default: FALSE). Whether to return densities on log scale.
##' @return Numeric vector.
##' @author Stephen R. Martin
##' @rdname tbeta
dtbeta_diff <- function(x, mus, psis, log = FALSE) {
    out <- .tbeta_prior_diff_exact(x, mus, psis)
    if(log) {
        out <- log(out)
    }
    return(out)
}
