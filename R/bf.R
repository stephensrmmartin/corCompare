##' @title Estimates (logspline) density of vector.
##' @param x Numeric vector.
##' @param lb Numeric. Lower bound.
##' @param ub Numeric. Upper bound.
##' @return Logspline object.
##' @author Stephen R. Martin
##' @import logspline
##' @keywords internal
.estimate_dist <- function(x, lb = -1, ub = 1) {
    ls <- logspline(x, lb, ub)
    return(ls)
}

##' @title Density wrapper.
##' @param x Numeric vector.
##' @param dist Density function or logspline object.
##' @param ... Arguments passed onto density function.
##' @return Numeric.
##' @author Stephen R. Martin
##' @keywords internal
.estimate_dens <- function(x, dist, ...) {
    dots <- list(...)
    if(class(dist) == "logspline") {
        do.call(dlogspline, list(x, dist))
    } else {
        args <- c(list(x), dots)
        do.call(dist, args)
    }
}

.tbeta_prior_diff <- function(N = 100000, mus, psis, dens = TRUE) {
    x1 <- rtbeta(N, mus[1], psis[1])
    x2 <- rtbeta(N, mus[2], psis[2])
    x <- x1 - x2
    if(dens) {
        return(.estimate_dist(x, -2, 2))
    } else {
        return(x)
    }
}

.tbeta_prior_diff_exact <- function(diff, mus, psis) {
    # [diff, r2] = f(r1, r2) = [r1 - r2, r2]
    # f_inv(diff, r2) = [diff + r2, r2]
    # p(diff, r2) = p(diff + r2, r2)
    # p(diff) = int p(diff + r2)p(r2)dr2
    f_r1 <- function(r2, diff, mu1, psi1) {
        dtbeta(r2 + diff, mu1, psi1)
    }
    f_r2 <- function(r2, mu2, psi2) {
        dtbeta(r2, mu2, psi2)
    }
    f <- function(r2, diff, mus, psis) {
        f_r1(r2, diff, mus[1], psis[1]) * f_r2(r2, mus[2], psis[2])
    }

    # Integrate
    ## out <- integrate(f, -1, 1, diff = diff, mus = mus, psis = psis)
    out <- sapply(diff, function(x) {
        integrate(f, -1, 1, diff = x, mus = mus, psis = psis)$value
    })

    return(out)
}

.bf_rho <- function(object) {
    pars <- "rho"
    mus <- object$meta$stan_args$tbeta_mu
    psis <- object$meta$stan_args$tbeta_psi
    denoms <- dtbeta(0, mus, psis)
    samps <- as.matrix(object$fit, pars = pars)
    posts <- apply(samps, 2, .estimate_dist)
    nums <- sapply(seq_len(length(posts)), function(p) {
        .estimate_dens(0, posts[[p]])
    })
    return(nums / denoms)
}


.bf_rho_diff <- function(object) {
    pars <- "rho_diff"
    mus <- object$meta$stan_args$tbeta_mu
    psis <- object$meta$stan_args$tbeta_psi
    samps <- as.matrix(object$fit, pars = pars)

    # Don't get diagonals
    inds <- .stan_to_inds(colnames(samps))
    inds_diag <- which(apply(inds, 1, function(r) {
        ref <- r[1]
        eq <- r == ref
        all(eq)
    }))
    inds <- inds[-inds_diag,]
    samps <- samps[, -inds_diag]

    posts <- apply(samps, 2, .estimate_dist, lb = -2, ub = 2)
    nums <- sapply(seq_len(length(posts)), function(p) {
        .estimate_dens(0, posts[[p]])
    })

    denoms <- apply(seq_len(nrow(inds)), 1, function(r) {
        g <- inds[r,]
        dtbeta_diff(0, mus[g], psis[g])
    })

    return(nums / denoms)
}

.stan_to_inds <- function(x) {

    rex.inner <- r"(.*\[([[:digit:]](?:,[[:digit:]]).*)\])"
    inner <- gsub(rex.inner, "\\1", x)
    split <- lapply(strsplit(inner, ","), as.numeric)
    inds <- do.call(rbind, split)
    return(inds)
    
}

