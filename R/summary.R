##' @title Summary method for corCompare objects.
##' @param object corCompare object.
##' @param prob Numeric (0, 1). Amount of probability mass to contain in the credible interval.
##' @param ... Not used.
##' @return summary.corCompare object. List containing the summary tables (summary) and meta-data (meta).
##' @author Stephen R. Martin
##' @export
summary.corCompare <- function(object, prob = .95, ...) {
    dots <- list(...)
    meta <- object$meta
    digits <- dots$digits %IfNull% 3
    meta$digits <- digits
    map <- object$meta$group_spec$levels

    rho <- .summarize(object, "rho", prob)
    rho_diff <- .summarize(object, "rho_diff", prob)

    inds.rho <- .stan_to_inds(rownames(rho))
    rownames(rho) <- map[inds.rho[,1]]

    inds.rho_diff <- .stan_to_inds(rownames(rho_diff))
    rownames(rho_diff) <- paste0(map[inds.rho_diff[,1]],"-",map[inds.rho_diff[,2]])

    out <- list(summary = list(rho = rho,
                               rho_diff = rho_diff),
                meta = meta)
    class(out) <- "summary.corCompare"
    return(out)
}
##' @title Print method for summary.corCompare
##' @param x summary.corCompare object.
##' @param ... Not used.
##' @return x (Invisibly).
##' @author Stephen R. Martin
##' @export
print.summary.corCompare <- function(x, ...) {
    dots <- list(...)
    digits <- dots$digits %IfNull% x$meta$digits

    .sep()
    .newline()
    cat("Correlations per group:")
    .newline()

    print(x$summary$rho, digits = digits)

    .sep()
    .newline()
    cat("All differences in correlations:")
    .newline()
    # Lower-diag only
    G <- x$meta$G
    inds <- which(lower.tri(diag(1, G, G)))
    print(x$summary$rho_diff[inds,], digits = digits)

    invisible(x)
}

.newline <- function(n = 1) {
    cat(rep("\n", n))
}

.sep <- function(sep = "---") {
    cat(sep)
}

.summarize <- function(object, pars, prob) {
    samps <- as.matrix(object$fit, pars)
    probs <- .prob_to_probs(prob)
    out <- apply(samps, 2, function(p) {
        m <- mean(p)
        s <- sd(p)
        mdn <- quantile(p, .5)
        cri <- quantile(p, probs)
        out <- c(m, s, mdn, cri)
        names(out) <- c("mean", "sd", "mdn", paste0("Q", probs * 100))
        out
    })
    out <- t(out)
    out <- cbind(out, rstan::summary(object$fit, rownames(out))$summary[,"Rhat"])

    colnames(out)[ncol(out)] <- "Rhat"

    # Add BFs
    if(pars == "rho") {
        out <- cbind(out, BF01 = .bf_rho(object))
    }
    if(pars == "rho_diff") {
        out <- cbind(out, BF01 = .bf_rho_diff(object))
    }

    return(out)
}

.prob_to_probs <- function(x) {
    lower <- (1 - x) / 2
    upper <- 1 - lower
    return(c(lower, upper))
}