
##' Fits the correlation-compare model.
##'
##' The formula is multipart: y ~ x | groups.
##' I.e., outcome ~ predictor | between + within + between2
##'
##' The model effectively estimates a covariance/correlation between outcome and predictor, for each possible combination of the specified grouping variables.
##' @title Fit cor_compare model.
##' @param formula Formula. See details.
##' @param data Data.frame.
##' @param ... Options to pass onto sampling.
##' @return Stanfit.
##' @author Stephen R. Martin
##' @import Formula
##' @importFrom rstan sampling
##' @export
cor_compare <- function(formula, data, ...) {
    dots <- list(...)
    formula <- as.Formula(formula)
    d <- .parse_formula(formula, data)
    meta <- d$meta

    stan_data <- d$data
    stan_data$prior_only <- dots$prior_only %IfNull% 0
    stan_data$tbeta_mu <- dots$tbeta_mu %IfNull% rep(0, d$meta$G)
    stan_data$tbeta_psi <- dots$tbeta_psi %IfNull% rep(2, d$meta$G)
    if(length(stan_data$tbeta_mu) == 1 & d$meta$G > 1) {
        stan_data$tbeta_mu <- rep(tbeta_mu, d$meta$G)
    }
    if(length(stan_data$tbeta_psi) == 1 & d$meta$G > 1) {
        stan_data$tbeta_psi <- rep(tbeta_psi, d$meta$G)
    }

    stan_args <- list(
        iter = dots$iter %IfNull% 2000,
        cores = dots$cores %IfNull% parallel::detectCores()
    )
    stan_args$control <- dots$control %IfNull% list(adapt_delta = .95)
    stan_args$control$adapt_delta <- stan_args$control$adapt_delta %IfNull% .95
    stan_args$object <- stanmodels$mixed

    meta$stan_args <- stan_args
    stan_args$data <- stan_data

    sOut <- do.call(sampling, stan_args)
    out <- list(fit = sOut,
                meta = meta)
    class(out) <- "corCompare"
    return(out)
}

.parse_formula <- function(formula, data) {
    formula <- Formula::as.Formula(formula)
    mf <- model.frame(formula, data, na.action = na.pass)

    # Remove missings.
    mf.c <- mf[complete.cases(mf), ]
    n_missings <- nrow(mf) - nrow(mf.c)
    if(n_missings > 1) {
        warning("Removed ", n_missings, "incomplete rows.")
    }
    mf <- mf.c

    # Check formula structure
    fLength <- length(formula)
    if(!all(fLength == c(1, 2))) {
        stop("Formula should contain one LHS, and two RHS terms [separated by a |].")
    }
    y <- as.vector(model.frame(formula, mf, lhs = 1, rhs = 0)[,1])
    x <- as.vector(model.frame(formula, mf, rhs = 1, lhs = 0)[,1])
    groups <- model.frame(formula, mf, rhs = 2, lhs = 0)
    # Create group-numerics
    group_char <- interaction(groups)
    group_numeric <- as.numeric(group_char)
    group_spec <- list(data = groups,
                       char = group_char,
                       numeric = group_numeric,
                       levels = levels(group_char))

    out <- list(meta = list(),
                data = list())
    out$meta$group_spec <- group_spec
    out$meta$n_missings <- n_missings
    out$meta$formula <- formula
    out$meta$G <- max(group_numeric)
    out$data <- list(
        N = nrow(mf),
        group = group_numeric,
        y = y,
        x = x
    )
    return(out)
}

