
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
