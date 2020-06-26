##' Default setting helper.
##'
##' If x is null, replace with expr.
##' @title Helper for defaulting NULL values.
##' @param x Object.
##' @param expr Expression.
##' @return 
##' @author Stephen Martin
##' @keywords internal
`%IfNull%` <- function(x, expr) {
    if(is.null(x)) {
        return(eval(expr))
    } else {
       return(x) 
    }
}
