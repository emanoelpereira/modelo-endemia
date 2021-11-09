#' end.of.epiweek
#'
#' @param x a date
#' @param end weekday that ends the week: 0 is Sunday (default is 6, Saturday)
#'
#' @return date of the last day of the input date's epiweek
#'
#' @examples
end.of.epiweek <- function(x, end = 6) {
    offset <- (end - 4) %% 7
    num.x <- as.numeric(x)
    return(x - (num.x %% 7) + offset + ifelse(num.x %% 7 > offset, 7, 0))
}
