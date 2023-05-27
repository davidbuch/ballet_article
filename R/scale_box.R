scale_box <- function(x, pad = TRUE){
  x <- as.matrix(x)
  nc <- ncol(x)
  if (pad) {
    f_low <- function(v) {
      v <- v[!is.na(v)]
      min(v) - sqrt(sum(v^2)/max(1, length(v) - 1L))
    }
    f_range <- function(v) {
      v <- v[!is.na(v)]
      max(v) - min(v) + 2 * sqrt(sum(v^2)/max(1, length(v) - 1L))
    }
  } else {
    f_low <- function(v) {
      v <- v[!is.na(v)]
      min(v)
    }
    f_range <- function(v) {
      v <- v[!is.na(v)]
      max(v) - min(v)
    }
  }
  
  lower <- apply(x, 2L, f_low)
  range <- apply(x, 2L, f_range)
  
  x <- sweep(x, 2L, lower, '-', check.margin = FALSE)
  x <- sweep(x, 2L, range, '/', check.margin = FALSE)
  
  attr(x, "scaled:lower") <- lower
  attr(x, "scaled:range") <- range
  
  x
}