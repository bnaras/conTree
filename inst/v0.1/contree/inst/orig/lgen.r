lgen <- function(n, f, s, dist) {
  if (dist == "logistic") {
    r <- runif(n)
    y <- f + s * log(r / (1 - r))
  }
  else if (dist == "normal") {
    y <- rnorm(n, f, 1.63107 * s)
  }
  else if (dist == "laplace") {
    y <- laplace(n, f, 1.585 * s)
  }
  else if (dist == "slash") {
    y <- f + 0.7522519 * s * rnorm(n) / runif(n)
  }
  else {
    stop("unsupported error distribution")
  }
  invisible(y)
}
