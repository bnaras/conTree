datgen <-
  function(p = 10, no = 25000, scl = 0.5, sclf = 0.75, er = 2,
           dist = "logistic", smax = 50, mseed = 2, xseed = 32, yseed = 97) {
    svs <- .Random.seed
    set.seed(xseed)
    x <- matrix(rnorm(p * no), ncol = p)
    set.seed(mseed)
    b <- rnorm(p)
    b <- b / sqrt(sum(b**2))
    z <- er * runif(p)
    f <- rep(0, no)
    for (j in 1:p) {
      u <- sign(x[, j]) * abs(x[, j])**z[j]
      v <- sqrt(var(u))
      f <- f + b[j] * u / v
    }
    z <- er * runif(p)
    s <- rep(0, no)
    b <- rnorm(p)
    b <- b * sclf / sqrt(sum(b**2))
    for (j in 1:p) {
      u <- sign(x[, j]) * abs(x[, j])**z[j]
      v <- sqrt(var(u))
      s <- s + b[j] * u / v
    }
    s <- pmin(smax, scl * exp(s))
    set.seed(yseed)
    y <- lgen(no, f, s, dist)
    set.seed(svs)
    invisible(list(x = x, y = y, f = f, s = s))
  }
