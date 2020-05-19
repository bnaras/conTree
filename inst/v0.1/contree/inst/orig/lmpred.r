lmpred <- function(x, lmm) {
  lmm$coefficients[1] + x %*% lmm$coefficients[2:length(lmm$coefficients)]
}
