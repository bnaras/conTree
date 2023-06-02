library('conTree')
data(census, package = "conTree")
attach(census)
dx <- 1:10000
dxt <- 10001:16281
tree <- contrast(xt[dx,], yt[dx], gblt[dx], type = 'prob')
## Compare
nodesum(tree, xt[dxt,], yt[dxt], gblt[dxt])

## with 
dfun <- function(y, z, w) {
  w  <- w / sum(w)
  abs(sum(w * (y - z)))
}

tree2 <- contrast(xt[dx,], yt[dx], gblt[dx], type = dfun)
nodesum(tree2, xt[dxt,], yt[dxt], gblt[dxt])
