% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/contrast.R
\name{prune.seq}
\alias{prune.seq}
\title{Show all possible pruned subtrees}
\usage{
prune.seq(tree, eps = 0.01, plot.it = TRUE)
}
\arguments{
\item{tree}{a tree model object output from contrast}

\item{eps}{small increment defining grid of threshold values}

\item{plot.it}{a logical flag indicating plot/don't plot of number of nodes versus threshold value for all pruned subtrees, default \code{TRUE}}
}
\value{
a named list of two items:
\itemize{
\item \code{thr} a set of threshold values that sequentially reduce tree size
\item \code{nodes} the corresponding tree sizes (number of terminal nodes)
}
}
\description{
Show all possible pruned subtrees
}
