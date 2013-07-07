\name{defaultNodeOrder}
\alias{defaultNodeOrder}
\title{Default layout function for BioFabric}
\usage{
  defaultNodeOrder(bfGraph, firstLabel = NULL)
}
\arguments{
  \item{inGraph}{The graph to order, in igraph format (see
  \code{\link{igraph}}). Node names MUST already be set
  using the vertex attribute "name".}

  \item{firstLabel}{(optional) The label of the node to
  start at.}
}
\value{
  The provided graph with the nodes reordered per the
  BioFabric default layout algorithm.
}
\description{
  Does breadth-first layout of a network for BioFabric
}
\details{
  Default node layout based on breadth-first search from
  highest degree node, searching in order of decreasing
  degree. Takes a network provided in igraph format (see
  \code{\link{igraph}}) and returns the same network, with
  the nodes reordered to the default breadth-first layout
  order
}
\examples{
\dontrun{
# Gotta have igraph!

library(igraph)

# defaultNodeOrder requires we name the graph, so we
# do that first. Typically, there is no reason to provide
# this function to bioFabric, since it uses it internally.
# However, we can gain extra functionality by providing the
# firstLabel argument, which starts the search at the specified
# node instead od the highest degree node. So wrapping the
# defaultNodeOrder as shown allows us to start at the top of
# the given tree (degree = 2), which is lower degree than
# its child nodes.

bfGraph = graph.tree(20, children=2, mode="out")
bfGraph <- autoNameForFabric(bfGraph)
startAtBF1 <- function(bfGraph) {
  return (defaultNodeOrder(bfGraph, firstLabel=V(bfGraph)[1]$name))
}
bioFabric(bfGraph, orderFunc=startAtBF1)
}
}
\author{
  Bill Longabaugh <wlongabaugh@systemsbiology.org>
}

