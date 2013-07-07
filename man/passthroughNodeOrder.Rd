\name{passthroughNodeOrder}
\alias{passthroughNodeOrder}
\title{Pass-through layout function for BioFabric}
\usage{
  passthroughNodeOrder(bfGraph)
}
\arguments{
  \item{inGraph}{The graph to order, in igraph format (see
  \code{\link{igraph}}). Node names MUST already be set
  using the vertex attribute "name".}
}
\value{
  The provided graph with no reordering.
}
\description{
  Does nothing
}
\details{
  Returns the provided graph with no changes. Use this
  function if you are providing a graph to bioFabric that
  is already node-ordered, so no reordering will occur.
}
\examples{
\dontrun{
library(igraph)
bfGraph = graph.tree(20, children=2, mode="out")
bfGraph <- autoNameForFabric(bfGraph)
bfGraph <- defaultNodeOrder(bfGraph)
bioFabric(bfGraph, orderFunc=passthroughNodeOrder)
}
}
\author{
  Bill Longabaugh <wlongabaugh@systemsbiology.org>
}

