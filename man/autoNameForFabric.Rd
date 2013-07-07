\name{autoNameForFabric}
\alias{autoNameForFabric}
\title{Node naming utility function for BioFabric}
\usage{
  autoNameForFabric(inGraph)
}
\arguments{
  \item{inGraph}{The graph to name, in igraph format (see
  \code{\link{igraph}}).}
}
\value{
  The provided graph with nodes named
}
\description{
  BioFabric requires nodes have names. This will do it!
}
\details{
  BioFabric requires that all nodes have names. If network
  is named, there is no change. If nodes are unnamed but
  have the attribute "label", that attribute is copied into
  the "name" attribute. Otherwise, names of the form "BFn"
  are formed, where n is the index of the node.
}
\examples{
\dontrun{

library(igraph)

bfGraph = graph.tree(20, children=2, mode="out")
bfGraph <- autoNameForFabric(bfGraph)
}
}
\author{
  Bill Longabaugh <wlongabaugh@systemsbiology.org>
}

