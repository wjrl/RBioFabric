\name{bioFabric}
\alias{bioFabric}
\title{BioFabric for R}
\usage{
  bioFabric(inGraph)
}
\arguments{
  \item{inGraph}{The graph to plot, in igraph format (see
  \code{\link{igraph}}). Node names should already be set
  using the vertex attribute "label". If not set, the
  function will assign number names.}
}
\description{
  An R implementation of the BioFabric network
  visualization tool
}
\details{
  Plots a network, provided in igraph format (see
  \code{\link{igraph}}), using BioFabric. This is a
  "pre-alpha" first release, and only handles the default
  layout algorithm. NOTE 1: Best results are obtained using
  a PDF display target, as shown in the examples. BioFabric
  is very sensitive to correct scaling, and the PDF output
  has been checked as being decent. NOTE 2: If you are an
  RStudio user (yay!), you won't be able to see much in the
  tiny Plot Frame. BUT, due to a known RStudio issue, the
  Zoom window will NOT give you correct label sizing or
  line width sizing. Use the PDF output!
}
\examples{
\dontrun{
# Gotta have igraph!

library(igraph)

# Generate graph

bfGraph = barabasi.game(100, m=6, directed=FALSE)

# Plot it up! For best results, make the PDF in the same
# aspect ratio as the network, though a little extra height
# covers the top labels

height <- vcount(bfGraph)
width <- ecount(bfGraph)
aspect <- height / width;
plotWidth <- 10.0
plotHeight <- plotWidth * (aspect * 1.2)
pdf("myBioFabricOutput.pdf", width=plotWidth, height=plotHeight)
bioFabric(bfGraph)
dev.off()
}
}
\author{
  Bill Longabaugh <wlongabaugh@systemsbiology.org>
}
\references{
  \url{http://www.BioFabric.org}

  \url{http://www.biomedcentral.com/1471-2105/13/275}
}

