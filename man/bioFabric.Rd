\name{bioFabric}
\alias{bioFabric}
\title{BioFabric for R}
\usage{
  bioFabric(inGraph, userOrder = NULL, orderFunc = NULL,
    shadowLinks = FALSE, dropNodeLabels = FALSE,
    dropZoneLabels = FALSE)
}
\arguments{
  \item{inGraph}{The graph to plot, in igraph format (see
  \code{\link{igraph}}). Node names should already be set
  using the vertex attribute "name", though the attribute
  "label" will work as well. If not set, the function will
  assign names of the form "BFn".}

  \item{userOrder}{(optional) An ordered list specifying
  the ordering of the node rows in the plot; nodes are
  identified using node name. All nodes must be present in
  the list. If this list is provided, the default node
  ordering routine will not be used.}

  \item{orderFunc}{(optional) A function that takes a graph
  to plot and returns the graph with the nodes reordered in
  the desired node order. Providing the defaultNodeOrder
  function for this parameter gives the default behavior.
  Note that if this parameter is provided, the userOrder
  parameter is ignored.}

  \item{shadowLinks}{(optional) Default value is FALSE. If
  set to TRUE, BioFabric will display shadow links. See the
  BioFabric paper or blog for more info.}

  \item{dropNodeLabels}{(optional) Default value is FALSE.
  If set to TRUE, no node labels are displayed. With large
  networks, some PDF readers cannot correctly handle tiny
  text, and labels may need to be dropped.}

  \item{dropZoneLabels}{(optional) Default value is FALSE.
  If set to TRUE, no node zone labels are displayed. With
  large networks, some PDF readers cannot correctly handle
  tiny text, and labels may need to be dropped.}
}
\description{
  An R implementation of the BioFabric network
  visualization tool
}
\details{
  Plots a network, provided in igraph format (see
  \code{\link{igraph}}), using BioFabric. This release
  supports user-specified node orders and the display of
  shadow links.

  Best results are obtained using a PDF display target, as
  shown in the example. BioFabric is very sensitive to
  correct relative scaling of linewidth, label size, and
  grid size. The PDF output has been targeted for correctly
  handling this scaling. PNG output should also be
  acceptable.

  SIZING: As the network grows, you will need to increase
  the size of the output to be able to get per-link
  resolution. PDF files that are not big enough (e.g. only
  10 inches wide) will create links below the hairline
  limit and may not print well. PNG outputs are the better
  bet to create low-resolution images of large networks in
  page-size dimensions. To be able to view an enormous
  network and zoom to see each link, create an enormous PDF
  file (e.g. 40, 80, 100+ inches wide!) or PNG file (e.g.
  30,000 pixels wide if the network has 5000 links).

  VIEWERS: PDF or PNG viewers that cannot do antialiasing
  of line art are a poor choice for BioFabric networks. The
  closely spaced parallel lines of a BioFabric plot MUST be
  antialiased to get acceptable results. Some specific
  viewer notes:

  Evince on Linux for PDF (Document Viewer 2.30.3 tested):
  Antialiased is fixed on, and visuals are good. Problems:
  Very tiny text below some size threshold explodes to a
  huge size. Cannot zoom above 400%, which is insufficient
  to explore network. No drag hand cursor to navigate,
  which is essential. Note that for Postscript output,
  Evince does not antialias, giving very poor results.

  Preview on Mac for PDF (Version 4.2 tested): Be sure that
  "Anti-alias text and line art" is checked on the PDF tab
  for Preferences, which gives fair visuals. Good maximum
  zoom level, and the Move cursor provides convenient
  draggable navigation.

  Adobe Reader on Mac for PDF (Version 9.5.5 tested): Be
  sure that "Smooth line art" is checked, and "Enhance thin
  lines" is NOT checked, on the Page Display Preferences.
  The Hand tool is available via Tools->Select & Zoom->Hand
  Tool, and the maximum 6400 percent zoom level is good for
  exploration.

  RStudio (0.97 tested): You won't be able to see much in
  the tiny Plot Frame. Due to a known RStudio issue, the
  Zoom window does not provide correct label sizing or line
  width sizing, since the plots are scaled based on the
  Plot Frame dimensions. Use the PDF or PNG output routes
  for now.

  General note on image viewers: If you have a plot with
  10,000 links, you are quickly heading towards PNG images
  that are over 32,000 pixels wide for high resolution.
  Lots of image viewers start to get cranky when viewing
  images that are more than 32,000 pixels wide.

  If you want to explore large networks with BioFabric,
  best to export the network as a tab-delimited SIF file
  and use the Java-based BioFabric application, which is
  designed for interactive large network visualization and
  exploration.
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
  Home page: \url{http://www.BioFabric.org}

  Blog: \url{http://biofabric.blogspot.com}

  Paper:
  \url{http://www.biomedcentral.com/1471-2105/13/275}
}

