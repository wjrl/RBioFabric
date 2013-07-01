# =================================================================================
#
# Released under "The MIT License":
# 
# Copyright (c) 2013 William J.R. Longabaugh
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# ================================================================================  

#'@title BioFabric for R
#'
#'@description
#'An R implementation of the BioFabric network visualization tool
#'
#'@details
#'Plots a network, provided in igraph format (see \code{\link{graph}}),
#'using BioFabric. This is a "pre-alpha" first release, and only handles
#'the default layout algorithm. NOTE 1: Best results are obtained using
#'a PDF display target, as shown in the examples. BioFabric is very 
#'sensitive to correct scaling, and the PDF output has been checked
#'as being decent. NOTE 2: If you are an RStudio user (yay!), you
#'won't be able to see much in the tiny Plot Frame. BUT, due to a known
#'RStudio issue, the Zoom window will NOT give you correct label
#'sizing or line width sizing. Use the PDF output!
#'
#'@param inGraph The graph to plot, in igraph format (see 
#'\code{\link{graph}}). Node names should already be set using the vertex
#'attribute "label". If not set, the function will assign number names.
#'@author Bill Longabaugh <wlongabaugh@@systemsbiology.org>
#'@export
#'@examples
#'
#' \dontrun{
#' # Gotta have igraph!
#' 
#' library(igraph)
#' 
#' # Generate graph
#' 
#' bfGraph = barabasi.game(100, m=6, directed=FALSE)
#'
#' # Plot it up! For best results, make the PDF in the same
#' # aspect ratio as the network, though a little extra height
#' # covers the top labels
#' 
#' height <- vcount(bfGraph)
#' width <- ecount(bfGraph)
#' aspect <- height / width;
#' plotWidth <- 10.0 
#' plotHeight <- plotWidth * (aspect 1.2)
#' pdf("myBioFabricOutput.pdf", width=plotWidth, height=plotHeight)
#' bioFabric(bfGraph)
#' dev.off()
#' }
#' 
#'@references \url{http://www.BioFabric.org} 
#'            \url{http://www.biomedcentral.com/1471-2105/13/275}
#'

bioFabric <- function(inGraph) {

  # ======================================================
  # Utility function to create the color cycle used in 
  # BioFabric. Nodes are drawn lighter than the base color, 
  # while edges are drawn darker than the base. Thus, we
  # take the multiplier and hit the rgb values with it;
  # capped at 255. This is also where we install the order
  # that will be used in the color cycle
  # ======================================================
  
  colCyc <- function(cList, order, mult) {
    oLen <- length(order);
    retval <- vector(length=oLen)
    for (i in 1:oLen) {
      taggedCol <- cList[order[[i]]]
      transCol <- col2rgb(taggedCol)
      myR <- min(255, transCol[1] * mult);
      myG <- min(255, transCol[2] * mult);
      myB <- min(255, transCol[3] * mult);
      retval[i] <- rgb(myR, myG, myB, maxColorValue = 255);
    }
    return (retval)
  }
  
  # ======================================================
  # Utility functions to create the separate coordinate
  # vectors for the plot: duplicate or alternate
  # ======================================================
  
  alter <- function(origA, origB) {
    oLen <- length(origA);
    retval <- vector()
    for (i in seq(1, oLen)) {
      retval[(2 * i) - 1] <- origA[[i]]
      retval[(2 * i)] <- origB[[i]]
    }
    return (retval)
  }
  
  duper <- function(orig) {
    return (alter(orig, orig))
  }
  
  # ======================================================
  # Define the colors used. These colors were originally
  # derived for use in BioTapestry (www.BioTapestry.org)
  # to create a set of colors that could be distinguished
  # at link crossings.
  # ======================================================
  
  baseColorNames <- c("EX-cyan",
                      "EX-dark-cyan",
                      "EX-yellow-orange",
                      "EX-pale-green",
                      "EX-dark-green",
                      "EX-pale-red-orange",
                      "EX-yellow-green",
                      "EX-yellow",
                      "EX-dark-gray-purple",
                      "EX-pale-magenta",
                      "EX-pale-purple",
                      "EX-purple",
                      "EX-dark-red",
                      "EX-red",
                      "EX-pale-yellow-green",
                      "EX-dark-purple",
                      "EX-pale-cyan",
                      "EX-pure-blue",
                      "EX-dark-yellow-green",
                      "EX-magenta",
                      "EX-dark-tan",
                      "EX-pale-blue",
                      "EX-orange",
                      "EX-medium-magenta",
                      "EX-blue-magenta",
                      "EX-green",
                      "EX-dark-magenta",
                      "EX-pale-blue-magenta",
                      "EX-pale-yellow orange",
                      "EX-dark-orange",
                      "EX-blue",
                      "EX-pale-red")
  
  baseColorColors <- c(rgb(0, 255, 255, maxColorValue = 255),
                       rgb(0, 100, 128, maxColorValue = 255),
                       rgb(255, 153, 0, maxColorValue = 255),
                       rgb(133, 205, 102, maxColorValue = 255),
                       rgb(39, 128, 0, maxColorValue = 255),
                       rgb(230, 156, 138, maxColorValue = 255),
                       rgb(154, 255, 0, maxColorValue = 255),
                       rgb(255, 203, 0, maxColorValue = 255),
                       rgb(0, 25, 128, maxColorValue = 255),
                       rgb(212, 138, 230, maxColorValue = 255),
                       rgb(149, 165, 230, maxColorValue = 255),
                       rgb(102, 51, 255, maxColorValue = 255),
                       rgb(140, 56, 56, maxColorValue = 255),
                       rgb(255, 0, 0, maxColorValue = 255),
                       rgb(222, 230, 138, maxColorValue = 255),
                       rgb(77, 56, 140, maxColorValue = 255),
                       rgb(138, 230, 181, maxColorValue = 255),
                       rgb(0, 0, 255, maxColorValue = 255),
                       rgb(114, 128, 0, maxColorValue = 255),
                       rgb(255, 0, 255, maxColorValue = 255),
                       rgb(166, 133, 83, maxColorValue = 255),
                       rgb(102, 183, 205, maxColorValue = 255),
                       rgb(255, 103, 0, maxColorValue = 255),
                       rgb(166, 83, 166, maxColorValue = 255),
                       rgb(155, 0, 255, maxColorValue = 255),
                       rgb(0, 255, 0, maxColorValue = 255),
                       rgb(102, 0, 128, maxColorValue = 255),
                       rgb(146, 102, 205, maxColorValue = 255),
                       rgb(205, 175, 102, maxColorValue = 255),
                       rgb(128, 92, 0, maxColorValue = 255),
                       rgb(0, 152, 255, maxColorValue = 255),
                       rgb(205, 102, 153, maxColorValue = 255))
  
  names(baseColorColors) <- baseColorNames
  
  # ======================================================
  # Define the color order used in the color cycle. This
  # order was originally derived for use in BioTapestry 
  # (www.BioTapestry.org) such that adjacent colors would
  # not be "too close", making it hard to distinguish when
  # they were coloring adjacent parallel links.
  # ======================================================
  
  tagOrder <- c("EX-blue",
                "EX-orange",
                "EX-dark-cyan",
                "EX-red",
                "EX-dark-orange",
                "EX-dark-gray-purple",
                "EX-cyan",
                "EX-yellow-orange",
                "EX-pure-blue",
                "EX-dark-yellow-green",
                "EX-dark-magenta",
                "EX-dark-green",
                "EX-blue-magenta",
                "EX-yellow-green",
                "EX-magenta",
                "EX-green",
                "EX-yellow",
                "EX-purple",
                "EX-dark-purple",
                "EX-dark-red",
                "EX-pale-green",
                "EX-pale-blue",
                "EX-dark-tan",
                "EX-pale-blue-magenta",
                "EX-pale-yellow orange",
                "EX-medium-magenta",
                "EX-pale-red",
                "EX-pale-cyan",
                "EX-pale-yellow-green",
                "EX-pale-purple",
                "EX-pale-magenta",
                "EX-pale-red-orange")
  
  # ======================================================
  # Define the node, link, and glyph color vectors. Again,
  # nodes are made lighter, links are make darker. Glyphs
  # have the same color as the links, but we need two in
  # a row to do the start and end glyphs.
  # ======================================================
  
  mcolN <- colCyc(baseColorColors, tagOrder, 1.4)
  mcolL <- colCyc(baseColorColors, tagOrder, 1.0 / 1.4)
  mcolG <- duper(mcolL)
  
  # ======================================================
  # Time to go to work on the graph. The graph HAS to 
  # have labeled vertices! If the igraph does not have
  # already, we glue them on. Also get the degree and
  # the vertex count!
  # ======================================================
  
  numV <- vcount(inGraph)
  vlabels <- get.vertex.attribute(inGraph, "label")
  if (is.null(vlabels)) {
    vlabels <- paste(1:vcount(inGraph))
    bfGraph <- set.vertex.attribute(inGraph, "label", value=vlabels);
  } else {
    bfGraph <- inGraph
  }
  degrees <- degree(bfGraph)
  isDirected <-is.directed(bfGraph)
  
  # ======================================================
  # Need to set up the breadth-first search of the network,
  # in order to get the default layout. Since graph.bfs searches
  # nodes using the existing node order in the graph, we will need
  # to permute the node order. The order we want is by decreasing
  # node degree, with ties broken using node label lexicographic
  # order. So we build a data frame with the required fields,
  # sort it, and create a permutation vector we can send into
  # permute.vertices. Then run the search!
  # ======================================================
  
  # build the frame, sort it:
  nodesToOrder <- data.frame(degrees, vlabels)
  ordered <- order(-nodesToOrder$degrees, nodesToOrder$vlabels)
  
  # Generate the permutation vector:
  nodePerm <- vector(length=numV)
  count <- 1
  for (i in ordered) {
    nodePerm[[i]] <- count
    count <- count + 1
  }
  
  # Permute the order, run the breath-first search. Note
  # that we ignore direction when doing this search on
  # directed graphs
  bfGraphPass1 <- permute.vertices(bfGraph, nodePerm)
  bfGraphBFS <- graph.bfs(bfGraphPass1, 1, neimode="all")
  
  # ======================================================
  # Use the results of the search to permute the graph
  # into the final order we need to lay it out.
  # ======================================================
  
  # Generate the final permutation vector that will
  # create the node order we need for the default
  # BioFabric layout:
  nodePerm2 <- vector(length=numV)
  count <- 1
  for (i in bfGraphBFS$order) {
    nodePerm2[[i]] <- count
    count <- count + 1
  }
  
  # ======================================================
  # Permute the order, thus creating the final
  # graph we will layout. Extract the vector of
  # labels, which we use to tag the node lines:
  # ======================================================
  
  bfGraph <- permute.vertices(bfGraphPass1, nodePerm2)
  bfGraphLabels <- get.vertex.attribute(bfGraph, "label")
  
  # ======================================================
  # Yay! Nodes are sorted. Now we need to sort the edges.
  # Are edges in the graph always organized so that the
  # lower-numbered node is in column 1? Maybe for undirected
  # graphs, but not for directed graphs. So make
  # sure it is the case, and build a dataFrame that we
  # can sort
  # ======================================================
  
  edgelist <- get.edgelist(bfGraph)
  numE <- ecount(bfGraph)
  
  minNodes <- vector(length=numE)
  maxNodes <- vector(length=numE)
  targNode <- vector(length=numE)
  for (i in 1:numE) {
    node1 <- edgelist[i,1]
    node2 <- edgelist[i,2]
    minNodes[[i]] <- min(node1, node2)
    maxNodes[[i]] <- max(node1, node2)
    targNode[[i]] <- node2
  }
  edgesForSort <- data.frame(top=minNodes, bot=maxNodes, targ=targNode)
  
  # ======================================================
  # Sort the edges so that edges with the topmost node go
  # first, with ties broken so the shorter link comes before
  # the longer link. Result: edge wedges!
  # ======================================================
  
  edgeOrder <- order(edgesForSort$top, edgesForSort$bot)
  orderedEdges <- edgesForSort[edgeOrder,]
  
  # ======================================================
  # We want the node lines to only be as long as they have
  # to be. How to do that? We want to create two vectors.
  # One lists the min column for each node, the other the
  # max column. To do this, we augment our edgeOrder frame
  # with an index to create a new combined frame. Then, for
  # each node, we extract all the rows in the frame that 
  # contain the node and create a vector of the edge indices
  # (which are the same as BioFabric edge column assignments).
  # Then we just need to add the min and the max of that
  # vector to the nodeMin and nodeMax we are building.
  #
  # Second, we need to find the "node zone" for each node,
  # which is the contiguous run of edges at the right end
  # of the node line. So we sort the vector of edge indices,
  # and walk through it backwards, stopping when we no
  # longer have a contiguous run.
  # ======================================================
  
  indexFrame <- data.frame(indx=1:numE)
  frameForNodeMinMax <- cbind(orderedEdges, indexFrame)
  
  nodeMin <- vector(length=numV)
  nodeMax <- vector(length=numV)
  zoneBoundsMin <- vector(length=numV)
  zoneBoundsMax <- vector(length=numV)
  zoneMid <- vector(length=numV)
  gotSingleton <- FALSE;
  for (i in 1:numV) {
    # Note the use of |, not ||, to OR in a vector:
    match <- frameForNodeMinMax[((frameForNodeMinMax$top == i) | (frameForNodeMinMax$bot == i)),]$indx
    # Singleton nodes will not have any match. But we still
    # want to show them out on the far right edge
    if (length(match) == 0) {
      match <- c(numE + 1)
      gotSingleton <- TRUE;
    }  
    nodeMin[[i]] <- min(match)
    maxCol <- max(match)
    nodeMax[[i]] <- maxCol
    # Node zone bounds calculated here:
    zoneBoundsMin[[i]] <- maxCol
    zoneBoundsMax[[i]] <- maxCol
    zoneMid[[i]] <- maxCol
    smatch <- sort(match)  
    matLen <- length(smatch)
    if (matLen > 1) {      
      for (j in matLen:2) {
        # Decrement the min if it is contiguous:
        if ((smatch[[j]] - 1) == smatch[[j - 1]]) {
          edgeAtCol <- orderedEdges[smatch[[j - 1]],]
          # if still contig, but now not at top, we are done too:
          if (is.na(edgeAtCol$top) || (edgeAtCol$top != i)) {
            break;
          }       
          zoneBoundsMin[[i]] <- smatch[[j]]
          zoneMid[[i]] <- (zoneBoundsMax[[i]] + smatch[[j - 1]]) / 2.0
        } else {
          break;
        }
      }
    }
  }
  
  # ======================================================
  # Graph width slightly wider if we need to show a column
  # of singletons
  # ======================================================
  
  graphWidth <- numE
  if (gotSingleton) {
    graphWidth <- graphWidth + 1;
  }
  
  # ======================================================
  # We do not label all the node zones. The zone has to exist
  # along the top edge of the network. 
  # ======================================================
  
  for (i in 1:numV) {
    edgeAtMax <- orderedEdges[nodeMax[[i]],]
    if (is.na(edgeAtMax$top) || (edgeAtMax$top != i)) {
      zoneBoundsMin[[i]] <- NA
      zoneBoundsMax[[i]] <- NA
      zoneMid[[i]] <- NA
    }
  }
  
  # ======================================================
  # Time to start drawing. Note, these "magic numbers" are
  # derived from tests of the best rendering behavior for
  # the BioFabric metaphor on Java2D. 
  # ======================================================
  
  grid <- 18
  halfBox <- 5
  strokeSize <- 3
  
  # ======================================================
  # These font-specific magic numbers come from trial 
  # and error in R
  # ======================================================
  
  labelRot <- 0.075
  nodeLabelFrac <- 0.45
  nodeLabelOffset <- 0.4
  nodeZoneLabelFrac <- 0.75
  maxTextHeightGrid <- 3
  minimumTextCexWithoutFailure <- 0.05 ## Labels too small? BOOM!
  
  # ======================================================
  # This number is the most nutso. What line width, in inches
  # corresponds to lwd = 1? This is the number I found in 
  # a 12/17/2012 blog posting by Xianjun Dong at the URL
  # http://www.r-bloggers.com/line-width-in-r-and-in-illustrator/
  # ======================================================
  
  baseLineWidth <- (1 / 96)
  
  # ======================================================
  # Set up the empty plot space. Important: the aspect
  # ratio is forced to 1! This is because BioFabric REQUIRES
  # that node and edge lines are equally spaced! 
  # ======================================================
  
  par(mai=c(0.25, 0.25, 0.25, 0.25))
  plot(0, 0, type="n", yaxt="n", xaxt="n", ylab="", xlab="",
       xlim=c(0.0, graphWidth * grid), ylim=c(0.0, numV * grid), asp=1)
  
  # ======================================================
  # Figure out the size, in inches, of our basic grid unit.
  # This is CRUCIAL, since we need to get our line width
  # in correct proportion to the gridding. Also need to
  # be able to size fonts in correct proportion to grid.
  # ======================================================
  
  gridInH <- par("pin")[1] / graphWidth
  gridInV <- par("pin")[2] / numV
  gridInMin <- min(gridInH, gridInV)
  inPerGrid <- gridInMin / grid
  lineIn <- strokeSize * inPerGrid
  cFac <- lineIn / baseLineWidth
  
  # ======================================================
  # Draw the node lines first. Note they only extend from
  # the min to the max. Also note we draw the first node
  # at the maximum Y, going top down, so we are working 
  # backwards in Y coordinates
  # ======================================================
  
  x0 <- nodeMin * grid
  y0 <- (numV * grid) - (1:numV * grid)
  x1 <- nodeMax * grid
  y1 <- (numV * grid) - (1:numV * grid)
  # We still want a tiny node line for singleton nodes
  for (i in 1:numV) {
    if (x0[[i]] == x1[[i]]) {
      x0[[i]] <- x0[[i]] - halfBox + 1;
      x1[[i]] <- x1[[i]] + halfBox - 1;
    }
  }
  segments(x0, y0, x1, y1, col = mcolN, lwd = cFac)
  
  # ======================================================
  # Next draw the edges over the node lines. 
  # ======================================================
  
  x0 <- 1:numE * grid
  x1 <- 1:numE * grid
  y0 <- (numV * grid) - (orderedEdges$top * grid)
  y1 <- (numV * grid) - (orderedEdges$bot * grid)
  segments(x0, y0, x1, y1, col = mcolL, lwd = cFac)
  
  # ======================================================
  # Next draw the endpoint glyphs. These come in matched
  # pairs for the top and bottom endpoints. Thus the 
  # construction of the alternating vector.
  # ======================================================
  
  avec <- alter(orderedEdges$top, orderedEdges$bot)
  xleft <- (duper(1:numE) * grid) - halfBox
  xright <- (duper(1:numE) * grid) + halfBox
  ybottom <- (numV * grid) - (avec * grid) - halfBox
  ytop <- (numV * grid) - (avec * grid) + halfBox
  rect(xleft, ybottom, xright, ytop, col = mcolG, border = "black", lwd = cFac)
  
  # ======================================================
  # Directed graphs get the arrow glyph added to the target
  # glyph. Looks like the polygon call likes to do one at
  # a time, so we loop. We have carried the target node along
  # in the orderedEdges frame until now so we can use it to
  # place the glyph, and we need to decide if the glyph goes
  # on top or on the bottom
  # ======================================================
  
  if (isDirected) {
    for (i in 1:numE) {
      arrowNode <- orderedEdges[i,]$targ 
      avLeftX <- (i * grid) - halfBox
      avRightX <- (i * grid) + halfBox
      avMidX <- (i * grid)
      # Need to decide if target arrow goes on the top or the bottom:
      botNode <- orderedEdges[i,]$bot
      onTop <- (arrowNode == botNode)
      sign <- -1
      if (onTop) {
        sign <- 1
      }
      avLeftY <- (numV * grid) - (arrowNode * grid) + (sign * (2 * halfBox))
      avRightY <- (numV * grid) - (arrowNode * grid) + (sign * (2 * halfBox))
      avMidY <- (numV * grid) - (arrowNode * grid) + (sign * halfBox)
      xVals <- c(avLeftX, avRightX, avMidX)
      yVals <- c(avLeftY, avRightY, avMidY)
      aColIndx <- ((i - 1) %% length(mcolL)) + 1
      polygon(xVals, yVals, col = mcolL[[aColIndx]], border = "black", lwd = cFac)
    }
  }
  
  # ======================================================
  # Node labels
  # ======================================================
  
  she <- strheight("XXQQ", units="in", cex=1.0)
  fracH <- (gridInMin / she)
  nodeLabelCex <- fracH * nodeLabelFrac
  if (nodeLabelCex < minimumTextCexWithoutFailure) {
    nodeLabelCex <- minimumTextCexWithoutFailure;
  }
  x0 <- nodeMin * grid
  y0 <- (numV * grid) - (1:numV * grid)
  text(x0, y0, bfGraphLabels, pos=2, offset=(nodeLabelOffset * fracH), cex=nodeLabelCex)

  # ======================================================
  # Zone labels
  # ======================================================
  
  numlab <- length(bfGraphLabels)
  for (i in 1:numlab) {
    zMid <- zoneMid[[i]]
    # No valid zone? Skip it!
    if (is.na(zMid)) {
      next
    }
    # Figure out the length and location of the zone
    currLab <- bfGraphLabels[[i]]
    slen <- strwidth(currLab, units="in", cex=1.0)
    currX <- zMid * grid
    currY <- (numV * grid) - (i * grid)
    # A zone with max == min still has one column!
    zoneLen <- (zoneBoundsMax[i] - zoneBoundsMin[i] + 1) * gridInMin;
    
    # Minimum label size is the node line label. If we need to
    # shrink the node zone label below that size to make it fit,
    # we don't. We rotate the label to vertical instead. We also cap
    # the size to keep from getting monster labels
    
    frac <- (nodeZoneLabelFrac * zoneLen) / slen
    nzh <- strheight("XXQQ", units="in", cex=frac)
    gridNZH <- nzh / gridInMin
    if (gridNZH > maxTextHeightGrid) {
      frac <- frac * (maxTextHeightGrid / gridNZH);
    }
    
    if (frac > nodeLabelCex) {
      currY <- currY + halfBox + 5 # bump up a bit...
      text(currX, currY, currLab, offset=0, pos=3, cex=frac)
    } else {
      halfStrAsGrid <- (((slen) * nodeLabelCex) / 2) / gridInMin
      currY <- currY + (halfStrAsGrid * grid) + grid # bump up a bit...
      text(currX, currY, currLab, offset=0, srt=90, cex=nodeLabelCex)
    }  
  }
  # Done! Thanks for using BioFabric. Questions? Email wlongabaugh@systemsbiology.org
  return
}
