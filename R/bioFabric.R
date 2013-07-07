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
# ================================================================================
#'@title Permutation utility function for BioFabric
#'
#'@description
#' Internal BioFabric utility function for generating a permutation vector
#'
#'@details
#' Given a vector of integers H, returns a vector V where V[n] gives
#' the index of integer n in H
#'
#'@param ordered A vector of integers. A vector of length n must contain
#'all integers from 1 to n (currently not checked)
#'
#'@return A vector V where V[n] gives the index of integer n in the input vector
#'
#'@keywords internal
#'
#'@noRd
#'
#'@author Bill Longabaugh <wlongabaugh@@systemsbiology.org>
#'
#'@examples
#'
#' \dontrun{
#' myord <- c(10, 3, 5, 2, 1, 9, 4, 7, 6, 8)
#' po <- permer(myord)
#' }
#' 

permer <- function(ordered) {
  numV <- length(ordered)
  nodePerm <- vector(length=numV)
  count <- 1
  for (i in ordered) {
    nodePerm[i] <- count
    count <- count + 1
  }
  return (nodePerm)
}

# ================================================================================
#'@title Default layout function for BioFabric
#'
#'@description
#'Does breadth-first layout of a network for BioFabric
#'
#'@details
#' Default node layout based on breadth-first search from
#' highest degree node, searching in order of decreasing
#' degree. Takes a network provided in igraph format (see 
#' \code{\link{igraph}}) and returns the same network, with
#' the nodes reordered to the default breadth-first layout
#' order
#'
#'@param inGraph The graph to order, in igraph format (see 
#'\code{\link{igraph}}). Node names MUST already be set using 
#'the vertex attribute "name".
#'@param firstLabel (optional) The label of the node to start at.
#'
#'@author Bill Longabaugh <wlongabaugh@@systemsbiology.org>
#'
#'@return The provided graph with the nodes reordered per the
#'BioFabric default layout algorithm.
#'
#'@export
#'@import igraph
#'
#'@examples
#'
#' \dontrun{
#' # Gotta have igraph!
#' 
#' library(igraph)
#' 
#' # defaultNodeOrder requires we name the graph, so we
#' # do that first. Typically, there is no reason to provide
#' # this function to bioFabric, since it uses it internally.
#' # However, we can gain extra functionality by providing the
#' # firstLabel argument, which starts the search at the specified
#' # node instead od the highest degree node. So wrapping the
#' # defaultNodeOrder as shown allows us to start at the top of
#' # the given tree (degree = 2), which is lower degree than
#' # its child nodes.
#' 
#' bfGraph = graph.tree(20, children=2, mode="out")
#' bfGraph <- autoNameForFabric(bfGraph)
#' startAtBF1 <- function(bfGraph) {
#'   return (defaultNodeOrder(bfGraph, firstLabel=V(bfGraph)[1]$name))
#' } 
#' bioFabric(bfGraph, orderFunc=startAtBF1)  
#' }

defaultNodeOrder <- function(bfGraph, firstLabel = NULL) {
  
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
  
  vlabels <- get.vertex.attribute(bfGraph, "name")
  if (is.null(vlabels)) {
    stop("Graph must have named nodes")
  }
  
  # build the frame, sort it:
  degrees <- degree(bfGraph)
  nodesToOrder <- data.frame(degrees, vlabels)
  ordered <- order(-nodesToOrder$degrees, nodesToOrder$vlabels)
  
  # Generate the permutation vector:
  nodePerm <- permer(ordered)
  
  # Permute the order, run the breath-first search. Note
  # that we ignore direction when doing this search on
  # directed graphs
  
  bfGraphPass1 <- permute.vertices(bfGraph, nodePerm)
  
  if (is.null(firstLabel)) {
    startIndex = 1
  } else {
    startIndex = which(V(bfGraphPass1)$name == firstLabel)
  }
  bfGraphBFS <- graph.bfs(bfGraphPass1, startIndex, neimode="all")
  
  # ======================================================
  # Use the results of the search to permute the graph
  # into the final order we need to lay it out.
  # ======================================================
  
  # Generate the final permutation vector that will
  # create the node order we need for the default
  # BioFabric layout:
  
  nodePerm2 <- permer(bfGraphBFS$order)
  
  # ======================================================
  # Permute the order, thus creating the final
  # graph we will layout. Extract the vector of
  # labels, which we use to tag the node lines:
  # ======================================================
  
  bfGraph <- permute.vertices(bfGraphPass1, nodePerm2)
  return (bfGraph)
}

# ================================================================================
#'@title Pass-through layout function for BioFabric
#'
#'@description
#'Does nothing
#'
#'@details
#' Returns the provided graph with no changes. Use this function
#' if you are providing a graph to bioFabric that is already
#' node-ordered, so no reordering will occur.
#'
#'@param inGraph The graph to order, in igraph format (see 
#'\code{\link{igraph}}). Node names MUST already be set using 
#'the vertex attribute "name".
#'
#'@return The provided graph with no reordering.
#'
#'@author Bill Longabaugh <wlongabaugh@@systemsbiology.org>
#'
#'@export
#'@import igraph
#'
#'@examples
#'
#' \dontrun{
#' library(igraph) 
#' bfGraph = graph.tree(20, children=2, mode="out")
#' bfGraph <- autoNameForFabric(bfGraph)
#' bfGraph <- defaultNodeOrder(bfGraph)
#' bioFabric(bfGraph, orderFunc=passthroughNodeOrder)  
#' }

passthroughNodeOrder <- function(bfGraph) {
  return (bfGraph)
}

# ================================================================================
#'@title Node naming utility function for BioFabric
#'
#'@description
#' BioFabric requires nodes have names. This will do it!
#'
#'@details
#' BioFabric requires that all nodes have names. If network is
#' named, there is no change. If nodes are unnamed but have the 
#' attribute "label", that attribute is copied into the "name"
#' attribute. Otherwise, names of the form "BFn" are formed,
#' where n is the index of the node.
#' 
#'@param inGraph The graph to name, in igraph format (see 
#'\code{\link{igraph}}).
#'
#'@return The provided graph with nodes named
#'
#'@author Bill Longabaugh <wlongabaugh@@systemsbiology.org>
#'
#'@examples
#'
#' \dontrun{
# Gotta have igraph!
#' 
#' library(igraph)
#' 
#' bfGraph = graph.tree(20, children=2, mode="out")
#' bfGraph <- autoNameForFabric(bfGraph)
#' }
#' 

autoNameForFabric <- function(inGraph) {
  
  vlabels <- get.vertex.attribute(inGraph, "name")
  if (is.null(vlabels)) {
    vlabels <- get.vertex.attribute(inGraph, "label")
    if (is.null(vlabels)) {
      vlabels <- paste("BF", 1:vcount(inGraph), sep="")
    } 
    bfGraph <- set.vertex.attribute(inGraph, "name", value=vlabels); 
  } else {
    bfGraph <- inGraph
  }
  return (bfGraph)
}

# ================================================================================
#'@title BioFabric for R
#'
#'@description
#'An R implementation of the BioFabric network visualization tool
#'
#'@details
#'Plots a network, provided in igraph format (see \code{\link{igraph}}),
#'using BioFabric. This release supports user-specified node orders and
#'the display of shadow links. 
#'
#'Best results are obtained using a PDF display target, as shown 
#'in the example. BioFabric is very sensitive to correct relative scaling of 
#'linewidth, label size, and grid size. The PDF output has been targeted 
#'for correctly handling this scaling. PNG output should also be acceptable.
#'
#'SIZING: As the network grows, you will need to increase the size of the
#'output to be able to get per-link resolution. PDF files that are not
#'big enough (e.g. only 10 inches wide) will create links below the hairline 
#'limit and may not print well. PNG outputs are the better bet to create
#'low-resolution images of large networks in page-size dimensions. To be
#'able to view an enormous network and zoom to see each link, create an 
#'enormous PDF file (e.g. 40, 80, 100+ inches wide!) or PNG file (e.g.
#'30,000 pixels wide if the network has 5000 links).
#'
#'VIEWERS: PDF or PNG viewers that cannot do antialiasing of line art are
#'a poor choice for BioFabric networks. The closely spaced parallel lines
#'of a BioFabric plot MUST be antialiased to get acceptable results. Some
#'specific viewer notes:
#'
#' Evince on Linux for PDF (Document Viewer 2.30.3 tested): Antialiased 
#' is fixed on, and visuals are good. Problems: Very tiny text below
#' some size threshold explodes to a huge size. Cannot zoom above 400%,
#' which is insufficient to explore network. No drag hand cursor to
#' navigate, which is essential. Note that for Postscript output,
#' Evince does not antialias, giving very poor results. 
#' 
#' Preview on Mac for PDF (Version 4.2 tested): Be sure that "Anti-alias 
#' text and line art" is checked on the PDF tab for Preferences, which
#' gives fair visuals. Good maximum zoom level, and the Move cursor
#' provides convenient draggable navigation.
#' 
#' Adobe Reader on Mac for PDF (Version 9.5.5 tested): Be sure that 
#' "Smooth line art" is checked, and "Enhance thin lines" is NOT checked,
#' on the Page Display Preferences. The Hand tool is available via
#' Tools->Select & Zoom->Hand Tool, and the maximum 6400 percent zoom level
#' is good for exploration.
#' 
#' RStudio (0.97 tested): You won't be able to see much in the tiny 
#' Plot Frame. Due to a known RStudio issue, the Zoom window does
#' not provide correct label sizing or line width sizing, since the
#' plots are scaled based on the Plot Frame dimensions. Use the PDF 
#' or PNG output routes for now.
#' 
#' General note on image viewers: If you have a plot with 10,000 links,
#' you are quickly heading towards PNG images that are over 32,000
#' pixels wide for high resolution. Lots of image viewers start to
#' get cranky when viewing images that are more than 32,000 pixels
#' wide.
#' 
#' If you want to explore large networks with BioFabric, best to 
#' export the network as a tab-delimited SIF file and use the 
#' Java-based BioFabric application, which is designed for 
#' interactive large network visualization and exploration.
#'
#'@param inGraph The graph to plot, in igraph format (see 
#'\code{\link{igraph}}). Node names should already be set using the vertex
#'attribute "name", though the attribute "label" will work as well.
#'If not set, the function will assign names of the form "BFn".
#'
#'@param userOrder (optional) An ordered list specifying the ordering of the
#'node rows in the plot; nodes are identified using node name. All nodes must 
#'be present in the list. If this list is provided, the default node ordering 
#'routine will not be used. 
#'
#'@param orderFunc (optional) A function that takes a graph to plot and returns
#'the graph with the nodes reordered in the desired node order. Providing the
#'defaultNodeOrder function for this parameter gives the default behavior. Note
#'that if this parameter is provided, the userOrder parameter is ignored.
#'
#'@param shadowLinks (optional) Default value is FALSE. If set to TRUE, BioFabric
#'will display shadow links. See the BioFabric paper or blog for more info.
#'
#'@param dropNodeLabels (optional) Default value is FALSE. If set to TRUE, no
#'node labels are displayed. With large networks, some PDF readers cannot
#'correctly handle tiny text, and labels may need to be dropped.
#'
#'@param dropZoneLabels (optional) Default value is FALSE. If set to TRUE, no
#'node zone labels are displayed. With large networks, some PDF readers cannot
#'correctly handle tiny text, and labels may need to be dropped.
#'
#'@author Bill Longabaugh <wlongabaugh@@systemsbiology.org>
#'
#'@export
#'@import igraph
#'
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
#' plotHeight <- plotWidth * (aspect * 1.2)
#' pdf("myBioFabricOutput.pdf", width=plotWidth, height=plotHeight)
#' bioFabric(bfGraph)
#' dev.off()
#' }
#' 
#'@references Home page: \url{http://www.BioFabric.org}
#'@references Blog: \url{http://biofabric.blogspot.com}  
#'@references Paper: \url{http://www.biomedcentral.com/1471-2105/13/275}
#'

bioFabric <- function(inGraph, userOrder=NULL, orderFunc=NULL, 
                      shadowLinks=FALSE, dropNodeLabels=FALSE,
                      dropZoneLabels=FALSE) {

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
      taggedCol <- cList[[order[i]]]
      transCol <- col2rgb(taggedCol$col)
      myR <- min(255, transCol[1] * mult);
      myG <- min(255, transCol[2] * mult);
      myB <- min(255, transCol[3] * mult);
      retval[i] <- rgb(myR, myG, myB, maxColorValue = 255);
    }
    return (retval)
  }
  
  # ======================================================
  # Create a named color from integer RGB values
  # ======================================================
  
  buildCol <- function(useName, r, g, b) { 
    return (list(name = useName, col = rgb(r, g, b, maxColorValue = 255)))
  }
  
  # ======================================================
  # Utility function to create permutation vector given 
  # existing name order and new name order
  # ======================================================
  
  namePermer <- function(currNameOrder, newNameOrder) {
    numV <- length(currNameOrder)
    nodePerm <- vector(length=numV)
    for (i in 1:numV) {
      newIndex = which(newNameOrder == currNameOrder[i])
      nodePerm[i] <- newIndex
    }
    return (nodePerm)
  }
  
  # ======================================================
  # Default node layout based on breadth-first search from
  # highest degree node, searching in order of decreasing
  # degree
  # ======================================================
  
  orderFromList <- function(bfGraph, nameList) {
    nodePerm <- namePermer(V(bfGraph)$name, nameList)
    return (permute.vertices(bfGraph, nodePerm))
  }

  # ======================================================
  # Default (non-shadow) node zone extraction
  # ======================================================
  
  findNodeZone <- function(orderedEdges, match) {
    smatch <- sort(match)
    maxCol <- max(match)
    matLen <- length(smatch)
    retval = list(zbmax=maxCol, zbmin=maxCol, zbmid=maxCol)
    if (matLen > 1) {      
      for (j in matLen:2) {
        # Decrement the min if it is contiguous:
        if ((smatch[j] - 1) == smatch[j - 1]) {
          edgeAtCol <- orderedEdges[smatch[j - 1],]
          # if still contig, but now not at top, we are done too:
          if (is.na(edgeAtCol$top) || (edgeAtCol$top != i)) {
            break;
          }       
          retval$zbmin <- smatch[j]
          retval$zbmid <- (maxCol + smatch[j - 1]) / 2.0
        } else {
          break;
        }
      }
    }

    # ======================================================
    # We do not label all the node zones. The zone has to exist
    # along the top edge of the network. 
    # ======================================================
    
    edgeAtMax <- orderedEdges[maxCol,]
    if (is.na(edgeAtMax$top) || (edgeAtMax$top != i)) {
      retval$zbmin <- NA
      retval$zbmax <- NA
      retval$zbmid <- NA
    }
    return (retval)
  }
  
  # ======================================================
  # Shadow links node zone extraction
  # ======================================================
  
  findNodeZoneShadow <- function(frameForNodeMinMax, nodeNum, orderedEdges) {
    match <- frameForNodeMinMax[(frameForNodeMinMax$zone == nodeNum),]$indx
    maxCol <- max(match)
    smatch <- sort(match)  
    matLen <- length(smatch)
    retval = list(zbmax=maxCol, zbmin=maxCol, zbmid=maxCol)
    if (matLen > 1) {      
      for (j in matLen:2) {
        # Decrement the min if it is contiguous:
        if ((smatch[j] - 1) == smatch[j - 1]) {
          retval$zbmin <- smatch[j]
          retval$zbmid <- (maxCol + smatch[j - 1]) / 2.0
        } else {
          break;
        }
      }
    }
    return (retval)
  }
    
  # ======================================================
  # Define the colors used. These colors were originally
  # derived for use in BioTapestry (www.BioTapestry.org)
  # to create a set of colors that could be distinguished
  # at link crossings.
  # ======================================================
  
  baseColorColors <- list(buildCol("EX-cyan", 0, 255, 255),
                          buildCol("EX-dark-cyan", 0, 100, 128),
                          buildCol("EX-yellow-orange", 255, 153, 0),
                          buildCol("EX-pale-green", 133, 205, 102),
                          buildCol("EX-dark-green", 39, 128, 0),
                          buildCol("EX-pale-red-orange", 230, 156, 138),
                          buildCol("EX-yellow-green", 154, 255, 0),
                          buildCol("EX-yellow", 255, 203, 0),
                          buildCol("EX-dark-gray-purple", 0, 25, 128),
                          buildCol("EX-pale-magenta", 212, 138, 230),
                          buildCol("EX-pale-purple", 149, 165, 230),
                          buildCol("EX-purple", 102, 51, 255),
                          buildCol("EX-dark-red", 140, 56, 56),
                          buildCol("EX-red", 255, 0, 0),
                          buildCol("EX-pale-yellow-green", 222, 230, 138),
                          buildCol("EX-dark-purple", 77, 56, 140),
                          buildCol("EX-pale-cyan", 138, 230, 181),
                          buildCol("EX-pure-blue", 0, 0, 255),
                          buildCol("EX-dark-yellow-green", 114, 128, 0),
                          buildCol("EX-magenta", 255, 0, 255),
                          buildCol("EX-dark-tan", 166, 133, 83),
                          buildCol("EX-pale-blue", 102, 183, 205),
                          buildCol("EX-orange", 255, 103, 0),
                          buildCol("EX-medium-magenta", 166, 83, 166),
                          buildCol("EX-blue-magenta", 155, 0, 255),
                          buildCol("EX-green", 0, 255, 0),
                          buildCol("EX-dark-magenta", 102, 0, 128),
                          buildCol("EX-pale-blue-magenta", 146, 102, 205),
                          buildCol("EX-pale-yellow orange", 205, 175, 102),
                          buildCol("EX-dark-orange", 128, 92, 0),
                          buildCol("EX-blue", 0, 152, 255),
                          buildCol("EX-pale-red", 205, 102, 153))
  
  # ======================================================
  # Create a map of color name to color def:
  # ======================================================
  
  numCol <-length(baseColorColors)
  baseColorNames <- vector(length = numCol)
  for (i in 1:numCol) {
    baseColorNames[i] <- baseColorColors[[i]]$name
  }
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
  mcolG <- as.vector(rbind(mcolL, mcolL))
    
  # ======================================================
  # Time to go to work on the graph. The graph HAS to 
  # have labeled vertices! If the igraph does not have
  # them already, we glue them on.
  # ======================================================
  
  bfGraph <- autoNameForFabric(inGraph)
  
  # ======================================================
  # We also require tags on the edges, and regular edges
  # are tagged as not being shadow links
  # ======================================================
  
  elabels <- get.edge.attribute(bfGraph, "name")
  if (is.null(elabels)) {
    elabels <- c("co") # gets recycled
    bfGraph <- set.edge.attribute(bfGraph, "name", value=elabels); 
  }
  eshads <- get.edge.attribute(bfGraph, "isShadow")
  if (is.null(eshads)) {
    eshads <- c(FALSE) # gets recycled
    bfGraph <- set.edge.attribute(bfGraph, "isShadow", value=eshads); 
  }
 
  # ======================================================
  # If the user specifies shadow links, we need to duplicate
  # the edges, retag them with the shadow qualifier, and mark
  # them as shadow links.
  # ======================================================
  
  if (shadowLinks) {
    numEOrig <- ecount(bfGraph)
    baseEdges <- get.edgelist(bfGraph, names=FALSE)
    enVecShad <- paste(c("shdw("), get.edge.attribute(bfGraph, "name"), c(")"), sep="")
    # Note the transpose! Else add.edges reads DOWN columns
    bfGraph <- add.edges(bfGraph, t(baseEdges))
    shadowIndices <- (numEOrig + 1):(2 * numEOrig)
    bfGraph <- set.edge.attribute(bfGraph, "name", index=shadowIndices, value=enVecShad)
    bfGraph <- set.edge.attribute(bfGraph, "isShadow", index=shadowIndices, value=c(TRUE))
  }
  
  numV <- vcount(inGraph)
  isDirected <- is.directed(bfGraph)
 
  # ======================================================
  # Sort the nodes. We may be given a user function that
  # takes an igraph and returns an igraph where the node 
  # indices have been reordered to the desired BioFabric
  # node order (node index 1 in top row). If we have the
  # function, use it; such a function overrides a node
  # list. If we have a node list, use that to assign node
  # rows. Finally, apply the default ordering.
  # ======================================================
  
  if (!is.null(orderFunc)) {
    if (!is.null(userOrder)) {
      warning("orderFunc overrides non-null userOrder argument")
    }
    bfGraph <- orderFunc(bfGraph)
  } else if (!is.null(userOrder)) {
    bfGraph <- orderFromList(bfGraph, userOrder)
  } else {
    bfGraph <- defaultNodeOrder(bfGraph)
  }  

  bfGraphLabels <- get.vertex.attribute(bfGraph, "name")
  
  # ======================================================
  # Yay! Nodes are sorted. Now we need to sort the edges.
  # Are edges in the graph always organized so that the
  # lower-numbered node is in column 1? Maybe for undirected
  # graphs, but not for directed graphs. So make
  # sure it is the case, and build a dataFrame that we
  # can sort
  # ======================================================
  
  # DANGER! If names not set to FALSE, edgelist comes back
  # with names, not indices, if graph is named. And we MUST
  # be working with indices to get the sort to work. With
  # node names, the order operation appears to be trying
  # to create factors from the node names to sort on, and
  # chaos ensues
 
  edgelist <- get.edgelist(bfGraph, names=FALSE)
  edgeNames <- get.edge.attribute(bfGraph, "name")
  edgeShadow <- get.edge.attribute(bfGraph, "isShadow")
  numE <- ecount(bfGraph)
  
  minNodes <- vector(length=numE)
  maxNodes <- vector(length=numE)
  targNode <- vector(length=numE)
  nodeZone <- vector(length=numE)
  for (i in 1:numE) {
    node1 <- edgelist[i,1]
    node2 <- edgelist[i,2]
    minNodes[[i]] <- min(node1, node2)
    maxNodes[[i]] <- max(node1, node2)
    targNode[[i]] <- node2
    if (edgeShadow[[i]] == FALSE) {
      nodeZone[[i]] <- minNodes[[i]]
    } else {
      nodeZone[[i]] <- maxNodes[[i]]
    }
  }
  edgesForSort <- data.frame(top=minNodes, bot=maxNodes, 
                             targ=targNode, tag=edgeNames, zone=nodeZone)
  
  # ======================================================
  # Sort the edges so that edges with the topmost node go
  # first, with ties broken so the shorter link comes before
  # the longer link. Result: edge wedges!
  # ======================================================
  
  edgeOrder <- order(edgesForSort$zone, edgesForSort$top, 
                     edgesForSort$bot, edgesForSort$tag)
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
    nodeMin[i] <- min(match)
    maxCol <- max(match)
    nodeMax[i] <- maxCol
    # Node zone bounds calculated here:
    if (shadowLinks) {
      results <- findNodeZoneShadow(frameForNodeMinMax, i, orderedEdges)    
    } else {  
      results <- findNodeZone(orderedEdges, match)
    }
    zoneBoundsMin[i] <- results$zbmin
    zoneBoundsMax[i] <- results$zbmax
    zoneMid[i] <- results$zbmid
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
  # Saw poor tiny-text behavior with Evince on Linux. Now
  # allowing user to suppress labels to deal with problem:
  # minimumTextCexWithoutFailure <- 0.05 ## Labels too small? BOOM!
  
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
    if (x0[i] == x1[i]) {
      x0[i] <- x0[i] - halfBox + 1;
      x1[i] <- x1[i] + halfBox - 1;
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
  
  # This call creates a vector (a[1], b[1], a[2], b[2]...)
  avec <- as.vector(rbind(orderedEdges$top, orderedEdges$bot))
  duperE <- as.vector(rbind(1:numE, 1:numE))
  xleft <- (duperE * grid) - halfBox
  xright <- (duperE * grid) + halfBox
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
      arrowNode <- orderedEdges[i, "targ"]
      avLeftX <- (i * grid) - halfBox
      avRightX <- (i * grid) + halfBox
      avMidX <- (i * grid)
      # Need to decide if target arrow goes on the top or the bottom:
      botNode <- orderedEdges[i, "bot"]
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
      polygon(xVals, yVals, col = mcolL[aColIndx], border = "black", lwd = cFac)
    }
  }
  
  # ======================================================
  # Node labels
  # ======================================================
 
  if (!dropNodeLabels) {
    she <- strheight("XXQQ", units="in", cex=1.0)
    fracH <- (gridInMin / she)
    nodeLabelCex <- fracH * nodeLabelFrac
    x0 <- nodeMin * grid
    y0 <- (numV * grid) - (1:numV * grid)
    text(x0, y0, bfGraphLabels, pos=2, offset=(nodeLabelOffset * fracH), cex=nodeLabelCex)
  }
  
  # ======================================================
  # Zone labels
  # ======================================================
  
  if (!dropZoneLabels) {
    numlab <- length(bfGraphLabels)
    for (i in 1:numlab) {
      zMid <- zoneMid[i]
      # No valid zone? Skip it!
      if (is.na(zMid)) {
        next
      }
      # Figure out the length and location of the zone
      currLab <- bfGraphLabels[i]
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
  }
  # Done! Thanks for using BioFabric. Questions? Email wlongabaugh@systemsbiology.org
}
