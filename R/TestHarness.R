#
# Test Script for BioFabric
#

library(igraph)
#bfGraph = read.graph("/home/bill/lesmiserables.txt", format="gml")
bfGraph = barabasi.game(100, m=6, directed=FALSE)
#bfGraph = erdos.renyi.game(100, 600, type="gnm", directed=FALSE)
height <- vcount(bfGraph)
width <- ecount(bfGraph)
aspect <- height / width;
plotWidth <- 10.0 
plotHeight <- plotWidth * (aspect * 1.2) 
pdf("/home/bill/testOut.pdf", width=plotWidth, height=plotHeight)
bioFabric(bfGraph)
dev.off()
