## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(RScelestial)
# We load igraph for drawing trees. If you do not want to draw,
# there is no need to import igraph.
library(igraph)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("RScelestial")

## -----------------------------------------------------------------------------
# Following command generates ten samples with 20 loci. 
# Rate of mutations on each edge of the evolutionary tree is 1.5. 
D = synthesis(10, 20, 5, seed = 7)
D

## -----------------------------------------------------------------------------
seq = as.ten.state.matrix(D$seqeunce)
SP = scelestial(seq, return.graph = TRUE)
SP

## ----fig.width=5, fig.height=5------------------------------------------------
tree.plot(SP, vertex.size = 30)

## ----fig.width=5, fig.height=5------------------------------------------------
SP = scelestial(seq, root.assign.method = "fix", root = "C8", return.graph = TRUE)
tree.plot(SP, vertex.size = 30)

## ----fig.width=5, fig.height=5------------------------------------------------
SP = scelestial(seq, root.assign.method = "balance", return.graph = TRUE)
tree.plot(SP, vertex.size = 30)

## -----------------------------------------------------------------------------
D.distance.matrix <- distance.matrix.true.tree(D)
D.distance.matrix
SP.distance.matrix <- distance.matrix.scelestial(SP)
SP.distance.matrix
## Difference between normalized distance matrices
vertices <- rownames(SP.distance.matrix)
sum(abs(D.distance.matrix[vertices,vertices] - SP.distance.matrix))

## ----load-libraries-----------------------------------------------------------
library(stringr)
library(seqinr)

## ----load-data----------------------------------------------------------------
data(phylip, package = "seqinr")


## ----data-cleaning------------------------------------------------------------
# Removing non-informative columns and duplicate rows.
mcb <-  toupper(t(sapply(seq(phylip$seq), function(i) unlist(strsplit(phylip$seq[[i]], '')))))
ccb <- as.character(phylip$seq)
occb <- order(ccb)
cbColMask <- sapply(seq(ncol(mcb)), function(j) length(levels(as.factor(mcb[,j]))) == 1)
cbRowMask <- rep(TRUE, length(ccb))
for (i in seq(length(ccb))) {
    if (i == 1 || ccb[occb[i]] != ccb[occb[i-1]]) {
        cbRowMask[occb[i]] <- FALSE
    }
}
mcbRows <- apply(mcb[!cbRowMask, !cbColMask], MARGIN = 1, FUN = function(a) paste0(str_replace(a, "-", "X"), collapse = ""))


## ----run-scelestial-----------------------------------------------------------
n.seq <- data.frame(nodes = phylip$nam[!cbRowMask], seq = mcbRows)
seq2 <- data.frame(t(as.ten.state.matrix.from.node.seq(n.seq)), stringsAsFactors = TRUE)
# Running Scelestial
SP = scelestial(seq2, return.graph = TRUE)

## ----plot---------------------------------------------------------------------
tree.plot(SP, vertex.size=20, vertex.label.dist=0, asp = 0, vertex.label.cex = 1)

