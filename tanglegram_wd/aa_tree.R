library("dendextend")
library("ape")
library("BiocManager")
library("DECIPHER")

tree <- read.tree("aa_tree.support.nwk")
tree_mid <- midpoint(tree)
tree_lad <- ladderize(tree_mid)

par(mar = c(5, 4, 4, 3) + 0.1)
plotBreakLongEdges(
  tree_lad,
  n = 2, #  2 long edge breaks
  show.tip.label = TRUE,
  cex = 0.5, # label size
  main = "EspY3 Homologous Proteins",
  no.margin = TRUE,
  label.offset = 0.001,
)

# Add dots on nodes with boostraps > 70
boots <- as.numeric(tree$node.label)
nodelabels(
  pch = ifelse(boots >= 70, 19, NA),
  col = "brown",
  cex = 0.7,
  )
add.scale.bar(ask = TRUE)


# Unrooted tree (midpoint root)
plotBreakLongEdges(tree_mid, n = 2, "u", show.tip.label = FALSE)

## Annotate Species and Sero

serovar <- read.csv("annotation_spp_type.csv")
rownames(serovar) <- serovar$tip.id

tree_mid$tip.label <- serovar[tree_mid$tip.label, "species"]

# Prune albertii clade
# no_albert <- tree_mid %>% prune(c("HEB1223389.1",
"HEB1528362.1",
"WMV65795.1",
"MCE7719977.1",
"MCZ8639243.1",
"EPF4820813.1",
"EPF4174552.1",
"EEW3329860.1",
"MCZ8625746.1",
"WDC28796.1",
"MCQ8939922.1"))