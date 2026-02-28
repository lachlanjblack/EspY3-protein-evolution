library(dendextend)
library(ape)
library(BiocManager)
library(DECIPHER)
library(phylogram)
library(phangorn)

tree <- read.tree("nt_tree_clean.nwk")
tree_mid <- midpoint(tree)

serovar <- read.csv("annotation_spp_type.csv", stringsAsFactors = FALSE)
rownames(serovar) <- serovar$tip.id   # tip.id must match ORIGINAL tip labels in the tree

## 1. Keep original protein IDs
protein_ids <- tree_mid$tip.label

## 2. Look up species for each protein id
species <- serovar[protein_ids, "species"]

## 3. Build new labels: [protein_id] | Species
new_labels <- paste0(protein_ids, " | ", species)

## 4. Put these on the tree
tree_mid$tip.label <- new_labels

## 5. Ladderise
tree_lad <- ladderize(tree_mid)

## 6. Plot
par(mar = c(5, 4, 4, 3) + 0.1)
plotBreakLongEdges(
  tree_lad,
  n = 2,
  show.tip.label = TRUE,
  cex = 0.5,
  main = "Nucleotide Tree",
  no.margin = FALSE,
  label.offset = 0.001
)

# Bootstrap dots
boots <- as.numeric(tree$node.label)
nodelabels(
  pch = ifelse(boots >= 70, 19, NA),
  col = "brown",
  cex = 0.7
)

legend(
  "topright",
  inset  = 0.01,
  title  = "Node support",
  legend = "≥ 70% bootstrap",
  pch    = 19,
  col    = "brown",
  pt.cex = 0.9,
  cex    = 0.9,
  bty    = "n"
)

# Scale bar position
usr <- par("usr")
x_left   <- usr[1] + 0.02 * diff(usr[1:2])    # 2% from left
y_bottom <- usr[3] + 0.04 * diff(usr[3:4])    # 2% from bottom

scale_len <- 0.05   # <-- adjust to match your branch lengths

# Draw scale bar
add.scale.bar(
  x      = x_left,
  y      = y_bottom,
  length = scale_len,
  lwd    = 2,
  col    = "black",
  ask    = FALSE
)