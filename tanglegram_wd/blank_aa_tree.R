library(dendextend)
library(ape)
library(BiocManager)
library(DECIPHER)
library(phylogram)
library(phangorn)

# ==============================================================================
# 1. SETUP & DATA LOADING
# ==============================================================================

# Read tree and midpoint-root
tree <- read.tree("aa_tree.support.nwk")
tree_mid <- midpoint(tree)

# Read annotation CSV (Only strictly needed for species names now)
serovar  <- read.csv("annotation_spp_type.csv", stringsAsFactors = FALSE)

# ==============================================================================
# 2. DATA SAFEGUARDS & MATCHING
# ==============================================================================

# Helper to ensure 'tip.id' column exists
fix_id_column <- function(df, filename) {
  if("tip.id" %in% colnames(df)) return(df)
  colnames(df)[1] <- "tip.id"
  return(df)
}

serovar <- fix_id_column(serovar, "annotation_spp_type.csv")
serovar$tip.id <- trimws(as.character(serovar$tip.id))

# Get clean tree IDs
raw_tree_ids   <- tree_mid$tip.label
clean_tree_ids <- trimws(sub("/.*", "", raw_tree_ids))

# Match Species to Tree Order
match_idx_spp <- match(clean_tree_ids, serovar$tip.id)
species_vec   <- serovar$species[match_idx_spp]

# ==============================================================================
# 3. LABEL CONSTRUCTION
# ==============================================================================

final_labels <- paste0(raw_tree_ids, " | ", species_vec)
final_labels[is.na(species_vec)] <- raw_tree_ids[is.na(species_vec)] 

tree_mid$tip.label <- final_labels
tree_lad <- ladderize(tree_mid)

# ==============================================================================
# 4. PLOTTING (PLAIN STYLE)
# ==============================================================================

pdf("aa_tree_plain.pdf", width = 12, height = 15)

# Margins: Reduced right margin since we removed the large isolation legend
par(mar = c(2, 2, 3, 2) + 0.1)

# Plot Tree WITH labels directly
plotBreakLongEdges(
  tree_lad,
  n = 2, 
  show.tip.label = TRUE,  # <--- Plain labels ON
  cex = 0.7,              # Adjust font size if text is too small/large
  label.offset = 0.002,   # Small gap between branch and text
  main = "Amino Acid Tree",
  no.margin = FALSE
)

# Node Support (Bootstrap > 70)
boots <- as.numeric(tree$node.label)
nodelabels(
  pch = ifelse(boots >= 70, 19, NA),
  col = "brown",
  cex = 1.0
)

# Node Support Legend
legend(
  "topright",
  inset  = 0.02,
  title  = "Node support",
  legend = "≥ 70% bootstrap",
  pch    = 19,
  col    = "brown",
  pt.cex = 1.2,
  cex    = 1.0,
  bty    = "n"
)

# Scale Bar
add.scale.bar(lwd = 3, length = 0.05, cex = 1.0)

dev.off()

print("Plot saved as 'aa_tree_plain.pdf'")