library(dendextend)
library(ape)
library(BiocManager)
library(DECIPHER)
library(phylogram)
library(phangorn)
library(RColorBrewer)

# ==============================================================================
# 1. SETUP & DATA LOADING
# ==============================================================================

# Read tree and midpoint-root
tree <- read.tree("aa_tree.support.nwk")
tree_mid <- midpoint(tree)

# Read annotation CSVs
serovar  <- read.csv("annotation_spp_type.csv", stringsAsFactors = FALSE)
iso_meta <- read.csv("anno_id_spp.csv", stringsAsFactors = FALSE)

# ==============================================================================
# 2. DATA SAFEGUARDS (ROBUST MATCHING)
# ==============================================================================

# Helper to ensure 'tip.id' column exists
fix_id_column <- function(df, filename) {
  if("tip.id" %in% colnames(df)) return(df)
  colnames(df)[1] <- "tip.id"
  return(df)
}

# Apply column fix
serovar  <- fix_id_column(serovar, "annotation_spp_type.csv")
iso_meta <- fix_id_column(iso_meta, "anno_id_spp.csv")

# Clean whitespace
serovar$tip.id  <- trimws(as.character(serovar$tip.id))
iso_meta$tip.id <- trimws(as.character(iso_meta$tip.id))

# Get clean tree IDs
raw_tree_ids   <- tree_mid$tip.label
clean_tree_ids <- trimws(sub("/.*", "", raw_tree_ids))

# Match Species
match_idx_spp <- match(clean_tree_ids, serovar$tip.id)
species_vec   <- serovar$species[match_idx_spp]

# Match Isolation Source (DIRECTLY from your edited CSV)
match_idx_iso <- match(clean_tree_ids, iso_meta$tip.id)
iso_vec       <- iso_meta$iso_source[match_idx_iso]

# Handle NAs
iso_vec[is.na(iso_vec)] <- "Unknown"

# ==============================================================================
# 3. COLOR GENERATION
# ==============================================================================

iso_levels <- sort(unique(iso_vec))
n_cats     <- length(iso_levels)

# Robust palette generation
if(n_cats < 3) {
  base_cols <- c("lightblue", "lightpink")[1:n_cats] 
} else if (n_cats <= 9) {
  # Set3 (Pastels) is best for background highlighting
  base_cols <- brewer.pal(n_cats, "Set3") 
} else {
  base_cols <- colorRampPalette(brewer.pal(9, "Set3"))(n_cats)
}

iso_cols <- setNames(base_cols, iso_levels)
tip_bg_colors <- iso_cols[iso_vec] 

# ==============================================================================
# 4. LABEL CONSTRUCTION
# ==============================================================================

final_labels <- paste0(raw_tree_ids, " | ", species_vec)
final_labels[is.na(species_vec)] <- raw_tree_ids[is.na(species_vec)] 

tree_mid$tip.label <- final_labels
tree_lad <- ladderize(tree_mid)

# ==============================================================================
# 5. PLOTTING (PDF OUTPUT)
# ==============================================================================

pdf("aa_tree_highlighted.pdf", width = 12, height = 15)

# Margins: Bottom, Left, Top, Right (Right=10 leaves room for legend)
par(mar = c(2, 2, 3, 10) + 0.1)

# A. Plot Tree (No labels initially)
plotBreakLongEdges(
  tree_lad,
  n = 2, 
  show.tip.label = FALSE,
  main = "Amino Acid Tree",
  no.margin = FALSE
)

# B. Add Highlighted Labels
tiplabels(
  text = tree_lad$tip.label,
  bg = tip_bg_colors,      
  col = "black",           
  frame = "rect",          
  cex = 0.7,               
  adj = 0,                 
  offset = 0.005           
)

# C. Node Support (Top Right)
boots <- as.numeric(tree$node.label)
nodelabels(
  pch = ifelse(boots >= 70, 19, NA),
  col = "brown",
  cex = 1.0
)

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

# D. Isolation Source Legend (Top Right, Stacked Below)
legend(
  "topright",
   #Inset: 0.02 from right, 0.10 from top (pushes it below first legend)
   #If overlapping, increase 0.10 to 0.12 or 0.15
  inset = c(0.02, 0.10), 
  title = "Host",
  legend = names(iso_cols),
  fill = iso_cols,       
  border = "black",
  cex = 1.0,
  )

# Scale Bar
add.scale.bar(lwd = 3, length = 0.05, cex = 1.0)

dev.off()

print("Plot saved as 'aa_tree.pdf'")