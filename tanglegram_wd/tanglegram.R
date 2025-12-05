############################################################
## All-in-one tanglegram script
############################################################

library(ape)
library(phangorn)
library(phylogram)
library(dendextend)
library(RColorBrewer)

## 1. Read trees ----
tree_aa_file <- "aa_tree.support.nwk"
tree_nt_file <- "nt_tree_clean.nwk"

tree_aa <- read.tree(tree_aa_file)
tree_nt <- read.tree(tree_nt_file)

## 2. Prune both trees to shared tips ----
shared_tips <- intersect(tree_aa$tip.label, tree_nt$tip.label)

tree_aa_pruned <- keep.tip(tree_aa, shared_tips)
tree_nt_pruned <- keep.tip(tree_nt, shared_tips)

## 3. Midpoint-root both trees ----
tree_aa_rooted <- midpoint(tree_aa_pruned)
tree_nt_rooted <- midpoint(tree_nt_pruned)

## 4 Prune E. albertii clade from both trees ----
ea_ids <- c(
  "HEB1223389.1",
  "HEB1528362.1",
  "WMV65795.1",
  "MCE7719977.1",
  "MCZ8639243.1",
  "EPF4820813.1",
  "EPF4174552.1",
  "EEW3329860.1",
  "MCZ8625746.1",
  "WDC28796.1",
  "MCQ8939922.1"
)

no_albert_aa <- tree_aa_rooted %>% prune(ea_ids)
no_albert_nt <- tree_nt_rooted %>% prune(ea_ids)

no_albert_aa <- ladderize(no_albert_aa)
no_albert_nt <- ladderize(no_albert_nt)



## 5. Convert rooted trees to dendrograms ----
dend_aa <- phylogram::as.dendrogram.phylo(no_albert_aa)
dend_nt <- phylogram::as.dendrogram.phylo(no_albert_nt)

## 6. Rotate NT dendrogram to match AA label order ----
target_order <- labels(dend_aa)
dend_nt_rot  <- dendextend::rotate(dend_nt, order = target_order)

## 6b. UNHANG / ALIGN LEAVES ----
dend_aa_flat <- dend_aa %>%
  hang.dendrogram %>%
  hang.dendrogram(hang = -1)

dend_nt_flat <- dend_nt %>%
  hang.dendrogram %>%
  hang.dendrogram(hang = -1)

## -- GRADIENT AA-longer to NT-longer

# root-to-tip distances for AA and NT trees
depth_aa  <- node.depth.edgelength(no_albert_aa)
depth_nt  <- node.depth.edgelength(no_albert_nt)

rt_aa <- depth_aa[1:Ntip(no_albert_aa)]
rt_nt <- depth_nt[1:Ntip(no_albert_nt)]

names(rt_aa) <- no_albert_aa$tip.label
names(rt_nt) <- no_albert_nt$tip.label

shared <- intersect(names(rt_aa), names(rt_nt))

df_rt <- data.frame(
  tip = shared,
  aa  = rt_aa[shared],
  nt  = rt_nt[shared],
  row.names = shared
)

# difference in normalized units:
df_rt$diff <- scale(df_rt$aa - df_rt$nt)[,1]  # z-score scaled

# continuous palette from blue → white → red
pal <- colorRampPalette(c("blue", "white", "red"))

# map each tip to a colour
cols_diff <- pal(100)[cut(df_rt$diff, breaks = 100)]
names(cols_diff) <- df_rt$tip

# apply gradient colours to tip labels
dend_aa_flat <- dend_aa_flat %>%
  set("labels_col", cols_diff[labels(dend_aa_flat)])

dend_nt_flat <- dend_nt_flat %>%
  set("labels_col", cols_diff[labels(dend_nt_flat)])

# --

## 7. Annotations: colour tip labels by isolation source ----

meta <- read.csv("anno_id_spp.csv", stringsAsFactors = FALSE)
rownames(meta) <- meta$tip_id   # protein IDs as rownames

labs_aa <- labels(dend_aa_flat)
labs_nt <- labels(dend_nt_flat)
all_labs <- unique(c(labs_aa, labs_nt))

# isolation source for ALL IDs in meta (named by protein ID)
iso_all <- meta$iso_source
names(iso_all) <- meta$tip_id
iso_all[is.na(iso_all)] <- "Unknown"

# define colours per isolation source based on actual levels
iso_levels <- sort(unique(iso_all))
n_iso <- length(iso_levels)
iso_cols <- setNames(
  brewer.pal(min(max(n_iso, 3), 8), "Set1")[seq_len(n_iso)],
  iso_levels
)

# get per-tree isolation source in label order
iso_aa <- iso_all[labs_aa]
iso_nt <- iso_all[labs_nt]

# get per-tree colour vectors in correct order
cols_aa <- iso_cols[iso_aa]
cols_nt <- iso_cols[iso_nt]

# apply label colours
dend_aa_flat <- dend_aa_flat %>%
  set("labels_col", cols_aa)

dend_nt_flat <- dend_nt_flat %>%
  set("labels_col", cols_nt)

## 8. Plot tanglegram in RStudio ----
par(mar = c(1, 4, 1, 1))
old_par <- par(font = 2)

dl <- dendlist(dend_aa_flat, dend_nt_flat) %>%
  set("highlight_branches_col")   # keep branch colouring you like

tanglegram(
  dl,
  sort = FALSE,
  main_left  = "AA tree (midpoint-rooted)",
  main_right = "NT tree (midpoint-rooted)",
  main       = "AA vs NT tanglegram (unhung leaves)",
  lab.cex    = 1,
  margin_inner = 6.5,
  lwd = 2,
  color_lines = "black",
  edge.lwd = 2,
  dLeaf_left  = -0.002,
  dLeaf_right = 0.001,
  highlight_distinct_edges = FALSE,
  highlight_branches_lwd = FALSE,
  common_subtrees_color_lines = TRUE,
  common_subtrees_color_branches = FALSE
)

grad_cols <- pal(3)

# gradient legend
legend(
  "bottomleft",                          # or "topright" etc.
  legend = c("NT longer", "similar", "AA longer"),
  col    = grad_cols,
  pch    = 15,
  pt.cex = 1.5,
  bty    = "n",
  title  = "AA vs NT branch length"
)

# iso source legend
legend(
  "topleft",
  legend = names(iso_cols),
  col    = iso_cols,
  pch    = 15,
  pt.cex = 1.5,
  bty    = "n",
  title  = "Isolation source"
)

cat("Done! Check the Plots pane.\n")