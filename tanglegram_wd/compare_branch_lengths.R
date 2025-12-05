library(ape)

# Root-to-tip distances (tips are 1:Ntip)
depth_aa_all <- node.depth.edgelength(no_albert_aa)
depth_nt_all <- node.depth.edgelength(no_albert_nt)

ntips_aa <- Ntip(no_albert_aa)
ntips_nt <- Ntip(no_albert_nt)

rt_aa <- depth_aa_all[1:ntips_aa]
names(rt_aa) <- no_albert_aa$tip.label

rt_nt <- depth_nt_all[1:ntips_nt]
names(rt_nt) <- no_albert_nt$tip.label

shared <- intersect(names(rt_aa), names(rt_nt))

df <- data.frame(
  tip = shared,
  aa = rt_aa[shared],
  nt = rt_nt[shared],
  row.names = shared
)

df$aa_norm <- df$aa / median(df$aa)
df$nt_norm <- df$nt / median(df$nt)

df$log2_aa_over_nt <- log2(df$aa_norm / df$nt_norm)

plot(
  df$nt_norm, df$aa_norm,
  xlab = "Normalised root-to-tip distances (NT tree)",
  ylab = "Normalised root-to-tip distances (AA tree)",
  pch = 19
)
abline(0, 1, lty = 2)

cor(df$aa_norm, df$nt_norm)         # overall correlation

# Which tips are most AA-long vs NT?
head(df[order(-df$log2_aa_over_nt), ], 10)  # AA relatively long
head(df[order(df$log2_aa_over_nt), ], 10)   # NT relatively long
