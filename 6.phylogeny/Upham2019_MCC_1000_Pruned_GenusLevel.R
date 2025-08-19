#this script prunes Upham Upham MCC-consensus tree with pupil data

if (!requireNamespace("ape", quietly=TRUE)) install.packages("ape")
if (!requireNamespace("dplyr", quietly=TRUE)) install.packages("dplyr")
if (!requireNamespace("stringr", quietly=TRUE)) install.packages("stringr")

library(ape)
library(dplyr)
library(stringr)

mcc_tree <- read.tree("Upham2019_MCC_1000.tre")
full_tips <- mcc_tree$tip.label

gen_df <- read.csv(
  file.path(getwd(), "Mammal_PupilTraits_GenusLevel.csv"),
  stringsAsFactors = FALSE
)
keep_bare <- unique(gen_df$Upham_TipLabel)

bare_tips_all <- str_remove(full_tips, "_[^_]+_[^_]+$")
tip_map_all <- setNames(full_tips, bare_tips_all)

keep_full <- na.omit(tip_map_all[keep_bare])

to_drop     <- setdiff(full_tips, keep_full)
pruned_tree <- drop.tip(mcc_tree, to_drop)

pruned_tree$tip.label <- str_remove(pruned_tree$tip.label, "_[^_]+_[^_]+$")

write.tree(pruned_tree, "Upham2019_MCC_1000_Pruned_GenusLevel.tre")


