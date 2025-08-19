
# Script that prunes the 4705-tip Álvarez tree with species from Mammal_PupilTraits_GenusLevel.csv

if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
library(ape)

csv_file <- "Mammal_PupilTraits_GenusLevel.csv"
df <- read.csv(csv_file, stringsAsFactors = FALSE)

df$Alvarez.Carretero_Tiplabel <- trimws(df$Alvarez.Carretero_Tiplabel)
species_list <- unique(na.omit(df$Alvarez.Carretero_Tiplabel))

tree_file <- "4705sp_mean.nwk"

tree <- read.tree(tree_file)

tree$tip.label <- trimws(tree$tip.label)

keep_tips <- intersect(tree$tip.label, species_list)

pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, keep_tips))

#Save pruned tree
out_file <- "4705sp_mean_pruned.nwk"
write.tree(pruned_tree, out_file)
