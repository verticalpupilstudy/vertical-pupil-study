# this script samples 1,000 trees from Upham 10,000 tree posterior
# and prunes it with genera that have pupil data
#requires "MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_all10k_v2_nexus.trees" can be downloaded from Upham et al. 2022


if (!requireNamespace("ape", quietly=TRUE))    install.packages("ape")
if (!requireNamespace("phangorn", quietly=TRUE)) install.packages("phangorn")
if (!requireNamespace("readr", quietly=TRUE))  install.packages("readr")
if (!requireNamespace("stringr", quietly=TRUE))install.packages("stringr")

library(ape)
library(phangorn)
library(readr)
library(stringr)

posterior_all <- read.nexus("MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_all10k_v2_nexus.trees")

idx1000    <- seq(1, length(posterior_all), by = 10)[1:1000]
posterior1k <- posterior_all[idx1000]

write.tree(posterior1k, file = "sampled1000_posterior.tre")

genus_df     <- read_csv("Mammal_PupilTraits_GenusLevel.csv", show_col_types = FALSE)
desired_tips <- genus_df$Upham_TipLabel
bare_species <- str_remove(desired_tips, "_[^_]+_[^_]+$")

pruned1k <- lapply(posterior1k, function(tr) {
  tr$tip.label <- str_remove(tr$tip.label, "_[^_]+_[^_]+$")
  keep <- intersect(tr$tip.label, bare_species)
  drop <- setdiff(tr$tip.label, keep)
  drop.tip(tr, drop)
})

write.tree(pruned1k, file = "sampled1000_pruned_posterior.tre")
