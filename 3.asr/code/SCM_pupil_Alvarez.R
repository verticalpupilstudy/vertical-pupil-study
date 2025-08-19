## This script runs stochastic character mapping (SCM) on pupil shape
# using the Alvarez-Carretero phylogeny (4705sp_mean).
#nsim=5000, ARD, pi = fitzjohn

#required libraries 
for (pkg in c("ape","phytools","dplyr","stringr","tibble")) {
  if (!requireNamespace(pkg, quietly=TRUE)) install.packages(pkg)
}
library(ape)
library(phytools)
library(dplyr)
library(stringr)
library(tibble)

#required files 
tree_file   <- "4705sp_mean_pruned.nwk"
traits_file <- "Mammal_PupilTraits_GenusLevel.csv"
tip_col     <- "Alvarez.Carretero_Tiplabel"
outdir      <- "outputs_alvarez"

# output directory 
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive=TRUE)
}
# reads tree and trait data
tree <- read.tree(tree_file)
traits_df <- read.csv(traits_file, stringsAsFactors=FALSE)

# Filter dataset to include only species in tree with valid states
valid_states <- c("Round","Vertical","Horizontal")
traits_df <- traits_df %>%
  filter(Pupil %in% valid_states,
         .data[[tip_col]] %in% tree$tip.label)

pupil_states <- setNames(traits_df$Pupil, traits_df[[tip_col]])
to_drop      <- setdiff(tree$tip.label, names(pupil_states))
if (length(to_drop) > 0) {
  tree <- drop.tip(tree, to_drop)
}

# Generate SCM with ARD and nsim = 5000, pi = fitzjohn
nsim <- 5000
simmaps <- make.simmap(
  tree  = tree,
  x     = pupil_states,
  model = "ARD",
  nsim  = nsim,
  pi    = "fitzjohn"
)

## Summaries (ancestral probs, transition counts, time in states,q matrix)
sim_desc        <- describe.simmap(simmaps)
ancestral_probs <- sim_desc$ace
write.csv(
  ancestral_probs,
  file     = file.path(outdir, "Alvarez_5000_ARD_AncestralProbs.csv"),
  row.names=TRUE
)

transition_counts <- sim_desc$count
transition_stats  <- as.data.frame(t(apply(
  transition_counts, 2,
  function(x) {
    m <- mean(x); s <- sd(x)
    c(Mean     = m,
      SD       = s,
      CI_Lower = m - 1.96 * s,
      CI_Upper = m + 1.96 * s,
      Min      = min(x),
      Max      = max(x))
  }
)))
write.csv(
  transition_stats,
  file     = file.path(outdir, "Alvarez_5000_ARD_TransitionStats.csv"),
  row.names=TRUE
)

time_df    <- as.data.frame(sim_desc$times)
time_stats <- time_df %>%
  summarise(across(
    everything(),
    list(
      Mean    = ~ mean(.),
      SD      = ~ sd(.),
      CI_Low  = ~ mean(.) - 1.96 * sd(.),
      CI_High = ~ mean(.) + 1.96 * sd(.)
    ),
    .names = "{col}_{fn}"
  ))
write.csv(
  time_stats,
  file     = file.path(outdir, "Alvarez_5000_ARD_TimeInStateStats.csv"),
  row.names=FALSE
)

qmat <- simmaps[[1]]$model$Q
saveRDS(
  qmat,
  file = file.path(outdir, "Alvarez_5000_ARD_Qmatrix.rds")
)

save.image(file = file.path(outdir, "Alvarez_5000_ARD_Results.RData"))
