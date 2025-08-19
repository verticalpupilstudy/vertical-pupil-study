## This script runs stochastic character mapping (SCM) on pupil shape
# using the Upham_MCC phylogeny.
# nsim = 5000, model = ARD, root prior = "fitzjohn"

#required libraries
for (pkg in c("ape", "phytools", "dplyr", "stringr", "tibble")) {
  if (!requireNamespace(pkg, quietly=TRUE)) install.packages(pkg)
}
library(ape)
library(phytools)
library(dplyr)
library(stringr)
library(tibble)

#required files
tree_file   <- "Upham2019_MCC_1000_Pruned_GenusLevel.tre"
traits_file <- "Mammal_PupilTraits_GenusLevel.csv"


# output directory
outdir <- "outputs"
if (!dir.exists(outdir)) dir.create(outdir)

# read tree and trait data

tree <- read.tree(tree_file)
traits_df <- read.csv(traits_file, stringsAsFactors = FALSE)


# keep only species in the tree with valid diel states 
valid_states <- c("Round","Vertical","Horizontal")
traits_df <- traits_df %>%
  filter(Pupil %in% valid_states, Upham_TipLabel %in% tree$tip.label)

# match states to tips and prune tree

pupil_states <- setNames(traits_df$Pupil, traits_df$Upham_TipLabel)
to_drop <- setdiff(tree$tip.label, names(pupil_states))
if (length(to_drop)>0) tree <- drop.tip(tree, to_drop)

# generate SCM with ARD and nsim = 5000, pi = "fitzjohn"
nsim <- 5000
simmaps <- make.simmap(
  tree  = tree,
  x     = pupil_states,
  model = "ARD",
  nsim  = nsim,
  pi    = "fitzjohn"
)

# summaries: ancestral probabilities, transition counts, time in states, Q matrix
sim_desc        <- describe.simmap(simmaps)
ancestral_probs <- sim_desc$ace
write.csv(
  ancestral_probs,
  file = file.path(outdir, "Upham_MCC_5000_ARD_AncestralProbs.csv"),
  row.names = TRUE
)

transition_counts <- sim_desc$count
transition_stats  <- as.data.frame(t(apply(
  transition_counts, 2,
  function(x) {
    m   <- mean(x)
    s   <- sd(x)
    c(
      Mean     = m,
      SD       = s,
      CI_Lower = m - 1.96 * s,
      CI_Upper = m + 1.96 * s,
      Min      = min(x),
      Max      = max(x)
    )
  }
)))
write.csv(
  transition_stats,
  file = file.path(outdir, "Upham_MCC_5000_ARD_TransitionStats.csv"),
  row.names = TRUE
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
  file = file.path(outdir, "Upham_MCC_5000_ARD_TimeInStateStats.csv"),
  row.names = FALSE
)

qmat <- simmaps[[1]]$model$Q
saveRDS(
  qmat,
  file = file.path(outdir, "Upham_MCC_5000_ARD_Qmatrix.rds")
)
save.image(file = file.path(outdir, "Upham_MCC_5000_ARD_Results.RData"))
