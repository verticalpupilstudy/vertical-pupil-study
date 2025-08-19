## This script runs stochastic character mapping (SCM) on diel activity
# using the Alvarez–Carretero phylogeny (4705sp_mean).
# nsim = 5000, model = ARD, root prior = "fitzjohn"

# required libraries
for (pkg in c("ape","phytools","dplyr","stringr","tibble","readr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(ape)
library(phytools)
library(dplyr)
library(stringr)
library(tibble)
library(readr)

# required files
tree_file   <- "4705sp_mean_pruned.nwk"
traits_file <- "Mammal_PupilTraits_GenusLevel.csv"
tip_col     <- "Alvarez.Carretero_Tiplabel"
outdir      <- "outputs_alvarez"

# output directory
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# read tree and trait data
tree      <- read.tree(tree_file)
traits_df <- read_csv(traits_file, show_col_types = FALSE)

# group diel categories (Crepuscular + Cathemeral + Mixed -> "Polyphasic")
traits_df <- traits_df %>%
  mutate(
    Diel_grouped = case_when(
      Diel %in% c("Crepuscular","Cathemeral","Mixed") ~ "Polyphasic",
      TRUE                                            ~ Diel
    )
  )

# keep only species in the tree with valid diel states
valid_states <- c("Nocturnal","Diurnal","Polyphasic")
traits_df <- traits_df %>%
  filter(.data[[tip_col]] %in% tree$tip.label,
         Diel_grouped %in% valid_states)

# match states to tips and prune tree
diel_states <- setNames(traits_df$Diel_grouped, traits_df[[tip_col]])
to_drop     <- setdiff(tree$tip.label, names(diel_states))
if (length(to_drop) > 0) {
  tree <- drop.tip(tree, to_drop)
}

# generate SCM with ARD and nsim = 5000, pi = "fitzjohn"
nsim <- 5000
simmaps <- make.simmap(
  tree  = tree,
  x     = diel_states,
  model = "ARD",
  nsim  = nsim,
  pi    = "fitzjohn"
)

# summaries: ancestral probabilities, transition counts, time in states, Q matrix
sim_desc        <- describe.simmap(simmaps)

ace_probs <- sim_desc$ace
write.csv(
  ace_probs,
  file      = file.path(outdir, "ASR_Diel_5000_ARD_AncestralProbs.csv"),
  row.names = TRUE
)

transition_counts <- sim_desc$count
transition_stats  <- as.data.frame(t(apply(
  transition_counts, 2,
  function(x) {
    m <- mean(x); s <- sd(x)
    c(Mean = m,
      SD = s,
      CI_Lower = m - 1.96*s,
      CI_Upper = m + 1.96*s,
      Min = min(x),
      Max = max(x))
  }
)))
write.csv(
  transition_stats,
  file      = file.path(outdir, "ASR_Diel_5000_ARD_TransitionStats.csv"),
  row.names = TRUE
)

time_df <- as.data.frame(sim_desc$times)
time_stats <- time_df %>%
  summarise(across(
    everything(),
    list(
      Mean    = ~ mean(.),
      SD      = ~ sd(.),
      CI_Low  = ~ mean(.) - 1.96*sd(.),
      CI_High = ~ mean(.) + 1.96*sd(.)
    ),
    .names = "{col}_{fn}"
  ))
write.csv(
  time_stats,
  file      = file.path(outdir, "ASR_Diel_5000_ARD_TimeInStateStats.csv"),
  row.names = FALSE
)

saveRDS(
  simmaps[[1]]$model$Q,
  file = file.path(outdir, "ASR_Diel_5000_ARD_Qmatrix.rds")
)
saveRDS(
  simmaps,
  file = file.path(outdir, "ASR_Diel_5000_ARD_Simmaps.rds")
)

save.image(file = file.path(outdir, "ASR_Diel_5000_ARD_Results.RData"))