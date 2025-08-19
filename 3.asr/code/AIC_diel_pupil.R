# Mk model comparison (ER, SYM, ARD) on Alvarez & Upham trees for
#pupil (round, vertical, horizontal)
#diel (diurnal, nocturnal + polyphasic(crepuscular, cathemeral + mixed))

library(ape)
library(geiger)
library(dplyr)
library(readr)

# required files
traits_file <- "Mammal_PupilTraits_GenusLevel.csv"
alvarez_tree <- "4705sp_mean_pruned.nwk"
alvarez_tip  <- "Alvarez.Carretero_Tiplabel"
upham_tree   <- "Upham2019_MCC_1000_Pruned_GenusLevel.tre"
upham_tip    <- "Upham_TipLabel"

#Pupil and diel categories
valid_shapes <- c("Round","Vertical","Horizontal")
group_diel <- function(x) {
  x <- trimws(x)
  dplyr::case_when(
    x %in% c("Crepuscular","Cathemeral","Mixed") ~ "Polyphasic",
    x %in% c("Diurnal","Nocturnal")              ~ x,
    TRUE                                         ~ NA_character_
  )
}
valid_diel <- c("Nocturnal","Diurnal","Polyphasic")

# Helpers
read_traits <- function(path) readr::read_csv(path, show_col_types = FALSE)

make_states <- function(traits_df, tip_col, trait_col, valid_levels) {
  df <- traits_df %>%
    filter(!is.na(.data[[tip_col]])) %>%
    mutate(across(all_of(trait_col), ~ trimws(.))) %>%
    filter(.data[[trait_col]] %in% valid_levels)
  setNames(df[[trait_col]], df[[tip_col]])
}

align_tree_states <- function(tree, states_named) {
  keep <- intersect(tree$tip.label, names(states_named))
  tree2 <- drop.tip(tree, setdiff(tree$tip.label, keep))
  states2 <- states_named[tree2$tip.label]
  list(tree = tree2, states = states2)
}

fit_mk_set <- function(tree, states, run_id, trait_label) {
  fit_ER  <- fitDiscrete(tree, states, model = "ER")
  fit_SYM <- fitDiscrete(tree, states, model = "SYM")
  fit_ARD <- fitDiscrete(tree, states, model = "ARD")

  tibble::tibble(
    Run      = run_id,
    Trait    = trait_label,
    N_tips   = length(states),
    Model    = c("ER","SYM","ARD"),
    logLik   = c(fit_ER$opt$lnL,  fit_SYM$opt$lnL,  fit_ARD$opt$lnL),
    k        = c(fit_ER$opt$k,    fit_SYM$opt$k,    fit_ARD$opt$k),
    AIC      = c(fit_ER$opt$aic,  fit_SYM$opt$aic,  fit_ARD$opt$aic),
    AICc     = c(fit_ER$opt$aicc, fit_SYM$opt$aicc, fit_ARD$opt$aicc)
  ) %>%
    arrange(AICc) %>%
    mutate(DeltaAICc = AICc - min(AICc))
}

# Data
traits <- read_traits(traits_file) %>%
  mutate(Diel_grouped = group_diel(Diel))
alz_tree <- read.tree(alvarez_tree)
uph_tree <- read.tree(upham_tree)

# Alvarez + Pupil
states_pupil_alz <- make_states(traits, alvarez_tip, "Pupil", valid_shapes)
alzp <- align_tree_states(alz_tree, states_pupil_alz)
res1 <- fit_mk_set(alzp$tree, alzp$states, "Alvarez-Pupil",
                   "Pupil (Round/Vertical/Horizontal)")

# Alvarez + Diel (grouped)
states_diel_alz <- make_states(traits, alvarez_tip, "Diel_grouped", valid_diel)
alzd <- align_tree_states(alz_tree, states_diel_alz)
res2 <- fit_mk_set(alzd$tree, alzd$states, "Alvarez-Diel",
                   "Diel (Nocturnal/Diurnal/Polyphasic)")

# Upham + Pupil
states_pupil_uph <- make_states(traits, upham_tip, "Pupil", valid_shapes)
uphp <- align_tree_states(uph_tree, states_pupil_uph)
res3 <- fit_mk_set(uphp$tree, uphp$states, "Upham-Pupil",
                   "Pupil (Round/Vertical/Horizontal)")

# Upham + Diel (grouped)
states_diel_uph <- make_states(traits, upham_tip, "Diel_grouped", valid_diel)
uphd <- align_tree_states(uph_tree, states_diel_uph)
res4 <- fit_mk_set(uphd$tree, uphd$states, "Upham-Diel",
                   "Diel (Nocturnal/Diurnal/Polyphasic)")

# Output
all_results <- dplyr::bind_rows(res1, res2, res3, res4)
readr::write_csv(all_results, "Mk_Model_Comparisons_Alvarez_Upham.csv")