#this script runs 12 phylogenetic logistic regression tests 
#for vertically slit pupils and diel activity 


library(phytools); library(ape); library(dplyr); library(phylolm); library(tidyr)

# Loads data and tree
ecol <- read.csv("Mammal_PupilTraits_GenusLevel.csv", stringsAsFactors=FALSE)
tree <- read.tree("Upham2019_MCC_1000_Pruned_GenusLevel.tre")

# normalises data names with tree trips
normalize <- function(x) {
  x <- gsub("_"," ", x)
  tolower(trimws(x))
}
ecol$Species   <- normalize(ecol$Upham_TipLabel)
tree$tip.label <- normalize(tree$tip.label)

# make some new variables, factors etc so it’s easier to filter later
ecol <- ecol %>%
  mutate(
    PupilShape = factor(Pupil, levels = c("Vertical", "Round", "Horizontal")),
    SlitGroup = case_when(
      Pupil == "Vertical" & Slit == 1 ~ "Vertical_Slit",
      Pupil == "Vertical" & Slit == 0 ~ "Vertical_Subcirc",
      TRUE ~ "Other"
    ),
    is_slit = SlitGroup == "Vertical_Slit",
    is_subcirc_or_round = SlitGroup != "Vertical_Slit",
    is_round = Pupil == "Round",
    is_strict_noct = Diel == "Nocturnal" & (Flexible == 0 | is.na(Flexible)),
    is_flex_noct = Diel == "Nocturnal" & Flexible == 1,
    is_crep_cat = Diel %in% c("Crepuscular", "Cathemeral"),
    is_nocturnal = Diel == "Nocturnal"
  )
#wrapper to run phylo regression 
run_phylo_print <- function(formula, data, tree, label, number) {
  keep <- intersect(tree$tip.label, data$Species)
  pruned <- drop.tip(tree, setdiff(tree$tip.label, keep))
  if(length(unique(eval(formula[[2]], data))) < 2) {
    message("\n== Test ", number, ": ", label, " ==\nSkipped: response has no variation")
    return()
  }
  rownames(data) <- data$Species
  message("\n== Test ", number, ": ", label, " ==")
  cat("N =", nrow(data), "\n")
  fit <- phyloglm(formula, data=data, phy=pruned,
                  method="logistic_MPLE", btol=50, boot=100)
  print(summary(fit)$coefficients)
}

# Define all 12 tests
tests <- list(
  list("Slit vs Subcirc+Round ~ Strict vs Flexible", 
       quote(Response ~ Predictor),
       ecol %>% filter(is_slit | is_subcirc_or_round, is_strict_noct | is_flex_noct) %>%
         mutate(Response = as.integer(is_slit), Predictor = as.integer(is_flex_noct))),
  
  list("Slit vs Subcirc+Round ~ Strict vs Flex+Crep+Cat", 
       quote(Response ~ Predictor),
       ecol %>% filter(is_slit | is_subcirc_or_round, is_strict_noct | is_flex_noct | is_crep_cat) %>%
         mutate(Response = as.integer(is_slit), Predictor = as.integer(!is_strict_noct))),
  
  list("Slit vs Subcirc+Round ~ Noct vs Crep+Cat", 
       quote(Response ~ Predictor),
       ecol %>% filter(is_slit | is_subcirc_or_round, is_nocturnal | is_crep_cat) %>%
         mutate(Response = as.integer(is_slit), Predictor = as.integer(is_crep_cat))),
  
  list("Slit vs Round ~ Strict vs Flexible", 
       quote(Response ~ Predictor),
       ecol %>% filter(is_slit | is_round, is_strict_noct | is_flex_noct) %>%
         mutate(Response = as.integer(is_slit), Predictor = as.integer(is_flex_noct))),
  
  list("Slit vs Round ~ Strict vs Flex+Crep+Cat", 
       quote(Response ~ Predictor),
       ecol %>% filter(is_slit | is_round, is_strict_noct | is_flex_noct | is_crep_cat) %>%
         mutate(Response = as.integer(is_slit), Predictor = as.integer(!is_strict_noct))),
  
  list("Slit vs Round ~ Noct vs Crep+Cat", 
       quote(Response ~ Predictor),
       ecol %>% filter(is_slit | is_round, is_nocturnal | is_crep_cat) %>%
         mutate(Response = as.integer(is_slit), Predictor = as.integer(is_crep_cat))),
  
  list("Strict vs Flexible ~ Slit vs Subcirc+Round", 
       quote(Response ~ Predictor),
       ecol %>% filter(is_slit | is_subcirc_or_round, is_strict_noct | is_flex_noct) %>%
         mutate(Response = as.integer(is_strict_noct), Predictor = as.integer(is_slit))),
  
  list("Strict vs Flex+Crep+Cat ~ Slit vs Subcirc+Round", 
       quote(Response ~ Predictor),
       ecol %>% filter(is_slit | is_subcirc_or_round, is_strict_noct | is_flex_noct | is_crep_cat) %>%
         mutate(Response = as.integer(is_strict_noct), Predictor = as.integer(is_slit))),
  
  list("Noct vs Crep+Cat ~ Slit vs Subcirc+Round", 
       quote(Response ~ Predictor),
       ecol %>% filter(is_slit | is_subcirc_or_round, is_nocturnal | is_crep_cat) %>%
         mutate(Response = as.integer(is_nocturnal), Predictor = as.integer(is_slit))),
  
  list("Strict vs Flexible ~ Slit vs Round", 
       quote(Response ~ Predictor),
       ecol %>% filter(is_slit | is_round, is_strict_noct | is_flex_noct) %>%
         mutate(Response = as.integer(is_strict_noct), Predictor = as.integer(is_slit))),
  
  list("Strict vs Flex+Crep+Cat ~ Slit vs Round", 
       quote(Response ~ Predictor),
       ecol %>% filter(is_slit | is_round, is_strict_noct | is_flex_noct | is_crep_cat) %>%
         mutate(Response = as.integer(is_strict_noct), Predictor = as.integer(is_slit))),
  
  list("Noct vs Crep+Cat ~ Slit vs Round", 
       quote(Response ~ Predictor),
       ecol %>% filter(is_slit | is_round, is_nocturnal | is_crep_cat) %>%
         mutate(Response = as.integer(is_nocturnal), Predictor = as.integer(is_slit)))
)

# Run all tests in loop 
for (i in seq_along(tests)) {
  label <- tests[[i]][[1]]
  formula <- tests[[i]][[2]]
  df <- tests[[i]][[3]]
  run_phylo_print(formula, df, tree, label, i)
}