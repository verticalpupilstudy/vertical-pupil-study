#this script runs 8 phylogenetic logistic regression tests 
#for vertically slit pupils and diet (both 3 and 4 categories) 

library(phytools)
library(ape)
library(dplyr)
library(phylolm)
library(tidyr)

# Loads data and tree
ecol <- read.csv("Mammal_PupilTraits_GenusLevel.csv", stringsAsFactors = FALSE)
tree <- read.tree("Upham2019_MCC_1000_Pruned_GenusLevel.tre")

# normalises data names with tree trips
norm <- function(x) tolower(trimws(gsub("_"," ", x)))
ecol$Species   <- norm(ecol$Upham_TipLabel)
tree$tip.label <- norm(tree$tip.label)

# Exclude diurnal AND horizontal pupils 
base <- ecol %>%
  filter(Diel != "Diurnal", Pupil != "Horizontal")


# pupil contrasts
# Verticalliy slit v vertically subcircular + round 
dat271 <- base %>%
  mutate(SlitNRound = ifelse(Slit == 1, 1L, 0L))

# Keeps only vertically slit and round; drop vertical subcircular (Slit==0)
dat259 <- base %>%
  filter(Slit == 1 | Pupil == "Round") %>%
  mutate(SlitVsRound = ifelse(Slit == 1, 1L, 0L))

# Diet3: Carnivore vs Herbivore only
mk_d3 <- function(df) {
  df %>%
    filter(Genus_Diet_3 %in% c("Carnivore","Herbivore")) %>%
    mutate(Diet3Bin = ifelse(Genus_Diet_3 == "Carnivore", 1L, 0L))
}

# Diet4: Carnivore vs (Insectivore + Herbivore)
mk_d4 <- function(df) {
  df %>%
    filter(Genus_Diet_4 %in% c("Carnivore","Insectivore","Herbivore")) %>%
    mutate(Diet4Bin = ifelse(Genus_Diet_4 == "Carnivore", 1L, 0L))
}

d271_d3 <- mk_d3(dat271)
d271_d4 <- mk_d4(dat271)
d259_d3 <- mk_d3(dat259)
d259_d4 <- mk_d4(dat259)

# aligns trees
align_tree <- function(df, tree) {
  rownames(df) <- df$Species
  keep <- intersect(tree$tip.label, df$Species)
  list(
    df   = df[keep, , drop = FALSE],
    tree = drop.tip(tree, setdiff(tree$tip.label, keep))
  )
}

run_phyloglm <- function(formula, df, tree, label, slit_var) {
  al <- align_tree(df, tree)
  df_al   <- al$df
  tree_al <- al$tree
  
  vars <- all.vars(formula)
  df_al <- df_al %>% drop_na(all_of(vars))
  
  cat("\n", label, "\n", sep = "")
  cat("n =", nrow(df_al))
  if (slit_var %in% names(df_al)) {
    cat(" |", slit_var, "=1:", sum(df_al[[slit_var]] == 1, na.rm = TRUE),
        "|", slit_var, "=0:", sum(df_al[[slit_var]] == 0, na.rm = TRUE))
  }
  cat("\n")
  
  mod <- phyloglm(formula, data = df_al, phy = tree_al,
                  method = "logistic_MPLE", btol = 50, boot = 100)
  print(summary(mod)$coefficients)
  invisible(mod)
}

#Runs 8 diet tests
# 1 Diet3 ~ [Slit = 1 vs Slit = 0 + Round]  
mod_A <- run_phyloglm(Diet3Bin ~ SlitNRound, d271_d3, tree,
                      "A) Diet3 ~ SlitNRound | 271-set [Slit vs Sub+Round]",
                      "SlitNRound")

# 2 Diet3 ~ [Slit = 1 vs Round]            
mod_B <- run_phyloglm(Diet3Bin ~ SlitVsRound, d259_d3, tree,
                      "B) Diet3 ~ SlitVsRound | 259-set [Slit vs Round]",
                      "SlitVsRound")

# 3 Diet4 ~ [Slit = 1 vs Slit = 0 + Round] 
mod_C <- run_phyloglm(Diet4Bin ~ SlitNRound, d271_d4, tree,
                      "C) Diet4 ~ SlitNRound | 271-set [Slit vs Sub+Round]",
                      "SlitNRound")

# 4 Diet4 ~ [Slit = 1 vs Round]            
mod_D <- run_phyloglm(Diet4Bin ~ SlitVsRound, d259_d4, tree,
                      "D) Diet4 ~ SlitVsRound | 259-set [Slit vs Round]",
                      "SlitVsRound")

# 5 [Slit = 1 vs Slit = 0 + Round] ~ Diet3 
mod_E <- run_phyloglm(SlitNRound ~ Diet3Bin, d271_d3, tree,
                      "E) SlitNRound ~ Diet3 | 271-set [Slit vs Sub+Round]",
                      "SlitNRound")

# 6 [Slit = 1 vs Round] ~ Diet3            
mod_F <- run_phyloglm(SlitVsRound ~ Diet3Bin, d259_d3, tree,
                      "F) SlitVsRound ~ Diet3 | 259-set [Slit vs Round]",
                      "SlitVsRound")

# 7 [Slit = 1 vs Slit = 0 + Round] ~ Diet4 
mod_G <- run_phyloglm(SlitNRound ~ Diet4Bin, d271_d4, tree,
                      "G) SlitNRound ~ Diet4 | 271-set [Slit vs Sub+Round]",
                      "SlitNRound")

# 8 [Slit = 1 vs Round] ~ Diet4            
mod_H <- run_phyloglm(SlitVsRound ~ Diet4Bin, d259_d4, tree,
                      "H) SlitVsRound ~ Diet4 | 259-set [Slit vs Round]",
                      "SlitVsRound")