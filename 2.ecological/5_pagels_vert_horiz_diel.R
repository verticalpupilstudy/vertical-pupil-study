
#This script runs four Pagel's binary correlation tests for pupil shape and diel

if (!requireNamespace("ape",      quietly=TRUE)) install.packages("ape")
if (!requireNamespace("phytools", quietly=TRUE)) install.packages("phytools")
if (!requireNamespace("dplyr",    quietly=TRUE)) install.packages("dplyr")

library(ape)
library(phytools)
library(dplyr)


ecol <- read.csv("Mammal_PupilTraits_GenusLevel.csv", stringsAsFactors=FALSE)
tree <- read.tree("Upham2019_MCC_1000_Pruned_GenusLevel.tre")

# normalises data names with tree trips
norm <- function(x) tolower(gsub("_"," ", trimws(x)))
ecol$Species    <- norm(ecol$Upham_TipLabel)
tree$tip.label  <- norm(tree$tip.label)

ecol <- ecol %>%
  mutate(DielOrdered = factor(
    ifelse(Diel %in% c("Crepuscular","Cathemeral"), "Polyphasic", Diel),
    levels = c("Diurnal","Polyphasic","Nocturnal")
  ))

#Helper function to run 
run_test <- function(p_up, p_other, diel_pos, diel_neg, tag) {
  d <- ecol %>%
    filter(Pupil       %in% c(p_up, p_other),
           !is.na(DielOrdered),
           Species      %in% tree$tip.label) %>%
    mutate(
      Pbin = ifelse(Pupil        == p_up, 1, 0),
      Dbin = ifelse(DielOrdered %in% diel_pos, 1, 0)
    ) %>%
    filter(!is.na(Dbin))
  
  keep <- intersect(tree$tip.label, d$Species)
  tr   <- drop.tip(tree, setdiff(tree$tip.label, keep))
  
  x <- setNames(d$Pbin, d$Species)
  y <- setNames(d$Dbin, d$Species)
  
  cat("\n---", tag, "---\n")
  print( fitPagel(tr, x, y, method="fitDiscrete", model="ARD") )
}

# Diel category 1; nocturnal vs (diurnal + polyphasic)
D1_pos <- "Nocturnal";             D1_neg <- c("Diurnal","Polyphasic")
# Diel category 2;(nocturnal + polyphasic) vs diurnal
D2_pos <- c("Nocturnal","Polyphasic"); D2_neg <- "Diurnal"

# Runs four Pagel tests 

# Horizontal vs (Vertical + Round)
run_test("Horizontal", c("Vertical","Round"), D1_pos, D1_neg, "H1_HOR_vs_VERCIR_D1")
run_test("Horizontal", c("Vertical","Round"), D2_pos, D2_neg, "H2_HOR_vs_VERCIR_D2")

# Vertical vs Round
run_test("Vertical", "Round", D1_pos, D1_neg, "V1_VE_vs_ROU_D1")
run_test("Vertical", "Round", D2_pos, D2_neg, "V2_VE_vs_ROU_D2")
