
#this script runs 12 phylogenetic logistic regression tests 
#for horizontal pupils and diel activity

library(phytools)
library(ape)
library(dplyr)
library(phylolm)

#load traits and tree
ecol <- read.csv("Mammal_PupilTraits_GenusLevel.csv", stringsAsFactors=FALSE)
tree <- read.tree("Upham2019_MCC_1000_Pruned_GenusLevel.tre")

# normalises data names with tree trips
norm <- function(x) tolower(trimws(gsub("_"," ", x)))
ecol$Species   <- norm(ecol$Upham_TipLabel)
tree$tip.label <- norm(tree$tip.label)

# make two datasets: broad = includes horiz, narrow = just round + vertical
broadH  <- ecol %>% filter(Pupil %in% c("Horizontal","Round","Vertical"))
narrowH <- ecol %>% filter(Pupil %in% c("Horizontal","Round"))
# add binary columns for tests
broadH <- broadH %>%
  mutate(
    PupilH1 = ifelse(Pupil=="Horizontal", 1L, 0L),
    PupilH2 = ifelse(Pupil=="Horizontal", 1L, 0L),
    DielA   = ifelse(Diel=="Diurnal", 1L,
                     ifelse(Diel=="Nocturnal", 0L, NA_integer_)),                          
    DielB   = ifelse(Diel %in% c("Diurnal","Cathemeral","Crepuscular"), 1L,
                     ifelse(Diel=="Nocturnal", 0L, NA_integer_)),                          
    DielC   = ifelse(Diel=="Diurnal", 1L,
                     ifelse(Diel %in% c("Cathemeral","Crepuscular","Nocturnal"), 0L, NA_integer_))
  )

narrowH <- narrowH %>%
  mutate(
    PupilH2 = ifelse(Pupil=="Horizontal", 1L, 0L),
    DielA   = ifelse(Diel=="Diurnal", 1L,
                     ifelse(Diel=="Nocturnal", 0L, NA_integer_)),
    DielB   = ifelse(Diel %in% c("Diurnal","Cathemeral","Crepuscular"), 1L,
                     ifelse(Diel=="Nocturnal", 0L, NA_integer_)),
    DielC   = ifelse(Diel=="Diurnal", 1L,
                     ifelse(Diel %in% c("Cathemeral","Crepuscular","Nocturnal"), 0L, NA_integer_))
  )

# Helper to run PLR
run_test <- function(df, formula, tree, label) {
  df2 <- df %>% drop_na(all.vars(formula))
  rownames(df2) <- df2$Species
  pr <- drop.tip(tree, setdiff(tree$tip.label, rownames(df2)))
  cat("\n", label, "— n =", nrow(df2), "\n")
  fit <- phyloglm(formula, data=df2, phy=pr,
                  method="logistic_MPLE", btol=50, boot=100)
  print(summary(fit)$coefficients)
}

#run the 12 tests 


# 1 DielA ~ PupilH1
run_test(broadH,  DielA ~ PupilH1, tree,
         "1: Diel [Diurnal vs Nocturnal] ~ Pupil [Horizontal vs (Round+Vert)]")

# 2 DielA ~ PupilH2
run_test(narrowH, DielA ~ PupilH2, tree,
         "2: Diel [Diurnal vs Nocturnal] ~ Pupil [Horizontal vs Round]")

# 3 DielB ~ PupilH1
run_test(broadH,  DielB ~ PupilH1, tree,
         "3: Diel [(Diurnal+Cathemeral+Crepuscular) vs Nocturnal] ~ Pupil [Horizontal vs (Round+Vert)]")

# 4 DielB ~ PupilH2
run_test(narrowH, DielB ~ PupilH2, tree,
         "4: Diel [(Diurnal+Cathemeral+Crepuscular) vs Nocturnal] ~ Pupil [Horizontal vs Round]")

# 5 DielC ~ PupilH1
run_test(broadH,  DielC ~ PupilH1, tree,
         "5: Diel [Diurnal vs (Cathemeral+Crepuscular+Nocturnal)] ~ Pupil [Horizontal vs (Round+Vert)]")

# 6 DielC ~ PupilH2
run_test(narrowH, DielC ~ PupilH2, tree,
         "6: Diel [Diurnal vs (Cathemeral+Crepuscular+Nocturnal)] ~ Pupil [Horizontal vs Round]")


# 7 PupilH1 ~ DielA
run_test(broadH,  PupilH1 ~ DielA, tree, 
         "7: Pupil [Horizontal vs (Round+Vert)] ~ Diel [Diurnal vs Nocturnal]")

# 8 PupilH2 ~ DielA
run_test(narrowH, PupilH2 ~ DielA, tree,
         "8: Pupil [Horizontal vs Round] ~ Diel [Diurnal vs Nocturnal]")

# 9 PupilH1 ~ DielB
run_test(broadH,  PupilH1 ~ DielB, tree,
         "9: Pupil [Horizontal vs (Round+Vert)] ~ Diel [(Diurnal+Cathemeral+Crepuscular) vs Nocturnal]")

# 10 PupilH2 ~ DielB
run_test(narrowH, PupilH2 ~ DielB, tree,
         "10: Pupil [Horizontal vs Round] ~ Diel [(Diurnal+Cathemeral+Crepuscular) vs Nocturnal]")

# 11 PupilH1 ~ DielC
run_test(broadH,  PupilH1 ~ DielC, tree,
         "11: Pupil [Horizontal vs (Round+Vert)] ~ Diel [Diurnal vs (Cathemeral+Crepuscular+Nocturnal)]")

# 12 PupilH2 ~ DielC
run_test(narrowH, PupilH2 ~ DielC, tree,
         "12: Pupil [Horizontal vs Round] ~ Diel [Diurnal vs (Cathemeral+Crepuscular+Nocturnal)]")