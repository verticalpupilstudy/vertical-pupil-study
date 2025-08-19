#this script runs 10 phylogenetic logistic regression tests 
#for vertical pupils and diel activity 

library(phytools)
library(ape)
library(dplyr)
library(phylolm)
library(tidyr)

ecol <- read.csv("Mammal_PupilTraits_GenusLevel.csv", stringsAsFactors=FALSE)
tree <- read.tree("Upham2019_MCC_1000_Pruned_GenusLevel.tre")

# normalises data names with tree trips
norm <- function(x) tolower(trimws(gsub("_"," ", x)))
ecol$Species   <- norm(ecol$Upham_TipLabel)
tree$tip.label <- norm(tree$tip.label)

# split dataset into broader (all pupil shapes) and narrower (just round vs vertical)
broad  <- ecol %>% filter(Pupil %in% c("Vertical","Round","Horizontal"))
narrow <- ecol %>% filter(Pupil %in% c("Vertical","Round"))

# add binary variables for tests
# PupilBin is vertical vs not, DielBin has different cutoffs for nocturnality
broad <- broad %>%
  mutate(
    PupilBin1 = ifelse(Pupil=="Vertical", 1L, 0L),
    PupilBin2 = ifelse(Pupil=="Vertical", 1L, 0L),  
    DielBin1  = ifelse(Diel %in% c("Nocturnal","Crepuscular","Cathemeral"), 1L, 0L),  
    DielBin2  = ifelse(Diel=="Nocturnal", 1L, 0L), 
    DielBin3  = ifelse(Diel=="Nocturnal", 1L,  
                       ifelse(Diel %in% c("Crepuscular","Cathemeral"), 0L, NA_integer_))
  )
  # same but for the narrower dataset

narrow <- narrow %>%
  mutate(
    PupilBin2 = ifelse(Pupil=="Vertical", 1L, 0L),
    DielBin1  = ifelse(Diel %in% c("Nocturnal","Crepuscular","Cathemeral"), 1L, 0L),
    DielBin2  = ifelse(Diel=="Nocturnal",1L,0L),
    DielBin3  = ifelse(Diel=="Nocturnal",1L,
                       ifelse(Diel %in% c("Crepuscular","Cathemeral"),0L,NA_integer_))
  )

#Helper to run the test
run_test <- function(df, formula, tree, label) {
  df2 <- df %>% drop_na(all.vars(formula))
  rownames(df2) <- df2$Species
  pr <- drop.tip(tree, setdiff(tree$tip.label, rownames(df2)))
  cat("\n", label, "— n =", nrow(df2), "\n")
  fit <- phyloglm(formula, data=df2, phy=pr,
                  method="logistic_MPLE", btol=50, boot=100)
  print(summary(fit)$coefficients)
}

# Run the 10 tests:

# 1. PupilBin1 ~ DielBin1
run_test(broad,  PupilBin1 ~ DielBin1, tree,
         "1: Pupil [Vert vs (Rnd+H)] ~ Diel [NC+C vs Diurnal]")

# 2. PupilBin2 ~ DielBin1
run_test(narrow, PupilBin2 ~ DielBin1, tree,
         "2: Pupil [Vert vs Round] ~ Diel [NC+C vs Diurnal]")

# 3. PupilBin1 ~ DielBin2
run_test(broad,  PupilBin1 ~ DielBin2, tree,
         "3: Pupil [Vert vs (Rnd+H)] ~ Diel [Nocturnal vs Others]")

# 4. PupilBin2 ~ DielBin2
run_test(narrow, PupilBin2 ~ DielBin2, tree,
         "4: Pupil [Vert vs Round] ~ Diel [Nocturnal vs Others]")

# 5. DielBin1 ~ PupilBin1
run_test(broad,  DielBin1 ~ PupilBin1, tree,
         "5: Diel [NC+C vs Diurnal] ~ Pupil [Vert vs (Rnd+H)]")

# 6. DielBin1 ~ PupilBin2
run_test(narrow, DielBin1 ~ PupilBin2, tree,
         "6: Diel [NC+C vs Diurnal] ~ Pupil [Vert vs Round]")

# 7. DielBin2 ~ PupilBin1
run_test(broad,  DielBin2 ~ PupilBin1, tree,
         "7: Diel [Nocturnal vs Others] ~ Pupil [Vert vs (Rnd+H)]")

# 8. DielBin2 ~ PupilBin2
run_test(narrow, DielBin2 ~ PupilBin2, tree,
         "8: Diel [Nocturnal vs Others] ~ Pupil [Vert vs Round]")

# 9. PupilBin2 ~ DielBin3
run_test(narrow, PupilBin2 ~ DielBin3, tree,
         "9: Pupil [Vert vs Round] ~ Diel [Nocturnal vs (C+C)]")
# 10. Diel [Nocturnal vs (C+C)] ~ Pupil [Vert vs Round]
run_test(narrow, DielBin3 ~ PupilBin2, tree,
         "10: Diel [Nocturnal vs (C+C)] ~ Pupil [Vert vs Round]")