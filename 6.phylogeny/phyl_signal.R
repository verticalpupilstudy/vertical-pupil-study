
# Script that computes Pagel's  λ and Blomberg's K for pupil shapes 


library(ape)      
library(dplyr)    
library(geiger)   
library(phytools) 

tree_file <- "Upham2019_MCC_1000_Pruned_GenusLevel.tre"
tree      <- read.tree(tree_file)


traits_file <- "Mammal_PupilTraits_GenusLevel.csv"
traits      <- read.csv(traits_file, stringsAsFactors = FALSE)


tree$tip.label    <- tolower(gsub("_", " ", tree$tip.label))
traits$TipLabel   <- tolower(gsub("_", " ",
                                  trimws(traits$Upham_TipLabel)))
traits$PupilShape <- toupper(trimws(traits$Pupil))


overlap <- intersect(tree$tip.label, traits$TipLabel)



matched_traits <- traits %>%
  filter(TipLabel %in% overlap)


pupil_states <- setNames(matched_traits$PupilShape,
                         matched_traits$TipLabel)

keep        <- names(pupil_states)
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, keep))


stopifnot(all(tree_pruned$tip.label %in% names(pupil_states)))


fit_ARD_lambda <- fitDiscrete(tree_pruned,
                              pupil_states[tree_pruned$tip.label],
                              model     = "ARD",
                              transform = "lambda")
cat("Pagel's λ (three-state):",
    round(fit_ARD_lambda$opt$lambda, 3), "\n")

bin_vert <- ifelse(pupil_states[tree_pruned$tip.label] == "VERTICAL", 1, 0)
K_vert   <- phylosig(tree_pruned, bin_vert,
                     method = "K", test = TRUE, nsim = 1000)
cat("Vertical pupils - Blomberg's K:",
    round(K_vert$K, 3), "; P =", signif(K_vert$P, 3), "\n")

lambda_vert <- phylosig(tree_pruned, bin_vert,
                        method = "lambda", test = TRUE)
cat("Vertical pupils - Pagel's λ:",
    round(lambda_vert$lambda, 3),
    "; LRT P =", signif(lambda_vert$P, 3), "\n")

bin_horiz <- ifelse(pupil_states[tree_pruned$tip.label] == "HORIZONTAL", 1, 0)

K_horiz <- phylosig(tree_pruned, bin_horiz, method = "K", test = TRUE, nsim = 1000)
cat("Horizontal pupils - Blomberg's K:",
    round(K_horiz$K, 3), "; P =", signif(K_horiz$P, 3), "\n")

lambda_horiz <- phylosig(tree_pruned, bin_horiz, method = "lambda", test = TRUE)
cat("Horizontal pupils - Pagel's λ:",
    round(lambda_horiz$lambda, 3), "; LRT P =", signif(lambda_horiz$P, 3), "\n")
bin_round <- ifelse(pupil_states[tree_pruned$tip.label] == "ROUND", 1, 0)

K_round <- phylosig(tree_pruned, bin_round, method = "K", test = TRUE, nsim = 1000)
cat("Round pupils - Blomberg's K:",
    round(K_round$K, 3), "; P =", signif(K_round$P, 3), "\n")

lambda_round <- phylosig(tree_pruned, bin_round, method = "lambda", test = TRUE)
cat("Round pupils - Pagel's λ:",
    round(lambda_round$lambda, 3), "; LRT P =", signif(lambda_round$P, 3), "\n")

