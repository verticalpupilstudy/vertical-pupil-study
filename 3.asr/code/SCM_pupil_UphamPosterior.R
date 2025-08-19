## Stochastic character mapping (SCM) on a 1000 tree sample 
#from the Upham posterior distribution 
# runs 500 simulations on each tree
## Set the tree block you want to run depending on computational power
#(e.g., 1:1000, 1:250, 251:500, etc.). [we used 250 tree batches]
#matches posterior nodes to Upham MCC 

# load packages
if (!requireNamespace("ape", quietly=TRUE))     install.packages("ape")
if (!requireNamespace("phytools", quietly=TRUE))install.packages("phytools")
if (!requireNamespace("openxlsx", quietly=TRUE))install.packages("openxlsx")

library(ape)
library(phytools)
library(openxlsx)

# quick helper: count transitions in a single simmap
count_transitions <- function(smap) {
  states <- unique(unlist(lapply(smap$maps, names)))
  M <- matrix(0, nrow=length(states), ncol=length(states),
              dimnames=list(states, states))
  for (seg in smap$maps) {
    st <- names(seg)
    if (length(st) > 1) {
      for (i in 2:length(st)) {
        M[ st[i-1], st[i] ] <- M[ st[i-1], st[i] ] + 1
      }
    }
  }
  M
}

# input files (Upham posterior + MCC consensus to match )
posterior_trees <- read.tree("sampled1000_pruned_posterior.tre")
consensus_tree  <- read.tree("Upham2019_MCC_1000_Pruned_GenusLevel.tre")

# load + align trait data to consensus tips
cat("Loading & aligning trait data...\n")
td <- read.csv("Mammal_PupilTraits_GenusLevel.csv", stringsAsFactors=FALSE)
if (!all(c("Upham_TipLabel","Pupil") %in% names(td)))
  stop("CSV needs 'Upham_TipLabel' and 'Pupil' columns")

valid_states <- c("Vertical","Round","Horizontal")
td <- subset(td, Pupil %in% valid_states & !is.na(Pupil))
td <- subset(td, Upham_TipLabel %in% consensus_tree$tip.label)
aligned_traits <- setNames(td$Pupil, td$Upham_TipLabel)

# choose the block to run
start_tree <- 1      # <- adjust this
end_tree   <- 250    # <- adjust this
if (start_tree < 1 || end_tree > length(posterior_trees) || start_tree > end_tree)
  stop("tree range must be within 1..", length(posterior_trees), " and start <= end")

tree_range <- start_tree:end_tree

# sim 
nsim <- 500   # number of stochastic maps per tree

# output folder + filename based on range
outdir    <- sprintf("outputs_upham_%d-%d", start_tree, end_tree)
outprefix <- sprintf("UphamPosterior_Tree%d_to_%d", start_tree, end_tree)
if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)

# accumulators
combined_probs        <- NULL
all_transition_counts <- data.frame()
q_matrix_list         <- list()

# loop over chosen posterior trees
for (tree_idx in tree_range) {
  cat("Tree", tree_idx, "— running", nsim, "simulations...\n")
  tree <- posterior_trees[[tree_idx]]

  # keep only tips with valid pupil 
  keep <- intersect(tree$tip.label, names(aligned_traits))
  if (length(keep) < 2) next
  if (length(keep) < length(tree$tip.label))
    tree <- drop.tip(tree, setdiff(tree$tip.label, keep))

  # run SCM (ARD, FitzJohn root prior)
  smaps <- make.simmap(tree, aligned_traits[keep], model="ARD",
                       nsim=nsim, pi="fitzjohn")

  # q-matrix
  q_matrix_list[[as.character(tree_idx)]] <- smaps[[1]]$Q

  # posterior node probabilities
  desc_all <- describe.simmap(smaps)
  pp       <- as.data.frame(desc_all$ace)
  pp$Posterior_Node <- rownames(pp)

  # match posterior nodes to MCC consensus
  mn <- matchNodes(consensus_tree, tree, method="descendants")
  md <- setNames(as.data.frame(mn), c("Consensus_Node","Posterior_Node"))
  md$Posterior_Node <- as.character(md$Posterior_Node)

  merged <- merge(md, pp, by="Posterior_Node", all.x=TRUE)

  # rename state cols with tree suffix
  ace_cols <- intersect(valid_states, colnames(merged))
  new_cols <- paste0(ace_cols, "_Tree", tree_idx)
  names(merged)[ match(ace_cols, names(merged)) ] <- new_cols

  keep_cols <- c("Consensus_Node",
                 paste0("Vertical_Tree",   tree_idx),
                 paste0("Round_Tree",      tree_idx),
                 paste0("Horizontal_Tree", tree_idx))
  merged2 <- merged[, intersect(keep_cols, colnames(merged)), drop=FALSE]

  combined_probs <- if (is.null(combined_probs)) merged2 else
                    merge(combined_probs, merged2, by="Consensus_Node", all=TRUE)

  # per-sim transitions
  for (sim_i in seq_along(smaps)) {
    M  <- count_transitions(smaps[[sim_i]])
    df <- as.data.frame(as.table(M))
    names(df) <- c("From","To","Count")
    df <- subset(df, From != To)
    df$Tree       <- tree_idx
    df$Simulation <- sim_i
    all_transition_counts <- rbind(all_transition_counts, df)
  }
}

# Outputs saved and labelled sampled tree range
combined_probs <- combined_probs[order(as.numeric(combined_probs$Consensus_Node)), ]

write.xlsx(combined_probs,
           file.path(outdir, paste0(outprefix, "_Posterior_Probs.xlsx")),
           rowNames=FALSE)
save(combined_probs,
     file=file.path(outdir, paste0(outprefix, "_Posterior_Probs.RData")))

write.xlsx(all_transition_counts,
           file.path(outdir, paste0(outprefix, "_All_Transition_Counts.xlsx")),
           rowNames=FALSE)

save(q_matrix_list,
     file=file.path(outdir, paste0(outprefix, "_Q_Matrices.RData")))

save(combined_probs, all_transition_counts, q_matrix_list,
     file=file.path(outdir, paste0(outprefix, "_Full_ASR_Results.RData")))

cat("Done. Wrote outputs to: ", outdir, "\n", sep = "")