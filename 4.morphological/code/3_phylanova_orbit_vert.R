
#Install / load required packages

if(!require(dplyr))    install.packages("dplyr");    library(dplyr)
if(!require(ape))      install.packages("ape");      library(ape)
if(!require(phytools)) install.packages("phytools"); library(phytools)

#reads data + tree

df   <- read.csv("Mammal_PupilTraits_GenusLevel.csv", stringsAsFactors = FALSE)
tree <- read.tree("Upham2019_MCC_1000_Pruned_GenusLevel.tre")

tree$tip.label <- sapply(strsplit(tree$tip.label, "_"),
                         function(x) tolower(paste(x[1:2], collapse = "_")))
df$tip2 <- tolower(df$Upham_TipLabel)


run_vert <- function(data_sub, nsim = 10000) {
  d <- data_sub %>%
    filter(Pupil %in% c("Round","Vertical"),
           !is.na(Orbit_Verticality))
  
  keep_tips <- intersect(tree$tip.label, d$tip2)
  tree_sub  <- drop.tip(tree, setdiff(tree$tip.label, keep_tips))
  d         <- d[match(tree_sub$tip.label, d$tip2), ]
  
  grp <- factor(d$Pupil, levels = c("Round","Vertical"))
  y   <- d$Orbit_Verticality; names(y) <- d$tip2
  
  set.seed(42)
  phylANOVA(tree_sub, grp, y, nsim = nsim, posthoc = TRUE, p.adj = "holm")
}

# Runs  phyl anova tests

# (1) All species with Orbit_Verticality
res_all      <- run_vert(df)

# (2a) Exclude Height > 500 cm 
res_h500     <- run_vert(df %>% filter(is.na(Height) | Height <= 500))

# (2b) Exclude Height > 450 cm 
res_h450     <- run_vert(df %>% filter(is.na(Height) | Height <= 450))

# (3a) Primates only
res_primate    <- run_vert(df %>% filter(tolower(Order.1.2) == "primates"))

# (3b) Carnivora only
res_carnivora  <- run_vert(df %>% filter(tolower(Order.1.2) == "carnivora"))

# (3c) Diprotodontia + Dasyuromorphia only
res_dip_dasy   <- run_vert(df %>% filter(tolower(Order.1.2) %in% c("diprotodontia","dasyuromorphia")))

# (3d) Rodentia only
res_rodentia   <- run_vert(df %>% filter(tolower(Order.1.2) == "rodentia"))


# 5) Print results

cat("\n(1) All species:\n");                     print(res_all)
cat("\n(2a) Exclude Height > 500 cm:\n");       print(res_h500)
cat("\n(2b) Exclude Height > 450 cm:\n");       print(res_h450)
cat("\n(3a) Primates:\n");                      print(res_primate)
cat("\n(3b) Carnivora:\n");                     print(res_carnivora)
cat("\n(3c) Diprotodontia + Dasyuromorphia:\n"); print(res_dip_dasy)
cat("\n(3d) Rodentia:\n");                     print(res_rodentia)
