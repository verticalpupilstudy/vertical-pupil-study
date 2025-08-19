#this script performs phyl. anova tests on carnivora + primate height/weight.

# Install / load required packages
if(!require(ape))      install.packages("ape");      library(ape)
if(!require(phytools)) install.packages("phytools"); library(phytools)

#reads data and tree
ecol <- read.csv("Mammal_PupilTraits_GenusLevel.csv", stringsAsFactors = FALSE)
tree <- read.tree("Upham2019_MCC_1000_Pruned_GenusLevel.tre")

tree$tip.label <- sapply(strsplit(tree$tip.label, "_"),
                         function(x) tolower(paste(x[1:2], collapse = "_")))
ecol$tip2 <- tolower(ecol$Upham_TipLabel)

# Carnivora (excl. families with horizontal pupils) height 

run_height_phylanova <- function(data_sub, nsim = 10000) {
  # keep only Round vs Vertical, drop diurnal, drop missing Height
  d <- subset(data_sub,
              Pupil %in% c("Round","Vertical") &
              tolower(Diel) != "diurnal" &
              !is.na(Height))
  # log-transform
  trait <- log(d$Height); names(trait) <- d$tip2
  # prune tree & reorder
  tips_keep <- intersect(tree$tip.label, d$tip2)
  tree_sub  <- drop.tip(tree, setdiff(tree$tip.label, tips_keep))
  d         <- d[match(tree_sub$tip.label, d$tip2), ]
  # grouping factor
  grp <- factor(d$Pupil, levels = c("Round","Vertical"))
  # runs phylo‐ANOVA
  set.seed(42)
  phylANOVA(tree_sub,
            x       = grp,
            y       = trait,
            nsim    = nsim,
            posthoc = TRUE,
            p.adj   = "holm")
}

carn <- subset(ecol, tolower(Order.1.2) == "carnivora")
res_car_height_log <- run_height_phylanova(carn)  # uses nsim=10000
cat("\nCarnivora height (log scale, 10,000 sims):\n")
print(res_car_height_log)


# Carnivora (excl. families with horizontal pupils) body mass
carn_sub <- subset(ecol,
  tolower(Order.1.2)    == "carnivora"     &
  tolower(Diel)         != "diurnal"       &
  !(Family.1.2 %in% c("Mephitidae","Mustelidae","Procyonidae")) &
  Pupil                %in% c("Round","Vertical") &
  !is.na(avg_mass_g)
)
keep_car      <- intersect(tree$tip.label, carn_sub$tip2)
tree_car      <- drop.tip(tree, setdiff(tree$tip.label, keep_car))
carn_sub      <- carn_sub[match(tree_car$tip.label, carn_sub$tip2), ]
grp_car       <- factor(carn_sub$Pupil, levels = c("Round","Vertical"))
logMass_car   <- log(carn_sub$avg_mass_g); names(logMass_car) <- carn_sub$tip2
set.seed(42)
res_car_mass <- phylANOVA(tree_car,
                          x       = grp_car,
                          y       = logMass_car,
                          nsim    = 10000,
                          posthoc = TRUE,
                          p.adj   = "holm")
cat("\n=== Carnivora mass (log) Phylo-ANOVA ===\n")
print(res_car_mass)

# Primate body mass
pri_sub <- subset(ecol,
  tolower(Order.1.2) == "primates"     &
  tolower(Diel)     != "diurnal"       &
  Pupil            %in% c("Round","Vertical") &
  !is.na(avg_mass_g)
)
keep_pri     <- intersect(tree$tip.label, pri_sub$tip2)
tree_pri     <- drop.tip(tree, setdiff(tree$tip.label, keep_pri))
pri_sub      <- pri_sub[match(tree_pri$tip.label, pri_sub$tip2), ]
grp_pri      <- factor(pri_sub$Pupil, levels = c("Round","Vertical"))
logMass_pri  <- log(pri_sub$avg_mass_g); names(logMass_pri) <- pri_sub$tip2
set.seed(42)
res_pri_mass <- phylANOVA(tree_pri,
                          x       = grp_pri,
                          y       = logMass_pri,
                          nsim    = 10000,
                          posthoc = TRUE,
                          p.adj   = "holm")
cat("\n=== Primates mass (log) Phylo-ANOVA ===\n")
print(res_pri_mass)

#Carnivora (excl. families with horizontal pupils) + primates) body mass
cp_sub <- subset(ecol,
  tolower(Order.1.2) %in% c("carnivora","primates") &
  tolower(Diel)      != "diurnal"       &
  Pupil             %in% c("Round","Vertical") &
  !is.na(avg_mass_g) &
  !(tolower(Order.1.2) == "carnivora" &
    Family.1.2 %in% c("Mephitidae","Mustelidae","Procyonidae"))
)
keep_cp     <- intersect(tree$tip.label, cp_sub$tip2)
tree_cp     <- drop.tip(tree, setdiff(tree$tip.label, keep_cp))
cp_sub      <- cp_sub[match(tree_cp$tip.label, cp_sub$tip2), ]
grp_cp      <- factor(cp_sub$Pupil, levels = c("Round","Vertical"))
logMass_cp  <- log(cp_sub$avg_mass_g); names(logMass_cp) <- cp_sub$tip2
set.seed(42)
res_cp_mass <- phylANOVA(tree_cp,
                         x       = grp_cp,
                         y       = logMass_cp,
                         nsim    = 10000,
                         posthoc = TRUE,
                         p.adj   = "holm")
cat("\n=== Carnivora + Primates mass (log) Phylo-ANOVA ===\n")
print(res_cp_mass)
