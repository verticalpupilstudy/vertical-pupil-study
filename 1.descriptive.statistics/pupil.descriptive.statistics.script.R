# Computes descriptive statistics for pupil shape data:
# counts and proportions of species, genera, families, and orders,
# plus frequency summaries for Supplementary Tables.

library(dplyr)

data_path <- "Mammal_PupilTraits_SpeciesLevel.csv"
raw       <- read.csv(data_path, stringsAsFactors = FALSE)

# Keep only eligible taxa with a recorded pupil shape (no NAs or blanks)
df0 <- raw %>%
  filter(Excluded == 0, !is.na(Pupil), Pupil != "")

# Main analysis - make Vulpes pupil = vertical, drop blank pupil shapes
df <- raw %>%
  filter(Excluded == 0) %>%
  mutate(Pupil = if_else(Genus.1.2 == "Vulpes", "Vertical", Pupil)) %>%
  filter(!is.na(Pupil), Pupil != "")

# For main analysis (excluding rodentia)
df_nr  <- df     %>% filter(Order.1.2 != "Rodentia")
elig_nr<- raw    %>% filter(Excluded == 0, Order.1.2 != "Rodentia")


# summary of eligible vs. ineligible genera 

elig_sp <- raw %>% filter(Excluded == 0) %>% distinct(Binomial.1.2) %>% nrow()
elig_ge <- raw %>% filter(Excluded == 0) %>% distinct(Genus.1.2)    %>% nrow()
elig_fa <- raw %>% filter(Excluded == 0) %>% distinct(Family.1.2)   %>% nrow()
elig_or <- raw %>% filter(Excluded == 0) %>% distinct(Order.1.2)    %>% nrow()

inelig_sp <- raw %>% filter(Excluded == 1) %>% distinct(Binomial.1.2) %>% nrow()
inelig_ge <- raw %>% filter(Excluded == 1) %>% distinct(Genus.1.2)    %>% nrow()
inelig_fa <- raw %>% filter(Excluded == 1) %>% distinct(Family.1.2)   %>% nrow()
inelig_or <- raw %>% filter(Excluded == 1) %>% distinct(Order.1.2)    %>% nrow()

pct_sp <- elig_sp / (elig_sp + inelig_sp) * 100
pct_ge <- elig_ge / (elig_ge + inelig_ge) * 100
pct_fa <- elig_fa / (elig_fa + inelig_fa) * 100
pct_or <- elig_or / (elig_or + inelig_or) * 100

cat("2. Summary of Eligible vs Ineligible Taxa\n")
cat(sprintf("  • Eligible taxa:   %d species, %d genera, %d families, %d orders\n",
            elig_sp, elig_ge, elig_fa, elig_or))
cat(sprintf("  • Ineligible taxa: %d species, %d genera, %d families, %d orders\n",
            inelig_sp, inelig_ge, inelig_fa, inelig_or))
cat(sprintf("  • Percent eligible: %.1f%% species, %.1f%% genera, %.1f%% families, %.1f%% orders\n\n",
            pct_sp, pct_ge, pct_fa, pct_or))


# overview of sampling (species level)

sampled_sp     <- df %>% distinct(Binomial.1.2) %>% nrow()
pct_sampled_sp <- sampled_sp / elig_sp * 100

cat("Species Sampling Overview\n")
cat(sprintf("  • Species with recorded pupil shape: %d of %d (%.1f%%)\n\n",
            sampled_sp, elig_sp, pct_sampled_sp))


# overview of sampling (genus level)

sampled_ge     <- df %>% distinct(Genus.1.2)  %>% nrow()
sampled_fa     <- df %>% distinct(Family.1.2) %>% nrow()
sampled_or     <- df %>% distinct(Order.1.2)  %>% nrow()

pct_ge_samp    <- sampled_ge / elig_ge * 100
pct_fa_samp    <- sampled_fa / elig_fa * 100
pct_or_samp    <- sampled_or / elig_or * 100

cat("Taxonomic Sampling Overview\n")
cat(sprintf("  • Genera sampled:  %d of %d (%.1f%%)\n", sampled_ge, elig_ge, pct_ge_samp))
cat(sprintf("  • Families sampled: %d of %d (%.1f%%)\n", sampled_fa, elig_fa, pct_fa_samp))
cat(sprintf("  • Orders sampled:   %d of %d (%.1f%%)\n\n", sampled_or, elig_or, pct_or_samp))


# genus level summaries

# Number of genera with more than 1 pupil shape
genus_shapes0 <- df0 %>%
  group_by(Genus.1.2) %>%
  summarize(n_shapes = n_distinct(Pupil), .groups="drop")
mixed_ge0 <- sum(genus_shapes0$n_shapes > 1)
total_ge0 <- nrow(genus_shapes0)

# Number of genera with 3 or more pupil shapes + how many mixed shapes
genus_species0 <- df0 %>%
  group_by(Genus.1.2) %>%
  summarize(
    n_sp     = n_distinct(Binomial.1.2),
    n_shapes = n_distinct(Pupil),
    .groups="drop"
  )
g3plus      <- filter(genus_species0, n_sp >= 3)
total3plus  <- nrow(g3plus)
mixed3plus  <- sum(g3plus$n_shapes > 1)

# Number of genera with each pupil shape
genus_by_shape <- df %>%
  group_by(Genus.1.2) %>%
  summarize(shape = first(Pupil), .groups="drop") %>%
  count(shape) %>%
  mutate(pct = n / sum(n) * 100)

cat("Genus-Level Analyses\n")
cat(sprintf("   Genera with inconsistent shapes: %d of %d\n", mixed_ge0, total_ge0))
cat(sprintf("  Genera with ≥3 species: %d (mixed: %d)\n", total3plus, mixed3plus))
cat("   Genera by pupil shape:\n"); print(genus_by_shape); cat("\n")


# how many species have each pupil shape as a percentage of total
species_shape <- df0 %>%
  count(Pupil) %>%
  mutate(pct = n / sum(n) * 100)

cat(" Species-Level Pupil Shape Prevalence\n"); print(species_shape); cat("\n")


# how many families/orders have eaach pupil shape

shape_taxa <- df0 %>%
  group_by(Pupil) %>%
  summarize(
    n_orders   = n_distinct(Order.1.2),
    n_families = n_distinct(Family.1.2),
    .groups="drop"
  )

cat(" Pupil Shape Spread Across Ranks\n"); print(shape_taxa); cat("\n")


#how manuy families have more than one shape, how many have mixed shapes

fam_shapes   <- df0 %>%
  group_by(Family.1.2) %>%
  summarize(n_shapes = n_distinct(Pupil), .groups="drop")
mixed_fa     <- sum(fam_shapes$n_shapes > 1)
exclusive_fa <- fam_shapes %>%
  filter(n_shapes == 1) %>%
  left_join(df0 %>% distinct(Family.1.2, Pupil), by="Family.1.2") %>%
  count(Pupil, name = "n_families")

cat(" Family-Level Pupil Shape Diversity\n")
cat(sprintf("  • Families with >1 shape: %d\n", mixed_fa))
cat("  • Families exclusively by shape:\n"); print(exclusive_fa); cat("\n")


# how manuy orders have more than one shape, how many have mixed shapes

order_shapes <- df0 %>%
  group_by(Order.1.2) %>%
  summarize(n_shapes = n_distinct(Pupil), .groups="drop")
mixed_or     <- sum(order_shapes$n_shapes > 1)
total_or0    <- nrow(order_shapes)
exclusive_or <- order_shapes %>%
  filter(n_shapes == 1) %>%
  left_join(df0 %>% distinct(Order.1.2, Pupil), by="Order.1.2") %>%
  count(Pupil, name = "n_orders_exclusive")

cat("Order-Level Pupil Shape Diversity\n")
cat(sprintf("  • Orders with >1 shape: %d of %d\n", mixed_or, total_or0))
cat("  • Orders exclusively by shape:\n"); print(exclusive_or); cat("\n")


#Overall datset covergage 

final_sp <- n_distinct(df$Binomial.1.2)
final_ge <- n_distinct(df$Genus.1.2)
final_fa <- n_distinct(df$Family.1.2)
final_or <- n_distinct(df$Order.1.2)

pct_ge   <- final_ge / elig_ge * 100
pct_fa   <- final_fa / elig_fa * 100

cat(" Coverage of final dataset\n")
cat(sprintf("  • %d species; %d/%d genera (%.1f%%); %d/%d families (%.1f%%); %d orders\n\n",
            final_sp, final_ge, elig_ge, pct_ge,
            final_fa, elig_fa, pct_fa,
            final_or))


# overall dataset coverage (withoud rodentia)
gen_nr     <- n_distinct(df_nr$Genus.1.2)
fam_nr     <- n_distinct(df_nr$Family.1.2)
ord_nr     <- n_distinct(df_nr$Order.1.2)

elig_ge_nr <- n_distinct(elig_nr$Genus.1.2)
elig_fa_nr <- n_distinct(elig_nr$Family.1.2)
elig_or_nr <- n_distinct(elig_nr$Order.1.2)

pct_ge_nr  <- gen_nr / elig_ge_nr * 100
pct_fa_nr  <- fam_nr / elig_fa_nr * 100

cat("coverage excluding rodentia\n")
cat(sprintf("  • %d/%d genera (%.1f%%); %d/%d families (%.1f%%); %d/%d orders\n",
            gen_nr, elig_ge_nr, pct_ge_nr,
            fam_nr, elig_fa_nr, pct_fa_nr,
            ord_nr, elig_or_nr))


