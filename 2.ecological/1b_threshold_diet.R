# Genus-level diet classification with two thresholds (species and genus)
# Workflow (applies to both  below):
#  Species-level threshold: assign each species a diet category using a cutoff
#    on its diet proportions.
#  Genus-level threshold: assign each genus the majority category if that
#    category reaches the genus cutoff; otherwise label the genus "Mixed".
# Schemes:
# Four-category: Herbivore, Vertebrate carnivore, Insectivore,Mixed
#   (uses Diet.Plant, Diet.Vertebrate, Diet.Invertebrate separately)
# Three-category: Herbivore, Carnivore (vertebrate+invertebrate),Mixed
#   (combines Diet.Vertebrate + Diet.Invertebrate)
# Metrics:
#  Accuracy: share of species whose category matches their genus assignment.
#  Mixed rate: share of genera that could not be assigned a single category.

# Libraries
library(readr)
library(dplyr)
library(tidyr)

#reads species level data + drops any excluded species 
data <- read_csv("Mammal_PupilTraits_SpeciesLevel.csv", show_col_types = FALSE) %>%
  filter(is.na(Excluded) | Excluded == 0)

# Filters genera without sampled pupil
valid_pupils    <- c("Vertical", "Round", "Horizontal")
eligible_genera <- data %>%
  filter(Pupil %in% valid_pupils) %>%
  pull(Genus.1.2) %>%
  unique()
data_subset <- data %>% filter(Genus.1.2 %in% eligible_genera)

#Threshold 
species_thresholds <- c(50, 60, 70, 80)
genus_thresholds   <- c(0.5, 0.6, 0.7, 0.8)


# four_rate: 4 Diet category threshold function (Herbivore, Carnivore, Insectivore, Mixed)

classify_four <- function(df, thr) {
  df %>%
    mutate(
      Diet_Category = case_when(
        Diet.Plant        >= thr           ~ "Herbivore",
        Diet.Vertebrate   >= thr           ~ "Carnivore",
        Diet.Invertebrate >= thr           ~ "Insectivore",
        TRUE                              ~ "Mixed"
      )
    )
}

four_rate <- expand.grid(
  species_thresh = species_thresholds,
  genus_thresh   = genus_thresholds
) %>%
  rowwise() %>%
  mutate({
    sp   <- classify_four(data_subset, species_thresh)
    genus_props <- sp %>%
      group_by(Genus.1.2, Diet_Category) %>%
      summarise(n = n(), .groups="drop") %>%
      group_by(Genus.1.2) %>%
      mutate(prop = n / sum(n)) %>%
      summarise(
        Genus_Diet = if_else(
          max(prop) >= genus_thresh && sum(prop == max(prop)) == 1,
          Diet_Category[which.max(prop)],
          "Mixed"
        ),
        .groups="drop"
      )
    joined <- sp %>%
      inner_join(genus_props, by="Genus.1.2") %>%
      filter(!is.na(Diet_Category), !is.na(Genus_Diet)) %>%
      mutate(match = Diet_Category == Genus_Diet)
    n_genus <- n_distinct(genus_props$Genus.1.2)
    n_mixed <- sum(genus_props$Genus_Diet == "Mixed")
    tibble(
      match_rate = mean(joined$match, na.rm = TRUE),
      n_species  = nrow(joined),
      n_match    = sum(joined$match),
      n_genus    = n_genus,
      n_mixed    = n_mixed,
      mixed_rate = n_mixed / n_genus
    )
  }) %>%
  unnest(cols = everything())

write_csv(four_rate, "four_rate_threshold_matrix.csv")
print(four_rate)


# same but for three_rate: 3 diet categories  Insectivore + Carnivore grouped (Herbivore vs Others)

classify_three <- function(df, thr) {
  df %>%
    mutate(
      Diet_Category = case_when(
        Diet.Plant                             >= thr           ~ "Herbivore",
        (Diet.Vertebrate + Diet.Invertebrate)  >= thr           ~ "Carnivore",
        TRUE                                               ~ "Mixed"
      )
    )
}

three_rate <- expand.grid(
  species_thresh = species_thresholds,
  genus_thresh   = genus_thresholds
) %>%
  rowwise() %>%
  mutate({
    sp   <- classify_three(data_subset, species_thresh)
    genus_props <- sp %>%
      group_by(Genus.1.2, Diet_Category) %>%
      summarise(n = n(), .groups="drop") %>%
      group_by(Genus.1.2) %>%
      mutate(prop = n / sum(n)) %>%
      summarise(
        Genus_Diet = if_else(
          max(prop) >= genus_thresh && sum(prop == max(prop)) == 1,
          Diet_Category[which.max(prop)],
          "Mixed"
        ),
        .groups="drop"
      )
    joined <- sp %>%
      inner_join(genus_props, by="Genus.1.2") %>%
      filter(!is.na(Diet_Category), !is.na(Genus_Diet)) %>%
      mutate(match = Diet_Category == Genus_Diet)
    n_genus <- n_distinct(genus_props$Genus.1.2)
    n_mixed <- sum(genus_props$Genus_Diet == "Mixed")
    tibble(
      match_rate = mean(joined$match, na.rm = TRUE),
      n_species  = nrow(joined),
      n_match    = sum(joined$match),
      n_genus    = n_genus,
      n_mixed    = n_mixed,
      mixed_rate = n_mixed / n_genus
    )
  }) %>%
  unnest(cols = everything())

write_csv(three_rate, "three_rate_threshold_matrix.csv")
print(three_rate)
