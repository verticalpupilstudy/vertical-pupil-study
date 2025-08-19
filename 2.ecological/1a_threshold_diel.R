# Script for genus-level classification of diel activity and nocturnal flexibility
# using threshold-based rules. Genus-level states are assigned from the proportions
# of species in each category, and then compared against the observed states of
# individual species to evaluate accuracy.
#
# Accuracy: proportion of species correctly classified.
# Mixed rate: proportion of genera that could not be assigned a single category
#             at the given threshold.

library(readr)
library(dplyr)
library(purrr)

species_csv <- "Mammal_PupilTraits_SpeciesLevel.csv"

# Threshold 
diel_T1_grid <- seq(0.50, 1.00, by = 0.05)  # main threshold for diel
flex_T_grid  <- seq(0.50, 1.00, by = 0.05)  # threshold for nocturnal subcategory
diel_T2_fallback <- 0.50                    # fallback for assigning Cathemeral when 
                                              #both diunral + nocturnal exceed this proportion

# Outputs
diel_out_csv <- "diel_threshold_analysis_updated.csv"
flex_out_csv <- "genus_flex_threshold_analysis.csv"


dat <- read_csv(species_csv, show_col_types = FALSE) %>%
  filter(Excluded == 0)

# Exclude genera without a sampled pupil shape
valid_pupils <- c("Vertical", "Round", "Horizontal")
genera_with_pupil <- dat %>%
  filter(Pupil %in% valid_pupils) %>%
  pull(Genus.1.2) %>%
  unique()

dat_genus <- dat %>% filter(Genus.1.2 %in% genera_with_pupil)


# Genus diel proportions
genus_props <- dat_genus %>%
  group_by(Genus.1.2) %>%
  summarise(
    total_species = n(),
    p_nocturnal   = mean(Diel_Cox == "Nocturnal",   na.rm = TRUE),
    p_diurnal     = mean(Diel_Cox == "Diurnal",     na.rm = TRUE),
    p_crepuscular = mean(Diel_Cox == "Crepuscular", na.rm = TRUE),
    p_cathemeral  = mean(Diel_Cox == "Cathemeral",  na.rm = TRUE),
    .groups = "drop"
  )
#function to assign genus diel at a given Threshold and evaluate accuracy 
assign_diel_at_T1 <- function(T1, T2 = diel_T2_fallback) {

  genus_labels <- genus_props %>%
    mutate(
      assigned_diel = case_when(
        abs(p_diurnal - p_nocturnal) < 1e-6 & p_diurnal > 0 & p_nocturnal > 0 ~ "Cathemeral",
        p_nocturnal   >= T1 ~ "Nocturnal",
        p_diurnal     >= T1 ~ "Diurnal",
        p_crepuscular >= T1 ~ "Crepuscular",
        p_cathemeral  >= T1 ~ "Cathemeral",
        (p_nocturnal + p_diurnal) >= T2 ~ "Cathemeral",
        TRUE ~ "Mixed"
      )
    )

  joined <- dat_genus %>%
    left_join(genus_labels %>% select(Genus.1.2, assigned_diel), by = "Genus.1.2") %>%
    filter(assigned_diel != "Mixed")

  tibble(
    T1         = T1,
    accuracy   = mean(joined$Diel_Cox == joined$assigned_diel, na.rm = TRUE),
    mixed_rate = mean(genus_labels$assigned_diel == "Mixed")
  )
}
#RUN AND SAVE
diel_results <- map_dfr(diel_T1_grid, assign_diel_at_T1)
write_csv(diel_results, diel_out_csv)

#analysis on nocturnal subcategory
nocturnal_only <- dat_genus %>% filter(Diel_Cox == "Nocturnal")

genus_flex <- nocturnal_only %>%
  filter(!is.na(Flexible_Cox)) %>%
  group_by(Genus.1.2) %>%
  summarise(
    n_nocturnal   = n(),
    prop_flexible = mean(Flexible_Cox == 1),
    .groups = "drop"
  )
#function to assign genus diel (nocturnal flex) at a given Threshold and evaluate accuracy 

evaluate_flex_at_T <- function(T) {
  
  genus_labels <- genus_flex %>%
    mutate(genus_flex_label = as.integer(prop_flexible >= T))

  joined <- nocturnal_only %>%
    filter(!is.na(Flexible_Cox)) %>%
    left_join(genus_labels, by = "Genus.1.2")

  tibble(
    threshold  = T,
    accuracy   = mean(joined$Flexible_Cox == joined$genus_flex_label),
    mixed_rate = 0
  )
}

flex_results <- map_dfr(flex_T_grid, evaluate_flex_at_T)
write_csv(flex_results, flex_out_csv)

# Summary

best_diel <- diel_results %>% arrange(desc(accuracy), mixed_rate, T1) %>% slice(1)
best_flex <- flex_results %>% arrange(desc(accuracy), threshold) %>% slice(1)

print(diel_results)
print(flex_results)
message(
  sprintf("Best diel: T1=%.2f | accuracy=%.3f | mixed=%.3f",
          best_diel$T1, best_diel$accuracy, best_diel$mixed_rate)
)
message(
  sprintf("Best flex: T=%.2f | accuracy=%.3f",
          best_flex$threshold, best_flex$accuracy)
)