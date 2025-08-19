# Script that adds average genus‐level mass to Mammal_PupilTraits_GenusLevel.csv

pkgs <- c("readr", "dplyr", "openxlsx")
for(p in pkgs) if(!requireNamespace(p, quietly=TRUE)) install.packages(p)

library(readr)
library(dplyr)
library(openxlsx)

species <- read_csv("Mammal_PupilTraits_SpeciesLevel.csv", show_col_types = FALSE) %>%
  filter(Excluded == 0)

valid_pupils    <- c("Vertical", "Round", "Horizontal")
eligible_genera <- species %>%
  filter(Pupil %in% valid_pupils) %>%
  pull(Genus.1.2) %>%
  unique()

genus_mass <- species %>%
  filter(Genus.1.2 %in% eligible_genera) %>%
  group_by(Genus.1.2) %>%
  summarise(
    avg_mass_g = mean(Mass.g, na.rm = TRUE),
    .groups = "drop"
  )

genus_level <- read_csv("Mammal_PupilTraits_GenusLevel.csv", show_col_types = FALSE)

genus_level_with_mass <- genus_level %>%
  left_join(genus_mass, by = "Genus.1.2")

# Save to a excel file
write.xlsx(
  genus_level_with_mass,
  "Mammal_PupilTraits_GenusLevel_withMass.xlsx",
  rowNames = FALSE
)
