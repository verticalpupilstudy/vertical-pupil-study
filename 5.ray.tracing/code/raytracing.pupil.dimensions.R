
# This script computes the vertical and horizontal pupil diameters/areas 
#for both round and vertical pupils at 19 incident angles
#saves results to csv/excel.

# Loads libraries
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)

# 1 reads chief rays document, computes pperpendicular unit vector 

chief <- read_excel("Chief_Rays.xlsx") %>%
  mutate(
    dx       = X1 - X0,
    dz       = Z1 - Z0,
    length   = sqrt(dx^2 + dz^2),
    ux_perp  = -dz / length,  # x  unit vector perpendicular to chief ray
    uz_perp  =  dx / length,  # z 
    p0x      = X0,
    p0z      = Z0
  ) %>%
  select(Incident_Angle, ux_perp, uz_perp, p0x, p0z)


# 2 Uses ray-tracing data to compute projections

main <- read_excel("Ray_Tracing_Main.xlsx")

joined <- main %>%
  left_join(chief, by = "Incident_Angle") %>%
  mutate(
    dx   = X1 - p0x,
    dz   = Z1 - p0z,
    proj = dx * ux_perp + dz * uz_perp
  )


# 3 Calculates width, height and area 

metrics <- joined %>%
  group_by(Pupil_Shape, Constriction_Level, Incident_Angle) %>%
  summarise(
    width  = max(proj, na.rm=TRUE) - min(proj, na.rm=TRUE),  # Full span across rays
    height = abs(Y1[Ray == "Superior"]) + abs(Y1[Ray == "Inferior"]),
    .groups = "drop"
  ) %>%
  mutate(
    area = pi * (width / 2) * (height / 2)
  )


# 4 Saves data to excel/csv document

write.xlsx(metrics, "Lens_Metrics_by_Angle.xlsx", rowNames = FALSE)
write.csv(metrics,  "Lens_Metrics_by_Angle.csv", row.names = FALSE)


# 5. Creates excel/csv document that compares vertical vs round pupils at constriction levels and angles

comparison <- metrics %>%
  pivot_wider(
    id_cols     = c(Constriction_Level, Incident_Angle),
    names_from  = Pupil_Shape,
    values_from = c(width, height, area),
    names_sep   = "_"
  ) %>%
  mutate(
    width_pct_diff  = (width_Vertical  - width_Round ) / width_Round  * 100,
    height_pct_diff = (height_Vertical - height_Round) / height_Round * 100,
    area_pct_diff   = (area_Vertical   - area_Round  ) / area_Round   * 100
  )

write.xlsx(comparison, "Vertical_vs_Round_Differences.xlsx", rowNames = FALSE)
write.csv(comparison,  "Vertical_vs_Round_Differences.csv", row.names = FALSE)


