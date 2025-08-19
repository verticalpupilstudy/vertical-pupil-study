
#this script averages ancestral probs accross posterior distribution
# Load required packages
library(readxl)
library(dplyr)

# Set working directory to folder containing transition count Excel files
# Replace with the relative path from the project root if needed
setwd("3.asr/output/SCM_pupil_UphamPosterior/posterior_probs/xlsx")

# List only real Excel files, exclude temp/lock files like those starting with '~$'
files <- list.files(pattern = "^[^~].*\\.xlsx$")

# Read and bind all transition data
all_data <- lapply(files, read_excel) %>%
  bind_rows()

# Summarise by transition type
summary_stats <- all_data %>%
  group_by(From, To) %>%
  summarise(
    mean_count = mean(Count),
    sd_count = sd(Count),
    n = n(),
    ci_95_lower = mean_count - 1.96 * sd_count / sqrt(n),
    ci_95_upper = mean_count + 1.96 * sd_count / sqrt(n),
    min_count = min(Count),
    max_count = max(Count),
    .groups = "drop"
  )

# Print summary
print(summary_stats)

# write to file
write.csv(summary_stats, "Transition_Summary_Stats.csv", row.names = FALSE)
