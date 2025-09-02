library(here)
library(glue)
library(dplyr)

# Set useful env variables
data_dir <- here("pipeline_inputs/")
inputs_dir <- glue(data_dir, "data_cleaning_inputs/")

# Load filepath of inputs
coords_fp <- glue(inputs_dir, "240905_hpc_screen_guide_coordinates_hg19.tsv")
acoords_fp <- glue(inputs_dir, "240905_hpc_screen_alt_guide_coordinates_hg19.tsv")
annotations_fp <- glue(inputs_dir, "240905_hpc_screen_guide_annotations_hg19.tsv")
dctap_targets_fp <- glue(inputs_dir, "250401_dctap_random_screen_wtc11_guide_targets.tsv")

# Read as dataframe without header
g_coords <- read.table(coords_fp, sep = "\t", header = FALSE)
g_altcoords <- read.table(acoords_fp, sep = "\t", header = FALSE)
g_annotations <- read.table(annotations_fp, sep = "\t", header = FALSE)
g_dctap_targets <- read.table(dctap_targets_fp, sep = "\t", header = TRUE)

# Append alternative guide target coordinates to rest of guide target coords
g_allcoords <- rbind(g_coords, g_altcoords)

# Append safe_targeting & negative_controls guide annotations
g_ctrls <- data.frame(
  V1 = c("NA", "NA"),
  V2 = c("NA", "NA"),
  V3 = c("NA", "NA"),
  V4 = c("safe_targeting", "negative_control"),
  stringsAsFactors = FALSE
)
g_allcoords <- rbind(g_allcoords, g_ctrls)

# Reassign headers
colnames(g_allcoords) <- c(
  "target_chr", 
  "target_start", 
  "target_end", 
  "target_name"
)
colnames(g_annotations) <- c("name", "target_name")

#' Merge using many-to-many join via dplyr inner_join function
#' Note the select function reorganizes target_name column to be at end
#' 
#' Next, assign target_type based on target_name
#' We assume that all target_name data is an enhancer target,
#' all safe_targeting is safe_targeting, all negative_control is
#' negative_control, and those ending with "*_TSS" are tss_control
guide_targets <- inner_join(g_annotations, g_allcoords, by = "target_name") %>%
  select(-target_name, everything(), target_name) %>%
  mutate(
    target_type = case_when(
      startsWith(target_name, "chr")    ~ "enh",
      target_name == "safe_targeting"   ~ "safe_targeting",
      target_name == "negative_control" ~ "negative_control",
      TRUE                              ~ "tss_control"
    )
  )

# Append WTC11 Random Screen guide target coordinates with rest of guide targets
g_wtc11_randscreen_annotations <- g_annotations %>% 
  as_tibble() %>% 
  filter(grepl("^WTC11", name))

g_prevscreen_ctrls <- data.frame(
  name = g_wtc11_randscreen_annotations$name,
  target_chr = "NA",
  target_start = 0,
  target_end = 0,
  target_name = g_wtc11_randscreen_annotations$target_name,
  target_type = "NA",
  stringsAsFactors = FALSE
)

g_dctap_targets <- g_dctap_targets %>%
  select(name, target_chr, target_start, target_end, target_name, target_type)
g_prevscreen_ctrls <- g_prevscreen_ctrls %>%
  rows_update(g_dctap_targets, by = "name", unmatched = "ignore")
guide_targets <- rbind(guide_targets, g_prevscreen_ctrls)

# Write to file
guide_targets_output_fp <- inputs_dir <- glue(
  data_dir, 
  "data_cleaning_outputs/guide_targets.tsv"
)
write.table(
  guide_targets, 
  file = guide_targets_output_fp,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)