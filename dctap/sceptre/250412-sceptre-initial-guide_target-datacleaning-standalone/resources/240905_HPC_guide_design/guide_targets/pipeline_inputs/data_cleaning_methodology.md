# Guide Targets Data Cleaning Methodology
This file details the generation of the `guide_targets.tsv` file required
for further downstream analysis processing.

## Overview
Each sgRNA library design pipeline have unique circumstances that requires 
minor adjustments. For example, off-target guide designs identified via
visual interactive plots such as [IGV]() requires designing
alternative guides for correction, scattering the information
required in `guide_targets.tsv` among many files instead of one.

The `guide_targets.tsv` file needed is an input to generate a
`sceptre_object` via the `grna_targets_data_frame` and `discovery_pair`
arguments. The details of how the `guide_targets.tsv` is used
in the pipeline can be referenced by the 
[DC_TAP_Paper](https://github.com/jamesgalante/DC_TAP_Paper)
repo. See `workflow/rules/process_validation_datasets/sceptre_setup.smk` 
in the github repo.

The *minimum* information needed in the `guide_targets.tsv` file is detailed
in the table below.

| name         | target_chr   | target_start | target_end   | target_name  | target_type  |
| ------------ | ------------ | ------------ | ------------ | ------------ | ------------ |
| YYMMDD_H...1 | chrX         | YYYYYYYYYYYY | ZZZZZZZZZZZZ | chrX:YYY-ZZZ | enh          |
| YYMMDD_H...2 | chrX         | YYYYYYYYYYYY | ZZZZZZZZZZZZ | AAA_TSS      | tss_pos_ctrl |
| YYMMDD_H...3 | chrX         | YYYYYYYYYYYY | ZZZZZZZZZZZZ | AAA_alt_TSS  | tss_pos_ctrl |
| YYMMDD_H...4 | NA           | NA           | NA           | safe_targ... | safe_targ... |
| YYMMDD_H...5 | NA           | NA           | NA           | negative_... | negative_... |
| ...          | ...          | ...          | ...          | ...          | ...          |

### Inputs
This section details the expected input files for data cleaning rules regarding 
`240905 HPC Screen Design` and the expected information within each file. 

1. `240905_hpc_screen_guide_annotations_hg19.tsv`.
   Note, file contains `name`(OligoID) and `target_name`(guideSet).
   File does not contain TSS target coordinates.

2. `240905_hpc_screen_guide_coordinates_hg19.tsv`.
   Note, file contains `target_chr`, `target_start`, `target_end`, and 
   `target_name` which also includes coordinates for tss targets. 
   File does not contain _alternative_ TSS target coordinates.
   File also does not contain _safe-targeting_ and _negative-control_
   guide annotations.

3. `240905_hpc_screen_alt_guide_coordinates_hg19.tsv`.
   Note, file contains `target_chr`, `target_start`, `target_end`, and 
   `target_name` which also includes coordinates for alternative tss targets.

### Outputs
This section details the output files for further downstream analysis processing
based on the pipeline detailed in this repo: 
[DC_TAP_Paper](https://github.com/jamesgalante/DC_TAP_Paper).

1. `guide_targets.tsv`

## Methodology
This section describes the simple methodology used to merge the file
listed in the sections above to generate the `guide_targets.tsv` file.

1. First, load `240905_hpc_screen_guide_coordinates_hg19.tsv` and read as a 
   tsv without headers. Next, similarly load 
   `240905_hpc_screen_alt_guide_coordinates_hg19.tsv`, and then
   append to the end of `240905_hpc_screen_guide_coordinates.tsv`.

2. Create a two-row dataframe to append `safe_targeting` and 
   `negative_control` guide annotations to the end of
   `240905_hpc_screen_guide_coordinates.tsv`. Refer back to 
   [Overview](#overview) for the expected data points associated
   with the safe targeting and negative control annotations.

3. Second, load `240905_hpc_screen_guide_annotations_hg19.tsv` and read as a
   tsv without headers.

4. Reassign headers of `240905_hpc_screen_guide_coordinates.tsv` and
   `240905_hpc_screen_guide_annotations_hg19.tsv` according to 
   `guide_targets.tsv` header specifications. Refer back to 
   [Overview](#overview) for expected header names.

5. Third, merge `240905_hpc_screen_guide_coordinates.tsv` and 
   `240905_hpc_screen_guide_annotations_hg19.tsv` on `target_name`.
   The recommended merge method is with a _many-to-many_ join function.
   For example, the `dplyr` library has the *inner_join()* function
   that accomplishes the task above without lost of duplicate 
   `target_name` data rows.

6. Fourth, assign `target_type` header based on the data row in the
   `target_name` on a case-by-case bases described below. We assume
   that We assume that all target names with a genome coordinate 
   is an enhancer target and the rest are those ending with "*_TSS" 
   which are tss_control targets. 

```r
target_type <- case_when(
   startsWith(target_name, "chr")    ~ "enh",
   target_name == "safe_targeting"   ~ "safe_targeting",
   target_name == "negative_control" ~ "negative_control",
   TRUE                              ~ "tss_control"
)
```

7. Finally, export the final dataframe as `guide_targets.tsv`.