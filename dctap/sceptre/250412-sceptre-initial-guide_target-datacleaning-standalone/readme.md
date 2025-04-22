# HPC Enhancer Screen via DC-TAP-seq

We are applying DC-TAP-seq to test over 30,000 E-G pairs across 2 cell states 
within a differentiating iPSC to hepatocyte cell model. In particular,
we designed guides against the iPSC and Definitive Endoderm (DE) state.

## HPC Screen

Hepatocytes comprises up to 80% of the total population and volume of the 
human liver. (Schulze) And, the liver represents one of the most pivotal 
organs of the human body in regulating glucose homeostasis, lipid metabolism, 
detoxification and many other physiological processes. (Li) Thus, we are
interested in using an established induced pluripotent stem cells (iPSC)
to hepatic progenitor cells (HPC) differentiation protocol as a biologically
relevant cellular model to screen which genes communicates with putative 
enhancer elements. 

### sgRNA Library Design

Initial HPC Screen guide design was created using 
[CRISPRDesigner](https://github.com/EngreitzLab/CRISPRDesigner)
pipeline. See below for details.

In order to design sgRNAs targeting candidate elements and gene promoters, 
we first created standard target regions. For candidate elements, we took 
the summit of DHS peaks as the central point and extended 150-bp in either 
direction, creating a 300-bp window in which to design sgRNAs. To target 
promoters, we created a 500-bp region spanning -250 to +250 relative to the 
Transcription Start Site (TSS). We then designed 15 independent sgRNAs per 
region according to our established pipeline (Fulco 2019, Nasser 2021, 
https://github.com/EngreitzLab/CRISPRDesigner). As negative controls, we 
included 400 non-targeting sgRNAs that lack matching sequence in the genome 
and 400 safe-targeting sgRNAs that align non-genic regions with no known 
open/active chromatin marks. For positive control gene TSS and enhancers, 
we included the same sgRNAs that we used in our previous study (Fulco 2019) 
that fall within the 500-bp TSS and the 300-bp enhancer region. In cases with 
fewer than 15 sgRNAs per region, we designed additional sgRNAs to supplement. 
We also added 5 additional TSS targeting positive control guides into 
the guide pool.

*Note:* The above description was copied from the DC-TAP-seq paper 
methodology section. See reference below.

#### 240905 HPC Screen sgRNA Design Files

This section notes where the original files associated with generating
the `guide_targets` file are located.

1. `d0_d2_guide_annotations_hg19.csv: {jray_gdrive}/Projects/CRISPRi_experiments/Assays/Diff_screen/Guide-design-large-scale/Attempt-1/Pipeline_outputs/To_order/Final`

2. `ChosenGenes.AllRegions.bed: {jray_gdrive}/Projects/CRISPRi_experiments/Assays/Diff_screen/Guide-design-large-scale/Attempt-1/Pipeline_outputs/`

2. `alt_TSS_ChosenGenes.AllRegions.bed: {jray_gdrive}/Projects/CRISPRi_experiments/Assays/Diff_screen/Guide-design-large-scale/Alt_TSS_design/ouputs/`

Note that the `*.bed` files are annotations and CRISPRDesigner pipeline are 
based on `hg19`.

## References

1. DC-TAP-seq Paper (submitting for publication pending).