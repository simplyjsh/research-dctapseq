# DC-TAP-seq qPCR Analysis Tool v0.0.1

> [!CAUTION]
> This tool will soon be depreciated.

This tool provides helper scripts and functions to
analyze qPCR data to orthogonally validate DC-TAP-seq
dataset effect sizes. The script adapted to the standardized qPCR
experiment design optimized by the folks at the Lander Lab and the
Broad. In addition, generic helper methods were added to format the
unique biological conditions of hiPSC differentiated samples.

## Installation

### Clone repository

```bash
git clone https://github.com/simplyjsh/research-dctapseq.git
```

### Setup environment

> [!TIP]
> Setting up with python virtual environment is highly recommended.

```bash
cd dctaptoolkit/qpcr/qpcr-analysis-tool-v001/
python -m venv .venv
source .venv/bin/activate
pip install .
```

## Usage

For performing analyses the recommendation is to create a notebook
(e.g. ipynb or marimo) and importing the necessary functions from
the **qpcr-analysis-tool-v001** package.

Below is an example generic notebook of the steps used to parse the raw data,
set the conditions of your experiment, and perform the qPCR analysis.

Note this qPCR tool was specifically adapted for our experimental
setup. Thus the parsing of the raw data may not apply to your
raw data input. However, the `get_deltaCq_*()` might still be
useful to you.

**Please use at your own discretion. An updated
and more flexible version is currently being developed**.

```python
# Cell 1
##############################################################################
# Imports
from dctap.qpcr.libs.pandas import set_defaultoptions, displaydf_full
from dctap.qpcr.constants import QPCRPATHS
from dctap.qpcr.utils.core import *
##############################################################################

# Cell 2
##############################################################################
# Set pandas settings
set_defaultoptions(pd, supresscopywarning=None)
output_notebook()
##############################################################################

# Cell 3
##############################################################################
# Read and annotate data from different plates and combine them
experiment_id = "experimentID-250306"
plate_ids = ["experimentID-250306-plate1", "experimentID-250306-plate2"]
dfs = []
for plate_id in plate_ids:
    dfs.append(get_plate_data(experiment_id, plate_id))
df = pd.concat(dfs)
df.reset_index(inplace=True, drop=True)
display(df)
##############################################################################

# Cell 4
##############################################################################
# Get sample metadata
df_samples = get_sample_metadata(cast(pd.Series, df.Sample), sep="_")
with displaydf_full():
  display(df_samples)
##############################################################################

# Cell 5
##############################################################################
# Set condition labels based off your experiment
conditions = ["bio_reps", "ctrl_calibrator", "cond_diff"]
df = set_conditions(
    df,
    df_samples,
    conditions=conditions,
    merge_cols=["0123", "0123", "012"],
)
display(df)
##############################################################################

# Cell 6
##############################################################################
# Get deltaCq for all samples
df1 = get_deltaCq_expression_bulkdata(
    df,
    ref_primer="GAPDH",
    test_primers=get_primers(df),
    drop_customcols=conditions,
)
display(df1)
##############################################################################

# Cell 7
##############################################################################
# Assign biological controls to calculate ddCq
df2 = get_deltaCq_stats(df1, biorep_col="bio_reps")
df_calibrators = get_calibrators(
    df2,
    ctrl_col="ctrl_calibrator",
    condition_col="cond_diff",
    assign_ctrl_samples=[
        "Unperturbed_ctrl",
        "Perturbed_ctrl_cond1",
        "Perturbed_ctrl_cond2",
        "Perturbed_ctrl_cond3",
        "Perturbed_ctrl_cond4",
    ],
    assign_cond_group=[
        "Unperturbed_NaN_ctrl",
        "Perturbed_sample_cond1",
        "Perturbed_sample_cond2",
        "Perturbed_sample_cond3",
        "Perturbed_sample_cond4",
    ],
)
with displaydf_full():
  display(df_calibrators)
##############################################################################

# Cell 8
##############################################################################
# Calculate ddCq fold change
df3 = get_deltadeltaCqMethod_foldchange(
    df1, df_calibrators, biorep_col="bio_reps", condition_col="cond_diff"
)
with displaydf_full():
  display(df3)
##############################################################################

# Cell 8 and beyond
##############################################################################
# You should now be able to use df3 for plotting purposes or save as tsv.
# Plot or write to disk as needed...
##############################################################################
```
