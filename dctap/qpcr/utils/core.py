import re
import glob
import numpy as np
import pandas as pd
from pathlib import Path

from dctap.qpcr.constants import QPCRPATHS
from dctap.qpcr.constants import PLATES


# -----------------------------------------------------------------------
# INFO: Helper methods
def _layout_to_annotation(
    experiment_id: str,
    plate_id: str,
    data_dir: Path = QPCRPATHS.DATADIR,
    template_layout_file: Path | str | None = None,
    primer_layout_file: Path | str | None = None,
):
    """
    Given a plate id, template layouts, and primer layouts,
    create an annotation csv file.
    """

    if template_layout_file is None:
        template_pattern: str = str(
            data_dir / experiment_id / ("*" + plate_id + "*template*layout.csv")
        )
        template_layout_files = glob.glob(template_pattern)
        assert len(template_layout_files) == 1
        template_layout_file = template_layout_files[0]

    if primer_layout_file is None:
        primer_pattern: str = str(
            data_dir / experiment_id / ("*" + plate_id + "*primer*layout.csv")
        )
        primer_layout_files = glob.glob(primer_pattern)
        assert len(primer_layout_files) == 1
        primer_layout_file = primer_layout_files[0]

    # Handles Template Layout
    with open(template_layout_file, "r") as f:
        lines = f.readlines()
    assert lines[0].startswith("Well")

    wells = []
    template_list = []
    for line in lines[1:]:
        temp = line.strip().split(",")
        assert len(temp) == 2
        wells = wells + temp[:1]
        template_list = template_list + temp[1:]

    assert len(wells) == PLATES.P384.TOTALWELLS
    assert len(template_list) == PLATES.P384.TOTALWELLS

    # Handles Primer Layout
    with open(primer_layout_file, "r") as f:
        lines = f.readlines()
    assert lines[0].startswith("Well")

    primer_list = []
    for line in lines[1:]:
        temp = line.strip().split(",")
        assert len(temp) == 2
        primer_list = primer_list + temp[1:]

    assert len(primer_list) == PLATES.P384.TOTALWELLS

    annotation_file = data_dir / experiment_id / (plate_id + "-annotation.csv")
    with open(annotation_file, "w") as f:
        f.write(",".join(["Well", "Sample", "Primer"]) + "\n")
        for i in range(PLATES.P384.TOTALWELLS):
            f.write(",".join([wells[i], template_list[i], primer_list[i]]) + "\n")

    return [str(annotation_file)]


def extract_key(s):
    m = re.match(r"^(.*)_TR\d+_(.*)$", s)
    if m:
        # Combine the two parts to form the key
        return m.group(1) + "_" + m.group(2)
    else:
        return s


# -----------------------------------------------------------------------
# INFO: Core Functions
def get_plate_data(
    experiment_id: str, plate_id: str, data_dir: Path = QPCRPATHS.DATADIR
):
    """
    Given a plate id, locate the annotations csv file and raw data csv files,
    then extract useful data, merge, and return the dataframe.
    """
    annotation_pattern = str(
        data_dir / experiment_id / ("*" + plate_id + "*annotation.csv")
    )
    annotation_files = glob.glob(annotation_pattern)

    if len(annotation_files) == 0:
        annotation_files = _layout_to_annotation(experiment_id, plate_id)
    assert len(annotation_files) == 1
    annotation_file = annotation_files[0]

    data_pattern = str(data_dir / experiment_id / ("*" + plate_id + ".csv"))
    data_files = glob.glob(data_pattern)
    assert len(data_files) == 1
    data_file = data_files[0]

    df = pd.read_csv(annotation_file)
    df_annotations = df[["Well", "Sample", "Primer"]]

    df = pd.read_csv(data_file)
    df_data = df[["Well", "Cq"]]

    df = df_annotations.merge(df_data, on="Well")
    df.dropna(inplace=True)
    df["plate_id"] = [plate_id] * len(df)
    df["experiment_id"] = [experiment_id] * len(df)

    return df


def get_deltaCq_expression_data(
    df,
    test_primer,
    ref_primer,
):
    """
    From the annotated dataframe containing raw qPCR data,
    calculate average Cq values of replicates and deltaCq values
    of each sample against the reference primer.
    """
    for i in ["Well", "Sample", "Primer", "Cq", "experiment_id"]:
        assert i in df.columns

    # Get test and reference primers
    df_ref = df[df.Primer == ref_primer].copy()
    df_ref.rename(columns={"Cq": "Cq_ref"}, inplace=True)

    df_test = df[df.Primer == test_primer].copy()
    df_test.rename(columns={"Cq": "Cq_test"}, inplace=True)

    # Get average Cq values of technical replicates
    df_ref_mean = df_ref.groupby("Sample", as_index=False)[["Cq_ref"]].mean()
    df_test_mean = df_test.groupby("Sample", as_index=False)[["Cq_test"]].mean()
    df = df_ref_mean.merge(df_test_mean, on="Sample")

    # Calculate deltaCq
    df["deltaCq"] = df.Cq_test - df.Cq_ref

    return df


if __name__ == "__main__":
    experiment_id = "JR95-250129"
    plate_ids = ["JR95-250129-plate1", "JR95-250129-plate2", "JR95-250129-plate3"]

    dfs = []
    for plate_id in plate_ids:
        dfs.append(get_plate_data(experiment_id, plate_id))

    df = pd.concat(dfs)
    df.reset_index(inplace=True, drop=True)

    # Adding a few helpful columns
    df["group"] = [df.Sample[i] + "___" + df.Primer[i] for i in range(len(df))]
    df["well_id"] = [df.plate_id[i] + "_" + df.Well[i] for i in range(len(df))]
    df["relExp_25"] = [2 ** (25 - df.Cq[i]) for i in range(len(df))]

    df1 = get_deltaCq_expression_data(df, test_primer="OCT4", ref_primer="GAPDH")
    df2 = get_deltaCq_expression_data(df, test_primer="CER1", ref_primer="GAPDH")

    print(df1)
