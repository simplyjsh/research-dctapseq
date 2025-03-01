import re
import glob
import numpy as np
import pandas as pd
from pathlib import Path
from typing import cast

from dctap.qpcr.constants import QPCRPATHS
from dctap.qpcr.constants import PLATES
from dctap.qpcr.constants.plates import Plate96
from dctap.qpcr.libs import Incrementor


# -----------------------------------------------------------------------
# INFO: Helper methods
def _layout_to_annotation(
    experiment_id: str,
    plate_id: str,
    data_dir: Path = QPCRPATHS.DATADIR,
    code_offset: int = 1,
    code_file: Path | str | None = None,
    template_layout_file: Path | str | None = None,
    primer_layout_file: Path | str | None = None,
):
    """
    Given a plate id, code, template layouts, and primer layouts,
    create an annotation csv file.
    """

    if code_file is None:
        code_pattern: str = str(
            data_dir / experiment_id / (experiment_id + "*code.csv")
        )
        code_files = glob.glob(code_pattern)
        assert len(code_files) == 1
        code_file = code_files[0]

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

    # Handles reading Code Layout
    with open(code_file, "r") as f:
        lines = f.readlines()
    assert lines[0].startswith("code")

    code_map = []
    for line in lines[1:]:
        temp = line.strip().split(",")
        assert len(temp) == 2
        code_map = code_map + temp[1:]

    # Handles reading Template Layout
    with open(template_layout_file, "r") as f:
        lines = f.readlines()
    assert lines[0].startswith("template")

    template_map = []
    for line in lines[1:]:
        temp = line.strip().split(",")
        assert len(temp) == 13
        template_map = template_map + temp[1:]

    assert len(template_map) == PLATES.P96.TOTALWELLS

    # Handles reading Primer Layout
    with open(primer_layout_file, "r") as f:
        lines = f.readlines()
    assert lines[0].startswith("primer")

    primer_map = []
    for line in lines[1:]:
        temp = line.strip().split(",")
        assert len(temp) == 13
        primer_map = primer_map + temp[1:]

    assert len(primer_map) == PLATES.P96.TOTALWELLS

    # Converts to 384-well plate layouts
    template_list = [None] * PLATES.P384.TOTALWELLS
    primer_list = [None] * PLATES.P384.TOTALWELLS
    increment = Incrementor()
    increment_primer = Incrementor()

    for dozen in range(0, len(template_map), PLATES.P96.COL):
        row = template_map[dozen : dozen + PLATES.P96.COL]
        for well in row:
            sample = ""
            primer = ""

            try:
                sample_id: int = int(well) - code_offset
                sample = code_map[sample_id]
                primer = primer_map[increment_primer]
            except ValueError:
                sample = ""
                primer = ""

            template_list[increment] = sample  # type: ignore
            template_list[increment + PLATES.P384.COL] = sample  # type: ignore

            primer_list[increment] = primer  # type: ignore
            primer_list[increment + PLATES.P384.COL] = primer  # type: ignore

            increment.step()

            template_list[increment] = sample  # type: ignore
            template_list[increment + PLATES.P384.COL] = sample  # type: ignore

            primer_list[increment] = primer  # type: ignore
            primer_list[increment + PLATES.P384.COL] = primer  # type: ignore

            increment.step()
            increment_primer.step()

        increment.skip(PLATES.P384.COL)

    annotation_file = data_dir / experiment_id / (plate_id + "-annotation.csv")
    with open(annotation_file, "w") as f:
        f.write(",".join(["Well", "Sample", "Primer"]) + "\n")
        for i in range(PLATES.P384.TOTALWELLS):
            f.write(
                ",".join([PLATES.P384.WELLS[i], template_list[i], primer_list[i]])  # type: ignore
                + "\n"
            )

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
    annotation_file = str(annotation_files[0])

    data_pattern = str(data_dir / experiment_id / ("*" + plate_id + ".csv"))
    data_files = glob.glob(data_pattern)
    assert len(data_files) == 1
    data_file = data_files[0]

    df = cast(pd.DataFrame, pd.read_csv(annotation_file))
    df_annotations = cast(pd.DataFrame, df[["Well", "Sample", "Primer"]])

    df = cast(pd.DataFrame, pd.read_csv(data_file))
    df_data = cast(pd.DataFrame, df[["Well", "Cq"]])

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

    df1 = get_deltaCq_expression_data(df, test_primer="CER1", ref_primer="GAPDH")
