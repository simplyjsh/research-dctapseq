import re
import glob
import numpy as np
import pandas as pd
from pathlib import Path
from typing import cast

from dctap.qpcr.constants import QPCRPATHS
from dctap.qpcr.constants import PLATES
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

    # Asserts that required layout files exists
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

            # Sometimes wells can be empty or Cq values is N/A
            try:
                sample_id: int = int(well) - code_offset
                sample = code_map[sample_id]
                primer = primer_map[increment_primer]
            except ValueError:
                sample = ""
                primer = ""

            """
            Munson G. 96-well plate to 384-well plate layout
            Each well on 96-well represents 4 neigbhoring wells on 384-well plate
            E.g. 96-well plate A01 := the same wells A01, A02, B01, B02 on 384-well plate
                 and so on.
            """
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

    # Write annotation file
    annotation_file = data_dir / experiment_id / (plate_id + "-annotation.csv")
    with open(annotation_file, "w") as f:
        f.write(",".join(["Well", "Sample", "Primer"]) + "\n")
        for i in range(PLATES.P384.TOTALWELLS):
            f.write(
                ",".join([PLATES.P384.WELLS[i], template_list[i], primer_list[i]])  # type: ignore
                + "\n"
            )

    # return path to annotation file as list[str] consistent with glob query
    return [str(annotation_file)]


def get_primers(df):
    primers = []
    if "Primer" in df.columns:
        primers = df.Primer.unique().tolist()
    return primers


# WARNING:
# Depreciating soon: no longer setting conditions after getting relative Cq values
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

    # Locate annotation csv file. If it does not exist, create one
    annotation_pattern = str(
        data_dir / experiment_id / ("*" + plate_id + "*annotation.csv")
    )
    annotation_files = glob.glob(annotation_pattern)

    if len(annotation_files) == 0:
        annotation_files = _layout_to_annotation(experiment_id, plate_id)
    assert len(annotation_files) == 1
    annotation_file = str(annotation_files[0])

    # Locate and read data
    data_pattern = str(data_dir / experiment_id / ("*" + plate_id + ".csv"))
    data_files = glob.glob(data_pattern)
    assert len(data_files) == 1
    data_file = data_files[0]

    # Merging data and annotations
    df = cast(pd.DataFrame, pd.read_csv(annotation_file))
    df_annotations = cast(pd.DataFrame, df[["Well", "Sample", "Primer"]])

    df = cast(pd.DataFrame, pd.read_csv(data_file))
    df_data = cast(pd.DataFrame, df[["Well", "Cq"]])

    df = df_annotations.merge(df_data, on="Well")
    df.dropna(inplace=True)

    first_col_index = 0
    df.insert(first_col_index, "experiment_id", ([experiment_id] * len(df)))  # type: ignore
    df["plate_id"] = [plate_id] * len(df)

    return df


def get_sample_metadata(sample: pd.Series | list[str], sep: str = "_") -> pd.DataFrame:
    """
    Given a list of strings and a separator, return a dataframe the string
    split by separators into columns. The purpose for this method is to
    help assign conditions.
    """
    # Convert to Series if provided a sample list
    if sample is not pd.Series:
        sample = pd.Series(sample)

    # Get unique sample names and split metadata for modular labelling
    sample = sample.copy().drop_duplicates().reset_index(drop=True)
    metadata = sample.str.split(sep, expand=True)

    first_col_index = 0
    metadata.insert(
        first_col_index,
        "Sample",
        sample,
    )

    return metadata


def set_conditions(
    df: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    conditions: list[str],
    merge_cols: list[str],
    additions: list[str] | None = None,
    append: list[str] | None = None,
) -> pd.DataFrame:
    """
    Creates condition columns for the df DataFrame based on the given conditions
    list of list strings using the sample_metadata as a mapping.
    """
    sample_metadata = sample_metadata.copy()
    metadata_colnames = sample_metadata.columns[1:].tolist()

    # Create user defined condition columns
    for i, condition in enumerate(conditions):
        cols = list(merge_cols[i])
        sample_metadata[condition] = sample_metadata[int(cols[0])].str.cat(
            [sample_metadata[int(col)] for col in cols[1:]], sep="_"
        )

    if additions is not None and append is not None:
        # Ensure the addition columns exists
        for addition in additions:
            assert addition in sample_metadata.columns

        # Append additional comments to defined condition columns
        for i, addition in enumerate(additions):
            sample_metadata[addition] = sample_metadata[addition] + "_" + append[i]

    # Drop sample_metadata columns
    sample_metadata = sample_metadata.drop(metadata_colnames, axis="columns")

    # Merge condition columns to main df
    df = df.copy()
    df = df.merge(sample_metadata, on="Sample", how="left")

    return df


def get_deltaCq_expression_data(
    df,
    ref_primer: str,
    test_primer: str,
):
    """
    From the annotated dataframe containing raw qPCR data,
    calculate average Cq values of replicates and deltaCq values
    of each sample against the reference primer.
    """
    for i in ["Well", "Sample", "Primer", "Cq"]:
        assert i in df.columns
    df = df.copy()

    ref_col = f"Cq_ref_{ref_primer}"
    tes_col = f"Cq_test_{test_primer}"

    # Get test and reference primers
    df_ref = df[df.Primer == ref_primer]
    df_ref.rename(columns={"Cq": ref_col}, inplace=True)

    df_test = df[df.Primer == test_primer]
    df_test.rename(columns={"Cq": tes_col}, inplace=True)

    # Get average Cq values of technical replicates
    ref_dict = {
        ref_col + "_mean": (ref_col, "mean"),
        ref_col + "_std": (ref_col, "std"),
        ref_col + "_count": (ref_col, "count"),
    }
    tes_dict = {
        tes_col + "_mean": (tes_col, "mean"),
        tes_col + "_std": (tes_col, "std"),
        tes_col + "_count": (tes_col, "count"),
    }

    df_ref_stats = df_ref.groupby("Sample", as_index=False).agg(**ref_dict)
    df_test_stats = df_test.groupby("Sample", as_index=False).agg(**tes_dict)

    df_ref_stats[ref_col + "_ste"] = df_ref_stats[ref_col + "_std"] / np.sqrt(
        df_ref_stats[ref_col + "_count"]
    )
    df_test_stats[tes_col + "_ste"] = df_test_stats[tes_col + "_std"] / np.sqrt(
        df_test_stats[tes_col + "_count"]
    )

    # Merge and drop metadata results
    df_results = df_ref_stats.merge(df_test_stats, on="Sample")
    df_results["deltaCq"] = (
        df_results[tes_col + "_mean"] - df_results[ref_col + "_mean"]
    )
    df_results.rename(
        columns={"deltaCq": f"deltaCq_{test_primer}v{ref_primer}"}, inplace=True
    )
    df_results = df_results.drop(
        [ref_col + "_count", tes_col + "_count"], axis="columns"
    )

    # Drop rows for primers that are not test primer
    df = df[df.apply(lambda row: row["Primer"] == test_primer, axis="columns")]
    df = df.drop(
        ["experiment_id", "Well", "Primer", "Cq", "plate_id"],
        axis="columns",
    )

    # Collapse df
    df = df.drop_duplicates(subset="Sample")
    df = df.merge(df_results, on="Sample", how="right")

    return df


# NOTE:
# This method probably needs to be refactored as there are many
# unnecessary computations occuring with running it through
# the get_deltaCq_expression_data() method
def get_deltaCq_expression_bulkdata(
    df, ref_primer: str, test_primers: list[str], drop_customcols: list[str] = []
):
    """
    From the annotated dataframe containing raw qPCR data,
    calculate average Cq values of replicates and deltaCq values
    of each sample against the reference primer.
    """
    # Remove ref_primer if include in test_primer list
    if ref_primer in test_primers:
        test_primers.remove(ref_primer)

    df = df.copy()
    df_results = pd.DataFrame(df["Sample"].drop_duplicates().reset_index(drop=True))
    for i, test_primer in enumerate(test_primers):
        df1 = get_deltaCq_expression_data(df, ref_primer, test_primer)

        # Drop unnecessary duplicate metadata columns generated from
        # get_deltaCq_expression_data()
        if i > 0:
            df1 = df1.drop(
                (
                    [
                        f"Cq_ref_{ref_primer}_mean",
                        f"Cq_ref_{ref_primer}_std",
                        f"Cq_ref_{ref_primer}_ste",
                    ]
                    + drop_customcols
                ),
                axis="columns",
            )
        df_results = df_results.merge(df1, on="Sample", how="left")

    return df_results


def get_deltaCq_stats(df, biorep_col: str):
    """
    Calculate average deltaCq values based on user provided biorep assignment
    """
    # Ensure deltaCq and user defined ctrl columns exist
    for i in ["Sample", biorep_col]:
        assert i in df.columns, (
            f"The assigned biological replicate column ({biorep_col})",
            "does not exist in the DataFrame",
        )

    deltacq_pattern = r"deltaCq_"
    assert any(re.match(deltacq_pattern, col) for col in df.columns), (
        "At least one column starting with 'deltaCq_' must exist"
    )

    # Drop all uncessary rows for results df
    df1 = df.copy()
    cq_pattern = r"Cq_"
    cq_cols = (df1.filter(regex=cq_pattern).columns).append(pd.Index(["Sample"]))
    df1 = (
        df1.drop(cq_cols, axis="columns")
        .drop_duplicates(subset=biorep_col, keep="first")
        .reset_index(drop=True)
    )
    df_results = pd.DataFrame(df1)

    # get all deltaCq columns
    deltacq_pattern = r"deltaCq_"
    deltacq_cols = df.filter(regex=deltacq_pattern).columns

    # Caclulate deltaCq stats
    for col in deltacq_cols:
        stats_dict = {
            col + "_mean": (col, "mean"),
            col + "_std": (col, "std"),
            col + "_count": (col, "count"),
        }
        df_stats = df.groupby(biorep_col, as_index=False).agg(**stats_dict)
        df_stats[col + "_ste"] = df_stats[col + "_std"] / np.sqrt(
            df_stats[col + "_count"]
        )

        # Drop metadata cols and merge to results
        df_stats = df_stats.drop([col + "_count"], axis="columns")
        df_results = df_results.merge(df_stats, on=biorep_col, how="left")

    return df_results


def get_calibrators(
    df,
    ctrl_col: str,
    condition_col: str,
    assign_ctrl_samples: list[str],
    assign_cond_group: list[str],
):
    """
    Get calibrator based on the user assigned ctrl and condition columns
    for deltadeltaCq calculations.
    """
    # Ensure deltaCq and user defined ctrl columns exist
    for i in [ctrl_col, condition_col]:
        assert i in df.columns, (
            f"The assigned control ({ctrl_col}) or condition ({condition_col})",
            "column does not exist in the DataFrame",
        )

    deltacq_pattern = r"deltaCq_.*_mean"
    assert any(re.match(deltacq_pattern, col) for col in df.columns), (
        "At least one column starting with 'deltaCq_*_mean' must exist"
    )

    # Ensure the assigned ctrl sames and condition groups exists
    assert set((assign_ctrl_samples)).issubset(set(df[ctrl_col].unique())), (
        f"Missing expected items: {set(assign_ctrl_samples) - set(df[ctrl_col].unique())}"
    )
    assert set((assign_cond_group)).issubset(set(df[condition_col].unique())), (
        f"Missing expected items: {set(assign_cond_group) - set(df[condition_col].unique())}"
    )

    df = df.copy()
    # Get calibrators based on assigned ctrls
    df_rows = df[df[ctrl_col].isin(assign_ctrl_samples)]
    deltacq_cols = [col for col in df.columns if re.match(deltacq_pattern, col)]

    calibrators = {}
    for _, row in df_rows.iterrows():
        ctrl_sample = row[condition_col]
        calibrators[ctrl_sample] = {
            re.sub(r"_mean$", "", calibrator): row[calibrator]
            for calibrator in deltacq_cols
        }

    return pd.DataFrame.from_dict(calibrators, orient="index")


def get_deltadeltaCqMethod_foldchange(
    df, df_calibrators, biorep_col: str, condition_col: str, keep_metadata: bool = False
):
    # TODO:
    # Notes about this method
    # Ensure that required cols and rows exists

    df = df.copy()
    # drop unnecessary Cq metadata/stats
    cq_pattern = r"^Cq_"
    cq_cols = (df.filter(regex=cq_pattern).columns).append(pd.Index(["Sample"]))
    df = df.drop(cq_cols, axis="columns").reset_index(drop=True)

    # Assign calibrators to the proper samples
    deltacq_pattern = r"deltaCq_.*"
    deltacq_cols = [col for col in df.columns if re.match(deltacq_pattern, col)]

    for col in deltacq_cols:
        calibrator_col = f"{col}_calibrator"
        deltadeltacq_col = f"delta{col}"
        foldchange = f"2^({deltadeltacq_col})"

        # Assignment
        df[calibrator_col] = df[condition_col].apply(
            lambda row: df_calibrators.loc[row, col]
        )

        # Calculations
        df[deltadeltacq_col] = df[col] - df[calibrator_col]
        df[foldchange] = np.power(2, -df[deltadeltacq_col])

        # Stats
        stats_dict = {
            foldchange + "_mean": (foldchange, "mean"),
            foldchange + "_std": (foldchange, "std"),
            foldchange + "_count": (foldchange, "count"),
        }
        df_stats = df.groupby(biorep_col, as_index=False).agg(**stats_dict)
        df_stats[foldchange + "_ste"] = df_stats[foldchange + "_std"] / np.sqrt(
            df_stats[foldchange + "_count"]
        )

        # Drop metadata
        if not keep_metadata:
            df_stats = df_stats.drop([foldchange + "_count"], axis="columns")

        df = df.merge(df_stats, on=biorep_col)

    # Drop metadata after calculation
    if not keep_metadata:
        pattern = r"^(delta|2\^\([^)]*\)$)"
        metadata_cols = df.filter(regex=pattern).columns
        df = (
            df.drop(metadata_cols, axis="columns")
            .drop_duplicates(subset=biorep_col, keep="first")
            .reset_index(drop=True)
        )

    return df


if __name__ == "__main__":
    # TODO: import testing methods instead of writing directly into this script
    pass
