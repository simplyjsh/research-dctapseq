"""
SceptreFormatterToolsCLI
------------------------
Entry point into command line interface.
"""

# -----------------------------------------------------------------------
# INFO: Imports
import click
import logging
import numpy as np
import pandas as pd

import sys
from pathlib import Path

# -----------------------------------------------------------------------
# INFO: Constants

# TODO: COLMAP dict for single truth access of column names.


# -----------------------------------------------------------------------
# INFO: Helper methods


def extract_attr(gtf: pd.DataFrame, fields: list) -> pd.DataFrame:
    """
    Extracts specified fields from the GTF attribute column.

    Parameters
    ----------
        gtf (pd.DataFrame): DataFrame containing a column named 'attribute'.
        fields (list of str): Field column headers to extract (e.g. ["gene_id", "gene_name"]).

    Returns
    -------
        results (pd.DataFrame): A DataFrame with one column per requested extracted attribute.
    """
    # Intialize results variable, allocate just enough memory
    results = pd.DataFrame(index=gtf.index)

    # Use regex and pandas' str extract to pull requested field attributes
    for field in fields:
        pattern = rf'{field} "([^"]+)"'
        results[field] = gtf["attribute"].str.extract(pattern)

    return results


# TODO: Write Function Docs
def clean_gene_id(gtfgenes: pd.DataFrame):
    colmap = {
        "s": "stable_id",
        "t": "tail",
        "j": "major_v",
        "n": "minor_v",
        "r": "remainder",
        "p": "pseudoautosomal",
        "e": "ensembl_ver",
    }

    # Split gene_id into stable ENSG ID and trailing ensembl version
    split1 = [colmap.get(key) for key in ["s", "t"]]
    df_splits = (
        gtfgenes.loc[:, "gene_id"]
        .copy()
        .str.split(".", n=1, expand=True)
        .set_axis(split1, axis=1)
    )

    # Split tail into major and trailing remainder
    split2 = [colmap.get(key) for key in ["j", "r"]]
    df_splits.loc[:, split2] = (
        df_splits.loc[:, colmap.get("t")]
        .str.split("_", n=1, expand=True)
        .set_axis(split2, axis=1)
    )

    # Split remainder into minor and pseudoautosomal
    split3 = [colmap.get(key) for key in ["n", "p"]]
    df_splits.loc[:, split3] = (
        df_splits.loc[:, colmap.get("r")]
        .str.split("_", n=1, expand=True)
        .set_axis(split3, axis=1)
    )

    # Rebuild ensembl_ver with major and minor version, if minor is None then minor version is zero.
    e, j, n = [colmap.get(key) for key in ["e", "j", "n"]]
    df_splits.loc[:, e] = (
        df_splits.loc[:, j].astype(str)
        + "."
        + df_splits.loc[:, n].fillna(0).astype(int).astype(str)
    ).astype(float)

    # Format dataframe for cols of interest
    df_cleaned = df_splits.loc[:, [colmap.get(key) for key in ["s", "e"]]].copy()

    return df_cleaned


# TODO: Write Function Docs
def find_ambigious_genes(gtfgenes: pd.DataFrame):
    # Parse genes that have the same gene_name for different multiple gene_ids for the same ensembl_version
    mask = gtfgenes.duplicated(subset=["gene_name", "ensembl_ver"], keep=False)
    df_ambigious = gtfgenes.loc[mask].copy()
    df_filtered = gtfgenes.loc[~mask].copy()

    # Sanity Check
    t, a, f = [len(gtfgenes), len(df_ambigious), len(df_filtered)]
    assert t == a + f, (
        f"Finding ambigious genes failed unexpectedly: total={t}, "
        f"ambigious={a}, filtered={f}"
    )

    # Format ambigious genes dataframe
    df_ambigious = (
        df_ambigious.dropna()
        .sort_values(["seqname", "gene_name"])
        .reset_index(drop=True)
    )

    logging.info(
        f"Found {len(df_ambigious.loc[:, 'gene_name'].unique())} unique gene "
        f"symbols with ambigious GENCODE gene annotations."
    )
    return df_filtered, df_ambigious


# TODO: Update Function Docs
def get_gtfgenes(gtf: str):
    """
    Load a GTF annotation file and return a DataFrame of gene features.

    This function reads a compressed GTF file, extracts only rows where the
    feature type is "gene", and parses the attributes column to retrieve
    `gene_id` and `gene_name`. `gene_id` goes through further formatting
    to obtain the latest Ensembl gene identifiers.

    The resulting DataFrame includes genomic coordinates and
    the ensembl identifiers for each gene.

    Parameters
    ----------
    gtf (str or path-like): Path to the GTF annotation file (supports `.gtf` or `.gtf.gz` format).

    Returns
    -------
    df_gtfgenes (pd.DataFrame):
        A DataFrame containing one row per gene with the following columns:

        - ``chr``           (str): Chromosome name (e.g. "chr1", "chrX").
        - ``start``         (int): Start coordinate of the gene.
        - ``end``           (int): End coordinate of the gene.
        - ``gene_name``     (str): Gene symbol (e.g. "DDX11L1").
        - ``gene_id``       (str): Ensembl gene identifier (e.g. "ENSG00000223972").
        - ``ensembl_ver``   (str): Ensembl identifier version used (e.g. "5.2").

    Notes
    -----
    - This function expects the GTF to follow the standard 9-column format
      with attributes in the last column. See https://www.ensembl.org/info/website/upload/gff.html.
    - Only features annotated as ``gene`` are retained.
    - Attributes are parsed using :func:`extract_attr`.

    Examples
    --------
    >>> df_genes = get_gtfgenes("gencode.v32lift37.annotation.gtf.gz")
    >>> df_genes.head()
         chr   start     end            gene_id gene_name ensembl_ver
    0   chr1   11869   14409    ENSG00000223972   DDX11L1         5.2
    1   chr1   14362   29806    ENSG00000227232    WASH7P         5.2
    """
    logging.info("Processing GENCODE GTF file.")

    # GTF standard colnames
    gtfcolnames = [
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ]

    # Colmap for single truth access of scoped variables
    colmap = {"i": "gene_id", "n": "gene_name", "s": "stable_id", "e": "ensembl_ver"}

    logging.info("Reading GTF file...")
    df_gtf = pd.read_csv(
        gtf, sep="\t", comment="#", header=None, names=gtfcolnames, compression="gzip"
    )

    logging.info("Formatting GTF file...")

    # Filter GTF for gene features
    df_gtfgenes = df_gtf.loc[df_gtf["feature"] == "gene", :]

    # Filter out unplaced or scaffold annotations
    df_gtfgenes = df_gtfgenes.loc[df_gtfgenes["seqname"].str.startswith("chr"), :]

    # Filter out mitochondrial genes
    df_gtfgenes = df_gtfgenes.loc[~df_gtfgenes["seqname"].str.startswith("chrM"), :]

    # Extract gene attribute metadata, interested in gene_id and gene_name
    fields = [colmap.get(key) for key in ["i", "n"]]
    df_gtfgenes.loc[:, fields] = extract_attr(df_gtfgenes, fields)

    # Split ensembl id from ensembl version.
    ensg_splits = [colmap.get(key) for key in ["s", "e"]]
    df_gtfgenes.loc[:, ensg_splits] = clean_gene_id(df_gtfgenes)

    # Find ambigious gene_ids, these genes have the same gene_name for multiple gene_ids and ensembl_version
    df_gtfgenes, ambigious_genes = find_ambigious_genes(df_gtfgenes)

    # Select latest version of ensembl id for duplicate gene listings by group sorting and selecting first
    sort_subset = [colmap.get(key) for key in ["n", "e"]]
    df_gtfgenes = df_gtfgenes.sort_values(
        sort_subset, ascending=[True, False]
    ).drop_duplicates(colmap.get("n"), keep="first")

    # Format DataFrame, keep specified cols and rename seqname -> chr and stable_id -> gene_id
    keepcols = ["seqname", "start", "end"] + [
        colmap.get(key) for key in ["n", "s", "e"]
    ]
    df_gtfgenes = (
        df_gtfgenes.loc[:, keepcols]
        .rename(columns={"seqname": "chr", colmap.get("s"): colmap.get("i")})
        .dropna()
        .sort_values(["chr"] + [colmap.get("n")])
        .reset_index(drop=True)
        .copy()
    )

    logging.info("STEP DONE.")
    return df_gtfgenes, ambigious_genes


# TODO: Write Function Doc
def standardize_to_gene_id(
    gene_file: str, gtf_gene_refs: pd.DataFrame, gtf_ambigious_refs: np.ndarray
):
    logging.info("Processing gene file input.")

    logging.info("Reading gene file...")
    s_genes = pd.read_csv(gene_file)["gene"]

    logging.info("Formatting gene file...")

    # Split genes from gene_file into series that either begins with ENSG or not
    mask = s_genes.str.startswith("ENSG", na=False)
    gene_id = pd.Series(s_genes.where(mask).dropna().reset_index(drop=True))
    gene_name = pd.Series(s_genes.where(~mask).dropna().reset_index(drop=True))

    # TEST: Todo make a gene_file that has a fake gene id not found in gtf_gene_refs to test this
    # Checks if any gene_ids that are not found within the user provided GENCODE annotations
    mask_nomatch_id = ~gene_id.isin(gtf_gene_refs.loc[:, "gene_id"].unique())
    gene_nomatch_id = pd.Series(gene_id[mask_nomatch_id])
    if not gene_nomatch_id.empty:
        logging.warning(
            f"[WARNING] The following ENSG IDs were not found in the provided GENCODE gene annotations "
            f"(potential version mismatch): {', '.join(map(str, gene_nomatch_id.tolist()))}",
        )

    # TEST: Todo make a gene_file that has a fake gene found in ambiguous to test this
    # Checks if any gene_name will not be mapped due to ambigious GENCODE annotations
    mask_ambig = gene_name.isin(gtf_ambigious_refs)
    gene_ambig = pd.Series(gene_name[mask_ambig])
    if not gene_ambig.empty:
        logging.warning(
            f"[WARNING] The following gene symbols are ambiguous in GENCODE and will NOT be mapped: "
            f"{', '.join(map(str, gene_ambig.tolist()))}",
        )

    # Remove ambigious genes and gene_ids from gene_ids and gene_name
    gene_id = pd.Series(gene_id[~mask_nomatch_id])
    gene_name = pd.Series(gene_name[~mask_ambig])

    # Setup mapping for conversion to gene ids, then convert. Default follows ensembl id standard.
    name_to_id = gtf_gene_refs.set_index("gene_name")["gene_id"]
    mapped_ids = gene_name.map(
        lambda gene_name: name_to_id.to_dict().get(gene_name, pd.NA)
    )
    mask_mapped = mapped_ids.notna()
    matched_gene_names = pd.Series(gene_name[mask_mapped])
    leftout_gene_names = pd.Series(gene_name[~mask_mapped])

    if not leftout_gene_names.empty:
        logging.warning(
            f"[WARNING] The following gene symbols had no match with any gene features in the "
            f"provided GENCODE annotation GTF file: {', '.join(map(str, leftout_gene_names.tolist()))}"
        )

    # TEST: Better methodology for checking?
    # Check gene inputs (at this point we should have matched gene names and ids to convert or it existed already)
    if gene_id.empty and matched_gene_names.empty:
        logging.error(
            "[ERROR] No matching gene ids found for the given GENCODE annotations. "
            "Please check the GENCODE annotation version used or the expected "
            "response gene file input."
        )
        sys.exit(1)

    # Format mappped gene ids
    mapped_gene_ids = pd.DataFrame(
        {
            "gene_name": pd.Series(gene_name[mask_mapped]).reset_index(drop=True),
            "gene_id": pd.Series(mapped_ids[mask_mapped]).reset_index(drop=True),
        }
    )

    logging.info("STEP DONE.")
    return mapped_gene_ids, gene_nomatch_id, gene_ambig, leftout_gene_names


# TODO: Write Function Doc and methodology
def calculate_relative_dist(
    gtf_split: pd.DataFrame, guide: pd.Series, dist_cutoff: int
):
    # Cal relative coords, subtracting from the guide's coord start
    gtf_split.loc[:, ["rel_coord_start", "rel_coord_end"]] = (
        gtf_split.loc[:, "start":"end"]
        .sub(guide.start)
        .set_axis(["rel_coord_start", "rel_coord_end"], axis=1)
    )

    # Classify relative position
    cond_left = (gtf_split["rel_coord_start"] < 0) & (gtf_split["rel_coord_end"] < 0)
    cond_right = (gtf_split["rel_coord_start"] > 0) & (gtf_split["rel_coord_end"] > 0)
    cond_start = gtf_split["rel_coord_start"].eq(0)
    cond_end = gtf_split["rel_coord_end"].eq(0)

    # Assign classification
    gtf_split["rel_pos"] = np.select(
        [cond_left, cond_right, cond_start, cond_end],
        ["abs_left", "abs_right", "abs_right", "abs_left"],
        default="check_both",
    )

    # Determine if the guide is near a particular gene based on dist_cutoff
    cond_near_left = (gtf_split["rel_pos"] == "abs_left") & (
        gtf_split["rel_coord_start"] > -dist_cutoff
    )
    cond_near_right = (gtf_split["rel_pos"] == "abs_right") & (
        gtf_split["rel_coord_end"] < dist_cutoff + 20
    )
    cond_near_both = (gtf_split["rel_pos"] == "check_both") & (
        (gtf_split["rel_coord_start"] > -dist_cutoff)
        & (gtf_split["rel_coord_end"] < dist_cutoff + 20)
    )

    gtf_split["nearby"] = np.select(
        [cond_near_left, cond_near_right, cond_near_both],
        [True, True, True],
        default=False,
    )

    filtered = gtf_split[gtf_split.loc[:, "nearby"]].copy()
    all_results = gtf_split.copy()

    return filtered, all_results


# -----------------------------------------------------------------------
# INFO: MAIN METHOD


def get_guidegene_pairing_file(
    gene_file,
    guide_file,
    dist_cutoff,
    gtf,
    output_path=Path(__file__).resolve().parent.parent.parent / "results",
):
    # Intialize results variables
    leftouts = dict()  # TODO:
    guide_gene_pairing = pd.DataFrame()

    # Read files & format files
    df_gtfgenes, df_ambigious_genes = get_gtfgenes(gtf)
    df_genes, nomatch_ensg, unmapped_ambig, unmapped_genes = standardize_to_gene_id(
        gene_file,
        df_gtfgenes.loc[:, ["gene_id", "gene_name"]],
        df_ambigious_genes.loc[:, "gene_name"].unique(),
    )
    df_guides = pd.read_csv(guide_file)[["name", "chr", "start", "end"]]

    # Filter gtf annotations to genes of interest
    df_gtfgenes_filtered = df_gtfgenes.loc[
        df_gtfgenes.loc[:, "gene_id"].isin(df_genes.loc[:, "gene_id"])
    ].reset_index(drop=True)

    # Split annotations by chromosomes and store as dict()
    gtf_splits = {
        chr_splits: group.reset_index(drop=True)
        for chr_splits, group in df_gtfgenes_filtered.groupby("chr")
    }

    # Split guides by chromosomes and store as dict()
    guide_splits = {
        chr_splits: group.reset_index(drop=True)
        for chr_splits, group in df_guides.groupby("chr")
    }

    pass


# -----------------------------------------------------------------------
# INFO: CLI


@click.group(name="fmt-inputs")
def fmt_inputs() -> None:
    pass


# TODO: Allow csv and tsv formats for gene and guide inputs.
@fmt_inputs.command(name="guide-gene-pairing")
def guide_gene_pairing() -> None:
    pass
