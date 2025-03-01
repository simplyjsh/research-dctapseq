from typing import Literal

import pandas as pd
from IPython.display import display


def set_defaultoptions(
    pandas,
    max_cols: int | None = None,
    max_rows: int | None = 10,
    expand_frame_repr: bool = False,
    large_repr: Literal["truncate", "info"] = "truncate",
):
    """
    Set default options for panda DataFrames.
    """

    # Display options
    pandas.set_option("display.max_columns", max_cols)

    """
    Note:
        Pandas DataFrame invokes max_rows only when value exceeds
        the instantiated df's total rows, otherwise, it defaults to 
        min_rows.

    e.g.:
        df.shape # outputs: (20, 10), (rows, cols)

        if max_rows <= 20, then df will only display based on min_rows
        which by defeault is 10. So only 10 rows will display.
        max_rows will setting will be used when >20.

    details:
        https://stackoverflow.com/questions/57860775/pandas-pd-options-display-max-rows-not-working-as-expected
    """
    pandas.set_option("display.max_rows", max_rows)
    pandas.set_option("display.min_rows", max_rows)

    pd.set_option("display.expand_frame_repr", expand_frame_repr)
    pd.set_option("display.large_repr", large_repr)


def displaydf_full():
    """
    Helper method to display a DataFrame in full, call with "with" context.
    e.g.
    with displaydf_full():
        display(df)
    """
    return pd.option_context("display.max_rows", None, "display.max_columns", None)
