"""
SceptreFormatterToolsCLI
------------------------
Entry point into command line interface.
"""

# -----------------------------------------------------------------------
# INFO: Imports
import click

from sceptretoolscli.core.formatters import fmt_inputs

# -----------------------------------------------------------------------
# INFO: Entry method


@click.group()
def cli() -> None:
    """
    Sceptre Formatter Tools

    This CLI tool runs helpful for formating sceptre inputs.
    Currently in developement and only supports formatting sceptre inputs.
    """
    pass


# -----------------------------------------------------------------------
# INFO: Entry to all other commands via click command groups

cli.add_command(fmt_inputs)
