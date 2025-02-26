from pathlib import Path


class Filepaths:
    ROOT: Path = Path(__file__).parent.parent.resolve()
    DATADIR: Path = ROOT / "data"
