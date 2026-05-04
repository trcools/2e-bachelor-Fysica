from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    generic_script = repo_root / "pdf-naar-md-bestand" / "pdf_naar_md.py"
    input_dir = Path(__file__).resolve().parent / "Samenvatting-REM"

    if not generic_script.exists():
        raise SystemExit(f"Converter script not found: {generic_script}")

    cmd = [
        sys.executable,
        str(generic_script),
        "--input-dir",
        str(input_dir),
    ]
    raise SystemExit(subprocess.call(cmd))


if __name__ == "__main__":
    main()
