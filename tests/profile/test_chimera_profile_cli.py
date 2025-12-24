#!/usr/bin/env python3
import subprocess
import sys
from pathlib import Path


def main() -> int:
    repo_root = Path(__file__).resolve().parents[2]
    chimera_py = repo_root / "chimera.py"
    cmd = [sys.executable, str(chimera_py), "profile", "-h"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write("chimera.py profile -h should succeed but failed.\n")
        sys.stderr.write(proc.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
