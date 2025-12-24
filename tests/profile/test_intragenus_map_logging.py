#!/usr/bin/env python3
import contextlib
import importlib.util
import io
import sys
import tempfile
import types
from pathlib import Path


def load_profile_module(repo_root: Path):
    profile_dir = repo_root / "src" / "profile"
    profile_path = profile_dir / "profile.py"
    taxonomy_path = profile_dir / "taxonomy_utils.py"

    pkg_name = "chimera_profile_pkg2"
    pkg = types.ModuleType(pkg_name)
    pkg.__path__ = [str(profile_dir)]
    sys.modules[pkg_name] = pkg

    tax_spec = importlib.util.spec_from_file_location(
        f"{pkg_name}.taxonomy_utils", taxonomy_path
    )
    if tax_spec is None or tax_spec.loader is None:
        raise RuntimeError("failed to load taxonomy_utils.py")
    tax_module = importlib.util.module_from_spec(tax_spec)
    sys.modules[f"{pkg_name}.taxonomy_utils"] = tax_module
    tax_spec.loader.exec_module(tax_module)

    prof_spec = importlib.util.spec_from_file_location(
        f"{pkg_name}.profile", profile_path
    )
    if prof_spec is None or prof_spec.loader is None:
        raise RuntimeError("failed to load profile.py")
    prof_module = importlib.util.module_from_spec(prof_spec)
    sys.modules[f"{pkg_name}.profile"] = prof_module
    prof_spec.loader.exec_module(prof_module)
    return prof_module


def main() -> int:
    repo_root = Path(__file__).resolve().parents[2]
    profile = load_profile_module(repo_root)

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        taxonomy_info = tmp_path / "taxonomy_info.tsv"
        lines = [
            "d1\t\t domain\t d__Bacteria",
            "g1\td1\t genus\t g__Alpha",
            "s1\tg1\t species\t s__Alpha s1",
            "s2\tg1\t species\t s__Alpha s2",
        ]
        for i in range(3, 21):
            gid = f"g{i}"
            sid = f"s{i}"
            lines.append(f"{gid}\td1\t genus\t g__G{i}")
            lines.append(f"{sid}\t{gid}\t species\t s__G{i} s{i}")
        taxonomy_info.write_text("\n".join(lines) + "\n", encoding="utf-8")

        classify_path = tmp_path / "ChimeraClassify.tsv"
        with classify_path.open("w", encoding="utf-8") as fh:
            for i in range(10):
                fh.write(f"read_s2_{i}\ts2:1\tPOST_TOPK=s2:0.9,s1:0.1\n")
            fh.write("read_fix\ts1:1\tPOST_TOPK=s1:0.51,s2:0.49\n")
            for i in range(3, 21):
                fh.write(f"read_{i}\ts{i}:1\n")

        out_path = tmp_path / "ChimeraAbundance"
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            profile.process_file(
                input_files=[str(classify_path)],
                output_file=str(out_path),
                taxonomy_kind="gtdb",
                taxonomy_info=str(taxonomy_info),
                abundance_mode="soft_seq",
                hd_intragenus_map=True,
            )
        output = buf.getvalue()
        if "Intragenus MAP" not in output:
            sys.stderr.write("expected intragenus MAP log line but not found\n")
            return 1
        if "corrected=0" in output:
            sys.stderr.write("expected intragenus MAP to correct at least one read\n")
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
