#!/usr/bin/env python3
import importlib.util
import re
import sys
import tempfile
import types
from pathlib import Path


def load_profile_module(repo_root: Path):
    profile_dir = repo_root / "src" / "profile"
    profile_path = profile_dir / "profile.py"
    taxonomy_path = profile_dir / "taxonomy_utils.py"

    pkg_name = "chimera_profile_pkg"
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


def test_dump_post_topk_default(repo_root: Path) -> None:
    header = repo_root / "src" / "classify" / "classifyConfig.hpp"
    text = header.read_text(encoding="utf-8", errors="ignore")
    m = re.search(r"dump_post_topk\s*=\s*([0-9]+)", text)
    assert m, "dump_post_topk default not found"
    value = int(m.group(1))
    assert value == 256, f"dump_post_topk default expected 256, got {value}"


def test_low_diversity_uses_top1(repo_root: Path) -> None:
    profile = load_profile_module(repo_root)
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        taxonomy_info = tmp_path / "taxonomy_info.tsv"
        lines = [
            "d1\t\t domain\t d__Bacteria",
            "g1\td1\t genus\t g__Alpha",
        ]
        for i in range(1, 21):
            lines.append(f"s{i}\tg1\t species\t s__Alpha_{i}")
        taxonomy_info.write_text("\n".join(lines) + "\n", encoding="utf-8")

        classify_path = tmp_path / "ChimeraClassify.tsv"
        with classify_path.open("w", encoding="utf-8") as fh:
            for i in range(10):
                tokens = ["s1:1"] + [f"s{j}:1" for j in range(2, 21)]
                fh.write(f"read{i}\t" + "\t".join(tokens) + "\n")

        out_path = tmp_path / "ChimeraAbundance"
        profile.process_file(
            input_files=[str(classify_path)],
            output_file=str(out_path),
            taxonomy_kind="gtdb",
            taxonomy_info=str(taxonomy_info),
            abundance_mode="soft_seq",
        )

        out_file = out_path.with_suffix(".tsv")
        header = out_file.read_text(encoding="utf-8", errors="ignore").splitlines()[0]
        assert header.strip() == "# abundance_mode=top1", (
            f"expected low-diversity auto switch to top1, got: {header}"
        )


if __name__ == "__main__":
    repo_root = Path(__file__).resolve().parents[2]
    try:
        test_dump_post_topk_default(repo_root)
        test_low_diversity_uses_top1(repo_root)
    except AssertionError as exc:
        print(f"TEST FAILED: {exc}")
        sys.exit(1)
    print("tests/profile/test_profile_low_diversity.py: all passed")
