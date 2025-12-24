#!/usr/bin/env python3
import importlib.util
import sys
import tempfile
import types
from pathlib import Path


def load_profile_module(repo_root: Path):
    profile_dir = repo_root / "src" / "profile"
    profile_path = profile_dir / "profile.py"
    taxonomy_path = profile_dir / "taxonomy_utils.py"

    pkg_name = "chimera_profile_prefilter_pkg"
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


def parse_species_section(path: Path):
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    species = []
    in_section = False
    for line in lines:
        if line.startswith("## ") and "Species" in line:
            in_section = True
            continue
        if in_section:
            if line.startswith("## "):
                break
            if not line.strip() or line.startswith("Taxon\t"):
                continue
            species.append(line.split("\t", 1)[0])
    return species


def test_species_prefilter_uses_eff_total_when_unclassified_high(repo_root: Path) -> None:
    profile = load_profile_module(repo_root)
    orig = (
        profile.MIN_REL_ABUNDANCE,
        profile.MIN_EVIDENCE,
        profile.HIGH_DIVERSITY_MIN_REL_ABUNDANCE,
        profile.HIGH_DIVERSITY_FINAL_REL_CUTOFF,
    )
    try:
        profile.MIN_REL_ABUNDANCE = 50.0
        profile.MIN_EVIDENCE = 0
        profile.HIGH_DIVERSITY_MIN_REL_ABUNDANCE = 0.0
        profile.HIGH_DIVERSITY_FINAL_REL_CUTOFF = 0.0

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            taxonomy_info = tmp_path / "taxonomy_info.tsv"
            lines = [
                "d1\t\t domain\t d__Bacteria",
                "g1\td1\t genus\t g__Alpha",
                "s1\tg1\t species\t s__Alpha_1",
            ]
            taxonomy_info.write_text("\n".join(lines) + "\n", encoding="utf-8")

            classify_path = tmp_path / "ChimeraClassify.tsv"
            with classify_path.open("w", encoding="utf-8") as fh:
                fh.write("read0\ts1:1\n")
                for i in range(1, 10):
                    fh.write(f"read{i}\tunclassified\n")

            out_path = tmp_path / "ChimeraAbundance"
            profile.process_file(
                input_files=[str(classify_path)],
                output_file=str(out_path),
                taxonomy_kind="gtdb",
                taxonomy_info=str(taxonomy_info),
                abundance_mode="top1",
                hd_species_head_mass=None,
                hd_genus_mode="none",
                hd_shape_mode="none",
                hd_intragenus_map=False,
            )

            species = parse_species_section(out_path.with_suffix(".tsv"))
            kept = [taxon for taxon in species if taxon != "unclassified"]
            assert kept, "expected species retained after prefilter"
    finally:
        (
            profile.MIN_REL_ABUNDANCE,
            profile.MIN_EVIDENCE,
            profile.HIGH_DIVERSITY_MIN_REL_ABUNDANCE,
            profile.HIGH_DIVERSITY_FINAL_REL_CUTOFF,
        ) = orig


if __name__ == "__main__":
    repo_root = Path(__file__).resolve().parents[2]
    try:
        test_species_prefilter_uses_eff_total_when_unclassified_high(repo_root)
    except AssertionError as exc:
        print(f"TEST FAILED: {exc}")
        sys.exit(1)
    print("tests/profile/test_profile_prefilter_denom.py: all passed")
