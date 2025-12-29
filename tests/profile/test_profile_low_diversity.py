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


def test_low_diversity_prefers_top1_even_with_post_topk(repo_root: Path) -> None:
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
                post = ",".join([f"s{j}:0.05" for j in range(1, 6)])
                fh.write(f"read{i}\t" + "\t".join(tokens) + f"\tPOST_TOPK={post}\n")

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


def test_low_diversity_falls_back_to_top1_without_post_topk(repo_root: Path) -> None:
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
            f"expected low-diversity fallback to top1, got: {header}"
        )


def test_low_diversity_rel_threshold_keeps_tail(repo_root: Path) -> None:
    profile = load_profile_module(repo_root)
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        taxonomy_info = tmp_path / "taxonomy_info.tsv"
        lines = [
            "d1\t\t domain\t d__Bacteria",
            "g1\td1\t genus\t g__Alpha",
        ]
        for i in range(1, 11):
            lines.append(f"s{i}\tg1\t species\t s__Alpha_{i}")
        taxonomy_info.write_text("\n".join(lines) + "\n", encoding="utf-8")

        classify_path = tmp_path / "ChimeraClassify.tsv"
        with classify_path.open("w", encoding="utf-8") as fh:
            idx = 0
            for _ in range(40):
                fh.write(f"read{idx}\ts1:1\n")
                idx += 1
            for _ in range(38):
                fh.write(f"read{idx}\ts2:1\n")
                idx += 1
            for _ in range(15):
                fh.write(f"read{idx}\ts3:1\n")
                idx += 1
            for taxid in range(4, 11):
                fh.write(f"read{idx}\ts{taxid}:1\n")
                idx += 1

        out_path = tmp_path / "ChimeraAbundance"
        profile.process_file(
            input_files=[str(classify_path)],
            output_file=str(out_path),
            taxonomy_kind="gtdb",
            taxonomy_info=str(taxonomy_info),
            abundance_mode="soft_seq",
        )

        out_file = out_path.with_suffix(".tsv")
        lines = out_file.read_text(encoding="utf-8", errors="ignore").splitlines()
        species_lines = []
        in_species = False
        for line in lines:
            if line.startswith("## Species Level ##"):
                in_species = True
                continue
            if in_species and line.startswith("## ") and line.endswith(" ##"):
                break
            if in_species and line and not line.startswith("#") and not line.startswith("Taxon"):
                species_lines.append(line)

        taxon_names = [ln.split("\t", 1)[0] for ln in species_lines]
        expected = [f"Alpha_{i}" for i in range(1, 11)]
        assert taxon_names == expected, (
            f"expected rel threshold to keep tail, got: {taxon_names}"
        )


def test_low_diversity_post_topk_highprob_rescue_adds_new_genus_species(
    repo_root: Path,
) -> None:
    profile = load_profile_module(repo_root)
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        taxonomy_info = tmp_path / "taxonomy_info.tsv"
        lines = [
            "d1\t\t domain\t d__Bacteria",
            "g1\td1\t genus\t g__Alpha",
            "g2\td1\t genus\t g__Beta",
            "s1\tg1\t species\t s__Alpha_1",
            "t1\tg2\t species\t s__Beta_1",
        ]
        taxonomy_info.write_text("\n".join(lines) + "\n", encoding="utf-8")

        classify_path = tmp_path / "ChimeraClassify.tsv"
        with classify_path.open("w", encoding="utf-8") as fh:
            for i in range(60):
                fh.write(
                    f"read{i}\ts1:1\tPOST_TOPK=t1:0.30,s1:0.70\n"
                )

        out_path = tmp_path / "ChimeraAbundance"
        profile.process_file(
            input_files=[str(classify_path)],
            output_file=str(out_path),
            taxonomy_kind="gtdb",
            taxonomy_info=str(taxonomy_info),
            abundance_mode="soft_seq",
        )

        out_file = out_path.with_suffix(".tsv")
        lines = out_file.read_text(encoding="utf-8", errors="ignore").splitlines()
        species_lines = []
        in_species = False
        for line in lines:
            if line.startswith("## Species Level ##"):
                in_species = True
                continue
            if in_species and line.startswith("## ") and line.endswith(" ##"):
                break
            if in_species and line and not line.startswith("#") and not line.startswith("Taxon"):
                species_lines.append(line)

        by_name = {ln.split("\t", 1)[0]: ln for ln in species_lines}
        assert "Beta_1" in by_name, "expected high-prob POST_TOPK rescue to add Beta_1"
        parts = by_name["Beta_1"].split("\t")
        assert len(parts) >= 2
        assert int(float(parts[1])) >= int(profile.MIN_EVIDENCE)


def test_low_diversity_headmass_prunes_tiny_tail(repo_root: Path) -> None:
    profile = load_profile_module(repo_root)
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        taxonomy_info = tmp_path / "taxonomy_info.tsv"
        lines = [
            "d1\t\t domain\t d__Bacteria",
            "g1\td1\t genus\t g__Alpha",
        ]
        for i in range(1, 10):
            lines.append(f"s{i}\tg1\t species\t s__Alpha_{i}")
        taxonomy_info.write_text("\n".join(lines) + "\n", encoding="utf-8")

        classify_path = tmp_path / "ChimeraClassify.tsv"
        counts = {
            1: 5000,
            2: 3000,
            3: 1500,
            4: 300,
            5: 100,
            6: 40,
            7: 30,
            8: 20,
            9: 10,
        }
        with classify_path.open("w", encoding="utf-8") as fh:
            idx = 0
            for taxid, n in counts.items():
                for _ in range(n):
                    fh.write(f"read{idx}\ts{taxid}:1\n")
                    idx += 1
        assert idx == 10000

        out_path = tmp_path / "ChimeraAbundance"
        profile.process_file(
            input_files=[str(classify_path)],
            output_file=str(out_path),
            taxonomy_kind="gtdb",
            taxonomy_info=str(taxonomy_info),
            abundance_mode="soft_seq",
        )

        out_file = out_path.with_suffix(".tsv")
        lines = out_file.read_text(encoding="utf-8", errors="ignore").splitlines()
        species_lines = []
        in_species = False
        for line in lines:
            if line.startswith("## Species Level ##"):
                in_species = True
                continue
            if in_species and line.startswith("## ") and line.endswith(" ##"):
                break
            if in_species and line and not line.startswith("#") and not line.startswith("Taxon"):
                species_lines.append(line)

        taxon_names = [ln.split("\t", 1)[0] for ln in species_lines]
        assert taxon_names == [f"Alpha_{i}" for i in range(1, 9)], (
            f"expected low-div head-mass to prune tiny tail, got: {taxon_names}"
        )


def test_low_diversity_post_topk_rescue_skips_bacterium_names(repo_root: Path) -> None:
    profile = load_profile_module(repo_root)
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        taxonomy_info = tmp_path / "taxonomy_info.tsv"
        lines = [
            "d1\t\t domain\t d__Bacteria",
            "g1\td1\t genus\t g__Alpha",
            "g2\td1\t genus\t g__Hyphomicrobiales",
            "s1\tg1\t species\t s__Alpha_1",
            "t1\tg2\t species\t s__Hyphomicrobiales bacterium",
        ]
        taxonomy_info.write_text("\n".join(lines) + "\n", encoding="utf-8")

        classify_path = tmp_path / "ChimeraClassify.tsv"
        with classify_path.open("w", encoding="utf-8") as fh:
            for i in range(60):
                fh.write(f"read{i}\ts1:1\tPOST_TOPK=t1:0.30,s1:0.70\n")

        out_path = tmp_path / "ChimeraAbundance"
        profile.process_file(
            input_files=[str(classify_path)],
            output_file=str(out_path),
            taxonomy_kind="gtdb",
            taxonomy_info=str(taxonomy_info),
            abundance_mode="soft_seq",
        )

        out_file = out_path.with_suffix(".tsv")
        text = out_file.read_text(encoding="utf-8", errors="ignore")
        assert "Hyphomicrobiales bacterium" not in text


if __name__ == "__main__":
    repo_root = Path(__file__).resolve().parents[2]
    try:
        test_dump_post_topk_default(repo_root)
        test_low_diversity_prefers_top1_even_with_post_topk(repo_root)
        test_low_diversity_falls_back_to_top1_without_post_topk(repo_root)
        test_low_diversity_rel_threshold_keeps_tail(repo_root)
        test_low_diversity_post_topk_highprob_rescue_adds_new_genus_species(repo_root)
        test_low_diversity_headmass_prunes_tiny_tail(repo_root)
        test_low_diversity_post_topk_rescue_skips_bacterium_names(repo_root)
    except AssertionError as exc:
        print(f"TEST FAILED: {exc}")
        sys.exit(1)
    print("tests/profile/test_profile_low_diversity.py: all passed")
