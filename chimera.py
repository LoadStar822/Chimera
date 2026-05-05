import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

from src.profile import profile


def default_threads():
    cpu = os.cpu_count()
    if not cpu or cpu <= 0:
        return 1
    return min(cpu, 192)


def get_downloader():
    from src.download import download as downloader
    return downloader


def get_chimera_path():
    env_path = os.environ.get("CHIMERA_BIN")
    if env_path:
        chimera_path = Path(env_path)
        if chimera_path.is_file() and os.access(chimera_path, os.X_OK):
            return str(chimera_path)
        raise FileNotFoundError(
            f"CHIMERA_BIN points to a missing or non-executable file: {env_path}"
        )
    # Use which command to find the chimera program path
    chimera_path = shutil.which("Chimera")
    if chimera_path is None:
        raise FileNotFoundError(
            "Cannot find 'Chimera' executable. Please ensure it is installed and in your PATH."
        )
    return chimera_path


def _infer_profile_taxonomy_paths(database_path: str):
    db = Path(database_path).resolve()
    search_dirs = []
    if db.parent not in search_dirs:
        search_dirs.append(db.parent)
    stem_path = db.with_suffix("")
    if stem_path.is_dir() and stem_path not in search_dirs:
        search_dirs.append(stem_path)

    taxonomy_meta = None
    taxonomy_info = None
    for base in search_dirs:
        if taxonomy_meta is None:
            candidate = base / "taxonomy.meta"
            if candidate.exists():
                taxonomy_meta = str(candidate)
        if taxonomy_info is None:
            candidate = base / "tax.info"
            if candidate.exists():
                taxonomy_info = str(candidate)
    return taxonomy_meta, taxonomy_info


def min_length_type(value: str):
    lowered = value.lower()
    if lowered == "auto":
        return "auto"
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Minimum length must be an integer or 'auto'") from exc
    if parsed < 0:
        raise argparse.ArgumentTypeError("Minimum length must be >= 0 or 'auto'")
    return parsed


def parse_arguments():
    if len(sys.argv) > 1 and sys.argv[1] == "download":
        parser = argparse.ArgumentParser(
            description="Chimera - A versatile tool for metagenomic classification"
        )
        subparsers = parser.add_subparsers(dest="command", required=False)
        subparsers.add_parser(
            "download", help="Download NCBI database sequences and resources"
        )

        args = parser.parse_args(sys.argv[1:2])
        return args

    parser = argparse.ArgumentParser(
        description="Chimera - A versatile tool for metagenomic classification"
    )

    parser.add_argument(
        "-v", "--version", action="store_true", help="Show version information"
    )

    subparsers = parser.add_subparsers(dest="command", required=False)

    # Download subcommand
    subparsers.add_parser(
        "download", help="Download NCBI database sequences and resources"
    )

    # Build subcommand
    build_parser = subparsers.add_parser(
        "build", help="Build a taxonomic sequence database"
    )
    build_parser.add_argument(
        "-i", "--input", required=True, help="Input file for building"
    )
    build_parser.add_argument(
        "-o", "--output", default="ChimeraDB", help="Output file for building"
    )
    build_parser.add_argument(
        "--strobe-k",
        type=int,
        default=28,
        help="Strobemer k-mer length",
    )
    build_parser.add_argument(
        "--strobe-order",
        type=int,
        default=2,
        help="Strobemer order (currently only 2 is supported)",
    )
    build_parser.add_argument(
        "--strobe-w-min",
        type=int,
        default=12,
        help="Strobemer minimum window",
    )
    build_parser.add_argument(
        "--strobe-w-max",
        type=int,
        default=32,
        help="Strobemer maximum window",
    )
    build_parser.add_argument(
        "-l",
        "--min-length",
        type=min_length_type,
        default="auto",
        help="Minimum sequence length for building (auto => strobemer minimum span)",
    )
    build_parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=default_threads(),
        help="Number of threads for building",
    )
    build_parser.add_argument(
        "--load-factor", type=float, default=0.85, help="IMCF filter load factor"
    )
    build_parser.add_argument(
        "--presence-unique-deg",
        type=int,
        default=1,
        help="Degree cutoff (<=) treated as unique signature for coverage meta",
    )
    build_parser.add_argument(
        "--taxonomy-kind",
        dest="taxonomy_kind",
        default="auto",
        choices=["auto", "ncbi", "gtdb"],
        help="Taxonomy source identifier (auto/ncbi/gtdb)",
    )
    build_parser.add_argument(
        "--taxonomy-version",
        dest="taxonomy_version",
        default="auto",
        help="Taxonomy version label, for example ncbi-taxdump-2025-09-15 or gtdb-rs226",
    )

    # Download and Build combined subcommand
    download_build_parser = subparsers.add_parser(
        "download_and_build",
        help="Download NCBI database sequences and resources and build a taxonomic sequence database",
    )
    download_build_parser.add_argument(
        "-o", "--output", default="ChimeraDB", help="Output file for building"
    )
    download_build_parser.add_argument(
        "--strobe-k",
        type=int,
        default=28,
        help="Strobemer k-mer length",
    )
    download_build_parser.add_argument(
        "--strobe-order",
        type=int,
        default=2,
        help="Strobemer order (currently only 2 is supported)",
    )
    download_build_parser.add_argument(
        "--strobe-w-min",
        type=int,
        default=12,
        help="Strobemer minimum window",
    )
    download_build_parser.add_argument(
        "--strobe-w-max",
        type=int,
        default=32,
        help="Strobemer maximum window",
    )
    download_build_parser.add_argument(
        "-l",
        "--min-length",
        type=min_length_type,
        default="auto",
        help="Minimum sequence length for building (auto => strobemer minimum span)",
    )
    download_build_parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=default_threads(),
        help="Number of threads for building",
    )
    download_build_parser.add_argument(
        "--load-factor", type=float, default=0.85, help="IMCF filter load factor"
    )
    download_build_parser.add_argument(
        "--presence-unique-deg",
        type=int,
        default=1,
        help="Degree cutoff (<=) treated as unique signature for coverage meta",
    )
    download_build_parser.add_argument(
        "--taxonomy-kind",
        dest="taxonomy_kind",
        default="auto",
        choices=["auto", "ncbi", "gtdb"],
        help="Taxonomy source identifier (auto/ncbi/gtdb)",
    )
    download_build_parser.add_argument(
        "--taxonomy-version",
        dest="taxonomy_version",
        default="auto",
        help="Taxonomy version label, for example ncbi-taxdump-2025-09-15 or gtdb-rs226",
    )

    # Classify subcommand
    classify_parser = subparsers.add_parser("classify", help="Classify sequences")
    classify_parser.add_argument(
        "-i",
        "--single",
        nargs="+",
        help="Input files for classifying (supports multiple files)",
    )
    classify_parser.add_argument(
        "-p", "--paired", nargs="+", help="Paired input files for classifying"
    )
    classify_parser.add_argument(
        "-o",
        "--output",
        default="ChimeraOutput",
        help="Output directory for classify, evidence, and profile results",
    )
    classify_parser.add_argument(
        "-d", "--database", required=True, help="Database file for classifying"
    )
    classify_parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=default_threads(),
        help="Number of threads for classifying",
    )
    classify_parser.add_argument(
        "-b", "--batch-size", type=int, default=400, help="Batch size for classifying"
    )
    classify_parser.add_argument(
        "--read-evidence",
        action="store_true",
        help="Write ChimeraReadEvidence.tsv for read-resolved audits",
    )
    # Profile subcommand
    profile_parser = subparsers.add_parser("profile", help="Generate sequence profile")
    profile_parser.add_argument(
        "-i", "--input", nargs="+", required=True, help="Input file(s) for profiling"
    )
    profile_parser.add_argument(
        "-o", "--output", default="ChimeraProfile", help="Output file for profiling"
    )
    profile_parser.add_argument(
        "-k", "--krona", action="store_true", help="Generate Krona chart"
    )
    profile_parser.add_argument(
        "--taxonomy-kind",
        dest="taxonomy_kind",
        default="auto",
        choices=["auto", "ncbi", "gtdb"],
        help="Taxonomy source for profile generation (auto attempts inference)",
    )
    profile_parser.add_argument(
        "--taxonomy-info",
        dest="taxonomy_info",
        default=None,
        help="Path to the tax.info file produced during build, used by GTDB and other custom taxonomies",
    )
    profile_parser.add_argument(
        "--taxonomy-meta",
        dest="taxonomy_meta",
        default=None,
        help="Path to taxonomy.meta metadata; if omitted, the input directory will be probed",
    )

    extract_parser = subparsers.add_parser(
        "extract-reads", help="Extract taxon-specific read bags from ChimeraReadEvidence.tsv"
    )
    extract_parser.add_argument(
        "--ledger", required=True, help="Path to ChimeraReadEvidence.tsv"
    )
    extract_parser.add_argument(
        "--reads",
        nargs="+",
        required=True,
        help="Input FASTQ file, or paired FASTQ files",
    )
    extract_parser.add_argument(
        "--taxid",
        action="append",
        default=[],
        help="Target taxid to extract; repeat for multiple taxa",
    )
    extract_parser.add_argument(
        "--taxids",
        default=None,
        help="File with one target taxid per line",
    )
    extract_parser.add_argument(
        "--min-profile-mass",
        type=float,
        default=0.0,
        help="Minimum per-read profile mass for the target taxid",
    )
    extract_parser.add_argument(
        "--out", required=True, help="Output directory for extracted read bags"
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    return args


def run_chimera(args, chimera_path=None):
    if args.command == "download":
        downloader = get_downloader()
        cmd_args = sys.argv
        try:
            download_index = cmd_args.index("download")
            download_args = cmd_args[download_index + 1 :]

            if download_args:
                downloader.download(interactive=False, raw_args=download_args)
            else:
                downloader.download(interactive=True)
        except ValueError:
            downloader.download(interactive=True)
        return 0

    if args.command == "download_and_build":
        downloader = get_downloader()
        options = downloader.download(interactive=True)
        args.command = "build"
        args.input = os.path.join(options.output_dir, "target.tsv")
        resolved_kind = getattr(options, "taxonomy_kind_resolved", None)
        if resolved_kind is None:
            resolved_kind = getattr(options, "taxonomy_mode", None)
        if resolved_kind is None:
            resolved_kind = "auto"
        args.taxonomy_kind = resolved_kind
        resolved_version = getattr(options, "taxonomy_version_resolved", None)
        if resolved_version is None:
            resolved_version = getattr(options, "version_label", None)
        if not resolved_version:
            resolved_version = "auto"
        args.taxonomy_version = resolved_version

    if args.command == "profile":
        if args.taxonomy_meta is None and args.input:
            candidate_meta = Path(args.input[0]).resolve().parent / "taxonomy.meta"
            if candidate_meta.exists():
                args.taxonomy_meta = str(candidate_meta)
        if args.taxonomy_info is None:
            if args.taxonomy_meta:
                meta_dir = Path(args.taxonomy_meta).resolve().parent
                candidate_info = meta_dir / "tax.info"
                if candidate_info.exists():
                    args.taxonomy_info = str(candidate_info)
            if args.taxonomy_info is None and args.input:
                candidate_info = Path(args.input[0]).resolve().parent / "tax.info"
                if candidate_info.exists():
                    args.taxonomy_info = str(candidate_info)
        if isinstance(args.output, str) and args.output.lower().endswith(".tsv"):
            args.output = args.output[: -len(".tsv")]
        if args.krona:
            downloader = get_downloader()
            from src.profile import conversion2Krona
            print("Converting file to Krona format...")
            conversion2Krona.convert_multiple_files_to_krona_format(
                args.input,
                args.output,
                taxonomy_kind=args.taxonomy_kind,
                taxonomy_info=args.taxonomy_info,
                taxonomy_meta=args.taxonomy_meta,
            )
            print("Conversion completed.")
            print("Generating Krona chart...")
            downloader.run_command(
                ["ktImportText", args.output + ".tsv", "-o", args.output + ".html"]
            )
            print("Krona chart generated.")
        try:
            profile.process_file(
                args.input,
                args.output,
                taxonomy_kind=args.taxonomy_kind,
                taxonomy_info=args.taxonomy_info,
                taxonomy_meta=args.taxonomy_meta,
            )
        except ValueError as exc:
            print(f"Profile failed: {exc}")
            return 1
        return 0

    if args.command == "extract-reads":
        from src.profile.read_evidence import (
            ReadEvidenceError,
            extract_read_bags,
            load_target_taxids,
        )

        try:
            target_taxids = load_target_taxids(args.taxid, args.taxids)
            counts = extract_read_bags(
                ledger_path=args.ledger,
                read_paths=args.reads,
                target_taxids=target_taxids,
                output_dir=args.out,
                min_profile_mass=args.min_profile_mass,
            )
        except ReadEvidenceError as exc:
            print(f"Extract reads failed: {exc}")
            return 1
        print(f"Extracted read bags into: {args.out}")
        for taxid in sorted(target_taxids, key=lambda x: (0, int(x)) if x.isdigit() else (1, x)):
            print(f"{taxid}\t{counts.get(taxid, 0)}")
        return 0

    if chimera_path is None:
        chimera_path = get_chimera_path()
    command = [chimera_path]

    if args.version:
        command.append("--version")

    if args.command:
        command.append(args.command)

    if args.command == "build":
        command.extend(["-i", args.input])
        command.extend(["-o", args.output])
        command.extend(["--strobe-k", str(args.strobe_k)])
        command.extend(["--strobe-order", str(args.strobe_order)])
        command.extend(["--strobe-w-min", str(args.strobe_w_min)])
        command.extend(["--strobe-w-max", str(args.strobe_w_max)])
        if args.min_length != "auto":
            command.extend(["-l", str(args.min_length)])
        command.extend(["-t", str(args.threads)])
        command.extend(["--load-factor", str(args.load_factor)])
        command.extend(
            ["--presence-unique-deg", str(args.presence_unique_deg)]
        )
        command.extend(["--taxonomy-kind", str(args.taxonomy_kind)])
        command.extend(["--taxonomy-version", str(args.taxonomy_version)])

    elif args.command == "classify":
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        classify_output = output_dir / "ChimeraClassify.tsv"
        evidence_output = output_dir / "ChimeraEvidence.tsv"
        profile_output = output_dir / "ChimeraProfile"
        if args.single:
            command.extend(["-i"] + args.single)
        if args.paired:
            command.extend(["-p"] + args.paired)
        command.extend(["-o", str(classify_output)])
        command.extend(["-d", args.database])
        command.extend(["-t", str(args.threads)])
        command.extend(["-b", str(args.batch_size)])
        if args.read_evidence:
            command.append("--read-evidence")

    if args.command == "classify":
        subprocess.run(command, check=True)
        if not evidence_output.is_file():
            print(f"Classify did not produce abundance evidence: {evidence_output}")
            return 1
        taxonomy_meta, taxonomy_info = _infer_profile_taxonomy_paths(args.database)
        try:
            profile.process_file(
                [str(evidence_output)],
                str(profile_output),
                taxonomy_kind="auto",
                taxonomy_info=taxonomy_info,
                taxonomy_meta=taxonomy_meta,
            )
        except ValueError as exc:
            print(f"Profile failed: {exc}")
            return 1
        return 0
    downloader = get_downloader()
    downloader.run_command(command)
    return 0


def main():
    args = parse_arguments()
    raise SystemExit(run_chimera(args) or 0)


if __name__ == "__main__":
    main()
