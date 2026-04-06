import argparse
import os
import shutil
import sys
from pathlib import Path

from src.profile import profile


def default_threads():
    cpu = os.cpu_count()
    return cpu if cpu and cpu > 0 else 1


def get_downloader():
    from src.download import download as downloader
    return downloader


def get_chimera_path():
    # Use which command to find the chimera program path
    chimera_path = shutil.which("Chimera")
    if chimera_path is None:
        raise FileNotFoundError(
            "Cannot find 'Chimera' executable. Please ensure it is installed and in your PATH."
        )
    return chimera_path


def kmer_type(x):
    x = int(x)
    if x < 1 or x > 50:
        raise argparse.ArgumentTypeError("Kmer size must be between 1 and 50")
    return x


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
        "-k",
        "--kmer",
        type=kmer_type,
        default=31,
        help="Kmer size for building (must be between 1 and 50)",
    )
    build_parser.add_argument(
        "-s",
        "--syncmer-s",
        type=int,
        default=16,
        help="Syncmer s-mer size (must be positive and smaller than k)",
    )
    build_parser.add_argument(
        "-P",
        "--syncmer-pos",
        type=int,
        default=7,
        help="Syncmer minimal s-mer offset (0-based)",
    )
    build_parser.add_argument(
        "-l",
        "--min-length",
        type=min_length_type,
        default="auto",
        help="Minimum sequence length for building (auto => k-mer size)",
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
        "-k",
        "--kmer",
        type=kmer_type,
        default=31,
        help="Kmer size for building (must be between 1 and 50)",
    )
    download_build_parser.add_argument(
        "-s",
        "--syncmer-s",
        type=int,
        default=16,
        help="Syncmer s-mer size (must be positive and smaller than k)",
    )
    download_build_parser.add_argument(
        "-P",
        "--syncmer-pos",
        type=int,
        default=7,
        help="Syncmer minimal s-mer offset (0-based)",
    )
    download_build_parser.add_argument(
        "-l",
        "--min-length",
        type=min_length_type,
        default="auto",
        help="Minimum sequence length for building (auto => k-mer size)",
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
        "-o", "--output", default="ChimeraClassify", help="Output file for classifying"
    )
    classify_parser.add_argument(
        "-d", "--database", required=True, help="Database file for classifying"
    )
    classify_parser.add_argument(
        "-s",
        "--shot-threshold",
        type=float,
        default=None,
        help="Shot threshold for classifying (defaults to 0.70 if unset)",
    )
    classify_parser.add_argument(
        "--first-filter-beta",
        type=float,
        default=None,
        help="Beta multiplier for first-stage filter (default 0.75)",
    )
    classify_parser.add_argument(
        "--presence-pi",
        type=float,
        default=None,
        help="Presence prior probability for coverage model (0-1)",
    )
    classify_parser.add_argument(
        "--presence-tau",
        type=float,
        default=None,
        help="Log posterior odds threshold for coverage model",
    )
    classify_parser.add_argument(
        "--presence-breadth-bits",
        type=int,
        default=None,
        help="Sketch bits for presence breadth (power-of-two suggested)",
    )
    classify_parser.add_argument(
        "--post-pi-min",
        type=float,
        default=None,
        help="Minimum global class weight (default 0.0005)",
    )
    # NOTE: Deprecated LCA/post-processing knobs are intentionally not exposed.
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
    profile_parser.add_argument(
        "--abundance-mode",
        dest="abundance_mode",
        default="soft_seq",
        choices=["soft_seq", "hard_seq", "top1", "soft_tax", "soft", "hard"],
        help="Abundance mode: soft_seq (default, sum all taxid:count), hard_seq (top1 taxid:count only), or top1 (count each sequence once)",
    )
    profile_parser.add_argument(
        "--hd-species-head-mass",
        type=float,
        default=97.0,
        help="Head-mass percentage for species output in high-diversity samples (default 97.0)",
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if args.command in {"build", "download_and_build"}:
        if args.syncmer_s < 1:
            parser.error("--syncmer-s must be greater than 0")
        if args.syncmer_s >= args.kmer:
            parser.error("--syncmer-s must be smaller than k-mer size")
        window_span = args.kmer - args.syncmer_s + 1
        if args.syncmer_pos < 0 or args.syncmer_pos >= window_span:
            parser.error("--syncmer-pos must satisfy 0 <= pos < k - s + 1")

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
                abundance_mode=args.abundance_mode,
                hd_species_head_mass=args.hd_species_head_mass,
            )
        except ValueError as exc:
            print(f"Profile failed: {exc}")
            return 1
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
        command.extend(["-k", str(args.kmer)])
        command.extend(["-s", str(args.syncmer_s)])
        command.extend(["-P", str(args.syncmer_pos)])
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
        if args.single:
            command.extend(["-i"] + args.single)
        if args.paired:
            command.extend(["-p"] + args.paired)
        command.extend(["-o", args.output])
        command.extend(["-d", args.database])
        if args.shot_threshold is not None:
            command.extend(["-s", str(args.shot_threshold)])
        if args.first_filter_beta is not None:
            command.extend(["--first-filter-beta", str(args.first_filter_beta)])
        if args.presence_pi is not None:
            command.extend(["--presence-pi", str(args.presence_pi)])
        if args.presence_tau is not None:
            command.extend(["--presence-tau", str(args.presence_tau)])
        if args.presence_breadth_bits is not None:
            command.extend(
                ["--presence-breadth-bits", str(args.presence_breadth_bits)]
            )
        if args.post_pi_min is not None:
            command.extend(["--post-pi-min", str(args.post_pi_min)])
        command.extend(["-t", str(args.threads)])
        command.extend(["-b", str(args.batch_size)])

    # Execute the command using the provided run function
    downloader = get_downloader()
    downloader.run_command(command)


def main():
    args = parse_arguments()
    run_chimera(args)


if __name__ == "__main__":
    main()
