import argparse
import os
import shutil
import sys

from src.download import download as downloader
from src.profile import conversion2Krona, profile


def default_threads():
    cpu = os.cpu_count()
    return cpu if cpu and cpu > 0 else 1


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
    if x < 1 or x > 31:
        raise argparse.ArgumentTypeError("Kmer size must be between 1 and 31")
    return x


def parse_arguments():
    if len(sys.argv) > 1 and sys.argv[1] == "download":
        parser = argparse.ArgumentParser(
            description="Chimera - A versatile tool for metagenomic classification"
        )
        subparsers = parser.add_subparsers(dest="command", required=False)
        download_parser = subparsers.add_parser(
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
    download_parser = subparsers.add_parser(
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
        "-m",
        "--mode",
        default="normal",
        choices=["fast", "normal"],
        help="Mode for building (choices: 'fast', 'normal'). Default is 'normal'.",
    )
    build_parser.add_argument(
        "-k",
        "--kmer",
        type=kmer_type,
        default=19,
        help="Kmer size for building (must be between 1 and 31)",
    )
    build_parser.add_argument(
        "-w", "--window", type=int, default=31, help="Window size for building"
    )
    build_parser.add_argument(
        "-l",
        "--min-length",
        type=int,
        default=0,
        help="Minimum length sequence for building",
    )
    build_parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=default_threads(),
        help="Number of threads for building",
    )
    build_parser.add_argument(
        "--load-factor", type=float, default=0.9, help="Loading ratio of CF"
    )
    build_parser.add_argument(
        "-M",
        "--max-hashes",
        type=int,
        default=2000000,
        help="Maximum number of hashes per taxid",
    )
    build_parser.add_argument(
        "-a", "--alpha", type=float, default=1.2, help="Alpha value for HICF"
    )
    build_parser.add_argument(
        "--relaxed-load-factor",
        type=float,
        default=0.95,
        help="Relaxed loading ratio of HICF",
    )
    build_parser.add_argument(
        "-f",
        "--filter",
        default="imcf",
        choices=["icf", "hicf", "imcf"],
        help="Filter for building (choices: 'icf', 'hicf, 'imcf'). Default is 'imcf'.",
    )
    build_parser.add_argument(
        "-c", "--fixed-cutoff", type=int, help="Fixed cutoff for minimizer (0 - 255)"
    )
    build_parser.add_argument(
        "-q", "--quiet", action="store_true", default=False, help="Quiet output"
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
        "-m",
        "--mode",
        default="normal",
        choices=["fast", "normal"],
        help="Mode for building (choices: 'fast', 'normal'). Default is 'normal'.",
    )
    download_build_parser.add_argument(
        "-k",
        "--kmer",
        type=kmer_type,
        default=19,
        help="Kmer size for building (must be between 1 and 31)",
    )
    download_build_parser.add_argument(
        "-w", "--window", type=int, default=31, help="Window size for building"
    )
    download_build_parser.add_argument(
        "-l",
        "--min-length",
        type=int,
        default=0,
        help="Minimum length sequence for building",
    )
    download_build_parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=default_threads(),
        help="Number of threads for building",
    )
    download_build_parser.add_argument(
        "--load-factor", type=float, default=0.58, help="Loading ratio of CF"
    )
    download_build_parser.add_argument(
        "-M",
        "--max-hashes",
        type=int,
        default=2000000,
        help="Maximum number of hashes per taxid",
    )
    download_build_parser.add_argument(
        "-a", "--alpha", type=float, default=1.2, help="Alpha value for HICF"
    )
    download_build_parser.add_argument(
        "--relaxed-load-factor",
        type=float,
        default=0.95,
        help="Relaxed loading ratio of HICF",
    )
    download_build_parser.add_argument(
        "-f",
        "--filter",
        default="imcf",
        choices=["icf", "hicf", "imcf"],
        help="Filter for building (choices: 'icf', 'hicf, 'imcf'). Default is 'imcf'.",
    )
    download_build_parser.add_argument(
        "-c", "--fixed-cutoff", type=int, help="Fixed cutoff for minimizer (0 - 255)"
    )
    download_build_parser.add_argument(
        "-q", "--quiet", action="store_true", default=False, help="Quiet output"
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
        help="Shot threshold for classifying (defaults to 0.65 if unset)",
    )
    classify_parser.add_argument(
        "--adaptive-shot",
        dest="adaptive_shot",
        action="store_true",
        default=None,
        help="Scale thresholds by evaluated minimizers (default on)",
    )
    classify_parser.add_argument(
        "--no-adaptive-shot",
        dest="adaptive_shot",
        action="store_false",
        help="Disable threshold scaling by evaluated minimizers",
    )
    classify_parser.add_argument(
        "--first-filter-beta",
        type=float,
        default=None,
        help="Beta multiplier for first-stage filter (default 0.75)",
    )
    classify_parser.add_argument(
        "--pre-em-topk",
        type=int,
        default=None,
        help="Limit candidates before EM/VEM to top K (default 24)",
    )
    classify_parser.add_argument(
        "--adaptive-fdr",
        dest="adaptive_fdr",
        action="store_true",
        default=None,
        help="Enable per-read noise-calibrated threshold (default on)",
    )
    classify_parser.add_argument(
        "--no-adaptive-fdr",
        dest="adaptive_fdr",
        action="store_false",
        help="Disable per-read noise-calibrated threshold",
    )
    classify_parser.add_argument(
        "--fdr-z",
        type=float,
        default=None,
        help="Z score for Poisson(mu)+Z*sqrt(mu) threshold (default 3.0)",
    )
    classify_parser.add_argument(
        "--min-eval-count",
        type=int,
        default=None,
        help="Minimum evaluated minimizers required to classify (default 24)",
    )
    classify_parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=default_threads(),
        help="Number of threads for classifying",
    )
    classify_parser.add_argument(
        "-m",
        "--mode",
        default="normal",
        choices=["fast", "normal"],
        help="Mode for classifying (choices: 'fast', 'normal'). Default is 'normal'.",
    )
    classify_parser.add_argument(
        "-f",
        "--filter",
        default="imcf",
        choices=["icf", "hicf", "imcf"],
        help="Filter for classifying (choices: 'icf', 'hicf', 'imcf'). Default is 'imcf'.",
    )
    classify_parser.add_argument(
        "-b", "--batch-size", type=int, default=400, help="Batch size for classifying"
    )
    group = classify_parser.add_mutually_exclusive_group()

    group.add_argument(
        "-l", "--lca", action="store_true", help="Use LCA algorithm for classification"
    )
    classify_parser.add_argument(
        "-T", "--tax-file", help="Taxonomy file for LCA classification"
    )

    group.add_argument(
        "-e",
        "--em",
        action="store_true",
        help="Use EM algorithm for classification (default)",
    )
    group.add_argument(
        "-V", "--vem", action="store_true", help="Use VEM algorithm for classification"
    )
    classify_parser.add_argument(
        "--em-iter", type=int, default=80, help="Number of EM iterations"
    )
    classify_parser.add_argument(
        "--em-threshold", type=float, default=0.001, help="EM threshold"
    )

    group.add_argument(
        "--none", action="store_true", help="Do not use LCA or EM for classification"
    )

    classify_parser.add_argument(
        "-q", "--quiet", action="store_true", default=False, help="Quiet output"
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

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if args.command == "classify":
        algorithms = [
            args.em,
            args.vem,
            args.lca,
            getattr(args, "none", False),
        ]
        if not any(algorithms):
            args.em = True

    return args


def run_chimera(args, chimera_path):
    if args.command == "download":
        # 获取原始命令行参数
        cmd_args = sys.argv
        # 找到"download"在参数列表中的位置
        try:
            download_index = cmd_args.index("download")
            # 提取"download"之后的所有参数
            download_args = cmd_args[download_index + 1 :]

            # 如果有额外参数，则使用非交互模式
            if download_args:
                downloader.download(interactive=False, raw_args=download_args)
            else:
                # 无参数时使用交互模式
                downloader.download(interactive=True)
        except ValueError:
            # 找不到"download"参数，使用交互模式
            downloader.download(interactive=True)
        return 0

    if args.command == "download_and_build":
        options = downloader.download(interactive=True)
        args.command = "build"
        args.input = os.path.join(options.output_dir, "target.tsv")

    if args.command == "profile":
        if args.krona:
            print("Converting file to Krona format...")
            conversion2Krona.convert_multiple_files_to_krona_format(
                args.input, args.output
            )
            print("Conversion completed.")
            print("Generating Krona chart...")
            downloader.run_command(
                ["ktImportText", args.output + ".tsv", "-o", args.output + ".html"]
            )
            print("Krona chart generated.")
        profile.process_file(args.input, args.output)
        return 0

    command = [chimera_path]

    if args.version:
        command.append("--version")

    # Add subcommand
    if args.command:
        command.append(args.command)

    # Add options and parameters based on parsed arguments
    if args.command == "download":
        if args.build:
            command.append("-b")

    if args.command == "build":
        command.extend(["-i", args.input])
        command.extend(["-o", args.output])
        command.extend(["-m", args.mode])
        command.extend(["-k", str(args.kmer)])
        command.extend(["-w", str(args.window)])
        command.extend(["-l", str(args.min_length)])
        command.extend(["-t", str(args.threads)])
        command.extend(["--load-factor", str(args.load_factor)])
        command.extend(["-a", str(args.alpha)])
        command.extend(["--relaxed-load-factor", str(args.relaxed_load_factor)])
        command.extend(["-f", args.filter])
        if args.fixed_cutoff:
            command.extend(["-c", str(args.fixed_cutoff)])
        command.extend(["-M", str(args.max_hashes)])
        if args.quiet:
            command.append("-q")

    elif args.command == "classify":
        if args.single:
            command.extend(["-i"] + args.single)
        if args.paired:
            command.extend(["-p"] + args.paired)
        command.extend(["-o", args.output])
        command.extend(["-d", args.database])
        if args.shot_threshold is not None:
            command.extend(["-s", str(args.shot_threshold)])
        if args.adaptive_shot is True:
            command.append("--adaptive-shot")
        elif args.adaptive_shot is False:
            command.append("--no-adaptive-shot")
        if args.first_filter_beta is not None:
            command.extend(["--first-filter-beta", str(args.first_filter_beta)])
        if args.pre_em_topk is not None:
            command.extend(["--pre-em-topk", str(args.pre_em_topk)])
        if args.adaptive_fdr is True:
            command.append("--adaptive-fdr")
        elif args.adaptive_fdr is False:
            command.append("--no-adaptive-fdr")
        if args.fdr_z is not None:
            command.extend(["--fdr-z", str(args.fdr_z)])
        if args.min_eval_count is not None:
            command.extend(["--min-eval-count", str(args.min_eval_count)])
        command.extend(["-t", str(args.threads)])
        command.extend(["-m", args.mode])
        command.extend(["-b", str(args.batch_size)])
        command.extend(["-f", args.filter])
        if args.lca:
            command.append("--lca")
            if not args.tax_file:
                raise ValueError(
                    "Taxonomy file must be provided when using LCA algorithm"
                )
            command.extend(["--tax-file", args.tax_file])
        if args.em:
            command.append("-e")
            command.extend(["--em-iter", str(args.em_iter)])
            command.extend(["--em-threshold", str(args.em_threshold)])
        if args.vem:
            command.append("-V")
            command.extend(["--em-iter", str(args.em_iter)])
            command.extend(["--em-threshold", str(args.em_threshold)])
        if args.quiet:
            command.append("-q")

    # Execute the command using the provided run function
    downloader.run_command(command)


def main():
    chimera_path = get_chimera_path()
    args = parse_arguments()
    run_chimera(args, chimera_path)


if __name__ == "__main__":
    main()
