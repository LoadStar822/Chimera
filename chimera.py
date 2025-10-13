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
    if x < 1 or x > 50:
        raise argparse.ArgumentTypeError("Kmer size must be between 1 and 50")
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
        "--load-factor", type=float, default=0.85, help="IMCF 滤器的负载因子"
    )
    build_parser.add_argument(
        "--toxic-quantile",
        type=float,
        default=0.999,
        help="过滤高频 syncmer 的分位阈值 (0-1]",
    )
    build_parser.add_argument(
        "--toxic-topN",
        dest="toxic_topn",
        type=int,
        default=0,
        help="优先剔除全库计数前 N 的 hash",
    )
    build_parser.add_argument(
        "--toxic-min-fraction",
        type=float,
        default=0.0001,
        help="毒性候选的最低全库占比",
    )
    build_parser.add_argument(
        "--toxic-safety-frac",
        type=float,
        default=0.1,
        help="每个 taxid 最少保留比例",
    )
    build_parser.add_argument(
        "--toxic-safety-min",
        type=int,
        default=1024,
        help="每个 taxid 最少保留数量",
    )
    build_parser.add_argument(
        "-M",
        "--max-hashes",
        type=int,
        default=2000000,
        help="Maximum number of hashes per taxid",
    )
    build_parser.add_argument(
        "-f",
        "--filter",
        default="imcf",
        choices=["imcf"],
        help="构建使用的滤器 (仅 imcf)",
    )
    build_parser.add_argument(
        "--adaptive-cutoff",
        action="store_true",
        help="启用基于文件规模的自适应 cutoff",
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
        "--load-factor", type=float, default=0.85, help="IMCF 滤器的负载因子"
    )
    download_build_parser.add_argument(
        "--toxic-quantile",
        type=float,
        default=0.999,
        help="过滤高频 syncmer 的分位阈值 (0-1]",
    )
    download_build_parser.add_argument(
        "--toxic-topN",
        dest="toxic_topn",
        type=int,
        default=0,
        help="优先剔除全库计数前 N 的 hash",
    )
    download_build_parser.add_argument(
        "--toxic-min-fraction",
        type=float,
        default=0.0001,
        help="毒性候选的最低全库占比",
    )
    download_build_parser.add_argument(
        "--toxic-safety-frac",
        type=float,
        default=0.1,
        help="每个 taxid 最少保留比例",
    )
    download_build_parser.add_argument(
        "--toxic-safety-min",
        type=int,
        default=1024,
        help="每个 taxid 最少保留数量",
    )
    download_build_parser.add_argument(
        "-M",
        "--max-hashes",
        type=int,
        default=2000000,
        help="Maximum number of hashes per taxid",
    )
    download_build_parser.add_argument(
        "-f",
        "--filter",
        default="imcf",
        choices=["imcf"],
        help="构建使用的滤器 (仅 imcf)",
    )
    download_build_parser.add_argument(
        "--adaptive-cutoff",
        action="store_true",
        help="启用基于文件规模的自适应 cutoff",
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
        help="Shot threshold for classifying (defaults to 0.62 if unset)",
    )
    classify_parser.add_argument(
        "--adaptive-shot",
        dest="adaptive_shot",
        action="store_true",
        default=None,
        help="Scale thresholds by evaluated syncmers (default on)",
    )
    classify_parser.add_argument(
        "--no-adaptive-shot",
        dest="adaptive_shot",
        action="store_false",
        help="Disable threshold scaling by evaluated syncmers",
    )
    classify_parser.add_argument(
        "--first-filter-beta",
        type=float,
        default=None,
        help="Beta multiplier for first-stage filter (default 0.80)",
    )
    classify_parser.add_argument(
        "--pre-em-topk",
        type=int,
        default=None,
        help="Limit candidates before EM/VEM to top K (default 32)",
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
        help="Minimum evaluated syncmers required to classify (default auto)",
    )
    classify_parser.add_argument(
        "--post-thres",
        type=float,
        default=None,
        help="Posterior acceptance threshold (default 0.56)",
    )
    classify_parser.add_argument(
        "--post-margin",
        type=float,
        default=None,
        help="Minimum gap between top posteriors (default 0.03)",
    )
    classify_parser.add_argument(
        "--post-ratio",
        type=float,
        default=None,
        help="Minimum ratio between top1 and top2 posteriors (default 1.30)",
    )
    classify_parser.add_argument(
        "--post-pi-min",
        type=float,
        default=None,
        help="Minimum global class weight (default 0.0001)",
    )
    classify_parser.add_argument(
        "--evidence-override",
        dest="evidence_override",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Allow strong pre-EM evidence to bypass posterior filter",
    )
    classify_parser.add_argument(
        "--output-posterior",
        dest="output_posterior",
        action="store_true",
        default=None,
        help="Write posterior probabilities to the TSV output",
    )
    classify_parser.add_argument(
        "--no-output-posterior",
        dest="output_posterior",
        action="store_false",
        help="Do not write posterior probabilities to the TSV output",
    )
    classify_parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=default_threads(),
        help="Number of threads for classifying",
    )
    classify_parser.add_argument(
        "-f",
        "--filter",
        default="imcf",
        choices=["imcf"],
        help="分类使用的滤器 (仅 imcf)",
    )
    classify_parser.add_argument(
        "-b", "--batch-size", type=int, default=400, help="Batch size for classifying"
    )
    classify_parser.add_argument(
        "-e",
        "--em",
        action="store_true",
        help="Use EM algorithm for classification (default)",
    )
    classify_parser.add_argument(
        "-V", "--vem", action="store_true", help="Use VEM algorithm for classification"
    )
    classify_parser.add_argument(
        "--em-iter", type=int, default=80, help="Number of EM iterations"
    )
    classify_parser.add_argument(
        "--em-threshold", type=float, default=0.001, help="EM threshold"
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
        if not any([args.em, args.vem]):
            args.em = True

    if args.command in {"build", "download_and_build"}:
        if args.syncmer_s < 1:
            parser.error("--syncmer-s must be greater than 0")
        if args.syncmer_s >= args.kmer:
            parser.error("--syncmer-s must be smaller than k-mer size")
        window_span = args.kmer - args.syncmer_s + 1
        if args.syncmer_pos < 0 or args.syncmer_pos >= window_span:
            parser.error("--syncmer-pos must satisfy 0 <= pos < k - s + 1")
        if args.toxic_topn < 0:
            parser.error("--toxic-topN 必须为非负整数")
        if args.toxic_topn == 0:
            if not (0 < args.toxic_quantile <= 1.0):
                parser.error(
                    "当未指定 --toxic-topN 时，--toxic-quantile 必须落在 (0, 1] 区间"
                )
        else:
            if not (0 <= args.toxic_quantile <= 1.0):
                parser.error("--toxic-quantile 必须落在 [0, 1] 区间")
            if args.toxic_quantile == 0:
                args.toxic_quantile = 1.0
        for value, name in [
            (args.toxic_min_fraction, "--toxic-min-fraction"),
            (args.toxic_safety_frac, "--toxic-safety-frac"),
        ]:
            if not (0.0 <= value <= 1.0):
                parser.error(f"{name} 必须落在 [0, 1] 区间")
        if args.toxic_safety_min < 0:
            parser.error("--toxic-safety-min 必须为非负整数")

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
        command.extend(["-k", str(args.kmer)])
        command.extend(["-s", str(args.syncmer_s)])
        command.extend(["-P", str(args.syncmer_pos)])
        command.extend(["-l", str(args.min_length)])
        command.extend(["-t", str(args.threads)])
        command.extend(["--load-factor", str(args.load_factor)])
        command.extend(["--toxic-quantile", str(args.toxic_quantile)])
        command.extend(["--toxic-topN", str(args.toxic_topn)])
        command.extend(["--toxic-min-fraction", str(args.toxic_min_fraction)])
        command.extend(["--toxic-safety-frac", str(args.toxic_safety_frac)])
        command.extend(["--toxic-safety-min", str(args.toxic_safety_min)])
        command.extend(["-M", str(args.max_hashes)])
        command.extend(["-f", args.filter])
        if args.adaptive_cutoff:
            command.append("--adaptive-cutoff")
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
        if args.post_thres is not None:
            command.extend(["--post-thres", str(args.post_thres)])
        if args.post_margin is not None:
            command.extend(["--post-margin", str(args.post_margin)])
        if args.post_ratio is not None:
            command.extend(["--post-ratio", str(args.post_ratio)])
        if args.post_pi_min is not None:
            command.extend(["--post-pi-min", str(args.post_pi_min)])
        if args.evidence_override is True:
            command.append("--evidence-override")
        elif args.evidence_override is False:
            command.append("--no-evidence-override")
        if args.output_posterior is True:
            command.append("--output-posterior")
        elif args.output_posterior is False:
            command.append("--no-output-posterior")
        command.extend(["-t", str(args.threads)])
        command.extend(["-b", str(args.batch_size)])
        command.extend(["-f", args.filter])
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
