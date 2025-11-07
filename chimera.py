import argparse
import os
import shutil
import sys
from pathlib import Path

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
        "--load-factor", type=float, default=0.9, help="IMCF 滤器的负载因子"
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
        "--taxonomy-kind",
        dest="taxonomy_kind",
        default="auto",
        choices=["auto", "ncbi", "gtdb"],
        help="taxonomy 数据源标识 (auto/ncbi/gtdb)",
    )
    build_parser.add_argument(
        "--taxonomy-version",
        dest="taxonomy_version",
        default="auto",
        help="taxonomy 数据版本标识，例如 ncbi-taxdump-2025-09-15 或 gtdb-rs226",
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
        "--load-factor", type=float, default=0.9, help="IMCF 滤器的负载因子"
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
        "--taxonomy-kind",
        dest="taxonomy_kind",
        default="auto",
        choices=["auto", "ncbi", "gtdb"],
        help="taxonomy 数据源标识 (auto/ncbi/gtdb)",
    )
    download_build_parser.add_argument(
        "--taxonomy-version",
        dest="taxonomy_version",
        default="auto",
        help="taxonomy 数据版本标识，例如 ncbi-taxdump-2025-09-15 或 gtdb-rs226",
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
        help="Beta multiplier for first-stage filter (default 0.75)",
    )
    classify_parser.add_argument(
        "--pre-em-topk",
        type=int,
        default=None,
        help="Limit candidates before EM/VEM to top K (default 32)",
    )
    classify_parser.add_argument(
        "--presence-caller",
        choices=["tdfdr", "hard_cutoff"],
        default=None,
        help="Presence caller strategy (tdfdr|hard_cutoff)",
    )
    classify_parser.add_argument(
        "--presence-q",
        type=float,
        default=None,
        help="Target q-value cutoff for presence calling",
    )
    classify_parser.add_argument(
        "--auto-q-tune",
        dest="auto_q_tune",
        action="store_true",
        default=None,
        help="Enable auto q tuning (default on)",
    )
    classify_parser.add_argument(
        "--no-auto-q-tune",
        dest="auto_q_tune",
        action="store_false",
        help="Disable auto q tuning",
    )
    classify_parser.add_argument(
        "--decoy-mode",
        default=None,
        help="Decoy generation mode (imcf-edge-shuffle)",
    )
    classify_parser.add_argument(
        "--decoy-reps",
        type=int,
        default=None,
        help="Number of IMCF edge-shuffle decoy replicates",
    )
    classify_parser.add_argument(
        "--exclusive-gamma",
        type=float,
        default=None,
        help="Exclusive edge weighting gamma",
    )
    classify_parser.add_argument(
        "--min-unique-evidence",
        type=int,
        default=None,
        help="Minimum IMCF unique-edge hits required to test presence",
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
    # NOTE: LCA 及后处理相关参数暂时废弃，内部逻辑继续沿用默认值
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
        "--taxonomy-kind",
        dest="taxonomy_kind",
        default="auto",
        choices=["auto", "ncbi", "gtdb"],
        help="分类期望的 taxonomy 数据源标识",
    )
    classify_parser.add_argument(
        "--taxonomy-version",
        dest="taxonomy_version",
        default="auto",
        help="分类期望的 taxonomy 数据版本标识",
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
    profile_parser.add_argument(
        "--taxonomy-kind",
        dest="taxonomy_kind",
        default="auto",
        choices=["auto", "ncbi", "gtdb"],
        help="Profile 使用的 taxonomy 数据源标识（auto 将尝试自动推断）",
    )
    profile_parser.add_argument(
        "--taxonomy-version",
        dest="taxonomy_version",
        default="auto",
        help="Profile 使用的 taxonomy 数据版本标识（可选）",
    )
    profile_parser.add_argument(
        "--taxonomy-info",
        dest="taxonomy_info",
        default=None,
        help="指向构建阶段生成的 tax.info 文件，用于 GTDB 等自定义 taxonomy",
    )
    profile_parser.add_argument(
        "--taxonomy-meta",
        dest="taxonomy_meta",
        default=None,
        help="指向 taxonomy.meta 元数据文件（若未提供将尝试在输入文件同目录查找）",
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
        if args.krona:
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
                taxonomy_version=args.taxonomy_version,
                taxonomy_info=args.taxonomy_info,
                taxonomy_meta=args.taxonomy_meta,
            )
        except ValueError as exc:
            print(f"Profile failed: {exc}")
            return 1
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
        if args.min_length != "auto":
            command.extend(["-l", str(args.min_length)])
        command.extend(["-t", str(args.threads)])
        command.extend(["--load-factor", str(args.load_factor)])
        command.extend(["-f", args.filter])
        command.extend(["--taxonomy-kind", args.taxonomy_kind])
        command.extend(["--taxonomy-version", args.taxonomy_version])
        if args.adaptive_cutoff:
            command.append("--adaptive-cutoff")
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
        if args.presence_caller is not None:
            command.extend(["--presence-caller", args.presence_caller])
        if args.presence_q is not None:
            command.extend(["--presence-q", str(args.presence_q)])
        if args.auto_q_tune is True:
            command.append("--auto-q-tune")
        elif args.auto_q_tune is False:
            command.append("--no-auto-q-tune")
        if args.decoy_mode is not None:
            command.extend(["--decoy-mode", args.decoy_mode])
        if args.decoy_reps is not None:
            command.extend(["--decoy-reps", str(args.decoy_reps)])
        if args.exclusive_gamma is not None:
            command.extend(["--exclusive-gamma", str(args.exclusive_gamma)])
        if args.min_unique_evidence is not None:
            command.extend(["--min-unique-evidence", str(args.min_unique_evidence)])
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
        command.extend(["--taxonomy-kind", args.taxonomy_kind])
        command.extend(["--taxonomy-version", args.taxonomy_version])
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
