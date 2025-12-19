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
        default=0,
        help="Maximum number of hashes per taxid",
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
        default=0,
        help="Maximum number of hashes per taxid",
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
        help="taxonomy 数据源标识 (auto/ncbi/gtdb)",
    )
    download_build_parser.add_argument(
        "--taxonomy-version",
        dest="taxonomy_version",
        default="auto",
        help="taxonomy 数据版本标识，例如 ncbi-taxdump-2025-09-15 或 gtdb-rs226",
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
        help="Limit candidates before EM/VEM to top K (default 16)",
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
        "--presence-noise",
        type=float,
        default=None,
        help="Override noise μ for coverage model; <=0 auto",
    )
    # presence-unique-deg 固定随数据库，分类侧不再暴露参数
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
        "--post-thres",
        type=float,
        default=None,
        help="Posterior acceptance threshold (default 0.56)",
    )
    classify_parser.add_argument(
        "--post-pi-min",
        type=float,
        default=None,
        help="Minimum global class weight (default 0.0005)",
    )
    # NOTE: LCA 及后处理相关参数暂时废弃，内部逻辑继续沿用默认值
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
        "-e",
        "--em",
        action="store_true",
        help="Use EM algorithm for classification (default)",
    )
    classify_parser.add_argument(
        "--em-iter", type=int, default=80, help="Number of EM iterations"
    )
    classify_parser.add_argument(
        "--em-threshold", type=float, default=0.001, help="EM threshold"
    )
    classify_parser.add_argument(
        "--em-temp",
        type=float,
        default=None,
        help="Softmax temperature inside EM/VEM",
    )
    classify_parser.add_argument(
        "--em-prior-strength",
        type=float,
        default=None,
        help="Dirichlet mass for abundance prior",
    )
    classify_parser.add_argument(
        "--em-coexist-penalty",
        type=float,
        default=None,
        help="Penalty applied when reads have near-tied taxa",
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
    profile_parser.add_argument(
        "--abundance-mode",
        dest="abundance_mode",
        default="soft_seq",
        choices=["soft_seq", "hard_seq", "top1", "soft_tax", "soft", "hard", "post_topk", "post"],
        help="丰度统计模式：soft_seq(默认，累加所有 taxid:count) / hard_seq(仅 top1 的 taxid:count) / top1(仅 top1，每条序列计 1) / post_topk(用 POST_TOPK posterior 将每条序列的总质量软分配到候选 taxid)。",
    )
    profile_parser.add_argument(
        "--hd-species-head-mass",
        type=float,
        default=97.0,
        help="高多样性样本：species 输出的 head mass（百分比），用于抑制非零长尾物种数（默认 97.0）",
    )
    profile_parser.add_argument(
        "--hd-genus-mode",
        choices=["none", "headmass", "topn_ratio"],
        default="none",
        help="高多样性样本：species 层级的 genus 内收敛模式（none/headmass/topn_ratio；默认 none）",
    )
    profile_parser.add_argument(
        "--hd-genus-head-mass",
        type=float,
        default=90.0,
        help="genus_headmass 模式：每个 genus 内保留覆盖的 head mass（百分比；默认 90.0）",
    )
    profile_parser.add_argument(
        "--hd-genus-max-species",
        type=int,
        default=3,
        help="每个 genus 最多保留的 species 数（默认 3）",
    )
    profile_parser.add_argument(
        "--hd-genus-topn",
        type=int,
        default=3,
        help="topn_ratio 模式：每个 genus 至少保留 topN（默认 3）",
    )
    profile_parser.add_argument(
        "--hd-genus-rel-min-in-genus",
        type=float,
        default=10.0,
        help="topn_ratio 模式：保留 genus 内相对质量 >= R%% 的 species（默认 10.0）",
    )
    profile_parser.add_argument(
        "--hd-shape-mode",
        choices=["none", "global_headmass", "genus_headmass"],
        default="none",
        help="高多样性样本：species 层级的 support–mass 形状特征过滤（none/global_headmass/genus_headmass；默认 none）",
    )
    profile_parser.add_argument(
        "--hd-shape-alpha",
        type=float,
        default=0.0,
        help="shape score: α in Score=logM+αlogN−βlog(M/N)−γ*skew+δlog(share_support)（默认 0.0）",
    )
    profile_parser.add_argument(
        "--hd-shape-beta",
        type=float,
        default=0.0,
        help="shape score: β in Score=logM+αlogN−βlog(M/N)−γ*skew+δlog(share_support)（默认 0.0）",
    )
    profile_parser.add_argument(
        "--hd-shape-gamma",
        type=float,
        default=0.0,
        help="shape score: γ in Score=logM+αlogN−βlog(M/N)−γ*skew+δlog(share_support)（默认 0.0）",
    )
    profile_parser.add_argument(
        "--hd-shape-delta",
        type=float,
        default=0.0,
        help="shape score: δ in Score=logM+αlogN−βlog(M/N)−γ*skew+δlog(share_support)（默认 0.0）",
    )
    profile_parser.add_argument(
        "--hd-shape-topk",
        type=int,
        default=8,
        help="shape score: 计算 skew 时每个物种保留的 topK contig 权重（默认 8；0 表示禁用 skew）",
    )
    profile_parser.add_argument(
        "--hd-shape-global-head-mass",
        type=float,
        default=97.0,
        help="global_headmass 模式：按 shape score 排序后保留覆盖的 head mass（百分比；默认 97.0）",
    )
    profile_parser.add_argument(
        "--hd-shape-genus-head-mass",
        type=float,
        default=90.0,
        help="genus_headmass 模式：每个 genus 内按 shape score 排序后保留覆盖的 head mass（百分比；默认 90.0）",
    )
    profile_parser.add_argument(
        "--hd-shape-genus-max-species",
        type=int,
        default=3,
        help="genus_headmass 模式：每个 genus 最多保留的 species 数（默认 3）",
    )
    profile_parser.add_argument(
        "--hd-intragenus-map",
        action="store_true",
        help="高多样性样本：使用 classify 输出的 POST_TOPK 做属内 MAP 纠错（默认关闭；需要重新跑一次 classify 生成 POST_TOPK）",
    )
    profile_parser.add_argument(
        "--hd-intragenus-topk",
        type=int,
        default=5,
        help="属内 MAP：使用 POST_TOPK 的前 K 个候选（默认 5）",
    )
    profile_parser.add_argument(
        "--hd-intragenus-beta",
        type=float,
        default=1.0,
        help="属内 MAP：全局属内先验的权重 β（默认 1.0）",
    )
    profile_parser.add_argument(
        "--hd-intragenus-gap",
        type=float,
        default=0.10,
        help="属内 MAP：仅当 top1-top2 posterior gap < 该阈值时尝试纠错（默认 0.10；设 0 可更激进）",
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if args.command == "classify" and not args.em:
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
                abundance_mode=args.abundance_mode,
                hd_species_head_mass=args.hd_species_head_mass,
                hd_genus_mode=args.hd_genus_mode,
                hd_genus_head_mass=args.hd_genus_head_mass,
                hd_genus_max_species=args.hd_genus_max_species,
                hd_genus_topn=args.hd_genus_topn,
                hd_genus_rel_min_in_genus=args.hd_genus_rel_min_in_genus,
                hd_shape_mode=args.hd_shape_mode,
                hd_shape_alpha=args.hd_shape_alpha,
                hd_shape_beta=args.hd_shape_beta,
                hd_shape_gamma=args.hd_shape_gamma,
                hd_shape_delta=args.hd_shape_delta,
                hd_shape_topk=args.hd_shape_topk,
                hd_shape_global_head_mass=args.hd_shape_global_head_mass,
                hd_shape_genus_head_mass=args.hd_shape_genus_head_mass,
                hd_shape_genus_max_species=args.hd_shape_genus_max_species,
                hd_intragenus_map=args.hd_intragenus_map,
                hd_intragenus_topk=args.hd_intragenus_topk,
                hd_intragenus_beta=args.hd_intragenus_beta,
                hd_intragenus_gap=args.hd_intragenus_gap,
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
        command.extend(["-M", str(args.max_hashes)])
        command.extend(
            ["--presence-unique-deg", str(args.presence_unique_deg)]
        )

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
        if args.presence_pi is not None:
            command.extend(["--presence-pi", str(args.presence_pi)])
        if args.presence_tau is not None:
            command.extend(["--presence-tau", str(args.presence_tau)])
        if args.presence_noise is not None:
            command.extend(["--presence-noise", str(args.presence_noise)])
        if args.presence_unique_deg is not None:
            command.extend(
                ["--presence-unique-deg", str(args.presence_unique_deg)]
            )
        if args.decoy_mode is not None:
            command.extend(["--decoy-mode", args.decoy_mode])
        if args.decoy_reps is not None:
            command.extend(["--decoy-reps", str(args.decoy_reps)])
        if args.exclusive_gamma is not None:
            command.extend(["--exclusive-gamma", str(args.exclusive_gamma)])
        if args.post_thres is not None:
            command.extend(["--post-thres", str(args.post_thres)])
        # posterior gating now只用 post_thres + post_pi_min
        if args.post_pi_min is not None:
            command.extend(["--post-pi-min", str(args.post_pi_min)])
        command.extend(["-t", str(args.threads)])
        command.extend(["-b", str(args.batch_size)])
        if args.em:
            command.append("-e")
            command.extend(["--em-iter", str(args.em_iter)])
            command.extend(["--em-threshold", str(args.em_threshold)])
        # EM 温度/先验/penalty 固定内部默认

    # Execute the command using the provided run function
    downloader.run_command(command)


def main():
    chimera_path = get_chimera_path()
    args = parse_arguments()
    run_chimera(args, chimera_path)


if __name__ == "__main__":
    main()
