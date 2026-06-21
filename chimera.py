import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path


NCBI_TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
NCBI_TAXDUMP_MD5_URL = NCBI_TAXDUMP_URL + ".md5"


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


def _is_auto_or_empty(value) -> bool:
    return value is None or str(value).strip() == "" or str(value).strip().lower() == "auto"


def _read_metadata_file(path: Path):
    values = {}
    if not path.exists():
        return values
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def _apply_build_taxonomy_meta_defaults(args) -> None:
    if not getattr(args, "input", None):
        return
    meta_path = Path(args.input).resolve().parent / "taxonomy.meta"
    meta = _read_metadata_file(meta_path)
    if not meta:
        return
    if _is_auto_or_empty(getattr(args, "taxonomy_kind", None)) and meta.get("taxonomy_kind"):
        args.taxonomy_kind = meta["taxonomy_kind"].lower()
    if _is_auto_or_empty(getattr(args, "taxonomy_version", None)) and meta.get("taxonomy_version"):
        args.taxonomy_version = meta["taxonomy_version"]
    if not getattr(args, "taxonomy_dir", None) and meta.get("taxonomy_dir"):
        args.taxonomy_dir = meta["taxonomy_dir"]


def _directory_has_taxdump(path: Path) -> bool:
    return path.is_dir() and (path / "nodes.dmp").is_file()


def _find_local_taxdump_archive(search_root: Path):
    candidates = [
        search_root / "taxdump.tar.gz",
    ]
    if search_root.is_dir():
        for child in search_root.iterdir():
            if child.is_dir():
                candidates.append(child / "taxdump.tar.gz")
    for candidate in candidates:
        if candidate.is_file() and candidate.stat().st_size > 0:
            return candidate
    return None


def _extract_taxdump_archive(archive_path: Path, extract_dir: Path) -> Path:
    if _directory_has_taxdump(extract_dir):
        return extract_dir
    extract_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path, "r:gz") as archive:
        for member in archive.getmembers():
            if not member.isfile():
                continue
            name = Path(member.name).name
            if not name:
                continue
            target = extract_dir / name
            source = archive.extractfile(member)
            if source is None:
                continue
            with source, target.open("wb") as out:
                shutil.copyfileobj(source, out)
    if not _directory_has_taxdump(extract_dir):
        raise RuntimeError(f"Taxdump archive did not contain nodes.dmp: {archive_path}")
    return extract_dir


def _download_url(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".part")
    try:
        with urllib.request.urlopen(url, timeout=120) as response, partial.open("wb") as out:
            shutil.copyfileobj(response, out)
        partial.replace(destination)
    finally:
        if partial.exists():
            partial.unlink()


def _download_ncbi_taxdump_archive(cache_dir: Path) -> Path:
    archive_path = cache_dir / "taxdump.tar.gz"
    md5_path = cache_dir / "taxdump.tar.gz.md5"
    if not archive_path.is_file() or archive_path.stat().st_size == 0:
        print(f"Downloading NCBI taxdump to {archive_path}")
        _download_url(NCBI_TAXDUMP_URL, archive_path)
    if not md5_path.is_file() or md5_path.stat().st_size == 0:
        _download_url(NCBI_TAXDUMP_MD5_URL, md5_path)
    expected = md5_path.read_text(encoding="utf-8").strip().split()[0].lower()
    digest = hashlib.md5()
    with archive_path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    actual = digest.hexdigest().lower()
    if expected != actual:
        archive_path.unlink(missing_ok=True)
        raise RuntimeError(
            f"Downloaded NCBI taxdump md5 mismatch: expected {expected}, got {actual}"
        )
    return archive_path


def _prepare_build_taxonomy_dir(args) -> None:
    if getattr(args, "no_local_resolution", False):
        return
    _apply_build_taxonomy_meta_defaults(args)
    if getattr(args, "taxonomy_dir", None):
        taxonomy_dir = Path(args.taxonomy_dir).resolve()
        if not _directory_has_taxdump(taxonomy_dir):
            raise FileNotFoundError(
                f"--taxonomy-dir must contain nodes.dmp: {taxonomy_dir}"
            )
        args.taxonomy_dir = str(taxonomy_dir)
        return

    kind = str(getattr(args, "taxonomy_kind", "auto") or "auto").strip().lower()
    if kind == "gtdb":
        return

    if not getattr(args, "input", None):
        return
    input_dir = Path(args.input).resolve().parent
    extract_dir = input_dir / "taxdump"
    if _directory_has_taxdump(extract_dir):
        args.taxonomy_dir = str(extract_dir)
        return
    if _directory_has_taxdump(input_dir):
        args.taxonomy_dir = str(input_dir)
        return

    meta = _read_metadata_file(input_dir / "taxonomy.meta")
    archive_value = meta.get("taxonomy_archive")
    archive_path = None
    if archive_value:
        archive_candidate = Path(archive_value).expanduser()
        if not archive_candidate.is_absolute():
            archive_candidate = input_dir / archive_candidate
        if archive_candidate.is_file() and archive_candidate.stat().st_size > 0:
            archive_path = archive_candidate
    if archive_path is None:
        archive_path = _find_local_taxdump_archive(input_dir)
    if archive_path is None:
        archive_path = _download_ncbi_taxdump_archive(input_dir)
    args.taxonomy_dir = str(_extract_taxdump_archive(archive_path, extract_dir))


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


def add_build_arguments(parser, require_input: bool) -> None:
    if require_input:
        parser.add_argument(
            "-i", "--input", required=True, help="Input file for building"
        )
    parser.add_argument(
        "-o", "--output", default="ChimeraDB", help="Output file for building"
    )
    parser.add_argument(
        "--strobe-k",
        type=int,
        default=28,
        help="Strobemer k-mer length",
    )
    parser.add_argument(
        "--strobe-order",
        type=int,
        default=2,
        help="Strobemer order (currently only 2 is supported)",
    )
    parser.add_argument(
        "--strobe-w-min",
        type=int,
        default=12,
        help="Strobemer minimum window",
    )
    parser.add_argument(
        "--strobe-w-max",
        type=int,
        default=32,
        help="Strobemer maximum window",
    )
    parser.add_argument(
        "-l",
        "--min-length",
        type=min_length_type,
        default="auto",
        help="Minimum sequence length for building (auto => strobemer minimum span)",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=default_threads(),
        help="Number of threads for building",
    )
    parser.add_argument(
        "--load-factor", type=float, default=0.85, help="IMCF filter load factor"
    )
    parser.add_argument(
        "--presence-unique-deg",
        type=int,
        default=1,
        help="Degree cutoff (<=) treated as unique signature for coverage meta",
    )
    parser.add_argument(
        "--taxonomy-kind",
        dest="taxonomy_kind",
        default="auto",
        choices=["auto", "ncbi", "gtdb"],
        help="Taxonomy source identifier (auto/ncbi/gtdb)",
    )
    parser.add_argument(
        "--taxonomy-version",
        dest="taxonomy_version",
        default="auto",
        help="Taxonomy version label, for example ncbi-taxdump-2025-09-15 or gtdb-rs226",
    )
    parser.add_argument(
        "--taxonomy-dir",
        dest="taxonomy_dir",
        default=None,
        help="Directory containing NCBI taxdump nodes.dmp; downloaded automatically for NCBI builds when omitted",
    )
    parser.add_argument(
        "--no-local-resolution",
        action="store_true",
        help="Do not build local read resolution (LPC) data",
    )


def append_build_command_args(command, args) -> None:
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
    command.extend(["--presence-unique-deg", str(args.presence_unique_deg)])
    command.extend(["--taxonomy-kind", str(args.taxonomy_kind)])
    command.extend(["--taxonomy-version", str(args.taxonomy_version)])
    if getattr(args, "taxonomy_dir", None):
        command.extend(["--taxonomy-dir", str(args.taxonomy_dir)])
    if getattr(args, "no_local_resolution", False):
        command.append("--no-local-resolution")


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
    add_build_arguments(build_parser, require_input=True)

    # Download and Build combined subcommand
    download_build_parser = subparsers.add_parser(
        "download_and_build",
        help="Download NCBI database sequences and resources and build a taxonomic sequence database",
    )
    add_build_arguments(download_build_parser, require_input=False)

    # Classify subcommand
    classify_parser = subparsers.add_parser("classify", help="Classify sequences")
    classify_inputs = classify_parser.add_mutually_exclusive_group(required=True)
    classify_inputs.add_argument(
        "-i",
        "--single",
        nargs="+",
        help="Input files for classifying (supports multiple files)",
    )
    classify_inputs.add_argument(
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
        "--no-local-resolution",
        action="store_true",
        help="Disable local read resolution (LPC) at classify time",
    )
    # Auxiliary profile utilities. Native abundance profiles are written by
    # `classify`; this command is for legacy aggregate conversion and Krona.
    profile_parser = subparsers.add_parser(
        "profile",
        help="Auxiliary profile conversion; classify already writes ChimeraProfile.tsv",
    )
    profile_parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=True,
        help="Existing ChimeraEvidence.tsv or aggregate taxid/count table",
    )
    profile_parser.add_argument(
        "-o", "--output", default="ChimeraProfile", help="Output prefix"
    )
    profile_parser.add_argument(
        "-k", "--krona", action="store_true", help="Generate Krona chart"
    )
    profile_parser.add_argument(
        "--taxonomy-kind",
        dest="taxonomy_kind",
        default="auto",
        choices=["auto", "ncbi", "gtdb"],
        help="Taxonomy source for auxiliary profile conversion",
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
        from src.profile import profile

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

    if chimera_path is None:
        chimera_path = get_chimera_path()
    command = [chimera_path]

    if args.version:
        command.append("--version")

    if args.command:
        command.append(args.command)

    if args.command == "build":
        _prepare_build_taxonomy_dir(args)
        append_build_command_args(command, args)

    elif args.command == "classify":
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        classify_output = output_dir / "ChimeraClassify.tsv"
        profile_output = output_dir / "ChimeraProfile.tsv"
        if args.single:
            command.extend(["-i"] + args.single)
        if args.paired:
            command.extend(["-p"] + args.paired)
        command.extend(["-o", str(classify_output)])
        command.extend(["-d", args.database])
        command.extend(["-t", str(args.threads)])
        command.extend(["-b", str(args.batch_size)])
        if getattr(args, "no_local_resolution", False):
            command.append("--no-local-resolution")
    if args.command == "classify":
        result = subprocess.run(command)
        if result.returncode != 0:
            return result.returncode
        if not profile_output.is_file():
            print(f"Classify did not produce native profile: {profile_output}")
            return 1
        return 0
    downloader = get_downloader()
    downloader.run_command(command)
    return 0


def main():
    args = parse_arguments()
    try:
        raise SystemExit(run_chimera(args) or 0)
    except (FileNotFoundError, RuntimeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
