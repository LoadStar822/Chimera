import argparse
import re
import sys
import subprocess
import sys
import time
import os
import shlex
from pathlib import Path
import gzip
import io
import tarfile
import pandas as pd
from multitax import NcbiTx, GtdbTx
from rich.console import Console
from rich.prompt import Confirm
from rich.table import Table
from rich.panel import Panel
from rich.text import Text
from rich.tree import Tree
from rich import box
import urllib.request
import urllib.error
import shutil
try:
    import questionary
except ImportError:  # pragma: no cover - optional dependency
    questionary = None

# Constants definition
VALID_DATABASES = ["genbank", "refseq"]
VALID_ORGANISM_GROUPS = [
    "archaea", "bacteria", "fungi", "human", "invertebrate", "metagenomes",
    "other", "plant", "protozoa", "vertebrate_mammalian", "vertebrate_other", "viral"
]
VALID_ASSEMBLY_LEVELS = ["complete genome", "chromosome", "scaffold", "contig"]
VALID_REFSEQ_CATEGORIES = ["reference genome", "na"]
VALID_FILE_TYPES = ["genomic.fna.gz", "assembly_report.txt", "protein.faa.gz", "genomic.gbff.gz"]
VALID_TAXONOMY_MODES = ["ncbi", "gtdb"]
VALID_DOWNLOADERS = ["wget", "curl"]
VALID_TAXONOMY_RANKS = [
    "superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"
]

DEFAULT_DATABASE = "refseq"
DEFAULT_ASSEMBLY_LEVEL = "complete genome"
DEFAULT_REFSEQ_CATEGORY = "reference genome"
DEFAULT_FILE_TYPE = "genomic.fna.gz"
DEFAULT_OUTPUT_DIR = "./genome_output"
DEFAULT_THREADS = "1"
DEFAULT_TAXONOMY_MODE = "ncbi"
DEFAULT_TAXONOMY_RANK = "species"
DEFAULT_RETRY_ATTEMPTS = "3"
DEFAULT_DOWNLOADER = "wget"


GTDB_METADATA_ACCESSION_COLUMNS = [
    "ncbi_refseq_assembly_accession",
    "ncbi_assembly_accession",
    "ncbi_genbank_assembly_accession",
    "accession"
]
GTDB_STANDARD_RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]

GTDB_ALLOWED_ORGANISM_GROUPS = ["archaea", "bacteria"]

TAXONOMY_METADATA_FILENAME = "taxonomy.meta"

GTDB_RELEASE_PATTERN = re.compile(r"(?:[_-](r\d+(?:\.\d+)?))", re.IGNORECASE)

GTDB_METADATA_SOURCES = [
    ("ar53_metadata.tar.gz", "https://data.gtdb.ecogenomic.org/releases/latest/ar53_metadata.tar.gz"),
    ("bac120_metadata.tar.gz", "https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tar.gz"),
    ("ar53_metadata.tsv.gz", "https://data.gtdb.ecogenomic.org/releases/latest/ar53_metadata.tsv.gz"),
    ("bac120_metadata.tsv.gz", "https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tsv.gz"),
]


def _normalize_gtdb_release(token):
    if not token:
        return None
    token = token.strip().lower()
    if not token:
        return None
    if token.startswith("rs"):
        body = token[2:]
    elif token.startswith("r"):
        body = token[1:]
    else:
        body = token
    if not body:
        return None
    if all(ch.isdigit() or ch == "." for ch in body):
        return f"rs{body}"
    return None


def _extract_release_from_part(part):
    lower = part.lower()
    match = GTDB_RELEASE_PATTERN.search(lower)
    if match:
        normalized = _normalize_gtdb_release(match.group(1))
        if normalized:
            return normalized
    if lower.startswith("rs") and len(lower) > 2:
        normalized = _normalize_gtdb_release(lower)
        if normalized:
            return normalized
    if lower.startswith("r") and len(lower) > 1:
        normalized = _normalize_gtdb_release(lower)
        if normalized:
            return normalized
    return None


def extract_gtdb_release_from_path(path):
    if not path:
        return None
    parts = Path(path).parts if not isinstance(path, str) else Path(path).parts
    for part in parts:
        release = _extract_release_from_part(part)
        if release:
            return release
    return None


def write_taxonomy_metadata(output_folder, taxonomy_kind, taxonomy_version):
    try:
        meta_path = Path(output_folder) / TAXONOMY_METADATA_FILENAME
        with meta_path.open("w", encoding="utf-8") as meta_file:
            meta_file.write(f"taxonomy_kind={taxonomy_kind}\n")
            meta_file.write(f"taxonomy_version={taxonomy_version}\n")
    except Exception:
        pass


console = Console()


def _split_defaults(defaults):
    if defaults is None:
        return set()
    if isinstance(defaults, (list, tuple, set)):
        return {str(item).strip() for item in defaults if str(item).strip()}
    return {part.strip() for part in str(defaults).split(",") if part.strip()}


def _render_choice_table(choices, defaults=None, title=None):
    if not choices:
        return
    default_set = _split_defaults(defaults)
    table = Table(
        title=title,
        show_header=False,
        expand=False,
        box=box.MINIMAL_DOUBLE_HEAD,
        padding=(0, 1),
    )
    table.add_column("编号", style="magenta", justify="right", width=4, no_wrap=True)
    table.add_column("可选值", style="cyan", no_wrap=True)

    for idx, choice in enumerate(choices, 1):
        style = "bold bright_cyan" if choice in default_set else "cyan"
        table.add_row(f"{idx:>2}", f"[{style}]{choice}[/]")
    console.print(table)


def _print_section(title, description=None):
    text = Text(title, style="bold cyan")
    if description:
        text.append("\n")
        text.append(description, style="dim")
    console.print(Panel(text, border_style="cyan", box=box.ROUNDED, expand=False))


def _format_bool(flag):
    return "✅ 是" if flag else "❌ 否"


def _format_value(value, placeholder="—"):
    if value is None:
        return placeholder
    if isinstance(value, str):
        value = value.strip()
        return value if value else placeholder
    return str(value)


def _format_multi_display(value, empty_label="全部"):
    if value is None:
        return empty_label
    parts = [part.strip() for part in str(value).split(",") if part.strip()]
    if not parts:
        return empty_label
    return ", ".join(parts)


def _interactive_available():
    return bool(questionary) and sys.stdin.isatty() and sys.stdout.isatty()


def sanitize_organism_group(value, taxonomy_mode):
    message = None
    if not taxonomy_mode or taxonomy_mode.lower() != "gtdb":
        return value or "", message
    tokens = []
    if isinstance(value, str) and value.strip():
        tokens = [item.strip() for item in value.split(",") if item.strip()]
    allowed = GTDB_ALLOWED_ORGANISM_GROUPS
    if not tokens:
        tokens = allowed
        message = "GTDB 模式仅支持 organism group=archaea,bacteria，已自动选择全部。"
    else:
        cleaned = [token for token in tokens if token in allowed]
        if not cleaned:
            cleaned = allowed
            message = "GTDB 模式仅支持 organism group=archaea,bacteria，已自动重置。"
        elif len(cleaned) != len(tokens):
            message = "GTDB 模式仅支持 organism group=archaea,bacteria，已忽略不兼容选项。"
        tokens = cleaned
    return ",".join(tokens), message


def sanitize_database(value, taxonomy_mode):
    if not value:
        value = DEFAULT_DATABASE
    tokens = [item.strip() for item in str(value).split(",") if item.strip() and item.strip() in VALID_DATABASES]
    if not tokens:
        tokens = [DEFAULT_DATABASE]
    tokens = list(dict.fromkeys(tokens))
    message = None
    if taxonomy_mode.lower() == "gtdb":
        if "refseq" not in tokens:
            tokens.insert(0, "refseq")
            message = "GTDB 模式建议至少包含 RefSeq 数据源；已自动附加。"
    return ",".join(tokens), message


def ensure_gtdb_metadata(output_dir, metadata_files=None, quiet=False):
    dest_dir = Path(output_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    existing = set()
    normalized_files = []
    if metadata_files:
        for path in metadata_files:
            p = Path(path)
            existing.add(p.name)
            normalized_files.append(p)

    for filename, url in GTDB_METADATA_SOURCES:
        if filename in existing:
            continue
        target = dest_dir / filename
        if target.exists():
            normalized_files.append(target)
            existing.add(filename)
            continue
        if not quiet:
            print(f" - Downloading GTDB metadata: {url}")
        try:
            with urllib.request.urlopen(url) as response, open(target, "wb") as out_file:
                shutil.copyfileobj(response, out_file)
            normalized_files.append(target)
            existing.add(filename)
        except Exception as exc:
            if not quiet:
                print(f"Warning: 下载 GTDB metadata 失败 {url}: {exc}")

    return normalized_files


def _tokenize_entries(response):
    return [token.strip() for token in re.split(r"[\s,]+", response) if token.strip()]


def _expand_numeric_token(token, choices):
    if not choices:
        return None
    if token == "*" or token.lower() in {"all", "全部"}:
        return list(choices)
    if re.fullmatch(r"\d+-\d+", token):
        start, end = token.split("-", 1)
        try:
            start_idx = int(start)
            end_idx = int(end)
        except ValueError:
            return []
        if start_idx > end_idx:
            start_idx, end_idx = end_idx, start_idx
        selected = []
        for idx in range(start_idx, end_idx + 1):
            if 1 <= idx <= len(choices):
                selected.append(choices[idx - 1])
        return selected
    if token.isdigit():
        idx = int(token)
        if 1 <= idx <= len(choices):
            return [choices[idx - 1]]
        return []
    return None


def _interpret_tokens(tokens, choices, allow_partial_numeric=False):
    resolved = []
    for token in tokens:
        expanded = _expand_numeric_token(token, choices)
        if expanded is None:
            resolved.append(token)
        elif expanded:
            resolved.extend(expanded)
        else:
            if allow_partial_numeric:
                continue
            return None
    return resolved


def ask_multi(label, choices, default=None, allow_empty=False, note=None):
    if choices and _interactive_available():
        default_set = _split_defaults(default)
        prompt = Text(label, style="bold cyan")
        if note:
            prompt.append(f"\n{note}", style="dim")
        console.print(prompt)
        try:
            answers = questionary.checkbox(
                "请选择 (空格选中，回车确认)",
                choices=[questionary.Choice(title=choice, value=choice, checked=choice in default_set) for choice in choices],
            ).ask()
        except Exception:
            answers = None
        if answers is None:
            if allow_empty and default is None:
                return ""
            if default is not None:
                return default
        if answers:
            return ",".join(answers)
        if allow_empty:
            return ""
        console.print("[yellow]未选择任何项，请至少选择一个选项[/yellow]")
        # fall through to text mode if user deselected all without default

    console.print(Text(label, style="bold cyan"))
    if note:
        console.print(Text(note, style="dim"))
    default_hint = None
    if default:
        default_hint = f"默认: {default}"
    elif allow_empty:
        default_hint = "默认: 全部/跳过"
    if default_hint:
        console.print(Text(default_hint, style="dim"))
    if choices:
        _render_choice_table(choices, default)
        console.print(Text("提示：可输入名称或编号（例如 '1 3 5-7' 或 '*'）", style="dim"))

    while True:
        response = console.input("[cyan]➜ [/cyan]").strip()
        if not response:
            if default is not None:
                return default
            if allow_empty:
                return ""
            console.print("[yellow]请输入至少一个值或使用默认值[/yellow]")
            continue

        tokens = _tokenize_entries(response)
        if not tokens:
            if allow_empty:
                return ""
            console.print("[yellow]请输入有效内容[/yellow]")
            continue

        interpreted = _interpret_tokens(tokens, choices)
        if interpreted is None:
            console.print("[red]包含无效编号，请重新输入[/red]")
            continue
        entries = interpreted
        if not entries:
            if allow_empty:
                return ""
            console.print("[yellow]请输入有效内容[/yellow]")
            continue

        if choices:
            normalized = []
            invalid = [entry for entry in entries if entry not in choices]
            if invalid:
                console.print(f"[red]无效选项: {', '.join(invalid)}[/red]")
                continue
            seen = set()
            for entry in entries:
                if entry not in seen:
                    normalized.append(entry)
                    seen.add(entry)
            entries = normalized
        else:
            seen = set()
            normalized = []
            for entry in entries:
                if entry not in seen:
                    normalized.append(entry)
                    seen.add(entry)
            entries = normalized

        if not entries:
            if allow_empty:
                return ""
            console.print("[yellow]请输入至少一个选项[/yellow]")
            continue

        return ",".join(entries)


def ask_text(label, default=None, allow_empty=True, note=None):
    if note:
        console.print(Text(note, style="dim"))
    if default not in (None, ""):
        console.print(Text(f"默认: {default}", style="dim"))
    elif default == "":
        console.print(Text("默认: 空", style="dim"))

    while True:
        response = console.input(f"[bold]{label}[/bold]\n[cyan]➜ [/cyan]").strip()
        if not response:
            if default is not None:
                return default
            if allow_empty:
                return ""
            console.print("[yellow]该项不能为空，请重新输入[/yellow]")
            continue
        return response


def ask_bool(label, default=False, note=None):
    if note:
        console.print(Text(note, style="dim"))
    return Confirm.ask(f"[bold]{label}[/bold]", default=default, console=console)


def ask_choice(label, choices, default=None, note=None):
    if choices and _interactive_available():
        prompt = Text(label, style="bold cyan")
        if note:
            prompt.append(f"\n{note}", style="dim")
        console.print(prompt)
        try:
            answer = questionary.select(
                "请选择",
                choices=[questionary.Choice(title=choice, value=choice) for choice in choices],
                default=default if default in choices else None,
            ).ask()
        except Exception:
            answer = None
        if answer is not None:
            return answer
        if default is not None:
            return default
    # fallback textual
    console.print(Text(label, style="bold cyan"))
    if note:
        console.print(Text(note, style="dim"))
    _render_choice_table(choices, default)
    if default is not None:
        console.print(Text(f"默认: {default}", style="dim"))
    console.print(Text("提示：可输入名称或编号，例如 '2'", style="dim"))

    while True:
        response = console.input("[cyan]➜ [/cyan]").strip()
        if not response:
            if default is not None:
                return default
            console.print("[yellow]请选择一个选项[/yellow]")
            continue
        tokens = _tokenize_entries(response)
        interpreted = _interpret_tokens(tokens, choices, allow_partial_numeric=True)
        if interpreted:
            for entry in interpreted:
                if entry in choices:
                    return entry
        if tokens:
            candidate = tokens[0]
            expanded = _expand_numeric_token(candidate, choices)
            if expanded:
                return expanded[0]
            if candidate in choices:
                return candidate
        console.print(f"[red]无效选项: {response}[/red]")


def validate_input(prompt, valid_options, default=None, allow_empty=False):
    """
    Validate user input
    
    Parameters:
    - prompt (str): Message to prompt user
    - valid_options (list): List of valid options
    - default (str): Default value
    - allow_empty (bool): Whether to allow empty input
    
    Returns:
    - str: Validated input value
    """
    while True:
        user_input = input(prompt).strip()
        if not user_input and default is not None:
            return default
        if allow_empty and not user_input:
            return ""
        entries = [entry.strip() for entry in user_input.split(",")]
        if all(entry in valid_options for entry in entries):
            return ",".join(entries)
        print(f"\nInvalid input. Please choose from: {', '.join(valid_options)}")


def prompt_user(options):
    """
    Interactively prompt user for download parameters
    
    Parameters:
    - options (argparse.Namespace): Command line arguments
    
    Returns:
    - dict: User input parameters
    """
    while True:
        console.rule("[bold cyan]Chimera 下载向导[/bold cyan]")
        intro = Text("欢迎使用 Chimera 数据库下载向导。\n", style="bold white")
        intro.append("我们将一步步收集参数，自动拼装 genome_updater.sh 指令。", style="dim")
        console.print(Panel(intro, border_style="cyan", box=box.DOUBLE, expand=False))
        console.print()

        if questionary is None and sys.stdin.isatty() and sys.stdout.isatty():
            console.print("[yellow]安装 questionary (pip install questionary) 可启用光标选择菜单。当前使用文本输入模式。[/yellow]")
            console.print()

        # Step 1: database & taxonomy
        _print_section("步骤 1/4 · 数据与分类", "选择要下载的数据库、物种范围与分类体系。")

        database_raw = getattr(options, "database", None)
        if database_raw is None:
            if _interactive_available():
                console.print(Text("数据库 (可多选)", style="bold cyan"))
                console.print(Text("回车可默认选择 RefSeq；若需覆盖 GenBank，可额外勾选。", style="dim"))
                try:
                    selection = questionary.checkbox(
                        "请选择 (空格选中，回车确认)",
                        choices=[questionary.Choice(title=choice, value=choice, checked=(choice == DEFAULT_DATABASE)) for choice in VALID_DATABASES],
                    ).ask()
                except Exception:
                    selection = None
                if not selection:
                    database_raw = DEFAULT_DATABASE
                else:
                    database_raw = ",".join(selection)
            else:
                database_raw = ask_multi(
                    "数据库 (可多选)",
                    VALID_DATABASES,
                    default=DEFAULT_DATABASE,
                    note="输入多个值请用逗号分隔。回车即使用默认 RefSeq。"
                )
        database_tokens = [item.strip() for item in str(database_raw).split(",") if item.strip() and item.strip() in VALID_DATABASES]
        if not database_tokens:
            database_tokens = [DEFAULT_DATABASE]
        database = ",".join(dict.fromkeys(database_tokens))
        use_refseq = "refseq" in database_tokens
        use_genbank = "genbank" in database_tokens

        organism_group = getattr(options, "organism_group", None)
        if organism_group is None:
            organism_group = ask_multi(
                "物种群组 (可选)",
                VALID_ORGANISM_GROUPS,
                default="",
                allow_empty=True,
                note="默认包含全部群组，指定多个时以逗号分隔。"
            )

        taxid = getattr(options, "taxid", None)
        if taxid is None:
            taxid = ask_text(
                "taxonomy ID (可留空)",
                default="",
                note="支持多个 ID，使用逗号分隔；留空表示不过滤。"
            )

        taxonomy_mode = getattr(options, "taxonomy_mode", None)
        if taxonomy_mode is None:
            taxonomy_mode = ask_choice(
                "taxonomy 模式",
                VALID_TAXONOMY_MODES,
                default=DEFAULT_TAXONOMY_MODE,
                note="选择 GTDB 时会自动解析官方 metadata 并写出 tax.info/target.tsv。"
            )
        else:
            taxonomy_mode = taxonomy_mode.lower()

        original_group = organism_group
        organism_group, group_message = sanitize_organism_group(organism_group, taxonomy_mode)
        if group_message:
            console.print(f"[yellow]{group_message}[/yellow]")
            if taxonomy_mode == "gtdb" and _interactive_available():
                organism_group = ask_multi(
                    "GTDB 可选物种群组",
                    GTDB_ALLOWED_ORGANISM_GROUPS,
                    default=organism_group or ",".join(GTDB_ALLOWED_ORGANISM_GROUPS),
                    note="GTDB 仅提供 archaea/bacteria 数据，请重新选择。"
                )

        if not organism_group:
            organism_group = ",".join(GTDB_ALLOWED_ORGANISM_GROUPS if taxonomy_mode == "gtdb" else [])

        database, db_message = sanitize_database(database, taxonomy_mode)
        if db_message:
            console.print(f"[yellow]{db_message}[/yellow]")
        database_tokens = [item.strip() for item in database.split(",") if item.strip()]
        use_refseq = "refseq" in database_tokens
        use_genbank = "genbank" in database_tokens

        taxonomy_rank = getattr(options, "taxonomy_rank", None) or ask_choice(
            "taxonomy 层级",
            VALID_TAXONOMY_RANKS,
            default=DEFAULT_TAXONOMY_RANK,
            note="决定 target.tsv 中目标节点的层级。"
        )

        limit_assembly = getattr(options, "limit_assembly", None)
        if limit_assembly is None:
            limit_assembly = ask_text(
                "下载数量限制",
                default="0",
                note="0 表示全部；也可输入 rank:number，例如 genus:3。"
            )

        console.print()

        # Step 2: filters
        _print_section("步骤 2/4 · 文件与过滤", "根据需要调整文件类型、RefSeq 分类及时间过滤。")
        filter_arguments = [
            getattr(options, "file_types", None),
            getattr(options, "refseq_category", None),
            getattr(options, "assembly_level", None),
            getattr(options, "start_date", None),
            getattr(options, "end_date", None),
            getattr(options, "custom_filter", None),
        ]
        filters_requested = any(value not in (None, "") for value in filter_arguments)
        if filters_requested or Confirm.ask("需要配置文件类型或其他过滤条件吗？", default=False, console=console):
            file_types = getattr(options, "file_types", None) or ask_multi(
                "文件类型",
                VALID_FILE_TYPES,
                default=DEFAULT_FILE_TYPE,
                note="可选择多个文件类型，逗号分隔。"
            )
            if use_refseq:
                refseq_note = "可多选；默认选择参考基因组。"
                refseq_category = getattr(options, "refseq_category", None) or ask_multi(
                    "RefSeq 分类",
                    VALID_REFSEQ_CATEGORIES,
                    default=DEFAULT_REFSEQ_CATEGORY,
                    note=refseq_note
                )
            else:
                refseq_category = ""
            assembly_level = getattr(options, "assembly_level", None) or ask_multi(
                "组装层级",
                VALID_ASSEMBLY_LEVELS,
                default=DEFAULT_ASSEMBLY_LEVEL,
                note="可多选；默认下载完整基因组。"
            )
            start_date = getattr(options, "start_date", None)
            if start_date is None:
                start_date = ask_text(
                    "发布日期下限 (YYYYMMDD，可留空)",
                    default="",
                    note="输入 YYYYMMDD 格式；留空表示不过滤。"
                )
            end_date = getattr(options, "end_date", None)
            if end_date is None:
                end_date = ask_text(
                    "发布日期上限 (YYYYMMDD，可留空)",
                    default="",
                    note="输入 YYYYMMDD 格式；留空表示不过滤。"
                )
            custom_filter = getattr(options, "custom_filter", None)
            if custom_filter is None:
                custom_filter = ask_text(
                    "自定义过滤 (colA:val1|colB:valX,valY)",
                    default="",
                    note="与 genome_updater.sh 语法一致；留空不启用。"
                )
        else:
            file_types = getattr(options, "file_types", None) or DEFAULT_FILE_TYPE
            if use_refseq:
                refseq_category = getattr(options, "refseq_category", None) or DEFAULT_REFSEQ_CATEGORY
            else:
                refseq_category = ""
            assembly_level = getattr(options, "assembly_level", None) or DEFAULT_ASSEMBLY_LEVEL
            start_date = getattr(options, "start_date", None) or ""
            end_date = getattr(options, "end_date", None) or ""
            custom_filter = getattr(options, "custom_filter", None) or ""

        console.print()

        # Step 3: runtime
        _print_section("步骤 3/4 · 下载与执行", "设置输出路径、线程数与运行模式。")
        output_dir = getattr(options, "output_dir", None) or ask_text(
            "输出目录",
            default=DEFAULT_OUTPUT_DIR,
            allow_empty=False,
            note="若目录已存在，可选择覆盖或继续。"
        )
        threads = getattr(options, "threads", None) or ask_text(
            "下载线程数",
            default=DEFAULT_THREADS,
            allow_empty=False,
            note="建议根据带宽与 I/O 能力设置。"
        )

        dry_run = ask_bool(
            "仅模拟 (dry-run) 不执行下载？",
            default=bool(getattr(options, "dry_run", False))
        )
        fix_mode = ask_bool(
            "只修复缺失/失败文件 (fix mode)？",
            default=bool(getattr(options, "fix_mode", False))
        )
        md5_check = ask_bool(
            "下载后校验 MD5？",
            default=True
        )

        console.print()

        _print_section("步骤 4/4 · 报告与高级选项", "可选择生成分析报告，并调整高级参数。")
        assembly_report = ask_bool(
            "生成 assembly accession 报告？",
            default=bool(getattr(options, "assembly_report", False))
        )
        sequence_report = ask_bool(
            "生成 sequence accession 报告？",
            default=bool(getattr(options, "sequence_report", False))
        )
        url_report = ask_bool(
            "生成下载 URL 报告？",
            default=bool(getattr(options, "url_report", False))
        )

        console.print()

        advanced_presets = [
            getattr(options, "version_label", None),
            getattr(options, "external_assembly", None),
            getattr(options, "alt_version_label", None),
            getattr(options, "retry_attempts", None),
            getattr(options, "conditional_exit", None),
            getattr(options, "ncbi_folders", None),
            getattr(options, "downloader", None),
            getattr(options, "delete_extra", None),
            getattr(options, "silent", None),
            getattr(options, "progress_only", None),
            getattr(options, "verbose", None),
            getattr(options, "debug", None),
        ]
        advanced_requested = any(value not in (None, "", False) for value in advanced_presets)

        if advanced_requested or Confirm.ask("需要继续配置高级选项吗？", default=False, console=console):
            version_label = getattr(options, "version_label", None)
            if version_label is None:
                version_label = ask_text(
                    "版本标签 (可留空自动使用时间戳)",
                    default=""
                )
            external_assembly = getattr(options, "external_assembly", None)
            if external_assembly is None:
                external_assembly = ask_text(
                    "外部 assembly_summary.txt (可留空)",
                    default=""
                )
            alt_version_label = getattr(options, "alt_version_label", None)
            if alt_version_label is None:
                alt_version_label = ask_text(
                    "备用版本标签 (可留空)",
                    default=""
                )
            retry_attempts = getattr(options, "retry_attempts", None) or ask_text(
                "下载重试次数",
                default=DEFAULT_RETRY_ATTEMPTS,
                allow_empty=False
            )
            conditional_exit = getattr(options, "conditional_exit", None) or ask_text(
                "失败阈值 (0 表示关闭)",
                default="0",
                allow_empty=False
            )
            ncbi_folders = ask_bool(
                "按照 NCBI FTP 目录结构存放文件？",
                default=bool(getattr(options, "ncbi_folders", False))
            )
            downloader = getattr(options, "downloader", None) or ask_choice(
                "下载工具",
                VALID_DOWNLOADERS,
                default=DEFAULT_DOWNLOADER,
                note="wget 支持断点续传；curl 更适合受限环境。"
            )
            delete_extra = ask_bool(
                "允许删除输出目录中多余的常规文件？",
                default=bool(getattr(options, "delete_extra", False))
            )
            silent = ask_bool(
                "开启静默输出 (silent)？",
                default=bool(getattr(options, "silent", False))
            )
            progress_only = ask_bool(
                "仅显示下载进度 (progress-only)？",
                default=bool(getattr(options, "progress_only", False))
            )
            verbose = ask_bool(
                "开启详细日志 (verbose)？",
                default=bool(getattr(options, "verbose", False))
            )
            debug = ask_bool(
                "开启调试模式 (debug)？",
                default=bool(getattr(options, "debug", False))
            )
        else:
            version_label = getattr(options, "version_label", None) or ""
            external_assembly = getattr(options, "external_assembly", None) or ""
            alt_version_label = getattr(options, "alt_version_label", None) or ""
            retry_attempts = getattr(options, "retry_attempts", None) or DEFAULT_RETRY_ATTEMPTS
            conditional_exit = getattr(options, "conditional_exit", None) or "0"
            ncbi_folders = bool(getattr(options, "ncbi_folders", False))
            downloader = getattr(options, "downloader", None) or DEFAULT_DOWNLOADER
            delete_extra = bool(getattr(options, "delete_extra", False))
            silent = bool(getattr(options, "silent", False))
            progress_only = bool(getattr(options, "progress_only", False))
            verbose = bool(getattr(options, "verbose", False))
            debug = bool(getattr(options, "debug", False))

        console.print()

        # Assemble configuration
        config = {
            "database": database,
            "organism_group": organism_group,
            "taxid": taxid,
            "file_types": file_types,
            "refseq_category": refseq_category,
            "assembly_level": assembly_level,
            "start_date": start_date,
            "end_date": end_date,
            "custom_filter": custom_filter,
            "taxonomy_mode": taxonomy_mode,
            "taxonomy_rank": taxonomy_rank,
            "limit_assembly": limit_assembly,
            "output_dir": output_dir,
            "threads": threads,
            "dry_run": dry_run,
            "fix_mode": fix_mode,
            "md5_check": md5_check,
            "assembly_report": assembly_report,
            "sequence_report": sequence_report,
            "url_report": url_report,
            "version_label": version_label,
            "external_assembly": external_assembly,
            "alt_version_label": alt_version_label,
            "retry_attempts": retry_attempts,
            "conditional_exit": conditional_exit,
            "ncbi_folders": ncbi_folders,
            "downloader": downloader,
            "delete_extra": delete_extra,
            "silent": silent,
            "progress_only": progress_only,
            "verbose": verbose,
            "debug": debug,
        }

        summary = Tree("[bold cyan]参数确认[/bold cyan]")
        basic_branch = summary.add("[bold]基础[/bold]")
        basic_branch.add(f"数据库: {_format_multi_display(database)}")
        basic_branch.add(f"物种群组: {_format_multi_display(organism_group)}")
        basic_branch.add(f"Taxonomy 模式: {taxonomy_mode}")
        basic_branch.add(f"Taxonomy 层级: {taxonomy_rank}")
        basic_branch.add(f"Taxonomy ID: {_format_value(taxid, '全部')}")
        basic_branch.add(f"数量限制: {limit_assembly}")

        filter_branch = summary.add("[bold]过滤[/bold]")
        filter_branch.add(f"文件类型: {_format_multi_display(file_types, '默认')}")
        filter_branch.add(f"RefSeq 分类: {_format_multi_display(refseq_category)}")
        filter_branch.add(f"组装层级: {_format_multi_display(assembly_level)}")
        filter_branch.add(f"发布日期范围: {_format_value(start_date, '不限')} ~ {_format_value(end_date, '不限')}")
        filter_branch.add(f"自定义过滤: {_format_value(custom_filter, '未设置')}")

        runtime_branch = summary.add("[bold]运行[/bold]")
        runtime_branch.add(f"输出目录: {output_dir}")
        runtime_branch.add(f"线程数: {threads}")
        runtime_branch.add(f"Dry-run: {_format_bool(dry_run)}")
        runtime_branch.add(f"Fix mode: {_format_bool(fix_mode)}")
        runtime_branch.add(f"MD5 校验: {_format_bool(md5_check)}")

        report_branch = summary.add("[bold]报告[/bold]")
        report_branch.add(f"Assembly 报告: {_format_bool(assembly_report)}")
        report_branch.add(f"Sequence 报告: {_format_bool(sequence_report)}")
        report_branch.add(f"URL 报告: {_format_bool(url_report)}")

        advanced_branch = summary.add("[bold]高级[/bold]")
        advanced_items = []
        if version_label:
            advanced_items.append(f"版本标签: {version_label}")
        if external_assembly:
            advanced_items.append(f"外部 assembly_summary: {external_assembly}")
        if alt_version_label:
            advanced_items.append(f"备用版本标签: {alt_version_label}")
        if retry_attempts != DEFAULT_RETRY_ATTEMPTS:
            advanced_items.append(f"重试次数: {retry_attempts}")
        if conditional_exit != "0":
            advanced_items.append(f"失败阈值: {conditional_exit}")
        if ncbi_folders:
            advanced_items.append("NCBI 目录结构: 已开启")
        if downloader != DEFAULT_DOWNLOADER:
            advanced_items.append(f"下载工具: {downloader}")
        if delete_extra:
            advanced_items.append("删除多余文件: 开启")
        if silent:
            advanced_items.append("静默输出: 开启")
        if progress_only:
            advanced_items.append("仅显示进度: 开启")
        if verbose:
            advanced_items.append("详细日志: 开启")
        if debug:
            advanced_items.append("调试模式: 开启")

        if not advanced_items:
            advanced_branch.add("使用默认配置")
        else:
            for item in advanced_items:
                advanced_branch.add(item)

        console.print(summary)
        console.print()

        if Confirm.ask("确认以上设置并继续执行下载流程吗？", default=True, console=console):
            return config

        if not Confirm.ask("需要重新配置参数吗？", default=True, console=console):
            console.print("[red]已取消操作。[/red]")
            sys.exit(1)

        console.print("[yellow]重新开始参数配置...\n[/yellow]")


def build_command(options):
    """
    Build genome_updater.sh command based on options
    
    Parameters:
    - options (argparse.Namespace): Command line arguments
    
    Returns:
    - str: Built command string
    """
    cmd = ["genome_updater.sh"]

    # Database options
    if options.database:
        cmd.append(f"-d '{options.database}'")

    # Organism options
    if options.organism_group:
        cmd.append(f"-g '{options.organism_group}'")

    if options.taxid:
        cmd.append(f"-T '{options.taxid}'")

    # File options
    if options.file_types:
        cmd.append(f"-f '{options.file_types}'")

    # Filter options
    if options.refseq_category:
        cmd.append(f"-c '{options.refseq_category}'")
        
    if options.assembly_level:
        cmd.append(f"-l '{options.assembly_level}'")
    
    if options.start_date:
        cmd.append(f"-D '{options.start_date}'")
    
    if options.end_date:
        cmd.append(f"-E '{options.end_date}'")
    
    if options.custom_filter:
        cmd.append(f"-F '{options.custom_filter}'")
    
    # Taxonomy options
    if options.taxonomy_mode:
        cmd.append(f"-M '{options.taxonomy_mode}'")
    
    
    if options.limit_assembly:
        if str(options.limit_assembly) == "0":
            pass
        elif ':' in str(options.limit_assembly):
            cmd.append(f"-A '{options.limit_assembly}'")
        else:
            cmd.append(f"-A 'species:{options.limit_assembly}'")
    
    # Keep taxonomy database
    cmd.append("-a")
    
    # Run options
    if options.output_dir:
        cmd.append(f"-o '{options.output_dir}'")

    if options.threads:
        cmd.append(f"-t {options.threads}")

    if options.dry_run:
        cmd.append("-k")

    if options.fix_mode:
        cmd.append("-i")

    if options.md5_check:
        cmd.append("-m")
    
    # Report options
    if getattr(options, 'assembly_report', False):
        cmd.append("-u")
    
    if getattr(options, 'sequence_report', False):
        cmd.append("-r")
    
    if getattr(options, 'url_report', False):
        cmd.append("-p")
    
    # Misc options
    if options.version_label:
        cmd.append(f"-b '{options.version_label}'")
    
    if options.external_assembly:
        cmd.append(f"-e '{options.external_assembly}'")
    
    if options.alt_version_label:
        cmd.append(f"-B '{options.alt_version_label}'")
    
    if options.retry_attempts:
        cmd.append(f"-R {options.retry_attempts}")
    
    if options.conditional_exit:
        cmd.append(f"-n {options.conditional_exit}")
    
    # NCBI folder structure
    if getattr(options, 'ncbi_folders', False):
        cmd.append("-N")
    
    if options.downloader:
        cmd.append(f"-L '{options.downloader}'")
    
    if getattr(options, 'delete_extra', False):
        cmd.append("-x")
    
    # Output mode settings
    if getattr(options, 'silent', False):
        cmd.append("-s")
    
    if getattr(options, 'progress_only', False):
        cmd.append("-w")
    
    if getattr(options, 'verbose', False):
        cmd.append("-V")
    
    if getattr(options, 'debug', False):
        cmd.append("-Z")

    return " ".join(cmd)


def run_command(cmd, ret_stdout=False, shell=False, quiet=False, debug=False, timeout=None):
    """
    Run command and handle results
    
    Parameters:
    - cmd (str or list): Command to execute
    - ret_stdout (bool): Whether to capture and return stdout
    - shell (bool): Whether to run command in shell
    - quiet (bool): Whether to suppress stderr output
    - debug (bool): Whether to enable debug mode
    - timeout (int): Command execution timeout in seconds
    
    Returns:
    - str: Captured stdout if ret_stdout is True, otherwise None
    """
    errcode = 1
    stdout = None

    try:
        if debug:
            print(f"Executing command: {cmd}")

        # Split command string into list if not using shell
        if isinstance(cmd, str) and not shell:
            cmd = shlex.split(cmd)

        # Run command
        process = subprocess.Popen(
            cmd,
            shell=shell,
            universal_newlines=True,
            stdout=subprocess.PIPE if ret_stdout else sys.stderr,
            stderr=subprocess.PIPE if quiet else sys.stderr
        )

        if ret_stdout:
            stdout, stderr = process.communicate(timeout=timeout)
            if stderr and not quiet:
                print(stderr)
        else:
            process.wait(timeout=timeout)

        errcode = process.returncode
        if errcode != 0:
            raise subprocess.CalledProcessError(errcode, cmd)

    except subprocess.CalledProcessError as e:
        print(f"Command execution failed with code {e.returncode}: {e.cmd}")
        sys.exit(e.returncode)
    except subprocess.TimeoutExpired:
        print(f"Command timed out after {timeout} seconds: {cmd}")
        process.kill()
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error running command: {cmd}")
        print(str(e))
        sys.exit(errcode)

    return stdout


def validate_input_files(input_files_folder, file_types, quiet=False, input_recursive=False):
    """
    Validate input files
    
    Parameters:
    - input_files_folder (list): List of input files or folders
    - file_types (list): List of file types
    - quiet (bool): Whether to suppress output
    - input_recursive (bool): Whether to search folders recursively
    
    Returns:
    - set: Set of valid files
    """
    valid_input_files = set()

    def _to_absolute(path_str):
        try:
            return str(Path(path_str).resolve(strict=False))
        except OSError:
            return os.path.abspath(path_str)

    for raw_path in input_files_folder:
        path = Path(raw_path).resolve(strict=False)
        if path.is_file():
            valid_input_files.add(str(path))
        elif path.is_dir():
            files_in_dir = 0

            if input_recursive:
                for file_type in file_types:
                    pattern = f"*{file_type}"
                    for candidate in path.rglob(pattern):
                        if candidate.is_file():
                            files_in_dir += 1
                            valid_input_files.add(str(candidate.resolve(strict=False)))
            else:
                for candidate in path.iterdir():
                    if not candidate.is_file():
                        continue
                    candidate_name = candidate.name
                    for file_type in file_types:
                        if candidate_name.endswith(file_type):
                            files_in_dir += 1
                            valid_input_files.add(str(candidate.resolve(strict=False)))
                            break

            if not quiet:
                print(
                    f" - {files_in_dir} valid files [types: {', '.join(file_types)}"
                    + (", recursive" if input_recursive else "")
                    + f"] found in {path}"
                )
        else:
            if not quiet:
                print(f" - Skipping invalid file/folder: {raw_path}")

    # Normalize collected file paths to absolute paths
    normalized_files = {_to_absolute(f) for f in valid_input_files}

    if not quiet:
        print(f" - Total valid files: {len(normalized_files)}")

    return normalized_files


def load_assembly_accession(input_files):
    """
    Load assembly accessions from input files
    
    Parameters:
    - input_files (set): Set of input files
    
    Returns:
    - pandas.DataFrame: DataFrame containing assembly accession information
    """
    info_cols = ["file", "target", "node", "specialization", "specialization_name"]

    # Regular expression to match GenBank/RefSeq assembly accessions
    pattern = re.compile(r"GC[AF]_[0-9]+\.[0-9]+")

    # Extract assembly accession or use filename as target
    data = []
    for f in input_files:
        abs_file = str(Path(f).resolve(strict=False))
        match = pattern.search(abs_file)
        target = match.group() if match else os.path.basename(abs_file)
        data.append((target, abs_file))

    # Create DataFrame with info_cols columns
    info = pd.DataFrame(data, columns=["target", "file"])

    # Add other columns with default values
    for col in info_cols:
        if col not in info.columns:
            info[col] = None

    # Clean data
    info.dropna(how="all", inplace=True)
    info.dropna(subset=["target"], inplace=True)
    info.drop_duplicates(subset=["target"], inplace=True)

    # Set 'target' as index
    info.set_index('target', inplace=True)

    return info


def parse_assembly_summary(info, assembly_summary, level="species"):
    """
    Parse assembly summary file
    
    Parameters:
    - info (pandas.DataFrame): Assembly information DataFrame
    - assembly_summary (str): Path to assembly summary file
    - level (str): Taxonomy level
    
    Returns:
    - dict: Assembly summary counts
    """
    count_assembly_summary = {}
    unique_acc = set(info.index)

    # Detect header lines
    header_lines = 0
    with open(assembly_summary, 'r') as ass_sum:
        for line in ass_sum:
            if line[0] == "#":
                header_lines += 1
            else:
                break

    # Read assembly summary file
    tmp_acc_node = pd.read_csv(
        assembly_summary,
        sep="\t",
        header=None,
        skiprows=header_lines,
        usecols=[0, 5, 7, 8],
        names=["target", "node", "organism_name", "infraspecific_name"],
        index_col="target",
        converters={"target": lambda x: x if x in unique_acc else None, "node": str}
    )
    
    # Keep only used sequence IDs
    tmp_acc_node = tmp_acc_node[tmp_acc_node.index.notnull()]

    # Save count
    count_assembly_summary[assembly_summary] = tmp_acc_node.shape[0]

    # Create specialization
    if level == "assembly":
        # Handle infraspecific_name prefix
        tmp_acc_node["infraspecific_name"] = tmp_acc_node["infraspecific_name"].replace(
            "^[a-z]+=", "", regex=True).fillna("")

        # Build name
        def build_name(n):
            if n.organism_name.endswith(n.infraspecific_name):
                return n.organism_name
            else:
                return n.organism_name + " " + n.infraspecific_name

        # Add infraspecific_name suffix
        tmp_acc_node["specialization_name"] = tmp_acc_node[["organism_name", "infraspecific_name"]].apply(
            lambda n: build_name(n), axis=1)
        tmp_acc_node["specialization"] = tmp_acc_node.index

    # Merge nodes and specializations found by target (accession)
    if count_assembly_summary[assembly_summary]:
        info.update(tmp_acc_node)
    
    del tmp_acc_node
    return count_assembly_summary


def get_file_info(options, info, build_output_folder, assembly_summary):
    """
    Get file information
    
    Parameters:
    - options (argparse.Namespace): Command line arguments
    - info (pandas.DataFrame): Assembly information DataFrame
    - build_output_folder (str): Build output folder
    - assembly_summary (str): Path to assembly summary file
    """
    # Parse assembly summary file
    start_time = time.time()
    quiet = getattr(options, 'silent', False) or getattr(options, 'progress_only', False)
    
    if not quiet:
        print("Parsing assembly summary file")
    
    count_assembly_summary = parse_assembly_summary(info, assembly_summary)

    # Output parsing results
    if not quiet:
        for assembly_summary_file, count in count_assembly_summary.items():
            file_name = assembly_summary_file.split("/")[-1]
            print(f" - Found {count} entries in {file_name}")

        print(f" - Done in {time.time() - start_time:.2f} seconds.\n")


def normalize_taxonomy_rank(taxonomy_rank, taxonomy_mode):
    """
    Normalize requested taxonomy rank based on taxonomy mode.
    """
    rank = (taxonomy_rank or "").lower()
    if taxonomy_mode == "gtdb":
        if rank == "superkingdom":
            return "domain"
        if rank == "strain":
            return "species"
    return rank


def gtdb_lineage_to_node(lineage):
    """
    Return the deepest GTDB node with a valid name from a lineage string.
    """
    if not isinstance(lineage, str):
        return None, None
    for entry in reversed([part.strip() for part in lineage.split(";")]):
        if len(entry) >= 4 and entry[3:]:
            return entry, entry[3:]
    return None, None


def extract_gtdb_accession_candidates(row):
    """
    Extract possible assembly accession identifiers from a GTDB metadata row.
    """
    candidates = set()
    for column in GTDB_METADATA_ACCESSION_COLUMNS:
        value = row.get(column, "")
        if not isinstance(value, str):
            continue
        for token in re.split(r"[;,]", value):
            acc = token.strip()
            if not acc or acc.lower() == "na":
                continue
            candidates.add(acc)
            if acc.startswith(("RS_", "GB_")) and len(acc) > 3:
                candidates.add(acc[3:])
    return candidates


def parse_gtdb_metadata_stream(stream, targets, release_holder=None):
    """
    Parse a GTDB metadata TSV stream and populate mappings for the provided targets set.
    """
    mapping = {}
    names = {}
    for chunk in pd.read_csv(stream, sep="\t", dtype=str, chunksize=100000):
        if release_holder is not None and release_holder[0] is None:
            for column in chunk.columns:
                if "release" not in column.lower():
                    continue
                series = chunk[column].dropna()
                for value in series:
                    normalized = _normalize_gtdb_release(str(value))
                    if normalized:
                        release_holder[0] = normalized
                        break
                if release_holder[0] is not None:
                    break
        chunk = chunk.fillna("")
        for _, row in chunk.iterrows():
            lineage = row.get("gtdb_taxonomy", "")
            node, name = gtdb_lineage_to_node(lineage)
            if not node:
                continue
            for acc in extract_gtdb_accession_candidates(row):
                if acc in targets:
                    mapping[acc] = node
                    names[acc] = name
                    targets.remove(acc)
            if not targets:
                break
        if not targets:
            break
    return mapping, names


def parse_gtdb_metadata_file(metadata_path, targets, release_holder=None):
    """
    Parse GTDB metadata from a local file (.tsv, .tsv.gz, or .tar.gz).
    """
    mapping = {}
    names = {}
    lower_name = metadata_path.name.lower()

    if release_holder is not None and release_holder[0] is None:
        release_holder[0] = extract_gtdb_release_from_path(metadata_path)

    if lower_name.endswith(".tar.gz"):
        with tarfile.open(metadata_path, "r:gz") as tar:
            members = [
                member for member in tar.getmembers()
                if member.isfile() and member.name.lower().endswith(".tsv")
            ]
            for member in members:
                if release_holder is not None and release_holder[0] is None:
                    release_holder[0] = extract_gtdb_release_from_path(member.name)
                with tar.extractfile(member) as fh:
                    if fh is None:
                        continue
                    with io.TextIOWrapper(fh, encoding="utf-8") as text_stream:
                        local_map, local_names = parse_gtdb_metadata_stream(text_stream, targets, release_holder)
                        mapping.update(local_map)
                        names.update(local_names)
                if not targets:
                    break
    elif lower_name.endswith(".tsv.gz"):
        with gzip.open(metadata_path, "rt") as fh:
            mapping, names = parse_gtdb_metadata_stream(fh, targets, release_holder)
    else:
        with open(metadata_path, "r") as fh:
            mapping, names = parse_gtdb_metadata_stream(fh, targets, release_holder)

    return mapping, names


def find_gtdb_metadata_files(search_roots):
    """
    Locate GTDB metadata files underneath the provided directories.
    """
    metadata_files = []
    seen = set()
    for root in search_roots:
        if not root:
            continue
        root_path = Path(root)
        if not root_path.is_dir():
            continue
        for dirpath, dirnames, filenames in os.walk(root_path):
            if "files" in dirnames:
                dirnames.remove("files")
            for filename in filenames:
                lower = filename.lower()
                if lower.endswith(("metadata.tar.gz", "metadata.tsv.gz", "metadata.tsv")):
                    path = Path(dirpath) / filename
                    if path not in seen:
                        metadata_files.append(path)
                        seen.add(path)
    metadata_files.sort()
    return metadata_files


def load_gtdb_taxonomy(info, metadata_files, quiet=False):
    """
    Load GTDB taxonomy assignments using downloaded metadata files.
    """
    remaining = set(info.index.astype(str))
    node_map = {}
    name_map = {}
    release_holder = [None]

    for metadata_path in metadata_files:
        if not remaining:
            break
        if release_holder[0] is None:
            release_holder[0] = extract_gtdb_release_from_path(metadata_path)
        if not quiet:
            print(f" - Parsing GTDB metadata: {metadata_path}")
        mapping, names = parse_gtdb_metadata_file(metadata_path, remaining, release_holder)
        if mapping:
            node_map.update(mapping)
            name_map.update(names)
            remaining -= set(mapping.keys())

    return node_map, name_map, remaining, release_holder[0]


def validate_taxonomy(info, tax, taxonomy_rank="species", taxonomy_mode=DEFAULT_TAXONOMY_MODE):
    """
    Validate taxonomy: convert to user-specified taxonomy rank nodes
    
    Parameters:
    - info (pandas.DataFrame): Assembly information DataFrame
    - tax (multitax.MultiTax): Taxonomy object
    - taxonomy_rank (str): Taxonomy rank to use (default: species)
    """
    normalized_rank = normalize_taxonomy_rank(taxonomy_rank, taxonomy_mode)

    if taxonomy_mode == "gtdb":
        def to_rank(node):
            if not isinstance(node, str):
                return None
            lineage_nodes = tax.lineage(node, ranks=GTDB_STANDARD_RANKS)
            rank_lookup = dict(zip(GTDB_STANDARD_RANKS, lineage_nodes))
            target_node = rank_lookup.get(normalized_rank)
            if target_node and target_node != tax.undefined_node:
                return target_node
            # Fallback to the deepest available classified rank
            for rank in reversed(GTDB_STANDARD_RANKS):
                candidate = rank_lookup.get(rank)
                if candidate and candidate != tax.undefined_node:
                    return candidate
            return None

        info["node"] = info["node"].apply(to_rank)
    else:
        info["node"] = info["node"].apply(lambda n: tax.parent_rank(n, normalized_rank))

    # Remove entries that could not be resolved to a valid taxonomy node
    na_entries = info["node"].isna().sum()
    if na_entries > 0:
        info.dropna(subset=["node"], inplace=True)


def remove_files(files):
    """
    Remove files
    
    Parameters:
    - files (str or list): Files to remove
    """
    if isinstance(files, str):
        files = [files]
    
    for f in files:
        if os.path.isfile(f):
            os.remove(f)


def write_tax(tax_file, info, tax):
    """
    Write taxonomy file
    
    Parameters:
    - tax_file (str): Path to taxonomy file
    - info (pandas.DataFrame): Assembly information DataFrame
    - tax (multitax.MultiTax): Taxonomy object
    """
    # Remove existing file
    remove_files(tax_file)
    
    # Write taxonomy to file
    tax.write(tax_file)

    # Load written taxonomy file
    tax_df = pd.read_csv(
        tax_file, 
        names=["node", "parent", "rank", "name"], 
        delimiter='\t', 
        dtype=str
    )

    # Write back updated DataFrame
    tax_df.to_csv(tax_file, sep="\t", header=False, index=False)

    # Create target.tsv file
    target_tsv_path = os.path.join(os.path.dirname(tax_file), "target.tsv")

    # Keep only file and node columns
    info_subset = info[["file", "node"]]

    # Write info_subset content to target_tsv_path
    with open(target_tsv_path, 'w') as f:
        f.write(info_subset.to_csv(sep="\t", header=False, index=False))


def setup_output_directory(options):
    """
    Set up output directory
    
    Parameters:
    - options (argparse.Namespace): Command line arguments
    
    Returns:
    - tuple: (output_folder, assembly_summary, tmp_folder)
    """
    output_folder = options.output_dir
    
    # Handle output directory
    if options.fix_mode:
        print("Fix mode enabled. Will re-download existing files.")
    elif os.path.exists(output_folder):
        while True:
            confirmation = input(
                f"Output folder '{output_folder}' exists. Do you want to clear it [y], continue [c], or cancel [n]? [y/c/n]: ").lower()
            
            if confirmation == 'y':
                print("Clearing output folder...")
                if os.name == 'nt':  # Windows
                    os.system(f"rmdir /s /q {output_folder}")
                else:  # Unix/Linux
                    os.system(f"rm -rf {output_folder}")
                os.makedirs(output_folder)  # Recreate folder after clearing
                break
            elif confirmation == 'c':
                print("Continuing without clearing output folder.")
                break
            elif confirmation == 'n':
                print("Operation cancelled.")
                sys.exit(1)
            else:
                print("Invalid option. Please choose 'y', 'c', or 'n'.")
    else:
        print(f"Creating output directory: {output_folder}")
        os.makedirs(output_folder, exist_ok=True)

    # Set up paths
    assembly_summary = os.path.join(output_folder, "assembly_summary.txt")
    tmp_folder = os.path.join(output_folder, "tmp")
    os.makedirs(tmp_folder, exist_ok=True)

    # Check existing files
    if os.path.isfile(assembly_summary) and os.path.getsize(assembly_summary) > 0:
        print(f"Assembly summary file exists: {assembly_summary}")

    return output_folder, assembly_summary, tmp_folder


def download(interactive=False, raw_args=None):
    """
    Download NCBI genome data
    
    Parameters:
    - interactive (bool): Whether to use interactive mode
    
    Returns:
    - argparse.Namespace: Command line arguments
    """
    # Create argument parser
    parser = argparse.ArgumentParser(
        description='Download NCBI genome data',
        epilog='''Examples:
        python download.py -d refseq -g "archaea,bacteria" -l "complete genome" -f genomic.fna.gz -o ./output
        ''',
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Add arguments
    # Database options
    parser.add_argument("-d", "--database", 
                       help=f"Database to download (default: {DEFAULT_DATABASE})\nOptions: {', '.join(VALID_DATABASES)}")
    
    # Organism options
    parser.add_argument("-g", "--organism-group",
                       help=f"Organism groups to download (comma-separated)\nOptions: {', '.join(VALID_ORGANISM_GROUPS)}")
    
    parser.add_argument("-T", "--taxid",
                       help="Taxonomy IDs to download (comma-separated)\nExample: 562,623 (for -M ncbi) or s__Escherichia coli (for -M gtdb)")
    
    # File options
    parser.add_argument("-f", "--file-types",
                       help=f"File types to download (default: {DEFAULT_FILE_TYPE})\nOptions: {', '.join(VALID_FILE_TYPES)}")
    
    # Filter options
    parser.add_argument("-c", "--refseq-category",
                       help=f"RefSeq categories to download (default: {DEFAULT_REFSEQ_CATEGORY})\nOptions: {', '.join(VALID_REFSEQ_CATEGORIES)}")
    
    parser.add_argument("-l", "--assembly-level",
                       help=f"Assembly levels to download (default: {DEFAULT_ASSEMBLY_LEVEL})\nOptions: {', '.join(VALID_ASSEMBLY_LEVELS)}")
    
    parser.add_argument("-D", "--start-date",
                       help="Start date (>=), based on sequence release date. Format YYYYMMDD")
    
    parser.add_argument("-E", "--end-date",
                       help="End date (<=), based on sequence release date. Format YYYYMMDD")
    
    parser.add_argument("-F", "--custom-filter",
                       help="Custom filter for assembly summary in format colA:val1|colB:valX,valY (case insensitive)")
    
    # Taxonomy options
    parser.add_argument("-M", "--taxonomy-mode",
                       help=f"Taxonomy mode (default: {DEFAULT_TAXONOMY_MODE})\nOptions: {', '.join(VALID_TAXONOMY_MODES)}")
    
    parser.add_argument("-Y", "--taxonomy-rank",
                       help=f"Taxonomy rank to use (default: {DEFAULT_TAXONOMY_RANK})\nOptions: {', '.join(VALID_TAXONOMY_RANKS)}")
    
    parser.add_argument("-A", "--limit-assembly",
                       help="Limit number of assemblies for each selected taxa. 0 for all.\nSelection by ranks also supported with rank:number (e.g., genus:3)")
    
    # Run options
    parser.add_argument("-o", "--output-dir",
                       help=f"Output directory path (default: {DEFAULT_OUTPUT_DIR})\nWarning: Will be cleared if exists")
    
    parser.add_argument("-t", "--threads",
                       help=f"Number of download threads (default: {DEFAULT_THREADS})")
    
    parser.add_argument("-k", "--dry-run",
                       action="store_true",
                       help="Enable dry run mode (no actual downloads)")
    
    parser.add_argument("-i", "--fix-mode",
                       action="store_true", 
                       help="Only re-download incomplete or failed files")
    
    parser.add_argument("-m", "--md5-check",
                       action="store_true",
                       help="Check MD5 of downloaded files")
    
    # Report options
    parser.add_argument("-u", "--assembly-report",
                       action="store_true",
                       help="Generate updated assembly accessions report")
    
    parser.add_argument("-r", "--sequence-report",
                       action="store_true",
                       help="Generate updated sequence accessions report")
    
    parser.add_argument("-p", "--url-report",
                       action="store_true",
                       help="Generate URLs report (success/failed)")
    
    # Misc options
    parser.add_argument("-b", "--version-label",
                       help="Version label (default: current timestamp)")
    
    parser.add_argument("-e", "--external-assembly",
                       help="External assembly_summary.txt file to recover data from")
    
    parser.add_argument("-B", "--alt-version-label",
                       help="Alternative version label to use as current version")
    
    parser.add_argument("-R", "--retry-attempts",
                       help=f"Number of attempts to retry download files in batches (default: {DEFAULT_RETRY_ATTEMPTS})")

    parser.add_argument("-n", "--conditional-exit",
                        help="Conditional exit status based on number of failures accepted")

    parser.add_argument("-N", "--ncbi-folders",
                        action="store_true",
                        help="Output files in folders like NCBI ftp structure")

    parser.add_argument("-L", "--downloader",
                        help=f"Downloader to use (default: {DEFAULT_DOWNLOADER})\nOptions: {', '.join(VALID_DOWNLOADERS)}")

    parser.add_argument("-x", "--delete-extra",
                        action="store_true",
                        help="Allow deletion of regular extra files found in output folder")

    parser.add_argument("-s", "--silent",
                        action="store_true",
                        help="Silent output")

    parser.add_argument("-w", "--progress-only",
                        action="store_true",
                        help="Silent output with download progress only")

    parser.add_argument("-V", "--verbose",
                        action="store_true",
                        help="Verbose log")

    parser.add_argument("-Z", "--debug",
                        action="store_true",
                        help="Print debug information and run in debug mode")

    
    if raw_args is not None:
        options = parser.parse_args(raw_args)
        interactive = False
    elif interactive:
        parsed_args = parser.parse_args([])
        options = argparse.Namespace(**prompt_user(parsed_args))
    else:
        # 命令行模式
        options = parser.parse_args()

    if not getattr(options, 'taxonomy_mode', None):
        options.taxonomy_mode = DEFAULT_TAXONOMY_MODE
    else:
        options.taxonomy_mode = options.taxonomy_mode.lower()
    taxonomy_mode = options.taxonomy_mode

    sanitized_group, group_message = sanitize_organism_group(getattr(options, 'organism_group', None), taxonomy_mode)
    if group_message and not (getattr(options, 'silent', False) or getattr(options, 'progress_only', False)):
        print(group_message)
    options.organism_group = sanitized_group

    sanitized_database, db_message = sanitize_database(getattr(options, 'database', None), taxonomy_mode)
    if db_message and not (getattr(options, 'silent', False) or getattr(options, 'progress_only', False)):
        print(db_message)
    options.database = sanitized_database
    
    # Build command
    command = build_command(options)
    
    # Print command information if not quiet
    quiet = getattr(options, 'silent', False) or getattr(options, 'progress_only', False)
    if not quiet:
        print(f"\nGenerated command: {command}")
        print(f"\nSetting up output directory...")
    
    # Set up output directory
    output_folder, assembly_summary, tmp_folder = setup_output_directory(options)
    
    if not quiet:
        print(f"Output folder: {output_folder}")
        print(f"Assembly summary file: {assembly_summary}")
        print(f"Temporary folder: {tmp_folder}")

    # Execute download
    try:
        # Download phase
        quiet = getattr(options, 'silent', False) or getattr(options, 'progress_only', False)
        if not quiet:
            download_start = time.time()
            print("\nDownloading data...")
            print(f"Database: {options.database if options.database else DEFAULT_DATABASE}")
            print(f"Organism group: {options.organism_group if options.organism_group else 'All'}")
            print(f"File types: {options.file_types if options.file_types else DEFAULT_FILE_TYPE}")
            if options.taxid:
                print(f"Taxonomy IDs: {options.taxid}")
            print(f"Using {options.threads if options.threads else DEFAULT_THREADS} threads for download")
        
        # Execute download command
        if not quiet:
            print("\nExecuting genome_updater.sh command...")
        run_command(command, shell=True, quiet=quiet)
        
        if not quiet:
            download_end = time.time()
            print(f"\nDownload completed in {download_end - download_start:.2f} seconds.")
            print(f"Downloaded files will be processed now...")

        # Processing phase
        process_start = time.time()
        if not quiet:
            print("\nProcessing data...")
            print(f"Starting post-download processing at {time.strftime('%Y-%m-%d %H:%M:%S')}")
        
        # Get file paths
        if not quiet:
            print("\nLocating downloaded files...")
        assembly_root_rel = os.path.dirname(os.readlink(assembly_summary))
        file_folder = os.path.join(options.output_dir, assembly_root_rel, "files")
        if not quiet:
            print(f"File folder path: {file_folder}")
            print(f"Looking for files with extensions: {options.file_types}")
        
        input_file = validate_input_files([file_folder], options.file_types.split(","), quiet, input_recursive=True)

        if not input_file:
            if not quiet:
                print("No valid files found. Exiting.")
            sys.exit(1)
        else:
            if not quiet:
                print("Valid files found.")

        # Load taxonomy information
        if not quiet:
            print("\nLoading taxonomy information...")
            print(f"Taxonomy mode: {taxonomy_mode}")
            tax_start_time = time.time()
        
        tax = GtdbTx() if taxonomy_mode == "gtdb" else NcbiTx()
        if not quiet:
            print("Loading assembly accessions from input files...")
        
        info = load_assembly_accession(input_file)
        info.index = info.index.astype(str)
        if not quiet:
            print(f"Found {len(info)} assembly accessions")
            print(f"Taxonomy loading completed in {time.time() - tax_start_time:.2f} seconds")
        
        if info.empty:
            if not quiet:
                print("No valid assembly accessions found. Exiting.")
            sys.exit(1)
        
        # Get file information and validate taxonomy
        taxonomy_version_value = "ncbi-taxdump"
        taxonomy_kind_value = taxonomy_mode

        if taxonomy_mode == "gtdb":
            if not quiet:
                print("\nResolving GTDB metadata for taxonomy assignments...")
            assembly_root_path = Path(options.output_dir) / assembly_root_rel
            metadata_roots = [
                assembly_root_path,
                Path(options.output_dir),
                Path(tmp_folder)
            ]
            metadata_files = find_gtdb_metadata_files(metadata_roots)
            metadata_files = ensure_gtdb_metadata(options.output_dir, metadata_files, quiet=quiet)
            if not metadata_files:
                raise FileNotFoundError(
                    "未能找到 GTDB metadata 文件，无法建立 taxonomy 映射。"
                )
            node_map, name_map, remaining, release_token = load_gtdb_taxonomy(info, metadata_files, quiet=quiet)
            if not node_map:
                raise RuntimeError("GTDB metadata 解析失败，未取得任何 taxonomy 映射。")
            info["node"] = info.index.to_series().map(node_map)
            if name_map:
                info["organism_name"] = info.index.to_series().map(name_map)
            mapped_count = info["node"].notna().sum()
            if not quiet:
                print(f"GTDB metadata 匹配到 {mapped_count}/{len(info)} 条目")
            if remaining and not quiet:
                print(f"警告：{len(remaining)} 个 assembly 在 GTDB metadata 中缺失 taxonomy 信息")
            taxonomy_version_value = f"gtdb-{release_token}" if release_token else "gtdb-unknown"
            if not quiet:
                print(f"GTDB release: {taxonomy_version_value}")
        else:
            if not quiet:
                print("\nParsing assembly summary file and extracting metadata...")
            get_file_info(options, info, tmp_folder, assembly_summary)
        
        requested_rank = getattr(options, 'taxonomy_rank', DEFAULT_TAXONOMY_RANK) or DEFAULT_TAXONOMY_RANK
        normalized_rank = normalize_taxonomy_rank(requested_rank, taxonomy_mode)

        if not quiet:
            if normalized_rank != (requested_rank or "").lower():
                print(f"\nValidating taxonomy using rank: {requested_rank} (标准化为 {normalized_rank})")
            else:
                print(f"\nValidating taxonomy using rank: {requested_rank}")
            taxonomy_start_time = time.time()
        
        validate_taxonomy(info, tax, requested_rank, taxonomy_mode)
        
        if not quiet:
            valid_entries = len(info)
            print(f"Valid taxonomy entries after validation: {valid_entries}")
            print(f"Taxonomy validation completed in {time.time() - taxonomy_start_time:.2f} seconds")
        
        # Filter taxonomy
        if not quiet:
            print("\nFiltering taxonomy nodes...")
            filter_start_time = time.time()
        
        unique_nodes = info["node"].unique()
        if not quiet:
            print(f"Found {len(unique_nodes)} unique taxonomy nodes")
        
        tax.filter(unique_nodes)
        
        if not quiet:
            print(f"Taxonomy filtering completed in {time.time() - filter_start_time:.2f} seconds")
        
        # Write taxonomy file
        if not quiet:
            print("\nWriting taxonomy information to file...")

        tax_file_path = os.path.join(output_folder, "tax.info")
        write_tax(tax_file_path, info, tax)
        write_taxonomy_metadata(output_folder, taxonomy_kind_value, taxonomy_version_value)

        if not quiet:
            print(f"Taxonomy information written to: {tax_file_path}")
            print(f"Target information written to: {os.path.join(os.path.dirname(tax_file_path), 'target.tsv')}")
        
        if not quiet:
            process_end = time.time()
            print(f"\nProcessing completed in {process_end - process_start:.2f} seconds.")

        # Clean up temporary folder
        if not quiet:
            print("\nCleaning up temporary files...")
        
        if os.name == 'nt':  # Windows
            os.system(f"rmdir /s /q {tmp_folder}")
        else:  # Unix/Linux
            os.system(f"rm -rf {tmp_folder}")
            
        if not quiet:
            print(f"Removed temporary folder: {tmp_folder}")

    except Exception as e:
        print(f"Error during download process: {str(e)}")
        sys.exit(1)

    return options

def main():
    """
    Main function to execute when script is run directly
    """
    try:
        # Auto-detect mode: use argument mode if command line arguments exist, otherwise use interactive mode
        options = download(interactive=None)
        print("\nDownload completed!")
        print(f"Output directory: {options.output_dir}")
    except KeyboardInterrupt:
        print("\nOperation interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\nError occurred: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
