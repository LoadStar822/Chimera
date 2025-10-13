import csv
import json
import math
import warnings
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

try:  # pragma: no cover - 兼容未安装 ete3 的环境
    from ete3 import NCBITaxa  # type: ignore
except ImportError:  # pragma: no cover - 运行环境缺失 ete3
    NCBITaxa = None  # type: ignore


LEVEL_ALIASES = {
    "acellular root": ("acellular root",),
    "realm": ("realm",),
    "domain": ("domain", "superkingdom"),
    "clade": ("clade",),
    "phylum": ("phylum",),
    "class": ("class",),
    "order": ("order",),
    "family": ("family",),
    "genus": ("genus",),
    "species": ("species",),
    "strain": ("strain",),
}

UNCLASSIFIED = "unclassified"

GLOBAL_EM_MAX_ITER = 60
GLOBAL_EM_TOL = 3e-6
GLOBAL_DIRICHLET_ALPHA = 0.05
GLOBAL_LENGTH_NORMALIZATION = True
GLOBAL_MIN_LENGTH_BP = 1000
GLOBAL_MIN_VIRAL_LENGTH_BP = 200


@dataclass
class TaxonRecord:
    taxid: int
    parent: Optional[int]
    rank: str
    name: str


class TaxonomyResolver:
    def __init__(self, tax_info: Optional[Path] = None) -> None:
        self._records: Dict[int, TaxonRecord] = {}
        self._ncbi: Optional[Any] = None
        if tax_info and tax_info.exists():
            self._load_tax_info(tax_info)
        if NCBITaxa is not None:
            try:
                self._ncbi = NCBITaxa()
            except Exception:  # pragma: no cover - ete3 初始化失败
                self._ncbi = None

    def _load_tax_info(self, path: Path) -> None:
        with open(path, "r", encoding="utf-8", errors="ignore") as infile:
            for line in infile:
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 4:
                    continue
                try:
                    taxid = int(parts[0])
                    parent = int(parts[1]) if parts[1] else None
                except ValueError:
                    continue
                rank = parts[2].lower()
                name = parts[3]
                self._records[taxid] = TaxonRecord(
                    taxid=taxid,
                    parent=parent if parent and parent > 0 else None,
                    rank=rank,
                    name=name,
                )

    def get_lineage(self, taxid: int) -> List[TaxonRecord]:
        if taxid <= 0:
            return []
        if self._records:
            lineage = self._lineage_from_records(taxid)
            if lineage and lineage[-1].taxid == taxid:
                return lineage
        return self._lineage_from_ncbi(taxid)

    def _lineage_from_records(self, taxid: int) -> List[TaxonRecord]:
        lineage: List[TaxonRecord] = []
        current = taxid
        seen = set()
        while current and current not in seen:
            seen.add(current)
            node = self._records.get(current)
            if not node:
                break
            lineage.append(node)
            if not node.parent or node.parent == current:
                break
            current = node.parent
        return list(reversed(lineage))

    def _lineage_from_ncbi(self, taxid: int) -> List[TaxonRecord]:
        if self._ncbi is None:
            return []
        try:
            lineage = self._ncbi.get_lineage(taxid)
            rank_map = self._ncbi.get_rank(lineage)
            name_map = self._ncbi.get_taxid_translator(lineage)
        except Exception:
            return []
        records: List[TaxonRecord] = []
        parent: Optional[int] = None
        for node in lineage:
            rank = (rank_map.get(node) or "").lower()
            name = name_map.get(node) or ""
            records.append(TaxonRecord(taxid=node, parent=parent, rank=rank, name=name))
            parent = node
        return records

    def has_taxonomy(self) -> bool:
        return bool(self._records) or self._ncbi is not None


def _normalize_rank(rank: str) -> Optional[str]:
    rank_lower = rank.lower()
    for level, aliases in LEVEL_ALIASES.items():
        if rank_lower in aliases:
            return level
    return None


def _should_mark_unclassified(level: str, has_level: Dict[str, bool]) -> bool:
    if level == "domain":
        return not has_level.get("acellular root", False)
    if level == "realm":
        return has_level.get("acellular root", False)
    if level == "acellular root":
        return has_level.get("realm", False) or not has_level.get("domain", False)
    return True


def calculate_shannon_index(taxon_dict: Dict[str, float]) -> float:
    total_count = sum(taxon_dict.values())
    if total_count == 0:
        return 0.0
    shannon_index = 0.0
    for count in taxon_dict.values():
        p_i = count / total_count
        if p_i > 0:
            shannon_index -= p_i * math.log(p_i)
    return shannon_index


def calculate_simpson_index(taxon_dict: Dict[str, float]) -> float:
    total_count = sum(taxon_dict.values())
    if total_count == 0:
        return 0.0
    simpson_index = 0.0
    for count in taxon_dict.values():
        p_i = count / total_count
        simpson_index += p_i ** 2
    return 1.0 - simpson_index


def _parse_taxid_weights(line: str) -> Tuple[List[Tuple[int, float]], bool]:
    line = line.strip()
    if not line:
        return [], False
    parts = [token.strip() for token in line.split("\t") if token.strip()]
    if len(parts) < 2:
        return [], False
    items: List[Tuple[int, float]] = []
    has_unclassified = False
    for token in parts[1:]:
        if not token:
            continue
        lowered = token.lower()
        if lowered == UNCLASSIFIED:
            has_unclassified = True
            continue
        if token.startswith("POST_") and "=" in token:
            token = token.split("=", 1)[1]
        if ":" in token:
            taxid_token, weight_token = token.split(":", 1)
        else:
            taxid_token, weight_token = token, ""
        try:
            taxid = int(taxid_token)
        except ValueError:
            continue
        if weight_token:
            try:
                weight = float(weight_token)
            except ValueError:
                weight = 1.0
        else:
            weight = 1.0
        items.append((taxid, weight))
    return items, has_unclassified


def _collect_taxid_weights(
    input_files: Iterable[str],
    keep_per_read: bool = False,
) -> Tuple[Dict[int, float], float, int, Optional[List[List[Tuple[int, float]]]]]:
    taxid_weights: Dict[int, float] = defaultdict(float)
    unclassified_reads = 0.0
    total_lines = 0
    per_read: Optional[List[List[Tuple[int, float]]]] = [] if keep_per_read else None
    for input_file in input_files:
        with open(input_file, "r", errors="ignore") as infile:
            for line in infile:
                total_lines += 1
                items, has_unclassified = _parse_taxid_weights(line)
                cleaned: List[Tuple[int, float]] = []
                for taxid, weight in items:
                    if math.isnan(weight) or weight <= 0.0:
                        continue
                    cleaned.append((taxid, weight))
                if not cleaned:
                    if has_unclassified:
                        unclassified_reads += 1.0
                    continue
                sum_weight = sum(weight for _, weight in cleaned)
                if sum_weight > 0.0:
                    if sum_weight > 1.0:
                        factor = 1.0 / sum_weight
                        candidate_weights = [(taxid, weight * factor) for taxid, weight in cleaned]
                    else:
                        if not has_unclassified:
                            factor = 1.0 / sum_weight
                            candidate_weights = [(taxid, weight * factor) for taxid, weight in cleaned]
                        else:
                            candidate_weights = cleaned
                else:
                    candidate_weights = []
                if not candidate_weights:
                    if has_unclassified:
                        unclassified_reads += 1.0
                    continue
                line_weight = 0.0
                for taxid, weight in candidate_weights:
                    taxid_weights[taxid] += weight
                    line_weight += weight
                if has_unclassified:
                    residue = max(0.0, 1.0 - line_weight)
                    unclassified_reads += residue
                if per_read is not None:
                    per_read.append(candidate_weights)
    return taxid_weights, unclassified_reads, total_lines, per_read


def _effective_length(length_bp: Optional[int], *, min_length: int = GLOBAL_MIN_LENGTH_BP) -> Optional[float]:
    if not length_bp or length_bp <= 0:
        return None
    floor = float(min_length)
    return max(float(length_bp), floor) / 1000.0


def _is_viral_lineage(resolver: Optional[TaxonomyResolver], taxid: int) -> bool:
    if resolver is None or taxid <= 0:
        return False
    lineage = resolver.get_lineage(taxid)
    for record in lineage:
        name_lower = (record.name or "").lower()
        if record.taxid == 10239:
            return True
        if record.rank in ("realm", "clade") and name_lower.endswith("viria"):
            return True
        if "virus" in name_lower or "viroid" in name_lower:
            return True
    return False


def _run_global_em(
    per_read_candidates: Sequence[Sequence[Tuple[int, float]]],
    initial_weights: Dict[int, float],
    tax_lengths: Optional[Dict[int, int]] = None,
    resolver: Optional[TaxonomyResolver] = None,
    max_iter: int = GLOBAL_EM_MAX_ITER,
    tol: float = GLOBAL_EM_TOL,
) -> Dict[int, float]:
    total_weight = sum(initial_weights.values())
    if total_weight <= 0.0:
        return dict(initial_weights)
    pi: Dict[int, float] = {}
    epsilon = 1e-12
    for taxid, weight in initial_weights.items():
        if weight > 0.0:
            pi[taxid] = weight / total_weight
    # Ensure every candidate taxid has a prior
    missing = False
    for read in per_read_candidates:
        for taxid, _ in read:
            if taxid not in pi:
                pi[taxid] = epsilon
                missing = True
    if missing:
        norm = sum(pi.values())
        if norm > 0.0:
            inv_norm = 1.0 / norm
            for taxid in pi:
                pi[taxid] *= inv_norm

    if not pi:
        return dict(initial_weights)

    dirichlet_alpha = GLOBAL_DIRICHLET_ALPHA
    if dirichlet_alpha < 0.0:
        dirichlet_alpha = 0.0

    length_scale_cache: Dict[int, Optional[float]] = {}
    viral_cache: Dict[int, bool] = {}

    def _resolve_length_scale(taxid: int) -> Optional[float]:
        if not tax_lengths:
            return None
        if taxid in length_scale_cache:
            return length_scale_cache[taxid]
        length_bp = tax_lengths.get(taxid)
        if not length_bp or length_bp <= 0:
            length_scale_cache[taxid] = None
            return None
        min_floor = GLOBAL_MIN_LENGTH_BP
        if resolver and length_bp < GLOBAL_MIN_LENGTH_BP:
            if taxid not in viral_cache:
                viral_cache[taxid] = _is_viral_lineage(resolver, taxid)
            if viral_cache.get(taxid):
                min_floor = GLOBAL_MIN_VIRAL_LENGTH_BP
        scale = _effective_length(length_bp, min_length=min_floor)
        length_scale_cache[taxid] = scale
        return scale

    for _ in range(max_iter):
        alpha_sum = dirichlet_alpha * len(pi)
        new_pi: Dict[int, float] = defaultdict(float)
        total_responsibility = 0.0
        for read in per_read_candidates:
            tmp: List[Tuple[int, float]] = []
            denom = 0.0
            for taxid, weight in read:
                if weight <= 0.0 or math.isnan(weight):
                    continue
                prior = pi.get(taxid, 0.0)
                if prior <= 0.0:
                    continue
                value = prior * weight
                if GLOBAL_LENGTH_NORMALIZATION and tax_lengths:
                    scale = _resolve_length_scale(taxid)
                    if scale:
                        value /= scale
                if value <= 0.0:
                    continue
                tmp.append((taxid, value))
                denom += value
            if denom <= 0.0:
                denom = sum(max(weight, 0.0) for _, weight in read)
                if denom <= 0.0:
                    continue
                for taxid, weight in read:
                    if weight <= 0.0 or math.isnan(weight):
                        continue
                    resp = weight / denom
                    new_pi[taxid] += resp
                    total_responsibility += resp
                continue
            inv_denom = 1.0 / denom
            for taxid, value in tmp:
                resp = value * inv_denom
                new_pi[taxid] += resp
                total_responsibility += resp
        if total_responsibility <= 0.0:
            break

        inv_total = 1.0 / (total_responsibility + alpha_sum) if (total_responsibility + alpha_sum) > 0.0 else 0.0
        max_delta = 0.0
        for taxid, value in new_pi.items():
            value = (value + dirichlet_alpha) * inv_total
            max_delta = max(max_delta, abs(value - pi.get(taxid, 0.0)))
            new_pi[taxid] = value
        for taxid, value in pi.items():
            if taxid not in new_pi:
                smoothed = dirichlet_alpha * inv_total if inv_total > 0.0 else 0.0
                max_delta = max(max_delta, abs(value - smoothed))
                new_pi[taxid] = smoothed
        pi = dict(new_pi)

        if pi:
            tiny_threshold = 1e-12
            tiny_keys = [taxid for taxid, value in pi.items() if value < tiny_threshold]
            if len(tiny_keys) == len(pi):
                keep_taxid = max(pi.items(), key=lambda item: item[1])[0]
                tiny_keys = [taxid for taxid in tiny_keys if taxid != keep_taxid]
            removed_mass = 0.0
            for taxid in tiny_keys:
                removed_mass += pi.pop(taxid, 0.0)
            if pi:
                if removed_mass > 0.0:
                    norm = sum(pi.values())
                    if norm > 0.0:
                        inv_norm = 1.0 / norm
                        for taxid in list(pi.keys()):
                            pi[taxid] *= inv_norm
                max_delta = max(max_delta, removed_mass)

        if max_delta < tol:
            break

    refined = {taxid: pi_val * total_weight for taxid, pi_val in pi.items()}
    for taxid in initial_weights:
        refined.setdefault(taxid, 0.0)
    return refined


def _load_snapshot(snapshot_path: Path) -> Optional[Dict[str, Any]]:
    try:
        with open(snapshot_path, "r", encoding="utf-8") as handle:
            return json.load(handle)
    except Exception as exc:
        warnings.warn(f"无法解析快照文件 {snapshot_path}: {exc}")
        return None


def _resolve_snapshot_path(entry: Any, base_dir: Path) -> Optional[Path]:
    if not isinstance(entry, dict):
        return None
    candidates: List[Path] = []
    relative = entry.get("relative")
    if isinstance(relative, str) and relative:
        candidate = base_dir if relative == "." else base_dir / relative
        candidates.append(candidate)
    absolute = entry.get("absolute")
    if isinstance(absolute, str) and absolute:
        candidates.append(Path(absolute))
    for candidate in candidates:
        try:
            resolved = candidate.resolve(strict=False)
        except Exception:
            resolved = candidate
        if resolved.exists():
            return resolved
    if candidates:
        try:
            return candidates[0].resolve(strict=False)
        except Exception:
            return candidates[0]
    return None


def _append_unique(paths: List[Path], candidate: Path) -> None:
    resolved = candidate.resolve(strict=False)
    if all(existing != resolved for existing in paths):
        paths.append(resolved)


def _find_metadata_paths(
    input_files: Sequence[str],
) -> Tuple[Optional[Path], List[Path], List[Path], Dict[int, int]]:
    seen_dirs: Set[Path] = set()
    dir_queue: List[Path] = []

    def enqueue(path: Optional[Path]) -> None:
        if path is None:
            return
        try:
            resolved = path.resolve(strict=False)
        except Exception:
            resolved = path
        if resolved in seen_dirs:
            return
        seen_dirs.add(resolved)
        dir_queue.append(resolved)

    for path_str in input_files:
        path = Path(path_str).resolve(strict=False)
        enqueue(path.parent)
        parent_parent = path.parent.parent
        if parent_parent != path.parent:
            enqueue(parent_parent)

    cwd = Path.cwd().resolve(strict=False)
    enqueue(cwd)
    enqueue(cwd / "genome_output")

    tax_info_path: Optional[Path] = None
    assembly_paths: List[Path] = []
    target_paths: List[Path] = []
    length_overrides: Dict[int, int] = {}
    processed_snapshots: Set[Path] = set()

    index = 0
    while index < len(dir_queue):
        directory = dir_queue[index]
        index += 1
        if not directory.exists():
            continue
        if directory.is_file():
            enqueue(directory.parent)
            continue

        tax_candidate = directory / "tax.info"
        if tax_candidate.exists() and tax_info_path is None:
            tax_info_path = tax_candidate.resolve(strict=False)
        assembly_candidate = directory / "assembly_summary.txt"
        if assembly_candidate.exists():
            _append_unique(assembly_paths, assembly_candidate)
        target_candidate = directory / "target.tsv"
        if target_candidate.exists():
            _append_unique(target_paths, target_candidate)

        build_dir = directory / "build"
        if build_dir.exists() and build_dir.is_dir():
            enqueue(build_dir)

        for snapshot in directory.glob("*.profile.json"):
            snapshot_resolved = snapshot.resolve(strict=False)
            if snapshot_resolved in processed_snapshots:
                continue
            processed_snapshots.add(snapshot_resolved)
            data = _load_snapshot(snapshot_resolved)
            if not data:
                continue
            snapshot_dir = snapshot_resolved.parent
            enqueue(snapshot_dir)

            db_info = data.get("database")
            dataset_root = None
            if isinstance(db_info, dict):
                dataset_entry = db_info.get("dataset_root")
                dataset_root = _resolve_snapshot_path(dataset_entry, snapshot_dir) or snapshot_dir
            else:
                dataset_root = snapshot_dir
            enqueue(dataset_root)

            if isinstance(db_info, dict):
                target_entry = db_info.get("target")
                if target_entry:
                    target_resolved = _resolve_snapshot_path(target_entry, dataset_root)
                    if target_resolved and target_resolved.exists():
                        _append_unique(target_paths, target_resolved)
                tax_entry = db_info.get("tax_info")
                if tax_entry:
                    tax_resolved = _resolve_snapshot_path(tax_entry, dataset_root)
                    if tax_resolved and tax_resolved.exists() and tax_info_path is None:
                        tax_info_path = tax_resolved
                assembly_entry = db_info.get("assembly_summary")
                if assembly_entry:
                    assembly_resolved = _resolve_snapshot_path(assembly_entry, dataset_root)
                    if assembly_resolved and assembly_resolved.exists():
                        _append_unique(assembly_paths, assembly_resolved)

            lengths = data.get("taxid_lengths")
            if isinstance(lengths, dict):
                for key, value in lengths.items():
                    try:
                        taxid = int(key)
                        length_val = int(value)
                    except (TypeError, ValueError):
                        continue
                    if length_val > 0:
                        length_overrides[taxid] = length_val

    return tax_info_path, assembly_paths, target_paths, length_overrides


def _load_target_accessions(target_paths: Sequence[Path]) -> Dict[int, List[str]]:
    mapping: Dict[int, List[str]] = defaultdict(list)
    for path in target_paths:
        with open(path, "r", encoding="utf-8", errors="ignore") as infile:
            for line in infile:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                try:
                    taxid = int(parts[1])
                except ValueError:
                    continue
                accession = Path(parts[0]).name
                if accession.endswith(".gz"):
                    accession = accession[:-3]
                accession = accession.split(".fna", 1)[0]
                if accession.endswith(".fa") or accession.endswith(".fasta"):
                    accession = accession.rsplit(".", 1)[0]
                token_parts = accession.split("_")
                if len(token_parts) >= 2:
                    accession = "_".join(token_parts[:2])
                if accession not in mapping[taxid]:
                    mapping[taxid].append(accession)
    return mapping


def _load_genome_lengths(assembly_paths: Sequence[Path]) -> Tuple[Dict[int, int], Dict[str, int]]:
    length_by_taxid: Dict[int, int] = defaultdict(int)
    length_by_accession: Dict[str, int] = {}
    seen_accessions = set()
    for path in assembly_paths:
        with open(path, "r", encoding="utf-8", errors="ignore") as infile:
            reader = csv.reader(infile, delimiter="\t")
            for row in reader:
                if not row or row[0].startswith("#"):
                    continue
                if len(row) <= 25:
                    continue
                accession = row[0]
                if accession in seen_accessions:
                    continue
                seen_accessions.add(accession)
                try:
                    taxid = int(row[5])
                    length = int(float(row[25]))
                except (ValueError, IndexError):
                    continue
                if length <= 0:
                    continue
                length_by_taxid[taxid] += length
                length_by_accession[accession] = length
    return length_by_taxid, length_by_accession


def _resolve_length(
    taxid: int,
    tax_lengths: Dict[int, int],
    accession_map: Dict[int, List[str]],
    accession_lengths: Dict[str, int],
) -> Optional[int]:
    if accession_map:
        accessions = accession_map.get(taxid)
        if accessions:
            total = 0
            for acc in accessions:
                length = accession_lengths.get(acc)
                if length:
                    total += length
            if total > 0:
                return total
    length = tax_lengths.get(taxid)
    if length and length > 0:
        return length
    return None


def _format_count(value: float) -> str:
    rounded = round(value)
    if abs(value - rounded) < 1e-6:
        return str(int(rounded))
    return f"{value:.6f}"


def _aggregate_levels(
    resolver: TaxonomyResolver,
    taxid_weights: Dict[int, float],
    base_unclassified: float,
    tax_lengths: Dict[int, int],
    accession_map: Dict[int, List[str]],
    accession_lengths: Dict[str, int],
) -> Tuple[
    Dict[str, Dict[str, float]],
    Dict[str, Dict[str, float]],
    Dict[str, Dict[str, float]],
    Dict[int, TaxonRecord],
]:
    count_by_level: Dict[str, Dict[str, float]] = {
        level: defaultdict(float) for level in LEVEL_ALIASES
    }
    length_by_level: Dict[str, Dict[str, float]] = {
        level: defaultdict(float) for level in LEVEL_ALIASES
    }
    rpk_by_level: Dict[str, Dict[str, float]] = {
        level: defaultdict(float) for level in LEVEL_ALIASES
    }
    taxon_records: Dict[int, TaxonRecord] = {}

    if base_unclassified:
        for level in count_by_level:
            count_by_level[level][UNCLASSIFIED] += base_unclassified

    for taxid, count in taxid_weights.items():
        if count <= 0.0:
            continue
        lineage = resolver.get_lineage(taxid)
        if not lineage:
            for level in count_by_level:
                count_by_level[level][UNCLASSIFIED] += count
            continue

        has_level = {level: False for level in LEVEL_ALIASES}
        hits: Dict[str, List[Tuple[int, str]]] = {level: [] for level in LEVEL_ALIASES}

        for record in lineage:
            normalized = _normalize_rank(record.rank)
            if not normalized:
                continue
            if normalized != "clade" and has_level[normalized]:
                continue
            name = record.name or str(record.taxid)
            hits[normalized].append((record.taxid, name))
            has_level[normalized] = True

        leaf = next((rec for rec in reversed(lineage) if rec.taxid == taxid), lineage[-1])
        taxon_records[taxid] = leaf
        leaf_length = _resolve_length(taxid, tax_lengths, accession_map, accession_lengths)
        leaf_rpk = (count / (leaf_length / 1000.0)) if leaf_length and leaf_length > 0 else 0.0

        for level, nodes in hits.items():
            if not nodes:
                continue
            for node_taxid, name in nodes:
                count_by_level[level][name] += count
                if leaf_length:
                    length_by_level[level][name] += leaf_length
                if leaf_rpk > 0.0:
                    rpk_by_level[level][name] += leaf_rpk

        for level in LEVEL_ALIASES:
            if hits[level]:
                continue
            if _should_mark_unclassified(level, has_level):
                count_by_level[level][UNCLASSIFIED] += count

    return count_by_level, length_by_level, rpk_by_level, taxon_records


def process_file(input_files: Iterable[str], output_file: str) -> None:
    inputs = list(input_files)
    (
        taxid_weights,
        base_unclassified,
        total_lines,
        per_read_candidates,
    ) = _collect_taxid_weights(inputs, keep_per_read=True)
    global_em_applied = False

    tax_info_path, assembly_paths, target_paths, snapshot_lengths = _find_metadata_paths(inputs)
    resolver = TaxonomyResolver(tax_info_path)
    if not resolver.has_taxonomy():
        warnings.warn("无法解析 taxonomy，所有记录将视为未分类", RuntimeWarning)

    accession_map = _load_target_accessions(target_paths) if target_paths else {}
    tax_lengths, accession_lengths = _load_genome_lengths(assembly_paths) if assembly_paths else ({}, {})
    if snapshot_lengths:
        for taxid, length in snapshot_lengths.items():
            if length > 0:
                tax_lengths[taxid] = length
    if accession_map and not accession_lengths:
        warnings.warn("未找到 assembly_summary.txt，TPM 将无法计算")

    if per_read_candidates:
        if any(len(candidates) > 1 for candidates in per_read_candidates):
            taxid_weights = _run_global_em(
                per_read_candidates,
                taxid_weights,
                tax_lengths,
                resolver=resolver,
            )
            global_em_applied = True
        per_read_candidates = None

    weighted_total = sum(taxid_weights.values()) + base_unclassified
    percentage_total = float(total_lines) if total_lines > 0 else weighted_total
    rpm_total = total_lines if total_lines > 0 else percentage_total

    count_by_level, length_by_level, rpk_by_level, taxon_records = _aggregate_levels(
        resolver,
        taxid_weights,
        base_unclassified,
        tax_lengths,
        accession_map,
        accession_lengths,
    )

    level_totals = {level: sum(counter.values()) for level, counter in count_by_level.items()}

    taxa_rows: List[Tuple[int, str, str, float, float, float, Optional[int], float]] = []
    for taxid, count in taxid_weights.items():
        record = taxon_records.get(taxid)
        name = record.name if record and record.name else str(taxid)
        rank = record.rank if record else "na"
        length = _resolve_length(taxid, tax_lengths, accession_map, accession_lengths)
        relative = (count / percentage_total * 100) if percentage_total > 0 else 0.0
        rpm = (count / rpm_total * 1_000_000) if rpm_total > 0 else 0.0
        rpk = (count / (length / 1000)) if length and length > 0 else 0.0
        taxa_rows.append((taxid, name, rank, count, relative, rpm, length, rpk))

    total_rpk = sum(row[7] for row in taxa_rows if row[7] > 0)

    taxa_output_path = Path(output_file + ".taxa.tsv")
    with open(taxa_output_path, "w", encoding="utf-8") as taxa_file:
        taxa_file.write(
            "TaxID\tName\tRank\tAssigned Reads\tRelative Abundance (%)\tReads per Million\t"
            "Genome Length (bp)\tRPK\tTPM\n"
        )
        for taxid, name, rank, count, relative, rpm, length, rpk in sorted(
            taxa_rows, key=lambda item: item[3], reverse=True
        ):
            tpm = (rpk / total_rpk * 1_000_000) if total_rpk > 0 and rpk > 0 else 0.0
            length_str = str(length) if length else ""
            taxa_file.write(
                f"{taxid}\t{name}\t{rank}\t{_format_count(count)}\t{relative:.2f}\t{rpm:.2f}\t"
                f"{length_str}\t{rpk:.6f}\t{tpm:.2f}\n"
            )

        if base_unclassified:
            relative = (base_unclassified / percentage_total * 100) if percentage_total > 0 else 0.0
            rpm = (base_unclassified / rpm_total * 1_000_000) if rpm_total > 0 else 0.0
            taxa_file.write(
                f"{UNCLASSIFIED}\t{UNCLASSIFIED}\tna\t{_format_count(base_unclassified)}\t{relative:.2f}\t"
                f"{rpm:.2f}\t\t\t\n"
            )

    output_path = Path(output_file + ".tsv")
    with open(output_path, "w", encoding="utf-8") as outfile:
        outfile.write(
            f"# raw_reads={total_lines} weighted_total={weighted_total:.6f} normalization=tpm "
            f"global_em={'on' if global_em_applied else 'off'} "
            f"length_norm={'on' if GLOBAL_LENGTH_NORMALIZATION else 'off'}\n"
        )
        outfile.write("# note: TPM is computed within each table independently\n")
        outfile.write(
            "Taxon\tReads\tRelative Abundance (%)\tReads per Million\tGenome Length (bp)\tRPK\tTPM\t"
            "Shannon Index\tSimpson Index\n"
        )

        for level in LEVEL_ALIASES:
            taxon_counter = count_by_level[level]
            total_count = level_totals.get(level, 0.0)
            if total_count <= 0:
                continue
            shannon_index = calculate_shannon_index(taxon_counter)
            simpson_index = calculate_simpson_index(taxon_counter)
            display_level = " ".join(part.capitalize() for part in level.split())
            outfile.write(
                f"\n## {display_level} Level ##\tShannon={shannon_index:.4f}\tSimpson={simpson_index:.4f}\n"
            )

            items: List[Tuple[str, float]] = []
            unclassified_item: Optional[Tuple[str, float]] = None
            for taxon, count in taxon_counter.items():
                if taxon == UNCLASSIFIED:
                    unclassified_item = (taxon, count)
                else:
                    items.append((taxon, count))
            items.sort(key=lambda item: item[1], reverse=True)
            if unclassified_item:
                items.append(unclassified_item)

            rpk_counter = rpk_by_level[level]
            total_rpk_level = sum(value for value in rpk_counter.values() if value > 0.0)

            for taxon, count in items:
                relative_level = (count / total_count * 100) if total_count > 0 else 0.0
                rpm = (count / rpm_total * 1_000_000) if rpm_total > 0 else 0.0
                length_val = length_by_level[level].get(taxon)
                rpk = rpk_counter.get(taxon, 0.0)
                tpm = (rpk / total_rpk_level * 1_000_000) if total_rpk_level > 0 and rpk > 0 else 0.0
                length_str = str(int(length_val)) if length_val and length_val > 0 else ""
                outfile.write(
                    f"{taxon}\t{_format_count(count)}\t{relative_level:.2f}\t{rpm:.2f}\t{length_str}\t"
                    f"{rpk:.6f}\t{tpm:.2f}\t{shannon_index:.4f}\t{simpson_index:.4f}\n"
                )

    if taxid_weights and total_lines and weighted_total < total_lines:
        warnings.warn(
            "部分 reads 未被计入丰度统计，可能因权重缺失或被标记为未分类。",
            RuntimeWarning,
        )
