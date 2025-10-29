import math
import warnings
from collections import Counter
from typing import Dict, Iterable, List, Optional, Tuple

from ete3 import NCBITaxa

from .taxonomy_utils import (
    GtdbTaxonomy,
    load_gtdb_taxonomy,
    normalize_kind,
    read_taxonomy_meta,
)

# 目标层级及其在 NCBI taxonomy 中可能出现的 rank 名称
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


def _normalize_rank(rank: Optional[str]) -> Optional[str]:
    if rank is None:
        return None
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


def calculate_shannon_index(taxon_dict: Dict[str, int]) -> float:
    total_count = sum(taxon_dict.values())
    if total_count == 0:
        return 0.0
    shannon_index = 0.0
    for count in taxon_dict.values():
        p_i = count / total_count
        if p_i > 0:
            shannon_index -= p_i * math.log(p_i)
    return shannon_index


def calculate_simpson_index(taxon_dict: Dict[str, int]) -> float:
    total_count = sum(taxon_dict.values())
    if total_count == 0:
        return 0.0
    simpson_index = 0.0
    for count in taxon_dict.values():
        p_i = count / total_count
        simpson_index += p_i ** 2
    return 1.0 - simpson_index


def _parse_primary_taxid(line: str, expect_numeric: bool) -> Optional[str]:
    line = line.strip()
    if not line:
        return None
    parts = [token.strip() for token in line.split("\t") if token.strip()]
    if len(parts) < 2:
        return None
    for token in parts[1:]:
        if not token or token.startswith("POST_TOP2="):
            continue
        if token.lower() == UNCLASSIFIED:
            return None
        tid_token = token.split(":", 1)[0]
        if expect_numeric:
            try:
                tid = int(tid_token)
            except ValueError:
                continue
            if tid > 0:
                return str(tid)
        else:
            if tid_token:
                return tid_token
    return None


def _collect_taxids(
    input_files: Iterable[str], expect_numeric: bool
) -> Tuple[Counter[str], int]:
    taxid_counts: Counter[str] = Counter()
    unclassified_reads = 0
    for input_file in input_files:
        with open(input_file, "r", errors="ignore") as infile:
            for line in infile:
                taxid = _parse_primary_taxid(line, expect_numeric=expect_numeric)
                if taxid is None:
                    unclassified_reads += 1
                else:
                    taxid_counts[taxid] += 1
    return taxid_counts, unclassified_reads


def _aggregate_levels(
    tax: Optional[NCBITaxa],
    taxid_counts: Counter[int],
    base_unclassified: int,
) -> Dict[str, Counter]:
    count_by_level: Dict[str, Counter[str]] = {
        level: Counter() for level in LEVEL_ALIASES
    }
    if base_unclassified:
        for level in count_by_level:
            count_by_level[level][UNCLASSIFIED] += base_unclassified
    if tax is None:
        return count_by_level

    for taxid, count in taxid_counts.items():
        try:
            lineage = tax.get_lineage(taxid)
            rank_map = tax.get_rank(lineage)
            name_map = tax.get_taxid_translator(lineage)
        except Exception:
            for level in count_by_level:
                count_by_level[level][UNCLASSIFIED] += count
            continue

        has_level = {level: False for level in LEVEL_ALIASES}
        hits: Dict[str, List[str]] = {level: [] for level in LEVEL_ALIASES}

        for node in lineage:
            rank = rank_map.get(node)
            if not rank or rank == "no rank":
                continue
            normalized = _normalize_rank(rank)
            if not normalized:
                continue
            name = name_map.get(node)
            if not name:
                continue
            if normalized != "clade" and has_level[normalized]:
                continue
            hits[normalized].append(name)
            has_level[normalized] = True

        for level, names in hits.items():
            if not names:
                continue
            for name in names:
                count_by_level[level][name] += count

        for level in LEVEL_ALIASES:
            if hits[level]:
                continue
            if _should_mark_unclassified(level, has_level):
                count_by_level[level][UNCLASSIFIED] += count

    return count_by_level


def _aggregate_gtdb_levels(
    taxonomy: GtdbTaxonomy,
    taxid_counts: Counter[str],
    base_unclassified: int,
) -> Dict[str, Counter[str]]:
    count_by_level: Dict[str, Counter[str]] = {
        level: Counter() for level in LEVEL_ALIASES
    }
    if base_unclassified:
        for level in count_by_level:
            count_by_level[level][UNCLASSIFIED] += base_unclassified

    for node, count in taxid_counts.items():
        if node not in taxonomy.rank:
            for level in count_by_level:
                count_by_level[level][UNCLASSIFIED] += count
            continue
        has_level = {level: False for level in LEVEL_ALIASES}
        recorded_any = False
        for current in taxonomy.iter_lineage(node):
            node_rank = taxonomy.rank.get(current)
            normalized = _normalize_rank(node_rank)
            if not normalized:
                continue
            if normalized != "clade" and has_level[normalized]:
                continue
            display_name = taxonomy.display_name(current)
            count_by_level[normalized][display_name] += count
            has_level[normalized] = True
            recorded_any = True
        if not recorded_any:
            for level in count_by_level:
                count_by_level[level][UNCLASSIFIED] += count
            continue
        for level in LEVEL_ALIASES:
            if level == "clade":
                continue
            if not has_level[level] and _should_mark_unclassified(level, has_level):
                count_by_level[level][UNCLASSIFIED] += count
    return count_by_level


def process_file(
    input_files: Iterable[str],
    output_file: str,
    taxonomy_kind: str = "auto",
    taxonomy_version: str = "auto",
    taxonomy_info: Optional[str] = None,
    taxonomy_meta: Optional[str] = None,
) -> None:
    meta = read_taxonomy_meta(taxonomy_meta)
    resolved_kind = normalize_kind(taxonomy_kind)
    if resolved_kind == "auto":
        resolved_kind = normalize_kind(meta.kind)
    if resolved_kind == "auto":
        resolved_kind = "ncbi"
    resolved_version = taxonomy_version or "auto"
    if resolved_version in {"auto", "", None}:
        resolved_version = meta.version

    expect_numeric = resolved_kind != "gtdb"
    taxid_counts, base_unclassified = _collect_taxids(
        input_files, expect_numeric=expect_numeric
    )

    if resolved_kind == "gtdb":
        if not taxonomy_info:
            raise ValueError("GTDB 模式需要提供 --taxonomy-info (tax.info 文件路径)")
        gtdb_taxonomy = load_gtdb_taxonomy(taxonomy_info)
        count_by_level = _aggregate_gtdb_levels(
            gtdb_taxonomy, taxid_counts, base_unclassified
        )
    else:
        try:
            tax = NCBITaxa()
        except Exception as exc:  # pragma: no cover - 初始化失败极少发生
            warnings.warn(f"无法初始化 NCBITaxa，所有条目将标记为未分类: {exc}")
            tax = None

        numeric_counts: Counter[int] = Counter()
        for taxid_str, count in taxid_counts.items():
            try:
                tid = int(taxid_str)
            except ValueError:
                continue
            if tid > 0:
                numeric_counts[tid] += count
        count_by_level = _aggregate_levels(tax, numeric_counts, base_unclassified)

    total_counts_by_level = {
        level: sum(counter.values()) for level, counter in count_by_level.items()
    }

    with open(output_file + ".tsv", "w") as outfile:
        outfile.write("Taxon\tCount\tRelative Abundance (%)\tShannon Index\tSimpson Index\n")

        for level, taxon_counter in count_by_level.items():
            total_count = total_counts_by_level[level]
            if total_count == 0:
                continue

            shannon_index = calculate_shannon_index(taxon_counter)
            simpson_index = calculate_simpson_index(taxon_counter)
            display_level = " ".join(part.capitalize() for part in level.split())
            outfile.write(f"\n## {display_level} Level ##\n")

            items = []
            unclassified_item = None
            for taxon, count in taxon_counter.items():
                if taxon == UNCLASSIFIED:
                    unclassified_item = (taxon, count)
                else:
                    items.append((taxon, count))
            items.sort(key=lambda item: item[1], reverse=True)
            if unclassified_item:
                items.append(unclassified_item)

            for taxon, count in items:
                relative_abundance = (count / total_count) * 100 if total_count > 0 else 0.0
                outfile.write(
                    f"{taxon}\t{count}\t{relative_abundance:.2f}\t{shannon_index:.4f}\t{simpson_index:.4f}\n"
                )
