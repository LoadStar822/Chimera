import math
import warnings
from collections import Counter
from typing import Dict, Iterable, List, Optional, Tuple

from ete3 import NCBITaxa

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


def _parse_primary_taxid(line: str) -> Optional[int]:
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
        try:
            tid = int(tid_token)
        except ValueError:
            continue
        if tid > 0:
            return tid
    return None


def _collect_taxids(input_files: Iterable[str]) -> Tuple[Counter[int], int]:
    taxid_counts: Counter[int] = Counter()
    unclassified_reads = 0
    for input_file in input_files:
        with open(input_file, "r", errors="ignore") as infile:
            for line in infile:
                taxid = _parse_primary_taxid(line)
                if taxid is None:
                    unclassified_reads += 1
                else:
                    taxid_counts[taxid] += 1
    return taxid_counts, unclassified_reads


def _aggregate_levels(
    tax: Optional[NCBITaxa],
    taxid_counts: Counter,
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


def process_file(input_files: Iterable[str], output_file: str) -> None:
    taxid_counts, base_unclassified = _collect_taxids(input_files)
    try:
        tax = NCBITaxa()
    except Exception as exc:  # pragma: no cover - 初始化失败极少发生
        warnings.warn(f"无法初始化 NCBITaxa，所有条目将标记为未分类: {exc}")
        tax = None

    count_by_level = _aggregate_levels(tax, taxid_counts, base_unclassified)
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
