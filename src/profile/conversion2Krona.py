from collections import Counter
from typing import Iterable, Optional

from multitax import NcbiTx

from .taxonomy_utils import (
    GtdbTaxonomy,
    load_gtdb_taxonomy,
    normalize_kind,
    read_taxonomy_meta,
)

LEVEL_ORDER = [
    "acellular root",
    "realm",
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "strain",
]
UNCLASSIFIED = "unclassified"


def _normalize_rank(rank: Optional[str]) -> Optional[str]:
    if rank is None:
        return None
    rank_lower = rank.lower()
    for level in LEVEL_ORDER:
        if rank_lower == level:
            return level
    # 兼容 NCBI 的 superkingdom 等命名
    if rank_lower == "superkingdom":
        return "domain"
    return None


def _parse_taxid_from_line(line: str, expect_numeric: bool) -> Optional[str]:
    parts = [token.strip() for token in line.strip().split("\t") if token.strip()]
    if len(parts) < 2:
        return None
    for token in parts[1:]:
        if not token or token.startswith("POST_"):
            continue
        if token.lower() == UNCLASSIFIED:
            return UNCLASSIFIED
        tid_token = token.split(":", 1)[0]
        if expect_numeric:
            if not tid_token.isdigit():
                continue
            return tid_token
        if tid_token:
            return tid_token
    return None


def _gtdb_lineage(taxonomy: GtdbTaxonomy, node: str) -> Optional[str]:
    if node not in taxonomy.rank:
        return None
    rank_hits = {}
    for current in taxonomy.iter_lineage(node):
        normalized = _normalize_rank(taxonomy.rank.get(current))
        if not normalized:
            continue
        if normalized not in rank_hits:
            rank_hits[normalized] = taxonomy.display_name(current)
    lineage = [rank_hits[level] for level in LEVEL_ORDER if level in rank_hits]
    if not lineage:
        return None
    return "\t".join(lineage)


def convert_multiple_files_to_krona_format(
    input_files: Iterable[str],
    output_file: str,
    taxonomy_kind: str = "auto",
    taxonomy_info: Optional[str] = None,
    taxonomy_meta: Optional[str] = None,
) -> None:
    meta = read_taxonomy_meta(taxonomy_meta)
    resolved_kind = normalize_kind(taxonomy_kind)
    if resolved_kind == "auto":
        resolved_kind = normalize_kind(meta.kind)
    if resolved_kind == "auto":
        resolved_kind = "ncbi"

    expect_numeric = resolved_kind != "gtdb"
    taxid_count: Counter[str] = Counter()
    for input_file in input_files:
        with open(input_file, "r", errors="ignore") as infile:
            for line in infile:
                taxid = _parse_taxid_from_line(line, expect_numeric=expect_numeric)
                if not taxid:
                    taxid = UNCLASSIFIED
                taxid_count[taxid] += 1

    if resolved_kind == "gtdb":
        if not taxonomy_info:
            raise ValueError("GTDB Krona 转换需要提供 --taxonomy-info (tax.info)")
        gtdb_taxonomy = load_gtdb_taxonomy(taxonomy_info)

    ncbi_tax = None
    if resolved_kind == "ncbi":
        try:
            ncbi_tax = NcbiTx()
        except Exception:
            ncbi_tax = None

    with open(output_file + ".tsv", "w", encoding="utf-8") as outfile:
        for taxid, count in taxid_count.items():
            if taxid == UNCLASSIFIED:
                outfile.write(f"{count + 1}\t{UNCLASSIFIED}\n")
                continue
            if resolved_kind == "gtdb":
                lineage = _gtdb_lineage(gtdb_taxonomy, taxid)
                if lineage is None:
                    outfile.write(f"{count + 1}\t{UNCLASSIFIED}\n")
                else:
                    outfile.write(f"{count + 1}\t{lineage}\n")
            else:
                if ncbi_tax is None:
                    outfile.write(f"{count + 1}\t{UNCLASSIFIED}\n")
                    continue
                try:
                    levels = ncbi_tax.name_lineage(int(taxid))
                except Exception:
                    outfile.write(f"{count + 1}\t{UNCLASSIFIED}\n")
                    continue
                taxonomy_string = "\t".join(levels)
                outfile.write(f"{count + 1}\t{taxonomy_string}\n")

    print(f"Conversion completed. The combined output file is saved as {output_file}.")
