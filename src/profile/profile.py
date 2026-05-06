import warnings
from collections import Counter
from pathlib import Path
from typing import Callable, Dict, Iterable, Optional, Tuple

from ete3 import NCBITaxa

from .taxonomy_utils import (
    GtdbTaxonomy,
    load_gtdb_taxonomy,
    normalize_kind,
    read_taxonomy_meta,
)

UNCLASSIFIED = "unclassified"


def _ncbi_taxid_to_species(
    tax: Optional[NCBITaxa], taxid: int, cache: Dict[int, Optional[str]]
) -> Optional[str]:
    if tax is None:
        return None
    if taxid in cache:
        return cache[taxid]
    try:
        lineage = tax.get_lineage(taxid)
        rank_map = tax.get_rank(lineage)
        name_map = tax.get_taxid_translator(lineage)
    except Exception:
        cache[taxid] = None
        return None

    species_name: Optional[str] = None
    for node in lineage:
        if rank_map.get(node) == "species" and name_map.get(node):
            species_name = name_map[node]
            break
    cache[taxid] = species_name
    return species_name


def _gtdb_node_to_species(
    taxonomy: GtdbTaxonomy, node: str, cache: Dict[str, Optional[str]]
) -> Optional[str]:
    if node in cache:
        return cache[node]
    if node not in taxonomy.rank:
        cache[node] = None
        return None
    species_name: Optional[str] = None
    for current in taxonomy.iter_lineage(node):
        if taxonomy.rank.get(current) == "species":
            species_name = taxonomy.display_name(current)
            break
    cache[node] = species_name
    return species_name


def _raise_invalid_profile_input(input_file: str) -> None:
    raise ValueError(
        "Profile input must be ChimeraEvidence.tsv or an aggregate table with "
        "a 'taxid\\tcount' header. ChimeraClassify.tsv is read-level output "
        f"and is not a valid abundance input: {input_file}"
    )


def _collect_taxids(
    input_files: Iterable[str], expect_numeric: bool
) -> Tuple[Counter[str], float]:
    taxid_counts: Counter[str] = Counter()
    unclassified_reads = 0.0
    for input_file in input_files:
        with open(input_file, "r", errors="ignore") as infile:
            saw_header = False
            aggregate_columns: Dict[str, int] = {}
            for raw in infile:
                line = raw.rstrip("\n")
                if not line:
                    continue
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                if not saw_header:
                    header = [part.strip().lower() for part in parts]
                    if len(header) < 2 or header[0] != "taxid" or header[1] != "count":
                        _raise_invalid_profile_input(input_file)
                    aggregate_columns = {
                        part.strip().lower(): idx for idx, part in enumerate(parts)
                    }
                    saw_header = True
                    continue
                tid_token = parts[0].strip()
                value_idx = aggregate_columns.get("count", 1)
                if value_idx >= len(parts):
                    continue
                try:
                    count = float(parts[value_idx])
                except ValueError:
                    continue
                if count <= 0:
                    continue
                if tid_token.lower() == UNCLASSIFIED:
                    unclassified_reads += count
                    continue
                if expect_numeric and not tid_token.isdigit():
                    _raise_invalid_profile_input(input_file)
                taxid_counts[tid_token] += count
            if not saw_header:
                _raise_invalid_profile_input(input_file)
    return taxid_counts, unclassified_reads


def _write_species_profile(
    output_file: str,
    count_unit: str,
    taxon_counter: Counter[str],
) -> None:
    output_base = output_file
    if output_base.lower().endswith(".tsv"):
        output_base = output_base[: -len(".tsv")]
    total_count = sum(float(v) for v in taxon_counter.values())
    output_path = Path(output_base + ".tsv")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as outfile:
        outfile.write("# profile_input=aggregate_evidence\n")
        outfile.write(f"# count_unit={count_unit}\n")
        outfile.write("Taxon\tCount\tRelative Abundance (%)\n")
        if not taxon_counter or total_count == 0:
            return
        outfile.write("\n## Species Level ##\n")
        items = [
            (taxon, float(count))
            for taxon, count in taxon_counter.items()
            if taxon != UNCLASSIFIED
        ]
        items.sort(key=lambda item: item[1], reverse=True)
        if UNCLASSIFIED in taxon_counter:
            items.append((UNCLASSIFIED, float(taxon_counter[UNCLASSIFIED])))
        for taxon, count in items:
            relative_abundance = (
                (float(count) / float(total_count)) * 100.0
                if total_count > 0
                else 0.0
            )
            outfile.write(f"{taxon}\t{count}\t{relative_abundance:.2f}\n")


def _aggregate_species_counts(
    taxon_counts: Counter,
    base_unclassified: float,
    resolve_species: Callable[[object], Optional[str]],
) -> Dict[str, Counter[str]]:
    species: Counter[str] = Counter()
    if base_unclassified:
        species[UNCLASSIFIED] += base_unclassified
    for taxon_id, count in taxon_counts.items():
        sp = resolve_species(taxon_id)
        if sp:
            species[sp] += count
        else:
            species[UNCLASSIFIED] += count
    return {"species": species}


def _aggregate_levels(
    tax: Optional[NCBITaxa],
    taxid_counts: Counter[int],
    base_unclassified: float,
) -> Dict[str, Counter[str]]:
    cache: Dict[int, Optional[str]] = {}
    return _aggregate_species_counts(
        taxid_counts,
        base_unclassified,
        lambda taxid: _ncbi_taxid_to_species(tax, int(taxid), cache),
    )


def _aggregate_gtdb_levels(
    taxonomy: GtdbTaxonomy,
    taxid_counts: Counter[str],
    base_unclassified: float,
) -> Dict[str, Counter[str]]:
    cache: Dict[str, Optional[str]] = {}

    def resolve_species(node: object) -> Optional[str]:
        node = str(node)
        if node not in taxonomy.rank:
            return None
        return _gtdb_node_to_species(taxonomy, node, cache)

    return _aggregate_species_counts(
        taxid_counts, base_unclassified, resolve_species
    )


def process_file(
    input_files: Iterable[str],
    output_file: str,
    taxonomy_kind: str = "auto",
    taxonomy_info: Optional[str] = None,
    taxonomy_meta: Optional[str] = None,
) -> None:
    input_files = list(input_files)
    meta = read_taxonomy_meta(taxonomy_meta)
    resolved_kind = normalize_kind(taxonomy_kind)
    if resolved_kind == "auto":
        resolved_kind = normalize_kind(meta.kind)
    if resolved_kind == "auto":
        resolved_kind = "ncbi"

    expect_numeric = resolved_kind != "gtdb"

    taxid_counts, base_unclassified = _collect_taxids(
        input_files, expect_numeric=expect_numeric
    )

    tax: Optional[NCBITaxa] = None
    gtdb_taxonomy: Optional[GtdbTaxonomy] = None
    if resolved_kind == "gtdb":
        if not taxonomy_info:
            raise ValueError("GTDB mode requires --taxonomy-info (path to tax.info)")
        gtdb_taxonomy = load_gtdb_taxonomy(taxonomy_info)
    else:
        try:
            tax = NCBITaxa()
        except Exception as exc:  # pragma: no cover
            warnings.warn(
                "Failed to initialize NCBITaxa; all entries will be marked as "
                f"unclassified: {exc}"
            )
            tax = None

    def aggregate_levels_for_counts(
        counts: Counter[str], unclassified_reads: int
    ) -> Dict[str, Counter[str]]:
        if resolved_kind == "gtdb":
            if gtdb_taxonomy is None:  # pragma: no cover
                raise ValueError("GTDB taxonomy not loaded")
            return _aggregate_gtdb_levels(gtdb_taxonomy, counts, unclassified_reads)
        numeric_counts: Counter[int] = Counter()
        for taxid_str, count in counts.items():
            try:
                tid = int(taxid_str)
            except ValueError:
                continue
            if tid > 0:
                numeric_counts[tid] += count
        return _aggregate_levels(tax, numeric_counts, unclassified_reads)

    count_by_level = aggregate_levels_for_counts(taxid_counts, base_unclassified)
    count_unit = "posterior_evidence"
    _write_species_profile(
        output_file=output_file,
        count_unit=count_unit,
        taxon_counter=count_by_level.get("species") or Counter(),
    )
