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
# 最小 evidence：posterior_evidence 票数下限，过滤低支持的长尾假阳性
MIN_EVIDENCE = 10
# 低多样性样本（典型：MOCK）更激进的 species 截断：减少长尾假阳性。
# 低多样性判定使用 top1 物种分布：头部质量高且有效物种数小。
LOW_DIVERSITY_TOP_K = 10
LOW_DIVERSITY_TOP_MASS = 0.75
LOW_DIVERSITY_EFF_SPECIES_MAX = 32.0
# 低多样性样本：在 species 层进一步用 head mass 裁剪极小尾巴，减少 presence FP。
# 该裁剪在“已被判定为 low-div”后执行，并以 eff_total（排除 unclassified）为分母。
LOW_DIVERSITY_SPECIES_HEAD_MASS = 99.85
LOW_DIVERSITY_SPECIES_HEAD_MIN_KEEP = 8
# 高多样性样本：最终写出时的 species 相对丰度硬阈值（%）。
HIGH_DIVERSITY_FINAL_REL_CUTOFF = 0.02
# 高多样性样本：散射型 top1 FP 往往是 singleton/doubleton。
# 用 top1 的“出现次数”做第二信号：只要 support 足够（>=K）就保留，
# 以便在提高 head mass（救 recall）时仍能压住 FP 物种种类数。
HIGH_DIVERSITY_SPECIES_SUPPORT_MIN = 3
# 高多样性样本：genus-local split（两遍扫描，基于 POST_TOPK 约束属内物种集合）。
HIGH_DIVERSITY_GENUS_LOCAL_SPLIT = True
HIGH_DIVERSITY_GENUS_LOCAL_THRES = 0.70
HIGH_DIVERSITY_GENUS_LOCAL_ALPHA = 0.90
HIGH_DIVERSITY_GENUS_LOCAL_TOPK = 3
HIGH_DIVERSITY_GENUS_LOCAL_POST_TOPK = 16


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


def _infer_genus(species_name: str) -> Optional[str]:
    """Infer genus name from a species display string.

    Intended for output-only post-processing in high-diversity species profiles.
    If parsing fails, return None (caller should skip genus contraction).
    """
    if not species_name:
        return None
    name = species_name.strip()
    if not name:
        return None
    if name.lower() == UNCLASSIFIED:
        return None

    # Strip common GTDB-style prefixes if present (keep it conservative).
    for prefix in ("s__", "g__"):
        if name.startswith(prefix):
            name = name[len(prefix) :].strip()
            break
    if not name:
        return None

    tokens = [t for t in name.split() if t]
    if not tokens:
        return None

    # "Candidatus X ..." => genus is the second token.
    if tokens[0].lower() == "candidatus":
        if len(tokens) < 2:
            return None
        genus = tokens[1]
    else:
        genus = tokens[0]

    for prefix in ("s__", "g__"):
        if genus.startswith(prefix):
            genus = genus[len(prefix) :]
            break
    genus = genus.strip()
    return genus or None


def _is_lowdiv_rescue_species_name_ok(species_name: str) -> bool:
    """Heuristic filter for low-diversity POST_TOPK rescue candidates.

    Goal: avoid adding generic/low-quality nodes like "X bacterium" or "Y sp."
    that often act as FP sinks in low-div mixtures.
    """
    if not species_name:
        return False
    name = species_name.strip()
    if not name:
        return False
    if name.lower() == UNCLASSIFIED:
        return False

    for prefix in ("s__", "g__"):
        if name.startswith(prefix):
            name = name[len(prefix) :].strip()
            break
    lower = name.lower()

    if "unclassified" in lower or "uncultured" in lower or "metagenome" in lower:
        return False
    # common low-quality placeholders
    if " bacterium" in lower:
        return False
    if " sp." in lower or lower.endswith(" sp") or " sp " in lower:
        return False
    return True


def _ncbi_species_name_for_taxid(
    tax: NCBITaxa, taxid: int, cache: Dict[int, Optional[str]]
) -> Optional[str]:
    if taxid in cache:
        return cache[taxid]
    try:
        lineage = tax.get_lineage(taxid)
        rank_map = tax.get_rank(lineage)
        name_map = tax.get_taxid_translator(lineage)
    except Exception:
        cache[taxid] = None
        return None
    name: Optional[str] = None
    for node in lineage:
        if rank_map.get(node) == "species":
            name = name_map.get(node)
            break
    if not name:
        name = name_map.get(taxid)
    cache[taxid] = name
    return name


def _gtdb_species_name_for_node(
    taxonomy: GtdbTaxonomy, node: str, cache: Dict[str, Optional[str]]
) -> Optional[str]:
    if node in cache:
        return cache[node]
    name: Optional[str] = None
    try:
        for current in taxonomy.iter_lineage(node):
            if taxonomy.rank.get(current) == "species":
                name = taxonomy.display_name(current)
                break
        if not name:
            if node in taxonomy.rank:
                name = taxonomy.display_name(node)
    except Exception:
        name = None
    cache[node] = name
    return name


def _parse_post_topk_from_tokens(
    tokens: List[str], expect_numeric: bool, topk: int
) -> List[Tuple[str, float]]:
    for tok in tokens:
        tok = tok.strip()
        if not tok.startswith("POST_TOPK="):
            continue
        payload = tok[len("POST_TOPK=") :].strip()
        if not payload:
            return []
        out: List[Tuple[str, float]] = []
        for item in payload.split(","):
            item = item.strip()
            if not item or ":" not in item:
                continue
            tid, p = item.split(":", 1)
            tid = tid.strip()
            if expect_numeric and not tid.isdigit():
                continue
            try:
                prob = float(p)
            except ValueError:
                continue
            if prob <= 0:
                continue
            out.append((tid, prob))
            if topk > 0 and len(out) >= topk:
                break
        return out
    return []


def _has_post_topk_tokens(
    input_files: Iterable[str], max_lines: int = 2000
) -> bool:
    if max_lines <= 0:
        return False
    remaining = max_lines
    for input_file in input_files:
        try:
            with open(input_file, "r", errors="ignore") as infile:
                for raw in infile:
                    if "POST_TOPK=" in raw:
                        return True
                    remaining -= 1
                    if remaining <= 0:
                        return False
        except OSError:
            continue
    return False


def _ncbi_taxid_to_species(
    tax: NCBITaxa, taxid: int, cache: Dict[int, Optional[str]]
) -> Optional[str]:
    """Map a NCBI taxid to its species display name (scientific name).

    Returns None if lineage lookup fails or species cannot be determined.
    """
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
        rank = rank_map.get(node)
        normalized = _normalize_rank(rank)
        if normalized != "species":
            continue
        name = name_map.get(node)
        if name:
            species_name = name
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
        node_rank = taxonomy.rank.get(current)
        normalized = _normalize_rank(node_rank)
        if normalized != "species":
            continue
        species_name = taxonomy.display_name(current)
        break
    cache[node] = species_name
    return species_name


def _collect_post_topk_species_stats(
    *,
    input_files: Iterable[str],
    expect_numeric: bool,
    resolved_kind: str,
    tax: Optional[NCBITaxa],
    gtdb_taxonomy: Optional[GtdbTaxonomy],
    prob_hi: float,
    topk: int,
) -> Tuple[Counter[str], Counter[str], Dict[str, float]]:
    """Collect POST_TOPK per-species stats in a single pass.

    Returns (any_support, hi_support, prob_mass):
    - any_support[species] += 1 if the species appears with prob>0 in POST_TOPK for a read
    - hi_support[species]  += 1 if max prob for that species in the read >= prob_hi
    - prob_mass[species]   += sum(probabilities for that species in the read)

    Notes:
    - per-read dedup is performed at species level (multiple strains map to one species)
    - candidates with prob<=0 are ignored (consistent with _parse_post_topk_from_tokens)
    """
    any_support: Counter[str] = Counter()
    hi_support: Counter[str] = Counter()
    prob_mass: Dict[str, float] = {}

    ncbi_cache: Dict[int, Optional[str]] = {}
    gtdb_cache: Dict[str, Optional[str]] = {}

    for input_file in input_files:
        with open(input_file, "r", errors="ignore") as infile:
            for raw in infile:
                line = raw.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                tokens = parts[1:]
                cand = _parse_post_topk_from_tokens(
                    [t.strip() for t in tokens if t.strip()],
                    expect_numeric=expect_numeric,
                    topk=int(topk),
                )
                if not cand:
                    continue

                # per-read aggregation at species level
                read_mass: Dict[str, float] = {}
                read_max: Dict[str, float] = {}
                for tid, p in cand:
                    sp: Optional[str] = None
                    if resolved_kind == "gtdb":
                        if gtdb_taxonomy is None:
                            continue
                        sp = _gtdb_node_to_species(gtdb_taxonomy, tid, gtdb_cache)
                    else:
                        if tax is None:
                            continue
                        if not tid.isdigit():
                            continue
                        sp = _ncbi_taxid_to_species(tax, int(tid), ncbi_cache)
                    if not sp:
                        continue
                    read_mass[sp] = read_mass.get(sp, 0.0) + float(p)
                    prev = read_max.get(sp, 0.0)
                    if p > prev:
                        read_max[sp] = float(p)

                for sp, m in read_mass.items():
                    any_support[sp] += 1
                    prob_mass[sp] = prob_mass.get(sp, 0.0) + float(m)
                    if float(read_max.get(sp, 0.0)) >= float(prob_hi):
                        hi_support[sp] += 1

    return any_support, hi_support, prob_mass


def _genus_local_split_species_counter(
    *,
    input_files: Iterable[str],
    expect_numeric: bool,
    resolved_kind: str,
    tax: Optional[NCBITaxa],
    gtdb_taxonomy: Optional[GtdbTaxonomy],
    topk: int,
    genus_thres: float,
    alpha: float,
    top_k: int,
) -> Counter[str]:
    """Genus-local split using POST_TOPK (two-pass)."""
    ncbi_cache: Dict[int, Optional[str]] = {}
    gtdb_cache: Dict[str, Optional[str]] = {}

    def taxid_to_species(tid: str) -> Optional[str]:
        if resolved_kind == "gtdb":
            if gtdb_taxonomy is None:
                return None
            return _gtdb_species_name_for_node(gtdb_taxonomy, tid, gtdb_cache)
        if tax is None:
            return None
        try:
            itid = int(tid)
        except ValueError:
            return None
        return _ncbi_species_name_for_taxid(tax, itid, ncbi_cache)

    # First pass: accumulate support per-genus per-species.
    support: Dict[str, Dict[str, float]] = {}
    total_reads = 0
    used_reads = 0
    for input_file in input_files:
        with open(input_file, "r", errors="ignore") as infile:
            for raw in infile:
                line = raw.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                tokens = [t.strip() for t in parts[1:] if t.strip()]
                if not tokens:
                    continue
                total_reads += 1
                primary = _parse_primary_taxid_and_count(line, expect_numeric=expect_numeric)
                if primary is None:
                    continue
                _, mass = primary

                cand = _parse_post_topk_from_tokens(tokens, expect_numeric=expect_numeric, topk=topk)
                if not cand:
                    continue

                genus_post: Dict[str, float] = {}
                sp_post: Dict[str, float] = {}
                for tid, prob in cand:
                    sp = taxid_to_species(tid)
                    if not sp:
                        continue
                    gname = _infer_genus(sp)
                    if gname is None:
                        continue
                    genus_post[gname] = genus_post.get(gname, 0.0) + prob
                    prev = sp_post.get(sp)
                    if prev is None or prob > prev:
                        sp_post[sp] = prob
                if not genus_post:
                    continue
                top_genus, top_mass = max(genus_post.items(), key=lambda kv: kv[1])
                if top_mass < float(genus_thres):
                    continue
                used_reads += 1
                for sp, prob in sp_post.items():
                    if _infer_genus(sp) != top_genus:
                        continue
                    support.setdefault(top_genus, {})
                    support[top_genus][sp] = support[top_genus].get(sp, 0.0) + mass * prob

    # Determine keep sets per genus.
    keep: Dict[str, set[str]] = {}
    for gname, mp in support.items():
        items = sorted(mp.items(), key=lambda kv: kv[1], reverse=True)
        total = sum(v for _, v in items)
        kept: set[str] = set()
        cum = 0.0
        for sp, val in items:
            if len(kept) >= int(top_k):
                break
            kept.add(sp)
            cum += val
            if total > 0 and cum >= float(alpha) * total:
                break
        if not kept and items:
            kept.add(items[0][0])
        keep[gname] = kept

    # Second pass: redistribute within genus keep set.
    out: Counter[str] = Counter()
    corrected = 0
    for input_file in input_files:
        with open(input_file, "r", errors="ignore") as infile:
            for raw in infile:
                line = raw.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                tokens = [t.strip() for t in parts[1:] if t.strip()]
                if not tokens:
                    continue
                primary = _parse_primary_taxid_and_count(line, expect_numeric=expect_numeric)
                if primary is None:
                    out[UNCLASSIFIED] += 1
                    continue
                top_tid, mass = primary
                base_sp = taxid_to_species(top_tid) or UNCLASSIFIED

                cand = _parse_post_topk_from_tokens(tokens, expect_numeric=expect_numeric, topk=topk)
                if not cand:
                    out[base_sp] += mass
                    continue

                genus_post: Dict[str, float] = {}
                sp_post: Dict[str, float] = {}
                for tid, prob in cand:
                    sp = taxid_to_species(tid)
                    if not sp:
                        continue
                    gname = _infer_genus(sp)
                    if gname is None:
                        continue
                    genus_post[gname] = genus_post.get(gname, 0.0) + prob
                    prev = sp_post.get(sp)
                    if prev is None or prob > prev:
                        sp_post[sp] = prob
                if not genus_post:
                    out[base_sp] += mass
                    continue

                top_genus, top_mass = max(genus_post.items(), key=lambda kv: kv[1])
                allowed = keep.get(top_genus)
                if top_mass >= float(genus_thres) and allowed:
                    sum_prob = 0.0
                    for sp, prob in sp_post.items():
                        if sp in allowed and _infer_genus(sp) == top_genus:
                            sum_prob += prob
                    if sum_prob > 0.0:
                        for sp, prob in sp_post.items():
                            if sp in allowed and _infer_genus(sp) == top_genus:
                                out[sp] += mass * (prob / sum_prob)
                        if base_sp not in allowed:
                            corrected += 1
                        continue
                out[base_sp] += mass

    print(
        "Genus-local split: total="
        f"{total_reads} used={used_reads} genera={len(keep)} corrected={corrected}"
    )
    return out


def _parse_primary_taxid(line: str, expect_numeric: bool) -> Optional[str]:
    line = line.strip()
    if not line:
        return None
    parts = [token.strip() for token in line.split("\t") if token.strip()]
    if len(parts) < 2:
        return None
    for token in parts[1:]:
        if not token or token.startswith("POST_"):
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


def _parse_primary_taxid_and_count(
    line: str, expect_numeric: bool
) -> Optional[Tuple[str, float]]:
    line = line.strip()
    if not line:
        return None
    parts = [token.strip() for token in line.split("\t") if token.strip()]
    if len(parts) < 2:
        return None
    for token in parts[1:]:
        if not token or token.startswith("POST_"):
            continue
        if token.lower() == UNCLASSIFIED:
            return None
        if ":" not in token:
            continue
        tid_token, count_token = token.split(":", 1)
        tid_token = tid_token.strip()
        if expect_numeric and not tid_token.isdigit():
            continue
        try:
            c = float(count_token)
        except ValueError:
            continue
        if c <= 0:
            continue
        return tid_token, c
    return None


def _collect_taxids(
    input_files: Iterable[str], expect_numeric: bool, abundance_mode: str = "soft_seq"
) -> Tuple[Counter[str], int]:
    mode = (abundance_mode or "soft_seq").strip().lower()
    if mode in {"soft", "soft_seq", "soft_tax"}:
        mode = "soft_seq"
    elif mode in {"hard", "hard_seq"}:
        mode = "hard_seq"
    elif mode in {"top1"}:
        mode = "top1"
    else:
        raise ValueError(f"Unknown abundance_mode: {abundance_mode}")

    taxid_counts: Counter[str] = Counter()
    unclassified_reads = 0
    for input_file in input_files:
        with open(input_file, "r", errors="ignore") as infile:
            for raw in infile:
                line = raw.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue

                tokens = parts[1:]
                if mode == "soft_seq":
                    has_taxid = False
                    for token in tokens:
                        token = token.strip()
                        if not token or token.startswith("POST_"):
                            continue
                        if token.lower() == UNCLASSIFIED:
                            continue
                        if ":" not in token:
                            continue
                        tid_token, count_token = token.split(":", 1)
                        tid_token = tid_token.strip()
                        if expect_numeric and not tid_token.isdigit():
                            continue
                        try:
                            c = float(count_token)
                        except ValueError:
                            continue
                        if c <= 0:
                            continue
                        has_taxid = True
                        taxid_counts[tid_token] += c
                    if not has_taxid:
                        unclassified_reads += 1
                elif mode == "hard_seq":
                    parsed = _parse_primary_taxid_and_count(
                        line, expect_numeric=expect_numeric
                    )
                    if parsed is None:
                        unclassified_reads += 1
                    else:
                        tid, c = parsed
                        taxid_counts[tid] += c
                else:  # top1
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
    abundance_mode: str = "soft_seq",
    hd_species_head_mass: Optional[float] = 97.0,
) -> None:
    input_files = list(input_files)
    meta = read_taxonomy_meta(taxonomy_meta)
    resolved_kind = normalize_kind(taxonomy_kind)
    if resolved_kind == "auto":
        resolved_kind = normalize_kind(meta.kind)
    if resolved_kind == "auto":
        resolved_kind = "ncbi"

    expect_numeric = resolved_kind != "gtdb"
    if hd_species_head_mass is not None and not (0.0 <= float(hd_species_head_mass) <= 100.0):
        raise ValueError("--hd-species-head-mass must be within [0, 100]")

    taxid_counts, base_unclassified = _collect_taxids(
        input_files, expect_numeric=expect_numeric, abundance_mode=abundance_mode
    )

    tax: Optional[NCBITaxa] = None
    gtdb_taxonomy: Optional[GtdbTaxonomy] = None
    if resolved_kind == "gtdb":
        if not taxonomy_info:
            raise ValueError("GTDB 模式需要提供 --taxonomy-info (tax.info 文件路径)")
        gtdb_taxonomy = load_gtdb_taxonomy(taxonomy_info)
    else:
        try:
            tax = NCBITaxa()
        except Exception as exc:  # pragma: no cover - 初始化失败极少发生
            warnings.warn(f"无法初始化 NCBITaxa，所有条目将标记为未分类: {exc}")
            tax = None

    def aggregate_levels_for_counts(
        counts: Counter[str], unclassified_reads: int
    ) -> Dict[str, Counter[str]]:
        if resolved_kind == "gtdb":
            if gtdb_taxonomy is None:  # pragma: no cover - should not happen
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

    total_counts_by_level = {
        level: sum(counter.values()) for level, counter in count_by_level.items()
    }

    top1_taxid_counts: Optional[Counter[str]] = None
    top1_unclassified: Optional[int] = None
    top1_by_level: Optional[Dict[str, Counter[str]]] = None

    def get_top1_by_level() -> Dict[str, Counter[str]]:
        nonlocal top1_taxid_counts, top1_unclassified, top1_by_level
        if top1_by_level is not None:
            return top1_by_level
        top1_taxid_counts, top1_unclassified = _collect_taxids(
            input_files, expect_numeric=expect_numeric, abundance_mode="top1"
        )
        top1_by_level = aggregate_levels_for_counts(
            top1_taxid_counts, int(top1_unclassified)
        )
        return top1_by_level

    # Auto mode override for low-diversity samples:
    # soft_seq tends to spread evidence into a long tail, which hurts presence precision
    # on MOCK-style mixtures. For low diversity, fall back to top1 (per-read) abundance,
    # then apply the species-level low-diversity cutoff.
    normalized_mode = (abundance_mode or "soft_seq").strip().lower()
    if normalized_mode in {"soft", "soft_seq", "soft_tax"}:
        normalized_mode = "soft_seq"
    elif normalized_mode in {"hard", "hard_seq"}:
        normalized_mode = "hard_seq"
    elif normalized_mode in {"top1"}:
        normalized_mode = "top1"
    else:
        normalized_mode = "soft_seq"

    species_counter = count_by_level.get("species")
    species_total = total_counts_by_level.get("species", 0)
    is_low_diversity = False
    force_low_diversity = False
    if species_counter and species_total > 0 and sum(species_counter.values()) > 0:
        # Use top1 distribution for low-diversity detection to avoid soft_seq tail inflation.
        detect_counter = species_counter
        detect_total = species_total
        try:
            detect_by_level = get_top1_by_level()
            detect_counter = detect_by_level.get("species") or species_counter
            detect_total = sum(detect_counter.values()) if detect_counter else species_total
        except Exception:
            detect_counter = species_counter
            detect_total = species_total

        uncls = detect_counter.get(UNCLASSIFIED, 0)
        eff_total = max(0, detect_total - uncls)
        rels = []
        shannon_index = 0.0
        if eff_total > 0:
            for taxon, count in detect_counter.items():
                if taxon == UNCLASSIFIED:
                    continue
                p = count / eff_total
                if p > 0:
                    shannon_index -= p * math.log(p)
                rels.append(p)
        rels.sort(reverse=True)
        top_mass = sum(rels[:LOW_DIVERSITY_TOP_K])
        eff_species = math.exp(shannon_index) if shannon_index > 0 else 0.0
        is_low_diversity = (top_mass >= LOW_DIVERSITY_TOP_MASS) and (
            eff_species <= LOW_DIVERSITY_EFF_SPECIES_MAX
        )
        print(
            f"[profile] low-div detect: shannon={shannon_index:.3f} "
            f"eff_species={eff_species:.2f} top_mass={top_mass * 100.0:.2f} "
            f"species_total={detect_total} is_low_diversity={is_low_diversity}"
        )

        if is_low_diversity:
            force_low_diversity = True
            if normalized_mode != "top1":
                print("[profile] detected low-diversity: switching abundance_mode to top1")
            abundance_mode = "top1"
            count_by_level = get_top1_by_level()
            total_counts_by_level = {
                level: sum(counter.values()) for level, counter in count_by_level.items()
            }

    # Optional: compute top1 support counts (sequence-level) for high-diversity species trimming.
    support_by_level: Optional[Dict[str, Counter[str]]] = None
    if HIGH_DIVERSITY_SPECIES_SUPPORT_MIN and HIGH_DIVERSITY_SPECIES_SUPPORT_MIN > 1:
        # If we already run top1 mode (e.g., low-diversity auto override), support is redundant.
        mode_for_support = (abundance_mode or "soft_seq").strip().lower()
        if mode_for_support not in {"top1"}:
            support_by_level = get_top1_by_level()

    # Low-diversity status used for species-level post-processing.
    # NOTE: This uses the *current* species distribution (soft_seq or top1 after auto-switch),
    # matching the original logic inside the output loop.
    species_is_low_diversity = force_low_diversity
    if not species_is_low_diversity:
        final_species_counter = count_by_level.get("species")
        final_species_total = total_counts_by_level.get("species", 0)
        if (
            final_species_counter
            and final_species_total > 0
            and sum(final_species_counter.values()) > 0
        ):
            uncls_count = final_species_counter.get(UNCLASSIFIED, 0)
            eff_total = max(0, final_species_total - uncls_count)
            rels = []
            shannon_eff = 0.0
            if eff_total > 0:
                for taxon, count in final_species_counter.items():
                    if taxon == UNCLASSIFIED:
                        continue
                    p = count / eff_total
                    if p > 0:
                        shannon_eff -= p * math.log(p)
                    rels.append(p)
            rels.sort(reverse=True)
            top_mass = sum(rels[:LOW_DIVERSITY_TOP_K])
            eff_species = math.exp(shannon_eff) if shannon_eff > 0 else 0.0
            species_is_low_diversity = (top_mass >= LOW_DIVERSITY_TOP_MASS) and (
                eff_species <= LOW_DIVERSITY_EFF_SPECIES_MAX
            )

    with open(output_file + ".tsv", "w") as outfile:
        mode = (abundance_mode or "soft_seq").strip().lower()
        if mode in {"soft", "soft_seq", "soft_tax"}:
            mode = "soft_seq"
        elif mode in {"hard", "hard_seq"}:
            mode = "hard_seq"
        elif mode in {"top1"}:
            mode = "top1"
        else:
            mode = "soft_seq"
        outfile.write(f"# abundance_mode={mode}\n")
        count_unit = "sequence" if mode == "top1" else "posterior_evidence"
        outfile.write(f"# count_unit={count_unit}\n")
        outfile.write("Taxon\tCount\tRelative Abundance (%)\n")

        for level, taxon_counter in count_by_level.items():
            total_count = total_counts_by_level[level]
            if total_count == 0:
                continue

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

            filtered_items = list(items)

            if level == "species" and filtered_items:
                uncls_count = taxon_counter.get(UNCLASSIFIED, 0)
                eff_total = max(0, total_count - uncls_count)

                if species_is_low_diversity:
                    # Low-diversity: head-mass prune on species to drop tiny FP tails.
                    # (Exclude UNCLASSIFIED from the denominator and keep at least K taxa.)
                    if (
                        LOW_DIVERSITY_SPECIES_HEAD_MASS
                        and float(LOW_DIVERSITY_SPECIES_HEAD_MASS) > 0
                        and eff_total > 0
                        and filtered_items
                    ):
                        non_uncls: List[Tuple[str, float]] = []
                        uncls_item = None
                        for taxon, count in filtered_items:
                            if taxon == UNCLASSIFIED:
                                uncls_item = (taxon, float(count))
                            else:
                                non_uncls.append((taxon, float(count)))
                        non_uncls.sort(key=lambda it: it[1], reverse=True)

                        # Use the mass after species-level filtering as denominator.
                        # `eff_total` may include tiny tails already removed by prefilters,
                        # which can make head-mass pruning keep one extra FP taxon.
                        denom_total = sum(c for _, c in non_uncls)

                        before_n = len(non_uncls)
                        kept: List[Tuple[str, float]] = []
                        cum = 0.0
                        head_mass = float(LOW_DIVERSITY_SPECIES_HEAD_MASS)
                        min_keep = int(LOW_DIVERSITY_SPECIES_HEAD_MIN_KEEP)
                        for i, (taxon, count) in enumerate(non_uncls):
                            kept.append((taxon, count))
                            cum += float(count)
                            if (
                                (i + 1) >= min_keep
                                and denom_total > 0
                                and ((cum / float(denom_total)) * 100.0) >= head_mass
                            ):
                                break
                        filtered_items = kept + ([uncls_item] if uncls_item else [])
                        after_n = len(kept)
                        kept_pct = (
                            ((cum / float(denom_total)) * 100.0) if denom_total > 0 else 0.0
                        )
                        print(
                            "[profile][auto] lowdiv_species_headmass="
                            f"enabled=1 head_mass={head_mass:.2f} min_keep={min_keep} "
                            f"before={before_n} after={after_n} kept_pct={kept_pct:.2f}"
                        )

                    # Low-diversity: ultra-conservative presence-only rescue for NEW genera
                    # using high-probability POST_TOPK evidence.
                    #
                    # Rationale: Some true low-abundance species may never win top1 but can
                    # appear as a strong alternative (p>=prob_min) in a limited number of reads.
                    # We only allow adding NEW genera (to avoid intragenus long-tail explosion),
                    # and we add only a tiny pseudo-count (>=MIN_EVIDENCE) so abundance stays stable.
                    has_post_topk = _has_post_topk_tokens(input_files)
                    prob_hi = 0.20
                    support_min = max(50, int(eff_total * 0.00002))
                    # Guard against "hitchhiker" taxa that show up in many POST_TOPK lists
                    # with tiny probabilities (common in low-div samples, especially within
                    # Enterobacteriaceae-like families). This prevents spurious taxa from
                    # passing support_min and inflating presence FPs.
                    support_hi_frac_min = 0.01

                    max_added = 3
                    added: List[Tuple[str, int, int, float]] = []
                    if has_post_topk and prob_hi > 0 and support_min > 0:
                        existing_species = {
                            taxon
                            for taxon, _ in filtered_items
                            if taxon != UNCLASSIFIED
                        }
                        existing_genera: set[str] = set()
                        for sp in existing_species:
                            g = _infer_genus(sp)
                            if g:
                                existing_genera.add(g)

                        any_support, hi_support, prob_mass = _collect_post_topk_species_stats(
                            input_files=input_files,
                            expect_numeric=expect_numeric,
                            resolved_kind=resolved_kind,
                            tax=tax,
                            gtdb_taxonomy=(
                                gtdb_taxonomy if resolved_kind == "gtdb" else None
                            ),
                            prob_hi=float(prob_hi),
                            topk=16,
                        )
                        candidates: List[Tuple[str, int, int, float, float]] = []
                        for sp, hi in hi_support.items():
                            if not sp or sp == UNCLASSIFIED:
                                continue
                            if not _is_lowdiv_rescue_species_name_ok(sp):
                                continue
                            if sp in existing_species:
                                continue
                            any_c = int(any_support.get(sp, 0))
                            if any_c <= 0:
                                continue
                            mass = float(prob_mass.get(sp, 0.0))
                            frac = (float(hi) / float(any_c)) if any_c > 0 else 0.0
                            if int(hi) < int(support_min):
                                continue
                            if frac < float(support_hi_frac_min):
                                continue
                            g = _infer_genus(sp)
                            if g and g in existing_genera:
                                continue
                            candidates.append((sp, int(hi), any_c, mass, frac))
                        # Prefer stronger evidence: higher posterior mass first.
                        candidates.sort(
                            key=lambda it: (it[3], it[1], it[4]),
                            reverse=True,
                        )

                        if candidates:
                            pseudo = float(MIN_EVIDENCE)
                            keep = candidates[: int(max_added)]
                            uncls_item = None
                            kept_items: List[Tuple[str, float]] = []
                            for taxon, count in filtered_items:
                                if taxon == UNCLASSIFIED:
                                    uncls_item = (taxon, count)
                                else:
                                    kept_items.append((taxon, float(count)))
                            # append rescued taxa (tiny counts; keep them at tail)
                            for sp, hi, any_c, mass, _ in keep:
                                kept_items.append((sp, pseudo))
                                added.append((sp, int(hi), int(any_c), float(mass)))
                            # re-sort non-unclassified part by count desc to keep output stable
                            kept_items.sort(key=lambda it: it[1], reverse=True)
                            filtered_items = kept_items + ([uncls_item] if uncls_item else [])

                    if has_post_topk:
                        added_text = (
                            ",".join(
                                [
                                    f"{sp}:hi={hi}/any={any_c}:mass={mass:.2f}"
                                    for sp, hi, any_c, mass in added
                                ]
                            )
                            if added
                            else ""
                        )
                        print(
                            "[profile][auto] lowdiv_posttopk_newgenus_rescue="
                            f"enabled=1 prob_hi={prob_hi:.2f} "
                            f"support_min={support_min} "
                            f"support_hi_frac_min={float(support_hi_frac_min):.4f} "
                            f"max_added={max_added} "
                            f"added={len(added)} [{added_text}]"
                        )
                    else:
                        print(
                            "[profile][auto] lowdiv_posttopk_newgenus_rescue="
                            "enabled=0 reason=no_post_topk"
                        )
                else:
                    if HIGH_DIVERSITY_GENUS_LOCAL_SPLIT:
                        split = _genus_local_split_species_counter(
                            input_files=input_files,
                            expect_numeric=expect_numeric,
                            resolved_kind=resolved_kind,
                            tax=tax,
                            gtdb_taxonomy=(gtdb_taxonomy if resolved_kind == "gtdb" else None),
                            topk=int(HIGH_DIVERSITY_GENUS_LOCAL_POST_TOPK),
                            genus_thres=float(HIGH_DIVERSITY_GENUS_LOCAL_THRES),
                            alpha=float(HIGH_DIVERSITY_GENUS_LOCAL_ALPHA),
                            top_k=int(HIGH_DIVERSITY_GENUS_LOCAL_TOPK),
                        )
                        if split:
                            taxon_counter = split
                            total_count = sum(taxon_counter.values())
                            uncls_count = taxon_counter.get(UNCLASSIFIED, 0)
                            eff_total = max(0, total_count - uncls_count)
                            items = taxon_counter.most_common()
                            filtered_items = list(items)

                    # High-diversity: trim long tail by head mass (exclude UNCLASSIFIED).
                    # This suppresses hundreds of tiny non-zero species that destroy
                    # presence precision, without affecting per-read assignments.
                    if hd_species_head_mass is not None:
                        denom = eff_total if eff_total > 0 else total_count
                        kept = set()
                        cum = 0.0
                        uncls = None
                        for taxon, count in filtered_items:
                            if taxon == UNCLASSIFIED:
                                uncls = (taxon, count)
                                continue
                            rel = (count / denom) * 100 if denom > 0 else 0.0
                            kept.add(taxon)
                            cum += rel
                            if cum >= float(hd_species_head_mass):
                                break

                        # Union with top1 support>=K to keep repeated low-abundance taxa
                        # while filtering scatter-type singletons.
                        if (
                            HIGH_DIVERSITY_SPECIES_SUPPORT_MIN
                            and HIGH_DIVERSITY_SPECIES_SUPPORT_MIN > 1
                            and support_by_level
                        ):
                            support_species = support_by_level.get("species")
                            if support_species:
                                for taxon, sup in support_species.items():
                                    if taxon != UNCLASSIFIED and sup >= HIGH_DIVERSITY_SPECIES_SUPPORT_MIN:
                                        kept.add(taxon)

                        if kept:
                            rebuilt = []
                            for taxon, count in items:
                                if taxon == UNCLASSIFIED:
                                    if uncls is not None:
                                        rebuilt.append(uncls)
                                    continue
                                if taxon in kept:
                                    rebuilt.append((taxon, count))
                            if uncls is not None and not rebuilt:
                                rebuilt.append(uncls)
                            if rebuilt:
                                filtered_items = rebuilt

                    # Final hard cutoff on species relative abundance (high-diversity only).
                    if (
                        level == "species"
                        and not species_is_low_diversity
                        and HIGH_DIVERSITY_FINAL_REL_CUTOFF
                        and HIGH_DIVERSITY_FINAL_REL_CUTOFF > 0
                        and filtered_items
                    ):
                        denom = eff_total if eff_total > 0 else total_count
                        final_items: List[Tuple[str, float]] = []
                        for taxon, count in filtered_items:
                            if taxon == UNCLASSIFIED:
                                final_items.append((taxon, count))
                                continue
                            rel = (count / denom) * 100 if denom > 0 else 0.0
                            if rel >= float(HIGH_DIVERSITY_FINAL_REL_CUTOFF):
                                final_items.append((taxon, count))
                        if final_items:
                            filtered_items = final_items

            if not filtered_items:
                # 阈值过高导致清空时，保留原始未截断结果
                filtered_items = items

            total_filtered = sum(c for _, c in filtered_items)

            for taxon, count in filtered_items:
                relative_abundance = (count / total_filtered) * 100 if total_filtered > 0 else 0.0
                outfile.write(
                    f"{taxon}\t{count}\t{relative_abundance:.2f}\n"
                )
