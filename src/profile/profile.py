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

# Profile 的核心输出是 species 级别（可聚合得到丰度 profile）。
# 为减少无用复杂度与输出噪音，这里只输出 species 一个层级。

UNCLASSIFIED = "unclassified"
# 统一基底量：K（头部尺度）与 p（POST_TOPK 高置信概率门）
TAIL_RISK_TOP_K = 10
PROB_HI_BASE = 0.20

# 最小 evidence：由 p 派生，减少 magic number
MIN_EVIDENCE = int(math.ceil(2.0 / PROB_HI_BASE))
# Tail-risk 锚点：集中生态（head-heavy）下，top1 物种分布通常 head mass 高且有效物种数小。
TAIL_RISK_TOP_MASS_ANCHOR = 0.75
TAIL_RISK_EFF_SPECIES_ANCHOR = float((3 * TAIL_RISK_TOP_K) + 2)
TAIL_RISK_R_ANCHOR = max(0.0, 1.0 - TAIL_RISK_TOP_MASS_ANCHOR)
# 集中生态时更激进的 species 头部裁剪：减少长尾假阳性。
HEAD_HEAVY_SPECIES_HEAD_MASS_BASE = 100.0 - (
    TAIL_RISK_TOP_MASS_ANCHOR * PROB_HI_BASE
)
SPECIES_HEAD_MIN_KEEP = int(
    math.ceil(0.8 * float(TAIL_RISK_TOP_K))
)
# 尾部更丰富时：最终写出时的 species 相对丰度硬阈值上限（%）。
TAIL_RICH_FINAL_REL_CUTOFF_MAX = (
    PROB_HI_BASE * PROB_HI_BASE
) / 2.0
# 尾部更丰富时：散射型 top1 FP 往往是 singleton/doubleton。
# 用 top1 的“出现次数”做第二信号：只要 support 足够（>=K）就保留，
# 以便在提高 head mass（救 recall）时仍能压住 FP 物种种类数。
TAIL_RICH_SPECIES_SUPPORT_MIN = 3
# head-heavy rescue support gate constants derived from {K, p}
HEAD_HEAVY_RESCUE_SUPPORT_MIN_FRAC = (PROB_HI_BASE**4) / (8.0 * float(TAIL_RISK_TOP_K))
HEAD_HEAVY_RESCUE_SUPPORT_MIN_FLOOR = int(5 * MIN_EVIDENCE)


def _clip01(value: float) -> float:
    if value < 0.0:
        return 0.0
    if value > 1.0:
        return 1.0
    return value


def _smoothstep(value: float, edge0: float, edge1: float) -> float:
    if edge1 <= edge0:
        return 1.0 if value >= edge1 else 0.0
    t = _clip01((float(value) - float(edge0)) / (float(edge1) - float(edge0)))
    return t * t * (3.0 - (2.0 * t))


def _compute_top1_tail_metrics(
    detect_counter: Counter[str], detect_total: float
) -> Tuple[float, float, float, float]:
    """Compute tail-risk probe metrics from top1 species distribution.

    Returns:
    - shannon_index
    - eff_species
    - top_mass (sum of top-K relative masses)
    - tail_risk_r = clip(1 - top_mass, 0, 1)
    """
    uncls = detect_counter.get(UNCLASSIFIED, 0)
    eff_total = max(0.0, float(detect_total) - float(uncls))
    rels: List[float] = []
    shannon_index = 0.0
    if eff_total > 0:
        for taxon, count in detect_counter.items():
            if taxon == UNCLASSIFIED:
                continue
            p = float(count) / eff_total
            if p > 0:
                shannon_index -= p * math.log(p)
            rels.append(p)
    rels.sort(reverse=True)
    top_mass = sum(rels[:TAIL_RISK_TOP_K])
    eff_species = math.exp(shannon_index) if shannon_index > 0 else 0.0
    tail_risk_r = _clip01(1.0 - top_mass)
    return shannon_index, eff_species, top_mass, tail_risk_r


def _compute_tail_risk_u(tail_risk_r: float, eff_species: float) -> float:
    """Unified risk intensity in [0,1] from tail mass and effective species."""
    r = _clip01(tail_risk_r)
    anchor = max(float(TAIL_RISK_R_ANCHOR), 1e-9)
    u_r = _clip01(r / (r + anchor))

    e0 = max(float(TAIL_RISK_EFF_SPECIES_ANCHOR), 1e-9)
    e = max(0.0, float(eff_species))
    u_e = _clip01((e - e0) / (e + e0))
    # Soft-OR merge: high in either channel pushes risk up.
    return _clip01(1.0 - ((1.0 - u_r) * (1.0 - u_e)))


def _compute_tail_risk_s(tail_risk_u: float) -> float:
    """Endpoint-hard continuous saturation of risk u in [0,1]."""
    u = _clip01(float(tail_risk_u))
    beta = max(1, int(math.ceil(1.0 / max(float(PROB_HI_BASE), 1e-9))))
    up = u**beta
    down = (1.0 - u) ** beta
    return _clip01(up / (up + down + 1e-12))


def _normalize_abundance_mode(abundance_mode: Optional[str]) -> str:
    mode = (abundance_mode or "soft_seq").strip().lower()
    if mode in {"soft", "soft_seq", "soft_tax"}:
        return "soft_seq"
    if mode in {"hard", "hard_seq"}:
        return "hard_seq"
    if mode in {"top1"}:
        return "top1"
    raise ValueError(f"Unknown abundance_mode: {abundance_mode}")


def _infer_genus(species_name: str) -> Optional[str]:
    """Infer genus name from a species display string.

    Intended for output-only post-processing in tail-rich species profiles.
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


def _is_rescue_species_name_ok(species_name: str) -> bool:
    """Heuristic filter for head-heavy POST_TOPK rescue candidates.

    Goal: avoid adding generic/low-quality nodes like "X bacterium" or "Y sp."
    that often act as FP sinks in head-heavy mixtures.
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
        if rank != "species":
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
        if node_rank != "species":
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
) -> Tuple[
    Counter[str],
    Counter[str],
    Dict[str, float],
    Dict[str, float],
]:
    """Collect POST_TOPK per-species stats in a single pass."""
    any_support: Counter[str] = Counter()
    hi_support: Counter[str] = Counter()
    prob_mass: Dict[str, float] = {}
    hi_mass: Dict[str, float] = {}

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
                token_list = [t.strip() for t in tokens if t.strip()]
                cand = _parse_post_topk_from_tokens(
                    token_list,
                    expect_numeric=expect_numeric,
                    topk=int(topk),
                )
                if not cand:
                    continue

                # per-read aggregation at species level
                read_mass: Dict[str, float] = {}
                read_max: Dict[str, float] = {}
                read_hi_mass: Dict[str, float] = {}
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
                    if p >= float(prob_hi):
                        read_hi_mass[sp] = read_hi_mass.get(sp, 0.0) + float(p)

                for sp, m in read_mass.items():
                    any_support[sp] += 1
                    prob_mass[sp] = prob_mass.get(sp, 0.0) + float(m)
                    if float(read_max.get(sp, 0.0)) >= float(prob_hi):
                        hi_support[sp] += 1
                    if sp in read_hi_mass:
                        hi_mass[sp] = hi_mass.get(sp, 0.0) + float(read_hi_mass[sp])

    return (
        any_support,
        hi_support,
        prob_mass,
        hi_mass,
    )


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
    mode = _normalize_abundance_mode(abundance_mode)

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
    species: Counter[str] = Counter()
    if base_unclassified:
        species[UNCLASSIFIED] += base_unclassified
    if tax is None:
        return {"species": species}

    cache: Dict[int, Optional[str]] = {}
    for taxid, count in taxid_counts.items():
        sp = _ncbi_taxid_to_species(tax, int(taxid), cache)
        if sp:
            species[sp] += count
        else:
            species[UNCLASSIFIED] += count
    return {"species": species}


def _aggregate_gtdb_levels(
    taxonomy: GtdbTaxonomy,
    taxid_counts: Counter[str],
    base_unclassified: int,
) -> Dict[str, Counter[str]]:
    species: Counter[str] = Counter()
    if base_unclassified:
        species[UNCLASSIFIED] += base_unclassified

    cache: Dict[str, Optional[str]] = {}
    for node, count in taxid_counts.items():
        if node not in taxonomy.rank:
            species[UNCLASSIFIED] += count
            continue
        sp = _gtdb_node_to_species(taxonomy, node, cache)
        if sp:
            species[sp] += count
        else:
            species[UNCLASSIFIED] += count
    return {"species": species}


def process_file(
    input_files: Iterable[str],
    output_file: str,
    taxonomy_kind: str = "auto",
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
    normalized_mode = _normalize_abundance_mode(abundance_mode)

    taxid_counts, base_unclassified = _collect_taxids(
        input_files, expect_numeric=expect_numeric, abundance_mode=normalized_mode
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

    # Unified mode handling:
    # keep abundance_mode stable, and use risk-driven blending in soft_seq.

    species_counter = count_by_level.get("species")
    species_total = total_counts_by_level.get("species", 0)
    detected_tail_risk_r = 1.0
    detected_eff_species = float("inf")
    detected_top_mass = 0.0
    if species_counter and species_total > 0:
        # Use top1 distribution for tail-risk detection to avoid soft_seq tail inflation.
        detect_counter = species_counter
        detect_total = species_total
        try:
            detect_by_level = get_top1_by_level()
            detect_counter = detect_by_level.get("species") or species_counter
            detect_total = sum(detect_counter.values()) if detect_counter else species_total
        except Exception:
            detect_counter = species_counter
            detect_total = species_total

        shannon_index, eff_species, top_mass, tail_risk_r = _compute_top1_tail_metrics(
            detect_counter, float(detect_total)
        )
        detected_eff_species = eff_species
        detected_top_mass = top_mass
        detected_tail_risk_r = tail_risk_r

        if normalized_mode == "soft_seq":
            # Continuous blend in relative-abundance space:
            # u≈0 => near top1, u≈1 => near soft_seq.
            try:
                top1_species = get_top1_by_level().get("species")
            except Exception:
                top1_species = None
            if top1_species:
                soft_species = species_counter or Counter()
                soft_uncls = float(soft_species.get(UNCLASSIFIED, 0.0))
                top1_uncls = float(top1_species.get(UNCLASSIFIED, 0.0))
                soft_eff = sum(float(v) for k, v in soft_species.items() if k != UNCLASSIFIED)
                top1_eff = sum(float(v) for k, v in top1_species.items() if k != UNCLASSIFIED)
                if soft_eff > 0 and top1_eff > 0:
                    u_mix = _compute_tail_risk_u(
                        detected_tail_risk_r, detected_eff_species
                    )
                    s_mix = _compute_tail_risk_s(u_mix)
                    w_soft = float(s_mix)
                    w_top1 = 1.0 - float(w_soft)
                    species_union = {
                        k for k in soft_species.keys() if k != UNCLASSIFIED
                    } | {
                        k for k in top1_species.keys() if k != UNCLASSIFIED
                    }
                    mixed: Counter[str] = Counter()
                    base_eff = (
                        (float(w_top1) * float(top1_eff))
                        + (float(w_soft) * float(soft_eff))
                    )
                    for sp in species_union:
                        p_soft = float(soft_species.get(sp, 0.0)) / float(soft_eff)
                        p_top1 = float(top1_species.get(sp, 0.0)) / float(top1_eff)
                        p_mix = (float(w_top1) * p_top1) + (float(w_soft) * p_soft)
                        if p_mix > 0:
                            mixed[sp] = float(p_mix) * float(base_eff)
                    soft_total = float(soft_eff + soft_uncls)
                    top1_total = float(top1_eff + top1_uncls)
                    uncls_rel_soft = (soft_uncls / soft_total) if soft_total > 0 else 0.0
                    uncls_rel_top1 = (top1_uncls / top1_total) if top1_total > 0 else 0.0
                    uncls_rel_mix = (float(w_top1) * uncls_rel_top1) + (
                        float(w_soft) * uncls_rel_soft
                    )
                    uncls_rel_mix = max(0.0, min(0.999999, float(uncls_rel_mix)))
                    total_mix = float(base_eff) / max(1e-9, (1.0 - float(uncls_rel_mix)))
                    uncls_mix = max(0.0, float(total_mix) - float(base_eff))
                    if uncls_mix > 0:
                        mixed[UNCLASSIFIED] = float(uncls_mix)
                    count_by_level["species"] = mixed
                    total_counts_by_level["species"] = sum(mixed.values())
                    species_counter = count_by_level.get("species")
                    species_total = total_counts_by_level.get("species", 0)
                    print(
                        "[profile][auto] mode_blend="
                        f"enabled=1 u={u_mix:.4f} s={s_mix:.4f} w_top1={w_top1:.4f} "
                        f"soft_eff={soft_eff:.1f} top1_eff={top1_eff:.1f} "
                        f"species_union={len(species_union)}"
                    )

    # Optional: compute top1 support counts (sequence-level) for tail-rich species trimming.
    support_by_level: Optional[Dict[str, Counter[str]]] = None
    if TAIL_RICH_SPECIES_SUPPORT_MIN and TAIL_RICH_SPECIES_SUPPORT_MIN > 1:
        # If we already run top1 mode, support is redundant.
        if normalized_mode != "top1":
            support_by_level = get_top1_by_level()

    # NOTE: detection is performed above using the top1 distribution to avoid soft_seq tail inflation.
    tail_risk_u = _compute_tail_risk_u(detected_tail_risk_r, detected_eff_species)
    tail_risk_s = _compute_tail_risk_s(tail_risk_u)
    tail_rich_strength = _smoothstep(float(tail_risk_s), 0.60, 0.95)
    head_heavy_strength = _clip01(1.0 - float(tail_rich_strength))
    policy_species_head_mass = None
    if hd_species_head_mass is not None:
        low_head = float(HEAD_HEAVY_SPECIES_HEAD_MASS_BASE)
        high_head = float(hd_species_head_mass)
        policy_species_head_mass = (low_head * (1.0 - tail_risk_s)) + (
            high_head * tail_risk_s
        )
    species_rel_cutoff = (
        float(TAIL_RICH_FINAL_REL_CUTOFF_MAX) * float(tail_rich_strength)
        if TAIL_RICH_FINAL_REL_CUTOFF_MAX and TAIL_RICH_FINAL_REL_CUTOFF_MAX > 0
        else 0.0
    )
    b_max = int(max(SPECIES_HEAD_MIN_KEEP, TAIL_RISK_TOP_K - 1))
    unified_rescue_budget = max(
        0,
        min(
            b_max,
            int(round(float(b_max) * max(0.0, 1.0 - tail_risk_u))),
        ),
    )
    head_heavy_rescue_budget = max(
        0,
        min(
            int(unified_rescue_budget),
            int(round(float(unified_rescue_budget) * head_heavy_strength)),
        ),
    )
    # Cutoff-boundary mass swap budget (1-in-1-out), activated continuously by tail-rich strength.
    tail_rich_cutoff_swap_budget_max = min(4, int(unified_rescue_budget) * 2)
    tail_rich_cutoff_swap_budget = max(
        0,
        int(round(float(tail_rich_cutoff_swap_budget_max) * tail_rich_strength)),
    )
    if species_counter and species_total > 0:
        print(
            "[profile][auto] tail-risk-policy="
            f"r={detected_tail_risk_r:.4f} top_mass={detected_top_mass * 100.0:.2f} "
            f"u={tail_risk_u:.4f} s={tail_risk_s:.4f} "
            f"head_heavy_strength={head_heavy_strength:.4f} "
            f"tail_rich_strength={tail_rich_strength:.4f} "
            f"head_mass={float(policy_species_head_mass) if policy_species_head_mass is not None else -1.0:.2f} "
            f"rel_cutoff={species_rel_cutoff:.4f} "
            f"head_heavy_rescue_budget={head_heavy_rescue_budget} "
            f"tail_rich_cutoff_swap_budget={tail_rich_cutoff_swap_budget}"
        )

    output_base = output_file
    if output_base.lower().endswith(".tsv"):
        output_base = output_base[: -len(".tsv")]

    with open(output_base + ".tsv", "w") as outfile:
        mode = _normalize_abundance_mode(abundance_mode)
        outfile.write(f"# abundance_mode={mode}\n")
        count_unit = "sequence" if mode == "top1" else "posterior_evidence"
        outfile.write(f"# count_unit={count_unit}\n")
        outfile.write("Taxon\tCount\tRelative Abundance (%)\n")

        taxon_counter = count_by_level.get("species")
        total_count = total_counts_by_level.get("species", 0)
        if not taxon_counter or total_count == 0:
            return

        outfile.write("\n## Species Level ##\n")

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

        filtered_items = list(items)

        if filtered_items:
            uncls_count = taxon_counter.get(UNCLASSIFIED, 0)
            eff_total = max(0, total_count - uncls_count)

            if head_heavy_rescue_budget > 0:
                # Head-heavy rescue stage: prune tiny FP tails first.
                # (Exclude UNCLASSIFIED from the denominator and keep at least K taxa.)
                if (
                    policy_species_head_mass is not None
                    and float(policy_species_head_mass) > 0
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

                    kept: List[Tuple[str, float]] = []
                    cum = 0.0
                    head_mass = float(policy_species_head_mass)
                    min_keep = int(SPECIES_HEAD_MIN_KEEP)
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

                # Head-heavy: ultra-conservative presence-only rescue for NEW genera
                # using high-probability POST_TOPK evidence.
                #
                # Rationale: Some true low-abundance species may never win top1 but can
                # appear as a strong alternative (p>=prob_min) in a limited number of reads.
                # We only allow adding NEW genera (to avoid intragenus long-tail explosion),
                # and we add only a tiny pseudo-count (>=MIN_EVIDENCE) so abundance stays stable.
                prob_hi = float(PROB_HI_BASE)
                support_min = max(
                    int(HEAD_HEAVY_RESCUE_SUPPORT_MIN_FLOOR),
                    int(float(eff_total) * float(HEAD_HEAVY_RESCUE_SUPPORT_MIN_FRAC)),
                )
                # Guard against "hitchhiker" taxa that show up in many POST_TOPK lists
                # with tiny probabilities (common in head-heavy samples, especially within
                # Enterobacteriaceae-like families). This prevents spurious taxa from
                # passing support_min and inflating presence FPs.
                support_hi_frac_min = float(PROB_HI_BASE * PROB_HI_BASE) / 4.0

                max_added = int(head_heavy_rescue_budget)
                added: List[Tuple[str, int, int, float, float, float]] = []
                if prob_hi > 0 and support_min > 0:
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

                    (
                        any_support,
                        hi_support,
                        prob_mass,
                        hi_mass,
                    ) = _collect_post_topk_species_stats(
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
                    candidates: List[Tuple[str, int, int, float, float, float, str, float]] = []
                    concentrated_frac_max = 1.0 - (float(PROB_HI_BASE) / 4.0)
                    narrow_any_max = int(float(support_min) / max(float(prob_hi), 1e-9))
                    for sp, hi in hi_support.items():
                        if not sp or sp == UNCLASSIFIED:
                            continue
                        if not _is_rescue_species_name_ok(sp):
                            continue
                        if sp in existing_species:
                            continue
                        any_c = int(any_support.get(sp, 0))
                        if any_c <= 0:
                            continue
                        mass = float(prob_mass.get(sp, 0.0))
                        hi_mass_value = float(hi_mass.get(sp, 0.0))
                        frac = (float(hi) / float(any_c)) if any_c > 0 else 0.0
                        if int(hi) < int(support_min):
                            continue
                        if frac < float(support_hi_frac_min):
                            continue
                        # Guard: very concentrated (hi≈any) but narrow support often acts as
                        # hitchhiker FP in head-heavy tails; prefer broader repeat evidence.
                        if frac >= float(concentrated_frac_max) and any_c <= int(narrow_any_max):
                            continue
                        g = _infer_genus(sp)
                        if g and g in existing_genera:
                            continue
                        tail_span = max(1, int(any_c - int(hi) + 1))
                        score = (
                            float(hi_mass_value)
                            * math.sqrt(max(float(frac), 1e-9))
                            * math.log1p(float(tail_span))
                        )
                        candidates.append(
                            (
                                sp,
                                int(hi),
                                any_c,
                                mass,
                                hi_mass_value,
                                frac,
                                g or "",
                                float(score),
                            )
                        )
                    # Prefer stronger evidence: high-confidence posterior mass with
                    # broad repeat support (anti-hitchhiker).
                    candidates.sort(
                        key=lambda it: (it[7], it[4], it[3], it[1]),
                        reverse=True,
                    )

                    if candidates:
                        pseudo = float(MIN_EVIDENCE)
                        keep: List[
                            Tuple[str, int, int, float, float, float, str, float]
                        ] = []
                        selected_new_genera: set[str] = set()
                        # Keep the first 3 additions behavior-compatible with the
                        # proven conservative baseline; use stricter gates only for
                        # expansion slots (>3) to avoid head-heavy tail FP spread.
                        base_keep = min(int(max_added), 3)
                        first_score = max(float(candidates[0][7]), 1e-12)
                        knee_score_ratio = max(0.01, float(support_hi_frac_min))
                        # Expansion candidates should look like true rare tails:
                        # enough high-prob support, but not too many broad ambiguous hits.
                        expansion_any_max = int(
                            max(
                                float(support_min),
                                float(support_min) / max(float(prob_hi), 1e-9),
                            )
                        )
                        expansion_first_score: Optional[float] = None
                        expansion_knee_ratio = float(prob_hi)
                        for cand in candidates:
                            if len(keep) >= int(max_added):
                                break
                            sp, hi, any_c, mass, hi_mass_value, frac, g, score = cand
                            if g and g in selected_new_genera:
                                continue
                            # Budget is an upper bound; stop if tail evidence is too weak.
                            if len(keep) >= 3 and float(score) < (first_score * knee_score_ratio):
                                break
                            if len(keep) >= int(base_keep):
                                # Expansion stage: suppress medium-confidence broad tails.
                                if int(any_c) > int(expansion_any_max):
                                    continue
                                if expansion_first_score is None:
                                    expansion_first_score = max(float(score), 1e-12)
                                elif (
                                    len(keep) >= int(base_keep + 1)
                                    and float(score)
                                    < float(expansion_first_score * expansion_knee_ratio)
                                ):
                                    break
                            keep.append(cand)
                            if g:
                                selected_new_genera.add(g)

                        uncls_item = None
                        kept_items: List[Tuple[str, float]] = []
                        for taxon, count in filtered_items:
                            if taxon == UNCLASSIFIED:
                                uncls_item = (taxon, count)
                            else:
                                kept_items.append((taxon, float(count)))
                        # append rescued taxa (tiny counts; keep them at tail)
                        for sp, hi, any_c, mass, hi_mass_value, _, _, score in keep:
                            kept_items.append((sp, pseudo))
                            added.append(
                                (
                                    sp,
                                    int(hi),
                                    int(any_c),
                                    float(mass),
                                    float(hi_mass_value),
                                    float(score),
                                )
                            )
                        # re-sort non-unclassified part by count desc to keep output stable
                        kept_items.sort(key=lambda it: it[1], reverse=True)
                        filtered_items = kept_items + ([uncls_item] if uncls_item else [])

                print(
                    "[profile][auto] headheavy_posttopk_newgenus_rescue="
                    f"enabled=1 prob_hi={prob_hi:.2f} "
                    f"support_min={support_min} "
                    f"support_hi_frac_min={float(support_hi_frac_min):.4f} "
                    f"expansion_any_max={int(max(float(support_min), float(support_min) / max(float(prob_hi), 1e-9)))} "
                    f"max_added={max_added} "
                    f"added={len(added)}"
                )
            if tail_rich_strength > 1e-9:
                # Tail-rich stage: trim long tail by head mass (exclude UNCLASSIFIED).
                # This suppresses hundreds of tiny non-zero species that destroy
                # presence precision, without affecting per-read assignments.
                if policy_species_head_mass is not None:
                    denom = eff_total if eff_total > 0 else total_count
                    kept = set()
                    kept_head = set()
                    cum = 0.0
                    uncls = None
                    ranked_non_uncls: List[Tuple[str, float]] = []
                    for taxon, count in filtered_items:
                        if taxon == UNCLASSIFIED:
                            uncls = (taxon, count)
                            continue
                        ranked_non_uncls.append((taxon, float(count)))

                    for taxon, count in ranked_non_uncls:
                        rel = (count / denom) * 100 if denom > 0 else 0.0
                        kept_head.add(taxon)
                        kept.add(taxon)
                        cum += rel
                        if cum >= float(policy_species_head_mass):
                            break

                    count_by_species = {taxon: count for taxon, count in ranked_non_uncls}
                    support_species = support_by_level.get("species") if support_by_level else None
                    support_k = int(max(1, int(TAIL_RICH_SPECIES_SUPPORT_MIN)))

                    # Union with top1 support>=K to keep repeated low-abundance taxa
                    # while filtering scatter-type singletons.
                    if (
                        TAIL_RICH_SPECIES_SUPPORT_MIN
                        and TAIL_RICH_SPECIES_SUPPORT_MIN > 1
                        and support_species
                    ):
                        support_gate = max(0.0, float(support_k) - float(PROB_HI_BASE))
                        for taxon, sup in support_species.items():
                            if (
                                taxon != UNCLASSIFIED
                                and (float(sup) * float(tail_rich_strength))
                                >= support_gate
                            ):
                                kept.add(taxon)

                    prob_hi = float(PROB_HI_BASE)
                    post_any_support: Counter[str] = Counter()
                    post_hi_support: Counter[str] = Counter()
                    post_prob_mass: Dict[str, float] = {}
                    post_hi_mass: Dict[str, float] = {}
                    if (
                        species_rel_cutoff > 0
                        and ranked_non_uncls
                    ):
                        (
                            post_any_support,
                            post_hi_support,
                            post_prob_mass,
                            post_hi_mass,
                        ) = _collect_post_topk_species_stats(
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

                    def _post_score(sp: str) -> Tuple[float, int, int, float]:
                        post_hi = int(post_hi_support.get(sp, 0))
                        post_any = int(post_any_support.get(sp, 0))
                        hi_mass = float(post_hi_mass.get(sp, 0.0))
                        frac = float(post_hi) / float(max(1, post_any))
                        score = float(hi_mass) * math.sqrt(max(frac, 1e-9))
                        return score, post_hi, post_any, hi_mass

                    # Tail-rich cutoff-boundary mass swap (1-in-1-out):
                    # move a limited mass from weak >=cutoff taxa to strong <cutoff taxa.
                    # Invariants:
                    # - final rel cutoff is unchanged
                    # - no cutoff bypass
                    # - swap is mass-conservative and budgeted.
                    swap_pairs: List[Tuple[str, str, float, float, float, float, float, float]] = []
                    if (
                        tail_rich_cutoff_swap_budget > 0
                        and species_rel_cutoff > 0
                        and ranked_non_uncls
                    ):
                        cutoff = float(species_rel_cutoff)
                        boundary_ratio = float(prob_hi)
                        promote_lo = cutoff * max(0.0, 1.0 - boundary_ratio)
                        demote_hi = cutoff * (1.0 + boundary_ratio)
                        eps_util = max(cutoff * 1e-6, 1e-9)
                        eps_cross = max(cutoff * 1e-4, 1e-8)
                        promote_pool: List[Tuple[str, float, int, int, int, float, float, str, float, float]] = []
                        demote_pool: List[Tuple[str, float, int, int, int, float, float, str, float, float]] = []
                        for taxon, _ in ranked_non_uncls:
                            rel = (
                                (float(count_by_species.get(taxon, 0.0)) / float(denom)) * 100.0
                                if denom > 0
                                else 0.0
                            )
                            top1_sup = int(support_species.get(taxon, 0)) if support_species else 0
                            score, post_hi, post_any, hi_mass = _post_score(taxon)
                            genus_name = _infer_genus(taxon) or ""
                            win_rate = (
                                float(top1_sup) / float(max(1, int(post_any)))
                                if int(post_any) > 0
                                else 0.0
                            )
                            reliability = (
                                float(hi_mass) / float(max(1e-12, float(post_prob_mass.get(taxon, 0.0))))
                                if float(post_prob_mass.get(taxon, 0.0)) > 0
                                else 0.0
                            )
                            if (
                                rel >= promote_lo
                                and rel < cutoff
                                and post_hi >= support_k
                            ):
                                promote_pool.append(
                                    (
                                        taxon,
                                        rel,
                                        top1_sup,
                                        post_hi,
                                        post_any,
                                        hi_mass,
                                        score,
                                        genus_name,
                                        win_rate,
                                        reliability,
                                    )
                                )
                            if taxon in kept and rel >= cutoff and rel <= demote_hi:
                                demote_pool.append(
                                    (
                                        taxon,
                                        rel,
                                        top1_sup,
                                        post_hi,
                                        post_any,
                                        hi_mass,
                                        score,
                                        genus_name,
                                        win_rate,
                                        reliability,
                                    )
                                )
                        promote_pool.sort(
                            key=lambda it: (it[6], it[5], it[3], it[8], it[9], it[1]),
                            reverse=True,
                        )
                        demote_pool.sort(
                            key=lambda it: (it[6], it[5], it[3], it[8], it[9], it[1]),
                        )

                        feasible_pairs: List[Tuple[float, str, str, float, float, float, float, float, float, str]] = []
                        for (
                            p_taxon,
                            p_rel0,
                            _p_top1_sup,
                            _p_post_hi,
                            _p_post_any,
                            _p_hi_mass,
                            p_score,
                            p_genus,
                            _p_win_rate,
                            _p_reliability,
                        ) in promote_pool:
                            rel_p = float(p_rel0)
                            if rel_p < promote_lo or rel_p >= cutoff:
                                continue
                            need_rel_for_p = (cutoff - rel_p) + eps_cross
                            # Promote utility reward: stronger evidence + closer to cutoff.
                            e_promote = float(p_score) / (need_rel_for_p + eps_util)

                            for (
                                d_taxon,
                                d_rel0,
                                _d_top1_sup,
                                _d_post_hi,
                                _d_post_any,
                                _d_hi_mass,
                                d_score,
                                _d_genus,
                                _d_win_rate,
                                _d_reliability,
                            ) in demote_pool:
                                if d_taxon == p_taxon:
                                    continue
                                rel_d = float(d_rel0)
                                if rel_d < cutoff or rel_d > demote_hi:
                                    continue
                                # Donor cost: penalize stronger evidence and farther-above-cutoff donors.
                                donor_cost = float(d_score) * (
                                    1.0 + ((rel_d - cutoff) / max(cutoff, eps_util))
                                )
                                need_rel_for_d = max(0.0, rel_d - (cutoff - eps_cross))
                                delta_rel = max(need_rel_for_p, need_rel_for_d)
                                delta_count = (delta_rel / 100.0) * float(denom)
                                if float(count_by_species.get(d_taxon, 0.0)) <= delta_count:
                                    continue
                                utility = (e_promote - donor_cost) / (delta_rel + eps_util)
                                if utility <= 0:
                                    continue
                                feasible_pairs.append(
                                    (
                                        float(utility),
                                        p_taxon,
                                        d_taxon,
                                        float(delta_rel),
                                        float(delta_count),
                                        float(p_score),
                                        float(d_score),
                                        float(rel_p),
                                        float(rel_d),
                                        p_genus,
                                    )
                                )

                        feasible_pairs.sort(key=lambda it: it[0], reverse=True)

                        selected_pairs: List[Tuple[float, str, str, float, float, float, float, float, float, str]] = []
                        used_promote: set[str] = set()
                        used_donor: set[str] = set()
                        used_genus: set[str] = set()
                        for pair in feasible_pairs:
                            (
                                utility,
                                p_taxon,
                                d_taxon,
                                delta_rel,
                                delta_count,
                                p_score,
                                d_score,
                                rel_p,
                                rel_d,
                                p_genus,
                            ) = pair
                            if len(selected_pairs) >= int(tail_rich_cutoff_swap_budget):
                                break
                            if p_taxon in used_promote or d_taxon in used_donor:
                                continue
                            if p_genus and p_genus in used_genus:
                                continue
                            selected_pairs.append(pair)
                            used_promote.add(p_taxon)
                            used_donor.add(d_taxon)
                            if p_genus:
                                used_genus.add(p_genus)

                        # Knee stop on utility sequence: avoid low-yield trailing swaps.
                        keep_n = len(selected_pairs)
                        if len(selected_pairs) > 1:
                            best_drop = 0.0
                            best_idx = -1
                            for i in range(len(selected_pairs) - 1):
                                drop = float(selected_pairs[i][0]) - float(selected_pairs[i + 1][0])
                                if drop > best_drop:
                                    best_drop = drop
                                    best_idx = i
                            if best_idx >= 0 and (best_idx + 1) < keep_n:
                                keep_n = best_idx + 1

                        executed_pairs = selected_pairs[:keep_n]

                        for (
                            utility,
                            p_taxon,
                            d_taxon,
                            delta_rel,
                            delta_count,
                            p_score,
                            d_score,
                            rel_p,
                            rel_d,
                            _p_genus,
                        ) in executed_pairs:
                            count_by_species[p_taxon] = float(
                                count_by_species.get(p_taxon, 0.0)
                            ) + float(delta_count)
                            count_by_species[d_taxon] = float(
                                count_by_species.get(d_taxon, 0.0)
                            ) - float(delta_count)
                            kept.add(p_taxon)
                            swap_pairs.append(
                                (
                                    p_taxon,
                                    d_taxon,
                                    float(delta_rel),
                                    float(utility),
                                    float(p_score),
                                    float(d_score),
                                    float(rel_p),
                                    float(rel_d),
                                )
                            )

                        print(
                            "[profile][auto] tailrich_cutoff_swap="
                            f"enabled=1 budget={int(tail_rich_cutoff_swap_budget)} "
                            f"swaps={len(swap_pairs)}"
                        )
                    elif tail_rich_cutoff_swap_budget > 0 and species_rel_cutoff > 0:
                        print(
                            "[profile][auto] tailrich_cutoff_swap="
                            f"enabled=0 reason=no_candidates "
                            f"budget={int(tail_rich_cutoff_swap_budget)}"
                        )

                    if kept:
                        rebuilt = []
                        for taxon, count in items:
                            if taxon == UNCLASSIFIED:
                                if uncls is not None:
                                    rebuilt.append(uncls)
                                continue
                            if taxon in kept:
                                rebuilt.append(
                                    (
                                        taxon,
                                        float(count_by_species.get(taxon, float(count))),
                                    )
                                )
                        if uncls is not None and not rebuilt:
                            rebuilt.append(uncls)
                        if rebuilt:
                            filtered_items = rebuilt

                # Final hard cutoff on species relative abundance (tail-rich weighted).
                if species_rel_cutoff > 0 and filtered_items:
                    denom = eff_total if eff_total > 0 else total_count
                    final_items: List[Tuple[str, float]] = []
                    for taxon, count in filtered_items:
                        if taxon == UNCLASSIFIED:
                            final_items.append((taxon, count))
                            continue
                        rel = (count / denom) * 100 if denom > 0 else 0.0
                        if rel >= float(species_rel_cutoff):
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
