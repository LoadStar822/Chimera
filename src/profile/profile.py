import heapq
import math
import warnings
from collections import Counter
from dataclasses import dataclass, field
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
# Tail 截断：以百分比计，0.01 = 0.01%
MIN_REL_ABUNDANCE = 0.01
# 最小 evidence：posterior_evidence 票数下限，过滤低支持的长尾假阳性
MIN_EVIDENCE = 10
# 低多样性样本（典型：MOCK）更激进的 species 截断：减少长尾假阳性。
# 低多样性判定不应直接用全量 Shannon（长尾会抬高 Shannon），这里辅以 head mass 触发。
LOW_DIVERSITY_SHANNON = 2.2
LOW_DIVERSITY_TOP_K = 10
LOW_DIVERSITY_TOP_MASS = 90.0
LOW_DIVERSITY_SPECIES_MIN_REL_ABUNDANCE = 0.12
# 高多样性样本：presence 对“非零长尾物种数”极敏感。
# 用 head mass 约束输出物种集合规模（仅用于 species 层级）。
HIGH_DIVERSITY_SPECIES_HEAD_MASS = 99.5
# 高多样性样本：更严格的最小相对丰度阈值（%），用于直接切掉尾巴物种。
HIGH_DIVERSITY_MIN_REL_ABUNDANCE = 0.05
# 高多样性样本：最终写出时的 species 相对丰度硬阈值（%）。
HIGH_DIVERSITY_FINAL_REL_CUTOFF = 0.02
# 高多样性样本：species presence 的主要误差常常来自“属内多物种混淆”
# （oracle_genus 上界高），因此在 species 层级对每个 genus 做二次收敛。
HIGH_DIVERSITY_GENUS_MODE = "none"  # none/headmass/topn_ratio
HIGH_DIVERSITY_GENUS_HEAD_MASS = 90.0
HIGH_DIVERSITY_GENUS_MAX_SPECIES = 3
HIGH_DIVERSITY_GENUS_TOPN = 3
HIGH_DIVERSITY_GENUS_REL_MIN_IN_GENUS = 10.0  # 百分比
# 高多样性样本：散射型 top1 FP 往往是 singleton/doubleton。
# 用 top1 的“出现次数”做第二信号：只要 support 足够（>=K）就保留，
# 以便在提高 head mass（救 recall）时仍能压住 FP 物种种类数。
HIGH_DIVERSITY_SPECIES_SUPPORT_MIN = 3
# 高多样性样本：tail gate（仅在 species 层级生效）。
# 低丰度尾巴需要 top1 或高 posterior 支撑才能保留。
HIGH_DIVERSITY_TAIL_REL_KEEP = 0.0
HIGH_DIVERSITY_TAIL_TOP1_MIN = 0
HIGH_DIVERSITY_TAIL_HIP_PROB = 0.90
HIGH_DIVERSITY_TAIL_HIP_MIN = 0
HIGH_DIVERSITY_TAIL_HIP_TOPK = 0
# 高多样性样本：profile-only 的最后一搏信号。
# 目标：利用 support–mass 的“形状特征”区分 TP / FP，而不是仅靠 abundance 排序删尾。
# 注意：该模块只影响 high-diversity 的 species 层级，low-diversity 分支完全不受影响。
HIGH_DIVERSITY_SHAPE_MODE = "none"  # none/global_headmass/genus_headmass
HIGH_DIVERSITY_SHAPE_ALPHA = 0.0
HIGH_DIVERSITY_SHAPE_BETA = 0.0
HIGH_DIVERSITY_SHAPE_GAMMA = 0.0
HIGH_DIVERSITY_SHAPE_DELTA = 0.0
HIGH_DIVERSITY_SHAPE_TOPK = 8
HIGH_DIVERSITY_SHAPE_GLOBAL_HEAD_MASS = 97.0
HIGH_DIVERSITY_SHAPE_GENUS_HEAD_MASS = 90.0
HIGH_DIVERSITY_SHAPE_GENUS_MAX_SPECIES = 3
# 高多样性样本：属内 MAP 纠错（需要 classify 输出 POST_TOPK=...）。
# 默认开启；仅在 high-diversity 的 species 层使用，不影响 low-diversity 分支。
HIGH_DIVERSITY_INTRAGENUS_MAP = True
HIGH_DIVERSITY_INTRAGENUS_TOPK = 5
HIGH_DIVERSITY_INTRAGENUS_BETA = 1.0
HIGH_DIVERSITY_INTRAGENUS_GAP = 0.10
# 高多样性样本：genus-local split（两遍扫描，基于 POST_TOPK 约束属内物种集合）。
HIGH_DIVERSITY_GENUS_LOCAL_SPLIT = True
HIGH_DIVERSITY_GENUS_LOCAL_THRES = 0.70
HIGH_DIVERSITY_GENUS_LOCAL_ALPHA = 0.90
HIGH_DIVERSITY_GENUS_LOCAL_TOPK = 3
HIGH_DIVERSITY_GENUS_LOCAL_POST_TOPK = 16
# Toggle post_topk abundance for high-diversity samples.
HIGH_DIVERSITY_USE_POST_TOPK = True
# High-diversity: use only top-K posterior candidates for post_topk abundance
# to avoid long-tail explosion.
POST_TOPK_ABUND_K = 8
# post_topk: if top1-top2 gap is large, collapse to top1 to suppress tail.
POST_TOPK_GAP_HARD = 0.20


@dataclass
class _Top1WeightStats:
    support: int = 0
    mass: float = 0.0
    max_w: float = 0.0
    topk_heap: List[float] = field(default_factory=list)  # min-heap
    topk_sum: float = 0.0

    def _push_topk(self, w: float, topk: int) -> None:
        if topk <= 0:
            return
        if len(self.topk_heap) < topk:
            heapq.heappush(self.topk_heap, w)
            self.topk_sum += w
            return
        if self.topk_heap and w > self.topk_heap[0]:
            smallest = heapq.heapreplace(self.topk_heap, w)
            self.topk_sum += w - smallest

    def add(self, w: float, topk: int) -> None:
        self.support += 1
        self.mass += w
        if w > self.max_w:
            self.max_w = w
        self._push_topk(w, topk=topk)


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


def _collect_top1_weight_stats(
    input_files: Iterable[str], expect_numeric: bool, topk: int
) -> Tuple[Dict[str, _Top1WeightStats], int]:
    """Collect per-taxid top1 support and weight-shape stats from ChimeraClassify.tsv.

    This is designed to be robust to:
    - multiple tokens per line (we pick the max-count token as top1)
    - extra POST_* tokens (ignored)
    - 'unclassified' lines
    """
    stats: Dict[str, _Top1WeightStats] = {}
    unclassified = 0

    for input_file in input_files:
        with open(input_file, "r", errors="ignore") as infile:
            for raw in infile:
                line = raw.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue

                best_tid: Optional[str] = None
                best_c = -1.0
                for token in parts[1:]:
                    token = token.strip()
                    if not token:
                        continue
                    if token.startswith("POST_"):
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
                    if c > best_c:
                        best_c = c
                        best_tid = tid_token

                if best_tid is None:
                    unclassified += 1
                    continue

                st = stats.get(best_tid)
                if st is None:
                    st = _Top1WeightStats()
                    stats[best_tid] = st
                st.add(best_c, topk=topk)

    return stats, unclassified


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


def _build_species_shape_stats(
    *,
    resolved_kind: str,
    tax: Optional[NCBITaxa],
    gtdb_taxonomy: Optional[GtdbTaxonomy],
    top1_stats: Dict[str, _Top1WeightStats],
    topk: int,
) -> Dict[str, _Top1WeightStats]:
    """Aggregate per-taxid top1 stats into per-species stats (keyed by display name)."""
    species_stats: Dict[str, _Top1WeightStats] = {}
    if resolved_kind == "gtdb":
        if gtdb_taxonomy is None:
            return species_stats
        cache: Dict[str, Optional[str]] = {}
        for node, st in top1_stats.items():
            sp = _gtdb_species_name_for_node(gtdb_taxonomy, node, cache)
            if not sp:
                continue
            agg = species_stats.get(sp)
            if agg is None:
                agg = _Top1WeightStats()
                species_stats[sp] = agg
            agg.support += st.support
            agg.mass += st.mass
            if st.max_w > agg.max_w:
                agg.max_w = st.max_w
            # merge topk heaps
            for w in st.topk_heap:
                agg._push_topk(w, topk=topk)
        return species_stats

    # NCBI
    if tax is None:
        return species_stats
    cache_i: Dict[int, Optional[str]] = {}
    for tid_str, st in top1_stats.items():
        try:
            tid = int(tid_str)
        except ValueError:
            continue
        sp = _ncbi_species_name_for_taxid(tax, tid, cache_i)
        if not sp:
            continue
        agg = species_stats.get(sp)
        if agg is None:
            agg = _Top1WeightStats()
            species_stats[sp] = agg
        agg.support += st.support
        agg.mass += st.mass
        if st.max_w > agg.max_w:
            agg.max_w = st.max_w
        for w in st.topk_heap:
            agg._push_topk(w, topk=topk)
    return species_stats


def _shape_score(
    *,
    mass: float,
    support: int,
    topk_share: float,
    share_support_in_genus: float,
    alpha: float,
    beta: float,
    gamma: float,
    delta: float,
) -> float:
    eps = 1e-12
    m = max(mass, eps)
    n = max(float(support), eps)
    a = max(m / n, eps)
    ss = max(share_support_in_genus, eps)
    # Score = log M + α log N − β log(M/N) − γ*skew + δ log(share_support+eps)
    return (
        math.log(m)
        + float(alpha) * math.log(n)
        - float(beta) * math.log(a)
        - float(gamma) * float(topk_share)
        + float(delta) * math.log(ss)
    )


def _apply_hd_shape_filter(
    *,
    filtered_items: List[Tuple[str, float]],
    species_stats: Dict[str, _Top1WeightStats],
    mode: str,
    denom: float,
    alpha: float,
    beta: float,
    gamma: float,
    delta: float,
    global_head_mass: float,
    genus_head_mass: float,
    genus_max_species: int,
) -> List[Tuple[str, float]]:
    """Apply an extra high-diversity species filter based on support–mass shape score."""
    normalized = (mode or "none").strip().lower()
    if normalized in {"off", "no", "false"}:
        normalized = "none"
    if normalized == "none":
        return filtered_items
    if normalized not in {"global_headmass", "genus_headmass"}:
        raise ValueError(f"Unknown hd_shape_mode: {mode}")

    uncls = None
    candidates: List[Tuple[str, float]] = []
    for taxon, count in filtered_items:
        if taxon == UNCLASSIFIED:
            uncls = (taxon, count)
            continue
        candidates.append((taxon, float(count)))
    if not candidates:
        return filtered_items

    # genus totals (support-based) for share_support
    genus_support: Dict[str, int] = {}
    for taxon, _ in candidates:
        genus = _infer_genus(taxon)
        if genus is None:
            continue
        st = species_stats.get(taxon)
        sup = int(st.support) if st is not None else 0
        genus_support[genus] = genus_support.get(genus, 0) + sup

    # compute scores
    scores: Dict[str, float] = {}
    for taxon, count in candidates:
        st = species_stats.get(taxon)
        sup = int(st.support) if st is not None else 0
        # Fall back to 1 to avoid pathological log(0) and to keep behavior stable
        # when support stats are missing (e.g. taxonomy lookup failure).
        if sup <= 0:
            sup = 1
        topk_share = 0.0
        if st is not None and st.mass > 0:
            topk_share = float(st.topk_sum) / float(st.mass)
        genus = _infer_genus(taxon)
        if genus is None:
            ss = 1.0
        else:
            gs = genus_support.get(genus, 0)
            ss = (sup / gs) if gs > 0 else 1.0
        scores[taxon] = _shape_score(
            mass=float(count),
            support=sup,
            topk_share=topk_share,
            share_support_in_genus=ss,
            alpha=alpha,
            beta=beta,
            gamma=gamma,
            delta=delta,
        )

    kept: set[str] = set()
    if normalized == "global_headmass":
        ordered = sorted(candidates, key=lambda it: scores.get(it[0], float("-inf")), reverse=True)
        cum = 0.0
        for taxon, count in ordered:
            kept.add(taxon)
            if denom > 0:
                cum += (count / denom) * 100.0
            if cum >= float(global_head_mass):
                break
    else:  # genus_headmass
        max_n = int(genus_max_species)
        by_genus: Dict[Optional[str], List[Tuple[str, float]]] = {}
        for taxon, count in candidates:
            genus = _infer_genus(taxon)
            by_genus.setdefault(genus, []).append((taxon, count))

        for genus, items in by_genus.items():
            if genus is None:
                for taxon, _ in items:
                    kept.add(taxon)
                continue
            items.sort(key=lambda it: scores.get(it[0], float("-inf")), reverse=True)
            genus_total = sum(c for _, c in items)
            cum = 0.0
            selected: List[str] = []
            for taxon, count in items:
                selected.append(taxon)
                if genus_total > 0:
                    cum += (count / genus_total) * 100.0
                if len(selected) >= max_n:
                    break
                if cum >= float(genus_head_mass):
                    break
            if items and not selected:
                selected = [items[0][0]]
            kept.update(selected)

    rebuilt: List[Tuple[str, float]] = []
    for taxon, count in filtered_items:
        if taxon == UNCLASSIFIED:
            continue
        if taxon in kept:
            rebuilt.append((taxon, count))
    if uncls is not None:
        rebuilt.append(uncls)
    return rebuilt or filtered_items


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


def _collect_post_topk_species_support(
    *,
    input_files: Iterable[str],
    expect_numeric: bool,
    resolved_kind: str,
    tax: Optional[NCBITaxa],
    gtdb_taxonomy: Optional[GtdbTaxonomy],
    prob_min: float,
    topk: int,
) -> Counter[str]:
    """Count per-sequence high-confidence POST_TOPK support at species level.

    For each sequence, if a species appears among POST_TOPK candidates with
    posterior >= prob_min, that species gets +1 support for this sequence.
    """
    species_support: Counter[str] = Counter()
    if prob_min <= 0:
        return species_support

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
                seen: set[str] = set()
                for tid, p in cand:
                    if p < prob_min:
                        continue
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
                    if sp:
                        seen.add(sp)
                for sp in seen:
                    species_support[sp] += 1
    return species_support


def _intragenus_map_species_counter(
    *,
    input_files: Iterable[str],
    expect_numeric: bool,
    resolved_kind: str,
    tax: Optional[NCBITaxa],
    gtdb_taxonomy: Optional[GtdbTaxonomy],
    prior_species_counter: Counter[str],
    topk: int,
    beta: float,
    gap: float,
) -> Counter[str]:
    """MAP correct top1 within genus using POST_TOPK posterior candidates.

    This only uses information present in ChimeraClassify.tsv:
    - top1 taxid:count (used as the mass added per sequence)
    - POST_TOPK=tid:posterior,... (used as read-level posterior evidence)

    The correction is restricted to the genus of the current top1 species.
    """
    # genus prior π(species|genus) from current (uncorrected) profile
    genus_total: Dict[str, float] = {}
    for sp, c in prior_species_counter.items():
        if sp == UNCLASSIFIED:
            continue
        gname = _infer_genus(sp)
        if gname is None:
            continue
        genus_total[gname] = genus_total.get(gname, 0.0) + float(c)

    genus_prior: Dict[str, Dict[str, float]] = {}
    for sp, c in prior_species_counter.items():
        if sp == UNCLASSIFIED:
            continue
        gname = _infer_genus(sp)
        if gname is None:
            continue
        total = genus_total.get(gname, 0.0)
        if total <= 0:
            continue
        genus_prior.setdefault(gname, {})[sp] = float(c) / total

    # taxid -> species display name cache
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

    out: Counter[str] = Counter()
    total_reads = 0
    with_primary = 0
    with_post_topk = 0
    corrected = 0
    gap_kept = 0
    no_genus = 0
    no_post_topk = 0
    same_genus_zero = 0
    same_genus_one = 0
    same_genus_multi = 0
    gap_values: List[float] = []
    top1_post_values: List[float] = []
    top2_post_values: List[float] = []
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
                    out[UNCLASSIFIED] += 1
                    continue
                top_tid, mass = primary
                with_primary += 1

                base_sp = taxid_to_species(top_tid) or UNCLASSIFIED
                if base_sp == UNCLASSIFIED:
                    out[UNCLASSIFIED] += 1
                    continue

                base_genus = _infer_genus(base_sp)
                if base_genus is None:
                    no_genus += 1
                    out[base_sp] += mass
                    continue

                cand = _parse_post_topk_from_tokens(tokens, expect_numeric=expect_numeric, topk=topk)
                if not cand:
                    no_post_topk += 1
                    out[base_sp] += mass
                    continue
                with_post_topk += 1

                # collapse candidates to species within the same genus (keep max posterior)
                sp_post: Dict[str, float] = {}
                for tid, prob in cand:
                    sp = taxid_to_species(tid)
                    if not sp:
                        continue
                    if _infer_genus(sp) != base_genus:
                        continue
                    prev = sp_post.get(sp)
                    if prev is None or prob > prev:
                        sp_post[sp] = prob

                if not sp_post:
                    same_genus_zero += 1
                    out[base_sp] += mass
                    continue

                ordered = sorted(sp_post.items(), key=lambda kv: kv[1], reverse=True)
                if len(ordered) >= 2:
                    gval = ordered[0][1] - ordered[1][1]
                    same_genus_multi += 1
                    gap_values.append(gval)
                    top1_post_values.append(ordered[0][1])
                    top2_post_values.append(ordered[1][1])
                else:
                    gval = 1.0
                    same_genus_one += 1
                    top1_post_values.append(ordered[0][1])

                # Only correct low-confidence cases (gap small); otherwise keep original.
                if gval >= float(gap):
                    gap_kept += 1
                    out[base_sp] += mass
                    continue

                pri = genus_prior.get(base_genus, {})
                eps = 1e-12
                best_sp = base_sp
                best_score = float("-inf")
                for sp, prob in ordered:
                    prior = pri.get(sp, eps)
                    score = math.log(max(prob, eps)) + float(beta) * math.log(max(prior, eps))
                    if score > best_score:
                        best_score = score
                        best_sp = sp

                if best_sp != base_sp:
                    corrected += 1
                out[best_sp] += mass
    def pct(p: int, q: int) -> float:
        return (100.0 * p / q) if q > 0 else 0.0

    def summarize(vals: List[float]) -> str:
        if not vals:
            return "n=0"
        vs = sorted(vals)
        n = len(vs)
        def pick(q: float) -> float:
            idx = int(round(q * (n - 1)))
            return vs[idx]
        return (
            f"n={n} p50={pick(0.50):.4f} p90={pick(0.90):.4f} "
            f"p99={pick(0.99):.4f}"
        )

    print(
        "Intragenus MAP: total="
        f"{total_reads} primary={with_primary} post_topk={with_post_topk} "
        f"corrected={corrected} gap_kept={gap_kept} "
        f"no_post_topk={no_post_topk} no_genus={no_genus}"
    )
    print(
        "Intragenus MAP: same_genus="
        f"zero={same_genus_zero} ({pct(same_genus_zero, with_post_topk):.2f}%) "
        f"one={same_genus_one} ({pct(same_genus_one, with_post_topk):.2f}%) "
        f"multi={same_genus_multi} ({pct(same_genus_multi, with_post_topk):.2f}%)"
    )
    print(f"Intragenus MAP: gap {summarize(gap_values)}")
    print(f"Intragenus MAP: top1_post {summarize(top1_post_values)}")
    print(f"Intragenus MAP: top2_post {summarize(top2_post_values)}")
    return out


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
    elif mode in {"post_topk", "post"}:
        mode = "post_topk"
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
                elif mode == "post_topk":
                    # Allocate per-sequence mass (sum of taxidCount counts) to POST_TOPK posterior candidates.
                    # This uses the dumped posterior probabilities (after EM) and can retain
                    # low-abundance true taxa that would otherwise be dropped by hard top1.
                    mass_total = 0.0
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
                        mass_total += c
                    if not has_taxid or mass_total <= 0:
                        unclassified_reads += 1
                        continue

                    cand = _parse_post_topk_from_tokens(
                        [t.strip() for t in tokens if t.strip()],
                        expect_numeric=expect_numeric,
                        topk=POST_TOPK_ABUND_K,
                    )
                    if not cand:
                        # Fallback: behave like soft_seq when POST_TOPK is missing.
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
                            taxid_counts[tid_token] += c
                        continue

                    # Optionally collapse to top1 when posterior gap is large.
                    cand_sorted = sorted(cand, key=lambda x: x[1], reverse=True)
                    if cand_sorted:
                        top1 = cand_sorted[0][1]
                        top2 = cand_sorted[1][1] if len(cand_sorted) > 1 else 0.0
                        if POST_TOPK_GAP_HARD > 0 and (top1 - top2) >= POST_TOPK_GAP_HARD:
                            taxid_counts[cand_sorted[0][0]] += mass_total
                            continue

                    sum_prob = sum(p for _, p in cand_sorted)
                    if sum_prob <= 0:
                        unclassified_reads += 1
                        continue
                    for tid, p in cand_sorted:
                        taxid_counts[tid] += mass_total * (p / sum_prob)
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
    hd_species_head_mass: Optional[float] = HIGH_DIVERSITY_SPECIES_HEAD_MASS,
    hd_genus_mode: str = HIGH_DIVERSITY_GENUS_MODE,
    hd_genus_head_mass: float = HIGH_DIVERSITY_GENUS_HEAD_MASS,
    hd_genus_max_species: int = HIGH_DIVERSITY_GENUS_MAX_SPECIES,
    hd_genus_topn: int = HIGH_DIVERSITY_GENUS_TOPN,
    hd_genus_rel_min_in_genus: float = HIGH_DIVERSITY_GENUS_REL_MIN_IN_GENUS,
    hd_shape_mode: str = HIGH_DIVERSITY_SHAPE_MODE,
    hd_shape_alpha: float = HIGH_DIVERSITY_SHAPE_ALPHA,
    hd_shape_beta: float = HIGH_DIVERSITY_SHAPE_BETA,
    hd_shape_gamma: float = HIGH_DIVERSITY_SHAPE_GAMMA,
    hd_shape_delta: float = HIGH_DIVERSITY_SHAPE_DELTA,
    hd_shape_topk: int = HIGH_DIVERSITY_SHAPE_TOPK,
    hd_shape_global_head_mass: float = HIGH_DIVERSITY_SHAPE_GLOBAL_HEAD_MASS,
    hd_shape_genus_head_mass: float = HIGH_DIVERSITY_SHAPE_GENUS_HEAD_MASS,
    hd_shape_genus_max_species: int = HIGH_DIVERSITY_SHAPE_GENUS_MAX_SPECIES,
    hd_intragenus_map: bool = HIGH_DIVERSITY_INTRAGENUS_MAP,
    hd_intragenus_topk: int = HIGH_DIVERSITY_INTRAGENUS_TOPK,
    hd_intragenus_beta: float = HIGH_DIVERSITY_INTRAGENUS_BETA,
    hd_intragenus_gap: float = HIGH_DIVERSITY_INTRAGENUS_GAP,
    tail_rel_keep: float = HIGH_DIVERSITY_TAIL_REL_KEEP,
    tail_top1_min: int = HIGH_DIVERSITY_TAIL_TOP1_MIN,
    tail_hip_prob: float = HIGH_DIVERSITY_TAIL_HIP_PROB,
    tail_hip_min: int = HIGH_DIVERSITY_TAIL_HIP_MIN,
    tail_hip_topk: int = HIGH_DIVERSITY_TAIL_HIP_TOPK,
) -> None:
    input_files = list(input_files)
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
    normalized_genus_mode = (hd_genus_mode or "none").strip().lower()
    if normalized_genus_mode in {"off", "no", "false"}:
        normalized_genus_mode = "none"
    if normalized_genus_mode not in {"none", "headmass", "topn_ratio"}:
        raise ValueError(f"Unknown hd_genus_mode: {hd_genus_mode}")
    if hd_species_head_mass is not None and not (0.0 <= float(hd_species_head_mass) <= 100.0):
        raise ValueError("--hd-species-head-mass must be within [0, 100]")
    if not (0.0 <= float(hd_genus_head_mass) <= 100.0):
        raise ValueError("--hd-genus-head-mass must be within [0, 100]")
    if int(hd_genus_max_species) < 1:
        raise ValueError("--hd-genus-max-species must be >= 1")
    if int(hd_genus_topn) < 1:
        raise ValueError("--hd-genus-topn must be >= 1")
    if float(hd_genus_rel_min_in_genus) < 0.0:
        raise ValueError("--hd-genus-rel-min-in-genus must be >= 0")
    normalized_shape_mode = (hd_shape_mode or "none").strip().lower()
    if normalized_shape_mode in {"off", "no", "false"}:
        normalized_shape_mode = "none"
    if normalized_shape_mode not in {"none", "global_headmass", "genus_headmass"}:
        raise ValueError("--hd-shape-mode must be one of: none/global_headmass/genus_headmass")
    if not (0.0 <= float(hd_shape_global_head_mass) <= 100.0):
        raise ValueError("--hd-shape-global-head-mass must be within [0, 100]")
    if not (0.0 <= float(hd_shape_genus_head_mass) <= 100.0):
        raise ValueError("--hd-shape-genus-head-mass must be within [0, 100]")
    if int(hd_shape_genus_max_species) < 1:
        raise ValueError("--hd-shape-genus-max-species must be >= 1")
    if int(hd_shape_topk) < 0:
        raise ValueError("--hd-shape-topk must be >= 0")
    if int(hd_intragenus_topk) < 0:
        raise ValueError("--hd-intragenus-topk must be >= 0")
    if float(hd_intragenus_gap) < 0.0:
        raise ValueError("--hd-intragenus-gap must be >= 0")
    if float(tail_rel_keep) < 0.0:
        raise ValueError("--tail-rel-keep must be >= 0")
    if int(tail_top1_min) < 0:
        raise ValueError("--tail-top1-min must be >= 0")
    if not (0.0 <= float(tail_hip_prob) <= 1.0):
        raise ValueError("--tail-hip-prob must be within [0, 1]")
    if int(tail_hip_min) < 0:
        raise ValueError("--tail-hip-min must be >= 0")
    if int(tail_hip_topk) < 0:
        raise ValueError("--tail-hip-topk must be >= 0")

    taxid_counts, base_unclassified = _collect_taxids(
        input_files, expect_numeric=expect_numeric, abundance_mode=abundance_mode
    )

    tax: Optional[NCBITaxa] = None
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

    tail_gate_enabled = (float(tail_rel_keep) > 0.0) and (
        int(tail_top1_min) > 0 or int(tail_hip_min) > 0
    )

    # Optional: compute top1 support counts (sequence-level) for high-diversity species trimming.
    support_by_level: Optional[Dict[str, Counter[str]]] = None
    if HIGH_DIVERSITY_SPECIES_SUPPORT_MIN and HIGH_DIVERSITY_SPECIES_SUPPORT_MIN > 1:
        # If we already run top1 mode (e.g., low-diversity auto override), support is redundant.
        mode_for_support = (abundance_mode or "soft_seq").strip().lower()
        if mode_for_support not in {"top1"}:
            support_taxid_counts, support_uncls = _collect_taxids(
                input_files, expect_numeric=expect_numeric, abundance_mode="top1"
            )
            if resolved_kind == "gtdb":
                support_by_level = _aggregate_gtdb_levels(
                    gtdb_taxonomy, support_taxid_counts, support_uncls
                )
            else:
                support_numeric: Counter[int] = Counter()
                for taxid_str, count in support_taxid_counts.items():
                    try:
                        tid = int(taxid_str)
                    except ValueError:
                        continue
                    if tid > 0:
                        support_numeric[tid] += count
                support_by_level = _aggregate_levels(tax, support_numeric, support_uncls)

    # Optional: tail gate support signals (high-diversity species only; does not affect low-diversity path).
    tail_top1_support_species: Optional[Counter[str]] = None
    if tail_gate_enabled and int(tail_top1_min) > 0:
        support_taxid_counts, support_uncls = _collect_taxids(
            input_files, expect_numeric=expect_numeric, abundance_mode="top1"
        )
        if resolved_kind == "gtdb":
            tail_support_by_level = _aggregate_gtdb_levels(
                gtdb_taxonomy, support_taxid_counts, support_uncls
            )
        else:
            support_numeric2: Counter[int] = Counter()
            for taxid_str, count in support_taxid_counts.items():
                try:
                    tid = int(taxid_str)
                except ValueError:
                    continue
                if tid > 0:
                    support_numeric2[tid] += count
            tail_support_by_level = _aggregate_levels(tax, support_numeric2, support_uncls)
        tail_top1_support_species = tail_support_by_level.get("species")

    tail_hip_support_species: Optional[Counter[str]] = None
    if tail_gate_enabled and int(tail_hip_min) > 0:
        tail_hip_support_species = _collect_post_topk_species_support(
            input_files=input_files,
            expect_numeric=expect_numeric,
            resolved_kind=resolved_kind,
            tax=tax,
            gtdb_taxonomy=(gtdb_taxonomy if resolved_kind == "gtdb" else None),
            prob_min=float(tail_hip_prob),
            topk=int(tail_hip_topk),
        )

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
    force_low_diversity_shannon: Optional[float] = None
    has_post_topk: Optional[bool] = None
    if (
        normalized_mode == "soft_seq"
        and species_counter
        and species_total > 0
        and sum(species_counter.values()) > 0
    ):
        # Use top1 distribution for low-diversity detection to avoid soft_seq tail inflation.
        detect_counter = species_counter
        detect_total = species_total
        try:
            detect_taxid_counts, detect_uncls = _collect_taxids(
                input_files, expect_numeric=expect_numeric, abundance_mode="top1"
            )
            if resolved_kind == "gtdb":
                detect_by_level = _aggregate_gtdb_levels(
                    gtdb_taxonomy, detect_taxid_counts, detect_uncls
                )
            else:
                detect_numeric: Counter[int] = Counter()
                for taxid_str, count in detect_taxid_counts.items():
                    try:
                        tid = int(taxid_str)
                    except ValueError:
                        continue
                    if tid > 0:
                        detect_numeric[tid] += count
                detect_by_level = _aggregate_levels(tax, detect_numeric, detect_uncls)
            detect_counter = detect_by_level.get("species") or species_counter
            detect_total = sum(detect_counter.values()) if detect_counter else species_total
        except Exception:
            detect_counter = species_counter
            detect_total = species_total

        uncls = detect_counter.get(UNCLASSIFIED, 0)
        eff_total = max(0, detect_total - uncls)
        rels = []
        for taxon, count in detect_counter.items():
            if taxon == UNCLASSIFIED:
                continue
            denom = eff_total if eff_total > 0 else detect_total
            rels.append((count / denom) * 100.0)
        rels.sort(reverse=True)
        top_mass = sum(rels[:LOW_DIVERSITY_TOP_K])
        shannon_index = calculate_shannon_index(detect_counter)
        is_low_diversity = (shannon_index < LOW_DIVERSITY_SHANNON) or (top_mass >= LOW_DIVERSITY_TOP_MASS)
        print(
            f"[profile] low-div detect: shannon={shannon_index:.3f} top_mass={top_mass:.2f} "
            f"species_total={detect_total} is_low_diversity={is_low_diversity}"
        )

        if is_low_diversity:
            force_low_diversity = True
            force_low_diversity_shannon = shannon_index
            if has_post_topk is None:
                has_post_topk = _has_post_topk_tokens(input_files)
            if has_post_topk:
                print("[profile] detected low-diversity: switching abundance_mode to post_topk")
                abundance_mode = "post_topk"
            else:
                print("[profile] detected low-diversity: switching abundance_mode to top1")
                abundance_mode = "top1"
            taxid_counts, base_unclassified = _collect_taxids(
                input_files, expect_numeric=expect_numeric, abundance_mode=abundance_mode
            )
            if resolved_kind == "gtdb":
                count_by_level = _aggregate_gtdb_levels(
                    gtdb_taxonomy, taxid_counts, base_unclassified
                )
            else:
                numeric_counts = Counter()
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
        else:
            # High-diversity: switch to post_topk to preserve secondary evidence without
            # squaring/rounding it away, but limit to top-K to avoid long-tail explosion.
            if HIGH_DIVERSITY_USE_POST_TOPK:
                print(
                    f"[profile] high-diversity: switching abundance_mode to post_topk "
                    f"(topk={POST_TOPK_ABUND_K})"
                )
                abundance_mode = "post_topk"
                taxid_counts, base_unclassified = _collect_taxids(
                    input_files, expect_numeric=expect_numeric, abundance_mode=abundance_mode
                )
                if resolved_kind == "gtdb":
                    count_by_level = _aggregate_gtdb_levels(
                        gtdb_taxonomy, taxid_counts, base_unclassified
                    )
                else:
                    numeric_counts = Counter()
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
        mode = (abundance_mode or "soft_seq").strip().lower()
        if mode in {"soft", "soft_seq", "soft_tax"}:
            mode = "soft_seq"
        elif mode in {"post_topk", "post"}:
            mode = "post_topk"
        elif mode in {"hard", "hard_seq"}:
            mode = "hard_seq"
        elif mode in {"top1"}:
            mode = "top1"
        else:
            mode = "soft_seq"
        outfile.write(f"# abundance_mode={mode}\n")
        count_unit = "sequence" if mode == "top1" else "posterior_evidence"
        outfile.write(f"# count_unit={count_unit}\n")
        if tail_gate_enabled:
            outfile.write(
                f"# tail_gate=on rel_keep={float(tail_rel_keep):.4f}% top1_min={int(tail_top1_min)} hip_prob={float(tail_hip_prob):.4f} hip_min={int(tail_hip_min)} hip_topk={int(tail_hip_topk)}\n"
            )
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

            # Tail 截断：移除丰度极低且证据不足的物种，减少虚假 presence
            prefilter_use_eff = False
            prefilter_denom = total_count
            prefilter_eff_total = total_count
            prefilter_uncls_frac = 0.0
            prefilter_before = 0
            prefilter_rescued = 0
            if level == "species":
                prefilter_uncls = taxon_counter.get(UNCLASSIFIED, 0)
                prefilter_eff_total = max(0, total_count - prefilter_uncls)
                prefilter_uncls_frac = (
                    (prefilter_uncls / total_count) if total_count > 0 else 0.0
                )
                prefilter_use_eff = (not is_low_diversity) or (prefilter_uncls_frac >= 0.05)
                if prefilter_use_eff and prefilter_eff_total > 0:
                    prefilter_denom = prefilter_eff_total

            filtered_items = []
            for taxon, count in items:
                if taxon == UNCLASSIFIED:
                    rel = (count / total_count) * 100 if total_count > 0 else 0.0
                    if rel < MIN_REL_ABUNDANCE and count < MIN_EVIDENCE:
                        continue
                    filtered_items.append((taxon, count))
                    continue

                if level == "species":
                    prefilter_before += 1
                    rel_total = (count / total_count) * 100 if total_count > 0 else 0.0
                    rel_eff = (
                        (count / prefilter_eff_total) * 100
                        if prefilter_eff_total > 0
                        else rel_total
                    )
                    rel = (count / prefilter_denom) * 100 if prefilter_denom > 0 else 0.0
                    if rel < MIN_REL_ABUNDANCE and count < MIN_EVIDENCE:
                        if (
                            prefilter_use_eff
                            and prefilter_eff_total > 0
                            and count < MIN_EVIDENCE
                            and rel_total < MIN_REL_ABUNDANCE
                            and rel_eff >= MIN_REL_ABUNDANCE
                        ):
                            prefilter_rescued += 1
                        continue
                    filtered_items.append((taxon, count))
                    continue

                rel = (count / total_count) * 100 if total_count > 0 else 0.0
                if rel < MIN_REL_ABUNDANCE and count < MIN_EVIDENCE:
                    continue
                filtered_items.append((taxon, count))

            if level == "species" and prefilter_use_eff:
                reasons = []
                if not is_low_diversity:
                    reasons.append("high_diversity")
                if prefilter_uncls_frac >= 0.05:
                    reasons.append(f"uncl_frac({prefilter_uncls_frac:.4f})>=0.05")
                reason_text = "|".join(reasons) if reasons else "none"
                prefilter_after = sum(
                    1 for taxon, _ in filtered_items if taxon != UNCLASSIFIED
                )
                denom_label = (
                    "eff_total"
                    if (prefilter_eff_total > 0 and prefilter_denom == prefilter_eff_total)
                    else "total"
                )
                print(
                    "[profile][auto] species_prefilter_denom="
                    f"{denom_label} enabled=1 reason={reason_text} "
                    f"uncl_frac={prefilter_uncls_frac:.4f} "
                    f"min_rel={MIN_REL_ABUNDANCE} min_ev={MIN_EVIDENCE} "
                    f"denom_total={float(total_count):.4f} "
                    f"denom_eff={float(prefilter_eff_total):.4f} "
                    f"before={prefilter_before} after={prefilter_after} "
                    f"rescued={prefilter_rescued}"
                )

            if level == "species" and filtered_items:
                # Low-diversity detection (robust to long tail):
                # - Shannon on raw distribution (may be inflated by long tail)
                # - Top-K mass on filtered distribution (dominance signal)
                uncls_count = taxon_counter.get(UNCLASSIFIED, 0)
                eff_total = max(0, total_count - uncls_count)
                rels = []
                for taxon, count in filtered_items:
                    if taxon == UNCLASSIFIED:
                        continue
                    denom = eff_total if eff_total > 0 else total_count
                    rel = (count / denom) * 100 if denom > 0 else 0.0
                    rels.append(rel)
                rels.sort(reverse=True)
                top_mass = sum(rels[:LOW_DIVERSITY_TOP_K])
                if force_low_diversity:
                    is_low_diversity = True
                else:
                    is_low_diversity = (shannon_index < LOW_DIVERSITY_SHANNON) or (
                        top_mass >= LOW_DIVERSITY_TOP_MASS
                    )

                if is_low_diversity:
                    low_div_shannon = (
                        force_low_diversity_shannon
                        if force_low_diversity_shannon is not None
                        else shannon_index
                    )
                    cap_n = int(math.exp(low_div_shannon) * 4.0)
                    cap_n = max(16, min(64, cap_n))
                    ranked = [
                        (taxon, count)
                        for taxon, count in filtered_items
                        if taxon != UNCLASSIFIED
                    ]
                    ranked.sort(key=lambda it: it[1], reverse=True)
                    topn = {taxon for taxon, _ in ranked[:cap_n]}
                    aggressive = []
                    for taxon, count in filtered_items:
                        if taxon == UNCLASSIFIED:
                            aggressive.append((taxon, count))
                            continue
                        denom = eff_total if eff_total > 0 else total_count
                        rel = (count / denom) * 100 if denom > 0 else 0.0
                        if taxon in topn or rel >= LOW_DIVERSITY_SPECIES_MIN_REL_ABUNDANCE:
                            aggressive.append((taxon, count))
                    if aggressive:
                        filtered_items = aggressive
                else:
                    # Optional: Intragenus MAP correction (requires POST_TOPK in classify output).
                    # Only for high-diversity samples at species level; does not affect low-diversity path.
                    if hd_intragenus_map:
                        corrected = _intragenus_map_species_counter(
                            input_files=input_files,
                            expect_numeric=expect_numeric,
                            resolved_kind=resolved_kind,
                            tax=tax,
                            gtdb_taxonomy=(gtdb_taxonomy if resolved_kind == "gtdb" else None),
                            prior_species_counter=taxon_counter,
                            topk=int(hd_intragenus_topk),
                            beta=float(hd_intragenus_beta),
                            gap=float(hd_intragenus_gap),
                        )
                        if corrected:
                            taxon_counter = corrected
                            total_count = sum(taxon_counter.values())
                            uncls_count = taxon_counter.get(UNCLASSIFIED, 0)
                            eff_total = max(0, total_count - uncls_count)
                            items = taxon_counter.most_common()
                            filtered_items = []
                            for taxon, count in items:
                                rel = (count / total_count) * 100 if total_count > 0 else 0.0
                                if rel < MIN_REL_ABUNDANCE and count < MIN_EVIDENCE:
                                    continue
                                filtered_items.append((taxon, count))

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
                            filtered_items = []
                            for taxon, count in items:
                                rel = (count / total_count) * 100 if total_count > 0 else 0.0
                                if rel < MIN_REL_ABUNDANCE and count < MIN_EVIDENCE:
                                    continue
                                filtered_items.append((taxon, count))

                    # High-diversity: hard minimum relative abundance cutoff.
                    if HIGH_DIVERSITY_MIN_REL_ABUNDANCE and HIGH_DIVERSITY_MIN_REL_ABUNDANCE > 0:
                        denom = eff_total if eff_total > 0 else total_count
                        pruned_items: List[Tuple[str, float]] = []
                        for taxon, count in filtered_items:
                            if taxon == UNCLASSIFIED:
                                pruned_items.append((taxon, count))
                                continue
                            rel = (count / denom) * 100 if denom > 0 else 0.0
                            if rel >= float(HIGH_DIVERSITY_MIN_REL_ABUNDANCE):
                                pruned_items.append((taxon, count))
                        if pruned_items:
                            filtered_items = pruned_items

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

                        # High-diversity: tail gate based on support signals.
                        # If a species is in the low-abundance tail, require either:
                        # - enough top1 support (sequence-level occurrences), OR
                        # - enough high-posterior support in POST_TOPK (p >= tail_hip_prob).
                        if tail_gate_enabled and filtered_items:
                            denom = eff_total if eff_total > 0 else total_count
                            kept_items: List[Tuple[str, float]] = []
                            dropped = 0
                            for taxon, count in filtered_items:
                                if taxon == UNCLASSIFIED:
                                    kept_items.append((taxon, count))
                                    continue
                                rel = (count / denom) * 100.0 if denom > 0 else 0.0
                                if rel >= float(tail_rel_keep):
                                    kept_items.append((taxon, count))
                                    continue
                                top1_sup = (
                                    int(tail_top1_support_species.get(taxon, 0))
                                    if tail_top1_support_species
                                    else 0
                                )
                                hip_sup = (
                                    int(tail_hip_support_species.get(taxon, 0))
                                    if tail_hip_support_species
                                    else 0
                                )
                                keep = False
                                if int(tail_top1_min) > 0 and top1_sup >= int(tail_top1_min):
                                    keep = True
                                if int(tail_hip_min) > 0 and hip_sup >= int(tail_hip_min):
                                    keep = True
                                if keep:
                                    kept_items.append((taxon, count))
                                else:
                                    dropped += 1
                            if dropped and kept_items:
                                filtered_items = kept_items

                    # High-diversity: optional shape-based filtering (profile-only "last try").
                    if normalized_shape_mode != "none" and filtered_items:
                        # Collect top1 per-sequence shape stats once and map to species display names.
                        top1_stats, _ = _collect_top1_weight_stats(
                            input_files, expect_numeric=expect_numeric, topk=int(hd_shape_topk)
                        )
                        species_shape = _build_species_shape_stats(
                            resolved_kind=resolved_kind,
                            tax=tax,
                            gtdb_taxonomy=(gtdb_taxonomy if resolved_kind == "gtdb" else None),
                            top1_stats=top1_stats,
                            topk=int(hd_shape_topk),
                        )
                        denom = eff_total if eff_total > 0 else total_count
                        filtered_items = _apply_hd_shape_filter(
                            filtered_items=filtered_items,
                            species_stats=species_shape,
                            mode=normalized_shape_mode,
                            denom=float(denom),
                            alpha=float(hd_shape_alpha),
                            beta=float(hd_shape_beta),
                            gamma=float(hd_shape_gamma),
                            delta=float(hd_shape_delta),
                            global_head_mass=float(hd_shape_global_head_mass),
                            genus_head_mass=float(hd_shape_genus_head_mass),
                            genus_max_species=int(hd_shape_genus_max_species),
                        )

                    # Genus-level contraction (high-diversity species only):
                    # Reduce within-genus species confusion by keeping only a small,
                    # high-mass subset of species per genus.
                    if normalized_genus_mode != "none" and filtered_items:
                        max_n = int(hd_genus_max_species)
                        keep_species = set()
                        genus_groups: Dict[str, List[Tuple[str, float]]] = {}
                        for taxon, count in filtered_items:
                            if taxon == UNCLASSIFIED:
                                continue
                            genus = _infer_genus(taxon)
                            if genus is None:
                                keep_species.add(taxon)
                                continue
                            genus_groups.setdefault(genus, []).append((taxon, float(count)))

                        for genus, pairs in genus_groups.items():
                            pairs.sort(key=lambda it: it[1], reverse=True)
                            genus_total = sum(c for _, c in pairs)
                            selected: List[str] = []

                            if normalized_genus_mode == "headmass":
                                cum = 0.0
                                for taxon, count in pairs:
                                    selected.append(taxon)
                                    if genus_total > 0:
                                        cum += (count / genus_total) * 100.0
                                    if len(selected) >= max_n:
                                        break
                                    if cum >= float(hd_genus_head_mass):
                                        break
                            else:  # topn_ratio
                                topn = int(hd_genus_topn)
                                rel_min = float(hd_genus_rel_min_in_genus)
                                selected_set = set()
                                for idx, (taxon, _) in enumerate(pairs):
                                    if idx < topn:
                                        selected_set.add(taxon)
                                if genus_total > 0:
                                    for taxon, count in pairs:
                                        rel = (count / genus_total) * 100.0
                                        if rel >= rel_min:
                                            selected_set.add(taxon)
                                for taxon, _ in pairs:
                                    if taxon in selected_set:
                                        selected.append(taxon)
                                    if len(selected) >= max_n:
                                        break
                                # Ensure top1 is always kept.
                                if pairs:
                                    top1 = pairs[0][0]
                                    if top1 in selected:
                                        selected = [top1] + [t for t in selected if t != top1]
                                    else:
                                        selected = [top1] + selected
                                    selected = selected[:max_n]

                            if pairs and not selected:
                                selected = [pairs[0][0]]
                            keep_species.update(selected)

                        rebuilt = []
                        for taxon, count in filtered_items:
                            if taxon == UNCLASSIFIED or taxon in keep_species:
                                rebuilt.append((taxon, count))
                        if rebuilt:
                            filtered_items = rebuilt

                    # Final hard cutoff on species relative abundance (high-diversity only).
                    if (
                        level == "species"
                        and not is_low_diversity
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
            shannon_index = calculate_shannon_index(Counter(dict(filtered_items)))
            simpson_index = calculate_simpson_index(Counter(dict(filtered_items)))

            for taxon, count in filtered_items:
                relative_abundance = (count / total_filtered) * 100 if total_filtered > 0 else 0.0
                outfile.write(
                    f"{taxon}\t{count}\t{relative_abundance:.2f}\t{shannon_index:.4f}\t{simpson_index:.4f}\n"
                )
