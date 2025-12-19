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
HIGH_DIVERSITY_SPECIES_HEAD_MASS = 97.0
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
HIGH_DIVERSITY_SPECIES_SUPPORT_MIN = 0
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
# 默认关闭；开启后仅在 high-diversity 的 species 层使用，不影响 low-diversity 分支。
HIGH_DIVERSITY_INTRAGENUS_MAP = False
HIGH_DIVERSITY_INTRAGENUS_TOPK = 5
HIGH_DIVERSITY_INTRAGENUS_BETA = 1.0
HIGH_DIVERSITY_INTRAGENUS_GAP = 0.10


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
                if base_sp == UNCLASSIFIED:
                    out[UNCLASSIFIED] += 1
                    continue

                base_genus = _infer_genus(base_sp)
                if base_genus is None:
                    out[base_sp] += mass
                    continue

                cand = _parse_post_topk_from_tokens(tokens, expect_numeric=expect_numeric, topk=topk)
                if not cand:
                    out[base_sp] += mass
                    continue

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
                    out[base_sp] += mass
                    continue

                ordered = sorted(sp_post.items(), key=lambda kv: kv[1], reverse=True)
                if len(ordered) >= 2:
                    gval = ordered[0][1] - ordered[1][1]
                else:
                    gval = 1.0

                # Only correct low-confidence cases (gap small); otherwise keep original.
                if gval >= float(gap):
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

                out[best_sp] += mass
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
                        topk=0,  # consume all dumped candidates
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

                    sum_prob = sum(p for _, p in cand)
                    if sum_prob <= 0:
                        unclassified_reads += 1
                        continue
                    for tid, p in cand:
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
    if (
        normalized_mode == "soft_seq"
        and species_counter
        and species_total > 0
        and sum(species_counter.values()) > 0
    ):
        uncls = species_counter.get(UNCLASSIFIED, 0)
        eff_total = max(0, species_total - uncls)
        rels = []
        for taxon, count in species_counter.items():
            if taxon == UNCLASSIFIED:
                continue
            denom = eff_total if eff_total > 0 else species_total
            rels.append((count / denom) * 100.0)
        rels.sort(reverse=True)
        top_mass = sum(rels[:LOW_DIVERSITY_TOP_K])
        shannon_index = calculate_shannon_index(species_counter)
        is_low_diversity = (shannon_index < LOW_DIVERSITY_SHANNON) or (top_mass >= LOW_DIVERSITY_TOP_MASS)

        if is_low_diversity:
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
            filtered_items = []
            for taxon, count in items:
                rel = (count / total_count) * 100 if total_count > 0 else 0.0
                if rel < MIN_REL_ABUNDANCE and count < MIN_EVIDENCE:
                    continue
                filtered_items.append((taxon, count))

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
                is_low_diversity = (shannon_index < LOW_DIVERSITY_SHANNON) or (top_mass >= LOW_DIVERSITY_TOP_MASS)

                if is_low_diversity:
                    aggressive = []
                    for taxon, count in filtered_items:
                        if taxon == UNCLASSIFIED:
                            aggressive.append((taxon, count))
                            continue
                        denom = eff_total if eff_total > 0 else total_count
                        rel = (count / denom) * 100 if denom > 0 else 0.0
                        if rel >= LOW_DIVERSITY_SPECIES_MIN_REL_ABUNDANCE:
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
