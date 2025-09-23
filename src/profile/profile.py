import math
import warnings
from typing import Dict, Iterable, Optional

from ete3 import NCBITaxa


# 目标层级及其在 NCBI taxonomy 中可能出现的 rank 名称（统一转小写）
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
    """将 NCBI rank 名称归一化到我们关注的层级。"""
    rank_lower = rank.lower()
    for level, aliases in LEVEL_ALIASES.items():
        if rank_lower in aliases:
            return level
    return None


def _should_mark_unclassified(level: str, has_level: Dict[str, bool]) -> bool:
    """判断当前层级在没有命中时是否需要记为 unclassified。"""
    if level == "domain":
        return not has_level.get("acellular root", False)
    if level == "realm":
        return has_level.get("acellular root", False)
    if level == "acellular root":
        return has_level.get("realm", False) or not has_level.get("domain", False)
    return True


def calculate_shannon_index(taxon_dict):
    total_count = sum(taxon_dict.values())
    if total_count == 0:
        return 0
    shannon_index = 0
    for count in taxon_dict.values():
        p_i = count / total_count
        if p_i > 0:
            shannon_index -= p_i * math.log(p_i)
    return shannon_index


def calculate_simpson_index(taxon_dict):
    total_count = sum(taxon_dict.values())
    if total_count == 0:
        return 0
    simpson_index = 0
    for count in taxon_dict.values():
        p_i = count / total_count
        simpson_index += p_i ** 2
    return 1 - simpson_index


def process_file(input_files: Iterable[str], output_file: str) -> None:
    try:
        tax = NCBITaxa()
    except Exception as exc:  # pragma: no cover - 初始化失败极少发生
        warnings.warn(f"无法初始化 NCBITaxa，所有条目将标记为未分类: {exc}")
        tax = None

    count_by_level: Dict[str, Dict[str, int]] = {level: {} for level in LEVEL_ALIASES}

    # 处理每个输入文件
    for input_file in input_files:
        with open(input_file, 'r') as infile:
            for line in infile:
                # 按制表符分割
                parts = line.strip().split('\t')

                if len(parts) < 2:
                    continue

                # 获取第二列中的 taxid 信息
                taxid_info = parts[1].split('\t')

                for info in taxid_info:
                    if info == UNCLASSIFIED:
                        # 如果是未分类的，所有层级都加1
                        for level in count_by_level:
                            count_by_level[level][UNCLASSIFIED] = (
                                count_by_level[level].get(UNCLASSIFIED, 0) + 1
                            )
                    else:
                        taxid = info.split(':')[0]

                        try:
                            taxid_int = int(taxid)
                        except ValueError:
                            for level in count_by_level:
                                count_by_level[level][UNCLASSIFIED] = (
                                    count_by_level[level].get(UNCLASSIFIED, 0) + 1
                                )
                            continue

                        if tax is None:
                            for level in count_by_level:
                                count_by_level[level][UNCLASSIFIED] = (
                                    count_by_level[level].get(UNCLASSIFIED, 0) + 1
                                )
                            continue

                        try:
                            lineage = tax.get_lineage(taxid_int)
                            ranks = tax.get_rank(lineage)
                            names = tax.get_taxid_translator(lineage)
                        except Exception:
                            for level in count_by_level:
                                count_by_level[level][UNCLASSIFIED] = (
                                    count_by_level[level].get(UNCLASSIFIED, 0) + 1
                                )
                            continue

                        has_level = {level: False for level in count_by_level}

                        for ancestor_taxid in lineage:
                            rank = ranks.get(ancestor_taxid)
                            if not rank or rank == "no rank":
                                continue

                            normalized = _normalize_rank(rank)
                            if not normalized:
                                continue

                            name = names.get(ancestor_taxid)
                            if not name:
                                continue

                            if normalized != "clade" and has_level[normalized]:
                                continue

                            count_by_level[normalized][name] = (
                                count_by_level[normalized].get(name, 0) + 1
                            )
                            has_level[normalized] = True

                        for level, seen in has_level.items():
                            if not seen and _should_mark_unclassified(level, has_level):
                                count_by_level[level][UNCLASSIFIED] = (
                                    count_by_level[level].get(UNCLASSIFIED, 0) + 1
                                )

    # 计算多样性指标及丰度，并输出
    total_counts_by_level = {level: sum(taxon_dict.values()) for level, taxon_dict in count_by_level.items()}

    with open(output_file + '.tsv', 'w') as outfile:
        # 写入头部信息
        outfile.write("Taxon\tCount\tRelative Abundance (%)\tShannon Index\tSimpson Index\n")

        for level, taxon_dict in count_by_level.items():
            shannon_index = calculate_shannon_index(taxon_dict)
            simpson_index = calculate_simpson_index(taxon_dict)
            total_count = total_counts_by_level[level]

            if total_count == 0:
                continue

            # 输出层级标题
            display_level = " ".join(part.capitalize() for part in level.split())
            outfile.write(f"\n## {display_level} Level ##\n")

            # 按照计数从大到小排序，并将unclassified放在最后
            sorted_items = []
            unclassified_item = None
            
            # 分离unclassified和其他项
            for taxon, count in taxon_dict.items():
                if taxon == UNCLASSIFIED:
                    unclassified_item = (taxon, count)
                else:
                    sorted_items.append((taxon, count))
            
            # 对其他项按计数从大到小排序
            sorted_items.sort(key=lambda x: x[1], reverse=True)
            
            # 如果有unclassified项，添加到最后
            if unclassified_item:
                sorted_items.append(unclassified_item)
            
            # 输出排序后的结果
            for taxon, count in sorted_items:
                relative_abundance = (count / total_count) * 100 if total_count > 0 else 0
                # 输出每个分类单元的计数和相对丰度
                outfile.write(
                    f"{taxon}\t{count}\t{relative_abundance:.2f}\t{shannon_index:.4f}\t{simpson_index:.4f}\n")
