import argparse
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

from multitax import GtdbTx, NcbiTx

try:
    from .profile import (
        TaxonomyResolver,
        UNCLASSIFIED,
        _parse_taxid_weights,
    )
except ImportError:  # pragma: no cover - 模块作为脚本运行
    from profile import (  # type: ignore
        TaxonomyResolver,
        UNCLASSIFIED,
        _parse_taxid_weights,
    )


class LineageProvider:
    def __init__(
        self,
        taxonomy_kind: str = "auto",
        tax_info: Optional[Path] = None,
    ) -> None:
        self.taxonomy_kind = (taxonomy_kind or "auto").lower()
        self.resolver: Optional[TaxonomyResolver] = None
        if tax_info:
            self.resolver = TaxonomyResolver(tax_info)
        if self.resolver and not self.resolver.has_taxonomy():
            self.resolver = None

        self.fallback = None
        if self.resolver is None:
            if self.taxonomy_kind in ("auto", "gtdb"):
                try:
                    self.fallback = GtdbTx()
                    self.taxonomy_kind = "gtdb"
                except Exception:
                    self.fallback = None
            if self.fallback is None:
                self.fallback = NcbiTx()
                if self.taxonomy_kind == "auto":
                    self.taxonomy_kind = "ncbi"

    def lineage_names(self, taxid: str) -> List[str]:
        if not taxid or taxid == UNCLASSIFIED:
            return []

        if self.resolver:
            lineage = self.resolver.get_lineage(taxid)
            if lineage:
                return [record.name or record.taxid for record in lineage]

        if self.fallback is None:
            return []

        node = self.fallback.latest(taxid)
        if node == self.fallback.undefined_node:
            return []

        names: List[str] = []
        for ancestor in self.fallback.lineage(node):
            name = self.fallback.name(ancestor)
            if name:
                names.append(name)
        return names


def _accumulate_counts(
    input_files: Sequence[Path],
) -> Dict[str, float]:
    counts: Dict[str, float] = defaultdict(float)
    for path in input_files:
        with path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                items, has_unclassified = _parse_taxid_weights(line)
                if not items and has_unclassified:
                    counts[UNCLASSIFIED] += 1.0
                    continue
                for taxid, weight in items:
                    counts[taxid] += weight if weight > 0 else 0.0
                if has_unclassified:
                    residue = max(
                        0.0,
                        1.0 - sum(weight for _, weight in items if weight > 0),
                    )
                    counts[UNCLASSIFIED] += residue
    return counts


def convert_multiple_files_to_krona_format(
    input_files: Iterable[str],
    output_prefix: str,
    taxonomy_kind: str = "auto",
    tax_info: Optional[str] = None,
) -> Path:
    inputs = [Path(path) for path in input_files]
    if not inputs:
        raise ValueError("至少需要一个分类结果输入文件。")
    output_path = Path(output_prefix).with_suffix(".tsv")
    provider = LineageProvider(taxonomy_kind=taxonomy_kind, tax_info=Path(tax_info) if tax_info else None)
    counts = _accumulate_counts(inputs)

    with output_path.open("w", encoding="utf-8") as outfile:
        for taxid, count in sorted(counts.items(), key=lambda kv: kv[1], reverse=True):
            if taxid == UNCLASSIFIED:
                outfile.write(f"{int(round(count))}\t{UNCLASSIFIED}\n")
                continue
            lineage = provider.lineage_names(taxid)
            if not lineage:
                lineage = [taxid]
            outfile.write(f"{int(round(count))}\t" + "\t".join(lineage) + "\n")

    return output_path


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="将 Chimera 分类结果转换为 Krona 兼容格式。",
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input_files",
        nargs="+",
        required=True,
        help="分类结果 TSV 文件（可多次指定）。",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        help="输出前缀（自动附加 .tsv）。",
    )
    parser.add_argument(
        "--taxonomy-kind",
        choices=["auto", "ncbi", "gtdb"],
        default="auto",
        help="分类体系：ncbi、gtdb 或 auto（默认自动推断）。",
    )
    parser.add_argument(
        "--tax-info",
        type=str,
        help="可选的 tax.info 路径，若提供则优先使用其中的层级信息。",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    output = convert_multiple_files_to_krona_format(
        args.input_files,
        args.output,
        taxonomy_kind=args.taxonomy_kind,
        tax_info=args.tax_info,
    )
    print(f"Krona 转换完成：{output}")


if __name__ == "__main__":
    main()
