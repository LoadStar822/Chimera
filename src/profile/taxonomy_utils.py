"""
轻量级 taxonomy 工具，供 profile 与 Krona 转换模块复用。
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, Optional


@dataclass
class TaxonomyMeta:
    kind: str = "auto"


@dataclass
class GtdbTaxonomy:
    parent: Dict[str, Optional[str]]
    rank: Dict[str, str]
    name: Dict[str, str]

    def iter_lineage(self, node: str) -> Iterator[str]:
        seen: set[str] = set()
        current = node
        while current and current not in seen:
            yield current
            seen.add(current)
            parent = self.parent.get(current)
            if not parent or parent == current:
                break
            current = parent

    def display_name(self, node: str) -> str:
        raw = self.name.get(node, node)
        if "__" in raw:
            prefix, suffix = raw.split("__", 1)
            if suffix:
                return suffix
        return raw


def normalize_kind(value: Optional[str]) -> str:
    if value is None:
        return "auto"
    normalized = value.strip().lower()
    if normalized not in {"auto", "ncbi", "gtdb"}:
        return "auto"
    return normalized


def read_taxonomy_meta(meta_path: Optional[str]) -> TaxonomyMeta:
    if not meta_path:
        return TaxonomyMeta()
    path = Path(meta_path)
    if not path.exists():
        return TaxonomyMeta()
    meta = TaxonomyMeta()
    try:
        with path.open("r", encoding="utf-8") as fh:
            for line in fh:
                if "=" not in line:
                    continue
                key, value = line.split("=", 1)
                key = key.strip().lower()
                value = value.strip()
                if key == "taxonomy_kind":
                    meta.kind = value or "auto"
    except OSError:
        return TaxonomyMeta()
    return meta


def load_gtdb_taxonomy(info_path: str) -> GtdbTaxonomy:
    parent: Dict[str, Optional[str]] = {}
    rank: Dict[str, str] = {}
    name: Dict[str, str] = {}
    path = Path(info_path)
    if not path.exists():
        raise FileNotFoundError(f"无法找到 GTDB taxonomy 信息文件: {info_path}")
    with path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            node, parent_id, node_rank, node_name = parts[:4]
            parent[node] = parent_id if parent_id else None
            rank[node] = node_rank.strip().lower()
            name[node] = node_name.strip()
    if not parent:
        raise ValueError(f"GTDB taxonomy 信息文件中未解析到任何节点: {info_path}")
    return GtdbTaxonomy(parent=parent, rank=rank, name=name)
