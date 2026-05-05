import csv
import gzip
from collections import defaultdict
from pathlib import Path
from statistics import median
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple


class ReadEvidenceError(ValueError):
    pass


def _open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    if path.endswith(".zst"):
        raise ReadEvidenceError(
            "zstd-compressed read evidence is not supported by this Python path; "
            "use the plain ChimeraReadEvidence.tsv output"
        )
    return open(path, "r", encoding="utf-8", newline="")


def parse_profile_evidence(value: str) -> Iterator[Tuple[str, float, float]]:
    if not value:
        return
    for entry in value.split(";"):
        if not entry:
            continue
        parts = entry.split(":")
        if len(parts) != 3:
            raise ReadEvidenceError(f"Invalid profile_evidence entry: {entry}")
        taxid, q_text, mass_text = parts
        if not taxid:
            continue
        try:
            yield taxid, float(q_text), float(mass_text)
        except ValueError as exc:
            raise ReadEvidenceError(
                f"Invalid numeric profile_evidence entry: {entry}"
            ) from exc


def load_profile_masses(evidence_path: str) -> Dict[str, float]:
    masses: Dict[str, float] = {}
    with _open_text(evidence_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "taxid" not in reader.fieldnames:
            raise ReadEvidenceError("ChimeraEvidence.tsv is missing taxid column")
        if "abundance_mass" not in reader.fieldnames:
            raise ReadEvidenceError(
                "ChimeraEvidence.tsv is missing abundance_mass column"
            )
        for row in reader:
            taxid = row.get("taxid", "")
            if not taxid:
                continue
            masses[taxid] = masses.get(taxid, 0.0) + float(
                row.get("abundance_mass") or 0.0
            )
    return masses


def reconstruct_profile_masses(read_evidence_path: str) -> Dict[str, float]:
    masses: Dict[str, float] = defaultdict(float)
    with _open_text(read_evidence_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "profile_evidence" not in reader.fieldnames:
            raise ReadEvidenceError(
                "ChimeraReadEvidence.tsv must be the wide table with "
                "profile_evidence column"
            )
        for row in reader:
            for taxid, _q, mass in parse_profile_evidence(
                row.get("profile_evidence", "")
            ):
                masses[taxid] += mass
    return dict(masses)


def load_presence_rows(presence_path: Optional[str]) -> Dict[str, Dict[str, str]]:
    if not presence_path:
        return {}
    rows: Dict[str, Dict[str, str]] = {}
    with _open_text(presence_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "taxid" not in reader.fieldnames:
            raise ReadEvidenceError("ChimeraPresence.tsv is missing taxid column")
        for row in reader:
            taxid = row.get("taxid", "")
            if taxid:
                rows[taxid] = row
    return rows


def _float_field(row: Dict[str, str], name: str, default: float = 0.0) -> float:
    try:
        return float(row.get(name, "") or default)
    except ValueError:
        return default


def _int_field(row: Dict[str, str], name: str, default: int = 0) -> int:
    try:
        return int(float(row.get(name, "") or default))
    except ValueError:
        return default


def write_profile_reconstruction_table(
    evidence_path: str,
    read_evidence_path: str,
    output_path: str,
    presence_path: Optional[str] = None,
    presence_posterior_min: float = 0.99,
    support_min: int = 5,
    support_field: str = "unique_obs",
    breadth_ratio_min: float = 0.0,
) -> Dict[str, float]:
    profile_masses = load_profile_masses(evidence_path)
    read_masses = reconstruct_profile_masses(read_evidence_path)
    presence_rows = load_presence_rows(presence_path)
    all_taxids = sorted(
        set(profile_masses) | set(read_masses),
        key=lambda x: (0, int(x)) if x.isdigit() else (1, x),
    )
    errors: List[float] = []
    extractable_mass = 0.0
    total_mass = 0.0

    parent = Path(output_path).parent
    if str(parent):
        parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "taxid",
                "profile_abundance_mass",
                "read_reconstructed_mass",
                "abs_error",
                "rel_error",
                "presence_posterior",
                "unique_reads",
                "unique_obs",
                "read_hits",
                "support_field",
                "support_value",
                "breadth_ratio",
                "extractable",
            ]
        )
        for taxid in all_taxids:
            profile_mass = profile_masses.get(taxid, 0.0)
            read_mass = read_masses.get(taxid, 0.0)
            abs_error = abs(profile_mass - read_mass)
            rel_error = abs_error / profile_mass if profile_mass > 0.0 else 0.0
            errors.append(abs_error)
            total_mass += profile_mass

            presence = presence_rows.get(taxid, {})
            presence_posterior = _float_field(presence, "presence_posterior", 0.0)
            unique_reads = _int_field(presence, "unique_reads", 0)
            unique_obs = _int_field(presence, "unique_obs", 0)
            read_hits = _int_field(presence, "read_hits", 0)
            support_value = _int_field(presence, support_field, 0)
            breadth_ratio = _float_field(presence, "breadth_ratio", 0.0)
            extractable = (
                presence_posterior >= presence_posterior_min
                and support_value >= support_min
                and breadth_ratio >= breadth_ratio_min
            )
            if extractable:
                extractable_mass += profile_mass
            writer.writerow(
                [
                    taxid,
                    f"{profile_mass:.17g}",
                    f"{read_mass:.17g}",
                    f"{abs_error:.17g}",
                    f"{rel_error:.17g}",
                    f"{presence_posterior:.17g}",
                    str(unique_reads),
                    str(unique_obs),
                    str(read_hits),
                    support_field,
                    str(support_value),
                    f"{breadth_ratio:.17g}",
                    "1" if extractable else "0",
                ]
            )

    max_abs = max(errors) if errors else 0.0
    median_abs = median(errors) if errors else 0.0
    return {
        "taxa": float(len(all_taxids)),
        "max_abs_error": max_abs,
        "median_abs_error": median_abs,
        "total_profile_mass": total_mass,
        "extractable_mass": extractable_mass,
        "extractable_fraction": extractable_mass / total_mass
        if total_mass > 0.0
        else 0.0,
    }


def load_target_taxids(taxids: Sequence[str], taxids_file: Optional[str]) -> Set[str]:
    out: Set[str] = set()
    for taxid in taxids:
        if taxid:
            out.add(str(taxid))
    if taxids_file:
        with open(taxids_file, "r", encoding="utf-8") as handle:
            for line in handle:
                text = line.strip()
                if not text or text.startswith("#"):
                    continue
                out.add(text.split()[0])
    if not out:
        raise ReadEvidenceError("At least one --taxid or --taxids file is required")
    return out


def read_id_keys(read_id: str) -> Set[str]:
    text = read_id.strip()
    keys = {text}
    if text:
        keys.add(text.split()[0])
    expanded = set(keys)
    for key in keys:
        if key.endswith("/1") or key.endswith("/2"):
            expanded.add(key[:-2])
    return {key for key in expanded if key}


def load_read_selection(
    ledger_path: str, target_taxids: Set[str], min_profile_mass: float
) -> Tuple[Dict[str, Set[str]], Dict[str, int]]:
    key_to_taxids: Dict[str, Set[str]] = defaultdict(set)
    counts: Dict[str, int] = defaultdict(int)
    with _open_text(ledger_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "read_id" not in reader.fieldnames:
            raise ReadEvidenceError("ChimeraReadEvidence.tsv is missing read_id column")
        if "profile_evidence" not in reader.fieldnames:
            raise ReadEvidenceError(
                "ChimeraReadEvidence.tsv must be the wide table with "
                "profile_evidence column"
            )
        for row in reader:
            matched_taxids: Set[str] = set()
            for taxid, _q, mass in parse_profile_evidence(
                row.get("profile_evidence", "")
            ):
                if taxid in target_taxids and mass >= min_profile_mass:
                    matched_taxids.add(taxid)
            if not matched_taxids:
                continue
            read_id = row.get("read_id", "")
            for key in read_id_keys(read_id):
                key_to_taxids[key].update(matched_taxids)
            for taxid in matched_taxids:
                counts[taxid] += 1
    return dict(key_to_taxids), dict(counts)


def _iter_fastq(path: str) -> Iterator[Tuple[str, str, str, str]]:
    with _open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not qual:
                raise ReadEvidenceError(f"Truncated FASTQ record in {path}")
            if not header.startswith("@"):
                raise ReadEvidenceError(f"Invalid FASTQ header in {path}: {header[:80]}")
            yield header, seq, plus, qual


def _header_keys(header: str) -> Set[str]:
    return read_id_keys(header[1:].rstrip("\n\r"))


def _safe_taxid_name(taxid: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in taxid)


def extract_read_bags(
    ledger_path: str,
    read_paths: Sequence[str],
    target_taxids: Set[str],
    output_dir: str,
    min_profile_mass: float = 0.0,
) -> Dict[str, int]:
    if len(read_paths) not in (1, 2):
        raise ReadEvidenceError("--reads expects one FASTQ or a paired FASTQ pair")
    key_to_taxids, _ledger_counts = load_read_selection(
        ledger_path, target_taxids, min_profile_mass
    )
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    handles: Dict[Tuple[str, int], object] = {}
    written: Dict[str, int] = defaultdict(int)

    def get_handle(taxid: str, mate_idx: int):
        key = (taxid, mate_idx)
        if key in handles:
            return handles[key]
        safe = _safe_taxid_name(taxid)
        if len(read_paths) == 1:
            path = out_dir / f"taxid_{safe}.fastq"
        else:
            path = out_dir / f"taxid_{safe}_R{mate_idx + 1}.fastq"
        handle = open(path, "w", encoding="utf-8", newline="")
        handles[key] = handle
        return handle

    def taxids_for_headers(headers: Iterable[str]) -> Set[str]:
        found: Set[str] = set()
        for header in headers:
            for key in _header_keys(header):
                found.update(key_to_taxids.get(key, set()))
        return found

    try:
        if len(read_paths) == 1:
            for record in _iter_fastq(read_paths[0]):
                taxids = taxids_for_headers([record[0]])
                for taxid in taxids:
                    handle = get_handle(taxid, 0)
                    handle.writelines(record)
                    written[taxid] += 1
        else:
            iter1 = _iter_fastq(read_paths[0])
            iter2 = _iter_fastq(read_paths[1])
            while True:
                try:
                    rec1 = next(iter1)
                except StopIteration:
                    try:
                        next(iter2)
                    except StopIteration:
                        break
                    raise ReadEvidenceError("Paired FASTQ files have different lengths")
                try:
                    rec2 = next(iter2)
                except StopIteration as exc:
                    raise ReadEvidenceError(
                        "Paired FASTQ files have different lengths"
                    ) from exc
                taxids = taxids_for_headers([rec1[0], rec2[0]])
                for taxid in taxids:
                    get_handle(taxid, 0).writelines(rec1)
                    get_handle(taxid, 1).writelines(rec2)
                    written[taxid] += 1
    finally:
        for handle in handles.values():
            handle.close()

    manifest_path = out_dir / "extract_reads_manifest.tsv"
    with open(manifest_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["taxid", "read_count", "output"])
        for taxid in sorted(
            target_taxids, key=lambda x: (0, int(x)) if x.isdigit() else (1, x)
        ):
            safe = _safe_taxid_name(taxid)
            if len(read_paths) == 1:
                output = f"taxid_{safe}.fastq"
            else:
                output = f"taxid_{safe}_R1.fastq;taxid_{safe}_R2.fastq"
            writer.writerow([taxid, str(written.get(taxid, 0)), output])
    return dict(written)
