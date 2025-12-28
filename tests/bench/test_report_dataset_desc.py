from pathlib import Path
import importlib.util
import sys

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "bench" / "scripts" / "report_dataset.py"

spec = importlib.util.spec_from_file_location("report_dataset_module", REPORT)
report = importlib.util.module_from_spec(spec)
sys.modules["report_dataset_module"] = report
spec.loader.exec_module(report)


def test_abundance_report_sections_include_exact():
    latest = {
        "l1_distance_pct": 1.0,
        "bray_curtis": 0.1,
        "presence_precision": 0.5,
        "presence_recall": 0.6,
        "presence_f1": 0.55,
        "exact_l1_distance_pct": 2.0,
        "exact_bray_curtis": 0.2,
        "exact_presence_precision": 0.4,
        "exact_presence_recall": 0.5,
        "exact_presence_f1": 0.45,
    }
    best = latest

    lines = report.abundance_report_sections(latest, best)
    joined = "\n".join(lines)
    assert "Abundance (descendant-aware)" in joined
    assert "Abundance (exact species)" in joined
