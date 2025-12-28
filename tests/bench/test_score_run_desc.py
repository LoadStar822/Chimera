from pathlib import Path
import importlib.util
import sys


ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / "bench" / "bench.py"

spec = importlib.util.spec_from_file_location("bench_module", BENCH)
bench = importlib.util.module_from_spec(spec)
sys.modules["bench_module"] = bench
spec.loader.exec_module(bench)


def test_score_from_abundance_prefers_desc():
    data = {"dataset": "cami-long-0", "desc_presence_f1": 0.9, "presence_f1": 0.4}
    mode, score = bench.score_from_abundance(data)
    assert mode == "abundance-presence-f1"
    assert score == 0.9


def test_score_from_per_read_prefers_desc():
    data = {"desc_f1": 0.8, "f1": 0.3}
    mode, score = bench.score_from_per_read(data)
    assert mode == "per-read-f1"
    assert score == 0.8
