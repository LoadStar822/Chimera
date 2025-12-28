from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT / "bench" / "scripts"))

import eval_classify  # noqa: E402


def test_is_descendant_basic():
    lineage = {
        2: {1, 2},
        3: {1, 2, 3},
        4: {1, 2, 3, 4},
    }
    assert eval_classify.is_descendant(3, 2, lineage) is True
    assert eval_classify.is_descendant(4, 3, lineage) is True
    assert eval_classify.is_descendant(3, 4, lineage) is False
    assert eval_classify.is_descendant(2, 2, lineage) is True
    assert eval_classify.is_descendant(2, 3, lineage) is False


def test_descendant_counts_basic():
    truth = {"r1": "2", "r2": "2", "r3": "3", "r4": "2"}
    preds = {"r1": "3", "r2": "unclassified", "r3": "4", "r4": "5"}
    lineage = {
        3: {1, 2, 3},
        4: {1, 2, 3, 4},
        5: {1, 5},
    }

    tp, fp, fn = eval_classify.descendant_counts(truth, preds, lineage)
    assert (tp, fp, fn) == (2, 1, 2)


def test_compute_desc_metrics_basic():
    truth = {"r1": "2", "r2": "2", "r3": "3", "r4": "2"}
    preds = {"r1": "3", "r2": "unclassified", "r3": "4", "r4": "5"}
    lineage = {
        3: {1, 2, 3},
        4: {1, 2, 3, 4},
        5: {1, 5},
    }

    metrics = eval_classify.compute_desc_metrics(truth, preds, lineage)
    assert metrics["tp"] == 2
    assert metrics["fp"] == 1
    assert metrics["fn"] == 2
