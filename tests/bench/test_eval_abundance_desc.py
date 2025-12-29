from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT / "bench" / "scripts"))

import eval_abundance  # noqa: E402


def test_collapse_to_truth_ancestors():
    parent = {"4": "3", "3": "2", "2": "1", "1": "1"}
    sci_for_taxid = {"2": "GenusX", "3": "SpeciesX", "4": "StrainX"}
    name_to_taxid = {"GenusX": "2", "SpeciesX": "3", "StrainX": "4"}
    truth_taxids = {"2"}  # genus-level truth
    pred = {"StrainX": 60.0, "Other": 40.0}

    collapsed = eval_abundance.collapse_to_truth_ancestors(
        pred,
        truth_taxids,
        name_to_taxid,
        sci_for_taxid,
        parent,
    )

    assert collapsed["GenusX"] == 60.0
    assert collapsed["Other"] == 40.0


def test_taxid_for_name_with_synonym():
    sci_for_taxid = {"2": "GenusX", "3": "SpeciesX"}
    name_to_taxid = {"GenusX": "2", "SpeciesX": "3"}
    syn_to_sci = {"SpeciesX alt": "SpeciesX"}
    sci_names = set(sci_for_taxid.values())

    taxid = eval_abundance.taxid_for_name("SpeciesX alt", name_to_taxid, syn_to_sci, sci_names)
    assert taxid == "3"


def test_drop_unclassified_also_drops_unclassified_sequences():
    pred = {"unclassified sequences": 50.0, "A": 50.0}
    out = eval_abundance.drop_unclassified(pred)
    assert "unclassified sequences" not in out
    assert out["A"] == 100.0
