"""
JEOL universal IDs (jason_id) carried through simpleNMR's working tables.

Input: a real JASON export of ethylbenzene (13C 1D, HSQC, HMBC, COSY) whose
peaks carry jason_id. The expected IDs are the ones confirmed by hand in
JASON (simpleNMR_JASON_findings_ethylbenzene.ipynb, October 2026).

Run from the repository root:  python -m pytest tests -q
(needs the full server environment: RDKit, JPype + Java, ...)
"""

import contextlib
import copy
import io
import json
from pathlib import Path

import pytest

from core.html_from_assignments import NMRProblem
from app.pipelines import prediction_pipeline as pp

DATA = Path(__file__).parent / "data" / "ethylbenzene_jeol_input.json"


def _solution(json_data):
    """Run the pipeline up to carbon assignment and return the NMRsolution."""
    with contextlib.redirect_stdout(io.StringIO()):
        problem = NMRProblem.from_mnova_dict(json_data)
        solution = pp._build_solution(problem)
        pp._assign_proton_types(solution)
        pp._transfer_hsqc_info(solution)
        pp._assign_carbons(solution, pp._previous_atom_assignments_lookup(json_data))
    return solution


@pytest.fixture(scope="module")
def solution():
    return _solution(json.loads(DATA.read_text()))


def _ids(table, key_cols):
    """{tuple of rounded key values: jason_id} for a working table."""
    return {
        tuple(round(float(row[c]), 4) for c in key_cols): row["jason_id"]
        for _, row in table.iterrows()
    }


def test_carbons_carry_their_1d_13c_peak(solution):
    ids = _ids(solution.c13, ["ppm"])
    assert ids[(144.3054,)] == "{db5d5e8d-a444-4c9f-8652-052226a52d1a}"   # C2, quaternary
    assert ids[(28.9897,)] == "{60d8daa1-0725-4220-9bc5-0ed42dc2e803}"    # C7, CH2
    assert ids[(15.7155,)] == "{84fce28e-910a-4c13-bf36-b2839f6ad45a}"    # C8, CH3


def test_protons_carry_their_hsqc_peak(solution):
    ids = _ids(solution.h1, ["ppm"])
    assert ids[(2.7775,)] == "{bfb0db5f-1680-4585-aa52-955107815816}"     # H7
    assert ids[(1.3676,)] == "{49dad087-f74e-45da-ada3-1bd45abb1bbc}"     # H8
    assert set(solution.hsqc["jason_id"]) == set(solution.h1["jason_id"])


def test_hmbc_rows_carry_ids_and_snapped_positions(solution):
    # snapped (f1, f2) = the positions the peak must be moved to in JASON
    ids = _ids(solution.hmbc, ["f1_ppm", "f2_ppm"])
    assert ids[(144.3054, 7.4012)] == "{1774e688-a307-4f0c-92ad-b60b6da1c87b}"   # meta-H -> C2
    assert ids[(127.9292, 7.2956)] == "{00a416e8-712d-48b2-b645-73536f816579}"   # para-H -> ortho C
    assert ids[(144.3054, 2.7775)] == "{9cb4a1c8-fb58-4864-a418-a446a0803e6e}"   # CH2 -> C2
    assert len(ids) == 6


def test_cosy_rows_carry_ids_and_snapped_positions(solution):
    ids = _ids(solution.cosy, ["f1_ppm", "f2_ppm"])
    assert ids[(1.3676, 2.7775)] == "{c722cbc6-e190-4a44-b790-757bcb2691ab}"     # CH3 - CH2
    assert len(ids) == 3


def test_input_without_ids_still_works():
    """MNOVA / TopSpin input has no jason_id: tables simply have no ID (or None)."""
    data = json.loads(DATA.read_text())
    for block in data.values():
        if isinstance(block, dict) and block.get("datatype") == "nmrspectrum":
            for kind in ("peaks", "integrals", "multiplets"):
                for p in block.get(kind, {}).get("data", {}).values():
                    p.pop("jason_id", None)
    plain = _solution(data)
    with_ids = _solution(json.loads(DATA.read_text()))
    assert list(plain.c13["ppm"]) == list(with_ids.c13["ppm"])
    for table in (plain.c13, plain.hsqc, plain.hmbc, plain.cosy):
        assert "jason_id" not in table.columns or table["jason_id"].isna().all()


# ---------------------------------------------------------------- stage 2: the peak map

def _run_pipeline(json_data):
    with contextlib.redirect_stdout(io.StringIO()):
        return pp.run(NMRProblem.from_mnova_dict(json_data), json_data)


@pytest.fixture(scope="module")
def peak_map():
    context = _run_pipeline(json.loads(DATA.read_text()))
    return json.loads(context["peak_map"])


def test_peak_map_nodes(peak_map):
    """Each carbon group: its 1D 13C peak and HSQC peak (IDs confirmed in JASON)."""
    nodes = peak_map["nodes"]
    assert len(nodes) == 6                       # the 6 distinct carbon groups
    c2 = next(n for n in nodes.values() if round(n["c13_ppm"], 4) == 144.3054)
    assert c2["c13_1d"] == "{db5d5e8d-a444-4c9f-8652-052226a52d1a}" and c2["hsqc"] == []
    ch2 = next(n for n in nodes.values() if round(n["c13_ppm"], 4) == 28.9897)
    assert ch2["c13_1d"] == "{60d8daa1-0725-4220-9bc5-0ed42dc2e803}"
    assert [h["jason_id"] for h in ch2["hsqc"]] == ["{bfb0db5f-1680-4585-aa52-955107815816}"]
    assert round(ch2["hsqc"][0]["h_ppm"], 4) == 2.7775


# JEOL ID -> (delta2, delta1) to move the peak to: the MOVES table that gave the
# correct display in JASON (simpleNMR_JASON_findings_ethylbenzene.ipynb)
EXPECTED_MOVES = {
    "{1774e688-a307-4f0c-92ad-b60b6da1c87b}": (7.4012, 144.3054),   # HMBC meta-H -> C2
    "{bb84d936-9adf-44bd-8c46-ebd0a7aad069}": (1.3676, 144.3054),   # HMBC CH3 -> C2
    "{17aa3cdf-e8f8-478b-b07c-124e637c561c}": (1.3676, 28.9897),    # HMBC CH3 -> C7
    "{9cb4a1c8-fb58-4864-a418-a446a0803e6e}": (2.7775, 144.3054),   # HMBC CH2 -> C2
    "{82292e56-7132-4af1-a2df-269c3383b84b}": (2.7775, 127.9292),   # HMBC CH2 -> ortho
    "{00a416e8-712d-48b2-b645-73536f816579}": (7.2956, 127.9292),   # HMBC para-H -> ortho
    "{c722cbc6-e190-4a44-b790-757bcb2691ab}": (2.7775, 1.3676),     # COSY CH2 - CH3
    "{0ce9e917-3c5f-4d11-8297-93b930ab33fd}": (7.4012, 7.3226),     # COSY meta - ortho
    "{e48c637d-1420-4c76-941c-39c6c30dd5a8}": (7.4012, 7.2956),     # COSY meta - para
}


def test_peak_map_correlations_match_the_moves_that_worked(peak_map):
    got = {c["jason_id"]: (round(c["delta2"], 4), round(c["delta1"], 4))
           for c in peak_map["correlations"]}
    assert got == EXPECTED_MOVES


def test_peak_map_correlation_ends_are_nodes(peak_map):
    for c in peak_map["correlations"]:
        assert c["delta2_node"] in peak_map["nodes"] and c["delta1_node"] in peak_map["nodes"]


def test_no_peak_map_without_ids():
    """MNOVA / TopSpin / NMRium input carries no jason_id: no peak map."""
    data = json.loads(DATA.read_text())
    for block in data.values():
        if isinstance(block, dict) and block.get("datatype") == "nmrspectrum":
            for kind in ("peaks", "integrals", "multiplets"):
                for p in block.get(kind, {}).get("data", {}).values():
                    p.pop("jason_id", None)
    assert _run_pipeline(data)["peak_map"] == "null"


def test_peak_map_refuses_mismatched_nodes(solution):
    """If moved_nodes ever stopped lining up with nodes, no map rather than a wrong one."""
    graph = {"nodes": [{"id": 1, "ppm": 144.3}, {"id": 5, "ppm": 128.4}],
             "moved_nodes": [{"id": 5, "ppm": 128.4}, {"id": 1, "ppm": 144.3}]}
    assert pp._build_peak_map(solution, graph) is None
