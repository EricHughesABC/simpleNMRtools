"""Tests for core/tolerance.py (resolution-based matching tolerances)."""
import pytest

from core import tolerance as T

DOC = {
    "HSQC_0": {"nucleus": "[13C, 1H]", "specfrequency": "[175.97, 699.74]",
               "sw_Hz": "[35193.1, 11363.6]", "spectral_size": "[4096, 2048]"},
    "HMBC_0": {"nucleus": "[13C, 1H]", "specfrequency": "[175.97, 699.74]",
               "sw_Hz": "[41352.3, 9090.9]", "spectral_size": "[4096, 4096]"},
    "COSY_0": {"nucleus": "[1H, 1H]", "specfrequency": "[699.74, 699.74]",
               "sw_Hz": "[3424.7, 3424.7]", "spectral_size": "[512, 4096]"},
    "H1_1D_0": {"nucleus": "1H", "specfrequency": 699.74, "sw_Hz": 13888.9, "spectral_size": 131072},
}


def test_parse_list():
    assert T.parse_list("[1.5, 2]") == [1.5, 2.0]
    assert T.parse_list(3) == [3.0]
    assert T.parse_list("") == []
    assert T.parse_list(None) == []
    assert T.parse_list("[a, 2]") == []


def test_point_and_nucleus():
    assert T.point_hz(DOC["HSQC_0"], 1) == pytest.approx(11363.6 / 2048)
    assert T.point_hz(DOC["H1_1D_0"], 0) == pytest.approx(13888.9 / 131072)
    assert T.nucleus_dim(DOC["HSQC_0"], "1H") == 1
    assert T.nucleus_dim(DOC["COSY_0"], "1H") == 1      # last = direct dimension
    assert T.nucleus_dim(DOC["HSQC_0"], "15N") is None


def test_hmbc_f2_uses_coarser_spectrum():
    # HSQC f2 point 5.55 Hz is coarser than HMBC f2 2.22 Hz -> 2 x 5.55 Hz
    tol = T.match_tolerance_ppm(DOC, "HMBC", 1)
    assert tol * 699.74 == pytest.approx(2 * 11363.6 / 2048, rel=1e-6)


def test_cosy_f1_coarse_dimension():
    tol = T.match_tolerance_ppm(DOC, "COSY", 0)
    assert tol * 699.74 == pytest.approx(2 * 3424.7 / 512, rel=1e-6)


def test_floor():
    doc = {k: dict(v) for k, v in DOC.items()}
    doc["HSQC_0"]["spectral_size"] = "[4096, 65536]"
    doc["HMBC_0"]["spectral_size"] = "[4096, 65536]"
    assert T.match_tolerance_ppm(doc, "HMBC", 1) * 699.74 == pytest.approx(T.FLOOR_HZ)


def test_missing_metadata_returns_none():
    doc = {k: dict(v) for k, v in DOC.items()}
    del doc["HMBC_0"]["sw_Hz"]
    assert T.match_tolerance_ppm(doc, "HMBC", 1) is None
    assert T.match_tolerance_ppm(DOC, "NOESY", 1) is None
    assert T.match_tolerance_ppm(None, "HMBC", 1) is None


def test_resolvable():
    assert T.resolvable(3.4280, 3.2525, DOC, "HMBC", 1) is True     # OH vs OCH3: 123 Hz
    assert T.resolvable(1.5100, 1.5121, DOC, "HMBC", 1) is False    # 1.5 Hz
