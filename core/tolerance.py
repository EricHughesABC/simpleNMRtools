"""Resolution-based matching tolerances from the spectrum metadata.

Self-contained (standard library only) so the same file can be copied into
other projects (e.g. simpleNMRclaudePeakRefine). Keep the copies identical;
change TOLERANCE_VERSION when the rules change.

A peak position is only as accurate as the digital resolution allows, so the
tolerance for matching a value picked in one spectrum (e.g. the proton of an
HMBC cross-peak) against a reference picked in another (the HSQC proton) is a
few digital points of the coarser of the two spectra in that nucleus, with a
small floor in Hz:

    tol_Hz  = max(n_points * max(pt_spec, pt_ref), floor_Hz)
    tol_ppm = tol_Hz / spectrometer frequency (MHz)

where pt = sw_Hz / spectral_size, the point spacing after zero filling (the
same as Mnova's "Digital Resolution"). If any metadata is missing the
functions return None and the caller keeps its old behaviour.

Metadata per spectrum, as exported by the Mnova scripts (2D values are
strings such as "[f1, f2]"):
    nucleus, specfrequency, sw_Hz, spectral_size, acquired_size
"""

TOLERANCE_VERSION = "1"

N_POINTS = 2.0      # digital points of the coarser spectrum
FLOOR_HZ = 4.0      # never tighter than this


def parse_list(value):
    """'[a, b]', [a, b], 'a' or a  ->  list of floats ([] if empty or unreadable)."""
    if value is None:
        return []
    if isinstance(value, (int, float)):
        return [float(value)]
    if isinstance(value, (list, tuple)):
        items = value
    else:
        text = str(value).strip().strip("[]")
        if not text:
            return []
        items = text.split(",")
    out = []
    for item in items:
        try:
            out.append(float(str(item).strip().strip("'\"")))
        except ValueError:
            return []
    return out


def parse_nuclei(value):
    """'[13C, 1H]' or '1H' -> ['13C', '1H']."""
    if isinstance(value, (list, tuple)):
        return [str(v).strip() for v in value]
    text = str(value or "").strip().strip("[]")
    return [t.strip().strip("'\"") for t in text.split(",") if t.strip()]


def spectrum(doc, label, index=0):
    """The spectrum block for a dialog label (e.g. 'HMBC' -> doc['HMBC_0']), or None."""
    if not isinstance(doc, dict):
        return None
    block = doc.get(f"{label}_{index}")
    return block if isinstance(block, dict) else None


def point_hz(block, dim):
    """Digital point spacing in Hz of dimension dim (0 = f1, -1/last = f2), or None."""
    if not block:
        return None
    sw, si = parse_list(block.get("sw_Hz")), parse_list(block.get("spectral_size"))
    if not sw or len(sw) != len(si):
        return None
    try:
        s, n = sw[dim], si[dim]
    except IndexError:
        return None
    return s / n if n > 0 and s > 0 else None


def frequency_mhz(block, dim):
    """Spectrometer frequency in MHz of dimension dim, or None."""
    sf = parse_list(block.get("specfrequency")) if block else []
    try:
        return sf[dim] if sf and sf[dim] > 0 else None
    except IndexError:
        return None


def nucleus_dim(block, nucleus):
    """Index of the last dimension with this nucleus ('1H', '13C'), or None.

    The last one is the direct dimension (f2) when both are the same, e.g. COSY.
    """
    nuclei = parse_nuclei(block.get("nucleus")) if block else []
    hits = [i for i, n in enumerate(nuclei) if n == nucleus]
    return hits[-1] if hits else None


def match_tolerance_ppm(doc, label, dim, ref_label="HSQC", ref_nucleus="1H",
                        n_points=N_POINTS, floor_hz=FLOOR_HZ):
    """Tolerance in ppm for matching dimension dim of spectrum `label` against
    the `ref_nucleus` dimension of spectrum `ref_label`; None if metadata is missing.

    dim: 0 for f1, 1 (or -1) for f2.
    """
    spec, ref = spectrum(doc, label), spectrum(doc, ref_label)
    if spec is None or ref is None:
        return None
    ref_dim = nucleus_dim(ref, ref_nucleus)
    if ref_dim is None:
        return None
    pts = [point_hz(spec, dim), point_hz(ref, ref_dim)]
    sf = frequency_mhz(spec, dim)
    if None in pts or sf is None:
        return None
    tol_hz = max(n_points * max(pts), floor_hz)
    return tol_hz / sf


def resolvable(delta_a_ppm, delta_b_ppm, doc, label, dim, **kw):
    """True if two shifts are further apart than the matching tolerance of this
    spectrum dimension (None if metadata is missing)."""
    tol = match_tolerance_ppm(doc, label, dim, **kw)
    return None if tol is None else abs(delta_a_ppm - delta_b_ppm) > tol
