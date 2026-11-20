"""
benchmark_dataset.py
====================
Run Simulated Annealing and CP-SAT on a dataset of simpleNMR JSON export
files and produce a comparison table suitable for a report.

Usage
-----
    # Run on all JSON files in a folder:
    python benchmark_dataset.py /path/to/json/folder

    # Run on specific files:
    python benchmark_dataset.py file1.json file2.json file3.json

    # Save results table to CSV:
    python benchmark_dataset.py /path/to/folder --csv results.csv

    # Increase CP-SAT time limit for larger molecules:
    python benchmark_dataset.py /path/to/folder --time-limit 120

    # Show verbose CP-SAT search log:
    python benchmark_dataset.py /path/to/folder --verbose

Must be run from the simpleNMRtools repo root so that config and core
packages are importable:

    cd /path/to/simpleNMRtools
    python /path/to/benchmark_dataset.py /path/to/json/folder
"""

import copy
import json
import math
import sys
import time
import argparse
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ── CP-SAT solver ─────────────────────────────────────────────────────────────
from ortools.sat.python import cp_model

# ── simpleNMR imports — add repo root so config/core are importable ───────────
# Works whether the script is run from the repo root or from elsewhere.
_repo_candidates = [
    Path(__file__).parent,                    # script is inside the repo
    Path.cwd(),                               # run from repo root
    Path(__file__).parent.parent,             # script one level above repo
]
for _candidate in _repo_candidates:
    if (_candidate / "config" / "globals.py").exists():
        sys.path.insert(0, str(_candidate))
        break
from config.globals import (
    SIMULATED_ANNEALING_RUNS,
    SIMULATED_ANNEALING_COOLING_RATE,
)
from core.simulated_annealing_claude import (
    SimulatedAnnealing2,
    compute_total_weight,
)


# ─────────────────────────────────────────────────────────────────────────────
# JSON loading with auto-fixes for older export formats
# ─────────────────────────────────────────────────────────────────────────────

def load_json(path: Path) -> Tuple[dict, List[str]]:
    """
    Load a simpleNMR export JSON and apply compatibility fixes.

    Returns (data_dict, list_of_fixes_applied).
    An empty fix list means the file was fully clean.
    """
    d = json.load(open(path))
    fixes = []

    # Fix 1: stringify atomNumber (mixed int/str in older exports)
    n_fixed = 0
    for n in d.get("nodes_now", []):
        if "atomNumber" in n and not isinstance(n["atomNumber"], str):
            n["atomNumber"] = str(n["atomNumber"])
            n_fixed += 1
    if n_fixed:
        fixes.append(f"atomNumber stringified ({n_fixed} nodes)")

    # Fix 2: orig_nodes_displayed missing or empty
    if not d.get("orig_nodes_displayed"):
        d["orig_nodes_displayed"] = [n["id"] for n in d["nodes_now"]]
        fixes.append("orig_nodes_displayed set to all node IDs")

    # Fix 3: possible_symmetry missing
    if "possible_symmetry" not in d:
        d["possible_symmetry"] = False
        fixes.append("possible_symmetry defaulted to False")

    # Fix 4: pool_ch3_ch1 missing
    if "pool_ch3_ch1" not in d:
        d["pool_ch3_ch1"] = False
        fixes.append("pool_ch3_ch1 defaulted to False")

    return d, fixes


# ─────────────────────────────────────────────────────────────────────────────
# Shared solver construction
# ─────────────────────────────────────────────────────────────────────────────

def make_solver(d: dict, randomize: bool = True) -> SimulatedAnnealing2:
    molfile = d["oldjsondata"]["molfile"]["data"]["0"]
    s = SimulatedAnnealing2.from_params(
        copy.deepcopy(d["nodes_now"]),
        copy.deepcopy(d["links"]),
        molfile,
        d["oldjsondata"],
        d["possible_symmetry"],
        d["orig_nodes_displayed"],
        d["pool_ch3_ch1"],
    )
    s.setup_run(
        cooling_rate=SIMULATED_ANNEALING_COOLING_RATE,
        randomize_mapping=randomize,
    )
    return s


# ─────────────────────────────────────────────────────────────────────────────
# CP-SAT model builder
# ─────────────────────────────────────────────────────────────────────────────

def _build_cpsat_model(
    carbon_graph,
    shortest_paths,
    assignable_groups,
    nP_reduced,
    phase2: bool = False,
):
    model = cp_model.CpModel()
    assign = {}
    for nk, ng in assignable_groups.items():
        for n in ng:
            for m in ng:
                assign[(n, m)] = model.new_bool_var(f"x_{n}_{m}")
    for nk, ng in assignable_groups.items():
        for n in ng:
            model.add_exactly_one(assign[(n, m)] for m in ng)
        for m in ng:
            model.add_exactly_one(assign[(n, m)] for n in ng)

    n2np = {int(n): int(p) for n, p in nP_reduced.items()}

    def edge_penalty(ic, ih, path):
        """Identical to compute_total_weight in the production SA code."""
        pen = 0
        if ic:
            pen += (path - 1) ** 3
        if ih:
            eff = path if path >= 3 else 2
            pen += (eff - 2) ** 3
        return pen

    pv, ot = {}, []
    const = 0   # penalty from correlations between two FIXED nodes
    for u, v, d2 in carbon_graph.edges(data=True):
        ic, ih = bool(d2.get("cosy")), bool(d2.get("hmbc"))
        if not (ic or ih):
            continue
        u_free, v_free = u in n2np, v in n2np

        if u_free and v_free:
            # both ends can move: penalty depends on the pair of assignments
            ug, vg = assignable_groups[n2np[u]], assignable_groups[n2np[v]]
            for mu in ug:
                for mv in vg:
                    pen = edge_penalty(ic, ih, shortest_paths[mu][mv])
                    if pen == 0:
                        continue
                    key = (u, v, mu, mv)
                    if key not in pv:
                        p = model.new_bool_var(f"p_{u}_{v}_{mu}_{mv}")
                        pv[key] = p
                        au, av = assign[(u, mu)], assign[(v, mv)]
                        model.add_bool_and([au, av]).only_enforce_if(p)
                        model.add_bool_or([au.Not(), av.Not()]).only_enforce_if(p.Not())
                    ot.append(pen * pv[key])

        elif u_free or v_free:
            # one end is fixed (maps to itself): penalty depends on ONE
            # assignment, so no product variable is needed
            free, fixed = (u, v) if u_free else (v, u)
            for mf in assignable_groups[n2np[free]]:
                pen = edge_penalty(ic, ih, shortest_paths[mf][fixed])
                if pen:
                    ot.append(pen * assign[(free, mf)])

        else:
            # both ends fixed: a constant, but it belongs in the reported score
            const += edge_penalty(ic, ih, shortest_paths[u][v])

    oe = sum(ot) + const
    if not phase2:
        model.minimize(oe)
    return model, assign, oe


class _AllSolutionsCb(cp_model.CpSolverSolutionCallback):
    def __init__(self, av, ag, fi, mx):
        super().__init__()
        self._av, self._ag, self._fi, self._mx = av, ag, fi, mx
        self.solutions = []
        self._seen = set()
        self.duplicates = 0     # times CP-SAT reported an assignment already seen

    def on_solution_callback(self):
        m = self._fi.copy()
        for nk, ng in self._ag.items():
            for n in ng:
                for mv in ng:
                    if self.Value(self._av[(n, mv)]):
                        m[n] = mv
                        break
        key = tuple(sorted(m.items()))
        if key in self._seen:           # CP-SAT can report a solution twice
            self.duplicates += 1
            return
        if len(self.solutions) >= self._mx:
            self.StopSearch()
            return
        self._seen.add(key)
        self.solutions.append(m)


# ─────────────────────────────────────────────────────────────────────────────
# Run SA
# ─────────────────────────────────────────────────────────────────────────────

def run_sa(d: dict, catoms_df, jgd_tmpl: dict, n_runs: int) -> dict:
    sa = make_solver(d)
    t0 = time.perf_counter()
    sa.run_optimization(num_times=n_runs)
    t_sa = time.perf_counter() - t0
    _, best = sa.process_results(catoms_df.copy(), copy.deepcopy(jgd_tmpl))
    w = sa.bestest_weight
    unique = len({
        "".join(str(v) for v in sorted(m.items()))
        for m in sa.results[w]["best_mapping"]
    })
    return {
        "time": t_sa,
        "weight": w,
        "unique_solutions": unique,
        "mae": best["best_mae"],
        "lae": best["best_lae"],
        "lae_atom": best["best_lae_atomNumber"],
        "proven_optimal": False,
        "success_rate": sa.results[w]["num_times"] / n_runs,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Run CP-SAT
# ─────────────────────────────────────────────────────────────────────────────

def run_cpsat(
    d: dict,
    catoms_df,
    jgd_tmpl: dict,
    max_solutions: int,
    time_limit: float,
    enum_time_limit: float,
    verbose: bool,
) -> dict:
    solver_obj = make_solver(d, randomize=False)
    displayed = set(solver_obj.orig_nodes_displayed)
    ag = {
        k: [n for n in v if n in displayed]
        for k, v in solver_obj.carbon_grouped_nodes.items()
    }
    ag = {k: v for k, v in ag.items() if len(v) > 1}
    nP = {
        k: v for k, v in solver_obj.nProtons_to_nodes.items()
        if k in displayed and int(v) in ag
    }
    fi = {node: node for node in solver_obj.carbon_graph.nodes()}

    # Phase 1 — prove W*
    t0 = time.perf_counter()
    m1, av1, oe1 = _build_cpsat_model(
        solver_obj.carbon_graph, solver_obj.shortest_paths, ag, nP
    )
    s1 = cp_model.CpSolver()
    s1.parameters.max_time_in_seconds = time_limit
    s1.parameters.num_workers = 8
    s1.parameters.log_search_progress = verbose
    st1 = s1.solve(m1)
    t_phase1 = time.perf_counter() - t0

    is_opt = (st1 == cp_model.OPTIMAL)
    is_feasible = st1 in (cp_model.OPTIMAL, cp_model.FEASIBLE)

    if not is_feasible:
        return {
            "time": t_phase1, "weight": float("inf"),
            "unique_solutions": 0, "mae": float("inf"),
            "lae": float("inf"), "lae_atom": -1,
            "proven_optimal": False, "phase1_time": t_phase1, "phase2_time": 0,
            "gap": float("inf"),
        }

    W_star = int(round(s1.objective_value))
    bound  = int(round(s1.best_objective_bound))
    gap    = abs(W_star - bound) / max(W_star, 1) if W_star > 0 else 0.0

    # Phase 2 — enumerate all at W* (pure feasibility, no minimize)
    t1 = time.perf_counter()
    m2, av2, oe2 = _build_cpsat_model(
        solver_obj.carbon_graph, solver_obj.shortest_paths, ag, nP,
        phase2=True,
    )
    m2.add(oe2 == W_star)
    cb = _AllSolutionsCb(av2, ag, fi, max_solutions)
    s2 = cp_model.CpSolver()
    s2.parameters.max_time_in_seconds = enum_time_limit
    s2.parameters.num_workers = 1      # enumeration is cheap; 1 worker is deterministic
    s2.parameters.enumerate_all_solutions = True
    s2.parameters.log_search_progress = False
    st2 = s2.solve(m2, cb)
    t_phase2 = time.perf_counter() - t1

    solutions = cb.solutions
    # complete = CP-SAT exhausted the search (not cut off by the cap or time limit)
    enum_complete = (st2 == cp_model.OPTIMAL) and len(solutions) < max_solutions
    if not solutions:
        enum_complete = False
        # fallback to Phase 1 solution
        m_fb = fi.copy()
        for nk, ng in ag.items():
            for n in ng:
                for mv in ng:
                    if s1.value(av1[(n, mv)]):
                        m_fb[n] = mv
                        break
        solutions = [m_fb]

    # Score solutions and process_results via solver object
    solver_obj.results = {
        W_star: {"num_times": len(solutions), "results": [], "best_mapping": solutions}
    }
    for sol in solutions:
        solver_obj.results[W_star]["results"].append((sol, {}, []))
    solver_obj.bestest_weight  = W_star
    solver_obj.bestest_mapping = solutions[0]
    solver_obj.bestest_stats   = {}
    solver_obj.bestest_weights = [(0, W_star, W_star, 0.0)]

    _, best = solver_obj.process_results(catoms_df.copy(), copy.deepcopy(jgd_tmpl))

    # distinct mappings (compare the full assignment, not just its set of values)
    unique = len({tuple(sorted(m.items())) for m in solutions})
    # safety check: every CP-SAT solution must score W* under the production scorer
    mismatch = sum(
        1 for m in solutions
        if compute_total_weight(solver_obj.carbon_graph, m,
                                solver_obj.shortest_paths) != W_star
    )
    total_time = t_phase1 + t_phase2

    return {
        "time": total_time,
        "weight": W_star,
        "unique_solutions": unique,
        "mae": best["best_mae"],
        "lae": best["best_lae"],
        "lae_atom": best["best_lae_atomNumber"],
        "proven_optimal": is_opt,
        "phase1_time": t_phase1,
        "phase2_time": t_phase2,
        "gap": gap,
        "enum_complete": enum_complete,
        "solutions": solutions,
        "duplicates": cb.duplicates,
        "weight_mismatch": mismatch,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Process one file
# ─────────────────────────────────────────────────────────────────────────────

def shift_only_start(d: dict, carbon_df) -> dict:
    """Rebuild the closest-carbon-shift assignment from an export.

    The export may already be optimised (the pipeline copies the optimised
    positions back into the atom table), so the pre-optimisation state is
    reconstructed the way the pipeline builds it: a Hungarian match of the
    experimental peaks to the predicted shifts, minimising total |delta shift|
    within proton-count groups (CH3 and CH1 pooled when pool_ch3_ch1 is set).
    Validated on exports saved with SA switched off, where it reproduces the
    stored assignment exactly.

    Returns {atom id the peak sits at in the export: atom id it is matched to}.
    """
    import numpy as np
    from scipy.optimize import linear_sum_assignment

    peaks = {n["id"]: (float(n["ppm"]), int(n["numProtons"])) for n in d["nodes_now"]}
    calc = {i: (float(carbon_df.loc[i, "ppm_calculated"]), int(carbon_df.loc[i, "numProtons"]))
            for i in peaks}
    pooled = bool(d.get("pool_ch3_ch1"))
    key = (lambda nH: "CH3+CH1" if nH in (1, 3) else nH) if pooled else (lambda nH: nH)
    S = {}
    for g in {key(p[1]) for p in peaks.values()}:
        pk = sorted(i for i, p in peaks.items() if key(p[1]) == g)
        at = sorted(i for i, c in calc.items() if key(c[1]) == g)
        if len(pk) != len(at):          # not expected: every atom in the set has a peak
            for i in pk:
                S[i] = i
            continue
        cost = np.abs(np.array([peaks[i][0] for i in pk])[:, None]
                      - np.array([calc[j][0] for j in at])[None, :])
        r, c = linear_sum_assignment(cost)
        for a, b in zip(r, c):
            S[pk[a]] = at[b]
    return S


def assess_move(solver, S, solutions, start_w, w_star, cp_mae, max_dmae, max_dlae) -> dict:
    """Judge the move from the shift-only start to the correlation optimum.

    Uses only the shifts and correlations in the export. 'optimum' is the lowest
    (MAE, LAE) assignment among all CP-SAT optimal assignments.  A carbon's move is
    'forced' when EVERY optimal assignment puts it somewhere other than its
    shift-only atom (so the correlations alone require the move); otherwise it is a
    tie-break choice.  Verdict rule: accept when the score falls, the shifts get
    only slightly worse, and every moved carbon has at least one correlation.
    """
    import numpy as np
    def _f(v):                      # non-numeric shifts (empty strings) are skipped,
        try:                        # as in the production MAE/LAE calculation
            return float(v)
        except (TypeError, ValueError):
            return float("nan")

    exp = {n["id"]: _f(n["ppm"]) for n in solver.nmr_nodes}
    calc = {i: _f(v) for i, v in solver.carbon_df["ppm_calculated"].items()}

    def err(m):
        e = [abs(exp[n] - calc[m[n]]) for n in exp
             if not (np.isnan(exp[n]) or np.isnan(calc.get(m[n], float("nan"))))]
        return (float(np.mean(e)), float(np.max(e))) if e else (float("nan"), float("nan"))

    stats = [err(m) for m in solutions]
    bi = min(range(len(solutions)), key=lambda i: stats[i])
    best = solutions[bi]
    mae1, lae1 = stats[bi]
    mae0, lae0 = err(S)

    g = solver.carbon_graph
    deg = {n: sum(1 for _, _, dd in g.edges(n, data=True) if dd.get("cosy") or dd.get("hmbc"))
           for n in exp}
    moved = [n for n in exp if best[n] != S[n]]
    forced = [n for n in moved if all(m[n] != S[n] for m in solutions)]
    ambiguous = [n for n in exp if len({m[n] for m in solutions}) > 1]
    moved_no_corr = [n for n in moved if deg[n] == 0]
    d_mae, d_lae = mae1 - mae0, lae1 - lae0

    if start_w <= w_star:
        verdict = "tie-break only" if moved else "no change needed"
    else:
        flags = []
        if d_mae > max_dmae:
            flags.append(f"MAE +{d_mae:.2f} ppm")
        if d_lae > max_dlae:
            flags.append(f"LAE +{d_lae:.1f} ppm")
        if moved_no_corr:
            flags.append(f"{len(moved_no_corr)} moved carbon(s) with no correlation")
        verdict = "accepted" if not flags else "FLAG: " + "; ".join(flags)

    return {
        "mae0": mae0, "mae1": mae1, "lae0": lae0, "lae1": lae1,
        "d_mae": d_mae, "d_lae": d_lae,
        "n_moved": len(moved), "n_forced": len(forced),
        "n_moved_no_corr": len(moved_no_corr), "n_ambiguous": len(ambiguous),
        "cov1": sum(1 for n in exp if deg[n] >= 1) / len(exp),
        "cov2": sum(1 for n in exp if deg[n] >= 2) / len(exp),
        "verdict": verdict,
        "selection_rule_differs": abs(cp_mae - mae1) > 1e-4,
    }


def get_smiles(d: dict) -> str:
    """SMILES from a simpleNMR export: top-level key first, then oldjsondata."""
    s = d.get("smilesString")
    if not s:
        try:
            s = d["oldjsondata"]["smiles"]["data"]["0"]
        except (KeyError, TypeError):
            s = ""
    return str(s).strip() if s else "n/a"


def process_file(
    path: Path,
    n_runs: int,
    cp_time_limit: float,
    cp_enum_limit: float,
    verbose: bool,
    max_solutions: int = 2000,
    max_dmae: float = 1.0,
    max_dlae: float = 10.0,
) -> dict:
    d, fixes = load_json(path)
    nodes = d["nodes_now"]
    links = d["links"]
    pg = Counter(n["numProtons"] for n in nodes)
    n_cosy = sum(1 for l in links if l.get("cosy"))
    n_hmbc = sum(1 for l in links if l.get("hmbc"))

    base = make_solver(d, randomize=False)
    catoms_df = base.carbon_df.copy()
    jgd_tmpl  = {"nodes": copy.deepcopy(nodes), "links": copy.deepcopy(links)}
    predicted_w = base.predicted_weight

    S = shift_only_start(d, catoms_df)
    start_w = compute_total_weight(base.carbon_graph, S, base.shortest_paths)
    start_n_differ = sum(1 for u, a in S.items() if a != u)
    source = d["oldjsondata"].get("MNOVAcalcMethod", {}).get("data", {}).get("0", "n/a")

    sa_result = run_sa(d, catoms_df, jgd_tmpl, n_runs)
    cp_result  = run_cpsat(
        d, catoms_df, jgd_tmpl, max_solutions,
        cp_time_limit, cp_enum_limit, verbose,
    )

    solutions = cp_result.pop("solutions", [])
    if solutions:
        assess = assess_move(base, S, solutions, start_w, cp_result["weight"],
                             cp_result["mae"], max_dmae, max_dlae)
    else:
        nan = float("nan")
        assess = {k: nan for k in ("mae0", "mae1", "lae0", "lae1", "d_mae", "d_lae",
                                   "n_moved", "n_forced", "n_moved_no_corr", "n_ambiguous",
                                   "cov1", "cov2")}
        assess.update(verdict="no CP-SAT solution", selection_rule_differs=False)

    return {
        "assess":      assess,
        "name":        path.stem,
        "smiles":      get_smiles(d),
        "source":      source,
        "start_w":     start_w,
        "start_n_differ": start_n_differ,
        "n_carbons":   len(nodes),
        "ch3":         pg[3], "ch2": pg[2], "ch1": pg[1], "c": pg[0],
        "n_cosy":      n_cosy,
        "n_hmbc":      n_hmbc,
        "predicted_w": predicted_w,
        "symmetry":    d.get("possible_symmetry", False),
        "fixes":       fixes,
        "sa":          sa_result,
        "cp":          cp_result,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Table printing
# ─────────────────────────────────────────────────────────────────────────────

def print_results_table(results: List[dict]) -> None:
    """Print a formatted comparison table to stdout."""

    # ── Per-molecule detail ───────────────────────────────────────────────────
    print(f"\n{'═'*90}")
    print("  DETAILED RESULTS")
    print(f"{'═'*90}")

    for r in results:
        sa = r["sa"]; cp = r["cp"]
        diff = sa["weight"] - cp["weight"]
        agree = abs(sa["mae"] - cp["mae"]) < 0.001
        speedup = sa["time"] / cp["time"] if cp["time"] > 0 else float("inf")
        faster = f"CP-SAT {speedup:.1f}× faster" if speedup >= 1 else \
                 f"SA {1/speedup:.1f}× faster"

        print(f"\n  {r['name']}")
        print(f"  {'─'*70}")
        print(f"  SMILES   : {r['smiles']}")
        print(f"  Carbons  : {r['n_carbons']}  "
              f"(CH3={r['ch3']} CH2={r['ch2']} CH={r['ch1']} C={r['c']})  "
              f"COSY={r['n_cosy']}  HMBC={r['n_hmbc']}  "
              f"symmetry={'yes' if r['symmetry'] else 'no'}")
        if r["fixes"]:
            print(f"  Auto-fix : {'; '.join(r['fixes'])}")
        print(f"  Predicted weight: {r['predicted_w']:.0f}   (score of the assignment stored in the export)")
        need = r["start_w"] > cp["weight"]
        print(f"  Shift-only start ({r['source']}): score {r['start_w']:.0f} -> optimum {cp['weight']:.0f}   "
              f"{'optimisation REQUIRED' if need else 'already optimal'};  "
              f"{r['start_n_differ']} of {r['n_carbons']} peaks differ from the export")
        print(f"  SA restarts that reached the optimum: {sa['success_rate']*100:.0f}%")
        a = r["assess"]
        print(f"  Shift-only -> optimum : MAE {a['mae0']:.2f} -> {a['mae1']:.2f} ppm ({a['d_mae']:+.2f}),  "
              f"LAE {a['lae0']:.1f} -> {a['lae1']:.1f} ppm ({a['d_lae']:+.1f})")
        print(f"  Carbons moved         : {a['n_moved']:.0f} of {r['n_carbons']}  "
              f"(forced by correlations {a['n_forced']:.0f}, tie-break only "
              f"{a['n_moved'] - a['n_forced']:.0f});  moved with no correlation: {a['n_moved_no_corr']:.0f}")
        print(f"  Correlation coverage  : {100*a['cov1']:.0f}% of carbons have >=1, "
              f"{100*a['cov2']:.0f}% have >=2;  carbons that differ among tied optima: {a['n_ambiguous']:.0f}")
        print(f"  Verdict               : {a['verdict']}")
        print()
        print(f"  {'Metric':<36} {'SA':>12} {'CP-SAT':>12}")
        print(f"  {'─'*60}")
        print(f"  {'Wall time (s)':<36} {sa['time']:>12.2f} {cp['time']:>12.2f}")
        print(f"  {'  Phase 1 / Phase 2 (s)':<36} {'N/A':>12} "
              f"{cp.get('phase1_time',0):>6.2f} / {cp.get('phase2_time',0):>5.2f}")
        print(f"  {'Speed comparison':<36} {'':>12} {faster:>12}")
        print(f"  {'Best correlation score (W)':<36} {sa['weight']:>12.0f} "
              f"{cp['weight']:>12.0f}")
        if diff != 0:
            print(f"  {'  ⚠  SA missed global optimum by':<36} {diff:>+12.0f} {'':>12}")
        print(f"  {'Unique optimal solutions':<36} {sa['unique_solutions']:>12} "
              f"{cp['unique_solutions']:>12}")
        enum_str = "Yes ✓" if cp.get("enum_complete") else "No - truncated"
        print(f"  {'CP-SAT optimal set complete?':<36} {'n/a':>12} {enum_str:>12}")
        if cp.get("weight_mismatch"):
            print(f"  ⚠  {cp['weight_mismatch']} CP-SAT solution(s) do not re-score to W* "
                  f"under the production scorer - check the model")
        print(f"  {'MAE — best solution (ppm)':<36} {sa['mae']:>12.4f} "
              f"{cp['mae']:>12.4f}  {'✓' if agree else '✗'}")
        print(f"  {'LAE — best solution (ppm)':<36} {sa['lae']:>12.4f} "
              f"{cp['lae']:>12.4f}  (atom #{cp['lae_atom']})")
        opt_str = "Yes ✓" if cp["proven_optimal"] else f"No (gap={cp.get('gap',0)*100:.1f}%)"
        print(f"  {'Proven global optimum?':<36} {'No':>12} {opt_str:>12}")

    # ── Summary table ─────────────────────────────────────────────────────────
    print(f"\n\n{'═'*110}")
    print("  SUMMARY TABLE")
    print(f"{'═'*110}")
    hdr = (f"  {'Molecule':<28} {'C':>3} {'COSY':>5} {'HMBC':>5} "
           f"{'Pred W':>7} {'SA W':>7} {'CP W*':>7} {'Δ':>5} "
           f"{'SA time':>9} {'CP time':>9} {'Speedup':>9} "
           f"{'SA MAE':>8} {'CP MAE':>8} {'Opt?':>6}  {'SMILES'}")
    print(hdr)
    print(f"  {'─'*106}")

    for r in results:
        sa = r["sa"]; cp = r["cp"]
        diff = sa["weight"] - cp["weight"]
        speedup = sa["time"] / cp["time"] if cp["time"] > 0 else float("inf")
        opt = "Yes ✓" if cp["proven_optimal"] else f"No({cp['gap']*100:.0f}%)"
        if not cp.get("enum_complete", True):
            opt += "*"
        faster_str = f"{speedup:.1f}×" if speedup >= 1 else f"0.{int(1/speedup*10)}×SA"

        print(
            f"  {r['name'][:28]:<28} {r['n_carbons']:>3} "
            f"{r['n_cosy']:>5} {r['n_hmbc']:>5} "
            f"{r['predicted_w']:>7.0f} {sa['weight']:>7.0f} {cp['weight']:>7.0f} "
            f"{diff:>+5.0f} "
            f"{sa['time']:>8.2f}s {cp['time']:>8.2f}s {faster_str:>9} "
            f"{sa['mae']:>8.4f} {cp['mae']:>8.4f} {opt:>6}  {r['smiles']}"
        )

    print(f"  {'─'*106}")

    # Aggregate stats
    n = len(results)
    sa_faster = sum(1 for r in results if r["sa"]["time"] < r["cp"]["time"])
    cp_faster = n - sa_faster
    sa_missed = sum(1 for r in results if r["sa"]["weight"] > r["cp"]["weight"])
    all_agree = sum(
        1 for r in results
        if abs(r["sa"]["mae"] - r["cp"]["mae"]) < 0.001
    )
    avg_speedup = sum(
        r["sa"]["time"] / r["cp"]["time"]
        for r in results if r["cp"]["time"] > 0
    ) / n

    print(f"\n  Molecules tested       : {n}")
    print(f"  CP-SAT faster          : {cp_faster}/{n}")
    print(f"  SA faster              : {sa_faster}/{n}")
    print(f"  Avg speedup (SA/CP)    : {avg_speedup:.1f}×")
    print(f"  SA found wrong answer  : {sa_missed}/{n}")
    print(f"  MAE agreement          : {all_agree}/{n}")
    print(f"  CP-SAT proven optimal  : "
          f"{sum(1 for r in results if r['cp']['proven_optimal'])}/{n}")
    n_complete = sum(1 for r in results if r["cp"].get("enum_complete"))
    print(f"  CP-SAT set complete    : {n_complete}/{n}")
    n_bad = sum(1 for r in results if r["cp"].get("weight_mismatch"))
    print(f"  Re-score mismatches    : {n_bad}/{n}  (should be 0)")
    need_n = sum(1 for r in results if r["start_w"] > r["cp"]["weight"])
    print(f"\n  Shift-only start already optimal   : {n - need_n}/{n}")
    print(f"  Correlation-based optimisation needed: {need_n}/{n}")
    by_src = {}
    for r in results:
        t = by_src.setdefault(r["source"], [0, 0])
        t[0] += 1
        t[1] += int(r["start_w"] > r["cp"]["weight"])
    for src_name, (tot, nd) in sorted(by_src.items()):
        print(f"    {src_name:<24}: needed in {nd}/{tot}  ({100*nd/tot:.0f}%)")
    sr = [r["sa"]["success_rate"] for r in results]
    print(f"  SA restarts reaching the optimum    : mean {100*sum(sr)/len(sr):.0f}%, "
          f"min {100*min(sr):.0f}%")
    verdicts = Counter(r["assess"]["verdict"].split(":")[0] for r in results)
    print("\n  Verdicts on the move from the shift-only start to the optimum:")
    for v_name, v_n in sorted(verdicts.items(), key=lambda kv: -kv[1]):
        print(f"    {v_name:<20}: {v_n}/{n}")
    improved = [r for r in results if r["start_w"] > r["cp"]["weight"]]
    if improved:
        dm = sorted(r["assess"]["d_mae"] for r in improved)
        print(f"  Where the score improved (n={len(improved)}): change in MAE median "
              f"{dm[len(dm)//2]:+.2f} ppm, max {dm[-1]:+.2f} ppm")
    for r in results:
        if r["assess"]["verdict"].startswith("FLAG"):
            print(f"    ! {r['name']}: {r['assess']['verdict']}")
    n_rule = sum(1 for r in results if r["assess"].get("selection_rule_differs"))
    print(f"  Selection rule picked a different assignment from the lowest (MAE, LAE): {n_rule}/{n}")
    if n_complete < n:
        print("  * = optimal set truncated (cap or time limit): the CP-SAT MAE is the "
              "best among the solutions found,\n      not necessarily among all optimal "
              "solutions. Raise --max-solutions / --enum-limit.")
    print(f"{'═'*110}\n")


def save_csv(results: List[dict], path: Path) -> None:
    """Write summary results to a CSV file."""
    import csv
    rows = []
    for r in results:
        sa = r["sa"]; cp = r["cp"]
        speedup = sa["time"] / cp["time"] if cp["time"] > 0 else float("inf")
        rows.append({
            "molecule":              r["name"],
            "smiles":                r["smiles"],
            "n_carbons":             r["n_carbons"],
            "ch3":                   r["ch3"],
            "ch2":                   r["ch2"],
            "ch1":                   r["ch1"],
            "c_quaternary":          r["c"],
            "n_cosy":                r["n_cosy"],
            "n_hmbc":                r["n_hmbc"],
            "symmetry":              r["symmetry"],
            "predicted_weight":      r["predicted_w"],
            "sa_weight":             sa["weight"],
            "cpsat_weight":          cp["weight"],
            "score_difference":      sa["weight"] - cp["weight"],
            "sa_time_s":             round(sa["time"], 3),
            "cpsat_time_s":          round(cp["time"], 3),
            "cpsat_phase1_s":        round(cp.get("phase1_time", 0), 3),
            "cpsat_phase2_s":        round(cp.get("phase2_time", 0), 3),
            "speedup_sa_over_cpsat": round(speedup, 2),
            "sa_unique_solutions":   sa["unique_solutions"],
            "cpsat_unique_solutions":cp["unique_solutions"],
            "sa_mae_ppm":            round(sa["mae"], 4),
            "cpsat_mae_ppm":         round(cp["mae"], 4),
            "sa_lae_ppm":            round(sa["lae"], 4),
            "cpsat_lae_ppm":         round(cp["lae"], 4),
            "calc_source":           r["source"],
            "shift_only_start_weight": r["start_w"],
            "shift_only_needs_optimisation": r["start_w"] > cp["weight"],
            "peaks_differing_from_shift_only": r["start_n_differ"],
            "sa_restart_success_rate": round(sa["success_rate"], 3),
            "shift_only_mae_ppm":    round(r["assess"]["mae0"], 4),
            "shift_only_lae_ppm":    round(r["assess"]["lae0"], 4),
            "optimum_mae_ppm":       round(r["assess"]["mae1"], 4),
            "optimum_lae_ppm":       round(r["assess"]["lae1"], 4),
            "delta_mae_ppm":         round(r["assess"]["d_mae"], 4),
            "delta_lae_ppm":         round(r["assess"]["d_lae"], 4),
            "carbons_moved_to_optimum": r["assess"]["n_moved"],
            "carbons_moved_forced":  r["assess"]["n_forced"],
            "carbons_moved_no_correlation": r["assess"]["n_moved_no_corr"],
            "carbons_ambiguous_among_ties": r["assess"]["n_ambiguous"],
            "coverage_ge1_correlation": round(r["assess"]["cov1"], 3),
            "coverage_ge2_correlations": round(r["assess"]["cov2"], 3),
            "move_verdict":          r["assess"]["verdict"],
            "selection_rule_differs": r["assess"].get("selection_rule_differs"),
            "cpsat_proven_optimal":  cp["proven_optimal"],
            "cpsat_set_complete":    cp.get("enum_complete"),
            "cpsat_rescore_mismatch":cp.get("weight_mismatch"),
            "cpsat_optimality_gap":  round(cp.get("gap", 0) * 100, 2),
            "auto_fixes":            "; ".join(r["fixes"]) if r["fixes"] else "",
        })

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"  CSV saved to: {path}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        prog="benchmark_dataset",
        description="Compare SA and CP-SAT on a dataset of simpleNMR JSON files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "inputs", nargs="+", type=Path,
        help="JSON files or a folder containing JSON files.",
    )
    parser.add_argument(
        "--runs", type=int, default=SIMULATED_ANNEALING_RUNS,
        help="Number of SA restarts.",
    )
    parser.add_argument(
        "--max-solutions", type=int, default=2000,
        help="Cap on optimal CP-SAT solutions enumerated (rows marked * if hit).",
    )
    parser.add_argument(
        "--max-dmae", type=float, default=1.0,
        help="Verdict threshold: largest acceptable MAE increase (ppm) when moving "
             "from the shift-only start to the correlation optimum.",
    )
    parser.add_argument(
        "--max-dlae", type=float, default=10.0,
        help="Verdict threshold: largest acceptable increase (ppm) in the largest "
             "single-carbon error.",
    )
    parser.add_argument(
        "--time-limit", type=float, default=60.0,
        help="CP-SAT Phase 1 time limit per molecule (seconds).",
    )
    parser.add_argument(
        "--enum-limit", type=float, default=30.0,
        help="CP-SAT Phase 2 enumeration limit per molecule (seconds).",
    )
    parser.add_argument(
        "--csv", type=Path, default=None,
        help="Save summary results to this CSV file.",
    )
    parser.add_argument(
        "--verbose", action="store_true", default=False,
        help="Print CP-SAT solver log for each molecule.",
    )

    args = parser.parse_args()

    # Collect JSON files
    json_files: List[Path] = []
    for inp in args.inputs:
        if inp.is_dir():
            json_files.extend(sorted(inp.glob("*.json")))
        elif inp.suffix == ".json":
            json_files.append(inp)
        else:
            print(f"  Skipping {inp} (not a JSON file or directory)")

    if not json_files:
        print("No JSON files found. Exiting.")
        sys.exit(1)

    print(f"\nSimpleNMR SA vs CP-SAT Benchmark")
    print(f"{'─'*50}")
    print(f"  Molecules     : {len(json_files)}")
    print(f"  SA runs       : {args.runs}  "
          f"(cooling rate: {SIMULATED_ANNEALING_COOLING_RATE})")
    print(f"  CP-SAT Ph1    : {args.time_limit}s per molecule")
    print(f"  CP-SAT Ph2    : {args.enum_limit}s per molecule")
    print()

    results = []
    for i, path in enumerate(json_files, 1):
        print(f"[{i}/{len(json_files)}] {path.name} ...", end="  ", flush=True)
        try:
            r = process_file(
                path, args.runs, args.time_limit,
                args.enum_limit, args.verbose, args.max_solutions,
                args.max_dmae, args.max_dlae,
            )
            results.append(r)
            sa = r["sa"]; cp = r["cp"]
            speedup = sa["time"] / cp["time"] if cp["time"] > 0 else float("inf")
            diff = sa["weight"] - cp["weight"]
            print(
                f"SA={sa['weight']:.0f}/{sa['time']:.1f}s  "
                f"CP={cp['weight']:.0f}/{cp['time']:.2f}s  "
                f"{'CP '+f'{speedup:.0f}×' if speedup>=1 else 'SA faster'}  "
                f"{'⚠ SA wrong!' if diff > 0 else 'agree'}"
            )
        except Exception as e:
            print(f"ERROR: {e}")
            import traceback; traceback.print_exc()

    if not results:
        print("No results to display.")
        sys.exit(1)

    print_results_table(results)

    if args.csv:
        save_csv(results, args.csv)


if __name__ == "__main__":
    main()
