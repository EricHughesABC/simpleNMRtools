"""
heldout_validation.py
=====================
Does minimising the HMBC/COSY score give an assignment that predicts
correlations the optimiser never saw?

Method (no reference assignments needed)
----------------------------------------
For each molecule the observed correlations are split into a training part and
a held-out part.  CP-SAT finds the optimal assignments for the training
correlations only and, exactly as the pipeline does, the one with the lowest
shift error (MAE, then LAE) is chosen among the tied optima.  That assignment
and the closest-shift (shift-only) starting assignment are then both scored on
the held-out correlations.  A lower held-out score for the optimised assignment
means the optimisation improved the assignment with respect to independent
measurements from the same sample, so concentration, pH and solvent effects on
the chemical shifts cancel out.

Two kinds of split are run:
  * k-fold      each correlation is held out exactly once per repeat
                (default 5 folds x 4 repeats = 20 train/test splits per molecule)
  * HMBC->COSY  train on HMBC only, test on COSY (an independent experiment)

Scores use the same cubic path-length penalty as the production code.  The
number of violated held-out correlations is reported as well, because it is
easier to interpret.

Usage (run from the simpleNMRtools folder, next to benchmark_cpsat_SA_dataset.py)
-------
    python heldout_validation.py ../simplecasenmr/examples
    python heldout_validation.py ../simplecasenmr/examples --verdicts results.csv --out heldout
    python heldout_validation.py file1.json file2.json --folds 5 --repeats 2

Outputs
-------
    <out>_splits.csv    one row per train/test split
    <out>_summary.csv   one row per molecule
"""

import argparse
import csv
import random
import sys
import time
import warnings
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import networkx as nx
import numpy as np
from ortools.sat.python import cp_model

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))

import benchmark_cpsat_SA_dataset as B          # loaders, CP-SAT model, shift-only start
from core.simulated_annealing_claude import compute_total_weight
from loguru import logger

logger.remove()
logger.add(sys.stderr, level="WARNING")


# ─────────────────────────────────────────────────────────────────────────────
# Scoring helpers
# ─────────────────────────────────────────────────────────────────────────────

def _f(v) -> float:
    """Float conversion that turns empty/non-numeric shifts into NaN."""
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def score_edges(mapping: Dict[int, int], edges: list, sp) -> Tuple[int, int]:
    """Cubic path-length penalty of ``mapping`` over ``edges``.

    Identical to compute_total_weight in the production SA code.  Returns
    (total penalty, number of correlations with a non-zero penalty).
    """
    total, violated = 0, 0
    for u, v, dd in edges:
        a, b = mapping[u], mapping[v]
        if a == b:
            continue
        path = sp[a][b]
        pen = 0
        if dd.get("cosy"):
            pen += (path - 1) ** 3
        if dd.get("hmbc"):
            eff = path if path >= 3 else 2
            pen += (eff - 2) ** 3
        total += pen
        violated += pen > 0
    return total, violated


class Molecule:
    """Everything needed to optimise and score one molecule."""

    def __init__(self, path: Path):
        d, _ = B.load_json(path)
        s = B.make_solver(d, randomize=False)
        self.name = path.stem
        self.s = s
        self.sp = s.shortest_paths
        self.nodes = list(s.carbon_graph.nodes())
        self.start = B.shift_only_start(d, s.carbon_df)        # closest-shift assignment
        disp = set(s.orig_nodes_displayed)
        ag = {k: [n for n in v if n in disp] for k, v in s.carbon_grouped_nodes.items()}
        self.ag = {k: v for k, v in ag.items() if len(v) > 1}
        self.nP = {k: v for k, v in s.nProtons_to_nodes.items()
                   if k in disp and int(v) in self.ag}
        self.identity = {n: n for n in self.nodes}
        self.edges = [(u, v, dict(dd)) for u, v, dd in s.carbon_graph.edges(data=True)
                      if dd.get("cosy") or dd.get("hmbc")]
        self.exp = {n["id"]: _f(n["ppm"]) for n in s.nmr_nodes}
        self.calc = {i: _f(v) for i, v in s.carbon_df["ppm_calculated"].items()}
        # safety check: this scorer must agree with the production one
        full, _ = score_edges(self.start, self.edges, self.sp)
        prod = compute_total_weight(s.carbon_graph, self.start, self.sp)
        if full != prod:
            print(f"  WARNING {self.name}: scorer {full} != production {prod}")

    def shift_error(self, mapping) -> Tuple[float, float]:
        e = [abs(self.exp[n] - self.calc[mapping[n]]) for n in self.exp
             if not (np.isnan(self.exp[n]) or np.isnan(self.calc.get(mapping[n], float("nan"))))]
        return (float(np.mean(e)), float(np.max(e))) if e else (float("inf"), float("inf"))

    def graph(self, edges: list) -> nx.Graph:
        g = nx.Graph()
        g.add_nodes_from(self.nodes)
        for u, v, dd in edges:
            g.add_edge(u, v, **dd)
        return g

    def optimise(self, train_edges: list, max_tied: int, workers: int):
        """Optimal assignments for ``train_edges``.

        Returns (optimal score, list of tied optimal assignments, truncated flag)
        or None if CP-SAT finds no solution.
        """
        g = self.graph(train_edges)
        m1, av1, _ = B._build_cpsat_model(g, self.sp, self.ag, self.nP)
        s1 = cp_model.CpSolver()
        s1.parameters.num_workers = workers
        s1.parameters.max_time_in_seconds = 30
        st = s1.solve(m1)
        if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            return None
        W = int(round(s1.objective_value))
        m2, av2, oe2 = B._build_cpsat_model(g, self.sp, self.ag, self.nP, phase2=True)
        m2.add(oe2 == W)
        cb = B._AllSolutionsCb(av2, self.ag, self.identity, max_tied)
        s2 = cp_model.CpSolver()
        s2.parameters.num_workers = 1                  # deterministic enumeration
        s2.parameters.max_time_in_seconds = 10
        s2.parameters.enumerate_all_solutions = True
        st2 = s2.solve(m2, cb)
        sols = cb.solutions
        truncated = not (st2 == cp_model.OPTIMAL and len(sols) < max_tied)
        if not sols:                                   # fall back to the Phase 1 solution
            m = self.identity.copy()
            for nk, ng in self.ag.items():
                for n in ng:
                    for mv in ng:
                        if s1.value(av1[(n, mv)]):
                            m[n] = mv
                            break
            sols, truncated = [m], True
        return W, sols, truncated

    def pick(self, sols: list) -> dict:
        """The pipeline's tie-break: lowest shift error (MAE, then LAE)."""
        return min(sols, key=self.shift_error)


# ─────────────────────────────────────────────────────────────────────────────
# Splits
# ─────────────────────────────────────────────────────────────────────────────

def kfold_splits(n_edges: int, folds: int, repeats: int, seed: int):
    k = max(2, min(folds, n_edges // 2))               # at least 2 correlations per fold
    for r in range(repeats):
        idx = list(range(n_edges))
        random.Random(seed + r).shuffle(idx)
        parts = [idx[i::k] for i in range(k)]
        for f in range(k):
            test = set(parts[f])
            yield "kfold", r, f, [i for i in idx if i not in test], sorted(test)


def run_split(mol: Molecule, mode: str, rep: int, fold: int,
              train_idx: list, test_idx: list, max_tied: int, workers: int) -> Optional[dict]:
    train = [mol.edges[i] for i in train_idx]
    test = [mol.edges[i] for i in test_idx]
    res = mol.optimise(train, max_tied, workers)
    if res is None:
        return None
    _, sols, truncated = res
    best = mol.pick(sols)
    tr0, _ = score_edges(mol.start, train, mol.sp)
    tr1, _ = score_edges(best, train, mol.sp)
    te0, v0 = score_edges(mol.start, test, mol.sp)
    te1, v1 = score_edges(best, test, mol.sp)
    return dict(molecule=mol.name, mode=mode, repeat=rep, fold=fold,
                n_train=len(train), n_test=len(test),
                train_start=tr0, train_opt=tr1, test_start=te0, test_opt=te1,
                violated_start=v0, violated_opt=v1, tied=len(sols), truncated=truncated,
                outcome="improved" if te1 < te0 else "worse" if te1 > te0 else "equal")


# ─────────────────────────────────────────────────────────────────────────────
# One molecule
# ─────────────────────────────────────────────────────────────────────────────

def validate(path: Path, folds: int, repeats: int, max_tied: int, seed: int,
             workers: int, experiment: bool) -> Tuple[dict, List[dict]]:
    mol = Molecule(path)
    full_start, _ = score_edges(mol.start, mol.edges, mol.sp)
    full = mol.optimise(mol.edges, max_tied, workers)
    full_opt = full[0] if full else float("nan")

    rows = []
    for mode, rep, fold, tr, te in kfold_splits(len(mol.edges), folds, repeats, seed):
        if len(tr) < 2 or not te:
            continue
        r = run_split(mol, mode, rep, fold, tr, te, max_tied, workers)
        if r:
            rows.append(r)

    if experiment:                                      # train on HMBC only, test on COSY
        hmbc = [i for i, e in enumerate(mol.edges) if e[2].get("hmbc") and not e[2].get("cosy")]
        cosy = [i for i, e in enumerate(mol.edges) if e[2].get("cosy")]
        if len(cosy) >= 3 and len(hmbc) >= 5:
            r = run_split(mol, "HMBC->COSY", 0, 0, hmbc, cosy, max_tied, workers)
            if r:
                rows.append(r)

    kf = [r for r in rows if r["mode"] == "kfold"]
    c = Counter(r["outcome"] for r in kf)
    exp_row = next((r for r in rows if r["mode"] == "HMBC->COSY"), None)
    summary = dict(
        molecule=mol.name, n_corr=len(mol.edges), full_start=full_start, full_opt=full_opt,
        needs_optimisation=full_start > full_opt, splits=len(kf),
        improved=c["improved"], equal=c["equal"], worse=c["worse"],
        median_test_start=float(np.median([r["test_start"] for r in kf])) if kf else float("nan"),
        median_test_opt=float(np.median([r["test_opt"] for r in kf])) if kf else float("nan"),
        violated_start=sum(r["violated_start"] for r in kf),
        violated_opt=sum(r["violated_opt"] for r in kf),
        truncated_splits=sum(r["truncated"] for r in kf),
        experiment_outcome=exp_row["outcome"] if exp_row else "n/a",
        experiment_test_start=exp_row["test_start"] if exp_row else "",
        experiment_test_opt=exp_row["test_opt"] if exp_row else "",
        generalises=("n/a" if not full_start > full_opt else
                     "yes" if c["improved"] > c["worse"] else "no"),
    )
    return summary, rows


# ─────────────────────────────────────────────────────────────────────────────
# Reporting
# ─────────────────────────────────────────────────────────────────────────────

def pct(a, b):
    return f"{100 * a / b:.0f}%" if b else "n/a"


def report(summaries: List[dict], splits: List[dict], verdicts: Dict[str, str]) -> None:
    from scipy.stats import binomtest

    need = [s for s in summaries if s["needs_optimisation"]]
    done = [s for s in summaries if not s["needs_optimisation"]]

    print(f"\n{'═' * 100}\n  HELD-OUT VALIDATION SUMMARY\n{'═' * 100}")
    print(f"  Molecules: {len(summaries)}  (optimisation needed on the full data: {len(need)}, "
          f"start already optimal: {len(done)})")

    def tot(rows):
        c = Counter()
        for s in rows:
            c.update(improved=s["improved"], equal=s["equal"], worse=s["worse"])
        return c["improved"], c["equal"], c["worse"]

    i, e, w = tot(need)
    n_sp = i + e + w
    print(f"\n  Molecules where optimisation was needed (n={len(need)}), "
          f"{n_sp} k-fold splits:")
    print(f"    held-out score improved {i} ({pct(i, n_sp)}), equal {e} ({pct(e, n_sp)}), "
          f"worse {w} ({pct(w, n_sp)})")
    vs = sum(s["violated_start"] for s in need)
    vo = sum(s["violated_opt"] for s in need)
    print(f"    held-out correlations violated: shift-only start {vs} -> optimised {vo} "
          f"({100 * (vs - vo) / vs:.0f}% fewer)" if vs else "    no violated held-out correlations")
    gen = [s for s in need if s["generalises"] == "yes"]
    print(f"    molecules where the optimised assignment predicts held-out correlations better "
          f"(improved > worse): {len(gen)}/{len(need)} ({pct(len(gen), len(need))})")
    nogen = [s for s in need if s["generalises"] != "yes"]
    for s in nogen:
        print(f"      ! {s['molecule']}: improved {s['improved']}, equal {s['equal']}, worse {s['worse']}")
    decided = [s for s in need if s["improved"] != s["worse"]]
    if decided:
        k = sum(1 for s in decided if s["improved"] > s["worse"])
        p = binomtest(k, len(decided), 0.5, alternative="greater").pvalue
        print(f"    molecule-level sign test (ties excluded): {k}/{len(decided)} better, one-sided p = {p:.2g}")
    print("    (splits within a molecule are not independent, so judge by molecules, not by splits)")

    if done:
        i2, e2, w2 = tot(done)
        print(f"\n  Start already optimal (n={len(done)}): held-out improved {i2}, equal {e2}, worse {w2}")
        for s in done:
            if s["worse"]:
                print(f"      ! {s['molecule']}: worse in {s['worse']} split(s) (training data pulled "
                      f"the optimum away from a start that is optimal on all the data)")

    exp_rows = [s for s in summaries if s["experiment_outcome"] != "n/a"]
    if exp_rows:
        c = Counter(s["experiment_outcome"] for s in exp_rows)
        print(f"\n  Independent-experiment split (train HMBC, test COSY), {len(exp_rows)} molecules "
              f"with enough COSY and HMBC peaks:")
        print(f"    improved {c['improved']}, equal {c['equal']}, worse {c['worse']}")

    if verdicts:
        print("\n  By the earlier move verdict (benchmark CSV), molecules needing optimisation:")
        groups = defaultdict(list)
        for s in need:
            groups[verdicts.get(s["molecule"], "unknown").split(":")[0]].append(s)
        for g_name, rows in sorted(groups.items()):
            gi, ge, gw = tot(rows)
            ng = sum(1 for s in rows if s["generalises"] == "yes")
            print(f"    {g_name:<10} n={len(rows):<3} molecules better on held-out data: {ng}/{len(rows)}"
                  f"   splits: improved {pct(gi, gi + ge + gw)}, worse {pct(gw, gi + ge + gw)}")
        for s in need:
            if verdicts.get(s["molecule"], "").startswith("FLAG"):
                print(f"      FLAGGED {s['molecule']:<46} held-out: improved {s['improved']}, "
                      f"equal {s['equal']}, worse {s['worse']}   ->  "
                      f"{'optimum still beats the start (suspect the shift predictions)' if s['generalises'] == 'yes' else 'does NOT generalise (check correlations / structure)'}")

    trunc = sum(s["truncated_splits"] for s in summaries)
    print(f"\n  Splits where the tied set was truncated at --max-tied: {trunc}/{len(splits)} "
          f"(the lowest-MAE pick is then among a sample of the tied optima)")
    print("═" * 100)


def write_csv(rows: List[dict], path: Path) -> None:
    if not rows:
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {path}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description="Held-out correlation validation of the "
                                             "HMBC/COSY optimisation.")
    ap.add_argument("inputs", nargs="+", type=Path, help="JSON files or folders of JSON files.")
    ap.add_argument("--folds", type=int, default=5, help="k for k-fold splitting.")
    ap.add_argument("--repeats", type=int, default=4, help="Repeats of the k-fold split.")
    ap.add_argument("--max-tied", type=int, default=500,
                    help="Cap on tied optimal assignments enumerated per split.")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=8, help="CP-SAT threads for the optimisation.")
    ap.add_argument("--verdicts", type=Path, default=None,
                    help="Benchmark results CSV; its move_verdict column is joined on molecule name.")
    ap.add_argument("--no-experiment-split", action="store_true",
                    help="Skip the HMBC->COSY split.")
    ap.add_argument("--limit", type=int, default=None, help="Only the first N molecules.")
    ap.add_argument("--out", type=str, default="heldout", help="Prefix for the output CSV files.")
    args = ap.parse_args()

    files: List[Path] = []
    for p in args.inputs:
        files.extend(sorted(p.glob("*.json")) if p.is_dir() else [p])
    if args.limit:
        files = files[:args.limit]
    if not files:
        sys.exit("No JSON files found.")

    verdicts = {}
    if args.verdicts:
        with open(args.verdicts, newline="") as fh:
            verdicts = {r["molecule"]: r.get("move_verdict", "") for r in csv.DictReader(fh)}

    print(f"\nHeld-out validation: {len(files)} molecules, {args.folds}-fold x {args.repeats} "
          f"repeats{'' if args.no_experiment_split else ' + HMBC->COSY'}\n")
    summaries, splits = [], []
    for i, path in enumerate(files, 1):
        t0 = time.perf_counter()
        print(f"[{i}/{len(files)}] {path.stem[:48]:<48}", end=" ", flush=True)
        try:
            s, rows = validate(path, args.folds, args.repeats, args.max_tied, args.seed,
                               args.workers, not args.no_experiment_split)
        except Exception as exc:                        # keep going; report at the end
            print(f"ERROR: {exc}")
            continue
        summaries.append(s)
        splits.extend(rows)
        print(f"{'needs opt' if s['needs_optimisation'] else 'start optimal':<13} "
              f"improved {s['improved']:>2} / equal {s['equal']:>2} / worse {s['worse']:>2}   "
              f"({time.perf_counter() - t0:.0f}s)")

    if not summaries:
        sys.exit("Nothing to report.")
    report(summaries, splits, verdicts)
    write_csv(splits, Path(f"{args.out}_splits.csv"))
    write_csv(summaries, Path(f"{args.out}_summary.csv"))


if __name__ == "__main__":
    main()
