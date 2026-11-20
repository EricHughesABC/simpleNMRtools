"""
cpsat_nmr_v2.py
===============
CP-SAT exact solver for NMR ¹³C chemical-shift assignment.

Drop-in replacement for the production ``simulated_annealing_claude.py``.
Built against the current production codebase (v6 symmetry handling).

Design overview
---------------
The problem is a Quadratic Assignment Problem (QAP):

  Find a bijection  assign: nmr_node → mol_node
  (restricted to same-numProtons groups, restricted to orig_nodes_displayed)
  that minimises compute_total_weight(carbon_graph, mapping, shortest_paths).

CP-SAT replaces the 50-run SA sampling strategy with a provably-exact
two-phase approach:

  Phase 1 — Prove the global minimum weight W*.
      A CP-SAT model with binary assignment variables and McCormick-
      linearised penalty terms finds the optimal mapping and proves no
      mapping can score lower than W*.

  Phase 2 — Enumerate every mapping that achieves W*.
      The objective is fixed to W* and CP-SAT's enumerate_all_solutions
      mode collects every distinct optimal mapping via a solution callback.
      These are stored in self.results[W*]["best_mapping"], exactly as SA
      would after 50 runs.

Stage 2 (MAE/LAE tie-break) is handled by the *unchanged* process_results()
method inherited from SimulatedAnnealing2, which reads the same results
structure.

Key production-codebase compatibility points
---------------------------------------------
* Uses the same orig_nodes_displayed restriction as run_optimization in SA:
  only nodes in the displayed set, in groups with ≥ 2 members, are
  assignable.  All other nodes are fixed to the identity mapping.
* Populates self.results, self.bestest_weight, self.bestest_mapping,
  self.bestest_stats, self.bestest_weights, self.predicted_weight,
  self.predicted_mapping, self.initial_mapping with the exact types and
  shapes process_results() expects.
* verbose defaults to False — no CP-SAT search log printed unless requested,
  so timing benchmarks are not skewed.

Usage (identical API to SimulatedAnnealing2)
--------------------------------------------
    solver = CpSatNMR.from_json_file(
        json_file,
        orig_nodes_displayed=orig_nodes_displayed,
        pool_ch3_ch1=pool_ch3_ch1,
    )
    solver.setup_run(time_limit=60.0, num_workers=8)
    solver.run_optimization()          # num_times ignored — see docstring
    json_out, best = solver.process_results(catoms_df, jsonGraphData)

Command-line
------------
    python cpsat_nmr_v2.py molecule.json
    python cpsat_nmr_v2.py molecule.json --time-limit 120 --workers 8 --verbose
    python cpsat_nmr_v2.py molecule.json --output result.json
"""

import copy
import json
import math
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import networkx as nx
from loguru import logger
from ortools.sat.python import cp_model

# ── Import everything from the production SA module ───────────────────────────
# We subclass SimulatedAnnealing2 to inherit __init__, from_json_file,
# from_params, create_graph_df, calc_xy3_coords, initialize_mapping,
# map_protons_to_nodes, process_results, and process_results_SA_skipped
# completely unchanged.
from simulated_annealing_claude import (
    SimulatedAnnealing2,
    compute_total_weight,
)


# ─────────────────────────────────────────────────────────────────────────────
# Solution callback — collects every solution found during enumeration
# ─────────────────────────────────────────────────────────────────────────────

class _AllSolutionsCallback(cp_model.CpSolverSolutionCallback):
    """
    Collect every feasible solution found by the solver into a list.

    Used in Phase 2 (enumeration) when enumerate_all_solutions=True.
    Each solution is stored as a mapping dict {nmr_node: mol_node} using
    the full carbon node set (identity for non-assignable nodes, solver
    value for assignable ones).

    Parameters
    ----------
    assign_vars : dict
        {(n, m): BoolVar} — only for assignable nodes.
    assignable_groups : dict
        {numProtons: [node, ...]} — the reduced grouped_nodes.
    full_identity : dict
        {node: node} — identity mapping for ALL carbon nodes, including
        those not in the assignable set.  We copy this and overlay the
        solver's assignments, so every mapping returned is complete.
    max_solutions : int
        Stop collecting after this many solutions (safety cap).
    """

    def __init__(
        self,
        assign_vars: Dict[Tuple[int, int], cp_model.IntVar],
        assignable_groups: Dict[int, List[int]],
        full_identity: Dict[int, int],
        max_solutions: int = 500,
    ) -> None:
        super().__init__()
        self._assign_vars    = assign_vars
        self._groups         = assignable_groups
        self._full_identity  = full_identity
        self._max_solutions  = max_solutions
        self.solutions: List[Dict[int, int]] = []
        self._count = 0
        self._seen = set()
        self.duplicates = 0   # times CP-SAT reported an assignment already seen

    def on_solution_callback(self) -> None:
        # Build a full mapping: start from identity, overlay solver values
        mapping = self._full_identity.copy()
        for nprotons, nodes in self._groups.items():
            for n in nodes:
                for m in nodes:
                    if self.Value(self._assign_vars[(n, m)]):
                        mapping[n] = m
                        break

        key = tuple(sorted(mapping.items()))
        if key in self._seen:            # CP-SAT can report a solution twice
            self.duplicates += 1
            return
        if self._count >= self._max_solutions:
            self.StopSearch()
            return

        self._seen.add(key)
        self.solutions.append(mapping)
        self._count += 1

    @property
    def solution_count(self) -> int:
        return self._count


# ─────────────────────────────────────────────────────────────────────────────
# Internal model builder
# ─────────────────────────────────────────────────────────────────────────────

def _build_model(
    carbon_graph: nx.Graph,
    shortest_paths: Dict[int, Dict[int, int]],
    assignable_groups: Dict[int, List[int]],
    nProtons_reduced: Dict[int, int],
    phase2_enumeration: bool = False,
) -> Tuple[cp_model.CpModel, Dict[Tuple[int, int], cp_model.IntVar], cp_model.LinearExprT]:
    """
    Build the CP-SAT model for the NMR assignment QAP.

    Variables
    ---------
    assign[(n, m)] ∈ {0,1}
        1 iff NMR node n is assigned to mol node m.
        Only created for (n, m) pairs within the same assignable group.

    pair[(u, v, mu, mv)] ∈ {0,1}
        = assign[(u,mu)] AND assign[(v,mv)].
        Only created for correlation edges with non-zero penalty.
        Enforced by three McCormick implications (tight for binary vars).

    Objective
    ---------
    Minimise  Σ  penalty(mu, mv, edge_type) × pair[(u,v,mu,mv)]

    where penalty mirrors compute_total_weight exactly:
        COSY : (shortest_path − 1)³
        HMBC : (max(shortest_path, 3) − 2)³   [floor path at 3 for HMBC]

    Parameters
    ----------
    carbon_graph   : NMR correlation graph (COSY/HMBC edges).
    shortest_paths : mol-graph all-pairs shortest paths.
    assignable_groups : {numProtons: [node, ...]} — reduced pool.
    nProtons_reduced  : {node: numProtons} — reduced pool.

    Returns
    -------
    model, assign_vars, objective_expression
    """

    model = cp_model.CpModel()

    # ── 1. Assignment variables ───────────────────────────────────────────────
    assign: Dict[Tuple[int, int], cp_model.IntVar] = {}
    for nprotons, nodes in assignable_groups.items():
        for n in nodes:
            for m in nodes:
                assign[(n, m)] = model.new_bool_var(f"x_{n}_{m}")

    # ── 2. Bijection constraints within each group ────────────────────────────
    for nprotons, nodes in assignable_groups.items():
        # Each NMR node maps to exactly one mol node
        for n in nodes:
            model.add_exactly_one(assign[(n, m)] for m in nodes)
        # Each mol node receives exactly one NMR node
        for m in nodes:
            model.add_exactly_one(assign[(n, m)] for n in nodes)

    # ── 3. Pre-compute penalties for each (edge, assignment-pair) triple ──────
    node_to_np: Dict[int, int] = {
        int(n): int(p) for n, p in nProtons_reduced.items()
    }

    def edge_penalty(is_cosy, is_hmbc, path):
        """Identical to compute_total_weight in the production SA code."""
        penalty = 0
        if is_cosy:
            penalty += (path - 1) ** 3
        if is_hmbc:
            eff = path if path >= 3 else 2
            penalty += (eff - 2) ** 3
        return penalty

    # pair_key → BoolVar, built lazily below
    pair_vars: Dict[Tuple, cp_model.IntVar] = {}
    objective_terms = []
    constant = 0     # penalty from correlations between two FIXED nodes

    for u, v, d in carbon_graph.edges(data=True):
        is_cosy = bool(d.get("cosy"))
        is_hmbc = bool(d.get("hmbc"))
        if not (is_cosy or is_hmbc):
            continue

        u_free, v_free = u in node_to_np, v in node_to_np

        if u_free and v_free:
            # Both ends can move: penalty depends on the PAIR of assignments
            for mu in assignable_groups[node_to_np[u]]:
                for mv in assignable_groups[node_to_np[v]]:
                    penalty = edge_penalty(is_cosy, is_hmbc, shortest_paths[mu][mv])
                    if penalty == 0:
                        continue  # zero-penalty pair: no variable needed

                    # Auxiliary product variable: pair = assign[u,mu] AND assign[v,mv]
                    key = (u, v, mu, mv)
                    if key not in pair_vars:
                        pv = model.new_bool_var(f"p_{u}_{v}_{mu}_{mv}")
                        pair_vars[key] = pv
                        a_u = assign[(u, mu)]
                        a_v = assign[(v, mv)]
                        # Linearisation (exact for binary variables):
                        model.add_bool_and([a_u, a_v]).only_enforce_if(pv)
                        model.add_bool_or([a_u.Not(), a_v.Not()]).only_enforce_if(
                            pv.Not()
                        )
                    objective_terms.append(penalty * pair_vars[key])

        elif u_free or v_free:
            # One end is fixed (it maps to itself): the penalty depends on ONE
            # assignment, so no product variable is needed.  These correlations
            # were previously dropped, which let CP-SAT report assignments that
            # the production scorer rates worse than W*.
            free, fixed = (u, v) if u_free else (v, u)
            for mf in assignable_groups[node_to_np[free]]:
                penalty = edge_penalty(is_cosy, is_hmbc, shortest_paths[mf][fixed])
                if penalty:
                    objective_terms.append(penalty * assign[(free, mf)])

        else:
            # Both ends fixed: a constant, but part of the reported score
            constant += edge_penalty(is_cosy, is_hmbc, shortest_paths[u][v])

    # ── 4. Objective ──────────────────────────────────────────────────────────
    objective_expr = sum(objective_terms) + constant
    if not phase2_enumeration:
        model.minimize(objective_expr)

    return model, assign, objective_expr


# ─────────────────────────────────────────────────────────────────────────────
# CpSatNMR — drop-in class
# ─────────────────────────────────────────────────────────────────────────────

class CpSatNMR(SimulatedAnnealing2):
    """
    CP-SAT exact solver for NMR carbon assignment.

    Drop-in replacement for SimulatedAnnealing2 in production.

    All graph construction, data loading and result-processing methods are
    inherited unchanged.  Only setup_run() and run_optimization() differ.

    Two-phase strategy
    ------------------
    Phase 1 (prove the minimum):
        CP-SAT finds and PROVES the globally optimal COSY/HMBC penalty
        score W*.  No temperature schedule.  No restarts needed.

    Phase 2 (enumerate all optimal mappings):
        The model is re-solved with the objective fixed to W* and
        enumerate_all_solutions=True.  Every distinct mapping that achieves
        W* is collected and stored in self.results[W*]["best_mapping"].

    Stage 2 (chemical shift tie-break) is handled by process_results(),
    which is inherited from SimulatedAnnealing2 without modification.

    run_optimization(num_times)
    ---------------------------
    The num_times argument is accepted for API compatibility.  For CP-SAT
    it controls the maximum number of optimal solutions to enumerate in
    Phase 2 (default 500).  Passing num_times=1 finds only the single
    best mapping; passing num_times=50 (the SA default) enumerates up to
    50 optimal solutions for the MAE/LAE tie-break.

    Checking results
    ----------------
    After run_optimization(), inspect:
        solver.bestest_stats["is_optimal"]      — True = proven global min
        solver.bestest_stats["optimality_gap"]  — 0.0 when optimal
        solver.bestest_stats["n_optimal_solutions"] — solutions found in Ph2
        solver.bestest_stats["phase1_time_s"]   — time to prove W*
        solver.bestest_stats["phase2_time_s"]   — time to enumerate
    """

    def setup_run(
        self,
        time_limit: float = 60.0,
        num_workers: int = 8,
        enumerate_time_limit: float = 30.0,
        max_solutions: int = 500,
        verbose: bool = False,
    ) -> None:
        """
        Configure CP-SAT parameters and pre-compute shared structures.

        Parameters
        ----------
        time_limit : float
            Wall-clock budget (seconds) for Phase 1 (proving W*).
            60 s is sufficient for most molecules at stigmaSterol's scale.
            Increase for larger molecules.
        num_workers : int
            Parallel CP-SAT threads.  Set to your CPU core count.
        enumerate_time_limit : float
            Wall-clock budget (seconds) for Phase 2 (enumerating all
            optimal mappings).  30 s is usually ample — the hard work
            was done in Phase 1.
        max_solutions : int
            Maximum number of optimal mappings to collect in Phase 2.
            Equivalent to the num_times argument in SA's run_optimization.
            500 is a safe ceiling; very few molecules have more than a
            handful of truly distinct optimal mappings.
        verbose : bool
            Print the CP-SAT branch-and-bound search log.  Default False
            to avoid skewing timing benchmarks.
        """
        self.time_limit           = time_limit
        self.num_workers          = num_workers
        self.enumerate_time_limit = enumerate_time_limit
        self.max_solutions        = max_solutions
        self.verbose              = verbose

        # ── Structures shared with SA (same as SA's setup_run) ────────────────
        self.shortest_paths = dict(
            nx.all_pairs_shortest_path_length(self.mol_graph)
        )
        self.nProtons_to_nodes = self.map_protons_to_nodes(self.graph_df)

        self.predicted_mapping = self.initialize_mapping(
            self.carbon_graph,
            self.carbon_grouped_nodes,
            randomize_mapping=False,
            displayed_nodes=self.orig_nodes_displayed,
        )
        self.predicted_weight = compute_total_weight(
            self.carbon_graph, self.predicted_mapping, self.shortest_paths
        )

        # ── Build the same restricted pools SA uses in run_optimization ───────
        # Only nodes in orig_nodes_displayed, in groups with ≥ 2 members.
        displayed_set = set(self.orig_nodes_displayed)

        carbon_grouped_nodes_displayed = {
            nprotons: [n for n in nodes if n in displayed_set]
            for nprotons, nodes in self.carbon_grouped_nodes.items()
        }
        self._assignable_groups = {
            k: v
            for k, v in carbon_grouped_nodes_displayed.items()
            if len(v) > 1
        }
        self._nProtons_reduced = {
            k: v
            for k, v in self.nProtons_to_nodes.items()
            if (k in displayed_set) and (int(v) in self._assignable_groups)
        }

        # Full identity mapping over all carbon nodes — used as the base
        # for building complete mappings from partial CP-SAT solutions.
        self._full_identity = {
            node: node for node in self.carbon_graph.nodes()
        }

        n_assignable = sum(len(v) for v in self._assignable_groups.values())
        n_total      = len(self.carbon_graph.nodes())
        n_cosy = sum(1 for _, _, d in self.carbon_graph.edges(data=True) if d.get("cosy"))
        n_hmbc = sum(1 for _, _, d in self.carbon_graph.edges(data=True) if d.get("hmbc"))

        logger.info(
            f"CP-SAT setup: {n_total} carbons, {n_assignable} assignable "
            f"({n_total - n_assignable} fixed to identity) | "
            f"COSY={n_cosy} HMBC={n_hmbc} | "
            f"time_limit={time_limit}s workers={num_workers}"
        )
        logger.info(f"  Predicted (identity) weight: {self.predicted_weight:.4f}")

    def run_optimization(self, num_times: int = 50) -> None:
        """
        Run the two-phase CP-SAT optimisation.

        Parameters
        ----------
        num_times : int
            Maximum number of optimal solutions to enumerate in Phase 2.
            Matches the SA convention of passing the number of restarts
            (50 by default).  Pass 1 to retrieve only the single best
            mapping without enumeration.
        """
        self.initial_mapping = self._full_identity.copy()
        self.initial_weight  = compute_total_weight(
            self.carbon_graph, self.initial_mapping, self.shortest_paths
        )
        self.results        = {}
        self.bestest_weight = float("inf")

        # ── Phase 1: Prove the global minimum weight W* ───────────────────────
        logger.info("CP-SAT Phase 1: proving global minimum …")
        t0 = time.perf_counter()

        model, assign_vars, _ = _build_model(
            self.carbon_graph,
            self.shortest_paths,
            self._assignable_groups,
            self._nProtons_reduced,
        )

        solver1 = cp_model.CpSolver()
        solver1.parameters.max_time_in_seconds  = self.time_limit
        solver1.parameters.num_workers          = self.num_workers
        solver1.parameters.log_search_progress  = self.verbose

        status1 = solver1.solve(model)
        phase1_time = time.perf_counter() - t0

        is_optimal = (status1 == cp_model.OPTIMAL)
        is_feasible = status1 in (cp_model.OPTIMAL, cp_model.FEASIBLE)

        if not is_feasible:
            logger.warning(
                f"CP-SAT Phase 1 returned '{solver1.status_name(status1)}' — "
                "no feasible solution found.  Try increasing --time-limit."
            )
            # Fall back to identity mapping so process_results() doesn't crash
            w = self.predicted_weight
            self._store_result(self._full_identity.copy(), w, {
                "is_optimal": False,
                "status": solver1.status_name(status1),
                "phase1_time_s": phase1_time,
                "phase2_time_s": 0.0,
                "n_optimal_solutions": 0,
                "optimality_gap": float("inf"),
            })
            return

        optimal_weight = int(round(solver1.objective_value))
        obj_bound      = int(round(solver1.best_objective_bound))
        gap = (
            abs(optimal_weight - obj_bound) / max(abs(optimal_weight), 1)
            if optimal_weight > 0 else 0.0
        )

        logger.info(
            f"  Phase 1 complete in {phase1_time:.2f}s | "
            f"W*={optimal_weight} | "
            f"proven={'yes' if is_optimal else 'no (time limit)'} | "
            f"gap={gap*100:.2f}%"
        )

        # ── Phase 2: Enumerate all mappings at W* ─────────────────────────────
        if num_times <= 1:
            # Caller only wants the single best — extract from Phase 1 directly
            mapping = self._extract_mapping(solver1, assign_vars)
            self._store_result(mapping, optimal_weight, {
                "is_optimal": is_optimal,
                "status": solver1.status_name(status1),
                "phase1_time_s": phase1_time,
                "phase2_time_s": 0.0,
                "n_optimal_solutions": 1,
                "optimality_gap": gap,
            })
            return

        logger.info(
            f"CP-SAT Phase 2: enumerating all mappings at W*={optimal_weight} "
            f"(max {num_times}, budget {self.enumerate_time_limit}s) …"
        )
        t1 = time.perf_counter()

        # Rebuild model, add constraint obj == W*, enumerate all solutions
        model2, assign_vars2, obj_expr2 = _build_model(
            self.carbon_graph,
            self.shortest_paths,
            self._assignable_groups,
            self._nProtons_reduced,
            phase2_enumeration=True,   # no minimize — pure feasibility at W*
        )
        # Fix objective to proven optimal value
        model2.add(obj_expr2 == optimal_weight)

        callback = _AllSolutionsCallback(
            assign_vars2,
            self._assignable_groups,
            self._full_identity,
            max_solutions=num_times,
        )

        solver2 = cp_model.CpSolver()
        solver2.parameters.max_time_in_seconds   = self.enumerate_time_limit
        solver2.parameters.num_workers           = 1   # deterministic enumeration
        solver2.parameters.enumerate_all_solutions = True
        solver2.parameters.log_search_progress   = self.verbose

        status2 = solver2.solve(model2, callback)
        phase2_time = time.perf_counter() - t1

        solutions = callback.solutions
        enum_complete = (status2 == cp_model.OPTIMAL) and len(solutions) < num_times
        if not enum_complete:
            logger.warning(
                f"Optimal set may be INCOMPLETE ({len(solutions)} found; cap "
                f"{num_times}, time budget {self.enumerate_time_limit}s). The best "
                f"MAE is the best among these only. Raise --num-solutions / "
                f"--enum-time-limit."
            )
        if not solutions:
            # Callback collected nothing — fall back to Phase 1 solution
            solutions = [self._extract_mapping(solver1, assign_vars)]

        logger.info(
            f"  Phase 2 complete in {phase2_time:.2f}s | "
            f"{len(solutions)} optimal solution(s) found"
        )

        stats = {
            "is_optimal": is_optimal,
            "status": solver1.status_name(status1),
            "phase1_time_s": phase1_time,
            "phase2_time_s": phase2_time,
            "n_optimal_solutions": len(solutions),
            "optimality_gap": gap,
            "enum_complete": enum_complete,
            "duplicates_ignored": callback.duplicates,
        }

        # Store all optimal mappings under the proven weight
        for mapping in solutions:
            self._store_result(mapping, optimal_weight, stats)

        logger.info(
            f"── CP-SAT results ──────────────────────────────\n"
            f"  Predicted weight       : {self.predicted_weight:.4f}\n"
            f"  Optimal weight (W*)    : {self.bestest_weight:.4f}\n"
            f"  Proven optimal         : {is_optimal}\n"
            f"  Optimality gap         : {gap*100:.2f}%\n"
            f"  Optimal solutions found: {len(solutions)}\n"
            f"  Phase 1 time           : {phase1_time:.2f}s\n"
            f"  Phase 2 time           : {phase2_time:.2f}s"
        )

    # ── Internal helpers ──────────────────────────────────────────────────────

    def _extract_mapping(
        self,
        solver: cp_model.CpSolver,
        assign_vars: Dict[Tuple[int, int], cp_model.IntVar],
    ) -> Dict[int, int]:
        """
        Extract a complete mapping from a solved CP-SAT model.

        Starts from the full identity and overlays the solver's assignment
        for every node in the assignable groups.
        """
        mapping = self._full_identity.copy()
        for nprotons, nodes in self._assignable_groups.items():
            for n in nodes:
                for m in nodes:
                    if solver.value(assign_vars[(n, m)]):
                        mapping[n] = m
                        break
        return mapping

    def _store_result(
        self,
        mapping: Dict[int, int],
        weight: float,
        stats: Dict,
    ) -> None:
        """
        Store a (mapping, weight) pair in self.results using the exact
        structure expected by process_results().

        The 'weights' history list uses the same 4-tuple as SA:
            (iteration, best_weight, current_weight, temperature)
        with temperature=0.0 (not applicable to CP-SAT).
        """
        history = [(0, weight, weight, 0.0)]

        if weight not in self.results:
            self.results[weight] = {
                "num_times"   : 0,
                "results"     : [],
                "best_mapping": [],
            }

        self.results[weight]["num_times"]    += 1
        self.results[weight]["results"].append((mapping, stats, history))
        self.results[weight]["best_mapping"].append(mapping)

        if weight < self.bestest_weight:
            self.bestest_weight  = weight
            self.bestest_mapping = mapping
            self.bestest_stats   = stats
            self.bestest_weights = history


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    """
    Command-line entry point — mirrors simulated_annealing_claude.py CLI.

    Usage
    -----
    python cpsat_nmr_v2.py molecule.json
    python cpsat_nmr_v2.py molecule.json --time-limit 120 --workers 8
    python cpsat_nmr_v2.py molecule.json --num-solutions 50 --output result.json
    python cpsat_nmr_v2.py molecule.json --verbose
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog="cpsat_nmr_v2",
        description="CP-SAT exact NMR assignment solver.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("json_file", type=Path,
                        help="Input NMR JSON file.")
    parser.add_argument("--orig-nodes-displayed", type=str, default=None,
                        metavar="ID,ID,...",
                        help="Comma-separated atom ids for the minimal "
                             "symmetry-consistent node set.")
    parser.add_argument("--pool-ch3-ch1", action="store_true", default=None,
                        help="Provenance flag — does not affect swap grouping.")
    parser.add_argument("--time-limit", type=float, default=60.0,
                        help="Phase 1 wall-clock budget (seconds).")
    parser.add_argument("--enum-time-limit", type=float, default=30.0,
                        help="Phase 2 enumeration budget (seconds).")
    parser.add_argument("--workers", type=int, default=8,
                        help="Parallel CP-SAT threads.")
    parser.add_argument("--num-solutions", type=int, default=2000,
                        help="Cap on optimal solutions enumerated (Phase 2). A "
                             "warning is logged if the optimal set is truncated.")
    parser.add_argument("--verbose", action="store_true", default=False,
                        help="Print CP-SAT search log.")
    parser.add_argument("--output", type=Path, default=None,
                        help="Save optimised graph JSON to this file.")

    args = parser.parse_args()

    if not args.json_file.exists():
        parser.error(f"File not found: {args.json_file}")

    orig_nodes_displayed = (
        [int(x.strip()) for x in args.orig_nodes_displayed.split(",") if x.strip()]
        if args.orig_nodes_displayed is not None
        else None
    )

    logger.info(f"Loading NMR data from: {args.json_file}")
    solver = CpSatNMR.from_json_file(
        args.json_file,
        orig_nodes_displayed=orig_nodes_displayed,
        pool_ch3_ch1=args.pool_ch3_ch1,
    )

    solver.setup_run(
        time_limit           = args.time_limit,
        enumerate_time_limit = args.enum_time_limit,
        num_workers          = args.workers,
        max_solutions        = args.num_solutions,
        verbose              = args.verbose,
    )

    t0 = time.perf_counter()
    solver.run_optimization(num_times=args.num_solutions)
    wall = time.perf_counter() - t0

    s = solver.bestest_stats
    logger.info(
        f"Total wall time: {wall:.2f}s | "
        f"W*={solver.bestest_weight:.4f} | "
        f"optimal={s.get('is_optimal')} | "
        f"solutions={s.get('n_optimal_solutions')}"
    )

    if args.output is not None:
        with open(args.json_file) as f:
            json_graph_data = json.load(f)

        catoms_df = solver.carbon_df.copy()
        result_json, best_results = solver.process_results(
            catoms_df, json_graph_data
        )

        logger.info(
            f"MAE={best_results['best_mae']:.4f} ppm | "
            f"LAE={best_results['best_lae']:.4f} ppm "
            f"(worst atom #{best_results['best_lae_atomNumber']})"
        )

        with open(args.output, "w") as f:
            json.dump(result_json, f, indent=2)
        logger.info(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()
