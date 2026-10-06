"""Prediction pipeline for simpleNMR.

Runs the complete NMR assignment workflow:

  1. Build and validate the NMRsolution object
  2. Assign CH3 / CH2 / CH1 proton types
  3. Transfer HSQC information into C13 and H1 tables
  4. Assign carbons to the C13 table
  5. Build the NMR network graph and molecular graph
  6. Prepare the carbon-atoms DataFrame (catoms_df)
  7. Run (or skip) simulated annealing optimisation
  8. Assemble and return the Jinja2 template context dict

Public API
----------
``run(problemdata_json, json_data) -> dict``
    Execute all steps in order.
    Returns a Jinja2 context dict ready for ``render_template``.
    Raises ``PipelineError`` on any recoverable failure.
"""

import copy
import json

import networkx as nx
from networkx.readwrite import json_graph
from loguru import logger

import core.nmrsolution as nmrsolution
from config.globals import SVG_DIMENSIONS as svgDimensions
from core.simulated_annealing_claude import SimulatedAnnealing2
from app.pipelines import PipelineError
from config.globals import SIMULATED_ANNEALING_RUNS, SIMULATED_ANNEALING_COOLING_RATE


# ── Internal utilities ────────────────────────────────────────────────────────

def _node_link_data(G) -> dict:
    """Version-safe wrapper for json_graph.node_link_data.

    NX < 3.3  — no edges kwarg; default produces {"links": [...]}
    NX 3.3–3.5 — edges kwarg added; explicit "links" silences FutureWarning
    NX >= 3.6  — default changed to "edges"; explicit "links" preserves key

    Always produces {"links": [...]} so downstream code is stable across
    all versions. Revisit when upgrading server beyond NX 3.6.
    """
    try:
        return json_graph.node_link_data(G, edges="links")
    except TypeError:  # NX < 3.3 — edges kwarg not yet supported
        return json_graph.node_link_data(G)


# ── Pipeline steps ────────────────────────────────────────────────────────────

def _build_solution(problemdata_json):
    """Initialise NMRsolution and populate it from the JSON problem data.

    Returns
    -------
    NMRsolution
        A fully initialised solution object.

    Raises
    ------
    PipelineError
        If NMRsolution construction fails, the nmrshiftdb lookup fails, or
        ``init_class_from_json`` returns an error.
    """
    solution = nmrsolution.NMRsolution(problemdata_json)

    if solution.nmrsolution_failed:
        raise PipelineError(
            solution.nmrsolution_error_message,
            solution.nmrsolution_error_code,
        )

    if solution.expected_molecule.nmrshiftdb_failed:
        raise PipelineError(
            solution.nmrshiftdb_error_message,
            solution.nmrshiftdb_error_code,
        )

    ok, msg = solution.init_class_from_json()
    if not ok:
        raise PipelineError(msg, 400)

    return solution


def _assign_proton_types(solution) -> None:
    """Run the CH3 / CH2 / CH1 proton-type assignment.

    Raises
    ------
    PipelineError
        If the assignment returns a non-ok status.
    """
    rtn_msg, rtn_num = solution.assign_CH3_CH2_CH1_overall()
    if rtn_msg != "ok":
        raise PipelineError(rtn_msg, rtn_num)


def _transfer_hsqc_info(solution) -> None:
    """Copy HSQC peak information into the C13 and H1 tables."""
    solution.transfer_hsqc_info_to_c13()
    solution.transfer_hsqc_info_to_h1()


def _assign_carbons(solution, previous_positions: dict) -> None:
    """Initialise and run the C13 carbon-assignment step.

    ``previous_positions`` (``{atomNumber: ppm}``, see
    ``_previous_atom_assignments_lookup``) is threaded through to
    ``attempt_assignment_CH3_CH2_CH1_to_C13_table``, which -- if the user
    opted to reuse a previous iteration's positions -- force-assigns each
    matching experimental peak back onto the same fixed atom slot it had
    last time, ahead of the normal Hungarian/nearest-ppm matching.

    Raises
    ------
    PipelineError
        If either the initialisation or the assignment returns a non-ok status.
    """
    rtn_msg, rtn_num = solution.initialise_prior_to_carbon_assignment()
    if rtn_msg != "ok":
        raise PipelineError(rtn_msg, rtn_num)

    rtn_msg, rtn_num = solution.attempt_assignment_CH3_CH2_CH1_to_C13_table(
        previous_positions=previous_positions
    )
    if rtn_msg != "ok":
        raise PipelineError(rtn_msg, rtn_num)

    solution.update_assignments_expt_dataframes()


def _previous_atom_assignments_lookup(json_data: dict) -> dict:
    """Parse the optional ``previousNodePositions`` block into a lookup.

    ``json_data`` may carry this block in the same ``{datatype, count,
    data}`` shape used elsewhere in the mnova payload, where each ``data``
    entry is ``{"atomNumber": ..., "ppm": ...}``. This is the "which
    experimental peak is currently attached to which fixed atom slot"
    state exported from a previous iteration's HTML (``exportToMnova`` in
    export.js) and fed back in by ``simplePREDICT_eeh.qs``
    (``readPreviousNodePositions``).

    Why atomNumber, not ppm, is the key here (the reverse of an earlier,
    incorrect attempt at this feature): atomNumber/x/y/ppm_calculated are
    fixed structural properties of the molecule, straight from the molfile
    and its RDKit depiction -- the same every run. What varies is which
    experimental peak (ppm) gets *attached* to a given atom slot, decided
    fresh each run by attempt_assignment_CH3_CH2_CH1_to_C13_table's
    Hungarian/nearest-ppm matching, with no memory of a previous run. This
    mirrors exactly what dragging a node onto another and clicking
    "Update" does client-side (see updateMovedAtoms/commitMoves in
    templates/partials/js/moves.js): it moves ppm/H1_ppm/jCouplingVals
    onto the destination slot, while that slot's atomNumber/x/y/
    ppm_calculated never change. So "atom slot K had experimental peak P
    attached" is the stable, meaningful thing to replay -- not "peak P was
    drawn at position (x, y)".

    Returns
    -------
    dict
        ``{atomNumber: ppm}``, empty if the block is absent/empty or
        malformed. atomNumber is left as whatever type it arrives as
        (typically a str); matching against it is done downstream.
    """
    block = json_data.get("previousNodePositions")
    if not block or not block.get("data"):
        return {}

    positions = {}
    for entry in block["data"].values():
        atom_number = entry.get("atomNumber")
        if atom_number is None:
            continue
        try:
            positions[atom_number] = float(entry["ppm"])
        except (KeyError, TypeError, ValueError):
            continue

    return positions


def _build_graph(solution, json_data) -> tuple:
    """Build the NMR network graph and the molecular connectivity graph.

    Returns
    -------
    tuple[dict, dict, dict, int]
        ``(jsonGraphData, jsonGraphData_mol, shortest_paths, number_of_hmbc_cosy_subgraphs)``
        where ``jsonGraphData`` is the node-link dict for the NMR graph,
        ``jsonGraphData_mol`` is the node-link dict for the molecule graph,
        and ``shortest_paths`` is the all-pairs Dijkstra path-length dict.
    """
    G2 = nmrsolution.create_network_graph(solution.c13, solution.h1)
    G2 = nmrsolution.add_all_cosy_edges_to_graph(
        G2, solution.cosy, solution.hsqc_clipcosy, solution.h1
    )
    G2 = nmrsolution.add_all_hmbc_edges_to_graph(
        G2, solution.hmbc, solution.h1, solution.c13
    )

    jsonGraphData = _node_link_data(G2)
    jsonGraphData["moved_nodes"] = jsonGraphData["nodes"]

    solution.initiate_molgraph(json_data, G2)
    jsonGraphData_mol = _node_link_data(solution.molgraph)
    shortest_paths = dict(nx.all_pairs_dijkstra_path_length(solution.molgraph))

    if nx.is_connected(G2):
        number_of_hmbc_cosy_subgraphs = 1
    else:
        # find the number of subgraphs
        number_of_hmbc_cosy_subgraphs = nx.number_connected_components(G2)

    components = list(nx.connected_components(G2))
    number_of_hmbc_cosy_subgraphs = len(components)

    for i, comp in enumerate(sorted(components, key=len, reverse=True)):
        print(f"Sub-network {i+1} ({len(comp)} nodes): {sorted(comp)}")

    return jsonGraphData, jsonGraphData_mol, shortest_paths, number_of_hmbc_cosy_subgraphs


def _build_catoms_df(solution):
    """Slice and prepare the carbon-atoms DataFrame from molprops_df.

    Returns
    -------
    pandas.DataFrame
        A copy of the relevant molprops columns, renamed and padded with
        empty assignment columns ready for simulated annealing.
    """
    catoms_df = solution.expected_molecule.molprops_df[
        [
            "ppm",
            "atom_idx",
            "atomNumber",
            "numProtons",
            "x",
            "y",
            "sym_atom_idx",
            "sym_atomNumber",
        ]
    ].copy()

    catoms_df.rename(
        columns={"ppm": "ppm_calculated", "atom_idx": "id"},
        inplace=True,
    )

    catoms_df["ppm"] = -1000000.0
    catoms_df["jCouplingVals"] = ""
    catoms_df["jCouplingClass"] = ""
    catoms_df["visible"] = "false"
    catoms_df["symbol"] = "C"
    catoms_df["iupacLabel"] = ""

    return catoms_df


def _run_simulated_annealing(solution, json_data, jsonGraphData, catoms_df) -> tuple:
    """Run (or skip) simulated annealing and copy optimised values back.

    Simulated annealing is skipped when either the SA flag in the JSON data
    is False or the predicted weight is zero. Per Eric: this always runs
    from whatever assignment it's handed (including one seeded by
    _assign_carbons's previous_positions) and is free to change any node --
    no special pinning of previously-restored peaks against SA reassignment.

    Returns
    -------
    tuple[dict, dict]
        ``(jsonGraphData, best_results)`` where ``jsonGraphData`` now contains
        a ``"moved_nodes"`` list with optimised positions, and ``best_results``
        is the scoring dict from the SA run (or the zero-weight sentinel from
        ``process_results_SA_skipped``).
    """
    possible_symmetry = (
        solution.c13_df.shape[0] < solution.expected_molecule.num_carbon_atoms
    )

    # Same source of truth as the rest of the pipeline: sym_molprops_df was
    # already rebuilt by _reconcile_molecule_symmetry (inside
    # initialise_prior_to_carbon_assignment, called by _assign_carbons above),
    # and ch3_ch1_evidence_based was set by _assign_proton_types -- neither
    # is recomputed here.
    orig_nodes_displayed = solution.expected_molecule.sym_molprops_df["atom_idx"].tolist()
    pool_ch3_ch1 = not getattr(solution, "ch3_ch1_evidence_based", False)

    simAnneal = SimulatedAnnealing2.from_params(
        copy.deepcopy(jsonGraphData["nodes"]),
        copy.deepcopy(jsonGraphData["links"]),
        json_data["molfile"]["data"]["0"],
        json_data,
        possible_symmetry,
        orig_nodes_displayed,
        pool_ch3_ch1,
    )
    simAnneal.setup_run(cooling_rate=SIMULATED_ANNEALING_COOLING_RATE)

    if json_data["simulatedAnnealing"]["data"]["0"] and simAnneal.predicted_weight > 0:
        simAnneal.run_optimization(SIMULATED_ANNEALING_RUNS)
        jsonGraphData, best_results = simAnneal.process_results(catoms_df, jsonGraphData)
    else:
        best_results = simAnneal.process_results_SA_skipped()

    # Copy optimised node positions back into catoms_df
    for idx, row in catoms_df.iterrows():
        node_id = row["id"]
        for node in jsonGraphData["moved_nodes"]:
            if node["id"] == node_id:
                catoms_df.loc[idx, "ppm"] = node["ppm"]
                catoms_df.loc[idx, "x"] = node["x"]
                catoms_df.loc[idx, "y"] = node["y"]
                catoms_df.loc[idx, "ppm_calculated"] = node["ppm_calculated"]
                catoms_df.loc[idx, "atomNumber"] = node["atomNumber"]

    return jsonGraphData, best_results


def _build_jinja_context(
    solution,
    problemdata_json,
    json_data,
    jsonGraphData,
    jsonGraphData_mol,
    shortest_paths,
    catoms_df,
    best_results,
    number_of_hmbc_cosy_subgraphs,
) -> dict:
    """Assemble the Jinja2 template context dict from all pipeline outputs.

    Returns
    -------
    dict
        A flat dict whose keys map directly to template variables in
        ``d3molplotmnova_template.html``.
    """
    catoms_str = json.dumps(catoms_df.to_dict(orient="records"), indent=4)
    exclusion_summary = solution.build_exclusion_summary()
    json_data_str = json.dumps({**json_data, "exclusion_summary": exclusion_summary}, indent=4)
    dataFrom = (
        "nmrshiftdb2" if problemdata_json.prediction_from_nmrshiftdb2() else "mnova"
    )

    return {
        "svg_container": solution.expected_molecule.svg_str,
        "graph_edges": jsonGraphData["links"],
        "graph_nodes": jsonGraphData["moved_nodes"],
        "orig_nodes": jsonGraphData["nodes"],
        "molgraph": jsonGraphData_mol,
        "shortest_paths": json.dumps(shortest_paths),
        "catoms": catoms_str,
        "title": "dummy_title",
        "smilesString": solution.smilesstr,
        "molFile": solution.molstr,
        "workingDirectory": problemdata_json.dataframes["workingDirectory"]
        .loc[0, "workingDirectory"]
        .strip(),
        "workingFilename": problemdata_json.dataframes["workingFilename"]
        .loc[0, "workingFilename"]
        .strip(),
        "dataFrom": dataFrom,
        "oldjsondata": json_data_str,
        "best_results": best_results,
        "number_of_hmbc_cosy_subgraphs": number_of_hmbc_cosy_subgraphs,
        "exclusion_summary": exclusion_summary,
    }


# ── Public API ────────────────────────────────────────────────────────────────

def run(problemdata_json, json_data: dict) -> dict:
    """Execute the full prediction + simulated annealing pipeline.

    Parameters
    ----------
    problemdata_json:
        ``NMRProblem`` instance created from the incoming POST data.
    json_data:
        The raw JSON dict from the request (passed through to several steps).

    Returns
    -------
    dict
        Jinja2 template context dict ready for ``render_template``.

    Raises
    ------
    PipelineError
        On any recoverable failure; carries ``.message`` and ``.code``
        for the Flask response.
    """
    logger.info("prediction_pipeline: starting")

    solution = _build_solution(problemdata_json)
    logger.debug("prediction_pipeline: solution built")

    _assign_proton_types(solution)
    logger.debug("prediction_pipeline: proton types assigned")

    _transfer_hsqc_info(solution)
    logger.debug("prediction_pipeline: HSQC info transferred")

    previous_positions = _previous_atom_assignments_lookup(json_data)

    _assign_carbons(solution, previous_positions)
    logger.debug("prediction_pipeline: carbons assigned")

    jsonGraphData, jsonGraphData_mol, shortest_paths, number_of_hmbc_cosy_subgraphs = _build_graph(solution, json_data)
    logger.debug(f"prediction_pipeline: graphs built, number_of_hmbc_cosy_subgraphs={number_of_hmbc_cosy_subgraphs}")

    catoms_df = _build_catoms_df(solution)
    logger.debug("prediction_pipeline: catoms_df prepared")

    jsonGraphData, best_results = _run_simulated_annealing(
        solution, json_data, jsonGraphData, catoms_df
    )
    logger.debug(f"prediction_pipeline: SA complete, best_weight={best_results.get('best_weight')}")

    context = _build_jinja_context(
        solution,
        problemdata_json,
        json_data,
        jsonGraphData,
        jsonGraphData_mol,
        shortest_paths,
        catoms_df,
        best_results,
        number_of_hmbc_cosy_subgraphs,
    )
    logger.info("prediction_pipeline: complete")
    return context
