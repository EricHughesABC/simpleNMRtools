"""Sync pipeline for simpleNMR.

Rebuilds the Jinja2 template context from previously computed NMR
assignments without re-running the solver.  Invoked when the incoming
JSON data represents a re-render / layout-sync rather than a fresh
prediction request (i.e. ``NMRProblem.is_prediction()`` is False).

Public API
----------
``run(problemdata_json, json_data) -> dict``
    Prepare the network graph and return a Jinja2 context dict.
    Does not raise ``PipelineError``; failures propagate as-is.
"""

import json

from loguru import logger


# ── Public API ────────────────────────────────────────────────────────────────

def run(problemdata_json, json_data: dict) -> dict:
    """Rebuild the Jinja2 context from previously computed assignments.

    Parameters
    ----------
    problemdata_json:
        ``NMRProblem`` instance created from the incoming POST data.
    json_data:
        The raw JSON dict from the request.

    Returns
    -------
    dict
        Jinja2 template context dict ready for ``render_template``.
    """
    logger.info("sync_pipeline: starting")

    problemdata_json.prepare_network_graph()

    json_data_str = json.dumps(json_data, indent=4)
    dataFrom = (
        "nmrshiftdb2" if problemdata_json.prediction_from_nmrshiftdb2() else "mnova"
    )

    jinja_template = problemdata_json.jinjadata
    jinja_template["smilesString"] = problemdata_json.dataframes["smiles"].loc[
        0, "smiles"
    ]
    jinja_template["molFile"] = problemdata_json.dataframes["molfile"].loc[
        0, "molfile"
    ]
    jinja_template["workingDirectory"] = (
        problemdata_json.dataframes["workingDirectory"]
        .loc[0, "workingDirectory"]
        .strip()
    )
    jinja_template["workingFilename"] = (
        problemdata_json.dataframes["workingFilename"]
        .loc[0, "workingFilename"]
        .strip()
    )
    jinja_template["oldjsondata"] = json_data_str
    jinja_template["dataFrom"] = dataFrom
    jinja_template["best_results"] = {
        "best_weight": 0,
        "best_mae": 0.0,
        "best_lae": 0.0,
    }

    logger.info("sync_pipeline: complete")
    return jinja_template
