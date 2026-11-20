"""Candidates pipeline: run the prediction pipeline on several structures.

Input (POST body of /simpleMNOVAcandidates), produced by simplecasenmr's
simplenmr_export.candidates_request():

    {
      "simplecasenmr": {"confidence": ..., "stage": ..., "reader_warnings": [...], ...},
      "candidates": [
        {"rank": 1, "smiles": "...", "combined": 0.03, "shift_mae": 0.11,
         "corr_penalty": 0.0, "n_long_correlations": 0,
         "input": { ...a complete simpleNMR input JSON for this structure... }},
        ...
      ]
    }

Each candidate's "input" goes through the normal prediction pipeline
(assignment + simulated annealing). A failure in one candidate does not stop
the others; it is reported in that candidate's tab.
"""

from __future__ import annotations

from loguru import logger

from app.pipelines import PipelineError, prediction_pipeline
from core.html_from_assignments import NMRProblem


def run(request_json: dict, max_candidates: int = 5) -> list[dict]:
    """Return one entry per candidate:
    {"meta": {...scores from simplecasenmr...}, "context": jinja context or None,
     "error": message or None}"""
    cands = request_json.get("candidates", [])[:max_candidates]
    if not cands:
        raise PipelineError("No candidates in request", 400)
    out = []
    for c in cands:
        meta = {k: v for k, v in c.items() if k != "input"}
        try:
            problem = NMRProblem.from_mnova_dict(c["input"])
            ctx = prediction_pipeline.run(problem, c["input"])
            out.append({"meta": meta, "context": ctx, "error": None})
        except PipelineError as e:
            logger.warning(f"candidate {meta.get('rank')}: {e.message}")
            out.append({"meta": meta, "context": None, "error": e.message})
        except Exception as e:  # keep the other candidates going
            logger.exception(f"candidate {meta.get('rank')} failed")
            out.append({"meta": meta, "context": None, "error": str(e)})
    return out
