"""NMR processing pipelines.

Two pipelines are available:

* ``prediction_pipeline`` — full NMR assignment + simulated annealing.
  Invoked when ``NMRProblem.is_prediction()`` is True.

* ``sync_pipeline`` — rebuilds the Jinja2 context from previously computed
  assignments without re-running the solver.
  Invoked when ``NMRProblem.is_prediction()`` is False.

Both pipelines expose a single ``run()`` function that returns a Jinja2
template context dict or raises ``PipelineError`` on failure.
"""


class PipelineError(Exception):
    """Raised by a pipeline when processing cannot continue.

    Parameters
    ----------
    message:
        Human-readable error description (returned verbatim to the client).
    code:
        HTTP status code to use in the Flask response (default: 400).
    """

    def __init__(self, message: str, code: int = 400) -> None:
        super().__init__(message)
        self.message = message
        self.code = code
