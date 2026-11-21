"""vsbench: a bias-audited virtual-screening benchmark harness.

Convention used throughout the package: a higher score means "more likely
active". Docking engines report energies (lower is better), so callers convert
with :func:`vsbench.metrics.to_higher_better` before computing any metric.
"""

__version__ = "0.1.0"
