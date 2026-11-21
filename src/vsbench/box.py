"""Docking search box derived from a reference ligand."""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

# AutoGrid / AutoDock-GPU accept at most 256 grid points per axis.
MAX_GRID_POINTS = 256


@dataclass(frozen=True)
class Box:
    center: tuple[float, float, float]
    size: tuple[float, float, float]
    spacing: float = 0.375

    @property
    def npts(self) -> tuple[int, int, int]:
        """Grid points per axis, rounded up to the even number AutoGrid requires."""
        out = []
        for edge in self.size:
            n = math.ceil(edge / self.spacing)
            out.append(n + (n % 2))
        return tuple(out)  # type: ignore[return-value]


def box_from_coords(coords, buffer: float = 6.0, min_edge: float = 24.0,
                    spacing: float = 0.375, cubic: bool = True) -> Box:
    """Box centred on the heavy-atom centroid of a reference ligand.

    The edge is the largest ligand span plus a buffer on each side (cubic) or
    the span along each axis plus the buffer (non-cubic), never below
    ``min_edge``. Raises if the box would exceed the grid-point limit.

    Defaults are deliberately generous (6 A buffer, 24 A minimum). In a pilot
    with a 4 A buffer, a large screening compound received a +61 kcal/mol
    gnina energy from steric clashes: a box sized to the template ligand
    penalises larger ligands, which is a size bias, not binding signal.
    """
    xyz = np.asarray(coords, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or len(xyz) < 3:
        raise ValueError("need at least three 3-D coordinates")
    center = tuple(float(v) for v in xyz.mean(axis=0))
    span = xyz.max(axis=0) - xyz.min(axis=0)
    if cubic:
        edge = max(float(span.max()) + 2 * buffer, min_edge)
        size = (edge, edge, edge)
    else:
        size = tuple(max(float(s) + 2 * buffer, min_edge) for s in span)
    box = Box(center=center, size=size, spacing=spacing)
    if max(box.npts) > MAX_GRID_POINTS:
        raise ValueError(
            f"box {size} at spacing {spacing} needs {box.npts} grid points; "
            f"limit is {MAX_GRID_POINTS} per axis"
        )
    return box
