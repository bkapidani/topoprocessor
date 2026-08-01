"""Python interface for Topoprocessor's mesher-neutral mesh model."""

from ._core import (
    Cell,
    CellKind,
    Facet,
    Mesh,
    facet_node_count,
    node_count,
    read_gmsh,
    read_netgen,
)
from .adapters import from_gmsh, from_netgen, from_ngsolve

__all__ = [
    "Cell",
    "CellKind",
    "Facet",
    "Mesh",
    "facet_node_count",
    "from_gmsh",
    "from_netgen",
    "from_ngsolve",
    "node_count",
    "read_gmsh",
    "read_netgen",
]
