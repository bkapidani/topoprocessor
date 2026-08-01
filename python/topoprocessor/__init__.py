"""Python interface for Topoprocessor's mesher-neutral mesh model."""

from ._core import (
    Cell,
    CellKind,
    CohomologyOptions,
    CohomologyResult,
    Facet,
    MaterialSelection,
    Mesh,
    OrientedEdge,
    compute_cohomology,
    facet_node_count,
    node_count,
    read_gmsh,
    read_netgen,
    write_h1,
)
from .adapters import from_gmsh, from_netgen, from_ngsolve
from .loaders import load_gmsh, load_netgen

__all__ = [
    "Cell",
    "CellKind",
    "CohomologyOptions",
    "CohomologyResult",
    "Facet",
    "MaterialSelection",
    "Mesh",
    "OrientedEdge",
    "compute_cohomology",
    "facet_node_count",
    "from_gmsh",
    "from_netgen",
    "from_ngsolve",
    "load_gmsh",
    "load_netgen",
    "node_count",
    "read_gmsh",
    "read_netgen",
    "write_h1",
]
