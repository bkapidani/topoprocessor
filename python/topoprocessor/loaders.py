"""File loaders backed by the current mesher Python libraries."""

import os

from .adapters import from_gmsh, from_netgen


def load_gmsh(filename, high_order="primary"):
    """Load a mesh through the installed Gmsh library.

    This is the supported route for modern MSH 4.1 files. Gmsh owns file
    parsing; Topoprocessor extracts the primary-node cochain complex.
    """

    try:
        import gmsh
    except ImportError as error:
        raise ImportError(
            "load_gmsh requires the optional 'gmsh' Python package"
        ) from error

    is_initialized = getattr(gmsh, "isInitialized", None)
    owns_session = not bool(is_initialized()) if callable(is_initialized) else True
    if owns_session:
        gmsh.initialize()
    try:
        gmsh.open(os.fspath(filename))
        return from_gmsh(gmsh.model, high_order=high_order)
    finally:
        if owns_session:
            gmsh.finalize()


def load_netgen(filename):
    """Load a Netgen ``.vol`` mesh through the installed Netgen library."""

    try:
        from netgen.meshing import Mesh as NetgenMesh
    except ImportError as error:
        raise ImportError(
            "load_netgen requires the optional 'netgen-mesher' Python package"
        ) from error

    try:
        mesh = NetgenMesh(dim=3)
    except TypeError:
        mesh = NetgenMesh()
    load = getattr(mesh, "Load", None)
    if not callable(load):
        raise TypeError("the installed Netgen Mesh object has no Load method")
    load(os.fspath(filename))
    return from_netgen(mesh)
