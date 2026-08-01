# In-memory Python mesh workflows

Topoprocessor's topology core consumes a canonical, mesher-neutral mesh. The
Python adapters translate live Netgen/NGSolve and Gmsh models into that mesh;
neither mesher nor a CAD kernel is linked into the C++ library.

## Netgen and NGSolve

Both a `netgen.meshing.Mesh` and its `ngsolve.Mesh` view are supported:

```python
from netgen.occ import Box, OCCGeometry
from topoprocessor.adapters import from_netgen, from_ngsolve

shape = Box((0, 0, 0), (1, 1, 1))
ngmesh = OCCGeometry(shape).GenerateMesh(maxh=0.25)

complex_from_netgen = from_netgen(ngmesh)

# If the application already uses NGSolve:
import ngsolve

mesh = ngsolve.Mesh(ngmesh)
complex_from_ngsolve = from_ngsolve(mesh)
```

The adapters normalize Netgen's one-based point identifiers and its
hexahedral local numbering. Volume material indices become cell labels and
surface indices become facet labels. Tetrahedral and hexahedral volume meshes
are accepted. Curved/high-order Netgen meshes expose their corner vertices
through this same low-order topology view.

## Gmsh and OpenCASCADE

The Gmsh adapter reads the active model through the official Python API, so it
works for both the built-in geometry kernel and `gmsh.model.occ`:

```python
import gmsh
from topoprocessor.adapters import from_gmsh

gmsh.initialize()
try:
    volume = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
    gmsh.model.occ.synchronize()
    gmsh.model.addPhysicalGroup(3, [volume], 7)
    gmsh.model.mesh.generate(3)

    complex = from_gmsh(gmsh.model)
finally:
    gmsh.finalize()
```

The physical-group tag of each surface or volume entity becomes its canonical
label. Ungrouped entities receive label zero. An entity in multiple physical
groups is rejected because selecting one tag implicitly would make material
classification ambiguous. Supported shapes are triangles, quadrangles,
tetrahedra, and hexahedra. Higher-order Gmsh elements use the primary nodes
reported by `getElementProperties`, identifying the embedded low-order
cochain complex. Pass `high_order="reject"` to `from_gmsh` for strict
first-order input. Prisms, pyramids, and hybrid meshes are rejected.

## Computing cohomology

The returned mesh can be passed directly to the C++ topology core:

```python
import topoprocessor as topo

materials = topo.MaterialSelection(
    conductor_labels=[2],
    insulator_labels=[1],
)
result = topo.compute_cohomology(complex_from_netgen, materials)

print(result.betti_number)
for edge, coefficient in zip(result.edges, result.generators[0]):
    if coefficient:
        print(edge.start, edge.end, coefficient)

# Optional compatibility output for existing solvers:
topo.write_h1(result, "h1.txt")
```

If only conductor labels are supplied, every other material is treated as an
insulator, matching the historical executable. Supplying insulator labels
enables strict classification: every mesh cell label must occur in exactly one
of the two lists.

The result contains integral edge cochains. The core checks the cocycle
condition before returning. A complex whose low-order incidence system cannot
be reduced integrally is rejected rather than silently rounded.

## Current and legacy files

Use the mesher-backed loaders for current native files:

```python
import topoprocessor as topo

gmsh_mesh = topo.load_gmsh("model.msh")  # current MSH 4.1 and older
netgen_mesh = topo.load_netgen("model.vol")
```

These functions require the optional `gmsh` and `netgen-mesher` packages,
respectively. They delegate parsing to the installed mesher release and then
use the same in-memory adapters.

`topoprocessor.read_netgen(path)` and `topoprocessor.read_gmsh(path)` remain
available for the recovered neutral and Gmsh 2.2 formats. New applications
should prefer `load_netgen`, `load_gmsh`, or the in-memory adapters so mesh
parsing remains the mesher's responsibility.
