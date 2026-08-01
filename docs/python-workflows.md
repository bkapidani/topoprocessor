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

The adapters normalize Netgen's one-based point identifiers. Volume material
indices become cell labels and surface indices become facet labels. Only
linear tetrahedral and hexahedral volume meshes are accepted.

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
classification ambiguous. Supported elements are linear triangles,
quadrangles, tetrahedra, and hexahedra.

## Legacy files

`topoprocessor.read_netgen(path)` and `topoprocessor.read_gmsh(path)` remain
available for the recovered neutral and Gmsh 2.2 formats. New applications
should prefer the in-memory adapters so mesh parsing remains the mesher's
responsibility.
