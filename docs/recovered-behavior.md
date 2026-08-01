# Recovered legacy behavior

This document freezes the evidence used before modernizing Topoprocessor. It is
intentionally descriptive: later layers must preserve these contracts unless a
change is explicitly documented.

## Baselines

* **Tetrahedra:** commit `14e4febb^` (the parent of the fixture-removal commit)
  is the last repository state containing both the tetrahedral implementation
  and its `furch.mesh` and `torus.mesh` inputs. Its cell and face types have four
  and three vertices respectively.
* **Hexahedra:** commit `0f63e379` introduces the eight-vertex volume and
  four-vertex surface implementation. Commit `20e0b8d` is the final source
  change and is therefore the recovered hexahedral baseline. The boundary
  reader's lingering triangular call was inconsistent with those types; it is
  restored to the quadrilateral parser in this recovery change.

The two implementations were never mesher-neutral variants of one core. The
modernization must first represent cell type explicitly rather than silently
interpreting every Netgen volume as a hexahedron.

## Netgen neutral-file contract

The recovered reader consumes three counted, whitespace-delimited sections:

1. points (`x y z`),
2. labelled volumes (`material node...`), and
3. labelled surface facets (`boundary node...`), including exterior boundaries
   and material interfaces.

Node identifiers are one-based in the file and are converted to zero-based
indices internally. Tetrahedra have four nodes and triangular facets;
hexahedra have eight nodes and quadrilateral facets. Hybrid files are not part
of the recovered contract.

## Frozen fixtures and invariants

`furch.mesh` and `torus.mesh` are byte-for-byte restorations from
`14e4febb^`. Their section counts and truncated SHA-256 digests are regression
oracles. `unit_hex.mesh` is the smallest closed hexahedral complex.

The automated recovery test checks:

* section counts and complete parsing;
* fixture digests for accidental edits;
* valid, non-repeated, in-range node identifiers;
* uniqueness of boundary facets;
* manifold facet incidence (one boundary or two incident volumes); and
* every declared surface is a facet of at least one volume.

## Modernization sequence

Follow-on changes should remain reviewable and proceed in this dependency
order: mesher-neutral C++ core; Python bindings; Netgen and Gmsh adapters; OCC
workflows; user documentation; CI; and packaging. Each layer must run these
recovery tests so adapters cannot change the established topology unnoticed.

The first layer is now represented by `mesh.hpp`: it stores zero-based points,
labelled cells, and labelled boundary facets without depending on a file format.
Cell kind is explicit, tetrahedral and hexahedral arities are checked, and the
current no-hybrid contract is enforced before algorithms consume the mesh.

The second layer is provided by the optional pybind11 module in
`python_bindings.cpp`. It exposes the same read-only mesh entities, explicit
cell kinds, arity helpers, and validation behavior to Python without adding a
file-format dependency. Build it with `-DTOPOPROCESSOR_BUILD_PYTHON=ON`; C++
validation failures are translated to the corresponding Python exceptions.

The third layer's legacy-file compatibility is implemented by
`mesh_adapters.hpp`. The Netgen adapter reads
the recovered counted neutral format and infers tetrahedral or hexahedral shape
from its volume records. The Gmsh adapter reads ASCII legacy 2.2 files, retains
the first (physical) element tag as the label, and supports linear triangles,
quadrilaterals, tetrahedra, and hexahedra. Both adapters normalize external node
identifiers to zero-based mesh indices and pass their result through the common
mesh validation. Their filename-based entry points are also available from the
Python module as `read_netgen` and `read_gmsh`.

The preferred adapters are Python-side `from_netgen`, `from_ngsolve`, and
`from_gmsh` functions. They consume live mesher objects through their public
APIs, normalize node identifiers, and construct the same C++ `Mesh`. OCC stays
a geometry source: `netgen.occ` and `gmsh.model.occ` generate meshes that flow
through these existing adapters without adding an OCC dependency to the core.
