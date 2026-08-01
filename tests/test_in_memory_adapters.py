import importlib
import sys
import types
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class FakeCellKind:
    tetrahedron = "tetrahedron"
    hexahedron = "hexahedron"


class FakeCell:
    def __init__(self, kind, label, nodes):
        self.kind = kind
        self.label = label
        self.nodes = list(nodes)


class FakeFacet:
    def __init__(self, label, nodes):
        self.label = label
        self.nodes = list(nodes)


class FakeMesh:
    def __init__(self, points, cells, boundary_facets):
        self.points = [tuple(point) for point in points]
        self.cells = list(cells)
        self.boundary_facets = list(boundary_facets)
        kinds = {cell.kind for cell in self.cells}
        if len(kinds) > 1:
            raise ValueError("hybrid meshes are not supported")


def install_fake_core():
    core = types.ModuleType("topoprocessor._core")
    core.Cell = FakeCell
    core.CellKind = FakeCellKind
    core.Facet = FakeFacet
    core.Mesh = FakeMesh
    core.facet_node_count = (
        lambda kind: 3 if kind == FakeCellKind.tetrahedron else 4
    )
    core.node_count = lambda kind: 4 if kind == FakeCellKind.tetrahedron else 8
    core.read_gmsh = lambda filename: filename
    core.read_netgen = lambda filename: filename
    sys.modules[core.__name__] = core


try:
    import topoprocessor as topo
except ImportError:
    sys.modules.pop("topoprocessor", None)
    install_fake_core()
    sys.path.insert(0, str(ROOT / "python"))
    import topoprocessor as topo

adapters = importlib.import_module("topoprocessor.adapters")
CellKind = topo.CellKind


class PointId:
    def __init__(self, number):
        self.nr = number


class NetgenPoint:
    def __init__(self, coordinates):
        self.p = coordinates


class NetgenElement:
    def __init__(self, index, vertices):
        self.index = index
        self.vertices = [PointId(vertex) for vertex in vertices]


class NetgenMesh:
    def __init__(self):
        self._points = [
            NetgenPoint((0, 0, 0)),
            NetgenPoint((1, 0, 0)),
            NetgenPoint((0, 1, 0)),
            NetgenPoint((0, 0, 1)),
        ]
        self._cells = [NetgenElement(7, [1, 2, 3, 4])]
        self._facets = [
            NetgenElement(11, [1, 2, 3]),
            NetgenElement(12, [1, 2, 4]),
            NetgenElement(13, [1, 3, 4]),
            NetgenElement(14, [2, 3, 4]),
        ]

    def Points(self):
        return self._points

    def Elements3D(self):
        return self._cells

    def Elements2D(self):
        return self._facets


class NGSolveMesh:
    def __init__(self, ngmesh):
        self.ngmesh = ngmesh


class GmshMeshApi:
    def __init__(self, volume_type=4):
        self.volume_type = volume_type

    def getNodes(self):
        return (
            [10, 20, 30, 40],
            [0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1],
            [],
        )

    def getElements(self, dimension, entity_tag):
        if dimension == 3:
            return [self.volume_type], [[100]], [[10, 20, 30, 40]]
        return [2], [[201, 202, 203, 204]], [
            [10, 20, 30, 10, 20, 40, 10, 30, 40, 20, 30, 40]
        ]


class GmshModel:
    def __init__(self, volume_type=4, duplicate_groups=False):
        self.mesh = GmshMeshApi(volume_type)
        self.duplicate_groups = duplicate_groups

    def getEntities(self, dimension):
        return [(dimension, 70 if dimension == 3 else 110)]

    def getPhysicalGroupsForEntity(self, dimension, entity_tag):
        if self.duplicate_groups and dimension == 3:
            return [7, 8]
        return [7 if dimension == 3 else 11]


class InMemoryAdaptersTest(unittest.TestCase):
    def test_netgen_mesh_is_normalized(self):
        mesh = adapters.from_netgen(NetgenMesh())

        self.assertEqual(list(mesh.points[1]), [1.0, 0.0, 0.0])
        self.assertEqual(mesh.cells[0].kind, CellKind.tetrahedron)
        self.assertEqual(mesh.cells[0].label, 7)
        self.assertEqual(mesh.cells[0].nodes, [0, 1, 2, 3])
        self.assertEqual(mesh.boundary_facets[0].label, 11)

    def test_ngsolve_mesh_unwraps_netgen_mesh(self):
        mesh = adapters.from_ngsolve(NGSolveMesh(NetgenMesh()))

        self.assertEqual(len(mesh.cells), 1)
        self.assertEqual(len(mesh.boundary_facets), 4)

    def test_gmsh_model_uses_physical_tags(self):
        mesh = adapters.from_gmsh(GmshModel())

        self.assertEqual(mesh.cells[0].label, 7)
        self.assertEqual(mesh.cells[0].nodes, [0, 1, 2, 3])
        self.assertEqual(len(mesh.boundary_facets), 4)
        self.assertEqual(mesh.boundary_facets[0].label, 11)

    def test_initialized_gmsh_module_is_accepted(self):
        module = types.SimpleNamespace(model=GmshModel())

        self.assertEqual(adapters.from_gmsh(module).cells[0].label, 7)

    def test_ambiguous_gmsh_physical_groups_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "multiple physical groups"):
            adapters.from_gmsh(GmshModel(duplicate_groups=True))

    def test_higher_order_gmsh_cells_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "element type 11"):
            adapters.from_gmsh(GmshModel(volume_type=11))


if __name__ == "__main__":
    unittest.main()
