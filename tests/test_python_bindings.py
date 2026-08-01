import unittest
from pathlib import Path

import topoprocessor as topo


FIXTURES = Path(__file__).parent / "fixtures"


class PythonBindingsTest(unittest.TestCase):
    def test_tetrahedral_mesh_round_trip(self):
        points = [
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        ]
        cell = topo.Cell(topo.CellKind.tetrahedron, 7, [0, 1, 2, 3])
        facets = [
            topo.Facet(11, [0, 1, 2]),
            topo.Facet(12, [0, 1, 3]),
            topo.Facet(13, [0, 2, 3]),
            topo.Facet(14, [1, 2, 3]),
        ]

        mesh = topo.Mesh(points, [cell], facets)

        self.assertEqual(mesh.cell_kind, topo.CellKind.tetrahedron)
        self.assertEqual(mesh.points, [list(point) for point in points])
        self.assertEqual(mesh.cells[0].label, 7)
        self.assertEqual(mesh.cells[0].nodes, [0, 1, 2, 3])
        self.assertEqual(mesh.boundary_facets[0].label, 11)
        self.assertIsNone(mesh.validate())

    def test_hexahedral_mesh_and_arity_helpers(self):
        points = [
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (1.0, 1.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
            (1.0, 0.0, 1.0),
            (1.0, 1.0, 1.0),
            (0.0, 1.0, 1.0),
        ]
        cell = topo.Cell(topo.CellKind.hexahedron, 3, list(range(8)))
        facet = topo.Facet(4, [0, 1, 2, 3])

        mesh = topo.Mesh(points, [cell], [facet])

        self.assertEqual(mesh.cell_kind, topo.CellKind.hexahedron)
        self.assertEqual(topo.node_count(mesh.cell_kind), 8)
        self.assertEqual(topo.facet_node_count(mesh.cell_kind), 4)

    def test_cpp_validation_errors_become_python_exceptions(self):
        with self.assertRaisesRegex(ValueError, "cell has 3 nodes; expected 4"):
            topo.Cell(topo.CellKind.tetrahedron, 1, [0, 1, 2])

        points = [(0.0, 0.0, 0.0)] * 4
        valid_cell = topo.Cell(topo.CellKind.tetrahedron, 1, [0, 1, 2, 3])
        with self.assertRaisesRegex(IndexError, "outside the point array"):
            topo.Mesh(points, [valid_cell], [topo.Facet(2, [0, 1, 4])])

        hex_cell = topo.Cell(topo.CellKind.hexahedron, 2, list(range(8)))
        with self.assertRaisesRegex(ValueError, "hybrid meshes"):
            topo.Mesh(points * 2, [valid_cell, hex_cell], [])

    def test_exposed_state_is_read_only(self):
        cell = topo.Cell(topo.CellKind.tetrahedron, 7, [0, 1, 2, 3])

        with self.assertRaises(AttributeError):
            cell.label = 8
        with self.assertRaises(AttributeError):
            cell.nodes = [3, 2, 1, 0]

    def test_netgen_adapter(self):
        mesh = topo.read_netgen(str(FIXTURES / "netgen" / "unit_hex.mesh"))

        self.assertEqual(mesh.cell_kind, topo.CellKind.hexahedron)
        self.assertEqual(len(mesh.points), 8)
        self.assertEqual(mesh.cells[0].nodes, [0, 1, 3, 2, 4, 5, 7, 6])
        self.assertEqual(len(mesh.boundary_facets), 6)

    def test_gmsh_adapter(self):
        mesh = topo.read_gmsh(str(FIXTURES / "gmsh" / "unit_tetra.msh"))

        self.assertEqual(mesh.cell_kind, topo.CellKind.tetrahedron)
        self.assertEqual(mesh.cells[0].label, 7)
        self.assertEqual(mesh.cells[0].nodes, [0, 1, 2, 3])
        self.assertEqual(mesh.boundary_facets[0].label, 11)

    def test_cohomology_api_for_contractible_mesh(self):
        points = [
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        ]
        mesh = topo.Mesh(
            points,
            [topo.Cell(topo.CellKind.tetrahedron, 1, [0, 1, 2, 3])],
            [],
        )

        result = topo.compute_cohomology(mesh)

        self.assertEqual(result.betti_number, 0)
        self.assertEqual(result.selected_cell_count, 1)
        self.assertEqual(len(result.edges), 6)


if __name__ == "__main__":
    unittest.main()
