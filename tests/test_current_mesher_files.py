import importlib.util
import tempfile
import unittest
from pathlib import Path

import topoprocessor as topo


class CurrentMesherFilesTest(unittest.TestCase):
    @unittest.skipUnless(importlib.util.find_spec("gmsh"), "gmsh is not installed")
    def test_current_gmsh_msh41_and_higher_order_tetrahedra(self):
        import gmsh

        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "second_order.msh"
            gmsh.initialize()
            try:
                gmsh.option.setNumber("General.Verbosity", 0)
                gmsh.model.add("second_order_tetrahedra")
                volume = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
                gmsh.model.occ.synchronize()
                gmsh.model.addPhysicalGroup(3, [volume], 7)
                gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
                gmsh.option.setNumber("Mesh.ElementOrder", 2)
                gmsh.model.mesh.generate(3)
                gmsh.write(str(filename))
            finally:
                gmsh.finalize()

            mesh = topo.load_gmsh(filename)

        self.assertEqual(mesh.cell_kind, topo.CellKind.tetrahedron)
        self.assertTrue(mesh.cells)
        self.assertEqual({cell.label for cell in mesh.cells}, {7})
        self.assertTrue(all(len(cell.nodes) == 4 for cell in mesh.cells))
        self.assertEqual(topo.compute_cohomology(mesh).betti_number, 0)

    @unittest.skipUnless(
        importlib.util.find_spec("netgen"), "netgen-mesher is not installed"
    )
    def test_current_netgen_vol(self):
        from netgen.csg import unit_cube

        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "unit_cube.vol"
            generated = unit_cube.GenerateMesh(maxh=0.7)
            generated.Save(str(filename))

            mesh = topo.load_netgen(filename)

        self.assertEqual(mesh.cell_kind, topo.CellKind.tetrahedron)
        self.assertTrue(mesh.cells)
        self.assertTrue(all(len(cell.nodes) == 4 for cell in mesh.cells))
        self.assertEqual(topo.compute_cohomology(mesh).betti_number, 0)


if __name__ == "__main__":
    unittest.main()
