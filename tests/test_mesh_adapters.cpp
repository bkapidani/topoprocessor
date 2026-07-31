#include "mesh_adapters.hpp"

#include <functional>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

template <typename Exception>
void expect_failure(const std::function<void()>& operation)
{
    try {
        operation();
    } catch (const Exception&) {
        return;
    }
    throw std::runtime_error("operation did not throw the expected exception");
}

std::string fixture(const std::string& root, const std::string& relative)
{
    return root + "/" + relative;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 2) {
        throw std::runtime_error("fixture root argument is required");
    }
    const std::string root = argv[1];

    const topo::Mesh tetra =
        topo::read_netgen(fixture(root, "netgen/furch.mesh"));
    if (tetra.cell_kind() != topo::CellKind::tetrahedron ||
        tetra.points().size() != 1132U || tetra.cells().size() != 5691U ||
        tetra.boundary_facets().size() != 846U) {
        throw std::runtime_error("Netgen tetrahedral fixture changed in transit");
    }

    const topo::Mesh hex =
        topo::read_netgen(fixture(root, "netgen/unit_hex.mesh"));
    if (hex.cell_kind() != topo::CellKind::hexahedron ||
        hex.points().size() != 8U || hex.cells().size() != 1U ||
        hex.boundary_facets().size() != 6U ||
        hex.cells().front().nodes().front() != 0U) {
        throw std::runtime_error("Netgen hexahedral fixture changed in transit");
    }

    const topo::Mesh gmsh_tetra =
        topo::read_gmsh(fixture(root, "gmsh/unit_tetra.msh"));
    if (gmsh_tetra.cell_kind() != topo::CellKind::tetrahedron ||
        gmsh_tetra.cells().front().label() != 7U ||
        gmsh_tetra.cells().front().nodes() !=
            std::vector<topo::NodeId>({0U, 1U, 2U, 3U}) ||
        gmsh_tetra.boundary_facets().front().label() != 11U) {
        throw std::runtime_error("Gmsh tetrahedral metadata was not preserved");
    }

    const topo::Mesh gmsh_hex =
        topo::read_gmsh(fixture(root, "gmsh/unit_hex.msh"));
    if (gmsh_hex.cell_kind() != topo::CellKind::hexahedron ||
        gmsh_hex.points().size() != 8U || gmsh_hex.cells().size() != 1U ||
        gmsh_hex.boundary_facets().size() != 6U ||
        gmsh_hex.cells().front().label() != 3U) {
        throw std::runtime_error("Gmsh hexahedral metadata was not preserved");
    }

    expect_failure<std::invalid_argument>([] {
        std::istringstream invalid(
            "4\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n"
            "1\n1 0 2 3 4\n0\n");
        topo::read_netgen(invalid);
    });

    expect_failure<std::invalid_argument>([] {
        std::istringstream invalid(
            "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"
            "$Nodes\n4\n1 0 0 0\n2 1 0 0\n3 0 1 0\n4 0 0 1\n$EndNodes\n"
            "$Elements\n1\n1 4 1 7 1 2 3 99\n$EndElements\n");
        topo::read_gmsh(invalid);
    });

    expect_failure<std::invalid_argument>([] {
        std::istringstream unsupported(
            "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"
            "$Nodes\n6\n1 0 0 0\n2 1 0 0\n3 0 1 0\n"
            "4 0 0 1\n5 1 0 1\n6 0 1 1\n$EndNodes\n"
            "$Elements\n1\n1 6 1 7 1 2 3 4 5 6\n$EndElements\n");
        topo::read_gmsh(unsupported);
    });
}
