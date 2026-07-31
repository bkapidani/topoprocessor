#include "mesh.hpp"

#include <functional>
#include <stdexcept>
#include <vector>

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

std::vector<topo::Point> points(std::size_t count)
{
    return std::vector<topo::Point>(count, topo::Point{{0.0, 0.0, 0.0}});
}

} // namespace

int main()
{
    using topo::Cell;
    using topo::CellKind;
    using topo::Facet;
    using topo::Mesh;

    const Mesh tetra(points(4), {Cell(CellKind::tetrahedron, 7, {0, 1, 2, 3})},
                     {Facet(1, {0, 1, 2})});
    if (tetra.cell_kind() != CellKind::tetrahedron ||
        tetra.cells().front().label() != 7) {
        throw std::runtime_error("tetrahedron metadata was not preserved");
    }

    const Mesh hex(points(8),
                   {Cell(CellKind::hexahedron, 2, {0, 1, 2, 3, 4, 5, 6, 7})},
                   {Facet(3, {0, 1, 2, 3})});
    if (hex.cell_kind() != CellKind::hexahedron) {
        throw std::runtime_error("hexahedron kind was not preserved");
    }

    expect_failure<std::invalid_argument>([] {
        Cell(CellKind::tetrahedron, 1, {0, 1, 2, 2});
    });
    expect_failure<std::out_of_range>([] {
        Mesh(points(4), {Cell(CellKind::tetrahedron, 1, {0, 1, 2, 4})}, {});
    });
    expect_failure<std::invalid_argument>([] {
        Mesh(points(8),
             {Cell(CellKind::tetrahedron, 1, {0, 1, 2, 3}),
              Cell(CellKind::hexahedron, 1, {0, 1, 2, 3, 4, 5, 6, 7})},
             {});
    });
    expect_failure<std::invalid_argument>([] {
        Mesh(points(4), {Cell(CellKind::tetrahedron, 1, {0, 1, 2, 3})},
             {Facet(1, {0, 1, 2}), Facet(2, {2, 1, 0})});
    });
    expect_failure<std::invalid_argument>([] {
        Mesh(points(8),
             {Cell(CellKind::hexahedron, 1, {0, 1, 2, 3, 4, 5, 6, 7})},
             {Facet(1, {0, 1, 2})});
    });
}
