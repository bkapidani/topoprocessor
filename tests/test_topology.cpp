#include "topology.hpp"

#include <functional>
#include <sstream>
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

topo::NodeId point_id(std::size_t x, std::size_t y, std::size_t z)
{
    return static_cast<topo::NodeId>(x + 4U * y + 16U * z);
}

std::vector<topo::Point> grid_points()
{
    std::vector<topo::Point> points;
    for (std::size_t z = 0; z < 2; ++z) {
        for (std::size_t y = 0; y < 4; ++y) {
            for (std::size_t x = 0; x < 4; ++x) {
                points.push_back(topo::Point{{static_cast<double>(x),
                                              static_cast<double>(y),
                                              static_cast<double>(z)}});
            }
        }
    }
    return points;
}

std::vector<topo::NodeId> cube_nodes(std::size_t x, std::size_t y)
{
    return {point_id(x, y, 0), point_id(x + 1U, y, 0),
            point_id(x + 1U, y + 1U, 0), point_id(x, y + 1U, 0),
            point_id(x, y, 1), point_id(x + 1U, y, 1),
            point_id(x + 1U, y + 1U, 1), point_id(x, y + 1U, 1)};
}

topo::Mesh cubical_grid()
{
    std::vector<topo::Cell> cells;
    for (std::size_t y = 0; y < 3; ++y) {
        for (std::size_t x = 0; x < 3; ++x) {
            const topo::Label label = (x == 1U && y == 1U) ? 2U : 1U;
            cells.emplace_back(topo::CellKind::hexahedron, label,
                               cube_nodes(x, y));
        }
    }
    return topo::Mesh(grid_points(), std::move(cells), {});
}

topo::Mesh tetrahedral_grid()
{
    std::vector<topo::Cell> cells;
    for (std::size_t y = 0; y < 3; ++y) {
        for (std::size_t x = 0; x < 3; ++x) {
            const topo::Label label = (x == 1U && y == 1U) ? 2U : 1U;
            const std::vector<topo::NodeId> n = cube_nodes(x, y);
            const std::vector<std::vector<topo::NodeId>> tetrahedra = {
                {n[0], n[1], n[2], n[6]}, {n[0], n[2], n[3], n[6]},
                {n[0], n[3], n[7], n[6]}, {n[0], n[7], n[4], n[6]},
                {n[0], n[4], n[5], n[6]}, {n[0], n[5], n[1], n[6]}};
            for (const auto& tetrahedron : tetrahedra) {
                cells.emplace_back(topo::CellKind::tetrahedron, label,
                                   tetrahedron);
            }
        }
    }
    return topo::Mesh(grid_points(), std::move(cells), {});
}

topo::Mesh double_hole_cubical_grid()
{
    constexpr std::size_t width = 6U;
    constexpr std::size_t height = 4U;
    const auto node = [=](std::size_t x, std::size_t y, std::size_t z) {
        return static_cast<topo::NodeId>(x + width * y + width * height * z);
    };

    std::vector<topo::Point> points;
    for (std::size_t z = 0; z < 2; ++z) {
        for (std::size_t y = 0; y < height; ++y) {
            for (std::size_t x = 0; x < width; ++x) {
                points.push_back(topo::Point{{static_cast<double>(x),
                                              static_cast<double>(y),
                                              static_cast<double>(z)}});
            }
        }
    }

    std::vector<topo::Cell> cells;
    for (std::size_t y = 0; y < 3; ++y) {
        for (std::size_t x = 0; x < 5; ++x) {
            const bool hole = y == 1U && (x == 1U || x == 3U);
            const topo::Label label = hole ? 2U : 1U;
            cells.emplace_back(
                topo::CellKind::hexahedron, label,
                std::vector<topo::NodeId>{
                    node(x, y, 0), node(x + 1U, y, 0),
                    node(x + 1U, y + 1U, 0), node(x, y + 1U, 0),
                    node(x, y, 1), node(x + 1U, y, 1),
                    node(x + 1U, y + 1U, 1), node(x, y + 1U, 1)});
        }
    }
    return topo::Mesh(std::move(points), std::move(cells), {});
}

void require_ring_result(const topo::CohomologyResult& result,
                         std::size_t selected_cells)
{
    if (result.betti_number() != 1U ||
        result.selected_cell_count != selected_cells || result.edges.empty() ||
        result.generators.front().size() != result.edges.size()) {
        throw std::runtime_error("annular mesh did not produce one generator");
    }
}

} // namespace

int main()
{
    const topo::Mesh hexes = cubical_grid();
    const topo::CohomologyResult solid = topo::compute_cohomology(hexes);
    if (solid.betti_number() != 0U || solid.selected_cell_count != 9U) {
        throw std::runtime_error("solid cubical grid should have trivial H^1");
    }

    const topo::MaterialSelection materials({2U}, {1U});
    const topo::CohomologyResult hex_ring =
        topo::compute_cohomology(hexes, materials);
    require_ring_result(hex_ring, 8U);

    const topo::CohomologyResult tet_ring =
        topo::compute_cohomology(tetrahedral_grid(), materials);
    require_ring_result(tet_ring, 48U);

    const topo::CohomologyResult double_hole =
        topo::compute_cohomology(double_hole_cubical_grid(), materials);
    if (double_hole.betti_number() != 2U ||
        double_hole.selected_cell_count != 13U) {
        throw std::runtime_error(
            "two-hole cubical mesh did not produce two generators");
    }

    std::ostringstream legacy_output;
    topo::write_h1(hex_ring, legacy_output);
    if (legacy_output.str().empty() || legacy_output.str().front() != '1') {
        throw std::runtime_error("h1 exporter did not write the generator count");
    }

    expect_failure<std::invalid_argument>([&] {
        topo::compute_cohomology(hexes,
                                 topo::MaterialSelection({1U}, {1U}));
    });
    expect_failure<std::invalid_argument>([&] {
        topo::compute_cohomology(hexes,
                                 topo::MaterialSelection({}, {1U}));
    });
}
