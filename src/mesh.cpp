#include "mesh.hpp"

#include <algorithm>
#include <set>
#include <string>
#include <utility>

namespace topo {
namespace {

void validate_nodes(const std::vector<NodeId>& nodes, std::size_t expected,
                    const char* entity)
{
    if (nodes.size() != expected) {
        throw std::invalid_argument(std::string(entity) + " has " +
                                    std::to_string(nodes.size()) +
                                    " nodes; expected " +
                                    std::to_string(expected));
    }

    const std::set<NodeId> unique(nodes.begin(), nodes.end());
    if (unique.size() != nodes.size()) {
        throw std::invalid_argument(std::string(entity) +
                                    " contains repeated nodes");
    }
}

void validate_range(const std::vector<NodeId>& nodes, std::size_t point_count,
                    const char* entity)
{
    for (const NodeId node : nodes) {
        if (node >= point_count) {
            throw std::out_of_range(std::string(entity) + " node " +
                                    std::to_string(node) +
                                    " is outside the point array");
        }
    }
}

} // namespace

std::size_t node_count(CellKind kind) noexcept
{
    return kind == CellKind::tetrahedron ? 4U : 8U;
}

std::size_t facet_node_count(CellKind kind) noexcept
{
    return kind == CellKind::tetrahedron ? 3U : 4U;
}

Cell::Cell(CellKind kind, Label label, std::vector<NodeId> nodes)
    : kind_(kind), label_(label), nodes_(std::move(nodes))
{
    validate_nodes(nodes_, node_count(kind_), "cell");
}

Facet::Facet(Label label, std::vector<NodeId> nodes)
    : label_(label), nodes_(std::move(nodes))
{
    if (nodes_.size() != 3U && nodes_.size() != 4U) {
        throw std::invalid_argument("facet must have three or four nodes");
    }
    validate_nodes(nodes_, nodes_.size(), "facet");
}

Mesh::Mesh(std::vector<Point> points, std::vector<Cell> cells,
           std::vector<Facet> boundary_facets)
    : points_(std::move(points)), cells_(std::move(cells)),
      boundary_facets_(std::move(boundary_facets))
{
    validate();
}

CellKind Mesh::cell_kind() const
{
    if (cells_.empty()) {
        throw std::logic_error("a mesh without cells has no cell kind");
    }
    return cells_.front().kind();
}

void Mesh::validate() const
{
    if (points_.empty()) {
        throw std::invalid_argument("mesh has no points");
    }
    if (cells_.empty()) {
        throw std::invalid_argument("mesh has no cells");
    }

    const CellKind kind = cells_.front().kind();
    for (const Cell& cell : cells_) {
        if (cell.kind() != kind) {
            throw std::invalid_argument("hybrid meshes are not supported");
        }
        validate_range(cell.nodes(), points_.size(), "cell");
    }

    const std::size_t expected_facet_nodes = facet_node_count(kind);
    std::set<std::vector<NodeId>> seen_facets;
    for (const Facet& facet : boundary_facets_) {
        if (facet.nodes().size() != expected_facet_nodes) {
            throw std::invalid_argument(
                "facet shape does not match the mesh cell kind");
        }
        validate_range(facet.nodes(), points_.size(), "facet");

        std::vector<NodeId> canonical = facet.nodes();
        std::sort(canonical.begin(), canonical.end());
        if (!seen_facets.insert(canonical).second) {
            throw std::invalid_argument("mesh contains a duplicate facet");
        }
    }
}

} // namespace topo
