#ifndef TOPOPROCESSOR_MESH_HPP
#define TOPOPROCESSOR_MESH_HPP

#include <array>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace topo {

using NodeId = std::uint32_t;
using Label = std::uint32_t;
using Point = std::array<double, 3>;

enum class CellKind {
    tetrahedron,
    hexahedron,
};

// Hexahedra use the cyclic Gmsh/VTK convention: nodes 0..3 form one
// quadrilateral face and nodes 4..7 form the corresponding opposite face.

class Cell {
public:
    Cell(CellKind kind, Label label, std::vector<NodeId> nodes);

    CellKind kind() const noexcept { return kind_; }
    Label label() const noexcept { return label_; }
    const std::vector<NodeId>& nodes() const noexcept { return nodes_; }

private:
    CellKind kind_;
    Label label_;
    std::vector<NodeId> nodes_;
};

class Facet {
public:
    Facet(Label label, std::vector<NodeId> nodes);

    Label label() const noexcept { return label_; }
    const std::vector<NodeId>& nodes() const noexcept { return nodes_; }

private:
    Label label_;
    std::vector<NodeId> nodes_;
};

class Mesh {
public:
    Mesh(std::vector<Point> points, std::vector<Cell> cells,
         std::vector<Facet> boundary_facets);

    const std::vector<Point>& points() const noexcept { return points_; }
    const std::vector<Cell>& cells() const noexcept { return cells_; }
    const std::vector<Facet>& boundary_facets() const noexcept
    {
        return boundary_facets_;
    }

    CellKind cell_kind() const;
    void validate() const;

private:
    std::vector<Point> points_;
    std::vector<Cell> cells_;
    std::vector<Facet> boundary_facets_;
};

std::size_t node_count(CellKind kind) noexcept;
std::size_t facet_node_count(CellKind kind) noexcept;

} // namespace topo

#endif
