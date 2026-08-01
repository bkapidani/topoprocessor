#include "mesh_adapters.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace topo {
namespace {

std::string trim(std::string value)
{
    const auto not_space = [](unsigned char character) {
        return !std::isspace(character);
    };
    value.erase(value.begin(),
                std::find_if(value.begin(), value.end(), not_space));
    value.erase(std::find_if(value.rbegin(), value.rend(), not_space).base(),
                value.end());
    return value;
}

std::string required_line(std::istream& input, const std::string& context)
{
    std::string line;
    while (std::getline(input, line)) {
        line = trim(line);
        if (!line.empty()) {
            return line;
        }
    }
    throw std::invalid_argument("unexpected end of file while reading " +
                                context);
}

std::vector<std::string> tokens(const std::string& line)
{
    std::istringstream stream(line);
    std::vector<std::string> result;
    std::string token;
    while (stream >> token) {
        result.push_back(token);
    }
    return result;
}

std::uint64_t unsigned_value(const std::string& token,
                             const std::string& context)
{
    if (token.empty() || token.front() == '-') {
        throw std::invalid_argument(context + " must be a non-negative integer");
    }

    std::size_t consumed = 0;
    unsigned long long value = 0;
    try {
        value = std::stoull(token, &consumed);
    } catch (const std::exception&) {
        throw std::invalid_argument("invalid integer for " + context);
    }
    if (consumed != token.size()) {
        throw std::invalid_argument("invalid integer for " + context);
    }
    return static_cast<std::uint64_t>(value);
}

std::size_t count_value(const std::string& line, const std::string& context)
{
    const std::vector<std::string> fields = tokens(line);
    if (fields.size() != 1U) {
        throw std::invalid_argument(context + " count must be on its own line");
    }
    const std::uint64_t value = unsigned_value(fields.front(), context + " count");
    if (value > std::numeric_limits<std::size_t>::max()) {
        throw std::invalid_argument(context + " count is too large");
    }
    return static_cast<std::size_t>(value);
}

std::uint32_t uint32_value(const std::string& token,
                           const std::string& context)
{
    const std::uint64_t value = unsigned_value(token, context);
    if (value > std::numeric_limits<std::uint32_t>::max()) {
        throw std::invalid_argument(context + " is too large");
    }
    return static_cast<std::uint32_t>(value);
}

double coordinate_value(const std::string& token, const std::string& context)
{
    std::size_t consumed = 0;
    double value = 0.0;
    try {
        value = std::stod(token, &consumed);
    } catch (const std::exception&) {
        throw std::invalid_argument("invalid coordinate for " + context);
    }
    if (consumed != token.size()) {
        throw std::invalid_argument("invalid coordinate for " + context);
    }
    return value;
}

void require_field_count(const std::vector<std::string>& fields,
                         std::size_t expected, const std::string& context)
{
    if (fields.size() != expected) {
        throw std::invalid_argument(context + " has " +
                                    std::to_string(fields.size()) +
                                    " fields; expected " +
                                    std::to_string(expected));
    }
}

NodeId one_based_node(const std::string& token, const std::string& context)
{
    const std::uint32_t value = uint32_value(token, context);
    if (value == 0U) {
        throw std::invalid_argument(context + " must be one-based");
    }
    return value - 1U;
}

std::vector<std::string> nonempty_lines(std::istream& input)
{
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(input, line)) {
        line = trim(line);
        if (!line.empty()) {
            lines.push_back(line);
        }
    }
    return lines;
}

std::size_t section_line(const std::vector<std::string>& lines,
                         const std::string& marker)
{
    const auto found = std::find(lines.begin(), lines.end(), marker);
    if (found == lines.end()) {
        throw std::invalid_argument("Gmsh file has no " + marker + " section");
    }
    return static_cast<std::size_t>(std::distance(lines.begin(), found));
}

void require_section_end(const std::vector<std::string>& lines,
                         std::size_t position, const std::string& marker)
{
    if (position >= lines.size() || lines[position] != marker) {
        throw std::invalid_argument("Gmsh section is missing " + marker);
    }
}

struct RawGmshEntity {
    Label label;
    std::vector<std::uint32_t> node_tags;
};

std::vector<NodeId> map_gmsh_nodes(
    const RawGmshEntity& entity,
    const std::unordered_map<std::uint32_t, NodeId>& node_indices,
    const char* context)
{
    std::vector<NodeId> result;
    result.reserve(entity.node_tags.size());
    for (const std::uint32_t tag : entity.node_tags) {
        const auto found = node_indices.find(tag);
        if (found == node_indices.end()) {
            throw std::invalid_argument(std::string("Gmsh ") + context +
                                        " references unknown node " +
                                        std::to_string(tag));
        }
        result.push_back(found->second);
    }
    return result;
}

bool unsupported_gmsh_facet(std::uint32_t type)
{
    switch (type) {
    case 9U:
    case 10U:
    case 16U:
    case 20U:
    case 21U:
    case 22U:
    case 23U:
    case 24U:
    case 25U:
        return true;
    default:
        return false;
    }
}

bool unsupported_gmsh_cell(std::uint32_t type)
{
    switch (type) {
    case 6U:
    case 7U:
    case 11U:
    case 12U:
    case 13U:
    case 14U:
    case 17U:
    case 18U:
    case 19U:
    case 29U:
    case 30U:
    case 31U:
    case 92U:
    case 93U:
        return true;
    default:
        return false;
    }
}

} // namespace

Mesh read_netgen(std::istream& input)
{
    const std::size_t point_count =
        count_value(required_line(input, "Netgen point count"), "point");
    std::vector<Point> points;
    points.reserve(point_count);
    for (std::size_t index = 0; index < point_count; ++index) {
        const std::vector<std::string> fields =
            tokens(required_line(input, "Netgen point"));
        require_field_count(fields, 3U, "Netgen point");
        points.push_back(Point{{coordinate_value(fields[0], "Netgen point"),
                                coordinate_value(fields[1], "Netgen point"),
                                coordinate_value(fields[2], "Netgen point")}});
    }

    const std::size_t cell_count =
        count_value(required_line(input, "Netgen cell count"), "cell");
    if (cell_count == 0U) {
        throw std::invalid_argument("Netgen mesh has no cells");
    }

    std::vector<Cell> cells;
    cells.reserve(cell_count);
    CellKind kind = CellKind::tetrahedron;
    for (std::size_t index = 0; index < cell_count; ++index) {
        const std::vector<std::string> fields =
            tokens(required_line(input, "Netgen cell"));
        if (index == 0U) {
            if (fields.size() == 5U) {
                kind = CellKind::tetrahedron;
            } else if (fields.size() == 9U) {
                kind = CellKind::hexahedron;
            } else {
                throw std::invalid_argument(
                    "Netgen cells must contain a label and four or eight nodes");
            }
        }
        require_field_count(fields, node_count(kind) + 1U, "Netgen cell");

        std::vector<NodeId> nodes;
        nodes.reserve(node_count(kind));
        for (std::size_t field = 1; field < fields.size(); ++field) {
            nodes.push_back(one_based_node(fields[field], "Netgen cell node"));
        }
        if (kind == CellKind::hexahedron) {
            // Netgen uses tensor-product vertex numbering. Normalize it to
            // the cyclic Gmsh/VTK convention used by the topology core.
            nodes = {nodes[0], nodes[1], nodes[3], nodes[2],
                     nodes[4], nodes[5], nodes[7], nodes[6]};
        }
        cells.emplace_back(kind, uint32_value(fields[0], "Netgen cell label"),
                           std::move(nodes));
    }

    const std::size_t facet_count =
        count_value(required_line(input, "Netgen facet count"), "facet");
    std::vector<Facet> facets;
    facets.reserve(facet_count);
    for (std::size_t index = 0; index < facet_count; ++index) {
        const std::vector<std::string> fields =
            tokens(required_line(input, "Netgen facet"));
        require_field_count(fields, facet_node_count(kind) + 1U,
                            "Netgen facet");

        std::vector<NodeId> nodes;
        nodes.reserve(facet_node_count(kind));
        for (std::size_t field = 1; field < fields.size(); ++field) {
            nodes.push_back(one_based_node(fields[field], "Netgen facet node"));
        }
        facets.emplace_back(uint32_value(fields[0], "Netgen facet label"),
                            std::move(nodes));
    }

    std::string trailing;
    while (std::getline(input, trailing)) {
        if (!trim(trailing).empty()) {
            throw std::invalid_argument("Netgen file has trailing data");
        }
    }
    return Mesh(std::move(points), std::move(cells), std::move(facets));
}

Mesh read_netgen(const std::string& filename)
{
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("could not open Netgen mesh: " + filename);
    }
    return read_netgen(input);
}

Mesh read_gmsh(std::istream& input)
{
    const std::vector<std::string> lines = nonempty_lines(input);

    const std::size_t format = section_line(lines, "$MeshFormat");
    if (format + 2U >= lines.size()) {
        throw std::invalid_argument("incomplete Gmsh mesh format section");
    }
    const std::vector<std::string> format_fields = tokens(lines[format + 1U]);
    require_field_count(format_fields, 3U, "Gmsh mesh format");
    if (format_fields[0] != "2.2" || format_fields[1] != "0") {
        throw std::invalid_argument("only ASCII Gmsh 2.2 files are supported");
    }
    require_section_end(lines, format + 2U, "$EndMeshFormat");

    const std::size_t node_section = section_line(lines, "$Nodes");
    if (node_section + 1U >= lines.size()) {
        throw std::invalid_argument("incomplete Gmsh node section");
    }
    const std::size_t point_count =
        count_value(lines[node_section + 1U], "Gmsh node");
    if (point_count > lines.size() - std::min(lines.size(), node_section + 2U)) {
        throw std::invalid_argument("incomplete Gmsh node section");
    }

    std::vector<Point> points;
    points.reserve(point_count);
    std::unordered_map<std::uint32_t, NodeId> node_indices;
    node_indices.reserve(point_count);
    for (std::size_t index = 0; index < point_count; ++index) {
        const std::vector<std::string> fields =
            tokens(lines[node_section + 2U + index]);
        require_field_count(fields, 4U, "Gmsh node");
        const std::uint32_t tag = uint32_value(fields[0], "Gmsh node tag");
        if (tag == 0U) {
            throw std::invalid_argument("Gmsh node tags must be positive");
        }
        if (!node_indices.emplace(tag, static_cast<NodeId>(index)).second) {
            throw std::invalid_argument("Gmsh file contains duplicate node tags");
        }
        points.push_back(Point{{coordinate_value(fields[1], "Gmsh node"),
                                coordinate_value(fields[2], "Gmsh node"),
                                coordinate_value(fields[3], "Gmsh node")}});
    }
    require_section_end(lines, node_section + 2U + point_count, "$EndNodes");

    const std::size_t element_section = section_line(lines, "$Elements");
    if (element_section + 1U >= lines.size()) {
        throw std::invalid_argument("incomplete Gmsh element section");
    }
    const std::size_t element_count =
        count_value(lines[element_section + 1U], "Gmsh element");
    if (element_count >
        lines.size() - std::min(lines.size(), element_section + 2U)) {
        throw std::invalid_argument("incomplete Gmsh element section");
    }

    std::vector<std::pair<CellKind, RawGmshEntity>> raw_cells;
    std::vector<RawGmshEntity> raw_facets;
    for (std::size_t index = 0; index < element_count; ++index) {
        const std::vector<std::string> fields =
            tokens(lines[element_section + 2U + index]);
        if (fields.size() < 3U) {
            throw std::invalid_argument("Gmsh element has fewer than three fields");
        }
        const std::uint32_t type = uint32_value(fields[1], "Gmsh element type");
        const std::size_t tag_count = static_cast<std::size_t>(
            uint32_value(fields[2], "Gmsh element tag count"));
        if (fields.size() < 3U + tag_count) {
            throw std::invalid_argument("Gmsh element has fewer tags than declared");
        }
        if (unsupported_gmsh_facet(type)) {
            throw std::invalid_argument(
                "higher-order Gmsh facets are not supported");
        }
        if (unsupported_gmsh_cell(type)) {
            throw std::invalid_argument(
                "Gmsh cells must be linear tetrahedra or hexahedra");
        }

        std::size_t entity_node_count = 0U;
        CellKind cell_kind = CellKind::tetrahedron;
        bool is_cell = false;
        bool is_facet = false;
        switch (type) {
        case 2U:
            entity_node_count = 3U;
            is_facet = true;
            break;
        case 3U:
            entity_node_count = 4U;
            is_facet = true;
            break;
        case 4U:
            entity_node_count = 4U;
            cell_kind = CellKind::tetrahedron;
            is_cell = true;
            break;
        case 5U:
            entity_node_count = 8U;
            cell_kind = CellKind::hexahedron;
            is_cell = true;
            break;
        default:
            continue;
        }
        require_field_count(fields, 3U + tag_count + entity_node_count,
                            "Gmsh element");

        RawGmshEntity entity;
        entity.label = tag_count == 0U
                           ? 0U
                           : uint32_value(fields[3], "Gmsh physical tag");
        entity.node_tags.reserve(entity_node_count);
        for (std::size_t field = 3U + tag_count; field < fields.size(); ++field) {
            const std::uint32_t tag =
                uint32_value(fields[field], "Gmsh element node tag");
            if (tag == 0U) {
                throw std::invalid_argument("Gmsh node tags must be positive");
            }
            entity.node_tags.push_back(tag);
        }

        if (is_cell) {
            raw_cells.emplace_back(cell_kind, std::move(entity));
        } else if (is_facet) {
            raw_facets.push_back(std::move(entity));
        }
    }
    require_section_end(lines, element_section + 2U + element_count,
                        "$EndElements");

    if (raw_cells.empty()) {
        throw std::invalid_argument(
            "Gmsh mesh contains no tetrahedral or hexahedral cells");
    }

    std::vector<Cell> cells;
    cells.reserve(raw_cells.size());
    for (const auto& raw : raw_cells) {
        cells.emplace_back(raw.first, raw.second.label,
                           map_gmsh_nodes(raw.second, node_indices, "cell"));
    }
    std::vector<Facet> facets;
    facets.reserve(raw_facets.size());
    for (const RawGmshEntity& raw : raw_facets) {
        facets.emplace_back(raw.label,
                            map_gmsh_nodes(raw, node_indices, "facet"));
    }
    return Mesh(std::move(points), std::move(cells), std::move(facets));
}

Mesh read_gmsh(const std::string& filename)
{
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("could not open Gmsh mesh: " + filename);
    }
    return read_gmsh(input);
}

} // namespace topo
