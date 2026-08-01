#include "topology.hpp"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>

namespace topo {
namespace {

using Integer = std::int64_t;

std::uint64_t magnitude(Integer value)
{
    if (value >= 0) {
        return static_cast<std::uint64_t>(value);
    }
    return static_cast<std::uint64_t>(-(value + 1)) + 1U;
}

std::uint64_t gcd_unsigned(std::uint64_t left, std::uint64_t right)
{
    while (right != 0U) {
        const std::uint64_t remainder = left % right;
        left = right;
        right = remainder;
    }
    return left;
}

Integer checked_add(Integer left, Integer right)
{
    if ((right > 0 && left > std::numeric_limits<Integer>::max() - right) ||
        (right < 0 && left < std::numeric_limits<Integer>::min() - right)) {
        throw std::overflow_error("integer overflow during cohomology reduction");
    }
    return left + right;
}

Integer checked_multiply(Integer left, Integer right)
{
    if (left == 0 || right == 0) {
        return 0;
    }
    if ((left == -1 && right == std::numeric_limits<Integer>::min()) ||
        (right == -1 && left == std::numeric_limits<Integer>::min())) {
        throw std::overflow_error("integer overflow during cohomology reduction");
    }
    if (left > 0) {
        if ((right > 0 && left > std::numeric_limits<Integer>::max() / right) ||
            (right < 0 && right < std::numeric_limits<Integer>::min() / left)) {
            throw std::overflow_error("integer overflow during cohomology reduction");
        }
    } else if ((right > 0 && left < std::numeric_limits<Integer>::min() / right) ||
               (right < 0 && left < std::numeric_limits<Integer>::max() / right)) {
        throw std::overflow_error("integer overflow during cohomology reduction");
    }
    return left * right;
}

Integer checked_negate(Integer value)
{
    if (value == std::numeric_limits<Integer>::min()) {
        throw std::overflow_error("integer overflow during cohomology reduction");
    }
    return -value;
}

class Fraction {
public:
    Fraction(Integer numerator = 0, Integer denominator = 1)
        : numerator_(numerator), denominator_(denominator)
    {
        normalize();
    }

    Integer numerator() const noexcept { return numerator_; }
    Integer denominator() const noexcept { return denominator_; }
    bool zero() const noexcept { return numerator_ == 0; }

    Fraction operator-() const
    {
        return Fraction(checked_negate(numerator_), denominator_);
    }

    Fraction& operator+=(const Fraction& other)
    {
        const std::uint64_t common = gcd_unsigned(
            static_cast<std::uint64_t>(denominator_),
            static_cast<std::uint64_t>(other.denominator_));
        const Integer left_scale = other.denominator_ /
                                   static_cast<Integer>(common);
        const Integer right_scale = denominator_ /
                                    static_cast<Integer>(common);
        numerator_ = checked_add(checked_multiply(numerator_, left_scale),
                                 checked_multiply(other.numerator_, right_scale));
        denominator_ = checked_multiply(denominator_, left_scale);
        normalize();
        return *this;
    }

    Fraction& operator-=(const Fraction& other)
    {
        return *this += -other;
    }

    Fraction& operator*=(const Fraction& other)
    {
        Integer left_numerator = numerator_;
        Integer right_numerator = other.numerator_;
        Integer left_denominator = denominator_;
        Integer right_denominator = other.denominator_;

        const std::uint64_t cross_left = gcd_unsigned(
            magnitude(left_numerator),
            static_cast<std::uint64_t>(right_denominator));
        const std::uint64_t cross_right = gcd_unsigned(
            magnitude(right_numerator),
            static_cast<std::uint64_t>(left_denominator));
        left_numerator /= static_cast<Integer>(cross_left);
        right_denominator /= static_cast<Integer>(cross_left);
        right_numerator /= static_cast<Integer>(cross_right);
        left_denominator /= static_cast<Integer>(cross_right);

        numerator_ = checked_multiply(left_numerator, right_numerator);
        denominator_ = checked_multiply(left_denominator, right_denominator);
        normalize();
        return *this;
    }

    Fraction& operator/=(const Fraction& other)
    {
        if (other.zero()) {
            throw std::logic_error("division by zero during cohomology reduction");
        }
        return *this *= Fraction(other.denominator(), other.numerator());
    }

private:
    void normalize()
    {
        if (denominator_ == 0) {
            throw std::logic_error("zero rational denominator");
        }
        if (numerator_ == 0) {
            denominator_ = 1;
            return;
        }
        if (denominator_ < 0) {
            numerator_ = checked_negate(numerator_);
            denominator_ = checked_negate(denominator_);
        }
        const std::uint64_t divisor = gcd_unsigned(
            magnitude(numerator_), static_cast<std::uint64_t>(denominator_));
        numerator_ /= static_cast<Integer>(divisor);
        denominator_ /= static_cast<Integer>(divisor);
    }

    Integer numerator_;
    Integer denominator_;
};

Fraction operator+(Fraction left, const Fraction& right)
{
    left += right;
    return left;
}

Fraction operator-(Fraction left, const Fraction& right)
{
    left -= right;
    return left;
}

Fraction operator*(Fraction left, const Fraction& right)
{
    left *= right;
    return left;
}

Fraction operator/(Fraction left, const Fraction& right)
{
    left /= right;
    return left;
}

using SparseRow = std::map<std::size_t, Fraction>;

void add_to_row(SparseRow& row, std::size_t column, const Fraction& value)
{
    auto found = row.find(column);
    if (found == row.end()) {
        if (!value.zero()) {
            row.emplace(column, value);
        }
        return;
    }
    found->second += value;
    if (found->second.zero()) {
        row.erase(found);
    }
}

void subtract_multiple(SparseRow& row, const SparseRow& pivot,
                       Fraction factor)
{
    for (const auto& entry : pivot) {
        add_to_row(row, entry.first, -(factor * entry.second));
    }
}

class DisjointSet {
public:
    explicit DisjointSet(std::size_t size) : parent_(size), rank_(size, 0)
    {
        std::iota(parent_.begin(), parent_.end(), 0U);
    }

    std::size_t find(std::size_t value)
    {
        if (parent_[value] != value) {
            parent_[value] = find(parent_[value]);
        }
        return parent_[value];
    }

    bool join(std::size_t left, std::size_t right)
    {
        left = find(left);
        right = find(right);
        if (left == right) {
            return false;
        }
        if (rank_[left] < rank_[right]) {
            std::swap(left, right);
        }
        parent_[right] = left;
        if (rank_[left] == rank_[right]) {
            ++rank_[left];
        }
        return true;
    }

private:
    std::vector<std::size_t> parent_;
    std::vector<unsigned char> rank_;
};

using Face = std::vector<NodeId>;

std::vector<Face> cell_faces(const Cell& cell)
{
    const std::vector<NodeId>& n = cell.nodes();
    if (cell.kind() == CellKind::tetrahedron) {
        return {{n[0], n[2], n[1]}, {n[0], n[1], n[3]},
                {n[1], n[2], n[3]}, {n[2], n[0], n[3]}};
    }
    // Canonical hexahedron ordering follows Gmsh/VTK: a cyclic bottom face,
    // followed by the corresponding cyclic top face.
    return {{n[0], n[3], n[2], n[1]}, {n[4], n[5], n[6], n[7]},
            {n[0], n[1], n[5], n[4]}, {n[1], n[2], n[6], n[5]},
            {n[2], n[3], n[7], n[6]}, {n[3], n[0], n[4], n[7]}};
}

bool contains(const std::set<Label>& labels, Label label)
{
    return labels.find(label) != labels.end();
}

std::vector<bool> selected_cells(const Mesh& mesh,
                                 const MaterialSelection& materials)
{
    const std::set<Label> conductors(materials.conductor_labels.begin(),
                                     materials.conductor_labels.end());
    const std::set<Label> insulators(materials.insulator_labels.begin(),
                                     materials.insulator_labels.end());
    for (const Label label : conductors) {
        if (contains(insulators, label)) {
            throw std::invalid_argument(
                "a material label cannot be both conductor and insulator");
        }
    }

    std::vector<bool> selected;
    selected.reserve(mesh.cells().size());
    for (const Cell& cell : mesh.cells()) {
        const Label label = cell.label();
        if (contains(conductors, label)) {
            selected.push_back(false);
        } else if (insulators.empty() || contains(insulators, label)) {
            selected.push_back(true);
        } else {
            throw std::invalid_argument("mesh material label " +
                                        std::to_string(label) +
                                        " is not classified");
        }
    }
    if (std::find(selected.begin(), selected.end(), true) == selected.end()) {
        throw std::invalid_argument("material selection contains no insulating cells");
    }
    return selected;
}

struct Complex {
    std::vector<OrientedEdge> edges;
    std::vector<std::vector<std::pair<std::size_t, Integer>>> face_edges;
    std::size_t selected_cells = 0;
};

Complex build_complex(const Mesh& mesh, const MaterialSelection& materials)
{
    const std::vector<bool> selected = selected_cells(mesh, materials);
    std::map<Face, std::pair<Face, std::size_t>> faces;
    std::set<std::vector<NodeId>> cell_keys;

    for (std::size_t index = 0; index < mesh.cells().size(); ++index) {
        if (!selected[index]) {
            continue;
        }
        std::vector<NodeId> cell_key = mesh.cells()[index].nodes();
        std::sort(cell_key.begin(), cell_key.end());
        if (!cell_keys.insert(cell_key).second) {
            throw std::invalid_argument("insulating submesh contains a duplicate cell");
        }
        for (const Face& face : cell_faces(mesh.cells()[index])) {
            Face key = face;
            std::sort(key.begin(), key.end());
            auto found = faces.find(key);
            if (found == faces.end()) {
                faces.emplace(std::move(key), std::make_pair(face, 1U));
            } else if (++found->second.second > 2U) {
                throw std::invalid_argument(
                    "insulating submesh contains a non-manifold face");
            }
        }
    }

    std::map<std::pair<NodeId, NodeId>, std::size_t> edge_indices;
    for (const auto& face_entry : faces) {
        const Face& face = face_entry.second.first;
        for (std::size_t i = 0; i < face.size(); ++i) {
            NodeId start = face[i];
            NodeId end = face[(i + 1U) % face.size()];
            if (end < start) {
                std::swap(start, end);
            }
            edge_indices.emplace(std::make_pair(start, end), 0U);
        }
    }

    Complex complex;
    complex.selected_cells = static_cast<std::size_t>(
        std::count(selected.begin(), selected.end(), true));
    complex.edges.reserve(edge_indices.size());
    std::size_t edge_index = 0;
    for (auto& entry : edge_indices) {
        entry.second = edge_index++;
        complex.edges.push_back({entry.first.first, entry.first.second});
    }

    complex.face_edges.reserve(faces.size());
    for (const auto& face_entry : faces) {
        const Face& face = face_entry.second.first;
        std::vector<std::pair<std::size_t, Integer>> boundary;
        boundary.reserve(face.size());
        for (std::size_t i = 0; i < face.size(); ++i) {
            const NodeId start = face[i];
            const NodeId end = face[(i + 1U) % face.size()];
            const auto key = std::make_pair(std::min(start, end),
                                            std::max(start, end));
            boundary.emplace_back(edge_indices.at(key), start < end ? 1 : -1);
        }
        complex.face_edges.push_back(std::move(boundary));
    }
    return complex;
}

std::vector<std::vector<Integer>> nullspace(const Complex& complex,
                                            std::size_t point_count)
{
    DisjointSet components(point_count);
    std::vector<bool> tree_edge(complex.edges.size(), false);
    for (std::size_t index = 0; index < complex.edges.size(); ++index) {
        const OrientedEdge& edge = complex.edges[index];
        tree_edge[index] = components.join(edge.start, edge.end);
    }

    std::vector<std::size_t> chord_column(complex.edges.size(),
                                          std::numeric_limits<std::size_t>::max());
    std::vector<std::size_t> chord_edge;
    for (std::size_t index = 0; index < complex.edges.size(); ++index) {
        if (!tree_edge[index]) {
            chord_column[index] = chord_edge.size();
            chord_edge.push_back(index);
        }
    }

    std::map<std::size_t, SparseRow> pivots;
    for (const auto& face : complex.face_edges) {
        SparseRow row;
        for (const auto& entry : face) {
            const std::size_t column = chord_column[entry.first];
            if (column != std::numeric_limits<std::size_t>::max()) {
                add_to_row(row, column, Fraction(entry.second));
            }
        }
        while (!row.empty()) {
            const std::size_t column = row.begin()->first;
            const auto pivot = pivots.find(column);
            if (pivot != pivots.end()) {
                subtract_multiple(row, pivot->second, row.begin()->second);
                continue;
            }
            const Fraction scale = row.begin()->second;
            for (auto& item : row) {
                item.second /= scale;
            }
            pivots.emplace(column, std::move(row));
            break;
        }
    }

    std::set<std::size_t> pivot_columns;
    for (const auto& pivot : pivots) {
        pivot_columns.insert(pivot.first);
    }

    std::vector<std::vector<Integer>> basis;
    for (std::size_t free_column = 0; free_column < chord_edge.size();
         ++free_column) {
        if (pivot_columns.find(free_column) != pivot_columns.end()) {
            continue;
        }
        std::vector<Fraction> values(chord_edge.size());
        values[free_column] = Fraction(1);
        for (auto pivot = pivots.rbegin(); pivot != pivots.rend(); ++pivot) {
            Fraction sum;
            for (const auto& entry : pivot->second) {
                if (entry.first != pivot->first) {
                    sum += entry.second * values[entry.first];
                }
            }
            values[pivot->first] = -sum;
        }

        for (const Fraction& value : values) {
            if (value.denominator() != 1) {
                throw std::invalid_argument(
                    "the low-order cell complex is not integrally unimodular; "
                    "an integral H^1 basis cannot be identified safely");
            }
        }

        std::vector<Integer> generator(complex.edges.size(), 0);
        Integer common_factor = 0;
        for (std::size_t column = 0; column < values.size(); ++column) {
            const Integer coefficient = values[column].numerator();
            generator[chord_edge[column]] = coefficient;
            common_factor = static_cast<Integer>(gcd_unsigned(
                magnitude(common_factor), magnitude(coefficient)));
        }
        if (common_factor > 1) {
            for (Integer& coefficient : generator) {
                coefficient /= common_factor;
            }
        }
        const auto first = std::find_if(generator.begin(), generator.end(),
                                        [](Integer value) { return value != 0; });
        if (first != generator.end() && *first < 0) {
            for (Integer& coefficient : generator) {
                coefficient = checked_negate(coefficient);
            }
        }
        basis.push_back(std::move(generator));
    }
    return basis;
}

void verify(const Complex& complex,
            const std::vector<std::vector<Integer>>& generators)
{
    for (const auto& generator : generators) {
        if (generator.size() != complex.edges.size()) {
            throw std::logic_error("cohomology generator has the wrong size");
        }
        for (const auto& face : complex.face_edges) {
            Integer sum = 0;
            for (const auto& entry : face) {
                sum = checked_add(sum,
                                  checked_multiply(entry.second,
                                                   generator[entry.first]));
            }
            if (sum != 0) {
                throw std::logic_error("computed generator is not a cocycle");
            }
        }
    }
}

} // namespace

MaterialSelection::MaterialSelection(std::vector<Label> conductors,
                                     std::vector<Label> insulators)
    : conductor_labels(std::move(conductors)),
      insulator_labels(std::move(insulators))
{
}

CohomologyResult compute_cohomology(const Mesh& mesh,
                                    const MaterialSelection& materials,
                                    const CohomologyOptions& options)
{
    mesh.validate();
    const Complex complex = build_complex(mesh, materials);
    CohomologyResult result;
    result.edges = complex.edges;
    result.generators = nullspace(complex, mesh.points().size());
    result.selected_cell_count = complex.selected_cells;
    result.face_count = complex.face_edges.size();
    if (options.verify_generators) {
        verify(complex, result.generators);
    }
    return result;
}

void write_h1(const CohomologyResult& result, std::ostream& output,
              bool one_based_nodes)
{
    output << result.generators.size() << '\n';
    for (const auto& generator : result.generators) {
        const std::size_t support = static_cast<std::size_t>(std::count_if(
            generator.begin(), generator.end(),
            [](Integer coefficient) { return coefficient != 0; }));
        output << support << '\n';
        for (std::size_t edge = 0; edge < generator.size(); ++edge) {
            if (generator[edge] == 0) {
                continue;
            }
            const NodeId offset = one_based_nodes ? 1U : 0U;
            output << result.edges[edge].start + offset << '\t'
                   << result.edges[edge].end + offset << '\t'
                   << generator[edge] << '\n';
        }
    }
    if (!output) {
        throw std::runtime_error("failed to write cohomology result");
    }
}

void write_h1(const CohomologyResult& result, const std::string& filename,
              bool one_based_nodes)
{
    std::ofstream output(filename);
    if (!output) {
        throw std::runtime_error("could not open output file: " + filename);
    }
    write_h1(result, output, one_based_nodes);
}

} // namespace topo
