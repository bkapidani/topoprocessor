#ifndef TOPOPROCESSOR_TOPOLOGY_HPP
#define TOPOPROCESSOR_TOPOLOGY_HPP

#include "mesh.hpp"

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

namespace topo {

struct OrientedEdge {
    NodeId start;
    NodeId end;
};

struct MaterialSelection {
    std::vector<Label> conductor_labels;
    std::vector<Label> insulator_labels;

    MaterialSelection(std::vector<Label> conductors = {},
                      std::vector<Label> insulators = {});
};

struct CohomologyOptions {
    bool verify_generators = true;
};

struct CohomologyResult {
    std::vector<OrientedEdge> edges;
    std::vector<std::vector<std::int64_t>> generators;
    std::size_t selected_cell_count = 0;
    std::size_t face_count = 0;

    std::size_t betti_number() const noexcept { return generators.size(); }
};

// Compute an integral basis of H^1 of the insulating cell subcomplex.
//
// With an empty selection every cell is insulating. If only conductor labels
// are supplied, all other labels are insulating (the historical CLI rule). If
// insulator labels are supplied, every mesh label must be listed explicitly as
// either insulating or conducting.
CohomologyResult compute_cohomology(
    const Mesh& mesh, const MaterialSelection& materials = {},
    const CohomologyOptions& options = {});

void write_h1(const CohomologyResult& result, std::ostream& output,
              bool one_based_nodes = true);
void write_h1(const CohomologyResult& result, const std::string& filename,
              bool one_based_nodes = true);

} // namespace topo

#endif
