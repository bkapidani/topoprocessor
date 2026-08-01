#include "mesh_adapters.hpp"
#include "topology.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv)
{
    if (argc != 2) {
        throw std::runtime_error("fixture root argument is required");
    }
    const std::string root = argv[1];

    const topo::Mesh furch = topo::read_netgen(root + "/netgen/furch.mesh");
    const topo::CohomologyResult furch_result = topo::compute_cohomology(
        furch, topo::MaterialSelection({1U}));
    std::cout << "furch b1=" << furch_result.betti_number() << '\n';
    if (furch_result.betti_number() != 1U) {
        throw std::runtime_error("furch complement should have first Betti number one");
    }

    const topo::Mesh torus = topo::read_netgen(root + "/netgen/torus.mesh");
    const topo::CohomologyResult torus_result = topo::compute_cohomology(
        torus, topo::MaterialSelection({2U}));
    std::cout << "torus b1=" << torus_result.betti_number() << '\n';

    if (torus_result.betti_number() != 1U) {
        throw std::runtime_error("torus complement should have first Betti number one");
    }
}
