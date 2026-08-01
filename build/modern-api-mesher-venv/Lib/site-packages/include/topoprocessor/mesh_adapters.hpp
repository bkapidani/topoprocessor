#ifndef TOPOPROCESSOR_MESH_ADAPTERS_HPP
#define TOPOPROCESSOR_MESH_ADAPTERS_HPP

#include "mesh.hpp"

#include <iosfwd>
#include <string>

namespace topo {

Mesh read_netgen(std::istream& input);
Mesh read_netgen(const std::string& filename);

Mesh read_gmsh(std::istream& input);
Mesh read_gmsh(const std::string& filename);

} // namespace topo

#endif
