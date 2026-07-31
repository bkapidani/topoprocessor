#include "mesh.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(topoprocessor, module)
{
    module.doc() = "Mesher-neutral Topoprocessor mesh model";

    py::enum_<topo::CellKind>(module, "CellKind")
        .value("tetrahedron", topo::CellKind::tetrahedron)
        .value("hexahedron", topo::CellKind::hexahedron)
        .export_values();

    py::class_<topo::Cell>(module, "Cell")
        .def(py::init<topo::CellKind, topo::Label,
                      std::vector<topo::NodeId>>(),
             py::arg("kind"), py::arg("label"), py::arg("nodes"))
        .def_property_readonly("kind", &topo::Cell::kind)
        .def_property_readonly("label", &topo::Cell::label)
        .def_property_readonly("nodes", &topo::Cell::nodes);

    py::class_<topo::Facet>(module, "Facet")
        .def(py::init<topo::Label, std::vector<topo::NodeId>>(),
             py::arg("label"), py::arg("nodes"))
        .def_property_readonly("label", &topo::Facet::label)
        .def_property_readonly("nodes", &topo::Facet::nodes);

    py::class_<topo::Mesh>(module, "Mesh")
        .def(py::init<std::vector<topo::Point>, std::vector<topo::Cell>,
                      std::vector<topo::Facet>>(),
             py::arg("points"), py::arg("cells"),
             py::arg("boundary_facets"))
        .def_property_readonly("points", &topo::Mesh::points)
        .def_property_readonly("cells", &topo::Mesh::cells)
        .def_property_readonly("boundary_facets",
                               &topo::Mesh::boundary_facets)
        .def_property_readonly("cell_kind", &topo::Mesh::cell_kind)
        .def("validate", &topo::Mesh::validate);

    module.def("node_count", &topo::node_count, py::arg("kind"));
    module.def("facet_node_count", &topo::facet_node_count,
               py::arg("kind"));
}
