#include "mesh.hpp"
#include "mesh_adapters.hpp"
#include "topology.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(_core, module)
{
    module.doc() = "Mesher-neutral Topoprocessor mesh and cohomology core";

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

    py::class_<topo::OrientedEdge>(module, "OrientedEdge")
        .def_property_readonly("start", [](const topo::OrientedEdge& edge) {
            return edge.start;
        })
        .def_property_readonly("end", [](const topo::OrientedEdge& edge) {
            return edge.end;
        });

    py::class_<topo::MaterialSelection>(module, "MaterialSelection")
        .def(py::init<std::vector<topo::Label>,
                      std::vector<topo::Label>>(),
             py::arg("conductor_labels") = std::vector<topo::Label>{},
             py::arg("insulator_labels") = std::vector<topo::Label>{})
        .def_readonly("conductor_labels",
                      &topo::MaterialSelection::conductor_labels)
        .def_readonly("insulator_labels",
                      &topo::MaterialSelection::insulator_labels);

    py::class_<topo::CohomologyOptions>(module, "CohomologyOptions")
        .def(py::init<>())
        .def_readwrite("verify_generators",
                       &topo::CohomologyOptions::verify_generators);

    py::class_<topo::CohomologyResult>(module, "CohomologyResult")
        .def_property_readonly("edges", [](const topo::CohomologyResult& result) {
            return result.edges;
        })
        .def_property_readonly(
            "generators", [](const topo::CohomologyResult& result) {
                return result.generators;
            })
        .def_property_readonly("selected_cell_count",
                               [](const topo::CohomologyResult& result) {
                                   return result.selected_cell_count;
                               })
        .def_property_readonly("face_count",
                               [](const topo::CohomologyResult& result) {
                                   return result.face_count;
                               })
        .def_property_readonly("betti_number",
                               &topo::CohomologyResult::betti_number);

    module.def("node_count", &topo::node_count, py::arg("kind"));
    module.def("facet_node_count", &topo::facet_node_count,
               py::arg("kind"));
    module.def(
        "read_netgen",
        static_cast<topo::Mesh (*)(const std::string&)>(&topo::read_netgen),
        py::arg("filename"));
    module.def(
        "read_gmsh",
        static_cast<topo::Mesh (*)(const std::string&)>(&topo::read_gmsh),
        py::arg("filename"));
    module.def("compute_cohomology", &topo::compute_cohomology,
               py::arg("mesh"),
               py::arg("materials") = topo::MaterialSelection{},
               py::arg("options") = topo::CohomologyOptions{},
               py::call_guard<py::gil_scoped_release>());
    module.def(
        "write_h1",
        static_cast<void (*)(const topo::CohomologyResult&,
                             const std::string&, bool)>(&topo::write_h1),
        py::arg("result"), py::arg("filename"),
        py::arg("one_based_nodes") = true);
}
