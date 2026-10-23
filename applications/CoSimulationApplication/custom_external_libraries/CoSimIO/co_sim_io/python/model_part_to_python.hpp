//     ______     _____ _           ________
//    / ____/___ / ___/(_)___ ___  /  _/ __ |
//   / /   / __ \\__ \/ / __ `__ \ / // / / /
//  / /___/ /_/ /__/ / / / / / / // // /_/ /
//  \____/\____/____/_/_/ /_/ /_/___/\____/
//  Kratos CoSimulationApplication
//
//  License:         BSD License, see license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#ifndef CO_SIM_IO_MODEL_PART_TO_PYHON_INCLUDED
#define CO_SIM_IO_MODEL_PART_TO_PYHON_INCLUDED

// Exposure of the CoSimIO to Python

// System includes
#include <string>
#include <sstream>

// pybind includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// CoSimIO include
#include "includes/model_part.hpp"

namespace CoSimIO {

void AddCoSimIOModelPartToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<CoSimIO::Node, CoSimIO::ModelPart::NodePointerType>(m,"Node")
        .def(py::init<const CoSimIO::IdType, const double, const double, const double>())
        .def(py::init<const CoSimIO::IdType, const CoSimIO::CoordinatesType&>())
        .def("Id", &CoSimIO::Node::Id)
        .def("X",  &CoSimIO::Node::X)
        .def("Y",  &CoSimIO::Node::Y)
        .def("Z",  &CoSimIO::Node::Z)
        .def("Coordinates",  &CoSimIO::Node::Coordinates)
        .def("__str__",   [](const CoSimIO::Node& I_Node)
            { std::stringstream ss; ss << I_Node; return ss.str(); } )
        ;

    py::class_<CoSimIO::Element, CoSimIO::ModelPart::ElementPointerType>(m,"Element")
        .def(py::init<const CoSimIO::IdType, const CoSimIO::ElementType, const CoSimIO::Element::NodesContainerType&>())
        .def("Id", &CoSimIO::Element::Id)
        .def("Type", &CoSimIO::Element::Type)
        .def("NumberOfNodes", &CoSimIO::Element::NumberOfNodes)
        .def_property_readonly("Nodes", [](CoSimIO::Element& I_Element) {
            return py::make_iterator(I_Element.NodesBegin(), I_Element.NodesEnd());
            }, py::keep_alive<0, 1>()) /* Keep vector alive while iterator is used */
        .def("__str__",   [](const CoSimIO::Element& I_Element)
            { std::stringstream ss; ss << I_Element; return ss.str(); } )
        ;

    py::class_<CoSimIO::ModelPart>(m,"ModelPart")
        .def(py::init<const std::string&>())
        .def("Name",                  &CoSimIO::ModelPart::Name)
        .def("NumberOfNodes",         &CoSimIO::ModelPart::NumberOfNodes)
        .def("NumberOfLocalNodes",    &CoSimIO::ModelPart::NumberOfLocalNodes)
        .def("NumberOfGhostNodes",    &CoSimIO::ModelPart::NumberOfGhostNodes)
        .def("NumberOfElements",      &CoSimIO::ModelPart::NumberOfElements)
        .def("CreateNewNode",         &CoSimIO::ModelPart::CreateNewNode, py::return_value_policy::reference_internal)
        .def("CreateNewNodes",        &CoSimIO::ModelPart::CreateNewNodes<std::vector<CoSimIO::IdType>,
                                                                          std::vector<double>>)
        .def("CreateNewGhostNode",    &CoSimIO::ModelPart::CreateNewGhostNode, py::return_value_policy::reference_internal)
        .def("CreateNewGhostNodes",   &CoSimIO::ModelPart::CreateNewGhostNodes<std::vector<CoSimIO::IdType>,
                                                                               std::vector<double>,
                                                                               std::vector<int>>)
        .def("CreateNewElement",      &CoSimIO::ModelPart::CreateNewElement, py::return_value_policy::reference_internal)
        .def("CreateNewElements",     &CoSimIO::ModelPart::CreateNewElements<std::vector<CoSimIO::IdType>,
                                                                             std::vector<CoSimIO::ElementType>,
                                                                             std::vector<CoSimIO::IdType>>)
        .def("GetNode",               [](CoSimIO::ModelPart& I_ModelPart, const CoSimIO::IdType I_Id){
            return I_ModelPart.pGetNode(I_Id);}, py::return_value_policy::reference_internal)
        .def("GetElement",            [](CoSimIO::ModelPart& I_ModelPart, const CoSimIO::IdType I_Id){
            return I_ModelPart.pGetElement(I_Id);}, py::return_value_policy::reference_internal)
        .def_property_readonly("Nodes", [](CoSimIO::ModelPart& I_ModelPart) {
            return py::make_iterator(I_ModelPart.NodesBegin(), I_ModelPart.NodesEnd());
            }, py::keep_alive<0, 1>()) /* Keep vector alive while iterator is used */
        .def_property_readonly("LocalNodes", [](CoSimIO::ModelPart& I_ModelPart) {
            return py::make_iterator(I_ModelPart.LocalNodesBegin(), I_ModelPart.LocalNodesEnd());
            }, py::keep_alive<0, 1>()) /* Keep vector alive while iterator is used */
        .def_property_readonly("GhostNodes", [](CoSimIO::ModelPart& I_ModelPart) {
            return py::make_iterator(I_ModelPart.GhostNodesBegin(), I_ModelPart.GhostNodesEnd());
            }, py::keep_alive<0, 1>()) /* Keep vector alive while iterator is used */
        .def_property_readonly("Elements", [](CoSimIO::ModelPart& I_ModelPart) {
            return py::make_iterator(I_ModelPart.ElementsBegin(), I_ModelPart.ElementsEnd());
            }, py::keep_alive<0, 1>()) /* Keep vector alive while iterator is used */
        .def("Clear", &CoSimIO::ModelPart::Clear)
        .def("__str__",   [](const CoSimIO::ModelPart& I_ModelPart)
            { std::stringstream ss; ss << I_ModelPart; return ss.str(); } )
        ;

    py::enum_<CoSimIO::ElementType>(m,"ElementType")
        .value("Hexahedra3D20",        CoSimIO::ElementType::Hexahedra3D20)
        .value("Hexahedra3D27",        CoSimIO::ElementType::Hexahedra3D27)
        .value("Hexahedra3D8",         CoSimIO::ElementType::Hexahedra3D8)
        .value("Prism3D15",            CoSimIO::ElementType::Prism3D15)
        .value("Prism3D6",             CoSimIO::ElementType::Prism3D6)
        .value("Pyramid3D13",          CoSimIO::ElementType::Pyramid3D13)
        .value("Pyramid3D5",           CoSimIO::ElementType::Pyramid3D5)
        .value("Quadrilateral2D4",     CoSimIO::ElementType::Quadrilateral2D4)
        .value("Quadrilateral2D8",     CoSimIO::ElementType::Quadrilateral2D8)
        .value("Quadrilateral2D9",     CoSimIO::ElementType::Quadrilateral2D9)
        .value("Quadrilateral3D4",     CoSimIO::ElementType::Quadrilateral3D4)
        .value("Quadrilateral3D8",     CoSimIO::ElementType::Quadrilateral3D8)
        .value("Quadrilateral3D9",     CoSimIO::ElementType::Quadrilateral3D9)
        .value("Tetrahedra3D10",       CoSimIO::ElementType::Tetrahedra3D10)
        .value("Tetrahedra3D4",        CoSimIO::ElementType::Tetrahedra3D4)
        .value("Triangle2D3",          CoSimIO::ElementType::Triangle2D3)
        .value("Triangle2D6",          CoSimIO::ElementType::Triangle2D6)
        .value("Triangle3D3",          CoSimIO::ElementType::Triangle3D3)
        .value("Triangle3D6",          CoSimIO::ElementType::Triangle3D6)
        .value("Line2D2",              CoSimIO::ElementType::Line2D2)
        .value("Line2D3",              CoSimIO::ElementType::Line2D3)
        .value("Line3D2",              CoSimIO::ElementType::Line3D2)
        .value("Line3D3",              CoSimIO::ElementType::Line3D3)
        .value("Point2D",              CoSimIO::ElementType::Point2D)
        .value("Point3D",              CoSimIO::ElementType::Point3D)
        ;
}

} // namespace CoSimIO

#endif // CO_SIM_IO_MODEL_PART_TO_PYHON_INCLUDED
