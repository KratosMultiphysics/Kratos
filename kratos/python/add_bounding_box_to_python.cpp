//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes

// External includes

// Project includes
#include "python/add_bounding_box_to_python.h"
#include "geometries/bounding_box.h"
#include "includes/model_part.h"

namespace Kratos::Python
{

using BoundingBoxType = BoundingBox<Point>;
using NodesContainerType = ModelPart::NodesContainerType;

void AddBoundingBoxToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<BoundingBoxType, BoundingBoxType::Pointer >(m, "BoundingBox")
    .def(py::init<Point, Point>())
    .def(py::init([](const NodesContainerType& rNodes){
        return BoundingBoxType(rNodes.begin(), rNodes.end());
    }))
    .def("Set", [](BoundingBoxType& rBoundingBox, const NodesContainerType& rNodes){
        rBoundingBox.Set(rNodes.begin(), rNodes.end());
    })
    .def("Extend", [](BoundingBoxType& rBoundingBox, double Margin){
        rBoundingBox.Extend(Margin);
    })
    .def("Extend", [](BoundingBoxType& rBoundingBox, const NodesContainerType& rNodes){
        rBoundingBox.Extend(rNodes.begin(), rNodes.end());
    })
    .def("GetMinPoint", [](BoundingBoxType& rBoundingBox){
        return rBoundingBox.GetMinPoint();
    })
    .def("GetMaxPoint", [](BoundingBoxType& rBoundingBox){
        return rBoundingBox.GetMaxPoint();
    })
    .def("GetPoints", [](BoundingBoxType& rBoundingBox){
        return std::make_tuple(rBoundingBox.GetMinPoint(), rBoundingBox.GetMaxPoint());
    })
    ;
}

}  // namespace Kratos::Python.
