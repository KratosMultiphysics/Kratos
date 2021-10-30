//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "python/add_bounding_box_to_python.h"
#include "geometries/bounding_box.h"
#include "includes/model_part.h"

namespace Kratos
{

namespace Python
{

typedef BoundingBox<Point> BoundingBoxType;
typedef ModelPart::NodesContainerType NodesContainerType;

BoundingBoxType BoundingBoxConstructorWithNodes(NodesContainerType& rNodes)
{
    return BoundingBoxType(rNodes.begin(), rNodes.end());
}

void Set(BoundingBoxType& rBoundingBox, NodesContainerType& rNodes)
{
    rBoundingBox.Set(rNodes.begin(), rNodes.end());
}

void Extend1(BoundingBoxType& rBoundingBox, NodesContainerType& rNodes)
{
    rBoundingBox.Extend(rNodes.begin(), rNodes.end());
}

void Extend2(BoundingBoxType& rBoundingBox, double Margin)
{
    rBoundingBox.Extend(Margin);
}

Point& GetMinPoint(BoundingBoxType& rBoundingBox)
{
    return rBoundingBox.GetMinPoint();
}

Point& GetMaxPoint(BoundingBoxType& rBoundingBox)
{
    return rBoundingBox.GetMaxPoint();
}

void AddBoundingBoxToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<BoundingBoxType, BoundingBoxType::Pointer >(m, "BoundingBox")
    // .def(py::init<>(BoundingBoxConstructorWithNodes))
    .def(py::init([](NodesContainerType& rNodes){
        return BoundingBoxType(rNodes.begin(), rNodes.end());
    }))
    .def(py::init<Point, Point>())
    .def("Set", Set)
    .def("Extend", Extend1)
    .def("Extend", Extend2)
    .def("GetMinPoint", GetMinPoint)
    .def("GetMaxPoint", GetMaxPoint)
    ;

}

}  // namespace Python.

} // Namespace Kratos
