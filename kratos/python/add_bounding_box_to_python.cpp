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
    std::vector<Point> points;
    for (auto it_node = rNodes.begin(); it_node < rNodes.end(); ++it_node) {
        points.push_back(*it_node);
    }
    return BoundingBoxType(points.begin(), points.end());
}

void Extend1(BoundingBoxType& rBoundingBox, NodesContainerType& rNodes)
{
    std::vector<Point> points;
    for (auto it_node = rNodes.begin(); it_node < rNodes.end(); ++it_node) {
        points.push_back(*it_node);
    }
    rBoundingBox.Extend(points.begin(), points.end());
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
        std::vector<Point> points;
        for (auto it_node = rNodes.begin(); it_node < rNodes.end(); ++it_node) {
            points.push_back(*it_node);
        }
        return BoundingBoxType(points.begin(), points.end());
        // return Kratos::make_shared<BoundingBoxType>(new BoundingBoxType(rNodes.begin(), rNodes.end()));
    }))
    .def(py::init<Point, Point>())
    // .def("Set", &BoundingBoxType::Set)
    .def("Extend", Extend1)
    .def("Extend", Extend2)
    .def("GetMinPoint", GetMinPoint)
    .def("GetMaxPoint", GetMaxPoint)
    ;

}

}  // namespace Python.

} // Namespace Kratos
