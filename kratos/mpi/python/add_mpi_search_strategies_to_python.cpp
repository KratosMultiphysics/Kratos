//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/model_part.h"
#include "mpi/python/add_mpi_search_strategies_to_python.h"
#include "mpi/spatial_containers/geometrical_objects_bins_mpi.h"

namespace Kratos::Python
{

void AddMPISearchStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using NodesContainerType = ModelPart::NodesContainerType;
    using ElementsContainerType = ModelPart::ElementsContainerType;
    using ConditionsContainerType = ModelPart::ConditionsContainerType;
    using ResultTypeContainerMap = SpatialSearchResultContainerMap<GeometricalObject>;

    py::class_<GeometricalObjectsBinsMPI, GeometricalObjectsBinsMPI::Pointer>(m, "GeometricalObjectsBinsMPI")
    .def(py::init<ElementsContainerType&, const DataCommunicator&>())
    .def(py::init<ConditionsContainerType&, const DataCommunicator&>())
    .def("GetBoundingBox", &GeometricalObjectsBinsMPI::GetBoundingBox)
    .def("SearchInRadius", [&](GeometricalObjectsBinsMPI& self, const NodesContainerType& rNodes, const double Radius) {
        // Perform the search
        ResultTypeContainerMap results; 
        self.SearchInRadius(rNodes.begin(), rNodes.end(), Radius, results);
        return results;
    })
    .def("SearchNearestInRadius", [&](GeometricalObjectsBinsMPI& self, const NodesContainerType& rNodes, const double Radius) {
        // Perform the search
        ResultTypeContainerMap results;
        self.SearchNearestInRadius(rNodes.begin(), rNodes.end(), Radius, results);
        return results;
    })
    .def("SearchNearest", [&](GeometricalObjectsBinsMPI& self, const NodesContainerType& rNodes) {
        // Perform the search
        ResultTypeContainerMap results;
        self.SearchNearest(rNodes.begin(), rNodes.end(), results);
        return results;
    })
    .def("SearchIsInside", [&](GeometricalObjectsBinsMPI& self, const NodesContainerType& rNodes) {
        // Perform the search
        ResultTypeContainerMap results;
        self.SearchIsInside(rNodes.begin(), rNodes.end(), results);
        return results;
    })
    ;
}

}