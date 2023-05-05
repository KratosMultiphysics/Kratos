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

    using ResultType = SpatialSearchResult<GeometricalObject>;
    using NodesContainerType = ModelPart::NodesContainerType;
    using ElementsContainerType = ModelPart::ElementsContainerType;
    using ConditionsContainerType = ModelPart::ConditionsContainerType;

    py::class_<GeometricalObjectsBinsMPI, GeometricalObjectsBinsMPI::Pointer>(m, "GeometricalObjectsBinsMPI")
    .def(py::init<ElementsContainerType&, const DataCommunicator&>())
    .def(py::init<ConditionsContainerType&, const DataCommunicator&>())
    .def("GetBoundingBox", &GeometricalObjectsBinsMPI::GetBoundingBox)
    .def("SearchInRadius", [&](GeometricalObjectsBinsMPI& self, const NodesContainerType& rNodes, const double Radius) {
        // Perform the search
        std::vector<std::vector<ResultType>> results;
        self.SearchInRadius(rNodes.begin(), rNodes.end(), Radius, results);

        // Copy the results to the python list
        py::list list_results;
        for (auto& r_result : results) {
            py::list sub_list_results;
            for (auto& r_sub_result : r_result) {
                sub_list_results.append(r_sub_result);
            }
            list_results.append(sub_list_results);
        }
        return list_results;
    })
    .def("SearchNearestInRadius", [&](GeometricalObjectsBinsMPI& self, const NodesContainerType& rNodes, const double Radius) {
        // Perform the search
        std::vector<ResultType> results = self.SearchNearestInRadius(rNodes.begin(), rNodes.end(), Radius);

        // Copy the results to the python list
        py::list list_results;
        for (auto& r_result : results) {
            list_results.append(r_result);
        }
        return list_results;
    })
    .def("SearchNearest", [&](GeometricalObjectsBinsMPI& self, const NodesContainerType& rNodes) {
        // Perform the search
        std::vector<ResultType> results = self.SearchNearest(rNodes.begin(), rNodes.end());

        // Copy the results to the python list
        py::list list_results;
        for (auto& r_result : results) {
            list_results.append(r_result);
        }
        return list_results;
    })
    .def("SearchIsInside", [&](GeometricalObjectsBinsMPI& self, const NodesContainerType& rNodes) {
        // Perform the search
        std::vector<ResultType> results = self.SearchIsInside(rNodes.begin(), rNodes.end());

        // Copy the results to the python list
        py::list list_results;
        for (auto& r_result : results) {
            list_results.append(r_result);
        }
        return list_results;
    })
    ;
}

}