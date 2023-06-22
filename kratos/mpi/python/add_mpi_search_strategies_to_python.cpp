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
#include "mpi/spatial_containers/specialized_spatial_search_mpi.h"
#include "mpi/spatial_containers/specialized_spatial_search_mpi_factory.h"

namespace Kratos::Python
{

/**
 * @brief Copies a list of doubles to a radius array
 * @param rListOfRadius The list of doubles to be copied 
 * @param rRadiusArray The radius array to be filled
 */
void CopyRadiusArrayToPython(
    const pybind11::list& rListOfRadius,
    SpatialSearch::RadiusArrayType& rRadiusArray
    )
{
    // Get the size of the radius array
    const std::size_t size_array = rListOfRadius.size();

    // Create the radius array
    rRadiusArray.resize(size_array);
    IndexPartition<std::size_t>(size_array).for_each([&](std::size_t i) {
        rRadiusArray[i] = rListOfRadius[i].cast<double>();
        //rRadiusArray[i] = rListOfRadius[i];
    });
}

/**
 * @brief Copies a Python list of radius to a C++ radius array.
 * @param rListOfRadius list of radius to copy
 * @return Radius array with copied radius
 */
SpatialSearch::RadiusArrayType CopyRadiusArrayToPython(const pybind11::list& rListOfRadius)
{
    // Create the radius array
    SpatialSearch::RadiusArrayType radius_array(rListOfRadius.size());
    CopyRadiusArrayToPython(rListOfRadius, radius_array);
    return radius_array;
}

template <SpatialContainer TSearchBackend>
void DefineSpecializedSpatialSearchMPI(pybind11::module& m, const std::string& rClassName)
{
    namespace py = pybind11;

    using SpatialSearchType = SpecializedSpatialSearchMPI<TSearchBackend>;
    using SpatialSearchPointerType = typename SpecializedSpatialSearchMPI<TSearchBackend>::Pointer;
    using BaseSpatialSearchType = SpecializedSpatialSearch<TSearchBackend>;

    using NodesContainerType = SpatialSearch::NodesContainerType;
    using ElementsContainerType = SpatialSearch::ElementsContainerType;
    using ConditionsContainerType = SpatialSearch::ConditionsContainerType;

    py::class_<SpatialSearchType, SpatialSearchPointerType, BaseSpatialSearchType>(m, rClassName.c_str())
    .def(py::init<>())
    .def(py::init<Parameters>())
    .def("SearchNodesOverPointsInRadius", [&](SpatialSearchType& self, const NodesContainerType& rStructureNodes, const NodesContainerType& rInputNodes, py::list& rListOfRadius, const DataCommunicator& rDataCommunicator) {
        return self.SearchNodesOverPointsInRadius(rStructureNodes, rInputNodes.begin(), rInputNodes.end(), CopyRadiusArrayToPython(rListOfRadius), rDataCommunicator);
    })
    .def("SearchElementsOverPointsInRadius", [&](SpatialSearchType& self, const ElementsContainerType& rStructureElements, const NodesContainerType& rInputNodes, py::list& rListOfRadius, const DataCommunicator& rDataCommunicator) {
        return self.SearchElementsOverPointsInRadius(rStructureElements, rInputNodes.begin(), rInputNodes.end(), CopyRadiusArrayToPython(rListOfRadius), rDataCommunicator);
    })
    .def("SearchConditionsOverPointsInRadius", [&](SpatialSearchType& self, const ConditionsContainerType& rStructureConditions, const NodesContainerType& rInputNodes, py::list& rListOfRadius, const DataCommunicator& rDataCommunicator) {
        return self.SearchConditionsOverPointsInRadius(rStructureConditions, rInputNodes.begin(), rInputNodes.end(), CopyRadiusArrayToPython(rListOfRadius), rDataCommunicator);
    })
    ;
}

void AddMPISearchStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using NodesContainerType = SpatialSearch::NodesContainerType;
    using ElementsContainerType = SpatialSearch::ElementsContainerType;
    using ConditionsContainerType = SpatialSearch::ConditionsContainerType;

    // The factory of the MPI search strategies
    py::class_<SpecializedSpatialSearchMPIFactory, SpecializedSpatialSearchMPIFactory::Pointer, SpatialSearch>(m, "SpecializedSpatialSearchMPI")
    .def(py::init< >())
    .def(py::init<Parameters>())
    .def("SearchNodesOverPointsInRadius", [&](SpecializedSpatialSearchMPIFactory& self, const NodesContainerType& rStructureNodes, const NodesContainerType& rInputNodes, py::list& rListOfRadius, const DataCommunicator& rDataCommunicator) {
        return self.SearchNodesOverPointsInRadius(rStructureNodes, rInputNodes.begin(), rInputNodes.end(), CopyRadiusArrayToPython(rListOfRadius), rDataCommunicator);
    })
    .def("SearchElementsOverPointsInRadius", [&](SpecializedSpatialSearchMPIFactory& self, const ElementsContainerType& rStructureElements, const NodesContainerType& rInputNodes, py::list& rListOfRadius, const DataCommunicator& rDataCommunicator) {
        return self.SearchElementsOverPointsInRadius(rStructureElements, rInputNodes.begin(), rInputNodes.end(), CopyRadiusArrayToPython(rListOfRadius), rDataCommunicator);
    })
    .def("SearchConditionsOverPointsInRadius", [&](SpecializedSpatialSearchMPIFactory& self, const ConditionsContainerType& rStructureConditions, const NodesContainerType& rInputNodes, py::list& rListOfRadius, const DataCommunicator& rDataCommunicator) {
        return self.SearchConditionsOverPointsInRadius(rStructureConditions, rInputNodes.begin(), rInputNodes.end(), CopyRadiusArrayToPython(rListOfRadius), rDataCommunicator);
    })
    ;

    // Register the specializations
    DefineSpecializedSpatialSearchMPI<SpatialContainer::KDTree>(m, "SpatialSearchKDTreeMPI");
    DefineSpecializedSpatialSearchMPI<SpatialContainer::Octree>(m, "SpatialSearchOctreeMPI");
    DefineSpecializedSpatialSearchMPI<SpatialContainer::BinsStatic>(m, "SpatialSearchBinsStaticMPI");
    DefineSpecializedSpatialSearchMPI<SpatialContainer::BinsDynamic>(m, "SpatialSearchBinsDynamicMPI");
}

}  // namespace Kratos::Python