//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "python/add_search_strategies_to_python.h"
#include "utilities/parallel_utilities.h"
#include "spatial_containers/spatial_search.h"
#include "spatial_containers/specialized_spatial_search.h"
#include "spatial_containers/specialized_spatial_search_factory.h"
#include "spatial_containers/geometrical_objects_bins.h"
#include "spatial_containers/spatial_search_result.h"

namespace Kratos::Python
{

/**
 * @brief Generates a list of lists from a vector of vectors
 * @param rList The list to be filled
 * @param rVector The vector of vectors to be copied
 * @tparam TClass The type of the vector of vectors
 */
template<class TClass>
void GenerateListFromVectorOfVector(
    pybind11::list& rList,
    const TClass& rVector
    )
{
    // Clear the result list
    rList.attr("clear")();

    // Copy to the result list
    for (std::size_t i = 0; i < rVector.size(); ++i) {
        pybind11::list results_list_i;
        for (auto& result : rVector[i]) {
            results_list_i.append(result);
        }
        rList.append(results_list_i);
    }
}

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

void AddSearchStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    /// The variables of the spatial search
    using RadiusArrayType = SpatialSearch::RadiusArrayType;
    using VectorResultElementsContainerType = SpatialSearch::VectorResultElementsContainerType;
    using VectorResultNodesContainerType = SpatialSearch::VectorResultNodesContainerType;
    using VectorResultConditionsContainerType = SpatialSearch::VectorResultConditionsContainerType;
    using VectorDistanceType = SpatialSearch::VectorDistanceType;
    using ElementsContainerType = SpatialSearch::ElementsContainerType;
    using NodesContainerType = SpatialSearch::NodesContainerType;
    using ConditionsContainerType = SpatialSearch::ConditionsContainerType;

    py::class_<SpatialSearch, SpatialSearch::Pointer>(m, "SpatialSearch")
    .def(py::init< >())
    .def("SearchElementsInRadiusExclusive", [&](SpatialSearch& self, ModelPart& rModelPart, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        KRATOS_DEBUG_ERROR_IF(size_array != rModelPart.NumberOfElements()) << "The size of the radius array must be equal to the size of the input elements array" << std::endl;

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultElementsContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchElementsInRadiusExclusive(rModelPart, radius_array, results, distances);

        // Copy the results to the python list
        std::tuple<py::list, py::list> results_tuple;
        GenerateListFromVectorOfVector(std::get<0>(results_tuple), results);
        GenerateListFromVectorOfVector(std::get<1>(results_tuple), distances);
        return results_tuple;
    })
    .def("SearchElementsInRadiusExclusive", [&](SpatialSearch& self, ModelPart& rModelPart,
    const ElementsContainerType& rInputElements, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        KRATOS_DEBUG_ERROR_IF(size_array != rInputElements.size()) << "The size of the radius array must be equal to the size of the input elements array" << std::endl;

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultElementsContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchElementsInRadiusExclusive(rModelPart, rInputElements, radius_array, results, distances);

        // Copy the results to the python list
        std::tuple<py::list, py::list> results_tuple;
        GenerateListFromVectorOfVector(std::get<0>(results_tuple), results);
        GenerateListFromVectorOfVector(std::get<1>(results_tuple), distances);
        return results_tuple;
    })
    .def("SearchElementsInRadiusExclusive", [&](SpatialSearch& self,
    const ElementsContainerType& rStructureElements, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        KRATOS_DEBUG_ERROR_IF(size_array != rStructureElements.size()) << "The size of the radius array must be equal to the size of the input elements array" << std::endl;

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultElementsContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchElementsInRadiusExclusive(rStructureElements, radius_array, results, distances);

        // Copy the results to the python list
        std::tuple<py::list, py::list> results_tuple;
        GenerateListFromVectorOfVector(std::get<0>(results_tuple), results);
        GenerateListFromVectorOfVector(std::get<1>(results_tuple), distances);
        return results_tuple;
    })
    .def("SearchElementsInRadiusExclusive", [&](SpatialSearch& self,
    const ElementsContainerType& rStructureElements,
    const ElementsContainerType& rInputElements, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        KRATOS_DEBUG_ERROR_IF(size_array != rInputElements.size()) << "The size of the radius array must be equal to the size of the input elements array" << std::endl;

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultElementsContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchElementsInRadiusExclusive(rStructureElements, rInputElements, radius_array, results, distances);

        // Copy the results to the python list
        std::tuple<py::list, py::list> results_tuple;
        GenerateListFromVectorOfVector(std::get<0>(results_tuple), results);
        GenerateListFromVectorOfVector(std::get<1>(results_tuple), distances);
        return results_tuple;
    })
    .def("SearchElementsInRadiusInclusive", [&](SpatialSearch& self, ModelPart& rModelPart, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        KRATOS_DEBUG_ERROR_IF(size_array != rModelPart.NumberOfElements()) << "The size of the radius array must be equal to the size of the input elements array" << std::endl;

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchElementsInRadiusInclusive(rModelPart, radius_array, results, distances);

        // Copy the results to the python list
        std::tuple<py::list, py::list> results_tuple;
        GenerateListFromVectorOfVector(std::get<0>(results_tuple), results);
        GenerateListFromVectorOfVector(std::get<1>(results_tuple), distances);
        return results_tuple;
    })
    .def("SearchElementsInRadiusInclusive", [&](SpatialSearch& self, ModelPart& rModelPart,
    const ElementsContainerType& rInputElements, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        KRATOS_DEBUG_ERROR_IF(size_array != rInputElements.size()) << "The size of the radius array must be equal to the size of the input elements array" << std::endl;

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchElementsInRadiusInclusive(rModelPart, rInputElements, radius_array, results, distances);

        // Copy the results to the python list
        std::tuple<py::list, py::list> results_tuple;
        GenerateListFromVectorOfVector(std::get<0>(results_tuple), results);
        GenerateListFromVectorOfVector(std::get<1>(results_tuple), distances);
        return results_tuple;
    })
    .def("SearchElementsInRadiusInclusive", [&](SpatialSearch& self,
    const ElementsContainerType& rStructureElements, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        KRATOS_DEBUG_ERROR_IF(size_array != rStructureElements.size()) << "The size of the radius array must be equal to the size of the input elements array" << std::endl;

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchElementsInRadiusInclusive(rStructureElements, radius_array, results, distances);

        // Copy the results to the python list
        std::tuple<py::list, py::list> results_tuple;
        GenerateListFromVectorOfVector(std::get<0>(results_tuple), results);
        GenerateListFromVectorOfVector(std::get<1>(results_tuple), distances);
        return results_tuple;
    })
    .def("SearchElementsInRadiusInclusive", [&](SpatialSearch& self,
    const ElementsContainerType& rStructureElements,
    const ElementsContainerType& rInputElements, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        KRATOS_DEBUG_ERROR_IF(size_array != rInputElements.size()) << "The size of the radius array must be equal to the size of the input elements array" << std::endl;

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchElementsInRadiusInclusive(rStructureElements, rInputElements, radius_array, results, distances);

        // Copy the results to the python list
        std::tuple<py::list, py::list> results_tuple;
        GenerateListFromVectorOfVector(std::get<0>(results_tuple), results);
        GenerateListFromVectorOfVector(std::get<1>(results_tuple), distances);
        return results_tuple;
    })
    .def("SearchNodesInRadiusExclusive", [&](SpatialSearch& self, ModelPart& rModelPart, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        KRATOS_DEBUG_ERROR_IF(size_array != rModelPart.NumberOfNodes()) << "The size of the radius array must be equal to the size of the input nodes array" << std::endl;

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchNodesInRadiusExclusive(rModelPart, radius_array, results, distances);

        // Copy the results to the python list
        std::tuple<py::list, py::list> results_tuple;
        GenerateListFromVectorOfVector(std::get<0>(results_tuple), results);
        GenerateListFromVectorOfVector(std::get<1>(results_tuple), distances);
        return results_tuple;
    })
    .def("SearchNodesInRadiusExclusive", [&](SpatialSearch& self, ModelPart& rModelPart,
    const NodesContainerType& rInputNodes, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        KRATOS_DEBUG_ERROR_IF(size_array != rInputNodes.size()) << "The size of the radius array must be equal to the size of the input nodes array" << std::endl;

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchNodesInRadiusExclusive(rModelPart, rInputNodes, radius_array, results, distances);

        // Copy the results to the python list
        std::tuple<py::list, py::list> results_tuple;
        GenerateListFromVectorOfVector(std::get<0>(results_tuple), results);
        GenerateListFromVectorOfVector(std::get<1>(results_tuple), distances);
        return results_tuple;
    })
    .def("SearchNodesInRadiusExclusive", [&](SpatialSearch& self,
    const NodesContainerType& rStructureNodes, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        KRATOS_DEBUG_ERROR_IF(size_array != rStructureNodes.size()) << "The size of the radius array must be equal to the size of the input nodes array" << std::endl;

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchNodesInRadiusExclusive(rStructureNodes, radius_array, results, distances);

        // Copy the results to the python list
        std::tuple<py::list, py::list> results_tuple;
        GenerateListFromVectorOfVector(std::get<0>(results_tuple), results);
        GenerateListFromVectorOfVector(std::get<1>(results_tuple), distances);
        return results_tuple;
    })
    .def("SearchNodesInRadiusExclusive", [&](SpatialSearch& self,
    const NodesContainerType& rStructureNodes,
    const NodesContainerType& rInputNodes, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        KRATOS_DEBUG_ERROR_IF(size_array != rInputNodes.size()) << "The size of the radius array must be equal to the size of the input nodes array" << std::endl;

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchNodesInRadiusExclusive(rStructureNodes, rInputNodes, radius_array, results, distances);

        // Copy the results to the python list
        std::tuple<py::list, py::list> results_tuple;
        GenerateListFromVectorOfVector(std::get<0>(results_tuple), results);
        GenerateListFromVectorOfVector(std::get<1>(results_tuple), distances);
        return results_tuple;
    })
    .def("SearchNodesInRadiusInclusive", [&](SpatialSearch& self, ModelPart& rModelPart, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        KRATOS_DEBUG_ERROR_IF(size_array != rModelPart.NumberOfNodes()) << "The size of the radius array must be equal to the size of the input nodes array" << std::endl;

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchNodesInRadiusInclusive(rModelPart, radius_array, results, distances);

        // Copy the results to the python list
        std::tuple<py::list, py::list> results_tuple;
        GenerateListFromVectorOfVector(std::get<0>(results_tuple), results);
        GenerateListFromVectorOfVector(std::get<1>(results_tuple), distances);
        return results_tuple;
    })
    .def("SearchNodesInRadiusInclusive", [&](SpatialSearch& self, ModelPart& rModelPart,
    const NodesContainerType& rInputNodes, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        KRATOS_DEBUG_ERROR_IF(size_array != rInputNodes.size()) << "The size of the radius array must be equal to the size of the input nodes array" << std::endl;

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchNodesInRadiusInclusive(rModelPart, rInputNodes, radius_array, results, distances);

        // Copy the results to the python list
        std::tuple<py::list, py::list> results_tuple;
        GenerateListFromVectorOfVector(std::get<0>(results_tuple), results);
        GenerateListFromVectorOfVector(std::get<1>(results_tuple), distances);
        return results_tuple;
    })
    .def("SearchNodesInRadiusInclusive", [&](SpatialSearch& self,
    const NodesContainerType& rStructureNodes, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        KRATOS_DEBUG_ERROR_IF(size_array != rStructureNodes.size()) << "The size of the radius array must be equal to the size of the input nodes array" << std::endl;

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchNodesInRadiusInclusive(rStructureNodes, radius_array, results, distances);

        // Copy the results to the python list
        std::tuple<py::list, py::list> results_tuple;
        GenerateListFromVectorOfVector(std::get<0>(results_tuple), results);
        GenerateListFromVectorOfVector(std::get<1>(results_tuple), distances);
        return results_tuple;
    })
    .def("SearchNodesInRadiusInclusive", [&](SpatialSearch& self,
    const NodesContainerType& rStructureNodes,
    const NodesContainerType& rInputNodes, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        KRATOS_DEBUG_ERROR_IF(size_array != rInputNodes.size()) << "The size of the radius array must be equal to the size of the input nodes array" << std::endl;

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchNodesInRadiusInclusive(rStructureNodes, rInputNodes, radius_array, results, distances);

        // Copy the results to the python list
        std::tuple<py::list, py::list> results_tuple;
        GenerateListFromVectorOfVector(std::get<0>(results_tuple), results);
        GenerateListFromVectorOfVector(std::get<1>(results_tuple), distances);
        return results_tuple;
    })
    .def("SearchConditionsInRadiusExclusive", [&](SpatialSearch& self, ModelPart& rModelPart, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        KRATOS_DEBUG_ERROR_IF(size_array != rModelPart.NumberOfConditions()) << "The size of the radius array must be equal to the size of the input conditions array" << std::endl;

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultConditionsContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchConditionsInRadiusExclusive(rModelPart, radius_array, results, distances);

        // Copy the results to the python list
        std::tuple<py::list, py::list> results_tuple;
        GenerateListFromVectorOfVector(std::get<0>(results_tuple), results);
        GenerateListFromVectorOfVector(std::get<1>(results_tuple), distances);
        return results_tuple;
    })
    .def("SearchConditionsInRadiusExclusive", [&](SpatialSearch& self, ModelPart& rModelPart,
    const ConditionsContainerType& rInputConditions, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        KRATOS_DEBUG_ERROR_IF(size_array != rInputConditions.size()) << "The size of the radius array must be equal to the size of the input conditions array" << std::endl;

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultConditionsContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchConditionsInRadiusExclusive(rModelPart, rInputConditions, radius_array, results, distances);

        // Copy the results to the python list
        std::tuple<py::list, py::list> results_tuple;
        GenerateListFromVectorOfVector(std::get<0>(results_tuple), results);
        GenerateListFromVectorOfVector(std::get<1>(results_tuple), distances);
        return results_tuple;
    })
    .def("SearchConditionsInRadiusExclusive", [&](SpatialSearch& self,
    const ConditionsContainerType& rStructureConditions, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        KRATOS_DEBUG_ERROR_IF(size_array != rStructureConditions.size()) << "The size of the radius array must be equal to the size of the input conditions array" << std::endl;

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultConditionsContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchConditionsInRadiusExclusive(rStructureConditions, radius_array, results, distances);

        // Copy the results to the python list
        std::tuple<py::list, py::list> results_tuple;
        GenerateListFromVectorOfVector(std::get<0>(results_tuple), results);
        GenerateListFromVectorOfVector(std::get<1>(results_tuple), distances);
        return results_tuple;
    })
    .def("SearchConditionsInRadiusExclusive", [&](SpatialSearch& self,
    const ConditionsContainerType& rStructureConditions,
    const ConditionsContainerType& rInputConditions, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        KRATOS_DEBUG_ERROR_IF(size_array != rInputConditions.size()) << "The size of the radius array must be equal to the size of the input conditions array" << std::endl;

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultConditionsContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchConditionsInRadiusExclusive(rStructureConditions, rInputConditions, radius_array, results, distances);

        // Copy the results to the python list
        std::tuple<py::list, py::list> results_tuple;
        GenerateListFromVectorOfVector(std::get<0>(results_tuple), results);
        GenerateListFromVectorOfVector(std::get<1>(results_tuple), distances);
        return results_tuple;
    })
    .def("SearchConditionsInRadiusInclusive", [&](SpatialSearch& self, ModelPart& rModelPart, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        KRATOS_DEBUG_ERROR_IF(size_array != rModelPart.NumberOfConditions()) << "The size of the radius array must be equal to the size of the input conditions array" << std::endl;

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchConditionsInRadiusInclusive(rModelPart, radius_array, results, distances);

        // Copy the results to the python list
        std::tuple<py::list, py::list> results_tuple;
        GenerateListFromVectorOfVector(std::get<0>(results_tuple), results);
        GenerateListFromVectorOfVector(std::get<1>(results_tuple), distances);
        return results_tuple;
    })
    .def("SearchConditionsInRadiusInclusive", [&](SpatialSearch& self, ModelPart& rModelPart,
    const ConditionsContainerType& rInputConditions, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        KRATOS_DEBUG_ERROR_IF(size_array != rInputConditions.size()) << "The size of the radius array must be equal to the size of the input conditions array" << std::endl;

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchConditionsInRadiusInclusive(rModelPart, rInputConditions, radius_array, results, distances);

        // Copy the results to the python list
        std::tuple<py::list, py::list> results_tuple;
        GenerateListFromVectorOfVector(std::get<0>(results_tuple), results);
        GenerateListFromVectorOfVector(std::get<1>(results_tuple), distances);
        return results_tuple;
    })
    .def("SearchConditionsInRadiusInclusive", [&](SpatialSearch& self,
    const ConditionsContainerType& rStructureConditions, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        KRATOS_DEBUG_ERROR_IF(size_array != rStructureConditions.size()) << "The size of the radius array must be equal to the size of the input conditions array" << std::endl;

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchConditionsInRadiusInclusive(rStructureConditions, radius_array, results, distances);

        // Copy the results to the python list
        std::tuple<py::list, py::list> results_tuple;
        GenerateListFromVectorOfVector(std::get<0>(results_tuple), results);
        GenerateListFromVectorOfVector(std::get<1>(results_tuple), distances);
        return results_tuple;
    })
    .def("SearchConditionsInRadiusInclusive", [&](SpatialSearch& self,
    const ConditionsContainerType& rStructureConditions,
    const ConditionsContainerType& rInputConditions, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        KRATOS_DEBUG_ERROR_IF(size_array != rInputConditions.size()) << "The size of the radius array must be equal to the size of the input conditions array" << std::endl;

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchConditionsInRadiusInclusive(rStructureConditions, rInputConditions, radius_array, results, distances);

        // Copy the results to the python list
        std::tuple<py::list, py::list> results_tuple;
        GenerateListFromVectorOfVector(std::get<0>(results_tuple), results);
        GenerateListFromVectorOfVector(std::get<1>(results_tuple), distances);
        return results_tuple;
    })
    ;

    py::class_<SpecializedSpatialSearchFactory, SpecializedSpatialSearchFactory::Pointer, SpatialSearch>(m, "SpecializedSpatialSearch")
    .def(py::init< >())
    .def(py::init<Parameters>())
    ;

    py::class_<SpecializedSpatialSearch<SpatialContainer::KDTree>, SpecializedSpatialSearch<SpatialContainer::KDTree>::Pointer, SpatialSearch>(m, "SpatialSearchKDTree")
    .def(py::init< >())
    .def(py::init<Parameters>())
    ;

    py::class_<SpecializedSpatialSearch<SpatialContainer::Octree>, SpecializedSpatialSearch<SpatialContainer::Octree>::Pointer, SpatialSearch>(m, "SpatialSearchOctree")
    .def(py::init< >())
    .def(py::init<Parameters>())
    ;

    py::class_<SpecializedSpatialSearch<SpatialContainer::BinsStatic>, SpecializedSpatialSearch<SpatialContainer::BinsStatic>::Pointer, SpatialSearch>(m, "SpatialSearchBinsStatic")
    .def(py::init< >())
    .def(py::init<Parameters>())
    ;

    py::class_<SpecializedSpatialSearch<SpatialContainer::BinsDynamic>, SpecializedSpatialSearch<SpatialContainer::BinsDynamic>::Pointer, SpatialSearch>(m, "SpatialSearchBinsDynamic")
    .def(py::init< >())
    .def(py::init<Parameters>())
    ;

    using ResultType = SpatialSearchResult<GeometricalObject>;

    py::class_<ResultType, ResultType::Pointer>(m, "ResultType")
    .def(py::init< >())
    .def(py::init<GeometricalObject*>())
    .def("Reset", &ResultType::Reset)
    .def("Get", [&](ResultType& self) {return self.Get().get();})
    .def("Set", &ResultType::Set)
    .def("GetDistance", &ResultType::GetDistance)
    .def("SetDistance", &ResultType::SetDistance)
    .def("IsObjectFound", &ResultType::IsObjectFound)
    .def("IsDistanceCalculated", &ResultType::IsDistanceCalculated)
    ;

    using NodesContainerType = ModelPart::NodesContainerType;
    using ElementsContainerType = ModelPart::ElementsContainerType;
    using ConditionsContainerType = ModelPart::ConditionsContainerType;

    py::class_<GeometricalObjectsBins, GeometricalObjectsBins::Pointer>(m, "GeometricalObjectsBins")
    .def(py::init<ElementsContainerType&>())
    .def(py::init<ConditionsContainerType&>())
    .def("GetBoundingBox", &GeometricalObjectsBins::GetBoundingBox)
    .def("GetCellSizes", &GeometricalObjectsBins::GetCellSizes)
    .def("GetNumberOfCells", &GeometricalObjectsBins::GetNumberOfCells)
    .def("GetTotalNumberOfCells", &GeometricalObjectsBins::GetTotalNumberOfCells)
    .def("SearchInRadius", [&](GeometricalObjectsBins& self, const Point& rPoint, const double Radius) {
        // Perform the search
        std::vector<ResultType> results;
        self.SearchInRadius(rPoint, Radius, results);

        // Copy the results to the python list
        py::list list_results;
        for (auto& r_result : results) {
            list_results.append(r_result);
        }
        return list_results;
    })
    .def("SearchInRadius", [&](GeometricalObjectsBins& self, const NodesContainerType& rNodes, const double Radius) {
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
    .def("SearchNearestInRadius", [&](GeometricalObjectsBins& self, const Point& rPoint, const double Radius) {
        // Perform the search
        return self.SearchNearestInRadius(rPoint, Radius);
    })
    .def("SearchNearestInRadius", [&](GeometricalObjectsBins& self, const NodesContainerType& rNodes, const double Radius) {
        // Perform the search
        std::vector<ResultType> results = self.SearchNearestInRadius(rNodes.begin(), rNodes.end(), Radius);

        // Copy the results to the python list
        py::list list_results;
        for (auto& r_result : results) {
            list_results.append(r_result);
        }
        return list_results;
    })
    .def("SearchNearest", [&](GeometricalObjectsBins& self, const Point& rPoint) {
        // Perform the search
        return self.SearchNearest(rPoint);
    })
    .def("SearchNearest", [&](GeometricalObjectsBins& self, const NodesContainerType& rNodes) {
        // Perform the search
        std::vector<ResultType> results = self.SearchNearest(rNodes.begin(), rNodes.end());

        // Copy the results to the python list
        py::list list_results;
        for (auto& r_result : results) {
            list_results.append(r_result);
        }
        return list_results;
    })
    .def("SearchIsInside", [&](GeometricalObjectsBins& self, const Point& rPoint) {
        // Perform the search
        return self.SearchIsInside(rPoint);
    })
    .def("SearchIsInside", [&](GeometricalObjectsBins& self, const NodesContainerType& rNodes) {
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

}  // namespace Kratos::Python.