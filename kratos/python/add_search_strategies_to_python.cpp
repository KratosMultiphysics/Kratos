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
    using VectorDistanceType = SpatialSearch::VectorDistanceType;
    using ElementsContainerType = SpatialSearch::ElementsContainerType;
    using NodesContainerType = SpatialSearch::NodesContainerType;

    py::class_<SpatialSearch, SpatialSearch::Pointer>(m, "SpatialSearch")
    .def(py::init< >())
    .def("SearchElementsInRadiusExclusive", [&](SpatialSearch& self, ModelPart& rModelPart, py::list& rListOfRadius, py::list& rResultsList, py::list& rDistancesList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultElementsContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchElementsInRadiusExclusive(rModelPart, radius_array, results, distances);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
        GenerateListFromVectorOfVector(rDistancesList, distances);
    })
    .def("SearchElementsInRadiusExclusive", [&](SpatialSearch& self, ModelPart& rModelPart,
    const ElementsContainerType& rInputElements, py::list& rListOfRadius, py::list& rResultsList, py::list& rDistancesList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultElementsContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchElementsInRadiusExclusive(rModelPart, rInputElements, radius_array, results, distances);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
        GenerateListFromVectorOfVector(rDistancesList, distances);
    })
    .def("SearchElementsInRadiusExclusive", [&](SpatialSearch& self,
    const ElementsContainerType& rStructureElements, py::list& rListOfRadius, py::list& rResultsList, py::list& rDistancesList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultElementsContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchElementsInRadiusExclusive(rStructureElements, radius_array, results, distances);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
        GenerateListFromVectorOfVector(rDistancesList, distances);
    })
    .def("SearchElementsInRadiusExclusive", [&](SpatialSearch& self,
    const ElementsContainerType& rStructureElements,
    const ElementsContainerType& rInputElements, py::list& rListOfRadius, py::list& rResultsList, py::list& rDistancesList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultElementsContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchElementsInRadiusExclusive(rStructureElements, rInputElements, radius_array, results, distances);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
        GenerateListFromVectorOfVector(rDistancesList, distances);
    })
    .def("SearchElementsInRadiusInclusive", [&](SpatialSearch& self, ModelPart& rModelPart, py::list& rListOfRadius, py::list& rResultsList, py::list& rDistancesList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchElementsInRadiusInclusive(rModelPart, radius_array, results, distances);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
        GenerateListFromVectorOfVector(rDistancesList, distances);
    })
    .def("SearchElementsInRadiusInclusive", [&](SpatialSearch& self, ModelPart& rModelPart,
    const ElementsContainerType& rInputElements, py::list& rListOfRadius, py::list& rResultsList, py::list& rDistancesList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchElementsInRadiusInclusive(rModelPart, rInputElements, radius_array, results, distances);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
        GenerateListFromVectorOfVector(rDistancesList, distances);
    })
    .def("SearchElementsInRadiusInclusive", [&](SpatialSearch& self,
    const ElementsContainerType& rStructureElements, py::list& rListOfRadius, py::list& rResultsList, py::list& rDistancesList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchElementsInRadiusInclusive(rStructureElements, radius_array, results, distances);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
        GenerateListFromVectorOfVector(rDistancesList, distances);
    })
    .def("SearchElementsInRadiusInclusive", [&](SpatialSearch& self,
    const ElementsContainerType& rStructureElements,
    const ElementsContainerType& rInputElements, py::list& rListOfRadius, py::list& rResultsList, py::list& rDistancesList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchElementsInRadiusInclusive(rStructureElements, rInputElements, radius_array, results, distances);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
        GenerateListFromVectorOfVector(rDistancesList, distances);
    })
    .def("SearchElementsInRadiusExclusive", [&](SpatialSearch& self, ModelPart& rModelPart, py::list& rListOfRadius, py::list& rResultsList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultElementsContainerType results(size_array);

        // Perform the search
        self.SearchElementsInRadiusExclusive(rModelPart, radius_array, results);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
    })
    .def("SearchElementsInRadiusExclusive", [&](SpatialSearch& self, ModelPart& rModelPart,
    const ElementsContainerType& rInputElements, py::list& rListOfRadius, py::list& rResultsList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultElementsContainerType results(size_array);

        // Perform the search
        self.SearchElementsInRadiusExclusive(rModelPart, rInputElements, radius_array, results);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
    })
    .def("SearchElementsInRadiusExclusive", [&](SpatialSearch& self,
    const ElementsContainerType& rStructureElements, py::list& rListOfRadius, py::list& rResultsList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultElementsContainerType results(size_array);

        // Perform the search
        self.SearchElementsInRadiusExclusive(rStructureElements, radius_array, results);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
    })
    .def("SearchElementsInRadiusExclusive", [&](SpatialSearch& self,
    const ElementsContainerType& rStructureElements,
    const ElementsContainerType& rInputElements, py::list& rListOfRadius, py::list& rResultsList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultElementsContainerType results(size_array);

        // Perform the search
        self.SearchElementsInRadiusExclusive(rStructureElements, rInputElements, radius_array, results);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
    })
    .def("SearchElementsInRadiusInclusive", [&](SpatialSearch& self, ModelPart& rModelPart, py::list& rListOfRadius, py::list& rResultsList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);

        // Perform the search
        self.SearchElementsInRadiusInclusive(rModelPart, radius_array, results);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
    })
    .def("SearchElementsInRadiusInclusive", [&](SpatialSearch& self, ModelPart& rModelPart,
    const ElementsContainerType& rInputElements, py::list& rListOfRadius, py::list& rResultsList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);

        // Perform the search
        self.SearchElementsInRadiusInclusive(rModelPart, rInputElements, radius_array, results);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
    })
    .def("SearchElementsInRadiusInclusive", [&](SpatialSearch& self,
    const ElementsContainerType& rStructureElements, py::list& rListOfRadius, py::list& rResultsList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);

        // Perform the search
        self.SearchElementsInRadiusInclusive(rStructureElements, radius_array, results);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
    })
    .def("SearchElementsInRadiusInclusive", [&](SpatialSearch& self,
    const ElementsContainerType& rStructureElements,
    const ElementsContainerType& rInputElements, py::list& rListOfRadius, py::list& rResultsList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);

        // Perform the search
        self.SearchElementsInRadiusInclusive(rStructureElements, rInputElements, radius_array, results);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
    })
    .def("SearchNodesInRadiusExclusive", [&](SpatialSearch& self, ModelPart& rModelPart, py::list& rListOfRadius, py::list& rResultsList, py::list& rDistancesList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchNodesInRadiusExclusive(rModelPart, radius_array, results, distances);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
        GenerateListFromVectorOfVector(rDistancesList, distances);
    })
    .def("SearchNodesInRadiusExclusive", [&](SpatialSearch& self, ModelPart& rModelPart,
    const NodesContainerType& rInputNodes, py::list& rListOfRadius, py::list& rResultsList, py::list& rDistancesList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchNodesInRadiusExclusive(rModelPart, rInputNodes, radius_array, results, distances);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
        GenerateListFromVectorOfVector(rDistancesList, distances);
    })
    .def("SearchNodesInRadiusExclusive", [&](SpatialSearch& self,
    const NodesContainerType& rStructureNodes, py::list& rListOfRadius, py::list& rResultsList, py::list& rDistancesList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchNodesInRadiusExclusive(rStructureNodes, radius_array, results, distances);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
        GenerateListFromVectorOfVector(rDistancesList, distances);
    })
    .def("SearchNodesInRadiusExclusive", [&](SpatialSearch& self,
    const NodesContainerType& rStructureNodes,
    const NodesContainerType& rInputNodes, py::list& rListOfRadius, py::list& rResultsList, py::list& rDistancesList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchNodesInRadiusExclusive(rStructureNodes, rInputNodes, radius_array, results, distances);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
        GenerateListFromVectorOfVector(rDistancesList, distances);
    })
    .def("SearchNodesInRadiusInclusive", [&](SpatialSearch& self, ModelPart& rModelPart, py::list& rListOfRadius, py::list& rResultsList, py::list& rDistancesList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchNodesInRadiusInclusive(rModelPart, radius_array, results, distances);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
        GenerateListFromVectorOfVector(rDistancesList, distances);
    })
    .def("SearchNodesInRadiusInclusive", [&](SpatialSearch& self, ModelPart& rModelPart,
    const NodesContainerType& rInputNodes, py::list& rListOfRadius, py::list& rResultsList, py::list& rDistancesList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchNodesInRadiusInclusive(rModelPart, rInputNodes, radius_array, results, distances);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
        GenerateListFromVectorOfVector(rDistancesList, distances);
    })
    .def("SearchNodesInRadiusInclusive", [&](SpatialSearch& self,
    const NodesContainerType& rStructureNodes, py::list& rListOfRadius, py::list& rResultsList, py::list& rDistancesList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchNodesInRadiusInclusive(rStructureNodes, radius_array, results, distances);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
        GenerateListFromVectorOfVector(rDistancesList, distances);
    })
    .def("SearchNodesInRadiusInclusive", [&](SpatialSearch& self,
    const NodesContainerType& rStructureNodes,
    const NodesContainerType& rInputNodes, py::list& rListOfRadius, py::list& rResultsList, py::list& rDistancesList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);

        // Perform the search
        self.SearchNodesInRadiusInclusive(rStructureNodes, rInputNodes, radius_array, results, distances);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
        GenerateListFromVectorOfVector(rDistancesList, distances);
    })
    .def("SearchNodesInRadiusExclusive", [&](SpatialSearch& self, ModelPart& rModelPart, py::list& rListOfRadius, py::list& rResultsList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);

        // Perform the search
        self.SearchNodesInRadiusExclusive(rModelPart, radius_array, results);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
    })
    .def("SearchNodesInRadiusExclusive", [&](SpatialSearch& self, ModelPart& rModelPart,
    const NodesContainerType& rInputNodes, py::list& rListOfRadius, py::list& rResultsList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);

        // Perform the search
        self.SearchNodesInRadiusExclusive(rModelPart, rInputNodes, radius_array, results);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
    })
    .def("SearchNodesInRadiusExclusive", [&](SpatialSearch& self,
    const NodesContainerType& rStructureNodes, py::list& rListOfRadius, py::list& rResultsList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);

        // Perform the search
        self.SearchNodesInRadiusExclusive(rStructureNodes, radius_array, results);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
    })
    .def("SearchNodesInRadiusExclusive", [&](SpatialSearch& self,
    const NodesContainerType& rStructureNodes,
    const NodesContainerType& rInputNodes, py::list& rListOfRadius, py::list& rResultsList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);

        // Perform the search
        self.SearchNodesInRadiusExclusive(rStructureNodes, rInputNodes, radius_array, results);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
    })
    .def("SearchNodesInRadiusInclusive", [&](SpatialSearch& self, ModelPart& rModelPart, py::list& rListOfRadius, py::list& rResultsList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);

        // Perform the search
        self.SearchNodesInRadiusInclusive(rModelPart, radius_array, results);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
    })
    .def("SearchNodesInRadiusInclusive", [&](SpatialSearch& self, ModelPart& rModelPart,
    const NodesContainerType& rInputNodes, py::list& rListOfRadius, py::list& rResultsList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);

        // Perform the search
        self.SearchNodesInRadiusInclusive(rModelPart, rInputNodes, radius_array, results);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
    })
    .def("SearchNodesInRadiusInclusive", [&](SpatialSearch& self,
    const NodesContainerType& rStructureNodes, py::list& rListOfRadius, py::list& rResultsList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);

        // Perform the search
        self.SearchNodesInRadiusInclusive(rStructureNodes, radius_array, results);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
    })
    .def("SearchNodesInRadiusInclusive", [&](SpatialSearch& self,
    const NodesContainerType& rStructureNodes,
    const NodesContainerType& rInputNodes, py::list& rListOfRadius, py::list& rResultsList) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();

        // Create the radius array
        RadiusArrayType radius_array(size_array);
        CopyRadiusArrayToPython(rListOfRadius, radius_array);

        // Create the results and distances arrays
        VectorResultNodesContainerType results(size_array);

        // Perform the search
        self.SearchNodesInRadiusInclusive(rStructureNodes, rInputNodes, radius_array, results);

        // Copy the results to the python list
        GenerateListFromVectorOfVector(rResultsList, results);
    })
    ;
}

}  // namespace Kratos::Python.

