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
#include "spatial_containers/spatial_search_result_container.h"

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

/**
 * @brief Binds a SpatialSearchResultContainer class to Python using pybind11.
 * @tparam TObjectType The type of object stored in the container.
 * @param m The pybind11 module to bind the class to.
 * @param rClassName The name of the class.
 */
template<typename TObjectType>
void BindSpatialSearchResultContainer(pybind11::module& m, const std::string& rClassName) {
    using ContainerType = SpatialSearchResultContainer<TObjectType>;
    auto cls = pybind11::class_<ContainerType, typename ContainerType::Pointer>(m, rClassName.c_str())
    .def(pybind11::init<>())
    .def("IsObjectFound", &ContainerType::IsObjectFound)
    .def("NumberOfLocalResults", &ContainerType::NumberOfLocalResults)
    .def("NumberOfGlobalResults", &ContainerType::NumberOfGlobalResults)
    .def("Reserve", &ContainerType::Reserve)
    .def("AddResult", [](ContainerType& self, TObjectType* pObject) {
        self.AddResult(pObject);
    })
    .def("AddResult", [](ContainerType& self, TObjectType* pObject, const double Distance) {
        self.AddResult(pObject, Distance);
    })
    .def("Clear", &ContainerType::Clear)
    .def("SynchronizeAll", &ContainerType::SynchronizeAll)
    .def("GetDistances", &ContainerType::GetDistances)
    .def("GetResultShapeFunctions", &ContainerType::GetResultShapeFunctions)
    .def("GetResultIndices", &ContainerType::GetResultIndices)
    .def("GetResultCoordinates", &ContainerType::GetResultCoordinates)
    .def("__getitem__", [](ContainerType& self, const std::size_t Index) {
        return *(self.GetLocalPointers().GetContainer().begin() + Index);
    })
    .def("__call__", [](ContainerType& self, const std::size_t Index) {
        // Check if the communicator has been created
        KRATOS_ERROR_IF(self.GetGlobalPointerCommunicator() == nullptr) << "The communicator has not been created. Therefore is not synchronized" << std::endl;
        return *(self.GetGlobalPointers().GetContainer().begin() + Index);
    })
    .def("__str__", PrintObject<ContainerType>)
    .def("__iter__", [](ContainerType& self) {
        return pybind11::make_iterator(self.begin(), self.end());
    }, pybind11::keep_alive<0, 1>()); /* Keep object alive while iterator is used */

    // Add the specific methods for the GeometricalObject
    if constexpr (std::is_same<TObjectType, GeometricalObject>::value) {  
        cls.def("AddResult", [](ContainerType& self, SpatialSearchResult<TObjectType>& rObject) {
            self.AddResult(rObject);
        });
    }
}

/**
 * @brief Binds a SpatialSearchResultContainerMap to a Python module.
 * @tparam T The type parameter of the SpatialSearchResultContainerMap.
 * @param m The Python module to bind the class to.
 * @param rClassName The name of the class in Python.
 */
template<typename T>
void BindSpatialSearchResultContainerMap(pybind11::module& m, const std::string& rClassName) {
    using ContainerMapType = SpatialSearchResultContainerMap<T>;
    pybind11::class_<ContainerMapType, typename ContainerMapType::Pointer>(m, rClassName.c_str())
    .def(pybind11::init<>())
    .def("NumberOfSearchResults", &ContainerMapType::NumberOfSearchResults)
    .def("InitializeResult", [](ContainerMapType& self, const std::size_t Index) {
        return self.InitializeResult(Index);
    })
    .def("InitializeResult", [](ContainerMapType& self, const array_1d<double, 3>& rCoordinates) {
        return self.InitializeResult(rCoordinates);
    })
    .def("HasResult", [](ContainerMapType& self, const std::size_t Index) {
        return self.HasResult(Index);
    })
    .def("HasResult", [](ContainerMapType& self, const array_1d<double, 3>& rCoordinates) {
        return self.HasResult(rCoordinates);
    })
    .def("Clear", &ContainerMapType::Clear)
    .def("__getitem__", [](ContainerMapType& self, const std::size_t Index) {
        return self[Index];
    })
    .def("__getitem__", [](ContainerMapType& self, const array_1d<double, 3>& rCoordinates) {
        return self[rCoordinates];
    })
    .def("__call__", [](ContainerMapType& self, const std::size_t Index) {
        return self(Index);
    })
    .def("__call__", [](ContainerMapType& self, const array_1d<double, 3>& rCoordinates) {
        return self(rCoordinates);
    })
    .def("__str__", PrintObject<ContainerMapType>)
    .def("__iter__", [](ContainerMapType& self) {
        return pybind11::make_iterator(self.begin(), self.end());
    }, pybind11::keep_alive<0, 1>()); /* Keep object alive while iterator is used */
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

/**
 * @brief Generates a tuple of two python lists containing the elements and distances of a spatial search result.
 * @param rResults a reference to a vector containing the result elements to be copied to the python list
 * @param rDistances a reference to a vector containing the distances to be copied to the python list
 * @return A tuple containing two python lists, the first one containing the result elements and the second one containing the distances
 * @tparam TClass The type of the class.
 */
template<class TClass>
std::tuple<pybind11::list, pybind11::list> GenerateSpatialSearchSolutionTuple(
    const TClass& rResults, 
    SpatialSearch::VectorDistanceType& rDistances
    )
{
    std::tuple<pybind11::list, pybind11::list> results_tuple;
    GenerateListFromVectorOfVector(std::get<0>(results_tuple), rResults);
    GenerateListFromVectorOfVector(std::get<1>(results_tuple), rDistances);
    return results_tuple;
}

/**
 * @brief Defines a specialized spatial search module in Pybind11.
 * @param m The Pybind11 module to define the specialized spatial search in.
 * @param rClassName The name of the specialized spatial search class.
 */
template <SpatialContainer TSearchBackend>
void DefineSpecializedSpatialSearch(pybind11::module& m, const std::string& rClassName)
{
    using SpatialSearchType = SpecializedSpatialSearch<TSearchBackend>;
    using SpatialSearchPointerType = typename SpecializedSpatialSearch<TSearchBackend>::Pointer;
    using BaseSpatialSearchType = SpatialSearch;

    pybind11::class_<SpatialSearchType, SpatialSearchPointerType, BaseSpatialSearchType>(m, rClassName.c_str())
    .def(pybind11::init<>())
    .def(pybind11::init<Parameters>())
    ;
}

void AddSearchStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    /// The variables of the spatial search
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

        // Perform the search
        VectorResultElementsContainerType results(size_array);
        VectorDistanceType distances(size_array);
        self.SearchElementsInRadiusExclusive(rModelPart, CopyRadiusArrayToPython(rListOfRadius), results, distances);

        // Copy the results to the python list
        return GenerateSpatialSearchSolutionTuple(results, distances);
    })
    .def("SearchElementsInRadiusExclusive", [&](SpatialSearch& self, ModelPart& rModelPart,
    const ElementsContainerType& rInputElements, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();
        KRATOS_DEBUG_ERROR_IF(size_array != rInputElements.size()) << "The size of the radius array must be equal to the size of the input elements array" << std::endl;

        // Perform the search
        VectorResultElementsContainerType results(size_array);
        VectorDistanceType distances(size_array);
        self.SearchElementsInRadiusExclusive(rModelPart, rInputElements, CopyRadiusArrayToPython(rListOfRadius), results, distances);

        // Copy the results to the python list
        return GenerateSpatialSearchSolutionTuple(results, distances);
    })
    .def("SearchElementsInRadiusExclusive", [&](SpatialSearch& self,
    const ElementsContainerType& rStructureElements, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();
        KRATOS_DEBUG_ERROR_IF(size_array != rStructureElements.size()) << "The size of the radius array must be equal to the size of the input elements array" << std::endl;

        // Perform the search
        VectorResultElementsContainerType results(size_array);
        VectorDistanceType distances(size_array);
        self.SearchElementsInRadiusExclusive(rStructureElements, CopyRadiusArrayToPython(rListOfRadius), results, distances);

        // Copy the results to the python list
        return GenerateSpatialSearchSolutionTuple(results, distances);
    })
    .def("SearchElementsInRadiusExclusive", [&](SpatialSearch& self,
    const ElementsContainerType& rStructureElements,
    const ElementsContainerType& rInputElements, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();
        KRATOS_DEBUG_ERROR_IF(size_array != rInputElements.size()) << "The size of the radius array must be equal to the size of the input elements array" << std::endl;

        // Perform the search
        VectorResultElementsContainerType results(size_array);
        VectorDistanceType distances(size_array);
        self.SearchElementsInRadiusExclusive(rStructureElements, rInputElements, CopyRadiusArrayToPython(rListOfRadius), results, distances);

        // Copy the results to the python list
        return GenerateSpatialSearchSolutionTuple(results, distances);
    })
    .def("SearchElementsInRadiusInclusive", [&](SpatialSearch& self, ModelPart& rModelPart, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();
        KRATOS_DEBUG_ERROR_IF(size_array != rModelPart.NumberOfElements()) << "The size of the radius array must be equal to the size of the input elements array" << std::endl;

        // Perform the search
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);
        self.SearchElementsInRadiusInclusive(rModelPart, CopyRadiusArrayToPython(rListOfRadius), results, distances);

        // Copy the results to the python list
        return GenerateSpatialSearchSolutionTuple(results, distances);
    })
    .def("SearchElementsInRadiusInclusive", [&](SpatialSearch& self, ModelPart& rModelPart,
    const ElementsContainerType& rInputElements, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();
        KRATOS_DEBUG_ERROR_IF(size_array != rInputElements.size()) << "The size of the radius array must be equal to the size of the input elements array" << std::endl;

        // Perform the search
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);
        self.SearchElementsInRadiusInclusive(rModelPart, rInputElements, CopyRadiusArrayToPython(rListOfRadius), results, distances);

        // Copy the results to the python list
        return GenerateSpatialSearchSolutionTuple(results, distances);
    })
    .def("SearchElementsInRadiusInclusive", [&](SpatialSearch& self,
    const ElementsContainerType& rStructureElements, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();
        KRATOS_DEBUG_ERROR_IF(size_array != rStructureElements.size()) << "The size of the radius array must be equal to the size of the input elements array" << std::endl;

        // Perform the search
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);
        self.SearchElementsInRadiusInclusive(rStructureElements, CopyRadiusArrayToPython(rListOfRadius), results, distances);

        // Copy the results to the python list
        return GenerateSpatialSearchSolutionTuple(results, distances);
    })
    .def("SearchElementsInRadiusInclusive", [&](SpatialSearch& self,
    const ElementsContainerType& rStructureElements,
    const ElementsContainerType& rInputElements, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();
        KRATOS_DEBUG_ERROR_IF(size_array != rInputElements.size()) << "The size of the radius array must be equal to the size of the input elements array" << std::endl;

        // Perform the search
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);
        self.SearchElementsInRadiusInclusive(rStructureElements, rInputElements, CopyRadiusArrayToPython(rListOfRadius), results, distances);

        // Copy the results to the python list
        return GenerateSpatialSearchSolutionTuple(results, distances);
    })
    .def("SearchNodesInRadiusExclusive", [&](SpatialSearch& self, ModelPart& rModelPart, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();
        KRATOS_DEBUG_ERROR_IF(size_array != rModelPart.NumberOfNodes()) << "The size of the radius array must be equal to the size of the input nodes array" << std::endl;

        // Perform the search
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);
        self.SearchNodesInRadiusExclusive(rModelPart, CopyRadiusArrayToPython(rListOfRadius), results, distances);

        // Copy the results to the python list
        return GenerateSpatialSearchSolutionTuple(results, distances);
    })
    .def("SearchNodesInRadiusExclusive", [&](SpatialSearch& self, ModelPart& rModelPart,
    const NodesContainerType& rInputNodes, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();
        KRATOS_DEBUG_ERROR_IF(size_array != rInputNodes.size()) << "The size of the radius array must be equal to the size of the input nodes array" << std::endl;

        // Perform the search
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);
        self.SearchNodesInRadiusExclusive(rModelPart, rInputNodes, CopyRadiusArrayToPython(rListOfRadius), results, distances);

        // Copy the results to the python list
        return GenerateSpatialSearchSolutionTuple(results, distances);
    })
    .def("SearchNodesInRadiusExclusive", [&](SpatialSearch& self,
    const NodesContainerType& rStructureNodes, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();
        KRATOS_DEBUG_ERROR_IF(size_array != rStructureNodes.size()) << "The size of the radius array must be equal to the size of the input nodes array" << std::endl;

        // Perform the search
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);
        self.SearchNodesInRadiusExclusive(rStructureNodes, CopyRadiusArrayToPython(rListOfRadius), results, distances);

        // Copy the results to the python list
        return GenerateSpatialSearchSolutionTuple(results, distances);
    })
    .def("SearchNodesInRadiusExclusive", [&](SpatialSearch& self,
    const NodesContainerType& rStructureNodes,
    const NodesContainerType& rInputNodes, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();
        KRATOS_DEBUG_ERROR_IF(size_array != rInputNodes.size()) << "The size of the radius array must be equal to the size of the input nodes array" << std::endl;

        // Perform the search
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);
        self.SearchNodesInRadiusExclusive(rStructureNodes, rInputNodes, CopyRadiusArrayToPython(rListOfRadius), results, distances);

        // Copy the results to the python list
        return GenerateSpatialSearchSolutionTuple(results, distances);
    })
    .def("SearchNodesInRadiusInclusive", [&](SpatialSearch& self, ModelPart& rModelPart, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();
        KRATOS_DEBUG_ERROR_IF(size_array != rModelPart.NumberOfNodes()) << "The size of the radius array must be equal to the size of the input nodes array" << std::endl;

        // Perform the search
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);
        self.SearchNodesInRadiusInclusive(rModelPart, CopyRadiusArrayToPython(rListOfRadius), results, distances);

        // Copy the results to the python list
        return GenerateSpatialSearchSolutionTuple(results, distances);
    })
    .def("SearchNodesInRadiusInclusive", [&](SpatialSearch& self, ModelPart& rModelPart,
    const NodesContainerType& rInputNodes, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();
        KRATOS_DEBUG_ERROR_IF(size_array != rInputNodes.size()) << "The size of the radius array must be equal to the size of the input nodes array" << std::endl;

        // Perform the search
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);
        self.SearchNodesInRadiusInclusive(rModelPart, rInputNodes, CopyRadiusArrayToPython(rListOfRadius), results, distances);

        // Copy the results to the python list
        return GenerateSpatialSearchSolutionTuple(results, distances);
    })
    .def("SearchNodesInRadiusInclusive", [&](SpatialSearch& self,
    const NodesContainerType& rStructureNodes, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();
        KRATOS_DEBUG_ERROR_IF(size_array != rStructureNodes.size()) << "The size of the radius array must be equal to the size of the input nodes array" << std::endl;

        // Perform the search
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);
        self.SearchNodesInRadiusInclusive(rStructureNodes, CopyRadiusArrayToPython(rListOfRadius), results, distances);

        // Copy the results to the python list
        return GenerateSpatialSearchSolutionTuple(results, distances);
    })
    .def("SearchNodesInRadiusInclusive", [&](SpatialSearch& self,
    const NodesContainerType& rStructureNodes,
    const NodesContainerType& rInputNodes, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();
        KRATOS_DEBUG_ERROR_IF(size_array != rInputNodes.size()) << "The size of the radius array must be equal to the size of the input nodes array" << std::endl;

        // Perform the search
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);
        self.SearchNodesInRadiusInclusive(rStructureNodes, rInputNodes, CopyRadiusArrayToPython(rListOfRadius), results, distances);

        // Copy the results to the python list
        return GenerateSpatialSearchSolutionTuple(results, distances);
    })
    .def("SearchConditionsInRadiusExclusive", [&](SpatialSearch& self, ModelPart& rModelPart, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();
        KRATOS_DEBUG_ERROR_IF(size_array != rModelPart.NumberOfConditions()) << "The size of the radius array must be equal to the size of the input conditions array" << std::endl;

        // Perform the search
        VectorResultConditionsContainerType results(size_array);
        VectorDistanceType distances(size_array);
        self.SearchConditionsInRadiusExclusive(rModelPart, CopyRadiusArrayToPython(rListOfRadius), results, distances);

        // Copy the results to the python list
        return GenerateSpatialSearchSolutionTuple(results, distances);
    })
    .def("SearchConditionsInRadiusExclusive", [&](SpatialSearch& self, ModelPart& rModelPart,
    const ConditionsContainerType& rInputConditions, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();
        KRATOS_DEBUG_ERROR_IF(size_array != rInputConditions.size()) << "The size of the radius array must be equal to the size of the input conditions array" << std::endl;

        // Perform the search
        VectorResultConditionsContainerType results(size_array);
        VectorDistanceType distances(size_array);
        self.SearchConditionsInRadiusExclusive(rModelPart, rInputConditions, CopyRadiusArrayToPython(rListOfRadius), results, distances);

        // Copy the results to the python list
        return GenerateSpatialSearchSolutionTuple(results, distances);
    })
    .def("SearchConditionsInRadiusExclusive", [&](SpatialSearch& self,
    const ConditionsContainerType& rStructureConditions, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();
        KRATOS_DEBUG_ERROR_IF(size_array != rStructureConditions.size()) << "The size of the radius array must be equal to the size of the input conditions array" << std::endl;

        // Perform the search
        VectorResultConditionsContainerType results(size_array);
        VectorDistanceType distances(size_array);
        self.SearchConditionsInRadiusExclusive(rStructureConditions, CopyRadiusArrayToPython(rListOfRadius), results, distances);

        // Copy the results to the python list
        return GenerateSpatialSearchSolutionTuple(results, distances);
    })
    .def("SearchConditionsInRadiusExclusive", [&](SpatialSearch& self,
    const ConditionsContainerType& rStructureConditions,
    const ConditionsContainerType& rInputConditions, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();
        KRATOS_DEBUG_ERROR_IF(size_array != rInputConditions.size()) << "The size of the radius array must be equal to the size of the input conditions array" << std::endl;

        // Perform the search
        VectorResultConditionsContainerType results(size_array);
        VectorDistanceType distances(size_array);
        self.SearchConditionsInRadiusExclusive(rStructureConditions, rInputConditions, CopyRadiusArrayToPython(rListOfRadius), results, distances);

        // Copy the results to the python list
        return GenerateSpatialSearchSolutionTuple(results, distances);
    })
    .def("SearchConditionsInRadiusInclusive", [&](SpatialSearch& self, ModelPart& rModelPart, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();
        KRATOS_DEBUG_ERROR_IF(size_array != rModelPart.NumberOfConditions()) << "The size of the radius array must be equal to the size of the input conditions array" << std::endl;

        // Perform the search
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);
        self.SearchConditionsInRadiusInclusive(rModelPart, CopyRadiusArrayToPython(rListOfRadius), results, distances);

        // Copy the results to the python list
        return GenerateSpatialSearchSolutionTuple(results, distances);
    })
    .def("SearchConditionsInRadiusInclusive", [&](SpatialSearch& self, ModelPart& rModelPart,
    const ConditionsContainerType& rInputConditions, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();
        KRATOS_DEBUG_ERROR_IF(size_array != rInputConditions.size()) << "The size of the radius array must be equal to the size of the input conditions array" << std::endl;

        // Perform the search
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);
        self.SearchConditionsInRadiusInclusive(rModelPart, rInputConditions, CopyRadiusArrayToPython(rListOfRadius), results, distances);

        // Copy the results to the python list
        return GenerateSpatialSearchSolutionTuple(results, distances);
    })
    .def("SearchConditionsInRadiusInclusive", [&](SpatialSearch& self,
    const ConditionsContainerType& rStructureConditions, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();
        KRATOS_DEBUG_ERROR_IF(size_array != rStructureConditions.size()) << "The size of the radius array must be equal to the size of the input conditions array" << std::endl;

        // Perform the search
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);
        self.SearchConditionsInRadiusInclusive(rStructureConditions, CopyRadiusArrayToPython(rListOfRadius), results, distances);

        // Copy the results to the python list
        return GenerateSpatialSearchSolutionTuple(results, distances);
    })
    .def("SearchConditionsInRadiusInclusive", [&](SpatialSearch& self,
    const ConditionsContainerType& rStructureConditions,
    const ConditionsContainerType& rInputConditions, py::list& rListOfRadius) {
        // Get the size of the radius array
        const std::size_t size_array = rListOfRadius.size();
        KRATOS_DEBUG_ERROR_IF(size_array != rInputConditions.size()) << "The size of the radius array must be equal to the size of the input conditions array" << std::endl;

        // Perform the search
        VectorResultNodesContainerType results(size_array);
        VectorDistanceType distances(size_array);
        self.SearchConditionsInRadiusInclusive(rStructureConditions, rInputConditions, CopyRadiusArrayToPython(rListOfRadius), results, distances);

        // Copy the results to the python list
        return GenerateSpatialSearchSolutionTuple(results, distances);
    })
    ;

    // The factory of the search strategies
    py::class_<SpecializedSpatialSearchFactory, SpecializedSpatialSearchFactory::Pointer, SpatialSearch>(m, "SpecializedSpatialSearch")
    .def(py::init< >())
    .def(py::init<Parameters>())
    ;

    // Register the specializations
    DefineSpecializedSpatialSearch<SpatialContainer::KDTree>(m, "SpatialSearchKDTree");
    DefineSpecializedSpatialSearch<SpatialContainer::Octree>(m, "SpatialSearchOctree");
    DefineSpecializedSpatialSearch<SpatialContainer::BinsStatic>(m, "SpatialSearchBinsStatic");
    DefineSpecializedSpatialSearch<SpatialContainer::BinsDynamic>(m, "SpatialSearchBinsDynamic");

    using ResultTypeGeometricalObject = SpatialSearchResult<GeometricalObject>;

    py::class_<ResultTypeGeometricalObject, ResultTypeGeometricalObject::Pointer>(m, "ResultTypeGeometricalObject")
    .def(py::init< >())
    .def(py::init<GeometricalObject*>())
    .def("Reset", &ResultTypeGeometricalObject::Reset)
    .def("Get", [&](ResultTypeGeometricalObject& self) {return self.Get().get();})
    .def("Set", &ResultTypeGeometricalObject::Set)
    .def("GetDistance", &ResultTypeGeometricalObject::GetDistance)
    .def("SetDistance", &ResultTypeGeometricalObject::SetDistance)
    .def("IsObjectFound", &ResultTypeGeometricalObject::IsObjectFound)
    .def("IsDistanceCalculated", &ResultTypeGeometricalObject::IsDistanceCalculated)
    ;

    // Containers
    BindSpatialSearchResultContainer<Node>(m, "ResultTypeContainerNode");
    BindSpatialSearchResultContainer<GeometricalObject>(m, "ResultTypeContainerGeometricalObject");

    // Containers map
    BindSpatialSearchResultContainerMap<Node>(m, "ResultTypeContainerMapNode");
    BindSpatialSearchResultContainerMap<GeometricalObject>(m, "ResultTypeContainerMapGeometricalObject");

    using ResultTypeContainerGeometricalObject = SpatialSearchResultContainer<GeometricalObject>;
    using ResultTypeContainerMapGeometricalObject = SpatialSearchResultContainerMap<GeometricalObject>;
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
        ResultTypeContainerGeometricalObject results;
        self.SearchInRadius(rPoint, Radius, results);
        return results;
    })
    .def("SearchInRadius", [&](GeometricalObjectsBins& self, const NodesContainerType& rNodes, const double Radius) {
        // Perform the search
        ResultTypeContainerMapGeometricalObject results; 
        self.SearchInRadius(rNodes.begin(), rNodes.end(), Radius, results);
        return results;
    })
    .def("SearchNearestInRadius", [&](GeometricalObjectsBins& self, const Point& rPoint, const double Radius) {
        // Perform the search
        ResultTypeContainerGeometricalObject results;
        self.SearchNearestInRadius(rPoint, Radius, results);
        return results;
    })
    .def("SearchNearestInRadius", [&](GeometricalObjectsBins& self, const NodesContainerType& rNodes, const double Radius) {
        // Perform the search
        ResultTypeContainerMapGeometricalObject results;
        self.SearchNearestInRadius(rNodes.begin(), rNodes.end(), Radius, results);
        return results;
    })
    .def("SearchNearest", [&](GeometricalObjectsBins& self, const Point& rPoint) {
        // Perform the search
        ResultTypeContainerGeometricalObject results;
        self.SearchNearest(rPoint, results);
        return results;
    })
    .def("SearchNearest", [&](GeometricalObjectsBins& self, const NodesContainerType& rNodes) {
        // Perform the search
        ResultTypeContainerMapGeometricalObject results;
        self.SearchNearest(rNodes.begin(), rNodes.end(), results);
        return results;
    })
    .def("SearchIsInside", [&](GeometricalObjectsBins& self, const Point& rPoint) {
        // Perform the search
        ResultTypeContainerGeometricalObject results;
        self.SearchIsInside(rPoint, results);
        return results;
    })
    .def("SearchIsInside", [&](GeometricalObjectsBins& self, const NodesContainerType& rNodes) {
        // Perform the search
        ResultTypeContainerMapGeometricalObject results;
        self.SearchIsInside(rNodes.begin(), rNodes.end(), results);
        return results;
    })
    ;
}

}  // namespace Kratos::Python.