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
#include "spatial_containers/spatial_search_result_container_vector.h"
#include "spatial_containers/parallel_spatial_search.h"

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
 * @brief Binds a SpatialSearchResult class to Python using pybind11.
 * @tparam TObjectType The type of object stored in the container.
 * @param m The pybind11 module to bind the class to.
 * @param rClassName The name of the class.
 */
template<typename TObjectType>
void BindSpatialSearchResult(pybind11::module& m, const std::string& rClassName) {
    using ResultType = SpatialSearchResult<TObjectType>;
    auto cls = pybind11::class_<ResultType, typename ResultType::Pointer>(m, rClassName.c_str())
    .def(pybind11::init< >())
    .def(pybind11::init<TObjectType*>())
    .def("Reset", &ResultType::Reset)
    .def("Get", [&](ResultType& self) {return self.Get().get();})
    .def("Set", &ResultType::Set)
    .def("GetDistance", &ResultType::GetDistance)
    .def("SetDistance", &ResultType::SetDistance)
    .def("IsObjectFound", &ResultType::IsObjectFound)
    .def("IsDistanceCalculated", &ResultType::IsDistanceCalculated)
    ;
}

/**
 * @brief Binds a SpatialSearchResultContainer class to Python using pybind11.
 * @tparam TObjectType The type of object stored in the container.
 * @tparam TObjectType The type of the object.
 * @tparam TSpatialSearchCommunication The type of spatial search communication considered.
 * @param rClassName The name of the class.
 */
template<class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication = SpatialSearchCommunication::SYNCHRONOUS>
void BindSpatialSearchResultContainer(pybind11::module& m, const std::string& rClassName) {
    using ContainerType = SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>;
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
    .def("GetGlobalIndex", &ContainerType::GetGlobalIndex)
    .def("SetGlobalIndex", &ContainerType::SetGlobalIndex)
    .def("GetLocalIndex", &ContainerType::GetLocalIndex)
    .def("SetLocalIndex", &ContainerType::SetLocalIndex)
    .def("GetLocalResults", &ContainerType::GetLocalResults)
    .def("GetGlobalResults", &ContainerType::GetGlobalResults)
    .def("__getitem__", [](ContainerType& self, const std::size_t Index) {
        return self.GetLocalResults()[Index];
    })
    .def("__call__", [](ContainerType& self, const std::size_t Index) {
        // Check if the communicator has been created
        KRATOS_ERROR_IF_NOT(self.IsSynchronized()) << "The data has not been synchronized" << std::endl;
        return *(self.GetGlobalResults().GetContainer().begin() + Index);
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
 * @brief Binds a SpatialSearchResultContainerVector to a Python module.
 * @tparam TObjectType The type of the object.
 * @tparam TSpatialSearchCommunication The type of spatial search communication considered.
 * @param m The Python module to bind the class to.
 * @param rClassName The name of the class in Python.
 */
template<class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication = SpatialSearchCommunication::SYNCHRONOUS>
void BindSpatialSearchResultContainerVector(pybind11::module& m, const std::string& rClassName) {
    using ContainerVectorType = SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>;
    pybind11::class_<ContainerVectorType, typename ContainerVectorType::Pointer>(m, rClassName.c_str())
    .def(pybind11::init<>())
    .def("NumberOfSearchResults", &ContainerVectorType::NumberOfSearchResults)
    .def("InitializeResult", &ContainerVectorType::InitializeResult)
    .def("HasResult",  &ContainerVectorType::HasResult)
    .def("Clear", &ContainerVectorType::Clear)
    .def("SynchronizeAll", &ContainerVectorType::SynchronizeAll)
    .def("GetDistances", [](ContainerVectorType& self) {
        std::vector<std::vector<double>> results;
        self.GetDistances(results);
        return results;
    })
    .def("GetResultIsLocal", [](
        ContainerVectorType& self,
        const DataCommunicator& rDataCommunicator
        ) {
        std::vector<std::vector<bool>> results;
        self.GetResultIsLocal(results, rDataCommunicator);
        return results;
    })
    .def("GetResultRank", [](ContainerVectorType& self) {
        std::vector<std::vector<int>> results;
        self.GetResultRank(results);
        return results;
    })
    .def("GetResultIsActive", [](
        ContainerVectorType& self,
        const DataCommunicator& rDataCommunicator
        ) {
        std::vector<std::vector<bool>> results;
        self.GetResultIsActive(results, rDataCommunicator);
        return results;
    })
    .def("GetResultIsInside", [](
        ContainerVectorType& self,
        ModelPart::NodesContainerType& rPoints,
        const DataCommunicator& rDataCommunicator,
        const double Tolerance = std::numeric_limits<double>::epsilon()
        ) {
        std::vector<std::vector<bool>> results;
        self.GetResultIsInside(results, rPoints, rDataCommunicator, Tolerance);
        return results;
    })
    .def("GetResultShapeFunctions", [](
        ContainerVectorType& self,
        ModelPart::NodesContainerType& rPoints,
        const DataCommunicator& rDataCommunicator
        ) {
        std::vector<std::vector<Vector>> results;
        self.GetResultShapeFunctions(results, rPoints, rDataCommunicator);
        return results;
    })
    .def("GetResultIndices", [](ContainerVectorType& self) {
        std::vector<std::vector<std::size_t>> results;
        self.GetResultIndices(results);
        return results;
    })
    .def("GetResultNodeIndices", [](ContainerVectorType& self) {
        std::vector<std::vector<std::vector<std::size_t>>> results;
        self.GetResultNodeIndices(results);
        return results;
    })
    .def("GetResultPartitionIndices", [](ContainerVectorType& self) {
        std::vector<std::vector<std::vector<int>>> results;
        self.GetResultPartitionIndices(results);
        return results;
    })
    .def("GetResultCoordinates", [](ContainerVectorType& self) {
        std::vector<std::vector<std::vector<array_1d<double, 3>>>> results;
        self.GetResultCoordinates(results);
        return results;
    })
    .def("GetGlobalPointerCommunicator", &ContainerVectorType::GetGlobalPointerCommunicator)
    .def("GetContainer", &ContainerVectorType::GetContainer)
    .def("__getitem__", [](ContainerVectorType& self, const std::size_t Index) {
        return self[Index];
    })
    .def("__call__", [](ContainerVectorType& self, const std::size_t Index) {
        return self(Index);
    })
    .def("__str__", PrintObject<ContainerVectorType>)
    .def("__iter__", [](ContainerVectorType& self) {
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

/**
 * @brief Defines a search wrapper module in Pybind11.
 * @param m The Pybind11 module to define the search wrapper search in.
 * @param rClassName The name of the search wrapper search class.
 * @tparam TSearchObject The seach object considered
 * @tparam TSpatialSearchCommunication The type of spatial search communication considered.
 */
template<class TSearchObject, SpatialSearchCommunication TSpatialSearchCommunication>
void DefineParallelSpatialSearch(pybind11::module& m, const std::string& rClassName)
{
    using NodesContainerType = ModelPart::NodesContainerType;
    using ElementsContainerType = ModelPart::ElementsContainerType;
    using ConditionsContainerType = ModelPart::ConditionsContainerType;
    using ObjectType = typename TSearchObject::ObjectType;
    using ResultTypeContainerVector = SpatialSearchResultContainerVector<ObjectType, TSpatialSearchCommunication>;
    using ParallelSpatialSearchType = ParallelSpatialSearch<TSearchObject, TSpatialSearchCommunication>;
    using ParallelSpatialSearchPointerType = typename ParallelSpatialSearch<TSearchObject, TSpatialSearchCommunication>::Pointer;

    /// Some constexpr flags
    static constexpr bool IsGeometricalObjectBins = std::is_same_v<TSearchObject, GeometricalObjectsBins>;

    auto parallel_spatial_search = pybind11::class_<ParallelSpatialSearchType, ParallelSpatialSearchPointerType>(m, rClassName.c_str());
    if constexpr (std::is_same<ObjectType, Node>::value) {
        parallel_spatial_search.def(pybind11::init<NodesContainerType&, const DataCommunicator&>());
        parallel_spatial_search.def(pybind11::init<NodesContainerType&, const DataCommunicator&, Parameters>());
    }
    if constexpr (std::is_same<ObjectType, Element>::value || std::is_same<ObjectType, GeometricalObject>::value) {
        parallel_spatial_search.def(pybind11::init<ElementsContainerType&, const DataCommunicator&>());
        parallel_spatial_search.def(pybind11::init<ElementsContainerType&, const DataCommunicator&, Parameters>());
    }
    if constexpr (std::is_same<ObjectType, Condition>::value || std::is_same<ObjectType, GeometricalObject>::value) {
        parallel_spatial_search.def(pybind11::init<ConditionsContainerType&, const DataCommunicator&>());
        parallel_spatial_search.def(pybind11::init<ConditionsContainerType&, const DataCommunicator&, Parameters>());
    }
    parallel_spatial_search.def("GetGlobalBoundingBox", &ParallelSpatialSearchType::GetGlobalBoundingBox);
    parallel_spatial_search.def("SearchInRadius", [&](ParallelSpatialSearchType& self, const NodesContainerType& rNodes, const double Radius) {
        // Perform the search
        auto p_results = Kratos::make_shared<ResultTypeContainerVector>();
        self.SearchInRadius(rNodes.begin(), rNodes.end(), Radius, *p_results);
        return p_results;
    });
    parallel_spatial_search.def("SearchNearestInRadius", [&](ParallelSpatialSearchType& self, const NodesContainerType& rNodes, const double Radius) {
        // Perform the search
        auto p_results = Kratos::make_shared<ResultTypeContainerVector>();
        self.SearchNearestInRadius(rNodes.begin(), rNodes.end(), Radius, *p_results);
        return p_results;
    });
    parallel_spatial_search.def("SearchNearest", [&](ParallelSpatialSearchType& self, const NodesContainerType& rNodes) {
        // Perform the search
        auto p_results = Kratos::make_shared<ResultTypeContainerVector>();
        self.SearchNearest(rNodes.begin(), rNodes.end(), *p_results);
        return p_results;
    });
    if constexpr (IsGeometricalObjectBins) {
        parallel_spatial_search.def("SearchIsInside", [&](ParallelSpatialSearchType& self, const NodesContainerType& rNodes) {
            // Perform the search
            auto p_results = Kratos::make_shared<ResultTypeContainerVector>();
            self.SearchIsInside(rNodes.begin(), rNodes.end(), *p_results);
            return p_results;
        });
    }
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

    // Result types
    BindSpatialSearchResult<Node>(m, "SpatialSearchResultNode");
    BindSpatialSearchResult<GeometricalObject>(m, "SpatialSearchResultGeometricalObject");
    BindSpatialSearchResult<Element>(m, "SpatialSearchResultElement");
    BindSpatialSearchResult<Condition>(m, "SpatialSearchResultCondition");

    // Containers
    BindSpatialSearchResultContainer<Node, SpatialSearchCommunication::SYNCHRONOUS>(m, "SpatialSearchResultContainerNode");
    BindSpatialSearchResultContainer<GeometricalObject, SpatialSearchCommunication::SYNCHRONOUS>(m, "SpatialSearchResultContainerGeometricalObject");
    BindSpatialSearchResultContainer<Element, SpatialSearchCommunication::SYNCHRONOUS>(m, "SpatialSearchResultContainerElement");
    BindSpatialSearchResultContainer<Condition, SpatialSearchCommunication::SYNCHRONOUS>(m, "SpatialSearchResultContainerCondition");

    // Containers vector
    BindSpatialSearchResultContainerVector<Node, SpatialSearchCommunication::SYNCHRONOUS>(m, "SpatialSearchResultContainerVectorNode");
    BindSpatialSearchResultContainerVector<GeometricalObject, SpatialSearchCommunication::SYNCHRONOUS>(m, "SpatialSearchResultContainerVectorGeometricalObject");
    BindSpatialSearchResultContainerVector<Element, SpatialSearchCommunication::SYNCHRONOUS>(m, "SpatialSearchResultContainerVectorElement");
    BindSpatialSearchResultContainerVector<Condition, SpatialSearchCommunication::SYNCHRONOUS>(m, "SpatialSearchResultContainerVectorCondition");

    // GeometricalObjectsBins
    using ResultTypeGeometricalObject = SpatialSearchResult<GeometricalObject>;
    using NodesContainerType = ModelPart::NodesContainerType;
    using ElementsContainerType = ModelPart::ElementsContainerType;
    using ConditionsContainerType = ModelPart::ConditionsContainerType;

    py::class_<GeometricalObjectsBins, GeometricalObjectsBins::Pointer>(m, "GeometricalObjectsBins")
    .def(py::init<ElementsContainerType&>())
    .def(py::init<ConditionsContainerType&>())
    .def(py::init<ElementsContainerType&, double>())
    .def(py::init<ConditionsContainerType&, double>())
    .def("GetBoundingBox", &GeometricalObjectsBins::GetBoundingBox)
    .def("GetCellSizes", &GeometricalObjectsBins::GetCellSizes)
    .def("GetNumberOfCells", &GeometricalObjectsBins::GetNumberOfCells)
    .def("GetTotalNumberOfCells", &GeometricalObjectsBins::GetTotalNumberOfCells)
    .def("SearchInRadius", [&](GeometricalObjectsBins& self, const Point& rPoint, const double Radius) {
        // Perform the search
        std::vector<ResultTypeGeometricalObject> results;
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
        std::vector<std::vector<ResultTypeGeometricalObject>> results;
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
        std::vector<ResultTypeGeometricalObject> results = self.SearchNearestInRadius(rNodes.begin(), rNodes.end(), Radius);

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
        std::vector<ResultTypeGeometricalObject> results = self.SearchNearest(rNodes.begin(), rNodes.end());

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
        std::vector<ResultTypeGeometricalObject> results = self.SearchIsInside(rNodes.begin(), rNodes.end());

        // Copy the results to the python list
        py::list list_results;
        for (auto& r_result : results) {
            list_results.append(r_result);
        }
        return list_results;
    })
    ;

    /* Define the search wrappers */
    // GeometricalObjectsBins
    DefineParallelSpatialSearch<GeometricalObjectsBins, SpatialSearchCommunication::SYNCHRONOUS>(m, "ParallelSpatialSearchGeometricalObjectBins");

    // KDTree
    DefineParallelSpatialSearch<Tree<KDTreePartition<Bucket<3ul, PointObject<Node>, std::vector<PointObject<Node>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS>(m, "ParallelSpatialSearchKDTreeNode");
    DefineParallelSpatialSearch<Tree<KDTreePartition<Bucket<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS>(m, "ParallelSpatialSearchKDTreeElement");
    DefineParallelSpatialSearch<Tree<KDTreePartition<Bucket<3ul, PointObject<Condition>, std::vector<PointObject<Condition>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS>(m, "ParallelSpatialSearchKDTreeCondition");

    // OCTree
    DefineParallelSpatialSearch<Tree<OCTreePartition<Bucket<3ul, PointObject<Node>, std::vector<PointObject<Node>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS>(m, "ParallelSpatialSearchOCTreeNode");
    DefineParallelSpatialSearch<Tree<OCTreePartition<Bucket<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS>(m, "ParallelSpatialSearchOCTreeElement");
    DefineParallelSpatialSearch<Tree<OCTreePartition<Bucket<3ul, PointObject<Condition>, std::vector<PointObject<Condition>::Pointer>>>>, SpatialSearchCommunication::SYNCHRONOUS>(m, "ParallelSpatialSearchOCTreeCondition");

    // StaticBinsTree
    DefineParallelSpatialSearch<Tree<Bins<3ul, PointObject<Node>, std::vector<PointObject<Node>::Pointer>>>, SpatialSearchCommunication::SYNCHRONOUS>(m, "ParallelSpatialSearchStaticBinsTreeNode");
    DefineParallelSpatialSearch<Tree<Bins<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>>, SpatialSearchCommunication::SYNCHRONOUS>(m, "ParallelSpatialSearchStaticBinsTreeElement");
    DefineParallelSpatialSearch<Tree<Bins<3ul, PointObject<Condition>, std::vector<PointObject<Condition>::Pointer>>>, SpatialSearchCommunication::SYNCHRONOUS>(m, "ParallelSpatialSearchStaticBinsTreeCondition");

    // DynamicBins
    DefineParallelSpatialSearch<BinsDynamic<3ul, PointObject<Node>, std::vector<PointObject<Node>::Pointer>>, SpatialSearchCommunication::SYNCHRONOUS>(m, "ParallelSpatialSearchDynamicBinsNode");
    DefineParallelSpatialSearch<BinsDynamic<3ul, PointObject<Element>, std::vector<PointObject<Element>::Pointer>>, SpatialSearchCommunication::SYNCHRONOUS>(m, "ParallelSpatialSearchDynamicBinsElement");
    DefineParallelSpatialSearch<BinsDynamic<3ul, PointObject<Condition>, std::vector<PointObject<Condition>::Pointer>>, SpatialSearchCommunication::SYNCHRONOUS>(m, "ParallelSpatialSearchDynamicBinsCondition");
}

}  // namespace Kratos::Python.