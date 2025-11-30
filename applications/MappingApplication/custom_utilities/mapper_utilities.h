//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#pragma once

// System includes
#include <array>
#include <vector>

// External includes

// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/math_utils.h"
#include "mapping_application_variables.h"
#include "mappers/mapper_flags.h"
#include "custom_utilities/mapper_local_system.h"

namespace Kratos {
namespace MapperUtilities {

typedef std::size_t SizeType;
typedef std::size_t IndexType;

typedef Node NodeType;

typedef Kratos::unique_ptr<MapperInterfaceInfo> MapperInterfaceInfoUniquePointerType;

typedef Kratos::shared_ptr<MapperInterfaceInfo> MapperInterfaceInfoPointerType;
typedef std::vector<std::vector<MapperInterfaceInfoPointerType>> MapperInterfaceInfoPointerVectorType;

typedef Kratos::unique_ptr<MapperLocalSystem> MapperLocalSystemPointer;
typedef std::vector<MapperLocalSystemPointer> MapperLocalSystemPointerVector;
typedef Kratos::shared_ptr<MapperLocalSystemPointerVector> MapperLocalSystemPointerVectorPointer;

using BoundingBoxType = std::array<double, 6>;

static void FillFunctionNodes(const NodeType& rNode,
                         const Variable<double>& rVariable,
                         double& rValue)
{
    rValue = rNode.FastGetSolutionStepValue(rVariable);
}

static void FillFunctionNodesNonHist(const NodeType& rNode,
                                const Variable<double>& rVariable,
                                double& rValue)
{
    rValue = rNode.GetValue(rVariable);
}

static inline std::function<void(const NodeType&, const Variable<double>&, double&)>
GetFillFunctionForNodes(const Kratos::Flags& rMappingOptions)
{
    if (rMappingOptions.Is(MapperFlags::FROM_NON_HISTORICAL))
        return &FillFunctionNodesNonHist;
    return &FillFunctionNodes;
}

static void UpdateFunctionNodes(NodeType& rNode,
                           const Variable<double>& rVariable,
                           const double Value,
                           const double Factor)
{
    rNode.FastGetSolutionStepValue(rVariable) = Value * Factor;
}

static void UpdateFunctionNodesWithAdd(NodeType& rNode,
                            const Variable<double>& rVariable,
                            const double Value,
                            const double Factor)
{
    rNode.FastGetSolutionStepValue(rVariable) += Value * Factor;
}

static void UpdateFunctionNodesNonHist(NodeType& rNode,
                            const Variable<double>& rVariable,
                            const double Value,
                            const double Factor)
{
    rNode.SetValue(rVariable, Value * Factor);
}

static void UpdateFunctionNodesNonHistWithAdd(NodeType& rNode,
                            const Variable<double>& rVariable,
                            const double Value,
                            const double Factor)
{
    rNode.GetValue(rVariable) += Value * Factor;
}

static inline std::function<void(NodeType&, const Variable<double>&, const double, const double)>
GetUpdateFunctionForNodes(const Kratos::Flags& rMappingOptions)
{
    if (rMappingOptions.Is(MapperFlags::ADD_VALUES) && rMappingOptions.Is(MapperFlags::TO_NON_HISTORICAL))
        return &UpdateFunctionNodesNonHistWithAdd;
    if (rMappingOptions.Is(MapperFlags::ADD_VALUES))
        return &UpdateFunctionNodesWithAdd;
    if (rMappingOptions.Is(MapperFlags::TO_NON_HISTORICAL))
        return &UpdateFunctionNodesNonHist;
    return &UpdateFunctionNodes;
}

template<class TEntityType>
static inline std::function<void(TEntityType&, const Variable<double>&, double, double)>
GetUpdateFunctionForEntities(const Kratos::Flags& rMappingOptions)
{
    if (rMappingOptions.Is(MapperFlags::ADD_VALUES))
        return [](TEntityType& rEntity, const Variable<double>& rVariable, double value, double factor){
            rEntity.SetValue(rVariable, rEntity.GetValue(rVariable) + factor * value);
        };

    return [](TEntityType& rEntity, const Variable<double>& rVariable, double value, double factor){
        rEntity.SetValue(rVariable, factor * value);
    };
}

template<class TVectorType, bool TParallel=true>
void UpdateSystemVectorFromModelPartNodes(
    TVectorType& rVector,
    const ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const Kratos::Flags& rMappingOptions,
    const bool InParallel=true)
{
    KRATOS_TRY;

    if (!rModelPart.GetCommunicator().GetDataCommunicator().IsDefinedOnThisRank()) return;

    // Here we construct a function pointer to not have the if all the time inside the loop
    const auto fill_fct = MapperUtilities::GetFillFunctionForNodes(rMappingOptions);

    const int num_local_nodes = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    const auto nodes_begin = rModelPart.GetCommunicator().LocalMesh().NodesBegin();

    // necessary bcs the Trilinos Vector is not threadsafe in the default configuration
    const int num_threads = InParallel ? ParallelUtilities::GetNumThreads() : 1;

    KRATOS_ERROR_IF(!rMappingOptions.Is(MapperFlags::FROM_NON_HISTORICAL) && !rModelPart.HasNodalSolutionStepVariable(rVariable)) << "Solution step variable \"" << rVariable.Name() << "\" missing in ModelPart \"" << rModelPart.FullName() << "\"!" << std::endl;

    IndexPartition<std::size_t>(num_local_nodes, num_threads).for_each([&](const std::size_t i){
        fill_fct(*(nodes_begin + i), rVariable, rVector[i]);
    });

    KRATOS_CATCH("");
}

template<class TVectorType, class TContainerType>
void UpdateSystemVectorFromGeometricalEntities(
    TVectorType& rVector,
    const TContainerType& rContainer,
    const Variable<double>& rVariable,
    const Kratos::Flags& rMappingOptions,
    const bool InParallel = true)
{
    KRATOS_TRY;

    using EntityType = typename TContainerType::value_type;

    const std::size_t n_entities = rContainer.size();
    auto it_begin = rContainer.begin();

    const int num_threads = InParallel ? ParallelUtilities::GetNumThreads() : 1;

    IndexPartition<std::size_t>(n_entities, num_threads).for_each([&](const std::size_t i){
            const EntityType& r_entity = *(it_begin + i);
            const double value = r_entity.GetValue(rVariable);
            rVector[i] = value;  // <-- always assign, no add
        });

    KRATOS_CATCH("");
}

template<class TVectorType>
void UpdateSystemVectorFromModelPartGeometries(
    TVectorType& rVector,
    const ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const Kratos::Flags& rMappingOptions,
    const bool InParallel = true)
{
    KRATOS_TRY;

    if (!rModelPart.GetCommunicator().GetDataCommunicator().IsDefinedOnThisRank())
        return;

    const auto& r_geometries = rModelPart.Geometries();

    UpdateSystemVectorFromGeometricalEntities(
        rVector,
        r_geometries,
        rVariable,
        rMappingOptions,
        InParallel);

    KRATOS_CATCH("");
}

template<class TVectorType>
void UpdateSystemVectorFromModelPartConditions(
    TVectorType& rVector,
    const ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const Kratos::Flags& rMappingOptions,
    const bool InParallel = true)
{
    KRATOS_TRY;

    if (!rModelPart.GetCommunicator().GetDataCommunicator().IsDefinedOnThisRank())
        return;

    const auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();

    UpdateSystemVectorFromGeometricalEntities(
        rVector,
        r_local_mesh.Conditions(),
        rVariable,
        rMappingOptions,
        InParallel);

    KRATOS_CATCH("");
}

template<class TVectorType>
void UpdateSystemVectorFromModelPartElements(
    TVectorType& rVector,
    const ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const Kratos::Flags& rMappingOptions,
    const bool InParallel = true)
{
    KRATOS_TRY;

    if (!rModelPart.GetCommunicator().GetDataCommunicator().IsDefinedOnThisRank())
        return;

    const auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();

    UpdateSystemVectorFromGeometricalEntities(
        rVector,
        r_local_mesh.Elements(),
        rVariable,
        rMappingOptions,
        InParallel);

    KRATOS_CATCH("");
}

template<class TVectorType>
void UpdateModelPartNodesFromSystemVector(
    const TVectorType& rVector,
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const Kratos::Flags& rMappingOptions,
    const bool InParallel=true)
{
    KRATOS_TRY;

    if (!rModelPart.GetCommunicator().GetDataCommunicator().IsDefinedOnThisRank()) return;

    const double factor = rMappingOptions.Is(MapperFlags::SWAP_SIGN) ? -1.0 : 1.0;

    // Here we construct a function pointer to not have the if all the time inside the loop
    const auto update_fct = std::bind(MapperUtilities::GetUpdateFunctionForNodes(rMappingOptions),
                                        std::placeholders::_1,
                                        std::placeholders::_2,
                                        std::placeholders::_3,
                                        factor);
    const int num_local_nodes = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    const auto nodes_begin = rModelPart.GetCommunicator().LocalMesh().NodesBegin();

    // necessary bcs the Trilinos Vector is not threadsafe in the default configuration
    const int num_threads = InParallel ? ParallelUtilities::GetNumThreads() : 1;

    KRATOS_ERROR_IF(!rMappingOptions.Is(MapperFlags::TO_NON_HISTORICAL) && !rModelPart.HasNodalSolutionStepVariable(rVariable)) << "Solution step variable \"" << rVariable.Name() << "\" missing in ModelPart \"" << rModelPart.FullName() << "\"!" << std::endl;

    IndexPartition<std::size_t>(num_local_nodes, num_threads).for_each([&](const std::size_t i){
        update_fct(*(nodes_begin + i), rVariable, rVector[i]);
    });

    if (rMappingOptions.Is(MapperFlags::TO_NON_HISTORICAL)) {
        rModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(rVariable);
    } else {
        rModelPart.GetCommunicator().SynchronizeVariable(rVariable);
    }

    KRATOS_CATCH("");
}

template<class TVectorType, class TContainerType>
void UpdateGeometricalEntitiesFromSystemVector(
    const TVectorType& rVector,
    TContainerType& rContainer,
    const Variable<double>& rVariable,
    const Kratos::Flags& rMappingOptions,
    const bool InParallel=true)
{
    KRATOS_TRY;

    using EntityType = typename TContainerType::value_type;

    const double factor = rMappingOptions.Is(MapperFlags::SWAP_SIGN) ? -1.0 : 1.0;

    const auto update_fct = MapperUtilities::GetUpdateFunctionForEntities<EntityType>(rMappingOptions);

    const std::size_t n_entities = rContainer.size();
    auto it_begin = rContainer.begin();

    const int num_threads = InParallel ? ParallelUtilities::GetNumThreads() : 1;

    IndexPartition<std::size_t>(n_entities, num_threads).for_each(
        [&](const std::size_t i){
            EntityType& r_obj = *(it_begin + i);
            update_fct(r_obj, rVariable, rVector[i], factor);
        });

    KRATOS_CATCH("");
}

template<class TVectorType>
void UpdateModelPartConditionsFromSystemVector(
    const TVectorType& rVector,
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const Kratos::Flags& rMappingOptions,
    const bool InParallel=true)
{
    if (!rModelPart.GetCommunicator().GetDataCommunicator().IsDefinedOnThisRank()) return;
    auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();

    UpdateGeometricalEntitiesFromSystemVector(
        rVector,
        r_local_mesh.Conditions(),
        rVariable,
        rMappingOptions,
        InParallel);
}

template<class TVectorType>
void UpdateModelPartElementsFromSystemVector(
    const TVectorType& rVector,
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const Kratos::Flags& rMappingOptions,
    const bool InParallel=true)
{
    if (!rModelPart.GetCommunicator().GetDataCommunicator().IsDefinedOnThisRank()) return;
    auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();

    UpdateGeometricalEntitiesFromSystemVector(
        rVector,
        r_local_mesh.Elements(),
        rVariable,
        rMappingOptions,
        InParallel);
}

template<class TVectorType>
void UpdateModelPartGeometriesFromSystemVector(
    const TVectorType& rVector,
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const Kratos::Flags& rMappingOptions,
    const bool InParallel = true)
{
    KRATOS_TRY;

    auto& r_geometries = rModelPart.Geometries();

    UpdateGeometricalEntitiesFromSystemVector(
        rVector,
        r_geometries,
        rVariable,
        rMappingOptions,
        InParallel);

    KRATOS_CATCH("");
}


/**
* @brief Assigning INTERFACE_EQUATION_IDs to the nodes, with and without MPI
* This function assigns the INTERFACE_EQUATION_IDs to the nodes, which
* act as EquationIds for the MappingMatrix. This work with and without MPI,
* in MPI a ScanSum is performed with the local number of nodes
* @param rModelPartCommunicator The Modelpart-Communicator to be used
* @author Philipp Bucher
*/
void KRATOS_API(MAPPING_APPLICATION) AssignInterfaceEquationIds(Communicator& rModelPartCommunicator);

void KRATOS_API(MAPPING_APPLICATION) CreateMapperLocalSystemsFromNodes(const MapperLocalSystem& rMapperLocalSystemPrototype,
                                       const Communicator& rModelPartCommunicator,
                                       std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rLocalSystems);

void CreateMapperLocalSystemsFromGeometries(const MapperLocalSystem& rMapperLocalSystemPrototype,
                                            const Communicator& rModelPartCommunicator,
                                            std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rLocalSystems);

template <class T1, class T2>
inline double ComputeDistance(const T1& rCoords1,
                              const T2& rCoords2)
{
    return std::sqrt( std::pow(rCoords1[0] - rCoords2[0] , 2) +
                      std::pow(rCoords1[1] - rCoords2[1] , 2) +
                      std::pow(rCoords1[2] - rCoords2[2] , 2) );
}

template <class T1, class T2, class T3>
bool PointsAreCollinear(
    const T1& rP1,
    const T2& rP2,
    const T3& rP3)
{
    // TODO this can probably be optimized
    const double a = MathUtils<double>::Norm3(rP1-rP2);
    const double b = MathUtils<double>::Norm3(rP2-rP3);
    const double c = MathUtils<double>::Norm3(rP3-rP1);

    const double s = (a+b+c) / 2.0;

    return (std::sqrt(s*(s-a)*(s-b)*(s-c))) < 1e-12;
}

template <typename TContainer>
double ComputeMaxEdgeLengthLocal(const TContainer& rEntityContainer);

double ComputeSearchRadius(const ModelPart& rModelPart, int EchoLevel);

double ComputeSearchRadius(const ModelPart& rModelPart1, const ModelPart& rModelPart2, const int EchoLevel);

void CheckInterfaceModelParts(const int CommRank);

BoundingBoxType KRATOS_API(MAPPING_APPLICATION) ComputeLocalBoundingBox(const ModelPart& rModelPart);

BoundingBoxType KRATOS_API(MAPPING_APPLICATION) ComputeGlobalBoundingBox(const ModelPart& rModelPart);

std::string BoundingBoxStringStream(const BoundingBoxType& rBoundingBox);

void KRATOS_API(MAPPING_APPLICATION) SaveCurrentConfiguration(ModelPart& rModelPart);
void KRATOS_API(MAPPING_APPLICATION) RestoreCurrentConfiguration(ModelPart& rModelPart);

template<class TDataType>
void EraseNodalVariable(ModelPart& rModelPart, const Variable<TDataType>& rVariable)
{
    KRATOS_TRY;

    block_for_each(rModelPart.Nodes(), [&](Node& rNode){
        rNode.GetData().Erase(rVariable);
    });

    KRATOS_CATCH("");
}

void KRATOS_API(MAPPING_APPLICATION) FillBufferBeforeLocalSearch(const MapperLocalSystemPointerVector& rMapperLocalSystems,
                                 const std::vector<double>& rBoundingBoxes,
                                 const SizeType BufferSizeEstimate,
                                 std::vector<std::vector<double>>& rSendBuffer,
                                 std::vector<int>& rSendSizes);

void KRATOS_API(MAPPING_APPLICATION) CreateMapperInterfaceInfosFromBuffer(const std::vector<std::vector<double>>& rRecvBuffer,
                                          const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                                          const int CommRank,
                                          MapperInterfaceInfoPointerVectorType& rMapperInterfaceInfosContainer);

void FillBufferAfterLocalSearch(MapperInterfaceInfoPointerVectorType& rMapperInterfaceInfosContainer,
                                const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                                const int CommRank,
                                std::vector<std::vector<char>>& rSendBuffer,
                                std::vector<int>& rSendSizes);

void AssignInterfaceInfosAfterRemoteSearch(const MapperInterfaceInfoPointerVectorType& rMapperInterfaceInfosContainer,
                                           MapperLocalSystemPointerVectorPointer& rpMapperLocalSystems);

void DeserializeMapperInterfaceInfosFromBuffer(
    const std::vector<std::vector<char>>& rSendBuffer,
    const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
    const int CommRank,
    MapperInterfaceInfoPointerVectorType& rMapperInterfaceInfosContainer);

/**
 * @class MapperInterfaceInfoSerializer
 * @ingroup MappingApplication
 * @brief Helper class to serialize/deserialize a vector containing MapperInterfaceInfos
 * @details This class serializes the vector containing the MapperInterfaceInfos (Shared Ptrs)
 * The goal of this class is to have a more efficient/faster implementation than the
 * one of the Serializer by avoiding the casting that is done in the serializer when pointers
 * are serialized
 * @TODO test the performance against the Serializer
 * @author Philipp Bucher
 */
class KRATOS_API(MAPPING_APPLICATION) MapperInterfaceInfoSerializer
{
public:

    MapperInterfaceInfoSerializer(std::vector<MapperInterfaceInfoPointerType>& rMapperInterfaceInfosContainer,
                                  const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo)
        : mrInterfaceInfos(rMapperInterfaceInfosContainer)
        , mrpRefInterfaceInfo(rpRefInterfaceInfo->Create())
        { }

private:

    std::vector<MapperInterfaceInfoPointerType>& mrInterfaceInfos;
    MapperInterfaceInfoPointerType mrpRefInterfaceInfo;

    friend class Kratos::Serializer; // Adding "Kratos::" is needed bcs of the "MapperUtilities"-namespace

    virtual void save(Kratos::Serializer& rSerializer) const;
    virtual void load(Kratos::Serializer& rSerializer);
};

}  // namespace MapperUtilities.

}  // namespace Kratos.
