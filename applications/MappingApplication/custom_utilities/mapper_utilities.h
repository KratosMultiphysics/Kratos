//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_MAPPER_UTILITIES_H_INCLUDED)
#define  KRATOS_MAPPER_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "custom_utilities/mapper_flags.h"
#include "custom_utilities/mapper_local_system.h"

namespace Kratos
{
namespace MapperUtilities
{

typedef std::size_t SizeType;
typedef std::size_t IndexType;

typedef Node<3> NodeType;

typedef Kratos::unique_ptr<MapperInterfaceInfo> MapperInterfaceInfoUniquePointerType;

typedef Kratos::shared_ptr<MapperInterfaceInfo> MapperInterfaceInfoPointerType;
typedef std::vector<std::vector<MapperInterfaceInfoPointerType>> MapperInterfaceInfoPointerVectorType;
typedef Kratos::unique_ptr<MapperInterfaceInfoPointerVectorType> MapperInterfaceInfoPointerVectorPointerType;

typedef Kratos::unique_ptr<MapperLocalSystem> MapperLocalSystemPointer;
typedef std::vector<MapperLocalSystemPointer> MapperLocalSystemPointerVector;
typedef Kratos::shared_ptr<MapperLocalSystemPointerVector> MapperLocalSystemPointerVectorPointer;


template< class TVarType >
static void FillFunction(const NodeType& rNode,
                         const TVarType& rVariable,
                         double& rValue)
{
    rValue = rNode.FastGetSolutionStepValue(rVariable);
}

template< class TVarType >
static void FillFunctionNonHist(const NodeType& rNode,
                                const TVarType& rVariable,
                                double& rValue)
{
    rValue = rNode.GetValue(rVariable);
}

template< class TVarType >
static std::function<void(const NodeType&, const TVarType&, double&)>
GetFillFunction(const Kratos::Flags& rMappingOptions)
{
    if (rMappingOptions.Is(MapperFlags::FROM_NON_HISTORICAL))
        return &FillFunctionNonHist<TVarType>;
    return &FillFunction<TVarType>;
}

template< class TVarType >
static void UpdateFunction(NodeType& rNode,
                           const TVarType& rVariable,
                           const double Value,
                           const double Factor)
{
    rNode.FastGetSolutionStepValue(rVariable) = Value * Factor;
}

template< class TVarType >
static void UpdateFunctionWithAdd(NodeType& rNode,
                            const TVarType& rVariable,
                            const double Value,
                            const double Factor)
{
    rNode.FastGetSolutionStepValue(rVariable) += Value * Factor;
}

template< class TVarType >
static void UpdateFunctionNonHist(NodeType& rNode,
                            const TVarType& rVariable,
                            const double Value,
                            const double Factor)
{
    rNode.GetValue(rVariable) = Value * Factor;
}

template< class TVarType >
static void UpdateFunctionNonHistWithAdd(NodeType& rNode,
                            const TVarType& rVariable,
                            const double Value,
                            const double Factor)
{
    rNode.GetValue(rVariable) += Value * Factor;
}

template< class TVarType >
static std::function<void(NodeType&, const TVarType&, const double, const double)>
GetUpdateFunction(const Kratos::Flags& rMappingOptions)
{
    if (rMappingOptions.Is(MapperFlags::ADD_VALUES) && rMappingOptions.Is(MapperFlags::TO_NON_HISTORICAL))
        return &UpdateFunctionNonHistWithAdd<TVarType>;
    if (rMappingOptions.Is(MapperFlags::ADD_VALUES))
        return &UpdateFunctionWithAdd<TVarType>;
    if (rMappingOptions.Is(MapperFlags::TO_NON_HISTORICAL))
        return &UpdateFunctionNonHist<TVarType>;
    return &UpdateFunction<TVarType>;
}

template< class TVectorType, class TVarType >
void UpdateSystemVectorFromModelPart(TVectorType& rVector,
                        ModelPart& rModelPart,
                        const TVarType& rVariable,
                        const Kratos::Flags& rMappingOptions)
{
    // Here we construct a function pointer to not have the if all the time inside the loop
    const auto fill_fct = MapperUtilities::GetFillFunction<TVarType>(rMappingOptions);

    const int num_local_nodes = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    const auto nodes_begin = rModelPart.GetCommunicator().LocalMesh().NodesBegin();

    #pragma omp parallel for
    for (int i=0; i<num_local_nodes; i++) {
        fill_fct(*(nodes_begin + i), rVariable, rVector[i]);
    }
}

template< class TVectorType, class TVarType >
void UpdateModelPartFromSystemVector(const TVectorType& rVector,
            ModelPart& rModelPart,
            const TVarType& rVariable,
            const Kratos::Flags& rMappingOptions)
{
    const double factor = rMappingOptions.Is(MapperFlags::SWAP_SIGN) ? -1.0 : 1.0;

    // Here we construct a function pointer to not have the if all the time inside the loop
    const auto update_fct = std::bind(MapperUtilities::GetUpdateFunction<TVarType>(rMappingOptions),
                                        std::placeholders::_1,
                                        std::placeholders::_2,
                                        std::placeholders::_3,
                                        factor);
    const int num_local_nodes = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    const auto nodes_begin = rModelPart.GetCommunicator().LocalMesh().NodesBegin();

    #pragma omp parallel for
    for (int i=0; i<num_local_nodes; i++) {
        update_fct(*(nodes_begin + i), rVariable, rVector[i]);
    }
}

/**
* @brief Assigning INTERFACE_EQUATION_IDs to the nodes, with and without MPI
* This function assigns the INTERFACE_EQUATION_IDs to the nodes, which
* act as EquationIds for the MappingMatrix. This work with and without MPI,
* in MPI a ScanSum is performed with the local number of nodes
* @param rModelPartCommunicator The Modelpart-Communicator to be used
* @author Philipp Bucher
*/
void AssignInterfaceEquationIds(Communicator& rModelPartCommunicator);

inline int ComputeNumberOfNodes(ModelPart& rModelPart)
{
    int num_nodes = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    rModelPart.GetCommunicator().SumAll(num_nodes); // Compute the sum among the partitions
    return num_nodes;
}

inline int ComputeNumberOfConditions(ModelPart& rModelPart)
{
    int num_conditions = rModelPart.GetCommunicator().LocalMesh().NumberOfConditions();
    rModelPart.GetCommunicator().SumAll(num_conditions); // Compute the sum among the partitions
    return num_conditions;
}

inline int ComputeNumberOfElements(ModelPart& rModelPart)
{
    int num_elements = rModelPart.GetCommunicator().LocalMesh().NumberOfElements();
    rModelPart.GetCommunicator().SumAll(num_elements); // Compute the sum among the partitions
    return num_elements;
}

inline double ComputeDistance(const array_1d<double, 3>& rCoords1,
                              const array_1d<double, 3>& rCoords2)
{
    return std::sqrt( std::pow(rCoords1[0] - rCoords2[0] , 2) +
                      std::pow(rCoords1[1] - rCoords2[1] , 2) +
                      std::pow(rCoords1[2] - rCoords2[2] , 2) );
}

template <typename T>
inline double ComputeMaxEdgeLengthLocal(const T& rEntityContainer)
{
    double max_element_size = 0.0f;
    // Loop through each edge of a geometrical entity ONCE
    for (auto& r_entity : rEntityContainer) {
        for (std::size_t i = 0; i < (r_entity.GetGeometry().size() - 1); ++i) {
            for (std::size_t j = i + 1; j < r_entity.GetGeometry().size(); ++j) {
                double edge_length = ComputeDistance(r_entity.GetGeometry()[i].Coordinates(),
                                                        r_entity.GetGeometry()[j].Coordinates());
                max_element_size = std::max(max_element_size, edge_length);
            }
        }
    }
    return max_element_size;
}

inline double ComputeMaxEdgeLengthLocal(const ModelPart::NodesContainerType& rNodes)
{
    double max_element_size = 0.0f;
    // TODO modify loop such that it loop only once over the nodes
    for (auto& r_node_1 : rNodes) {
        for (auto& r_node_2 : rNodes) {
            double edge_length = ComputeDistance(r_node_1.Coordinates(),
                                                    r_node_2.Coordinates());
            max_element_size = std::max(max_element_size, edge_length);
        }
    }
    return max_element_size;
}

double ComputeSearchRadius(ModelPart& rModelPart, int EchoLevel);

inline double ComputeSearchRadius(ModelPart& rModelPart1, ModelPart& rModelPart2, int EchoLevel)
{
    double search_radius = std::max(ComputeSearchRadius(rModelPart1, EchoLevel),
                                    ComputeSearchRadius(rModelPart2, EchoLevel));
    return search_radius;
}

void CheckInterfaceModelParts(const int CommRank);

std::vector<double> ComputeLocalBoundingBox(ModelPart& rModelPart);

void ComputeBoundingBoxesWithTolerance(const std::vector<double>& rBoundingBoxes,
                                       const double Tolerance,
                                       std::vector<double>& rBoundingBoxesWithTolerance);

std::string BoundingBoxStringStream(const std::vector<double>& rBoundingBox);

bool PointIsInsideBoundingBox(const std::vector<double>& rBoundingBox,
                              const Point& rPoint);

void FillBufferBeforeLocalSearch(const MapperLocalSystemPointerVectorPointer& rpMapperLocalSystems,
                                 const std::vector<double>& rBoundingBoxes,
                                 const SizeType BufferSizeEstimate,
                                 std::vector<std::vector<double>>& rSendBuffer,
                                 std::vector<int>& rSendSizes);

void CreateMapperInterfaceInfosFromBuffer(const std::vector<std::vector<double>>& rRecvBuffer,
                                          const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                                          const int CommRank,
                                          MapperInterfaceInfoPointerVectorPointerType& rpMapperInterfaceInfosContainer);

void FillBufferAfterLocalSearch(const MapperInterfaceInfoPointerVectorPointerType& rpMapperInterfaceInfosContainer,
                                const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                                const int CommRank,
                                std::vector<std::vector<char>>& rSendBuffer,
                                std::vector<int>& rSendSizes);

void AssignInterfaceInfosAfterRemoteSearch(const MapperInterfaceInfoPointerVectorPointerType& rpMapperInterfaceInfosContainer,
                                           MapperLocalSystemPointerVectorPointer& rpMapperLocalSystems);

void DeserializeMapperInterfaceInfosFromBuffer(
    const std::vector<std::vector<char>>& rSendBuffer,
    const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
    const int CommRank,
    MapperInterfaceInfoPointerVectorPointerType& rpMapperInterfaceInfosContainer);

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
class MapperInterfaceInfoSerializer
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

    friend class Kratos::Serializer; // Adding "Kratos::" is nedded bcs of the "MapperUtilities"-namespace

    virtual void save(Kratos::Serializer& rSerializer) const;
    virtual void load(Kratos::Serializer& rSerializer);
};

}  // namespace MapperUtilities.

}  // namespace Kratos.

#endif // KRATOS_MAPPER_UTILITIES_H_INCLUDED  defined
