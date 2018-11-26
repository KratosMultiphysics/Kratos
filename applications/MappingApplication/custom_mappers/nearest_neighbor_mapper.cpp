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

// System includes

// External includes

// Project includes
#include "nearest_neighbor_mapper.h"
#include "mapping_application_variables.h"
#include "custom_utilities/mapper_typedefs.h"

namespace Kratos
{

typedef std::size_t IndexType;
typedef std::size_t SizeType;

/***********************************************************************************/
/* PUBLIC Methods */
/***********************************************************************************/
void NearestNeighborInterfaceInfo::ProcessSearchResult(const InterfaceObject::Pointer& rpInterfaceObject,
                                                      const double NeighborDistance)
{
    SetLocalSearchWasSuccessful();

    if (NeighborDistance < mNearestNeighborDistance) {
        mNearestNeighborDistance = NeighborDistance;
        mNearestNeighborId = rpInterfaceObject->pGetBaseNode()->GetValue(INTERFACE_EQUATION_ID);
    }
}

void NearestNeighborLocalSystem::CalculateAll(MatrixType& rLocalMappingMatrix,
                    EquationIdVectorType& rOriginIds,
                    EquationIdVectorType& rDestinationIds,
                    MapperLocalSystem::PairingStatus& rPairingStatus) const
{
    if (mInterfaceInfos.size() > 0) {
        rPairingStatus = MapperLocalSystem::PairingStatus::InterfaceInfoFound;

        if (rLocalMappingMatrix.size1() != 1 || rLocalMappingMatrix.size2() != 1) {
            rLocalMappingMatrix.resize(1, 1, false);
        }
        if (rOriginIds.size()      != 1) rOriginIds.resize(1);
        if (rDestinationIds.size() != 1) rDestinationIds.resize(1);

        int nearest_neighbor_id;
        double nearest_neighbor_distance;
        mInterfaceInfos[0]->GetValue(nearest_neighbor_id);
        mInterfaceInfos[0]->GetValue(nearest_neighbor_distance);

        for (SizeType i=1; i<mInterfaceInfos.size(); ++i) {
            // no check if this InterfaceInfo is an approximation is necessary
            // bcs this does not exist for NearestNeighbor
            double distance;
            mInterfaceInfos[i]->GetValue(distance);

            if (distance < nearest_neighbor_distance) {
                nearest_neighbor_distance = distance;
                mInterfaceInfos[i]->GetValue(nearest_neighbor_id);
            }
        }

        rLocalMappingMatrix(0,0) = 1.0;
        rOriginIds[0] = nearest_neighbor_id;
        KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;
        rDestinationIds[0] = mpNode->GetValue(INTERFACE_EQUATION_ID);
    }
    else ResizeToZero(rLocalMappingMatrix, rOriginIds, rDestinationIds, rPairingStatus);
}

std::string NearestNeighborLocalSystem::PairingInfo(const int EchoLevel, const int CommRank) const
{
    KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

    std::stringstream buffer;
    buffer << "NearestNeighborLocalSystem based on " << mpNode->Info();
    if (EchoLevel > 1) // TODO leave here?
        buffer << " at Coodinates " << Coordinates()[0] << " | " << Coordinates()[1] << " | " << Coordinates()[2];
    buffer << " in rank " << CommRank;
    return buffer.str();
}

/* Performs operations that are needed for Initialization and when the interface is updated (=> Remeshed)
I.e. Operations that can be performed several times in the livetime of the mapper
*/
template<class TSparseSpace, class TDenseSpace>
void NearestNeighborMapper<TSparseSpace, TDenseSpace>::InitializeInterface(Kratos::Flags MappingOptions)
{
    mpMapperLocalSystems->clear();

    const MapperLocalSystemPointer p_ref_local_system = Kratos::make_unique<NearestNeighborLocalSystem>();;

    KRATOS_ERROR_IF_NOT(mpInterfacePreprocessor) << "mpInterfacePreprocessor is a nullptr!" << std::endl;
    mpInterfacePreprocessor->CreateMapperLocalSystems(p_ref_local_system);

    BuildMappingMatrix(MappingOptions);
}

/* Performs operations that are needed for Initialization and when the interface is updated (All cases)
I.e. Operations that can be performed several times in the livetime of the mapper
*/
template<class TSparseSpace, class TDenseSpace>
void NearestNeighborMapper<TSparseSpace, TDenseSpace>::BuildMappingMatrix(Kratos::Flags MappingOptions)
{
    AssignInterfaceEquationIds(); // Has to be done ever time in case of overlapping interfaces!

    KRATOS_ERROR_IF_NOT(mpIntefaceCommunicator) << "mpIntefaceCommunicator is a nullptr!" << std::endl;

    const MapperInterfaceInfoUniquePointerType p_ref_interface_info = Kratos::make_unique<NearestNeighborInterfaceInfo>();
    const auto interface_object_construction_type_origin = InterfaceObject::Node_Coords;

    mpIntefaceCommunicator->ExchangeInterfaceData(mrModelPartDestination.GetCommunicator(),
                                                  MappingOptions,
                                                  p_ref_interface_info,
                                                  interface_object_construction_type_origin);

    KRATOS_ERROR_IF_NOT(mpMappingMatrixBuilder) << "mpMappingMatrixBuilder is a nullptr!" << std::endl;

    mpMappingMatrixBuilder->BuildMappingMatrix(mpInterfaceVectorContainerOrigin,
                                                  mpInterfaceVectorContainerDestination,
                                                  *mpMapperLocalSystems);

    // if (mEchoLevel > 0) PrintPairingInfo();
}

template<class TSparseSpace, class TDenseSpace>
void NearestNeighborMapper<TSparseSpace, TDenseSpace>::InitializeMappingMatrixBuilder()
{
    mpMappingMatrixBuilder = Kratos::make_unique<MappingMatrixBuilder<TSparseSpace, TDenseSpace>>();
}

template<>
void NearestNeighborMapper<MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType>::InitializeInterfaceCommunicator()
{
    Parameters search_settings(R"({})"); // TODO fill this
    search_settings.ValidateAndAssignDefaults(mMapperSettings);
    mpIntefaceCommunicator = Kratos::make_unique<InterfaceCommunicator>(mrModelPartOrigin,
                                                                      mpMapperLocalSystems,
                                                                      search_settings);
}

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
template<>
void NearestNeighborMapper<MapperDefinitions::MPISparseSpaceType, MapperDefinitions::DenseSpaceType>::InitializeInterfaceCommunicator()
{
    Parameters search_settings(R"({})"); // TODO fill this
    search_settings.ValidateAndAssignDefaults(mMapperSettings);
    mpIntefaceCommunicator = Kratos::make_unique<InterfaceCommunicatorMPI>(mrModelPartOrigin,
                                                                         mpMapperLocalSystems,
                                                                         search_settings);
}
#endif

template<class TSparseSpace, class TDenseSpace>
void NearestNeighborMapper<TSparseSpace, TDenseSpace>::MapInternal(
    const Variable<double>& rOriginVariable,
    const Variable<double>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    mpInterfaceVectorContainerOrigin->UpdateSystemVectorFromModelPart(rOriginVariable, MappingOptions);

    TSparseSpace::Mult(
        mpMappingMatrixBuilder->GetMappingMatrix(),
        mpInterfaceVectorContainerOrigin->GetVector(),
        mpInterfaceVectorContainerDestination->GetVector()); // rQd = rMdo * rQo

    mpInterfaceVectorContainerDestination->UpdateModelPartFromSystemVector(rDestinationVariable, MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void NearestNeighborMapper<TSparseSpace, TDenseSpace>::MapInternalTranspose(
    const Variable<double>& rOriginVariable,
    const Variable<double>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    mpInterfaceVectorContainerDestination->UpdateSystemVectorFromModelPart(rDestinationVariable, MappingOptions);

    TSparseSpace::TransposeMult(
        mpMappingMatrixBuilder->GetMappingMatrix(),
        mpInterfaceVectorContainerDestination->GetVector(),
        mpInterfaceVectorContainerOrigin->GetVector()); // rQo = rMdo^T * rQd

    mpInterfaceVectorContainerOrigin->UpdateModelPartFromSystemVector(rOriginVariable, MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void NearestNeighborMapper<TSparseSpace, TDenseSpace>::MapInternal(
    const Variable<array_1d<double, 3>>& rOriginVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    const auto& var_x_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + "_X");
    const auto& var_y_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + "_Y");
    const auto& var_z_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + "_Z");

    const auto& var_x_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + "_X");
    const auto& var_y_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + "_Y");
    const auto& var_z_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + "_Z");

    // X-Component
    mpInterfaceVectorContainerOrigin->UpdateSystemVectorFromModelPart(var_x_origin, MappingOptions);

    TSparseSpace::Mult(
        mpMappingMatrixBuilder->GetMappingMatrix(),
        mpInterfaceVectorContainerOrigin->GetVector(),
        mpInterfaceVectorContainerDestination->GetVector()); // rQd = rMdo * rQo

    mpInterfaceVectorContainerDestination->UpdateModelPartFromSystemVector(var_x_destination, MappingOptions);

    // Y-Component
    mpInterfaceVectorContainerOrigin->UpdateSystemVectorFromModelPart(var_y_origin, MappingOptions);

    TSparseSpace::Mult(
        mpMappingMatrixBuilder->GetMappingMatrix(),
        mpInterfaceVectorContainerOrigin->GetVector(),
        mpInterfaceVectorContainerDestination->GetVector()); // rQd = rMdo * rQo

    mpInterfaceVectorContainerDestination->UpdateModelPartFromSystemVector(var_y_destination, MappingOptions);

    // Z-Component
    mpInterfaceVectorContainerOrigin->UpdateSystemVectorFromModelPart(var_z_origin, MappingOptions);

    TSparseSpace::Mult(
        mpMappingMatrixBuilder->GetMappingMatrix(),
        mpInterfaceVectorContainerOrigin->GetVector(),
        mpInterfaceVectorContainerDestination->GetVector()); // rQd = rMdo * rQo

    mpInterfaceVectorContainerDestination->UpdateModelPartFromSystemVector(var_z_destination, MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void NearestNeighborMapper<TSparseSpace, TDenseSpace>::MapInternalTranspose(
    const Variable<array_1d<double, 3>>& rOriginVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    const auto& var_x_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + "_X");
    const auto& var_y_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + "_Y");
    const auto& var_z_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + "_Z");

    const auto& var_x_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + "_X");
    const auto& var_y_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + "_Y");
    const auto& var_z_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + "_Z");

    // X-Component
    mpInterfaceVectorContainerDestination->UpdateSystemVectorFromModelPart(var_x_destination, MappingOptions);

    TSparseSpace::TransposeMult(
        mpMappingMatrixBuilder->GetMappingMatrix(),
        mpInterfaceVectorContainerOrigin->GetVector(),
        mpInterfaceVectorContainerDestination->GetVector()); // rQo = rMdo^T * rQd

    mpInterfaceVectorContainerOrigin->UpdateModelPartFromSystemVector(var_x_origin, MappingOptions);

    // Y-Component
    mpInterfaceVectorContainerDestination->UpdateSystemVectorFromModelPart(var_y_destination, MappingOptions);

    TSparseSpace::TransposeMult(
        mpMappingMatrixBuilder->GetMappingMatrix(),
        mpInterfaceVectorContainerOrigin->GetVector(),
        mpInterfaceVectorContainerDestination->GetVector()); // rQo = rMdo^T * rQd

    mpInterfaceVectorContainerOrigin->UpdateModelPartFromSystemVector(var_y_origin, MappingOptions);

    // Z-Component
    mpInterfaceVectorContainerDestination->UpdateSystemVectorFromModelPart(var_z_destination, MappingOptions);

    TSparseSpace::TransposeMult(
        mpMappingMatrixBuilder->GetMappingMatrix(),
        mpInterfaceVectorContainerOrigin->GetVector(),
        mpInterfaceVectorContainerDestination->GetVector()); // rQo = rMdo^T * rQd

    mpInterfaceVectorContainerOrigin->UpdateModelPartFromSystemVector(var_z_origin, MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void NearestNeighborMapper<TSparseSpace, TDenseSpace>::ValidateInput(Parameters MapperSettings)
{
    MapperUtilities::CheckInterfaceModelParts(0);
    ValidateParameters(MapperSettings);

    // mEchoLevel = MapperSettings["echo_level"].GetInt();

    if (mMapperSettings["search_radius"].GetDouble() < 0.0) {
        const double search_radius = MapperUtilities::ComputeSearchRadius(mrModelPartOrigin,
                                        mrModelPartDestination,
                                        0);
        mMapperSettings["search_radius"].SetDouble(search_radius);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class NearestNeighborMapper< MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType >;
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
template class NearestNeighborMapper< MapperDefinitions::MPISparseSpaceType, MapperDefinitions::DenseSpaceType >;
#endif

}  // namespace Kratos.
