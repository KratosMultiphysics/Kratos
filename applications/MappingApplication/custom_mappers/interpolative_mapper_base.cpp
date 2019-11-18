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
#include "interpolative_mapper_base.h"
#include "custom_utilities/mapper_typedefs.h"
#include "custom_utilities/mapping_matrix_utilities.h"
#include "mapping_application_variables.h"
#include "custom_utilities/mapper_utilities.h"
#include "input_output/vtk_output.h"
#include "utilities/variable_utils.h"
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
#include "custom_searching/interface_communicator_mpi.h"
#endif

namespace Kratos {

typedef std::size_t IndexType;
typedef std::size_t SizeType;


template<class TSparseSpace, class TDenseSpace>
void InterpolativeMapperBase<TSparseSpace, TDenseSpace>::InitializeInterface(Kratos::Flags MappingOptions)
{
    CreateMapperLocalSystems(mrModelPartDestination.GetCommunicator(),
                             mMapperLocalSystems);

    BuildMappingMatrix(MappingOptions);
}

/* Performs operations that are needed for Initialization and when the interface is updated (All cases)
I.e. Operations that can be performed several times in the livetime of the mapper
*/
template<class TSparseSpace, class TDenseSpace>
void InterpolativeMapperBase<TSparseSpace, TDenseSpace>::BuildMappingMatrix(Kratos::Flags MappingOptions)
{
    AssignInterfaceEquationIds(); // Has to be done ever time in case of overlapping interfaces!

    KRATOS_ERROR_IF_NOT(mpIntefaceCommunicator) << "mpIntefaceCommunicator is a nullptr!" << std::endl;

    const MapperInterfaceInfoUniquePointerType p_ref_interface_info = GetMapperInterfaceInfo();

    mpIntefaceCommunicator->ExchangeInterfaceData(mrModelPartDestination.GetCommunicator(),
                                                  MappingOptions,
                                                  p_ref_interface_info);

    const int echo_level = mMapperSettings["echo_level"].GetInt();

    MappingMatrixUtilities::BuildMappingMatrix<TSparseSpace, TDenseSpace>(
        mpMappingMatrix,
        mpInterfaceVectorContainerOrigin->pGetVector(),
        mpInterfaceVectorContainerDestination->pGetVector(),
        mpInterfaceVectorContainerOrigin->GetModelPart(),
        mpInterfaceVectorContainerDestination->GetModelPart(),
        mMapperLocalSystems,
        echo_level);

    if (echo_level > 0) {
        PrintPairingInfo(echo_level);
    }
}

template<>
void InterpolativeMapperBase<MapperDefinitions::SparseSpaceType,
    MapperDefinitions::DenseSpaceType>::InitializeInterfaceCommunicator()
{
    mpIntefaceCommunicator = Kratos::make_unique<InterfaceCommunicator>(mrModelPartOrigin,
                                                                        mMapperLocalSystems,
                                                                        mMapperSettings);
}

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
template<>
void InterpolativeMapperBase<MapperDefinitions::MPISparseSpaceType,
    MapperDefinitions::DenseSpaceType>::InitializeInterfaceCommunicator()
{
    mpIntefaceCommunicator = Kratos::make_unique<InterfaceCommunicatorMPI>(mrModelPartOrigin,
                                                                           mMapperLocalSystems,
                                                                           mMapperSettings);
}
#endif

template<class TSparseSpace, class TDenseSpace>
void InterpolativeMapperBase<TSparseSpace, TDenseSpace>::MapInternal(
    const Variable<double>& rOriginVariable,
    const Variable<double>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    mpInterfaceVectorContainerOrigin->UpdateSystemVectorFromModelPart(rOriginVariable, MappingOptions);

    TSparseSpace::Mult(
        *mpMappingMatrix,
        mpInterfaceVectorContainerOrigin->GetVector(),
        mpInterfaceVectorContainerDestination->GetVector()); // rQd = rMdo * rQo

    mpInterfaceVectorContainerDestination->UpdateModelPartFromSystemVector(rDestinationVariable, MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void InterpolativeMapperBase<TSparseSpace, TDenseSpace>::MapInternalTranspose(
    const Variable<double>& rOriginVariable,
    const Variable<double>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    mpInterfaceVectorContainerDestination->UpdateSystemVectorFromModelPart(rDestinationVariable, MappingOptions);

    TSparseSpace::TransposeMult(
        *mpMappingMatrix,
        mpInterfaceVectorContainerDestination->GetVector(),
        mpInterfaceVectorContainerOrigin->GetVector()); // rQo = rMdo^T * rQd

    mpInterfaceVectorContainerOrigin->UpdateModelPartFromSystemVector(rOriginVariable, MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void InterpolativeMapperBase<TSparseSpace, TDenseSpace>::MapInternal(
    const Variable<array_1d<double, 3>>& rOriginVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    const std::vector<std::string> var_comps{"_X", "_Y", "_Z"};

    for (const auto& var_ext : var_comps) {
        const auto& var_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + var_ext);
        const auto& var_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + var_ext);

        mpInterfaceVectorContainerOrigin->UpdateSystemVectorFromModelPart(var_origin, MappingOptions);

        TSparseSpace::Mult(
            *mpMappingMatrix,
            mpInterfaceVectorContainerOrigin->GetVector(),
            mpInterfaceVectorContainerDestination->GetVector()); // rQd = rMdo * rQo

        mpInterfaceVectorContainerDestination->UpdateModelPartFromSystemVector(var_destination, MappingOptions);
    }
}

template<class TSparseSpace, class TDenseSpace>
void InterpolativeMapperBase<TSparseSpace, TDenseSpace>::MapInternalTranspose(
    const Variable<array_1d<double, 3>>& rOriginVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    const std::vector<std::string> var_comps{"_X", "_Y", "_Z"};

    for (const auto& var_ext : var_comps) {
        const auto& var_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + var_ext);
        const auto& var_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + var_ext);

        mpInterfaceVectorContainerDestination->UpdateSystemVectorFromModelPart(var_destination, MappingOptions);

        TSparseSpace::TransposeMult(
            *mpMappingMatrix,
            mpInterfaceVectorContainerDestination->GetVector(),
            mpInterfaceVectorContainerOrigin->GetVector()); // rQo = rMdo^T * rQd

        mpInterfaceVectorContainerOrigin->UpdateModelPartFromSystemVector(var_origin, MappingOptions);
    }
}

template<class TSparseSpace, class TDenseSpace>
void InterpolativeMapperBase<TSparseSpace, TDenseSpace>::ValidateInput()
{
    MapperUtilities::CheckInterfaceModelParts(0);

    Parameters mapper_default_settings(GetMapperDefaultSettings());
    mMapperSettings.ValidateAndAssignDefaults(mapper_default_settings);

    if (mMapperSettings["search_radius"].GetDouble() < 0.0) {
        const double search_radius = MapperUtilities::ComputeSearchRadius(
                                        mrModelPartOrigin,
                                        mrModelPartDestination,
                                        mMapperSettings["echo_level"].GetInt());
        mMapperSettings["search_radius"].SetDouble(search_radius);
    }
}

template<class TSparseSpace, class TDenseSpace>
void InterpolativeMapperBase<TSparseSpace, TDenseSpace>::PrintPairingInfo(const int EchoLevel)
{
    std::stringstream warning_msg;

    if (EchoLevel > 1) {
        // Initialize the values for printing later
        VariableUtils().SetNonHistoricalVariable(PAIRING_STATUS, 1, mrModelPartDestination.Nodes());
    }

    for (const auto& rp_local_sys : mMapperLocalSystems) {
        const auto pairing_status = rp_local_sys->GetPairingStatus();

        if (pairing_status != MapperLocalSystem::PairingStatus::InterfaceInfoFound) {
            warning_msg << rp_local_sys->PairingInfo(EchoLevel);

            if (pairing_status == MapperLocalSystem::PairingStatus::Approximation)
                warning_msg << " is using an approximation";
            else if (pairing_status == MapperLocalSystem::PairingStatus::NoInterfaceInfo)
                warning_msg << " has not found a neighbor";

            KRATOS_WARNING_ALL_RANKS("Mapper") << warning_msg.str() << std::endl; // TODO use data-comm of the destination-MP

            // reset the stringstream
            warning_msg.str( std::string() );
            warning_msg.clear();
        }
    }

    if (EchoLevel > 1) {
        // print a debug ModelPart to check the pairing

        std::string prefix = Info() + "_PairingStatus_";

        Parameters vtk_params( R"({
            "file_format"                        : "binary",
            "output_precision"                   : 7,
            "output_control_type"                : "step",
            "custom_name_prefix"                 : "",
            "save_output_files_in_folder"        : false,
            "nodal_data_value_variables"         : ["PAIRING_STATUS"]
        })");

        vtk_params["custom_name_prefix"].SetString(prefix);
        VtkOutput(mrModelPartDestination, vtk_params).PrintOutput();
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class InterpolativeMapperBase< MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType >;
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
template class InterpolativeMapperBase< MapperDefinitions::MPISparseSpaceType, MapperDefinitions::DenseSpaceType >;
#endif

}  // namespace Kratos.
