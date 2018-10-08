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
#include "mapper.h"
#include "custom_utilities/mapper_typedefs.h"
#include "custom_utilities/mapper_utilities.h"
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
#include "custom_searching/interface_communicator_mpi.h"
#endif

namespace Kratos
{
/***********************************************************************************/
/* PUBLIC Methods */
/***********************************************************************************/

/***********************************************************************************/
/* PROTECTED Methods */
/***********************************************************************************/
template<class TSparseSpace, class TDenseSpace>
Mapper<TSparseSpace, TDenseSpace>::Mapper(ModelPart& rModelPartOrigin,
                ModelPart& rModelPartDestination,
                Parameters MapperSettings) :
                    mrModelPartOrigin(rModelPartOrigin),
                    mrModelPartDestination(rModelPartDestination),
                    mGeneralMapperSettings(MapperSettings)
{
    ValidateInput(MapperSettings);

    // TODO throw error in case of MPI-execution with one core


    mEchoLevel = MapperSettings["echo_level"].GetInt();
}

/* This function initializes the Mapper
I.e. Operations that should be performed ONLY ONCE in the livetime of the mapper
*/
template<class TSparseSpace, class TDenseSpace>
void Mapper<TSparseSpace, TDenseSpace>::Initialize()
{
    mpMapperLocalSystems = Kratos::make_shared<MapperLocalSystemPointerVector>();

    mpInterfacePreprocessor = Kratos::make_unique<InterfacePreprocessor>(mrModelPartDestination,
                                                                         mpMapperLocalSystems);
    InitializeMappingOperationUtility();

    InitializeInterfaceCommunicator();

    InitializeInterface();
}

/* Performs operations that are needed for Initialization and when the interface is updated (=> Remeshed)
I.e. Operations that can be performed several times in the livetime of the mapper
*/
template<class TSparseSpace, class TDenseSpace>
void Mapper<TSparseSpace, TDenseSpace>::InitializeInterface(Kratos::Flags MappingOptions)
{
    mpMapperLocalSystems->clear();

    const MapperLocalSystemPointer p_ref_local_system = GetMapperLocalSystem();

    KRATOS_ERROR_IF_NOT(mpInterfacePreprocessor) << "mpInterfacePreprocessor is a nullptr!" << std::endl;
    mpInterfacePreprocessor->CreateMapperLocalSystems(p_ref_local_system);

    BuildMappingMatrix(MappingOptions);
}

/* Performs operations that are needed for Initialization and when the interface is updated (All cases)
I.e. Operations that can be performed several times in the livetime of the mapper
*/
template<class TSparseSpace, class TDenseSpace>
void Mapper<TSparseSpace, TDenseSpace>::BuildMappingMatrix(Kratos::Flags MappingOptions)
{
    AssignInterfaceEquationIds(); // Has to be done ever time in case of overlapping interfaces!

    KRATOS_ERROR_IF_NOT(mpIntefaceCommunicator) << "mpIntefaceCommunicator is a nullptr!" << std::endl;

    const MapperInterfaceInfoUniquePointerType p_ref_interface_info = GetMapperInterfaceInfo();
    const InterfaceObject::ConstructionType interface_object_construction_type_origin =
        GetInterfaceObjectConstructionTypeOrigin();

    mpIntefaceCommunicator->ExchangeInterfaceData(mrModelPartDestination.GetCommunicator(),
                                             MappingOptions,
                                             p_ref_interface_info,
                                             interface_object_construction_type_origin);

    KRATOS_ERROR_IF_NOT(mpMappingOperationUtility) << "mpMappingOperationUtility is a nullptr!" << std::endl;

    // this function can always be called, it won't do anything if the sizes are correct
    mpMappingOperationUtility->BuildMappingSystem(mpMdo, mpQo, mpQd,
                                                  mrModelPartOrigin,
                                                  mrModelPartDestination,
                                                  *mpMapperLocalSystems);

    if (mEchoLevel > 0) PrintPairingInfo();
}

template<class TSparseSpace, class TDenseSpace>
void Mapper<TSparseSpace, TDenseSpace>::InitializeMappingOperationUtility()
{
    // here we could return the MatrixFree variant in the future
    Parameters utility_settings(R"({})"); // TODO fill this
    utility_settings.ValidateAndAssignDefaults(mGeneralMapperSettings);
    mpMappingOperationUtility = Kratos::make_unique<MatrixBasedMappingOperationUtility<TSparseSpace, TDenseSpace>>(utility_settings);
}

template<>
void Mapper<MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType>::InitializeInterfaceCommunicator()
{
    // The current solution is just a workaround since cloning doesn't work
    // Parameters search_settings(mGeneralMapperSettings.Clone());
    // Parameters search_settings = mGeneralMapperSettings.Clone(); //TODO I think this would be the correct solution ...
    Parameters search_settings(R"({})"); // TODO fill this
    search_settings.ValidateAndAssignDefaults(mGeneralMapperSettings);
    mpIntefaceCommunicator = Kratos::make_unique<InterfaceCommunicator>(mrModelPartOrigin,
                                                                      mpMapperLocalSystems,
                                                                      search_settings);
}

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
template<>
void Mapper<MapperDefinitions::MPISparseSpaceType, MapperDefinitions::DenseSpaceType>::InitializeInterfaceCommunicator()
{
    // The current solution is just a workaround since cloning doesn't work
    // Parameters search_settings = mGeneralMapperSettings.Clone();
    Parameters search_settings(R"({})"); // TODO fill this
    search_settings.ValidateAndAssignDefaults(mGeneralMapperSettings);
    mpIntefaceCommunicator = Kratos::make_unique<InterfaceCommunicatorMPI>(mrModelPartOrigin,
                                                                         mpMapperLocalSystems,
                                                                         search_settings);
}
#endif

/*
This function contains the actual Implementation Of the UpdateInterface function
It is done like this bcs in this way the same operation can be called on the InverseMapper (if it is initialized)
*/
template<class TSparseSpace, class TDenseSpace>
void Mapper<TSparseSpace, TDenseSpace>::UpdateInterfaceInternal(Kratos::Flags MappingOptions, double SearchRadius)
{
    // Set the Flags according to the type of remeshing
    if (MappingOptions.Is(MapperFlags::REMESHED))
    {
        if (MappingOptions.IsDefined(MapperFlags::ORIGIN_ONLY))
        {
            KRATOS_INFO("Mapper-UpdateInterface") << "If the domain is remeshed then "
                << "setting \"ORIGIN_ONLY\" has no effect" << std::endl;
            MappingOptions.Reset(MapperFlags::ORIGIN_ONLY);
        }
        if (MappingOptions.IsDefined(MapperFlags::DESTINATION_ONLY))
        {
            KRATOS_INFO("Mapper-UpdateInterface") << "If the domain is remeshed then "
                << "setting \"DESTINATION_ONLY\" has no effect" << std::endl;
            MappingOptions.Reset(MapperFlags::DESTINATION_ONLY);
        }
        InitializeInterface(MappingOptions);
    }
    else BuildMappingMatrix(MappingOptions);
}


// template<class TSparseSpace, class TDenseSpace> template<typename T>
// void Mapper<TSparseSpace, TDenseSpace>::TestFunction(T someParam)
// {

// }
/***********************************************************************************/
/* PRIVATE Methods */
/***********************************************************************************/

template<class TSparseSpace, class TDenseSpace>
void Mapper<TSparseSpace, TDenseSpace>::ValidateInput(Parameters MapperSettings)
{
    MapperUtilities::CheckInterfaceModelParts(0);
    ValidateParameters(MapperSettings);

    mEchoLevel = MapperSettings["echo_level"].GetInt();

    if (mGeneralMapperSettings["search_radius"].GetDouble() < 0.0)
    {
        const double search_radius = MapperUtilities::ComputeSearchRadius(mrModelPartOrigin,
                                        mrModelPartDestination,
                                        0);
        mGeneralMapperSettings["search_radius"].SetDouble(search_radius);
    }
}

template<class TSparseSpace, class TDenseSpace>
void Mapper<TSparseSpace, TDenseSpace>::AssignInterfaceEquationIds()
{
    MapperUtilities::AssignInterfaceEquationIds(mrModelPartOrigin.GetCommunicator());
    MapperUtilities::AssignInterfaceEquationIds(mrModelPartDestination.GetCommunicator());
}

template<class TSparseSpace, class TDenseSpace>
void Mapper<TSparseSpace, TDenseSpace>::PrintPairingInfo()
{
    const int comm_rank = mrModelPartDestination.GetCommunicator().MyPID();
    std::stringstream warning_msg;

    for (const auto& rp_local_sys : *mpMapperLocalSystems)
    {
        const auto pairing_status = rp_local_sys->GetPairingStatus();

        if (pairing_status != MapperLocalSystem::PairingStatus::InterfaceInfoFound)
        {
            warning_msg << rp_local_sys->PairingInfo(mEchoLevel, comm_rank);

            if (pairing_status == MapperLocalSystem::PairingStatus::Approximation)
                warning_msg << " is using an approximation";
            else if (pairing_status == MapperLocalSystem::PairingStatus::NoInterfaceInfo)
                warning_msg << " has not found a neighbor";

            KRATOS_WARNING("Mapper") << warning_msg.str() << std::endl;

            // reset the stringstream
            warning_msg.str( std::string() );
            warning_msg.clear();
        }
    }
}

// /// input stream function
// inline std::istream & operator >> (std::istream& rIStream, Mapper& rThis);

// /// output stream function
// inline std::ostream & operator << (std::ostream& rOStream, const Mapper& rThis) {
//   rThis.PrintInfo(rOStream);
//   rOStream << " : " << std::endl;
//   rThis.PrintData(rOStream);
//   return rOStream;
// }

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class Mapper< MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType >;
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
template class Mapper< MapperDefinitions::MPISparseSpaceType, MapperDefinitions::DenseSpaceType >;
#endif

}  // namespace Kratos.
