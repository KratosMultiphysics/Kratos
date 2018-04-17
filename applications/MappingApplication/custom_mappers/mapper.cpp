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

namespace Kratos
{
/***********************************************************************************/
/* PUBLIC Methods */
/***********************************************************************************/
template<class TSparseSpace, class TDenseSpace>
void Mapper<TSparseSpace, TDenseSpace>::UpdateInterface(Kratos::Flags MappingOptions, double SearchRadius)
{

}

template<class TSparseSpace, class TDenseSpace>
void Mapper<TSparseSpace, TDenseSpace>::Map(const Variable<double>& rOriginVariable,
                    const Variable<double>& rDestinationVariable,
                    Kratos::Flags MappingOptions)
{
    TMap(rOriginVariable, rDestinationVariable, MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void Mapper<TSparseSpace, TDenseSpace>::Map(const Variable< array_1d<double, 3> >& rOriginVariable,
                    const Variable< array_1d<double, 3> >& rDestinationVariable,
                    Kratos::Flags MappingOptions)
{
    TMap(rOriginVariable, rDestinationVariable, MappingOptions);
    // if (MappingOptions.Is(MapperFlags::CONSERVATIVE))
    // {
    //     MappingOptions.Reset(MapperFlags::CONSERVATIVE); // TODO test this!!!
    //     MappingOptions.Set(MapperFlags::USE_TRANSPOSE); // TODO test this!!!

    //     InverseMap(rOriginVariable, rDestinationVariable, MappingOptions);
    // }
    // else
    // {
    //     VectorComponentType var_component_x_origin =
    //         KratosComponents< VectorComponentType >::Get(rOriginVariable.Name()+std::string("_X"));
    //     VectorComponentType var_component_y_origin =
    //         KratosComponents< VectorComponentType >::Get(rOriginVariable.Name()+std::string("_Y"));
    //     VectorComponentType var_component_z_origin =
    //         KratosComponents< VectorComponentType >::Get(rOriginVariable.Name()+std::string("_Z"));

    //     VectorComponentType var_component_x_destination =
    //         KratosComponents< VectorComponentType >::Get(rDestinationVariable.Name()+std::string("_X"));
    //     VectorComponentType var_component_y_destination =
    //         KratosComponents< VectorComponentType >::Get(rDestinationVariable.Name()+std::string("_Y"));
    //     VectorComponentType var_component_z_destination =
    //         KratosComponents< VectorComponentType >::Get(rDestinationVariable.Name()+std::string("_Z"));

    //     // X-Component
    //     InitializeMappingStep<VectorComponentType>(var_component_x_origin,
    //                                                var_component_x_destination,
    //                                                MappingOptions);

    //     ExecuteMappingStep(MappingOptions);

    //     FinalizeMappingStep<VectorComponentType>(var_component_x_origin,
    //                                              var_component_x_destination,
    //                                              MappingOptions);

    //     // Y-Component
    //     InitializeMappingStep<VectorComponentType>(var_component_y_origin,
    //                                                var_component_y_destination,
    //                                                MappingOptions);

    //     ExecuteMappingStep(MappingOptions);

    //     FinalizeMappingStep<VectorComponentType>(var_component_y_origin,
    //                                              var_component_y_destination,
    //                                              MappingOptions);

    //     // Z-Component
    //     InitializeMappingStep<VectorComponentType>(var_component_z_origin,
    //                                                var_component_z_destination,
    //                                                MappingOptions);

    //     ExecuteMappingStep(MappingOptions);

    //     FinalizeMappingStep<VectorComponentType>(var_component_z_origin,
    //                                              var_component_z_destination,
    //                                              MappingOptions);
    // }
}

template<class TSparseSpace, class TDenseSpace>
void Mapper<TSparseSpace, TDenseSpace>::InverseMap(const Variable<double>& rOriginVariable,
                        const Variable<double>& rDestinationVariable,
                        Kratos::Flags MappingOptions)
{
    TInverseMap(rOriginVariable, rDestinationVariable, MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void Mapper<TSparseSpace, TDenseSpace>::InverseMap(const Variable< array_1d<double, 3> >& rOriginVariable,
                        const Variable< array_1d<double, 3> >& rDestinationVariable,
                        Kratos::Flags MappingOptions)
{
    TInverseMap(rOriginVariable, rDestinationVariable, MappingOptions);
}



/***********************************************************************************/
/* PROTECTED Methods */
/***********************************************************************************/
template<class TSparseSpace, class TDenseSpace>
Mapper<TSparseSpace, TDenseSpace>::Mapper(ModelPart& rModelPartOrigin,
                ModelPart& rModelPartDestination,
                Parameters MapperSettings,
                const bool IsMPIExecution) :
                    mrModelPartOrigin(rModelPartOrigin),
                    mrModelPartDestination(rModelPartDestination),
                    mIsMPIExecution(IsMPIExecution)
{
#ifndef KRATOS_USING_MPI // serial compilation
    KRATOS_ERROR_IF(IsMPIExecution) << "Trying to construct an MPI-Mapper "
        << "in a serial compilation!" << std::endl;
#endif

    // ValidateParameters(MapperSettings);
    // mEchoLevel = MapperSettings["echo_level"].GetInt();

    // mpInterfaceModelPart = Kratos::make_shared<ModelPart>("Mapper-Interface");
    // mpInterfacePreprocessor = Kratos::make_shared<InterfacePreprocessor>(mrModelPartDestination,
    //                                                                      mpInterfaceModelPart);

    // GenerateInterfaceModelPart();
    Initialize();
}

template<class TSparseSpace, class TDenseSpace>
void Mapper<TSparseSpace, TDenseSpace>::Initialize()
{
    InitializeInterfaceCommunicator();
    InitializeMappingOperationUtility();
}

/***********************************************************************************/
/* PRIVATE Methods */
/***********************************************************************************/
template<class TSparseSpace, class TDenseSpace>
void Mapper<TSparseSpace, TDenseSpace>::InitializeInterfaceCommunicator()
{
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
    int mpi_initialized;
    MPI_Initialized(&mpi_initialized);
    if (mpi_initialized) // parallel execution, i.e. mpi imported in python
    {
        int comm_size;
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
        if (comm_size > 1)
            mpInterfaceCommunicator = CreateMPIInterfaceCommunicator(mrModelPartOrigin, mpInterfaceModelPart);
    }
#else // serial compilation
    mpInterfaceCommunicator = CreateInterfaceCommunicator(mrModelPartOrigin, mpInterfaceModelPart);
#endif
}

template<class TSparseSpace, class TDenseSpace>
void Mapper<TSparseSpace, TDenseSpace>::InitializeMappingOperationUtility()
{
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
    int mpi_initialized;
    MPI_Initialized(&mpi_initialized);
    if (mpi_initialized) // parallel execution, i.e. mpi imported in python
    {
        int comm_size;
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
        if (comm_size > 1)
            mpMappingOperationUtility = CreateMPIMappingOperationUtility(mpInterfaceModelPart);
    }
#else // serial compilation
    mpMappingOperationUtility = CreateMappingOperationUtility(mpInterfaceModelPart);
#endif
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
template class Mapper< MapperDefinitions::UblasSparseSpaceType, MapperDefinitions::DenseSpaceType >;
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
template class Mapper< MapperDefinitions::TrilinosSparseSpaceType, MapperDefinitions::DenseSpaceType >;
#endif

}  // namespace Kratos.
