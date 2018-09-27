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
#include "matrix_based_mapping_operation_utility.h"
#include "custom_utilities/mapper_typedefs.h"
#include "custom_utilities/mapper_utilities.h"
#include "mapping_application_variables.h"

namespace Kratos
{
using SparseSpaceType = MapperDefinitions::MPISparseSpaceType;
using DenseSpaceType = MapperDefinitions::DenseSpaceType;

using UtilityType = MatrixBasedMappingOperationUtility<SparseSpaceType, DenseSpaceType>;

using EquationIdVectorType = typename MapperLocalSystem::EquationIdVectorType;
using MappingWeightsVector = typename MapperLocalSystem::MappingWeightsVector;

/***********************************************************************************/
/* PUBLIC Methods */
/***********************************************************************************/
template<>
UtilityType::MatrixBasedMappingOperationUtility(Parameters Settings)
    : MappingOperationUtility<SparseSpaceType, DenseSpaceType>(Settings)
{
    KRATOS_ERROR_IF_NOT(SparseSpaceType::IsDistributed())
        << "Using a non-distributed Space!" << std::endl;
}

template<>
void UtilityType::ResizeAndInitializeVectors(
    TSystemMatrixUniquePointerType& rpMdo,
    TSystemVectorUniquePointerType& rpQo,
    TSystemVectorUniquePointerType& rpQd,
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    MapperLocalSystemPointerVector& rMapperLocalSystems) const
{
    // big TODO what if the rank doesn't have local nodes ... ?

    const int num_local_nodes_orig = rModelPartOrigin.GetCommunicator().LocalMesh().NumberOfNodes();
    const int num_local_nodes_dest = rModelPartDestination.GetCommunicator().LocalMesh().NumberOfNodes();

    const Epetra_MpiComm epetra_comm(MPI_COMM_WORLD); // TODO do I have to save this as member???

    std::vector<int> global_elements_orig(num_local_nodes_orig);
    std::vector<int> global_elements_dest(num_local_nodes_dest);

    const auto nodes_begin_orig = rModelPartOrigin.GetCommunicator().LocalMesh().NodesBegin();
    #pragma omp parallel for
    for (int i=0; i<num_local_nodes_orig; ++i) {
        global_elements_orig[i] = ( nodes_begin_orig + i )->GetValue(INTERFACE_EQUATION_ID);
    }

    const auto nodes_begin_dest = rModelPartDestination.GetCommunicator().LocalMesh().NodesBegin();
    #pragma omp parallel for
    for (int i=0; i<num_local_nodes_dest; ++i) {
        global_elements_dest[i] = ( nodes_begin_dest + i )->GetValue(INTERFACE_EQUATION_ID);
    }

    // Construct vectors containing all the indices this processor contributes to
    std::vector<int> row_indices(num_local_nodes_dest*2); // using number of nodes as size estimation
    std::vector<int> col_indices(num_local_nodes_orig*2);

    EquationIdVectorType origin_ids;
    EquationIdVectorType destination_ids;

    for (auto& rp_local_sys : rMapperLocalSystems)
    {
        rp_local_sys->EquationIdVectors(origin_ids, destination_ids);

        KRATOS_DEBUG_ERROR_IF(origin_ids.size() != destination_ids.size())
            << "EquationID vectors have size mismatch!" << std::endl;

        col_indices.reserve( col_indices.size() + origin_ids.size() );
        col_indices.insert( col_indices.end(), origin_ids.begin(), origin_ids.end() );
        row_indices.reserve( row_indices.size() + destination_ids.size() );
        row_indices.insert( row_indices.end(), destination_ids.begin(), destination_ids.end() );
    }

    // "Uniqueify" the vectors
    std::sort( col_indices.begin(), col_indices.end() );
    col_indices.erase( std::unique( col_indices.begin(), col_indices.end() ), col_indices.end() );
    std::sort( row_indices.begin(), row_indices.end() );
    row_indices.erase( std::unique( row_indices.begin(), row_indices.end() ), row_indices.end() );


    // Epetra_Map (long long NumGlobalElements, int NumMyElements, const long long *MyGlobalElements, int IndexBase, const Epetra_Comm &Comm)

    const int num_global_elements = -1; // this means its gonna be computed by Epetra_Map
    const int index_base = 0; // for C/C++

    Epetra_Map epetra_col_map(num_global_elements,
                              col_indices.size(),
                              col_indices.data(), // taken as const
                              index_base,
                              epetra_comm);

    Epetra_Map epetra_row_map(num_global_elements,
                              row_indices.size(),
                              row_indices.data(), // taken as const
                              index_base,
                              epetra_comm);

    std::cout << epetra_row_map << std::endl;
    std::cout << epetra_col_map << std::endl;

    Epetra_Map epetra_domain_map(num_global_elements,
                                 num_local_nodes_orig,
                                 global_elements_orig.data(), // taken as const
                                 index_base,
                                 epetra_comm);

    Epetra_Map epetra_range_map(num_global_elements,
                                num_local_nodes_dest,
                                global_elements_dest.data(), // taken as const
                                index_base,
                                epetra_comm);


    // explanation in here: https://trilinos.org/docs/dev/packages/epetra/doc/html/classEpetra__CrsGraph.html
    const int num_indices_per_row = 5; // TODO this is to be tested => set to zero maybe ...

    // TODO do I even need the graph? I think I could directly use the Matrix and perform the same operations ...
    // Performance optimization see https://trilinos.org/docs/dev/packages/epetra/doc/html/classEpetra__CrsGraph.html
    Epetra_FECrsGraph epetra_graph(Epetra_DataAccess::Copy,
                                   epetra_row_map,
                                   epetra_col_map,
                                   num_indices_per_row);


    for (auto& rp_local_sys : rMapperLocalSystems)
    {
        rp_local_sys->EquationIdVectors(origin_ids, destination_ids);

        KRATOS_DEBUG_ERROR_IF(origin_ids.size() != destination_ids.size())
            << "EquationID vectors have size mismatch!" << std::endl;

        if (origin_ids.size() > 0)
        {
            const int ierr = epetra_graph.InsertGlobalIndices(
                destination_ids.size(), destination_ids.data(),
                origin_ids.size(),      origin_ids.data());

            KRATOS_ERROR_IF( ierr < 0 ) << "Epetra failure in Epetra_FECrsGraph.InsertGlobalIndices. "
                << "Error code: " << ierr << std::endl;
        }
    }

    int ierr = epetra_graph.GlobalAssemble(epetra_domain_map, epetra_range_map); // TODO check if it should call "FillComplete"
    KRATOS_ERROR_IF( ierr != 0 ) << "Epetra failure in Epetra_FECrsGraph.GlobalAssemble. "
        << "Error code: " << ierr << std::endl;

    epetra_graph.OptimizeStorage();

    // // TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(Copy,Agraph) );
    // // https://trilinos.org/docs/dev/packages/epetra/doc/html/Epetra__DataAccess_8h.html#ad1a985e79f94ad63030815a0d7d90928
    TSystemMatrixUniquePointerType p_Mdo = Kratos::make_unique<TSystemMatrixType>(Epetra_DataAccess::Copy, epetra_graph);
    rpMdo.swap(p_Mdo);

    TSystemVectorUniquePointerType p_new_vector_destination = Kratos::make_unique<TSystemVectorType>(epetra_range_map);
    rpQd.swap(p_new_vector_destination);

    TSystemVectorUniquePointerType p_new_vector_origin = Kratos::make_unique<TSystemVectorType>(epetra_domain_map);
    rpQo.swap(p_new_vector_origin);

    // rpQo->GlobalAssemble();
    // rpQd->GlobalAssemble();
    MappingWeightsVector mapping_weights;

    std::cout << "Before Assembly" << std::endl;

    for (auto& rp_local_sys : rMapperLocalSystems)
    {
        rp_local_sys->CalculateLocalSystem(mapping_weights, origin_ids, destination_ids);

        KRATOS_DEBUG_ERROR_IF(mapping_weights.size() != origin_ids.size())
            << "OriginID vector size mismatch" << std::endl;
        KRATOS_DEBUG_ERROR_IF(mapping_weights.size() != destination_ids.size())
            << "DestinationID vector size mismatch" << std::endl;

        KRATOS_WATCH(mapping_weights)
        KRATOS_WATCH(origin_ids)
        KRATOS_WATCH(destination_ids)
        std::cout << std::endl;

        if (mapping_weights.size() > 0)
        {
            const int ierr = rpMdo->SumIntoGlobalValues(
                destination_ids.size(), destination_ids.data(),
                origin_ids.size(),      origin_ids.data(),
                mapping_weights.data(),
                Epetra_FECrsMatrix::ROW_MAJOR );

            KRATOS_ERROR_IF( ierr < 0 ) << "Epetra failure in Epetra_FECrsMatrix.SumIntoGlobalValues. "
                << "Error code: " << ierr << std::endl;
        }

        rp_local_sys->Clear();
    }

    rModelPartOrigin.GetCommunicator().Barrier();

    std::cout << "After Assembly\n" << *rpMdo << std::endl;

    rpMdo->GlobalAssemble(epetra_domain_map, epetra_range_map);

    rModelPartOrigin.GetCommunicator().Barrier();

    std::cout << "After GlobalAssemble\n" << *rpMdo << std::endl;

    if (GetEchoLevel() > 2)
        SparseSpaceType::WriteMatrixMarketMatrix("TrilinosMappingMatrix", *rpMdo, false);
}

// The "Build" function
template<>
void UtilityType::BuildMappingMatrix(
    const MapperLocalSystemPointerVector& rMapperLocalSystems,
    TSystemMatrixType& rMdo) const
{
    // MappingWeightsVector mapping_weights;

    // EquationIdVectorType origin_ids;
    // EquationIdVectorType destination_ids;

    // std::cout << "Before Assembly" << std::endl;

    // for (auto& rp_local_sys : rMapperLocalSystems)
    // {
    //     rp_local_sys->CalculateLocalSystem(mapping_weights, origin_ids, destination_ids);

    //     KRATOS_DEBUG_ERROR_IF(mapping_weights.size() != origin_ids.size())
    //         << "OriginID vector size mismatch" << std::endl;
    //     KRATOS_DEBUG_ERROR_IF(mapping_weights.size() != destination_ids.size())
    //         << "DestinationID vector size mismatch" << std::endl;

    //     KRATOS_WATCH(mapping_weights)
    //     KRATOS_WATCH(origin_ids)
    //     KRATOS_WATCH(destination_ids)
    //     std::cout << std::endl;

    //     if (mapping_weights.size() > 0)
    //     {
    //         const int ierr = rMdo.SumIntoGlobalValues(
    //             destination_ids.size(), destination_ids.data(),
    //             origin_ids.size(),      origin_ids.data(),
    //             mapping_weights.data());

    //         KRATOS_ERROR_IF( ierr < 0 ) << "Epetra failure in Epetra_FECrsMatrix.SumIntoGlobalValues. "
    //             << "Error code: " << ierr << std::endl;
    //     }

    //     rp_local_sys->Clear();
    // }

    // std::cout << "After Assembly" << std::endl;

    // // rMdo.GlobalAssemble(epetra_domain_map, epetra_range_map);

    // if (GetEchoLevel() > 2)
    //     SparseSpaceType::WriteMatrixMarketMatrix("TrilinosMappingMatrix", rMdo, false);
}

template< class TVarType >
void FillSystemVector(UtilityType::TSystemVectorType& rVector,
                        ModelPart& rModelPart,
                        const TVarType& rVariable,
                        const Kratos::Flags& rMappingOptions,
                        const int EchoLevel)
{
    // Here we construct a function pointer to not have the if all the time inside the loop
    const auto fill_fct = MapperUtilities::GetFillFunction<TVarType>(rMappingOptions);

    for (const auto& r_node : rModelPart.GetCommunicator().LocalMesh().Nodes())
    {
        const int global_index = r_node.GetValue(INTERFACE_EQUATION_ID); // TODO find a better solution, this might have been overwritten by a different mapper!!! (if parts of the interface are shared)
        const int local_index = rVector.Map().LID(global_index);
        fill_fct(r_node, rVariable, rVector[0][local_index]);
    }

    // rVector.GlobalAssemble(); // I am quite sure this is not needed, since one node is one entry ...
    if (EchoLevel > 2)
        SparseSpaceType::WriteMatrixMarketVector("TrilinosFillSystemVector", rVector);
}

template< class TVarType >
void Update(UtilityType::TSystemVectorType& rVector,
            ModelPart& rModelPart,
            const TVarType& rVariable,
            const Kratos::Flags& rMappingOptions,
            const int EchoLevel)
{
    const double factor = rMappingOptions.Is(MapperFlags::SWAP_SIGN) ? -1.0 : 1.0;

    // Here we construct a function pointer to not have the if all the time inside the loop
    const auto update_fct = std::bind(MapperUtilities::GetUpdateFunction<TVarType>(rMappingOptions),
                                        std::placeholders::_1,
                                        std::placeholders::_2,
                                        std::placeholders::_3,
                                        factor);

    for (auto& r_node : rModelPart.GetCommunicator().LocalMesh().Nodes())
    {
        const int global_index = r_node.GetValue(INTERFACE_EQUATION_ID);
        const int local_index = rVector.Map().LID(global_index);
        update_fct(r_node, rVariable, rVector[0][local_index]);
    }

    // for (int localIndex = 0; localIndex < rVector.MyLength(); ++localIndex)
    //     std::cout << "Updates | rVector[localIndex]: " << rVector[0][localIndex] << std::endl;
    // rVector.GlobalAssemble(); // I am quite sure this is not needed, since one node is one entry ...
    if (EchoLevel > 2)
        SparseSpaceType::WriteMatrixMarketVector("TrilinosUpdate", rVector);
}


template<class TVarType>
void TInitializeMappingStep(UtilityType::TSystemMatrixType& rMdo,
    UtilityType::TSystemVectorType& rQo,
    UtilityType::TSystemVectorType& rQd,
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    const TVarType& rOriginVariable,
    const TVarType& rDestinationVariable,
    const Kratos::Flags MappingOptions,
    const bool UseTranspose,
    const int EchoLevel)
{
    if (UseTranspose)
        FillSystemVector(rQd, rModelPartDestination, rDestinationVariable, MappingOptions, EchoLevel);
    else
        FillSystemVector(rQo, rModelPartOrigin, rOriginVariable, MappingOptions, EchoLevel);
}

template<>
void UtilityType::InitializeMappingStep(
    TSystemMatrixType& rMdo,
    TSystemVectorType& rQo,
    TSystemVectorType& rQd,
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    const DoubleVariableType& rOriginVariable,
    const DoubleVariableType& rDestinationVariable,
    const Kratos::Flags MappingOptions,
    const bool UseTranspose) const
{
    TInitializeMappingStep(rMdo, rQo, rQd,
                            rModelPartOrigin, rModelPartDestination,
                            rOriginVariable, rDestinationVariable,
                            MappingOptions, UseTranspose, GetEchoLevel());
}

template<>
void UtilityType::InitializeMappingStep(
    TSystemMatrixType& rMdo,
    TSystemVectorType& rQo,
    TSystemVectorType& rQd,
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    const ComponentVariableType& rOriginVariable,
    const ComponentVariableType& rDestinationVariable,
    const Kratos::Flags MappingOptions,
    const bool UseTranspose) const
{
    TInitializeMappingStep(rMdo, rQo, rQd,
                            rModelPartOrigin, rModelPartDestination,
                            rOriginVariable, rDestinationVariable,
                            MappingOptions, UseTranspose, GetEchoLevel());
}

void ExecuteMapping(UtilityType::TSystemMatrixType& rMdo,
                    UtilityType::TSystemVectorType& rQo,
                    UtilityType::TSystemVectorType& rQd,
                    const bool UseTranspose)
{
    if (UseTranspose)
        SparseSpaceType::TransposeMult(rMdo, rQd, rQo); // rQo = rMdo^T * rQo
    else
        SparseSpaceType::Mult(rMdo, rQo, rQd); // rQd = rMdo * rQo
}

template<>
void UtilityType::ExecuteMappingStep(
    TSystemMatrixType& rMdo,
    TSystemVectorType& rQo,
    TSystemVectorType& rQd,
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    const DoubleVariableType& rOriginVariable,
    const DoubleVariableType& rDestinationVariable,
    const Kratos::Flags MappingOptions,
    const bool UseTranspose) const
{
    ExecuteMapping(rMdo, rQo, rQd, UseTranspose);
}

template<>
void UtilityType::ExecuteMappingStep(
    TSystemMatrixType& rMdo,
    TSystemVectorType& rQo,
    TSystemVectorType& rQd,
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    const ComponentVariableType& rOriginVariable,
    const ComponentVariableType& rDestinationVariable,
    const Kratos::Flags MappingOptions,
    const bool UseTranspose) const
{
    ExecuteMapping(rMdo, rQo, rQd, UseTranspose);
}

template<class TVarType>
void TFinalizeMappingStep(UtilityType::TSystemMatrixType& rMdo,
    UtilityType::TSystemVectorType& rQo,
    UtilityType::TSystemVectorType& rQd,
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    const TVarType& rOriginVariable,
    const TVarType& rDestinationVariable,
    const Kratos::Flags MappingOptions,
    const bool UseTranspose,
    const int EchoLevel)
{
    if (UseTranspose)
        Update(rQo, rModelPartOrigin, rOriginVariable, MappingOptions, EchoLevel);
    else
        Update(rQd, rModelPartDestination, rDestinationVariable, MappingOptions, EchoLevel);
}

template<>
void UtilityType::FinalizeMappingStep(
    TSystemMatrixType& rMdo,
    TSystemVectorType& rQo,
    TSystemVectorType& rQd,
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    const DoubleVariableType& rOriginVariable,
    const DoubleVariableType& rDestinationVariable,
    const Kratos::Flags MappingOptions,
    const bool UseTranspose) const
{
    TFinalizeMappingStep(rMdo, rQo, rQd,
                            rModelPartOrigin, rModelPartDestination,
                            rOriginVariable, rDestinationVariable,
                            MappingOptions, UseTranspose, GetEchoLevel());
}

template<>
void UtilityType::FinalizeMappingStep(
    TSystemMatrixType& rMdo,
    TSystemVectorType& rQo,
    TSystemVectorType& rQd,
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    const ComponentVariableType& rOriginVariable,
    const ComponentVariableType& rDestinationVariable,
    const Kratos::Flags MappingOptions,
    const bool UseTranspose) const
{
    TFinalizeMappingStep(rMdo, rQo, rQd,
                            rModelPartOrigin, rModelPartDestination,
                            rOriginVariable, rDestinationVariable,
                            MappingOptions, UseTranspose, GetEchoLevel());
}

/***********************************************************************************/
/* PROTECTED Methods */
/***********************************************************************************/
// template<>
// void UtilityType::ConstructMatrixStructure(MapperLocalSystemPointerVector& rMapperLocalSystems,
//                                            TSystemMatrixType& rMdo) const
// {

// }

/***********************************************************************************/
/* PRIVATE Methods */
/***********************************************************************************/

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class MatrixBasedMappingOperationUtility< SparseSpaceType, DenseSpaceType >;

}  // namespace Kratos.
