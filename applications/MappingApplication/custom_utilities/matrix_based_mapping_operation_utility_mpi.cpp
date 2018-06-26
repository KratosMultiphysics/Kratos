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


    typedef MapperDefinitions::MPISparseSpaceType SparseSpaceType;

    typedef typename SparseSpaceType::MatrixType SystemMatrixType; // Epetra_FECrsMatrix
    typedef typename SparseSpaceType::MatrixPointerType SystemMatrixPointerType;

    typedef typename SparseSpaceType::VectorType SystemVectorType; // Epetra_FEVector
    typedef typename SparseSpaceType::VectorPointerType SystemVectorPointerType;

    SystemMatrixPointerType mpMappingMatrix;
    SystemVectorPointerType mpVectorOrigin;
    SystemVectorPointerType mpVectorDestination;

    // big TODO what if the rank doesn't have local nodes ... ?

    const int num_local_nodes_orig = rModelPartOrigin.GetCommunicator().LocalMesh().NumberOfNodes();
    const int num_local_nodes_dest = rModelPartDestination.GetCommunicator().LocalMesh().NumberOfNodes();

    const Epetra_MpiComm epetra_comm(MPI_COMM_WORLD); // TODO do I have to save this as member???

    std::vector<int> global_elements_orig(num_local_nodes_orig);
    std::vector<int> global_elements_dest(num_local_nodes_dest);

    const auto nodes_begin_orig = rModelPartOrigin.GetCommunicator().LocalMesh().NodesBegin();
    #pragma omp parallel for
    for (int i=0; i<num_local_nodes_orig; ++i)
        global_elements_orig[i] = ( nodes_begin_orig + i )->GetValue(INTERFACE_EQUATION_ID);

    const auto nodes_begin_dest = rModelPartDestination.GetCommunicator().LocalMesh().NodesBegin();
    #pragma omp parallel for
    for (int i=0; i<num_local_nodes_dest; ++i)
        global_elements_dest[i] = ( nodes_begin_dest + i )->GetValue(INTERFACE_EQUATION_ID);

    // Epetra_Map (long long NumGlobalElements, int NumMyElements, const long long *MyGlobalElements, int IndexBase, const Epetra_Comm &Comm)

    const int num_global_elements = -1; // this means its gonna be computed by Epetra_Map
    const int index_base = 0; // for C/C++

    Epetra_Map epetra_map_cols(num_global_elements,
                                num_local_nodes_orig,
                                global_elements_orig.data(), // taken as const
                                index_base,
                                epetra_comm);

    Epetra_Map epetra_map_rows(num_global_elements,
                                num_local_nodes_dest,
                                global_elements_dest.data(), // taken as const
                                index_base,
                                epetra_comm);


    // explanation in here: https://trilinos.org/docs/dev/packages/epetra/doc/html/classEpetra__CrsGraph.html
    const int num_indices_per_row = 25; // TODO this is to be tested => set to zero maybe ...

    // TODO do I even need the graph? I think I could directly use the Matrix and perform the same operations ...
    // TODO I should construct the graph with two maps! => one for row and one for column
    // Performance optimization see https://trilinos.org/docs/dev/packages/epetra/doc/html/classEpetra__CrsGraph.html
    Epetra_FECrsGraph epetra_graph(Epetra_DataAccess::Copy,
                                   epetra_map_rows,
                                   epetra_map_cols,
                                   num_indices_per_row);

    EquationIdVectorType origin_ids;
    EquationIdVectorType destination_ids;

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

    int ierr = epetra_graph.GlobalAssemble(); // TODO check if it should call "FillComplete"
    KRATOS_ERROR_IF( ierr != 0 ) << "Epetra failure in Epetra_FECrsGraph.GlobalAssemble. "
        << "Error code: " << ierr << std::endl;

    // // TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(Copy,Agraph) );
    // // https://trilinos.org/docs/dev/packages/epetra/doc/html/Epetra__DataAccess_8h.html#ad1a985e79f94ad63030815a0d7d90928
    SystemMatrixPointerType p_new_mapping_matrix = Kratos::make_shared<SystemMatrixType>(Epetra_DataAccess::Copy, epetra_graph);
    mpMappingMatrix.swap(p_new_mapping_matrix);

    SystemVectorPointerType p_new_vector_destination = Kratos::make_shared<SystemVectorType>(epetra_map_rows);
    mpVectorDestination.swap(p_new_vector_destination);

    SystemVectorPointerType p_new_vector_origin = Kratos::make_shared<SystemVectorType>(epetra_map_cols);
    mpVectorOrigin.swap(p_new_vector_origin);

    std::cout << "AFTER Trilinos" << std::cout;

    Philipp check if you assigned the rows/colums correctly!!!

}

// The "Build" function
template<>
void UtilityType::BuildMappingMatrix(
    const MapperLocalSystemPointerVector& rMapperLocalSystems,
    TSystemMatrixType& rMdo) const
{
    MappingWeightsVector mapping_weights;

    EquationIdVectorType origin_ids;
    EquationIdVectorType destination_ids;

    for (auto& rp_local_sys : rMapperLocalSystems)
    {
        rp_local_sys->CalculateLocalSystem(mapping_weights, origin_ids, destination_ids);

        KRATOS_DEBUG_ERROR_IF(mapping_weights.size() != origin_ids.size())
            << "OriginID vector size mismatch" << std::endl;
        KRATOS_DEBUG_ERROR_IF(mapping_weights.size() != destination_ids.size())
            << "DestinationID vector size mismatch" << std::endl;

        if (mapping_weights.size() > 0)
        {
            const int ierr = rMdo.SumIntoGlobalValues(
                destination_ids.size(), destination_ids.data(),
                origin_ids.size(),      origin_ids.data(),
                mapping_weights.data());

            KRATOS_ERROR_IF( ierr < 0 ) << "Epetra failure in Epetra_FECrsMatrix.SumIntoGlobalValues. "
                << "Error code: " << ierr << std::endl;
        }

        rp_local_sys->Clear();
    }
}

template< class TVarType >
void FillSystemVector(UtilityType::TSystemVectorType& rVector,
                        ModelPart& rModelPart,
                        const TVarType& rVariable,
                        const Kratos::Flags& rMappingOptions)
{
    // Here we construct a function pointer to not have the if all the time inside the loop
    const auto fill_fct = MapperUtilities::GetFillFunction<TVarType>(rMappingOptions);

    for (const auto& r_node : rModelPart.GetCommunicator().LocalMesh().Nodes())
    {
        const int global_index = r_node.GetValue(INTERFACE_EQUATION_ID);
        const int local_index = rVector.Map().LID(global_index);
        fill_fct(r_node, rVariable, *rVector[local_index]);
    }
}

template< class TVarType >
void Update(UtilityType::TSystemVectorType& rVector,
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

    for (auto& r_node : rModelPart.GetCommunicator().LocalMesh().Nodes())
    {
        const int global_index = r_node.GetValue(INTERFACE_EQUATION_ID);
        const int local_index = rVector.Map().LID(global_index);
        update_fct(r_node, rVariable, *rVector[local_index]);
    }
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
    const bool UseTranspose)
{
    if (UseTranspose)
        FillSystemVector(rQd, rModelPartDestination, rDestinationVariable, MappingOptions);
    else
        FillSystemVector(rQo, rModelPartOrigin, rOriginVariable, MappingOptions);
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
                            MappingOptions, UseTranspose);
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
                            MappingOptions, UseTranspose);
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
    const bool UseTranspose)
{
    if (UseTranspose)
        Update(rQo, rModelPartOrigin, rOriginVariable, MappingOptions);
    else
        Update(rQd, rModelPartDestination, rDestinationVariable, MappingOptions);
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
                            MappingOptions, UseTranspose);
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
                            MappingOptions, UseTranspose);
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
