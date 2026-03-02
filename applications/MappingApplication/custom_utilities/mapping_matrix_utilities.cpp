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

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "mapping_matrix_utilities.h"
#include "mappers/mapper_define.h"
#include "custom_utilities/mapper_utilities.h"

namespace Kratos {

namespace {

typedef typename MapperDefinitions::SparseSpaceType MappingSparseSpaceType;
typedef typename MapperDefinitions::DenseSpaceType  DenseSpaceType;

typedef MappingMatrixUtilities<MappingSparseSpaceType, DenseSpaceType> MappingMatrixUtilitiesType;

typedef typename MapperLocalSystem::MatrixType MatrixType;
typedef typename MapperLocalSystem::EquationIdVectorType EquationIdVectorType;

typedef std::size_t IndexType;
typedef std::size_t SizeType;

/***********************************************************************************/
/* Functions for internal use in this file */
/***********************************************************************************/

void ConstructMatrixStructure(Kratos::unique_ptr<typename MappingSparseSpaceType::MatrixType>& rpMdo,
                              std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rMapperLocalSystems,
                              const SizeType NumNodesOrigin,
                              const SizeType NumNodesDestination)
{
    // one set for each row storing the corresponding col-IDs
    std::vector<std::unordered_set<IndexType> > indices(NumNodesDestination);

    // preallocate memory for the column indices
    for (IndexType i=0; i<NumNodesDestination; ++i) {
        // TODO I guess this can be optimized...
        // this highly depends on the used mapper => same goes for the Graph in Trilinos
        indices[i].reserve(3);
    }

    EquationIdVectorType origin_ids;
    EquationIdVectorType destination_ids;

    // Looping the local-systems to get the entries for the matrix
    // TODO omp
    for (/*const*/auto& r_local_sys : rMapperLocalSystems) { // TODO I think this can be const bcs it is the ptr
        r_local_sys->EquationIdVectors(origin_ids, destination_ids);
        for (const auto dest_idx : destination_ids) {
            indices[dest_idx].insert(origin_ids.begin(), origin_ids.end());
        }
    }

    // computing the number of non-zero entries
    SizeType num_non_zero_entries = 0;
    for (const auto& r_row_indices : indices) { // looping the indices per row
        num_non_zero_entries += r_row_indices.size(); // adding the number of col-indices
    }

    auto p_Mdo = Kratos::make_unique<typename MappingSparseSpaceType::MatrixType>(
        NumNodesDestination,
        NumNodesOrigin,
        num_non_zero_entries);

    double* p_matrix_values = p_Mdo->value_data().begin();
    IndexType* p_matrix_row_indices = p_Mdo->index1_data().begin();
    IndexType* p_matrix_col_indices = p_Mdo->index2_data().begin();

    // filling the index1 vector - do NOT make parallel!
    p_matrix_row_indices[0] = 0;
    for (IndexType i=0; i<NumNodesDestination; ++i) {
        p_matrix_row_indices[i+1] = p_matrix_row_indices[i] + indices[i].size();
    }

    for (IndexType i=0; i<NumNodesDestination; ++i) {
        const IndexType row_begin = p_matrix_row_indices[i];
        const IndexType row_end = p_matrix_row_indices[i+1];
        IndexType j = row_begin;
        for (const auto index : indices[i]) {
            p_matrix_col_indices[j] = index;
            p_matrix_values[j] = 0.0;
            ++j;
        }

        indices[i].clear(); //deallocating the memory // TODO necessary?

        std::sort(&p_matrix_col_indices[row_begin], &p_matrix_col_indices[row_end]);
    }

    p_Mdo->set_filled(indices.size()+1, num_non_zero_entries);

    rpMdo.swap(p_Mdo);
}

void BuildMatrix(Kratos::unique_ptr<typename MappingSparseSpaceType::MatrixType>& rpMdo,
                 std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rMapperLocalSystems)
{
    MatrixType local_mapping_matrix;
    EquationIdVectorType origin_ids;
    EquationIdVectorType destination_ids;

    for (auto& r_local_sys : rMapperLocalSystems) { // TODO omp

        r_local_sys->CalculateLocalSystem(local_mapping_matrix, origin_ids, destination_ids);

        KRATOS_DEBUG_ERROR_IF(local_mapping_matrix.size1() != destination_ids.size()) << "MappingMatrixAssembly: DestinationID vector size mismatch: LocalMappingMatrix-Size1: " << local_mapping_matrix.size1() << " | DestinationIDs-size: " << destination_ids.size() << std::endl;
        KRATOS_DEBUG_ERROR_IF(local_mapping_matrix.size2() != origin_ids.size()) << "MappingMatrixAssembly: OriginID vector size mismatch: LocalMappingMatrix-Size2: " << local_mapping_matrix.size2() << " | OriginIDs-size: " << origin_ids.size() << std::endl;

        for (IndexType i=0; i<destination_ids.size(); ++i) {
            for (IndexType j=0; j<origin_ids.size(); ++j) {
                // #pragma omp atomic
                (*rpMdo)(destination_ids[i], origin_ids[j]) += local_mapping_matrix(i,j);
            }
        }

        r_local_sys->Clear();
    }
}

} // anonymous namespace

template<>
void MappingMatrixUtilitiesType::CheckRowSum(
    const typename MappingSparseSpaceType::MatrixType& rM,
    const std::string& rBaseFileName,
    const bool ThrowError,
    const double Tolerance)
{
    typename MappingSparseSpaceType::VectorType unit_vector(MappingSparseSpaceType::Size2(rM));
    MappingSparseSpaceType::Set(unit_vector, 1.0);

    typename MappingSparseSpaceType::VectorType row_sums_vector(MappingSparseSpaceType::Size1(rM));

    MappingSparseSpaceType::Mult(rM, unit_vector, row_sums_vector);

    bool write_mm_file = false;
    for (std::size_t i = 0; i < MappingSparseSpaceType::Size(row_sums_vector); ++i) {
        if (std::abs(row_sums_vector[i] - 1.0) > Tolerance) {
            KRATOS_WARNING("MappingMatrixAssembly") << "The row sum in row " << i << " is unequal 1.0: " << row_sums_vector[i] << std::endl;
            write_mm_file = true;
        }
    }

    if (write_mm_file) {
        MappingSparseSpaceType::WriteMatrixMarketVector(("RowSumVector_" + rBaseFileName).c_str(), row_sums_vector);
        KRATOS_ERROR_IF(ThrowError) << "Mapping matrix does not sum to unity. Please check file " << rBaseFileName << " in your project directory for row sums\n";
    }
}

template<>
void MappingMatrixUtilitiesType::InitializeSystemVector(
    Kratos::unique_ptr<typename MappingSparseSpaceType::VectorType>& rpVector,
    const std::size_t VectorSize)
{
    // The vectors dont have graphs, that why we don't always have to reinitialize them
    if (rpVector == nullptr || rpVector->size() != VectorSize) { //if the pointer is not initialized initialize it to an empty vector
        Kratos::unique_ptr<typename MappingSparseSpaceType::VectorType> p_new_vector = Kratos::make_unique<typename MappingSparseSpaceType::VectorType>(VectorSize);
        rpVector.swap(p_new_vector);

        // TODO do I also have to set to zero the contents?
    }
    else {
        MappingSparseSpaceType::SetToZero(*rpVector);
    }
}

template<>
void MappingMatrixUtilitiesType::BuildMappingMatrix(
    Kratos::unique_ptr<typename MappingSparseSpaceType::MatrixType>& rpMappingMatrix,
    Kratos::unique_ptr<typename MappingSparseSpaceType::VectorType>& rpInterfaceVectorOrigin,
    Kratos::unique_ptr<typename MappingSparseSpaceType::VectorType>& rpInterfaceVectorDestination,
    const ModelPart& rModelPartOrigin,
    const ModelPart& rModelPartDestination,
    std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rMapperLocalSystems,
    const int EchoLevel)
{
    KRATOS_TRY

    static_assert(!MappingSparseSpaceType::IsDistributed(), "Using a distributed Space!");

    const SizeType num_nodes_origin = rModelPartOrigin.NumberOfNodes();
    const SizeType num_nodes_destination = rModelPartDestination.NumberOfNodes();

    // Initialize the Matrix
    // This has to be done always since the Graph has changed if the Interface is updated!
    ConstructMatrixStructure(rpMappingMatrix, rMapperLocalSystems,
                             num_nodes_origin, num_nodes_destination);

    BuildMatrix(rpMappingMatrix, rMapperLocalSystems);

    MappingMatrixUtilitiesType::InitializeSystemVector(rpInterfaceVectorOrigin, num_nodes_origin);
    MappingMatrixUtilitiesType::InitializeSystemVector(rpInterfaceVectorDestination, num_nodes_destination);

    KRATOS_CATCH("")
}

template<>
void MappingMatrixUtilitiesType::BuildMappingMatrixRBFMapper(
    Kratos::unique_ptr<typename MappingSparseSpaceType::MatrixType>& rpMappingMatrix,
    Kratos::unique_ptr<typename MappingSparseSpaceType::VectorType>& rpInterfaceVectorOrigin,
    Kratos::unique_ptr<typename MappingSparseSpaceType::VectorType>& rpInterfaceVectorDestination,
    const ModelPart& rModelPartOrigin,
    const ModelPart& rModelPartDestination,
    std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rMapperLocalSystems,
    const IndexType NumberOfPolynomialTerms,
    const bool BuildOriginInterpolationMatrix,
    const bool OriginIsIga,
    const int EchoLevel)
{
    KRATOS_TRY

    static_assert(!MappingSparseSpaceType::IsDistributed(), "Using a distributed Space!");

    const SizeType num_nodes_origin = rModelPartOrigin.NumberOfNodes();
    const SizeType num_conditions_origin = rModelPartOrigin.NumberOfConditions();
    const SizeType num_nodes_destination = rModelPartDestination.NumberOfNodes();
    const SizeType num_conditions_destination = rModelPartDestination.NumberOfConditions();

    IndexType origin_size;
    IndexType destination_size;

    if (!OriginIsIga){
        origin_size = num_nodes_origin;
        destination_size = num_nodes_destination;
    } else if (OriginIsIga && BuildOriginInterpolationMatrix) {
        origin_size = num_conditions_origin;
        destination_size = num_conditions_destination;
    } else if (OriginIsIga && !BuildOriginInterpolationMatrix) {
        origin_size = num_conditions_origin;
        destination_size = num_nodes_destination;
    }
    
    // Initialize the Matrix
    // This has to be done always since the Graph has changed if the Interface is updated!
    if (BuildOriginInterpolationMatrix){
        ConstructMatrixStructure(rpMappingMatrix, rMapperLocalSystems,
                                origin_size + NumberOfPolynomialTerms, destination_size + NumberOfPolynomialTerms);
    } else {
        ConstructMatrixStructure(rpMappingMatrix, rMapperLocalSystems,
                                origin_size + NumberOfPolynomialTerms, destination_size);
    }

    BuildMatrix(rpMappingMatrix, rMapperLocalSystems);

    MappingMatrixUtilitiesType::InitializeSystemVector(rpInterfaceVectorOrigin, num_nodes_origin);
    MappingMatrixUtilitiesType::InitializeSystemVector(rpInterfaceVectorDestination, num_nodes_destination);

    KRATOS_CATCH("")
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class MappingMatrixUtilities< MappingSparseSpaceType, DenseSpaceType >;

}  // namespace Kratos.
