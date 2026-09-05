// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Manuel Messmer
//

// System includes
#include <algorithm> // std::sort, std::unique

// External includes

// Project includes
#include "custom_utilities/perturb_geometry_sparse_utility.h"
#include "custom_utilities/node_search_utility.h"
#include "utilities/builtin_timer.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{

int PerturbGeometrySparseUtility::CreateRandomFieldVectors(){
    KRATOS_TRY;

    int number_of_nodes = mrInitialModelPart.NumberOfNodes();

    BuiltinTimer assemble_cl_matrix_timer;

    // Construct and initialize required containers
    ModelPart::NodesContainerType nodes = mrInitialModelPart.Nodes();
    ResultNodesContainerType  results;

    // Construct and initialize searcher
    NodeSearchUtility searcher(nodes);

    // Collect the correlation entries per row, then write the CSR arrays
    // directly: element-by-element operator() insertion is a uBLAS-only idiom
    // and is not valid to build up the Eigen backend matrix.
    int counter = 0;
    const auto it_node_begin = mrInitialModelPart.NodesBegin();

    std::vector<std::vector<std::pair<int, double>>> row_entries(number_of_nodes);
    for( int i = 0; i < number_of_nodes; i++){
        auto it_node = it_node_begin + i;
        // Add current nodes to results
        results = { mrInitialModelPart.Nodes().GetContainer()[i] };
        // Find all neighbouring nodes
        searcher.SearchNodesInRadius(&*it_node, 3.0*mCorrelationLength, results);
        for( size_t j = 0; j < results.size(); j++){
            counter++;
            int index = results[j]->GetId() - 1;
            auto it_node_results = it_node_begin + index;

            row_entries[i].emplace_back(index, CorrelationFunction( *it_node, *it_node_results, mCorrelationLength));
        }
        auto& r_row = row_entries[i];
        std::sort(r_row.begin(), r_row.end());
        // The search may report a node more than once (e.g. the row node is
        // seeded and found again); the assignment in the old element-wise
        // build silently collapsed those duplicates.
        r_row.erase(std::unique(r_row.begin(), r_row.end(),
                                [](const auto& rA, const auto& rB){ return rA.first == rB.first; }),
                    r_row.end());
    }

    std::size_t nnz = 0;
    for( const auto& r_row : row_entries)
        nnz += r_row.size();

    SparseMatrixType correlation_matrix(number_of_nodes, number_of_nodes, nnz);
    {
        auto row_ptr = correlation_matrix.index1_data().begin();
        auto col_ptr = correlation_matrix.index2_data().begin();
        auto val_ptr = correlation_matrix.value_data().begin();
        row_ptr[0] = 0;
        std::size_t k = 0;
        for( int i = 0; i < number_of_nodes; i++){
            for( const auto& r_entry : row_entries[i]){
                col_ptr[k] = r_entry.first;
                val_ptr[k] = r_entry.second;
                ++k;
            }
            row_ptr[i + 1] = k;
        }
        correlation_matrix.set_filled(number_of_nodes + 1, nnz);
    }

    KRATOS_INFO_IF("PerturbGeometrySparseUtility: Build Correlation Matrix Time", mEchoLevel > 0)
            << assemble_cl_matrix_timer.ElapsedSeconds() << std::endl
            << (double) counter / ((double)number_of_nodes*(double)number_of_nodes) * 100.0 << "% of Matrix is populated" << std::endl;

    DenseMatrixType eigenvectors;
    DenseVectorType eigenvalues;

    // Build the identity by writing the CSR arrays directly (valid for both backends)
    SparseMatrixType identity_matrix(number_of_nodes, number_of_nodes, number_of_nodes);
    {
        auto row_ptr = identity_matrix.index1_data().begin();
        auto col_ptr = identity_matrix.index2_data().begin();
        auto val_ptr = identity_matrix.value_data().begin();
        row_ptr[0] = 0;
        for( int i = 0; i < number_of_nodes; i++){
            col_ptr[i] = i;
            val_ptr[i] = 1.0;
            row_ptr[i + 1] = i + 1;
        }
        identity_matrix.set_filled(number_of_nodes + 1, number_of_nodes);
    }

    mpEigenSolver->Solve( identity_matrix,
                          correlation_matrix,
                          eigenvalues,
                          eigenvectors);

    // Find number of required eigenvalues to statisfy convergence criterion
    double sum_eigenvalues =  1 / eigenvalues(0);
    int number_required_eigenvalues = 0;

    for( size_t i = 1; i < eigenvalues.size(); i++){
        double epsilon = 1 - sum_eigenvalues / ( sum_eigenvalues + 1 / eigenvalues(i) );
        if( epsilon < mTruncationError){
            number_required_eigenvalues = i + 1;
            KRATOS_INFO_IF("PerturbGeometrySparseUtility", mEchoLevel > 0)
                << "Truncation Error  (" << mTruncationError
                << ") is achieved with " << number_required_eigenvalues << " Eigenvalues" << std::endl;

            break;
        }
        else if( i == eigenvalues.size() - 1){
            number_required_eigenvalues = eigenvalues.size();
            KRATOS_INFO_IF("PerturbGeometrySparseUtility", mEchoLevel > 0)
                << "Truncation Error is NOT achieved:  " << epsilon << " / " << mTruncationError << std::endl
                << "Maximum number of computed eigenvalues is used: " << number_required_eigenvalues << std::endl;
        }
        sum_eigenvalues += 1/eigenvalues(i);
    }

    // Normalize required eigenvectors
    for( int i =  0; i < number_required_eigenvalues; i++){
        double eucledian_norm =  norm_2( row(eigenvectors,i) );
        row(eigenvectors,i) = 1.0 / eucledian_norm * row(eigenvectors,i);
    }

    int number_of_random_variables = number_required_eigenvalues;

    // Construct and initialize perturbation matrix
    DenseMatrixType& rPerturbationMatrix = *mpPerturbationMatrix;
    rPerturbationMatrix.resize(number_of_nodes,number_of_random_variables);

    // Assemble perturbation matrix. The product with a row of the sparse
    // correlation matrix is accumulated over its CSR arrays (a sparse uBLAS
    // row proxy does not exist for the Eigen backend matrix).
    BuiltinTimer assemble_random_field_time;
    const auto& corr_row_ptr = correlation_matrix.index1_data();
    const auto& corr_col_idx = correlation_matrix.index2_data();
    const auto& corr_values = correlation_matrix.value_data();
    IndexPartition<unsigned int>(number_of_nodes).for_each(
        [number_of_random_variables, &eigenvalues, &eigenvectors, &corr_row_ptr, &corr_col_idx, &corr_values, &rPerturbationMatrix](unsigned int i){
            for( int j = 0; j < number_of_random_variables; j++){
                double dot = 0.0;
                for( auto k = corr_row_ptr[i]; k < corr_row_ptr[i + 1]; ++k){
                    dot += corr_values[k] * eigenvectors(j, corr_col_idx[k]);
                }
                rPerturbationMatrix(i,j) = std::sqrt(eigenvalues(j)) * dot;
            }
        }
    );

    KRATOS_INFO_IF("PerturbGeometrySparseUtility: Assemble Random Field Time", mEchoLevel > 0)
        << assemble_random_field_time.ElapsedSeconds() << std::endl;

    return number_of_random_variables;

    KRATOS_CATCH("")

}

} // namespace Kratos