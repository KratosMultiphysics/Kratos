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

// External includes

// Project includes
#include "custom_utilities/perturb_geometry_subgrid_utility.h"
#include "custom_utilities/node_search_utility.h"
#include "utilities/builtin_timer.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"

namespace Kratos
{

int PerturbGeometrySubgridUtility::CreateRandomFieldVectors(){
    KRATOS_TRY;

    int num_of_nodes = mrInitialModelPart.NumberOfNodes();
    // Search radius. Defines the minimum distance among the nodes in the subgrid
    double radius = mMinDistanceSubgrid;

    // Get all nodes
    ModelPart::NodesContainerType nodes = mrInitialModelPart.Nodes();
    // Construct and initialize searcher class
    NodeSearchUtility searcher(nodes);

    // Define the reduced space vector
    std::vector<ModelPart::NodeType::Pointer> reduced_space_nodes;
    ResultNodesContainerType  results;

    BuiltinTimer reduced_space_timer;
    // Mark all nodes as unvisited
    VariableUtils().SetFlag(VISITED, false, mrInitialModelPart.Nodes());

    // Generate reduced space
    // It is assumed that the reduced space contains less the 50% of all ndoes
    reduced_space_nodes.reserve((int)0.5*num_of_nodes);
    const auto it_node_begin = mrInitialModelPart.NodesBegin();
    for(int i = 0; i < num_of_nodes; ++i){
        auto it_node = it_node_begin + i;
        if(!it_node->Is(VISITED)) {
            it_node->Set(VISITED,true);
            reduced_space_nodes.push_back(&*it_node);
            results = {};
            searcher.SearchNodesInRadius(&*it_node, radius, results);
            for(std::size_t j = 0; j < results.size(); ++j){
                results[j]->Set(VISITED,true);
            }
        }
    }
    const int num_nodes_reduced_space = reduced_space_nodes.size();

    KRATOS_INFO_IF("PerturbGeometrySubgridUtility: Find Reduced Space Time", mEchoLevel > 0)
        << reduced_space_timer.ElapsedSeconds() << std::endl
        << "Number of Nodes in Reduced Space: " << num_nodes_reduced_space << " / " << num_of_nodes << std::endl;

    // Construct and initialize correlation matrix
    BuiltinTimer build_cl_matrix_timer;
    DenseMatrixType correlation_matrix;
    correlation_matrix.resize(num_nodes_reduced_space,num_nodes_reduced_space);
    // Assemble correlation matrix
    const auto it_node_reduced_begin = reduced_space_nodes.begin();
    IndexPartition<unsigned int>(num_nodes_reduced_space).for_each(
        [this, num_nodes_reduced_space, it_node_reduced_begin, &correlation_matrix](unsigned int row_counter){
            auto it_node = it_node_reduced_begin + row_counter;
            for( int column_counter = 0; column_counter < num_nodes_reduced_space; column_counter++){
                auto it_node_2 = it_node_reduced_begin + column_counter;
                correlation_matrix(row_counter, column_counter) = CorrelationFunction(**it_node, **it_node_2, this->mCorrelationLength);
            }
        }
    );

    KRATOS_INFO_IF("PerturbGeometrySubgridUtility: Build Correlation Matrix Time", mEchoLevel > 0)
        << build_cl_matrix_timer.ElapsedSeconds() << std::endl;

    // Construct eigensolver and solve eigenproblem
    DenseVectorType eigenvalues; // Vector is resized inside the solver
    DenseMatrixType eigenvectors; // Matrix is resized inside the solver
    DenseMatrixType dummy;
    mpEigenSolver->Solve(correlation_matrix, dummy, eigenvalues, eigenvectors);

    // Find number of required eigenvalues to statisfy convergence criterion
    // Eigenvalues are sorted in ascending order and are normalized to an euclidean length of one!!
    double total_sum_eigenvalues = 0.0;
    double reduced_sum_eigenvalues = 0.0;
    int num_eigenvalues_required = 0;

    for( std::size_t i = 0; i < eigenvalues.size(); ++i){
        total_sum_eigenvalues += eigenvalues(i);
    }
    for( std::size_t i = 0; i < eigenvalues.size(); ++i){
        reduced_sum_eigenvalues += eigenvalues(i);
        num_eigenvalues_required++;
        if( reduced_sum_eigenvalues > (1 - mTruncationError)*total_sum_eigenvalues){
            KRATOS_INFO_IF("PerturbGeometrySubgridUtility", mEchoLevel > 0)
                << "Truncation Error (" <<  mTruncationError
                << ") is achieved with " << num_eigenvalues_required << " Eigenvalues" << std::endl;
            break;
        }
    }

    int num_random_variables = num_eigenvalues_required;

    // Get and resize final perturbation matrix
    DenseMatrixType& rPerturbationMatrix = *mpPerturbationMatrix;
    rPerturbationMatrix.resize(num_of_nodes, num_random_variables);

    // Construct Correlation vector
    DenseVectorType correlation_vector;
    correlation_vector.resize(num_nodes_reduced_space);
    // Generate random field over full domain
    BuiltinTimer assemble_random_field_time;
    #pragma omp parallel for firstprivate(correlation_vector)
    for( int i = 0; i < num_of_nodes; i++){
        auto it_node = it_node_begin + i;
        // Assemble correlation vector
        for( int j = 0; j < num_nodes_reduced_space; j++){
            auto it_node_reduced = it_node_reduced_begin + j;
            correlation_vector(j) = CorrelationFunction( *it_node, **it_node_reduced, mCorrelationLength);
        }
        // Assemble perturbation field
        for( int j = 0; j < num_random_variables; j++){
            rPerturbationMatrix(i,j) = std::sqrt( 1.0/eigenvalues(j) ) * inner_prod(correlation_vector, column(eigenvectors,j) );
        }
    }

    KRATOS_INFO_IF("PerturbGeometrySubgridUtility: Assemble Random Field Time", mEchoLevel > 0)
            << assemble_random_field_time.ElapsedSeconds() << std::endl;

    return num_random_variables;

    KRATOS_CATCH("");
}

} // namespace Kratos