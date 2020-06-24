/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Author: Manuel Messmer
*/

// System includes

// External includes
#include <Eigen/Core>
#include <Eigen/Dense>

// Project includes
#include "custom_processes/perturb_geometry_subgrid_process.h"
#include "utilities/builtin_timer.h"

#include "utilities/openmp_utils.h"

namespace Kratos
{

typedef ModelPart::NodesContainerType::ContainerType                 ResultNodesContainerType;

/**
 * @brief Creates Eigenvectors of correlation matrix in a subgrid
 * @details Finds a subgrid (coarser mesh). Generates correlation matrix in subgrid. Decomposes correlation matrix.
 * @param correlation_matrix Correlation matrix. Stores correlation value for all nodes in the subgrid.
 * @param rPerturbationMatrix Perturbation matrix. Stores eigenvectors of correlation matrix.
 */
int PerturbGeometrySubgridProcess::CreateEigenvectors( ModelPart& rThisModelPart, double minDistance, double correlationLength, double truncationTolerance ){
    KRATOS_TRY;

    int num_of_nodes = rThisModelPart.NumberOfNodes();
    // Search radius. Defines the minimum distance among the nodes in the subgrid
    double radius = minDistance;

    // Get all nodes
    ModelPart::NodesContainerType nodes = rThisModelPart.Nodes();
    // Construct and initialize searcher classs
    OMP_NodeSearch searcher;
    searcher.InitializeSearch(nodes);

    // Define the reduced space vector
    std::vector<ModelPart::NodeIterator> reduced_space_nodes;
    ResultNodesContainerType  results;

    BuiltinTimer reduced_space_timer;
    // Mark all nodes as unvisited
    #pragma omp parallel
    {
        const auto it_node_begin = rThisModelPart.NodesBegin();
        #pragma omp for
        for (int i = 0; i < num_of_nodes; i++){
            auto it_node = it_node_begin + i;
            it_node->Set(VISITED,false);
        }
    }

    // Generate reduced space
    for (ModelPart::NodeIterator it_node =rThisModelPart.NodesBegin(); it_node != rThisModelPart.NodesEnd(); it_node++)
    {
        if( !it_node->Is(VISITED) ) {
            it_node->Set(VISITED,true);
            reduced_space_nodes.push_back(it_node);
            results = {};
            searcher.SearchNodesInRadiusExclusiveImplementation(nodes,it_node->GetId()-1,radius,results);
            for( size_t i = 0; i < results.size(); i++ ){
                results[i]->Set(VISITED,true);
            }
        }
    }
    const int num_nodes_reduced_space = reduced_space_nodes.size();

    KRATOS_INFO_IF("PerturbGeometrySubgridProcess: Find Reduced Space Time", mEchoLevel > 0)
        << reduced_space_timer.ElapsedSeconds() << std::endl
        << "Number of Nodes in Reduced Space: " << num_nodes_reduced_space << " / " << num_of_nodes << std::endl;

    // Construct and initialize correlation matrix
    BuiltinTimer build_cl_matrix_timer;
    Eigen::MatrixXd correlation_matrix;
    correlation_matrix.resize(num_nodes_reduced_space,num_nodes_reduced_space);
    // Assemble correlation matrix
    #pragma omp parallel
    {
        const auto it_node_begin = reduced_space_nodes.begin();
        #pragma omp for
        for( int row_counter = 0; row_counter < num_nodes_reduced_space; row_counter++)
        {
            auto it_node = it_node_begin + row_counter;
            for( int column_counter = 0; column_counter < num_nodes_reduced_space; column_counter++)
            {
                auto it_node_2 = it_node_begin + column_counter;
                correlation_matrix(row_counter ,column_counter ) = CorrelationFunction( *it_node, *it_node_2, correlationLength);
            }
        }
    }
    KRATOS_INFO_IF("PerturbGeometrySubgridProcess: Build Correlation Matrix Time", mEchoLevel > 0)
        << build_cl_matrix_timer.ElapsedSeconds() << std::endl;

    // Construct eigensolver and solve eigenproblem
    BuiltinTimer eigensolver_time;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(num_nodes_reduced_space);
    es.compute(correlation_matrix,Eigen::ComputeEigenvectors);
    KRATOS_INFO_IF("PerturbGeometrySubgridProcess: Eigensolver Time", mEchoLevel > 0)
            << eigensolver_time.ElapsedSeconds() << std::endl;

    if( es.info() != Eigen::Success)
    {
        KRATOS_INFO("PerturbGeometrySubgridProcess: Eigensolver")
        << "Decomposition failed" << std::endl;
    }

    // Find number of required eigenvalues to statisfy convergence criterion
    // Eigenvalues are sorted in ascending order and are normalized to an euclidean length of one!!
    double total_sum_eigenvalues = 0.0;
    double reduced_sum_eigenvalues = 0.0;
    int num_eigenvalues_required = 0;

    for( int i = es.eigenvalues().size() - 1; i >= 0; i--)
    {
        total_sum_eigenvalues += es.eigenvalues()(i);
    }
    for( int i = es.eigenvalues().size() - 1; i >= 0; i--)
    {
        reduced_sum_eigenvalues += es.eigenvalues()(i);
        num_eigenvalues_required++;
        if( reduced_sum_eigenvalues > (1 - truncationTolerance)*total_sum_eigenvalues)
        {
            KRATOS_INFO_IF("PerturbGeometrySubgridProcess", mEchoLevel > 0)
                << "Truncation Error (" <<  truncationTolerance
                << ") is achieved with " << num_eigenvalues_required << " Eigenvalues" << std::endl;
            break;
        }
    }

    int num_random_variables = num_eigenvalues_required;

    // Retrieve eigenvector matrix
    Eigen::MatrixXd eigenvectors = es.eigenvectors();
    // Truncate eigenvector matrix
    // Attention: eigenvectors are sorted according to ascending order of eigenvalues
    int cut_off = num_nodes_reduced_space - num_random_variables;
    Eigen::MatrixXd eigenvectors_truncated = Eigen::Map<Eigen::MatrixXd,0,Eigen::OuterStride<> >(
    eigenvectors.real().data() + cut_off , es.eigenvectors().real().rows(), num_random_variables, Eigen::OuterStride<>(es.eigenvectors().outerStride()) ) ;

    // Retrieve and truncate eigenvalue vector
    Eigen::VectorXd eigenvalues = es.eigenvalues().real().tail(num_random_variables);

    // Get and resize final displacement matrix
    DenseMatrixType& rPerturbationMatrix = *mpPerturbationMatrix;
    rPerturbationMatrix.resize(num_of_nodes,num_random_variables);

    // Generate random field over full domain
    BuiltinTimer assemble_random_field_time;
    #pragma omp parallel
    {
        const auto it_node_begin = rThisModelPart.NodesBegin();

        #pragma omp for
        for( int i = 0; i < num_of_nodes; i++){

            Eigen::VectorXd correlation_vector = Eigen::RowVectorXd::Zero(num_nodes_reduced_space);
            auto it_node = it_node_begin + i;
            const auto it_node_reduced_begin = reduced_space_nodes.begin();
            // Assemble correlation vector
            for( int j = 0; j < num_nodes_reduced_space; j++){
                auto it_node_reduced = it_node_reduced_begin + j;
                correlation_vector(j) = CorrelationFunction( it_node, *it_node_reduced, correlationLength);
            }
            // Assemble perturbation field
            for( int j = 0; j < num_random_variables; j++){
                rPerturbationMatrix(i,j) = sqrt( 1/eigenvalues(j) ) * (correlation_vector).dot(eigenvectors_truncated.col(j));
            }
        }
    }
    KRATOS_INFO_IF("PerturbGeometrySubgridProcess: Assemble Random Field Time", mEchoLevel > 0)
            << assemble_random_field_time.ElapsedSeconds() << std::endl;


    return num_random_variables;

    KRATOS_CATCH("");
}

} // namespace Kratos
