// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Manuel Messmer
//
//

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

int PerturbGeometrySubgridProcess::CreateEigenvectors( ModelPart& rThisModelPart, double minDistance, double correlationLength, double truncationTolerance ){
    KRATOS_TRY;

    int num_of_nodes = rThisModelPart.NumberOfNodes();
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
    for (ModelPart::NodeIterator itNode =rThisModelPart.NodesBegin(); itNode != rThisModelPart.NodesEnd(); itNode++)
    {
        if( !itNode->Is(VISITED) ) {
            itNode->Set(VISITED,true);
            reduced_space_nodes.push_back(itNode);
            results = {};
            searcher.SearchNodesInRadiusExclusiveImplementation(nodes,itNode->GetId()-1,radius,results);
            for( int i = 0; i < results.size(); i++ ){
                results[i]->Set(VISITED,true);
            }
        }
    }
    const int num_nodes_reduced_space = reduced_space_nodes.size();

    KRATOS_INFO("PerturbGeometrySubgridProcess: Find Reduced Space Time")
        << reduced_space_timer.ElapsedSeconds() << std::endl
        << "Number of Nodes in Reduced Space: " << num_nodes_reduced_space << " / " << num_of_nodes << std::endl;

    // Construct and initialize correlation matrix
    BuiltinTimer build_cl_matrix_timer;
    Eigen::MatrixXd CorrelationMatrix;
    CorrelationMatrix.resize(num_nodes_reduced_space,num_nodes_reduced_space);
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
                CorrelationMatrix(row_counter ,column_counter ) = CorrelationFunction( *it_node, *it_node_2, correlationLength);
            }
        }
    }
    KRATOS_INFO("PerturbGeometrySubgridProcess: Build Correlation Matrix Time")
        << build_cl_matrix_timer.ElapsedSeconds() << std::endl;

    // Construct eigensolver and solve eigenproblem
    BuiltinTimer eigensolver_time;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(num_nodes_reduced_space);
    es.compute(CorrelationMatrix,Eigen::ComputeEigenvectors);
    KRATOS_INFO("PerturbGeometrySubgridProcess: Eigensolver Time")
            << eigensolver_time.ElapsedSeconds() << std::endl;

    if( es.info() != Eigen::Success)
    {
        KRATOS_INFO("PerturbGeometrySubgridProcess")
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
            KRATOS_INFO("PerturbGeometrySubgridProcess")
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
    Eigen::MatrixXd Eigenvectors = Eigen::Map<Eigen::MatrixXd,0,Eigen::OuterStride<> >(
    eigenvectors.real().data() + cut_off , es.eigenvectors().real().rows(), num_random_variables, Eigen::OuterStride<>(es.eigenvectors().outerStride()) ) ;

    // Retrieve and truncate eigenvalue vector
    Eigen::VectorXd Eigenvalues = es.eigenvalues().real().tail(num_random_variables);

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

            Eigen::VectorXd CorrelationVector = Eigen::RowVectorXd::Zero(num_nodes_reduced_space);
            auto it_node = it_node_begin + i;
            const auto it_node_reduced_begin = reduced_space_nodes.begin();
            // Assemble correlation vector
            for( int j = 0; j < num_nodes_reduced_space; j++){
                auto it_node_reduced = it_node_reduced_begin + j;
                CorrelationVector(j) = CorrelationFunction( it_node, *it_node_reduced, correlationLength);
            }
            // Assemble perturbation field
            for( int j = 0; j < num_random_variables; j++){
                rPerturbationMatrix(i,j) = sqrt( 1/Eigenvalues(j) ) * (CorrelationVector).dot(Eigenvectors.col(j));
            }
        }
    }
    KRATOS_INFO("PerturbGeometrySubgridProcess: Assemble Random Field Time")
            << assemble_random_field_time.ElapsedSeconds() << std::endl;


    return num_random_variables;
    KRATOS_CATCH("");
}

} // namespace Kratos
