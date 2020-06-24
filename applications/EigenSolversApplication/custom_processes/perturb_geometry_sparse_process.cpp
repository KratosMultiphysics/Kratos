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
#include "custom_processes/perturb_geometry_sparse_process.h"
#include "custom_utilities/omp_node_search.h"
#include "utilities/builtin_timer.h"


namespace Kratos
{
    typedef TUblasSparseSpace<double> TSparseSpaceType;

    typedef ModelPart::NodesContainerType::ContainerType                ResultNodesContainerType;

    typedef typename TSparseSpaceType::MatrixType                       SparseMatrixType;


/**
 * @brief Creates Eigenvectors of a sparse correlation matrix
 * @details Generates sparse correlation matrix. Decomposes correlation matrix.
 * @param correlation_matrix Correlation matrix. Stores correlation value for all nodes.
 * @param rPerturbationMatrix Perturbation matrix. Stores eigenvectors of correlation matrix.
 */
int PerturbGeometrySparseProcess::CreateEigenvectors( ModelPart& rModelPart, double correlationLength, double truncationTolerance, Parameters eigenvalueSettings ){
    KRATOS_TRY;

    int number_of_nodes = rModelPart.NumberOfNodes();

    BuiltinTimer assemble_cl_matrix_timer;

    // Construct and initialize required containers
    const ModelPart::NodesContainerType nodes = rModelPart.Nodes();
    ResultNodesContainerType  results;
    SparseMatrixType correlation_matrix;
    correlation_matrix.resize(number_of_nodes,number_of_nodes);
    TSparseSpaceType::SetToZero(correlation_matrix);
    // Construct and initialize searcher
    OMP_NodeSearch searcher;
    searcher.InitializeSearch(nodes);

    // Assemble correlation matrix
    int counter = 0;
    const auto it_node_begin = rModelPart.NodesBegin();

    for( int i = 0; i < number_of_nodes; i++){
        auto it_node = it_node_begin + i;
        // Add current nodes to results
        results = { rModelPart.Nodes().GetContainer()[i] };
        // Find all neighbouring nodes
        searcher.SearchNodesInRadiusExclusiveImplementation(nodes,it_node->GetId()-1,3*correlationLength,results);
        for( size_t j = 0; j < results.size(); j++){
            counter++;
            int index = results[j]->GetId() - 1;
            auto it_node_results = it_node_begin + index;

            correlation_matrix(i , index ) =  CorrelationFunction( it_node, it_node_results , correlationLength);
        }
    }

    KRATOS_INFO_IF("PerturbGeometrySparseProcess: Build Correlation Matrix Time", mEchoLevel > 0)
            << assemble_cl_matrix_timer.ElapsedSeconds() << std::endl
            << (double) counter / (double)(number_of_nodes*number_of_nodes) * 100 << "% of Matrix is populated" << std::endl;

    DenseMatrixType eigenvectors;
    DenseVectorType eigenvalues;

    SparseMatrixType identity_matrix =  IdentityMatrix(number_of_nodes,number_of_nodes);

    EigensystemSolver<> eigensolver(eigenvalueSettings);
    eigensolver.Solve(  identity_matrix,
                        correlation_matrix,
                        eigenvalues,
                        eigenvectors);

    // Find number of required eigenvalues to statisfy convergence criterion
    double sum_eigenvalues =  1 / eigenvalues(0);
    int number_required_eigenvalues = 0;
    double epsilon;
    for( size_t i = 1; i < eigenvalues.size(); i++)
    {
        epsilon = 1 - sum_eigenvalues / ( sum_eigenvalues + 1 / eigenvalues(i) );
        if( epsilon < truncationTolerance)
        {
            number_required_eigenvalues = i + 1;
            KRATOS_INFO_IF("PerturbGeometrySparseProcess", mEchoLevel > 0)
                << "Truncation Error  (" << truncationTolerance
                << ") is achieved with " << number_required_eigenvalues << " Eigenvalues" << std::endl;

            break;
        }
        else if( i == eigenvalues.size() - 1)
        {
            number_required_eigenvalues = eigenvalues.size();
            KRATOS_INFO_IF("PerturbGeometrySparseProcess", mEchoLevel > 0)
                << "Truncation Error is NOT achieved:  " << epsilon << " / " << truncationTolerance << std::endl
                << "Maximum number of computed eigenvalues is used: " << number_required_eigenvalues << std::endl;
        }
        sum_eigenvalues += 1/eigenvalues(i);
    }

    // Normalize required eigenvectors
    double eucledian_norm = 0;
    for( int i =  0; i < number_required_eigenvalues; i++)
    {
        eucledian_norm =  norm_2( row(eigenvectors,i) );
        row(eigenvectors,i) = 1.0 / eucledian_norm * row(eigenvectors,i);
    }

    int number_of_random_variables = number_required_eigenvalues;

    // Construct and initialize perturbation matrix
    DenseMatrixType& rPerturbationMatrix = *mpPerturbationMatrix;
    rPerturbationMatrix.resize(number_of_nodes,number_of_random_variables);

    // Assemble perturbation matrix
    BuiltinTimer assemble_random_field_time;
    #pragma omp parallel for
    for( int i = 0; i < number_of_nodes; i++){
        for( int j = 0; j < number_of_random_variables; j++)
        {
            rPerturbationMatrix(i,j) = sqrt( eigenvalues(j) ) * inner_prod( row(eigenvectors,j),  row( correlation_matrix, i));
        }
    }

    KRATOS_INFO_IF("PerturbGeometrySparseProcess: Assemble Random Field Time", mEchoLevel > 0)
        << assemble_random_field_time.ElapsedSeconds() << std::endl;

    return number_of_random_variables;

    KRATOS_CATCH("")

}

} // namespace Kratos
