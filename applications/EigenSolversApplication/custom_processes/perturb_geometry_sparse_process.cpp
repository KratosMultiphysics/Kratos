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
#include "custom_processes/perturb_geometry_sparse_process.h"
#include "custom_utilities/omp_node_search.h"
#include "utilities/builtin_timer.h"
#include "custom_utilities/ublas_wrapper.h"



namespace Kratos
{
    typedef TUblasSparseSpace<double> TSparseSpaceType;
    typedef TUblasDenseSpace<double> TDenseSpaceType;

    typedef ModelPart::NodesContainerType::ContainerType                ResultNodesContainerType;

    typedef typename TSparseSpaceType::MatrixType                       SparseMatrixType;

    typedef typename TDenseSpaceType::VectorType                        DenseVectorType;

    typedef typename TDenseSpaceType::MatrixType                        DenseMatrixType;

int PerturbGeometrySparseProcess::CreateEigenvectors( ModelPart& rModelPart, double correlationLength, double truncationTolerance ){
    KRATOS_TRY;

    int number_of_nodes = rModelPart.NumberOfNodes();

    BuiltinTimer assemble_cl_matrix_timer;

    // Construct and initialize required containers
    const ModelPart::NodesContainerType nodes = rModelPart.Nodes();
    ResultNodesContainerType  results;
    SparseMatrixType CorrelationMatrix;
    CorrelationMatrix.resize(number_of_nodes,number_of_nodes);
    TSparseSpaceType::SetToZero(CorrelationMatrix);
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
        for( int j = 0; j < results.size(); j++){
            counter++;
            int index = results[j]->GetId() - 1;
            auto it_node_results = it_node_begin + index;
            CorrelationMatrix(i , index ) =  CorrelationFunction( it_node, it_node_results , correlationLength);
        }
    }

    KRATOS_INFO("PerturbGeometrySparseProcess: Assemble Correlation Matrix Time")
            << assemble_cl_matrix_timer.ElapsedSeconds() << std::endl
            << (double) counter / (double)(number_of_nodes*number_of_nodes) * 100 << "% of Matrix is populated" << std::endl;

    // Solve Eigenvalue Problem
    Parameters params(R"(
        {
            "solver_type": "eigen_eigensystem",
            "number_of_eigenvalues": 100,
            "normalize_eigenvectors": false,
            "max_iteration": 1000,
            "tolerance": 1e-4,
            "echo_level": 1
        })");
    DenseMatrixType Eigenvectors;
    DenseVectorType Eigenvalues;

    SparseMatrixType IdentityMatrix_ =  IdentityMatrix(number_of_nodes,number_of_nodes);

    EigensystemSolver<> eigensolver(params);
    eigensolver.Solve(  IdentityMatrix_,
                        CorrelationMatrix,
                        Eigenvalues,
                        Eigenvectors);

    // Find number of required eigenvalues to statisfy convergence criterion
    double sum_eigenvalues =  1 / Eigenvalues(0);
    int number_required_eigevalues = 0;
    double epsilon;
    for( int i = 1; i < Eigenvalues.size(); i++)
    {
        epsilon = 1 - sum_eigenvalues / ( sum_eigenvalues + 1 / Eigenvalues(i) );
        if( epsilon < truncationTolerance)
        {
            number_required_eigevalues = i + 1;
            KRATOS_INFO("PerturbGeometrySparseProcess")
                << "Truncation Error  (" <<  truncationTolerance << std::endl
                << ") is achieved with : " << number_required_eigevalues << " Eigenvalues" << std::endl;
            break;
        }
        else if( i == Eigenvalues.size() - 1)
        {
            number_required_eigevalues = Eigenvalues.size();
            KRATOS_INFO("PerturbGeometrySparseProcess")
                << "Truncation NOT converged, achieved tolerance:  " <<  epsilon << " / " << truncationTolerance << std::endl
                << "Maximum number of computed eigenvalues is used: " << number_required_eigevalues << std::endl;
        }
        sum_eigenvalues += 1/Eigenvalues(i);
    }

    int NumOfRandomVariables = number_required_eigevalues;
    double eucledian_norm;

    // Construct and initialize perturbation matrix
    DenseMatrixType& rPerturbationMatrix = *mpPerturbationMatrix;
    rPerturbationMatrix.resize(number_of_nodes,NumOfRandomVariables);
    DenseVectorType CorrelationVector;

    // Assemble perturbation matrix
    for( int i = 0; i < number_of_nodes; i++){
        for( int j = 0; j < NumOfRandomVariables; j++)
        {
            rPerturbationMatrix(i,j) = sqrt( Eigenvalues(j) ) * inner_prod( row(Eigenvectors,j),  row( CorrelationMatrix, i));
        }
    }
    // int j = -1;
    // for (ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); itNode++)
    // {
    //     j++;
    //     for( int i = 0; i < NumOfRandomVariables; i++)
    //     {
    //         rPerturbationMatrix(j,i) = sqrt( Eigenvalues(i) ) * inner_prod( row(Eigenvectors,i),  row( CorrelationMatrix, j));
    //     }
    // }

    return NumOfRandomVariables;

    KRATOS_CATCH("")

}

} // namespace Kratos
