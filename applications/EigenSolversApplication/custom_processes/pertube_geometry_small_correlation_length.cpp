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

#include "utilities/mortar_utilities.h"
// Project includes
#include "custom_processes/pertube_geometry_small_correlation_length.h"
#include "utilities/builtin_timer.h"
#include "custom_utilities/ublas_wrapper.h"



namespace Kratos
{ 

    typedef ModelPart::NodesContainerType::ContainerType                 ResultNodesContainerType;
                                 
    typedef Kratos::intrusive_ptr<Kratos::Node<3UL> >                     NodesType;
    typedef TUblasSparseSpace<double> TSparseSpaceType;
    typedef TUblasDenseSpace<double> TDenseSpaceType;

    typedef typename TSparseSpaceType::MatrixType                       SparseMatrixType;

    typedef typename TDenseSpaceType::VectorType DenseVectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

int PertubeGeometrySmallCorrelationLength::CreateEigenvectors( double correlationLength, double truncationTolerance ){
    KRATOS_TRY;
    
    int NumOfNodes = mrThisModelPart.NumberOfNodes();
    KRATOS_WATCH( NumOfNodes );

    int NumOfRandomVariables = 6;

    BuiltinTimer reduceModel;

    searcher = new OMP_NodeSearch;
    ModelPart::NodesContainerType nodes = mrThisModelPart.Nodes();

    ResultNodesContainerType  results;
    std::vector<std::vector<double>> resutlsDistance;
    
    searcher->InitializeSearch(nodes);
    
    SparseMatrixType CorrelationMatrix;
    CorrelationMatrix.resize(NumOfNodes,NumOfNodes);  
    TSparseSpaceType::SetToZero(CorrelationMatrix);

    int row_counter = -1;
    int counter = 0;
    for ( ModelPart::NodeIterator itNode = mrThisModelPart.NodesBegin(); itNode != mrThisModelPart.NodesEnd(); itNode++ )
    {
        row_counter++;
        results = { mrThisModelPart.Nodes().GetContainer()[ itNode->GetId() - 1 ] };
        searcher->SearchNodesInRadiusExclusiveImplementation(nodes,itNode->GetId()-1,3*correlationLength,results);

        for( int i = 0; i < results.size(); i++ ){
            counter++;
            CorrelationMatrix(row_counter , results[i]->GetId() - 1 ) =  CorrelationFunction( itNode, results[i] , correlationLength);
        }
    }
    KRATOS_WATCH( counter );
    
    // Solve Eigenvalue Problem
    Parameters params(R"(
        {
            "solver_type": "eigen_eigensystem",
            "number_of_eigenvalues": 150,
            "normalize_eigenvectors": false,
            "max_iteration": 1000,
            "tolerance": 1e-3,
            "echo_level": 1
        })");
    DenseMatrixType Eigenvectors;
    DenseVectorType Eigenvalues;

    SparseMatrixType IdentityMatrix_ =  IdentityMatrix(NumOfNodes,NumOfNodes);

    mpEigenSolver = new EigensystemSolver<>(params);
    //CorrelationMatrix = - CorrelationMatrix;
    mpEigenSolver->Solve( IdentityMatrix_,
                          CorrelationMatrix,
                          Eigenvalues,
                          Eigenvectors  );

    double sum_eig =  1 / Eigenvalues(0);
    int numb = 0;
    double espilon = 0.0;
    for( int i = 1; i < Eigenvalues.size(); i++)
    {
        espilon = 1 - sum_eig / ( sum_eig + 1 / Eigenvalues(i) );
        if( espilon < truncationTolerance)
        {
            numb = i + 1;
            std::cout << "Truncation converged with a tolerance of:  " <<  truncationTolerance << std::endl;
            std::cout << "Needed Eigenvectors: " << numb << std::endl;
            break;
        }
        else if( i == Eigenvalues.size() - 1)
        {
            numb = Eigenvalues.size();
            std::cout << "Truncation NOT voncerged, achieved tolerance:  " <<  espilon << " / " << truncationTolerance << std::endl;
            std::cout << "Maximum number of computed eigenvalues is used: " << numb << std::endl;
        }
        sum_eig += 1/Eigenvalues(i);
    }

    NumOfRandomVariables = numb;
    double eucledian_norm;
    
    for( int i = 0; i < NumOfRandomVariables; i++ )
    {
        double eucledian_norm =  norm_2( row(Eigenvectors,i) );
        noalias( row(Eigenvectors,i) )= 1.0 / eucledian_norm * row(Eigenvectors,i);
    }
    KRATOS_WATCH( Eigenvalues );
    KRATOS_WATCH( Eigenvectors.size1() );
    KRATOS_WATCH( Eigenvectors.size2() );
    KRATOS_WATCH( row(Eigenvectors,0).size() );
    KRATOS_WATCH( CorrelationMatrix.size1() );
    KRATOS_WATCH( CorrelationMatrix.size2() );
    //KRATOS_WATCH( Eigenvectors );
    
    mDisplacement.resize(NumOfNodes,NumOfRandomVariables);
    int j = -1;
    DenseVectorType CorrelationVector;
    for (ModelPart::NodeIterator itNode = mrThisModelPart.NodesBegin(); itNode != mrThisModelPart.NodesEnd(); itNode++)
    {
        j++;
        for( int i = 0; i < NumOfRandomVariables; i++)
        {
            mDisplacement(j,i) = sqrt( Eigenvalues(i) ) * inner_prod( row(Eigenvectors,i),  row( CorrelationMatrix, j));
        }
    }
    //mrThisModelPart.AddNodalSolutionStepVariable(NORMAL);
    MortarUtilities::ComputeNodesMeanNormalModelPart( mrThisModelPart, false );
    return NumOfRandomVariables;

    KRATOS_CATCH("")

}

void PertubeGeometrySmallCorrelationLength::AssembleEigenvectors( const std::vector<double>& variables, double maxDisplacement )
{
    int NumOfRandomVariables = variables.size();
    KRATOS_WATCH( variables );
    int j = -1;
    double max = 0.0;
    array_1d<double, 3> normal;
    for (ModelPart::NodeIterator itNode = mrThisModelPart.NodesBegin(); itNode != mrThisModelPart.NodesEnd(); itNode++)
    {
        j++;
        double tmp = 0.0;
        normal =  itNode->FastGetSolutionStepValue(NORMAL);
        for( int i = 0; i < NumOfRandomVariables; i++)
        {
            itNode->GetInitialPosition().Coordinates() += normal*maxDisplacement*variables[i]*mDisplacement(j,i);
            tmp += maxDisplacement*variables[i]*mDisplacement(j,i);
        } 
        if( std::abs(tmp) > std::abs(max) )
        {
            max = tmp;
        }     
    }
    std::cout << "Maximal Displacement: " << max << std::endl;
}

double PertubeGeometrySmallCorrelationLength::CorrelationFunction( ModelPart::NodeIterator itNode1,NodesType itNode2, double CorrelationLenth)
{
    array_1d<double, 3> coorrdinate;
    coorrdinate = itNode1->GetInitialPosition().Coordinates() - itNode2->GetInitialPosition().Coordinates();
    
    double norm = sqrt( coorrdinate(0)*coorrdinate(0) + coorrdinate(1)*coorrdinate(1) + coorrdinate(2)*coorrdinate(2) );

    return( exp( - norm*norm / (CorrelationLenth*CorrelationLenth) ) );  
}


double PertubeGeometrySmallCorrelationLength::Kernel( double x, double sigma )
{
    return exp( - (x*x) / (0.45*0.45*sigma*sigma) );
}

} // namespace Kratos
