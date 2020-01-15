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

int PertubeGeometrySmallCorrelationLength::CreateEigenvectors( ModelPart& rModelPart, double correlationLength, double truncationTolerance ){
    KRATOS_TRY;
    
    int NumOfNodes = rModelPart.NumberOfNodes();
    KRATOS_WATCH( NumOfNodes );

    int NumOfRandomVariables = 6;

    BuiltinTimer reduceModel;

    searcher = new OMP_NodeSearch;
    ModelPart::NodesContainerType nodes = rModelPart.Nodes();

    ResultNodesContainerType  results;
    std::vector<std::vector<double>> resutlsDistance;
    
    searcher->InitializeSearch(nodes);
    
    SparseMatrixType CorrelationMatrix;
    CorrelationMatrix.resize(NumOfNodes,NumOfNodes);  
    TSparseSpaceType::SetToZero(CorrelationMatrix);

    int row_counter = -1;
    int counter = 0;
    for ( ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); itNode++ )
    {
        row_counter++;
        results = { rModelPart.Nodes().GetContainer()[ itNode->GetId() - 1 ] };
        searcher->SearchNodesInRadiusExclusiveImplementation(nodes,itNode->GetId()-1,3*correlationLength,results);

        for( int i = 0; i < results.size(); i++ ){
            counter++;
            CorrelationMatrix(row_counter , results[i]->GetId() - 1 ) =  CorrelationFunction( itNode, results[i] , correlationLength);
        }
    }
    KRATOS_WATCH( counter );
    for(int i = 0; i < 10; i++)
    {
        for( int j = 0; j < 10; j++)
        {
            std::cout << CorrelationMatrix(i,j) << ", ";
            
        }
        std::cout << std::endl;
    }
    for( int i = 0; i < 100; i++)
    {
        for( int j = 0; j < 100; j++)
        {
            CorrelationMatrix_check_orig(i,j) = CorrelationMatrix(i,j);
        }
    }
    
    // Solve Eigenvalue Problem
    Parameters params(R"(
        {
            "solver_type": "eigen_eigensystem",
            "number_of_eigenvalues": 90,
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
    // KRATOS_WATCH( Eigenvalues );
    // KRATOS_WATCH( Eigenvectors.size1() );
    // KRATOS_WATCH( Eigenvectors.size2() );
    // KRATOS_WATCH( row(Eigenvectors,0).size() );
    // KRATOS_WATCH( CorrelationMatrix.size1() );
    // KRATOS_WATCH( CorrelationMatrix.size2() );
    //KRATOS_WATCH( Eigenvectors );
    
    mDisplacement.resize(NumOfNodes,NumOfRandomVariables);
    int j = -1;
    DenseVectorType CorrelationVector;
    for (ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); itNode++)
    {
        j++;
        for( int i = 0; i < NumOfRandomVariables; i++)
        {
            mDisplacement(j,i) = sqrt( Eigenvalues(i) ) * inner_prod( row(Eigenvectors,i),  row( CorrelationMatrix, j));
        }
    }
    //rModelPart.AddNodalSolutionStepVariable(NORMAL);
    //MortarUtilities::ComputeNodesMeanNormalModelPart( rModelPart, false );
    return NumOfRandomVariables;

    KRATOS_CATCH("")

}

void PertubeGeometrySmallCorrelationLength::AssembleEigenvectors( ModelPart& rModelPart, const std::vector<double>& variables )
{
    int NumOfRandomVariables = variables.size();
    //KRATOS_WATCH( variables );
    int j = -1;
    double max = 0.0;
    array_1d<double, 3> normal;
    ModelPart::NodeIterator itNode_initial = mrInitialModelPart.NodesBegin();
    for (ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); itNode++)
    {
        j++;
        double tmp = 0.0;
        normal =  itNode_initial->FastGetSolutionStepValue(NORMAL);
        itNode_initial = itNode_initial + 1;
        for( int i = 0; i < NumOfRandomVariables; i++)
        {
            itNode->GetInitialPosition().Coordinates() += normal*mMaximalDisplacement*variables[i]*mDisplacement(j,i);
            tmp += mMaximalDisplacement*variables[i]*mDisplacement(j,i);
        } 
        if( std::abs(tmp) > std::abs(max) )
        {
            max = tmp;
        }     
    }
    //std::cout << "Maximal Displacement: " << max << std::endl;
    Eigen::MatrixXd CorrelationMatrix_tmp;
    CorrelationMatrix_tmp = Eigen::MatrixXd::Zero(rModelPart.NumberOfNodes(), rModelPart.NumberOfNodes());
    // Remove this later again!
    // ################################################################################
    int row_counter = -1;
    for ( ModelPart::NodeIterator it = rModelPart.NodesBegin(); it != rModelPart.NodesEnd(); it++ )
    {
        row_counter++;
        int column_counter = -1;
        for( ModelPart::NodeIterator it_inner = rModelPart.NodesBegin(); it_inner != rModelPart.NodesEnd(); it_inner++ )
        {
            column_counter++;
            array_1d<double, 3> coorrdinate1;
            array_1d<double, 3> coorrdinate2;
            //coorrdinate1 = it->GetInitialPosition().Coordinates();
            //coorrdinate2 = it_inner->GetInitialPosition().Coordinates();
            double c1 = it->GetInitialPosition().Coordinates()(2);
            double c2 = it_inner->GetInitialPosition().Coordinates()(2);
            double norm1 = sqrt( coorrdinate1(0)*coorrdinate1(0) + coorrdinate1(1)*coorrdinate1(1) + coorrdinate1(2)*coorrdinate1(2) );
            double norm2 = sqrt( coorrdinate2(0)*coorrdinate2(0) + coorrdinate2(1)*coorrdinate2(1) + coorrdinate2(2)*coorrdinate2(2) );
            CorrelationMatrix_tmp(row_counter ,column_counter ) =  c1*c2;
        }
    }
    CorrelationMatrix_check += CorrelationMatrix_tmp;
    //################################################################################
}

void PertubeGeometrySmallCorrelationLength::Average(int number)
{
    for(int i = 0; i < CorrelationMatrix_check.rows(); i++)
    {
        for( int j = 0; j < CorrelationMatrix_check.cols(); j++)
        {
            CorrelationMatrix_check(i,j) = CorrelationMatrix_check(i,j)/ (double)number;
        }
    }
    for(int i = 0; i < 10; i++)
    {
        for( int j = 0; j < 10; j++)
        {
            std::cout << CorrelationMatrix_check(i,j) << ", ";
        }
        std::cout << std::endl;
    }
    // KRATOS_WATCH( CorrelationMatrix_check_orig - CorrelationMatrix_check );
    // KRATOS_WATCH( CorrelationMatrix_check_orig(0,4) );
    // KRATOS_WATCH( CorrelationMatrix_check(0,4) );
    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(CorrelationMatrix_check.rows());
    // es.compute(CorrelationMatrix_check,Eigen::ComputeEigenvectors);
    // for( int i = 0; i < es.eigenvectors().col(0).size(); i++)
    // {
    //     std::cout << i << ": " << es.eigenvectors().col(0)(i) << "\t " << CorrelationMatrix_check_orig.col(0)(i) << std::endl;
    // }
    KRATOS_WATCH( ( CorrelationMatrix_check  ).rows() );
    //KRATOS_WATCH( ( CorrelationMatrix_check.cwiseAbs() - CorrelationMatrix_check_orig.cwiseAbs() ) );
    KRATOS_WATCH( ( CorrelationMatrix_check_orig - CorrelationMatrix_check ).norm() / CorrelationMatrix_check_orig.norm() );
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
