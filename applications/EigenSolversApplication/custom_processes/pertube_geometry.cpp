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
#include "custom_processes/pertube_geometry.h"
#include "utilities/builtin_timer.h"



namespace Kratos
{ 

typedef ModelPart::NodesContainerType::ContainerType                 ResultNodesContainerType;

int PertubeGeometryProcess::CreateEigenvectors( ModelPart& rThisModelPart, double minDistance, double correlationLength, double truncationTolerance ){
    KRATOS_TRY;
    
    int NumOfNodes = rThisModelPart.NumberOfNodes();
    KRATOS_WATCH( NumOfNodes );
    // double minDistance = 0.05;
    // double correlationLength = 0.5;

    BuiltinTimer reduceModel;

    searcher = new OMP_NodeSearch;
    ModelPart::NodesContainerType nodes = rThisModelPart.Nodes();
    ModelPart::NodesContainerType nodes2 = {};
    double radius = minDistance;
    ResultNodesContainerType  results;
    std::vector<std::vector<double>> resutlsDistance;
    
    searcher->InitializeSearch(nodes);
    std::vector<ModelPart::NodeIterator> reduced_space_nodes;
    int counter = 1;
    ModelPart::NodeIterator it_n = rThisModelPart.NodesBegin();
    
    for (ModelPart::NodeIterator itNode = rThisModelPart.NodesBegin(); itNode != rThisModelPart.NodesEnd(); itNode++)
    {   
        itNode->Set(VISITED,false);
    }
    for (ModelPart::NodeIterator itNode =rThisModelPart.NodesBegin(); itNode != rThisModelPart.NodesEnd(); itNode++)
    {   
        if( !itNode->Is(VISITED) ) {
            itNode->Set(VISITED,true);
            reduced_space_nodes.push_back(itNode);
            results = {};
            searcher->SearchNodesInRadiusExclusiveImplementation(nodes,itNode->GetId()-1,radius,results);
            for( int i = 0; i < results.size(); i++ ){
                results[i]->Set(VISITED,true);
            }
        }
    }
    KRATOS_INFO("System Construction Time")
            << reduceModel.ElapsedSeconds() << std::endl;
    
    Eigen::MatrixXd CorrelationMatrix; 
    double redSpace = reduced_space_nodes.size();
    KRATOS_WATCH( redSpace );
    std::cout << "Size: " <<  reduced_space_nodes.size() << std::endl;
    CorrelationMatrix.resize(redSpace,redSpace);
    int row_counter = -1;
    for ( std::vector<ModelPart::NodeIterator>::iterator it = reduced_space_nodes.begin(); it != reduced_space_nodes.end(); it++ )
    {
        row_counter++;
        int column_counter = -1;
        for( std::vector<ModelPart::NodeIterator>::iterator it_inner = reduced_space_nodes.begin(); it_inner != reduced_space_nodes.end(); it_inner++ )
        {
            column_counter++;
            CorrelationMatrix(row_counter ,column_counter ) =  CorrelationFunction( (*it), (*it_inner), correlationLength);  
        }
    }
    // Remove this
    row_counter = -1;
    for ( ModelPart::NodeIterator it = rThisModelPart.NodesBegin(); it != rThisModelPart.NodesEnd(); it++ )
    {
        row_counter++;
        int column_counter = -1;
        for( ModelPart::NodeIterator it_inner = rThisModelPart.NodesBegin(); it_inner != rThisModelPart.NodesEnd(); it_inner++ )
        {
            column_counter++;
            CorrelationMatrix_check_orig(row_counter ,column_counter ) =  CorrelationFunction( (it), (it_inner), correlationLength);  
        }
    }
    for(int i = 0; i < 10; i++)
    {
        for( int j = 0; j < 10; j++)
        {
            std::cout << CorrelationMatrix_check_orig(i,j) << ", ";
        }
        std::cout << std::endl;
    }
    //#######################
    // Remove this later
    //CorrelationMatrix_check_orig = CorrelationMatrix;
    //KRATOS_WATCH( CorrelationMatrix );
    // Solve Eigenvalue Problem
    Parameters params(R"(
        {
            "solver_type": "eigen_eigensystem",
            "number_of_eigenvalues": 5,
            "normalize_eigenvectors": true,
            "max_iteration": 100,
            "tolerance": 1e-6,
            "echo_level": 1
        })");
    //KRATOS_WATCH(CorrelationMatrix);
    std::cout << "Eigensolver start: " << std::endl;
    BuiltinTimer Eigensolver_time;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(redSpace);
    es.compute(CorrelationMatrix,Eigen::ComputeEigenvectors);
    KRATOS_INFO("Eigensolver Time")
            << Eigensolver_time.ElapsedSeconds() << std::endl;
    std::cout << "Eigensolver done: " << std::endl;
    //es.compute(CorrelationMatrix,Eigen::ComputeEigenvectors);
    if( es.info() != Eigen::Success)
    {
        std::cout << "Decomposition failed!" << std::endl;  
    }

    // Find number of neccessary eigenvectors
    // KRATOS_WATCH( es.eigenvectors().col(0) );


    double sum_eig =  es.eigenvalues()( es.eigenvalues().size()-1 );
    int numb = 0;
    double espilon = 0.0;
    for( int i = es.eigenvalues().size() - 2; i >= 0; i--)
    {
        espilon = 1 - sum_eig / ( sum_eig + es.eigenvalues()(i) );
        numb++;
        KRATOS_WATCH( espilon );
        if( espilon < truncationTolerance)
        {
            std::cout << "Truncation converged with a tolerance of:  " <<  truncationTolerance << std::endl;
            std::cout << "Needed Eigenvectors: " << numb << std::endl;
            break;
        }
        else if( i == es.eigenvalues().size() - 1)
        {
            numb = es.eigenvalues().size();
            std::cout << "Truncation NOT concerged, achieved tolerance:  " <<  espilon << " / " << 0.01 << std::endl;
            std::cout << "Maximum number of computed eigenvalues is used: " << numb << std::endl;
        }
        sum_eig += es.eigenvalues()(i);
    }
    
    // Different criteria
    // ###################################################################################
    double sum_a = 0.0;
    double sum_b = 0.0;
    KRATOS_WATCH( es.eigenvalues().size() );
    for( int i = es.eigenvalues().size() - 1; i > es.eigenvalues().size() - numb - 1; i--)
    {
        sum_a += es.eigenvalues()(i);
        //std::cout << "a: " << sum_a << std::endl;
    }
    for( int i = es.eigenvalues().size() - 1; i >= 0; i--)
    {
        sum_b += es.eigenvalues()(i);
    }
    std::cout << "espilone2: " << 1 - sum_a/sum_b << std::endl;
    // ######################################################################################
    
    int NumOfRandomVariables = numb;
    // Truncate Matrix
    Eigen::MatrixXd Eigenvectors2 = es.eigenvectors();

    int cutOff = redSpace - NumOfRandomVariables;
    Eigen::MatrixXd Eigenvectors = Eigen::Map<Eigen::MatrixXd,0,Eigen::OuterStride<> >(
    Eigenvectors2.real().data() + cutOff , es.eigenvectors().real().rows(), NumOfRandomVariables, Eigen::OuterStride<>(es.eigenvectors().outerStride()) ) ;
    //KRATOS_WATCH( Eigenvectors );
    // Truncate Vector
    //KRATOS_WATCH( es.eigenvalues() );
    Eigen::VectorXd Eigenvalues = es.eigenvalues().real().tail(NumOfRandomVariables);
    //KRATOS_WATCH( Eigenvalues );

    int j = -1;

    Displacement = Eigen::MatrixXd::Zero(NumOfNodes,NumOfRandomVariables); 
    Eigen::VectorXd CorrelationVector = Eigen::VectorXd::Zero(redSpace);

    for (ModelPart::NodeIterator itNode = rThisModelPart.NodesBegin(); itNode != rThisModelPart.NodesEnd(); itNode++)
    {
        j++;
        double tmp = 0;
        CorrelationVector = Eigen::RowVectorXd::Zero(redSpace);
        int k = -1;
        for( std::vector<ModelPart::NodeIterator>::iterator it = reduced_space_nodes.begin(); it != reduced_space_nodes.end(); it++){
            k++;
            // Remove if statement
            if( (*it) == itNode )
            {
                CorrelationVector(k) = CorrelationFunction( itNode, (*it), correlationLength);
            }
            else {
                CorrelationVector(k) = CorrelationFunction( itNode, (*it), correlationLength);
            }
        }

        for( int i = 0; i < NumOfRandomVariables; i++)
        {
            Displacement(j,i) = sqrt( 1/Eigenvalues(i) ) * (CorrelationVector).dot(Eigenvectors.col(i));

        }
        
    }
    
    return NumOfRandomVariables;
    KRATOS_CATCH("");
}

void PertubeGeometryProcess::AssembleEigenvectors( ModelPart& rThisModelPart, const std::vector<double>& variables, double correlationLength )
{
    //Todo: Check if NumOfRandomVariables==NumberOfEigenvectors
    int NumOfRandomVariables = variables.size();
    // KRATOS_WATCH( NumOfRandomVariables );
    //KRATOS_WATCH( variables );
    int j = -1;
    double max = 0.0;
    array_1d<double, 3> normal;
    ModelPart::NodeIterator itNode_intial = mrInitialModelPart.NodesBegin();
    for (ModelPart::NodeIterator itNode = rThisModelPart.NodesBegin(); itNode != rThisModelPart.NodesEnd(); itNode++)
    {
        j++;
        double tmp = 0.0;
        normal =  itNode_intial->FastGetSolutionStepValue(NORMAL);
        itNode_intial = itNode_intial + 1;
        for( int i = 0; i < NumOfRandomVariables; i++)
        {
            itNode->GetInitialPosition().Coordinates() += normal*mMaximalDisplacement*variables[i]*Displacement(j,i);
            //itNode->GetInitialPosition().Coordinates()(1) += maxDisplacement*variables[i]*Displacement(j,i);
            tmp += mMaximalDisplacement*variables[i]*Displacement(j,i);
        } 
        if( std::abs(tmp) > std::abs(max) )
        {
            max = tmp;
        }     
    }
    Eigen::MatrixXd CorrelationMatrix_tmp;
    CorrelationMatrix_tmp = Eigen::MatrixXd::Zero(rThisModelPart.NumberOfNodes(), rThisModelPart.NumberOfNodes());
    // Remove this later again!
    // ################################################################################
    int row_counter = -1;
    for ( ModelPart::NodeIterator it = rThisModelPart.NodesBegin(); it != rThisModelPart.NodesEnd(); it++ )
    {
        row_counter++;
        int column_counter = -1;
        for( ModelPart::NodeIterator it_inner = rThisModelPart.NodesBegin(); it_inner != rThisModelPart.NodesEnd(); it_inner++ )
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
    //std::cout << "Maximal Displacement: " << max << std::endl;
}

void PertubeGeometryProcess::Average(int number)
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
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(CorrelationMatrix_check.rows());
    es.compute(CorrelationMatrix_check,Eigen::ComputeEigenvectors);
    // for( int i = 0; i < es.eigenvectors().col(0).size(); i++)
    // {
    //     std::cout << i << ": " << es.eigenvectors().col(0)(i) << "\t " << CorrelationMatrix_check_orig.col(0)(i) << std::endl;
    // }

    KRATOS_WATCH( (CorrelationMatrix_check_orig - CorrelationMatrix_check ).norm() / CorrelationMatrix_check_orig.norm()  );
}


double PertubeGeometryProcess::CorrelationFunction( ModelPart::NodeIterator itNode1, ModelPart::NodeIterator itNode2, double CorrelationLenth)
{
    array_1d<double, 3> coorrdinate;
    coorrdinate = itNode1->GetInitialPosition().Coordinates() - itNode2->GetInitialPosition().Coordinates();

    double norm = sqrt( coorrdinate(0)*coorrdinate(0) + coorrdinate(1)*coorrdinate(1) + coorrdinate(2)*coorrdinate(2) );

    return( exp( - norm*norm / (CorrelationLenth*CorrelationLenth) ) );  
}


double PertubeGeometryProcess::Kernel( double x, double sigma )
{
    return exp( - (x*x) / (0.45*0.45*sigma*sigma) );
}

} // namespace Kratos
