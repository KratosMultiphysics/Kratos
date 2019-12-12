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
#include "custom_processes/pertube_geometry.h"
#include "utilities/builtin_timer.h"



namespace Kratos
{ 

typedef ModelPart::NodesContainerType::ContainerType                 ResultNodesContainerType;

int PertubeGeometryProcess::CreateEigenvectors( double minDistance, double correlationLength, double truncationTolerance ){
    KRATOS_TRY;
    
    int NumOfNodes = mrThisModelPart.NumberOfNodes();
    KRATOS_WATCH( NumOfNodes );
    // double minDistance = 0.05;
    // double correlationLength = 0.5;

    BuiltinTimer reduceModel;

    searcher = new OMP_NodeSearch;
    ModelPart::NodesContainerType nodes = mrThisModelPart.Nodes();
    ModelPart::NodesContainerType nodes2 = {};
    double radius = minDistance;
    ResultNodesContainerType  results;
    std::vector<std::vector<double>> resutlsDistance;
    
    searcher->InitializeSearch(nodes);
    std::vector<ModelPart::NodeIterator> reduced_space_nodes;
    int counter = 1;
    ModelPart::NodeIterator it_n = mrThisModelPart.NodesBegin();
    
    for (ModelPart::NodeIterator itNode = mrThisModelPart.NodesBegin(); itNode != mrThisModelPart.NodesEnd(); itNode++)
    {   
        itNode->Set(VISITED,false);
    }
    for (ModelPart::NodeIterator itNode = mrThisModelPart.NodesBegin(); itNode != mrThisModelPart.NodesEnd(); itNode++)
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
    // KRATOS_WATCH( es.eigenvalues() );
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
            std::cout << "Truncation NOT voncerged, achieved tolerance:  " <<  espilon << " / " << 0.01 << std::endl;
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

    for (ModelPart::NodeIterator itNode = mrThisModelPart.NodesBegin(); itNode != mrThisModelPart.NodesEnd(); itNode++)
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
            // if( Displacement(j,i) != 0.0 )
            // {
            //     Displacement(j,i) = 1.0;
            // }
        }
        
    }
    MortarUtilities::ComputeNodesMeanNormalModelPart( mrThisModelPart, false );
    return NumOfRandomVariables;
    KRATOS_CATCH("");
}

void PertubeGeometryProcess::AssembleEigenvectors( const std::vector<double>& variables, double maxDisplacement )
{
    int NumOfRandomVariables = variables.size();
    KRATOS_WATCH( variables );
    int j = -1;
    double max = 0.0;
    for (ModelPart::NodeIterator itNode = mrThisModelPart.NodesBegin(); itNode != mrThisModelPart.NodesEnd(); itNode++)
    {
        j++;
        double tmp = 0.0;
        array_1d<double, 3> normal;
        normal =  itNode->FastGetSolutionStepValue(NORMAL);
        for( int i = 0; i < NumOfRandomVariables; i++)
        {
            itNode->GetInitialPosition().Coordinates() += normal*maxDisplacement*variables[i]*Displacement(j,i);
            //itNode->GetInitialPosition().Coordinates()(1) += maxDisplacement*variables[i]*Displacement(j,i);
            tmp += maxDisplacement*variables[i]*Displacement(j,i);
        } 
        if( std::abs(tmp) > std::abs(max) )
        {
            max = tmp;
        }     
    }
    std::cout << "Maximal Displacement: " << max << std::endl;
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
