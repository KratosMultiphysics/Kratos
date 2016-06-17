// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: structural_mechanics_application/license.txt
//
//  Main authors:    Mohamed Khalil
//                   Vicente Mataix Ferrándiz
//

// System includes

// External includes

// Project includes
#include "custom_conditions/mortar_contact_2D_condition.hpp"
#include "custom_utilities/structural_mechanics_math_utilities.hpp"
#include "utilities/math_utils.h"
#include "structural_mechanics_application.h"
#include "structural_mechanics_application_variables.h"
// #include "custom_utilities/tree_contact_search.h"

#include <algorithm>

namespace Kratos 
{
/**
 * Flags related to the condition computation
 */
KRATOS_CREATE_LOCAL_FLAG( MortarContact2DCondition, COMPUTE_RHS_VECTOR,                 0 );
KRATOS_CREATE_LOCAL_FLAG( MortarContact2DCondition, COMPUTE_LHS_MATRIX,                 1 );
KRATOS_CREATE_LOCAL_FLAG( MortarContact2DCondition, COMPUTE_RHS_VECTOR_WITH_COMPONENTS, 2 );
KRATOS_CREATE_LOCAL_FLAG( MortarContact2DCondition, COMPUTE_LHS_MATRIX_WITH_COMPONENTS, 3 );

/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

MortarContact2DCondition::MortarContact2DCondition(
    IndexType NewId,
    GeometryType::Pointer pGeometry ) :
    Condition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

MortarContact2DCondition::MortarContact2DCondition( 
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties ) :
    Condition( NewId, pGeometry, pProperties )
{
    
}

/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

MortarContact2DCondition::MortarContact2DCondition( 
    MortarContact2DCondition const& rOther ) :
    Condition( rOther )
{
    
}

/************************************* OPERATIONS **********************************/
/***********************************************************************************/

Condition::Pointer MortarContact2DCondition::Create( 
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties ) const
{
    return boost::make_shared<MortarContact2DCondition>( NewId, GetGeometry( ).Create( ThisNodes ), pProperties );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

MortarContact2DCondition::~MortarContact2DCondition( )
{
    
}

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::Initialize( )
{
    KRATOS_TRY;

    // TODO Add content
    
    //std::vector<contact_container> * ContactContainer = this->GetValue(CONTACT_CONTAINERS);
    
    KRATOS_WATCH("Hola");

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    // TODO Add content
    // Populate the vector of master elements
    std::vector<contact_container> * all_containers = this->GetValue(CONTACT_CONTAINERS);
    mThisMasterElements.resize( all_containers->size( ) );
    
    for ( unsigned int i_cond = 0; i_cond < all_containers->size(); ++i_cond )
    {
        mThisMasterElements[i_cond] = (*all_containers)[i_cond].condition;
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    // TODO Add content
    // TODO GetMasterElements
    
    // NOTE: Maybe is a good idea to add the tree_contact_search initialization here instead of the calling from python
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

double MortarContact2DCondition::LagrangeMultiplierShapeFunctionValue( 
    const IndexType& rPointNumber,
    const IndexType& rShapeFunctionIndex 
    )
{
    // NOTE: For more information look at the thesis of Popp page 93/236
    
    const unsigned number_nodes = GetGeometry( ).PointsNumber();
    double phi = 0.0;

    if (number_nodes == 2) // First order
    {
        const GeometryType::IntegrationPointType& pt = GetGeometry().IntegrationPoints( mThisIntegrationMethod )[rPointNumber];
        GeometryType::CoordinatesArrayType local_coordinates;
        double eta = GetGeometry( ).PointLocalCoordinates( local_coordinates, pt.Coordinates( ) )(0);
        if (rShapeFunctionIndex == 1 )
        {
            phi = ( 0.5 * ( 1.0 - 3.0 * eta ) );
        }
        else if (rShapeFunctionIndex == 2 )
        {
            phi = ( 0.5 * ( 1.0 + 3.0 * eta ) );
        }
        else
        {
            KRATOS_THROW_ERROR( std::logic_error, " The rShapeFunctionIndex is wrong: ", rShapeFunctionIndex );
        }
    }
    else if (number_nodes == 3) // Second order
    {
        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);
        if (rShapeFunctionIndex == 1 )
        {
            phi = Ncontainer(rPointNumber, 0) -  3.0/4.0 * Ncontainer(rPointNumber, 2) + 0.5;
        }
        else if (rShapeFunctionIndex == 2 )
        {
            phi = Ncontainer(rPointNumber, 1) -  3.0/4.0 * Ncontainer(rPointNumber, 2) + 0.5;
        }
        else if (rShapeFunctionIndex == 3 )
        {
            phi = 5.0/2.0 * Ncontainer(rPointNumber, 2) - 1.0;
        }
        else
        {
            KRATOS_THROW_ERROR( std::logic_error, " The rShapeFunctionIndex is wrong: ", rShapeFunctionIndex );
        }
    }
        
    return phi;
}

/***********************************************************************************/
/***********************************************************************************/

const Matrix MortarContact2DCondition::LagrangeMultiplierShapeFunctionLocalGradient( const IndexType& rPointNumber )
{
    // For 2D2N elements as presented in Popp's and Gitterle's work
    const unsigned int num_slave_nodes = GetGeometry( ).PointsNumber( );
    const unsigned int local_dimension_slave = GetGeometry( ).LocalSpaceDimension( );
    Matrix DPhi_De = Matrix( num_slave_nodes, local_dimension_slave );

    DPhi_De( 0, 0 ) = -3.0 / 2.0;
    DPhi_De( 0, 1 ) = +3.0 / 2.0;

    return DPhi_De;
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateLocalSystem( 
    std::vector<MatrixType>& rLeftHandSideMatrices,
    const std::vector<Variable<MatrixType> >& rLHSVariables,
    std::vector<VectorType>& rRightHandSideVectors,
    const std::vector<Variable<VectorType> >& rRHSVariables,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set(MortarContact2DCondition::COMPUTE_LHS_MATRIX_WITH_COMPONENTS);
    LocalSystem.CalculationFlags.Set(MortarContact2DCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);

    //Initialize sizes for the system components:
    if ( rLHSVariables.size( ) != rLeftHandSideMatrices.size( ) )
    {
        rLeftHandSideMatrices.resize( rLHSVariables.size( ) );
    }

    if ( rRHSVariables.size( ) != rRightHandSideVectors.size( ) )
    {
        rRightHandSideVectors.resize( rRHSVariables.size( ) );
    }

    LocalSystem.CalculationFlags.Set(MortarContact2DCondition::COMPUTE_LHS_MATRIX);
    for ( unsigned int i = 0; i < rLeftHandSideMatrices.size( ); i++ )
    {
        // Note: rRightHandSideVectors.size() > 0
        this->InitializeSystemMatrices( rLeftHandSideMatrices[i],
                                        rRightHandSideVectors[0],
                                        LocalSystem.CalculationFlags 
                                        );
    }

    LocalSystem.CalculationFlags.Set( MortarContact2DCondition::COMPUTE_RHS_VECTOR );
    LocalSystem.CalculationFlags.Set( MortarContact2DCondition::COMPUTE_LHS_MATRIX, false ); // Temporarily only
    for ( unsigned int i = 0; i < rRightHandSideVectors.size( ); i++ )
    {
        // Note: rLeftHandSideMatrices.size() > 0
        this->InitializeSystemMatrices( rLeftHandSideMatrices[0],
                                        rRightHandSideVectors[i],
                                        LocalSystem.CalculationFlags 
                                        );
    }
    LocalSystem.CalculationFlags.Set( MortarContact2DCondition::COMPUTE_LHS_MATRIX, true ); // Reactivated again


    // Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrices( rLeftHandSideMatrices );
    LocalSystem.SetRightHandSideVectors( rRightHandSideVectors );

    LocalSystem.SetLeftHandSideVariables( rLHSVariables );
    LocalSystem.SetRightHandSideVariables( rRHSVariables );

    // Calculate condition system
    this->CalculateConditionSystem( LocalSystem, rCurrentProcessInfo );

}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;

    // TODO: ADD CONTENT!!!!!
    
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set( MortarContact2DCondition::COMPUTE_LHS_MATRIX );
    LocalSystem.CalculationFlags.Set( MortarContact2DCondition::COMPUTE_RHS_VECTOR );

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix,
                                    rRightHandSideVector,
                                    LocalSystem.CalculationFlags 
                                  );
    
    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix( rLeftHandSideMatrix );
    LocalSystem.SetRightHandSideVector( rRightHandSideVector );

    //Calculate condition system
    this->CalculateConditionSystem( LocalSystem, rCurrentProcessInfo );

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateLeftHandSide( 
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContact2DCondition::COMPUTE_LHS_MATRIX );

    VectorType RightHandSideVector = Vector( );

    // Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix,
                                    RightHandSideVector,
                                    LocalSystem.CalculationFlags 
                                    );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix( rLeftHandSideMatrix );
    LocalSystem.SetRightHandSideVector( RightHandSideVector );

    //Calculate condition system
    this->CalculateConditionSystem( LocalSystem, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateLeftHandSide( 
    std::vector< MatrixType >& rLeftHandSideMatrices,
    const std::vector< Variable< MatrixType > >& rLHSVariables,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    LocalSystemComponents LocalSystem;

    LocalSystem.CalculationFlags.Set( MortarContact2DCondition::COMPUTE_LHS_MATRIX );
    LocalSystem.CalculationFlags.Set( MortarContact2DCondition::COMPUTE_LHS_MATRIX_WITH_COMPONENTS );

    VectorType RightHandSideVector = Vector( );

    for( unsigned int i = 0; i<rLeftHandSideMatrices.size( ); i++ )
    {
        this->InitializeSystemMatrices( rLeftHandSideMatrices[i],
                                        RightHandSideVector,
                                        LocalSystem.CalculationFlags 
                                        );
    }

    LocalSystem.SetLeftHandSideMatrices( rLeftHandSideMatrices );
    LocalSystem.SetRightHandSideVector( RightHandSideVector );

    this->CalculateConditionSystem( LocalSystem, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateRightHandSide( 
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContact2DCondition::COMPUTE_RHS_VECTOR);

    MatrixType LeftHandSideMatrix = Matrix( );

    // Initialize size for the system components
    this->InitializeSystemMatrices( LeftHandSideMatrix,
                                    rRightHandSideVector,
                                    LocalSystem.CalculationFlags 
                                  );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix( LeftHandSideMatrix );
    LocalSystem.SetRightHandSideVector( rRightHandSideVector );

    //Calculate condition system
    this->CalculateConditionSystem( LocalSystem, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateRightHandSide( 
    std::vector< VectorType >& rRightHandSideVectors,
    const std::vector< Variable< VectorType > >& rRHSVariables,
    ProcessInfo& rCurrentProcessInfo )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContact2DCondition::COMPUTE_RHS_VECTOR );
    LocalSystem.CalculationFlags.Set( MortarContact2DCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS );

    MatrixType LeftHandSideMatrix = Matrix( );

    // Initialize size for the system components
    for( unsigned int i=0; i<rRightHandSideVectors.size(); i++ )
    {
        this->InitializeSystemMatrices( LeftHandSideMatrix,
                                        rRightHandSideVectors[i],
                                        LocalSystem.CalculationFlags 
                                      );
    }

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix( LeftHandSideMatrix );
    LocalSystem.SetRightHandSideVectors( rRightHandSideVectors );

    //Calculate condition system
    this->CalculateConditionSystem( LocalSystem, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::InitializeSystemMatrices( 
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    Flags& rCalculationFlags 
    )
{
    const unsigned int dimension = 2;
    const unsigned int num_slave_nodes = GetGeometry( ).PointsNumber( );
    const unsigned int num_master_nodes_total = CalculateTotalNumberOfMasterNodes( );

    const unsigned int condition_size = this->CalculateConditionSize( );

    if ( rCalculationFlags.Is( MortarContact2DCondition::COMPUTE_LHS_MATRIX ) ) // Calculation of the matrix is required
    {
        noalias( rLeftHandSideMatrix ) = ZeroMatrix( condition_size, condition_size ); // Resetting LHS
    }

    // Resizing as needed the RHS
    if ( rCalculationFlags.Is( MortarContact2DCondition::COMPUTE_RHS_VECTOR ) ) // Calculation of the matrix is required
    {
        rRightHandSideVector = ZeroVector( condition_size ); // Resetting RHS
    }
}

/***********************************************************************************/
/***********************************************************************************/

const unsigned int& MortarContact2DCondition::CalculateConditionSize( )
{
    const unsigned int dimension = 2;

    unsigned int condition_size = 0;
        
    // Master and slave displacement nodes
    condition_size += this->CalculateTotalNumberOfMasterNodes( );
    condition_size += mThisMasterElements.size( ) * GetGeometry().PointsNumber( );
        
    for ( unsigned int i_master_elem = 0; i_master_elem < mThisMasterElements.size( ); ++i_master_elem )
    {
        condition_size += CalculateNumberOfActiveNodesInContactPair( i_master_elem );
    }
  
    condition_size *= dimension;
    
    return condition_size;
}

/***********************************************************************************/
/***********************************************************************************/

const unsigned int& MortarContact2DCondition::CalculateTotalNumberOfMasterNodes( )
{
    unsigned int num_master_nodes_total = 0;
    for ( unsigned int i_master_elem = 0; i_master_elem < mThisMasterElements.size( ); ++i_master_elem )
    {
        num_master_nodes_total += mThisMasterElements[i_master_elem]->GetGeometry().PointsNumber();
    }

    return num_master_nodes_total;
}

/***********************************************************************************/
/***********************************************************************************/

const unsigned int& MortarContact2DCondition::CalculateNumberOfActiveNodesInContactPair( const unsigned int& rPairIndex )
{
    const unsigned int num_slave_nodes = GetGeometry().PointsNumber();
    const contact_container& current_container = ( *( this->GetValue( CONTACT_CONTAINERS ) ) )[rPairIndex];
    
    unsigned int num_active_nodes = 0;
    for ( unsigned int iSlave = 0; iSlave < num_slave_nodes; ++iSlave )
    {
        if( current_container.active_nodes[iSlave] )
        {
            ++num_active_nodes;
        }
    }

    return num_active_nodes;
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateConditionSystem( 
    LocalSystemComponents& rLocalSystem,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
    
    // Create and initialize condition variables:
    GeneralVariables Variables;

    // Reading integration points
    bool calculate_slave_contributions = true;
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry( ).IntegrationPoints( mThisIntegrationMethod );

    for (unsigned int i_master_elem = 0; i_master_elem < mThisMasterElements.size( ); ++i_master_elem)
    {
        // Initialize general variables for the current master element
        this->InitializeGeneralVariables( Variables, i_master_elem );
        this->InitializeActiveSet( Variables );

        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size( ); PointNumber++ )
        {
            double integration_weight = integration_points[PointNumber].Weight( );

            this->CalculateKinematics( Variables, PointNumber );

            this->InitializeConditionMatrices( Variables,
                                               integration_weight,
                                               PointNumber,
                                               calculate_slave_contributions);

            // To prevent recalculation of slave contributions
            calculate_slave_contributions = false;

            // Calculation of the matrix is required
            if ( rLocalSystem.CalculationFlags.Is( MortarContact2DCondition::COMPUTE_LHS_MATRIX ) ||
                 rLocalSystem.CalculationFlags.Is( MortarContact2DCondition::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
            {
                // Contributions to stiffness matrix calculated on the reference config
                this->CalculateAndAddLHS( rLocalSystem, Variables, integration_weight );
            }

            // Calculation of the vector is required
            if ( rLocalSystem.CalculationFlags.Is( MortarContact2DCondition::COMPUTE_RHS_VECTOR ) ||
                 rLocalSystem.CalculationFlags.Is( MortarContact2DCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
            {
                // Contribution to previous step contact force and residuals vector
                this->CalculateAndAddRHS( rLocalSystem, Variables, integration_weight );
            }
        }
    }
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::InitializeGeneralVariables(
    GeneralVariables& rVariables,
    const unsigned int& rMasterElementIndex 
    )
{
    // Working space dimension
    const unsigned int dimension = 2;

    // Slave segment info
    GeometryType& current_slave_element = this->GetGeometry();
    const unsigned int number_of_nodes_slave = current_slave_element.PointsNumber();
    const unsigned int local_dimension_slave = current_slave_element.LocalSpaceDimension();

    // Master segment info
    GeometryType& current_master_element = mThisMasterElements[rMasterElementIndex]->GetGeometry();
    const unsigned int number_of_nodes_master = current_master_element.PointsNumber();
    const unsigned int local_dimension_master = current_master_element.LocalSpaceDimension();

    // Slave element info
    rVariables.Initialize(local_dimension_master,
                          number_of_nodes_master,
                          local_dimension_slave,
                          number_of_nodes_slave,
                          dimension 
                         );

    rVariables.SetMasterElement( current_master_element );
    rVariables.SetMasterElementIndex( rMasterElementIndex );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::InitializeActiveSet( GeneralVariables& rVariables )
{
    std::vector< unsigned int > active_nodes;
    std::vector< unsigned int > inactive_nodes;
  
    this->DetermineActiveAndInactiveSets( active_nodes, inactive_nodes, rVariables.GetMasterElementIndex( ) );
  
    rVariables.SetActiveSet( active_nodes );
    rVariables.SetInactiveSet( inactive_nodes );
}

/***********************************************************************************/
/***********************************************************************************/
  
void MortarContact2DCondition::DetermineActiveAndInactiveSets(
        std::vector<unsigned int>& rActiveNodes,
        std::vector<unsigned int>& rInactiveNodes,
        const unsigned int& rCondIndex
        )
{
    contact_container& current_container = ( *( this->GetValue( CONTACT_CONTAINERS ) ) )[rCondIndex];
    const unsigned int num_slave_nodes = GetGeometry( ).PointsNumber( );
  
    rActiveNodes.resize( 0 );
    rInactiveNodes.resize( 0 );
    
    for ( unsigned int iSlave = 0; iSlave < num_slave_nodes; ++iSlave )
    {
        if( current_container.active_nodes[iSlave] )
        {
            rActiveNodes.push_back( iSlave );
        }
        else
        {
            rInactiveNodes.push_back( iSlave );
        }
    }
}

//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************

void MortarContact2DCondition::CalculateKinematics( 
    GeneralVariables& rVariables,
    const double& rPointNumber 
    )
{
    KRATOS_TRY;
    
    /* DEFINITIONS */
    GeometryType& slave_nodes = GetGeometry( );
    GeometryType& master_nodes = rVariables.GetMasterElement( );
    const unsigned int number_of_master_nodes = master_nodes.PointsNumber( );
    const unsigned int number_of_slave_nodes = slave_nodes.PointsNumber( );

    /* RESIZE MATRICES AND VECTORS */
    rVariables.Phi_LagrangeMultipliers.resize( number_of_slave_nodes );
    rVariables.N_Master.resize( number_of_master_nodes );
    rVariables.N_Slave.resize( number_of_slave_nodes );

    rVariables.DPhi_De_LagrangeMultipliers.resize( number_of_slave_nodes, slave_nodes.LocalSpaceDimension( ) );
    rVariables.DN_De_Master.resize( number_of_master_nodes, master_nodes.LocalSpaceDimension( ) );
    rVariables.DN_De_Slave.resize( number_of_slave_nodes, slave_nodes.LocalSpaceDimension( ) );

    /*  POPULATE MATRICES AND VECTORS */
    for( unsigned int iNode = 0; iNode < number_of_master_nodes; ++iNode )
    {
        rVariables.N_Master[iNode] = master_nodes.ShapeFunctionValue( rPointNumber, iNode, mThisIntegrationMethod );
    }

    for( unsigned int iNode = 0; iNode < number_of_slave_nodes; ++iNode )
    {
        rVariables.N_Slave[iNode] = slave_nodes.ShapeFunctionValue( rPointNumber, iNode, mThisIntegrationMethod );
        rVariables.Phi_LagrangeMultipliers[iNode] = LagrangeMultiplierShapeFunctionValue( rPointNumber, iNode );
    }

    rVariables.DN_De_Master = master_nodes.ShapeFunctionLocalGradient( rPointNumber );
    rVariables.DN_De_Slave  = slave_nodes.ShapeFunctionLocalGradient( rPointNumber );
    rVariables.DPhi_De_LagrangeMultipliers = LagrangeMultiplierShapeFunctionLocalGradient( rPointNumber );

    slave_nodes.Jacobian( rVariables.j_Slave, mThisIntegrationMethod );

    rVariables.DetJSlave = slave_nodes.DeterminantOfJacobian( rPointNumber, mThisIntegrationMethod );

    this->CalculateNormalGapAtIntegrationPoint( rVariables, rPointNumber );

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateNormalGapAtIntegrationPoint(
    GeneralVariables& rVariables,
    const unsigned int& rPointNumber
    )
{
    // We interpolate the gap bewtween the nodes
    const std::vector<double> current_contact_gap = (*( this->GetValue( CONTACT_CONTAINERS ) ) )[rVariables.GetMasterElementIndex()].contact_gap;
        rVariables.IntegrationPointNormalGap[rPointNumber] = 0.0; 
  
    for ( unsigned int int_node = 0; int_node < rVariables.N_Master.size(); ++int_node )
    {
        rVariables.IntegrationPointNormalGap[rPointNumber] += rVariables.N_Master[int_node] * current_contact_gap[int_node]; 
    }
}
  
/***********************************************************************************/
/***********************************************************************************/
  
void MortarContact2DCondition::InitializeConditionMatrices(
    GeneralVariables& rVariables,
    double& rIntegrationWeight,
    const unsigned int& rPointNumber,
    const bool& rCalculateSlaveContributions 
    )
{
    // Working space dimension
    const unsigned int dimension = 2;

    // Slave segment info
    GeometryType& current_slave_element = this->GetGeometry( );
    const unsigned int number_of_nodes_slave = current_slave_element.PointsNumber( );
    const unsigned int local_dimension_slave = current_slave_element.LocalSpaceDimension( );

    // Master segment info
    const GeometryType& current_master_element = rVariables.GetMasterElement( );
    const unsigned int number_of_nodes_master = current_master_element.PointsNumber( );
    const unsigned int local_dimension_master = current_master_element.LocalSpaceDimension( );

    // Compute mortar condition matrices
    mThisMortarConditionMatrices.Initialize( rCalculateSlaveContributions,
                                             local_dimension_master,
                                             local_dimension_slave,
                                             number_of_nodes_master,
                                             number_of_nodes_slave,
                                             dimension );

    this->CalculateDAndM( rVariables, rIntegrationWeight, rCalculateSlaveContributions );
}

/***********************************************************************************/
/*************** METHODS TO CALCULATE THE CONTACT CONDITION MATRICES ***************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateDAndM( 
    GeneralVariables& rVariables,
    double& rIntegrationWeight,
    const bool& rCalculateSlaveContributions 
    )
{
    const unsigned int dimension = 2;
    const unsigned int num_master_nodes = rVariables.GetMasterElement( ).PointsNumber( );

    const Vector& N_s = rVariables.N_Slave;
    const Vector& N_m = rVariables.N_Master;
    const Vector& Phi = rVariables.Phi_LagrangeMultipliers;
    const double& J_s = rVariables.DetJSlave;
    
    MatrixType& D = mThisMortarConditionMatrices.D;
    MatrixType& M = mThisMortarConditionMatrices.M;

    // Inactive Slave DOFs
    for ( unsigned int i_inactive = 0; i_inactive < rVariables.GetInactiveSet().size( ); ++i_inactive )
    {
        const unsigned int iNode = rVariables.GetInactiveSet()[i_inactive];
        const unsigned int i = i_inactive * dimension;
        
        // tTo avoid recalculating the D matrix
        if( rCalculateSlaveContributions )
        {
            D( i    , i     ) = rIntegrationWeight * N_s( iNode ) * J_s;
            D( i + 1, i + 1 ) = D( i, i );
        }

        for ( unsigned int jNode = 0; jNode < num_master_nodes; ++jNode )
        {
            const unsigned int j = jNode * dimension;
            M( i  , j   ) = rIntegrationWeight * Phi( iNode ) * N_m( jNode ) * J_s;
            M( i+1, j+1 ) = M( i, j );
        }
     }

     // Active Slave DOFs
     for ( unsigned int i_active = 0; i_active < rVariables.GetActiveSet().size( ); ++i_active )
     {
         const unsigned int iNode = rVariables.GetActiveSet()[i_active];
         const unsigned int i = ( i_active + rVariables.GetInactiveSet().size( ) ) * dimension;    // active follows slave

         // To avoid recalculating the D matrix
         if( rCalculateSlaveContributions )
         {
             D( i    , i     ) = rIntegrationWeight * N_s( iNode ) * J_s;
             D( i + 1, i + 1 ) = D( i, i );
         }

         for ( unsigned int jNode = 0; jNode < num_master_nodes; ++jNode )
         {
             const unsigned int j = jNode * dimension;
             M( i    , j     ) = rIntegrationWeight * Phi( iNode ) * N_m( jNode ) * J_s;
             M( i + 1, j + 1 ) = M( i, j );
         }
     }
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateAndAddLHS(
    LocalSystemComponents& rLocalSystem,
    GeneralVariables& rVariables,
    double& rIntegrationWeight
    )
{
    // TODO Check LHS Variables and discuss with Vicente whether those are the final ones or not !!

    // Contact pair variables
    const unsigned int dimension = 2;
    const unsigned int num_master_nodes = rVariables.GetMasterElement( ).PointsNumber( );
    const unsigned int num_slave_nodes = GetGeometry( ).PointsNumber( );
    const unsigned int num_total_nodes = num_slave_nodes + num_master_nodes;

    if ( rLocalSystem.CalculationFlags.Is( MortarContact2DCondition::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
    {
        /* COMPONENT-WISE LHS MATRIX */
        const std::vector<Variable<MatrixType> >& rLeftHandSideVariables = rLocalSystem.GetLeftHandSideVariables( );
        bool calculated;

        for ( unsigned int i = 0; i < rLeftHandSideVariables.size( ); i++ )
        {
            calculated = false;

            if ( rLeftHandSideVariables[i] == MORTAR_CONTACT_OPERATOR )
            {
                MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrices( )[i];
                                MatrixType LHS_contact_pair = ZeroMatrix( rLeftHandSideMatrix.size1( ), rLeftHandSideMatrix.size2( ) );
                
                // Calculate
                this->CalculateAndAddMortarContactOperator( LHS_contact_pair, rVariables, rIntegrationWeight );

                // Assemble
                this->AssembleContactPairLHSToConditionSystem( rVariables.GetMasterElementIndex( ),
                                                               LHS_contact_pair,
                                                               rLeftHandSideMatrix );
                calculated = true;
            }

            if ( calculated == false )
            {
                KRATOS_THROW_ERROR( std::logic_error,  " CONDITION can not supply the required local system variable: ", rLeftHandSideVariables[i] );
            }
        }
    }
    else
    {
        /* SINGLE LHS MATRIX */
        MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix( );
        MatrixType LHS_contact_pair = ZeroMatrix( rLeftHandSideMatrix.size1( ), rLeftHandSideMatrix.size2( ) );
                
        // Calculate
        this->CalculateAndAddMortarContactOperator( LHS_contact_pair, rVariables, rIntegrationWeight );

        // Assemble
        this->AssembleContactPairLHSToConditionSystem( rVariables.GetMasterElementIndex( ),
                                                       LHS_contact_pair,
                                                       rLeftHandSideMatrix );

//      KRATOS_WATCH( rLeftHandSideMatrix )
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::AssembleContactPairLHSToConditionSystem(
    const unsigned int rPairIndex,
    MatrixType& rPairLHS,
    MatrixType& rConditionLHS 
    )
{
    const unsigned int dimension = 2;
    const unsigned int num_slave_nodes = GetGeometry( ).PointsNumber( );
    const unsigned int num_current_master_nodes = mThisMasterElements[rPairIndex]->GetGeometry( ).PointsNumber( );
        const unsigned int num_current_active_nodes = CalculateNumberOfActiveNodesInContactPair( rPairIndex );
        const unsigned int current_pair_size = dimension * ( num_current_master_nodes + num_slave_nodes +  num_current_active_nodes );
  
    // Find location of the piar's master DOFs in ConditionLHS
    unsigned int index_begin = 0;

    for ( unsigned int i_master_elem = 0; i_master_elem < rPairIndex - 1; ++i_master_elem )
    {
        index_begin += mThisMasterElements[i_master_elem]->GetGeometry( ).PointsNumber( );
        index_begin += num_slave_nodes;
        index_begin += CalculateNumberOfActiveNodesInContactPair( i_master_elem );
    }

        index_begin *= dimension;
    
    Matrix pair_block_in_rConditionLHS = subrange( rConditionLHS,
                                                   index_begin,
                                                   index_begin + current_pair_size,
                                                   index_begin,
                                                   index_begin + current_pair_size
                                                 );

    pair_block_in_rConditionLHS = rPairLHS;
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateAndAddRHS(
    LocalSystemComponents& rLocalSystem,
    GeneralVariables& rVariables,
    double& rIntegrationWeight 
    )
{
    // NOTE: For a standard LM method the RHS added is zero!!!
    if( rLocalSystem.CalculationFlags.Is( MortarContact2DCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
    {
        std::vector<VectorType>& rRightHandSideVectors = rLocalSystem.GetRightHandSideVectors();
        const std::vector< Variable< VectorType > >& rRightHandSideVariables = rLocalSystem.GetRightHandSideVariables();
        for( unsigned int i=0; i<rRightHandSideVariables.size(); i++ )
        {
                    this->CalculateAndAddInternalForces( rRightHandSideVectors[i], rVariables, rIntegrationWeight );
        }
    }
    else
    {
        VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();
        this->CalculateAndAddInternalForces( rRightHandSideVector, rVariables, rIntegrationWeight );
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateAndAddInternalForces(
        VectorType& rRightHandSideVector,
    GeneralVariables & rVariables,
    double& rIntegrationWeight
        )
{
    KRATOS_TRY;

    VectorType InternalForces = ZeroVector(this->CalculateConditionSize( ));
    noalias( rRightHandSideVector ) -= InternalForces;

    KRATOS_CATCH( "" );
}
  
/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::AssembleContactPairRHSToConditionSystem(
    const unsigned int rPairIndex,
    VectorType& rPairRHS,
    VectorType& rConditionRHS 
    )
{
    const unsigned int dimension = 2;
    const unsigned int num_slave_nodes = GetGeometry( ).PointsNumber( );
    const unsigned int num_current_master_nodes = mThisMasterElements[rPairIndex]->GetGeometry( ).PointsNumber( );
    const unsigned int num_current_active_nodes = CalculateNumberOfActiveNodesInContactPair( rPairIndex );
    const unsigned int current_pair_size = dimension * ( num_current_master_nodes + num_slave_nodes +  num_current_active_nodes );
  
    // Find location of the piar's master DOFs in ConditionLHS
    unsigned int index_begin = 0;

    for ( unsigned int i_master_elem = 0; i_master_elem < rPairIndex - 1; ++i_master_elem )
    {
        index_begin += mThisMasterElements[i_master_elem]->GetGeometry( ).PointsNumber( );
        index_begin += num_slave_nodes;
        index_begin += CalculateNumberOfActiveNodesInContactPair( i_master_elem );
    }

        index_begin *= dimension;
    
    Vector pair_block_in_rConditionRHS = subrange( rConditionRHS,
                                                   index_begin,
                                                   index_begin + current_pair_size
                                                 );

    pair_block_in_rConditionRHS = rPairRHS;
}

/***********************************************************************************/
/**************** AUXILLIARY METHODS FOR CONDITION LHS CONTRIBUTION ****************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateAndAddMortarContactOperator( 
    MatrixType& rLeftHandSideMatrix,
    GeneralVariables& rVariables,
    double& rIntegrationWeight 
    )
{
    KRATOS_TRY;
  
    // ********* DONE **********
  
    // Contact pair variables
    const unsigned int dimension = 2;
    const unsigned int num_slave_nodes = GetGeometry( ).PointsNumber( );
    const unsigned int num_master_nodes = rVariables.GetMasterElement( ).PointsNumber( );  
    const unsigned int num_total_nodes = num_slave_nodes + num_master_nodes;
    const unsigned int num_active_nodes = rVariables.GetActiveSet( ).size( );
    const unsigned int num_inactive_nodes = rVariables.GetInactiveSet( ).size( );

    const Matrix& D = mThisMortarConditionMatrices.D;
    const Matrix& M = mThisMortarConditionMatrices.M;
  
    // Calculate the blocks of B_co and B_co_transpose
    unsigned int i = 0, j = 0, j_node = 0;
    for ( unsigned int i_active = 0; i_active < num_active_nodes; ++i_active )
    {
        i = ( i_active + num_total_nodes ) * dimension;
        
        for ( unsigned int j_master = 0; j_master < num_master_nodes; ++j_master )
        {
                // fill the -M and -M' parts
            j = j_master * dimension;
                rLeftHandSideMatrix( i,     j     ) = -M( i, j );
                rLeftHandSideMatrix( i + 1, j + 1 ) = -M( i, j );
          
                rLeftHandSideMatrix( j,     i     ) = -M( i, j );
                rLeftHandSideMatrix( j + 1, i + 1 ) = -M( i, j );
        }
      
        // fill the D and D' parts
        j = ( i_active + num_master_nodes + num_inactive_nodes ) * dimension;
        rLeftHandSideMatrix( i,     j     ) = D( i_active, i_active );
        rLeftHandSideMatrix( i + 1, j + 1 ) = D( i_active, i_active );

        rLeftHandSideMatrix( j,     i     ) = D( i_active, i_active );
        rLeftHandSideMatrix( j + 1, i + 1 ) = D( i_active, i_active );
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY;

        const unsigned int num_slave_nodes = GetGeometry().PointsNumber();
        const std::vector<contact_container> all_conditions = *( this->GetValue( CONTACT_CONTAINERS ) );
        
    rResult.resize( 0, false );

    /* ORDER - [ MASTER, SLAVE_I, SLAVE_A, LAMBDA ] */
  
    for ( unsigned int i_cond = 0; all_conditions.size( ); ++i_cond )
    {
        const contact_container& current_container = all_conditions[i_cond];
        GeometryType& current_master = current_container.condition->GetGeometry( );
        
        std::vector< unsigned int > active_nodes;
        std::vector< unsigned int > inactive_nodes;
        this->DetermineActiveAndInactiveSets( active_nodes, inactive_nodes, i_cond );
      
        // Master Nodes Displacement Equation IDs
        for ( unsigned int iNode = 0; iNode < current_master.PointsNumber( ); iNode++ )
        {
            NodeType& master_node = current_master[iNode];
            rResult.push_back( master_node.GetDof( DISPLACEMENT_X ).EquationId( ) );
            rResult.push_back( master_node.GetDof( DISPLACEMENT_Y ).EquationId( ) );
        }
        
        // Inactive Nodes Displacement Equation IDs
        for ( unsigned int i_slave = 0; i_slave < inactive_nodes.size( ); ++i_slave )
        {
            NodeType& slave_node = GetGeometry()[ inactive_nodes[i_slave] ];
            rResult.push_back( slave_node.GetDof( DISPLACEMENT_X ).EquationId( ) );
            rResult.push_back( slave_node.GetDof( DISPLACEMENT_Y ).EquationId( ) );
        }

        // Active Nodes Displacement Equation IDs
        for ( unsigned int i_slave = 0; i_slave < active_nodes.size( ); ++i_slave )
        {
            NodeType& slave_node = GetGeometry()[ active_nodes[i_slave] ];
            rResult.push_back( slave_node.GetDof( DISPLACEMENT_X ).EquationId( ) );
            rResult.push_back( slave_node.GetDof( DISPLACEMENT_Y ).EquationId( ) );
        }

        // Active Nodes Lambda Equation IDs
        for ( unsigned int i_slave = 0; i_slave < active_nodes.size( ); ++i_slave )
        {
            NodeType& slave_node = GetGeometry()[ active_nodes[i_slave] ];
            rResult.push_back( slave_node.GetDof( LAGRANGE_MULTIPLIER_X ).EquationId( ) );
            rResult.push_back( slave_node.GetDof( LAGRANGE_MULTIPLIER_Y ).EquationId( ) );
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::GetDofList(
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    const unsigned int num_slave_nodes = GetGeometry().PointsNumber();
    const std::vector<contact_container> all_conditions = *( this->GetValue( CONTACT_CONTAINERS ) );
        
    rConditionalDofList.resize( 0 );

    /* ORDER - [ MASTER, SLAVE_I, SLAVE_A, LAMBDA ] */
  
    for ( unsigned int i_cond = 0; all_conditions.size( ); ++i_cond )
    {
        const contact_container& current_container = all_conditions[i_cond];
        GeometryType& current_master = current_container.condition->GetGeometry( );
        
        std::vector< unsigned int > active_nodes;
        std::vector< unsigned int > inactive_nodes;
        this->DetermineActiveAndInactiveSets( active_nodes, inactive_nodes, i_cond );
      
        // Master Nodes Displacement Equation IDs
        for ( unsigned int iNode = 0; iNode < current_master.PointsNumber( ); iNode++ )
        {
            NodeType& master_node = current_master[iNode];
            rConditionalDofList.push_back( master_node.pGetDof( DISPLACEMENT_X ) );
            rConditionalDofList.push_back( master_node.pGetDof( DISPLACEMENT_Y ) );
        }
        
        // Inactive Nodes Displacement Equation IDs
        for ( unsigned int i_slave = 0; i_slave < inactive_nodes.size( ); ++i_slave )
        {
            NodeType& slave_node = GetGeometry()[ inactive_nodes[i_slave] ];
            rConditionalDofList.push_back( slave_node.pGetDof( DISPLACEMENT_X ) );
            rConditionalDofList.push_back( slave_node.pGetDof( DISPLACEMENT_Y ) );
        }

        // Active Nodes Displacement Equation IDs
        for ( unsigned int i_slave = 0; i_slave < active_nodes.size( ); ++i_slave )
        {
            NodeType& slave_node = GetGeometry()[ active_nodes[i_slave] ];
            rConditionalDofList.push_back( slave_node.pGetDof( DISPLACEMENT_X ) );
            rConditionalDofList.push_back( slave_node.pGetDof( DISPLACEMENT_Y ) );
        }

        // Active Nodes Lambda Equation IDs
        for ( unsigned int i_slave = 0; i_slave < active_nodes.size( ); ++i_slave )
        {
            NodeType& slave_node = GetGeometry()[ active_nodes[i_slave] ];
            rConditionalDofList.push_back( slave_node.pGetDof( LAGRANGE_MULTIPLIER_X ) );
            rConditionalDofList.push_back( slave_node.pGetDof( LAGRANGE_MULTIPLIER_Y ) );
        }
    }

    KRATOS_CATCH( "" );
}

//******************************* GET DOUBLE VALUE *********************************/
/***********************************************************************************/

void MortarContact2DCondition::GetValueOnIntegrationPoints( 
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateOnIntegrationPoints( 
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;

    unsigned int number_of_integration_pts = GetGeometry( ).IntegrationPoints( mThisIntegrationMethod ).size( );
    if ( rOutput.size( ) != number_of_integration_pts )
    {
        rOutput.resize( number_of_integration_pts );
    }

    // TODO: ADD CONTENT!!!!!

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::GetNodalDeltaMovements( 
    Vector& rValues,
    const int& rNode 
    )
{
    unsigned int dimension = 2;

    if( rValues.size() != dimension )
    {
        rValues.resize( dimension );
    }

    rValues = ZeroVector( dimension );

    Vector CurrentValueVector = ZeroVector(3);
    CurrentValueVector = GetCurrentValue( DISPLACEMENT, CurrentValueVector, rNode );

    Vector PreviousValueVector = ZeroVector(3);
    CurrentValueVector = GetPreviousValue( DISPLACEMENT, CurrentValueVector, rNode );

    rValues[0] = CurrentValueVector[0] - PreviousValueVector[0];
    rValues[1] = CurrentValueVector[1] - PreviousValueVector[1];

    if( dimension == 3 ) // We are working in 2D, this is not supposed to do anything
    {
        rValues[2] = CurrentValueVector[2] - PreviousValueVector[2];
    }
}

//************************************************************************************
//************************************************************************************

Vector& MortarContact2DCondition::GetCurrentValue( 
    const Variable<array_1d<double,3> >& rVariable,
    Vector& rValue,
    const unsigned int& rNode
    )
{
    KRATOS_TRY;

    const unsigned int dimension = 2;

    array_1d<double, 3> ArrayValue;
    ArrayValue = GetGeometry( )[rNode].FastGetSolutionStepValue( rVariable );

    if( rValue.size() != dimension )
    {
        rValue.resize( dimension, false );
    }

    for( unsigned int i = 0; i < dimension; i++ )
    {
        rValue[i] = ArrayValue[i];
    }

    return rValue;

    KRATOS_CATCH( "" );
}

//************************************************************************************
//************************************************************************************

Vector& MortarContact2DCondition::GetPreviousValue( 
    const Variable<array_1d<double,3> >& rVariable,
    Vector& rValue,
    const unsigned int& rNode
    )
{
    KRATOS_TRY;

    const unsigned int dimension = 2;

    array_1d<double, 3> ArrayValue;
    ArrayValue = GetGeometry( )[rNode].FastGetSolutionStepValue( rVariable, 1 );

    if( rValue.size() != dimension )
    {
        rValue.resize( dimension, false );
    }

    for( unsigned int i = 0; i < dimension; i++ )
    {
        rValue[i] = ArrayValue[i];
    }

    return rValue;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
}

void MortarContact2DCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
}

} // Namespace Kratos