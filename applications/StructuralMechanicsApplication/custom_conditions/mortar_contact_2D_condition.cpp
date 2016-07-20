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
#include "custom_utilities/contact_utilities.h"
#include "utilities/math_utils.h"
#include "structural_mechanics_application.h"
#include "structural_mechanics_application_variables.h"

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
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();

    //DO NOT ADD DOFS HERE!!!
}

/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

MortarContact2DCondition::MortarContact2DCondition( 
    MortarContact2DCondition const& rOther ) :
    Condition( rOther )
{
    //DO NOT ADD DOFS HERE!!!
}

/************************************* OPERATIONS **********************************/
/***********************************************************************************/

Condition::Pointer MortarContact2DCondition::Create( 
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new MortarContact2DCondition( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer MortarContact2DCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer( new MortarContact2DCondition( NewId, pGeom, pProperties ) );
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

    if( this->Has(INTEGRATION_ORDER_CONTACT) )
    {
        if (this->GetValue(INTEGRATION_ORDER_CONTACT) == 1)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
        }
        else if (this->GetValue(INTEGRATION_ORDER_CONTACT) == 2)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
        }
        else if (this->GetValue(INTEGRATION_ORDER_CONTACT) == 3)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
        }
        else if (this->GetValue(INTEGRATION_ORDER_CONTACT) == 4)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_4;
        }
        else if (this->GetValue(INTEGRATION_ORDER_CONTACT) == 5)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_5;
        }
        else
        {
            std::cout << "The number of integration points is not defined.  INTEGRATION_ORDER_CONTACT: "<< this->GetValue(INTEGRATION_ORDER_CONTACT) << std::endl;
            std::cout << "Options are: 1, 2, 3, 4, 5  " << std::endl;
            std::cout << "Taking default number of integration points (INTEGRATION_ORDER_CONTACT = 1)  " << std::endl;
            mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
        }
    }
    else
    {
        mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // Populate the vector of master elements
    std::vector<contact_container> * all_containers = this->GetValue(CONTACT_CONTAINERS);
    mThisMasterElements.clear();
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
    
    // Populate the vector of master elements
    std::vector<contact_container> *& all_containers = this->GetValue(CONTACT_CONTAINERS);
    
    double ActiveCheckFactor = 0.2;
    if( GetProperties().Has(ACTIVE_CHECK_FACTOR) )
    {
        ActiveCheckFactor = GetProperties().GetValue(ACTIVE_CHECK_FACTOR);
    }
    
    for ( unsigned int i_cond = 0; i_cond < all_containers->size(); ++i_cond )
    {
        Condition::Pointer pCond = (*all_containers)[i_cond].condition;
    
        // Fill the condition
        ContactUtilities::ContactContainerFiller((*all_containers)[i_cond], pCond->GetGeometry().Center(), GetGeometry(), pCond->GetGeometry(), 
                                                 this->GetValue(NORMAL), pCond->GetValue(NORMAL), ActiveCheckFactor, mThisIntegrationMethod);
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    // TODO: Add things if needed
}

/***********************************************************************************/
/***********************************************************************************/

const Vector MortarContact2DCondition::LagrangeMultiplierShapeFunctionValue(const double xi_local)
{
    // NOTE: For more information look at the thesis of Popp page 93/236
    
    const unsigned int num_slave_nodes = GetGeometry( ).PointsNumber();
    Vector Phi = ZeroVector( num_slave_nodes );

    if (num_slave_nodes == 2) // First order
    {
        Phi[0] = ( 0.5 * ( 1.0 - 3.0 * xi_local ) );
        Phi[1] = ( 0.5 * ( 1.0 + 3.0 * xi_local ) );
    }
    else if (num_slave_nodes == 3) // Second order
    {
        array_1d<double,3> aux_coordinates = ZeroVector(3);
        aux_coordinates[0] = xi_local;
        Vector Ncontainer = ZeroVector(3);
        Ncontainer = GetGeometry().ShapeFunctionsValues(Ncontainer, aux_coordinates);

        Phi[0] = Ncontainer(0) -  3.0/4.0 * Ncontainer(2) + 0.5;
        Phi[1] = Ncontainer(1) -  3.0/4.0 * Ncontainer(2) + 0.5;
        Phi[2] = 5.0/2.0 * Ncontainer(2) - 1.0;
    }
        
    return Phi;
}

/***********************************************************************************/
/***********************************************************************************/

const Matrix MortarContact2DCondition::LagrangeMultiplierShapeFunctionLocalGradient( const double xi_local )
{
    // For 2D2N elements as presented in Popp's and Gitterle's work
    const unsigned int num_slave_nodes = GetGeometry( ).PointsNumber( );
    const unsigned int local_dimension_slave = GetGeometry( ).LocalSpaceDimension( );
    Matrix DPhi_De = ZeroMatrix( num_slave_nodes, local_dimension_slave );

    if (num_slave_nodes == 2) // First order
    {
        DPhi_De( 0, 0 ) = - 3.0 / 2.0;
        DPhi_De( 1, 0 ) = + 3.0 / 2.0;
    }
    else if (num_slave_nodes == 3) // Second order
    {
        array_1d<double,3> aux_coordinates = ZeroVector(3);
        aux_coordinates[0] = xi_local;
        Matrix DNcontainer = ZeroMatrix(3, 1);
        DNcontainer = GetGeometry().ShapeFunctionsLocalGradients( DNcontainer , aux_coordinates );
        
        DPhi_De( 0, 0 ) = DNcontainer(0, 0) -  3.0/4.0 * DNcontainer(2, 0);
        DPhi_De( 1, 0 ) = DNcontainer(1, 0) -  3.0/4.0 * DNcontainer(2, 0);
        DPhi_De( 2, 0 ) = 5.0/2.0 * DNcontainer(2, 0);
    }

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
    LocalSystem.CalculationFlags.Set(MortarContact2DCondition::COMPUTE_LHS_MATRIX_WITH_COMPONENTS, true);
    LocalSystem.CalculationFlags.Set(MortarContact2DCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS, true);

    //Initialize sizes for the system components:
    if ( rLHSVariables.size( ) != rLeftHandSideMatrices.size( ) )
    {
        rLeftHandSideMatrices.resize( rLHSVariables.size( ) );
    }

    if ( rRHSVariables.size( ) != rRightHandSideVectors.size( ) )
    {
        rRightHandSideVectors.resize( rRHSVariables.size( ) );
    }

    LocalSystem.CalculationFlags.Set(MortarContact2DCondition::COMPUTE_LHS_MATRIX, true);
    for ( unsigned int i = 0; i < rLeftHandSideMatrices.size( ); i++ )
    {
        // Note: rRightHandSideVectors.size() > 0
        this->InitializeSystemMatrices( rLeftHandSideMatrices[i], rRightHandSideVectors[0],LocalSystem.CalculationFlags );
    }

    LocalSystem.CalculationFlags.Set( MortarContact2DCondition::COMPUTE_RHS_VECTOR, true );
    LocalSystem.CalculationFlags.Set( MortarContact2DCondition::COMPUTE_LHS_MATRIX, false ); // Temporarily only
    for ( unsigned int i = 0; i < rRightHandSideVectors.size( ); i++ )
    {
        // Note: rLeftHandSideMatrices.size() > 0
        this->InitializeSystemMatrices( rLeftHandSideMatrices[0], rRightHandSideVectors[i], LocalSystem.CalculationFlags  );
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
    
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContact2DCondition::COMPUTE_LHS_MATRIX, true );
    LocalSystem.CalculationFlags.Set( MortarContact2DCondition::COMPUTE_RHS_VECTOR, true );

    // Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );
    
    // Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix( rLeftHandSideMatrix );
    LocalSystem.SetRightHandSideVector( rRightHandSideVector );

    // Calculate condition system
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
    LocalSystem.CalculationFlags.Set( MortarContact2DCondition::COMPUTE_LHS_MATRIX, true );

    VectorType RightHandSideVector = Vector( );

    // Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, RightHandSideVector, LocalSystem.CalculationFlags );

    // Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix( rLeftHandSideMatrix );
    LocalSystem.SetRightHandSideVector( RightHandSideVector );

    // Calculate condition system
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
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContact2DCondition::COMPUTE_LHS_MATRIX, true );
    LocalSystem.CalculationFlags.Set( MortarContact2DCondition::COMPUTE_LHS_MATRIX_WITH_COMPONENTS, true );

    VectorType RightHandSideVector = Vector( );

    // Initialize size for the system components
    for( unsigned int i = 0; i < rLeftHandSideMatrices.size( ); i++ )
    {
        this->InitializeSystemMatrices( rLeftHandSideMatrices[i], RightHandSideVector, LocalSystem.CalculationFlags );
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
    LocalSystem.CalculationFlags.Set( MortarContact2DCondition::COMPUTE_RHS_VECTOR, true);

    MatrixType LeftHandSideMatrix = Matrix( );

    // Initialize size for the system components
    this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVector,LocalSystem.CalculationFlags);

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
    LocalSystem.CalculationFlags.Set( MortarContact2DCondition::COMPUTE_RHS_VECTOR, true );
    LocalSystem.CalculationFlags.Set( MortarContact2DCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS, true );

    MatrixType LeftHandSideMatrix = Matrix( );

    // Initialize size for the system components
    for( unsigned int i = 0; i < rRightHandSideVectors.size(); i++ )
    {
        this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVectors[i], LocalSystem.CalculationFlags );
    }

    // Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix( LeftHandSideMatrix );
    LocalSystem.SetRightHandSideVectors( rRightHandSideVectors );

    // Calculate condition system
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
    const unsigned int condition_size = this->CalculateConditionSize( );
    
    // Resizing as needed the LHS
    if ( rCalculationFlags.Is( MortarContact2DCondition::COMPUTE_LHS_MATRIX ) ) // Calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != condition_size )
        {
            rLeftHandSideMatrix.resize( condition_size, condition_size, false );
        }
        noalias( rLeftHandSideMatrix ) = ZeroMatrix( condition_size, condition_size ); // Resetting LHS
    }

    // Resizing as needed the RHS
    if ( rCalculationFlags.Is( MortarContact2DCondition::COMPUTE_RHS_VECTOR ) ) // Calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != condition_size )
        {
            rRightHandSideVector.resize( condition_size, false );
        }
        rRightHandSideVector = ZeroVector( condition_size ); // Resetting RHS
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateMassMatrix( 
    MatrixType& rMassMatrix, 
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
    
    rMassMatrix.resize(0, 0, false);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateDampingMatrix( 
    MatrixType& rDampingMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    rDampingMatrix.resize(0, 0, false);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

const unsigned int MortarContact2DCondition::CalculateConditionSize( )
{
    const unsigned int dimension = 2;

    unsigned int condition_size = 0;
        
    // Master and slave displacement nodes
    condition_size += this->CalculateTotalNumberOfMasterNodes( );
        
    for ( unsigned int i_master_elem = 0; i_master_elem < mThisMasterElements.size( ); ++i_master_elem )
    {
        condition_size += mThisMasterElements[i_master_elem]->GetGeometry( ).PointsNumber( );
        condition_size += GetGeometry( ).PointsNumber( );
    }
  
    condition_size *= dimension;
    
    return condition_size;
}

/***********************************************************************************/
/***********************************************************************************/

const unsigned int MortarContact2DCondition::CalculateTotalNumberOfMasterNodes( )
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

const unsigned int MortarContact2DCondition::CalculateNumberOfActiveNodesInContactPair( const unsigned int& rPairIndex )
{
    const unsigned int num_slave_nodes = GetGeometry().PointsNumber();
    
    unsigned int num_active_nodes = 0;
    for ( unsigned int iSlave = 0; iSlave < num_slave_nodes; ++iSlave )
    {
        if( GetGeometry()[iSlave].Is(ACTIVE) == true )
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

    // Working space dimension
    const unsigned int dimension = 2;

    // Slave segment info
    const unsigned int number_of_nodes_slave     = GetGeometry().PointsNumber( );
    const unsigned int number_of_elements_master = mThisMasterElements.size( );
    const unsigned int local_dimension_slave     = GetGeometry().LocalSpaceDimension( );
    
    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry( ).IntegrationPoints( mThisIntegrationMethod );
    
    MortarWeightedGaps WeightedGaps;
    for (unsigned int iNode = 0; iNode < number_of_nodes_slave; iNode++)
    {
        GetGeometry()[iNode].GetValue(WEIGHTED_GAP) = 1.0e9; // Value extremely large 
    }
    
    for (unsigned int i_master_elem = 0; i_master_elem < number_of_elements_master; ++i_master_elem)
    {
        // Initialize general variables for the current master element
        this->InitializeGeneralVariables( Variables, i_master_elem );
        this->InitializeActiveSet( Variables );

        // Master segment info
        const GeometryType& current_master_element = Variables.GetMasterElement( );
        const unsigned int number_of_nodes_master = current_master_element.PointsNumber( );
        const unsigned int local_dimension_master = current_master_element.LocalSpaceDimension( );

        // Compute mortar condition matrices
        MortarConditionMatrices ThisMortarConditionMatrices;
        ThisMortarConditionMatrices.Initialize( local_dimension_master, local_dimension_slave,number_of_nodes_master, number_of_nodes_slave, dimension );
        
        // Weighted gaps
        WeightedGaps.Initialize();
        VectorType& gn = WeightedGaps.gn;
        
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            const double IntegrationWeight = integration_points[PointNumber].Weight( );
            this->CalculateKinematics( Variables, PointNumber, i_master_elem );
            this->CalculateDAndM( Variables, IntegrationWeight, ThisMortarConditionMatrices, gn );
        }
        
        for (unsigned int iNode = 0; iNode < number_of_nodes_slave; iNode++)
        {
            if (GetGeometry()[iNode].GetValue(WEIGHTED_GAP) > gn[iNode])
            {
                GetGeometry()[iNode].GetValue(WEIGHTED_GAP) = gn[iNode]; // Saving the smallest one
            }
        }
        
        // Calculation of the matrix is required
        if ( rLocalSystem.CalculationFlags.Is( MortarContact2DCondition::COMPUTE_LHS_MATRIX ) ||
                rLocalSystem.CalculationFlags.Is( MortarContact2DCondition::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
        {
            // Contributions to stiffness matrix calculated on the reference config
            this->CalculateAndAddLHS( rLocalSystem, Variables, ThisMortarConditionMatrices );
        }

        // Calculation of the vector is required
        if ( rLocalSystem.CalculationFlags.Is( MortarContact2DCondition::COMPUTE_RHS_VECTOR ) ||
                rLocalSystem.CalculationFlags.Is( MortarContact2DCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
        {
            // Contribution to previous step contact force and residuals vector
            this->CalculateAndAddRHS( rLocalSystem, Variables, ThisMortarConditionMatrices, gn );
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

    // Resize according to the integration order
    const unsigned int& integration_point_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );
    
    // Slave segment info
    GeometryType& current_slave_element = this->GetGeometry();
    const unsigned int number_of_nodes_slave = current_slave_element.PointsNumber();
    const unsigned int local_dimension_slave = current_slave_element.LocalSpaceDimension();

    // Master segment info
    GeometryType& current_master_element = mThisMasterElements[rMasterElementIndex]->GetGeometry();
    const unsigned int number_of_nodes_master = current_master_element.PointsNumber();
    const unsigned int local_dimension_master = current_master_element.LocalSpaceDimension();

    // Slave element info
    rVariables.Initialize(local_dimension_master, number_of_nodes_master, local_dimension_slave, number_of_nodes_slave, dimension, integration_point_number );

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
        const unsigned int& rPairIndex
        )
{
    const unsigned int num_slave_nodes = GetGeometry( ).PointsNumber( );
  
    rActiveNodes.resize( 0, false );
    rInactiveNodes.resize( 0, false );
    
    for ( unsigned int iSlave = 0; iSlave < num_slave_nodes; ++iSlave )
    {
        if( GetGeometry()[iSlave].Is(ACTIVE) == true )
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
    const double& rPointNumber,
    const unsigned int& rPairIndex
    )
{
    KRATOS_TRY;
    
    /* DEFINITIONS */
    GeometryType& slave_nodes  = GetGeometry( );
    GeometryType& master_nodes = rVariables.GetMasterElement( );
    const unsigned int number_of_master_nodes = master_nodes.PointsNumber( );
    const unsigned int number_of_slave_nodes  =  slave_nodes.PointsNumber( );

    /* LOCAL COORDINATES */ // TODO: Think to move this, because we repeat it every time
    const contact_container& current_container = ( *( this->GetValue( CONTACT_CONTAINERS ) ) )[rPairIndex];
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry( ).IntegrationPoints(mThisIntegrationMethod);
    const PointType& eta = integration_points[rPointNumber].Coordinates();
    const double xi_local_slave  =  0.5 * (1.0 - eta[0]) * current_container.local_coordinates_slave[0]  + 0.5 * (1.0 + eta[0]) * current_container.local_coordinates_slave[1];
    Point<3> local_point;
    local_point.Coordinates() = ZeroVector(3);
    local_point.Coordinate(1) = xi_local_slave;
    
    /* RESIZE MATRICES AND VECTORS */
    rVariables.Phi_LagrangeMultipliers.resize( number_of_slave_nodes, false );
    rVariables.N_Master.resize( number_of_master_nodes, false );
    rVariables.N_Slave.resize( number_of_slave_nodes, false );

    rVariables.DN_De_Master.resize( number_of_master_nodes, master_nodes.LocalSpaceDimension( ), false );
    rVariables.DN_De_Slave.resize( number_of_slave_nodes, slave_nodes.LocalSpaceDimension( ), false );
    rVariables.DPhi_De_LagrangeMultipliers.resize( number_of_slave_nodes, slave_nodes.LocalSpaceDimension( ), false );
    
    /*  POPULATE MATRICES AND VECTORS */
    
    /// SLAVE CONDITION ///
    // SHAPE FUNCTIONS 
    rVariables.N_Slave = slave_nodes.ShapeFunctionsValues( rVariables.N_Slave, local_point.Coordinates() );
    rVariables.Phi_LagrangeMultipliers= LagrangeMultiplierShapeFunctionValue( xi_local_slave );
    // SHAPE FUNCTION DERIVATIVES
    rVariables.DN_De_Slave  =  slave_nodes.ShapeFunctionsLocalGradients( rVariables.DN_De_Slave , local_point );
    rVariables.DPhi_De_LagrangeMultipliers = LagrangeMultiplierShapeFunctionLocalGradient( xi_local_slave );
    
    // MASTER CONDITION
    this->MasterShapeFunctionValue( rVariables, eta, local_point);

    /* CALCULATE JACOBIAN AND JACOBIAN DETERMINANT */
    slave_nodes.Jacobian( rVariables.j_Slave, local_point.Coordinates() );
    rVariables.DetJSlave = (current_container.local_coordinates_slave[1] - current_container.local_coordinates_slave[0])/2.0 * rVariables.j_Slave(0, 0); // TODO: Check if it is correct

    KRATOS_CATCH( "" );
}
 
/***********************************************************************************/
/*************** METHODS TO CALCULATE THE CONTACT CONDITION MATRICES ***************/
/***********************************************************************************/

void MortarContact2DCondition::MasterShapeFunctionValue(
    GeneralVariables& rVariables,
    const PointType& rSlaveIntegrationPoint,
    const PointType& local_point 
    )
{
    rVariables.N_Master = ZeroVector(2);
    
    rVariables.IntegrationPointNormalGap = 0.0;
    PointType projected_gp_global;
    GeometryType& master_seg = rVariables.GetMasterElement( );

    array_1d<double,3> normal = ZeroVector(3);
    for( unsigned int iNode = 0; iNode < GetGeometry().PointsNumber(); ++iNode )
    {
        normal -= rVariables.N_Slave[iNode] * GetGeometry()[iNode].GetValue(NORMAL); // The opposite direction
    }
    
//     normal = normal/norm_2(normal); // It is suppossed to be already unitary (just in case)
    
    // Doing calculations with eta
    GeometryType::CoordinatesArrayType slave_gp_global;
//     this->GetGeometry( ).GlobalCoordinates( slave_gp_global, rSlaveIntegrationPoint );
    this->GetGeometry( ).GlobalCoordinates( slave_gp_global, local_point );
    
    double aux_dist = 0.0;
    ContactUtilities::ProjectDirection( master_seg, slave_gp_global, projected_gp_global, aux_dist, normal );
    
    GeometryType::CoordinatesArrayType projected_gp_local;
    if( master_seg.IsInside( projected_gp_global.Coordinates( ), projected_gp_local ) )
    {
        // SHAPE FUNCTIONS 
        rVariables.N_Master     = master_seg.ShapeFunctionsValues(         rVariables.N_Master,     projected_gp_local );         
        rVariables.DN_De_Master = master_seg.ShapeFunctionsLocalGradients( rVariables.DN_De_Master, projected_gp_local );         
    }
    
    // Doing calculations with local xi
    this->GetGeometry( ).GlobalCoordinates( slave_gp_global, local_point );
    ContactUtilities::ProjectDirection( master_seg, slave_gp_global, projected_gp_global, rVariables.IntegrationPointNormalGap, normal );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateDAndM( 
    GeneralVariables& rVariables,
    const double& rIntegrationWeight,
    MortarConditionMatrices& ThisMortarConditionMatrices,
    VectorType& gn
    )
{
    const unsigned int dimension = 2;
    const unsigned int num_slave_nodes  = GetGeometry().PointsNumber( );
    const unsigned int num_master_nodes = rVariables.GetMasterElement( ).PointsNumber( );

    const Vector& N_s = rVariables.N_Slave;
    const Vector& N_m = rVariables.N_Master;
    const Vector& Phi = rVariables.Phi_LagrangeMultipliers;
    const double& gap = rVariables.IntegrationPointNormalGap;
    const double& J_s = rVariables.DetJSlave;
    
    MatrixType& D  = ThisMortarConditionMatrices.D;
    MatrixType& M  = ThisMortarConditionMatrices.M;

    // For all the nodes
    for ( unsigned int iNode = 0; iNode < num_slave_nodes; ++iNode )
    {
        gn[iNode] += rIntegrationWeight * Phi( iNode ) * gap * J_s;
        
        const unsigned int i = iNode * dimension;
        
        D( i    , i     ) += rIntegrationWeight * Phi( iNode ) * N_s( iNode ) * J_s;
        D( i + 1, i + 1 )  = D( i, i ); // NOTE: We can avoid this, because we are storing twice the same value
         
        for ( unsigned int jNode = 0; jNode < num_master_nodes; ++jNode )
        {
            const unsigned int j = jNode * dimension;
            M( i    , j     ) += rIntegrationWeight * Phi( iNode ) * N_m( jNode ) * J_s;
            M( i + 1, j + 1 )  = M( i, j ); // NOTE: We can avoid this, because we are storing twice the same value
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateAndAddLHS(
    LocalSystemComponents& rLocalSystem,
    GeneralVariables& rVariables,
    const MortarConditionMatrices& ThisMortarConditionMatrices
    )
{
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
                MatrixType LHS_contact_pair = ZeroMatrix(rLeftHandSideMatrix.size1(), rLeftHandSideMatrix.size2());
                
                // Calculate
                this->CalculateAndAddMortarContactOperator( LHS_contact_pair, rVariables, ThisMortarConditionMatrices );

                // Assemble
                this->AssembleContactPairLHSToConditionSystem( rVariables.GetMasterElementIndex( ), LHS_contact_pair, rLeftHandSideMatrix );
                
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
        const unsigned int pair_size = 2 * ( rVariables.GetMasterElement( ).PointsNumber( ) + 2 * GetGeometry( ).PointsNumber( ) );       
        MatrixType LHS_contact_pair = ZeroMatrix(pair_size, pair_size);
        
        // Calculate
        this->CalculateAndAddMortarContactOperator( LHS_contact_pair, rVariables, ThisMortarConditionMatrices );
        
        // Assemble
        this->AssembleContactPairLHSToConditionSystem( rVariables.GetMasterElementIndex( ), LHS_contact_pair, rLeftHandSideMatrix );

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
    const unsigned int num_master_nodes = mThisMasterElements[rPairIndex]->GetGeometry( ).PointsNumber( );
    const unsigned int current_pair_size = dimension * ( num_master_nodes + 2 * num_slave_nodes );
  
    // Find location of the piar's master DOFs in ConditionLHS
    unsigned int index_begin = 0;

    for ( unsigned int i_master_elem = 0; i_master_elem < rPairIndex ; ++i_master_elem ) 
    {
        index_begin += mThisMasterElements[i_master_elem]->GetGeometry( ).PointsNumber( );
        index_begin += 2 * num_slave_nodes;
    }

    index_begin *= dimension;
  
    const unsigned int index_end = index_begin + current_pair_size;
    
    subrange( rConditionLHS, index_begin, index_end, index_begin, index_end) += rPairLHS;
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateAndAddRHS(
    LocalSystemComponents& rLocalSystem,
    GeneralVariables& rVariables,
    const MortarConditionMatrices& ThisMortarConditionMatrices,
    const VectorType& gn
    )
{
    if ( rLocalSystem.CalculationFlags.Is( MortarContact2DCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
    {
        /* COMPONENT-WISE RHS MATRIX */
        const std::vector<Variable<VectorType> >& rRightHandSideVariables = rLocalSystem.GetRightHandSideVariables( );
        bool calculated;

        for ( unsigned int i = 0; i < rRightHandSideVariables.size( ); i++ )
        {
            calculated = false;

            if ( rRightHandSideVariables[i] == MORTAR_CONTACT_OPERATOR )
            {
                VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVectors()[i];
                VectorType RHS_contact_pair = ZeroVector(rRightHandSideVector.size());
                
                // Calculate
                this->CalculateAndAddMortarContactOperator( RHS_contact_pair, rVariables, ThisMortarConditionMatrices, gn );

                // Assemble
                this->AssembleContactPairRHSToConditionSystem( rVariables.GetMasterElementIndex( ), rVariables.GetActiveSet( ), RHS_contact_pair, rRightHandSideVector );
                
                calculated = true;
            }

            if ( calculated == false )
            {
                KRATOS_THROW_ERROR( std::logic_error,  " CONDITION can not supply the required local system variable: ", rRightHandSideVariables[i] );
            }
        }
    }
    else 
    {
        /* SINGLE LHS MATRIX */
        VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();
        const unsigned int pair_size = 2 * ( rVariables.GetMasterElement( ).PointsNumber( ) + 2 * GetGeometry( ).PointsNumber( ) ); 
        VectorType RHS_contact_pair = ZeroVector(pair_size);
        
        // Calculate
        this->CalculateAndAddMortarContactOperator( RHS_contact_pair, rVariables, ThisMortarConditionMatrices, gn );
        
        // Assemble
        this->AssembleContactPairRHSToConditionSystem( rVariables.GetMasterElementIndex( ), rVariables.GetActiveSet( ), RHS_contact_pair, rRightHandSideVector );
    }
    
}
  
/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::AssembleContactPairRHSToConditionSystem(
    const unsigned int rPairIndex,
    const std::vector< unsigned int > active_set,
    VectorType& rPairRHS,
    VectorType& rConditionRHS 
    )
{
    const unsigned int dimension = 2;
    const unsigned int num_slave_nodes = GetGeometry( ).PointsNumber( );
    const unsigned int num_master_nodes = mThisMasterElements[rPairIndex]->GetGeometry( ).PointsNumber( );
    const unsigned int current_pair_size = dimension * ( num_master_nodes + 2 * num_slave_nodes );
  
    // Find location of the piar's master DOFs in ConditionRHS
    unsigned int index_begin = 0;

    for ( unsigned int i_master_elem = 0; i_master_elem < rPairIndex; ++i_master_elem )
    {
        index_begin += mThisMasterElements[i_master_elem]->GetGeometry( ).PointsNumber( );
        index_begin += 2 * num_slave_nodes;
    }

    index_begin *= dimension;
    const unsigned int index_end = index_begin + current_pair_size;
    
    // Computing subrange
    subrange( rConditionRHS, index_begin, index_end ) += rPairRHS;
}

/***********************************************************************************/
/**************** AUXILLIARY METHODS FOR CONDITION LHS CONTRIBUTION ****************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateAndAddMortarContactOperator( 
    MatrixType& rLeftHandSideMatrix,
    GeneralVariables& rVariables,
    const MortarConditionMatrices& ThisMortarConditionMatrices
    )
{
    KRATOS_TRY;
  
    // Contact pair variables
    const unsigned int dimension = 2;
    const unsigned int num_slave_nodes    = GetGeometry( ).PointsNumber( );
    const unsigned int num_master_nodes   = rVariables.GetMasterElement( ).PointsNumber( );  
    const unsigned int num_total_nodes    = num_slave_nodes + num_master_nodes;
    const unsigned int num_active_nodes   = rVariables.GetActiveSet( ).size( );
    const unsigned int num_inactive_nodes = rVariables.GetInactiveSet( ).size( );

    const Matrix& D = ThisMortarConditionMatrices.D;
    const Matrix& M = ThisMortarConditionMatrices.M;
    
    unsigned int i = 0, j = 0;
    for ( unsigned int i_inactive = 0; i_inactive < num_inactive_nodes; ++i_inactive )
    {
        i = ( i_inactive + num_total_nodes ) * dimension;
        rLeftHandSideMatrix( i    , i )    = 1.0;
        rLeftHandSideMatrix( i + 1, i + 1) = 1.0;
    }
    
    // Calculate the blocks of B_co and B_co_transpose
    i = 0;
    for ( unsigned int i_active = 0; i_active < num_active_nodes; ++i_active )
    {
        i = ( i_active + num_inactive_nodes + num_total_nodes ) * dimension;
        const unsigned int iNode = rVariables.GetActiveSet()[i_active];
        
        for ( unsigned int j_master = 0; j_master < num_master_nodes; ++j_master )
        {
            // Fill the -M and -M' parts
            j = j_master * dimension;
            double minus_M_ij = - M( iNode * dimension, j_master * dimension );
            rLeftHandSideMatrix( i,     j     ) += minus_M_ij;
            rLeftHandSideMatrix( i + 1, j + 1 ) += minus_M_ij;
                                                  
            rLeftHandSideMatrix( j    , i     ) += minus_M_ij;
            rLeftHandSideMatrix( j + 1, i + 1 ) += minus_M_ij;
        }
      
        // Fill the D and D' parts
        j = ( i_active + num_master_nodes + num_inactive_nodes ) * dimension;
        double D_ii = D( iNode * dimension, iNode * dimension );
        rLeftHandSideMatrix( i,     j     ) += D_ii;
        rLeftHandSideMatrix( i + 1, j + 1 ) += D_ii;
                                              
        rLeftHandSideMatrix( j,     i     ) += D_ii;
        rLeftHandSideMatrix( j + 1, i + 1 ) += D_ii;
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateAndAddMortarContactOperator( 
    VectorType& rRightHandSideVector,
    GeneralVariables& rVariables,
    const MortarConditionMatrices& ThisMortarConditionMatrices,
    const VectorType& gn
    )
{
    KRATOS_TRY;
  
    // Contact pair variables
    const unsigned int dimension = 2;
    const unsigned int num_slave_nodes    = GetGeometry( ).PointsNumber( );
    const unsigned int num_master_nodes   = rVariables.GetMasterElement( ).PointsNumber( );  
    const unsigned int num_total_nodes    = num_slave_nodes + num_master_nodes;
    const unsigned int num_active_nodes   = rVariables.GetActiveSet( ).size( );
    const unsigned int num_inactive_nodes = rVariables.GetInactiveSet( ).size( );

    const MatrixType& D  = ThisMortarConditionMatrices.D;
    const MatrixType& M  = ThisMortarConditionMatrices.M;
    
    // Calculate the block of r_co
    unsigned int j = 0;
    array_1d<double, 3> lagrange_multiplier; 
    for ( unsigned int i_active = 0; i_active < num_active_nodes; ++i_active )
    {
        const unsigned int iNode = rVariables.GetActiveSet()[i_active];
        noalias(lagrange_multiplier)  = GetGeometry( )[iNode].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER, 0); 
    
        // Fill the - lambda * - M part
        for ( unsigned int j_master = 0; j_master < num_master_nodes; ++j_master )
        {
            j = j_master * dimension;
            const double minus_M_ij = - M( iNode * dimension, j );
            rRightHandSideVector[ j     ] -= lagrange_multiplier[0] * minus_M_ij;
            rRightHandSideVector[ j + 1 ] -= lagrange_multiplier[1] * minus_M_ij;
        }
        
        // Fill the - lambda *  D part
        j = ( i_active + num_master_nodes + num_inactive_nodes ) * dimension;
        const double D_ii = D( iNode * dimension, iNode * dimension );
        rRightHandSideVector[ j     ] -= lagrange_multiplier[0] * D_ii;
        rRightHandSideVector[ j + 1 ] -= lagrange_multiplier[1] * D_ii;
        
        // Adding th ezero equality to the RHS
        const array_1d<double, 3> gap_decomp = gn[iNode] * GetGeometry()[iNode].GetValue(NORMAL); 
        j = (num_total_nodes + i_active) * dimension; 
        rRightHandSideVector[j    ] = 0.0;
        rRightHandSideVector[j + 1] = 0.0;
        
        // Adding the gap to the RHS
        const array_1d<double, 3> gap_decomp = gn[iNode] * GetGeometry()[iNode].GetValue(NORMAL); 
        j = (num_total_nodes + num_inactive_nodes + i_active) * dimension; 
        rRightHandSideVector[j    ] -= gap_decomp[0]; 
        rRightHandSideVector[j + 1] -= gap_decomp[1];
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
        
    rResult.resize( 0, false );

    const std::vector<contact_container> all_conditions = *( this->GetValue( CONTACT_CONTAINERS ) );
        
    /* ORDER - [ MASTER, SLAVE_I, SLAVE_A, LAMBDA ] */
    
    for ( unsigned int i_cond = 0;  i_cond < all_conditions.size( ); ++i_cond )
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
        
        // Inactive Nodes  Lambda Equation IDs
        for ( unsigned int i_slave = 0; i_slave < inactive_nodes.size( ); ++i_slave )
        {
            NodeType& slave_node = GetGeometry()[ inactive_nodes[i_slave] ];
            rResult.push_back( slave_node.GetDof( LAGRANGE_MULTIPLIER_X ).EquationId( ) );
            rResult.push_back( slave_node.GetDof( LAGRANGE_MULTIPLIER_Y ).EquationId( ) );
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
    
    rConditionalDofList.resize( 0 );

    const std::vector<contact_container> all_conditions = *( this->GetValue( CONTACT_CONTAINERS ) );
    
    /* ORDER - [ MASTER, SLAVE_I, SLAVE_A, LAMBDA ] */
    
    for ( unsigned int i_cond = 0; i_cond < all_conditions.size( ); ++i_cond )
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
        
        // Inactive Nodes Lambda Equation IDs
        for ( unsigned int i_slave = 0; i_slave < inactive_nodes.size( ); ++i_slave )
        {
            NodeType& slave_node = GetGeometry()[ inactive_nodes[i_slave] ];
            rConditionalDofList.push_back( slave_node.pGetDof( LAGRANGE_MULTIPLIER_X ) );
            rConditionalDofList.push_back( slave_node.pGetDof( LAGRANGE_MULTIPLIER_Y ) );
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

    const unsigned int number_of_integration_pts = GetGeometry( ).IntegrationPointsNumber( mThisIntegrationMethod );
    if ( rOutput.size( ) != number_of_integration_pts )
    {
        rOutput.resize( number_of_integration_pts, false );
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
