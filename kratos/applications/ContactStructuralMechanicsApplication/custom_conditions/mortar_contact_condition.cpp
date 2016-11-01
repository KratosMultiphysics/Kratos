// KRATOS  ___|  |       |       |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//           | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License: BSD License
//   license: structural_mechanics_application/license.txt
//
//  Main authors:  Mohamed Khalil
//                 Vicente Mataix Ferrándiz
//

// System includes

// External includes

// Project includes
#include "custom_conditions/mortar_contact_condition.h"
#include "custom_utilities/contact_utilities.h"
#include "utilities/math_utils.h"
#include "structural_mechanics_application.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_math_utilities.hpp"
#include <algorithm>

/* Includes of particular contact conditions */
#include "contact_2D_2N_2N.hpp"
#include "contact_3D_3N_3N.hpp"
#include "contact_3D_4N_4N.hpp"

// Logging format include
#include "custom_utilities/logging_settings.hpp"

namespace Kratos 
{
/**
 * Flags related to the condition computation
 */
KRATOS_CREATE_LOCAL_FLAG( MortarContactCondition, COMPUTE_RHS_VECTOR,     0 );
KRATOS_CREATE_LOCAL_FLAG( MortarContactCondition, COMPUTE_LHS_MATRIX,     1 );
KRATOS_CREATE_LOCAL_FLAG( MortarContactCondition, COMPUTE_RHS_VECTOR_WITH_COMPONENTS, 2 );
KRATOS_CREATE_LOCAL_FLAG( MortarContactCondition, COMPUTE_LHS_MATRIX_WITH_COMPONENTS, 3 );

/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

MortarContactCondition::MortarContactCondition(
    IndexType NewId,
    GeometryType::Pointer pGeometry ) :
    Condition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

MortarContactCondition::MortarContactCondition( 
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties ) :
    Condition( NewId, pGeometry, pProperties )
{
//     mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod(); // NOTE: If this is considered just one GP is considered

    const unsigned int dimension = GetGeometry( ).WorkingSpaceDimension( );
    const unsigned int local_dimension_slave = GetGeometry( ).LocalSpaceDimension( );
    const unsigned int number_of_slave_nodes = GetGeometry( ).PointsNumber( );
    
    InitializeIntegrationMethod(dimension, local_dimension_slave, number_of_slave_nodes); // NOTE: The integration method here defined depends of the properties of the parents elements, that's why I consider again this in the Initialize() 
    
    //DO NOT ADD DOFS HERE!!!
}

/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

MortarContactCondition::MortarContactCondition( 
    MortarContactCondition const& rOther ) :
    Condition( rOther )
{
    //DO NOT ADD DOFS HERE!!!
}

/************************************* OPERATIONS **********************************/
/***********************************************************************************/

Condition::Pointer MortarContactCondition::Create( 
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new MortarContactCondition( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer MortarContactCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer( new MortarContactCondition( NewId, pGeom, pProperties ) );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

MortarContactCondition::~MortarContactCondition( )
{
}

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

void MortarContactCondition::Initialize( ) 
{
    KRATOS_TRY;
    
    const unsigned int dimension = GetGeometry( ).WorkingSpaceDimension( );
    const unsigned int local_dimension_slave = GetGeometry( ).LocalSpaceDimension( );
    const unsigned int number_of_slave_nodes = GetGeometry( ).PointsNumber( );
    
    InitializeIntegrationMethod(dimension, local_dimension_slave, number_of_slave_nodes);
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContactCondition::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // Populate the vector of master elements
    std::vector<contact_container> * all_containers = this->GetValue(CONTACT_CONTAINERS);
    mThisMasterElements.clear();
    mThisMasterElements.resize( all_containers->size( ) );
    
    const double ActiveCheckFactor = GetProperties().GetValue(ACTIVE_CHECK_FACTOR); 
    
    for ( unsigned int i_cond = 0; i_cond < all_containers->size(); ++i_cond )
    {
        mThisMasterElements[i_cond] = (*all_containers)[i_cond].condition;

        // Fill the condition
        Condition::Pointer& pCond = mThisMasterElements[i_cond];
        
        // Initializes only via gap tolerance for the first step in the solution
        // Do the mortar segmentation in all time steps
        ContactUtilities::ContactContainerFiller((*all_containers)[i_cond], pCond->GetGeometry().Center(), GetGeometry(), pCond->GetGeometry(), 
                                                   this->GetValue(NORMAL), pCond->GetValue(NORMAL), ActiveCheckFactor);
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContactCondition::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContactCondition::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    // TODO: Add things if needed
}

/***********************************************************************************/
/***********************************************************************************/

IntegrationMethod MortarContactCondition::GetIntegrationMethod()
{   
    return mThisIntegrationMethod;
}

/***********************************************************************************/
/***********************************************************************************/

const Vector MortarContactCondition::LagrangeMultiplierShapeFunctionValue(
    const double xi_local,
    const double eta_local
    )
{
    // NOTE: For more information look at the thesis of Popp page 93/236
    
    // Working space dimension
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    
    // Number of slave nodes
    const unsigned int number_of_slave_nodes = GetGeometry( ).PointsNumber();
    
    Vector Phi = ZeroVector( number_of_slave_nodes );
            
    if (dimension == 2)
    {
        if (number_of_slave_nodes == 2) // First order
        {
            Phi[0] = ( 0.5 * ( 1.0 - 3.0 * xi_local ) );
            Phi[1] = ( 0.5 * ( 1.0 + 3.0 * xi_local ) );
        }
        else if (number_of_slave_nodes == 3) // Second order
        {
            array_1d<double,3> aux_coordinates = ZeroVector(3);
            aux_coordinates[0] = xi_local;
            Vector Ncontainer = ZeroVector(3);
            Ncontainer = GetGeometry().ShapeFunctionsValues(Ncontainer, aux_coordinates);

            Phi[0] = Ncontainer(0) -  3.0/4.0 * Ncontainer(2) + 0.5;
            Phi[1] = Ncontainer(1) -  3.0/4.0 * Ncontainer(2) + 0.5;
            Phi[2] = 5.0/2.0 * Ncontainer(2) - 1.0;
        }
    }
    else if (dimension == 3)
    {
        if (number_of_slave_nodes == 3) // Triangle
        {
            Phi[0] = 3.0 - 4.0 * xi_local - 4.0 * eta_local;
            Phi[1] = 4.0 * xi_local  - 1.0;
            Phi[2] = 4.0 * eta_local - 1.0;
        }
        else if (number_of_slave_nodes == 4) // Quadrilateral
        {
            array_1d<double,3> aux_coordinates = ZeroVector(3);
            aux_coordinates[0] =  xi_local;
            aux_coordinates[1] = eta_local;
            Vector Ncontainer = ZeroVector(3);
            Ncontainer = GetGeometry().ShapeFunctionsValues(Ncontainer, aux_coordinates);

            Phi[0] =   4.0 * Ncontainer(0) - 2.0 * Ncontainer(1) + 1.0 * Ncontainer(2) - 2.0 * Ncontainer(3);
            Phi[1] = - 2.0 * Ncontainer(0) + 4.0 * Ncontainer(1) - 2.0 * Ncontainer(2) + 1.0 * Ncontainer(3);
            Phi[2] =   1.0 * Ncontainer(0) - 2.0 * Ncontainer(1) + 4.0 * Ncontainer(2) - 2.0 * Ncontainer(3);
            Phi[3] = - 2.0 * Ncontainer(0) + 1.0 * Ncontainer(1) - 2.0 * Ncontainer(2) + 4.0 * Ncontainer(3);
        }
        else
        {
            KRATOS_THROW_ERROR( std::logic_error,  " CONDITION can not supply the required number_of_slave_nodes: ", number_of_slave_nodes );
        }
    }
        
    return Phi;
}

/***********************************************************************************/
/***********************************************************************************/

const Matrix MortarContactCondition::LagrangeMultiplierShapeFunctionLocalGradient( 
    const double xi_local, 
    const double eta_local 
    )
{
    // Working space dimension
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    
    // Number of slave nodes
    const unsigned int number_of_slave_nodes = GetGeometry( ).PointsNumber( );
    const unsigned int local_dimension_slave = GetGeometry( ).LocalSpaceDimension( );
    
    Matrix DPhi_De = ZeroMatrix( number_of_slave_nodes, local_dimension_slave );
            
    if (dimension == 2)
    {
        if (number_of_slave_nodes == 2) // First order
        {
            DPhi_De( 0, 0 ) = - 3.0 / 2.0;
            DPhi_De( 1, 0 ) = + 3.0 / 2.0;
        }
        else if (number_of_slave_nodes == 3) // Second order
        {
            array_1d<double,3> aux_coordinates = ZeroVector(3);
            aux_coordinates[0] = xi_local;
            Matrix DNcontainer = ZeroMatrix(number_of_slave_nodes, local_dimension_slave);
            DNcontainer = GetGeometry().ShapeFunctionsLocalGradients( DNcontainer , aux_coordinates );
            
            DPhi_De( 0, 0 ) = DNcontainer(0, 0) -  3.0/4.0 * DNcontainer(2, 0);
            DPhi_De( 1, 0 ) = DNcontainer(1, 0) -  3.0/4.0 * DNcontainer(2, 0);
            DPhi_De( 2, 0 ) = 5.0/2.0 * DNcontainer(2, 0);
        }
    }
    else if (dimension == 3)
    {
        if (number_of_slave_nodes == 3) // Triangle
        {
            DPhi_De( 0, 0 ) = - 4.0;
            DPhi_De( 1, 0 ) =   4.0;
            DPhi_De( 2, 0 ) =   4.0;
            
            DPhi_De( 0, 1 ) = - 4.0;
            DPhi_De( 1, 1 ) =   0.0;
            DPhi_De( 2, 1 ) =   4.0;
        }
        else if (number_of_slave_nodes == 4) // Quadrilateral
        {
            array_1d<double,3> aux_coordinates = ZeroVector(3);
            aux_coordinates[0] =  xi_local;
            aux_coordinates[1] = eta_local;
            Matrix DNcontainer = ZeroMatrix(number_of_slave_nodes, local_dimension_slave);
            DNcontainer = GetGeometry().ShapeFunctionsLocalGradients(DNcontainer, aux_coordinates);

            DPhi_De( 0, 0 ) =   4.0 * DNcontainer( 0, 0 ) - 2.0 * DNcontainer( 1, 0 ) + 1.0 * DNcontainer( 2, 0 ) - 2.0 * DNcontainer( 3, 0 );
            DPhi_De( 1, 0 ) = - 2.0 * DNcontainer( 0, 0 ) + 4.0 * DNcontainer( 1, 0 ) - 2.0 * DNcontainer( 2, 0 ) + 1.0 * DNcontainer( 3, 0 );
            DPhi_De( 2, 0 ) =   1.0 * DNcontainer( 0, 0 ) - 2.0 * DNcontainer( 1, 0 ) + 4.0 * DNcontainer( 2, 0 ) - 2.0 * DNcontainer( 3, 0 );
            DPhi_De( 3, 0 ) = - 2.0 * DNcontainer( 0, 0 ) + 1.0 * DNcontainer( 1, 0 ) - 2.0 * DNcontainer( 2, 0 ) + 4.0 * DNcontainer( 3, 0 );
            DPhi_De( 0, 1 ) =   4.0 * DNcontainer( 0, 1 ) - 2.0 * DNcontainer( 1, 1 ) + 1.0 * DNcontainer( 2, 1 ) - 2.0 * DNcontainer( 3, 1 );
            DPhi_De( 1, 1 ) = - 2.0 * DNcontainer( 0, 1 ) + 4.0 * DNcontainer( 1, 1 ) - 2.0 * DNcontainer( 2, 1 ) + 1.0 * DNcontainer( 3, 1 );
            DPhi_De( 2, 1 ) =   1.0 * DNcontainer( 0, 1 ) - 2.0 * DNcontainer( 1, 1 ) + 4.0 * DNcontainer( 2, 1 ) - 2.0 * DNcontainer( 3, 1 );
            DPhi_De( 3, 1 ) = - 2.0 * DNcontainer( 0, 1 ) + 1.0 * DNcontainer( 1, 1 ) - 2.0 * DNcontainer( 2, 1 ) + 4.0 * DNcontainer( 3, 1 );
        }
        else
        {
            KRATOS_THROW_ERROR( std::logic_error,  " CONDITION can not supply the required number_of_slave_nodes: ", number_of_slave_nodes );
        }
    }

    return DPhi_De;
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContactCondition::CalculateLocalSystem( 
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
    LocalSystem.CalculationFlags.Set(MortarContactCondition::COMPUTE_LHS_MATRIX_WITH_COMPONENTS, true);
    LocalSystem.CalculationFlags.Set(MortarContactCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS, true);

    //Initialize sizes for the system components:
    if ( rLHSVariables.size( ) != rLeftHandSideMatrices.size( ) )
    {
        rLeftHandSideMatrices.resize( rLHSVariables.size( ) );
    }

    if ( rRHSVariables.size( ) != rRightHandSideVectors.size( ) )
    {
        rRightHandSideVectors.resize( rRHSVariables.size( ) );
    }

    LocalSystem.CalculationFlags.Set(MortarContactCondition::COMPUTE_LHS_MATRIX, true);
    for ( unsigned int i = 0; i < rLeftHandSideMatrices.size( ); i++ )
    {
        // Note: rRightHandSideVectors.size() > 0
        this->InitializeSystemMatrices( rLeftHandSideMatrices[i], rRightHandSideVectors[0],LocalSystem.CalculationFlags );
    }

    LocalSystem.CalculationFlags.Set( MortarContactCondition::COMPUTE_RHS_VECTOR, true );
    LocalSystem.CalculationFlags.Set( MortarContactCondition::COMPUTE_LHS_MATRIX, false ); // Temporarily only
    for ( unsigned int i = 0; i < rRightHandSideVectors.size( ); i++ )
    {
        // Note: rLeftHandSideMatrices.size() > 0
        this->InitializeSystemMatrices( rLeftHandSideMatrices[0], rRightHandSideVectors[i], LocalSystem.CalculationFlags  );
    }
    LocalSystem.CalculationFlags.Set( MortarContactCondition::COMPUTE_LHS_MATRIX, true ); // Reactivated again

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

void MortarContactCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;
    
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContactCondition::COMPUTE_LHS_MATRIX, true );
    LocalSystem.CalculationFlags.Set( MortarContactCondition::COMPUTE_RHS_VECTOR, true );

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

void MortarContactCondition::CalculateLeftHandSide( 
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContactCondition::COMPUTE_LHS_MATRIX, true );

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

void MortarContactCondition::CalculateLeftHandSide( 
    std::vector< MatrixType >& rLeftHandSideMatrices,
    const std::vector< Variable< MatrixType > >& rLHSVariables,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContactCondition::COMPUTE_LHS_MATRIX, true );
    LocalSystem.CalculationFlags.Set( MortarContactCondition::COMPUTE_LHS_MATRIX_WITH_COMPONENTS, true );

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

void MortarContactCondition::CalculateRightHandSide( 
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContactCondition::COMPUTE_RHS_VECTOR, true);

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

void MortarContactCondition::CalculateRightHandSide( 
    std::vector< VectorType >& rRightHandSideVectors,
    const std::vector< Variable< VectorType > >& rRHSVariables,
    ProcessInfo& rCurrentProcessInfo )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContactCondition::COMPUTE_RHS_VECTOR, true );
    LocalSystem.CalculationFlags.Set( MortarContactCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS, true );

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

void MortarContactCondition::InitializeSystemMatrices( 
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    Flags& rCalculationFlags 
    )
{
    const unsigned int condition_size = this->CalculateConditionSize( );
    
    // Resizing as needed the LHS
    if ( rCalculationFlags.Is( MortarContactCondition::COMPUTE_LHS_MATRIX ) ) // Calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != condition_size )
        {
            rLeftHandSideMatrix.resize( condition_size, condition_size, false );
        }
        noalias( rLeftHandSideMatrix ) = ZeroMatrix( condition_size, condition_size ); // Resetting LHS
    }

    // Resizing as needed the RHS
    if ( rCalculationFlags.Is( MortarContactCondition::COMPUTE_RHS_VECTOR ) ) // Calculation of the matrix is required
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

void MortarContactCondition::CalculateMassMatrix( 
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

void MortarContactCondition::CalculateDampingMatrix( 
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

const unsigned int MortarContactCondition::CalculateConditionSize( )
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    unsigned int condition_size = 0;
        
    for ( unsigned int i_master_elem = 0; i_master_elem < mThisMasterElements.size( ); ++i_master_elem )
    {
        condition_size += mThisMasterElements[i_master_elem]->GetGeometry( ).PointsNumber( );
        condition_size += 2 * GetGeometry( ).PointsNumber( );
    }
  
    condition_size *= dimension;
    
    return condition_size;
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContactCondition::CalculateConditionSystem( 
    LocalSystemComponents& rLocalSystem,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
    
    // Create and initialize condition variables:
    GeneralVariables Variables;
    
    // Initialize the current contact data
    ContactData rContactData;

    // Slave segment info
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_elements_master = mThisMasterElements.size( );
    const unsigned int number_of_nodes_slave     = GetGeometry().PointsNumber();
    
    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = mUseManualColocationIntegration ?
                                                                         mColocationIntegration.IntegrationPoints( ) :
                                                                         GetGeometry( ).IntegrationPoints( mThisIntegrationMethod );
    this->InitializeContactData(rContactData, rCurrentProcessInfo);
    
    for (unsigned int PairIndex = 0; PairIndex < number_of_elements_master; ++PairIndex)
    {   
        // Initialize general variables for the current master element
        this->InitializeGeneralVariables( Variables, rCurrentProcessInfo, PairIndex );
        
        // Update the contact data
        this->UpdateContactData(rContactData, PairIndex);
         
        // Master segment info
        const GeometryType& current_master_element = Variables.GetMasterElement( );
        
        Vector aux_int_gap      = ZeroVector(number_of_nodes_slave);
        Vector aux_int_slip     = ZeroVector(number_of_nodes_slave);
        Vector aux_int_friction = ZeroVector(number_of_nodes_slave);
        
        const unsigned int pair_size = dimension * ( current_master_element.PointsNumber( ) + 2 * GetGeometry( ).PointsNumber( ) ); 
        MatrixType LHS_contact_pair = ZeroMatrix(pair_size, pair_size);
        VectorType RHS_contact_pair = ZeroVector(pair_size);
        
        double total_weight = 0.0;
        
        // Integrating the LHS and RHS
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            // Calculate the kinematic variables
            const bool inside = this->CalculateKinematics( Variables, PointNumber, PairIndex, integration_points );
        
//             Variables.print();
            
            // Calculate the gap in the integration node and check tolerance
            const double integration_point_gap = inner_prod(rContactData.Gaps, Variables.N_Slave);
            
            if (inside == true)
            {
                const double IntegrationWeight = integration_points[PointNumber].Weight();
                total_weight += IntegrationWeight;
               
                // The normal LM augmented
                const double augmented_normal_lm = this->AugmentedNormalLM(Variables, rContactData, integration_point_gap, dimension);
                
                // The slip of th GP
                double integration_point_slip;
                
                // The tangent LM augmented
                const double augmented_tangent_lm = this->AugmentedTangentLM(Variables, rContactData, current_master_element, integration_point_slip, dimension);
                
                // Calculation of the matrix is required
                if ( rLocalSystem.CalculationFlags.Is( MortarContactCondition::COMPUTE_LHS_MATRIX ) ||
                        rLocalSystem.CalculationFlags.Is( MortarContactCondition::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
                {
                    this->CalculateLocalLHS( LHS_contact_pair, Variables, rContactData, IntegrationWeight, augmented_normal_lm, augmented_tangent_lm, integration_point_gap,  integration_point_slip );
                }

                // Calculation of the vector is required
                if ( rLocalSystem.CalculationFlags.Is( MortarContactCondition::COMPUTE_RHS_VECTOR ) ||
                        rLocalSystem.CalculationFlags.Is( MortarContactCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
                {
                    this->CalculateLocalRHS( RHS_contact_pair, Variables, rContactData, IntegrationWeight, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip );
                }
                
                for (unsigned int iNode = 0; iNode < number_of_nodes_slave; iNode++)
                {
                    aux_int_gap[iNode]      += IntegrationWeight * integration_point_gap  * Variables.DetJSlave * Variables.Phi_LagrangeMultipliers[iNode];
                    aux_int_slip[iNode]     += IntegrationWeight * integration_point_slip * Variables.DetJSlave * Variables.Phi_LagrangeMultipliers[iNode];
                    aux_int_friction[iNode] += IntegrationWeight * Variables.mu           * Variables.DetJSlave * Variables.Phi_LagrangeMultipliers[iNode];
                }
            }
        }
        
        // We can consider the pair if at least one of the collocation point is inside 
        if (total_weight > 0.0)
        {
            // Assemble of the matrix is required
            if ( rLocalSystem.CalculationFlags.Is( MortarContactCondition::COMPUTE_LHS_MATRIX ) ||
                    rLocalSystem.CalculationFlags.Is( MortarContactCondition::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
            {
                // Contributions to stiffness matrix calculated on the reference config
                this->CalculateAndAddLHS( rLocalSystem, LHS_contact_pair, PairIndex, current_master_element );
            }

            // Assemble of the vector is required
            if ( rLocalSystem.CalculationFlags.Is( MortarContactCondition::COMPUTE_RHS_VECTOR ) ||
                    rLocalSystem.CalculationFlags.Is( MortarContactCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
            {
                // Contribution to previous step contact force and residuals vector
                this->CalculateAndAddRHS( rLocalSystem, RHS_contact_pair, PairIndex, current_master_element );
            }
            
//             std::cout << "--------------------------------------------------" << std::endl;
//             KRATOS_WATCH(this->Id());
//             KRATOS_WATCH(PairIndex);
// //             KRATOS_WATCH(LHS_contact_pair);
//             LOG_MATRIX_PRETTY( LHS_contact_pair );
// //             KRATOS_WATCH(RHS_contact_pair);
//             LOG_VECTOR_PRETTY( RHS_contact_pair );
            
            for (unsigned int iNode = 0; iNode < number_of_nodes_slave; iNode++)
            {
                GetGeometry()[iNode].GetValue(WEIGHTED_GAP)      += aux_int_gap[iNode]; 
                GetGeometry()[iNode].GetValue(WEIGHTED_SLIP)     += aux_int_slip[iNode]; 
                GetGeometry()[iNode].GetValue(WEIGHTED_FRICTION) += aux_int_friction[iNode]; 
            }
        }
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContactCondition::InitializeGeneralVariables(
    GeneralVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo,
    const unsigned int& rMasterElementIndex
    )
{
    // Working space dimension
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Resize according to the integration order
    const unsigned int& integration_point_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );
    
    // Slave segment info
    GeometryType& current_slave_element = this->GetGeometry();
    const unsigned int number_of_slave_nodes = current_slave_element.PointsNumber();
    const unsigned int local_dimension_slave = current_slave_element.LocalSpaceDimension();

    // Master segment info
    GeometryType& current_master_element = mThisMasterElements[rMasterElementIndex]->GetGeometry();
    const unsigned int number_nodes_master = current_master_element.PointsNumber();
    const unsigned int local_dimension_master = current_master_element.LocalSpaceDimension();
    
    // Slave element info
    rVariables.Initialize(local_dimension_master, local_dimension_slave, number_nodes_master, number_of_slave_nodes, dimension, integration_point_number );

    rVariables.SetMasterElement( current_master_element );
    rVariables.SetMasterElementIndex( rMasterElementIndex );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContactCondition::InitializeContactData(
    ContactData& rContactData,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // Working space dimension
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Slave segment info
    GeometryType& current_slave_element = this->GetGeometry();
    const unsigned int number_of_slave_nodes = current_slave_element.PointsNumber();
    
    // Slave element info
    rContactData.Initialize(GetGeometry(), number_of_slave_nodes, dimension );
    
    /* Set Delta time */
    rContactData.Dt = rCurrentProcessInfo[DELTA_TIME];
    
    /* LM */
    rContactData.LagrangeMultipliers = ContactUtilities::GetVariableMatrix(GetGeometry(), VECTOR_LAGRANGE_MULTIPLIER, 0); 
    
    /* NORMALS AND TANGENTS */ 
    rContactData.NormalsSlave = ContactUtilities::GetVariableMatrix(GetGeometry(),  NORMAL); 
    rContactData.Tangent1Slave = ContactUtilities::GetVariableMatrix(GetGeometry(), TANGENT_XI); 
    if (dimension == 3)
    {
        rContactData.Tangent2Slave = ContactUtilities::GetVariableMatrix(GetGeometry(), TANGENT_ETA); 
    }
    
    if (GetProperties().Has(NORMAL_AUGMENTATION_FACTOR) == true)
    {
        rContactData.epsilon_normal  = GetProperties().GetValue(NORMAL_AUGMENTATION_FACTOR);
    }
    if (GetProperties().Has(TANGENT_AUGMENTATION_FACTOR) == true)
    {
        rContactData.epsilon_tangent = GetProperties().GetValue(TANGENT_AUGMENTATION_FACTOR);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContactCondition::UpdateContactData(
    ContactData& rContactData,
    const unsigned int& rMasterElementIndex
    )
{
    // Working space dimension
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Slave segment info
    GeometryType& current_slave_element = this->GetGeometry();
    const unsigned int number_of_slave_nodes = current_slave_element.PointsNumber();
    
    // Master segment info
    GeometryType& current_master_element = mThisMasterElements[rMasterElementIndex]->GetGeometry();
    const unsigned int number_nodes_master = current_master_element.PointsNumber();
    
    // Slave element info
    rContactData.UpdateMasterPair(current_master_element, number_nodes_master, dimension );
    
    /* NORMALS AND GAPS */
    for (unsigned int iNode = 0; iNode < number_of_slave_nodes; iNode++)
    {
        array_1d<double,3> normal = GetGeometry()[iNode].GetValue(NORMAL);
        
        PointType projected_global;
        ContactUtilities::ProjectDirection(current_master_element, GetGeometry()[iNode], projected_global, rContactData.Gaps(iNode), normal ); // NOTE: This is not the CPP, so the solution can be wrong
    }
}

/*********************************COMPUTE KINEMATICS*********************************/
/************************************************************************************/

bool MortarContactCondition::CalculateKinematics( 
    GeneralVariables& rVariables,
    const double& rPointNumber,
    const unsigned int& rPairIndex,
    const GeometryType::IntegrationPointsArrayType& integration_points
    )
{
    /* DEFINITIONS */
    GeometryType& slave_nodes  = GetGeometry( );
    GeometryType& master_nodes = rVariables.GetMasterElement( );
    const unsigned int number_of_master_nodes = master_nodes.PointsNumber( );
    const unsigned int number_of_slave_nodes  =  slave_nodes.PointsNumber( );

    /* LOCAL COORDINATES */
    const PointType& local_point = integration_points[rPointNumber].Coordinates();
    
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
//     rVariables.Phi_LagrangeMultipliers = rVariables.N_Slave; // TODO: This could be needed in the future to be different than the standart shape functions 
    rVariables.Phi_LagrangeMultipliers = LagrangeMultiplierShapeFunctionValue( local_point.Coordinate(1), local_point.Coordinate(2) );
    
    // SHAPE FUNCTION DERIVATIVES
    rVariables.DN_De_Slave  =  slave_nodes.ShapeFunctionsLocalGradients( rVariables.DN_De_Slave , local_point );
//     rVariables.DPhi_De_LagrangeMultipliers = slave_nodes.ShapeFunctionsLocalGradients( rVariables.DN_De_Slave , local_point );// TODO: This could be needed in the future to be different than the standart shape functions
   rVariables.DPhi_De_LagrangeMultipliers = LagrangeMultiplierShapeFunctionLocalGradient( local_point.Coordinate(1), local_point.Coordinate(2) );
    
    // MASTER CONDITION
    const bool inside = this->MasterShapeFunctionValue( rVariables, local_point);
    
    /* CALCULATE JACOBIAN AND JACOBIAN DETERMINANT */
    slave_nodes.Jacobian( rVariables.j_Slave, local_point.Coordinates() );
    rVariables.DetJSlave = ContactUtilities::ContactElementDetJacobian( rVariables.j_Slave );
    
    /* FRICTION COEFFICIENT */
    if (GetProperties().Has(FRICTION_COEFFICIENT) == true) // NOTE: In function of the friction law
    {
        rVariables.mu = GetProperties().GetValue(FRICTION_COEFFICIENT);
    }        
    
    return inside;
}
 
/***********************************************************************************/
/*************** METHODS TO CALCULATE THE CONTACT CONDITION MATRICES ***************/
/***********************************************************************************/

bool MortarContactCondition::MasterShapeFunctionValue(
    GeneralVariables& rVariables,
    const PointType& local_point 
    )
{    
    GeometryType& master_seg = rVariables.GetMasterElement( );
    rVariables.N_Master     = ZeroVector( master_seg.PointsNumber( ) );
    rVariables.DN_De_Master = ZeroMatrix( master_seg.PointsNumber( ), master_seg.LocalSpaceDimension( ) );

    PointType projected_gp_global;
    const VectorType normal = ContactUtilities::GaussPointNormal(rVariables.N_Slave, GetGeometry());
    
    GeometryType::CoordinatesArrayType slave_gp_global;
    double aux_dist = 0.0;
    this->GetGeometry( ).GlobalCoordinates( slave_gp_global, local_point );
    ContactUtilities::ProjectDirection( master_seg, slave_gp_global, projected_gp_global, aux_dist, -normal ); // The opposite direction
    
    GeometryType::CoordinatesArrayType projected_gp_local;
    
    const bool inside = master_seg.IsInside( projected_gp_global.Coordinates( ), projected_gp_local ) ;
    
    if( inside == true )
    {
        // SHAPE FUNCTIONS 
        rVariables.N_Master     = master_seg.ShapeFunctionsValues(         rVariables.N_Master,     projected_gp_local );         
        rVariables.DN_De_Master = master_seg.ShapeFunctionsLocalGradients( rVariables.DN_De_Master, projected_gp_local );
    }
    
    return inside;
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContactCondition::CalculateAndAddLHS(
    LocalSystemComponents& rLocalSystem,
    MatrixType& LHS_contact_pair, 
    const unsigned int rPairIndex,
    const GeometryType& current_master_element
    )
{
    /* DEFINITIONS */
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_slave_nodes = GetGeometry( ).PointsNumber( );
    const unsigned int number_of_master_nodes = current_master_element.PointsNumber( );
    const unsigned int number_of_total_nodes = number_of_slave_nodes + number_of_master_nodes;
        
    unsigned int index_1 = 0;
    for (unsigned int i_slave = 0; i_slave < number_of_slave_nodes; ++i_slave )
    {
//         // We impose a 0 zero LM in the inactive nodes
//         if (GetGeometry( )[i_slave].Is(ACTIVE) == false)
//         {
//             subrange(LHS_contact_pair, index_1, index_1 + dimension, index_1, index_1 + dimension) = IdentityMatrix(dimension, dimension);
//         }
//         // And the tangent direction (sliding)
//         else
//         {            
            Matrix tangent_matrix;
            tangent_matrix.resize(dimension - 1, dimension);

            const array_1d<double, 3> tangent1 = GetGeometry()[i_slave].GetValue(TANGENT_XI);
            const array_1d<double, 3> tangent2 = GetGeometry()[i_slave].GetValue(TANGENT_ETA);
            for (unsigned int i = 0; i < dimension; i++)
            {
                tangent_matrix(0, i) = tangent1[i];
                if (dimension == 3)
                {
                    tangent_matrix(1, i) = tangent2[i];
                }
            }
            
            index_1 = (number_of_total_nodes + i_slave) * dimension;
            subrange(LHS_contact_pair, index_1 + 1, index_1 + dimension, index_1 , index_1 + dimension) = tangent_matrix;
        }
//     }
    
    if ( rLocalSystem.CalculationFlags.Is( MortarContactCondition::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
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
                
                // Assemble in the correct position
                this->AssembleContactPairLHSToConditionSystem(LHS_contact_pair, rLeftHandSideMatrix, rPairIndex);
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
        
        // Assemble in the correct position
        this->AssembleContactPairLHSToConditionSystem(LHS_contact_pair, rLeftHandSideMatrix, rPairIndex);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContactCondition::AssembleContactPairLHSToConditionSystem(
    MatrixType& rPairLHS,
    MatrixType& rConditionLHS,
    const unsigned int rPairIndex
    )
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_slave_nodes = GetGeometry( ).PointsNumber( );
    const unsigned int number_of_master_nodes = mThisMasterElements[rPairIndex]->GetGeometry( ).PointsNumber( );
    const unsigned int current_pair_size = dimension * ( number_of_master_nodes + 2 * number_of_slave_nodes );
  
    // Find location of the piar's master DOFs in ConditionLHS
    unsigned int index_begin = 0;

    for ( unsigned int i_master_elem = 0; i_master_elem < rPairIndex ; ++i_master_elem ) 
    {
        index_begin += mThisMasterElements[i_master_elem]->GetGeometry( ).PointsNumber( );
        index_begin += 2 * number_of_slave_nodes;
    }

    index_begin *= dimension;
  
    const unsigned int index_end = index_begin + current_pair_size;
    
    subrange( rConditionLHS, index_begin, index_end, index_begin, index_end) += rPairLHS;
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContactCondition::CalculateLocalLHS(
    MatrixType& rPairLHS,
    GeneralVariables& rVariables,
    const ContactData& rContactData,
    const double& rIntegrationWeight,
    const double& augmented_normal_lm,
    const double& augmented_tangent_lm,
    const double& integration_point_gap,
    const double& integration_point_slip
    )
{
    /* DEFINITIONS */
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    GeometryType& slave_nodes  = GetGeometry( );
    GeometryType& master_nodes = rVariables.GetMasterElement( );
    const unsigned int number_of_master_nodes = master_nodes.PointsNumber( );
    const unsigned int number_of_slave_nodes  =  slave_nodes.PointsNumber( );
    
    const Vector N1           = rVariables.N_Slave;
    const Vector N2           = rVariables.N_Master;
    const Vector Phi          = rVariables.Phi_LagrangeMultipliers;
//     const Matrix DPhi         = rVariables.DPhi_De_LagrangeMultipliers;
    const double detJ         = rVariables.DetJSlave; 
    
    
    if (dimension == 2)
    {
        if (number_of_master_nodes == 2 &&  number_of_slave_nodes == 2)
        {
            if (augmented_normal_lm <= 0.0)
            {
                // Contact active
                rPairLHS += rIntegrationWeight * Contact2D2N2N::ComputeGaussPointActiveLHS( N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
                
                if (std::abs(augmented_tangent_lm) - rVariables.mu * augmented_normal_lm >= 0.0)
                {
                    // Slip
//                     rPairLHS += rIntegrationWeight * Contact2D2N2N::ComputeGaussPointStickLHS( N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
                    rPairLHS += rIntegrationWeight * Contact2D2N2N::ComputeGaussPointSlipLHS( N1, N2, Phi, detJ, rVariables.mu, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
                }
                else
                {
                    // Stick
                    rPairLHS += rIntegrationWeight * Contact2D2N2N::ComputeGaussPointStickLHS( N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
                }
            }
//             else
//             {        
//                 // Contact inactive
//                 rPairLHS += rIntegrationWeight * Contact2D2N2N::ComputeGaussPointInactiveLHS( N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//             }
        }
        else
        {
            KRATOS_WATCH(number_of_master_nodes);
            KRATOS_THROW_ERROR( std::logic_error,  " COMBINATION OF GEOMETRIES NOT IMPLEMENTED. Number of slave elements: ", number_of_slave_nodes );
        }
    }
    else
    {
        if (number_of_master_nodes == 3 &&  number_of_slave_nodes == 3)
        {
//             if (augmented_normal_lm <= 0.0)
//             {
//                 // Contact active
//                 rPairLHS += rIntegrationWeight * Contact3D3N3N::ComputeGaussPointActiveLHS( N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//                 
//                 if (std::abs(augmented_tangent_lm) - rVariables.mu * augmented_normal_lm >= 0.0)
//                 {
//                     // Slip
//                     rPairLHS += rIntegrationWeight * Contact3D3N3N::ComputeGaussPointSlipLHS( N1, N2, Phi, detJ, rVariables.mu, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//                 }
//                 else
//                 {
//                     // Stick
//                     rPairLHS += rIntegrationWeight * Contact3D3N3N::ComputeGaussPointStickLHS( N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//                 }
//             }
//             else
//             {        
//                 // Contact inactive
//                 rPairLHS += rIntegrationWeight * Contact3D3N3N::ComputeGaussPointInactiveLHS( N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//             }                       
        }
        else if (number_of_master_nodes == 4 &&  number_of_slave_nodes == 4)
        {
            const Matrix DN1 = rVariables.DN_De_Slave;
            
//             if (augmented_normal_lm <= 0.0)
//             {
//                 // Contact active
//                 rPairLHS += rIntegrationWeight * Contact3D4N4N::ComputeGaussPointActiveLHS( N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//                 
//                 if (std::abs(augmented_tangent_lm) - rVariables.mu * augmented_normal_lm >= 0.0)
//                 {
//                     // Slip
//                     rPairLHS += rIntegrationWeight * Contact3D4N4N::ComputeGaussPointSlipLHS( N1, N2, Phi, detJ, rVariables.mu, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//                 }
//                 else
//                 {
//                     // Stick
//                     rPairLHS += rIntegrationWeight * Contact3D4N4N::ComputeGaussPointStickLHS( N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//                 }
//             }
//             else
//             {        
//                 // Contact inactive
//                 rPairLHS += rIntegrationWeight * Contact3D4N4N::ComputeGaussPointInactiveLHS( N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//             }      
        }
        else
        {
            KRATOS_WATCH(number_of_master_nodes);
            KRATOS_THROW_ERROR( std::logic_error,  " COMBINATION OF GEOMETRIES NOT IMPLEMENTED. Number of slave elements: ", number_of_slave_nodes );
        }
    }
    
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContactCondition::CalculateAndAddRHS(
    LocalSystemComponents& rLocalSystem,
    VectorType& RHS_contact_pair,
    const unsigned int rPairIndex,
    const GeometryType& current_master_element
    )
{   
    
//     /* DEFINITIONS */
//     const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
//     const unsigned int number_of_slave_nodes = GetGeometry( ).PointsNumber( );
//     const unsigned int number_of_master_nodes = current_master_element.PointsNumber( );
//     const unsigned int number_of_total_nodes = number_of_slave_nodes + number_of_master_nodes;
//         
//     unsigned int index = 0;
//     
//     for (unsigned int i_slave = 0; i_slave < number_of_slave_nodes; ++i_slave )
//     {
//         index = (number_of_total_nodes + i_slave) * dimension; 
// 
//         if (GetGeometry( )[i_slave].Is(ACTIVE) == false)
//         {
//             subrange(RHS_contact_pair, index, index + dimension) = ZeroVector(dimension);
//         }
//     }
        
    if ( rLocalSystem.CalculationFlags.Is( MortarContactCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
    {
        /* COMPONENT-WISE RHS VECTOR */
        const std::vector<Variable<VectorType> >& rRightHandSideVariables = rLocalSystem.GetRightHandSideVariables( );
        bool calculated;

        for ( unsigned int i = 0; i < rRightHandSideVariables.size( ); i++ )
        {
            calculated = false;

            if ( rRightHandSideVariables[i] == MORTAR_CONTACT_OPERATOR )
            {
                VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVectors()[i];

                // Assemble
                this->AssembleContactPairRHSToConditionSystem( RHS_contact_pair, rRightHandSideVector, rPairIndex );
                
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
        /* SINGLE RHS VECTOR */
        VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();
        
        // Assemble
        this->AssembleContactPairRHSToConditionSystem( RHS_contact_pair, rRightHandSideVector, rPairIndex );
    }
    
}
  
/***********************************************************************************/
/***********************************************************************************/

void MortarContactCondition::AssembleContactPairRHSToConditionSystem(
    VectorType& rPairRHS,
    VectorType& rConditionRHS,
    const unsigned int rPairIndex
    )
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_slave_nodes = GetGeometry( ).PointsNumber( );
    const unsigned int number_of_master_nodes = mThisMasterElements[rPairIndex]->GetGeometry( ).PointsNumber( );
    const unsigned int current_pair_size = dimension * ( number_of_master_nodes + 2 * number_of_slave_nodes );
  
    // Find location of the pair's master DOFs in ConditionRHS
    unsigned int index_begin = 0;

    for ( unsigned int i_master_elem = 0; i_master_elem < rPairIndex; ++i_master_elem )
    {
        index_begin += mThisMasterElements[i_master_elem]->GetGeometry( ).PointsNumber( );
        index_begin += 2 * number_of_slave_nodes;
    }

    index_begin *= dimension;
    const unsigned int index_end = index_begin + current_pair_size;
    
    // Computing subrange
    subrange( rConditionRHS, index_begin, index_end ) += rPairRHS;
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContactCondition::CalculateLocalRHS(
    VectorType& rPairRHS,
    GeneralVariables& rVariables,
    const ContactData& rContactData,
    const double& rIntegrationWeight,
    const double& augmented_normal_lm,
    const double& augmented_tangent_lm,
    const double& integration_point_gap,
    const double& integration_point_slip
    )
{
    /* DEFINITIONS */
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    GeometryType& slave_nodes  = GetGeometry( );
    GeometryType& master_nodes = rVariables.GetMasterElement( );
    const unsigned int number_of_master_nodes = master_nodes.PointsNumber( );
    const unsigned int number_of_slave_nodes  =  slave_nodes.PointsNumber( );
    
    const Vector N1           = rVariables.N_Slave;
    const Vector N2           = rVariables.N_Master;
    const Vector Phi          = rVariables.Phi_LagrangeMultipliers;
    const double detJ         = rVariables.DetJSlave;
    
    if (dimension == 2)
    {
        if (number_of_master_nodes == 2 &&  number_of_slave_nodes == 2)
        {
            if (augmented_normal_lm <= 0.0) 
            {
                // Contact active
                rPairRHS += rIntegrationWeight * Contact2D2N2N::ComputeGaussPointActiveRHS(N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
                
                if (std::abs(augmented_tangent_lm) - rVariables.mu * augmented_normal_lm >= 0.0)
                {
                    // Slip
//                     rPairRHS += rIntegrationWeight * Contact2D2N2N::ComputeGaussPointStickRHS(N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
                    rPairRHS += rIntegrationWeight * Contact2D2N2N::ComputeGaussPointSlipRHS(N1, N2, Phi, detJ, rVariables.mu, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
                }
                else
                {
                    // Stick
                    rPairRHS += rIntegrationWeight * Contact2D2N2N::ComputeGaussPointStickRHS(N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
                }
            }
//             else
//             {
//                 // Contact inactive
//                 rPairRHS += rIntegrationWeight * Contact2D2N2N::ComputeGaussPointInactiveRHS(N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//             }
        }
        else
        {
            KRATOS_WATCH(number_of_master_nodes);
            KRATOS_THROW_ERROR( std::logic_error,  " COMBINATION OF GEOMETRIES NOT IMPLEMENTED. Number of slave elements: ", number_of_slave_nodes );
        }
    }
    else
    {
        if (number_of_master_nodes == 3 &&  number_of_slave_nodes == 3)
        {
//             if (augmented_normal_lm <= 0.0)
//             {
//                 // Contact active
//                 rPairRHS += rIntegrationWeight *  Contact3D3N3N::ComputeGaussPointRHS(N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//                 if (std::abs(augmented_tangent_lm) - rVariables.mu * augmented_normal_lm >= 0.0)
//                 {
//                     // Slip
//                     rPairRHS += rIntegrationWeight * Contact3D3N3N::ComputeGaussPointSlipRHS(N1, N2, Phi, detJ, rVariables.mu, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//                 }
//                 else
//                 {
//                     // Stick
//                     rPairRHS += rIntegrationWeight * Contact3D3N3N::ComputeGaussPointStickRHS(N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip); 
//                 }
//             }
//             else
//             {
//                 // Contact inactive
//                 rPairRHS += rIntegrationWeight * Contact3D3N3N::ComputeGaussPointInactiveRHS(N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//             }
        }
        if (number_of_master_nodes == 4 &&  number_of_slave_nodes == 4)
        {
            const Matrix DN1 = rVariables.DN_De_Slave;
                        
//             if (augmented_normal_lm <= 0.0)
//             {
//                 // Contact active
//                 rPairRHS += rIntegrationWeight *  Contact3D4N4N::ComputeGaussPointRHS(N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//                 if (std::abs(augmented_tangent_lm) - rVariables.mu * augmented_normal_lm >= 0.0)
//                 {
//                     // Slip
//                     rPairRHS += rIntegrationWeight * Contact3D4N4N::ComputeGaussPointSlipRHS(N1, N2, Phi, detJ, rVariables.mu, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//                 }
//                 else
//                 {
//                     // Stick
//                     rPairRHS += rIntegrationWeight * Contact3D4N4N::ComputeGaussPointStickRHS(N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip); 
//                 }
//             }
//             else
//             {
//                 // Contact inactive
//                 rPairRHS += rIntegrationWeight * Contact3D4N4N::ComputeGaussPointInactiveRHS(N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//             }
        }
        else
        {
            KRATOS_WATCH(number_of_master_nodes);
            KRATOS_THROW_ERROR( std::logic_error,  " COMBINATION OF GEOMETRIES NOT IMPLEMENTED. Number of slave elements: ", number_of_slave_nodes );
        }
    }
    
    
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContactCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY;
        
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    
    rResult.resize( 0, false );

    const std::vector<contact_container> all_conditions = *( this->GetValue( CONTACT_CONTAINERS ) );
        
    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    
    const unsigned int number_of_slave_nodes = GetGeometry().PointsNumber( );
    
    for ( unsigned int i_cond = 0;  i_cond < all_conditions.size( ); ++i_cond )
    {   
        GeometryType& current_master = all_conditions[i_cond].condition->GetGeometry( );
            
        // Master Nodes Displacement Equation IDs
        for ( unsigned int iNode = 0; iNode < current_master.PointsNumber( ); iNode++ )
        {
            NodeType& master_node = current_master[iNode];
            rResult.push_back( master_node.GetDof( DISPLACEMENT_X ).EquationId( ) );
            rResult.push_back( master_node.GetDof( DISPLACEMENT_Y ).EquationId( ) );
            if (dimension == 3)
            {
                rResult.push_back( master_node.GetDof( DISPLACEMENT_Z ).EquationId( ) );
            }
        }
        
        // Slave Nodes Displacement Equation IDs
        for ( unsigned int i_slave = 0; i_slave < number_of_slave_nodes; ++i_slave )
        {
            NodeType& slave_node = GetGeometry()[ i_slave ];
            rResult.push_back( slave_node.GetDof( DISPLACEMENT_X ).EquationId( ) );
            rResult.push_back( slave_node.GetDof( DISPLACEMENT_Y ).EquationId( ) );
            if (dimension == 3)
            {
                rResult.push_back( slave_node.GetDof( DISPLACEMENT_Z ).EquationId( ) );
            }
        }
        
        // Slave Nodes  Lambda Equation IDs
        for ( unsigned int i_slave = 0; i_slave < number_of_slave_nodes; ++i_slave )
        {
            NodeType& slave_node = GetGeometry()[ i_slave ];
            rResult.push_back( slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_X ).EquationId( ) );
            rResult.push_back( slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Y ).EquationId( ) );
            if (dimension == 3)
            {
                rResult.push_back( slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Z ).EquationId( ) );
            }
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContactCondition::GetDofList(
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    
    rConditionalDofList.resize( 0 );

    const std::vector<contact_container> all_conditions = *( this->GetValue( CONTACT_CONTAINERS ) );
    
    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    
    const unsigned int number_of_slave_nodes = GetGeometry().PointsNumber( );
        
    for ( unsigned int i_cond = 0; i_cond < all_conditions.size( ); ++i_cond )
    {
        GeometryType& current_master = all_conditions[i_cond].condition->GetGeometry( );
        
        // Master Nodes Displacement Equation IDs
        for ( unsigned int iNode = 0; iNode < current_master.PointsNumber( ); iNode++ )
        {
            NodeType& master_node = current_master[iNode];
            rConditionalDofList.push_back( master_node.pGetDof( DISPLACEMENT_X ) );
            rConditionalDofList.push_back( master_node.pGetDof( DISPLACEMENT_Y ) );
            if (dimension == 3)
            {
                rConditionalDofList.push_back( master_node.pGetDof( DISPLACEMENT_Z ) );
            }
        }
        
        // Slave Nodes Displacement Equation IDs
        for ( unsigned int i_slave = 0; i_slave < number_of_slave_nodes; ++i_slave )
        {
            NodeType& slave_node = GetGeometry()[ i_slave ];
            rConditionalDofList.push_back( slave_node.pGetDof( DISPLACEMENT_X ) );
            rConditionalDofList.push_back( slave_node.pGetDof( DISPLACEMENT_Y ) );
            if (dimension == 3)
            {
                rConditionalDofList.push_back( slave_node.pGetDof( DISPLACEMENT_Z ) );
            }
        }
        
        // Slave Nodes Lambda Equation IDs
        for ( unsigned int i_slave = 0; i_slave < number_of_slave_nodes; ++i_slave )
        {
            NodeType& slave_node = GetGeometry()[ i_slave ];
            rConditionalDofList.push_back( slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_X ) );
            rConditionalDofList.push_back( slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Y ) );
            if (dimension == 3)
            {
                rConditionalDofList.push_back( slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Z ) );
            }
        }
    }

    KRATOS_CATCH( "" );
}

//******************************* GET DOUBLE VALUE *********************************/
/***********************************************************************************/

void MortarContactCondition::GetValueOnIntegrationPoints( 
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContactCondition::CalculateOnIntegrationPoints( 
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;
    
    // TODO: Add the FRICTION_COEFFICIENT, and maybe if it is ACTIVE or SLIPPING the GP

    // Create and initialize condition variables:
    GeneralVariables Variables;
    
    // Initialize the current contact data
    ContactData rContactData;

    // Slave segment info
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_elements_master = mThisMasterElements.size( );
    const unsigned int number_of_nodes_slave     = GetGeometry().PointsNumber();
    
    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = mUseManualColocationIntegration ?
                                                                         mColocationIntegration.IntegrationPoints( ) :
                                                                         GetGeometry( ).IntegrationPoints( mThisIntegrationMethod );
                                               
    const unsigned int number_of_integration_pts =integration_points.size();
    if ( rOutput.size( ) != number_of_integration_pts )
    {
        rOutput.resize( number_of_integration_pts, false );
    }
    
    const std::vector<double> zero_vector (number_of_integration_pts, 0.0);
    rOutput = zero_vector;
    
    this->InitializeContactData(rContactData, rCurrentProcessInfo);
    
    for (unsigned int PairIndex = 0; PairIndex < number_of_elements_master; ++PairIndex)
    {   
        // Initialize general variables for the current master element
        this->InitializeGeneralVariables( Variables, rCurrentProcessInfo, PairIndex );
        
        // Update the contact data
        this->UpdateContactData(rContactData, PairIndex);
        
        // Master segment info
        const GeometryType& current_master_element = Variables.GetMasterElement( );
        
        // Calculating the values
        for ( unsigned int PointNumber = 0; PointNumber < number_of_integration_pts; PointNumber++ )
        {
            // Calculate the kinematic variables
            const bool inside = this->CalculateKinematics( Variables, PointNumber, PairIndex, integration_points );
            
            // Calculate the gap in the integration node and check tolerance
            const double integration_point_gap = inner_prod(rContactData.Gaps, Variables.N_Slave);

            if (inside == true)
            {
                if (rVariable == NORMAL_CONTACT_STRESS || rVariable == TANGENTIAL_CONTACT_STRESS || rVariable == SLIP_GP)
                {
                    // The normal LM augmented
                    const double augmented_normal_lm = this->AugmentedNormalLM(Variables, rContactData, integration_point_gap, dimension);
                    
                    if ( rVariable == NORMAL_CONTACT_STRESS )
                    {
                        rOutput[PointNumber] += augmented_normal_lm;
                    }
                    else
                    {
                        // The slip of th GP
                        double integration_point_slip;
                    
                        // The tangent LM augmented
                        const double augmented_tangent_lm = this->AugmentedTangentLM(Variables, rContactData, current_master_element, integration_point_slip, dimension);
                        
                        if ( rVariable == TANGENTIAL_CONTACT_STRESS )
                        {
                            rOutput[PointNumber] += augmented_tangent_lm;
                        }
                        else if (rVariable == SLIP_GP )
                        {
                            rOutput[PointNumber] += integration_point_slip; // TODO: Maybe it is necessary to check if it is the smallest of all the gaps 
                        }
                    }
                }
                else if (rVariable == GAP_GP )
                {
                    rOutput[PointNumber] += integration_point_gap; // TODO: Maybe it is necessary to check if it is the smallest of all the gaps 
                }
            }
        }
    } 
    
    KRATOS_CATCH( "" );
}

// TODO: Add the vectorial components, this way we can get the normal, tangent, etc...

/******************* AUXILLIARY METHODS FOR GENERAL CALCULATIONS *******************/
/***********************************************************************************/

double MortarContactCondition::AugmentedNormalLM(
    const GeneralVariables& rVariables,
    const ContactData& rContactData,
    const double& integration_point_gap,
    const unsigned int& dimension
    )
{
    double augmented_normal_lm; 
    
    if (dimension == 2)
    {
        // The LM of the GP 
        const array_1d<double, 2> lm_gp     = prod(trans(rContactData.LagrangeMultipliers), rVariables.Phi_LagrangeMultipliers);
    
        // The normal of the GP
        const array_1d<double, 2> normal_gp = prod(trans(rContactData.NormalsSlave), rVariables.N_Slave);
        
        // The LM in the tangent direction
        augmented_normal_lm  = inner_prod(normal_gp, lm_gp) - rContactData.epsilon_normal * integration_point_gap;
    }
    else // If it is not 2D it is 3D
    {
        // The LM of the GP 
        const array_1d<double, 3> lm_gp     = prod(trans(rContactData.LagrangeMultipliers), rVariables.Phi_LagrangeMultipliers);
    
        // The normal of the GP
        const array_1d<double, 3> normal_gp = prod(trans(rContactData.NormalsSlave), rVariables.N_Slave);
            
        // The LM in the tangent direction
        augmented_normal_lm  = inner_prod(normal_gp, lm_gp) - rContactData.epsilon_normal * integration_point_gap;
    }
                
    return augmented_normal_lm;
}

/***********************************************************************************/
/***********************************************************************************/

double MortarContactCondition::AugmentedTangentLM(
    const GeneralVariables& rVariables,
    const ContactData& rContactData,
    const GeometryType& current_master_element, 
    double& integration_point_slip,
    const unsigned int& dimension
    )
{    
    double augmented_tangent_lm; 
    if (dimension == 2)
    {
        // The LM of the GP 
        const array_1d<double, 2> lm_gp = prod(trans(rContactData.LagrangeMultipliers), rVariables.Phi_LagrangeMultipliers);
        
        // The tangents of the GP
        const array_1d<double, 2> tangent_gp = prod(trans(rContactData.Tangent1Slave), rVariables.N_Slave);
        
        // The tangential LM
        const double tangent_lm  = inner_prod(tangent_gp, lm_gp);
        
        // The velocities matrices
        const Matrix v1 = ContactUtilities::GetVariableMatrix(GetGeometry(),          VELOCITY, 0); 
        const Matrix v2 = ContactUtilities::GetVariableMatrix(current_master_element, VELOCITY, 0);
        
        // The slip of the LM
        const array_1d<double, 2> vector_integration_point_slip = rContactData.Dt * (prod(trans(v1), rVariables.N_Slave) - prod(trans(v2), rVariables.N_Master));
        integration_point_slip = inner_prod(vector_integration_point_slip, tangent_gp);

        // The augmented LM in the tangent direction
        augmented_tangent_lm = tangent_lm + rContactData.epsilon_tangent * integration_point_slip; 
    }
    else
    {
        // The LM of the GP 
        const array_1d<double, 3> lm_gp     = prod(trans(rContactData.LagrangeMultipliers), rVariables.Phi_LagrangeMultipliers);
        
        // The tangents of the GP
        const array_1d<double, 3> tangent_xi_gp  = prod(trans(rContactData.Tangent1Slave), rVariables.N_Slave);
        const array_1d<double, 3> tangent_eta_gp = prod(trans(rContactData.Tangent2Slave), rVariables.N_Slave);
        
        // The tangential LM
        const double tangent_xi_lm  = inner_prod(tangent_xi_gp, lm_gp);
        const double tangent_eta_lm = inner_prod(tangent_eta_gp, lm_gp);
        const double tangent_lm = std::sqrt(tangent_xi_lm * tangent_xi_lm + tangent_eta_lm * tangent_eta_lm); 
        
        // The resultant direction of the LM
        const array_1d<double, 3>& tangent_gp =  (tangent_xi_lm * tangent_xi_gp + tangent_eta_lm * tangent_eta_gp)/tangent_lm; // NOTE: This is the direction of the slip (using the LM as reference)
        
        // The velocities matrices
        const Matrix v1 = ContactUtilities::GetVariableMatrix(GetGeometry(),          VELOCITY, 0); 
        const Matrix v2 = ContactUtilities::GetVariableMatrix(current_master_element, VELOCITY, 0);
        
        // The slip of the LM
        const array_1d<double, 3> vector_integration_point_slip = rContactData.Dt * (prod(trans(v1), rVariables.N_Slave) - prod(trans(v2), rVariables.N_Master));
        integration_point_slip = inner_prod(vector_integration_point_slip, tangent_gp);

        // The augmented LM in the tangent direction
        augmented_tangent_lm = tangent_lm + rContactData.epsilon_tangent * integration_point_slip; 
    }
    
    return augmented_tangent_lm;
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContactCondition::InitializeIntegrationMethod( 
    const unsigned int dimension,
    const unsigned int local_dimension_slave,
    const unsigned int number_of_slave_nodes
    )
{
    mUseManualColocationIntegration = false;
    if( GetProperties().Has(INTEGRATION_ORDER_CONTACT) )
    {
        const double integration_order = GetProperties().GetValue(INTEGRATION_ORDER_CONTACT);
        if (dimension == 2)
        {
            if (integration_order == 3)
            {
                mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_1;
            }
            else if (integration_order == 5)
            {
                mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_2;
            }
            else if (integration_order == 7)
            {
                mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_3;
            }
            else if (integration_order == 9)
            {
                mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_4;
            }
            else if (integration_order == 11)
            {
                mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_5;
            }
            else
            {
                mUseManualColocationIntegration = true;
                mColocationIntegration.Initialize( integration_order,local_dimension_slave, number_of_slave_nodes );
            }
        }
        else
        {
            if (number_of_slave_nodes == 3) // TODO: Complete
            {
//                 if (integration_order == 3)
//                 {
//                 }
//                 else
//                 {
                    mUseManualColocationIntegration = true;
                    mColocationIntegration.Initialize( integration_order,local_dimension_slave, number_of_slave_nodes );
//                 }
            }
            else if (number_of_slave_nodes == 4) // TODO: Complete
            {
//                 if (integration_order == )
//                 {
//                 }
//                 else
//                 {
                    mUseManualColocationIntegration = true;
                    mColocationIntegration.Initialize( integration_order,local_dimension_slave, number_of_slave_nodes );
//                 }
            }
        }
    }
    else
    {
        mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_5;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContactCondition::GetNodalDeltaMovements( 
    Vector& rValues,
    const int& rNode 
    )
{
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();;

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

    if( dimension == 3 ) 
    {
        rValues[2] = CurrentValueVector[2] - PreviousValueVector[2];
    }
}

//************************************************************************************
//************************************************************************************

Vector& MortarContactCondition::GetCurrentValue( 
    const Variable<array_1d<double,3> >& rVariable,
    Vector& rValue,
    const unsigned int& rNode
    )
{
    KRATOS_TRY;

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();;

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

Vector& MortarContactCondition::GetPreviousValue( 
    const Variable<array_1d<double,3> >& rVariable,
    Vector& rValue,
    const unsigned int& rNode
    )
{
    KRATOS_TRY;

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();;

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

void MortarContactCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
}

void MortarContactCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
}

} // Namespace Kratos
