// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Mohamed Khalil
//                   Vicente Mataix Ferrándiz
//

// System includes

// External includes

// Project includes
#include "custom_conditions/mortar_contact_3D_condition.hpp"
#include "custom_utilities/contact_utilities.h"
#include "utilities/math_utils.h"
#include "contact_structural_mechanics_application.h"
#include "contact_structural_mechanics_application_variables.h"

namespace Kratos
{
/**
 * Flags related to the condition computation
 */
KRATOS_CREATE_LOCAL_FLAG( MortarContact3DCondition, COMPUTE_RHS_VECTOR,                 0 );
KRATOS_CREATE_LOCAL_FLAG( MortarContact3DCondition, COMPUTE_LHS_MATRIX,                 1 );
KRATOS_CREATE_LOCAL_FLAG( MortarContact3DCondition, COMPUTE_RHS_VECTOR_WITH_COMPONENTS, 2 );
KRATOS_CREATE_LOCAL_FLAG( MortarContact3DCondition, COMPUTE_LHS_MATRIX_WITH_COMPONENTS, 3 );

/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

MortarContact3DCondition::MortarContact3DCondition(
    IndexType NewId,
    GeometryType::Pointer pGeometry) :
        Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

MortarContact3DCondition::MortarContact3DCondition(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties ) :
    Condition( NewId, pGeometry, pProperties )
{
}
/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

MortarContact3DCondition::MortarContact3DCondition( 
    MortarContact3DCondition const& rOther ) :
    Condition( rOther )
{
    //DO NOT ADD DOFS HERE!!!
}

/************************************* OPERATIONS **********************************/
/***********************************************************************************/

Condition::Pointer MortarContact3DCondition::Create( 
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties ) const
{
    return ( Condition::Pointer( new MortarContact3DCondition( NewId, GetGeometry().Create( rThisNodes ), pProperties ) ) );
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer MortarContact3DCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer( new MortarContact3DCondition( NewId, pGeom, pProperties ) );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

MortarContact3DCondition::~MortarContact3DCondition( )
{
}

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::Initialize( ) 
{
//    KRATOS_TRY;
//    
//    mUseColocationIntegration = true;
//    if( this->Has(INTEGRATION_ORDER_CONTACT) )
//    {
//        mColocationIntegration.Initialize( this->GetValue(INTEGRATION_ORDER_CONTACT) );
//    }
//    else
//    {
//        std::cout << RED << "Guys. There is a bug !!" << std::endl;
//    }
//
//    KRATOS_CATCH( "" );
    KRATOS_TRY;
    mUseColocationIntegration = false;

    if( GetProperties().Has(INTEGRATION_ORDER_CONTACT) )
    {
        if (GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) == 1)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
        }
        else if (GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) == 2)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
        }
        else if (GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) == 3)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
        }
        else if (GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) == 4)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_4;
        }
        else if (GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) == 5)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_5;
        }
        else
        {
//            std::cout << "The number of integration points is not defined for Gauss quadrature. Using colocation integration with "
//                << GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) << " points." << std::endl;
            mUseColocationIntegration = true;
            mColocationIntegration.Initialize( GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) );
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

void MortarContact3DCondition::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // Populate the vector of master elements
    std::vector<contact_container> * all_containers = this->GetValue(CONTACT_CONTAINERS);
    mThisMasterElements.clear();
    mThisMasterElements.resize( all_containers->size( ) );
    
//     const double ActiveCheckFactor = GetProperties().GetValue(ACTIVE_CHECK_FACTOR); 
    
    for ( unsigned int i_cond = 0; i_cond < all_containers->size(); ++i_cond )
    {
        mThisMasterElements[i_cond] = (*all_containers)[i_cond].condition;
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    // TODO: Add things if needed

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // TODO: Add things if needed
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::CalculateLocalSystem( 
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
    LocalSystem.CalculationFlags.Set(MortarContact3DCondition::COMPUTE_LHS_MATRIX_WITH_COMPONENTS, true);
    LocalSystem.CalculationFlags.Set(MortarContact3DCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS, true);

    //Initialize sizes for the system components:
    if ( rLHSVariables.size( ) != rLeftHandSideMatrices.size( ) )
    {
        rLeftHandSideMatrices.resize( rLHSVariables.size( ) );
    }

    if ( rRHSVariables.size( ) != rRightHandSideVectors.size( ) )
    {
        rRightHandSideVectors.resize( rRHSVariables.size( ) );
    }

    LocalSystem.CalculationFlags.Set(MortarContact3DCondition::COMPUTE_LHS_MATRIX, true);
    for ( unsigned int i = 0; i < rLeftHandSideMatrices.size( ); i++ )
    {
        // Note: rRightHandSideVectors.size() > 0
        this->InitializeSystemMatrices( rLeftHandSideMatrices[i], rRightHandSideVectors[0],LocalSystem.CalculationFlags );
    }

    LocalSystem.CalculationFlags.Set( MortarContact3DCondition::COMPUTE_RHS_VECTOR, true );
    LocalSystem.CalculationFlags.Set( MortarContact3DCondition::COMPUTE_LHS_MATRIX, false ); // Temporarily only
    for ( unsigned int i = 0; i < rRightHandSideVectors.size( ); i++ )
    {
        // Note: rLeftHandSideMatrices.size() > 0
        this->InitializeSystemMatrices( rLeftHandSideMatrices[0], rRightHandSideVectors[i], LocalSystem.CalculationFlags  );
    }
    LocalSystem.CalculationFlags.Set( MortarContact3DCondition::COMPUTE_LHS_MATRIX, true ); // Reactivated again

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

void MortarContact3DCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;
    
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContact3DCondition::COMPUTE_LHS_MATRIX, true );
    LocalSystem.CalculationFlags.Set( MortarContact3DCondition::COMPUTE_RHS_VECTOR, true );

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

void MortarContact3DCondition::CalculateLeftHandSide( 
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContact3DCondition::COMPUTE_LHS_MATRIX, true );

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

void MortarContact3DCondition::CalculateLeftHandSide( 
    std::vector< MatrixType >& rLeftHandSideMatrices,
    const std::vector< Variable< MatrixType > >& rLHSVariables,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContact3DCondition::COMPUTE_LHS_MATRIX, true );
    LocalSystem.CalculationFlags.Set( MortarContact3DCondition::COMPUTE_LHS_MATRIX_WITH_COMPONENTS, true );

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

void MortarContact3DCondition::CalculateRightHandSide( 
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContact3DCondition::COMPUTE_RHS_VECTOR, true);

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

void MortarContact3DCondition::CalculateRightHandSide( 
    std::vector< VectorType >& rRightHandSideVectors,
    const std::vector< Variable< VectorType > >& rRHSVariables,
    ProcessInfo& rCurrentProcessInfo )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContact3DCondition::COMPUTE_RHS_VECTOR, true );
    LocalSystem.CalculationFlags.Set( MortarContact3DCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS, true );

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

void MortarContact3DCondition::InitializeSystemMatrices( 
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    Flags& rCalculationFlags 
    )
{
    const unsigned int condition_size = this->CalculateConditionSize( );
    
    // Resizing as needed the LHS
    if ( rCalculationFlags.Is( MortarContact3DCondition::COMPUTE_LHS_MATRIX ) ) // Calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != condition_size )
        {
            rLeftHandSideMatrix.resize( condition_size, condition_size, false );
        }
        noalias( rLeftHandSideMatrix ) = ZeroMatrix( condition_size, condition_size ); // Resetting LHS
    }

    // Resizing as needed the RHS
    if ( rCalculationFlags.Is( MortarContact3DCondition::COMPUTE_RHS_VECTOR ) ) // Calculation of the matrix is required
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

void MortarContact3DCondition::CalculateMassMatrix( 
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

void MortarContact3DCondition::CalculateDampingMatrix( 
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

const unsigned int MortarContact3DCondition::CalculateConditionSize( )
{
    const unsigned int dimension = 3;

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


void MortarContact3DCondition::CalculateConditionSystem( 
    LocalSystemComponents& rLocalSystem,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
    
    // // LOG_GENERAL( CYAN, "Function: ", __FUNCTION__ )
    
    // Create and initialize condition variables:
    GeneralVariables Variables;

    // Working space dimension
    const unsigned int dimension = 3;

    // Slave segment info
    const unsigned int number_of_nodes_slave     = GetGeometry().PointsNumber( );
    const unsigned int number_of_elements_master = mThisMasterElements.size( );
    
    // Reading integration points
//    const GeometryType::IntegrationPointsArrayType& integration_points = mUseColocationIntegration ?
//                                                                         mColocationIntegration.IntegrationPoints( ):
//                                                                         GetGeometry( ).IntegrationPoints( mThisIntegrationMethod );
    
    const GeometryType::IntegrationPointsArrayType& integration_points = mColocationIntegration.IntegrationPoints( );
    
    for (unsigned int iNode = 0; iNode < number_of_nodes_slave; iNode++)
    {
        GetGeometry()[iNode].GetValue(WEIGHTED_GAP) = 1.0e9; // Value extremely large 
    }
    
    for (unsigned int i_master_elem = 0; i_master_elem < number_of_elements_master; ++i_master_elem)
    {
        // Initialize general variables for the current master element
        this->InitializeGeneralVariables( Variables, i_master_elem );
        
        // Master segment info
        const GeometryType& current_master_element = Variables.GetMasterElement( );
        const unsigned int& number_of_nodes_master = current_master_element.PointsNumber( );

        // Compute mortar condition matrices
        MortarConditionMatrices ThisMortarConditionMatrices;
        ThisMortarConditionMatrices.Initialize( number_of_nodes_master, number_of_nodes_slave, dimension );
        
        // Weighted gaps
        VectorType& gn = ThisMortarConditionMatrices.NodalWeightedGaps;
        
      //  LOG_CONDITION_HEADER( Variables.GetMasterElement( ), GetGeometry( ) )
        
        double aux_integ_area = 0;
        
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            // LOG_INFO( PointNumber )
            this->CalculateKinematics( Variables, PointNumber, i_master_elem );
            const double IntegrationWeight = integration_points[PointNumber].Weight( ) * Variables.ColocationWeightCoeff;
            
            aux_integ_area += IntegrationWeight;
            
            if( IntegrationWeight > 1e-10 ) 
            {
                this->CalculateDAndM( Variables, IntegrationWeight, ThisMortarConditionMatrices );

                // DIRECTIONAL DERIVATIVES
                this->CalculateDeltaDetJSlave(                     Variables, ThisMortarConditionMatrices );
                //            this->CalculateDeltaPhi(        Variables, IntegrationWeight, ThisMortarConditionMatrices );
                this->CalculateDeltaIntegrationSegmentPoint(       Variables, ThisMortarConditionMatrices );
                this->CalculateLambdaTDeltaBco( Variables, IntegrationWeight, ThisMortarConditionMatrices );
                this->CalculateN(               Variables, IntegrationWeight, ThisMortarConditionMatrices );
            }
        }
        if( aux_integ_area > 1e-9 ) // if no integration_pt is projectable on the master then the system matix and RHS vector are zeros
        {
            for (unsigned int iNode = 0; iNode < number_of_nodes_slave; iNode++)
            {
                if (GetGeometry()[iNode].GetValue(WEIGHTED_GAP) > gn[iNode])
                {
                    GetGeometry()[iNode].GetValue(WEIGHTED_GAP) = gn[iNode]; // Saving the smallest one
                }
            }

            // Calculation of the matrix is required
            if ( rLocalSystem.CalculationFlags.Is( MortarContact3DCondition::COMPUTE_LHS_MATRIX ) ||
                rLocalSystem.CalculationFlags.Is( MortarContact3DCondition::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
            {
                // Contributions to stiffness matrix calculated on the reference config

                this->CalculateAndAddLHS( rLocalSystem, Variables, ThisMortarConditionMatrices );
            }

            // Calculation of the vector is required
            if ( rLocalSystem.CalculationFlags.Is( MortarContact3DCondition::COMPUTE_RHS_VECTOR ) ||
                rLocalSystem.CalculationFlags.Is( MortarContact3DCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
            {
                // Contribution to previous step contact force and residuals vector
                this->CalculateAndAddRHS( rLocalSystem, Variables, ThisMortarConditionMatrices );
            }
        }
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::InitializeGeneralVariables(
    GeneralVariables& rVariables,
    const unsigned int& rMasterElementIndex 
    )
{
    // Working space dimension
    const unsigned int dimension = 3;

    // Slave segment info
    GeometryType& current_slave_element = this->GetGeometry();
    const unsigned int number_of_nodes_slave = current_slave_element.PointsNumber();

    // Master segment info
    GeometryType& current_master_element = mThisMasterElements[rMasterElementIndex]->GetGeometry();
    const unsigned int number_of_nodes_master = current_master_element.PointsNumber();

    // Slave element info
    rVariables.Initialize(number_of_nodes_master, number_of_nodes_slave, dimension );
    
    rVariables.SetMasterElement( current_master_element );
    rVariables.SetMasterElementIndex( rMasterElementIndex );
}

//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************

void MortarContact3DCondition::CalculateKinematics( 
    GeneralVariables& rVariables,
    const double& rPointNumber,
    const unsigned int& rPairIndex
    )
{
    KRATOS_TRY;
    
    // // LOG_GENERAL( CYAN, "Function: ", __FUNCTION__ )
    
    /* DEFINITIONS */
    GeometryType& slave_nodes  = GetGeometry( );
    GeometryType& master_nodes = rVariables.GetMasterElement( );
    const unsigned int number_of_master_nodes = master_nodes.PointsNumber( );
    const unsigned int number_of_slave_nodes  =  slave_nodes.PointsNumber( );

    const GeometryType::IntegrationPointsArrayType& integration_points = mUseColocationIntegration ?
                                                                         mColocationIntegration.IntegrationPoints( ):
                                                                         GetGeometry( ).IntegrationPoints( mThisIntegrationMethod );
    
//    const GeometryType::IntegrationPointsArrayType& integration_points = mColocationIntegration.IntegrationPoints( );
    
    const PointType& local_point = integration_points[rPointNumber].Coordinates();
    
    /* RESIZE MATRICES AND VECTORS */
    rVariables.Phi_LagrangeMultipliers.resize( number_of_slave_nodes, false );
    rVariables.N_Master.resize( number_of_master_nodes, false );
    rVariables.N_Slave.resize( number_of_slave_nodes, false );

    rVariables.DPhi_De_LagrangeMultipliers.resize( number_of_slave_nodes, slave_nodes.LocalSpaceDimension( ), false );
    rVariables.DN_De_Master.resize( number_of_master_nodes, master_nodes.LocalSpaceDimension( ), false );
    rVariables.DN_De_Slave.resize( number_of_slave_nodes, slave_nodes.LocalSpaceDimension( ), false );
    
    /*  POPULATE MATRICES AND VECTORS */
    
    /// SLAVE CONDITION ///
    rVariables.N_Slave = slave_nodes.ShapeFunctionsValues( rVariables.N_Slave, local_point.Coordinates() );
    rVariables.DN_De_Slave  =  slave_nodes.ShapeFunctionsLocalGradients( rVariables.DN_De_Slave , local_point );
    
    /// LM SPACE CHOICE - STANDARD OR DUAL ///
    rVariables.Phi_LagrangeMultipliers = rVariables.N_Slave;
    rVariables.DPhi_De_LagrangeMultipliers = rVariables.DN_De_Slave;
//    rVariables.DPhi_De_LagrangeMultipliers = LagrangeMultiplierShapeFunctionLocalGradient( local_point[0], local_point[1] );
//    rVariables.Phi_LagrangeMultipliers= LagrangeMultiplierShapeFunctionValue( local_point[0], local_point[1] );
    
    /// MASTER CONDITION ///
    this->MasterShapeFunctionValue( rVariables, local_point);

    /* CALCULATE JACOBIAN AND JACOBIAN DETERMINANT */
    /*
     * This is the jacobian of the contact element.
     */
    slave_nodes.Jacobian( rVariables.j_Slave, local_point.Coordinates( ) );
    rVariables.DetJSlave = ContactUtilities::ContactElementDetJacobian( rVariables.j_Slave );
    
    KRATOS_CATCH( "" );
}

/**********************************************************************************/
/*********************** AUXILLIARY SHAPE FUNCTIONS METHODS ***********************/
/**********************************************************************************/

const Vector MortarContact3DCondition::LagrangeMultiplierShapeFunctionValue(
    const double xi_local,
    const double eta_local
)
{
    // NOTE: For more information look at the thesis of Popp page 93/236
    
    const unsigned int num_slave_nodes = GetGeometry( ).PointsNumber();
    Vector Phi = ZeroVector( num_slave_nodes );

    if (num_slave_nodes == 3) // triangular element
    {
        Phi[0] = 3 - 4 * xi_local - 4 * eta_local;
        Phi[1] = 4 * xi_local  - 1;
        Phi[3] = 4 * eta_local - 1;
    }
    else if (num_slave_nodes == 4) // quad element
    {
        array_1d<double,3> aux_coordinates = ZeroVector(3);
        aux_coordinates[0] =  xi_local;
        aux_coordinates[1] = eta_local;
        Vector N = ZeroVector(4);
        GetGeometry().ShapeFunctionsValues(N, aux_coordinates);
        
        Phi( 0 ) =  4*N(0) - 2*N(1) +   N(2) - 2*N(3); 
        Phi( 1 ) = -2*N(0) + 4*N(1) - 2*N(2) +   N(3);
        Phi( 2 ) =    N(0) - 2*N(1) + 4*N(2) - 2*N(3);
        Phi( 3 ) = -2*N(0) +   N(1) - 2*N(2) + 4*N(3);
    }
    else
    {
        KRATOS_THROW_ERROR( std::logic_error, "Higher order contact elements are not implemented in 3D.", "" )
    }
        
    return Phi;
}

/***********************************************************************************/
/***********************************************************************************/

const Matrix MortarContact3DCondition::LagrangeMultiplierShapeFunctionLocalGradient(
    const double xi_local,
    const double eta_local
)
{
    
    // // LOG_GENERAL( CYAN, "Function: ", __FUNCTION__ )
    
    const unsigned int num_slave_nodes = GetGeometry( ).PointsNumber( );
    const unsigned int local_dimension_slave = GetGeometry( ).LocalSpaceDimension( );
    Matrix DPhi_De = ZeroMatrix( num_slave_nodes, local_dimension_slave );

    if (num_slave_nodes == 3) // triangular element
    {
        // DPhi/DXi
        DPhi_De( 0, 0 ) = -4.0; 
        DPhi_De( 1, 0 ) =  4.0;
        DPhi_De( 2, 0 ) =  0.0;
        
        // DPhi/DEta
        DPhi_De( 0, 1 ) = -4.0; 
        DPhi_De( 1, 1 ) =  0.0;
        DPhi_De( 2, 1 ) =  4.0;
    }
    else if (num_slave_nodes == 4) // quad element
    {
        array_1d<double,3> aux_coordinates = ZeroVector(3);
        aux_coordinates[0] =  xi_local;
        aux_coordinates[1] = eta_local;
        Matrix DN_De = ZeroMatrix(4, local_dimension_slave);
        GetGeometry().ShapeFunctionsLocalGradients(DN_De, aux_coordinates);

        // DPhi/DXi
        DPhi_De( 0, 0 ) =  4*DN_De(0,0) - 2*DN_De(1,0) +   DN_De(2,0) - 2*DN_De(3,0); 
        DPhi_De( 1, 0 ) = -2*DN_De(0,0) + 4*DN_De(1,0) - 2*DN_De(2,0) +   DN_De(3,0);
        DPhi_De( 2, 0 ) =    DN_De(0,0) - 2*DN_De(1,0) + 4*DN_De(2,0) - 2*DN_De(3,0);
        DPhi_De( 3, 0 ) = -2*DN_De(0,0) +   DN_De(1,0) - 2*DN_De(2,0) + 4*DN_De(3,0);

        // DPhi/DEta
        DPhi_De( 0, 1 ) =  4*DN_De(0,1) - 2*DN_De(1,1) +   DN_De(2,1) - 2*DN_De(3,1); 
        DPhi_De( 1, 1 ) = -2*DN_De(0,1) + 4*DN_De(1,1) - 2*DN_De(2,1) +   DN_De(3,1);
        DPhi_De( 2, 1 ) =    DN_De(0,1) - 2*DN_De(1,1) + 4*DN_De(2,1) - 2*DN_De(3,1);
        DPhi_De( 3, 1 ) = -2*DN_De(0,1) +   DN_De(1,1) - 2*DN_De(2,1) + 4*DN_De(3,1);
    }
    else
    {
        KRATOS_THROW_ERROR( std::logic_error, "Higher order contact elements are not implemented in 3D.", "" )
    }

    return DPhi_De;
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::MasterShapeFunctionValue(
    GeneralVariables& rVariables,
    const PointType& local_point 
    )
{
    // // LOG_GENERAL( CYAN, "Function: ", __FUNCTION__ )
    
    GeometryType& master_seg = rVariables.GetMasterElement( );
    rVariables.N_Master     = ZeroVector( master_seg.PointsNumber( ) );
    rVariables.DN_De_Master = ZeroMatrix( master_seg.PointsNumber( ), master_seg.LocalSpaceDimension( ) );

    PointType projected_gp_global;
    rVariables.IntegrationPointNormalGap = 0.0;
    VectorType& normal = rVariables.IntegrationPointNormalVector;
    normal = ContactUtilities::GaussPointNormal(rVariables.N_Slave, GetGeometry());
    
    GeometryType::CoordinatesArrayType slave_gp_global;
    this->GetGeometry( ).GlobalCoordinates( slave_gp_global, local_point );
    ContactUtilities::ProjectDirection( master_seg, slave_gp_global, projected_gp_global, rVariables.IntegrationPointNormalGap, -normal ); // The opposite direction
    
    GeometryType::CoordinatesArrayType projected_gp_local;
    if( master_seg.IsInside( projected_gp_global.Coordinates( ), projected_gp_local ) )
    {
        rVariables.N_Master     = master_seg.ShapeFunctionsValues(         rVariables.N_Master,     projected_gp_local );         
        rVariables.DN_De_Master = master_seg.ShapeFunctionsLocalGradients( rVariables.DN_De_Master, projected_gp_local );
        master_seg.Jacobian( rVariables.j_Master, projected_gp_local );
        rVariables.ColocationWeightCoeff = 1;
    }
    else
    {
        if( mUseColocationIntegration == true )
        {
//            std::cout << LT_RED << "Integration pt is outside. Coords: " << local_point.Coordinates( ) << RESET << std::endl;
            rVariables.ColocationWeightCoeff = 0;
        }
    }
}

/***********************************************************************************/
/*************** METHODS TO CALCULATE THE CONTACT CONDITION MATRICES ***************/
/***********************************************************************************/

void MortarContact3DCondition::CalculateDAndM( 
    GeneralVariables& rVariables,
    const double& rIntegrationWeight,
    MortarConditionMatrices& ThisMortarConditionMatrices
    )
{
    // // LOG_GENERAL( CYAN, "Function: ", __FUNCTION__ )
    
    const unsigned int num_slave_nodes  = GetGeometry().PointsNumber( );
    const unsigned int num_master_nodes = rVariables.GetMasterElement( ).PointsNumber( );

    const Vector& N_s = rVariables.N_Slave;                   
    const Vector& N_m = rVariables.N_Master;                  
    const Vector& Phi = rVariables.Phi_LagrangeMultipliers;  
    const double& gap = rVariables.IntegrationPointNormalGap;
    const double& J_s = rVariables.DetJSlave;                
    
    MatrixType& D  = ThisMortarConditionMatrices.D;
    MatrixType& M  = ThisMortarConditionMatrices.M;
    VectorType& gn = ThisMortarConditionMatrices.NodalWeightedGaps;

    // For all the nodes
    for ( unsigned int i_slave = 0; i_slave < num_slave_nodes; ++i_slave )
    {
        gn[i_slave] += rIntegrationWeight * Phi( i_slave ) * gap * J_s;

//        D( i_slave, i_slave ) += rIntegrationWeight * Phi( i_slave ) * N_s( i_slave ) * J_s;
        for ( unsigned int j_slave = 0; j_slave < num_slave_nodes; ++j_slave )
        {
            D( i_slave, j_slave  ) += rIntegrationWeight * Phi( i_slave ) * N_s( j_slave ) * J_s;
        }

        for ( unsigned int j_master = 0; j_master < num_master_nodes; ++j_master )
        {
            M( i_slave, j_master ) += rIntegrationWeight * Phi( i_slave ) * N_m( j_master ) * J_s;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::CalculateDeltaDetJSlave(
   GeneralVariables& rVariables,
   MortarConditionMatrices& ThisMortarConditionMatrices
   )
{
   // LOG_GENERAL( LT_YELLOW, "|.... ", __FUNCTION__ );

   // prerequisites
   const unsigned int dimension = 3;
   const unsigned int num_slave_nodes  = GetGeometry( ).PointsNumber( );

   const VectorType& DN_Dxi  = column( rVariables.DN_De_Slave, 0 );
   const VectorType& DN_Deta = column( rVariables.DN_De_Slave, 1 );
   const VectorType& J_xi    = column( rVariables.j_Slave, 0 );
   const VectorType& J_eta   = column( rVariables.j_Slave, 1 );
   
   const VectorType& n = rVariables.IntegrationPointNormalVector;

   MatrixType& DeltaJ_s = ThisMortarConditionMatrices.DeltaDetJSlave;
   
   MatrixType Delta_Jxi_x_Jeta = Matrix(dimension, dimension);

   /*
    * In 3D, Delta ||J|| is more complex than in 2D, because we have two local directions
    * Also, in Popp's thesis, a different approach (mortar segments) is followed in 3D, so we don't have a reference for this
    * but we can rely on our intuition and engineering sense.. here is the rational
    * 
    * RULE 1: Delta ||a|| = (a.Delta_a) / ||a||
    * RULE 2: From Popp's thesis (eqn. A.24), the Jacobian of the integration cell is defined by the cross product of the two non-unit local vectors
    * Remember, the integration cell is a linear tri-element.. xi_cell = x2-x1 and eta_cell = x3-x1
    * 
    * From that, we can generalize that the Jacobian of a 2D element in 3D space is ( vector_xi x vector_eta )
    * vector_xi is the first column of the Jacobian matrix in rVariables.j_Slave, let's call it J_xi
    * vector_eta is the second column of the Jacobian matrix in rVariables.j_Slave, let's call it J_eta
    * 
    * therefore, J = J_xi x J_eta, which is the non-unit normal vector at the current integration point => J = J_xi x J_eta = n_hat
    * substitute back into RULE 1 => Delta ||J|| = (n_hat/||n_hat||) . Delta(J_xi x J_eta) = n . Delta(J_xi x J_eta)
    *                                            = n . ( Delta_J_xi x J_eta - Delta_J_eta x J_xi ) 
    * 
    * J_xi and J_eta define the tangent vectors in the element's plane
    * the derivation of how this bracket is calculated can be found in ContactUtilities::ComputeDeltaNodesMeanNormalModelPart
    * it is similar to how Delta_ne_adj 
    */
   for ( unsigned int i_slave = 0, i = 0; i_slave < num_slave_nodes; ++i_slave, i += dimension )
   {
       if ( GetGeometry( )[i_slave].Is( ACTIVE ) )
       {
           Delta_Jxi_x_Jeta(0,0) = 0.0;
           Delta_Jxi_x_Jeta(0,1) =  J_eta(2) * DN_Dxi(i_slave) - J_xi(2) * DN_Deta(i_slave); 
           Delta_Jxi_x_Jeta(0,2) = -J_eta(1) * DN_Dxi(i_slave) + J_xi(1) * DN_Deta(i_slave); 
           Delta_Jxi_x_Jeta(1,0) = -J_eta(2) * DN_Dxi(i_slave) + J_xi(2) * DN_Deta(i_slave); 
           Delta_Jxi_x_Jeta(1,1) = 0.0;
           Delta_Jxi_x_Jeta(1,2) =  J_eta(0) * DN_Dxi(i_slave) - J_xi(0) * DN_Deta(i_slave);
           Delta_Jxi_x_Jeta(2,0) =  J_eta(1) * DN_Dxi(i_slave) - J_xi(1) * DN_Deta(i_slave); 
           Delta_Jxi_x_Jeta(2,1) = -J_eta(0) * DN_Dxi(i_slave) + J_xi(0) * DN_Deta(i_slave); 
           Delta_Jxi_x_Jeta(2,2) = 0.0;
           
           DeltaJ_s( 0, i     ) = inner_prod( n, column( Delta_Jxi_x_Jeta, 0 ) );
           DeltaJ_s( 0, i + 1 ) = inner_prod( n, column( Delta_Jxi_x_Jeta, 1 ) );
           DeltaJ_s( 0, i + 2 ) = inner_prod( n, column( Delta_Jxi_x_Jeta, 2 ) );
       }
   }
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::CalculateDeltaIntegrationSegmentPoint(
    GeneralVariables& rVariables,
    MortarConditionMatrices& ThisMortarConditionMatrices
    )
{
    // LOG_GENERAL( LT_YELLOW, "|.... ", __FUNCTION__ );

    // Geometries
    const unsigned int& dimension = 3;
    const GeometryType& slave_nodes  = GetGeometry( );
    const GeometryType& master_nodes = rVariables.GetMasterElement( );
    const unsigned int& num_slave_nodes  =  slave_nodes.PointsNumber( );
    const unsigned int& num_master_nodes = master_nodes.PointsNumber( );
    const unsigned int num_total_nodes = num_master_nodes + num_slave_nodes;
    
    // Shape fucntions and their dertivatives
    const VectorType& N_m               = rVariables.N_Master;
    const VectorType& N_s               = rVariables.N_Slave;
    const array_1d<double, 3>& ng       = rVariables.IntegrationPointNormalVector;
    const array_1d<double, 3>& x2hat_x1 = rVariables.IntegrationPointNormalGap * -ng;
    const VectorType& J_xi_s      = column( rVariables.j_Slave,  0 );      // J_xi_s
    const VectorType& J_eta_s     = column( rVariables.j_Slave,  1 );      // J_eta_s
    const VectorType& J_xi_m      = column( rVariables.j_Master, 0 );      // J_xi_m
    const VectorType& J_eta_m     = column( rVariables.j_Master, 1 );      // J_eta_m
    
    /*
     * Master side integration point (the projection of the slave integration point)                             
     * Even using collocation, this won't evaluate to zero because of the linearization of the projection operator
     * Note: using collocation, this is calculated only when the slave integration point is contact the master segment
     * ( i.e. when  the weight coefficient is not zero )
     * This is again complex in 3D because we have two local direction to account for
     * 
     * DERIVATION OF THIS EQUATION - because it is not in Popp's thesis so we use our engineering sense and math intuition
     * 
     * Looking at the definition of the projection in eq. 3.30 in Popp's thesis
     * let's define: gap = x2hat - x1 
     * 
     * taking the derivative of 3.30 (as he suggested in 2D case in the A.1.4)
     * Delta_gap x ng + gap x Delta_ng = D1 + D2 = 0
     * 
     * REMARK: DeltaXi2 and DeltaEta2 are inside Delta_gap
     * ng = -( t_eta x t_xi )
     * 
     ******** SIMPLIFICATION OF D1 
     * 
     * Cross Product Property: a x (b x c) = b.a.c - c.a.b
     * 
     * therefore, D1 = t_eta.( Delta_gap . t_xi ) + t_xi.( Delta_gap . t_eta )
     *               = ( outer( t_eta, t_xi ) - outer( t_eta, t_xi ) ) . Delta_gap
     *               = T_tau . Delta_gap
     *               
     * Delta_gap = [ N_m,xi . x_m * DeltaXi2 + N_m,eta . x_m * DeltaEta2 ] + [ N_m . Delta_xm ] - [ N_s . Delta_xs ]
     *           = [    J_xi' . I . DeltaXi2 +     J_eta . I . DeltaEta2 ] + [ N_m . Delta_xm ] - [ N_s . Delta_xs ]
     *           
     * substitute back in D1
     * D1 = T_tau . [    J_xi' . I . DeltaXi2 +    J_eta' . I . DeltaEta2 ]
     *    + T_tau . ( [ N_m . Delta_xm ] - [ N_s . Delta_xs ] )
     *    
     * from t_xi = xi and t_eta = eta and t_xi . t_eta = 0 => Delta( t_xi . t_eta ) = 0
     * DeltaEta2 = ( (t_eta ./ t_xi)' * I ) . DeltaXi2 = R . Delta_Xi2
     * 
     * substitute back in D1
     * D1 = T_tau . [ J_xi' . I + J_eta' . I . R ] . DeltaXi2 + T_tau . ( [ N_m . Delta_xm ] - [ N_s . Delta_xs ] )
     *    =                             coeff_DXi2 . DeltaXi2 + T_tau . ( [ N_m . Delta_xm ] - [ N_s . Delta_xs ] )
     *    
     *    
     ******** SIMPLIFICATION OF D2 
     * 
     * D2 = gap x Delta_ng = T_gap . Delta_ng => both are clearly defined in the code below
     * 
     * 
     ******** FINDING DeltaXi2
     * 
     * coeff_DXi2 . DeltaXi2 + T_tau . ( [ N_m . Delta_xm ] - [ N_s . Delta_xs ] ) + T_gap . Delta_ng = 0
     * 
     * DeltaXi2 = -inv( coeff_DXi2 ) . ( T_gap . Delta_ng + T_tau . ( [ N_m . Delta_xm ] - [ N_s . Delta_xs ] ) )
     * 
     * -inv( coeff_DXi2 ) . T_tau = -inv( T_J ) => this is an inverse of a diagonal matrix
     * 
     * DeltaXi2 = -inv( coeff_DXi2 ) . T_gap . Delta_ng - inv( T_J ). ( [ N_m . Delta_xm ] - [ N_s . Delta_xs ] )
     *
     */

    // Integration Point Directional Derivatives
    MatrixType& DeltaXi_gp_2  = ThisMortarConditionMatrices.DeltaEtaMasterIntegrationPoint;
    MatrixType& DeltaEta_gp_2 = ThisMortarConditionMatrices.DeltaXiMasterIntegrationPoint;
    
    
    // DEFINITIONS FOR AUXILLIARY TENSORS //
    // tangent vectors
    const double tol = 1e-10;
    array_1d<double, 3> t_xi_s  = J_xi_s  / norm_2(J_xi_s );
    array_1d<double, 3> t_eta_s = J_eta_s / norm_2(J_eta_s);
    array_1d<double, 3> t_xi_m  = J_xi_m  / norm_2(J_xi_m );
    array_1d<double, 3> t_eta_m = J_eta_m / norm_2(J_eta_m);
    
    // T_J = N_m,xi . x_m + N_eta . x_m . ( t_eta./t_xi * I ) = J_xi + J_eta .* (t_eta ./ t_xi)
    array_1d<double, 3> T_J = element_prod( t_eta_m, J_eta_m );
    if( t_xi_m[0] > tol ) T_J[0] /= t_xi_m[0]; else T_J[0] = 0;
    if( t_xi_m[1] > tol ) T_J[1] /= t_xi_m[1]; else T_J[1] = 0;
    if( t_xi_m[2] > tol ) T_J[2] /= t_xi_m[2]; else T_J[2] = 0;
    T_J += J_xi_m;
    
    // T_tau
    const MatrixType T_tau = -( outer_prod( t_eta_s, t_xi_s ) - outer_prod( t_xi_s, t_eta_s ) );
    
    // denom_matrix -> an inverted matrix analogous to "denom" in 2D
    MatrixType coeff_DXi2 = T_tau;   // Coeff_DXi2 = T_tau * ( T_J' * I )
    row(coeff_DXi2, 0) *= T_J[0]; 
    row(coeff_DXi2, 1) *= T_J[1]; 
    row(coeff_DXi2, 2) *= T_J[2]; 
    MatrixType coeff_DXi2_inv = ZeroMatrix(dimension, dimension);
    double det_T_tau_T_J = MathUtils<double>::Det3( coeff_DXi2 );
    
    if( det_T_tau_T_J > tol )
    {
        MathUtils<double>::InvertMatrix3( coeff_DXi2, coeff_DXi2_inv, det_T_tau_T_J );

        // T_g -> ( x2hat - x1 ) x Delta_ng = T_g . Delta_ng
        MatrixType T_g = Matrix(dimension, dimension);
        T_g(0,0) = 0.0;
        T_g(0,1) = -x2hat_x1[2];
        T_g(0,2) =  x2hat_x1[1];
        T_g(1,0) =  x2hat_x1[2];
        T_g(1,1) = 0.0;
        T_g(1,2) = -x2hat_x1[0];
        T_g(2,0) = -x2hat_x1[1];
        T_g(2,1) =  x2hat_x1[0];
        T_g(2,2) = 0.0;

        // T_g * Delta_ng
        MatrixType Delta_ng = ZeroMatrix( dimension, dimension * num_slave_nodes );
        for ( unsigned int i_slave = 0, j = 0; i_slave < num_slave_nodes; ++i_slave, j += dimension )
        {
            subrange( Delta_ng,
                      0,
                      dimension,
                      j,
                      j + dimension ) = N_s( i_slave ) * slave_nodes[i_slave].GetValue( DELTA_NORMAL );
        }
        const MatrixType I = IdentityMatrix( dimension, dimension );
        const MatrixType norm_Delta_ng = prod( ( I - outer_prod( ng, ng ) ), Delta_ng ) / rVariables.DetJSlave;
        Delta_ng = norm_Delta_ng;

        // ASSEMBLY OF MATRIX //
        // add Delta_ng contribution
        const MatrixType Tg_Dng = prod( T_g, Delta_ng );
        subrange( DeltaXi_gp_2,
            0, 
            dimension, 
            dimension * num_master_nodes, 
            dimension * num_total_nodes ) += prod( coeff_DXi2_inv, Tg_Dng );

        // add master shape functions contributions -> -Nm / T_J
        for ( unsigned int i_master = 0, j = 0; i_master < num_master_nodes; ++i_master, j += dimension )
        {
            if( T_J[0] > tol ) DeltaXi_gp_2( 0, j     ) -= N_m(i_master) / T_J[0];
            if( T_J[1] > tol ) DeltaXi_gp_2( 1, j + 1 ) -= N_m(i_master) / T_J[1];
            if( T_J[2] > tol ) DeltaXi_gp_2( 2, j + 2 ) -= N_m(i_master) / T_J[2];
        }

        // add slave shape functions contributions -> Ns / T_J
        for ( unsigned int i_slave = 0, j = num_master_nodes * dimension; i_slave < num_slave_nodes; ++i_slave, j += dimension )
        {
            if( T_J[0] > tol ) DeltaXi_gp_2( 0, j     ) += N_s(i_slave) / T_J[0];
            if( T_J[1] > tol ) DeltaXi_gp_2( 1, j + 1 ) += N_s(i_slave) / T_J[1];
            if( T_J[2] > tol ) DeltaXi_gp_2( 2, j + 2 ) += N_s(i_slave) / T_J[2];
        }

        // CALCULATE DeltaEta_gp_2 //
        if( t_xi_m[0] > tol ) row( DeltaEta_gp_2, 0 ) += ( t_eta_m[0]/t_xi_m[0] ) * row( DeltaXi_gp_2, 0 );
        if( t_xi_m[1] > tol ) row( DeltaEta_gp_2, 1 ) += ( t_eta_m[1]/t_xi_m[1] ) * row( DeltaXi_gp_2, 1 );
        if( t_xi_m[2] > tol ) row( DeltaEta_gp_2, 2 ) += ( t_eta_m[2]/t_xi_m[2] ) * row( DeltaXi_gp_2, 2 );
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::CalculateLambdaTDeltaBco( 
    GeneralVariables& rVariables,
    const double& rIntegrationWeight,
    MortarConditionMatrices& ThisMortarConditionMatrices
    )
{
    // LOG_GENERAL( LT_YELLOW, "|.... ", __FUNCTION__ );

    const unsigned int dimension = 3;
    const unsigned int num_slave_nodes  = GetGeometry( ).PointsNumber( );
    const unsigned int num_master_nodes = rVariables.GetMasterElement( ).PointsNumber( );
//     const unsigned int num_total_nodes  = num_slave_nodes + num_master_nodes;
    const unsigned int master_dofs      = num_master_nodes * dimension;

    const Vector& N_s      = rVariables.N_Slave;
    const Vector& N_m      = rVariables.N_Master;
    const Vector& Phi      = rVariables.Phi_LagrangeMultipliers;
    const Vector& DN_Dxi_m  = column( rVariables.DN_De_Master, 0 );  // DN_master/Dxi
    const Vector& DN_Deta_m = column( rVariables.DN_De_Master, 1 );  // DN_master/Deta
    const double& det_J_s  = rVariables.DetJSlave;

    const MatrixType& DeltaDetJ     = ThisMortarConditionMatrices.DeltaDetJSlave;
//     const MatrixType& DeltaPhi      = ThisMortarConditionMatrices.DeltaPhiLagrangeMultipliers;
    const MatrixType& DeltaXi_gp_2  = ThisMortarConditionMatrices.DeltaXiMasterIntegrationPoint;
    const MatrixType& DeltaEta_gp_2 = ThisMortarConditionMatrices.DeltaEtaMasterIntegrationPoint;
    
    MatrixType& DeltaD = ThisMortarConditionMatrices.LambdaTDeltaD;
    MatrixType& DeltaM = ThisMortarConditionMatrices.LambdaTDeltaM;
    
    VectorType lambda = ZeroVector(dimension);
    for ( unsigned int i_slave = 0; i_slave < num_slave_nodes; ++i_slave )
    {
        lambda[0] += Phi[i_slave] * GetGeometry( )[i_slave].GetValue( VECTOR_LAGRANGE_MULTIPLIER_X );
        lambda[1] += Phi[i_slave] * GetGeometry( )[i_slave].GetValue( VECTOR_LAGRANGE_MULTIPLIER_Y );
        lambda[2] += Phi[i_slave] * GetGeometry( )[i_slave].GetValue( VECTOR_LAGRANGE_MULTIPLIER_Z );
    }
    
    // DeltaD' * lambda
    for ( unsigned int k_slave = 0,  k = 0; k_slave  < num_slave_nodes;  ++k_slave,  k += dimension )
    {
        row( DeltaD, k     ) += rIntegrationWeight * N_s[k_slave] * lambda[0] * row( DeltaDetJ, 0 );
        row( DeltaD, k + 1 ) += rIntegrationWeight * N_s[k_slave] * lambda[1] * row( DeltaDetJ, 0 );
        row( DeltaD, k + 2 ) += rIntegrationWeight * N_s[k_slave] * lambda[2] * row( DeltaDetJ, 0 );
    }

    // DeltaM' * lambda
    for ( unsigned int l_master = 0, l = 0; l_master < num_master_nodes; ++l_master, l += dimension )
    {
        // the line with Delta_Phi
        subrange( DeltaM, l,     l + 1, master_dofs, DeltaM.size2( ) ) += rIntegrationWeight * lambda[0] * N_m[l_master] * DeltaDetJ;
        subrange( DeltaM, l + 1, l + 2, master_dofs, DeltaM.size2( ) ) += rIntegrationWeight * lambda[1] * N_m[l_master] * DeltaDetJ;
        subrange( DeltaM, l + 2, l + 3, master_dofs, DeltaM.size2( ) ) += rIntegrationWeight * lambda[2] * N_m[l_master] * DeltaDetJ;
        
        // the line with Delta_Nm
        row( DeltaM, l     ) += rIntegrationWeight * lambda[0] * det_J_s * ( row(DeltaXi_gp_2,0)  * DN_Dxi_m[l_master] + row(DeltaEta_gp_2,0) * DN_Deta_m[l_master] );
        row( DeltaM, l + 1 ) += rIntegrationWeight * lambda[1] * det_J_s * ( row(DeltaXi_gp_2,1)  * DN_Dxi_m[l_master] + row(DeltaEta_gp_2,1) * DN_Deta_m[l_master] );
        row( DeltaM, l + 2 ) += rIntegrationWeight * lambda[2] * det_J_s * ( row(DeltaXi_gp_2,2)  * DN_Dxi_m[l_master] + row(DeltaEta_gp_2,2) * DN_Deta_m[l_master] );
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::CalculateN(
    GeneralVariables& rVariables,
    const double& rIntegrationWeight,
    MortarConditionMatrices& ThisMortarConditionMatrices
    )
{
    // LOG_GENERAL( LT_YELLOW, "|.... ", __FUNCTION__ );

    // Contact pair variables
    const unsigned int& dimension = 3;
    const unsigned int& num_slave_nodes  = GetGeometry( ).PointsNumber( );
    const unsigned int& num_master_nodes = rVariables.GetMasterElement( ).PointsNumber( );  
    const unsigned int& num_total_nodes  = num_slave_nodes + num_master_nodes;

    // Variables
    const double& gap_nh  = rVariables.IntegrationPointNormalGap;
    const double& det_J_s = rVariables.DetJSlave; 
    const VectorType& phi = rVariables.Phi_LagrangeMultipliers; 

    // Directional Derivatives
    const MatrixType& Delta_J    = ThisMortarConditionMatrices.DeltaDetJSlave;
    
    //Calculation
    MatrixType& N = ThisMortarConditionMatrices.N;
                  
    MatrixType const_partial = ZeroMatrix( 1, dimension * num_total_nodes );
    this->CalculateDeltaDiscreteGap( rVariables, ThisMortarConditionMatrices, const_partial );
    
    // const_partial = Delta_gap_nh * || J || + gap_nh * Delta_J
    const_partial *= det_J_s;
    subrange( const_partial,
              0,
              1,
              dimension * num_master_nodes,
              dimension * num_total_nodes ) += gap_nh * Delta_J;

    for ( unsigned int i_slave = 0; i_slave < num_slave_nodes; ++i_slave )
    {
        if ( GetGeometry( )[i_slave].Is( ACTIVE ) )
        {
            subrange( N, i_slave, i_slave + 1, 0, dimension * num_total_nodes ) += const_partial * phi( i_slave ) * rIntegrationWeight;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::CalculateDeltaDiscreteGap( 
         GeneralVariables& rVariables,
         MortarConditionMatrices& ThisMortarConditionMatrices,
         MatrixType& rDeltaDiscreteGap
         )
{
    // LOG_GENERAL( YELLOW, "|........ ", __FUNCTION__ );


    // Geometries
    const unsigned int& dimension = 3;
    const GeometryType& slave_nodes  = GetGeometry( );
    const GeometryType& master_nodes = rVariables.GetMasterElement( );
    const unsigned int& num_slave_nodes  =  slave_nodes.PointsNumber( );
    const unsigned int& num_master_nodes = master_nodes.PointsNumber( );
    const unsigned int& num_total_nodes = num_slave_nodes + num_master_nodes;
    
    // Shape fucntions and their dertivatives
    const VectorType& N_s       = rVariables.N_Slave;
    const VectorType& N_m       = rVariables.N_Master;
    const VectorType& Jxi_m  = column( rVariables.j_Master, 0 );
//     const VectorType& Jeta_m = column( rVariables.j_Master, 1 );
    const array_1d<double, 3>& ng       = rVariables.IntegrationPointNormalVector;
    const array_1d<double, 3>& x2hat_x1 = rVariables.IntegrationPointNormalGap * -ng;
    
    // Directional derivative of Xi2 and Eta2 - needed even for collocation
    const MatrixType& DeltaXi2 = ThisMortarConditionMatrices.DeltaXiMasterIntegrationPoint;
    const MatrixType& DeltaXi2_x = subrange( DeltaXi2, 0, 1, 0, dimension * num_total_nodes );
    const MatrixType& DeltaXi2_y = subrange( DeltaXi2, 1, 2, 0, dimension * num_total_nodes );
    const MatrixType& DeltaXi2_z = subrange( DeltaXi2, 2, 3, 0, dimension * num_total_nodes );
    const MatrixType& DeltaEta2 = ThisMortarConditionMatrices.DeltaEtaMasterIntegrationPoint;
    const MatrixType& DeltaEta2_x = subrange( DeltaEta2, 0, 1, 0, dimension * num_total_nodes );
    const MatrixType& DeltaEta2_y = subrange( DeltaEta2, 1, 2, 0, dimension * num_total_nodes );
    const MatrixType& DeltaEta2_z = subrange( DeltaEta2, 2, 3, 0, dimension * num_total_nodes );
    
    // ng' . ( N2 . Dx2 )
    for ( unsigned int j_master = 0, j = 0; j_master < num_master_nodes; ++j_master, j += dimension )
    {
        rDeltaDiscreteGap( 0, j     ) += ng[0] * N_m(j_master); 
        rDeltaDiscreteGap( 0, j + 1 ) += ng[1] * N_m(j_master);
        rDeltaDiscreteGap( 0, j + 2 ) += ng[2] * N_m(j_master);
    }
    
    // - ng' . ( N1 . Dx1 )
    for ( unsigned int j_slave = 0, j = num_master_nodes * dimension; j_slave < num_slave_nodes; ++j_slave, j += dimension )
    {
        if ( slave_nodes[j_slave].Is( ACTIVE ) )
        {
        rDeltaDiscreteGap( 0, j     ) -= ng[0] * N_s(j_slave); 
        rDeltaDiscreteGap( 0, j + 1 ) -= ng[1] * N_s(j_slave);
        rDeltaDiscreteGap( 0, j + 2 ) -= ng[2] * N_s(j_slave);
        }
    }
    
    // ( x2hat - x1 )'. DNg 
    const MatrixType I = IdentityMatrix( dimension, dimension );
    const MatrixType Delta_ng_norm = ( I - outer_prod( ng, ng ) ) / rVariables.DetJSlave;
    for ( unsigned int j_slave = 0, j = num_master_nodes * dimension; j_slave < num_slave_nodes; ++j_slave, j += dimension )
    {
        const MatrixType Dn_j = prod( Delta_ng_norm, slave_nodes[j_slave].GetValue( DELTA_NORMAL ) );
        rDeltaDiscreteGap( 0, j     ) += N_s(j_slave) * ( x2hat_x1[0] * Dn_j( 0, 0 ) + x2hat_x1[1] * Dn_j( 1, 0 ) + x2hat_x1[2] * Dn_j( 2, 0 ) );
        rDeltaDiscreteGap( 0, j + 1 ) += N_s(j_slave) * ( x2hat_x1[0] * Dn_j( 0, 1 ) + x2hat_x1[1] * Dn_j( 1, 1 ) + x2hat_x1[2] * Dn_j( 2, 1 ) );
        rDeltaDiscreteGap( 0, j + 2 ) += N_s(j_slave) * ( x2hat_x1[0] * Dn_j( 0, 2 ) + x2hat_x1[1] * Dn_j( 1, 2 ) + x2hat_x1[2] * Dn_j( 2, 2 ) );
    }
    
    // Directional derivative of master segment
    const MatrixType DeltaNm_contribution = Jxi_m[0] * ng[0] * ( DeltaXi2_x + DeltaEta2_x )
                                          + Jxi_m[1] * ng[1] * ( DeltaXi2_y + DeltaEta2_y )
                                          + Jxi_m[2] * ng[2] * ( DeltaXi2_z + DeltaEta2_z );
    
    rDeltaDiscreteGap += DeltaNm_contribution;
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::CalculateF( MatrixType& rF )
{
    // In 3D, we calculate two rows for each node, Fsxi and Fseta
    
    // LOG_GENERAL( LT_YELLOW, "|.. ", __FUNCTION__ );

    const unsigned int dimension = 3;
    const GeometryType& slave_nodes = GetGeometry(); 
    const unsigned int num_slave_nodes  = slave_nodes.PointsNumber( );

    const MatrixType I = IdentityMatrix( dimension, dimension );
    
    unsigned int i, j = 0;
    VectorType DN_Dxi_j  = ZeroVector(num_slave_nodes);
    VectorType DN_Deta_j = ZeroVector(num_slave_nodes);
    
    if( num_slave_nodes == 3 )  // triangle
    {
        DN_Dxi_j[0]  = -1.0;
        DN_Dxi_j[1]  =  1.0;
        DN_Dxi_j[2]  =  0.0;
        
        DN_Deta_j[0] = -1.0;
        DN_Deta_j[1] =  0.0;
        DN_Deta_j[2] =  1.0;
    }
    else if( num_slave_nodes == 4 ) //quad - assuming no fancy users with 4-node tri elements
    {
        DN_Dxi_j[0] = -0.5;
        DN_Dxi_j[1] =  0.5;
        DN_Dxi_j[2] =  0.5;
        DN_Dxi_j[3] = -0.5;
        
        DN_Dxi_j[0] = -0.5;
        DN_Dxi_j[1] = -0.5;
        DN_Dxi_j[2] =  0.5;
        DN_Dxi_j[3] =  0.5;
    }
    else
    {
        KRATOS_THROW_ERROR( std::logic_error, "Mortar condition is not implemented for higher order 2D contact elements. Number of nodes: ", num_slave_nodes );
    }
    
    for (unsigned int iSlave = 0; iSlave < num_slave_nodes; iSlave++)
    {
        if( GetGeometry( )[iSlave].Is( ACTIVE ) )
        {
            i = iSlave;
            j = iSlave * dimension;

            const array_1d<double, 3>& LM = GetGeometry( )[iSlave].FastGetSolutionStepValue( VECTOR_LAGRANGE_MULTIPLIER );

            // Calculate nodal tangents
            array_1d<double, 3> t_xi  = ZeroVector(3);
            array_1d<double, 3> t_eta = ZeroVector(3);
            ContactUtilities::NodalTangents( t_xi, t_eta, slave_nodes, iSlave );
            
            const double norm_t_xi  = norm_2(t_xi );
            MatrixType DeltaTangentXi = ( I/norm_t_xi - outer_prod( t_xi, t_xi )/(norm_t_xi*norm_t_xi*norm_t_xi) ) * DN_Dxi_j[iSlave];
            
            rF( i, j     ) += DeltaTangentXi(0, 0) * LM[0] + DeltaTangentXi(0, 1) * LM[1] + DeltaTangentXi(0, 2) * LM[2];
            rF( i, j + 1 ) += DeltaTangentXi(1, 0) * LM[0] + DeltaTangentXi(1, 1) * LM[1] + DeltaTangentXi(1, 2) * LM[2];
            rF( i, j + 2 ) += DeltaTangentXi(2, 0) * LM[0] + DeltaTangentXi(2, 1) * LM[1] + DeltaTangentXi(2, 2) * LM[2];

            const double norm_t_eta = norm_2(t_eta);
            MatrixType DeltaTangentEta = ( I/norm_t_eta - outer_prod( t_eta, t_eta )/(norm_t_eta*norm_t_eta*norm_t_eta) ) * DN_Deta_j[iSlave];
            
            rF( i + 1, j     ) += DeltaTangentEta(0, 0) * LM[0] + DeltaTangentEta(0, 1) * LM[1] + DeltaTangentEta(0, 2) * LM[2];
            rF( i + 1, j + 1 ) += DeltaTangentEta(1, 0) * LM[0] + DeltaTangentEta(1, 1) * LM[1] + DeltaTangentEta(1, 2) * LM[2];
            rF( i + 1, j + 2 ) += DeltaTangentEta(2, 0) * LM[0] + DeltaTangentEta(2, 1) * LM[1] + DeltaTangentEta(2, 2) * LM[2];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::CalculateAndAddLHS(
    LocalSystemComponents& rLocalSystem,
    GeneralVariables& rVariables,
    const MortarConditionMatrices& ThisMortarConditionMatrices
    )
{
    const unsigned int dimension = 3;

    if ( rLocalSystem.CalculationFlags.Is( MortarContact3DCondition::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
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
        const unsigned int pair_size = dimension * ( rVariables.GetMasterElement( ).PointsNumber( ) + 2 * GetGeometry( ).PointsNumber( ) );       
        MatrixType LHS_contact_pair = ZeroMatrix(pair_size, pair_size);
        
        // Calculate
        this->CalculateAndAddMortarContactOperator(  LHS_contact_pair, rVariables, ThisMortarConditionMatrices );
        this->CalculateAndAddContactStiffnessMatrix( LHS_contact_pair, rVariables, ThisMortarConditionMatrices );
        this->CalculateAndAddNormalGapLinearization( LHS_contact_pair, rVariables, ThisMortarConditionMatrices );
        
        // Assemble
        this->AssembleContactPairLHSToConditionSystem( rVariables.GetMasterElementIndex( ), LHS_contact_pair, rLeftHandSideMatrix );
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::AssembleContactPairLHSToConditionSystem(
    const unsigned int rPairIndex,
    MatrixType& rPairLHS,
    MatrixType& rConditionLHS 
    )
{
    // // LOG_GENERAL( CYAN, "Function: ", __FUNCTION__ )
        
    const unsigned int dimension = 3;
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

void MortarContact3DCondition::CalculateAndAddRHS(
    LocalSystemComponents& rLocalSystem,
    GeneralVariables& rVariables,
    const MortarConditionMatrices& ThisMortarConditionMatrices
    )
{
    const unsigned int dimension = 3;
    
    if ( rLocalSystem.CalculationFlags.Is( MortarContact3DCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
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
                this->CalculateAndAddMortarContactOperator( RHS_contact_pair, rVariables, ThisMortarConditionMatrices );

                // Assemble
                this->AssembleContactPairRHSToConditionSystem( rVariables.GetMasterElementIndex( ), RHS_contact_pair, rRightHandSideVector );
                
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
        const unsigned int pair_size = dimension * ( rVariables.GetMasterElement( ).PointsNumber( ) + 2 * GetGeometry( ).PointsNumber( ) ); 
        VectorType RHS_contact_pair = ZeroVector(pair_size);
        
        // Calculate
        this->CalculateAndAddMortarContactOperator( RHS_contact_pair, rVariables, ThisMortarConditionMatrices );
        this->CalculateAndAddGap(                   RHS_contact_pair, rVariables, ThisMortarConditionMatrices );
        
        // Assemble
        this->AssembleContactPairRHSToConditionSystem( rVariables.GetMasterElementIndex( ), RHS_contact_pair, rRightHandSideVector );
    }
    
}
  
/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::AssembleContactPairRHSToConditionSystem(
    const unsigned int rPairIndex,
    VectorType& rPairRHS,
    VectorType& rConditionRHS 
    )
{
    // // LOG_GENERAL( CYAN, "Function: ", __FUNCTION__ )
    
    const unsigned int dimension = 3;
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

void MortarContact3DCondition::CalculateAndAddMortarContactOperator( 
    MatrixType& rLeftHandSideMatrix,
    GeneralVariables& rVariables,
    const MortarConditionMatrices& ThisMortarConditionMatrices
    )
{
    KRATOS_TRY;
  
    // // LOG_GENERAL( CYAN, "Function: ", __FUNCTION__ )
    
    // Contact pair variables
    const unsigned int dimension = 3;
    const unsigned int num_slave_nodes  = GetGeometry( ).PointsNumber( );
    const unsigned int num_master_nodes = rVariables.GetMasterElement( ).PointsNumber( );  
    const unsigned int num_total_nodes  = num_slave_nodes + num_master_nodes;

    const Matrix& D = ThisMortarConditionMatrices.D;
    const Matrix& M = ThisMortarConditionMatrices.M;
    
    unsigned int i = 0, j = 0;
    for ( unsigned int i_slave = 0; i_slave < num_slave_nodes; ++i_slave )
    {
        i = ( i_slave + num_total_nodes ) * dimension;

        if (GetGeometry( )[i_slave].Is(ACTIVE) == true)
        {
            // Calculate the blocks of B_co and B_co_transpose
            for ( unsigned int j_master = 0; j_master < num_master_nodes; ++j_master )
            {
                // Fill the -M and -M' parts
                j = j_master * dimension;
                double minus_M_ij = - M( i_slave, j_master );
                // uncomment for mesh tying case only
//                rLeftHandSideMatrix( i,     j     ) += minus_M_ij;
//                rLeftHandSideMatrix( i + 1, j + 1 ) += minus_M_ij;
//                rLeftHandSideMatrix( i + 2, j + 2 ) += minus_M_ij;
                                                    
                rLeftHandSideMatrix( j    , i     ) += minus_M_ij;
                rLeftHandSideMatrix( j + 1, i + 1 ) += minus_M_ij;
                rLeftHandSideMatrix( j + 2, i + 2 ) += minus_M_ij;
            }
      
            // Fill the D and D' parts
            for ( unsigned int j_slave = 0; j_slave < num_slave_nodes; ++j_slave )
            {
                j = ( j_slave + num_master_nodes ) * dimension;
                double D_ii = D( i_slave, j_slave );
                // uncomment for mesh tying case only
//                rLeftHandSideMatrix( i,     j     ) += D_ii;
//                rLeftHandSideMatrix( i + 1, j + 1 ) += D_ii;
//                rLeftHandSideMatrix( i + 2, j + 2 ) += D_ii;

                rLeftHandSideMatrix( j,     i     ) += D_ii;
                rLeftHandSideMatrix( j + 1, i + 1 ) += D_ii;
                rLeftHandSideMatrix( j + 2, i + 2 ) += D_ii;
            }
        }
        else
        {
            // We impose a 0 zero LM in the inactive nodes
            rLeftHandSideMatrix( i    , i     ) = 1.0;
            rLeftHandSideMatrix( i + 1, i + 1 ) = 1.0;
            rLeftHandSideMatrix( i + 2, i + 2 ) = 1.0;
        }
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::CalculateAndAddContactStiffnessMatrix(
        MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const MortarConditionMatrices& ThisMortarConditionMatrices)
{
    // LOG_GENERAL( MAGENTA, "|.... ", __FUNCTION__ );
    
    KRATOS_TRY
    
    // Contact pair variables
    const unsigned int dimension = 3;
    const unsigned int num_slave_nodes  = GetGeometry( ).PointsNumber( );
    const unsigned int num_master_nodes = rVariables.GetMasterElement( ).PointsNumber( );  
    const unsigned int num_total_nodes  = num_slave_nodes + num_master_nodes;
    
    subrange( rLeftHandSideMatrix,
              0,
              dimension * num_master_nodes,
              0,
              dimension * num_total_nodes ) = -ThisMortarConditionMatrices.LambdaTDeltaM;
    
    subrange( rLeftHandSideMatrix,
              dimension * num_master_nodes,
              dimension * num_total_nodes,
              dimension * num_master_nodes,
              dimension * num_total_nodes )  =  ThisMortarConditionMatrices.LambdaTDeltaD;
    
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::CalculateAndAddNormalGapLinearization(
    MatrixType& rLeftHandSideMatrix,
    GeneralVariables& rVariables,
    const MortarConditionMatrices& ThisMortarConditionMatrices )
{
    // LOG_GENERAL( MAGENTA, "|.... ", __FUNCTION__ );
    
    KRATOS_TRY
    
    // Contact pair variables
    const unsigned int dimension = 3;
    const GeometryType& slave_nodes = GetGeometry(); 
    const unsigned int num_slave_nodes  = slave_nodes.PointsNumber( );
    const unsigned int num_master_nodes = rVariables.GetMasterElement( ).PointsNumber( );  
    const unsigned int num_total_nodes  = num_slave_nodes + num_master_nodes;
  
    // Assemble N block
    const MatrixType& N = ThisMortarConditionMatrices.N;
    subslice( rLeftHandSideMatrix,
              dimension * num_total_nodes,
              dimension,
              N.size1( ),
              0,
              1,
              N.size2( ) ) += N;
    
    // Assemble F block - in 3D Fsxi and Fseta 
    MatrixType F = ZeroMatrix( num_slave_nodes * 2, num_slave_nodes * dimension );  // [ (dim-1) * slave_nodes x dim * slave_nodes ]
    this->CalculateF( F );
    subslice( rLeftHandSideMatrix,
              dimension * num_total_nodes + 1,
              dimension,
              F.size1( ),
              dimension * num_master_nodes,
              1,
              F.size2( ) ) += F; 

    
    // Assemble T Block - in 3D T_xi' and T_eta'
    for ( unsigned int i_slave = 0, i = num_total_nodes * dimension; i_slave < num_slave_nodes; ++i_slave, i+=dimension )
    {
        if (GetGeometry( )[i_slave].Is(ACTIVE) == true)
        {
            // Calculate nodal tangents
            array_1d<double, 3> t_xi  = ZeroVector(3);
            array_1d<double, 3> t_eta = ZeroVector(3);
            ContactUtilities::NodalTangents( t_xi, t_eta, slave_nodes, i_slave );
            t_xi  /= norm_2( t_xi  );
            t_eta /= norm_2( t_eta );
            
            rLeftHandSideMatrix( i + 1, i     ) = t_xi[0];
            rLeftHandSideMatrix( i + 1, i + 1 ) = t_xi[1];
            rLeftHandSideMatrix( i + 1, i + 2 ) = t_xi[2];
            
            rLeftHandSideMatrix( i + 2, i     ) = t_eta[0];
            rLeftHandSideMatrix( i + 2, i + 1 ) = t_eta[1];
            rLeftHandSideMatrix( i + 2, i + 2 ) = t_eta[2];
        }
    }
    
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/**************** AUXILLIARY METHODS FOR CONDITION RHS CONTRIBUTION ****************/
/***********************************************************************************/

void MortarContact3DCondition::CalculateAndAddMortarContactOperator( 
    VectorType& rRightHandSideVector,
    GeneralVariables& rVariables,
    const MortarConditionMatrices& ThisMortarConditionMatrices
    )
{
    KRATOS_TRY;
  
    // // LOG_GENERAL( CYAN, "Function: ", __FUNCTION__ )
    
    // Contact pair variables
    const unsigned int dimension = 3;
    const unsigned int num_slave_nodes  = GetGeometry( ).PointsNumber( );
    const unsigned int num_master_nodes = rVariables.GetMasterElement( ).PointsNumber( );  
//     const unsigned int num_total_nodes  = num_slave_nodes + num_master_nodes;

    const MatrixType& D  = ThisMortarConditionMatrices.D;
    const MatrixType& M  = ThisMortarConditionMatrices.M;
    
    // Calculate the block of r_co
    unsigned int j = 0;
    array_1d<double, 3> lagrange_multiplier = ZeroVector(3); 
    for ( unsigned int i_slave = 0; i_slave < num_slave_nodes; ++i_slave )
    {
        noalias(lagrange_multiplier) = GetGeometry( )[i_slave].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER, 0); 

        // Fill the - lambda * - M part
        for ( unsigned int j_master = 0; j_master < num_master_nodes; ++j_master )
        {
            j = j_master * dimension;
            const double minus_M_ij = - M( i_slave, j_master );
            rRightHandSideVector[ j     ] -= lagrange_multiplier[0] * minus_M_ij;
            rRightHandSideVector[ j + 1 ] -= lagrange_multiplier[1] * minus_M_ij;
            rRightHandSideVector[ j + 2 ] -= lagrange_multiplier[2] * minus_M_ij;
        }

        // Fill the - lambda *  D part
        for ( unsigned int j_slave = 0; j_slave < num_slave_nodes; ++j_slave )
        {
            j = ( num_master_nodes + j_slave  ) * dimension;
            const double D_ii = D( i_slave, j_slave );
            rRightHandSideVector[ j     ] -= lagrange_multiplier[0] * D_ii;
            rRightHandSideVector[ j + 1 ] -= lagrange_multiplier[1] * D_ii;
            rRightHandSideVector[ j + 2 ] -= lagrange_multiplier[2] * D_ii;
        }
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::CalculateAndAddGap( 
    VectorType& rRightHandSideVector,
    GeneralVariables& rVariables,
    const MortarConditionMatrices& ThisMortarConditionMatrices
    )
{
    // LOG_GENERAL( MAGENTA, "|.... ", __FUNCTION__ );

    KRATOS_TRY;
  
    // Contact pair variables
    const unsigned int dimension = 3;
    const unsigned int num_slave_nodes  = GetGeometry( ).PointsNumber( );
    const unsigned int num_master_nodes = rVariables.GetMasterElement( ).PointsNumber( );  
    const unsigned int num_total_nodes  = num_slave_nodes + num_master_nodes;

//     const MatrixType& D  = ThisMortarConditionMatrices.D;
//     const MatrixType& M  = ThisMortarConditionMatrices.M;
    const VectorType& gn = ThisMortarConditionMatrices.NodalWeightedGaps;
    
    unsigned int j = 0;
    array_1d<double, 3> lagrange_multiplier = ZeroVector(3); 
    array_1d<double, 3> t_xi  = ZeroVector(3); 
    array_1d<double, 3> t_eta = ZeroVector(3); 
    for ( unsigned int i_slave = 0; i_slave < num_slave_nodes; ++i_slave )
    {
        noalias(lagrange_multiplier) = GetGeometry( )[i_slave].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER, 0); 
            
        if (GetGeometry( )[i_slave].Is(ACTIVE) == true)
        {
            // Adding the gap to the RHS
            j = (num_total_nodes + i_slave) * dimension;
            
            // for mesh tying
//            const array_1d<double, 3> gap_decomp = gn[i_slave] * GetGeometry()[i_slave].GetValue(NORMAL);
//            rRightHandSideVector[ j     ] -= gap_decomp[0]; 
//            rRightHandSideVector[ j + 1 ] -= gap_decomp[1];
//            rRightHandSideVector[ j + 2 ] -= gap_decomp[2];
            
            // for unilateral contact
            ContactUtilities::NodalTangents( t_xi, t_eta, GetGeometry( ), i_slave );
            rRightHandSideVector[ j     ] -= -gn[i_slave];
            rRightHandSideVector[ j + 1 ] -= inner_prod( t_xi,  lagrange_multiplier ); 
            rRightHandSideVector[ j + 2 ] -= inner_prod( t_eta, lagrange_multiplier ); 
        }
        else
        {
            // Adding the zero equality to the RHS
            j = (num_total_nodes + i_slave) * dimension; 
            rRightHandSideVector[ j     ] = 0.0;
            rRightHandSideVector[ j + 1 ] = 0.0;
            rRightHandSideVector[ j + 2 ] = 0.0;
        }
    }
    
    KRATOS_CATCH( "" );
}

/******************************************************************/
/********** AUXILLIARY METHODS FOR GENERAL CALCULATIONS ***********/
/******************************************************************/

void MortarContact3DCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY;
        
    rResult.resize( 0, false );

    const std::vector<contact_container> all_conditions = *( this->GetValue( CONTACT_CONTAINERS ) );
        
    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    
    const unsigned int num_slave_nodes = GetGeometry().PointsNumber( );
    
    for ( unsigned int i_cond = 0;  i_cond < all_conditions.size( ); ++i_cond )
    {   
        GeometryType& current_master = all_conditions[i_cond].condition->GetGeometry( );
            
        // Master Nodes Displacement Equation IDs
        for ( unsigned int iNode = 0; iNode < current_master.PointsNumber( ); iNode++ )
        {
            NodeType& master_node = current_master[iNode];
            rResult.push_back( master_node.GetDof( DISPLACEMENT_X ).EquationId( ) );
            rResult.push_back( master_node.GetDof( DISPLACEMENT_Y ).EquationId( ) );
            rResult.push_back( master_node.GetDof( DISPLACEMENT_Z ).EquationId( ) );
        }
        
        // Slave Nodes Displacement Equation IDs
        for ( unsigned int i_slave = 0; i_slave < num_slave_nodes; ++i_slave )
        {
            NodeType& slave_node = GetGeometry()[ i_slave ];
            rResult.push_back( slave_node.GetDof( DISPLACEMENT_X ).EquationId( ) );
            rResult.push_back( slave_node.GetDof( DISPLACEMENT_Y ).EquationId( ) );
            rResult.push_back( slave_node.GetDof( DISPLACEMENT_Z ).EquationId( ) );
        }
        
        // Slave Nodes  Lambda Equation IDs
        for ( unsigned int i_slave = 0; i_slave < num_slave_nodes; ++i_slave )
        {
            NodeType& slave_node = GetGeometry()[ i_slave ];
            rResult.push_back( slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_X ).EquationId( ) );
            rResult.push_back( slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Y ).EquationId( ) );
            rResult.push_back( slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Z ).EquationId( ) );
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::GetDofList(
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    rConditionalDofList.resize( 0 );

    const std::vector<contact_container> all_conditions = *( this->GetValue( CONTACT_CONTAINERS ) );
    
    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    
    const unsigned int num_slave_nodes = GetGeometry().PointsNumber( );
        
    for ( unsigned int i_cond = 0; i_cond < all_conditions.size( ); ++i_cond )
    {
        GeometryType& current_master = all_conditions[i_cond].condition->GetGeometry( );
        
        // Master Nodes Displacement Equation IDs
        for ( unsigned int iNode = 0; iNode < current_master.PointsNumber( ); iNode++ )
        {
            NodeType& master_node = current_master[iNode];
            rConditionalDofList.push_back( master_node.pGetDof( DISPLACEMENT_X ) );
            rConditionalDofList.push_back( master_node.pGetDof( DISPLACEMENT_Y ) );
            rConditionalDofList.push_back( master_node.pGetDof( DISPLACEMENT_Z ) );
        }
        
        // Slave Nodes Displacement Equation IDs
        for ( unsigned int i_slave = 0; i_slave < num_slave_nodes; ++i_slave )
        {
            NodeType& slave_node = GetGeometry()[ i_slave ];
            rConditionalDofList.push_back( slave_node.pGetDof( DISPLACEMENT_X ) );
            rConditionalDofList.push_back( slave_node.pGetDof( DISPLACEMENT_Y ) );
            rConditionalDofList.push_back( slave_node.pGetDof( DISPLACEMENT_Z ) );
        }
        
        // Slave Nodes Lambda Equation IDs
        for ( unsigned int i_slave = 0; i_slave < num_slave_nodes; ++i_slave )
        {
            NodeType& slave_node = GetGeometry()[ i_slave ];
            rConditionalDofList.push_back( slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_X ) );
            rConditionalDofList.push_back( slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Y ) );
            rConditionalDofList.push_back( slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Z ) );
        }
    }

    KRATOS_CATCH( "" );
}

//******************************* GET DOUBLE VALUE *********************************/
/***********************************************************************************/

void MortarContact3DCondition::GetValueOnIntegrationPoints( 
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::CalculateOnIntegrationPoints( 
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;

    const unsigned int number_of_integration_pts = mColocationIntegration.IntegrationPoints( ).size( );
    if ( rOutput.size( ) != number_of_integration_pts )
    {
        rOutput.resize( number_of_integration_pts, false );
    }

    // TODO: ADD CONTENT!!!!!

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact3DCondition::GetNodalDeltaMovements( 
    Vector& rValues,
    const int& rNode 
    )
{
    unsigned int dimension = 3;

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
    rValues[2] = CurrentValueVector[2] - PreviousValueVector[2];
}

//************************************************************************************
//************************************************************************************

Vector& MortarContact3DCondition::GetCurrentValue( 
    const Variable<array_1d<double,3> >& rVariable,
    Vector& rValue,
    const unsigned int& rNode
    )
{
    KRATOS_TRY;

    const unsigned int dimension = 3;

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

Vector& MortarContact3DCondition::GetPreviousValue( 
    const Variable<array_1d<double,3> >& rVariable,
    Vector& rValue,
    const unsigned int& rNode
    )
{
    KRATOS_TRY;

    const unsigned int dimension = 3;

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

void MortarContact3DCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
}

void MortarContact3DCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
}

} // Namespace Kratos
