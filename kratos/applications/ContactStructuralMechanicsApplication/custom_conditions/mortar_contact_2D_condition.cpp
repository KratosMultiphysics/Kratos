// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
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
#include "contact_structural_mechanics_application.h"
#include "contact_structural_mechanics_application_variables.h"

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
    mUseColocationIntegration = false;

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
        else if (this->GetValue(INTEGRATION_ORDER_CONTACT) == 11)
        {
            mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_5;
        }
        else
        {
//            std::cout << "The number of integration points is not defined.  INTEGRATION_ORDER_CONTACT: "<< this->GetValue(INTEGRATION_ORDER_CONTACT) << std::endl;
//            std::cout << "Options are: 1, 2, 3, 4, 5  " << std::endl;
//            std::cout << "Taking default number of integration points (INTEGRATION_ORDER_CONTACT = 1)  " << std::endl;
//            mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
//            std::cout << "The number of integration points is not defined for Gauss quadrature. Using colocation integration with "
//                << this->GetValue(INTEGRATION_ORDER_CONTACT) << " points." << std::endl;
            mUseColocationIntegration = true;
            mColocationIntegration.Initialize( this->GetValue(INTEGRATION_ORDER_CONTACT) );
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
    
    const double ActiveCheckFactor = GetProperties().GetValue(ACTIVE_CHECK_FACTOR);
    
    for ( unsigned int i_cond = 0; i_cond < all_containers->size(); ++i_cond )
    {
        mThisMasterElements[i_cond] = (*all_containers)[i_cond].condition;

        // Fill the condition
        Condition::Pointer& pCond = mThisMasterElements[i_cond];
        
        ContactUtilities::ContactContainerFiller((*all_containers)[i_cond], pCond->GetGeometry().Center(), GetGeometry(), pCond->GetGeometry(), 
                                                this->GetValue(NORMAL), pCond->GetValue(NORMAL), ActiveCheckFactor);
        
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    
/*
 * 
 * commented out because we should use the PDASS instead of the gap tolerance to update active and incactive during the iterations
 * Gap tolerance should be used only at the beginning of the solution step to add additional "close" nodes to the active set
 * the following lines have been merged with this->InitializeSolutionStep()
 * 
 */
    
//    KRATOS_TRY;
//    
//    // Populate the vector of master elements
//    std::vector<contact_container> *& all_containers = this->GetValue(CONTACT_CONTAINERS);
//    
//    double ActiveCheckFactor = 0.005;
//    if( GetProperties().Has(ACTIVE_CHECK_FACTOR) )
//    {
//        ActiveCheckFactor = GetProperties().GetValue(ACTIVE_CHECK_FACTOR);
//    }
//    
//    for ( unsigned int i_cond = 0; i_cond < all_containers->size(); ++i_cond )
//    {
//        Condition::Pointer pCond = (*all_containers)[i_cond].condition;
//    
//        // Fill the condition
//        ContactUtilities::ContactContainerFiller((*all_containers)[i_cond], pCond->GetGeometry().Center(), GetGeometry(), pCond->GetGeometry(), 
//                                                 this->GetValue(NORMAL), pCond->GetValue(NORMAL), ActiveCheckFactor);
//    }
//    
//    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    // TODO: Add things if needed
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
    
    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = mUseColocationIntegration ?
                                                                         mColocationIntegration.IntegrationPoints( ) :
                                                                         GetGeometry( ).IntegrationPoints( mThisIntegrationMethod );
    
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
        const unsigned int number_of_nodes_master = current_master_element.PointsNumber( );

        // Compute mortar condition matrices
        MortarConditionMatrices ThisMortarConditionMatrices;
        ThisMortarConditionMatrices.Initialize( number_of_nodes_master, number_of_nodes_slave, dimension );
        
        // Weighted gaps
        VectorType& gn = ThisMortarConditionMatrices.NodalWeightedGaps;
        
        // LOG_CONDITION_HEADER( Variables.GetMasterElement( ), GetGeometry( ) )
        
        // LOG_GENERAL( LT_CYAN, "   |_ L_X : ", GetGeometry( )[0].FastGetSolutionStepValue( VECTOR_LAGRANGE_MULTIPLIER_X ) << ", " << GetGeometry( )[1].FastGetSolutionStepValue( VECTOR_LAGRANGE_MULTIPLIER_X ) )

        double aux_integ_area = 0.0;
        
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            // // LOG_INFO( PointNumber )
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
                if (GetGeometry()[iNode].GetValue(WEIGHTED_GAP) > gn[iNode] )
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
                this->CalculateAndAddRHS( rLocalSystem, Variables, ThisMortarConditionMatrices );
            }
        }
        else
        {
//            std::cout << RED << "All integration pts projections are outside the master surface. Ignoring this condition !!" << RESET << std::endl;
        }
    }
    
//     std::cout << "--------------------------------------------------" << std::endl;
//     KRATOS_WATCH(this->Id())
//     MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix( );   
//     VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector( );   
//     KRATOS_WATCH(rLeftHandSideMatrix);
//     KRATOS_WATCH(rRightHandSideVector);

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

    /* LOCAL COORDINATES */
    const GeometryType::IntegrationPointsArrayType& integration_points = mUseColocationIntegration ?
                                                                         mColocationIntegration.IntegrationPoints( ) :
                                                                         GetGeometry( ).IntegrationPoints( mThisIntegrationMethod );

    const PointType& eta = integration_points[rPointNumber].Coordinates();
    Point<3> local_point;
    local_point.Coordinates() = ZeroVector(3);
    
    local_point.Coordinate(1) = eta[0];
    rVariables.SegmentProportion = 1.0;
    
//    if (mUseColocationIntegration == true)
//    {
//        local_point.Coordinate(1) = eta[0];
//        rVariables.SegmentProportion = 1.0;
//    }
//    else
//    {
//        local_point.Coordinate(1) = 0.5 * (1.0 - eta[0]) * current_container.local_coordinates_slave[0]
//                                  + 0.5 * (1.0 + eta[0]) * current_container.local_coordinates_slave[1];
//        rVariables.SegmentProportion = (current_container.local_coordinates_slave[1] - current_container.local_coordinates_slave[0])/2.0;
//    }
    
    /* RESIZE MATRICES AND VECTORS */
    rVariables.Phi_LagrangeMultipliers.resize( number_of_slave_nodes, false );
    rVariables.N_Master.resize( number_of_master_nodes, false );
    rVariables.N_Slave.resize( number_of_slave_nodes, false );

    rVariables.DN_De_Master.resize( number_of_master_nodes, master_nodes.LocalSpaceDimension( ), false );
    rVariables.DN_De_Slave.resize( number_of_slave_nodes, slave_nodes.LocalSpaceDimension( ), false );
    rVariables.DPhi_De_LagrangeMultipliers.resize( number_of_slave_nodes, slave_nodes.LocalSpaceDimension( ), false );
    
    /*  POPULATE MATRICES AND VECTORS */
    
    /// SLAVE CONDITION ///
    rVariables.N_Slave = slave_nodes.ShapeFunctionsValues( rVariables.N_Slave, local_point.Coordinates() );
    rVariables.DN_De_Slave  =  slave_nodes.ShapeFunctionsLocalGradients( rVariables.DN_De_Slave , local_point );

    /// LM SPACE CHOICE - STANDARD OR DUAL ///
    rVariables.Phi_LagrangeMultipliers = rVariables.N_Slave;
    rVariables.DPhi_De_LagrangeMultipliers = rVariables.DN_De_Slave;
//    rVariables.Phi_LagrangeMultipliers = LagrangeMultiplierShapeFunctionValue( local_point[0] );
//    rVariables.DPhi_De_LagrangeMultipliers = LagrangeMultiplierShapeFunctionLocalGradient( local_point[0] );
    
    /// MASTER CONDITION ///
    this->MasterShapeFunctionValue( rVariables, local_point);

    /* CALCULATE JACOBIAN AND JACOBIAN DETERMINANT */
    slave_nodes.Jacobian( rVariables.j_Slave, local_point.Coordinates() );
    rVariables.DetJSlave = ContactUtilities::ContactElementDetJacobian( rVariables.j_Slave );

    KRATOS_CATCH( "" );
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

void MortarContact2DCondition::MasterShapeFunctionValue(
    GeneralVariables& rVariables,
    const PointType& local_point 
    )
{
    GeometryType& master_seg = rVariables.GetMasterElement( );
    rVariables.N_Master     = ZeroVector( master_seg.PointsNumber( ) );
    rVariables.DN_De_Master = ZeroMatrix( master_seg.PointsNumber( ), master_seg.LocalSpaceDimension( ) );

    PointType projected_gp_global;
    rVariables.IntegrationPointNormalGap = 0.0;
    VectorType& normal = rVariables.IntegrationPointNormalVector;
    normal = ContactUtilities::GaussPointNormal(rVariables.N_Slave, GetGeometry());
    
    GeometryType::CoordinatesArrayType slave_gp_global;
    double aux_dist = 0.0;
    this->GetGeometry( ).GlobalCoordinates( slave_gp_global, local_point );
    ContactUtilities::ProjectDirection( master_seg, slave_gp_global, projected_gp_global, aux_dist, -normal ); // The opposite direction
    rVariables.IntegrationPointNormalGap = -aux_dist;
    
    GeometryType::CoordinatesArrayType projected_gp_local;
    if( master_seg.IsInside( projected_gp_global.Coordinates( ), projected_gp_local ) )
    {
        // SHAPE FUNCTIONS 
        rVariables.N_Master     = master_seg.ShapeFunctionsValues(         rVariables.N_Master,     projected_gp_local );         
        rVariables.DN_De_Master = master_seg.ShapeFunctionsLocalGradients( rVariables.DN_De_Master, projected_gp_local );
        master_seg.Jacobian( rVariables.j_Master, projected_gp_local );
        rVariables.ColocationWeightCoeff = 1;
    }
    else
    {
//        if( mUseColocationIntegration )
//        {
////            std::cout << RED << "Integration pt is outside. Coords: " << local_point.Coordinates( ) << RESET << std::endl;
//            rVariables.ColocationWeightCoeff = 0;
//        }
            rVariables.ColocationWeightCoeff = 0;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateDAndM( 
    GeneralVariables& rVariables,
    const double& rIntegrationWeight,
    MortarConditionMatrices& ThisMortarConditionMatrices
    )
{
    // LOG_GENERAL( LT_YELLOW, "|.... ", __FUNCTION__ );
    
    const unsigned int num_slave_nodes  = GetGeometry().PointsNumber( );
    const unsigned int num_master_nodes = rVariables.GetMasterElement( ).PointsNumber( );

    const Vector& N_s = rVariables.N_Slave;
    const Vector& N_m = rVariables.N_Master;
    const Vector& Phi = rVariables.Phi_LagrangeMultipliers;
    const double& gap = rVariables.IntegrationPointNormalGap;
    const double& J_s = rVariables.DetJSlave;
    const double& SegProp = rVariables.SegmentProportion;
    
    MatrixType& D  = ThisMortarConditionMatrices.D;
    MatrixType& M  = ThisMortarConditionMatrices.M;
    VectorType& gn  = ThisMortarConditionMatrices.NodalWeightedGaps;

    // For all the nodes
    for ( unsigned int i_slave = 0; i_slave < num_slave_nodes; ++i_slave )
    {
        gn[i_slave] += rIntegrationWeight * Phi( i_slave ) * gap * J_s * SegProp;

//        D( i_slave, i_slave ) += rIntegrationWeight * Phi( i_slave ) * N_s( i_slave ) * J_s * SegProp;
        for ( unsigned int j_slave = 0; j_slave < num_slave_nodes; ++j_slave )
        {
            D( i_slave, j_slave  ) += rIntegrationWeight * Phi( i_slave ) * N_s( j_slave ) * J_s;
        }

        for ( unsigned int j_master = 0; j_master < num_master_nodes; ++j_master )
        {
            M( i_slave, j_master ) += rIntegrationWeight * Phi( i_slave ) * N_m( j_master ) * J_s * SegProp;
        }
    }
}

 /***********************************************************************************/
 /***********************************************************************************/
 
void MortarContact2DCondition::CalculateDeltaDetJSlave(
    GeneralVariables& rVariables,
    MortarConditionMatrices& ThisMortarConditionMatrices
    )
{
    // LOG_GENERAL( LT_YELLOW, "|.... ", __FUNCTION__ );

    // prerequisites
    const unsigned int dimension = 2;
    const unsigned int num_slave_nodes  = GetGeometry( ).PointsNumber( );

    const double& det_J_s     = rVariables.DetJSlave;
    const VectorType& DN_De_s = column( rVariables.DN_De_Slave, 0 );
    const VectorType& J_gp    = column( rVariables.j_Slave, 0 );

    MatrixType& DeltaJ_s = ThisMortarConditionMatrices.DeltaDetJSlave; 

    // Fill up the elements corresponding to the slave DOFs - the rest remains zero
    for ( unsigned int i_slave = 0, i = 0; i_slave < num_slave_nodes; ++i_slave, i += dimension )
    {
        if ( GetGeometry( )[i_slave].Is( ACTIVE ) )
        {
            DeltaJ_s( 0, i     ) = J_gp( 0 ) * DN_De_s( i_slave ) / det_J_s;
            DeltaJ_s( 0, i + 1 ) = J_gp( 1 ) * DN_De_s( i_slave ) / det_J_s;
        }
    }
}
 
 /***********************************************************************************/
 /***********************************************************************************/
 
 void MortarContact2DCondition::CalculateDeltaPhi( 
     GeneralVariables& rVariables,
     const double& rIntegrationWeight,
     MortarConditionMatrices& ThisMortarConditionMatrices
     )
 {
//     // LOG_GENERAL( LT_YELLOW, "|.... ", __FUNCTION__ );
//
//     const unsigned int dimension = 2;
//     const unsigned int num_slave_nodes  = GetGeometry( ).PointsNumber( );
//     const unsigned int num_master_nodes = rVariables.GetMasterElement( ).PointsNumber( );
//
// 
//     const Vector& N_s      = rVariables.N_Slave;
//     const double& det_J_s  = rVariables.DetJSlave;
// 
//     const MatrixType& DeltaJ       = ThisMortarConditionMatrices.DeltaDetJSlave;
//     const MatrixType& DeltaXi_gp_s = ThisMortarConditionMatrices.DeltaIntegrationPoint[0];
//     MatrixType&       DeltaPhi     = ThisMortarConditionMatrices.DeltaPhiLagrangeMultipliers;
// 
//     MatrixType A_e    = ZeroMatrix( num_slave_nodes, num_slave_nodes );
//     MatrixType Minv_e =     Matrix( num_slave_nodes, num_slave_nodes );
//     this->CalculateA( A_e, Minv_e );
//
//     for ( unsigned int i = 0; i < num_slave_nodes; ++i )
//     {
//         for ( unsigned int j = 0; j < num_slave_nodes; ++j )
//         {
//             Matrix DeltaA_ij = ZeroMatrix( dimension, dimension );
//             
//             DeltaA_ij += N_s[i] * Minv_e(i, j) * subrange( DeltaJ,
//                                                            0,
//                                                            dimension,
//                                                            i * dimension,
//                                                            i * dimension + dimension );
//             
//             for ( unsigned int k = 0; k < num_slave_nodes; ++k )
//             {
//                 for ( unsigned int l = 0; l < num_slave_nodes; ++l )
//                 {
//                     DeltaA_ij += A_e(i,k) * N_s[k] * N_s[l] * subrange( DeltaJ,
//                                                                         0,
//                                                                         dimension,
//                                                                         l * dimension,
//                                                                         l * dimension + dimension );
//                 }
//             }
//             
//             subrange( DeltaPhi,
//                       i * dimension,
//                       i * dimension + dimension,
//                       j * dimension,
//                       j * dimension ) += DeltaA_ij * N_s[j] * rIntegrationWeight;
//         }
//     }
 }

 /***********************************************************************************/
 /***********************************************************************************/

 void MortarContact2DCondition::CalculateA( MatrixType& rA, MatrixType& rMinv_e )
 {
//     const unsigned int dimension = 2;
//     const unsigned int num_slave_nodes = GetGeometry( ).PointsNumber( );
//     const GeometryType::IntegrationPointsArrayType& integration_points = mUseColocationIntegration ?
//                                                                          mColocationIntegration.IntegrationPoints( ) :
//                                                                          GetGeometry( ).IntegrationPoints( mThisIntegrationMethod );
//
//     /*
//      * Calculate the modified D and M matrices - M_e and D_e ( refer to Popp's thesis eq. 4.57 )
//      * and Invert M_e
//      * Note: D_e and M_e are stored as non-blocks
//      * This will actually make life easier to calculate DeltaPhi
//      */
//     Matrix D_e = ZeroMatrix( num_slave_nodes, num_slave_nodes );
//     Matrix M_e = ZeroMatrix( num_slave_nodes, num_slave_nodes );
//     
//     for ( unsigned int i_point = 0; i_point < integration_points.size( ); ++i_point )
//     {
//         const double& integration_weight          = integration_points[i_point].Weight( );
//         const array_1d<double,3>& local_pt_coords = integration_points[i_point].Coordinates( );
//         MatrixType J_gp = ZeroMatrix( num_slave_nodes, 1 );
//         GetGeometry( ).Jacobian( J_gp, i_point, mThisIntegrationMethod );
//         const double& det_J_s = norm_2( column( J_gp , 0 ) ); 
//         
//         // TODO check if this is correct - not clear in Popp's thesis
//         for ( unsigned int i = 0; i < num_slave_nodes; i++ )
//         {
//             D_e( i, i ) += integration_weight * det_J_s * GetGeometry( ).ShapeFunctionValue( i, local_pt_coords );
//             
//             for ( unsigned int j = 0; j < num_slave_nodes; j++ )
//                 M_e( j, i ) += D_e( i, i ) * GetGeometry( ).ShapeFunctionValue( j, local_pt_coords );
//         }
//     }
// 
//     double det_M_e = 0.0;
//     MathUtils<double>::InvertMatrix( M_e, rMinv_e, det_M_e );
// 
//     // Calculate matrix A = D * Minv
//     rA = prod( D_e, rMinv_e );
 }
 
 /***********************************************************************************/
 /***********************************************************************************/

 void MortarContact2DCondition::CalculateDeltaIntegrationSegmentPoint(
     GeneralVariables& rVariables,
     MortarConditionMatrices& ThisMortarConditionMatrices
     )
 {
     // LOG_GENERAL( LT_YELLOW, "|.... ", __FUNCTION__ );

     // Geometries
     const unsigned int& dimension = 2;
     const GeometryType& slave_nodes  = GetGeometry( );
     const GeometryType& master_nodes = rVariables.GetMasterElement( );
     const unsigned int& num_slave_nodes  =  slave_nodes.PointsNumber( );
     const unsigned int& num_master_nodes = master_nodes.PointsNumber( );
     const unsigned int num_total_nodes = num_master_nodes + num_slave_nodes;
     
     // Shape fucntions and their dertivatives
     const VectorType& N_m      = rVariables.N_Master;
     const VectorType& N_s      = rVariables.N_Slave;
     const VectorType& DN_Dxi_m = column( rVariables.DN_De_Master, 0 );  // DN_master/Dxi
     const array_1d<double, 3>& ng        = rVariables.IntegrationPointNormalVector;
     const array_1d<double, 3>& x2hat_x1  = rVariables.IntegrationPointNormalGap * ng;

     // Integration Point Directional Derivatives
     MatrixType& DeltaXi_gp_2 = ThisMortarConditionMatrices.DeltaMasterIntegrationPoint;
     
     /*
      * Master side integration point (the projection of the slave integration point)                             
      * Even using collocation, this won't evaluate to zero because of the linearization of the projection operator
      * Note: using collocation, this is calculated only when the slave integration point is contact the mastet segment
      * ( i.e. when  the weight coefficient is not zero )
      */
     
     double denom = 0.0;
     for ( unsigned int i_master = 0; i_master < num_master_nodes; ++i_master )
     {
         denom += DN_Dxi_m( i_master ) * ( master_nodes[i_master].Y( ) * ng[0] - master_nodes[i_master].X( ) * ng[1] );
     }

     // calculating Dn_g
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
     const MatrixType norm_Delta_ng = prod( I - outer_prod( ng, ng ), Delta_ng ) / rVariables.DetJSlave;
     Delta_ng = norm_Delta_ng;
     
     // [ -( y2hat - y1 ) * Delta_ng_x,   ( x2hat - x1 ) * Delta_ng_y ]' / denom  
     row( Delta_ng, 0 ) *= -x2hat_x1[1] / denom;
     row( Delta_ng, 1 ) *=  x2hat_x1[0] / denom;

     // now ready for addition to DeltaXi_gp
     subrange( DeltaXi_gp_2,
         0, 
         dimension, 
         dimension * num_master_nodes, 
         dimension * num_total_nodes ) += Delta_ng;

     // add master shape functions contributions
     // [ N_m * ngy / denom,    0
     //   0,                    -N_m * ngx / denom ] 
     for ( unsigned int i_master = 0, j = 0; i_master < num_master_nodes; ++i_master, j += dimension )
     {
         DeltaXi_gp_2( 0, j     ) += N_m(i_master) * ng[1] / denom;
         DeltaXi_gp_2( 1, j + 1 ) -= N_m(i_master) * ng[0] / denom;
     }

     // add slave shape functions contributions
     // [ -N_s * ngy / denom,   0
     //   0,                    N_s * ngx / denom ] 
     for ( unsigned int i_slave = 0, j = num_master_nodes * dimension; i_slave < num_slave_nodes; ++i_slave, j += dimension )
     {
         DeltaXi_gp_2( 0, j     ) -= N_s(i_slave) * ng[1] / denom;
         DeltaXi_gp_2( 1, j + 1 ) += N_s(i_slave) * ng[0] / denom;
     }
 }

 /***********************************************************************************/
 /***********************************************************************************/

 void MortarContact2DCondition::CalculateLambdaTDeltaBco( 
     GeneralVariables& rVariables,
     const double& rIntegrationWeight,
     MortarConditionMatrices& ThisMortarConditionMatrices
     )
 {
     // LOG_GENERAL( LT_YELLOW, "|.... ", __FUNCTION__ );

     const unsigned int dimension = 2;
     const unsigned int num_slave_nodes  = GetGeometry( ).PointsNumber( );
     const unsigned int num_master_nodes = rVariables.GetMasterElement( ).PointsNumber( );
     const unsigned int master_dofs      = num_master_nodes * dimension;
 
     const Vector& N_s      = rVariables.N_Slave;
     const Vector& N_m      = rVariables.N_Master;
     const Vector& Phi      = rVariables.Phi_LagrangeMultipliers;
     const Vector& DN_Dxi_m = column( rVariables.DN_De_Master, 0 );  // DN_master/Dxi
     const double& det_J_s  = rVariables.DetJSlave;
 
     const MatrixType& DeltaDetJ    = ThisMortarConditionMatrices.DeltaDetJSlave;
     const MatrixType& DeltaXi_gp_2 = ThisMortarConditionMatrices.DeltaMasterIntegrationPoint;
     
     MatrixType& DeltaD = ThisMortarConditionMatrices.LambdaTDeltaD;
     MatrixType& DeltaM = ThisMortarConditionMatrices.LambdaTDeltaM;
     
     VectorType lambda = ZeroVector(dimension);
     for ( unsigned int i_slave = 0; i_slave < num_slave_nodes; ++i_slave )
     {
         lambda[0] += Phi[i_slave] * GetGeometry( )[i_slave].FastGetSolutionStepValue( VECTOR_LAGRANGE_MULTIPLIER_X );
         lambda[1] += Phi[i_slave] * GetGeometry( )[i_slave].FastGetSolutionStepValue( VECTOR_LAGRANGE_MULTIPLIER_Y );
     }
     
     // DeltaD' * lambda
     for ( unsigned int k_slave = 0,  k = 0; k_slave  < num_slave_nodes;  ++k_slave,  k += dimension )
     {
         row( DeltaD, k     ) += rIntegrationWeight * N_s[k_slave] * lambda[0] * row( DeltaDetJ, 0 );
         row( DeltaD, k + 1 ) += rIntegrationWeight * N_s[k_slave] * lambda[1] * row( DeltaDetJ, 0 );
     }

     // DeltaM' * lambda
     for ( unsigned int l_master = 0, l = 0; l_master < num_master_nodes; ++l_master, l += dimension )
     {
         subrange( DeltaM, l,     l + 1, master_dofs, DeltaM.size2( ) ) += rIntegrationWeight * lambda[0] * DeltaDetJ            *      N_m[l_master];
         subrange( DeltaM, l + 1, l + 2, master_dofs, DeltaM.size2( ) ) += rIntegrationWeight * lambda[1] * DeltaDetJ            *      N_m[l_master];
              row( DeltaM, l                                          ) += rIntegrationWeight * lambda[0] * row(DeltaXi_gp_2, 0) * DN_Dxi_m[l_master] * det_J_s;
              row( DeltaM, l + 1                                      ) += rIntegrationWeight * lambda[1] * row(DeltaXi_gp_2, 1) * DN_Dxi_m[l_master] * det_J_s;
     }
 }

 /***********************************************************************************/
 /***********************************************************************************/

 void MortarContact2DCondition::CalculateN(
     GeneralVariables& rVariables,
     const double& rIntegrationWeight,
     MortarConditionMatrices& ThisMortarConditionMatrices
     )
 {
     // LOG_GENERAL( LT_YELLOW, "|.... ", __FUNCTION__ );

     // Contact pair variables
     const unsigned int& dimension = 2;
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

 void MortarContact2DCondition::CalculateDeltaDiscreteGap( 
          GeneralVariables& rVariables,
          MortarConditionMatrices& ThisMortarConditionMatrices,
          MatrixType& rDeltaDiscreteGap
          )
 {
     // LOG_GENERAL( YELLOW, "|........ ", __FUNCTION__ );


     // Geometries
     const unsigned int& dimension = 2;
     const GeometryType& slave_nodes  = GetGeometry( );
     const GeometryType& master_nodes = rVariables.GetMasterElement( );
     const unsigned int& num_slave_nodes  =  slave_nodes.PointsNumber( );
     const unsigned int& num_master_nodes = master_nodes.PointsNumber( );
     const unsigned int& num_total_nodes = num_slave_nodes + num_master_nodes;
     
     // Shape fucntions and their dertivatives
     const VectorType& N_s      = rVariables.N_Slave;
     const VectorType& N_m      = rVariables.N_Master;
     const VectorType& Jxi_m    = column( rVariables.j_Master, 0 );
     const array_1d<double, 3>& ng       = rVariables.IntegrationPointNormalVector;
     const array_1d<double, 3>& x2hat_x1 = rVariables.IntegrationPointNormalGap * ng;
     
     // Directional derivative of Xi2 - needed even for collocation
     const MatrixType& DeltaXi2 = ThisMortarConditionMatrices.DeltaMasterIntegrationPoint;
     const MatrixType& DeltaXi2_x = subrange( DeltaXi2, 0, 1, 0, dimension * num_total_nodes );
     const MatrixType& DeltaXi2_y = subrange( DeltaXi2, 1, 2, 0, dimension * num_total_nodes );
     
     // ng' . ( N2 . Dx2 )
     for ( unsigned int j_master = 0, j = 0; j_master < num_master_nodes; ++j_master, j += dimension )
     {
         rDeltaDiscreteGap( 0, j     ) += ng[0] * N_m(j_master); 
         rDeltaDiscreteGap( 0, j + 1 ) += ng[1] * N_m(j_master);
     }
     
     // - ng' . ( N1 . Dx1 )
     for ( unsigned int j_slave = 0, j = num_master_nodes * dimension; j_slave < num_slave_nodes; ++j_slave, j += dimension )
     {
         if ( slave_nodes[j_slave].Is( ACTIVE ) )
         {
         rDeltaDiscreteGap( 0, j     ) -= ng[0] * N_s(j_slave); 
         rDeltaDiscreteGap( 0, j + 1 ) -= ng[1] * N_s(j_slave);
         }
     }
     
     // ( x2hat - x1 )'. DNg 
     const MatrixType I = IdentityMatrix( dimension, dimension );
     const MatrixType Delta_ng_norm = ( I - outer_prod( ng, ng ) ) / rVariables.DetJSlave;
     for ( unsigned int j_slave = 0, j = num_master_nodes * dimension; j_slave < num_slave_nodes; ++j_slave, j += dimension )
     {
         const MatrixType Dn_j = prod( Delta_ng_norm, slave_nodes[j_slave].GetValue( DELTA_NORMAL ) );
         rDeltaDiscreteGap( 0, j     ) += N_s(j_slave) * ( x2hat_x1[0] * Dn_j( 0, 0 ) + x2hat_x1[1] * Dn_j( 1, 0 ) );
         rDeltaDiscreteGap( 0, j + 1 ) += N_s(j_slave) * ( x2hat_x1[0] * Dn_j( 0, 1 ) + x2hat_x1[1] * Dn_j( 1, 1 ) );
     }
     
     // Directional derivative of master segment
//     MatrixType DeltaNm_contribution = ZeroMatrix( 1, dimension * num_total_nodes );
//     for ( unsigned int i_master = 0; i_master < num_master_nodes; ++i_master )
//     {
//         DeltaNm_contribution += DN_Dxi_m[i_master] * ( ng[0] * master_nodes[i_master].X( ) * DeltaXi2_x + ng[1] * master_nodes[i_master].Y( ) * DeltaXi2_y ); 
//     }
//     
//     rDeltaDiscreteGap += DeltaNm_contribution;
     rDeltaDiscreteGap += ( Jxi_m[0] * ng[0] * DeltaXi2_x + Jxi_m[1] * ng[1] * DeltaXi2_y );
 }

/***********************************************************************************/
/***********************************************************************************/

 void MortarContact2DCondition::CalculateF( MatrixType& rF )
 {
     // In 2D this only calculates Fs_xi

     // LOG_GENERAL( LT_YELLOW, "|.. ", __FUNCTION__ );

     const unsigned int dimension = 2;
     const unsigned int num_slave_nodes  = GetGeometry().PointsNumber( );

     /*                 
      *  DT = DN x e3 = | N12,  -N11 | 
      *                 | N22,  -N21 | 
      */

     unsigned int i, j = 0;
     for (unsigned int iSlave = 0; iSlave < num_slave_nodes; iSlave++)
     {
         if( GetGeometry( )[iSlave].Is( ACTIVE ) )
         {
             i = iSlave;
             j = iSlave * dimension;

             const array_1d<double, 3>& LM = GetGeometry( )[iSlave].FastGetSolutionStepValue( VECTOR_LAGRANGE_MULTIPLIER, 0 );
             const MatrixType& DeltaNormal = GetGeometry( )[iSlave].GetValue(DELTA_NORMAL);
             MatrixType DeltaTangentXi = MatrixType( dimension, dimension );
             column( DeltaTangentXi, 0 ) =  column( DeltaNormal, 1 );
             column( DeltaTangentXi, 1 ) = -column( DeltaNormal, 0 );
             
             rF( i, j     ) += DeltaTangentXi(0, 0) * LM[0] + DeltaTangentXi(0, 1) * LM[1];
             rF( i, j + 1 ) += DeltaTangentXi(1, 0) * LM[0] + DeltaTangentXi(1, 1) * LM[1];
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
    
   // LOG_GENERAL( LT_MAGENTA << BOLD, "\n", __FUNCTION__ );
    
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
        this->CalculateAndAddMortarContactOperator(  LHS_contact_pair, rVariables, ThisMortarConditionMatrices );
        this->CalculateAndAddContactStiffnessMatrix( LHS_contact_pair, rVariables, ThisMortarConditionMatrices );
        this->CalculateAndAddNormalGapLinearization( LHS_contact_pair, rVariables, ThisMortarConditionMatrices );
        
      // LOG_MATRIX_PRETTY( LHS_contact_pair )
        
//        /********** DEBUG **********/
//        std::cout << CYAN;
//        std::cout.precision(3);
//        std::cout << "LHS_contact_pair " << "[ " << LHS_contact_pair.size1( ) << " x " << LHS_contact_pair.size2( ) << " ] :" << std::endl;
//        for ( unsigned int i = 0; i < LHS_contact_pair.size1( ); ++i )
//        {
//            for ( unsigned int j = 0; j < LHS_contact_pair.size2( ); ++j )
//                std::cout << LHS_contact_pair(i,j) << "\t";
//            std::cout << std::endl;
//        }
//        std::cout << std::endl;
//        std::cout << RESET;
//        /********** DEBUG **********/
        
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
    // LOG_GENERAL( MAGENTA, "|.... ", __FUNCTION__ );

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
    const MortarConditionMatrices& ThisMortarConditionMatrices
    )
{
   // LOG_GENERAL( LT_MAGENTA << BOLD, "\n", __FUNCTION__ );
    
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
        /* SINGLE LHS MATRIX */
        VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();
        const unsigned int pair_size = 2 * ( rVariables.GetMasterElement( ).PointsNumber( ) + 2 * GetGeometry( ).PointsNumber( ) ); 
        VectorType RHS_contact_pair = ZeroVector(pair_size);
        
        // Calculate
        this->CalculateAndAddMortarContactOperator( RHS_contact_pair, rVariables, ThisMortarConditionMatrices );
        this->CalculateAndAddGap(                   RHS_contact_pair, rVariables, ThisMortarConditionMatrices );

      // LOG_VECTOR_PRETTY( RHS_contact_pair )
        
        // Assemble
        this->AssembleContactPairRHSToConditionSystem( rVariables.GetMasterElementIndex( ), RHS_contact_pair, rRightHandSideVector );
    }
    
}
  
/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::AssembleContactPairRHSToConditionSystem(
    const unsigned int rPairIndex,
    VectorType& rPairRHS,
    VectorType& rConditionRHS 
    )
{
    // LOG_GENERAL( MAGENTA, "|.... ", __FUNCTION__ );

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
    // LOG_GENERAL( MAGENTA, "|.... ", __FUNCTION__ );
    
    KRATOS_TRY;
  
    // Contact pair variables
    const unsigned int dimension = 2;
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
                // uncomment for mesh tying
//                rLeftHandSideMatrix( i,     j     ) += minus_M_ij;
//                rLeftHandSideMatrix( i + 1, j + 1 ) += minus_M_ij;
                                                    
                rLeftHandSideMatrix( j    , i     ) += minus_M_ij;
                rLeftHandSideMatrix( j + 1, i + 1 ) += minus_M_ij;
            }
      
            // Fill the D and D' parts
            for ( unsigned int j_slave = 0; j_slave < num_slave_nodes; ++j_slave )
            {
                j = ( j_slave + num_master_nodes ) * dimension;
                double D_ii = D( i_slave, j_slave );
                // uncomment for mesh tying
//                rLeftHandSideMatrix( i,     j     ) += D_ii;
//                rLeftHandSideMatrix( i + 1, j + 1 ) += D_ii;

                rLeftHandSideMatrix( j,     i     ) += D_ii;
                rLeftHandSideMatrix( j + 1, i + 1 ) += D_ii;
            }
        }
        else
        {
            // We impose a 0 zero LM in the inactive nodes
            rLeftHandSideMatrix( i    , i )    = 1.0;
            rLeftHandSideMatrix( i + 1, i + 1) = 1.0;
        }
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateAndAddMortarContactOperator( 
    VectorType& rRightHandSideVector,
    GeneralVariables& rVariables,
    const MortarConditionMatrices& ThisMortarConditionMatrices
    )
{
    // LOG_GENERAL( MAGENTA, "|.... ", __FUNCTION__ );

    KRATOS_TRY;
  
    // Contact pair variables
    const unsigned int dimension = 2;
    const unsigned int num_slave_nodes  = GetGeometry( ).PointsNumber( );
    const unsigned int num_master_nodes = rVariables.GetMasterElement( ).PointsNumber( );  

    const MatrixType& D   = ThisMortarConditionMatrices.D;
    const MatrixType& M   = ThisMortarConditionMatrices.M;
    
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
        }

        // Fill the - lambda *  D part
        for ( unsigned int j_slave = 0; j_slave < num_slave_nodes; ++j_slave )
        {
            j = ( num_master_nodes + j_slave  ) * dimension;
            const double D_ii = D( i_slave, j_slave );
            rRightHandSideVector[ j     ] -= lagrange_multiplier[0] * D_ii;
            rRightHandSideVector[ j + 1 ] -= lagrange_multiplier[1] * D_ii;
        }
            
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateAndAddGap( 
    VectorType& rRightHandSideVector,
    GeneralVariables& rVariables,
    const MortarConditionMatrices& ThisMortarConditionMatrices
    )
{
    // LOG_GENERAL( MAGENTA, "|.... ", __FUNCTION__ );

    KRATOS_TRY;
  
    // Contact pair variables
    const unsigned int dimension = 2;
    const unsigned int num_slave_nodes  = GetGeometry( ).PointsNumber( );
    const unsigned int num_master_nodes = rVariables.GetMasterElement( ).PointsNumber( );  
    const unsigned int num_total_nodes  = num_slave_nodes + num_master_nodes;

    const VectorType& gn  = ThisMortarConditionMatrices.NodalWeightedGaps;
    
    unsigned int j = 0;
    array_1d<double, 3> lagrange_multiplier = ZeroVector(3); 
    for ( unsigned int i_slave = 0; i_slave < num_slave_nodes; ++i_slave )
    {
        noalias(lagrange_multiplier) = GetGeometry( )[i_slave].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER, 0); 
            
        if (GetGeometry( )[i_slave].Is(ACTIVE) == true)
        {
            // Adding the gap to the RHS
            j = (num_total_nodes + i_slave) * dimension;
            
            const array_1d<double, 3>& normal = GetGeometry( )[i_slave].GetValue( NORMAL );
            // for mesh tying
//            const array_1d<double, 3> gap_decomp = gn[i_slave] * normal; 
//            rRightHandSideVector[ j     ] -= gap_decomp[0]; 
//            rRightHandSideVector[ j + 1 ] -= gap_decomp[1];
            
            // for unilateral contact
            rRightHandSideVector[ j     ] -= gn[i_slave];
            rRightHandSideVector[ j + 1 ] -= ( +normal[1]*lagrange_multiplier[0] - normal[0]*lagrange_multiplier[1] ); // tx*LMx + ty*LMy 
        }
        else
        {
            // Adding the zero equality to the RHS
            j = (num_total_nodes + i_slave) * dimension; 
            rRightHandSideVector[ j     ] = 0.0;
            rRightHandSideVector[ j + 1 ] = 0.0;
        }
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateAndAddContactStiffnessMatrix(
        MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const MortarConditionMatrices& ThisMortarConditionMatrices)
{
    // LOG_GENERAL( MAGENTA, "|.... ", __FUNCTION__ );
    
    KRATOS_TRY
    
    // Contact pair variables
    const unsigned int dimension = 2;
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

void MortarContact2DCondition::CalculateAndAddNormalGapLinearization(
    MatrixType& rLeftHandSideMatrix,
    GeneralVariables& rVariables,
    const MortarConditionMatrices& ThisMortarConditionMatrices )
{
    // LOG_GENERAL( MAGENTA, "|.... ", __FUNCTION__ );
    
    KRATOS_TRY
    
    // Contact pair variables
    const unsigned int dimension = 2;
    const unsigned int num_slave_nodes  = GetGeometry( ).PointsNumber( );
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
    
    // Assemble F block - in 2D only Fsxi 
    MatrixType F = ZeroMatrix( num_slave_nodes, num_slave_nodes * dimension );
    this->CalculateF( F );
    subslice( rLeftHandSideMatrix,
              dimension * num_total_nodes + 1,
              dimension,
              F.size1( ),
              dimension * num_master_nodes,
              1,
              F.size2( ) ) += F; 

    
    // Assemble T Block - in 2D only T_xi; T_xi = [ +ny -nx ]
    for ( unsigned int i_slave = 0, i = num_total_nodes * dimension; i_slave < num_slave_nodes; ++i_slave, i+=dimension )
    {
        if (GetGeometry( )[i_slave].Is(ACTIVE) == true)
        {
            const array_1d<double, 3>& normal = GetGeometry( )[i_slave].GetValue( NORMAL );

            rLeftHandSideMatrix( i + 1, i     ) = +normal[1];
            rLeftHandSideMatrix( i + 1, i + 1 ) = -normal[0];
        }
    }
    
    KRATOS_CATCH( "" )
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
        }
        
        // Slave Nodes Displacement Equation IDs
        for ( unsigned int i_slave = 0; i_slave < num_slave_nodes; ++i_slave )
        {
            NodeType& slave_node = GetGeometry()[ i_slave ];
            rResult.push_back( slave_node.GetDof( DISPLACEMENT_X ).EquationId( ) );
            rResult.push_back( slave_node.GetDof( DISPLACEMENT_Y ).EquationId( ) );
        }
        
        // Slave Nodes  Lambda Equation IDs
        for ( unsigned int i_slave = 0; i_slave < num_slave_nodes; ++i_slave )
        {
            NodeType& slave_node = GetGeometry()[ i_slave ];
            rResult.push_back( slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_X ).EquationId( ) );
            rResult.push_back( slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Y ).EquationId( ) );
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
        }
        
        // Slave Nodes Displacement Equation IDs
        for ( unsigned int i_slave = 0; i_slave < num_slave_nodes; ++i_slave )
        {
            NodeType& slave_node = GetGeometry()[ i_slave ];
            rConditionalDofList.push_back( slave_node.pGetDof( DISPLACEMENT_X ) );
            rConditionalDofList.push_back( slave_node.pGetDof( DISPLACEMENT_Y ) );
        }
        
        // Slave Nodes Lambda Equation IDs
        for ( unsigned int i_slave = 0; i_slave < num_slave_nodes; ++i_slave )
        {
            NodeType& slave_node = GetGeometry()[ i_slave ];
            rConditionalDofList.push_back( slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_X ) );
            rConditionalDofList.push_back( slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Y ) );
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
