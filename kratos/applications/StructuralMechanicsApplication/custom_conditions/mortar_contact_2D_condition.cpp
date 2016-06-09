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
    // NOTE: Maybe is a good idea to add the tree_contact_search initialization here instead of the calling from python

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    // TODO Add content

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    // TODO Add content
    // TODO GetMasterElements

    this->DetermineActiveNodes( );

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::DetermineActiveNodes( )
{
    const unsigned int num_slave_nodes = GetGeometry( ).PointsNumber( );

    mThisInactiveSlaveNodes.resize( 0 );
    mThisActiveSlaveNodes.resize( 0 );

    for ( unsigned int iSlave = 0; iSlave < num_slave_nodes; ++iSlave )
    {
        if( GetGeometry( )[iSlave].Is( ACTIVE ) )
        {
            mThisActiveSlaveNodes.push_back( iSlave );
        }
        else
        {
            mThisInactiveSlaveNodes.push_back( iSlave );
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

double MortarContact2DCondition::LagrangeMultiplierShapeFunctionValue( 
    const IndexType& rPointNumber,
    const IndexType& rShapeFunctionIndex 
    )
{
    // NOTE: For more information look at the thesis of Popp page 93/236
    
    const GeometryType::IntegrationPointType& pt = GetGeometry( ).IntegrationPoints( mThisIntegrationMethod )[rPointNumber];
    GeometryType::CoordinatesArrayType local_coordinates;

    double eta = GetGeometry( ).PointLocalCoordinates( local_coordinates, pt.Coordinates( ) )( 0 );
    double phi = 0;

    if (rShapeFunctionIndex == 1 )
    {
        phi = ( 0.5 * ( 1 - 3 * eta ) );
    }
    else if (rShapeFunctionIndex == 2 )
    {
        phi = ( 0.5 * ( 1 + 3 * eta ) );
    }
    else
    {
        KRATOS_THROW_ERROR( std::logic_error, " The rShapeFunctionIndex is wrong: ", rShapeFunctionIndex );
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

    DPhi_De( 0, 0 ) = -3.0 / 2;
    DPhi_De( 0, 1 ) = +3.0 / 2;

    return DPhi_De;
}

/***********************************************************************************/
/***********************************************************************************/

const unsigned int MortarContact2DCondition::CalculateTotalNumberOfMasterNodes( )
{
    unsigned int num_master_nodes_total = 0;
    for ( unsigned int i_master_elem = 0; i_master_elem < mThisMasterElements.size( ); ++i_master_elem )
    {
        num_master_nodes_total += mThisMasterElements[i_master_elem]->PointsNumber( );
    }

    return num_master_nodes_total;
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

    for( unsigned int i=0; i<rLeftHandSideMatrices.size( ); i++ )
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

    const unsigned int condition_size = dimension * ( num_slave_nodes + num_master_nodes_total );

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
    GeometryType& current_slave_element = this->GetGeometry( );
    const unsigned int number_of_nodes_slave = current_slave_element.PointsNumber( );
    const unsigned int local_dimension_slave = current_slave_element.LocalSpaceDimension( );

    // Master segment info
    GeometryType& current_master_element = *( this->mThisMasterElements[rMasterElementIndex] );
    const unsigned int number_of_nodes_master = current_master_element.PointsNumber( );
    const unsigned int local_dimension_master = current_master_element.LocalSpaceDimension( );

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
    this->CalculateAndAppendMortarProjectionMatrixTranspose( rVariables );
    this->CalculateDeltaDAndDeltaM( rVariables, rIntegrationWeight, rCalculateSlaveContributions );
    this->CalculateDeltaJSlave( rVariables, rIntegrationWeight, rPointNumber );
    this->CalculateDeltaPhi( rVariables, rIntegrationWeight );
    // TODO: CalculateDiscreteGap goes here
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
    const unsigned int num_slave_nodes  = GetGeometry( ).PointsNumber( );
    const unsigned int num_master_nodes = rVariables.GetMasterElement( ).PointsNumber( );

    const Vector& N_s = rVariables.N_Slave;
    const Vector& N_m = rVariables.N_Master;
    const Vector& Phi = rVariables.Phi_LagrangeMultipliers;
    const double& J_s = rVariables.DetJSlave;

    MatrixType& D = mThisMortarConditionMatrices.D;
    MatrixType& M = mThisMortarConditionMatrices.M;

    unsigned int i = 0, j = 0;
    for ( unsigned int iNode = 0; iNode < num_slave_nodes; ++iNode )
    {
        i = iNode * dimension;
        // To avoid recalculating the D matrix
        if( rCalculateSlaveContributions )
        {
            D( i, i ) = rIntegrationWeight * N_s( iNode ) * J_s;
            D( i + 1, i + 1 ) = D( iNode, iNode );
        }

        for ( unsigned int jNode = 0; jNode < num_master_nodes; ++jNode )
        {
            j = jNode * dimension;
            M( i    , j      ) = rIntegrationWeight * Phi( iNode ) * N_m( jNode ) * J_s;
            M( i + 1,  j + 1 ) = M( i, j );
        }
    }

}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateDeltaDAndDeltaM( 
    GeneralVariables& rVariables,
    double& rIntegrationWeight,
    const bool& rCalculateSlaveContributions 
    )
{
    const unsigned int dimension = 2;
    const unsigned int num_slave_nodes  = GetGeometry( ).PointsNumber( );
    const unsigned int num_master_nodes = rVariables.GetMasterElement( ).PointsNumber( );

    const Vector& N_s = rVariables.N_Slave;
    const Vector& N_m = rVariables.N_Master;
    const Vector& Phi = rVariables.Phi_LagrangeMultipliers;

    MatrixType& DeltaD   = mThisMortarConditionMatrices.DeltaD;
    MatrixType& DeltaM   = mThisMortarConditionMatrices.DeltaM;
    MatrixType& DeltaJ_s = mThisMortarConditionMatrices.DeltaJSlave;

    // Calculate directional derivatives of mortar condition matrices
    for ( unsigned int iNode = 0; iNode < num_slave_nodes; ++iNode )
    {
        Matrix DeltaJ_iNode = subrange( DeltaJ_s, iNode, iNode + dimension, 0, num_slave_nodes * dimension );
        Matrix DeltaD_iNode = subrange( DeltaD,   iNode, iNode + dimension, 0, num_slave_nodes * dimension );
        Matrix DeltaM_iNode = subrange( DeltaM,   iNode, iNode + dimension, 0, num_slave_nodes * dimension );

        if( rCalculateSlaveContributions )
        {
            DeltaD_iNode = rIntegrationWeight * N_s( iNode ) * DeltaJ_iNode;
        }
        for ( unsigned int jNode = 0; jNode < num_master_nodes; ++jNode )
        {
            DeltaM_iNode = rIntegrationWeight * Phi( iNode ) * N_m( jNode ) * DeltaJ_iNode;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateDeltaJSlave( 
    GeneralVariables& rVariables,
    double& rIntegrationWeight,
    const IndexType& rPointNumber 
    )
{
    // Prerequisites
    const unsigned int dimension = 2;
    const unsigned int num_active_nodes = mThisActiveSlaveNodes.size( );

    const double& J_s         = rVariables.DetJSlave;
    const MatrixType& DN_De_s = rVariables.DN_De_Slave;
    const MatrixType& J_gp    = rVariables.j_Slave[rPointNumber];

    MatrixType& DeltaJ_s = mThisMortarConditionMatrices.DeltaJSlave;

    // Fill up the elements corresponding to the slave DOFs - the rest remains zero
    unsigned int i = 0, iNode = 0;
    for ( unsigned int iActive = 0; iActive < num_active_nodes; ++iActive )
    {
        iNode = mThisActiveSlaveNodes[iActive];
        i = iNode * dimension;

        DeltaJ_s( iNode, i     ) = J_gp( 0, 0 ) * ( 1 / J_s ) * DN_De_s( iNode, 0 );
        DeltaJ_s( iNode, i + 1 ) = J_gp( 1, 0 ) * ( 1 / J_s ) * DN_De_s( iNode, 0 );
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateDeltaPhi( 
    GeneralVariables& rVariables,
    double& rIntegrationWeight 
    )
{
    const unsigned int dimension = 2;
    const unsigned int num_slave_nodes  = GetGeometry( ).PointsNumber( );

    const Vector& N_s = rVariables.N_Slave;
    const double& J_s = rVariables.DetJSlave;

    MatrixType& D = mThisMortarConditionMatrices.D;
    MatrixType& M = mThisMortarConditionMatrices.M;

    MatrixType& DeltaD = mThisMortarConditionMatrices.DeltaD;
    MatrixType& DeltaM = mThisMortarConditionMatrices.DeltaM;
    MatrixType& DeltaPhi = mThisMortarConditionMatrices.DeltaPhiLagrangeMultipliers;

   /*
    * Invert M matrix - this is a modified version of Mortar's M matrix
    * refer to Popp's thesis eq. 4.57
    */
    // First calculate the modified M in reduced size to make the inversion explicit
    Matrix M_reduced = Matrix( num_slave_nodes, num_slave_nodes );
    Matrix Minv_reduced = Matrix( num_slave_nodes, num_slave_nodes );
    for ( unsigned int i = 0; i < num_slave_nodes; ++i )
    {
        for ( unsigned int j = 0; j < num_slave_nodes; ++j )
        {
            M_reduced( i, j ) = rIntegrationWeight * N_s( i ) * N_s( j ) * J_s;
        }
    }

    double det_M_reduced = MathUtils<double>::Det( M_reduced );
    MathUtils<double>::InvertMatrix( M_reduced, Minv_reduced, det_M_reduced );

    // Then fill up the blocks again to maintain sizes consistency
    MatrixType Minv = ZeroMatrix( M.size1( ), M.size2( ) );
    unsigned int i = 0, j = 0;
    for ( unsigned int iBlock = 0; iBlock < num_slave_nodes; ++iBlock )
    {
        for ( unsigned int jBlock = 0; jBlock < num_slave_nodes; ++jBlock )
        {
            i = iBlock * dimension;
            j = jBlock * dimension;
            Minv( i    , j     ) = Minv_reduced( iBlock, jBlock );
            Minv( i + 1, j + 1 ) = Minv( i, j );
        }
    }

    // Calculate matrix A = D * Minv
    MatrixType A = ZeroMatrix( num_slave_nodes * dimension, num_slave_nodes * dimension );

    for ( unsigned int i = 0; i < num_slave_nodes; i += dimension )
    {
        for ( unsigned int j = 0; j < num_slave_nodes; j += dimension )
        {
            A( i, j ) = D( i, i ) * Minv ( i, j );
            A( i+1, j+1 ) = A( i, j );
        }
    }

    // Finally calculate DeltaPhi
    DeltaPhi = prod( ( DeltaD + prod( A, DeltaM ) ), Minv );
    for ( unsigned int i = 0; i < num_slave_nodes; ++i )
    {
        Matrix DeltaPhi_i = subrange( DeltaPhi, i, i + dimension, 0, num_slave_nodes * dimension );
        DeltaPhi_i *= N_s( i );
    }
}


/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateNormalGapAtIntegrationPoint( 
    GeneralVariables& rVariables,
    const IndexType& rPointNumber
    )
{
    const unsigned int dimension = 2;
    const unsigned int num_slave_nodes = GetGeometry( ).PointsNumber( );

    // Calculate the normal at the integration point
    for ( unsigned int i = 0; i < num_slave_nodes; ++i )
    {
        const array_1d<double, 3> & Normal = GetGeometry()[i].FastGetSolutionStepValue(NORMAL);
        rVariables.IntegrationPointNormalVector( 0 ) += rVariables.N_Slave( i ) * Normal[0];
        rVariables.IntegrationPointNormalVector( 1 ) += rVariables.N_Slave( i ) * Normal[1];
    }

    // Calculate the vector between the gauss point and its projection on the master surface
    Vector x_gp = Vector( dimension );
    x_gp( 0 ) = GetGeometry( ).IntegrationPoints( mThisIntegrationMethod )[rPointNumber].Coordinate( 0 );
    x_gp( 1 ) = GetGeometry( ).IntegrationPoints( mThisIntegrationMethod )[rPointNumber].Coordinate( 1 );

    Vector x_gp_proj = ZeroVector( dimension ); // FIXME use the mapper that vicente wrote

    rVariables.IntegrationPointNormalGap = inner_prod( rVariables.IntegrationPointNormalVector, ( x_gp_proj - x_gp ) );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateDeltaNodalNormal( const IndexType& rPointNumber )
{
    // TODO: Add content
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateDeltaDiscreteGap( const IndexType& rPointNumber )
{
    // TODO: Add content
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateAndAppendMortarProjectionMatrixTranspose( GeneralVariables& rVariables )
{
    const unsigned int dimension = 2;
    const unsigned int num_active_nodes = mThisActiveSlaveNodes.size( );
    const unsigned int num_master_nodes = rVariables.GetMasterElement( ).size( );

    const MatrixType& D = mThisMortarConditionMatrices.D;
    const MatrixType& M = mThisMortarConditionMatrices.M;

    MatrixType P_transpose = ZeroMatrix( num_active_nodes * dimension, num_master_nodes * dimension );
    // P_transpose = inv( D_AA_tranpose ) * M_A_transpose
    // D_AA is always a diagonal matrix so it is inverted by simply dividing the diagonal elements
    unsigned int i = 0, j = 0;
    for ( unsigned int iNode = 0; iNode < num_master_nodes; ++iNode )
    {
        i = dimension * iNode;
        for ( unsigned int jActive = 0; jActive < num_active_nodes; ++jActive )
        {
            j = dimension * mThisActiveSlaveNodes[jActive];
            P_transpose( i    , j     ) = M( j, i ) / D( i, i );
            P_transpose( i + 1, j + 1 ) = P_transpose( i, j );
        }
    }

    mThisMortarConditionMatrices.AppendMortarProjectionMatrixTranspose( P_transpose );
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
    const unsigned int number_of_slave_nodes  =  slave_nodes.PointsNumber( );

    /* RESIZE MATRICES AND VECTORS */
    rVariables.Phi_LagrangeMultipliers.resize( number_of_slave_nodes );
    rVariables.N_Master.resize( number_of_master_nodes );
    rVariables.N_Slave.resize( number_of_slave_nodes );

    rVariables.DPhi_De_LagrangeMultipliers.resize( number_of_slave_nodes, slave_nodes.LocalSpaceDimension( ) );
    rVariables.DN_De_Master.resize( number_of_master_nodes, master_nodes.LocalSpaceDimension( ) );
    rVariables.DN_De_Slave.resize( number_of_slave_nodes, slave_nodes.LocalSpaceDimension( ) );

    /* POPULATE MATRICES AND VECTORS */
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

    KRATOS_CATCH( "" );
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

            if ( rLeftHandSideVariables[i] == CONTACT_STIFFNESS_MATRIX )
            {
                // Calculate
                MatrixType LHS_contact_pair = ZeroMatrix( dimension * num_total_nodes, dimension * num_total_nodes );
                this->CalculateAndAddKco( LHS_contact_pair, rVariables, rIntegrationWeight );

                // Assemble
                MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrices( )[i];
                this->AssembleContactPairLHSToConditionSystem( rVariables.GetMasterElementIndex( ),
                                                               LHS_contact_pair,
                                                               rLeftHandSideMatrix );
                calculated = true;
            }

            if ( rLeftHandSideVariables[i] == GAP_DERIVATIVES_MATRIX )
            {
                // Calculate
                MatrixType LHS_contact_pair = ZeroMatrix( dimension * num_total_nodes, dimension * num_total_nodes );
                this->CalculateAndAddGapDerivatives( LHS_contact_pair, rVariables, rIntegrationWeight );

                // Assemble
                MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrices( )[i];
                this->AssembleContactPairLHSToConditionSystem( rVariables.GetMasterElementIndex( ),
                                                               LHS_contact_pair,
                                                               rLeftHandSideMatrix );
                calculated = true;
            }

            if ( calculated == false )
            {
                KRATOS_THROW_ERROR( std::logic_error, " CONDITION can not supply the required local system variable: ", rLeftHandSideVariables[i] );
            }
        }
    }
    else
    {
        /* SINGLE LHS MATRIX */
        MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();
        MatrixType LHS_contact_pair = ZeroMatrix( dimension * num_total_nodes, dimension * num_total_nodes );

        // Calculate
        this->CalculateAndAddKco( LHS_contact_pair, rVariables, rIntegrationWeight );
        this->CalculateAndAddGapDerivatives( LHS_contact_pair, rVariables, rIntegrationWeight );

        // Assemble
        this->AssembleContactPairLHSToConditionSystem( rVariables.GetMasterElementIndex( ),
                                                       LHS_contact_pair,
                                                       rLeftHandSideMatrix );

//         KRATOS_WATCH( rLeftHandSideMatrix );
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::AssembleContactPairLHSToConditionSystem(
    const unsigned int rMasterElementIndex,
    MatrixType& rPairLHS,
    MatrixType& rConditionLHS 
    )
{
    const unsigned int dimension = 2;
    const unsigned int num_slave_nodes = GetGeometry( ).PointsNumber( );
    const unsigned int num_current_master_nodes = mThisMasterElements[rMasterElementIndex]->PointsNumber( );

    // Find location of the piar's master DOFs in ConditionLHS
    unsigned int index_begin = 0;
    unsigned int index_end = 0;

    for (unsigned int i_master_elem = 0; i_master_elem < rMasterElementIndex - 1; ++i_master_elem)
    {
        index_begin += mThisMasterElements[i_master_elem]->PointsNumber( );
    }

    index_end = index_begin + num_current_master_nodes;

    index_begin *= dimension;
    index_end   *= dimension;

    Matrix pair_LHS_master = subrange( rPairLHS, 0, num_current_master_nodes * dimension, 0, rPairLHS.size2( ) );
    Matrix pair_LHS_slave  = subrange( rPairLHS, num_current_master_nodes * dimension , rPairLHS.size1( ), 0, rPairLHS.size2( )  );

    Matrix LHS_master = subrange( rConditionLHS, index_begin, index_end, 0, rPairLHS.size2( )  );
    Matrix LHS_slave  = subrange( rConditionLHS, rConditionLHS.size1( ) - num_slave_nodes, rConditionLHS.size1( ), 0, rPairLHS.size2( )  );

    LHS_master = pair_LHS_master;
    LHS_slave += pair_LHS_slave;
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateAndAddRHS( 
    LocalSystemComponents& rLocalSystem,
    GeneralVariables& rVariables,
    double& rIntegrationWeight 
    )
{
    // TODO Check RHS Variables and discuss with Vicente whether those are the final ones or not !!

    // Contact pair variables
    const unsigned int dimension = 2;
    const unsigned int num_master_nodes = rVariables.GetMasterElement( ).PointsNumber( );
    const unsigned int num_slave_nodes = GetGeometry( ).PointsNumber( );
    const unsigned int num_total_nodes = num_slave_nodes + num_master_nodes;

    if ( rLocalSystem.CalculationFlags.Is( MortarContact2DCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
    {
        /* COMPONENT-WISE LHS MATRIX */
        const std::vector<Variable<VectorType>>& rRightHandSideVariables = rLocalSystem.GetRightHandSideVariables( );
        bool calculated; // TODO: Change by a local flag
//         std::vector<VectorType>& rRightHandSideVectors = rLocalSystem.GetRightHandSideVectors( );

        for ( unsigned int i = 0; i < rRightHandSideVariables.size( ); i++ )
        {
            calculated = false;

            if ( rRightHandSideVariables[i] == CONTACT_FORCES_VECTOR )
            {
                // Calculate
                Vector RHS_contact_pair = ZeroVector( dimension * num_total_nodes );
                this->CalculateAndAddContactForces( RHS_contact_pair, rVariables, rIntegrationWeight );

                // Assemble
                VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVectors( )[i];
                this->AssembleContactPairRHSToConditionSystem( rVariables.GetMasterElementIndex( ),
                                                                RHS_contact_pair,
                                                                rRightHandSideVector );
                calculated = true;
            }
            
            if ( rRightHandSideVariables[i] == NORMAL_GAPS_VECTOR )
            {
                // Calculate
                Vector RHS_contact_pair = ZeroVector( dimension * num_total_nodes );
                this->CalculateAndAddGapsVector( RHS_contact_pair, rVariables, rIntegrationWeight );

                // Assemble
                VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVectors( )[i];
                this->AssembleContactPairRHSToConditionSystem( rVariables.GetMasterElementIndex( ),
                                                                RHS_contact_pair,
                                                                rRightHandSideVector );
                calculated = true;
            }

            if ( calculated == false )
            {
                KRATOS_THROW_ERROR( std::logic_error, " CONDITION can not supply the required local system variable: ", rRightHandSideVariables[i] );
            }
        }
    }
    else
    {
        /* SINGLE LHS MATRIX */
        VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector( );
        Vector RHS_contact_pair = ZeroVector( dimension * num_total_nodes );

        this->CalculateAndAddContactForces( RHS_contact_pair, rVariables, rIntegrationWeight );
        this->CalculateAndAddGapsVector(    RHS_contact_pair, rVariables, rIntegrationWeight );

        this->AssembleContactPairRHSToConditionSystem( rVariables.GetMasterElementIndex( ),
                                                       RHS_contact_pair,
                                                       rRightHandSideVector );
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::AssembleContactPairRHSToConditionSystem( 
    const unsigned int rMasterElementIndex,
    VectorType& rPairRHS,
    VectorType& rConditionRHS 
    )
{
    const unsigned int dimension = 2;
    const unsigned int num_slave_nodes = GetGeometry( ).PointsNumber( );
    const unsigned int num_current_master_nodes = mThisMasterElements[rMasterElementIndex]->PointsNumber( );

    // Find location of the piar's master DOFs in ConditionRHS
    unsigned int index_begin = 0;
    unsigned int index_end = 0;

    for (unsigned int i_master_elem = 0; i_master_elem < rMasterElementIndex - 1; ++i_master_elem)
    {
        index_begin += mThisMasterElements[i_master_elem]->PointsNumber( );
    }

    index_end = index_begin + num_current_master_nodes;

    index_begin *= dimension;
    index_end *= dimension;

    Vector pair_RHS_master = subrange( rPairRHS, 0, num_current_master_nodes * dimension );
    Vector pair_RHS_slave  = subrange( rPairRHS, num_current_master_nodes * dimension , rPairRHS.size( ) );

    Vector RHS_master = subrange( rConditionRHS, index_begin, index_end );
    Vector RHS_slave  = subrange( rConditionRHS, rConditionRHS.size( ) - num_slave_nodes, rConditionRHS.size( ) );

    RHS_master = pair_RHS_master;
    RHS_slave += pair_RHS_slave;
}

/***********************************************************************************/
/**************** AUXILLIARY METHODS FOR CONDITION LHS CONTRIBUTION ****************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateAndAddKco( 
    MatrixType& rLeftHandSideMatrix,
    GeneralVariables& rVariables,
    double& rIntegrationWeight )
{
    KRATOS_TRY;

    // Contact pair variables
    const unsigned int dimension = 2;
    const unsigned int num_master_nodes = rVariables.GetMasterElement( ).PointsNumber( );
    const unsigned int num_active_nodes = mThisActiveSlaveNodes.size( );

    // Contact stiffness matrix blocks
    MatrixType K_co_master = ZeroMatrix( dimension * num_master_nodes, dimension * num_active_nodes );
    MatrixType K_co_active = ZeroMatrix( dimension * num_active_nodes, dimension * num_active_nodes );

    // Calculate first the projection operator - for the sake of performance it is already calculated tranposed
    const Matrix& P_transpose = mThisMortarConditionMatrices.GetMortarProjectionMatrixTranspose( rVariables.GetMasterElementIndex( ) );
    const Matrix& DeltaD      = mThisMortarConditionMatrices.DeltaD;
    const Matrix& DeltaM      = mThisMortarConditionMatrices.DeltaM;

    // Calculate the blocks K_co_master and K_co_active
    unsigned int i = 0, j = 0, jNode = 0;
    for ( unsigned int iMaster = 0; iMaster < num_master_nodes; ++iMaster )
    {
        i = iMaster * dimension;
        for ( unsigned int jActive = 0; jActive < num_active_nodes; ++jActive )
        {
            jNode = mThisActiveSlaveNodes[jActive];
            j = jNode * dimension;

            const array_1d<double, 3> & lagrange_mult = GetGeometry()[jNode].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER);

            K_co_master( i,     j     ) -= DeltaM( j    , i     ) * lagrange_mult[0];
            K_co_master( i + 1, j + 1 ) -= DeltaM( j + 1, i + 1 ) * lagrange_mult[1];

            K_co_active( j    , j     ) += DeltaD( j    , j     ) * lagrange_mult[0];
            K_co_active( j + 1, j + 1 ) += DeltaD( j + 1, j + 1 ) * lagrange_mult[1];
        }
    }

    // Fill up the rows of the LHS associated with the master nodes' residuals
    Matrix rLHS_master = subrange( rLeftHandSideMatrix, 0, num_master_nodes * dimension, 0, rLeftHandSideMatrix.size2( )  );
    rLHS_master += K_co_master + prod( P_transpose, K_co_active );

    KRATOS_CATCH( "" );
}


/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateAndAddGapDerivatives( 
    MatrixType& rLeftHandSideMatrix,
    GeneralVariables& rVariables,
    double& rIntegrationWeight 
    )
{
    KRATOS_TRY;

    // Contact pair variables
    const unsigned int dimension = 2;
    const unsigned int num_master_nodes = rVariables.GetMasterElement( ).PointsNumber( );
    const unsigned int num_slave_nodes = GetGeometry( ).PointsNumber( );
    const unsigned int num_active_nodes = mThisActiveSlaveNodes.size( );
    const unsigned int num_total_nodes = num_slave_nodes + num_master_nodes;

    // Condition matrices
    const Matrix& DeltaJ_s = mThisMortarConditionMatrices.DeltaJSlave;
    const Matrix& DeltaPhi = mThisMortarConditionMatrices.DeltaPhiLagrangeMultipliers;

    const Vector& Phi = rVariables.Phi_LagrangeMultipliers;
    const double& J_s = rVariables.DetJSlave;
    const double& g_n = rVariables.IntegrationPointNormalGap;

    // Gaps lie in the last block of the system of equations of the condition pair
    unsigned int iLHS = num_total_nodes - num_active_nodes ;
    unsigned int iNode = 0;
    for ( unsigned int iActive = 0; iActive < num_active_nodes; ++iActive, ++iLHS )
    {
        iNode = mThisActiveSlaveNodes[iActive];
        Matrix N_im = subrange( rLeftHandSideMatrix, iLHS, iLHS+1, 0, num_master_nodes * dimension );
        Matrix N_is = subrange( rLeftHandSideMatrix, iLHS, iLHS+1, num_master_nodes * dimension, num_total_nodes * dimension );

        Matrix DeltaDiscreteGap_im = subrange( rVariables.DeltaDiscreteGap, iNode, iNode+1, 0, num_master_nodes * dimension );
        Matrix DeltaDiscreteGap_is = subrange( rVariables.DeltaDiscreteGap, iNode, iNode+1, num_master_nodes * dimension, num_total_nodes * dimension );
        Matrix DeltaPhi_i = subrange( DeltaPhi, iNode, iNode+1, 0, num_slave_nodes );
        Matrix DeltaJ_i   = subrange( DeltaJ_s, iNode, iNode+1, 0, num_slave_nodes );

        N_im += rIntegrationWeight * ( Phi( iNode ) * J_s * DeltaDiscreteGap_im );

        N_is += rIntegrationWeight * ( g_n * J_s * DeltaPhi_i +
                                       Phi( iNode ) * g_n * DeltaJ_i +
                                       Phi( iNode ) * J_s * DeltaDiscreteGap_is );
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/**************** AUXILLIARY METHODS FOR CONDITION RHS CONTRIBUTION ****************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateAndAddGapsVector( 
    VectorType& rRightHandSideVector,
    GeneralVariables& rVariables,
    double& rIntegrationWeight 
    )
{
    const unsigned int num_active_nodes = mThisActiveSlaveNodes.size( );
    const unsigned int num_master_nodes = rVariables.GetMasterElement( ).PointsNumber( );
    const unsigned int num_slave_nodes  = GetGeometry( ).PointsNumber( );
    const unsigned int num_total_nodes  = num_master_nodes + num_slave_nodes;

    const Vector& Phi = rVariables.Phi_LagrangeMultipliers;
    const double& J_s = rVariables.DetJSlave;
    const double& g_n = rVariables.IntegrationPointNormalGap;

    // Gaps lie in the last block of the system of equations of the condition
    unsigned int iRHS = num_total_nodes - num_active_nodes ;
    unsigned int iNode = 0;
    for ( unsigned int iActive = 0; iActive < num_active_nodes; ++iActive, ++iRHS )
    {
        iNode = mThisActiveSlaveNodes[iActive];
        rRightHandSideVector( iRHS ) += rIntegrationWeight * Phi(iNode) * g_n * J_s;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::CalculateAndAddContactForces( 
    VectorType& rRightHandSideVector,
    GeneralVariables& rVariables,
    double& rIntegrationWeight 
    )
{
    KRATOS_TRY;

    // Contact pair variables
    const unsigned int dimension = 2;
    const unsigned int num_active_nodes = mThisActiveSlaveNodes.size( );
    const unsigned int num_master_nodes = rVariables.GetMasterElement( ).PointsNumber( );
    const unsigned int num_slave_nodes  = GetGeometry( ).PointsNumber( );
    const unsigned int num_total_nodes  = num_master_nodes + num_slave_nodes;

    const Matrix& D = mThisMortarConditionMatrices.D;
    const Matrix& M = mThisMortarConditionMatrices.M;

    unsigned int i = 0, j = 0, jNode = 0;
    for ( unsigned int iMaster = 0; iMaster < num_master_nodes; ++iMaster )
    {
        i = iMaster * dimension;
        for ( unsigned int jActive = 0; jActive < num_active_nodes; ++jActive )
        {
            jNode = mThisActiveSlaveNodes[jActive];
            j = jNode * dimension;
            const array_1d<double, 3> & lagrange_mult = GetGeometry()[jNode].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER);

            rRightHandSideVector( i )     -= M( j    , i     ) * lagrange_mult[0];
            rRightHandSideVector( i + 1 ) -= M( j + 1, i + 1 ) * lagrange_mult[1];

            rRightHandSideVector( num_total_nodes - num_active_nodes + j )   += D( j    , j     ) * lagrange_mult[0];
            rRightHandSideVector( num_total_nodes - num_active_nodes + j+1 ) += D( j + 1, j + 1 ) * lagrange_mult[1];
        }
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::EquationIdVector( 
    EquationIdVectorType& rResult,
    ProcessInfo& CurrentProcessInfo 
    )
{
    KRATOS_TRY;

    rResult.resize( 0, false );

    /* ORDER - [ MASTER, SLAVE_I, SLAVE_A ] */

    // Master Nodes Equation IDs
    for (unsigned int i_master_elem = 0; i_master_elem < mThisMasterElements.size( ); ++i_master_elem)
    {
        const unsigned num_master_nodes = mThisMasterElements[i_master_elem]->PointsNumber( );
        for ( unsigned int iNode = 0; iNode < num_master_nodes; iNode++ )
        {
            NodeType& master_node = ( *( mThisMasterElements[i_master_elem] ) )[iNode];
            rResult.push_back( master_node.GetDof( DISPLACEMENT_X ).EquationId( ) );
            rResult.push_back( master_node.GetDof( DISPLACEMENT_Y ).EquationId( ) );
        }
    }

    // Inactive Nodes Equation IDs
    for (unsigned int i_inactive = 0; i_inactive < mThisInactiveSlaveNodes.size( ); ++i_inactive)
    {
        const unsigned int iNode = mThisInactiveSlaveNodes[i_inactive];
        NodeType& slave_node = GetGeometry( )[iNode];
        if( slave_node.IsNot( ACTIVE ) )
        {
            rResult.push_back( slave_node.GetDof( DISPLACEMENT_X ).EquationId( ) );
            rResult.push_back( slave_node.GetDof( DISPLACEMENT_Y ).EquationId( ) );
        }
    }

    // Active Nodes Equation IDs
    for (unsigned int i_active = 0; i_active < mThisInactiveSlaveNodes.size( ); ++i_active)
    {
        const unsigned int iNode = mThisInactiveSlaveNodes[i_active];
        NodeType& slave_node = GetGeometry( )[iNode];
        if( slave_node.Is( ACTIVE ) )
        {
            rResult.push_back( slave_node.GetDof( DISPLACEMENT_X ).EquationId( ) );
            rResult.push_back( slave_node.GetDof( DISPLACEMENT_Y ).EquationId( ) );
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::GetDofList( 
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    rConditionalDofList.resize( 0 );

    /* ORDER - [ MASTER, SLAVE_I, SLAVE_A ] */

    // Master Nodes Equation IDs
    for ( unsigned int i_master_elem = 0; i_master_elem < mThisMasterElements.size( ); ++i_master_elem )
    {
        const unsigned num_master_nodes = mThisMasterElements[i_master_elem]->PointsNumber( );
        for ( unsigned int iNode = 0; iNode < num_master_nodes; iNode++ )
        {
            NodeType& master_node = ( *( mThisMasterElements[i_master_elem] ) )[iNode];
            rConditionalDofList.push_back( master_node.pGetDof( DISPLACEMENT_X ) );
            rConditionalDofList.push_back( master_node.pGetDof( DISPLACEMENT_Y ) );
        }
    }

    // Inactive Nodes Equation IDs
    for ( unsigned int i_inactive = 0; i_inactive < mThisInactiveSlaveNodes.size( ); ++i_inactive )
    {
        const unsigned int iNode = mThisInactiveSlaveNodes[i_inactive];
        NodeType& slave_node = GetGeometry( )[iNode];
        if( slave_node.IsNot( ACTIVE ) )
        {
            rConditionalDofList.push_back( slave_node.pGetDof( DISPLACEMENT_X ) );
            rConditionalDofList.push_back( slave_node.pGetDof( DISPLACEMENT_Y ) );
        }
    }

    // Active Nodes Equation IDs
    for (unsigned int i_active = 0; i_active < mThisInactiveSlaveNodes.size( ); ++i_active)
    {
        const unsigned int iNode = mThisInactiveSlaveNodes[i_active];
        NodeType& slave_node = GetGeometry( )[iNode];
        if( slave_node.Is( ACTIVE ) )
        {
            rConditionalDofList.push_back( slave_node.pGetDof( DISPLACEMENT_X ) );
            rConditionalDofList.push_back( slave_node.pGetDof( DISPLACEMENT_Y ) );
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void MortarContact2DCondition::GetMortarProjectionMatricesTranspose( std::vector<MatrixType>& rMortarProjectionMatricesTranspose )
{
    rMortarProjectionMatricesTranspose.resize( 0 );

    for (unsigned int i_master_elem = 0; i_master_elem < mThisMasterElements.size( ); ++i_master_elem)
    {
        rMortarProjectionMatricesTranspose.push_back(mThisMortarConditionMatrices.GetMortarProjectionMatrixTranspose( i_master_elem ) );
    }
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

    // TODO: ADD CONT ENT!!!!!
//     switch (rVariable)
//     {
//             case SOME_VARIABLE:
//                     for (int g = 0; g < number_of_integration_pts; ++g)
//                     {
//                     }
//                     break;
// 
//             default:
//                     break;
//     }

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
