// KRATOS  ___|  |       |       |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//           | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License: BSD License
//   license: structural_mechanics_application/license.txt
//
//  Main authors:  Vicente Mataix Ferrándiz
//

// System includes

// External includes

// Project includes
/* Mortar includes */
#include "custom_conditions/ALM_frictionless_mortar_contact_condition.h"

/* Additional includes */
#include <algorithm>

/* Utilities */
#include "custom_utilities/contact_utilities.h"
#include "utilities/math_utils.h"
#include "custom_utilities/structural_mechanics_math_utilities.hpp" // NOTE: Change for a more performant solver
// #include "../FSIapplication/custom_utilities/qr_utility.h" // QR decomposition utility used in matrix inversion.

namespace Kratos 
{
/**
 * Flags related to the condition computation 
 */
// Avoiding using the macro since this has a template parameter. If there was no template plase use the KRATOS_CREATE_LOCAL_FLAG macro
template< unsigned int TDim, unsigned int TNumNodes>
const Kratos::Flags AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR(Kratos::Flags::Create(0));
template< unsigned int TDim, unsigned int TNumNodes>
const Kratos::Flags AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX(Kratos::Flags::Create(1));
template< unsigned int TDim, unsigned int TNumNodes>
const Kratos::Flags AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS(Kratos::Flags::Create(2));
template< unsigned int TDim, unsigned int TNumNodes>
const Kratos::Flags AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS(Kratos::Flags::Create(3));

/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::Create( 
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties ) const
{
    return boost::make_shared< AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes> >( NewId, this->GetGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return boost::make_shared< AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes> >( NewId, pGeom, pProperties );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::~AugmentedLagrangianMethodFrictionlessMortarContactCondition( )
{
}

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::Initialize( ) 
{
    KRATOS_TRY;
    
    InitializeIntegrationMethod();
    
    // First populate of the vector of master elements
    const std::vector<contact_container> * all_containers = this->GetValue( CONTACT_CONTAINERS );
    mPairSize = all_containers->size();

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    // First populate of the vector of master elements
    const std::vector<contact_container> * all_containers = this->GetValue( CONTACT_CONTAINERS );
    mPairSize = all_containers->size();
    mThisMasterElements.resize( mPairSize );
    
    for ( unsigned int i_cond = 0; i_cond < mPairSize; ++i_cond )
    {
        mThisMasterElements[i_cond] = (*all_containers)[i_cond].condition;
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    // NOTE: Add things if necessary
        
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // NOTE: Add things if necessary
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // Create and initialize condition variables:#pragma omp critical
    GeneralVariables rVariables;
    
    // Create the current contact data
    DerivativeData rDerivativeData;
    
    // Create the mortar operators
    MortarConditionMatrices rThisMortarConditionMatrices;

//     // Reading integration points 
//     const GeometryType::IntegrationPointsArrayType& integration_points = mUseManualColocationIntegration ?
//                                                                          mColocationIntegration.IntegrationPoints( ) :
//                                                                          GetGeometry( ).IntegrationPoints( mThisIntegrationMethod );
                                                                                  
    this->InitializeDerivativeData(rDerivativeData, rCurrentProcessInfo);
    
    // Compute Ae and its derivative
    this->CalculateAeAndDeltaAe(rDerivativeData, rVariables, rCurrentProcessInfo); 
    
    // Iterate over the master segments
    for (unsigned int PairIndex = 0; PairIndex < mPairSize; ++PairIndex)
    {   
        // Reading integration points
        this->ComputeSelectiveIntegrationMethod(PairIndex);
        const GeometryType::IntegrationPointsArrayType& integration_points = mUseManualColocationIntegration ?
                                                                         mColocationIntegration.IntegrationPoints( ) :
                                                                         GetGeometry( ).IntegrationPoints( mThisIntegrationMethod );
                                                                         
        // Initialize general variables for the current master element
        this->InitializeGeneralVariables( rVariables, rCurrentProcessInfo, PairIndex );
        
        // Update the contact data
        this->UpdateDerivativeData(rDerivativeData, PairIndex);
        
        // Initialize the mortar operators
        rThisMortarConditionMatrices.Initialize();
         
        // Compute the normal derivatives of the master
        this->CalculateDeltaNormalMaster(rDerivativeData);
        
        // Initialize the integration weight
        double total_weight = 0.0;
        
        // Integrating the mortar operators
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            // Calculate the kinematic variables
            bool inside = this->CalculateKinematics( rVariables, rDerivativeData, PointNumber, integration_points );
            
            if (inside == true)
            {                    
                const double IntegrationWeight = integration_points[PointNumber].Weight();
                total_weight += IntegrationWeight;
                
                this->CalculateMortarOperators(rThisMortarConditionMatrices, rVariables, IntegrationWeight);
            }
        }
        
        // We can consider the pair if at least one of the collocation point is inside 
        if (total_weight > 0.0)
        {            
            // Setting the weighted gap
            // Mortar condition matrices - DOperator and MOperator
            const bounded_matrix<double, TNumNodes, TNumNodes>& DOperator = rThisMortarConditionMatrices.DOperator;
            const bounded_matrix<double, TNumNodes, TNumNodes>& MOperator = rThisMortarConditionMatrices.MOperator;
            
            // Current coordinates 
            const bounded_matrix<double, TNumNodes, TDim> x1 = GetCoordinates(this->GetGeometry());
            const bounded_matrix<double, TNumNodes, TDim> x2 = GetCoordinates(mThisMasterElements[PairIndex]->GetGeometry());
    
            const bounded_matrix<double, TNumNodes, TDim> Dx1Mx2 = prod(DOperator, x1) - prod(MOperator, x2); 
            
            for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
            {
                const array_1d<double, TDim> normal = subrange(GetGeometry()[iNode].GetValue(NORMAL), 0, TDim);
                const array_1d<double, TDim> aux_array = row(Dx1Mx2, iNode);
                
                #pragma omp atomic 
                GetGeometry()[iNode].GetValue(WEIGHTED_GAP) += inner_prod(aux_array, normal); 
            }
        }
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
IntegrationMethod AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::GetIntegrationMethod()
{   
    return mThisIntegrationMethod;
}

template< unsigned int TDim, unsigned int TNumNodes>
IntegrationMethod AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::GetIntegrationMethod(
    const unsigned int integration_order, 
    const bool collocation
    )
{
    if (collocation == false)
    {
        if (integration_order == 1)
        {
            return GeometryData::GI_GAUSS_1;
        }
        else if (integration_order == 2)
        {
            return GeometryData::GI_GAUSS_2;
        }
        else if (integration_order == 3)
        {
            return GeometryData::GI_GAUSS_3;
        }
        else if (integration_order == 4)
        {
            return GeometryData::GI_GAUSS_4;
        }
        else if (integration_order == 5)
        {
            return GeometryData::GI_GAUSS_5;
        }
        else
        {
            return GeometryData::GI_GAUSS_5; // NOTE: Maximium by default
//             return GetGeometry().GetDefaultIntegrationMethod();
        }
    }
    else
    {
        if (integration_order == 1)
        {
            return GeometryData::GI_EXTENDED_GAUSS_1;
        }
        else if (integration_order == 2)
        {
            return GeometryData::GI_EXTENDED_GAUSS_2;
        }
        else if (integration_order == 3)
        {
            return GeometryData::GI_EXTENDED_GAUSS_3;
        }
        else if (integration_order == 4)
        {
            return GeometryData::GI_EXTENDED_GAUSS_4;
        }
        else if (integration_order == 5)
        {
            return GeometryData::GI_EXTENDED_GAUSS_5;
        }
        else
        {
            return this->GetIntegrationMethod();
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateLocalSystem( 
    std::vector<MatrixType>& rLeftHandSideMatrices,
    const std::vector<Variable<MatrixType> >& rLHSVariables,
    std::vector<VectorType>& rRightHandSideVectors,
    const std::vector<Variable<VectorType> >& rRHSVariables,
    ProcessInfo& rCurrentProcessInfo 
    )
{    
    // Calculates the size of the system
    constexpr unsigned int TMatrixSize = TDim * (TNumNodes + TNumNodes) + TNumNodes;
    
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set(AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS, true);
    LocalSystem.CalculationFlags.Set(AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS, true);

    //Initialize sizes for the system components
    if ( rLHSVariables.size( ) != rLeftHandSideMatrices.size( ) )
    {
        rLeftHandSideMatrices.resize( rLHSVariables.size( ) );
    }

    if ( rRHSVariables.size( ) != rRightHandSideVectors.size( ) )
    {
        rRightHandSideVectors.resize( rRHSVariables.size( ) );
    }

    LocalSystem.CalculationFlags.Set(AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX, true);
    for ( unsigned int i = 0; i < rLeftHandSideMatrices.size( ); i++ )
    {
        // Note: rRightHandSideVectors.size() > 0
        this->InitializeSystemMatrices<TMatrixSize>( rLeftHandSideMatrices[i], rRightHandSideVectors[0],LocalSystem.CalculationFlags );
    }

    LocalSystem.CalculationFlags.Set( AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR, true );
    LocalSystem.CalculationFlags.Set( AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX, false ); // Temporarily only
    for ( unsigned int i = 0; i < rRightHandSideVectors.size( ); i++ )
    {
        // Note: rLeftHandSideMatrices.size() > 0
        this->InitializeSystemMatrices<TMatrixSize>( rLeftHandSideMatrices[0], rRightHandSideVectors[i], LocalSystem.CalculationFlags  );
    }
    LocalSystem.CalculationFlags.Set( AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX, true ); // Reactivated again

    // Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrices( rLeftHandSideMatrices );
    LocalSystem.SetRightHandSideVectors( rRightHandSideVectors );

    LocalSystem.SetLeftHandSideVariables( rLHSVariables );
    LocalSystem.SetRightHandSideVariables( rRHSVariables );

    // Calculate condition system
    this->CalculateConditionSystem<TMatrixSize>( LocalSystem, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;
    
    // Calculates the size of the system
    constexpr unsigned int TMatrixSize = TDim * (TNumNodes + TNumNodes) + TNumNodes;
    
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX, true );
    LocalSystem.CalculationFlags.Set( AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR, true );

    // Initialize sizes for the system components:
    this->InitializeSystemMatrices<TMatrixSize>( rLeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );
    
    // Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix( rLeftHandSideMatrix );
    LocalSystem.SetRightHandSideVector( rRightHandSideVector );

    // Calculate condition system
    this->CalculateConditionSystem<TMatrixSize>( LocalSystem, rCurrentProcessInfo );

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateLeftHandSide( 
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Calculates the size of the system
    constexpr unsigned int TMatrixSize = TDim * (TNumNodes + TNumNodes) + TNumNodes;
    
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX, true );

    VectorType RightHandSideVector = Vector( );

    // Initialize sizes for the system components:
    this->InitializeSystemMatrices<TMatrixSize>( rLeftHandSideMatrix, RightHandSideVector, LocalSystem.CalculationFlags );

    // Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix( rLeftHandSideMatrix );
    LocalSystem.SetRightHandSideVector( RightHandSideVector );

    // Calculate condition system
    this->CalculateConditionSystem<TMatrixSize>( LocalSystem, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateLeftHandSide( 
    std::vector< MatrixType >& rLeftHandSideMatrices,
    const std::vector< Variable< MatrixType > >& rLHSVariables,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Calculates the size of the system
    constexpr unsigned int TMatrixSize = TDim * (TNumNodes + TNumNodes) + TNumNodes;
    
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX, true );
    LocalSystem.CalculationFlags.Set( AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS, true );

    VectorType RightHandSideVector = Vector( );

    // Initialize size for the system components
    for( unsigned int i = 0; i < rLeftHandSideMatrices.size( ); i++ )
    {
        this->InitializeSystemMatrices<TMatrixSize>( rLeftHandSideMatrices[i], RightHandSideVector, LocalSystem.CalculationFlags );
    }

    LocalSystem.SetLeftHandSideMatrices( rLeftHandSideMatrices );
    LocalSystem.SetRightHandSideVector( RightHandSideVector );

    // Calculate condition system
    this->CalculateConditionSystem<TMatrixSize>( LocalSystem, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateRightHandSide( 
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Calculates the size of the system
    constexpr unsigned int TMatrixSize = TDim * (TNumNodes + TNumNodes) + TNumNodes;
    
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR, true);

    MatrixType LeftHandSideMatrix = Matrix( );

    // Initialize size for the system components
    this->InitializeSystemMatrices<TMatrixSize>( LeftHandSideMatrix, rRightHandSideVector,LocalSystem.CalculationFlags);

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix( LeftHandSideMatrix );
    LocalSystem.SetRightHandSideVector( rRightHandSideVector );

    // Calculate condition system
    this->CalculateConditionSystem<TMatrixSize>( LocalSystem, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateRightHandSide( 
    std::vector< VectorType >& rRightHandSideVectors,
    const std::vector< Variable< VectorType > >& rRHSVariables,
    ProcessInfo& rCurrentProcessInfo )
{
    // Calculates the size of the system
    constexpr unsigned int TMatrixSize = TDim * (TNumNodes + TNumNodes) + TNumNodes;
        
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR, true );
    LocalSystem.CalculationFlags.Set( AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS, true );

    MatrixType LeftHandSideMatrix = Matrix( );

    // Initialize size for the system components
    for( unsigned int i = 0; i < rRightHandSideVectors.size(); i++ )
    {
        this->InitializeSystemMatrices<TMatrixSize>( LeftHandSideMatrix, rRightHandSideVectors[i], LocalSystem.CalculationFlags );
    }

    // Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix( LeftHandSideMatrix );
    LocalSystem.SetRightHandSideVectors( rRightHandSideVectors );

    // Calculate condition system
    this->CalculateConditionSystem<TMatrixSize>( LocalSystem, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
template< unsigned int TMatrixSize >
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::InitializeSystemMatrices( 
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    Flags& rCalculationFlags 
    )
{
    const unsigned int condition_size = this->CalculateConditionSize<TMatrixSize>( );
    
    // Resizing as needed the LHS
    if ( rCalculationFlags.Is( AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX ) ) // Calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != condition_size )
        {
            rLeftHandSideMatrix.resize( condition_size, condition_size, false );
        }
        noalias( rLeftHandSideMatrix ) = ZeroMatrix( condition_size, condition_size ); // Resetting LHS
    }

    // Resizing as needed the RHS
    if ( rCalculationFlags.Is( AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR ) ) // Calculation of the matrix is required
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

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateMassMatrix( 
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

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateDampingMatrix( 
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

template< unsigned int TDim, unsigned int TNumNodes>
template< unsigned int TMatrixSize >
const unsigned int AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateConditionSize( )
{
    const unsigned int condition_size = mPairSize * TMatrixSize;
    
    return condition_size;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
template< unsigned int TMatrixSize>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim, TNumNodes>::CalculateConditionSystem( 
    LocalSystemComponents& rLocalSystem,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
    
    // Create and initialize condition variables:#pragma omp critical
    GeneralVariables rVariables;
    
    // Create the current contact data
    DerivativeData rDerivativeData;
    
    // Create the mortar operators
    MortarConditionMatrices rThisMortarConditionMatrices;

//     // Reading integration points 
//     const GeometryType::IntegrationPointsArrayType& integration_points = mUseManualColocationIntegration ?
//                                                                          mColocationIntegration.IntegrationPoints( ) :
//                                                                          GetGeometry( ).IntegrationPoints( mThisIntegrationMethod );
                                                                                  
    this->InitializeDerivativeData(rDerivativeData, rCurrentProcessInfo);
    
    // Compute Ae and its derivative
    this->CalculateAeAndDeltaAe(rDerivativeData, rVariables, rCurrentProcessInfo); 
    
    // Iterate over the master segments
    for (unsigned int PairIndex = 0; PairIndex < mPairSize; ++PairIndex)
    {   
        // Reading integration points
        this->ComputeSelectiveIntegrationMethod(PairIndex);
        const GeometryType::IntegrationPointsArrayType& integration_points = mUseManualColocationIntegration ?
                                                                         mColocationIntegration.IntegrationPoints( ) :
                                                                         GetGeometry( ).IntegrationPoints( mThisIntegrationMethod );
                                                                         
        // Initialize general variables for the current master element
        this->InitializeGeneralVariables( rVariables, rCurrentProcessInfo, PairIndex );
        
        // Update the contact data
        this->UpdateDerivativeData(rDerivativeData, PairIndex);
        
        // Initialize the mortar operators
        rThisMortarConditionMatrices.Initialize();
         
        // Compute the normal derivatives of the master
        this->CalculateDeltaNormalMaster(rDerivativeData);
        
        // Initialize the integration weight
        double total_weight = 0.0;
        
        // Integrating the mortar operators
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            // Calculate the kinematic variables
            bool inside = this->CalculateKinematics( rVariables, rDerivativeData, PointNumber, integration_points );
            
            if (inside == true)
            {   
                /* Update the derivatives */
                // Update the derivative of DetJ
                this->CalculateDeltaDetJSlave(rVariables, rDerivativeData);
                // Update the derivatives of the shape functions and the gap
                this->CalculateDeltaN(rVariables, rDerivativeData); // FIXME: This is the old version!!!!
                // The derivatives of the dual shape function
                this->CalculateDeltaPhi(rVariables, rDerivativeData);
                
                const double IntegrationWeight = integration_points[PointNumber].Weight();
                total_weight += IntegrationWeight;
                
                this->CalculateMortarOperators(rThisMortarConditionMatrices, rVariables, IntegrationWeight);
//                 this->CalculateMortarOperators(rThisMortarConditionMatrices, rVariables, rDerivativeData, IntegrationWeight);
            }
        }
        
        // We can consider the pair if at least one of the collocation point is inside 
        if (total_weight > 0.0)
        {
            std::cout << "--------------------------------------------------" << std::endl;
            rThisMortarConditionMatrices.print();
            KRATOS_WATCH(this->Id());
            KRATOS_WATCH(PairIndex);
            
            // Calculates the active/inactive combination pair
            const unsigned int ActiveInactive = GetActiveInactiveValue(this->GetGeometry());
            
            // Assemble of the matrix is required
            if ( rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX ) ||
                    rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
            {
                // Calculate the local contribution
                bounded_matrix<double, TMatrixSize, TMatrixSize> LHS_contact_pair = this->CalculateLocalLHS<TMatrixSize>( rThisMortarConditionMatrices, PairIndex, ActiveInactive);
                
                // Contributions to stiffness matrix calculated on the reference config
                this->CalculateAndAddLHS<TMatrixSize>( rLocalSystem, LHS_contact_pair, PairIndex );
                
//                 KRATOS_WATCH(LHS_contact_pair);
                LOG_MATRIX_PRETTY( LHS_contact_pair );
            }

            // Assemble of the vector is required
            if ( rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR ) ||
                    rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
            {
                // Calculate the local contribution
                const array_1d<double, TMatrixSize> RHS_contact_pair = this->CalculateLocalRHS<TMatrixSize>( rThisMortarConditionMatrices, PairIndex, ActiveInactive);
                
                // Contribution to previous step contact force and residuals vector
                this->CalculateAndAddRHS<TMatrixSize>( rLocalSystem, RHS_contact_pair, PairIndex );
                
//                 KRATOS_WATCH(RHS_contact_pair);
                LOG_VECTOR_PRETTY( RHS_contact_pair );
            }
            
//             // Setting the weighted gap
//             // Mortar condition matrices - DOperator and MOperator
//             const bounded_matrix<double, TNumNodes, TNumNodes>& DOperator = rThisMortarConditionMatrices.DOperator;
//             const bounded_matrix<double, TNumNodes, TNumNodes>& MOperator = rThisMortarConditionMatrices.MOperator;
//             
//             // Current coordinates 
//             const bounded_matrix<double, TNumNodes, TDim> x1 = GetCoordinates(this->GetGeometry());
//             const bounded_matrix<double, TNumNodes, TDim> x2 = GetCoordinates(mThisMasterElements[PairIndex]->GetGeometry());
//     
//             const bounded_matrix<double, TNumNodes, TDim> Dx1Mx2 = prod(DOperator, x1) - prod(MOperator, x2); 
//             
//             for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
//             {
//                 const array_1d<double, TDim> normal = subrange(GetGeometry()[iNode].GetValue(NORMAL), 0, TDim);
//                 const array_1d<double, TDim> aux_array = row(Dx1Mx2, iNode);
//                 
//                 #pragma omp atomic 
//                 GetGeometry()[iNode].GetValue(WEIGHTED_GAP) += inner_prod(aux_array, normal); 
//             }
        }
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::InitializeGeneralVariables(
    GeneralVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo,
    const unsigned int& rMasterElementIndex
    )
{
    // Master segment info
    GeometryType& CurrentMasterElement = mThisMasterElements[rMasterElementIndex]->GetGeometry();
    
    // Slave element info
    rVariables.Initialize();

    rVariables.SetMasterElement( CurrentMasterElement );
    rVariables.SetMasterElementIndex( rMasterElementIndex );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateAeAndDeltaAe(
    DerivativeData& rDerivativeData,
    GeneralVariables& rVariables,
//     const GeometryType::IntegrationPointsArrayType& integration_points,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    double total_weight = 0.0; // NOTE: The integral is supposed to be in the domain partially integrated, I don't know if consider any additional thing for the partial integration
    
    // We initilize the Ae components
    AeData rAeData;
    rAeData.Initialize();
    
    rDerivativeData.InitializeDeltaAeComponents();
    for (unsigned int PairIndex = 0; PairIndex < mPairSize; ++PairIndex)
    {   
        // Reading integration points
        this->ComputeSelectiveIntegrationMethod(PairIndex);
        const GeometryType::IntegrationPointsArrayType& integration_points = mUseManualColocationIntegration ?
                                                                         mColocationIntegration.IntegrationPoints( ) :
                                                                         GetGeometry( ).IntegrationPoints( mThisIntegrationMethod );

        // Initialize general variables for the current master element
        this->InitializeGeneralVariables( rVariables, rCurrentProcessInfo, PairIndex );
        
        // Update the contact data
        this->UpdateDerivativeData(rDerivativeData, PairIndex);
            
        // Calculating the proportion between the integrated area and segment area
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            // Calculate the kinematic variables
            bool inside = this->CalculateKinematics( rVariables, rDerivativeData, PointNumber, integration_points );
            
            if (inside == true)
            {   
                const double IntegrationWeight = integration_points[PointNumber].Weight();
                total_weight += IntegrationWeight;
                this->CalculateDeltaAeComponents(rVariables, rDerivativeData, rAeData, IntegrationWeight);
            }
        }
    }
    
    // We can consider the pair if at least one of the collocation point is inside (TODO: Change this if collocation is not used)
    if (total_weight > 0.0)
    {
        this->CalculateDeltaAe(rDerivativeData, rAeData);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::InitializeDerivativeData(
    DerivativeData& rDerivativeData,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // Slave element info
    rDerivativeData.Initialize(GetGeometry());
    
    /* LM */
    rDerivativeData.LagrangeMultipliers = GetVariableVector(GetGeometry(), NORMAL_CONTACT_STRESS, 0); 
    
    /* NORMALS */
    rDerivativeData.Normal_s = GetVariableMatrix(GetGeometry(),  NORMAL); 
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::UpdateDerivativeData(
    DerivativeData& rDerivativeData,
    const unsigned int& rMasterElementIndex
    )
{    
    // Slave element info
    rDerivativeData.UpdateMasterPair(mThisMasterElements[rMasterElementIndex] );
    
    /* NORMALS */
    for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
    {
//         array_1d<double,3> normal = this->GetValue(NORMAL);
        array_1d<double,3> normal = GetGeometry()[iNode].GetValue(NORMAL);
    }
}

/*********************************COMPUTE KINEMATICS*********************************/
/************************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
bool AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateKinematics( 
    GeneralVariables& rVariables,
    const DerivativeData rDerivativeData,
    const double& rPointNumber,
    const GeometryType::IntegrationPointsArrayType& integration_points
    )
{
    /* DEFINITIONS */
    GeometryType& slave_nodes  = GetGeometry( );

    /* LOCAL COORDINATES */
    const PointType& local_point = integration_points[rPointNumber].Coordinates();
    
    /*  POPULATE MATRICES AND VECTORS */
    
    /// SLAVE CONDITION ///
    
    // SHAPE FUNCTIONS 
    slave_nodes.ShapeFunctionsValues( rVariables.N_Slave, local_point.Coordinates() );
    rVariables.Phi_LagrangeMultipliers = prod(rDerivativeData.Ae, rVariables.N_Slave);
//     rVariables.Phi_LagrangeMultipliers = rVariables.N_Slave; // TODO: This could be needed in the future to be different than the standart shape functions 
    
    // SHAPE FUNCTION DERIVATIVES
    slave_nodes.ShapeFunctionsLocalGradients( rVariables.DN_De_Slave, local_point );
//     slave_nodes.ShapeFunctionsLocalGradients( rVariables.DN_De_Slave , local_point );// TODO: This could be needed in the future to be different than the standart shape functions
    
    // MASTER CONDITION
    const bool inside = this->MasterShapeFunctionValue( rVariables, local_point);
    
    /* CALCULATE JACOBIAN AND JACOBIAN DETERMINANT */
    slave_nodes.Jacobian( rVariables.j_Slave, local_point.Coordinates() );
    rVariables.DetJSlave = ContactUtilities::ContactElementDetJacobian( rVariables.j_Slave );
    
    return inside;
}
 
/***********************************************************************************/
/*************** METHODS TO CALCULATE THE CONTACT CONDITION MATRICES ***************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
bool AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::MasterShapeFunctionValue(
    GeneralVariables& rVariables,
    const PointType& local_point 
    )
{    
    GeometryType& master_seg = rVariables.GetMasterElement( );

    PointType projected_gp_global;
//     const array_1d<double,3> normal = this->GetValue(NORMAL);
    const array_1d<double,3> normal = ContactUtilities::GaussPointNormal(rVariables.N_Slave, GetGeometry());
    
    GeometryType::CoordinatesArrayType slave_gp_global;
    double aux_dist = 0.0;
    this->GetGeometry( ).GlobalCoordinates( slave_gp_global, local_point );
    ContactUtilities::ProjectDirection( master_seg, slave_gp_global, projected_gp_global, aux_dist, -normal ); // The opposite direction
    
    GeometryType::CoordinatesArrayType projected_gp_local;
    
    const bool inside = master_seg.IsInside( projected_gp_global.Coordinates( ), projected_gp_local ) ;
    
    if( inside == true )
    {
        // SHAPE FUNCTIONS 
        master_seg.ShapeFunctionsValues(         rVariables.N_Master,     projected_gp_local );         
        master_seg.ShapeFunctionsLocalGradients( rVariables.DN_De_Master, projected_gp_local );
    }
    
    return inside;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateMortarOperators(
    MortarConditionMatrices& rThisMortarConditionMatrices,
    GeneralVariables& rVariables,
    DerivativeData& rDerivativeData,
    const double& rIntegrationWeight
    )
{
    /* DEFINITIONS */
    const double J_s = rVariables.DetJSlave; 
    const VectorType Phi = rVariables.Phi_LagrangeMultipliers;
    const VectorType N1  = rVariables.N_Slave;
    const VectorType N2  = rVariables.N_Master;
    
    // Derivatives
    static const unsigned int Size1 = (TNumNodes * TDim);
    static const unsigned int Size2 = 2 * (TNumNodes * TDim);

    const array_1d<double, Size1> DeltaJ_s  = rDerivativeData.DeltaJ_s;
    const array_1d<array_1d<double, TNumNodes >, Size1> DeltaPhi = rDerivativeData.DeltaPhi;
    const array_1d<array_1d<double, TNumNodes >, Size2> DeltaN1  = rDerivativeData.DeltaN1;
    const array_1d<array_1d<double, TNumNodes >, Size2> DeltaN2  = rDerivativeData.DeltaN2;
    
    // Mortar condition matrices - DOperator and MOperator
    bounded_matrix<double, TNumNodes, TNumNodes>& DOperator = rThisMortarConditionMatrices.DOperator;
    bounded_matrix<double, TNumNodes, TNumNodes>& MOperator = rThisMortarConditionMatrices.MOperator;
    
    // D and M directional derivatives
    array_1d<bounded_matrix<double, TNumNodes, TNumNodes>, Size2>& DeltaDOperator = rThisMortarConditionMatrices.DeltaDOperator;
    array_1d<bounded_matrix<double, TNumNodes, TNumNodes>, Size2>& DeltaMOperator = rThisMortarConditionMatrices.DeltaMOperator;
    
    for (unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++)
    {
        for (unsigned int j_slave = 0; j_slave < TNumNodes; j_slave++)
        {
            const double phi = Phi[j_slave];
            const double n1  = N1[i_slave];
            const double n2  = N2[i_slave];
            
            DOperator(i_slave, j_slave) += J_s * rIntegrationWeight * phi * n1;
            MOperator(i_slave, j_slave) += J_s * rIntegrationWeight * phi * n2;
            
            for (unsigned int i = 0; i < TDim * TNumNodes; i++)
            {
                DeltaDOperator[i](i_slave, j_slave) += DeltaJ_s[i] * rIntegrationWeight * phi* n1           \
                                                     + J_s * rIntegrationWeight * DeltaPhi[i][j_slave] * n1 \
                                                     + J_s * rIntegrationWeight * phi* DeltaN1[i][i_slave];
                                                                            
                DeltaMOperator[i](i_slave, j_slave) += DeltaJ_s[i] * rIntegrationWeight * phi* n2           \
                                                     + J_s * rIntegrationWeight * DeltaPhi[i][j_slave] * n2 \
                                                     + J_s * rIntegrationWeight * phi* DeltaN2[i][i_slave];
            }
            for (unsigned int i = TDim * TNumNodes; i < 2 * TDim * TNumNodes; i++)
            {
                DeltaDOperator[i](i_slave, j_slave) += J_s * rIntegrationWeight * phi* DeltaN1[i][i_slave];
                                                                            
                DeltaMOperator[i](i_slave, j_slave) += J_s * rIntegrationWeight * phi* DeltaN2[i][i_slave];
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateMortarOperators(
    MortarConditionMatrices& rThisMortarConditionMatrices,
    GeneralVariables& rVariables,
    const double& rIntegrationWeight
    )
{
    /* DEFINITIONS */
    const double J_s = rVariables.DetJSlave; 
    const VectorType Phi = rVariables.Phi_LagrangeMultipliers;
    const VectorType N1  = rVariables.N_Slave;
    const VectorType N2  = rVariables.N_Master;

    // Mortar condition matrices - DOperator and MOperator
    bounded_matrix<double, TNumNodes, TNumNodes>& DOperator = rThisMortarConditionMatrices.DOperator;
    bounded_matrix<double, TNumNodes, TNumNodes>& MOperator = rThisMortarConditionMatrices.MOperator;
    
    for (unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++)
    {
        for (unsigned int j_slave = 0; j_slave < TNumNodes; j_slave++)
        {
            const double phi = Phi[j_slave];
            
            DOperator(i_slave, j_slave) += J_s * rIntegrationWeight * phi * N1[i_slave];
            MOperator(i_slave, j_slave) += J_s * rIntegrationWeight * phi * N2[i_slave];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateDeltaAeComponents(
    GeneralVariables& rVariables,
    DerivativeData& rDerivativeData,
    AeData& rAeData,
    const double& rIntegrationWeight
    )
{
    /* DEFINITIONS */
    const VectorType N1 = rVariables.N_Slave;
    const double detJ = rVariables.DetJSlave; 
     
    rAeData.De += rIntegrationWeight * this->ComputeDe( N1, detJ);
    rAeData.Me += rIntegrationWeight * this->ComputeMe( N1, detJ);
    
    for (unsigned int i = 0; i < TDim * TNumNodes; i++)
    {
        const double DeltaDetJ = rDerivativeData.DeltaJ_s[i];
        
        bounded_matrix<double, TNumNodes, TNumNodes> DeltaDe;
        const bounded_matrix<double, TNumNodes, TNumNodes> DeltaMe  = DeltaDetJ * outer_prod(N1, N1);
        
        for (unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++)
        {
            for (unsigned int j_slave = 0; j_slave < TNumNodes; j_slave++)
            {
                if (i_slave == j_slave)
                {
                    DeltaDe(i_slave, i_slave) = DeltaDetJ * N1[i_slave];
                }
                else
                {
                    DeltaDe(i_slave, j_slave) = 0.0;
                }
            }
        }
        
        rAeData.DeltaDe[i] += rIntegrationWeight * DeltaDe;
        rAeData.DeltaMe[i] += rIntegrationWeight * DeltaMe;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
bounded_matrix<double, TNumNodes, TNumNodes> AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::ComputeDe(
        const array_1d<double, TNumNodes> N1, 
        const double detJ 
        )
{
    bounded_matrix<double, TNumNodes, TNumNodes> De;
    
    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        for (unsigned int j = 0; j < TNumNodes; j++)
        {
            if (i == j)
            {
                De(i,i) = detJ * N1[i];
            }
        }
    }
    
    return De;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
bounded_matrix<double, TNumNodes, TNumNodes> AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::ComputeMe(
        const array_1d<double, TNumNodes> N1, 
        const double detJ 
        )
{
    bounded_matrix<double, TNumNodes, TNumNodes>  Me;
    
    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        for (unsigned int j = 0; j < TNumNodes; j++)
        {
            Me(i,j) = detJ * N1[i] * N1[j];
        }
    }
    
    return Me;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateDeltaAe(
    DerivativeData& rDerivativeData,
    AeData& rAeData
    )
{        
    Matrix InvMe = ZeroMatrix(TNumNodes, TNumNodes);
    // NOTE: Legacy inversion. In case Me is almost singular or singular (few GP integrated), will be considered as ZeroMatrix 
    if (TNumNodes == 2)
    {
        StructuralMechanicsMathUtilities::InvMat2x2(rAeData.Me, InvMe); // TODO: Change this for something more performant
    }
    else if (TNumNodes == 3)
    {
        StructuralMechanicsMathUtilities::InvMat3x3(rAeData.Me, InvMe);
    }   
    else
    {
        StructuralMechanicsMathUtilities::InvertMatrix(rAeData.Me, InvMe);
    }   
    
//     // Inversion using the QR decompisition // NOTE: Giving problems in the cases of almost singular matrix
//     QR<double, row_major> QR_decomposition;     
//     QR_decomposition.compute(TNumNodes, TNumNodes, &rDerivativeData.Me(0, 0));
//     
//     Matrix aux_I = identity_matrix<double>( TNumNodes ); // NOTE: Simplify this code¿?
//     for (unsigned int col = 0; col < TNumNodes; col++)
//     {   Vector aux_I_col = column(aux_I, col);
//         Vector aux_InvMe_col = ZeroVector(TNumNodes);
//         QR_decomposition.solve(&aux_I_col[0], &aux_InvMe_col[0]);
//         column(InvMe, col) = aux_InvMe_col;
//     }
    
    noalias(rDerivativeData.Ae) = prod(rAeData.De, InvMe);
    
    static const unsigned int Size1 = (TNumNodes * TDim);
    array_1d<bounded_matrix<double, TNumNodes, TNumNodes> , Size1>& DeltaAe = rDerivativeData.DeltaAe;
    
    for (unsigned int i = 0; i < TDim * TNumNodes; i++)
    {
        DeltaAe[i] = rAeData.DeltaDe[i] - prod(rDerivativeData.Ae, rAeData.DeltaMe[i]);
        DeltaAe[i] = prod(rDerivativeData.DeltaAe[i], InvMe);
//         DeltaAe[i] = ZeroMatrix(TNumNodes, TNumNodes); // NOTE: Test with zero derivative
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
template< unsigned int TMatrixSize >
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateAndAddLHS(
    LocalSystemComponents& rLocalSystem,
    const bounded_matrix<double, TMatrixSize, TMatrixSize>& LHS_contact_pair, 
    const unsigned int rPairIndex
    )
{
    /* DEFINITIONS */

    if ( rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
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
                this->AssembleContactPairLHSToConditionSystem<TMatrixSize>(LHS_contact_pair, rLeftHandSideMatrix, rPairIndex);
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
        this->AssembleContactPairLHSToConditionSystem<TMatrixSize>(LHS_contact_pair, rLeftHandSideMatrix, rPairIndex);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
template< unsigned int TMatrixSize>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::AssembleContactPairLHSToConditionSystem(
    const bounded_matrix<double, TMatrixSize, TMatrixSize>& rPairLHS,
    MatrixType& rConditionLHS,
    const unsigned int rPairIndex
    )
{
    // Find location of the pair's master DOFs in ConditionRHS
    const unsigned int index_begin = rPairIndex * TMatrixSize;
    const unsigned int index_end  = index_begin + TMatrixSize;
    
    subrange( rConditionLHS, index_begin, index_end, index_begin, index_end) += rPairLHS;
}

/***************************** BEGIN AD REPLACEMENT ********************************/
/***********************************************************************************/

template<>
template<>
bounded_matrix<double, 10, 10> AugmentedLagrangianMethodFrictionlessMortarContactCondition<2,2>::CalculateLocalLHS<10>(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const unsigned int& rMasterElementIndex,
        const unsigned int& rActiveInactive
        )
{
    bounded_matrix<double,10,10> lhs;
    
    // Master segment info
    GeometryType& CurrentMasterElement = mThisMasterElements[rMasterElementIndex]->GetGeometry();

    // Initialize values
    const bounded_matrix<double, 2, 2> u1 = GetVariableMatrix(this->GetGeometry(), DISPLACEMENT, 0);
    const bounded_matrix<double, 2, 2> u2 = GetVariableMatrix(CurrentMasterElement, DISPLACEMENT, 0);
    const bounded_matrix<double, 2, 2> X1 = GetCoordinates(this->GetGeometry(), false);
    const bounded_matrix<double, 2, 2> X2 = GetCoordinates(CurrentMasterElement, false);
    
    const array_1d<double, 2> lmnormal = GetVariableVector(this->GetGeometry(), NORMAL_CONTACT_STRESS, 0); 
    
    const bounded_matrix<double, 2, 2> normalslave = GetVariableMatrix(this->GetGeometry(),  NORMAL); 
    
    // Augmentation parameters
    double scale_factor = 1.0;
    double penalty_parameter = 0.0;
    if (GetProperties().Has(SCALE_FACTOR) == true)
    {
        scale_factor  = GetProperties().GetValue(SCALE_FACTOR);
    }
    if (GetProperties().Has(PENALTY_FACTOR) == true)
    {
        penalty_parameter = GetProperties().GetValue(PENALTY_FACTOR);
    }
    
    // Mortar operators
    const bounded_matrix<double, 2, 2> MOperator = rMortarConditionMatrices.MOperator;
    const bounded_matrix<double, 2, 2> DOperator = rMortarConditionMatrices.DOperator;
    // Mortar operators derivatives
    const array_1d<bounded_matrix<double, 2, 2>, 8> DeltaMOperator = rMortarConditionMatrices.DeltaMOperator;
    const array_1d<bounded_matrix<double, 2, 2>, 8> DeltaDOperator = rMortarConditionMatrices.DeltaDOperator;

    if (rActiveInactive == 0 )
    {
        const double clhs0 =     0.5*std::pow(scale_factor, 2.0)/penalty_parameter;
    
        lhs(0,0)=0;
        lhs(0,1)=0;
        lhs(0,2)=0;
        lhs(0,3)=0;
        lhs(0,4)=0;
        lhs(0,5)=0;
        lhs(0,6)=0;
        lhs(0,7)=0;
        lhs(0,8)=0;
        lhs(0,9)=0;
        lhs(1,0)=0;
        lhs(1,1)=0;
        lhs(1,2)=0;
        lhs(1,3)=0;
        lhs(1,4)=0;
        lhs(1,5)=0;
        lhs(1,6)=0;
        lhs(1,7)=0;
        lhs(1,8)=0;
        lhs(1,9)=0;
        lhs(2,0)=0;
        lhs(2,1)=0;
        lhs(2,2)=0;
        lhs(2,3)=0;
        lhs(2,4)=0;
        lhs(2,5)=0;
        lhs(2,6)=0;
        lhs(2,7)=0;
        lhs(2,8)=0;
        lhs(2,9)=0;
        lhs(3,0)=0;
        lhs(3,1)=0;
        lhs(3,2)=0;
        lhs(3,3)=0;
        lhs(3,4)=0;
        lhs(3,5)=0;
        lhs(3,6)=0;
        lhs(3,7)=0;
        lhs(3,8)=0;
        lhs(3,9)=0;
        lhs(4,0)=0;
        lhs(4,1)=0;
        lhs(4,2)=0;
        lhs(4,3)=0;
        lhs(4,4)=0;
        lhs(4,5)=0;
        lhs(4,6)=0;
        lhs(4,7)=0;
        lhs(4,8)=0;
        lhs(4,9)=0;
        lhs(5,0)=0;
        lhs(5,1)=0;
        lhs(5,2)=0;
        lhs(5,3)=0;
        lhs(5,4)=0;
        lhs(5,5)=0;
        lhs(5,6)=0;
        lhs(5,7)=0;
        lhs(5,8)=0;
        lhs(5,9)=0;
        lhs(6,0)=0;
        lhs(6,1)=0;
        lhs(6,2)=0;
        lhs(6,3)=0;
        lhs(6,4)=0;
        lhs(6,5)=0;
        lhs(6,6)=0;
        lhs(6,7)=0;
        lhs(6,8)=0;
        lhs(6,9)=0;
        lhs(7,0)=0;
        lhs(7,1)=0;
        lhs(7,2)=0;
        lhs(7,3)=0;
        lhs(7,4)=0;
        lhs(7,5)=0;
        lhs(7,6)=0;
        lhs(7,7)=0;
        lhs(7,8)=0;
        lhs(7,9)=0;
        lhs(8,0)=0;
        lhs(8,1)=0;
        lhs(8,2)=0;
        lhs(8,3)=0;
        lhs(8,4)=0;
        lhs(8,5)=0;
        lhs(8,6)=0;
        lhs(8,7)=0;
        lhs(8,8)=clhs0;
        lhs(8,9)=0;
        lhs(9,0)=0;
        lhs(9,1)=0;
        lhs(9,2)=0;
        lhs(9,3)=0;
        lhs(9,4)=0;
        lhs(9,5)=0;
        lhs(9,6)=0;
        lhs(9,7)=0;
        lhs(9,8)=0;
        lhs(9,9)=clhs0;
    }
    else if (rActiveInactive == 2 )
    {
        const double clhs0 =     MOperator(1,0); // MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs1 =     DeltaMOperator[4](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs2 =     X1(0,0) + u1(0,0);
        const double clhs3 =     DOperator(1,0); // DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs4 =     X1(1,0) + u1(1,0);
        const double clhs5 =     DOperator(1,1); // DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs6 =     X2(0,0) + u2(0,0);
        const double clhs7 =     X2(1,0) + u2(1,0);
        const double clhs8 =     MOperator(1,1); // MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs9 =     X1(0,1) + u1(0,1);
        const double clhs10 =     X1(1,1) + u1(1,1);
        const double clhs11 =     X2(0,1) + u2(0,1);
        const double clhs12 =     X2(1,1) + u2(1,1);
        const double clhs13 =     lmnormal[1]*scale_factor + penalty_parameter*(normalslave(1,0)*(-clhs0*clhs6 + clhs2*clhs3 + clhs4*clhs5 - clhs7*clhs8) + normalslave(1,1)*(-clhs0*clhs11 + clhs10*clhs5 - clhs12*clhs8 + clhs3*clhs9));
        const double clhs14 =     DeltaDOperator[4](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs15 =     DeltaDOperator[4](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs16 =     DeltaMOperator[4](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs17 =     normalslave(1,1)*(-clhs1*clhs11 + clhs10*clhs15 - clhs12*clhs16 + clhs14*clhs9);
        const double clhs18 =     normalslave(1,0)*(clhs0 + clhs1*clhs6 - clhs14*clhs2 - clhs15*clhs4 + clhs16*clhs7);
        const double clhs19 =     -clhs17 + clhs18;
        const double clhs20 =     clhs19*penalty_parameter;
        const double clhs21 =     -clhs0*clhs20 + clhs1*clhs13;
        const double clhs22 =     DeltaMOperator[5](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs23 =     DeltaDOperator[5](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs24 =     DeltaDOperator[5](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs25 =     DeltaMOperator[5](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs26 =     normalslave(1,0)*(clhs2*clhs23 - clhs22*clhs6 + clhs24*clhs4 - clhs25*clhs7) - normalslave(1,1)*(clhs0 - clhs10*clhs24 + clhs11*clhs22 + clhs12*clhs25 - clhs23*clhs9);
        const double clhs27 =     clhs26*penalty_parameter;
        const double clhs28 =     clhs0*clhs27 + clhs13*clhs22;
        const double clhs29 =     DeltaMOperator[6](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs30 =     DeltaDOperator[6](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs31 =     DeltaDOperator[6](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs32 =     DeltaMOperator[6](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs33 =     normalslave(1,1)*(clhs10*clhs31 - clhs11*clhs29 - clhs12*clhs32 + clhs30*clhs9);
        const double clhs34 =     normalslave(1,0)*(-clhs2*clhs30 + clhs29*clhs6 - clhs31*clhs4 + clhs32*clhs7 + clhs8);
        const double clhs35 =     -clhs33 + clhs34;
        const double clhs36 =     clhs35*penalty_parameter;
        const double clhs37 =     -clhs0*clhs36 + clhs13*clhs29;
        const double clhs38 =     DeltaMOperator[7](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs39 =     DeltaDOperator[7](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs40 =     DeltaDOperator[7](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs41 =     DeltaMOperator[7](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs42 =     normalslave(1,0)*(clhs2*clhs39 - clhs38*clhs6 + clhs4*clhs40 - clhs41*clhs7) - normalslave(1,1)*(-clhs10*clhs40 + clhs11*clhs38 + clhs12*clhs41 - clhs39*clhs9 + clhs8);
        const double clhs43 =     clhs42*penalty_parameter;
        const double clhs44 =     clhs0*clhs43 + clhs13*clhs38;
        const double clhs45 =     DeltaMOperator[0](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs46 =     DeltaDOperator[0](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs47 =     DeltaDOperator[0](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs48 =     DeltaMOperator[0](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs49 =     normalslave(1,0)*(clhs2*clhs46 + clhs3 + clhs4*clhs47 - clhs45*clhs6 - clhs48*clhs7) + normalslave(1,1)*(clhs10*clhs47 - clhs11*clhs45 - clhs12*clhs48 + clhs46*clhs9);
        const double clhs50 =     clhs49*penalty_parameter;
        const double clhs51 =     clhs0*clhs50 + clhs13*clhs45;
        const double clhs52 =     DeltaMOperator[1](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs53 =     DeltaDOperator[1](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs54 =     DeltaDOperator[1](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs55 =     DeltaMOperator[1](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs56 =     normalslave(1,0)*(clhs2*clhs53 + clhs4*clhs54 - clhs52*clhs6 - clhs55*clhs7) + normalslave(1,1)*(clhs10*clhs54 - clhs11*clhs52 - clhs12*clhs55 + clhs3 + clhs53*clhs9);
        const double clhs57 =     clhs56*penalty_parameter;
        const double clhs58 =     clhs0*clhs57 + clhs13*clhs52;
        const double clhs59 =     DeltaMOperator[2](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs60 =     DeltaDOperator[2](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs61 =     DeltaDOperator[2](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs62 =     DeltaMOperator[2](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs63 =     normalslave(1,0)*(clhs2*clhs60 + clhs4*clhs61 + clhs5 - clhs59*clhs6 - clhs62*clhs7) + normalslave(1,1)*(clhs10*clhs61 - clhs11*clhs59 - clhs12*clhs62 + clhs60*clhs9);
        const double clhs64 =     clhs63*penalty_parameter;
        const double clhs65 =     clhs0*clhs64 + clhs13*clhs59;
        const double clhs66 =     DeltaMOperator[3](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs67 =     DeltaDOperator[3](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs68 =     DeltaDOperator[3](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs69 =     DeltaMOperator[3](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs70 =     normalslave(1,0)*(clhs2*clhs67 + clhs4*clhs68 - clhs6*clhs66 - clhs69*clhs7) + normalslave(1,1)*(clhs10*clhs68 - clhs11*clhs66 - clhs12*clhs69 + clhs5 + clhs67*clhs9);
        const double clhs71 =     clhs70*penalty_parameter;
        const double clhs72 =     clhs0*clhs71 + clhs13*clhs66;
        const double clhs73 =     normalslave(1,0)*scale_factor;
        const double clhs74 =     normalslave(1,1)*scale_factor;
        const double clhs75 =     clhs13*clhs16 - clhs20*clhs8;
        const double clhs76 =     clhs13*clhs25 + clhs27*clhs8;
        const double clhs77 =     clhs13*clhs32 - clhs36*clhs8;
        const double clhs78 =     clhs13*clhs41 + clhs43*clhs8;
        const double clhs79 =     clhs13*clhs48 + clhs50*clhs8;
        const double clhs80 =     clhs13*clhs55 + clhs57*clhs8;
        const double clhs81 =     clhs13*clhs62 + clhs64*clhs8;
        const double clhs82 =     clhs13*clhs69 + clhs71*clhs8;
        const double clhs83 =     penalty_parameter*(clhs17 - clhs18);
        const double clhs84 =     clhs13*clhs14 + clhs3*clhs83;
        const double clhs85 =     clhs13*clhs23 + clhs27*clhs3;
        const double clhs86 =     penalty_parameter*(clhs33 - clhs34);
        const double clhs87 =     clhs13*clhs30 + clhs3*clhs86;
        const double clhs88 =     clhs13*clhs39 + clhs3*clhs43;
        const double clhs89 =     clhs13*clhs46 + clhs3*clhs50;
        const double clhs90 =     clhs13*clhs53 + clhs3*clhs57;
        const double clhs91 =     clhs13*clhs60 + clhs3*clhs64;
        const double clhs92 =     clhs13*clhs67 + clhs3*clhs71;
        const double clhs93 =     clhs13*clhs15 + clhs5*clhs83;
        const double clhs94 =     clhs13*clhs24 + clhs27*clhs5;
        const double clhs95 =     clhs13*clhs31 + clhs5*clhs86;
        const double clhs96 =     clhs13*clhs40 + clhs43*clhs5;
        const double clhs97 =     clhs13*clhs47 + clhs5*clhs50;
        const double clhs98 =     clhs13*clhs54 + clhs5*clhs57;
        const double clhs99 =     clhs13*clhs61 + clhs5*clhs64;
        const double clhs100 =     clhs13*clhs68 + clhs5*clhs71;
    
        lhs(0,0)=clhs21*normalslave(1,0);
        lhs(0,1)=clhs28*normalslave(1,0);
        lhs(0,2)=clhs37*normalslave(1,0);
        lhs(0,3)=clhs44*normalslave(1,0);
        lhs(0,4)=clhs51*normalslave(1,0);
        lhs(0,5)=clhs58*normalslave(1,0);
        lhs(0,6)=clhs65*normalslave(1,0);
        lhs(0,7)=clhs72*normalslave(1,0);
        lhs(0,8)=0;
        lhs(0,9)=clhs0*clhs73;
        lhs(1,0)=clhs21*normalslave(1,1);
        lhs(1,1)=clhs28*normalslave(1,1);
        lhs(1,2)=clhs37*normalslave(1,1);
        lhs(1,3)=clhs44*normalslave(1,1);
        lhs(1,4)=clhs51*normalslave(1,1);
        lhs(1,5)=clhs58*normalslave(1,1);
        lhs(1,6)=clhs65*normalslave(1,1);
        lhs(1,7)=clhs72*normalslave(1,1);
        lhs(1,8)=0;
        lhs(1,9)=clhs0*clhs74;
        lhs(2,0)=clhs75*normalslave(1,0);
        lhs(2,1)=clhs76*normalslave(1,0);
        lhs(2,2)=clhs77*normalslave(1,0);
        lhs(2,3)=clhs78*normalslave(1,0);
        lhs(2,4)=clhs79*normalslave(1,0);
        lhs(2,5)=clhs80*normalslave(1,0);
        lhs(2,6)=clhs81*normalslave(1,0);
        lhs(2,7)=clhs82*normalslave(1,0);
        lhs(2,8)=0;
        lhs(2,9)=clhs73*clhs8;
        lhs(3,0)=clhs75*normalslave(1,1);
        lhs(3,1)=clhs76*normalslave(1,1);
        lhs(3,2)=clhs77*normalslave(1,1);
        lhs(3,3)=clhs78*normalslave(1,1);
        lhs(3,4)=clhs79*normalslave(1,1);
        lhs(3,5)=clhs80*normalslave(1,1);
        lhs(3,6)=clhs81*normalslave(1,1);
        lhs(3,7)=clhs82*normalslave(1,1);
        lhs(3,8)=0;
        lhs(3,9)=clhs74*clhs8;
        lhs(4,0)=-clhs84*normalslave(1,0);
        lhs(4,1)=-clhs85*normalslave(1,0);
        lhs(4,2)=-clhs87*normalslave(1,0);
        lhs(4,3)=-clhs88*normalslave(1,0);
        lhs(4,4)=-clhs89*normalslave(1,0);
        lhs(4,5)=-clhs90*normalslave(1,0);
        lhs(4,6)=-clhs91*normalslave(1,0);
        lhs(4,7)=-clhs92*normalslave(1,0);
        lhs(4,8)=0;
        lhs(4,9)=-clhs3*clhs73;
        lhs(5,0)=-clhs84*normalslave(1,1);
        lhs(5,1)=-clhs85*normalslave(1,1);
        lhs(5,2)=-clhs87*normalslave(1,1);
        lhs(5,3)=-clhs88*normalslave(1,1);
        lhs(5,4)=-clhs89*normalslave(1,1);
        lhs(5,5)=-clhs90*normalslave(1,1);
        lhs(5,6)=-clhs91*normalslave(1,1);
        lhs(5,7)=-clhs92*normalslave(1,1);
        lhs(5,8)=0;
        lhs(5,9)=-clhs3*clhs74;
        lhs(6,0)=-clhs93*normalslave(1,0);
        lhs(6,1)=-clhs94*normalslave(1,0);
        lhs(6,2)=-clhs95*normalslave(1,0);
        lhs(6,3)=-clhs96*normalslave(1,0);
        lhs(6,4)=-clhs97*normalslave(1,0);
        lhs(6,5)=-clhs98*normalslave(1,0);
        lhs(6,6)=-clhs99*normalslave(1,0);
        lhs(6,7)=-clhs100*normalslave(1,0);
        lhs(6,8)=0;
        lhs(6,9)=-clhs5*clhs73;
        lhs(7,0)=-clhs93*normalslave(1,1);
        lhs(7,1)=-clhs94*normalslave(1,1);
        lhs(7,2)=-clhs95*normalslave(1,1);
        lhs(7,3)=-clhs96*normalslave(1,1);
        lhs(7,4)=-clhs97*normalslave(1,1);
        lhs(7,5)=-clhs98*normalslave(1,1);
        lhs(7,6)=-clhs99*normalslave(1,1);
        lhs(7,7)=-clhs100*normalslave(1,1);
        lhs(7,8)=0;
        lhs(7,9)=-clhs5*clhs74;
        lhs(8,0)=0;
        lhs(8,1)=0;
        lhs(8,2)=0;
        lhs(8,3)=0;
        lhs(8,4)=0;
        lhs(8,5)=0;
        lhs(8,6)=0;
        lhs(8,7)=0;
        lhs(8,8)=0.5*std::pow(scale_factor, 2.0)/penalty_parameter;
        lhs(8,9)=0;
        lhs(9,0)=clhs19*scale_factor;
        lhs(9,1)=-clhs26*scale_factor;
        lhs(9,2)=clhs35*scale_factor;
        lhs(9,3)=-clhs42*scale_factor;
        lhs(9,4)=-clhs49*scale_factor;
        lhs(9,5)=-clhs56*scale_factor;
        lhs(9,6)=-clhs63*scale_factor;
        lhs(9,7)=-clhs70*scale_factor;
        lhs(9,8)=0;
        lhs(9,9)=0;
    }
    else if (rActiveInactive == 1 )
    {
        const double clhs0 =     MOperator(0,0); // MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs1 =     DeltaMOperator[4](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs2 =     X1(0,0) + u1(0,0);
        const double clhs3 =     DOperator(0,0); // DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs4 =     X1(1,0) + u1(1,0);
        const double clhs5 =     DOperator(0,1); // DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs6 =     X2(0,0) + u2(0,0);
        const double clhs7 =     X2(1,0) + u2(1,0);
        const double clhs8 =     MOperator(0,1); // MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs9 =     X1(0,1) + u1(0,1);
        const double clhs10 =     X1(1,1) + u1(1,1);
        const double clhs11 =     X2(0,1) + u2(0,1);
        const double clhs12 =     X2(1,1) + u2(1,1);
        const double clhs13 =     lmnormal[0]*scale_factor + penalty_parameter*(normalslave(0,0)*(-clhs0*clhs6 + clhs2*clhs3 + clhs4*clhs5 - clhs7*clhs8) + normalslave(0,1)*(-clhs0*clhs11 + clhs10*clhs5 - clhs12*clhs8 + clhs3*clhs9));
        const double clhs14 =     DeltaDOperator[4](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs15 =     DeltaDOperator[4](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs16 =     DeltaMOperator[4](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs17 =     normalslave(0,1)*(-clhs1*clhs11 + clhs10*clhs15 - clhs12*clhs16 + clhs14*clhs9);
        const double clhs18 =     normalslave(0,0)*(clhs0 + clhs1*clhs6 - clhs14*clhs2 - clhs15*clhs4 + clhs16*clhs7);
        const double clhs19 =     -clhs17 + clhs18;
        const double clhs20 =     clhs19*penalty_parameter;
        const double clhs21 =     -clhs0*clhs20 + clhs1*clhs13;
        const double clhs22 =     DeltaMOperator[5](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs23 =     DeltaDOperator[5](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs24 =     DeltaDOperator[5](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs25 =     DeltaMOperator[5](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs26 =     normalslave(0,0)*(clhs2*clhs23 - clhs22*clhs6 + clhs24*clhs4 - clhs25*clhs7) - normalslave(0,1)*(clhs0 - clhs10*clhs24 + clhs11*clhs22 + clhs12*clhs25 - clhs23*clhs9);
        const double clhs27 =     clhs26*penalty_parameter;
        const double clhs28 =     clhs0*clhs27 + clhs13*clhs22;
        const double clhs29 =     DeltaMOperator[6](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs30 =     DeltaDOperator[6](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs31 =     DeltaDOperator[6](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs32 =     DeltaMOperator[6](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs33 =     normalslave(0,1)*(clhs10*clhs31 - clhs11*clhs29 - clhs12*clhs32 + clhs30*clhs9);
        const double clhs34 =     normalslave(0,0)*(-clhs2*clhs30 + clhs29*clhs6 - clhs31*clhs4 + clhs32*clhs7 + clhs8);
        const double clhs35 =     -clhs33 + clhs34;
        const double clhs36 =     clhs35*penalty_parameter;
        const double clhs37 =     -clhs0*clhs36 + clhs13*clhs29;
        const double clhs38 =     DeltaMOperator[7](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs39 =     DeltaDOperator[7](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs40 =     DeltaDOperator[7](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs41 =     DeltaMOperator[7](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs42 =     normalslave(0,0)*(clhs2*clhs39 - clhs38*clhs6 + clhs4*clhs40 - clhs41*clhs7) - normalslave(0,1)*(-clhs10*clhs40 + clhs11*clhs38 + clhs12*clhs41 - clhs39*clhs9 + clhs8);
        const double clhs43 =     clhs42*penalty_parameter;
        const double clhs44 =     clhs0*clhs43 + clhs13*clhs38;
        const double clhs45 =     DeltaMOperator[0](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs46 =     DeltaDOperator[0](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs47 =     DeltaDOperator[0](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs48 =     DeltaMOperator[0](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs49 =     normalslave(0,0)*(clhs2*clhs46 + clhs3 + clhs4*clhs47 - clhs45*clhs6 - clhs48*clhs7) + normalslave(0,1)*(clhs10*clhs47 - clhs11*clhs45 - clhs12*clhs48 + clhs46*clhs9);
        const double clhs50 =     clhs49*penalty_parameter;
        const double clhs51 =     clhs0*clhs50 + clhs13*clhs45;
        const double clhs52 =     DeltaMOperator[1](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs53 =     DeltaDOperator[1](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs54 =     DeltaDOperator[1](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs55 =     DeltaMOperator[1](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs56 =     normalslave(0,0)*(clhs2*clhs53 + clhs4*clhs54 - clhs52*clhs6 - clhs55*clhs7) + normalslave(0,1)*(clhs10*clhs54 - clhs11*clhs52 - clhs12*clhs55 + clhs3 + clhs53*clhs9);
        const double clhs57 =     clhs56*penalty_parameter;
        const double clhs58 =     clhs0*clhs57 + clhs13*clhs52;
        const double clhs59 =     DeltaMOperator[2](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs60 =     DeltaDOperator[2](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs61 =     DeltaDOperator[2](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs62 =     DeltaMOperator[2](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs63 =     normalslave(0,0)*(clhs2*clhs60 + clhs4*clhs61 + clhs5 - clhs59*clhs6 - clhs62*clhs7) + normalslave(0,1)*(clhs10*clhs61 - clhs11*clhs59 - clhs12*clhs62 + clhs60*clhs9);
        const double clhs64 =     clhs63*penalty_parameter;
        const double clhs65 =     clhs0*clhs64 + clhs13*clhs59;
        const double clhs66 =     DeltaMOperator[3](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs67 =     DeltaDOperator[3](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs68 =     DeltaDOperator[3](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs69 =     DeltaMOperator[3](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs70 =     normalslave(0,0)*(clhs2*clhs67 + clhs4*clhs68 - clhs6*clhs66 - clhs69*clhs7) + normalslave(0,1)*(clhs10*clhs68 - clhs11*clhs66 - clhs12*clhs69 + clhs5 + clhs67*clhs9);
        const double clhs71 =     clhs70*penalty_parameter;
        const double clhs72 =     clhs0*clhs71 + clhs13*clhs66;
        const double clhs73 =     normalslave(0,0)*scale_factor;
        const double clhs74 =     normalslave(0,1)*scale_factor;
        const double clhs75 =     clhs13*clhs16 - clhs20*clhs8;
        const double clhs76 =     clhs13*clhs25 + clhs27*clhs8;
        const double clhs77 =     clhs13*clhs32 - clhs36*clhs8;
        const double clhs78 =     clhs13*clhs41 + clhs43*clhs8;
        const double clhs79 =     clhs13*clhs48 + clhs50*clhs8;
        const double clhs80 =     clhs13*clhs55 + clhs57*clhs8;
        const double clhs81 =     clhs13*clhs62 + clhs64*clhs8;
        const double clhs82 =     clhs13*clhs69 + clhs71*clhs8;
        const double clhs83 =     penalty_parameter*(clhs17 - clhs18);
        const double clhs84 =     clhs13*clhs14 + clhs3*clhs83;
        const double clhs85 =     clhs13*clhs23 + clhs27*clhs3;
        const double clhs86 =     penalty_parameter*(clhs33 - clhs34);
        const double clhs87 =     clhs13*clhs30 + clhs3*clhs86;
        const double clhs88 =     clhs13*clhs39 + clhs3*clhs43;
        const double clhs89 =     clhs13*clhs46 + clhs3*clhs50;
        const double clhs90 =     clhs13*clhs53 + clhs3*clhs57;
        const double clhs91 =     clhs13*clhs60 + clhs3*clhs64;
        const double clhs92 =     clhs13*clhs67 + clhs3*clhs71;
        const double clhs93 =     clhs13*clhs15 + clhs5*clhs83;
        const double clhs94 =     clhs13*clhs24 + clhs27*clhs5;
        const double clhs95 =     clhs13*clhs31 + clhs5*clhs86;
        const double clhs96 =     clhs13*clhs40 + clhs43*clhs5;
        const double clhs97 =     clhs13*clhs47 + clhs5*clhs50;
        const double clhs98 =     clhs13*clhs54 + clhs5*clhs57;
        const double clhs99 =     clhs13*clhs61 + clhs5*clhs64;
        const double clhs100 =     clhs13*clhs68 + clhs5*clhs71;
    
        lhs(0,0)=clhs21*normalslave(0,0);
        lhs(0,1)=clhs28*normalslave(0,0);
        lhs(0,2)=clhs37*normalslave(0,0);
        lhs(0,3)=clhs44*normalslave(0,0);
        lhs(0,4)=clhs51*normalslave(0,0);
        lhs(0,5)=clhs58*normalslave(0,0);
        lhs(0,6)=clhs65*normalslave(0,0);
        lhs(0,7)=clhs72*normalslave(0,0);
        lhs(0,8)=clhs0*clhs73;
        lhs(0,9)=0;
        lhs(1,0)=clhs21*normalslave(0,1);
        lhs(1,1)=clhs28*normalslave(0,1);
        lhs(1,2)=clhs37*normalslave(0,1);
        lhs(1,3)=clhs44*normalslave(0,1);
        lhs(1,4)=clhs51*normalslave(0,1);
        lhs(1,5)=clhs58*normalslave(0,1);
        lhs(1,6)=clhs65*normalslave(0,1);
        lhs(1,7)=clhs72*normalslave(0,1);
        lhs(1,8)=clhs0*clhs74;
        lhs(1,9)=0;
        lhs(2,0)=clhs75*normalslave(0,0);
        lhs(2,1)=clhs76*normalslave(0,0);
        lhs(2,2)=clhs77*normalslave(0,0);
        lhs(2,3)=clhs78*normalslave(0,0);
        lhs(2,4)=clhs79*normalslave(0,0);
        lhs(2,5)=clhs80*normalslave(0,0);
        lhs(2,6)=clhs81*normalslave(0,0);
        lhs(2,7)=clhs82*normalslave(0,0);
        lhs(2,8)=clhs73*clhs8;
        lhs(2,9)=0;
        lhs(3,0)=clhs75*normalslave(0,1);
        lhs(3,1)=clhs76*normalslave(0,1);
        lhs(3,2)=clhs77*normalslave(0,1);
        lhs(3,3)=clhs78*normalslave(0,1);
        lhs(3,4)=clhs79*normalslave(0,1);
        lhs(3,5)=clhs80*normalslave(0,1);
        lhs(3,6)=clhs81*normalslave(0,1);
        lhs(3,7)=clhs82*normalslave(0,1);
        lhs(3,8)=clhs74*clhs8;
        lhs(3,9)=0;
        lhs(4,0)=-clhs84*normalslave(0,0);
        lhs(4,1)=-clhs85*normalslave(0,0);
        lhs(4,2)=-clhs87*normalslave(0,0);
        lhs(4,3)=-clhs88*normalslave(0,0);
        lhs(4,4)=-clhs89*normalslave(0,0);
        lhs(4,5)=-clhs90*normalslave(0,0);
        lhs(4,6)=-clhs91*normalslave(0,0);
        lhs(4,7)=-clhs92*normalslave(0,0);
        lhs(4,8)=-clhs3*clhs73;
        lhs(4,9)=0;
        lhs(5,0)=-clhs84*normalslave(0,1);
        lhs(5,1)=-clhs85*normalslave(0,1);
        lhs(5,2)=-clhs87*normalslave(0,1);
        lhs(5,3)=-clhs88*normalslave(0,1);
        lhs(5,4)=-clhs89*normalslave(0,1);
        lhs(5,5)=-clhs90*normalslave(0,1);
        lhs(5,6)=-clhs91*normalslave(0,1);
        lhs(5,7)=-clhs92*normalslave(0,1);
        lhs(5,8)=-clhs3*clhs74;
        lhs(5,9)=0;
        lhs(6,0)=-clhs93*normalslave(0,0);
        lhs(6,1)=-clhs94*normalslave(0,0);
        lhs(6,2)=-clhs95*normalslave(0,0);
        lhs(6,3)=-clhs96*normalslave(0,0);
        lhs(6,4)=-clhs97*normalslave(0,0);
        lhs(6,5)=-clhs98*normalslave(0,0);
        lhs(6,6)=-clhs99*normalslave(0,0);
        lhs(6,7)=-clhs100*normalslave(0,0);
        lhs(6,8)=-clhs5*clhs73;
        lhs(6,9)=0;
        lhs(7,0)=-clhs93*normalslave(0,1);
        lhs(7,1)=-clhs94*normalslave(0,1);
        lhs(7,2)=-clhs95*normalslave(0,1);
        lhs(7,3)=-clhs96*normalslave(0,1);
        lhs(7,4)=-clhs97*normalslave(0,1);
        lhs(7,5)=-clhs98*normalslave(0,1);
        lhs(7,6)=-clhs99*normalslave(0,1);
        lhs(7,7)=-clhs100*normalslave(0,1);
        lhs(7,8)=-clhs5*clhs74;
        lhs(7,9)=0;
        lhs(8,0)=clhs19*scale_factor;
        lhs(8,1)=-clhs26*scale_factor;
        lhs(8,2)=clhs35*scale_factor;
        lhs(8,3)=-clhs42*scale_factor;
        lhs(8,4)=-clhs49*scale_factor;
        lhs(8,5)=-clhs56*scale_factor;
        lhs(8,6)=-clhs63*scale_factor;
        lhs(8,7)=-clhs70*scale_factor;
        lhs(8,8)=0;
        lhs(8,9)=0;
        lhs(9,0)=0;
        lhs(9,1)=0;
        lhs(9,2)=0;
        lhs(9,3)=0;
        lhs(9,4)=0;
        lhs(9,5)=0;
        lhs(9,6)=0;
        lhs(9,7)=0;
        lhs(9,8)=0;
        lhs(9,9)=0.5*std::pow(scale_factor, 2.0)/penalty_parameter;
    }
    else if (rActiveInactive == 3 )
    {
        const double clhs0 =     MOperator(0,0); // MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs1 =     DeltaMOperator[4](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs2 =     X1(0,0) + u1(0,0);
        const double clhs3 =     DOperator(0,0); // DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs4 =     X1(1,0) + u1(1,0);
        const double clhs5 =     DOperator(0,1); // DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs6 =     X2(0,0) + u2(0,0);
        const double clhs7 =     X2(1,0) + u2(1,0);
        const double clhs8 =     MOperator(0,1); // MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs9 =     X1(0,1) + u1(0,1);
        const double clhs10 =     X1(1,1) + u1(1,1);
        const double clhs11 =     X2(0,1) + u2(0,1);
        const double clhs12 =     X2(1,1) + u2(1,1);
        const double clhs13 =     lmnormal[0]*scale_factor + penalty_parameter*(normalslave(0,0)*(-clhs0*clhs6 + clhs2*clhs3 + clhs4*clhs5 - clhs7*clhs8) + normalslave(0,1)*(-clhs0*clhs11 + clhs10*clhs5 - clhs12*clhs8 + clhs3*clhs9));
        const double clhs14 =     clhs13*normalslave(0,0);
        const double clhs15 =     MOperator(1,0); // MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs16 =     DeltaMOperator[4](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs17 =     DOperator(1,0); // DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs18 =     DOperator(1,1); // DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs19 =     MOperator(1,1); // MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs20 =     lmnormal[1]*scale_factor + penalty_parameter*(normalslave(1,0)*(-clhs15*clhs6 + clhs17*clhs2 + clhs18*clhs4 - clhs19*clhs7) + normalslave(1,1)*(clhs10*clhs18 - clhs11*clhs15 - clhs12*clhs19 + clhs17*clhs9));
        const double clhs21 =     clhs20*normalslave(1,0);
        const double clhs22 =     DeltaDOperator[4](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs23 =     DeltaDOperator[4](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs24 =     DeltaMOperator[4](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs25 =     normalslave(0,1)*(-clhs1*clhs11 + clhs10*clhs23 - clhs12*clhs24 + clhs22*clhs9);
        const double clhs26 =     normalslave(0,0)*(clhs0 + clhs1*clhs6 - clhs2*clhs22 - clhs23*clhs4 + clhs24*clhs7);
        const double clhs27 =     -clhs25 + clhs26;
        const double clhs28 =     clhs27*normalslave(0,0)*penalty_parameter;
        const double clhs29 =     DeltaDOperator[4](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs30 =     DeltaDOperator[4](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs31 =     DeltaMOperator[4](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs32 =     normalslave(1,1)*(clhs10*clhs30 - clhs11*clhs16 - clhs12*clhs31 + clhs29*clhs9);
        const double clhs33 =     normalslave(1,0)*(clhs15 + clhs16*clhs6 - clhs2*clhs29 - clhs30*clhs4 + clhs31*clhs7);
        const double clhs34 =     -clhs32 + clhs33;
        const double clhs35 =     clhs34*normalslave(1,0)*penalty_parameter;
        const double clhs36 =     DeltaMOperator[5](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs37 =     DeltaMOperator[5](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs38 =     DeltaDOperator[5](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs39 =     DeltaDOperator[5](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs40 =     DeltaMOperator[5](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs41 =     normalslave(0,0)*(clhs2*clhs38 - clhs36*clhs6 + clhs39*clhs4 - clhs40*clhs7) - normalslave(0,1)*(clhs0 - clhs10*clhs39 + clhs11*clhs36 + clhs12*clhs40 - clhs38*clhs9);
        const double clhs42 =     clhs41*normalslave(0,0)*penalty_parameter;
        const double clhs43 =     DeltaDOperator[5](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs44 =     DeltaDOperator[5](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs45 =     DeltaMOperator[5](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs46 =     normalslave(1,0)*(clhs2*clhs43 - clhs37*clhs6 + clhs4*clhs44 - clhs45*clhs7) - normalslave(1,1)*(-clhs10*clhs44 + clhs11*clhs37 + clhs12*clhs45 + clhs15 - clhs43*clhs9);
        const double clhs47 =     clhs46*normalslave(1,0)*penalty_parameter;
        const double clhs48 =     DeltaMOperator[6](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs49 =     DeltaMOperator[6](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs50 =     DeltaDOperator[6](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs51 =     DeltaDOperator[6](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs52 =     DeltaMOperator[6](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs53 =     normalslave(0,1)*(clhs10*clhs51 - clhs11*clhs48 - clhs12*clhs52 + clhs50*clhs9);
        const double clhs54 =     normalslave(0,0)*(-clhs2*clhs50 - clhs4*clhs51 + clhs48*clhs6 + clhs52*clhs7 + clhs8);
        const double clhs55 =     -clhs53 + clhs54;
        const double clhs56 =     clhs55*normalslave(0,0)*penalty_parameter;
        const double clhs57 =     DeltaDOperator[6](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs58 =     DeltaDOperator[6](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs59 =     DeltaMOperator[6](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs60 =     normalslave(1,1)*(clhs10*clhs58 - clhs11*clhs49 - clhs12*clhs59 + clhs57*clhs9);
        const double clhs61 =     normalslave(1,0)*(clhs19 - clhs2*clhs57 - clhs4*clhs58 + clhs49*clhs6 + clhs59*clhs7);
        const double clhs62 =     -clhs60 + clhs61;
        const double clhs63 =     clhs62*normalslave(1,0)*penalty_parameter;
        const double clhs64 =     DeltaMOperator[7](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs65 =     DeltaMOperator[7](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs66 =     DeltaDOperator[7](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs67 =     DeltaDOperator[7](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs68 =     DeltaMOperator[7](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs69 =     normalslave(0,0)*(clhs2*clhs66 + clhs4*clhs67 - clhs6*clhs64 - clhs68*clhs7) - normalslave(0,1)*(-clhs10*clhs67 + clhs11*clhs64 + clhs12*clhs68 - clhs66*clhs9 + clhs8);
        const double clhs70 =     clhs69*normalslave(0,0)*penalty_parameter;
        const double clhs71 =     DeltaDOperator[7](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs72 =     DeltaDOperator[7](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs73 =     DeltaMOperator[7](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs74 =     normalslave(1,0)*(clhs2*clhs71 + clhs4*clhs72 - clhs6*clhs65 - clhs7*clhs73) - normalslave(1,1)*(-clhs10*clhs72 + clhs11*clhs65 + clhs12*clhs73 + clhs19 - clhs71*clhs9);
        const double clhs75 =     clhs74*normalslave(1,0)*penalty_parameter;
        const double clhs76 =     DeltaMOperator[0](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs77 =     DeltaMOperator[0](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs78 =     DeltaDOperator[0](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs79 =     DeltaDOperator[0](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs80 =     DeltaMOperator[0](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs81 =     normalslave(0,0)*(clhs2*clhs78 + clhs3 + clhs4*clhs79 - clhs6*clhs76 - clhs7*clhs80) + normalslave(0,1)*(clhs10*clhs79 - clhs11*clhs76 - clhs12*clhs80 + clhs78*clhs9);
        const double clhs82 =     clhs81*normalslave(0,0)*penalty_parameter;
        const double clhs83 =     DeltaDOperator[0](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs84 =     DeltaDOperator[0](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs85 =     DeltaMOperator[0](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs86 =     normalslave(1,0)*(clhs17 + clhs2*clhs83 + clhs4*clhs84 - clhs6*clhs77 - clhs7*clhs85) + normalslave(1,1)*(clhs10*clhs84 - clhs11*clhs77 - clhs12*clhs85 + clhs83*clhs9);
        const double clhs87 =     clhs86*normalslave(1,0)*penalty_parameter;
        const double clhs88 =     DeltaMOperator[1](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs89 =     DeltaMOperator[1](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs90 =     DeltaDOperator[1](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs91 =     DeltaDOperator[1](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs92 =     DeltaMOperator[1](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs93 =     normalslave(0,0)*(clhs2*clhs90 + clhs4*clhs91 - clhs6*clhs88 - clhs7*clhs92) + normalslave(0,1)*(clhs10*clhs91 - clhs11*clhs88 - clhs12*clhs92 + clhs3 + clhs9*clhs90);
        const double clhs94 =     clhs93*normalslave(0,0)*penalty_parameter;
        const double clhs95 =     DeltaDOperator[1](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs96 =     DeltaDOperator[1](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs97 =     DeltaMOperator[1](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs98 =     normalslave(1,0)*(clhs2*clhs95 + clhs4*clhs96 - clhs6*clhs89 - clhs7*clhs97) + normalslave(1,1)*(clhs10*clhs96 - clhs11*clhs89 - clhs12*clhs97 + clhs17 + clhs9*clhs95);
        const double clhs99 =     clhs98*normalslave(1,0)*penalty_parameter;
        const double clhs100 =     DeltaMOperator[2](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs101 =     DeltaMOperator[2](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs102 =     DeltaDOperator[2](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs103 =     DeltaDOperator[2](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs104 =     DeltaMOperator[2](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs105 =     normalslave(0,0)*(-clhs100*clhs6 + clhs102*clhs2 + clhs103*clhs4 - clhs104*clhs7 + clhs5) + normalslave(0,1)*(clhs10*clhs103 - clhs100*clhs11 + clhs102*clhs9 - clhs104*clhs12);
        const double clhs106 =     clhs105*normalslave(0,0)*penalty_parameter;
        const double clhs107 =     DeltaDOperator[2](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs108 =     DeltaDOperator[2](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs109 =     DeltaMOperator[2](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs110 =     normalslave(1,0)*(-clhs101*clhs6 + clhs107*clhs2 + clhs108*clhs4 - clhs109*clhs7 + clhs18) + normalslave(1,1)*(clhs10*clhs108 - clhs101*clhs11 + clhs107*clhs9 - clhs109*clhs12);
        const double clhs111 =     clhs110*normalslave(1,0)*penalty_parameter;
        const double clhs112 =     DeltaMOperator[3](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs113 =     DeltaMOperator[3](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs114 =     DeltaDOperator[3](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs115 =     DeltaDOperator[3](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs116 =     DeltaMOperator[3](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs117 =     normalslave(0,0)*(-clhs112*clhs6 + clhs114*clhs2 + clhs115*clhs4 - clhs116*clhs7) + normalslave(0,1)*(clhs10*clhs115 - clhs11*clhs112 + clhs114*clhs9 - clhs116*clhs12 + clhs5);
        const double clhs118 =     clhs117*normalslave(0,0)*penalty_parameter;
        const double clhs119 =     DeltaDOperator[3](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs120 =     DeltaDOperator[3](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs121 =     DeltaMOperator[3](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs122 =     normalslave(1,0)*(-clhs113*clhs6 + clhs119*clhs2 + clhs120*clhs4 - clhs121*clhs7) + normalslave(1,1)*(clhs10*clhs120 - clhs11*clhs113 + clhs119*clhs9 - clhs12*clhs121 + clhs18);
        const double clhs123 =     clhs122*normalslave(1,0)*penalty_parameter;
        const double clhs124 =     normalslave(0,0)*scale_factor;
        const double clhs125 =     normalslave(1,0)*scale_factor;
        const double clhs126 =     clhs13*normalslave(0,1);
        const double clhs127 =     clhs20*normalslave(1,1);
        const double clhs128 =     clhs27*normalslave(0,1)*penalty_parameter;
        const double clhs129 =     clhs34*normalslave(1,1)*penalty_parameter;
        const double clhs130 =     clhs41*normalslave(0,1)*penalty_parameter;
        const double clhs131 =     clhs46*normalslave(1,1)*penalty_parameter;
        const double clhs132 =     clhs55*normalslave(0,1)*penalty_parameter;
        const double clhs133 =     clhs62*normalslave(1,1)*penalty_parameter;
        const double clhs134 =     clhs69*normalslave(0,1)*penalty_parameter;
        const double clhs135 =     clhs74*normalslave(1,1)*penalty_parameter;
        const double clhs136 =     clhs81*normalslave(0,1)*penalty_parameter;
        const double clhs137 =     clhs86*normalslave(1,1)*penalty_parameter;
        const double clhs138 =     clhs93*normalslave(0,1)*penalty_parameter;
        const double clhs139 =     clhs98*normalslave(1,1)*penalty_parameter;
        const double clhs140 =     clhs105*normalslave(0,1)*penalty_parameter;
        const double clhs141 =     clhs110*normalslave(1,1)*penalty_parameter;
        const double clhs142 =     clhs117*normalslave(0,1)*penalty_parameter;
        const double clhs143 =     clhs122*normalslave(1,1)*penalty_parameter;
        const double clhs144 =     normalslave(0,1)*scale_factor;
        const double clhs145 =     normalslave(1,1)*scale_factor;
        const double clhs146 =     clhs25 - clhs26;
        const double clhs147 =     clhs146*normalslave(0,0)*penalty_parameter;
        const double clhs148 =     clhs32 - clhs33;
        const double clhs149 =     clhs148*normalslave(1,0)*penalty_parameter;
        const double clhs150 =     clhs53 - clhs54;
        const double clhs151 =     clhs150*normalslave(0,0)*penalty_parameter;
        const double clhs152 =     clhs60 - clhs61;
        const double clhs153 =     clhs152*normalslave(1,0)*penalty_parameter;
        const double clhs154 =     clhs146*normalslave(0,1)*penalty_parameter;
        const double clhs155 =     clhs148*normalslave(1,1)*penalty_parameter;
        const double clhs156 =     clhs150*normalslave(0,1)*penalty_parameter;
        const double clhs157 =     clhs152*normalslave(1,1)*penalty_parameter;
    
        lhs(0,0)=-clhs0*clhs28 + clhs1*clhs14 - clhs15*clhs35 + clhs16*clhs21;
        lhs(0,1)=clhs0*clhs42 + clhs14*clhs36 + clhs15*clhs47 + clhs21*clhs37;
        lhs(0,2)=-clhs0*clhs56 + clhs14*clhs48 - clhs15*clhs63 + clhs21*clhs49;
        lhs(0,3)=clhs0*clhs70 + clhs14*clhs64 + clhs15*clhs75 + clhs21*clhs65;
        lhs(0,4)=clhs0*clhs82 + clhs14*clhs76 + clhs15*clhs87 + clhs21*clhs77;
        lhs(0,5)=clhs0*clhs94 + clhs14*clhs88 + clhs15*clhs99 + clhs21*clhs89;
        lhs(0,6)=clhs0*clhs106 + clhs100*clhs14 + clhs101*clhs21 + clhs111*clhs15;
        lhs(0,7)=clhs0*clhs118 + clhs112*clhs14 + clhs113*clhs21 + clhs123*clhs15;
        lhs(0,8)=clhs0*clhs124;
        lhs(0,9)=clhs125*clhs15;
        lhs(1,0)=-clhs0*clhs128 + clhs1*clhs126 + clhs127*clhs16 - clhs129*clhs15;
        lhs(1,1)=clhs0*clhs130 + clhs126*clhs36 + clhs127*clhs37 + clhs131*clhs15;
        lhs(1,2)=-clhs0*clhs132 + clhs126*clhs48 + clhs127*clhs49 - clhs133*clhs15;
        lhs(1,3)=clhs0*clhs134 + clhs126*clhs64 + clhs127*clhs65 + clhs135*clhs15;
        lhs(1,4)=clhs0*clhs136 + clhs126*clhs76 + clhs127*clhs77 + clhs137*clhs15;
        lhs(1,5)=clhs0*clhs138 + clhs126*clhs88 + clhs127*clhs89 + clhs139*clhs15;
        lhs(1,6)=clhs0*clhs140 + clhs100*clhs126 + clhs101*clhs127 + clhs141*clhs15;
        lhs(1,7)=clhs0*clhs142 + clhs112*clhs126 + clhs113*clhs127 + clhs143*clhs15;
        lhs(1,8)=clhs0*clhs144;
        lhs(1,9)=clhs145*clhs15;
        lhs(2,0)=clhs14*clhs24 - clhs19*clhs35 + clhs21*clhs31 - clhs28*clhs8;
        lhs(2,1)=clhs14*clhs40 + clhs19*clhs47 + clhs21*clhs45 + clhs42*clhs8;
        lhs(2,2)=clhs14*clhs52 - clhs19*clhs63 + clhs21*clhs59 - clhs56*clhs8;
        lhs(2,3)=clhs14*clhs68 + clhs19*clhs75 + clhs21*clhs73 + clhs70*clhs8;
        lhs(2,4)=clhs14*clhs80 + clhs19*clhs87 + clhs21*clhs85 + clhs8*clhs82;
        lhs(2,5)=clhs14*clhs92 + clhs19*clhs99 + clhs21*clhs97 + clhs8*clhs94;
        lhs(2,6)=clhs104*clhs14 + clhs106*clhs8 + clhs109*clhs21 + clhs111*clhs19;
        lhs(2,7)=clhs116*clhs14 + clhs118*clhs8 + clhs121*clhs21 + clhs123*clhs19;
        lhs(2,8)=clhs124*clhs8;
        lhs(2,9)=clhs125*clhs19;
        lhs(3,0)=clhs126*clhs24 + clhs127*clhs31 - clhs128*clhs8 - clhs129*clhs19;
        lhs(3,1)=clhs126*clhs40 + clhs127*clhs45 + clhs130*clhs8 + clhs131*clhs19;
        lhs(3,2)=clhs126*clhs52 + clhs127*clhs59 - clhs132*clhs8 - clhs133*clhs19;
        lhs(3,3)=clhs126*clhs68 + clhs127*clhs73 + clhs134*clhs8 + clhs135*clhs19;
        lhs(3,4)=clhs126*clhs80 + clhs127*clhs85 + clhs136*clhs8 + clhs137*clhs19;
        lhs(3,5)=clhs126*clhs92 + clhs127*clhs97 + clhs138*clhs8 + clhs139*clhs19;
        lhs(3,6)=clhs104*clhs126 + clhs109*clhs127 + clhs140*clhs8 + clhs141*clhs19;
        lhs(3,7)=clhs116*clhs126 + clhs121*clhs127 + clhs142*clhs8 + clhs143*clhs19;
        lhs(3,8)=clhs144*clhs8;
        lhs(3,9)=clhs145*clhs19;
        lhs(4,0)=-clhs14*clhs22 - clhs147*clhs3 - clhs149*clhs17 - clhs21*clhs29;
        lhs(4,1)=-clhs14*clhs38 - clhs17*clhs47 - clhs21*clhs43 - clhs3*clhs42;
        lhs(4,2)=-clhs14*clhs50 - clhs151*clhs3 - clhs153*clhs17 - clhs21*clhs57;
        lhs(4,3)=-clhs14*clhs66 - clhs17*clhs75 - clhs21*clhs71 - clhs3*clhs70;
        lhs(4,4)=-clhs14*clhs78 - clhs17*clhs87 - clhs21*clhs83 - clhs3*clhs82;
        lhs(4,5)=-clhs14*clhs90 - clhs17*clhs99 - clhs21*clhs95 - clhs3*clhs94;
        lhs(4,6)=-clhs102*clhs14 - clhs106*clhs3 - clhs107*clhs21 - clhs111*clhs17;
        lhs(4,7)=-clhs114*clhs14 - clhs118*clhs3 - clhs119*clhs21 - clhs123*clhs17;
        lhs(4,8)=-clhs124*clhs3;
        lhs(4,9)=-clhs125*clhs17;
        lhs(5,0)=-clhs126*clhs22 - clhs127*clhs29 - clhs154*clhs3 - clhs155*clhs17;
        lhs(5,1)=-clhs126*clhs38 - clhs127*clhs43 - clhs130*clhs3 - clhs131*clhs17;
        lhs(5,2)=-clhs126*clhs50 - clhs127*clhs57 - clhs156*clhs3 - clhs157*clhs17;
        lhs(5,3)=-clhs126*clhs66 - clhs127*clhs71 - clhs134*clhs3 - clhs135*clhs17;
        lhs(5,4)=-clhs126*clhs78 - clhs127*clhs83 - clhs136*clhs3 - clhs137*clhs17;
        lhs(5,5)=-clhs126*clhs90 - clhs127*clhs95 - clhs138*clhs3 - clhs139*clhs17;
        lhs(5,6)=-clhs102*clhs126 - clhs107*clhs127 - clhs140*clhs3 - clhs141*clhs17;
        lhs(5,7)=-clhs114*clhs126 - clhs119*clhs127 - clhs142*clhs3 - clhs143*clhs17;
        lhs(5,8)=-clhs144*clhs3;
        lhs(5,9)=-clhs145*clhs17;
        lhs(6,0)=-clhs14*clhs23 - clhs147*clhs5 - clhs149*clhs18 - clhs21*clhs30;
        lhs(6,1)=-clhs14*clhs39 - clhs18*clhs47 - clhs21*clhs44 - clhs42*clhs5;
        lhs(6,2)=-clhs14*clhs51 - clhs151*clhs5 - clhs153*clhs18 - clhs21*clhs58;
        lhs(6,3)=-clhs14*clhs67 - clhs18*clhs75 - clhs21*clhs72 - clhs5*clhs70;
        lhs(6,4)=-clhs14*clhs79 - clhs18*clhs87 - clhs21*clhs84 - clhs5*clhs82;
        lhs(6,5)=-clhs14*clhs91 - clhs18*clhs99 - clhs21*clhs96 - clhs5*clhs94;
        lhs(6,6)=-clhs103*clhs14 - clhs106*clhs5 - clhs108*clhs21 - clhs111*clhs18;
        lhs(6,7)=-clhs115*clhs14 - clhs118*clhs5 - clhs120*clhs21 - clhs123*clhs18;
        lhs(6,8)=-clhs124*clhs5;
        lhs(6,9)=-clhs125*clhs18;
        lhs(7,0)=-clhs126*clhs23 - clhs127*clhs30 - clhs154*clhs5 - clhs155*clhs18;
        lhs(7,1)=-clhs126*clhs39 - clhs127*clhs44 - clhs130*clhs5 - clhs131*clhs18;
        lhs(7,2)=-clhs126*clhs51 - clhs127*clhs58 - clhs156*clhs5 - clhs157*clhs18;
        lhs(7,3)=-clhs126*clhs67 - clhs127*clhs72 - clhs134*clhs5 - clhs135*clhs18;
        lhs(7,4)=-clhs126*clhs79 - clhs127*clhs84 - clhs136*clhs5 - clhs137*clhs18;
        lhs(7,5)=-clhs126*clhs91 - clhs127*clhs96 - clhs138*clhs5 - clhs139*clhs18;
        lhs(7,6)=-clhs103*clhs126 - clhs108*clhs127 - clhs140*clhs5 - clhs141*clhs18;
        lhs(7,7)=-clhs115*clhs126 - clhs120*clhs127 - clhs142*clhs5 - clhs143*clhs18;
        lhs(7,8)=-clhs144*clhs5;
        lhs(7,9)=-clhs145*clhs18;
        lhs(8,0)=clhs27*scale_factor;
        lhs(8,1)=-clhs41*scale_factor;
        lhs(8,2)=clhs55*scale_factor;
        lhs(8,3)=-clhs69*scale_factor;
        lhs(8,4)=-clhs81*scale_factor;
        lhs(8,5)=-clhs93*scale_factor;
        lhs(8,6)=-clhs105*scale_factor;
        lhs(8,7)=-clhs117*scale_factor;
        lhs(8,8)=0;
        lhs(8,9)=0;
        lhs(9,0)=clhs34*scale_factor;
        lhs(9,1)=-clhs46*scale_factor;
        lhs(9,2)=clhs62*scale_factor;
        lhs(9,3)=-clhs74*scale_factor;
        lhs(9,4)=-clhs86*scale_factor;
        lhs(9,5)=-clhs98*scale_factor;
        lhs(9,6)=-clhs110*scale_factor;
        lhs(9,7)=-clhs122*scale_factor;
        lhs(9,8)=0;
        lhs(9,9)=0;
    }


    return lhs;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
template<>
bounded_matrix<double, 21, 21> AugmentedLagrangianMethodFrictionlessMortarContactCondition<3,3>::CalculateLocalLHS<21>(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const unsigned int& rMasterElementIndex,
        const unsigned int& rActiveInactive
        )
{
    bounded_matrix<double,21,21> lhs;
    
    // Master segment info
    GeometryType& CurrentMasterElement = mThisMasterElements[rMasterElementIndex]->GetGeometry();

    // Initialize values
    const bounded_matrix<double, 3, 3> u1 = GetVariableMatrix(this->GetGeometry(), DISPLACEMENT, 0);
    const bounded_matrix<double, 3, 3> u2 = GetVariableMatrix(CurrentMasterElement, DISPLACEMENT, 0);
    const bounded_matrix<double, 3, 3> X1 = GetCoordinates(this->GetGeometry(), false);
    const bounded_matrix<double, 3, 3> X2 = GetCoordinates(CurrentMasterElement, false);
    
    const array_1d<double, 3> lmnormal = GetVariableVector(this->GetGeometry(), NORMAL_CONTACT_STRESS, 0); 
    
    const bounded_matrix<double, 3, 3> normalslave = GetVariableMatrix(this->GetGeometry(),  NORMAL); 
    
    // Augmentation parameters
    double scale_factor = 1.0;
    double penalty_parameter = 0.0;
    if (GetProperties().Has(SCALE_FACTOR) == true)
    {
        scale_factor  = GetProperties().GetValue(SCALE_FACTOR);
    }
    if (GetProperties().Has(PENALTY_FACTOR) == true)
    {
        penalty_parameter = GetProperties().GetValue(PENALTY_FACTOR);
    }
    
    // Mortar operators
    const bounded_matrix<double, 3, 3> MOperator = rMortarConditionMatrices.MOperator;
    const bounded_matrix<double, 3, 3> DOperator = rMortarConditionMatrices.DOperator;
    // Mortar operators derivatives
    const array_1d<bounded_matrix<double, 3, 3>, 23> DeltaMOperator = rMortarConditionMatrices.DeltaMOperator;
    const array_1d<bounded_matrix<double, 3, 3>, 23> DeltaDOperator = rMortarConditionMatrices.DeltaDOperator;
    
    return lhs;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
template<>
bounded_matrix<double, 28, 28> AugmentedLagrangianMethodFrictionlessMortarContactCondition<3,4>::CalculateLocalLHS<28>(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const unsigned int& rMasterElementIndex,
        const unsigned int& rActiveInactive
        )
{
    bounded_matrix<double,28,28> lhs;
    
    // Master segment info
    GeometryType& CurrentMasterElement = mThisMasterElements[rMasterElementIndex]->GetGeometry();

    // Initialize values
    const bounded_matrix<double, 4, 3> u1 = GetVariableMatrix(this->GetGeometry(), DISPLACEMENT, 0);
    const bounded_matrix<double, 4, 3> u2 = GetVariableMatrix(CurrentMasterElement, DISPLACEMENT, 0);
    const bounded_matrix<double, 4, 3> X1 = GetCoordinates(this->GetGeometry(), false);
    const bounded_matrix<double, 4, 3> X2 = GetCoordinates(CurrentMasterElement, false);
    
    const array_1d<double, 4> lmnormal = GetVariableVector(this->GetGeometry(), NORMAL_CONTACT_STRESS, 0); 
    
    const bounded_matrix<double, 4, 3> normalslave = GetVariableMatrix(this->GetGeometry(),  NORMAL); 
    
    // Augmentation parameters
    double scale_factor = 1.0;
    double penalty_parameter = 0.0;
    if (GetProperties().Has(SCALE_FACTOR) == true)
    {
        scale_factor  = GetProperties().GetValue(SCALE_FACTOR);
    }
    if (GetProperties().Has(PENALTY_FACTOR) == true)
    {
        penalty_parameter = GetProperties().GetValue(PENALTY_FACTOR);
    }
    
    // Mortar operators
    const bounded_matrix<double, 4, 4> MOperator = rMortarConditionMatrices.MOperator;
    const bounded_matrix<double, 4, 4> DOperator = rMortarConditionMatrices.DOperator;
    // Mortar operators derivatives
    const array_1d<bounded_matrix<double, 4, 4>, 24> DeltaMOperator = rMortarConditionMatrices.DeltaMOperator;
    const array_1d<bounded_matrix<double, 4, 4>, 24> DeltaDOperator = rMortarConditionMatrices.DeltaDOperator;
    
    return lhs;
}

/****************************** END AD REPLACEMENT *********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
template< unsigned int TMatrixSize >
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateAndAddRHS(
    LocalSystemComponents& rLocalSystem,
    const array_1d<double, TMatrixSize>& RHS_contact_pair,
    const unsigned int rPairIndex
    )
{   
    if ( rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
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
                this->AssembleContactPairRHSToConditionSystem<TMatrixSize>( RHS_contact_pair, rRightHandSideVector, rPairIndex );
                
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
        this->AssembleContactPairRHSToConditionSystem<TMatrixSize>( RHS_contact_pair, rRightHandSideVector, rPairIndex );
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
template< unsigned int TMatrixSize>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::AssembleContactPairRHSToConditionSystem(
    const array_1d<double, TMatrixSize>& rPairRHS,
    VectorType& rConditionRHS,
    const unsigned int rPairIndex
    )
{
    // Find location of the pair's master DOFs in ConditionRHS
    const unsigned int index_begin = rPairIndex * TMatrixSize;
    const unsigned int index_end  = index_begin + TMatrixSize;
    
    subrange( rConditionRHS, index_begin, index_end) += rPairRHS;
}

/***************************** BEGIN AD REPLACEMENT ********************************/
/***********************************************************************************/

template<>
template<>
array_1d<double, 10> AugmentedLagrangianMethodFrictionlessMortarContactCondition<2,2>::CalculateLocalRHS<10>(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const unsigned int& rMasterElementIndex,
        const unsigned int& rActiveInactive
        )
{
    array_1d<double,10> rhs;

    // Master segment info
    GeometryType& CurrentMasterElement = mThisMasterElements[rMasterElementIndex]->GetGeometry();

    // Initialize values
    const bounded_matrix<double, 2, 2> u1 = GetVariableMatrix(this->GetGeometry(), DISPLACEMENT, 0);
    const bounded_matrix<double, 2, 2> u2 = GetVariableMatrix(CurrentMasterElement, DISPLACEMENT, 0);
    const bounded_matrix<double, 2, 2> X1 = GetCoordinates(this->GetGeometry(), false);
    const bounded_matrix<double, 2, 2> X2 = GetCoordinates(CurrentMasterElement, false);
    
    const array_1d<double, 2> lmnormal = GetVariableVector(this->GetGeometry(), NORMAL_CONTACT_STRESS, 0); 
    
    const bounded_matrix<double, 2, 2> normalslave = GetVariableMatrix(this->GetGeometry(),  NORMAL); 
    
    // Augmentation parameters
    double scale_factor = 1.0;
    double penalty_parameter = 0.0;
    if (GetProperties().Has(SCALE_FACTOR) == true)
    {
        scale_factor  = GetProperties().GetValue(SCALE_FACTOR);
    }
    if (GetProperties().Has(PENALTY_FACTOR) == true)
    {
        penalty_parameter = GetProperties().GetValue(PENALTY_FACTOR);
    }
    
    // Mortar operators
    const bounded_matrix<double, 2, 2> MOperator = rMortarConditionMatrices.MOperator;
    const bounded_matrix<double, 2, 2> DOperator = rMortarConditionMatrices.DOperator;

    if (rActiveInactive == 0 )
    {
        const double crhs0 =     0.5*std::pow(scale_factor, 2.0)/penalty_parameter;
    
        rhs[0]=0;
        rhs[1]=0;
        rhs[2]=0;
        rhs[3]=0;
        rhs[4]=0;
        rhs[5]=0;
        rhs[6]=0;
        rhs[7]=0;
        rhs[8]=-crhs0*lmnormal[0];
        rhs[9]=-crhs0*lmnormal[1];
    }
    else if (rActiveInactive == 2 )
    {
        const double crhs0 =     MOperator(1,0); // MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs1 =     DOperator(1,0); // DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs2 =     DOperator(1,1); // DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs3 =     MOperator(1,1); // MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs4 =     normalslave(1,0)*(-crhs0*(X2(0,0) + u2(0,0)) + crhs1*(X1(0,0) + u1(0,0)) + crhs2*(X1(1,0) + u1(1,0)) - crhs3*(X2(1,0) + u2(1,0))) + normalslave(1,1)*(-crhs0*(X2(0,1) + u2(0,1)) + crhs1*(X1(0,1) + u1(0,1)) + crhs2*(X1(1,1) + u1(1,1)) - crhs3*(X2(1,1) + u2(1,1)));
        const double crhs5 =     crhs4*penalty_parameter + lmnormal[1]*scale_factor;
        const double crhs6 =     crhs5*normalslave(1,0);
        const double crhs7 =     crhs5*normalslave(1,1);
    
        rhs[0]=-crhs0*crhs6;
        rhs[1]=-crhs0*crhs7;
        rhs[2]=-crhs3*crhs6;
        rhs[3]=-crhs3*crhs7;
        rhs[4]=crhs1*crhs6;
        rhs[5]=crhs1*crhs7;
        rhs[6]=crhs2*crhs6;
        rhs[7]=crhs2*crhs7;
        rhs[8]=-0.5*lmnormal[0]*std::pow(scale_factor, 2.0)/penalty_parameter;
        rhs[9]=crhs4*scale_factor;
    }
    else if (rActiveInactive == 1 )
    {
        const double crhs0 =     MOperator(0,0); // MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs1 =     DOperator(0,0); // DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs2 =     DOperator(0,1); // DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs3 =     MOperator(0,1); // MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs4 =     normalslave(0,0)*(-crhs0*(X2(0,0) + u2(0,0)) + crhs1*(X1(0,0) + u1(0,0)) + crhs2*(X1(1,0) + u1(1,0)) - crhs3*(X2(1,0) + u2(1,0))) + normalslave(0,1)*(-crhs0*(X2(0,1) + u2(0,1)) + crhs1*(X1(0,1) + u1(0,1)) + crhs2*(X1(1,1) + u1(1,1)) - crhs3*(X2(1,1) + u2(1,1)));
        const double crhs5 =     crhs4*penalty_parameter + lmnormal[0]*scale_factor;
        const double crhs6 =     crhs5*normalslave(0,0);
        const double crhs7 =     crhs5*normalslave(0,1);
    
        rhs[0]=-crhs0*crhs6;
        rhs[1]=-crhs0*crhs7;
        rhs[2]=-crhs3*crhs6;
        rhs[3]=-crhs3*crhs7;
        rhs[4]=crhs1*crhs6;
        rhs[5]=crhs1*crhs7;
        rhs[6]=crhs2*crhs6;
        rhs[7]=crhs2*crhs7;
        rhs[8]=crhs4*scale_factor;
        rhs[9]=-0.5*lmnormal[1]*std::pow(scale_factor, 2.0)/penalty_parameter;
    }
    else if (rActiveInactive == 3 )
    {
        const double crhs0 =     MOperator(0,0); // MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs1 =     X1(0,0) + u1(0,0);
        const double crhs2 =     DOperator(0,0); // DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs3 =     X1(1,0) + u1(1,0);
        const double crhs4 =     DOperator(0,1); // DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs5 =     X2(0,0) + u2(0,0);
        const double crhs6 =     X2(1,0) + u2(1,0);
        const double crhs7 =     MOperator(0,1); // MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs8 =     X1(0,1) + u1(0,1);
        const double crhs9 =     X1(1,1) + u1(1,1);
        const double crhs10 =     X2(0,1) + u2(0,1);
        const double crhs11 =     X2(1,1) + u2(1,1);
        const double crhs12 =     normalslave(0,0)*(-crhs0*crhs5 + crhs1*crhs2 + crhs3*crhs4 - crhs6*crhs7) + normalslave(0,1)*(-crhs0*crhs10 - crhs11*crhs7 + crhs2*crhs8 + crhs4*crhs9);
        const double crhs13 =     crhs12*penalty_parameter + lmnormal[0]*scale_factor;
        const double crhs14 =     crhs13*normalslave(0,0);
        const double crhs15 =     MOperator(1,0); // MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs16 =     DOperator(1,0); // DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs17 =     DOperator(1,1); // DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs18 =     MOperator(1,1); // MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs19 =     normalslave(1,0)*(crhs1*crhs16 - crhs15*crhs5 + crhs17*crhs3 - crhs18*crhs6) + normalslave(1,1)*(-crhs10*crhs15 - crhs11*crhs18 + crhs16*crhs8 + crhs17*crhs9);
        const double crhs20 =     crhs19*penalty_parameter + lmnormal[1]*scale_factor;
        const double crhs21 =     crhs20*normalslave(1,0);
        const double crhs22 =     crhs13*normalslave(0,1);
        const double crhs23 =     crhs20*normalslave(1,1);
    
        rhs[0]=-crhs0*crhs14 - crhs15*crhs21;
        rhs[1]=-crhs0*crhs22 - crhs15*crhs23;
        rhs[2]=-crhs14*crhs7 - crhs18*crhs21;
        rhs[3]=-crhs18*crhs23 - crhs22*crhs7;
        rhs[4]=crhs14*crhs2 + crhs16*crhs21;
        rhs[5]=crhs16*crhs23 + crhs2*crhs22;
        rhs[6]=crhs14*crhs4 + crhs17*crhs21;
        rhs[7]=crhs17*crhs23 + crhs22*crhs4;
        rhs[8]=crhs12*scale_factor;
        rhs[9]=crhs19*scale_factor;
    }


    return rhs;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
template<>
array_1d<double, 21> AugmentedLagrangianMethodFrictionlessMortarContactCondition<3,3>::CalculateLocalRHS<21>(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const unsigned int& rMasterElementIndex,
        const unsigned int& rActiveInactive
        )
{
    array_1d<double,21> rhs;

    // Master segment info
    GeometryType& CurrentMasterElement = mThisMasterElements[rMasterElementIndex]->GetGeometry();

    // Initialize values
    const bounded_matrix<double, 3, 3> u1 = GetVariableMatrix(this->GetGeometry(), DISPLACEMENT, 0);
    const bounded_matrix<double, 3, 3> u2 = GetVariableMatrix(CurrentMasterElement, DISPLACEMENT, 0);
    const bounded_matrix<double, 3, 3> X1 = GetCoordinates(this->GetGeometry(), false);
    const bounded_matrix<double, 3, 3> X2 = GetCoordinates(CurrentMasterElement, false);
    
    const array_1d<double, 3> lmnormal = GetVariableVector(this->GetGeometry(), NORMAL_CONTACT_STRESS, 0); 
    
    const bounded_matrix<double, 3, 3> normalslave = GetVariableMatrix(this->GetGeometry(),  NORMAL); 
    
    // Augmentation parameters
    double scale_factor = 1.0;
    double penalty_parameter = 0.0;
    if (GetProperties().Has(SCALE_FACTOR) == true)
    {
        scale_factor  = GetProperties().GetValue(SCALE_FACTOR);
    }
    if (GetProperties().Has(PENALTY_FACTOR) == true)
    {
        penalty_parameter = GetProperties().GetValue(PENALTY_FACTOR);
    }
    
    // Mortar operators
    const bounded_matrix<double, 3, 3> MOperator = rMortarConditionMatrices.MOperator;
    const bounded_matrix<double, 3, 3> DOperator = rMortarConditionMatrices.DOperator;

    return rhs;
}
/***********************************************************************************/
/***********************************************************************************/

template<>
template<>
array_1d<double, 28> AugmentedLagrangianMethodFrictionlessMortarContactCondition<3,4>::CalculateLocalRHS<28>(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const unsigned int& rMasterElementIndex,
        const unsigned int& rActiveInactive
        )
{
    array_1d<double,28> rhs;

    // Master segment info
    GeometryType& CurrentMasterElement = mThisMasterElements[rMasterElementIndex]->GetGeometry();

    // Initialize values
    const bounded_matrix<double, 4, 3> u1 = GetVariableMatrix(this->GetGeometry(), DISPLACEMENT, 0);
    const bounded_matrix<double, 4, 3> u2 = GetVariableMatrix(CurrentMasterElement, DISPLACEMENT, 0);
    const bounded_matrix<double, 4, 3> X1 = GetCoordinates(this->GetGeometry(), false);
    const bounded_matrix<double, 4, 3> X2 = GetCoordinates(CurrentMasterElement, false);
    
    const array_1d<double, 4> lmnormal = GetVariableVector(this->GetGeometry(), NORMAL_CONTACT_STRESS, 0); 
    
    const bounded_matrix<double, 4, 3> normalslave = GetVariableMatrix(this->GetGeometry(),  NORMAL); 
    
    // Augmentation parameters
    double scale_factor = 1.0;
    double penalty_parameter = 0.0;
    if (GetProperties().Has(SCALE_FACTOR) == true)
    {
        scale_factor  = GetProperties().GetValue(SCALE_FACTOR);
    }
    if (GetProperties().Has(PENALTY_FACTOR) == true)
    {
        penalty_parameter = GetProperties().GetValue(PENALTY_FACTOR);
    }
    
    // Mortar operators
    const bounded_matrix<double, 4, 4> MOperator = rMortarConditionMatrices.MOperator;
    const bounded_matrix<double, 4, 4> DOperator = rMortarConditionMatrices.DOperator;

    return rhs;
}

/****************************** END AD REPLACEMENT *********************************/
/***********************************************************************************/

// TODO: Do an explicit specialization!!!!

template< unsigned int TDim, unsigned int TNumNodes> // NOTE: Formulation taken from Mohamed Khalil work
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateDeltaDetJSlave(
   GeneralVariables& rVariables,
   DerivativeData& rDerivativeData
   )
{
    if (TDim == 2)
    {
        // Fill up the elements corresponding to the slave DOFs - the rest remains zero
        for ( unsigned int i_slave = 0, i = 0; i_slave < TNumNodes; ++i_slave, i += TDim )
        {
            rDerivativeData.DeltaJ_s[i    ] = rVariables.j_Slave( 0, 0 ) * rVariables.DN_De_Slave( i_slave, 0) / rVariables.DetJSlave;
            rDerivativeData.DeltaJ_s[i + 1] = rVariables.j_Slave( 1, 0 ) * rVariables.DN_De_Slave( i_slave, 0) / rVariables.DetJSlave;
        }
    }
    else
    {
        const array_1d<double,TNumNodes>& DN_Dxi  = column( rVariables.DN_De_Slave, 0 );
        const array_1d<double,TNumNodes>& DN_Deta = column( rVariables.DN_De_Slave, 1 );
        
        const array_1d<double,TDim>& J_xi    = column( rVariables.j_Slave, 0 );
        const array_1d<double,TDim>& J_eta   = column( rVariables.j_Slave, 1 );
        
        const array_1d<double,TDim>& n = prod(trans(rDerivativeData.Normal_m), rVariables.N_Slave); // FIXME: Check this!!!!
        
        bounded_matrix<double, TDim, TDim> Delta_Jxi_x_Jeta;
        
        for ( unsigned int i_slave = 0, i = 0; i_slave < TNumNodes; ++i_slave, i += TDim )
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
            
            rDerivativeData.DeltaJ_s[i    ] = inner_prod( n, column( Delta_Jxi_x_Jeta, 0 ) );
            rDerivativeData.DeltaJ_s[i + 1] = inner_prod( n, column( Delta_Jxi_x_Jeta, 1 ) );
            rDerivativeData.DeltaJ_s[i + 2] = inner_prod( n, column( Delta_Jxi_x_Jeta, 2 ) );
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes> // NOTE: Formulation taken from Mohamed Khalil work
bounded_matrix<double, TDim, TDim> AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::LocalDeltaNormal( // NOTE: Not the mean, look in the contact utilities 
    const GeometryType& CondGeometry,
    const unsigned int node_index
    )
{
    // Tolerance
    const double tol = 1.0e-14;
        
    bounded_matrix<double, TDim, TDim> Delta_ne_adj;
    bounded_matrix<double, TDim, TDim> Ce;
    
    const bounded_matrix<double, TDim, TDim> I = IdentityMatrix(TDim, TDim);
    
    bounded_matrix<double, TDim, TDim> DeltaNormal = ZeroMatrix(TDim,TDim);
    
    const array_1d<double, 3> ne = this->GetValue(NORMAL);   // normalized condition normal
    bounded_matrix<double, TDim, TDim> ne_o_ne = subrange( outer_prod( ne, ne ), 0, TDim, 0, TDim );
    
    if (TDim == 2)
    {
        const double ne_norm = this->GetGeometry( ).Length( ); // The norm of a geometry's normal is its characteristic dimension - length for 2D and area for 3D 
                        
        Delta_ne_adj( 0, 0 ) =  0.0;
        Delta_ne_adj( 0, 1 ) = -1.0;
        Delta_ne_adj( 1, 0 ) =  1.0;
        Delta_ne_adj( 1, 1 ) =  0.0;
        
        double DN_De_j = 0.0;
        if( node_index == 0 )
        {
            DN_De_j = - 0.5;
        }
        else if( node_index == 1 )
        {
            DN_De_j =   0.5;
        }
        
        Ce = prod( I - ne_o_ne, Delta_ne_adj ) / ne_norm;     // In 2D, Delta_ne_adj is node-independent => evaluated outside the nodes loop
        
        DeltaNormal = - 2.0 * Ce * DN_De_j; // NOTE: Check why - 2???!!!, it was the only wayto ensure the same value as the symbolic. You will need to repeat this in 3D            
//         DeltaNormal = Ce * DN_De_j;             
    }
    else
    {
        const double ne_norm = this->GetGeometry( ).Area( ); // The norm of a geometry's normal is its characteristic dimension - length for 2D and area for 3D 
        
        MatrixType J = ZeroMatrix( 3, 2 ); // Jacobian [ 3D global x 2D local ]
        array_1d<double, 2> DN_De_j;
        array_1d<double, 3> local_coords_j;
        
        if( TNumNodes == 3 )    // linear triangle element
        {
            if( node_index == 0 )
            {
                local_coords_j[0] = 0.0;
                local_coords_j[1] = 0.0;
                DN_De_j[0] = - 1.0;
                DN_De_j[1] = - 1.0;
            }
            else if( node_index == 1 )
            {
                local_coords_j[0] = 1.0;
                local_coords_j[1] = 0.0;
                DN_De_j[0] = 1.0;
                DN_De_j[1] = 0.0;
            }
            else // node_index == 2
            {
                local_coords_j[0] = 0.0;
                local_coords_j[1] = 1.0;
                DN_De_j[0] = 0.0;
                DN_De_j[1] = 1.0;
            }
        }
        else if( TNumNodes == 4 )    // linear quad element 
        {
            if( node_index == 0 )
            {
                local_coords_j[0] = - 1.0;
                local_coords_j[1] = - 1.0;
                DN_De_j[0] = - 0.5;
                DN_De_j[1] = - 0.5;
            }
            else if( node_index == 1 )
            {
                local_coords_j[0] =   1.0;
                local_coords_j[1] = - 1.0;
                DN_De_j[0] =   0.5;
                DN_De_j[1] = - 0.5;
            }
            else if( node_index == 2 )
            {
                local_coords_j[0] =  1.0;
                local_coords_j[1] =  1.0;
                DN_De_j[0] =  0.5;
                DN_De_j[1] =  0.5;
            }
            else // node_index == 3
            {
                local_coords_j[0] = - 1.0;
                local_coords_j[1] =   1.0;
                DN_De_j[0] = - 0.5;
                DN_De_j[1] =   0.5;
            }
        }
        
        this->GetGeometry( ).Jacobian( J, local_coords_j );
        
        Delta_ne_adj(0,0) = 0.0;
        Delta_ne_adj(0,1) = +J(2,1) * DN_De_j[0] - J(2,0) * DN_De_j[1]; 
        Delta_ne_adj(0,2) = -J(1,1) * DN_De_j[0] + J(1,0) * DN_De_j[1]; 
        Delta_ne_adj(1,0) = -J(2,1) * DN_De_j[0] + J(2,0) * DN_De_j[1]; 
        Delta_ne_adj(1,1) = 0.0;                   
        Delta_ne_adj(1,2) = +J(0,1) * DN_De_j[0] - J(0,0) * DN_De_j[1]; 
        Delta_ne_adj(2,0) = +J(1,1) * DN_De_j[0] - J(1,0) * DN_De_j[1]; 
        Delta_ne_adj(2,1) = -J(0,1) * DN_De_j[0] + J(0,0) * DN_De_j[1]; 
        Delta_ne_adj(2,2) = 0.0;
        
        Ce = prod( I - ne_o_ne, Delta_ne_adj ) / ne_norm;
        DeltaNormal = Ce;
    }
    
    const double ne_norm = norm_2( ne );
    const double ne_norm_3 = ne_norm * ne_norm * ne_norm;
    
    if ( ne_norm_3 > tol )
    {
        const bounded_matrix<double, TDim, TDim> Cj = I / ne_norm - ne_o_ne / ne_norm_3;
        DeltaNormal = prod( Cj, DeltaNormal );
    }
        
    return DeltaNormal; 
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes> // NOTE: Formulation taken from Mohamed Khalil work
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateDeltaNormalSlave(DerivativeData& rDerivativeData)
{
    if (TDim == 2) // TODO: Use explicit 
    {
        for ( unsigned int i_slave = 0, i = 0; i_slave < TNumNodes; ++i_slave, i += TDim )
        {
//             bounded_matrix<double, 2, 2> DeltaNormal = GetGeometry()[i_slave].GetValue(DELTA_NORMAL);
            bounded_matrix<double, 2, 2> DeltaNormal = this->LocalDeltaNormal(GetGeometry(), i_slave);
            for (unsigned i_dof = 0; i_dof < TDim; i_dof++) 
            {
                for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
                {
                    row(rDerivativeData.Delta_Normal_s[i_slave * TDim + i_dof], i_node)      =   trans(column(DeltaNormal, i_dof)); 
                }
            }
        }
    }
    else
    {
        for ( unsigned int i_slave = 0, i = 0; i_slave < TNumNodes; ++i_slave, i += TDim )
        {
//             bounded_matrix<double, 3, 3> DeltaNormal = GetGeometry()[i_slave]->GetValue(DELTA_NORMAL);
            const bounded_matrix<double, 3, 3> DeltaNormal = this->LocalDeltaNormal(this->GetGeometry(), i_slave);
            
            // Calculate nodal tangents
            
            const MatrixType I = IdentityMatrix( TDim, TDim );
            
            array_1d<double, 2> DN_De_j;
            if( TNumNodes == 3 )    // linear triangle element // TODO: Use an enum
            {
                if( i_slave == 0 )
                {
                    DN_De_j[0] = -1.0;
                    DN_De_j[1] = -1.0;
                }
                else if( i_slave == 1 )
                {
                    DN_De_j[0] = 1.0;
                    DN_De_j[1] = 0.0;
                }
                else // i_slave == 2
                {
                    DN_De_j[0] = 0.0;
                    DN_De_j[1] = 1.0;
                }
            }
            else if( TNumNodes == 4 )    // linear quad element 
            {
                if( i_slave == 0 )
                {
                    DN_De_j[0] = -0.5;
                    DN_De_j[1] = -0.5;
                }
                else if( i_slave == 1 )
                {
                    DN_De_j[0] =  0.5;
                    DN_De_j[1] = -0.5;
                }
                else if( i_slave == 2 )
                {
                    DN_De_j[0] =  0.5;
                    DN_De_j[1] =  0.5;
                }
                else // i_slave == 3
                {
                    DN_De_j[0] = -0.5;
                    DN_De_j[1] =  0.5;
                }
            }
            
            for (unsigned i_dof = 0; i_dof < TDim; i_dof++) 
            {
                for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
                {
                    row(rDerivativeData.Delta_Normal_s[i_slave * TDim + i_dof], i_node) = trans(column(DeltaNormal, i_dof)); 
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes> // NOTE: Formulation taken from Mohamed Khalil work
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateDeltaNormalMaster(DerivativeData& rDerivativeData)
{
    if (TDim == 2)
    {
        for ( unsigned int i_master = 0, i = 0; i_master < TNumNodes; ++i_master, i += TDim )
        {
//             bounded_matrix<double, 2, 2> DeltaNormal = GetGeometry[i_master].GetValue(DELTA_NORMAL);
            bounded_matrix<double, 2, 2> DeltaNormal = this->LocalDeltaNormal(rDerivativeData.MasterGeometry, i_master);
            for (unsigned i_dof = 0; i_dof < TDim; i_dof++) 
            {
                for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
                {
                    row(rDerivativeData.Delta_Normal_m[i_master * TDim + i_dof], i_node) = trans(column(DeltaNormal, i_dof)); 
                }
            }
        }
    }
    else
    {
        for ( unsigned int i_master = 0, i = 0; i_master < TNumNodes; ++i_master, i += TDim )
        {
//             bounded_matrix<double, 3, 3> DeltaNormal = GetGeometry[i_master]->GetValue(DELTA_NORMAL);
            const bounded_matrix<double, 3, 3> DeltaNormal = this->LocalDeltaNormal(rDerivativeData.MasterGeometry, i_master);
            for (unsigned i_dof = 0; i_dof < TDim; i_dof++) 
            {
                for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
                {
                    row(rDerivativeData.Delta_Normal_m[i_master * TDim + i_dof], i_node)  = trans(column(DeltaNormal, i_dof)); 
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateDeltaN(
   GeneralVariables& rVariables,
   DerivativeData& rDerivativeData
   )
{
    static const unsigned int Size1 = (TNumNodes * TDim);
    
    // Shape functions
    const VectorType  N1 = rVariables.N_Slave;
    const VectorType  N2 = rVariables.N_Master;
    
    // Coordinates
    const bounded_matrix<double, TNumNodes, TDim> u1 = rDerivativeData.u1;
    const bounded_matrix<double, TNumNodes, TDim> X1 = rDerivativeData.X1;
    const bounded_matrix<double, TNumNodes, TDim> u2 = rDerivativeData.u2;
    const bounded_matrix<double, TNumNodes, TDim> X2 = rDerivativeData.X2;
    
    // Normals
    const array_1d<double, TDim >  Normal_sg = prod(trans(rDerivativeData.Normal_s), N1);
    const array_1d<double, TDim >  Normal_mg = prod(trans(rDerivativeData.Normal_m), N2);
    
    const array_1d<bounded_matrix<double, TNumNodes, TDim>, Size1> DNormal_s = rDerivativeData.Delta_Normal_s;
    const array_1d<bounded_matrix<double, TNumNodes, TDim>, Size1> DNormal_m = rDerivativeData.Delta_Normal_m;
    
    bool compute = false;
    if (TDim == 2)
    {
       if (TNumNodes == 2)
       {
          compute = true;
       }
    }
    else
    {
      if (TNumNodes == 3)
       {
          compute = true;
       }
    }
    
   /* Calculate Delta N */ // TODO: Do the same for the N1
   if (compute == true)
   {
      const double tol = 1.0e-18;
      
//       const array_1d<double, TNumNodes > vector_nodes = trans(row(X2 + u2, 0)) - prod(trans(X1 + u1), N1); // NOTE: This is the way I considered in the symbolic
      const array_1d<double, TNumNodes > vector_nodes =  prod(trans(X2 + u2), N2) - prod(trans(X1 + u1), N1);
      const double Dist = inner_prod(vector_nodes, Normal_mg)/(inner_prod(Normal_sg, Normal_mg) + tol);
      
      double div1,div2 = 1.0;
      double mult1,mult2,mult3,mult4,mult5 = 0.0;
      if (TDim == 2)
      {
         if (TNumNodes == 2)
         {
//             div1 = (X1(0,0) + X1(0,1) - X1(1,0) - X1(1,1) + u1(0,0) + u1(0,1) - u1(1,0) - u1(1,1)) + tol; // NOTE: You will need to compute DeltaN1
            div2 = (X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + tol;
         }
      }
      else
      {
         if (TNumNodes == 3)
         {
             div1 = ((u2(0, 0) + u2(0, 2) - u2(1, 0) - u2(1, 2) + X2(0, 0) + 
                      X2(0, 2) - X2(1, 0) - X2(1, 2)) * (u2(0, 0) + u2(0, 1) - 
                      u2(2, 0) - u2(2, 1) + X2(0, 0) + X2(0, 1) - X2(2, 0) - 
                      X2(2, 1)) - (u2(0, 0) + u2(0, 1) - u2(1, 0) - u2(1, 1) + 
                      X2(0, 0) + X2(0, 1) - X2(1, 0) - X2(1, 1)) * (u2(0, 0) + 
                      u2(0, 2) - u2(2, 0) - u2(2, 2) + X2(0, 0) + X2(0, 2) - 
                      X2(2, 0) - X2(2, 2))) + tol;
             div2 = (-(u2(1, 1) + X2(1, 1)) * (u2(2, 0) + X2(2, 0)) + (u2(1, 2) + 
                       X2(1, 2)) * (u2(2, 0) + X2(2, 0)) + (u2(0, 2) + 
                       X2(0, 2)) * (u2(1, 0) + u2(1, 1) - u2(2, 0) - u2(2, 1) + 
                       X2(1, 0) + X2(1, 1) - X2(2, 0) - X2(2, 1)) + (u2(1, 0) + 
                       X2(1, 0)) * (u2(2, 1) + X2(2, 1)) + (u2(1, 2) + 
                       X2(1, 2)) * (u2(2, 1) + X2(2, 1)) - (u2(1, 0) + u2(1, 1) + 
                       X2(1, 0) + X2(1, 1)) * (u2(2, 2) + X2(2, 2)) + (u2(0, 1) + 
                       X2(0, 1)) * (-u2(1, 0) - u2(1, 2) + u2(2, 0) + u2(2, 2) - 
                       X2(1, 0) - X2(1, 2) + X2(2, 0) + X2(2, 2)) + (u2(0, 0) + 
                       X2(0, 0)) * (u2(1, 1) - u2(1, 2) - u2(2, 1) + u2(2, 2) + 
                       X2(1, 1) - X2(1, 2) - X2(2, 1) + X2(2, 2))) + tol;
                        
             mult1 = (-(u2(0, 0) + u2(0, 2) - u2(2, 0) - u2(2, 2) + X2(0, 0) + X2(0, 2) - X2(2, 0) - X2(2, 2)));
             mult2 = ( (u2(0, 0) + u2(0, 1) - u2(2, 0) - u2(2, 1) + X2(0, 0) + X2(0, 1) - X2(2, 0) - X2(2, 1)));
             mult3 = ( (u2(0, 1) - u2(0, 2) - u2(1, 1) + u2(1, 2) + X2(0, 1) - X2(0, 2) - X2(1, 1) + X2(1, 2)));
             mult4 = ((-u2(0, 0) - u2(0, 2) + u2(1, 0) + u2(1, 2) - X2(0, 0) - X2(0, 2) + X2(1, 0) + X2(1, 2)));
             mult5 = ( (u2(0, 0) + u2(0, 1) - u2(1, 0) - u2(1, 1) + X2(0, 0) + X2(0, 1) - X2(1, 0) - X2(1, 1)));
         }
      }
      
      // Derivatives slave
      for (unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++)
      {
         for (unsigned int i_dim = 0; i_dim < TDim; i_dim++)
         {
            const unsigned int i_dof = i_slave * TDim + i_dim;
            
            const array_1d<double, TNumNodes > DNormal_sg = prod(trans(DNormal_s[i_dof]), N1);
            array_1d<double, TNumNodes > aux_vector = ZeroVector(TNumNodes);
            aux_vector[i_dim] = 1.0;
            const array_1d<double, TNumNodes > Deltavector_nodes = - N1[i_slave] * aux_vector;
            const double DeltaDist = ((inner_prod(Deltavector_nodes, Normal_mg))* inner_prod(Normal_sg, Normal_mg) - inner_prod(vector_nodes, Normal_mg) * (inner_prod(DNormal_sg, Normal_mg)))/std::pow(inner_prod(Normal_sg, Normal_mg) + tol, 2);
            const array_1d<double, TNumNodes > Deltax1g = N1[i_slave] * aux_vector;
            const array_1d<double, TNumNodes > Deltax2g = Deltax1g + DeltaDist * Normal_sg + Dist * DNormal_sg;
            
            if (TDim == 2)
            {
               if (TNumNodes == 2)
               {
                  rDerivativeData.DeltaN2[i_dof][0] =  (Deltax2g[0] + Deltax2g[1])/div2;
                  rDerivativeData.DeltaN2[i_dof][1] =  - rDerivativeData.DeltaN2[i_dof][0];
               }
            }
            else
            {
               if (TNumNodes == 3)
               {
                  rDerivativeData.DeltaN2[i_dof][1] = - (mult1 * (Deltax2g[0] + Deltax2g[1]) + mult2 * (Deltax2g[0] + Deltax2g[2]))/div1;
                  rDerivativeData.DeltaN2[i_dof][2] =   (mult3 *  Deltax2g[0] + mult4 * Deltax2g[1] + mult5 * Deltax2g[2])/div2;
                  rDerivativeData.DeltaN2[i_dof][0] =  - rDerivativeData.DeltaN2[i_dof][1] - rDerivativeData.DeltaN2[i_dof][2];
               }
            }
         }
      }
      
      // Derivatives master
      for (unsigned int i_master = 0; i_master < TNumNodes; i_master++)
      {
         for (unsigned int i_dim = 0; i_dim < TDim; i_dim++)
         {
            const unsigned int i_dof = (TNumNodes + i_master) * TDim + i_dim;
            
            const array_1d<double, TNumNodes > DNormal_mg = prod(trans(DNormal_m[i_dof - TNumNodes * TDim]), N2);
            array_1d<double, TNumNodes > aux_vector = ZeroVector(TNumNodes);
//             if (i_master == 0) // NOTE: This is the way I considered in the symbolic
//             {
                aux_vector[i_dim] = 1.0;
//             }
            const array_1d<double, TNumNodes > Deltavector_nodes = N2[i_master] * aux_vector;
            const double DeltaDist = ((inner_prod(Deltavector_nodes, Normal_mg) + inner_prod(vector_nodes, DNormal_mg))* inner_prod(Normal_sg, Normal_mg) - inner_prod(vector_nodes, Normal_mg) * (inner_prod(Normal_sg, DNormal_mg)))/std::pow(inner_prod(Normal_sg, Normal_mg) + tol, 2);;
            const array_1d<double, TNumNodes > x2g = prod(trans(X2 + u2), N2);
            const array_1d<double, TNumNodes > Deltax2g = DeltaDist * Normal_sg;
            
            if (TDim == 2)
            {
               if (TNumNodes == 2)
               {
                   if (i_master == 0)
                   {
                       rDerivativeData.DeltaN2[i_dof][0] =  ( ((u2(1,0) + X2(1,0) + u2(1,1) + X2(1,1)) - x2g[0] - x2g[1]) + div2 * (Deltax2g[0] + Deltax2g[1]))/std::pow(div2, 2);
                   }
                   else 
                   {
                       rDerivativeData.DeltaN2[i_dof][0] =  ((-(u2(0,0) + X2(0,0) + u2(0,1) + X2(0,1)) + x2g[0] + x2g[1]) + div2 * (Deltax2g[0] + Deltax2g[1]))/std::pow(div2, 2);
                   }
                  
                  rDerivativeData.DeltaN2[i_dof][1] =  - rDerivativeData.DeltaN2[i_dof][0];
               }
            }
            else
            {
               if (TNumNodes == 3)
               {    
                   const double multmaster0 = ( (u2(0, 0) + u2(0, 2) - u2(2, 0) - u2(2, 2) + X2(0, 0) + X2(0, 2) - X2(2, 0) - X2(2, 2)));
                   const double multmaster1 = ((-u2(0, 0) - u2(0, 1) + u2(2, 0) + u2(2, 1) - X2(0, 0) - X2(0, 1) + X2(2, 0) + X2(2, 1)) );
                   const double multmaster2 = ( (u2(0, 0) + u2(0, 1) + X2(0, 0) + X2(0, 1) - x2g[0] - x2g[1]));
                   const double multmaster3 = ( (u2(0, 0) + u2(0, 2) + X2(0, 0) + X2(0, 2) - x2g[0] - x2g[2]));
                   const double multmaster4 = ( (u2(1, 0) + u2(1, 1) + X2(1, 0) + X2(1, 1)) * x2g[2] + (u2(0, 0) + X2(0, 0)));
                   const double multmaster5 = ((u2(1, 0) + u2(1, 1) + X2(1, 0) + X2(1, 1) - x2g[0] - x2g[1]));
                   const double multmaster6 = ( (u2(1, 1) - u2(1, 2) + X2(1, 1) - X2(1, 2) - x2g[1] + x2g[2]));
                   const double coefmaster1 = ( (u2(1, 0) + X2(1, 0)) * x2g[1] + (u2(1, 2) + X2(1, 2)) * x2g[1]);
                   const double coefmaster2 = ( -(u2(1, 1) + X2(1, 1)) * x2g[0] + (u2(1, 2) + X2(1, 2)) * x2g[0]);
                   const double coefmaster3 = ( u2(0, 0) + u2(0, 2) + X2(0, 0) + X2(0, 2) - x2g[0] - x2g[2]);
                   const double coefmaster4 = ( -u2(0, 0) - u2(0, 1) - X2(0, 0) - X2(0, 1) + x2g[0] + x2g[1]);
                   const double coefmaster5 = ( (u2(0, 1) + X2(0, 1)) * (u2(1, 0) + u2(1, 2) + X2(1, 0) + X2(1, 2) - x2g[0] - x2g[2]));
                   
                   if ((i_dof - TDim * TNumNodes) == 0)
                   {
                       rDerivativeData.DeltaN2[i_dof][1] =  ((u2(1, 1) - u2(1, 2) - u2(2, 1) + u2(2, 2) + X2(1, 1) - 
                                                           X2(1, 2) - X2(2, 1) + X2(2, 2)) * (multmaster0*multmaster2 + 
                                                           multmaster1*multmaster3) - (div1) * (u2(0, 1) - u2(0, 2) + 
                                                           X2(0, 1) - X2(0, 2) - x2g[1] + x2g[2] + 
                                                           mult1 * (-1 + Deltax2g[0] + Deltax2g[1]) + 
                                                           mult2 * (-1 + Deltax2g[0] + Deltax2g[2])) )/std::pow(div1, 2);
                                                           
                       rDerivativeData.DeltaN2[i_dof][2] =  ( -(u2(1, 1) - u2(1, 2) - u2(2, 1) + u2(2, 2) + X2(1, 1) - 
                                                            X2(1, 2) - X2(2, 1) + X2(2, 2)) * (coefmaster2 + (u2(0, 2) + X2(0, 2)) * 
                                                            multmaster5 + coefmaster1 - coefmaster5 - 
                                                            multmaster4 * multmaster6) + (div2) * (u2(1, 1) - u2(1, 2) + 
                                                            X2(1, 1) - X2(1, 2) - x2g[1] + x2g[2] + mult3 * Deltax2g[0] + 
                                                            mult4 * Deltax2g[1] + mult5 * Deltax2g[2]))/std::pow(div2, 2);
                   }
                   else if ((i_dof - TDim * TNumNodes) == 1)
                   {
                       rDerivativeData.DeltaN2[i_dof][1] =  ((-u2(1, 0) - u2(1, 2) + u2(2, 0) + u2(2, 2) - X2(1, 0) - 
                                                            X2(1, 2) + X2(2, 0) + X2(2, 2)) * (multmaster0*multmaster2 + 
                                                            multmaster1*multmaster3) - (div1) * (-u2(0, 0) - u2(0, 2) - 
                                                            X2(0, 0) - X2(0, 2) + x2g[0] + x2g[2] + 
                                                            mult1 * (-1 + Deltax2g[0] + Deltax2g[1]) + 
                                                            mult2 * (Deltax2g[0] + Deltax2g[2])) )/std::pow(div1, 2);
                                                            
                       rDerivativeData.DeltaN2[i_dof][2] =  (-(-u2(1, 0) - u2(1, 2) + u2(2, 0) + u2(2, 2) - X2(1, 0) - 
                                                             X2(1, 2) + X2(2, 0) + X2(2, 2)) * (coefmaster2 + (u2(0, 2) + X2(0, 2)) * 
                                                             multmaster5 + coefmaster1 - coefmaster5 - 
                                                             multmaster4 * multmaster6) + (div2) * (-u2(1, 0) - u2(1, 2) - 
                                                             X2(1, 0) - X2(1, 2) + x2g[0] + x2g[2] + mult3 * Deltax2g[0] + 
                                                             mult4 * Deltax2g[1] + mult5 * Deltax2g[2]) )/std::pow(div2, 2);
                   }
                   else if ((i_dof - TDim * TNumNodes) == 2)
                   {
                       rDerivativeData.DeltaN2[i_dof][1] =  ( (u2(1, 0) + u2(1, 1) - u2(2, 0) - u2(2, 1) + X2(1, 0) + 
                                                            X2(1, 1) - X2(2, 0) - X2(2, 1)) * (multmaster0*multmaster2 + 
                                                            multmaster1*multmaster3) - (div1) * (u2(0, 0) + u2(0, 1) + 
                                                            X2(0, 0) + X2(0, 1) - x2g[0] - x2g[1] + 
                                                            mult1 * (Deltax2g[0] + Deltax2g[1]) + 
                                                            mult2 * (-1 + Deltax2g[0] + Deltax2g[2])))/std::pow(div1, 2);
                        
                       rDerivativeData.DeltaN2[i_dof][2] =  (-(u2(1, 0) + u2(1, 1) - u2(2, 0) - u2(2, 1) + X2(1, 0) + 
                                                            X2(1, 1) - X2(2, 0) - X2(2, 1)) * (coefmaster2 + (u2(0, 2) + X2(0, 2)) * 
                                                            multmaster5 + coefmaster1 - coefmaster5 - 
                                                            multmaster4 * multmaster6) + (div2) * (u2(1, 0) + u2(1, 1) + 
                                                            X2(1, 0) + X2(1, 1) - x2g[0] - x2g[1] + mult3 * Deltax2g[0] + 
                                                            mult4 * Deltax2g[1] + mult5 * Deltax2g[2]) )/std::pow(div2, 2);
                   }
                   else if ((i_dof - TDim * TNumNodes) == 3)
                   {
                       rDerivativeData.DeltaN2[i_dof][1] =  ( (-u2(0, 1) + u2(0, 2) + u2(2, 1) - u2(2, 2) - X2(0, 1) + 
                                                            X2(0, 2) + X2(2, 1) - X2(2, 2)) * (multmaster0*multmaster2 + 
                                                            multmaster1 * multmaster3) - (div1) * (+mult1 * (Deltax2g[0] + Deltax2g[1]) + 
                                                            mult2 * (Deltax2g[0] + Deltax2g[2])))/std::pow(div1, 2);
                                                            
                       rDerivativeData.DeltaN2[i_dof][2] =  ( -(-u2(0, 1) + u2(0, 2) + u2(2, 1) - u2(2, 2) - X2(0, 1) + 
                                                              X2(0, 2) + X2(2, 1) - X2(2, 2)) * (coefmaster2 + (u2(0, 2) + X2(0, 2)) * 
                                                              multmaster5 + coefmaster1 - coefmaster5 - 
                                                              multmaster4 * multmaster6) + (div2) * (-u2(0, 1) + u2(0, 2) - 
                                                              X2(0, 1) + X2(0, 2) + x2g[1] - x2g[2] + mult3 * Deltax2g[0] + 
                                                              mult4 * Deltax2g[1] + mult5 * Deltax2g[2]))/std::pow(div2, 2);
                   }
                   else if ((i_dof - TDim * TNumNodes) == 4)
                   {
                       rDerivativeData.DeltaN2[i_dof][1] =  ( multmaster0*(multmaster0*multmaster2 + 
                                                           multmaster1 * multmaster3) - (div1) * (+mult1 * (Deltax2g[0] + Deltax2g[1]) + 
                                                           mult2 * (Deltax2g[0] + Deltax2g[2])))/std::pow(div1, 2);
                                                           
                       rDerivativeData.DeltaN2[i_dof][2] =  ( +mult1 * (coefmaster2 + (u2(0, 2) + X2(0, 2)) * multmaster5 + 
                                                            coefmaster1 - coefmaster5 - multmaster4 * multmaster6) + (div2) * (coefmaster3 + 
                                                            mult3 * Deltax2g[0] + mult4 * Deltax2g[1] + mult5 * Deltax2g[2]))/std::pow(div2, 2);
                   }
                   else if ((i_dof - TDim * TNumNodes) == 5)
                   {
                       rDerivativeData.DeltaN2[i_dof][1] =  ( multmaster1 *(multmaster0*multmaster2 + 
                                                           multmaster1 * multmaster3) - (div1) * (+mult1 * (Deltax2g[0] + Deltax2g[1]) + 
                                                           mult2 * (Deltax2g[0] + Deltax2g[2])))/std::pow(div1, 2);
                                                           
                       rDerivativeData.DeltaN2[i_dof][2] =  ( -multmaster1 *  (coefmaster2 + (u2(0, 2) + X2(0, 2)) * multmaster5 +
                                                            coefmaster1 - coefmaster5 -  multmaster4 * multmaster6) + (div2) * (coefmaster4 + 
                                                            mult3 * Deltax2g[0] + mult4 * Deltax2g[1] + mult5 * Deltax2g[2]))/std::pow(div2, 2);
                   }
                   else if ((i_dof - TDim * TNumNodes) == 6)
                   {
                       rDerivativeData.DeltaN2[i_dof][1] =  (mult3 *(multmaster0*multmaster2 + 
                                                          multmaster1*multmaster3) - (div1) * (-u2(0, 1) + u2(0, 2) - 
                                                          X2(0, 1) + X2(0, 2) + x2g[1] - x2g[2] + 
                                                          mult1 * (Deltax2g[0] + Deltax2g[1]) + 
                                                          mult2 * (Deltax2g[0] + Deltax2g[2])) )/std::pow(div1, 2);
                                                          
                       rDerivativeData.DeltaN2[i_dof][2] =  ( -mult3 * (coefmaster2 + (u2(0, 2) + X2(0, 2)) * multmaster5 + 
                                                            coefmaster1 - coefmaster5 - multmaster4 * multmaster6) + (div2) * (mult3 * Deltax2g[0] + 
                                                            mult4 * Deltax2g[1] + mult5 * Deltax2g[2]))/std::pow(div2, 2);
                   }
                   else if ((i_dof - TDim * TNumNodes) == 7)
                   {
                       rDerivativeData.DeltaN2[i_dof][1] =  (mult4 *(multmaster0*multmaster2 + 
                                                          multmaster1*multmaster3) - (div1) * (coefmaster3 + 
                                                          mult1 * (Deltax2g[0] + Deltax2g[1]) + 
                                                          mult2 * (Deltax2g[0] + Deltax2g[2])))/std::pow(div1, 2);
                                                          
                       rDerivativeData.DeltaN2[i_dof][2] =  ( -mult4 * (coefmaster2 + (u2(0, 2) + X2(0, 2)) * multmaster5 + 
                                                            coefmaster1 - coefmaster5 - multmaster4 * multmaster6) + (div2) * (mult3 * Deltax2g[0] + 
                                                            mult4 * Deltax2g[1] + mult5 * Deltax2g[2]))/std::pow(div2, 2);
                   }
                   else if ((i_dof - TDim * TNumNodes) == 8)
                   {
                       rDerivativeData.DeltaN2[i_dof][1] =  (mult5 *(multmaster0*multmaster2 + 
                                                          multmaster1*multmaster3) - (div1) * (coefmaster4 + 
                                                          mult1 * (Deltax2g[0] + Deltax2g[1]) + 
                                                          mult2 * (Deltax2g[0] + Deltax2g[2])))/std::pow(div1, 2);
                                                          
                       rDerivativeData.DeltaN2[i_dof][2] =  (-mult5 * (coefmaster2 + (u2(0, 2) + X2(0, 2)) * multmaster5 + 
                                                          coefmaster1 - coefmaster5 - multmaster4 * multmaster6) + (div2) * (mult3 * Deltax2g[0] + 
                                                          mult4 * Deltax2g[1] + mult5 * Deltax2g[2]) )/std::pow(div2, 2);
                   }

                  rDerivativeData.DeltaN2[i_dof][0] =  - rDerivativeData.DeltaN2[i_dof][1] - rDerivativeData.DeltaN2[i_dof][2];
               }
            }
         }
      }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateDeltaPhi(
   GeneralVariables& rVariables,
   DerivativeData& rDerivativeData
   )
{
    // Shape functions
    const VectorType N1 = rVariables.N_Slave;
    
    for (unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++)
    {
        for (unsigned int i_dim = 0; i_dim < TDim; i_dim++)
        {
            const unsigned int i_dof = i_slave * TDim + i_dim;
            
            rDerivativeData.DeltaPhi[i_dof] = prod(rDerivativeData.DeltaAe[i_dof], N1);;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& CurrentProcessInfo 
    )
{
    KRATOS_TRY;   
    
    const std::vector<contact_container> all_conditions = *( this->GetValue( CONTACT_CONTAINERS ) );
    
    // Calculates the size of the system
    const unsigned int condition_size = (TDim * ( TNumNodes + TNumNodes) + TNumNodes)* all_conditions.size(); 
    
    if (rResult.size() != condition_size)
    {
        rResult.resize( condition_size, false );
    }
    
    unsigned int index = 0;
    
    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    for ( unsigned int i_cond = 0;  i_cond < all_conditions.size(); ++i_cond )
    {   
        GeometryType& current_master = all_conditions[i_cond].condition->GetGeometry( );
        
        // Master Nodes Displacement Equation IDs
        for ( unsigned int i_master = 0; i_master < TNumNodes; i_master++ ) // NOTE: Assuming same number of nodes for master and slave
        {
            NodeType& master_node = current_master[i_master];
            rResult[index++] = master_node.GetDof( DISPLACEMENT_X ).EquationId( );
            rResult[index++] = master_node.GetDof( DISPLACEMENT_Y ).EquationId( );
            if (TDim == 3)
            {
                rResult[index++] = master_node.GetDof( DISPLACEMENT_Z ).EquationId( );
            }
        }

        // Slave Nodes Displacement Equation IDs
        for ( unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++ )
        {
            NodeType& slave_node = GetGeometry()[ i_slave ];
            rResult[index++] = slave_node.GetDof( DISPLACEMENT_X ).EquationId( );
            rResult[index++] = slave_node.GetDof( DISPLACEMENT_Y ).EquationId( );
            if (TDim == 3)
            {
                rResult[index++] = slave_node.GetDof( DISPLACEMENT_Z ).EquationId( );
            }
        }

        // Slave Nodes  Lambda Equation IDs
        for ( unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++ )
        {
            NodeType& slave_node = GetGeometry()[ i_slave ];
            rResult[index++] = slave_node.GetDof( NORMAL_CONTACT_STRESS ).EquationId( );
        }
        
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim, TNumNodes>::GetDofList(
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo 
)
{
    KRATOS_TRY;
    
    const std::vector<contact_container> all_conditions = *( this->GetValue( CONTACT_CONTAINERS ) );
    
    // Calculates the size of the system
    const unsigned int condition_size = (TDim * (TNumNodes + TNumNodes) + TNumNodes) * all_conditions.size(); 
    
    if (rConditionalDofList.size() != condition_size)
    {
        rConditionalDofList.resize( condition_size );
    }
    
    unsigned int index = 0;
    
    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    for ( unsigned int i_cond = 0; i_cond < all_conditions.size(); ++i_cond )
    {
        GeometryType& current_master = all_conditions[i_cond].condition->GetGeometry( );   

        // Master Nodes Displacement Equation IDs
        for ( unsigned int i_master = 0; i_master < TNumNodes; i_master++ ) // NOTE: Assuming same number of nodes for master and slave
        {
            NodeType& master_node = current_master[i_master];
            rConditionalDofList[index++] =master_node.pGetDof( DISPLACEMENT_X );
            rConditionalDofList[index++] =master_node.pGetDof( DISPLACEMENT_Y );
            if (TDim == 3)
            {
                rConditionalDofList[index++] =master_node.pGetDof( DISPLACEMENT_Z );
            }
        }

        // Slave Nodes Displacement Equation IDs
        for ( unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++ )
        {
            NodeType& slave_node = GetGeometry()[ i_slave ];
            rConditionalDofList[index++] =slave_node.pGetDof( DISPLACEMENT_X );
            rConditionalDofList[index++] =slave_node.pGetDof( DISPLACEMENT_Y );
            if (TDim == 3)
            {
                rConditionalDofList[index++] =slave_node.pGetDof( DISPLACEMENT_Z );
            }
        }

        // Slave Nodes Lambda Equation IDs
        for ( unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++ )
        {
            NodeType& slave_node = GetGeometry()[ i_slave ];
            rConditionalDofList[index++] =slave_node.pGetDof( NORMAL_CONTACT_STRESS );
        }
    }
    
    KRATOS_CATCH( "" );
}


//******************************* GET DOUBLE VALUE *********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::GetValueOnIntegrationPoints( 
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//******************************* GET ARRAY_1D VALUE *******************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::GetValueOnIntegrationPoints( 
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//******************************* GET VECTOR VALUE *********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::GetValueOnIntegrationPoints( 
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateOnIntegrationPoints( 
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;
    
    // TODO: Add the FRICTION_COEFFICIENT, and maybe if it is ACTIVE or SLIPPING the GP

    // Create and initialize condition variables:
    GeneralVariables rVariables;
    
    // Initialize the current contact data
    DerivativeData rDerivativeData;
    
    // Reading integration points
    const double integration_order = GetProperties().GetValue(INTEGRATION_ORDER_CONTACT);
    mColocationIntegration.Initialize( integration_order);
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

    // TODO: Add eventually
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateOnIntegrationPoints( 
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector< array_1d<double, 3 > >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
    
    // Create and initialize condition variables:
    GeneralVariables rVariables;
    
    // Initialize the current contact data
    DerivativeData rDerivativeData;
    
    // Reading integration points
    const double integration_order = GetProperties().GetValue(INTEGRATION_ORDER_CONTACT);
    mColocationIntegration.Initialize( integration_order);
    const GeometryType::IntegrationPointsArrayType& integration_points = mUseManualColocationIntegration ?
                                                                         mColocationIntegration.IntegrationPoints( ) :
                                                                         GetGeometry( ).IntegrationPoints( mThisIntegrationMethod );
                                                                                                                        
    const unsigned int number_of_integration_pts = integration_points.size();
    if ( rOutput.size() != number_of_integration_pts )
    {
        rOutput.resize( number_of_integration_pts );
    }
    
    const array_1d<double, 3> zero_vector = ZeroVector(3);
    for (unsigned int PointNumber = 0; PointNumber < number_of_integration_pts; PointNumber++)
    {
        rOutput[PointNumber] = zero_vector;
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateOnIntegrationPoints( 
    const Variable<Vector>& rVariable, 
    std::vector<Vector>& rOutput, 
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
    
    // TODO: Fill this!!!
    
    KRATOS_CATCH( "" );
}

/******************* AUXILLIARY METHODS FOR GENERAL CALCULATIONS *******************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::ComputeSelectiveIntegrationMethod(const unsigned int rPairIndex)
{
    const double integration_order = GetProperties().GetValue(INTEGRATION_ORDER_CONTACT);
    mUseManualColocationIntegration = true;
    
    if (TDim == 2)
    {
        if (TNumNodes == 2) // NOTE: Total weight of a  line is 2.0
        {
            // Using standart integration methods (I am using collocation)
            mColocationIntegration.Initialize( integration_order);
            
//             // Using exact integration
//             const double tol = 1.0e-4; 
//             const IntegrationMethod AuxIntegrationMethod = GetIntegrationMethod(integration_order, false);
//             GeometryType::IntegrationPointsArrayType IntegrationPointsConsidered;
//             
//             double total_weight = 0.0;
//             array_1d<double,2> coor_aux = ZeroVector(2);
//             
//             // Declaring auxiliar values
//             PointType projected_gp_global;
//             GeometryType::CoordinatesArrayType projected_gp_local;
//             const array_1d<double, 3> normal = this->GetValue(NORMAL);
//             double aux_dist = 0.0;
//             
//             // The master geometry
//             GeometryType& master_seg = mThisMasterElements[rPairIndex]->GetGeometry();
//             
//             // Declare the boolean of full integral
//             bool full_int = true;
//             
//             // First look if the edges of the slave are inside of the master, if not check if the opposite is true, if not then the element is not in contact
//             for (unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++)
//             {
//                 ContactUtilities::ProjectDirection(master_seg, GetGeometry()[i_slave].Coordinates(), projected_gp_global, aux_dist, -normal ); // The opposite direction
//                 
//                 const bool inside = master_seg.IsInside( projected_gp_global.Coordinates( ), projected_gp_local );
//                 
//                 if (inside == false)
//                 {
//                     full_int = false;
//                 }
//                 else
//                 {
//                     if (i_slave == 0)
//                     {
//                         coor_aux[0] = - 1.0;
//                     }
//                     else if (i_slave == 1)
//                     {
//                         coor_aux[1] =   1.0;
//                     }
//                 }
//             }
//             
//             if (full_int == true)
//             {
//                 total_weight = 2.0;
//             }
//             else
//             {
//                 std::vector<double> aux_xi;
//                 for (unsigned int i_master = 0; i_master < TNumNodes; i_master++)
//                 {
//                     ContactUtilities::ProjectDirection(GetGeometry(), master_seg[i_master].Coordinates(), projected_gp_global, aux_dist, normal );
// 
//                     const bool inside = GetGeometry().IsInside( projected_gp_global.Coordinates( ), projected_gp_local );
//                     
//                     if (inside == true)
//                     {
//                         aux_xi.push_back(projected_gp_local[0]);
//                     }
//                 }
//                 
//                 if (aux_xi.size() == 1)
//                 {
//                     if (coor_aux[0] == - 1.0)
//                     {
//                         coor_aux[1] = aux_xi[0];
//                     }
//                     else if (coor_aux[1] == 1.0)
//                     {
//                         coor_aux[0] = aux_xi[0];
//                     }
//                     else
//                     {
//                         KRATOS_WATCH("WARNING: THIS IS NOT SUPPOSED TO HAPPEN!!!!");
//                     }
//                 }
//                 else  if (aux_xi.size() == 2)
//                 {
//                     if (aux_xi[0] < aux_xi[1])
//                     {
//                         coor_aux[0] = aux_xi[0];
//                         coor_aux[1] = aux_xi[1];
//                     }
//                     else
//                     {
//                         coor_aux[1] = aux_xi[0];
//                         coor_aux[0] = aux_xi[1];
//                     }
//                 }
//                 
//                 total_weight = coor_aux[1] - coor_aux[0];
//             }
//             
//             if(total_weight < 0.0)
//             {
//                 KRATOS_THROW_ERROR( std::logic_error, "WAAAAAAAAAAAAARNING!!!!!!!!, wrong order of the coordinates", coor_aux);
//             }
//             
//             if (total_weight > tol)
// //             if (total_weight > 0.0)
//             {
//                 // With the proportion of the weigth you recalculate the integration weight. Change the coordinates of the integration to accomodate
//                 const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(AuxIntegrationMethod);
//                 IntegrationPointsConsidered.resize(integration_points.size());
//                 for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
//                 {
//                     const double weight = integration_points[PointNumber].Weight() * total_weight/2.0;
//                     const double xi = 0.5 * (1.0 - integration_points[PointNumber].Coordinate(1)) * coor_aux[0] 
//                                     + 0.5 * (1.0 + integration_points[PointNumber].Coordinate(1)) * coor_aux[1];
//                     
//                     IntegrationPointsConsidered[PointNumber] = IntegrationPoint<2>( xi, weight );
//                 }
//             }
//             else
//             {
// //                 IntegrationPointsConsidered.resize(0); // An empty std::vector
//                 IntegrationPointsConsidered.clear(); // An empty std::vector
//             }
//             
//             mColocationIntegration.SetIntegrationPoints(IntegrationPointsConsidered);
// //             if (IntegrationPointsConsidered.size() > 0)
// //             {
// //                 std::cout <<  GetGeometry()[0].X() << " " << GetGeometry()[0].Y() << " " << GetGeometry()[1].X() << " " << GetGeometry()[1].Y() << std::endl;
// //                 std::cout <<  master_seg[0].X() << " " << master_seg[0].Y() << " " << master_seg[1].X() << " " << master_seg[1].Y() << std::endl;
// //                 KRATOS_WATCH(coor_aux);
// //                 mColocationIntegration.print();
// //             }
        }
        else
        {
            // Using standart integration methods (I am using collocation)
            mColocationIntegration.Initialize( integration_order);
        }
    }
    else
    {
        if (TNumNodes == 3) // NOTE: Total weight of a triangle is 0.5
        {
            mColocationIntegration.Initialize( integration_order);
//             // TODO: Finish this
//             // Compute the local Coordinates of the master condition
//             PointType projected_gp_global;
//             const array_1d<double,3> normal = this->GetValue(NORMAL);
//             
//             GeometryType::CoordinatesArrayType slave_gp_global;
//             double aux_dist = 0.0;
//             
//             for (unsigned int i = 0; i < 3; i++)
//             {
//                 this->GetGeometry( ).GlobalCoordinates( slave_gp_global, local_point );
//                 ContactUtilities::ProjectDirection( master_seg, slave_gp_global, projected_gp_global, aux_dist, -normal ); // The opposite direction
//                 
//                 GeometryType::CoordinatesArrayType projected_gp_local;
//                 
//                 const bool inside = master_seg.IsInside( projected_gp_global.Coordinates( ), projected_gp_local ) ;
//             }
        }
        else
        {
            // Using standart integration methods (I am consideing collocation)
            mColocationIntegration.Initialize( integration_order);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< >
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<2, 2>::InitializeIntegrationMethod()
{
    mUseManualColocationIntegration = false;
    if( GetProperties().Has(INTEGRATION_ORDER_CONTACT) )
    {
        const double integration_order = GetProperties().GetValue(INTEGRATION_ORDER_CONTACT);

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
                mColocationIntegration.Initialize( integration_order);
            }
    }
    else
    {
        mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_5;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< >
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<2, 3>::InitializeIntegrationMethod()
{
    mUseManualColocationIntegration = false;
    if( GetProperties().Has(INTEGRATION_ORDER_CONTACT) )
    {
        const double integration_order = GetProperties().GetValue(INTEGRATION_ORDER_CONTACT);

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
                mColocationIntegration.Initialize( integration_order);
            }
    }
    else
    {
        mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_5;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< >
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 3>::InitializeIntegrationMethod()
{
    mUseManualColocationIntegration = false;
    if( GetProperties().Has(INTEGRATION_ORDER_CONTACT) )
    {
        const double integration_order = GetProperties().GetValue(INTEGRATION_ORDER_CONTACT);

//         if (integration_order == 3)
//         {
//         }
//         else
//         {
            mUseManualColocationIntegration = true;
            mColocationIntegration.Initialize( integration_order);
//         }
    }
    else
    {
        mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_5;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< >
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 4>::InitializeIntegrationMethod()
{
    mUseManualColocationIntegration = false;
    if( GetProperties().Has(INTEGRATION_ORDER_CONTACT) )
    {
        const double integration_order = GetProperties().GetValue(INTEGRATION_ORDER_CONTACT);

//         if (integration_order == )
//         {
//         }
//         else
//         {
            mUseManualColocationIntegration = true;
            mColocationIntegration.Initialize( integration_order);
//         }
    }
    else
    {
        mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_5;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class AugmentedLagrangianMethodFrictionlessMortarContactCondition<2, 2>;
template class AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 3>;
template class AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 4>;

} // Namespace Kratos
