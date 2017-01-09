// KRATOS  ___|  |       |       |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//           | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License: BSD License
//   license: structural_mechanics_application/license.txt
//
//  Main authors:  Vicente Mataix Ferrándiz
//                 Mohamed Khalil
//

// System includes

// External includes

// Project includes
/* Mortar includes */
#include "custom_conditions/mortar_contact_condition.h"
// #include "structural_mechanics_application.h"
// #include "structural_mechanics_application_variables.h"

/* Additional includes */
#include <algorithm>

/* Utilities */
#include "custom_utilities/contact_utilities.h"
#include "utilities/math_utils.h"
#include "custom_utilities/structural_mechanics_math_utilities.hpp"
// #include "../FSIapplication/custom_utilities/qr_utility.h" // QR decomposition utility used in matrix inversion.

/* Includes of particular contact conditions */
#include "contact_2D_2N_2N.hpp"
#include "contact_3D_3N_3N.hpp"
#include "contact_3D_4N_4N.hpp"
#include "contact_2D_2N_2N_double_lm.hpp"
// #include "contact_3D_3N_3N_double_lm.hpp"
// #include "contact_3D_4N_4N_double_lm.hpp"

/* Logging format include */
#include "custom_utilities/logging_settings.hpp"

namespace Kratos 
{
/**
 * Flags related to the condition computation 
 */
// Avoiding using the macro since this has a template parameter. If there was no template plase use the KRATOS_CREATE_LOCAL_FLAG macro
template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
const Kratos::Flags MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_RHS_VECTOR(Kratos::Flags::Create(0));
template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
const Kratos::Flags MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_LHS_MATRIX(Kratos::Flags::Create(1));
template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
const Kratos::Flags MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS(Kratos::Flags::Create(2));
template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
const Kratos::Flags MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS(Kratos::Flags::Create(3));

/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
Condition::Pointer MortarContactCondition<TDim,TNumNodes,TDoubleLM>::Create( 
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties ) const
{
    return boost::make_shared< MortarContactCondition<TDim,TNumNodes,TDoubleLM> >( NewId, this->GetGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
Condition::Pointer MortarContactCondition<TDim,TNumNodes,TDoubleLM>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return boost::make_shared< MortarContactCondition<TDim,TNumNodes,TDoubleLM> >( NewId, pGeom, pProperties );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
MortarContactCondition<TDim,TNumNodes,TDoubleLM>::~MortarContactCondition( )
{
}

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::Initialize( ) 
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

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    // First populate of the vector of master elements
    const std::vector<contact_container> * all_containers = this->GetValue( CONTACT_CONTAINERS );
    mPairSize = all_containers->size();
//     mThisMasterElements.clear(); // NOTE: Be careful, this calls the destructor
    mThisMasterElements.resize( mPairSize );
    
    for ( unsigned int i_cond = 0; i_cond < mPairSize; ++i_cond )
    {
        mThisMasterElements[i_cond] = (*all_containers)[i_cond].condition;
    }
   
//     if (mPairSize == 0) // NOTE: This is not supposed to be necessary
//     {
//         this->Set(ACTIVE, false);
//     }
//     else
//     {
//         this->Set(ACTIVE, true);
//     }

    // Set the initial value of the LM for the quasi-static problems
    if (rCurrentProcessInfo[TIME_STEPS] == 1) 
    {
        const double tol = 1.0e-16;
        for (unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++)
        {
            if (GetGeometry()[i_slave].Is(ACTIVE) == true)
            {
                const array_1d<double, 3> normal      = GetGeometry()[i_slave].GetValue(NORMAL);
                const array_1d<double, 3> tangent_xi  = GetGeometry()[i_slave].GetValue(TANGENT_XI);
                const array_1d<double, 3> tangent_eta = GetGeometry()[i_slave].GetValue(TANGENT_ETA);
//                 GetGeometry()[i_slave].FastGetSolutionStepValue( VECTOR_LAGRANGE_MULTIPLIER, 0) = tol * (- normal);
                GetGeometry()[i_slave].FastGetSolutionStepValue( VECTOR_LAGRANGE_MULTIPLIER, 0) = tol * (- normal + tangent_xi + tangent_eta);
            }
        }
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    // NOTE: Add things if necessary
        
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // NOTE: Add things if necessary
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // Create and initialize condition variables:
    GeneralVariables rVariables;
    
    // Initialize the current contact data
    ContactData<TDim, TNumNodes> rContactData;
    
    // Reading integration points
//     const GeometryType::IntegrationPointsArrayType& integration_points = mUseManualColocationIntegration ?
//                                                                          mColocationIntegration.IntegrationPoints( ) :
//                                                                          GetGeometry( ).IntegrationPoints( mThisIntegrationMethod );
    this->InitializeContactData(rContactData, rCurrentProcessInfo);
    
    this->CalculateAeAndDeltaAe(rContactData, rVariables, rCurrentProcessInfo);
    
    std::vector<contact_container> *& all_containers = this->GetValue(CONTACT_CONTAINERS);
    
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
        this->UpdateContactData(rContactData, PairIndex);
         
        // Master segment info
        const GeometryType& CurrentMasterElement = rVariables.GetMasterElement( );
        
        array_1d<double, TNumNodes> aux_int_gap      = ZeroVector(TNumNodes);
        array_1d<double, TNumNodes> aux_int_slip     = ZeroVector(TNumNodes);
        array_1d<double, TNumNodes> aux_int_friction = ZeroVector(TNumNodes);
        
        double total_weight = 0.0;
        
        // Integrating the LHS and RHS
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            // Calculate the kinematic variables
            const bool inside = this->CalculateKinematics( rVariables, rContactData, PointNumber, PairIndex, integration_points );
            
            if (inside == true)
            {
                // Calculate the gap in the integration node and check tolerance
                const double IntegrationPointGap = inner_prod(rContactData.Gaps, rVariables.N_Slave);
            
                const double IntegrationWeight = integration_points[PointNumber].Weight();
                total_weight += IntegrationWeight;

                // The slip of th GP
                double IntegrationPointSlip;
                this->CalculateAugmentedTangentLM(rVariables, rContactData, CurrentMasterElement, IntegrationPointSlip);
                
                for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
                {
                    const double int_aux = IntegrationWeight * rVariables.DetJSlave * rVariables.Phi_LagrangeMultipliers[iNode];
                    aux_int_gap[iNode]      +=  IntegrationPointGap  * int_aux;
                    aux_int_slip[iNode]     +=  IntegrationPointSlip * int_aux;
                    aux_int_friction[iNode] +=  rVariables.mu          * int_aux;
                }
            }
        }
        
        // We can consider the pair if at least one of the collocation point is inside 
        if (total_weight > 0.0)
        {
            (*all_containers)[PairIndex].active_pair = true;

            for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
            {
                #pragma omp atomic 
                GetGeometry()[iNode].GetValue(WEIGHTED_GAP)      += aux_int_gap[iNode]; 
                #pragma omp atomic 
                GetGeometry()[iNode].GetValue(WEIGHTED_SLIP)     += aux_int_slip[iNode]; 
                #pragma omp atomic 
                GetGeometry()[iNode].GetValue(WEIGHTED_FRICTION) += aux_int_friction[iNode]; 
            }
        }
        else
        {
            (*all_containers)[PairIndex].active_pair = false;
        }
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
IntegrationMethod MortarContactCondition<TDim,TNumNodes,TDoubleLM>::GetIntegrationMethod()
{   
    return mThisIntegrationMethod;
}

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
IntegrationMethod MortarContactCondition<TDim,TNumNodes,TDoubleLM>::GetIntegrationMethod(
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

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateLocalSystem( 
    std::vector<MatrixType>& rLeftHandSideMatrices,
    const std::vector<Variable<MatrixType> >& rLHSVariables,
    std::vector<VectorType>& rRightHandSideVectors,
    const std::vector<Variable<VectorType> >& rRHSVariables,
    ProcessInfo& rCurrentProcessInfo 
    )
{    
    // Calculates the size of the system
    constexpr unsigned int TMatrixSize = TDim * ((2 + TDoubleLM) * TNumNodes + TNumNodes);
    
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set(MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS, true);
    LocalSystem.CalculationFlags.Set(MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS, true);

    //Initialize sizes for the system components:
    if ( rLHSVariables.size( ) != rLeftHandSideMatrices.size( ) )
    {
        rLeftHandSideMatrices.resize( rLHSVariables.size( ) );
    }

    if ( rRHSVariables.size( ) != rRightHandSideVectors.size( ) )
    {
        rRightHandSideVectors.resize( rRHSVariables.size( ) );
    }

    LocalSystem.CalculationFlags.Set(MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_LHS_MATRIX, true);
    for ( unsigned int i = 0; i < rLeftHandSideMatrices.size( ); i++ )
    {
        // Note: rRightHandSideVectors.size() > 0
        this->InitializeSystemMatrices<TMatrixSize>( rLeftHandSideMatrices[i], rRightHandSideVectors[0],LocalSystem.CalculationFlags );
    }

    LocalSystem.CalculationFlags.Set( MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_RHS_VECTOR, true );
    LocalSystem.CalculationFlags.Set( MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_LHS_MATRIX, false ); // Temporarily only
    for ( unsigned int i = 0; i < rRightHandSideVectors.size( ); i++ )
    {
        // Note: rLeftHandSideMatrices.size() > 0
        this->InitializeSystemMatrices<TMatrixSize>( rLeftHandSideMatrices[0], rRightHandSideVectors[i], LocalSystem.CalculationFlags  );
    }
    LocalSystem.CalculationFlags.Set( MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_LHS_MATRIX, true ); // Reactivated again

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

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;
    
    // Calculates the size of the system
    constexpr unsigned int TMatrixSize = TDim * ((2 + TDoubleLM) * TNumNodes + TNumNodes);
    
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_LHS_MATRIX, true );
    LocalSystem.CalculationFlags.Set( MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_RHS_VECTOR, true );

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

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateLeftHandSide( 
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Calculates the size of the system
    constexpr unsigned int TMatrixSize = TDim * ((2 + TDoubleLM) * TNumNodes + TNumNodes);
    
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_LHS_MATRIX, true );

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

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateLeftHandSide( 
    std::vector< MatrixType >& rLeftHandSideMatrices,
    const std::vector< Variable< MatrixType > >& rLHSVariables,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Calculates the size of the system
    constexpr unsigned int TMatrixSize = TDim * ((2 + TDoubleLM) * TNumNodes + TNumNodes);
    
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_LHS_MATRIX, true );
    LocalSystem.CalculationFlags.Set( MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS, true );

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

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateRightHandSide( 
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Calculates the size of the system
    constexpr unsigned int TMatrixSize = TDim * ((2 + TDoubleLM) * TNumNodes + TNumNodes);
    
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_RHS_VECTOR, true);

    MatrixType LeftHandSideMatrix = Matrix( );

    // Initialize size for the system components
    this->InitializeSystemMatrices<TMatrixSize>( LeftHandSideMatrix, rRightHandSideVector,LocalSystem.CalculationFlags);

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix( LeftHandSideMatrix );
    LocalSystem.SetRightHandSideVector( rRightHandSideVector );

    //Calculate condition system
    this->CalculateConditionSystem<TMatrixSize>( LocalSystem, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateRightHandSide( 
    std::vector< VectorType >& rRightHandSideVectors,
    const std::vector< Variable< VectorType > >& rRHSVariables,
    ProcessInfo& rCurrentProcessInfo )
{
    // Calculates the size of the system
    constexpr unsigned int TMatrixSize = TDim * ((2 + TDoubleLM) * TNumNodes + TNumNodes);
        
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_RHS_VECTOR, true );
    LocalSystem.CalculationFlags.Set( MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS, true );

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

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
template< unsigned int TMatrixSize >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::InitializeSystemMatrices( 
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    Flags& rCalculationFlags 
    )
{
    const unsigned int condition_size = this->CalculateConditionSize<TMatrixSize>( );
    
    // Resizing as needed the LHS
    if ( rCalculationFlags.Is( MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_LHS_MATRIX ) ) // Calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != condition_size )
        {
            rLeftHandSideMatrix.resize( condition_size, condition_size, false );
        }
        noalias( rLeftHandSideMatrix ) = ZeroMatrix( condition_size, condition_size ); // Resetting LHS
    }

    // Resizing as needed the RHS
    if ( rCalculationFlags.Is( MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_RHS_VECTOR ) ) // Calculation of the matrix is required
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

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateMassMatrix( 
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

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateDampingMatrix( 
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

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
template< unsigned int TMatrixSize >
const unsigned int MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateConditionSize( )
{
    const unsigned int condition_size = mPairSize * TMatrixSize; // NOTE: Assuming same number of nodes for the master
    
    return condition_size;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
template< unsigned int TMatrixSize>
void MortarContactCondition<TDim, TNumNodes, TDoubleLM>::CalculateConditionSystem( 
    LocalSystemComponents& rLocalSystem,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
    
    // Create and initialize condition variables:#pragma omp critical
    GeneralVariables rVariables;
    
    // Initialize the current contact data
    ContactData<TDim, TNumNodes> rContactData;

//     // Reading integration points 
//     const GeometryType::IntegrationPointsArrayType& integration_points = mUseManualColocationIntegration ?
//                                                                          mColocationIntegration.IntegrationPoints( ) :
//                                                                          GetGeometry( ).IntegrationPoints( mThisIntegrationMethod );
                                                                                  
    this->InitializeContactData(rContactData, rCurrentProcessInfo);
    
    // Compute Ae and its derivative
    this->CalculateAeAndDeltaAe(rContactData, rVariables, rCurrentProcessInfo); 
    
    // Compute the normal and tangent derivatives of the slave
    this->CalculateDeltaNormalTangentSlave(rContactData);
    
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
        this->UpdateContactData(rContactData, PairIndex);
         
        // Compute the normal derivatives of the master
        this->CalculateDeltaNormalMaster(rContactData);
        
        // Master segment info
        const GeometryType& CurrentMasterElement = rVariables.GetMasterElement( );
        
        // Initialize the integration weight
        double total_weight = 0.0;
        
        // Integrating the LHS and RHS
        bounded_matrix<double, TMatrixSize, TMatrixSize> LHS_contact_pair = ZeroMatrix(TMatrixSize, TMatrixSize);
        array_1d<double, TMatrixSize> RHS_contact_pair = ZeroVector(TMatrixSize);
        
        // We check the case to compute
        const unsigned int compute = this->CaseToCompute(integration_points);
        
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            // Calculate the kinematic variables
            bool inside = this->CalculateKinematics( rVariables, rContactData, PointNumber, PairIndex, integration_points );
            
            if (inside == true)
            {   
                /* Update the derivatives */
                // Update the derivative of DetJ
                this->CalculateDeltaDetJSlave(rVariables, rContactData);
                // Update the derivatives of the shape functions and the gap
                this->CalculateDeltaNAndDeltaGap(rVariables, rContactData, compute);
                // The derivatives of the dual shape function
                this->CalculateDeltaPhi(rVariables, rContactData);
                // Compute in each node the tangent functional and their derivatives
                double AugmentedNormalLM = this->CalculateCtanAndDeltaCtan(rVariables, rContactData);
                
                if (AugmentedNormalLM < 0.0)
                {
                    const double IntegrationWeight = integration_points[PointNumber].Weight();
                    total_weight += IntegrationWeight;
                
                    // Calculation of the matrix is required
                    if ( rLocalSystem.CalculationFlags.Is( MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_LHS_MATRIX ) ||
                            rLocalSystem.CalculationFlags.Is( MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
                    {
                        this->CalculateLocalLHS<TMatrixSize>( LHS_contact_pair, rVariables, rContactData, IntegrationWeight);
                    }
                    
                    // Calculation of the vector is required
                    if ( rLocalSystem.CalculationFlags.Is( MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_RHS_VECTOR ) ||
                            rLocalSystem.CalculationFlags.Is( MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
                    {
                        this->CalculateLocalRHS<TMatrixSize>( RHS_contact_pair, rVariables, rContactData, IntegrationWeight);
                    }
                }
            }
        }
        
        // We can consider the pair if at least one of the collocation point is inside 
        if (total_weight > 0.0)
        {
            // Assemble of the matrix is required
            if ( rLocalSystem.CalculationFlags.Is( MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_LHS_MATRIX ) ||
                    rLocalSystem.CalculationFlags.Is( MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
            {
                // Contributions to stiffness matrix calculated on the reference config
                this->CalculateAndAddLHS<TMatrixSize>( rLocalSystem, LHS_contact_pair, PairIndex, CurrentMasterElement );
            }

            // Assemble of the vector is required
            if ( rLocalSystem.CalculationFlags.Is( MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_RHS_VECTOR ) ||
                    rLocalSystem.CalculationFlags.Is( MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
            {
                // Contribution to previous step contact force and residuals vector
                this->CalculateAndAddRHS<TMatrixSize>( rLocalSystem, RHS_contact_pair, PairIndex, CurrentMasterElement );
            }
            
//             std::cout << "--------------------------------------------------" << std::endl;
//             KRATOS_WATCH(this->Id());
//             KRATOS_WATCH(PairIndex);
// //             KRATOS_WATCH(rContactData.u1);
// //             KRATOS_WATCH(rContactData.LagrangeMultipliers);
// // //             KRATOS_WATCH(rContactData.DoubleLagrangeMultipliers);
// //             KRATOS_WATCH(LHS_contact_pair);
//             LOG_MATRIX_PRETTY( LHS_contact_pair );
// //             KRATOS_WATCH(RHS_contact_pair);
//             LOG_VECTOR_PRETTY( RHS_contact_pair );
        }
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::InitializeGeneralVariables(
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

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateAeAndDeltaAe(
    ContactData<TDim, TNumNodes>& rContactData,
    GeneralVariables& rVariables,
//     const GeometryType::IntegrationPointsArrayType& integration_points,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    double total_weight = 0.0; // NOTE: The integral is supposed to be in the domain partially integrated, I don't know if consider any additional thing for the partial integration
    
    rContactData.InitializeDeltaAeComponents();
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
        this->UpdateContactData(rContactData, PairIndex);
            
        // Calculating the proportion between the integrated area and segment area
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            // Calculate the kinematic variables
            bool inside = this->CalculateKinematics( rVariables, rContactData, PointNumber, PairIndex, integration_points );
            
            if (inside == true)
            {   
                const double IntegrationWeight = integration_points[PointNumber].Weight();
                total_weight += IntegrationWeight;
                this->CalculateDeltaAeComponents(rVariables, rContactData, IntegrationWeight);
            }
        }
    }
    
    // We can consider the pair if at least one of the collocation point is inside 
    if (total_weight > 0.0)
    {
        this->CalculateDeltaAe(rContactData);
        rContactData.ClearDeltaAeComponents();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::InitializeContactData(
    ContactData<TDim, TNumNodes>& rContactData,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // Slave element info
    rContactData.Initialize(GetGeometry());
    
    /* Set Delta time */
    rContactData.Dt = rCurrentProcessInfo[DELTA_TIME];
    
    /* LM */
    rContactData.LagrangeMultipliers = GetVariableMatrix(GetGeometry(), VECTOR_LAGRANGE_MULTIPLIER, 0); 
    if (TDoubleLM == true)
    {
        rContactData.DoubleLagrangeMultipliers = GetVariableMatrix(GetGeometry(), DOUBLE_LM, 0); 
    }
    
    /* NORMALS AND TANGENTS */ // TODO: To interpolate it is necessary to consider a smooth function
    const array_1d<double, 3> normal      = this->GetValue(NORMAL);
    const array_1d<double, 3> tangent_xi  = this->GetValue(TANGENT_XI);
    const array_1d<double, 3> tangent_eta = this->GetValue(TANGENT_ETA);
    for (unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++)
    {
        for (unsigned int i_dim = 0; i_dim < TDim; i_dim++)
        {
            rContactData.Normal_s(i_slave, i_dim)      = normal[i_dim];
            rContactData.Tangent_xi_s(i_slave, i_dim)  = tangent_xi[i_dim];
            rContactData.Tangent_eta_s(i_slave, i_dim) = tangent_eta[i_dim];
        }
    }
    
//     rContactData.Normal_s = GetVariableMatrix(GetGeometry(),  NORMAL); 
//     rContactData.Tangent_xi_s = GetVariableMatrix(GetGeometry(), TANGENT_XI); 
//     if (TDim == 3)
//     {
//         rContactData.Tangent_eta_s = GetVariableMatrix(GetGeometry(), TANGENT_ETA); 
//     }

    if (GetProperties().Has(NORMAL_AUGMENTATION_FACTOR) == true)
    {
        rContactData.epsilon_normal  = GetProperties().GetValue(NORMAL_AUGMENTATION_FACTOR);
    }
    if (GetProperties().Has(TANGENT_AUGMENTATION_FACTOR) == true)
    {
        rContactData.epsilon_tangent = GetProperties().GetValue(TANGENT_AUGMENTATION_FACTOR);
    }
    
    if (TDoubleLM == true)
    {
        if (GetProperties().Has(DOUBLE_LM_FACTOR) == true)
        {
            rContactData.epsilon = GetProperties().GetValue(DOUBLE_LM_FACTOR);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::UpdateContactData(
    ContactData<TDim, TNumNodes>& rContactData,
    const unsigned int& rMasterElementIndex
    )
{    
    // Master segment info
    GeometryType& CurrentMasterElement = mThisMasterElements[rMasterElementIndex]->GetGeometry();
    
    // Slave element info
    rContactData.UpdateMasterPair(mThisMasterElements[rMasterElementIndex] );
//     rContactData.UpdateMasterPair(CurrentMasterElement);
    
    /* NORMALS AND GAPS */
    for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
    {
        array_1d<double,3> normal = this->GetValue(NORMAL);
//         array_1d<double,3> normal = GetGeometry()[iNode].GetValue(NORMAL);
        
        PointType projected_global;
        ContactUtilities::ProjectDirection(CurrentMasterElement, GetGeometry()[iNode], projected_global, rContactData.Gaps(iNode), normal ); // NOTE: This is not the CPP, so the solution can be wrong
    }
}

/*********************************COMPUTE KINEMATICS*********************************/
/************************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
bool MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateKinematics( 
    GeneralVariables& rVariables,
    const ContactData<TDim, TNumNodes> rContactData,
    const double& rPointNumber,
    const unsigned int& rPairIndex,
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
    rVariables.Phi_LagrangeMultipliers = prod(rContactData.Ae, rVariables.N_Slave);
//     rVariables.Phi_LagrangeMultipliers = rVariables.N_Slave; // TODO: This could be needed in the future to be different than the standart shape functions 
    
    // SHAPE FUNCTION DERIVATIVES
    slave_nodes.ShapeFunctionsLocalGradients( rVariables.DN_De_Slave, local_point );
//     slave_nodes.ShapeFunctionsLocalGradients( rVariables.DN_De_Slave , local_point );// TODO: This could be needed in the future to be different than the standart shape functions
    
    // MASTER CONDITION
    const bool inside = this->MasterShapeFunctionValue( rVariables, local_point);
    
    /* CALCULATE JACOBIAN AND JACOBIAN DETERMINANT */
    slave_nodes.Jacobian( rVariables.j_Slave, local_point.Coordinates() );
    rVariables.DetJSlave = ContactUtilities::ContactElementDetJacobian( rVariables.j_Slave );
    
    /* FRICTION COEFFICIENT */
    if (GetProperties().Has(FRICTION_COEFFICIENT) == true) // NOTE: In function of the friction law (TODO: Remove this!!!)
    {
        rVariables.mu = GetProperties().GetValue(FRICTION_COEFFICIENT);
    }        
    
    return inside;
}
 
/***********************************************************************************/
/*************** METHODS TO CALCULATE THE CONTACT CONDITION MATRICES ***************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
bool MortarContactCondition<TDim,TNumNodes,TDoubleLM>::MasterShapeFunctionValue(
    GeneralVariables& rVariables,
    const PointType& local_point 
    )
{    
    GeometryType& master_seg = rVariables.GetMasterElement( );

    PointType projected_gp_global;
    const array_1d<double,3> normal = this->GetValue(NORMAL);
//     const array_1d<double,3> normal = ContactUtilities::GaussPointNormal(rVariables.N_Slave, GetGeometry());
    
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

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateDeltaAeComponents(
    GeneralVariables& rVariables,
    ContactData<TDim, TNumNodes>& rContactData,
    const double& rIntegrationWeight
    )
{
    /* DEFINITIONS */
    const array_1d<double, TNumNodes> N1 = rVariables.N_Slave;
    const double detJ = rVariables.DetJSlave; 
     
    rContactData.De += rIntegrationWeight * this->ComputeDe( N1, detJ);
    rContactData.Me += rIntegrationWeight * this->ComputeMe( N1, detJ);
    
    for (unsigned int i = 0; i < TDim * TNumNodes; i++)
    {
        const double DeltaDetJ = rContactData.DeltaJ_s[i];
        
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
        
        rContactData.DeltaDe[i] += rIntegrationWeight * DeltaDe;
        rContactData.DeltaMe[i] += rIntegrationWeight * DeltaMe;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
bounded_matrix<double, TNumNodes, TNumNodes> MortarContactCondition<TDim,TNumNodes,TDoubleLM>::ComputeDe(
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

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
bounded_matrix<double, TNumNodes, TNumNodes> MortarContactCondition<TDim,TNumNodes,TDoubleLM>::ComputeMe(
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

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateDeltaAe(ContactData<TDim, TNumNodes>& rContactData)
{        
    Matrix InvMe = ZeroMatrix(TNumNodes, TNumNodes);
    // NOTE: Legacy inversion. In case Me is almost singular or singular (few GP integrated), will be considered as ZeroMatrix 
    if (TNumNodes == 2)
    {
        StructuralMechanicsMathUtilities::InvMat2x2(rContactData.Me, InvMe);
    }
    else if (TNumNodes == 3)
    {
        StructuralMechanicsMathUtilities::InvMat3x3(rContactData.Me, InvMe);
    }   
    else
    {
        StructuralMechanicsMathUtilities::InvertMatrix(rContactData.Me, InvMe);
    }   
    
//     // Inversion using the QR decompisition // NOTE: Giving problems in the cases of almost singular matrix
//     QR<double, row_major> QR_decomposition;     
//     QR_decomposition.compute(TNumNodes, TNumNodes, &rContactData.Me(0, 0));
//     
//     Matrix aux_I = identity_matrix<double>( TNumNodes ); // NOTE: Simplify this code¿?
//     for (unsigned int col = 0; col < TNumNodes; col++)
//     {   Vector aux_I_col = column(aux_I, col);
//         Vector aux_InvMe_col = ZeroVector(TNumNodes);
//         QR_decomposition.solve(&aux_I_col[0], &aux_InvMe_col[0]);
//         column(InvMe, col) = aux_InvMe_col;
//     }
    
    noalias(rContactData.Ae) = prod(rContactData.De, InvMe);
    
    for (unsigned int i = 0; i < TDim * TNumNodes; i++)
    {
        rContactData.DeltaAe[i] = rContactData.DeltaDe[i] - prod(rContactData.Ae, rContactData.DeltaMe[i]);
        rContactData.DeltaAe[i] = prod(rContactData.DeltaAe[i], InvMe);
//         rContactData.DeltaAe[i] = ZeroMatrix(TNumNodes, TNumNodes); // NOTE: Test with zero derivative
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
template< unsigned int TMatrixSize >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateAndAddLHS(
    LocalSystemComponents& rLocalSystem,
    bounded_matrix<double, TMatrixSize, TMatrixSize>& LHS_contact_pair, 
    const unsigned int rPairIndex,
    const GeometryType& CurrentMasterElement
    )
{
    /* DEFINITIONS */
    const unsigned int number_of_total_nodes = TNumNodes + TNumNodes; // NOTE: Assuming same number of nodes for the master
        
    unsigned int index_1 = 0;
    for (unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave )
    {
        index_1 = (number_of_total_nodes + i_slave) * TDim;
        
        const bounded_matrix<double, TDim, TDim> RotationMatrix = GetRotationMatrixSlave(i_slave);
        for (unsigned int i_node = 0; i_node < number_of_total_nodes; i_node++)
//         for (unsigned int i_node = 0; i_node < (number_of_total_nodes + (1 + TDoubleLM) * TNumNodes); i_node++)
        {
            const unsigned index_rot = i_node * TDim;
            subrange(LHS_contact_pair, index_1, index_1 + TDim, index_rot, index_rot + TDim) = prod(RotationMatrix, subrange(LHS_contact_pair, index_1, index_1 + TDim, index_rot, index_rot + TDim));
        }
        if (TDoubleLM == true) // TODO: For a consistent Double LM calculate an independent rotationmatrix
        {
            const unsigned int index_2 = index_1 + TNumNodes * TDim;
            for (unsigned int i_node = 0; i_node < number_of_total_nodes; i_node++)
//             for (unsigned int i_node = 0; i_node < (number_of_total_nodes + 2 * TNumNodes); i_node++)
            {
                const unsigned index_rot = i_node * TDim;
                subrange(LHS_contact_pair, index_2, index_2 + TDim, index_rot, index_rot + TDim) = prod(RotationMatrix, subrange(LHS_contact_pair, index_2, index_2 + TDim, index_rot, index_rot + TDim));
            } 
        }
        
        // We impose a 0 zero LM in the inactive nodes
        if (GetGeometry( )[i_slave].Is(ACTIVE) == false)
        {
//             subrange(LHS_contact_pair, 0, number_of_total_nodes * TDim, index_1, index_1 + TDim) = ZeroMatrix(number_of_total_nodes * TDim, TDim); // TODO: Consider or not?
            subrange(LHS_contact_pair, index_1, index_1 + TDim, 0, number_of_total_nodes * TDim) = ZeroMatrix(TDim, number_of_total_nodes * TDim);
            subrange(LHS_contact_pair, index_1, index_1 + TDim, index_1, index_1 + TDim)         = IdentityMatrix(TDim, TDim);
            
            if (TDoubleLM == true)
            {
                const unsigned int index_2 = index_1 + TNumNodes * TDim;
//                 subrange(LHS_contact_pair, 0, number_of_total_nodes * TDim, index_1, index_1 + TDim) = ZeroMatrix(number_of_total_nodes * TDim, TDim); // TODO: Consider or not?
                subrange(LHS_contact_pair, index_2, index_2 + TDim, 0, number_of_total_nodes * TDim) = ZeroMatrix(TDim, number_of_total_nodes * TDim);
                subrange(LHS_contact_pair, index_2, index_2 + TDim, index_2, index_2 + TDim)         = IdentityMatrix(TDim, TDim); 
                subrange(LHS_contact_pair, index_1, index_1 + TDim, index_2, index_2 + TDim)         = ZeroMatrix(TDim, TDim); 
                subrange(LHS_contact_pair, index_2, index_2 + TDim, index_1, index_1 + TDim)         = ZeroMatrix(TDim, TDim); 
            }
        }
        // And the tangent direction (sliding)
        else
        {            
            if (TDoubleLM == false)
            {
                Matrix tangent_matrix;
                tangent_matrix.resize(TDim - 1, TDim);

                const array_1d<double, 3> tangent1 = GetGeometry()[i_slave].GetValue(TANGENT_XI);
                const array_1d<double, 3> tangent2 = GetGeometry()[i_slave].GetValue(TANGENT_ETA);
                for (unsigned int i = 0; i < TDim; i++)
                {
                    tangent_matrix(0, i) = tangent1[i];
                    if (TDim == 3)
                    {
                        tangent_matrix(1, i) = tangent2[i];
                    }
                }
                
                subrange(LHS_contact_pair, index_1 + 1, index_1 + TDim, index_1 , index_1 + TDim) = tangent_matrix;
            }
        }
    }
    
    if ( rLocalSystem.CalculationFlags.Is( MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
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

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
template< unsigned int TMatrixSize>
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::AssembleContactPairLHSToConditionSystem(
    bounded_matrix<double, TMatrixSize, TMatrixSize>& rPairLHS,
    MatrixType& rConditionLHS,
    const unsigned int rPairIndex
    )
{
    // Find location of the pair's master DOFs in ConditionRHS
    const unsigned int index_begin = rPairIndex * TMatrixSize;
    const unsigned int index_end  = index_begin + TMatrixSize;
    
    subrange( rConditionLHS, index_begin, index_end, index_begin, index_end) += rPairLHS;
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MortarContactCondition<2, 2, false>::CalculateLocalLHS<12>(
    bounded_matrix<double, 12, 12>& rPairLHS,
    GeneralVariables& rVariables,
    const ContactData<2, 2>& rContactData,
    const double& rIntegrationWeight
    )
{
    /* DEFINITIONS */
    const array_1d<double, 2> N1  = rVariables.N_Slave;
    const array_1d<double, 2> N2  = rVariables.N_Master;
    const array_1d<double, 2> Phi = rVariables.Phi_LagrangeMultipliers;
    const double detJ = rVariables.DetJSlave; 

    rPairLHS += rIntegrationWeight * Contact2D2N2N::ComputeGaussPointFrictionlessLHS( N1, N2, Phi, detJ, rContactData);
//     rPairLHS += rIntegrationWeight * Contact2D2N2N::ComputeGaussPointLHS( N1, N2, Phi, detJ, rVariables.mu, rContactData);
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MortarContactCondition<2, 3, false>::CalculateLocalLHS<18>(
    bounded_matrix<double, 18, 18>& rPairLHS,
    GeneralVariables& rVariables,
    const ContactData<2, 3>& rContactData,
    const double& rIntegrationWeight
    )
{
//     /* DEFINITIONS */
//     const array_1d<double, 3> N1  = rVariables.N_Slave;
//     const array_1d<double, 3> N2  = rVariables.N_Master;
//     const array_1d<double, 3> Phi = rVariables.Phi_LagrangeMultipliers;
//     const double detJ = rVariables.DetJSlave; 
// 
    // TODO: Finish this!!!!
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MortarContactCondition<3, 3, false>::CalculateLocalLHS<27>(
    bounded_matrix<double, 27, 27>& rPairLHS,
    GeneralVariables& rVariables,
    const ContactData<3, 3>& rContactData,
    const double& rIntegrationWeight
    )
{
    /* DEFINITIONS */
    const array_1d<double, 3> N1  = rVariables.N_Slave;
    const array_1d<double, 3> N2  = rVariables.N_Master;
    const array_1d<double, 3> Phi = rVariables.Phi_LagrangeMultipliers;
    const double detJ         = rVariables.DetJSlave; 

    rPairLHS += rIntegrationWeight * Contact3D3N3N::ComputeGaussPointFrictionlessLHS( N1, N2, Phi, detJ, rContactData);
//     rPairLHS += rIntegrationWeight * Contact3D3N3N::ComputeGaussPointLHS( N1, N2, Phi, detJ, rVariables.mu, rContactData);
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MortarContactCondition<3, 4, false>::CalculateLocalLHS<36>(
    bounded_matrix<double, 36, 36>& rPairLHS,
    GeneralVariables& rVariables,
    const ContactData<3, 4>& rContactData,
    const double& rIntegrationWeight
    )
{
    /* DEFINITIONS */
    const array_1d<double, 4> N1  = rVariables.N_Slave;
    const array_1d<double, 4> N2  = rVariables.N_Master;
    const array_1d<double, 4> Phi = rVariables.Phi_LagrangeMultipliers;
    const double detJ = rVariables.DetJSlave; 

    rPairLHS += rIntegrationWeight * Contact3D4N4N::ComputeGaussPointFrictionlessLHS( N1, N2, Phi, detJ, rContactData);
//     rPairLHS += rIntegrationWeight * Contact3D4N4N::ComputeGaussPointLHS( N1, N2, Phi, detJ, rVariables.mu, rContactData);
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MortarContactCondition<2, 2, true>::CalculateLocalLHS<16>(
    bounded_matrix<double, 16, 16>& rPairLHS,
    GeneralVariables& rVariables,
    const ContactData<2, 2>& rContactData,
    const double& rIntegrationWeight
    )
{
    /* DEFINITIONS */
    const array_1d<double, 2>  N1  = rVariables.N_Slave;
    const array_1d<double, 2>  N2  = rVariables.N_Master;
    const array_1d<double, 2>  Phi = rVariables.Phi_LagrangeMultipliers;
    const double detJ = rVariables.DetJSlave; 

    rPairLHS += rIntegrationWeight * Contact2D2N2NDLM::ComputeGaussPointFrictionlessLHS( N1, N2, Phi, detJ, rContactData);
//     rPairLHS += rIntegrationWeight * Contact2D2N2NDLM::ComputeGaussPointLHS( N1, N2, Phi, detJ, rVariables.mu, rContactData);
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MortarContactCondition<2, 3, true>::CalculateLocalLHS<24>(
    bounded_matrix<double, 24, 24>& rPairLHS,
    GeneralVariables& rVariables,
    const ContactData<2, 3>& rContactData,
    const double& rIntegrationWeight
    )
{
//     /* DEFINITIONS */
//     const array_1d<double, 3> N1  = rVariables.N_Slave;
//     const array_1d<double, 3> N2  = rVariables.N_Master;
//     const array_1d<double, 3> Phi = rVariables.Phi_LagrangeMultipliers;
//     const double detJ         = rVariables.DetJSlave; 
// 
    // TODO: Finish this!!!!
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MortarContactCondition<3, 3, true>::CalculateLocalLHS<36>(
    bounded_matrix<double, 36, 36>& rPairLHS,
    GeneralVariables& rVariables,
    const ContactData<3, 3>& rContactData,
    const double& rIntegrationWeight
    )
{
//     /* DEFINITIONS */
//     const array_1d<double, 3> N1  = rVariables.N_Slave;
//     const array_1d<double, 3> N2  = rVariables.N_Master;
//     const array_1d<double, 3> Phi = rVariables.Phi_LagrangeMultipliers;
//     const double detJ         = rVariables.DetJSlave; 
//                
//     rPairLHS += rIntegrationWeight * Contact3D3N3NDLM::ComputeGaussPointFrictionlessLHS( N1, N2, Phi, detJ, rContactData);
// //     rPairLHS += rIntegrationWeight * Contact3D3N3NDLM::ComputeGaussPointLHS( N1, N2, Phi, detJ, rVariables.mu, rContactData);
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MortarContactCondition<3, 4, true>::CalculateLocalLHS<48>(
    bounded_matrix<double, 48, 48>& rPairLHS,
    GeneralVariables& rVariables,
    const ContactData<3, 4>& rContactData,
    const double& rIntegrationWeight
    )
{
//     /* DEFINITIONS */
//     const array_1d<double, 4> N1  = rVariables.N_Slave;
//     const array_1d<double, 4> N2  = rVariables.N_Master;
//     const array_1d<double, 4> Phi = rVariables.Phi_LagrangeMultipliers;
//     const double detJ = rVariables.DetJSlave; 
// 
//     rPairLHS += rIntegrationWeight * Contact3D4N4NDLM::ComputeGaussPointFrictionlessLHS( N1, DN1, N2, DN2, Phi, detJ, rContactData);
// //     rPairLHS += rIntegrationWeight * Contact3D4N4NDLM::ComputeGaussPointLHS( N1, DN1, N2, DN2, Phi, detJ, rVariables.mu, rContactData);
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
template< unsigned int TMatrixSize >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateAndAddRHS(
    LocalSystemComponents& rLocalSystem,
    array_1d<double, TMatrixSize>& RHS_contact_pair,
    const unsigned int rPairIndex,
    const GeometryType& CurrentMasterElement
    )
{   
    /* DEFINITIONS */
    const unsigned int number_of_total_nodes = TNumNodes + TNumNodes; // NOTE: Assuming same number of nodes for the master
        
    unsigned int index = 0;
    
    for (unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave )
    {
        index = (number_of_total_nodes + i_slave) * TDim; 

        const bounded_matrix<double, TDim, TDim> RotationMatrix = GetRotationMatrixSlave(i_slave);
    
        array_1d<double, TDim> aux_rhs = subrange(RHS_contact_pair, index, index + TDim);
        array_1d<double, TDim> aux_product;
        axpy_prod(RotationMatrix, aux_rhs, aux_product, true);

        for (unsigned i_dim = 0; i_dim < TDim; i_dim++)
        {
            RHS_contact_pair[index + i_dim] = aux_product[i_dim];
        }

        if (TDoubleLM == true) // TODO: For a consistent Double LM calculate an independent rotationmatrix
        {
            const unsigned int index_dlm = index + TNumNodes * TDim;
            
            aux_rhs = subrange(RHS_contact_pair, index_dlm, index_dlm + TDim);
            axpy_prod(RotationMatrix, aux_rhs, aux_product, true);
            
            for (unsigned i_dim = 0; i_dim < TDim; i_dim++)
            {
                RHS_contact_pair[index_dlm + i_dim] = aux_product[i_dim];
            }
        }
        
        if (GetGeometry( )[i_slave].Is(ACTIVE) == false)
        {
            for (unsigned int i = index; i < index + TDim; i++)
            {
                RHS_contact_pair[i] = 0.0;
            }
            if (TDoubleLM == true)
            {
                index += TNumNodes * TDim;
                for (unsigned int i = index; i < index + TDim; i++)
                {
                    RHS_contact_pair[i] = 0.0;
                }
            }
        }
        else // NOTE: In the tangent direction we impose sliding
        {   
            for (unsigned int i = index + 1; i < index + TDim; i++)
            {
                RHS_contact_pair[i] = 0.0;
            }
            if (TDoubleLM == true)
            {
                index += TNumNodes * TDim;
                for (unsigned int i = index + 1; i < index + TDim; i++)
                {
                    RHS_contact_pair[i] = 0.0;
                }
            }
        }
    }
        
    if ( rLocalSystem.CalculationFlags.Is( MortarContactCondition<TDim,TNumNodes,TDoubleLM>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
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

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
template< unsigned int TMatrixSize >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::AssembleContactPairRHSToConditionSystem(
    array_1d<double, TMatrixSize>& rPairRHS,
    VectorType& rConditionRHS,
    const unsigned int rPairIndex
    )
{
    // Find location of the pair's master DOFs in ConditionRHS
    const unsigned int index_begin = rPairIndex * TMatrixSize;
    const unsigned int index_end  = index_begin + TMatrixSize;
    
    // Computing subrange
    subrange( rConditionRHS, index_begin, index_end ) += rPairRHS;
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MortarContactCondition<2, 2, false>::CalculateLocalRHS<12>(
    array_1d<double,12>& rPairRHS,
    GeneralVariables& rVariables,
    const ContactData<2, 2>& rContactData,
    const double& rIntegrationWeight
    )
{
    /* DEFINITIONS */    
    const array_1d<double, 2> N1  = rVariables.N_Slave;
    const array_1d<double, 2> N2  = rVariables.N_Master;
    const array_1d<double, 2> Phi = rVariables.Phi_LagrangeMultipliers;
    const double detJ = rVariables.DetJSlave;
    
    rPairRHS += rIntegrationWeight * Contact2D2N2N::ComputeGaussPointFrictionlessRHS(N1, N2, Phi, detJ, rContactData);
//     rPairRHS += rIntegrationWeight * Contact2D2N2N::ComputeGaussPointRHS(N1, N2, Phi, detJ, rVariables.mu, rContactData);
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MortarContactCondition<2, 3, false>::CalculateLocalRHS<18>(
    array_1d<double,18>& rPairRHS,
    GeneralVariables& rVariables,
    const ContactData<2, 3>& rContactData,
    const double& rIntegrationWeight
    )
{
//     /* DEFINITIONS */     
//     const array_1d<double, 3> N1  = rVariables.N_Slave;
//     const array_1d<double, 3> N2  = rVariables.N_Master;
//     const array_1d<double, 3> Phi = rVariables.Phi_LagrangeMultipliers;
//     const double detJ         = rVariables.DetJSlave;
//     
    // TODO: Finish this!!!
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MortarContactCondition<3, 3, false>::CalculateLocalRHS<27>(
    array_1d<double,27>& rPairRHS,
    GeneralVariables& rVariables,
    const ContactData<3, 3>& rContactData,
    const double& rIntegrationWeight
    )
{
    /* DEFINITIONS */
    const array_1d<double, 3> N1  = rVariables.N_Slave;
    const array_1d<double, 3> N2  = rVariables.N_Master;
    const array_1d<double, 3> Phi = rVariables.Phi_LagrangeMultipliers;
    const double detJ = rVariables.DetJSlave;
    
    rPairRHS += rIntegrationWeight * Contact3D3N3N::ComputeGaussPointFrictionlessRHS(N1, N2, Phi, detJ, rContactData);
//     rPairRHS += rIntegrationWeight * Contact3D3N3N::ComputeGaussPointRHS(N1, N2, Phi, detJ, rVariables.mu, rContactData);
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MortarContactCondition<3, 4, false>::CalculateLocalRHS<36>(
    array_1d<double,36>& rPairRHS,
    GeneralVariables& rVariables,
    const ContactData<3, 4>& rContactData,
    const double& rIntegrationWeight
    )
{
    /* DEFINITIONS */
    const array_1d<double, 4> N1  = rVariables.N_Slave;
    const array_1d<double, 4> N2  = rVariables.N_Master;
    const array_1d<double, 4> Phi = rVariables.Phi_LagrangeMultipliers;
    const double detJ = rVariables.DetJSlave;
    
    rPairRHS += rIntegrationWeight * Contact3D4N4N::ComputeGaussPointFrictionlessRHS(N1, N2, Phi, detJ, rContactData);
//     rPairRHS += rIntegrationWeight * Contact3D4N4N::ComputeGaussPointRHS(N1, N2, Phi, detJ, rVariables.mu, rContactData);
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MortarContactCondition<2, 2, true>::CalculateLocalRHS<16>(
    array_1d<double,16>& rPairRHS,
    GeneralVariables& rVariables,
    const ContactData<2, 2>& rContactData,
    const double& rIntegrationWeight
    )
{
    /* DEFINITIONS */    
    const array_1d<double, 2>  N1  = rVariables.N_Slave;
    const array_1d<double, 2>  N2  = rVariables.N_Master;
    const array_1d<double, 2>  Phi = rVariables.Phi_LagrangeMultipliers;
    const double detJ = rVariables.DetJSlave;
    
    rPairRHS += rIntegrationWeight * Contact2D2N2NDLM::ComputeGaussPointFrictionlessRHS(N1, N2, Phi, detJ, rContactData);
//     rPairRHS += rIntegrationWeight * Contact2D2N2NDLM::ComputeGaussPointRHS(N1, N2, Phi, detJ, rVariables.mu, rContactData);
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MortarContactCondition<2, 3, true>::CalculateLocalRHS<24>(
    array_1d<double,24>& rPairRHS,
    GeneralVariables& rVariables,
    const ContactData<2, 3>& rContactData,
    const double& rIntegrationWeight
    )
{
//     /* DEFINITIONS */     
//     const array_1d<double, 3> N1  = rVariables.N_Slave;
//     const array_1d<double, 3> N2  = rVariables.N_Master;
//     const array_1d<double, 3> Phi = rVariables.Phi_LagrangeMultipliers;
//     const double detJ         = rVariables.DetJSlave;
//     
    // TODO: Finish this!!!
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MortarContactCondition<3, 3, true>::CalculateLocalRHS<36>(
    array_1d<double,36>& rPairRHS,
    GeneralVariables& rVariables,
    const ContactData<3, 3>& rContactData,
    const double& rIntegrationWeight
    )
{
//     /* DEFINITIONS */
//     const array_1d<double, 3> N1  = rVariables.N_Slave;
//     const array_1d<double, 3> N2  = rVariables.N_Master;
//     const array_1d<double, 3> Phi = rVariables.Phi_LagrangeMultipliers;
//     const double detJ = rVariables.DetJSlave;
//     
//     rPairRHS += rIntegrationWeight * Contact3D3N3NDLM::ComputeGaussPointFrictionlessRHS(N1, N2, Phi, detJ, rContactData);
// //     rPairRHS += rIntegrationWeight * Contact3D3N3NDLM::ComputeGaussPointRHS(N1, N2, Phi, detJ, rVariables.mu, rContactData);
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MortarContactCondition<3, 4, true>::CalculateLocalRHS<48>(
    array_1d<double,48>& rPairRHS,
    GeneralVariables& rVariables,
    const ContactData<3,4>& rContactData,
    const double& rIntegrationWeight
    )
{
//     /* DEFINITIONS */
//     const array_1d<double, 4> N1  = rVariables.N_Slave;
//     const array_1d<double, 4> N2  = rVariables.N_Master;
//     const array_1d<double, 4> Phi = rVariables.Phi_LagrangeMultipliers;
//     const double detJ = rVariables.DetJSlave;
//     
//     rPairRHS += rIntegrationWeight * Contact3D4N4NDLM::ComputeGaussPointFrictionlessRHS(N1, N2, Phi, detJ, rContactData);
// //     rPairRHS += rIntegrationWeight * Contact3D4N4NDLM::ComputeGaussPointRHS(N1, N2, Phi, detJ, rVariables.mu, rContactData);
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM > // NOTE: Formulation taken from Mohamed Khalil work
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateDeltaDetJSlave(
   GeneralVariables& rVariables,
   ContactData<TDim, TNumNodes>& rContactData
   )
{
    if (TDim == 2)
    {
        // Fill up the elements corresponding to the slave DOFs - the rest remains zero
        for ( unsigned int i_slave = 0, i = 0; i_slave < TNumNodes; ++i_slave, i += TDim )
        {
            rContactData.DeltaJ_s[i    ] = rVariables.j_Slave( 0, 0 ) * rVariables.DN_De_Slave( i_slave, 0) / rVariables.DetJSlave;
            rContactData.DeltaJ_s[i + 1] = rVariables.j_Slave( 1, 0 ) * rVariables.DN_De_Slave( i_slave, 0) / rVariables.DetJSlave;
        }
    }
    else
    {
        const array_1d<double,TNumNodes>& DN_Dxi  = column( rVariables.DN_De_Slave, 0 );
        const array_1d<double,TNumNodes>& DN_Deta = column( rVariables.DN_De_Slave, 1 );
        
        const array_1d<double,TDim>& J_xi    = column( rVariables.j_Slave, 0 );
        const array_1d<double,TDim>& J_eta   = column( rVariables.j_Slave, 1 );
        
        const array_1d<double,TDim>& n = prod(trans(rContactData.Normal_m), rVariables.N_Slave);
        
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
            
            rContactData.DeltaJ_s[i    ] = inner_prod( n, column( Delta_Jxi_x_Jeta, 0 ) );
            rContactData.DeltaJ_s[i + 1] = inner_prod( n, column( Delta_Jxi_x_Jeta, 1 ) );
            rContactData.DeltaJ_s[i + 2] = inner_prod( n, column( Delta_Jxi_x_Jeta, 2 ) );
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM > // NOTE: Formulation taken from Mohamed Khalil work
bounded_matrix<double, TDim, TDim> MortarContactCondition<TDim,TNumNodes,TDoubleLM>::LocalDeltaNormal( // NOTE: Not the mean, look in the contact utilities 
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

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM > // NOTE: Formulation taken from Mohamed Khalil work
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateDeltaNormalTangentSlave(ContactData<TDim, TNumNodes>& rContactData)
{
    if (TDim == 2)
    {
        for ( unsigned int i_slave = 0, i = 0; i_slave < TNumNodes; ++i_slave, i += TDim )
        {
//             bounded_matrix<double, 2, 2> DeltaNormal = GetGeometry()[i_slave].GetValue(DELTA_NORMAL);
            bounded_matrix<double, 2, 2> DeltaNormal = this->LocalDeltaNormal(GetGeometry(), i_slave);
            for (unsigned i_dof = 0; i_dof < TDim; i_dof++) 
            {
                for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
                {
                    row(rContactData.Delta_Normal_s[i_slave * TDim + i_dof], i_node)      =   trans(column(DeltaNormal, i_dof)); 
                    rContactData.Delta_Tangent_xi_s[i_slave * TDim  + i_dof]( i_node, 0 ) =   DeltaNormal( i_dof, 1 );
                    rContactData.Delta_Tangent_xi_s[i_slave * TDim  + i_dof]( i_node, 1 ) = - DeltaNormal( i_dof, 0 );
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
            array_1d<double, 3> t_xi  = ZeroVector(3);
            array_1d<double, 3> t_eta = ZeroVector(3);
            ContactUtilities::NodalTangents( t_xi, t_eta, GetGeometry(), i_slave );
            
            const MatrixType I = IdentityMatrix( TDim, TDim );
            
            array_1d<double, 2> DN_De_j;
            if( TNumNodes == 3 )    // linear triangle element
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
                        
            const double norm_t_xi  = norm_2(t_xi );
            MatrixType DeltaTangentXi = ( I/norm_t_xi - outer_prod( t_xi, t_xi )/(norm_t_xi*norm_t_xi*norm_t_xi) ) * DN_De_j[0];
            
            const double norm_t_eta = norm_2(t_eta);
            MatrixType DeltaTangentEta = ( I/norm_t_eta - outer_prod( t_eta, t_eta )/(norm_t_eta*norm_t_eta*norm_t_eta) ) * DN_De_j[1];
            
            for (unsigned i_dof = 0; i_dof < TDim; i_dof++) 
            {
                for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
                {
                    row(rContactData.Delta_Normal_s[i_slave * TDim + i_dof],      i_node) = trans(column(DeltaNormal,     i_dof)); 
                    row(rContactData.Delta_Tangent_xi_s[i_slave * TDim + i_dof],  i_node) = trans(column(DeltaTangentXi,  i_dof)); 
                    row(rContactData.Delta_Tangent_eta_s[i_slave * TDim + i_dof], i_node) = trans(column(DeltaTangentEta, i_dof)); 
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM > // NOTE: Formulation taken from Mohamed Khalil work
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateDeltaNormalMaster(ContactData<TDim, TNumNodes>& rContactData)
{
    if (TDim == 2)
    {
        for ( unsigned int i_master = 0, i = 0; i_master < TNumNodes; ++i_master, i += TDim )
        {
//             bounded_matrix<double, 2, 2> DeltaNormal = GetGeometry[i_master].GetValue(DELTA_NORMAL);
            bounded_matrix<double, 2, 2> DeltaNormal = this->LocalDeltaNormal(rContactData.MasterGeometry, i_master);
            for (unsigned i_dof = 0; i_dof < TDim; i_dof++) 
            {
                for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
                {
                    row(rContactData.Delta_Normal_m[i_master * TDim + i_dof], i_node) = trans(column(DeltaNormal, i_dof)); 
                }
            }
        }
    }
    else
    {
        for ( unsigned int i_master = 0, i = 0; i_master < TNumNodes; ++i_master, i += TDim )
        {
//             bounded_matrix<double, 3, 3> DeltaNormal = GetGeometry[i_master]->GetValue(DELTA_NORMAL);
            const bounded_matrix<double, 3, 3> DeltaNormal = this->LocalDeltaNormal(rContactData.MasterGeometry, i_master);
            for (unsigned i_dof = 0; i_dof < TDim; i_dof++) 
            {
                for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
                {
                    row(rContactData.Delta_Normal_m[i_master * TDim + i_dof], i_node)  = trans(column(DeltaNormal, i_dof)); 
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateDeltaNAndDeltaGap(
   GeneralVariables& rVariables,
   ContactData<TDim, TNumNodes>& rContactData,
   const unsigned int CaseToCompute // NOTE: Just in 2D
   )
{
    /* IMPORTANT: The meaning of the CaseToCompute: 
     * CaseToCompute == 0: The slave is enterely integrated and just DeltaN2 is calculated
     * CaseToCompute == 1: Just the first node of the slave is inside, then we calculate DeltaN2 of the first and DeltaN1 of the second node
     * CaseToCompute == 2: Just the second node of the slave is inside, then we calculate DeltaN2 of the second and DeltaN1 of the first node
     * CaseToCompute == 3: The master is enterely integrated and just DeltaN1 is calculated
     */
    
    // Shape functions
    const array_1d<double, TNumNodes >  N1 = rVariables.N_Slave;
    const array_1d<double, TNumNodes >  N2 = rVariables.N_Master;
    
    // Coordinates
    const bounded_matrix<double, TNumNodes, TDim> u1 = rContactData.u1;
    const bounded_matrix<double, TNumNodes, TDim> X1 = rContactData.X1;
    const bounded_matrix<double, TNumNodes, TDim> u2 = rContactData.u2;
    const bounded_matrix<double, TNumNodes, TDim> X2 = rContactData.X2;
    
    // Normals
    const array_1d<double, TDim >  Normal_sg = prod(trans(rContactData.Normal_s), N1);
    const array_1d<double, TDim >  Normal_mg = prod(trans(rContactData.Normal_m), N2);
    
    const std::vector<bounded_matrix<double, TNumNodes, TDim>> DNormal_s = rContactData.Delta_Normal_s;
    const std::vector<bounded_matrix<double, TNumNodes, TDim>> DNormal_m = rContactData.Delta_Normal_m;
    
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
                  rContactData.DeltaN2[i_dof][0] =  (Deltax2g[0] + Deltax2g[1])/div2;
                  rContactData.DeltaN2[i_dof][1] =  - rContactData.DeltaN2[i_dof][0];
               }
            }
            else
            {
               if (TNumNodes == 3)
               {
                  rContactData.DeltaN2[i_dof][1] = - (mult1 * (Deltax2g[0] + Deltax2g[1]) + mult2 * (Deltax2g[0] + Deltax2g[2]))/div1;
                  rContactData.DeltaN2[i_dof][2] =   (mult3 *  Deltax2g[0] + mult4 * Deltax2g[1] + mult5 * Deltax2g[2])/div2;
                  rContactData.DeltaN2[i_dof][0] =  - rContactData.DeltaN2[i_dof][1] - rContactData.DeltaN2[i_dof][2];
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
                       rContactData.DeltaN2[i_dof][0] =  ( ((u2(1,0) + X2(1,0) + u2(1,1) + X2(1,1)) - x2g[0] - x2g[1]) + div2 * (Deltax2g[0] + Deltax2g[1]))/std::pow(div2, 2);
                   }
                   else 
                   {
                       rContactData.DeltaN2[i_dof][0] =  ((-(u2(0,0) + X2(0,0) + u2(0,1) + X2(0,1)) + x2g[0] + x2g[1]) + div2 * (Deltax2g[0] + Deltax2g[1]))/std::pow(div2, 2);
                   }
                  
                  rContactData.DeltaN2[i_dof][1] =  - rContactData.DeltaN2[i_dof][0];
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
                       rContactData.DeltaN2[i_dof][1] =  ((u2(1, 1) - u2(1, 2) - u2(2, 1) + u2(2, 2) + X2(1, 1) - 
                                                           X2(1, 2) - X2(2, 1) + X2(2, 2)) * (multmaster0*multmaster2 + 
                                                           multmaster1*multmaster3) - (div1) * (u2(0, 1) - u2(0, 2) + 
                                                           X2(0, 1) - X2(0, 2) - x2g[1] + x2g[2] + 
                                                           mult1 * (-1 + Deltax2g[0] + Deltax2g[1]) + 
                                                           mult2 * (-1 + Deltax2g[0] + Deltax2g[2])) )/std::pow(div1, 2);
                                                           
                       rContactData.DeltaN2[i_dof][2] =  ( -(u2(1, 1) - u2(1, 2) - u2(2, 1) + u2(2, 2) + X2(1, 1) - 
                                                            X2(1, 2) - X2(2, 1) + X2(2, 2)) * (coefmaster2 + (u2(0, 2) + X2(0, 2)) * 
                                                            multmaster5 + coefmaster1 - coefmaster5 - 
                                                            multmaster4 * multmaster6) + (div2) * (u2(1, 1) - u2(1, 2) + 
                                                            X2(1, 1) - X2(1, 2) - x2g[1] + x2g[2] + mult3 * Deltax2g[0] + 
                                                            mult4 * Deltax2g[1] + mult5 * Deltax2g[2]))/std::pow(div2, 2);
                   }
                   else if ((i_dof - TDim * TNumNodes) == 1)
                   {
                       rContactData.DeltaN2[i_dof][1] =  ((-u2(1, 0) - u2(1, 2) + u2(2, 0) + u2(2, 2) - X2(1, 0) - 
                                                            X2(1, 2) + X2(2, 0) + X2(2, 2)) * (multmaster0*multmaster2 + 
                                                            multmaster1*multmaster3) - (div1) * (-u2(0, 0) - u2(0, 2) - 
                                                            X2(0, 0) - X2(0, 2) + x2g[0] + x2g[2] + 
                                                            mult1 * (-1 + Deltax2g[0] + Deltax2g[1]) + 
                                                            mult2 * (Deltax2g[0] + Deltax2g[2])) )/std::pow(div1, 2);
                                                            
                       rContactData.DeltaN2[i_dof][2] =  (-(-u2(1, 0) - u2(1, 2) + u2(2, 0) + u2(2, 2) - X2(1, 0) - 
                                                             X2(1, 2) + X2(2, 0) + X2(2, 2)) * (coefmaster2 + (u2(0, 2) + X2(0, 2)) * 
                                                             multmaster5 + coefmaster1 - coefmaster5 - 
                                                             multmaster4 * multmaster6) + (div2) * (-u2(1, 0) - u2(1, 2) - 
                                                             X2(1, 0) - X2(1, 2) + x2g[0] + x2g[2] + mult3 * Deltax2g[0] + 
                                                             mult4 * Deltax2g[1] + mult5 * Deltax2g[2]) )/std::pow(div2, 2);
                   }
                   else if ((i_dof - TDim * TNumNodes) == 2)
                   {
                       rContactData.DeltaN2[i_dof][1] =  ( (u2(1, 0) + u2(1, 1) - u2(2, 0) - u2(2, 1) + X2(1, 0) + 
                                                            X2(1, 1) - X2(2, 0) - X2(2, 1)) * (multmaster0*multmaster2 + 
                                                            multmaster1*multmaster3) - (div1) * (u2(0, 0) + u2(0, 1) + 
                                                            X2(0, 0) + X2(0, 1) - x2g[0] - x2g[1] + 
                                                            mult1 * (Deltax2g[0] + Deltax2g[1]) + 
                                                            mult2 * (-1 + Deltax2g[0] + Deltax2g[2])))/std::pow(div1, 2);
                        
                       rContactData.DeltaN2[i_dof][2] =  (-(u2(1, 0) + u2(1, 1) - u2(2, 0) - u2(2, 1) + X2(1, 0) + 
                                                            X2(1, 1) - X2(2, 0) - X2(2, 1)) * (coefmaster2 + (u2(0, 2) + X2(0, 2)) * 
                                                            multmaster5 + coefmaster1 - coefmaster5 - 
                                                            multmaster4 * multmaster6) + (div2) * (u2(1, 0) + u2(1, 1) + 
                                                            X2(1, 0) + X2(1, 1) - x2g[0] - x2g[1] + mult3 * Deltax2g[0] + 
                                                            mult4 * Deltax2g[1] + mult5 * Deltax2g[2]) )/std::pow(div2, 2);
                   }
                   else if ((i_dof - TDim * TNumNodes) == 3)
                   {
                       rContactData.DeltaN2[i_dof][1] =  ( (-u2(0, 1) + u2(0, 2) + u2(2, 1) - u2(2, 2) - X2(0, 1) + 
                                                            X2(0, 2) + X2(2, 1) - X2(2, 2)) * (multmaster0*multmaster2 + 
                                                            multmaster1 * multmaster3) - (div1) * (+mult1 * (Deltax2g[0] + Deltax2g[1]) + 
                                                            mult2 * (Deltax2g[0] + Deltax2g[2])))/std::pow(div1, 2);
                                                            
                       rContactData.DeltaN2[i_dof][2] =  ( -(-u2(0, 1) + u2(0, 2) + u2(2, 1) - u2(2, 2) - X2(0, 1) + 
                                                              X2(0, 2) + X2(2, 1) - X2(2, 2)) * (coefmaster2 + (u2(0, 2) + X2(0, 2)) * 
                                                              multmaster5 + coefmaster1 - coefmaster5 - 
                                                              multmaster4 * multmaster6) + (div2) * (-u2(0, 1) + u2(0, 2) - 
                                                              X2(0, 1) + X2(0, 2) + x2g[1] - x2g[2] + mult3 * Deltax2g[0] + 
                                                              mult4 * Deltax2g[1] + mult5 * Deltax2g[2]))/std::pow(div2, 2);
                   }
                   else if ((i_dof - TDim * TNumNodes) == 4)
                   {
                       rContactData.DeltaN2[i_dof][1] =  ( multmaster0*(multmaster0*multmaster2 + 
                                                           multmaster1 * multmaster3) - (div1) * (+mult1 * (Deltax2g[0] + Deltax2g[1]) + 
                                                           mult2 * (Deltax2g[0] + Deltax2g[2])))/std::pow(div1, 2);
                                                           
                       rContactData.DeltaN2[i_dof][2] =  ( +mult1 * (coefmaster2 + (u2(0, 2) + X2(0, 2)) * multmaster5 + 
                                                            coefmaster1 - coefmaster5 - multmaster4 * multmaster6) + (div2) * (coefmaster3 + 
                                                            mult3 * Deltax2g[0] + mult4 * Deltax2g[1] + mult5 * Deltax2g[2]))/std::pow(div2, 2);
                   }
                   else if ((i_dof - TDim * TNumNodes) == 5)
                   {
                       rContactData.DeltaN2[i_dof][1] =  ( multmaster1 *(multmaster0*multmaster2 + 
                                                           multmaster1 * multmaster3) - (div1) * (+mult1 * (Deltax2g[0] + Deltax2g[1]) + 
                                                           mult2 * (Deltax2g[0] + Deltax2g[2])))/std::pow(div1, 2);
                                                           
                       rContactData.DeltaN2[i_dof][2] =  ( -multmaster1 *  (coefmaster2 + (u2(0, 2) + X2(0, 2)) * multmaster5 +
                                                            coefmaster1 - coefmaster5 -  multmaster4 * multmaster6) + (div2) * (coefmaster4 + 
                                                            mult3 * Deltax2g[0] + mult4 * Deltax2g[1] + mult5 * Deltax2g[2]))/std::pow(div2, 2);
                   }
                   else if ((i_dof - TDim * TNumNodes) == 6)
                   {
                       rContactData.DeltaN2[i_dof][1] =  (mult3 *(multmaster0*multmaster2 + 
                                                          multmaster1*multmaster3) - (div1) * (-u2(0, 1) + u2(0, 2) - 
                                                          X2(0, 1) + X2(0, 2) + x2g[1] - x2g[2] + 
                                                          mult1 * (Deltax2g[0] + Deltax2g[1]) + 
                                                          mult2 * (Deltax2g[0] + Deltax2g[2])) )/std::pow(div1, 2);
                                                          
                       rContactData.DeltaN2[i_dof][2] =  ( -mult3 * (coefmaster2 + (u2(0, 2) + X2(0, 2)) * multmaster5 + 
                                                            coefmaster1 - coefmaster5 - multmaster4 * multmaster6) + (div2) * (mult3 * Deltax2g[0] + 
                                                            mult4 * Deltax2g[1] + mult5 * Deltax2g[2]))/std::pow(div2, 2);
                   }
                   else if ((i_dof - TDim * TNumNodes) == 7)
                   {
                       rContactData.DeltaN2[i_dof][1] =  (mult4 *(multmaster0*multmaster2 + 
                                                          multmaster1*multmaster3) - (div1) * (coefmaster3 + 
                                                          mult1 * (Deltax2g[0] + Deltax2g[1]) + 
                                                          mult2 * (Deltax2g[0] + Deltax2g[2])))/std::pow(div1, 2);
                                                          
                       rContactData.DeltaN2[i_dof][2] =  ( -mult4 * (coefmaster2 + (u2(0, 2) + X2(0, 2)) * multmaster5 + 
                                                            coefmaster1 - coefmaster5 - multmaster4 * multmaster6) + (div2) * (mult3 * Deltax2g[0] + 
                                                            mult4 * Deltax2g[1] + mult5 * Deltax2g[2]))/std::pow(div2, 2);
                   }
                   else if ((i_dof - TDim * TNumNodes) == 8)
                   {
                       rContactData.DeltaN2[i_dof][1] =  (mult5 *(multmaster0*multmaster2 + 
                                                          multmaster1*multmaster3) - (div1) * (coefmaster4 + 
                                                          mult1 * (Deltax2g[0] + Deltax2g[1]) + 
                                                          mult2 * (Deltax2g[0] + Deltax2g[2])))/std::pow(div1, 2);
                                                          
                       rContactData.DeltaN2[i_dof][2] =  (-mult5 * (coefmaster2 + (u2(0, 2) + X2(0, 2)) * multmaster5 + 
                                                          coefmaster1 - coefmaster5 - multmaster4 * multmaster6) + (div2) * (mult3 * Deltax2g[0] + 
                                                          mult4 * Deltax2g[1] + mult5 * Deltax2g[2]) )/std::pow(div2, 2);
                   }

                  rContactData.DeltaN2[i_dof][0] =  - rContactData.DeltaN2[i_dof][1] - rContactData.DeltaN2[i_dof][2];
               }
            }
         }
      }
    }

    /* Calculate Delta Gap */ // TODO: Consider DeltaN1 
    // Derivatives slave
    for (unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++)
    {
        for (unsigned int i_dim = 0; i_dim < TDim; i_dim++)
        {
            const unsigned int i_dof = i_slave * TDim + i_dim;
            
            const array_1d<double, TNumNodes > DN2 = rContactData.DeltaN2[i_dof];
            
            const array_1d<double, TDim >  DeltaNormal = prod(trans(rContactData.Delta_Normal_s[i_dof]), N1);
            
            const array_1d<double, TDim >  RelPos = prod(trans(X1 + u1), N1) - prod(trans(X2 + u2), N2);
            
            rContactData.DeltaGap[i_dof] = - inner_prod(DeltaNormal, RelPos) - Normal_sg[i_dim] * N1[i_slave] + inner_prod(Normal_sg, prod(trans(X2 + u2), DN2));
        }
    }
    
    // Derivatives master
    for (unsigned int i_master = 0; i_master < TNumNodes; i_master++)
    {
        for (unsigned int i_dim = 0; i_dim < TDim; i_dim++)
        {
            const unsigned int i_dof = (TNumNodes + i_master) * TDim + i_dim;
            
            const array_1d<double, TNumNodes >  DN2 = rContactData.DeltaN2[i_dof];
            
            rContactData.DeltaGap[i_dof] = Normal_sg[i_dim] * N2[i_master] +  inner_prod(Normal_sg, prod(trans(X2 + u2), DN2));
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateDeltaPhi(
   GeneralVariables& rVariables,
   ContactData<TDim, TNumNodes>& rContactData
   )
{
    // Shape functions
    const Vector N1 = rVariables.N_Slave;
    
    for (unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++)
    {
        for (unsigned int i_dim = 0; i_dim < TDim; i_dim++)
        {
            const unsigned int i_dof = i_slave * TDim + i_dim;
            
            rContactData.DeltaPhi[i_dof] = prod(rContactData.DeltaAe[i_dof], N1);;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
double MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateCtanAndDeltaCtan(
   GeneralVariables& rVariables,
   ContactData<TDim, TNumNodes>& rContactData
   )
{   
    // NOTE: I calculate the increment of slip manually, Gitterle uses the Mortar operator, but I don't calculate them, technically is the same working with GP 
    
    // Master segment info
    const GeometryType& CurrentMasterElement = rVariables.GetMasterElement( );
    
    const double mu = rVariables.mu; 
    const double epsilon_normal  = rContactData.epsilon_normal;
    const double epsilon_tangent = rContactData.epsilon_tangent;

    /****** CALCULATING THE COMPONENTS OF THE COMPLEMENTARY FUNCTION ******/
    
    // Local kinematic update
    const array_1d<double, TNumNodes> N1  = rVariables.N_Slave;
    const array_1d<double, TNumNodes> N2  = rVariables.N_Master;
    const array_1d<double, TNumNodes> Phi = rVariables.Phi_LagrangeMultipliers;
    
    // The LM of the node
    const array_1d<double, TDim> LMGP = prod(trans(rContactData.LagrangeMultipliers), Phi);

    // Calculate the augmented normal LM
    const double AugmentedNormalLM = this->CalculateAugmentedNormalLM(rVariables, rContactData, inner_prod(rContactData.Gaps, N1));
    
    if (AugmentedNormalLM < 0.0) // Compression (active GP)
    {
        // The velocities matrices
        const bounded_matrix<double, TNumNodes, TDim> v1 = rContactData.v1; 
        const bounded_matrix<double, TNumNodes, TDim> v2 = rContactData.v2;
        
        // Calculating the tangent vectors
        const array_1d<double, TDim> NormalsGP     = prod(trans(rContactData.Normal_s),      N1);
        const array_1d<double, TDim> TangentXisGP  = prod(trans(rContactData.Tangent_xi_s),  N1);
        const array_1d<double, TDim> TangentEtasGP = prod(trans(rContactData.Tangent_eta_s), N1);
            
        // Calculate the augmented tangent LM
        array_1d<double, (TDim - 1)> node_slip;
        const array_1d<double, (TDim - 1)> AugmentedTangentLM = this->CalculateAugmentedTangentLM(rVariables, rContactData, CurrentMasterElement, node_slip);
        
        if ((std::abs(AugmentedTangentLM[0]) + mu * AugmentedNormalLM) < 0.0) // Stick
        {
            rContactData.Ctan[0] =  - mu * AugmentedNormalLM * epsilon_tangent * node_slip[0];
        }
        else
        {
            rContactData.Ctan[0] = std::abs(AugmentedTangentLM[0]) * inner_prod(TangentXisGP, LMGP) - mu * AugmentedNormalLM * AugmentedTangentLM[0];
        }
        
        if (TDim == 3)
        {
            if ((std::abs(AugmentedTangentLM[1]) + mu * AugmentedNormalLM) < 0.0) // Stick
            {
                rContactData.Ctan[1] =  - mu * AugmentedNormalLM * epsilon_tangent * node_slip[1];
            }
            else
            {
                rContactData.Ctan[1] = std::abs(AugmentedTangentLM[1]) * inner_prod(TangentEtasGP, LMGP) - mu * AugmentedNormalLM * AugmentedTangentLM[1];
            }
        }
        
        /****** CALCULATING THE COMPONENTS OF THE DERIVATIVE COMPLEMENTARY FUNCTION ******/

        // Calculate common components
        const array_1d<double, TDim> VectorIntegrationPointSlip = rContactData.Dt * (prod(trans(v1), N1) - prod(trans(v2), N2));
        array_1d<double, TDim> DeltaVectorIntegrationPointSlip;
        
        // Derivatives slave
        for (unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++)
        {
            for (unsigned int i_dim = 0; i_dim < TDim; i_dim++)
            {
                const unsigned int i_dof = i_slave * TDim + i_dim;
                
                array_1d<double, TDim> aux_vec = ZeroVector(TDim);
                aux_vec[i_dim] = N1[i_slave];
                
                const array_1d<double, TNumNodes> DeltaN1 = rContactData.DeltaN1[i_dof];
                const array_1d<double, TNumNodes> DeltaN2 = rContactData.DeltaN2[i_dof];
                DeltaVectorIntegrationPointSlip = aux_vec + rContactData.Dt * (prod(trans(v1), DeltaN1) - prod(trans(v2), DeltaN2));
                
                const array_1d<double, TDim> DeltaNormalsGP    = prod(trans(rContactData.Delta_Normal_s[i_dof]),     N1);
                const double DeltaLMNormal = inner_prod(DeltaNormalsGP, LMGP);
                const double DeltaGap = rContactData.DeltaGap[i_dof];
                // Delta slip
                const array_1d<double, TDim> DeltaTangentXisGP = prod(trans(rContactData.Delta_Tangent_xi_s[i_dof]), N1);
                array_1d<double, (TDim - 1)> DeltaIntegrationPointSlip;
                DeltaIntegrationPointSlip[0]  = inner_prod(VectorIntegrationPointSlip, DeltaTangentXisGP);
                DeltaIntegrationPointSlip[0] += inner_prod(DeltaVectorIntegrationPointSlip, TangentXisGP);
                
                if ((std::abs(AugmentedTangentLM[0]) + mu * AugmentedNormalLM) < 0.0) // Stick
                {
                    rContactData.DeltaCtan[i_dof][0] = - mu * AugmentedNormalLM * epsilon_tangent * DeltaIntegrationPointSlip[0] // NOTE: (Gitterle p.75)
                                                       - mu * (DeltaLMNormal + epsilon_normal * DeltaGap) * epsilon_tangent * node_slip[0]; 
                }
                else // Slip
                {
                    const double      LMXi = inner_prod(     TangentXisGP, LMGP);
                    const double DeltaLMXi = inner_prod(DeltaTangentXisGP, LMGP);
                    const double SignTangPress = boost::math::sign(AugmentedTangentLM[0]);
                    rContactData.DeltaCtan[i_dof][0] = std::abs(AugmentedTangentLM[0]) * DeltaLMXi 
                                                       + (SignTangPress * LMXi - mu * AugmentedNormalLM) * (DeltaLMXi + epsilon_tangent * DeltaIntegrationPointSlip[0])
                                                       - mu * (DeltaLMNormal + epsilon_normal * DeltaGap) * AugmentedTangentLM[0];
                }
                
                if (TDim == 3)
                {
                    const array_1d<double, TDim> DeltaTangentEtasGP = prod(trans(rContactData.Delta_Tangent_eta_s[i_dof]), N1);
                    DeltaIntegrationPointSlip[1]  = inner_prod(VectorIntegrationPointSlip, DeltaTangentEtasGP);
                    DeltaIntegrationPointSlip[1] += inner_prod(DeltaVectorIntegrationPointSlip, TangentEtasGP);
                    if ((std::abs(AugmentedTangentLM[1]) + mu * AugmentedNormalLM) < 0.0) // Stick
                    {
                        rContactData.DeltaCtan[i_dof][1] = - mu * AugmentedNormalLM * epsilon_tangent * DeltaIntegrationPointSlip[1] // NOTE: (Gitterle p.75)
                                                           - mu * (DeltaLMNormal + epsilon_normal * DeltaGap) * epsilon_tangent * node_slip[1]; 
                    }
                    else // Slip
                    {
                        const double      LMEta = inner_prod(     TangentEtasGP, LMGP);
                        const double DeltaLMEta = inner_prod(DeltaTangentEtasGP, LMGP);
                        const double SignTangPress = boost::math::sign(AugmentedTangentLM[1]);
                        rContactData.DeltaCtan[i_dof][1] = std::abs(AugmentedTangentLM[1]) * DeltaLMEta 
                                                       + (SignTangPress *  LMEta - mu * AugmentedNormalLM ) * (DeltaLMEta + epsilon_tangent * DeltaIntegrationPointSlip[1])
                                                       - mu * (DeltaLMNormal + epsilon_normal * DeltaGap) * AugmentedTangentLM[1];
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
                
                array_1d<double, TDim> aux_vec = ZeroVector(TDim);
                aux_vec[i_dim] = N2[i_master];
                
                const array_1d<double, TNumNodes> DeltaN1 = rContactData.DeltaN1[i_dof];
                const array_1d<double, TNumNodes> DeltaN2 = rContactData.DeltaN2[i_dof];
                DeltaVectorIntegrationPointSlip = - aux_vec + rContactData.Dt * (prod(trans(v1), DeltaN1) - prod(trans(v2), DeltaN2));
                
                const double DeltaGap = rContactData.DeltaGap[i_dof];
                array_1d<double, (TDim - 1)> DeltaIntegrationPointSlip;
                DeltaIntegrationPointSlip[0] = inner_prod(DeltaVectorIntegrationPointSlip, TangentXisGP);
                
                if ((std::abs(AugmentedTangentLM[0]) + mu * AugmentedNormalLM) < 0.0) // Stick
                {
                    rContactData.DeltaCtan[i_dof][0] = - mu * AugmentedNormalLM * epsilon_tangent * DeltaIntegrationPointSlip[0]// NOTE: (Gitterle p.75)
                                                       - mu * ( + epsilon_normal * DeltaGap) * epsilon_tangent * node_slip[0]; 
                }
                else // Slip
                {
                    const double      LMXi = inner_prod(     TangentXisGP, LMGP);
                    const double SignTangPress = boost::math::sign(AugmentedTangentLM[0]);
                    rContactData.DeltaCtan[i_dof][0] = (SignTangPress * LMXi - mu * AugmentedNormalLM) * ( + epsilon_tangent * DeltaIntegrationPointSlip[0])
                                                       - mu * ( + epsilon_normal * DeltaGap) * AugmentedTangentLM[0];
                }
                
                if (TDim == 3)
                {
                    DeltaIntegrationPointSlip[1] = inner_prod(DeltaVectorIntegrationPointSlip, TangentEtasGP);
                    if ((std::abs(AugmentedTangentLM[1]) + mu * AugmentedNormalLM) < 0.0) // Stick
                    {
                        rContactData.DeltaCtan[i_dof][1] = - mu * AugmentedNormalLM * epsilon_tangent * DeltaIntegrationPointSlip[1]// NOTE: (Gitterle p.75)
                                                           - mu * ( + epsilon_normal * DeltaGap) * epsilon_tangent * node_slip[1]; 
                    }
                    else // Slip
                    {
                        const double LMEta = inner_prod(TangentEtasGP, LMGP);
                        const double SignTangPress = boost::math::sign(AugmentedTangentLM[1]);
                        rContactData.DeltaCtan[i_dof][1] = (SignTangPress * LMEta - mu * AugmentedNormalLM) * ( + epsilon_tangent * DeltaIntegrationPointSlip[1])
                                                           - mu * ( + epsilon_normal * DeltaGap) * AugmentedTangentLM[1];
                    }
                }
            }
            
        }
        
        // Derivatives lagrange multiplier
        for (unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++)
        {
            for (unsigned int i_dim = 0; i_dim < TDim; i_dim++)
            {
                const unsigned int i_dof = (2 * TNumNodes + i_slave) * TDim + i_dim;
                
                array_1d<double, TDim> aux_vec = ZeroVector(TDim);
                aux_vec[i_dim] = Phi[i_slave];
                const double DeltaLMNormal = inner_prod(NormalsGP, aux_vec);

                if ((std::abs(AugmentedTangentLM[0]) + mu * AugmentedNormalLM) < 0.0) // Stick
                {
                    rContactData.DeltaCtan[i_dof][0] = mu * (DeltaLMNormal) * epsilon_tangent * node_slip[0]; 
                }
                else // Slip
                {
                    const double      LMXi = inner_prod(TangentXisGP, LMGP);
                    const double DeltaLMXi = inner_prod(TangentXisGP, aux_vec);
                    const double SignTangPress = boost::math::sign(AugmentedTangentLM[0]);
                    rContactData.DeltaCtan[i_dof][0] = DeltaLMXi * (std::abs(AugmentedTangentLM[0])
                                                      + SignTangPress * LMXi + mu * AugmentedNormalLM)
                                                      + mu * DeltaLMNormal * AugmentedTangentLM[0];
                }
                if (TDim == 3)
                {
                    if ((std::abs(AugmentedTangentLM[1]) + mu * AugmentedNormalLM) < 0.0) // Stick
                    {
                        rContactData.DeltaCtan[i_dof][1] = mu * (DeltaLMNormal) * epsilon_tangent * node_slip[1]; 
                    }
                    else // Slip
                    {
                        const double      LMEta = inner_prod(TangentEtasGP, LMGP);
                        const double DeltaLMEta = inner_prod(TangentEtasGP, aux_vec);
                        const double SignTangPress = boost::math::sign(AugmentedTangentLM[1]);
                        rContactData.DeltaCtan[i_dof][1] = DeltaLMEta * (std::abs(AugmentedTangentLM[1])  
                                                          + SignTangPress * LMEta + mu * AugmentedNormalLM) 
                                                          + mu * DeltaLMNormal * AugmentedTangentLM[1];
                    }
                }
            }
        }
    }
    
    return AugmentedNormalLM;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& CurrentProcessInfo 
    )
{
    KRATOS_TRY;   
    
    const std::vector<contact_container> all_conditions = *( this->GetValue( CONTACT_CONTAINERS ) );
    
    // Calculates the size of the system
    const unsigned int condition_size = TDim * ((2 + TDoubleLM) * TNumNodes + TNumNodes) * all_conditions.size(); 
    
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
            rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_X ).EquationId( );
            rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Y ).EquationId( );
            if (TDim == 3)
            {
                rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Z ).EquationId( );
            }
        }
        
        if (TDoubleLM == true)
        {
            // Slave Nodes  Lambda Equation IDs
            for ( unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++ )
            {
                NodeType& slave_node = GetGeometry()[ i_slave ];
                rResult[index++] = slave_node.GetDof( DOUBLE_LM_X ).EquationId( );
                rResult[index++] = slave_node.GetDof( DOUBLE_LM_Y ).EquationId( );
                if (TDim == 3)
                {
                    rResult[index++] = slave_node.GetDof( DOUBLE_LM_Z ).EquationId( );
                }
            }
        }
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim, TNumNodes, TDoubleLM>::GetDofList(
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo 
)
{
    KRATOS_TRY;
    
    const std::vector<contact_container> all_conditions = *( this->GetValue( CONTACT_CONTAINERS ) );
    
    // Calculates the size of the system
    const unsigned int condition_size = TDim * ((2 + TDoubleLM) * TNumNodes + TNumNodes) * all_conditions.size(); 
    
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
            rConditionalDofList[index++] =slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_X );
            rConditionalDofList[index++] =slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Y );
            if (TDim == 3)
            {
                rConditionalDofList[index++] =slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Z );
            }
        }

        if (TDoubleLM == true)
        {
            // Slave Nodes Lambda Equation IDs
            for ( unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++ )
            {
                NodeType& slave_node = GetGeometry()[ i_slave ];
                rConditionalDofList[index++] =slave_node.pGetDof( DOUBLE_LM_X );
                rConditionalDofList[index++] =slave_node.pGetDof( DOUBLE_LM_Y );
                if (TDim == 3)
                {
                    rConditionalDofList[index++] =slave_node.pGetDof( DOUBLE_LM_Z );
                }
            }
        }
    }
    
    KRATOS_CATCH( "" );
}


//******************************* GET DOUBLE VALUE *********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::GetValueOnIntegrationPoints( 
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//******************************* GET ARRAY_1D VALUE *******************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::GetValueOnIntegrationPoints( 
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//******************************* GET VECTOR VALUE *********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::GetValueOnIntegrationPoints( 
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateOnIntegrationPoints( 
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
    ContactData<TDim, TNumNodes> rContactData;
    
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
    
    this->InitializeContactData(rContactData, rCurrentProcessInfo);
    
    this->CalculateAeAndDeltaAe(rContactData, rVariables, rCurrentProcessInfo);
    
    for (unsigned int PairIndex = 0; PairIndex < mPairSize; ++PairIndex)
    {   
        // Initialize general variables for the current master element
        this->InitializeGeneralVariables( rVariables, rCurrentProcessInfo, PairIndex );
        
        // Update the contact data
        this->UpdateContactData(rContactData, PairIndex);
        
        // Master segment info
        const GeometryType& CurrentMasterElement = rVariables.GetMasterElement( );

        // Calculating the values
        for ( unsigned int PointNumber = 0; PointNumber < number_of_integration_pts; PointNumber++ )
        {
            // Calculate the kinematic variables
            const bool inside = this->CalculateKinematics( rVariables, rContactData, PointNumber, PairIndex, integration_points );

            if (inside == true)
            {
                // Calculate the gap in the integration node and check tolerance
                const double IntegrationPointGap = inner_prod(rContactData.Gaps, rVariables.N_Slave);
            
                if (rVariable == NORMAL_CONTACT_STRESS || rVariable == TANGENTIAL_CONTACT_STRESS || rVariable == SLIP_GP)
                {
                    // The normal LM augmented
                    const double AugmentedNormalLM = this->CalculateAugmentedNormalLM(rVariables, rContactData, IntegrationPointGap);
                    
                    if (AugmentedNormalLM < 0)
                    {
                        if ( rVariable == NORMAL_CONTACT_STRESS )
                        {
                            rOutput[PointNumber] += AugmentedNormalLM;
                        }
                        else
                        {
                            // The slip of th GP
                            double IntegrationPointSlip;
                        
                            // The tangent LM augmented
                            const double AugmentedTangentLM = this->CalculateAugmentedTangentLM(rVariables, rContactData, CurrentMasterElement, IntegrationPointSlip);
                            
                            if ( rVariable == TANGENTIAL_CONTACT_STRESS )
                            {
                                rOutput[PointNumber] += AugmentedTangentLM;
                            }
                            else if (rVariable == SLIP_GP )
                            {
                                rOutput[PointNumber] += IntegrationPointSlip; // TODO: Maybe it is necessary to check if it is the smallest of all the gaps 
                            }
                        }
                    }
                }
                else if (rVariable == GAP_GP )
                {
                    rOutput[PointNumber] += IntegrationPointGap; // TODO: Maybe it is necessary to check if it is the smallest of all the gaps 
                }
            }
        }
    } 
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateOnIntegrationPoints( 
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector< array_1d<double, 3 > >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
    
    // Create and initialize condition variables:
    GeneralVariables rVariables;
    
    // Initialize the current contact data
    ContactData<TDim, TNumNodes> rContactData;
    
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
    
    this->InitializeContactData(rContactData, rCurrentProcessInfo);
    
    this->CalculateAeAndDeltaAe(rContactData, rVariables, rCurrentProcessInfo);
    
    if (rVariable == NORMAL_GP || rVariable == TANGENT_GP)
    {
        if (mPairSize > 0)
        {
            // Initialize general variables for the current master element
            this->InitializeGeneralVariables( rVariables, rCurrentProcessInfo, 0 ); // NOTE: The pair does not matter
            // Initialize Ae
            rContactData.InitializeDeltaAeComponents(); // NOTE: We need at least a zero Ae matrix to avoid problems in the CalculateKinematics
            
            // Calculating the values
            for ( unsigned int PointNumber = 0; PointNumber < number_of_integration_pts; PointNumber++ )
            {
                // Calculate the kinematic variables
                const bool inside = this->CalculateKinematics( rVariables, rContactData, PointNumber, 0, integration_points ); // NOTE: The pair does not matter
                
                if (inside == true)
                {
                    if (rVariable == NORMAL_GP)
                    {   
                        const array_1d<double, 3> NormalGP  = prod(trans(rContactData.Normal_s), rVariables.N_Slave);
                        rOutput[PointNumber] = NormalGP;
                    }
                    else
                    {
                        const array_1d<double, 3> TangentGP = prod(trans(rContactData.Tangent_xi_s), rVariables.N_Slave); 
                        rOutput[PointNumber] = TangentGP;
                    }
                }
            }
        }
    }
    else
    {
        for (unsigned int PairIndex = 0; PairIndex < mPairSize; ++PairIndex)
        {   
            // Initialize general variables for the current master element
            this->InitializeGeneralVariables( rVariables, rCurrentProcessInfo, PairIndex );
            
            // Update the contact data
            this->UpdateContactData(rContactData, PairIndex);
            
            // Master segment info
            const GeometryType& CurrentMasterElement = rVariables.GetMasterElement( );
            
            // Calculating the values
            for ( unsigned int PointNumber = 0; PointNumber < number_of_integration_pts; PointNumber++ )
            {
                // Calculate the kinematic variables
                const bool inside = this->CalculateKinematics( rVariables, rContactData, PointNumber, PairIndex, integration_points );
            
                if (inside == true)
                {
                    // Calculate the gap in the integration node and check tolerance
                    const double IntegrationPointGap = inner_prod(rContactData.Gaps, rVariables.N_Slave);
                
                    if (rVariable == NORMAL_CONTACT_STRESS_GP || rVariable == TANGENTIAL_CONTACT_STRESS_GP)
                    {
                        // The normal LM augmented
                        const double AugmentedNormalLM = this->CalculateAugmentedNormalLM(rVariables, rContactData, IntegrationPointGap);
                        
                        if (AugmentedNormalLM < 0)
                        {
                            if ( rVariable == NORMAL_CONTACT_STRESS_GP )
                            {
                                const array_1d<double, 3> NormalGP  = prod(trans(rContactData.Normal_s), rVariables.N_Slave);
                                rOutput[PointNumber] += AugmentedNormalLM * NormalGP;
                            }
                            else
                            {
                                // The slip of th GP
                                double IntegrationPointSlip;
                            
                                // The tangent LM augmented
                                const double AugmentedTangentLM = this->CalculateAugmentedTangentLM(rVariables, rContactData, CurrentMasterElement, IntegrationPointSlip);
                                
                                const array_1d<double, 3> TangentGP = prod(trans(rContactData.Tangent_xi_s), rVariables.N_Slave); // TODO: Extend to 3D case
                                rOutput[PointNumber] += AugmentedTangentLM * TangentGP;
                            }
                        }
                    }
                }
            }
        }
    } 
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateOnIntegrationPoints( 
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

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
double MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateAugmentedNormalLM(
    const GeneralVariables& rVariables,
    const ContactData<TDim, TNumNodes>& rContactData,
    const double& IntegrationPointGap
    )
{
    // The LM of the GP 
    const array_1d<double, TDim> LMGP = prod(trans(rContactData.LagrangeMultipliers), rVariables.Phi_LagrangeMultipliers);
    
    // The normal of the GP
    const array_1d<double, TDim> NormalGP = prod(trans(rContactData.Normal_s), rVariables.N_Slave);
    
    // The LM in the tangent direction
    double AugmentedNormalLM  = inner_prod(NormalGP, LMGP);// + rContactData.epsilon_normal * IntegrationPointGap;
                
    if (IntegrationPointGap < 0.0) // Penetration (TODO: Think about this!!!)
    {
        AugmentedNormalLM += rContactData.epsilon_normal * IntegrationPointGap;
    }
                
    return AugmentedNormalLM;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
double MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateAugmentedTangentLM(
    const GeneralVariables& rVariables,
    const ContactData<TDim, TNumNodes>& rContactData,
    const GeometryType& CurrentMasterElement, 
    double& IntegrationPointSlip
    )
{    
    // The LM of the GP 
    const array_1d<double, TDim> LMGP = prod(trans(rContactData.LagrangeMultipliers), rVariables.Phi_LagrangeMultipliers);
    
    // The tangents of the GP
    const array_1d<double, TDim> TangentXisGP  = prod(trans(rContactData.Tangent_xi_s),  rVariables.N_Slave);
    const array_1d<double, TDim> TangentEtasGP = prod(trans(rContactData.Tangent_eta_s), rVariables.N_Slave);
    
    // The tangential LM
    const double TangentXiLM  = inner_prod(TangentXisGP, LMGP);
    const double TangentEtaLM = inner_prod(TangentEtasGP, LMGP);
    const double TangentLM = std::sqrt(TangentXiLM * TangentXiLM + TangentEtaLM * TangentEtaLM); 
    
    // The resultant direction of the LM
    const array_1d<double, TDim>& TangentGP =  (TangentXiLM * TangentXisGP + TangentEtaLM * TangentEtasGP)/TangentLM; // NOTE: This is the direction of the slip (using the LM as reference)
    
    // The velocities matrices
    const bounded_matrix<double, TNumNodes, TDim> v1 = rContactData.v1; 
    const bounded_matrix<double, TNumNodes, TDim> v2 = rContactData.v2;
    
    // The slip of the LM
    const array_1d<double, TDim> vector_IntegrationPointSlip = rContactData.Dt * (prod(trans(v1), rVariables.N_Slave) - prod(trans(v2), rVariables.N_Master));
    IntegrationPointSlip = inner_prod(vector_IntegrationPointSlip, TangentGP);

    // The augmented LM in the tangent direction
    const double AugmentedTangentLM = TangentLM + rContactData.epsilon_tangent * IntegrationPointSlip; 

    return AugmentedTangentLM;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
array_1d<double, (TDim - 1)> MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CalculateAugmentedTangentLM(
    const GeneralVariables& rVariables,
    const ContactData<TDim, TNumNodes>& rContactData,
    const GeometryType& CurrentMasterElement, 
    array_1d<double, (TDim - 1)>& IntegrationPointSlip
    )
{  
    array_1d<double, (TDim - 1)> AugmentedTangentLM;
    
    // The LM of the GP 
    const array_1d<double, TDim> LMGP = prod(trans(rContactData.LagrangeMultipliers), rVariables.Phi_LagrangeMultipliers);
    
    // The tangents of the GP
    const array_1d<double, TDim> TangentXisGP  = prod(trans(rContactData.Tangent_xi_s), rVariables.N_Slave);
    
    // The tangential LM
    const double TangentXiLM  = inner_prod(TangentXisGP, LMGP);
    
    // The velocities matrices
    const bounded_matrix<double, TNumNodes, TDim> v1 = rContactData.v1; 
    const bounded_matrix<double, TNumNodes, TDim> v2 = rContactData.v2;
    
    // The slip of the LM
    const array_1d<double, TDim> VectorIntegrationPointSlip = rContactData.Dt * (prod(trans(v1), rVariables.N_Slave) - prod(trans(v2), rVariables.N_Master));
    IntegrationPointSlip[0] = inner_prod(VectorIntegrationPointSlip, TangentXisGP);

    // The augmented LM in the tangent direction
    AugmentedTangentLM[0] = TangentXiLM + rContactData.epsilon_tangent * IntegrationPointSlip[0]; 
    
    if (TDim == 3)
    {
        // The tangents of the GP
        const array_1d<double, TDim> TangentEtasGP = prod(trans(rContactData.Tangent_eta_s), rVariables.N_Slave);
        
        // The tangential LM
        const double TangentEtaLM = inner_prod(TangentEtasGP, LMGP);
        
        // The slip of the LM
        IntegrationPointSlip[1] = inner_prod(VectorIntegrationPointSlip, TangentEtasGP);

        // The augmented LM in the tangent direction
        AugmentedTangentLM[1] = TangentEtaLM + rContactData.epsilon_tangent * IntegrationPointSlip[1]; 
    }

    return AugmentedTangentLM;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
unsigned int MortarContactCondition<TDim,TNumNodes,TDoubleLM>::CaseToCompute(const GeometryType::IntegrationPointsArrayType& integration_points)
{
    unsigned int compute = 0;
    
    if (TDim == 2)
    {
        if (TNumNodes == 2)
        {
            if (integration_points.size() > 0)
            {
                const double xi_0 = integration_points[0].Coordinate(1);
                const double xi_1 = (integration_points.back()).Coordinate(1);
                
                if (xi_0 == -1.0)
                {
                    if (xi_1 < 1.0)
                    {
                        compute = 1;
                    }
                }
                else if (xi_0 > -1.0)
                {
                    if (xi_1 == 1.0)
                    {
                        compute = 2;
                    }
                    else
                    {
                        compute = 3;
                    }
                }
            }
        }
    }

    return compute;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
bounded_matrix<double, TDim, TDim> MortarContactCondition<TDim,TNumNodes,TDoubleLM>::GetRotationMatrixSlave(unsigned int i_slave)
{
    bounded_matrix<double, TDim, TDim> RotationMatrix;
   
    const array_1d<double, 3> normal      = GetGeometry()[i_slave].GetValue(NORMAL);
    const array_1d<double, 3> tangent_xi  = GetGeometry()[i_slave].GetValue(TANGENT_XI);
    const array_1d<double, 3> tangent_eta = GetGeometry()[i_slave].GetValue(TANGENT_ETA);
    
    for (unsigned int i_dim = 0; i_dim < TDim; i_dim++)
    {
        RotationMatrix(0, i_dim) = normal[i_dim];
        RotationMatrix(1, i_dim) = tangent_xi[i_dim];
        if (TDim == 3)
        {
            RotationMatrix(2, i_dim) = tangent_eta[i_dim];
        }
    }
    
    return RotationMatrix;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes , bool TDoubleLM >
void MortarContactCondition<TDim,TNumNodes,TDoubleLM>::ComputeSelectiveIntegrationMethod(const unsigned int rPairIndex)
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
void MortarContactCondition<2, 2, false>::InitializeIntegrationMethod()
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
void MortarContactCondition<2, 3, false>::InitializeIntegrationMethod()
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
void MortarContactCondition<3, 3, false>::InitializeIntegrationMethod()
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
void MortarContactCondition<3, 4, false>::InitializeIntegrationMethod()
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

/***********************************************************************************/ // TODO: Look for an alternative way to do this
/***********************************************************************************/

template< >
void MortarContactCondition<2, 2, true>::InitializeIntegrationMethod()
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
void MortarContactCondition<2, 3, true>::InitializeIntegrationMethod()
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
void MortarContactCondition<3, 3, true>::InitializeIntegrationMethod()
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
void MortarContactCondition<3, 4, true>::InitializeIntegrationMethod()
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

template class MortarContactCondition<2, 2, false>;
template class MortarContactCondition<2, 3, false>;
template class MortarContactCondition<3, 3, false>;
template class MortarContactCondition<3, 4, false>;
template class MortarContactCondition<2, 2, true>;
template class MortarContactCondition<2, 3, true>;
template class MortarContactCondition<3, 3, true>;
template class MortarContactCondition<3, 4, true>;

} // Namespace Kratos
