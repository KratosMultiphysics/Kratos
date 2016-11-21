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
// Avoiding using the macro since this has a template parameter. If there was no template plase use the KRATOS_CREATE_LOCAL_FLAG macro
template< unsigned int TDim, unsigned int TNumNodes >
const Kratos::Flags MortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR(Kratos::Flags::Create(0));
template< unsigned int TDim, unsigned int TNumNodes >
const Kratos::Flags MortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX(Kratos::Flags::Create(1));
template< unsigned int TDim, unsigned int TNumNodes >
const Kratos::Flags MortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS(Kratos::Flags::Create(2));
template< unsigned int TDim, unsigned int TNumNodes >
const Kratos::Flags MortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS(Kratos::Flags::Create(3));

/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer MortarContactCondition<TDim,TNumNodes>::Create( 
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties ) const
{
    return boost::make_shared< MortarContactCondition<TDim,TNumNodes> >( NewId, this->GetGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer MortarContactCondition<TDim,TNumNodes>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return boost::make_shared< MortarContactCondition<TDim,TNumNodes> >( NewId, pGeom, pProperties );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
MortarContactCondition<TDim,TNumNodes>::~MortarContactCondition( )
{
}

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void MortarContactCondition<TDim,TNumNodes>::Initialize( ) 
{
    KRATOS_TRY;
    
    InitializeIntegrationMethod();
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void MortarContactCondition<TDim,TNumNodes>::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
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

template< unsigned int TDim, unsigned int TNumNodes >
void MortarContactCondition<TDim,TNumNodes>::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    // NOTE: Add things if necessary

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void MortarContactCondition<TDim,TNumNodes>::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // Create and initialize condition variables:
    GeneralVariables Variables;
    
    // Initialize the current contact data
    ContactData rContactData;
    
    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = mUseManualColocationIntegration ?
                                                                         mColocationIntegration.IntegrationPoints( ) :
                                                                         GetGeometry( ).IntegrationPoints( mThisIntegrationMethod );
    this->InitializeContactData(rContactData, rCurrentProcessInfo);
    
    std::vector<contact_container> * all_containers = this->GetValue(CONTACT_CONTAINERS);
    
    for (unsigned int PairIndex = 0; PairIndex < mThisMasterElements.size( ); ++PairIndex)
    {   
        // Initialize general variables for the current master element
        this->InitializeGeneralVariables( Variables, rCurrentProcessInfo, PairIndex );
        
        // Update the contact data
        this->UpdateContactData(rContactData, PairIndex);
         
        // Master segment info
        const GeometryType& current_master_element = Variables.GetMasterElement( );
        
        array_1d<double, TNumNodes> aux_int_gap      = ZeroVector(TNumNodes);
        array_1d<double, TNumNodes> aux_int_slip     = ZeroVector(TNumNodes);
        array_1d<double, TNumNodes> aux_int_friction = ZeroVector(TNumNodes);
        
        double total_weight = 0.0;
        
        // Integrating the LHS and RHS
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            // Calculate the kinematic variables
            const bool inside = this->CalculateKinematics( Variables, PointNumber, PairIndex, integration_points );
        
            // Calculate the gap in the integration node and check tolerance
            const double integration_point_gap = inner_prod(rContactData.Gaps, Variables.N_Slave);
            
            if (inside == true)
            {
                const double IntegrationWeight = integration_points[PointNumber].Weight();
                     
                total_weight += IntegrationWeight;
                
                // The slip of th GP
                double integration_point_slip;
                // The tangent LM augmented
                const double augmented_tangent_lm = this->AugmentedTangentLM(Variables, rContactData, current_master_element, integration_point_slip);
                
                for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
                {
                    const double int_aux = IntegrationWeight * Variables.DetJSlave * Variables.Phi_LagrangeMultipliers[iNode];
                    aux_int_gap[iNode]      +=  integration_point_gap  * int_aux;
                    aux_int_slip[iNode]     +=  integration_point_slip * int_aux;
                    aux_int_friction[iNode] +=  Variables.mu           * int_aux;
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

template< unsigned int TDim, unsigned int TNumNodes >
IntegrationMethod MortarContactCondition<TDim,TNumNodes>::GetIntegrationMethod()
{   
    return mThisIntegrationMethod;
}

/***********************************************************************************/
/***********************************************************************************/

template< >
const array_1d<double, 2> MortarContactCondition<2, 2>::LagrangeMultiplierShapeFunctionValue(
    const double xi_local,
    const double eta_local
    )
{
    // NOTE: For more information look at the thesis of Popp page 93/236
   
    array_1d<double, 2> Phi;
            
    Phi[0] = ( 0.5 * ( 1.0 - 3.0 * xi_local ) );
    Phi[1] = ( 0.5 * ( 1.0 + 3.0 * xi_local ) );
  
    return Phi;
}

/***********************************************************************************/
/***********************************************************************************/

template< >
const array_1d<double, 3> MortarContactCondition<2, 3>::LagrangeMultiplierShapeFunctionValue(
    const double xi_local,
    const double eta_local
    )
{
    // NOTE: For more information look at the thesis of Popp page 93/236
   
    array_1d<double, 3> Phi;
           
    array_1d<double,3> aux_coordinates = ZeroVector(3);
    aux_coordinates[0] = xi_local;
    VectorType Ncontainer;
    Ncontainer = GetGeometry().ShapeFunctionsValues(Ncontainer, aux_coordinates);

    Phi[0] = Ncontainer(0) -  3.0/4.0 * Ncontainer(2) + 0.5;
    Phi[1] = Ncontainer(1) -  3.0/4.0 * Ncontainer(2) + 0.5;
    Phi[2] = 5.0/2.0 * Ncontainer(2) - 1.0;
       
    return Phi;
}

/***********************************************************************************/
/***********************************************************************************/

template< >
const array_1d<double, 3> MortarContactCondition<3, 3>::LagrangeMultiplierShapeFunctionValue(
    const double xi_local,
    const double eta_local
    )
{
    // NOTE: For more information look at the thesis of Popp page 93/236
   
    array_1d<double, 3> Phi;
            
    Phi[0] = 3.0 - 4.0 * xi_local - 4.0 * eta_local;
    Phi[1] = 4.0 * xi_local  - 1.0;
    Phi[2] = 4.0 * eta_local - 1.0;    
        
    return Phi;
}


/***********************************************************************************/
/***********************************************************************************/

template< >
const array_1d<double, 4> MortarContactCondition<3, 4>::LagrangeMultiplierShapeFunctionValue(
    const double xi_local,
    const double eta_local
    )
{
    // NOTE: For more information look at the thesis of Popp page 93/236
   
    array_1d<double, 4> Phi;
            
    array_1d<double,3> aux_coordinates = ZeroVector(3);
    aux_coordinates[0] =  xi_local;
    aux_coordinates[1] = eta_local;
    VectorType Ncontainer;
    Ncontainer = GetGeometry().ShapeFunctionsValues(Ncontainer, aux_coordinates);

    Phi[0] =   4.0 * Ncontainer(0) - 2.0 * Ncontainer(1) + 1.0 * Ncontainer(2) - 2.0 * Ncontainer(3);
    Phi[1] = - 2.0 * Ncontainer(0) + 4.0 * Ncontainer(1) - 2.0 * Ncontainer(2) + 1.0 * Ncontainer(3);
    Phi[2] =   1.0 * Ncontainer(0) - 2.0 * Ncontainer(1) + 4.0 * Ncontainer(2) - 2.0 * Ncontainer(3);
    Phi[3] = - 2.0 * Ncontainer(0) + 1.0 * Ncontainer(1) - 2.0 * Ncontainer(2) + 4.0 * Ncontainer(3);

    return Phi;
}

/***********************************************************************************/
/***********************************************************************************/

template< >
const bounded_matrix<double, 2, 1> MortarContactCondition<2, 2>::LagrangeMultiplierShapeFunctionLocalGradient( 
    const double xi_local, 
    const double eta_local 
    )
{
    bounded_matrix<double, 2, 1> DPhi_De;
            
    DPhi_De( 0, 0 ) = - 3.0 / 2.0;
    DPhi_De( 1, 0 ) = + 3.0 / 2.0;
   
    return DPhi_De;
}

/***********************************************************************************/
/***********************************************************************************/

template< >
const bounded_matrix<double, 3, 1> MortarContactCondition<2, 3>::LagrangeMultiplierShapeFunctionLocalGradient( 
    const double xi_local, 
    const double eta_local 
    )
{
    bounded_matrix<double, 3, 1> DPhi_De;
            
    array_1d<double,3> aux_coordinates = ZeroVector(3);
    aux_coordinates[0] = xi_local;
    MatrixType DNcontainer;
    DNcontainer = GetGeometry().ShapeFunctionsLocalGradients( DNcontainer , aux_coordinates );
    
    DPhi_De( 0, 0 ) = DNcontainer(0, 0) -  3.0/4.0 * DNcontainer(2, 0);
    DPhi_De( 1, 0 ) = DNcontainer(1, 0) -  3.0/4.0 * DNcontainer(2, 0);
    DPhi_De( 2, 0 ) = 5.0/2.0 * DNcontainer(2, 0);
    
    return DPhi_De;
}

/***********************************************************************************/
/***********************************************************************************/

template< >
const bounded_matrix<double, 3, 2> MortarContactCondition<3 ,3>::LagrangeMultiplierShapeFunctionLocalGradient( 
    const double xi_local, 
    const double eta_local 
    )
{
    bounded_matrix<double, 3, 2> DPhi_De;
            
    DPhi_De( 0, 0 ) = - 4.0;
    DPhi_De( 1, 0 ) =   4.0;
    DPhi_De( 2, 0 ) =   4.0;
    
    DPhi_De( 0, 1 ) = - 4.0;
    DPhi_De( 1, 1 ) =   0.0;
    DPhi_De( 2, 1 ) =   4.0;

    return DPhi_De;
}

/***********************************************************************************/
/***********************************************************************************/

template< >
const bounded_matrix<double, 4, 2> MortarContactCondition<3, 4>::LagrangeMultiplierShapeFunctionLocalGradient( 
    const double xi_local, 
    const double eta_local 
    )
{
    bounded_matrix<double, 4, 2> DPhi_De;
    
    array_1d<double,3> aux_coordinates = ZeroVector(3);
    aux_coordinates[0] =  xi_local;
    aux_coordinates[1] = eta_local;
    MatrixType DNcontainer;
    DNcontainer = GetGeometry().ShapeFunctionsLocalGradients(DNcontainer, aux_coordinates);

    DPhi_De( 0, 0 ) =   4.0 * DNcontainer( 0, 0 ) - 2.0 * DNcontainer( 1, 0 ) + 1.0 * DNcontainer( 2, 0 ) - 2.0 * DNcontainer( 3, 0 );
    DPhi_De( 1, 0 ) = - 2.0 * DNcontainer( 0, 0 ) + 4.0 * DNcontainer( 1, 0 ) - 2.0 * DNcontainer( 2, 0 ) + 1.0 * DNcontainer( 3, 0 );
    DPhi_De( 2, 0 ) =   1.0 * DNcontainer( 0, 0 ) - 2.0 * DNcontainer( 1, 0 ) + 4.0 * DNcontainer( 2, 0 ) - 2.0 * DNcontainer( 3, 0 );
    DPhi_De( 3, 0 ) = - 2.0 * DNcontainer( 0, 0 ) + 1.0 * DNcontainer( 1, 0 ) - 2.0 * DNcontainer( 2, 0 ) + 4.0 * DNcontainer( 3, 0 );
    DPhi_De( 0, 1 ) =   4.0 * DNcontainer( 0, 1 ) - 2.0 * DNcontainer( 1, 1 ) + 1.0 * DNcontainer( 2, 1 ) - 2.0 * DNcontainer( 3, 1 );
    DPhi_De( 1, 1 ) = - 2.0 * DNcontainer( 0, 1 ) + 4.0 * DNcontainer( 1, 1 ) - 2.0 * DNcontainer( 2, 1 ) + 1.0 * DNcontainer( 3, 1 );
    DPhi_De( 2, 1 ) =   1.0 * DNcontainer( 0, 1 ) - 2.0 * DNcontainer( 1, 1 ) + 4.0 * DNcontainer( 2, 1 ) - 2.0 * DNcontainer( 3, 1 );
    DPhi_De( 3, 1 ) = - 2.0 * DNcontainer( 0, 1 ) + 1.0 * DNcontainer( 1, 1 ) - 2.0 * DNcontainer( 2, 1 ) + 4.0 * DNcontainer( 3, 1 );

    return DPhi_De;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void MortarContactCondition<TDim,TNumNodes>::CalculateLocalSystem( 
    std::vector<MatrixType>& rLeftHandSideMatrices,
    const std::vector<Variable<MatrixType> >& rLHSVariables,
    std::vector<VectorType>& rRightHandSideVectors,
    const std::vector<Variable<VectorType> >& rRHSVariables,
    ProcessInfo& rCurrentProcessInfo 
    )
{    
    // Calculates the size of the system
    constexpr unsigned int MatrixSize = TDim * (2 * TNumNodes + TNumNodes);
    
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set(MortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS, true);
    LocalSystem.CalculationFlags.Set(MortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS, true);

    //Initialize sizes for the system components:
    if ( rLHSVariables.size( ) != rLeftHandSideMatrices.size( ) )
    {
        rLeftHandSideMatrices.resize( rLHSVariables.size( ) );
    }

    if ( rRHSVariables.size( ) != rRightHandSideVectors.size( ) )
    {
        rRightHandSideVectors.resize( rRHSVariables.size( ) );
    }

    LocalSystem.CalculationFlags.Set(MortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX, true);
    for ( unsigned int i = 0; i < rLeftHandSideMatrices.size( ); i++ )
    {
        // Note: rRightHandSideVectors.size() > 0
        this->InitializeSystemMatrices<MatrixSize>( rLeftHandSideMatrices[i], rRightHandSideVectors[0],LocalSystem.CalculationFlags );
    }

    LocalSystem.CalculationFlags.Set( MortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR, true );
    LocalSystem.CalculationFlags.Set( MortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX, false ); // Temporarily only
    for ( unsigned int i = 0; i < rRightHandSideVectors.size( ); i++ )
    {
        // Note: rLeftHandSideMatrices.size() > 0
        this->InitializeSystemMatrices<MatrixSize>( rLeftHandSideMatrices[0], rRightHandSideVectors[i], LocalSystem.CalculationFlags  );
    }
    LocalSystem.CalculationFlags.Set( MortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX, true ); // Reactivated again

    // Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrices( rLeftHandSideMatrices );
    LocalSystem.SetRightHandSideVectors( rRightHandSideVectors );

    LocalSystem.SetLeftHandSideVariables( rLHSVariables );
    LocalSystem.SetRightHandSideVariables( rRHSVariables );

    // Calculate condition system
    this->CalculateConditionSystem<MatrixSize>( LocalSystem, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void MortarContactCondition<TDim,TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;
    
    // Calculates the size of the system
    constexpr unsigned int MatrixSize = TDim * (2 * TNumNodes + TNumNodes);
    
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX, true );
    LocalSystem.CalculationFlags.Set( MortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR, true );

    // Initialize sizes for the system components:
    this->InitializeSystemMatrices<MatrixSize>( rLeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );
    
    // Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix( rLeftHandSideMatrix );
    LocalSystem.SetRightHandSideVector( rRightHandSideVector );

    // Calculate condition system
    this->CalculateConditionSystem<MatrixSize>( LocalSystem, rCurrentProcessInfo );

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void MortarContactCondition<TDim,TNumNodes>::CalculateLeftHandSide( 
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Calculates the size of the system
    constexpr unsigned int MatrixSize = TDim * (2 * TNumNodes + TNumNodes);
    
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX, true );

    VectorType RightHandSideVector = Vector( );

    // Initialize sizes for the system components:
    this->InitializeSystemMatrices<MatrixSize>( rLeftHandSideMatrix, RightHandSideVector, LocalSystem.CalculationFlags );

    // Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix( rLeftHandSideMatrix );
    LocalSystem.SetRightHandSideVector( RightHandSideVector );

    // Calculate condition system
    this->CalculateConditionSystem<MatrixSize>( LocalSystem, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void MortarContactCondition<TDim,TNumNodes>::CalculateLeftHandSide( 
    std::vector< MatrixType >& rLeftHandSideMatrices,
    const std::vector< Variable< MatrixType > >& rLHSVariables,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Calculates the size of the system
    constexpr unsigned int MatrixSize = TDim * (2 * TNumNodes + TNumNodes);
    
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX, true );
    LocalSystem.CalculationFlags.Set( MortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS, true );

    VectorType RightHandSideVector = Vector( );

    // Initialize size for the system components
    for( unsigned int i = 0; i < rLeftHandSideMatrices.size( ); i++ )
    {
        this->InitializeSystemMatrices<MatrixSize>( rLeftHandSideMatrices[i], RightHandSideVector, LocalSystem.CalculationFlags );
    }

    LocalSystem.SetLeftHandSideMatrices( rLeftHandSideMatrices );
    LocalSystem.SetRightHandSideVector( RightHandSideVector );

    // Calculate condition system
    this->CalculateConditionSystem<MatrixSize>( LocalSystem, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void MortarContactCondition<TDim,TNumNodes>::CalculateRightHandSide( 
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Calculates the size of the system
    constexpr unsigned int MatrixSize = TDim * (2 * TNumNodes + TNumNodes);
    
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR, true);

    MatrixType LeftHandSideMatrix = Matrix( );

    // Initialize size for the system components
    this->InitializeSystemMatrices<MatrixSize>( LeftHandSideMatrix, rRightHandSideVector,LocalSystem.CalculationFlags);

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix( LeftHandSideMatrix );
    LocalSystem.SetRightHandSideVector( rRightHandSideVector );

    //Calculate condition system
    this->CalculateConditionSystem<MatrixSize>( LocalSystem, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void MortarContactCondition<TDim,TNumNodes>::CalculateRightHandSide( 
    std::vector< VectorType >& rRightHandSideVectors,
    const std::vector< Variable< VectorType > >& rRHSVariables,
    ProcessInfo& rCurrentProcessInfo )
{
    // Calculates the size of the system
    constexpr unsigned int MatrixSize = TDim * (2 * TNumNodes + TNumNodes);
        
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR, true );
    LocalSystem.CalculationFlags.Set( MortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS, true );

    MatrixType LeftHandSideMatrix = Matrix( );

    // Initialize size for the system components
    for( unsigned int i = 0; i < rRightHandSideVectors.size(); i++ )
    {
        this->InitializeSystemMatrices<MatrixSize>( LeftHandSideMatrix, rRightHandSideVectors[i], LocalSystem.CalculationFlags );
    }

    // Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix( LeftHandSideMatrix );
    LocalSystem.SetRightHandSideVectors( rRightHandSideVectors );

    // Calculate condition system
    this->CalculateConditionSystem<MatrixSize>( LocalSystem, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
template< unsigned int MatrixSize >
void MortarContactCondition<TDim,TNumNodes>::InitializeSystemMatrices( 
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    Flags& rCalculationFlags 
    )
{
    const unsigned int condition_size = this->CalculateConditionSize<MatrixSize>( );
    
    // Resizing as needed the LHS
    if ( rCalculationFlags.Is( MortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX ) ) // Calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != condition_size )
        {
            rLeftHandSideMatrix.resize( condition_size, condition_size, false );
        }
        noalias( rLeftHandSideMatrix ) = ZeroMatrix( condition_size, condition_size ); // Resetting LHS
    }

    // Resizing as needed the RHS
    if ( rCalculationFlags.Is( MortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR ) ) // Calculation of the matrix is required
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

template< unsigned int TDim, unsigned int TNumNodes >
void MortarContactCondition<TDim,TNumNodes>::CalculateMassMatrix( 
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

template< unsigned int TDim, unsigned int TNumNodes >
void MortarContactCondition<TDim,TNumNodes>::CalculateDampingMatrix( 
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

template< unsigned int TDim, unsigned int TNumNodes >
template< unsigned int MatrixSize >
const unsigned int MortarContactCondition<TDim,TNumNodes>::CalculateConditionSize( )
{
    const unsigned int condition_size = mThisMasterElements.size( ) * MatrixSize; // NOTE: Assuming same number of nodes for the master
    
    return condition_size;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
template< unsigned int MatrixSize>
void MortarContactCondition<TDim, TNumNodes>::CalculateConditionSystem( 
    LocalSystemComponents& rLocalSystem,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
    
    // Create and initialize condition variables:#pragma omp critical
    GeneralVariables Variables;
    
    // Initialize the current contact data
    ContactData rContactData;

    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = mUseManualColocationIntegration ?
                                                                         mColocationIntegration.IntegrationPoints( ) :
                                                                         GetGeometry( ).IntegrationPoints( mThisIntegrationMethod );
                                                                                  
    this->InitializeContactData(rContactData, rCurrentProcessInfo);
    
    for (unsigned int PairIndex = 0; PairIndex < mThisMasterElements.size( ); ++PairIndex)
    {   
        // Initialize general variables for the current master element
        this->InitializeGeneralVariables( Variables, rCurrentProcessInfo, PairIndex );
        
        // Update the contact data
        this->UpdateContactData(rContactData, PairIndex);
         
        // Master segment info
        const GeometryType& current_master_element = Variables.GetMasterElement( );

        bounded_matrix<double, MatrixSize, MatrixSize> LHS_contact_pair = ZeroMatrix(MatrixSize, MatrixSize);
        array_1d<double, MatrixSize> RHS_contact_pair = ZeroVector(MatrixSize);
        
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
                total_weight += IntegrationWeight; // NOTE: I keep this, because with ALM if it is inactive there is a penalty term in the LHS and the RHS
               
                // The normal LM augmented
                double augmented_normal_lm = this->AugmentedNormalLM(Variables, rContactData, integration_point_gap);
                  
                // The slip of th GP
                double integration_point_slip;
                // The tangent LM augmented
                const double augmented_tangent_lm = this->AugmentedTangentLM(Variables, rContactData, current_master_element, integration_point_slip);
                
                if (rCurrentProcessInfo[TIME_STEPS] == 1) // NOTE: To avoid problems the first iteration
                {
                    augmented_normal_lm -= 1.0e-16;
                }
                
                // Calculation of the matrix is required
                if ( rLocalSystem.CalculationFlags.Is( MortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX ) ||
                        rLocalSystem.CalculationFlags.Is( MortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
                {
                    this->CalculateLocalLHS<MatrixSize>( LHS_contact_pair, Variables, rContactData, IntegrationWeight, augmented_normal_lm, augmented_tangent_lm, integration_point_gap,  integration_point_slip );
                }
                
                // Calculation of the vector is required
                if ( rLocalSystem.CalculationFlags.Is( MortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR ) ||
                        rLocalSystem.CalculationFlags.Is( MortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
                {
                    this->CalculateLocalRHS<MatrixSize>( RHS_contact_pair, Variables, rContactData, IntegrationWeight, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip );
                }
            }
        }
        
        // We can consider the pair if at least one of the collocation point is inside 
        if (total_weight > 0.0)
        {
            // Assemble of the matrix is required
            if ( rLocalSystem.CalculationFlags.Is( MortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX ) ||
                    rLocalSystem.CalculationFlags.Is( MortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
            {
                // Contributions to stiffness matrix calculated on the reference config
                this->CalculateAndAddLHS<MatrixSize>( rLocalSystem, LHS_contact_pair, PairIndex, current_master_element );
            }

            // Assemble of the vector is required
            if ( rLocalSystem.CalculationFlags.Is( MortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR ) ||
                    rLocalSystem.CalculationFlags.Is( MortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
            {
                // Contribution to previous step contact force and residuals vector
                this->CalculateAndAddRHS<MatrixSize>( rLocalSystem, RHS_contact_pair, PairIndex, current_master_element );
            }
            
//             std::cout << "--------------------------------------------------" << std::endl;
//             KRATOS_WATCH(this->Id());
//             KRATOS_WATCH(PairIndex);
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

template< unsigned int TDim, unsigned int TNumNodes >
void MortarContactCondition<TDim,TNumNodes>::InitializeGeneralVariables(
    GeneralVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo,
    const unsigned int& rMasterElementIndex
    )
{
    // Master segment info
    GeometryType& current_master_element = mThisMasterElements[rMasterElementIndex]->GetGeometry();
    
    // Slave element info
    rVariables.Initialize();

    rVariables.SetMasterElement( current_master_element );
    rVariables.SetMasterElementIndex( rMasterElementIndex );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void MortarContactCondition<TDim,TNumNodes>::InitializeContactData(
    ContactData& rContactData,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // Slave element info
    rContactData.Initialize(GetGeometry(), TNumNodes, TDim );
    
    /* Set Delta time */
    rContactData.Dt = rCurrentProcessInfo[DELTA_TIME];
    
    /* LM */
    rContactData.LagrangeMultipliers = ContactUtilities::GetVariableMatrix(GetGeometry(), VECTOR_LAGRANGE_MULTIPLIER, 0); // NOTE: Valgrind says there is a problem here
    
    /* NORMALS AND TANGENTS */ // TODO: To interpolate it is necessary to consider a smooth function
    const array_1d<double, 3> normal      = this->GetValue(NORMAL);
    const array_1d<double, 3> tangent_xi  = this->GetValue(TANGENT_XI);
    const array_1d<double, 3> tangent_eta = this->GetValue(TANGENT_ETA);
    for (unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++)
    {
        for (unsigned int i_dim = 0; i_dim < TDim; i_dim++)
        {
            rContactData.NormalsSlave(i_slave, i_dim) = normal[i_dim];
            rContactData.Tangent1Slave(i_slave, i_dim) = tangent_xi[i_dim];
            rContactData.Tangent2Slave(i_slave, i_dim) = tangent_eta[i_dim];
        }
    }
    
//     rContactData.NormalsSlave = ContactUtilities::GetVariableMatrix(GetGeometry(),  NORMAL); 
//     rContactData.Tangent1Slave = ContactUtilities::GetVariableMatrix(GetGeometry(), TANGENT_XI); 
//     if (TDim == 3)
//     {
//         rContactData.Tangent2Slave = ContactUtilities::GetVariableMatrix(GetGeometry(), TANGENT_ETA); 
//     }
//     
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

template< unsigned int TDim, unsigned int TNumNodes >
void MortarContactCondition<TDim,TNumNodes>::UpdateContactData(
    ContactData& rContactData,
    const unsigned int& rMasterElementIndex
    )
{    
    // Master segment info
    GeometryType& current_master_element = mThisMasterElements[rMasterElementIndex]->GetGeometry();
    
    // Slave element info
    rContactData.UpdateMasterPair(mThisMasterElements[rMasterElementIndex], TNumNodes, TDim );
//     rContactData.UpdateMasterPair(current_master_element, TNumNodes, TDim );
    
    /* NORMALS AND GAPS */
    for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
    {
        array_1d<double,3> normal = this->GetValue(NORMAL);
//         array_1d<double,3> normal = GetGeometry()[iNode].GetValue(NORMAL);
        
        PointType projected_global;
        ContactUtilities::ProjectDirection(current_master_element, GetGeometry()[iNode], projected_global, rContactData.Gaps(iNode), normal ); // NOTE: This is not the CPP, so the solution can be wrong
    }
}

/*********************************COMPUTE KINEMATICS*********************************/
/************************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
bool MortarContactCondition<TDim,TNumNodes>::CalculateKinematics( 
    GeneralVariables& rVariables,
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
    if (GetProperties().Has(FRICTION_COEFFICIENT) == true) // NOTE: In function of the friction law (TODO: Remove this!!!)
    {
        rVariables.mu = GetProperties().GetValue(FRICTION_COEFFICIENT);
    }        
    
    return inside;
}
 
/***********************************************************************************/
/*************** METHODS TO CALCULATE THE CONTACT CONDITION MATRICES ***************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
bool MortarContactCondition<TDim,TNumNodes>::MasterShapeFunctionValue(
    GeneralVariables& rVariables,
    const PointType& local_point 
    )
{    
    GeometryType& master_seg = rVariables.GetMasterElement( );
    rVariables.N_Master.clear();
    rVariables.DN_De_Master.clear();

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

template< unsigned int TDim, unsigned int TNumNodes >
template< unsigned int MatrixSize >
void MortarContactCondition<TDim,TNumNodes>::CalculateAndAddLHS(
    LocalSystemComponents& rLocalSystem,
    bounded_matrix<double, MatrixSize, MatrixSize>& LHS_contact_pair, 
    const unsigned int rPairIndex,
    const GeometryType& current_master_element
    )
{
    /* DEFINITIONS */
    const unsigned int number_of_total_nodes = TNumNodes + TNumNodes; // NOTE: Assuming same number of nodes for the master
        
    unsigned int index_1 = 0;
    for (unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave )
    {
        index_1 = (number_of_total_nodes + i_slave) * TDim;
                    
        // We impose a 0 zero LM in the inactive nodes
        if (GetGeometry( )[i_slave].Is(ACTIVE) == false)
        {
//             subrange(LHS_contact_pair, 0, number_of_total_nodes * TDim, index_1, index_1 + TDim) = ZeroMatrix(number_of_total_nodes * TDim, TDim); // TODO: Consider or not?
            subrange(LHS_contact_pair, index_1, index_1 + TDim, 0, number_of_total_nodes * TDim) = ZeroMatrix(TDim, number_of_total_nodes * TDim);
            subrange(LHS_contact_pair, index_1, index_1 + TDim, index_1, index_1 + TDim)         = IdentityMatrix(TDim, TDim);
        }
        // And the tangent direction (sliding)
        else
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
    
    if ( rLocalSystem.CalculationFlags.Is( MortarContactCondition<TDim,TNumNodes>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
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
                this->AssembleContactPairLHSToConditionSystem<MatrixSize>(LHS_contact_pair, rLeftHandSideMatrix, rPairIndex);
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
        this->AssembleContactPairLHSToConditionSystem<MatrixSize>(LHS_contact_pair, rLeftHandSideMatrix, rPairIndex);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
template< unsigned int MatrixSize>
void MortarContactCondition<TDim,TNumNodes>::AssembleContactPairLHSToConditionSystem(
    bounded_matrix<double, MatrixSize, MatrixSize>& rPairLHS,
    MatrixType& rConditionLHS,
    const unsigned int rPairIndex
    )
{
    // Find location of the pair's master DOFs in ConditionRHS
    const unsigned int index_begin = rPairIndex * MatrixSize;
    const unsigned int index_end  = index_begin + MatrixSize;
    
    subrange( rConditionLHS, index_begin, index_end, index_begin, index_end) += rPairLHS;
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MortarContactCondition<2, 2>::CalculateLocalLHS<12>(
    bounded_matrix<double, 12, 12>& rPairLHS,
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
    const Vector N1           = rVariables.N_Slave;
    const Vector N2           = rVariables.N_Master;
    const Vector Phi          = rVariables.Phi_LagrangeMultipliers;
    const double detJ         = rVariables.DetJSlave; 

    if (augmented_normal_lm < 0.0)
    {                                
        // Contact active
        rPairLHS += rIntegrationWeight * Contact2D2N2N::ComputeGaussPointActiveLHS( N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
        
        if (std::abs(augmented_tangent_lm) - rVariables.mu * augmented_normal_lm >= 0.0)
        {
            // Slip
            rPairLHS += rIntegrationWeight * Contact2D2N2N::ComputeGaussPointSlipLHS( N1, N2, Phi, detJ, rVariables.mu, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
        }
        else
        {
            // Stick
            rPairLHS += rIntegrationWeight * Contact2D2N2N::ComputeGaussPointStickLHS( N1, N2, Phi, detJ, rVariables.mu, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
        }
    }
//     else
//     {        
//         // Contact inactive
//         rPairLHS += rIntegrationWeight * Contact2D2N2N::ComputeGaussPointInactiveLHS( N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//     }
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MortarContactCondition<2, 3>::CalculateLocalLHS<18>(
    bounded_matrix<double, 18, 18>& rPairLHS,
    GeneralVariables& rVariables,
    const ContactData& rContactData,
    const double& rIntegrationWeight,
    const double& augmented_normal_lm,
    const double& augmented_tangent_lm,
    const double& integration_point_gap,
    const double& integration_point_slip
    )
{
//     /* DEFINITIONS */
//     const Vector N1           = rVariables.N_Slave;
//     const Vector N2           = rVariables.N_Master;
//     const Vector Phi          = rVariables.Phi_LagrangeMultipliers;
//     const double detJ         = rVariables.DetJSlave; 
// 
    // TODO: Finish this!!!!
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MortarContactCondition<3, 3>::CalculateLocalLHS<27>(
    bounded_matrix<double, 27, 27>& rPairLHS,
    GeneralVariables& rVariables,
    const ContactData& rContactData,
    const double& rIntegrationWeight,
    const double& augmented_normal_lm,
    const double& augmented_tangent_lm,
    const double& integration_point_gap,
    const double& integration_point_slip
    )
{
//     /* DEFINITIONS */
//     const Vector N1           = rVariables.N_Slave;
//     const Vector N2           = rVariables.N_Master;
//     const Vector Phi          = rVariables.Phi_LagrangeMultipliers;
//     const double detJ         = rVariables.DetJSlave; 
// 
//     if (augmented_normal_lm < 0.0)
//     {                                
//         // Contact active
//         rPairLHS += rIntegrationWeight * Contact3D3N3N::ComputeGaussPointActiveLHS( N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//         
//         if (std::abs(augmented_tangent_lm) - rVariables.mu * augmented_normal_lm >= 0.0)
//         {
//             // Slip
//             rPairLHS += rIntegrationWeight * Contact3D3N3N::ComputeGaussPointSlipLHS( N1, N2, Phi, detJ, rVariables.mu, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//         }
//         else
//         {
//             // Stick
//             rPairLHS += rIntegrationWeight * Contact3D3N3N::ComputeGaussPointStickLHS( N1, N2, Phi, detJ, rVariables.mu, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//         }
//     }
// //     else
// //     {        
// //         // Contact inactive
// //         rPairLHS += rIntegrationWeight * Contact3D3N3N::ComputeGaussPointInactiveLHS( N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
// //     }
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MortarContactCondition<3, 4>::CalculateLocalLHS<36>(
    bounded_matrix<double, 36, 36>& rPairLHS,
    GeneralVariables& rVariables,
    const ContactData& rContactData,
    const double& rIntegrationWeight,
    const double& augmented_normal_lm,
    const double& augmented_tangent_lm,
    const double& integration_point_gap,
    const double& integration_point_slip
    )
{
//     /* DEFINITIONS */
//     const Vector N1           = rVariables.N_Slave;
//     const Vector N2           = rVariables.N_Master;
//     const Vector Phi          = rVariables.Phi_LagrangeMultipliers;
//     const Matrix DPhi         = rVariables.DPhi_De_LagrangeMultipliers;
//     const double detJ         = rVariables.DetJSlave; 
// 
//     if (augmented_normal_lm < 0.0)
//     {                                
//         // Contact active
//         rPairLHS += rIntegrationWeight * Contact3D4N4N::ComputeGaussPointActiveLHS( N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//         
//         if (std::abs(augmented_tangent_lm) - rVariables.mu * augmented_normal_lm >= 0.0)
//         {
//             // Slip
//             rPairLHS += rIntegrationWeight * Contact3D4N4N::ComputeGaussPointSlipLHS( N1, N2, Phi, detJ, rVariables.mu, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//         }
//         else
//         {
//             // Stick
//             rPairLHS += rIntegrationWeight * Contact3D4N4N::ComputeGaussPointStickLHS( N1, N2, Phi, detJ, rVariables.mu, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//         }
//     }
// //     else
// //     {        
// //         // Contact inactive
// //         rPairLHS += rIntegrationWeight * Contact3D4N4N::ComputeGaussPointInactiveLHS( N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
// //     }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
template< unsigned int MatrixSize >
void MortarContactCondition<TDim,TNumNodes>::CalculateAndAddRHS(
    LocalSystemComponents& rLocalSystem,
    array_1d<double, MatrixSize>& RHS_contact_pair,
    const unsigned int rPairIndex,
    const GeometryType& current_master_element
    )
{   
    /* DEFINITIONS */
    const unsigned int number_of_total_nodes = TNumNodes + TNumNodes; // NOTE: Assuming same number of nodes for the master
        
    unsigned int index = 0;
    
    for (unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave )
    {
        index = (number_of_total_nodes + i_slave) * TDim; 

        if (GetGeometry( )[i_slave].Is(ACTIVE) == false)
        {
            for (unsigned int i = index; i < index + TDim; i++)
            {
                RHS_contact_pair[i] = 0.0;
            }
        }
        else
        {   
            for (unsigned int i = index + 1; i < index + TDim; i++)
            {
                RHS_contact_pair[i] = 0.0;
            }
        }
    }
        
    if ( rLocalSystem.CalculationFlags.Is( MortarContactCondition<TDim,TNumNodes>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
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
                this->AssembleContactPairRHSToConditionSystem<MatrixSize>( RHS_contact_pair, rRightHandSideVector, rPairIndex );
                
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
        this->AssembleContactPairRHSToConditionSystem<MatrixSize>( RHS_contact_pair, rRightHandSideVector, rPairIndex );
    }
}
  
/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
template< unsigned int MatrixSize >
void MortarContactCondition<TDim,TNumNodes>::AssembleContactPairRHSToConditionSystem(
    array_1d<double, MatrixSize>& rPairRHS,
    VectorType& rConditionRHS,
    const unsigned int rPairIndex
    )
{
    // Find location of the pair's master DOFs in ConditionRHS
    const unsigned int index_begin = rPairIndex * MatrixSize;
    const unsigned int index_end  = index_begin + MatrixSize;
    
    // Computing subrange
    subrange( rConditionRHS, index_begin, index_end ) += rPairRHS;
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MortarContactCondition<2, 2>::CalculateLocalRHS<12>(
    array_1d<double,12>& rPairRHS,
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
    const Vector N1           = rVariables.N_Slave;
    const Vector N2           = rVariables.N_Master;
    const Vector Phi          = rVariables.Phi_LagrangeMultipliers;
    const double detJ         = rVariables.DetJSlave;
    
    if (augmented_normal_lm < 0.0)  // TODO: This is a conflict (< or <=¿?¿?¿?)
    {
        // Contact active
        rPairRHS += rIntegrationWeight * Contact2D2N2N::ComputeGaussPointActiveRHS(N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
        
        if (std::abs(augmented_tangent_lm) - rVariables.mu * augmented_normal_lm >= 0.0)
        {
            // Slip
            rPairRHS += rIntegrationWeight * Contact2D2N2N::ComputeGaussPointSlipRHS(N1, N2, Phi, detJ, rVariables.mu, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
        }
        else
        {
            // Stick
            rPairRHS += rIntegrationWeight * Contact2D2N2N::ComputeGaussPointStickRHS(N1, N2, Phi, detJ, rVariables.mu, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
        }
    }
//     else
//     {
//         // Contact inactive
//         rPairRHS += rIntegrationWeight * Contact2D2N2N::ComputeGaussPointInactiveRHS(N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//     }
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MortarContactCondition<2, 3>::CalculateLocalRHS<18>(
    array_1d<double,18>& rPairRHS,
    GeneralVariables& rVariables,
    const ContactData& rContactData,
    const double& rIntegrationWeight,
    const double& augmented_normal_lm,
    const double& augmented_tangent_lm,
    const double& integration_point_gap,
    const double& integration_point_slip
    )
{
//     /* DEFINITIONS */     
//     const Vector N1           = rVariables.N_Slave;
//     const Vector N2           = rVariables.N_Master;
//     const Vector Phi          = rVariables.Phi_LagrangeMultipliers;
//     const double detJ         = rVariables.DetJSlave;
//     
    // TODO: Finish this!!!
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MortarContactCondition<3, 3>::CalculateLocalRHS<27>(
    array_1d<double,27>& rPairRHS,
    GeneralVariables& rVariables,
    const ContactData& rContactData,
    const double& rIntegrationWeight,
    const double& augmented_normal_lm,
    const double& augmented_tangent_lm,
    const double& integration_point_gap,
    const double& integration_point_slip
    )
{
//     /* DEFINITIONS */
//     const Vector N1           = rVariables.N_Slave;
//     const Vector N2           = rVariables.N_Master;
//     const Vector Phi          = rVariables.Phi_LagrangeMultipliers;
//     const double detJ         = rVariables.DetJSlave;
//     
//     if (augmented_normal_lm < 0.0)  // TODO: This is a conflict (< or <=¿?¿?¿?)
//     {
//         // Contact active
//         rPairRHS += rIntegrationWeight * Contact3D3N3N::ComputeGaussPointActiveRHS(N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//         
//         if (std::abs(augmented_tangent_lm) - rVariables.mu * augmented_normal_lm >= 0.0)
//         {
//             // Slip
//             rPairRHS += rIntegrationWeight * Contact3D3N3N::ComputeGaussPointSlipRHS(N1, N2, Phi, detJ, rVariables.mu, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//         }
//         else
//         {
//             // Stick
//             rPairRHS += rIntegrationWeight * Contact3D3N3N::ComputeGaussPointStickRHS(N1, N2, Phi, detJ, rVariables.mu, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//         }
//     }
// //     else
// //     {
// //         // Contact inactive
// //         rPairRHS += rIntegrationWeight * Contact3D3N3N::ComputeGaussPointInactiveRHS(N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
// //     }
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
void MortarContactCondition<3, 4>::CalculateLocalRHS<36>(
    array_1d<double,36>& rPairRHS,
    GeneralVariables& rVariables,
    const ContactData& rContactData,
    const double& rIntegrationWeight,
    const double& augmented_normal_lm,
    const double& augmented_tangent_lm,
    const double& integration_point_gap,
    const double& integration_point_slip
    )
{
//     /* DEFINITIONS */
//     const Vector N1           = rVariables.N_Slave;
//     const Vector N2           = rVariables.N_Master;
//     const Vector Phi          = rVariables.Phi_LagrangeMultipliers;
//     const Matrix DPhi         = rVariables.DPhi_De_LagrangeMultipliers;
//     const double detJ         = rVariables.DetJSlave;
//     
//     if (augmented_normal_lm < 0.0)  // TODO: This is a conflict (< or <=¿?¿?¿?)
//     {
//         // Contact active
//         rPairRHS += rIntegrationWeight * Contact3D4N4N::ComputeGaussPointActiveRHS(N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//         
//         if (std::abs(augmented_tangent_lm) - rVariables.mu * augmented_normal_lm >= 0.0)
//         {
//             // Slip
//             rPairRHS += rIntegrationWeight * Contact3D4N4N::ComputeGaussPointSlipRHS(N1, N2, Phi, detJ, rVariables.mu, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//         }
//         else
//         {
//             // Stick
//             rPairRHS += rIntegrationWeight * Contact3D4N4N::ComputeGaussPointStickRHS(N1, N2, Phi, detJ, rVariables.mu, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
//         }
//     }
// //     else
// //     {
// //         // Contact inactive
// //         rPairRHS += rIntegrationWeight * Contact3D4N4N::ComputeGaussPointInactiveRHS(N1, N2, Phi, detJ, rContactData, augmented_normal_lm, augmented_tangent_lm, integration_point_gap, integration_point_slip);
// //     }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void MortarContactCondition<TDim,TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& CurrentProcessInfo 
    )
{
    KRATOS_TRY;   

    rResult.resize( 0, false );

    const std::vector<contact_container> all_conditions = *( this->GetValue( CONTACT_CONTAINERS ) );

    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    for ( unsigned int i_cond = 0;  i_cond < all_conditions.size( ); ++i_cond )
    {   
        GeometryType& current_master = all_conditions[i_cond].condition->GetGeometry( );

        // Master Nodes Displacement Equation IDs
        for ( unsigned int i_master = 0; i_master < TNumNodes; i_master++ ) // NOTE: Assuming same number of nodes for master and slave
        {
            NodeType& master_node = current_master[i_master];
            rResult.push_back( master_node.GetDof( DISPLACEMENT_X ).EquationId( ) );
            rResult.push_back( master_node.GetDof( DISPLACEMENT_Y ).EquationId( ) );
            if (TDim == 3)
            {
                rResult.push_back( master_node.GetDof( DISPLACEMENT_Z ).EquationId( ) );
            }
        }

        // Slave Nodes Displacement Equation IDs
        for ( unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave )
        {
            NodeType& slave_node = GetGeometry()[ i_slave ];
            rResult.push_back( slave_node.GetDof( DISPLACEMENT_X ).EquationId( ) );
            rResult.push_back( slave_node.GetDof( DISPLACEMENT_Y ).EquationId( ) );
            if (TDim == 3)
            {
                rResult.push_back( slave_node.GetDof( DISPLACEMENT_Z ).EquationId( ) );
            }
        }

        // Slave Nodes  Lambda Equation IDs
        for ( unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave )
        {
            NodeType& slave_node = GetGeometry()[ i_slave ];
            rResult.push_back( slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_X ).EquationId( ) );
            rResult.push_back( slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Y ).EquationId( ) );
            if (TDim == 3)
            {
                rResult.push_back( slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Z ).EquationId( ) );
            }
        }
    }

//     // Calculates the size of the system
//     constexpr unsigned int MatrixSize = TDim * (2 * TNumNodes + TNumNodes);
//     const unsigned int condition_size = this->CalculateConditionSize<MatrixSize>( );
//     
//     if (rResult.size() != condition_size)
//     {
//         rResult.resize( condition_size, false );
//     }
//     
//     unsigned int index = 0;
//         
//     /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
//     for (unsigned int PairIndex = 0; PairIndex < mThisMasterElements.size( ); ++PairIndex)
//     {
//         GeometryType& current_master = mThisMasterElements[PairIndex]->GetGeometry( );
//             
//         // Master Nodes Displacement Equation IDs
//         for ( unsigned int i_master = 0; i_master < TNumNodes; i_master++ ) //NOTE: Assuming the master has the same number of nodes
//         {
//             NodeType& master_node = current_master[i_master];
//             rResult[index++] =  master_node.GetDof( DISPLACEMENT_X ).EquationId( );
//             rResult[index++] =  master_node.GetDof( DISPLACEMENT_Y ).EquationId( );
//             if (TDim == 3)
//             {
//                 rResult[index++] =  master_node.GetDof( DISPLACEMENT_Z ).EquationId( );
//             }
//         }
//         
//         // Slave Nodes Displacement Equation IDs
//         for ( unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave )
//         {
//             NodeType& slave_node = GetGeometry()[ i_slave ];
//             rResult[index++] =  slave_node.GetDof( DISPLACEMENT_X ).EquationId( );
//             rResult[index++] =  slave_node.GetDof( DISPLACEMENT_Y ).EquationId( );
//             if (TDim == 3)
//             {
//                 rResult[index++] =  slave_node.GetDof( DISPLACEMENT_Z ).EquationId( );
//             }
//         }
//         
//         // Slave Nodes  Lambda Equation IDs
//         for ( unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave )
//         {
//             NodeType& slave_node = GetGeometry()[ i_slave ];
//             rResult[index++] =  slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_X ).EquationId( );
//             rResult[index++] =  slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Y ).EquationId( );
//             if (TDim == 3)
//             {
//                 rResult[index++] =  slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Z ).EquationId( );
//             }
//         }
//     }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void MortarContactCondition<TDim, TNumNodes>::GetDofList(
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo 
)
{
    KRATOS_TRY;

    rConditionalDofList.resize( 0 );
    const std::vector<contact_container> all_conditions = *( this->GetValue( CONTACT_CONTAINERS ) );

    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    for ( unsigned int i_cond = 0; i_cond < all_conditions.size( ); ++i_cond )
    {
        GeometryType& current_master = all_conditions[i_cond].condition->GetGeometry( );   

        // Master Nodes Displacement Equation IDs
        for ( unsigned int i_master = 0; i_master < TNumNodes; i_master++ ) // NOTE: Assuming same number of nodes for master and slave
        {
            NodeType& master_node = current_master[i_master];
            rConditionalDofList.push_back( master_node.pGetDof( DISPLACEMENT_X ) );
            rConditionalDofList.push_back( master_node.pGetDof( DISPLACEMENT_Y ) );
            if (TDim == 3)
            {
                rConditionalDofList.push_back( master_node.pGetDof( DISPLACEMENT_Z ) );
            }
        }

        // Slave Nodes Displacement Equation IDs
        for ( unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave )
        {
            NodeType& slave_node = GetGeometry()[ i_slave ];
            rConditionalDofList.push_back( slave_node.pGetDof( DISPLACEMENT_X ) );
            rConditionalDofList.push_back( slave_node.pGetDof( DISPLACEMENT_Y ) );
            if (TDim == 3)
            {
                rConditionalDofList.push_back( slave_node.pGetDof( DISPLACEMENT_Z ) );
            }
        }

        // Slave Nodes Lambda Equation IDs
        for ( unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave )
        {
            NodeType& slave_node = GetGeometry()[ i_slave ];
            rConditionalDofList.push_back( slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_X ) );
            rConditionalDofList.push_back( slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Y ) );
            if (TDim == 3)
            {
                rConditionalDofList.push_back( slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Z ) );
            }
        }
    }
    
//     // Calculates the size of the system
//     constexpr unsigned int MatrixSize = TDim * (2 * TNumNodes + TNumNodes);
//     const unsigned int condition_size = this->CalculateConditionSize<MatrixSize>( );
//     
//     if (rConditionalDofList.size() != condition_size)
//     {
//         rConditionalDofList.resize( condition_size );
//     }
//     
//     unsigned int index = 0;
//     
//     /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
//     for (unsigned int PairIndex = 0; PairIndex < mThisMasterElements.size( ); ++PairIndex)
//     {
//         GeometryType& current_master = mThisMasterElements[PairIndex]->GetGeometry( );
//         
//         // Master Nodes Displacement Equation IDs
//         for ( unsigned int i_master = 0; i_master < TNumNodes; i_master++ ) //NOTE: Assuming the master has the same number of nodes
//         {
//             NodeType& master_node = current_master[i_master];
//             rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_X );
//             rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_Y );
//             if (TDim == 3)
//             {
//                 rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_Z );
//             }
//         }
//         
//         // Slave Nodes Displacement Equation IDs
//         for ( unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave )
//         {
//             NodeType& slave_node = GetGeometry()[ i_slave ];
//             rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_X );
//             rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_Y );
//             if (TDim == 3)
//             {
//                 rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_Z );
//             }
//         }
//         
//         // Slave Nodes Lambda Equation IDs
//         for ( unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave )
//         {
//             NodeType& slave_node = GetGeometry()[ i_slave ];
//             rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_X );
//             rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Y );
//             if (TDim == 3)
//             {
//                 rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Z );
//             }
//         }
//     }

    KRATOS_CATCH( "" );
}


//******************************* GET DOUBLE VALUE *********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void MortarContactCondition<TDim,TNumNodes>::GetValueOnIntegrationPoints( 
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//******************************* GET ARRAY_1D VALUE *******************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void MortarContactCondition<TDim,TNumNodes>::GetValueOnIntegrationPoints( 
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//******************************* GET VECTOR VALUE *********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void MortarContactCondition<TDim,TNumNodes>::GetValueOnIntegrationPoints( 
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void MortarContactCondition<TDim,TNumNodes>::CalculateOnIntegrationPoints( 
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
    const unsigned int number_of_elements_master = mThisMasterElements.size( );
    
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
                    const double augmented_normal_lm = this->AugmentedNormalLM(Variables, rContactData, integration_point_gap);
                    
                    if (augmented_normal_lm < 0)
                    {
                        if ( rVariable == NORMAL_CONTACT_STRESS )
                        {
                            rOutput[PointNumber] += augmented_normal_lm;
                        }
                        else
                        {
                            // The slip of th GP
                            double integration_point_slip;
                        
                            // The tangent LM augmented
                            const double augmented_tangent_lm = this->AugmentedTangentLM(Variables, rContactData, current_master_element, integration_point_slip);
                            
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

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void MortarContactCondition<TDim,TNumNodes>::CalculateOnIntegrationPoints( 
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector< array_1d<double, 3 > >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
    
    // Create and initialize condition variables:
    GeneralVariables Variables;
    
    // Initialize the current contact data
    ContactData rContactData;

    // Slave segment info
    const unsigned int number_of_elements_master = mThisMasterElements.size( );
    
    // Reading integration points
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
    
    if (rVariable == NORMAL_GP || rVariable == TANGENT_GP)
    {
        if (number_of_elements_master > 0)
        {
            // Initialize general variables for the current master element
            this->InitializeGeneralVariables( Variables, rCurrentProcessInfo, 0 ); // NOTE: The pair does not matter
            
            // Calculating the values
            for ( unsigned int PointNumber = 0; PointNumber < number_of_integration_pts; PointNumber++ )
            {
                // Calculate the kinematic variables
                const bool inside = this->CalculateKinematics( Variables, PointNumber, 0, integration_points ); // NOTE: The pair does not matter
                
                if (inside == true)
                {
                    if (rVariable == NORMAL_GP)
                    {   
                        const array_1d<double, 3> normal_gp  = prod(trans(rContactData.NormalsSlave), Variables.N_Slave);
                        rOutput[PointNumber] = normal_gp;
                    }
                    else
                    {
                        const array_1d<double, 3> tangent_gp = prod(trans(rContactData.Tangent1Slave), Variables.N_Slave); 
                        rOutput[PointNumber] = tangent_gp;
                    }
                }
            }
        }
    }
    else
    {
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
                    if (rVariable == NORMAL_CONTACT_STRESS_GP || rVariable == TANGENTIAL_CONTACT_STRESS_GP)
                    {
                        // The normal LM augmented
                        const double augmented_normal_lm = this->AugmentedNormalLM(Variables, rContactData, integration_point_gap);
                        
                        if (augmented_normal_lm < 0)
                        {
                            if ( rVariable == NORMAL_CONTACT_STRESS_GP )
                            {
                                const array_1d<double, 3> normal_gp  = prod(trans(rContactData.NormalsSlave), Variables.N_Slave);
                                rOutput[PointNumber] += augmented_normal_lm * normal_gp;
                            }
                            else
                            {
                                // The slip of th GP
                                double integration_point_slip;
                            
                                // The tangent LM augmented
                                const double augmented_tangent_lm = this->AugmentedTangentLM(Variables, rContactData, current_master_element, integration_point_slip);
                                
                                const array_1d<double, 3> tangent_gp = prod(trans(rContactData.Tangent1Slave), Variables.N_Slave); // TODO: Extend to 3D case
                                rOutput[PointNumber] += augmented_tangent_lm * tangent_gp;
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

template< unsigned int TDim, unsigned int TNumNodes >
void MortarContactCondition<TDim,TNumNodes>::CalculateOnIntegrationPoints( 
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

template< unsigned int TDim, unsigned int TNumNodes >
double MortarContactCondition<TDim,TNumNodes>::AugmentedNormalLM(
    const GeneralVariables& rVariables,
    const ContactData& rContactData,
    const double& integration_point_gap
    )
{
    // The LM of the GP 
    const array_1d<double, TDim> lm_gp     = prod(trans(rContactData.LagrangeMultipliers), rVariables.Phi_LagrangeMultipliers);
    
    // The normal of the GP
    const array_1d<double, TDim> normal_gp = prod(trans(rContactData.NormalsSlave), rVariables.N_Slave);
    
    // The LM in the tangent direction
    double augmented_normal_lm  = inner_prod(normal_gp, lm_gp);// + rContactData.epsilon_normal * integration_point_gap;
                
    if (integration_point_gap < 0.0) // Penetration (TODO: Think about this!!!)
    {
        augmented_normal_lm += rContactData.epsilon_normal * integration_point_gap;
    }
                
    return augmented_normal_lm;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
double MortarContactCondition<TDim,TNumNodes>::AugmentedTangentLM(
    const GeneralVariables& rVariables,
    const ContactData& rContactData,
    const GeometryType& current_master_element, 
    double& integration_point_slip
    )
{    
    // The LM of the GP 
    const array_1d<double, TDim> lm_gp = prod(trans(rContactData.LagrangeMultipliers), rVariables.Phi_LagrangeMultipliers);
    
    // The tangents of the GP
    const array_1d<double, TDim> tangent_xi_gp  = prod(trans(rContactData.Tangent1Slave), rVariables.N_Slave);
    const array_1d<double, TDim> tangent_eta_gp = prod(trans(rContactData.Tangent2Slave), rVariables.N_Slave);
    
    // The tangential LM
    const double tangent_xi_lm  = inner_prod(tangent_xi_gp, lm_gp);
    const double tangent_eta_lm = inner_prod(tangent_eta_gp, lm_gp);
    const double tangent_lm = std::sqrt(tangent_xi_lm * tangent_xi_lm + tangent_eta_lm * tangent_eta_lm); 
    
    // The resultant direction of the LM
    const array_1d<double, TDim>& tangent_gp =  (tangent_xi_lm * tangent_xi_gp + tangent_eta_lm * tangent_eta_gp)/tangent_lm; // NOTE: This is the direction of the slip (using the LM as reference)
    
    // The velocities matrices
    const Matrix v1 = ContactUtilities::GetVariableMatrix(GetGeometry(),          VELOCITY, 0); 
    const Matrix v2 = ContactUtilities::GetVariableMatrix(current_master_element, VELOCITY, 0);
    
    // The slip of the LM
    const array_1d<double, TDim> vector_integration_point_slip = rContactData.Dt * (prod(trans(v1), rVariables.N_Slave) - prod(trans(v2), rVariables.N_Master));
    integration_point_slip = inner_prod(vector_integration_point_slip, tangent_gp);

    // The augmented LM in the tangent direction
    const double augmented_tangent_lm = tangent_lm + rContactData.epsilon_tangent * integration_point_slip; 

    return augmented_tangent_lm;
}

/***********************************************************************************/
/***********************************************************************************/

template< >
void MortarContactCondition<2, 2>::InitializeIntegrationMethod()
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
void MortarContactCondition<2, 3>::InitializeIntegrationMethod()
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
void MortarContactCondition<3, 3>::InitializeIntegrationMethod()
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
void MortarContactCondition<3, 4>::InitializeIntegrationMethod()
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

template class MortarContactCondition<2,2>;
template class MortarContactCondition<2,3>;
template class MortarContactCondition<3,3>;
template class MortarContactCondition<3,4>;

} // Namespace Kratos
