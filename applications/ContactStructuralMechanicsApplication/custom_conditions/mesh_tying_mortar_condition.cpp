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
#include "custom_conditions/mesh_tying_mortar_condition.h"

/* Additional includes */
#include <algorithm>

/* Utilities */
#include "custom_utilities/contact_utilities.h"
#include "custom_utilities/exact_mortar_segmentation_utility.h"
#include "utilities/math_utils.h"
#include "custom_utilities/structural_mechanics_math_utilities.hpp" // NOTE: Change for a more performant solver

namespace Kratos 
{
/**
 * Flags related to the condition computation 
 */
// Avoiding using the macro since this has a template parameter. If there was no template plase use the KRATOS_CREATE_LOCAL_FLAG macro
template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
const Kratos::Flags MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR(Kratos::Flags::Create(0));
template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
const Kratos::Flags MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX(Kratos::Flags::Create(1));
template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
const Kratos::Flags MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS(Kratos::Flags::Create(2));
template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
const Kratos::Flags MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS(Kratos::Flags::Create(3));

/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
Condition::Pointer MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::Create( 
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties ) const
{
    return boost::make_shared< MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor> >( NewId, this->GetGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
Condition::Pointer MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return boost::make_shared< MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor> >( NewId, pGeom, pProperties );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::~MeshTyingMortarCondition( )
{
}

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::Initialize( ) 
{
    KRATOS_TRY;
    
    // Pointer to the reference element 
    mThisSlaveElement = this->GetValue(ELEMENT_POINTER);
    
    // Populate of the vector of master elements (it is supposed to be constant)
    const std::vector<contact_container> * all_containers = this->GetValue( CONTACT_CONTAINERS );
    
    const unsigned int IntegrationOrder = GetProperties().GetValue(INTEGRATION_ORDER_CONTACT);
    
    ExactMortarIntegrationUtility<TDim,TNumNodes> IntegrationUtility = ExactMortarIntegrationUtility<TDim,TNumNodes>(GetGeometry(), IntegrationOrder);
//     mPairSize = all_containers->size();
//     mThisMasterConditions.resize( mPairSize );
//     
    mPairSize = 0;
    for ( unsigned int i_cond = 0; i_cond < all_containers->size(); ++i_cond )
    {
        Condition::Pointer pCond = (*all_containers)[i_cond].condition;
        
        IntegrationPointsType IntegrationPoints;
        const bool is_inside = IntegrationUtility.GetExactIntegration(pCond->GetGeometry(), IntegrationPoints);
        
        if (is_inside == true)
        {
            mPairSize += 1;
            mThisMasterConditions.push_back(pCond);
            mThisMasterElements.push_back(pCond->GetValue(ELEMENT_POINTER));
            mThisSlaveIntegrationPoints.push_back(IntegrationPoints);
        }
    }
    
    // Creation of the variables 
    if (TTensor == ScalarValue)
    {
        mTyingVarScalar = KratosComponents<Variable<double>>::Get(GetProperties().GetValue(TYING_VARIABLE));
    }
    else
    {
        mTyingVarVector[0] = KratosComponents<Variable<array_1d_component_type>>::Get(GetProperties().GetValue(TYING_VARIABLE)+"_X");
        mTyingVarVector[1] = KratosComponents<Variable<array_1d_component_type>>::Get(GetProperties().GetValue(TYING_VARIABLE)+"_Y");
        if (TDim == 3)
        {
            mTyingVarVector[2] = KratosComponents<Variable<array_1d_component_type>>::Get(GetProperties().GetValue(TYING_VARIABLE)+"_Z");
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    // NOTE: Add things if necessary
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    // NOTE: Add things if necessary
        
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // NOTE: Add things if necessary
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // TODO: Recalculate Lagrange multiplers
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::CalculateLocalSystem( 
    std::vector<MatrixType>& rLeftHandSideMatrices,
    const std::vector<Variable<MatrixType> >& rLHSVariables,
    std::vector<VectorType>& rRightHandSideVectors,
    const std::vector<Variable<VectorType> >& rRHSVariables,
    ProcessInfo& rCurrentProcessInfo 
    )
{    
    // Calculates the size of the system
    constexpr unsigned int TMatrixSize = TDim * TNumNodesElem;
    
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set(MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS, true);
    LocalSystem.CalculationFlags.Set(MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS, true);

    //Initialize sizes for the system components
    if ( rLHSVariables.size( ) != rLeftHandSideMatrices.size( ) )
    {
        rLeftHandSideMatrices.resize( rLHSVariables.size( ) );
    }

    if ( rRHSVariables.size( ) != rRightHandSideVectors.size( ) )
    {
        rRightHandSideVectors.resize( rRHSVariables.size( ) );
    }

    LocalSystem.CalculationFlags.Set(MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX, true);
    for ( unsigned int i = 0; i < rLeftHandSideMatrices.size( ); i++ )
    {
        // Note: rRightHandSideVectors.size() > 0
        this->InitializeSystemMatrices<TMatrixSize>( rLeftHandSideMatrices[i], rRightHandSideVectors[0],LocalSystem.CalculationFlags );
    }

    LocalSystem.CalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR, true );
    LocalSystem.CalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX, false ); // Temporarily only
    for ( unsigned int i = 0; i < rRightHandSideVectors.size( ); i++ )
    {
        // Note: rLeftHandSideMatrices.size() > 0
        this->InitializeSystemMatrices<TMatrixSize>( rLeftHandSideMatrices[0], rRightHandSideVectors[i], LocalSystem.CalculationFlags  );
    }
    LocalSystem.CalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX, true ); // Reactivated again

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

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;
    
    // Calculates the size of the system
    constexpr unsigned int TMatrixSize = TDim * TNumNodesElem;
    
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX, true );
    LocalSystem.CalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR, true );

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

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::CalculateLeftHandSide( 
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Calculates the size of the system
    constexpr unsigned int TMatrixSize = TDim * TNumNodesElem;
    
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX, true );

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

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::CalculateLeftHandSide( 
    std::vector< MatrixType >& rLeftHandSideMatrices,
    const std::vector< Variable< MatrixType > >& rLHSVariables,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Calculates the size of the system
    constexpr unsigned int TMatrixSize = TDim * TNumNodesElem;
    
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX, true );
    LocalSystem.CalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS, true );

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

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::CalculateRightHandSide( 
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Calculates the size of the system
    constexpr unsigned int TMatrixSize = TDim * TNumNodesElem;
    
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR, true);

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

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::CalculateRightHandSide( 
    std::vector< VectorType >& rRightHandSideVectors,
    const std::vector< Variable< VectorType > >& rRHSVariables,
    ProcessInfo& rCurrentProcessInfo )
{
    // Calculates the size of the system
    constexpr unsigned int TMatrixSize = TDim * TNumNodesElem;
        
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR, true );
    LocalSystem.CalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS, true );

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

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
template< unsigned int TMatrixSize >
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::InitializeSystemMatrices( 
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    Flags& rCalculationFlags 
    )
{
    const unsigned int condition_size = this->CalculateConditionSize<TMatrixSize>( );
    
    // Resizing as needed the LHS
    if ( rCalculationFlags.Is( MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX ) ) // Calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != condition_size )
        {
            rLeftHandSideMatrix.resize( condition_size, condition_size, false );
        }
        noalias( rLeftHandSideMatrix ) = ZeroMatrix( condition_size, condition_size ); // Resetting LHS
    }

    // Resizing as needed the RHS
    if ( rCalculationFlags.Is( MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR ) ) // Calculation of the matrix is required
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

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::CalculateMassMatrix( 
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

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::CalculateDampingMatrix( 
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

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
template< unsigned int TMatrixSize >
const unsigned int MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::CalculateConditionSize( )
{
    const unsigned int condition_size = mPairSize * TMatrixSize;
    
    return condition_size;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
template< unsigned int TMatrixSize>
void MeshTyingMortarCondition<TDim, TNumNodes>::CalculateConditionSystem( 
    LocalSystemComponents& rLocalSystem,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
    
//                 // DoF of the slave       
//             if (TTensor == ScalarValue)
//             {
//                 for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
//                 {
//                     const double var  = SlaveGeometry[iNode].FastGetSolutionStepValue(mTyingVarScalar);
//                     u1(iNode, 0) = var;
//                 }
//             }
//             else
//             {
//                 for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
//                 {
//                     for (unsigned int iDof = 0; iDof < TDim; iDof++)
//                     {
//                         const double var = SlaveGeometry[iNode].FastGetSolutionStepValue(mTyingVarVector[iDof]);
//                         u1(iNode, iDof) = var;
//                     }
//                 }
//             }
    
    // Create and initialize condition variables:#pragma omp critical
    GeneralVariables rVariables;
    
    // Create the current contact data
    DofData rDofData;
    
    // Create the mortar operators
    MortarConditionMatrices rThisMortarConditionMatrices;
                                                                                  
    this->InitializeDofData(rDofData, rCurrentProcessInfo);
    
    // Compute Ae and its derivative
    this->CalculateAe(rDofData, rVariables, rCurrentProcessInfo); 
    
    // Iterate over the master segments
    for (unsigned int PairIndex = 0; PairIndex < mPairSize; ++PairIndex)
    {   
        // Reading integration points
//         this->ComputeSelectiveIntegrationMethod(PairIndex);
//         const GeometryType::IntegrationPointsArrayType& integration_points = mIntegrationPoints.IntegrationPoints( );
                                                                         
        // Initialize general variables for the current master element
        this->InitializeGeneralVariables( rVariables, rCurrentProcessInfo, PairIndex );
        
        // Update the contact data
        this->UpdateDofData(rDofData, PairIndex);
        
        // Initialize the mortar operators
        rThisMortarConditionMatrices.Initialize();
        
        // Initialize the integration weight
        double total_weight = 0.0;
        
        // Integrating the mortar operators
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            // Calculate the kinematic variables
            this->CalculateKinematics( rVariables, rDofData, PointNumber, integration_points );
            
            this->CalculateMortarOperators(rThisMortarConditionMatrices, rVariables, IntegrationWeight);
        }
        
        // We can consider the pair if at least one of the collocation point is inside 
        if (total_weight > 0.0)
        {            
//             // Debug
//             std::cout << "--------------------------------------------------" << std::endl;
//             KRATOS_WATCH(this->Id());
//             KRATOS_WATCH(PairIndex);
//             rThisMortarConditionMatrices.print();
            
            // Calculates the active/inactive combination pair
            const unsigned int ActiveInactive = GetActiveInactiveValue(this->GetGeometry());
            
            // Assemble of the matrix is required
            if ( rLocalSystem.CalculationFlags.Is( MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX ) ||
                    rLocalSystem.CalculationFlags.Is( MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
            {
                // Calculate the local contribution
                bounded_matrix<double, TMatrixSize, TMatrixSize> LHS_contact_pair = this->CalculateLocalLHS<TMatrixSize>( rThisMortarConditionMatrices, PairIndex, ActiveInactive);
                
                // Contributions to stiffness matrix calculated on the reference config
                this->CalculateAndAddLHS<TMatrixSize>( rLocalSystem, LHS_contact_pair, PairIndex );
                
//                 // Debug
// //                 KRATOS_WATCH(LHS_contact_pair);
//                 LOG_MATRIX_PRETTY( LHS_contact_pair );
            }

            // Assemble of the vector is required
            if ( rLocalSystem.CalculationFlags.Is( MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR ) ||
                    rLocalSystem.CalculationFlags.Is( MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
            {
                // Calculate the local contribution
                const array_1d<double, TMatrixSize> RHS_contact_pair = this->CalculateLocalRHS<TMatrixSize>( rThisMortarConditionMatrices, PairIndex, ActiveInactive);
                
                // Contribution to previous step contact force and residuals vector
                this->CalculateAndAddRHS<TMatrixSize>( rLocalSystem, RHS_contact_pair, PairIndex );
                
//                 // Debug
// //                 KRATOS_WATCH(RHS_contact_pair);
//                 LOG_VECTOR_PRETTY( RHS_contact_pair );
            }
        }
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::InitializeGeneralVariables(
    GeneralVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo,
    const unsigned int& rMasterElementIndex
    )
{
    // Master segment info
    GeometryType& CurrentMasterElement = mThisMasterConditions[rMasterElementIndex]->GetGeometry();
    
    // Slave element info
    rVariables.Initialize();

    rVariables.SetMasterElement( CurrentMasterElement );
    rVariables.SetMasterElementIndex( rMasterElementIndex );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::CalculateAe(
    DofData& rDofData,
    GeneralVariables& rVariables,
//     const GeometryType::IntegrationPointsArrayType& integration_points,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    double total_weight = 0.0; // NOTE: The integral is supposed to be in the domain partially integrated, I don't know if consider any additional thing for the partial integration
    
    // We initilize the Ae components
    AeData rAeData;
    rAeData.Initialize();
    
    rDofData.InitializeAeComponents();
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
        this->UpdateDofData(rDofData, PairIndex);
            
        // Calculating the proportion between the integrated area and segment area
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            // Calculate the kinematic variables
            this->CalculateKinematics( rVariables, rDofData, PointNumber, integration_points );
            
            const double IntegrationWeight = integration_points[PointNumber].Weight();
            total_weight += IntegrationWeight;
            this->CalculateAeComponents(rVariables, rDofData, rAeData, IntegrationWeight);
        }
    }
    
    // We can consider the pair if at least one of the collocation point is inside (TODO: Change this if collocation is not used)
    if (total_weight > 0.0)
    {
        this->CalculateAe(rDofData, rAeData);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::InitializeDofData(
    DofData& rDofData,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // Slave element info
    rDofData.Initialize(GetGeometry());
    
    /* LM */
    rDofData.LagrangeMultipliers = ContactUtilities::GetVariableVector<TNumNodes>(GetGeometry(), NORMAL_CONTACT_STRESS, 0); 
    
    /* NORMALS */
    rDofData.Normal_s = ContactUtilities::GetVariableMatrix<TDim,TNumNodes,TNumNodesElem,TTensor>(GetGeometry(),  NORMAL); 
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::UpdateDofData(
    DofData& rDofData,
    const unsigned int& rMasterElementIndex
    )
{    
    // Slave element info
    rDofData.UpdateMasterPair(mThisMasterConditions[rMasterElementIndex] );
    
    /* NORMALS */
    for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
    {
//         array_1d<double,3> normal = this->GetValue(NORMAL);
        array_1d<double,3> normal = GetGeometry()[iNode].GetValue(NORMAL);
    }
}

/*********************************COMPUTE KINEMATICS*********************************/
/************************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::CalculateKinematics( 
    GeneralVariables& rVariables,
    const DofData rDofData,
    const double& rPointNumber,
    const GeometryType::IntegrationPointsArrayType& integration_points
    )
{
    /* LOCAL COORDINATES */
    const PointType& local_point = integration_points[rPointNumber].Coordinates();
    
    /*  POPULATE MATRICES AND VECTORS */
    
    /// SLAVE CONDITION ///
    
    // SHAPE FUNCTIONS 
    GetGeometry( ).ShapeFunctionsValues( rVariables.N_Slave, local_point.Coordinates() );
    rVariables.Phi_LagrangeMultipliers = prod(rDofData.Ae, rVariables.N_Slave);
//     rVariables.Phi_LagrangeMultipliers = rVariables.N_Slave; // TODO: This could be needed in the future to be different than the standart shape functions 
    
    /// MASTER CONDITION ///
    this->MasterShapeFunctionValue( rVariables, local_point);
    
    /* CALCULATE JACOBIAN AND JACOBIAN DETERMINANT */
    GetGeometry( ).Jacobian( rVariables.j_Slave, local_point.Coordinates() );
    rVariables.DetJSlave = ContactUtilities::ContactElementDetJacobian( rVariables.j_Slave );
}
 
/***********************************************************************************/
/*************** METHODS TO CALCULATE THE CONTACT CONDITION MATRICES ***************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
bool MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::MasterShapeFunctionValue(
    GeneralVariables& rVariables,
    const PointType& local_point 
    )
{    
    GeometryType& master_seg = rVariables.GetMasterElement( );

    PointType projected_gp_global;
    const array_1d<double,3> normal = ContactUtilities::GaussPointNormal(rVariables.N_Slave, GetGeometry());
    
    GeometryType::CoordinatesArrayType slave_gp_global;
    double aux_dist = 0.0;
    this->GetGeometry( ).GlobalCoordinates( slave_gp_global, local_point );
    ContactUtilities::ProjectDirection( master_seg, slave_gp_global, projected_gp_global, aux_dist, -normal ); // The opposite direction
    
    GeometryType::CoordinatesArrayType projected_gp_local;
    
    master_seg.PointLocalCoordinates(projected_gp_local, projected_gp_global.Coordinates( ) ) ;
    
    // SHAPE FUNCTIONS 
    master_seg.ShapeFunctionsValues(         rVariables.N_Master,     projected_gp_local );         
    
    return inside;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::CalculateMortarOperators(
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
            const double phi = Phi[i_slave];
            
            DOperator(i_slave, j_slave) += J_s * rIntegrationWeight * phi * N1[j_slave];
            MOperator(i_slave, j_slave) += J_s * rIntegrationWeight * phi * N2[j_slave];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
bounded_matrix<double, TNumNodes, TNumNodes> MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::ComputeDe(
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

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
bounded_matrix<double, TNumNodes, TNumNodes> MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::ComputeMe(
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

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::CalculateAe(
    DofData& rDofData,
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
    
    noalias(rDofData.Ae) = prod(rAeData.De, InvMe);
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
template< unsigned int TMatrixSize >
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::CalculateAndAddLHS(
    LocalSystemComponents& rLocalSystem,
    const bounded_matrix<double, TMatrixSize, TMatrixSize>& LHS_contact_pair, 
    const unsigned int rPairIndex
    )
{
    /* DEFINITIONS */

    if ( rLocalSystem.CalculationFlags.Is( MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
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
                KRATOS_ERROR <<  " CONDITION can not supply the required local system variable: " << rLeftHandSideVariables[i] << std::endl;
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

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
template< unsigned int TMatrixSize>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::AssembleContactPairLHSToConditionSystem(
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

/****************************** END AD REPLACEMENT *********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
template< unsigned int TMatrixSize >
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::CalculateAndAddRHS(
    LocalSystemComponents& rLocalSystem,
    const array_1d<double, TMatrixSize>& RHS_contact_pair,
    const unsigned int rPairIndex
    )
{   
    if ( rLocalSystem.CalculationFlags.Is( MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
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
                KRATOS_ERROR << " CONDITION can not supply the required local system variable: " << rRightHandSideVariables[i] << std::endl;
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

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
template< unsigned int TMatrixSize>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::AssembleContactPairRHSToConditionSystem(
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

/****************************** END AD REPLACEMENT *********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::EquationIdVector(
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
        
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim, TNumNodes>::GetDofList(
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
    }
    
    KRATOS_CATCH( "" );
}


//******************************* GET DOUBLE VALUE *********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::GetValueOnIntegrationPoints( 
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//******************************* GET ARRAY_1D VALUE *******************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::GetValueOnIntegrationPoints( 
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//******************************* GET VECTOR VALUE *********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::GetValueOnIntegrationPoints( 
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::CalculateOnIntegrationPoints( 
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;

    // Create and initialize condition variables:
    GeneralVariables rVariables;
    
    // Initialize the current contact data
    DofData rDofData;
    
        // Reading integration points
//         this->ComputeSelectiveIntegrationMethod(PairIndex);
//         const GeometryType::IntegrationPointsArrayType& integration_points = mIntegrationPoints.IntegrationPoints( );
                                               
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

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::CalculateOnIntegrationPoints( 
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector< array_1d<double, 3 > >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
    
    // Create and initialize condition variables:
    GeneralVariables rVariables;
    
    // Initialize the current contact data
    DofData rDofData;
    
        // Reading integration points
//         this->ComputeSelectiveIntegrationMethod(PairIndex);
//         const GeometryType::IntegrationPointsArrayType& integration_points = mIntegrationPoints.IntegrationPoints( );
                                                                                                                        
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

template< unsigned int TDim, unsigned int TNumNodes, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodes,TNumNodesElem,TTensor>::CalculateOnIntegrationPoints( 
    const Variable<Vector>& rVariable, 
    std::vector<Vector>& rOutput, 
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
    
    // TODO: Fill this!!!
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template class MeshTyingMortarCondition<2, 2, DISPLACEMENT>;
template class MeshTyingMortarCondition<3, 3, DISPLACEMENT>;
template class MeshTyingMortarCondition<3, 4, DISPLACEMENT>;
template class MeshTyingMortarCondition<2, 2, TEMPERATURE>;
template class MeshTyingMortarCondition<3, 3, TEMPERATURE>;
template class MeshTyingMortarCondition<3, 4, TEMPERATURE>;

} // Namespace Kratos
