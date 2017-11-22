// KRATOS  ___|  |       |       |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//           | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License: BSD License
//   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
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
#include "utilities/exact_mortar_segmentation_utility.h"
#include "utilities/math_utils.h"

namespace Kratos 
{
/**
 * Flags related to the condition computation 
 */
// Avoiding using the macro since this has a template parameter. If there was no template plase use the KRATOS_CREATE_LOCAL_FLAG macro
template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
const Kratos::Flags MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR(Kratos::Flags::Create(0));
template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
const Kratos::Flags MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX(Kratos::Flags::Create(1));
template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
const Kratos::Flags MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS(Kratos::Flags::Create(2));
template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
const Kratos::Flags MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS(Kratos::Flags::Create(3));

/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
Condition::Pointer MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::Create( 
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties ) const
{
    return boost::make_shared< MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor> >( NewId, this->GetGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
Condition::Pointer MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return boost::make_shared< MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor> >( NewId, pGeom, pProperties );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::~MeshTyingMortarCondition( )
{
}

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::Initialize( ) 
{
    KRATOS_TRY;
    
    // Pointer to the reference element 
//     mThisSlaveElement = this->GetValue(ELEMENT_POINTER);
    
    // Populate of the vector of master elements (it is supposed to be constant)    
    ConditionMap::Pointer& AllConditionMaps = this->GetValue( MAPPING_PAIRS );
    
    mIntegrationOrder = GetProperties().GetValue(INTEGRATION_ORDER_CONTACT);

    IntegrationMethod ThisIntegrationMethod = GetIntegrationMethod();
    
    mPairSize = 0;
    ExactMortarIntegrationUtility<TDim, NumNodes> IntUtil = ExactMortarIntegrationUtility<TDim, NumNodes>(mIntegrationOrder);
    
    // Create and initialize condition variables:
    GeneralVariables rVariables;
    
    for (auto itPair = AllConditionMaps->begin(); itPair != AllConditionMaps->end(); ++itPair )
    {
        Condition::Pointer pCond = *(itPair);
        
        // Reading integration points
        ConditionArrayListType ConditionsPointsSlave;
        const bool IsInside = IntUtil.GetExactIntegration(this->GetGeometry(), this->GetValue(NORMAL), pCond->GetGeometry(), pCond->GetValue(NORMAL), ConditionsPointsSlave);
        
        if (IsInside == true)
        {
            GeometryType::IntegrationPointsArrayType AllIntegrationPointsSlave;
            
            for (unsigned int i_geom = 0; i_geom < ConditionsPointsSlave.size(); i_geom++)
            {
                std::vector<PointType::Pointer> PointsArray (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
                for (unsigned int i_node = 0; i_node < TDim; i_node++)
                {
                    PointType GlobalPoint;
                    GetGeometry().GlobalCoordinates(GlobalPoint, ConditionsPointsSlave[i_geom][i_node]);
                    PointsArray[i_node] = boost::make_shared<PointType>(GlobalPoint);
                }
                
                DecompositionType DecompGeom( PointsArray );
                
                const bool BadShape = (TDim == 2) ? false : MortarUtilities::HeronCheck(DecompGeom);
                
                if (BadShape == false)
                {
                    const GeometryType::IntegrationPointsArrayType& IntegrationPointsSlave = DecompGeom.IntegrationPoints( ThisIntegrationMethod );
                    
                    // Integrating the mortar operators
                    for ( unsigned int PointNumber = 0; PointNumber < IntegrationPointsSlave.size(); PointNumber++ )
                    {
                        const double Weight = IntegrationPointsSlave[PointNumber].Weight();
                        const PointType LocalPointDecomp = IntegrationPointsSlave[PointNumber].Coordinates();
                        PointType LocalPointParent;
                        PointType GPGlobal;
                        DecompGeom.GlobalCoordinates(GPGlobal, LocalPointDecomp);
                        GetGeometry().PointLocalCoordinates(LocalPointParent, GPGlobal);
                        
                        const double DetJ = DecompGeom.DeterminantOfJacobian( LocalPointDecomp );
                        
                        AllIntegrationPointsSlave.push_back( IntegrationPointType( LocalPointParent.Coordinate(1), LocalPointParent.Coordinate(2), Weight * DetJ ));
                    }
                }
            }
 
            if (AllIntegrationPointsSlave.size() > 0)
            {
                mPairSize += 1;
                mThisMasterConditions.push_back(pCond);
                mIntegrationPointsVector.push_back(AllIntegrationPointsSlave);
//                 mIntegrationPointsVector.push_back(IntegrationPointsSlave);
//                 mThisMasterElements.push_back(pCond->GetValue(ELEMENT_POINTER));
            }
        }
    }
        
//     // Creation of the variables 
//     if (TTensor == ScalarValue)
//     {
//         mTyingVarScalar = KratosComponents<Variable<double>>::Get(GetProperties().GetValue(TYING_VARIABLE));
//     }
//     else
//     {
//         mTyingVarVector[0] = KratosComponents<Variable<array_1d_component_type>>::Get(GetProperties().GetValue(TYING_VARIABLE)+"_X");
//         mTyingVarVector[1] = KratosComponents<Variable<array_1d_component_type>>::Get(GetProperties().GetValue(TYING_VARIABLE)+"_Y");
//         if (TDim == 3)
//         {
//             mTyingVarVector[2] = KratosComponents<Variable<array_1d_component_type>>::Get(GetProperties().GetValue(TYING_VARIABLE)+"_Z");
//         }
//     }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // NOTE: Add things if necessary
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // NOTE: Add things if necessary
        
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // NOTE: Add things if necessary
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // TODO: Recalculate Lagrange multiplers and Slave DoF (remember is increment, so you need to add to the existing ones)
    
    // Equation: $$ invD (-RHSslave - KSN Delta uElements - KSS DeltauS) $$
    
//     // DoF of the slave       // TODO: Add the values of the elements
//     if (TTensor == ScalarValue)
//     {
//         for (unsigned int iNode = 0; iNode < NumNodes; iNode++)
//         {
//             const double var  = SlaveGeometry[iNode].FastGetSolutionStepValue(mTyingVarScalar);
//             u1(iNode, 0) = var;
//         }
//     }
//     else
//     {
//         for (unsigned int iNode = 0; iNode < NumNodes; iNode++)
//         {
//             for (unsigned int iDof = 0; iDof < TDim; iDof++)
//             {
//                 const double var = SlaveGeometry[iNode].FastGetSolutionStepValue(mTyingVarVector[iDof]);
//                 u1(iNode, iDof) = var;
//             }
//         }
//     }
    
//     if (TTensor == ScalarValue) // SCALAR_LAGRANGE_MULTIPLIER
//     {
// 
//     }
//     else // VECTOR_LAGRANGE_MULTIPLIER
//     {
//         
//     }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateLocalSystem( 
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
    LocalSystem.CalculationFlags.Set(MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS, true);
    LocalSystem.CalculationFlags.Set(MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS, true);

    //Initialize sizes for the system components
    if ( rLHSVariables.size( ) != rLeftHandSideMatrices.size( ) )
    {
        rLeftHandSideMatrices.resize( rLHSVariables.size( ) );
    }

    if ( rRHSVariables.size( ) != rRightHandSideVectors.size( ) )
    {
        rRightHandSideVectors.resize( rRHSVariables.size( ) );
    }

    LocalSystem.CalculationFlags.Set(MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX, true);
    for ( unsigned int i = 0; i < rLeftHandSideMatrices.size( ); i++ )
    {
        // Note: rRightHandSideVectors.size() > 0
        this->InitializeSystemMatrices( rLeftHandSideMatrices[i], rRightHandSideVectors[0],LocalSystem.CalculationFlags );
    }

    LocalSystem.CalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR, true );
    LocalSystem.CalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX, false ); // Temporarily only
    for ( unsigned int i = 0; i < rRightHandSideVectors.size( ); i++ )
    {
        // Note: rLeftHandSideMatrices.size() > 0
        this->InitializeSystemMatrices( rLeftHandSideMatrices[0], rRightHandSideVectors[i], LocalSystem.CalculationFlags  );
    }
    LocalSystem.CalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX, true ); // Reactivated again

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

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;

    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX, true );
    LocalSystem.CalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR, true );

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

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateLeftHandSide( 
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX, true );

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

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateLeftHandSide( 
    std::vector< MatrixType >& rLeftHandSideMatrices,
    const std::vector< Variable< MatrixType > >& rLHSVariables,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX, true );
    LocalSystem.CalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS, true );

    VectorType RightHandSideVector = Vector( );

    // Initialize size for the system components
    for( unsigned int i = 0; i < rLeftHandSideMatrices.size( ); i++ )
    {
        this->InitializeSystemMatrices( rLeftHandSideMatrices[i], RightHandSideVector, LocalSystem.CalculationFlags );
    }

    LocalSystem.SetLeftHandSideMatrices( rLeftHandSideMatrices );
    LocalSystem.SetRightHandSideVector( RightHandSideVector );

    // Calculate condition system
    this->CalculateConditionSystem( LocalSystem, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateRightHandSide( 
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR, true);

    MatrixType LeftHandSideMatrix = Matrix( );

    // Initialize size for the system components
    this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVector,LocalSystem.CalculationFlags);

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix( LeftHandSideMatrix );
    LocalSystem.SetRightHandSideVector( rRightHandSideVector );

    // Calculate condition system
    this->CalculateConditionSystem( LocalSystem, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateRightHandSide( 
    std::vector< VectorType >& rRightHandSideVectors,
    const std::vector< Variable< VectorType > >& rRHSVariables,
    ProcessInfo& rCurrentProcessInfo )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR, true );
    LocalSystem.CalculationFlags.Set( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS, true );

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

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::InitializeSystemMatrices( 
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    Flags& rCalculationFlags 
    )
{
    const unsigned int ConditionSize = this->CalculateConditionSize( );
    
    // Resizing as needed the LHS
    if ( rCalculationFlags.Is( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX ) ) // Calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != ConditionSize )
        {
            rLeftHandSideMatrix.resize( ConditionSize, ConditionSize, false );
        }
        noalias( rLeftHandSideMatrix ) = ZeroMatrix( ConditionSize, ConditionSize ); // Resetting LHS
    }

    // Resizing as needed the RHS
    if ( rCalculationFlags.Is( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR ) ) // Calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != ConditionSize )
        {
            rRightHandSideVector.resize( ConditionSize, false );
        }
        rRightHandSideVector = ZeroVector( ConditionSize ); // Resetting RHS
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateMassMatrix( 
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

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateDampingMatrix( 
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

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
const unsigned int MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateConditionSize( )
{
    const unsigned int ConditionSize = mPairSize * MatrixSize;
    
    return ConditionSize;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim, TNumNodesElem, TTensor>::CalculateConditionSystem( 
    LocalSystemComponents& rLocalSystem,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
    
    // Create and initialize condition variables:
    GeneralVariables rVariables;
    
    // Create the current DoF data
    DofData rDofData;
    
    // Create the mortar operators
    MortarConditionMatrices rThisMortarConditionMatrices;
                                            
    // Initialize the DoF data
    this->InitializeDofData(rDofData, rCurrentProcessInfo);
    
    // Compute Ae and its derivative
    this->CalculateAe(rDofData, rVariables, rCurrentProcessInfo); 
    
//     // We calculate the Equation ID, LHS and RHS of the slave parent element
//     boost::numeric::ublas::bounded_matrix<double, DimensionLocalElem, DimensionLocalElem> LHS_SlaveElem_Contribution;
//     array_1d<double, DimensionLocalElem> RHS_SlaveElem_Contribution;
//     Element::EquationIdVectorType& EquationIdSlaveElem;
//     
//     (mThisSlaveElement) -> EquationIdVector(EquationIdSlaveElem, rCurrentProcessInfo);
//     
//     if ( rLocalSystem.CalculationFlags.Is( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX ) ||
//                 rLocalSystem.CalculationFlags.Is( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
//     {
//         // We calculate the contributions of the local elements
//         (mThisSlaveElement) -> CalculateLeftHandSide(LHS_SlaveElem_Contribution,rCurrentProcessInfo);
//     }
//     if ( rLocalSystem.CalculationFlags.Is( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR ) ||
//             rLocalSystem.CalculationFlags.Is( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
//     {
//         // We calculate the contributions of the local elements
//         (mThisSlaveElement) -> CalculateRightHandSide(RHS_SlaveElem_Contribution,rCurrentProcessInfo);
//     }
    
    // Iterate over the master segments
//     ExactMortarIntegrationUtility<TDim, NumNodes> IntUtil = ExactMortarIntegrationUtility<TDim, NumNodes>(mIntegrationOrder);
    for (unsigned int PairIndex = 0; PairIndex < mPairSize; ++PairIndex)
    {
        // The normal of the master condition
        const array_1d<double, 3> MasterNormal = mThisMasterConditions[PairIndex]->GetValue(NORMAL);
        
        // Initialize general variables for the current master element
        this->InitializeGeneralVariables( rVariables, rCurrentProcessInfo, PairIndex );
        
        // Update master pair info
        rDofData.UpdateMasterPair(mThisMasterConditions[PairIndex]);
        
        // Initialize the mortar operators
        rThisMortarConditionMatrices.Initialize();
        
        // Get the integration points
        const IntegrationPointsType IntegrationPointsSlave = mIntegrationPointsVector[PairIndex];
//         IntUtil.GetExactIntegration(this->GetGeometry(), this->GetValue(NORMAL), mThisMasterConditions[PairIndex]->GetGeometry(), mThisMasterConditions[PairIndex]->GetValue(NORMAL), IntegrationPointsSlave);
        
        const unsigned int NumberOfIntegrationPoints = IntegrationPointsSlave.size();
        
        // Integrating the mortar operators
        for ( unsigned int PointNumber = 0; PointNumber < NumberOfIntegrationPoints; PointNumber++ )
        {            
            // Calculate the kinematic variables
            this->CalculateKinematics( rVariables, rDofData, MasterNormal, PointNumber, IntegrationPointsSlave );
            
            this->CalculateMortarOperators(rThisMortarConditionMatrices, rVariables, IntegrationPointsSlave[PointNumber].Weight());
        }
                
        if (NumberOfIntegrationPoints > 0)
        {
//             // Debug
//             std::cout << "--------------------------------------------------" << std::endl;
//             KRATOS_WATCH(this->Id());
//             KRATOS_WATCH(PairIndex);
//             rThisMortarConditionMatrices.print();
            
            // Assemble of the matrix is required
            if ( rLocalSystem.CalculationFlags.Is( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX ) ||
                    rLocalSystem.CalculationFlags.Is( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
            {        
                // Calculate the local contribution
                const boost::numeric::ublas::bounded_matrix<double, MatrixSize, MatrixSize> LHS_contact_pair = this->CalculateLocalLHS<MatrixSize>( rThisMortarConditionMatrices, rDofData, PairIndex, rCurrentProcessInfo);
    //             const boost::numeric::ublas::bounded_matrix<double, MatrixSize, MatrixSize> LHS_contact_pair = this->CalculateLocalLHS<MatrixSize>( rThisMortarConditionMatrices, rDofData, LHS_SlaveElem_Contribution, EquationIdSlaveElem, PairIndex, rCurrentProcessInfo);
                
//                 // Debug
// //                 KRATOS_WATCH(LHS_contact_pair);
//                 LOG_MATRIX_PRETTY( LHS_contact_pair );
                
                // Contributions to stiffness matrix calculated on the reference config
                this->CalculateAndAddLHS( rLocalSystem, LHS_contact_pair, PairIndex );
            }

            // Assemble of the vector is required
            if ( rLocalSystem.CalculationFlags.Is( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR ) ||
                    rLocalSystem.CalculationFlags.Is( MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
            {
                // Calculate the local contribution
                const array_1d<double, MatrixSize> RHS_contact_pair = this->CalculateLocalRHS<MatrixSize>( rThisMortarConditionMatrices, rDofData, PairIndex, rCurrentProcessInfo);
    //             const array_1d<double, MatrixSize> RHS_contact_pair = this->CalculateLocalRHS<MatrixSize>( rThisMortarConditionMatrices, rDofData, RHS_SlaveElem_Contribution, EquationIdSlaveElem, PairIndex, rCurrentProcessInfo);
                
//                 // Debug
// //                 KRATOS_WATCH(RHS_contact_pair);
//                 LOG_VECTOR_PRETTY( RHS_contact_pair );
                
                // Contribution to previous step contact force and residuals vector
                this->CalculateAndAddRHS( rLocalSystem, RHS_contact_pair, PairIndex );
            }
        }
    }
    
    // Debug
//     VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();
// //     KRATOS_WATCH(rRightHandSideVector);
//     LOG_VECTOR_PRETTY( rRightHandSideVector );
    
//     // Debug
//     MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix( );  
// //     KRATOS_WATCH(rLeftHandSideMatrix);
//     LOG_MATRIX_PRETTY( rLeftHandSideMatrix );
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::InitializeGeneralVariables(
    GeneralVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo,
    const unsigned int rMasterElementIndex
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

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateAe(
    DofData& rDofData,
    GeneralVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We initilize the Ae components
    AeData rAeData;
    rAeData.Initialize();
    
    rDofData.InitializeAeComponents();
    
    double TotalWeight = 0.0;
    
    for (unsigned int PairIndex = 0; PairIndex < mPairSize; ++PairIndex)
    {
        // The normal of the master condition
        const array_1d<double, 3> MasterNormal = mThisMasterConditions[PairIndex]->GetValue(NORMAL);
        
        // Get the integration points
        const IntegrationPointsType IntegrationPointsSlave = mIntegrationPointsVector[PairIndex];
        
        // Initialize general variables for the current master element
        this->InitializeGeneralVariables( rVariables, rCurrentProcessInfo, PairIndex );
            
        // Calculating the proportion between the integrated area and segment area
        for ( unsigned int PointNumber = 0; PointNumber < IntegrationPointsSlave.size(); PointNumber++ )
        {            
            // Calculate the kinematic variables
            this->CalculateKinematics( rVariables, rDofData, MasterNormal, PointNumber, IntegrationPointsSlave );
            
            const double& IntegrationWeight = IntegrationPointsSlave[PointNumber].Weight();
            TotalWeight += IntegrationWeight;
            
            this->CalculateAeComponents(rVariables, rAeData, IntegrationWeight);
        }
    }
    
    if (TotalWeight > 0.0)
    {
        this->CalculateAe(rDofData, rAeData);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::InitializeDofData(
    DofData& rDofData,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // Slave element info
    rDofData.Initialize(GetGeometry());
    
    if (TTensor == 1)
    {
        for (unsigned int iNode = 0; iNode < NumNodes; iNode++)
        {
            const double Value = GetGeometry()[iNode].FastGetSolutionStepValue(TEMPERATURE);
            const double LM = GetGeometry()[iNode].FastGetSolutionStepValue(SCALAR_LAGRANGE_MULTIPLIER);
            rDofData.u1(iNode, 0) = Value;
            rDofData.LagrangeMultipliers(iNode, 0) = LM;
        }
    }
    else
    {
        for (unsigned int iNode = 0; iNode < NumNodes; iNode++)
        {
            const array_1d<double, 3> Value = GetGeometry()[iNode].FastGetSolutionStepValue(DISPLACEMENT);
            const array_1d<double, 3> LM = GetGeometry()[iNode].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
            for (unsigned int iDof = 0; iDof < TTensor; iDof++)
            {
                rDofData.u1(iNode, iDof) = Value[iDof];
                rDofData.LagrangeMultipliers(iNode, iDof) = LM[iDof];
            }
        }
    }
}

/*********************************COMPUTE KINEMATICS*********************************/
/************************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateKinematics( 
    GeneralVariables& rVariables,
    const DofData rDofData,
    const array_1d<double, 3> MasterNormal,
    const double& rPointNumber,
    const IntegrationPointsType& IntegrationPointsSlave
    )
{       
    /* LOCAL COORDINATES */
    const PointType& LocalPoint = IntegrationPointsSlave[rPointNumber].Coordinates();
       
    /// SLAVE CONDITION ///
    GetGeometry( ).ShapeFunctionsValues( rVariables.N_Slave, LocalPoint.Coordinates() );
    rVariables.Phi_LagrangeMultipliers = prod(rDofData.Ae, rVariables.N_Slave);
//     rVariables.Phi_LagrangeMultipliers = rVariables.N_Slave; // TODO: This could be needed in the future to be different than the standart shape functions 
    
    /* CALCULATE JACOBIAN AND JACOBIAN DETERMINANT */
    rVariables.DetJSlave = 1.0;
//     rVariables.DetJSlave = GetGeometry( ).DeterminantOfJacobian( LocalPoint );
    
    /// MASTER CONDITION ///
    this->MasterShapeFunctionValue( rVariables, MasterNormal, LocalPoint);
}
 
/***********************************************************************************/
/*************** METHODS TO CALCULATE THE CONTACT CONDITION MATRICES ***************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::MasterShapeFunctionValue(
    GeneralVariables& rVariables,
    const array_1d<double, 3> MasterNormal,
    const PointType& LocalPoint 
    )
{    
    GeometryType& MasterSegment = rVariables.GetMasterElement( );

    PointType ProjectedGPGlobal;
    const array_1d<double,3> GPNormal = MortarUtilities::GaussPointUnitNormal(rVariables.N_Slave, GetGeometry());
    
    GeometryType::CoordinatesArrayType SlaveGPGlobal;
    this->GetGeometry( ).GlobalCoordinates( SlaveGPGlobal, LocalPoint );
    MortarUtilities::FastProjectDirection( MasterSegment, SlaveGPGlobal, ProjectedGPGlobal, MasterNormal, -GPNormal ); // The opposite direction
    
    GeometryType::CoordinatesArrayType ProjectedGPLocal;
    
    MasterSegment.PointLocalCoordinates(ProjectedGPLocal, ProjectedGPGlobal.Coordinates( ) ) ;
    
    // SHAPE FUNCTIONS 
    MasterSegment.ShapeFunctionsValues( rVariables.N_Master, ProjectedGPLocal );         
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateMortarOperators(
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
    boost::numeric::ublas::bounded_matrix<double, NumNodes, NumNodes>& DOperator = rThisMortarConditionMatrices.DOperator;
    boost::numeric::ublas::bounded_matrix<double, NumNodes, NumNodes>& MOperator = rThisMortarConditionMatrices.MOperator;
    
    for (unsigned int i_slave = 0; i_slave < NumNodes; i_slave++)
    {
        for (unsigned int j_slave = 0; j_slave < NumNodes; j_slave++)
        {
            const double phi = Phi[i_slave];
            
            DOperator(i_slave, j_slave) += J_s * rIntegrationWeight * phi * N1[j_slave];
            MOperator(i_slave, j_slave) += J_s * rIntegrationWeight * phi * N2[j_slave];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateAeComponents(
    GeneralVariables& rVariables,
    AeData& rAeData,
    const double& rIntegrationWeight
    )
{
    /* DEFINITIONS */
    const VectorType N1 = rVariables.N_Slave;
    const double DetJ = rVariables.DetJSlave; 
     
    rAeData.De += rIntegrationWeight * (ComputeDe(N1, DetJ));
    rAeData.Me += rIntegrationWeight * DetJ * outer_prod(N1, N1);
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateAe(
    DofData& rDofData,
    AeData& rAeData
    )
{        
    double auxdet;
    const double Tolerance = std::numeric_limits<double>::epsilon();

    // We compute the norm
    const double NormMe = norm_frobenius(rAeData.Me);
    
    // Now we normalize the matrix
    const boost::numeric::ublas::bounded_matrix<double, NumNodes, NumNodes> NormalizedMe = rAeData.Me/NormMe;
    
    // We compute the normalized inverse
    const boost::numeric::ublas::bounded_matrix<double, NumNodes, NumNodes> NormalizedInvMe = MathUtils<double>::InvertMatrix<NumNodes>(NormalizedMe, auxdet, Tolerance); 
    
    // Now we compute the inverse
    const boost::numeric::ublas::bounded_matrix<double, NumNodes, NumNodes> InvMe = NormalizedInvMe/NormMe;

    noalias(rDofData.Ae) = prod(rAeData.De, InvMe);
}

/***************************** BEGIN AD REPLACEMENT ********************************/
/***********************************************************************************/

template< >
template< >
boost::numeric::ublas::bounded_matrix<double, 6, 6> MeshTyingMortarCondition<2,3,ScalarValue>::CalculateLocalLHS<6>(
    const MortarConditionMatrices& rMortarConditionMatrices,
    DofData& rDofData,
    const unsigned int rMasterElementIndex,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    boost::numeric::ublas::bounded_matrix<double, 6, 6> lhs;

    // We get the mortar operators
    const boost::numeric::ublas::bounded_matrix<double, 2, 2>& MOperator = rMortarConditionMatrices.MOperator;
    const boost::numeric::ublas::bounded_matrix<double, 2, 2>& DOperator = rMortarConditionMatrices.DOperator;

    const double clhs0 =     -MOperator(0,0);
    const double clhs1 =     -MOperator(1,0);
    const double clhs2 =     -MOperator(0,1);
    const double clhs3 =     -MOperator(1,1);

    lhs(0,0)=0;
    lhs(0,1)=0;
    lhs(0,2)=0;
    lhs(0,3)=0;
    lhs(0,4)=clhs0;
    lhs(0,5)=clhs1;
    lhs(1,0)=0;
    lhs(1,1)=0;
    lhs(1,2)=0;
    lhs(1,3)=0;
    lhs(1,4)=clhs2;
    lhs(1,5)=clhs3;
    lhs(2,0)=0;
    lhs(2,1)=0;
    lhs(2,2)=0;
    lhs(2,3)=0;
    lhs(2,4)=DOperator(0,0);
    lhs(2,5)=DOperator(1,0);
    lhs(3,0)=0;
    lhs(3,1)=0;
    lhs(3,2)=0;
    lhs(3,3)=0;
    lhs(3,4)=DOperator(0,1);
    lhs(3,5)=DOperator(1,1);
    lhs(4,0)=clhs0;
    lhs(4,1)=clhs2;
    lhs(4,2)=DOperator(0,0);
    lhs(4,3)=DOperator(0,1);
    lhs(4,4)=0;
    lhs(4,5)=0;
    lhs(5,0)=clhs1;
    lhs(5,1)=clhs3;
    lhs(5,2)=DOperator(1,0);
    lhs(5,3)=DOperator(1,1);
    lhs(5,4)=0;
    lhs(5,5)=0;


    return lhs;
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
boost::numeric::ublas::bounded_matrix<double, 12, 12> MeshTyingMortarCondition<2,3,Vector2DValue>::CalculateLocalLHS<12>(
    const MortarConditionMatrices& rMortarConditionMatrices,
    DofData& rDofData,
    const unsigned int rMasterElementIndex,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    boost::numeric::ublas::bounded_matrix<double, 12, 12> lhs;

    // We get the mortar operators
    const boost::numeric::ublas::bounded_matrix<double, 2, 2>& MOperator = rMortarConditionMatrices.MOperator;
    const boost::numeric::ublas::bounded_matrix<double, 2, 2>& DOperator = rMortarConditionMatrices.DOperator;

    const double clhs0 =     -MOperator(0,0);
    const double clhs1 =     -MOperator(1,0);
    const double clhs2 =     -MOperator(0,1);
    const double clhs3 =     -MOperator(1,1);

    lhs(0,0)=0;
    lhs(0,1)=0;
    lhs(0,2)=0;
    lhs(0,3)=0;
    lhs(0,4)=0;
    lhs(0,5)=0;
    lhs(0,6)=0;
    lhs(0,7)=0;
    lhs(0,8)=clhs0;
    lhs(0,9)=0;
    lhs(0,10)=clhs1;
    lhs(0,11)=0;
    lhs(1,0)=0;
    lhs(1,1)=0;
    lhs(1,2)=0;
    lhs(1,3)=0;
    lhs(1,4)=0;
    lhs(1,5)=0;
    lhs(1,6)=0;
    lhs(1,7)=0;
    lhs(1,8)=0;
    lhs(1,9)=clhs0;
    lhs(1,10)=0;
    lhs(1,11)=clhs1;
    lhs(2,0)=0;
    lhs(2,1)=0;
    lhs(2,2)=0;
    lhs(2,3)=0;
    lhs(2,4)=0;
    lhs(2,5)=0;
    lhs(2,6)=0;
    lhs(2,7)=0;
    lhs(2,8)=clhs2;
    lhs(2,9)=0;
    lhs(2,10)=clhs3;
    lhs(2,11)=0;
    lhs(3,0)=0;
    lhs(3,1)=0;
    lhs(3,2)=0;
    lhs(3,3)=0;
    lhs(3,4)=0;
    lhs(3,5)=0;
    lhs(3,6)=0;
    lhs(3,7)=0;
    lhs(3,8)=0;
    lhs(3,9)=clhs2;
    lhs(3,10)=0;
    lhs(3,11)=clhs3;
    lhs(4,0)=0;
    lhs(4,1)=0;
    lhs(4,2)=0;
    lhs(4,3)=0;
    lhs(4,4)=0;
    lhs(4,5)=0;
    lhs(4,6)=0;
    lhs(4,7)=0;
    lhs(4,8)=DOperator(0,0);
    lhs(4,9)=0;
    lhs(4,10)=DOperator(1,0);
    lhs(4,11)=0;
    lhs(5,0)=0;
    lhs(5,1)=0;
    lhs(5,2)=0;
    lhs(5,3)=0;
    lhs(5,4)=0;
    lhs(5,5)=0;
    lhs(5,6)=0;
    lhs(5,7)=0;
    lhs(5,8)=0;
    lhs(5,9)=DOperator(0,0);
    lhs(5,10)=0;
    lhs(5,11)=DOperator(1,0);
    lhs(6,0)=0;
    lhs(6,1)=0;
    lhs(6,2)=0;
    lhs(6,3)=0;
    lhs(6,4)=0;
    lhs(6,5)=0;
    lhs(6,6)=0;
    lhs(6,7)=0;
    lhs(6,8)=DOperator(0,1);
    lhs(6,9)=0;
    lhs(6,10)=DOperator(1,1);
    lhs(6,11)=0;
    lhs(7,0)=0;
    lhs(7,1)=0;
    lhs(7,2)=0;
    lhs(7,3)=0;
    lhs(7,4)=0;
    lhs(7,5)=0;
    lhs(7,6)=0;
    lhs(7,7)=0;
    lhs(7,8)=0;
    lhs(7,9)=DOperator(0,1);
    lhs(7,10)=0;
    lhs(7,11)=DOperator(1,1);
    lhs(8,0)=clhs0;
    lhs(8,1)=0;
    lhs(8,2)=clhs2;
    lhs(8,3)=0;
    lhs(8,4)=DOperator(0,0);
    lhs(8,5)=0;
    lhs(8,6)=DOperator(0,1);
    lhs(8,7)=0;
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=0;
    lhs(9,1)=clhs0;
    lhs(9,2)=0;
    lhs(9,3)=clhs2;
    lhs(9,4)=0;
    lhs(9,5)=DOperator(0,0);
    lhs(9,6)=0;
    lhs(9,7)=DOperator(0,1);
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=clhs1;
    lhs(10,1)=0;
    lhs(10,2)=clhs3;
    lhs(10,3)=0;
    lhs(10,4)=DOperator(1,0);
    lhs(10,5)=0;
    lhs(10,6)=DOperator(1,1);
    lhs(10,7)=0;
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=0;
    lhs(11,1)=clhs1;
    lhs(11,2)=0;
    lhs(11,3)=clhs3;
    lhs(11,4)=0;
    lhs(11,5)=DOperator(1,0);
    lhs(11,6)=0;
    lhs(11,7)=DOperator(1,1);
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;


    return lhs;
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
boost::numeric::ublas::bounded_matrix<double, 6, 6> MeshTyingMortarCondition<2,4,ScalarValue>::CalculateLocalLHS<6>(
    const MortarConditionMatrices& rMortarConditionMatrices,
    DofData& rDofData,
    const unsigned int rMasterElementIndex,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    boost::numeric::ublas::bounded_matrix<double, 6, 6> lhs;

    // We get the mortar operators
    const boost::numeric::ublas::bounded_matrix<double, 2, 2>& MOperator = rMortarConditionMatrices.MOperator;
    const boost::numeric::ublas::bounded_matrix<double, 2, 2>& DOperator = rMortarConditionMatrices.DOperator;

    const double clhs0 =     -MOperator(0,0);
    const double clhs1 =     -MOperator(1,0);
    const double clhs2 =     -MOperator(0,1);
    const double clhs3 =     -MOperator(1,1);

    lhs(0,0)=0;
    lhs(0,1)=0;
    lhs(0,2)=0;
    lhs(0,3)=0;
    lhs(0,4)=clhs0;
    lhs(0,5)=clhs1;
    lhs(1,0)=0;
    lhs(1,1)=0;
    lhs(1,2)=0;
    lhs(1,3)=0;
    lhs(1,4)=clhs2;
    lhs(1,5)=clhs3;
    lhs(2,0)=0;
    lhs(2,1)=0;
    lhs(2,2)=0;
    lhs(2,3)=0;
    lhs(2,4)=DOperator(0,0);
    lhs(2,5)=DOperator(1,0);
    lhs(3,0)=0;
    lhs(3,1)=0;
    lhs(3,2)=0;
    lhs(3,3)=0;
    lhs(3,4)=DOperator(0,1);
    lhs(3,5)=DOperator(1,1);
    lhs(4,0)=clhs0;
    lhs(4,1)=clhs2;
    lhs(4,2)=DOperator(0,0);
    lhs(4,3)=DOperator(0,1);
    lhs(4,4)=0;
    lhs(4,5)=0;
    lhs(5,0)=clhs1;
    lhs(5,1)=clhs3;
    lhs(5,2)=DOperator(1,0);
    lhs(5,3)=DOperator(1,1);
    lhs(5,4)=0;
    lhs(5,5)=0;


    return lhs;
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
boost::numeric::ublas::bounded_matrix<double, 12, 12> MeshTyingMortarCondition<2,4,Vector2DValue>::CalculateLocalLHS<12>(
    const MortarConditionMatrices& rMortarConditionMatrices,
    DofData& rDofData,
    const unsigned int rMasterElementIndex,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    boost::numeric::ublas::bounded_matrix<double, 12, 12> lhs;

    // We get the mortar operators
    const boost::numeric::ublas::bounded_matrix<double, 2, 2>& MOperator = rMortarConditionMatrices.MOperator;
    const boost::numeric::ublas::bounded_matrix<double, 2, 2>& DOperator = rMortarConditionMatrices.DOperator;

    const double clhs0 =     -MOperator(0,0);
    const double clhs1 =     -MOperator(1,0);
    const double clhs2 =     -MOperator(0,1);
    const double clhs3 =     -MOperator(1,1);

    lhs(0,0)=0;
    lhs(0,1)=0;
    lhs(0,2)=0;
    lhs(0,3)=0;
    lhs(0,4)=0;
    lhs(0,5)=0;
    lhs(0,6)=0;
    lhs(0,7)=0;
    lhs(0,8)=clhs0;
    lhs(0,9)=0;
    lhs(0,10)=clhs1;
    lhs(0,11)=0;
    lhs(1,0)=0;
    lhs(1,1)=0;
    lhs(1,2)=0;
    lhs(1,3)=0;
    lhs(1,4)=0;
    lhs(1,5)=0;
    lhs(1,6)=0;
    lhs(1,7)=0;
    lhs(1,8)=0;
    lhs(1,9)=clhs0;
    lhs(1,10)=0;
    lhs(1,11)=clhs1;
    lhs(2,0)=0;
    lhs(2,1)=0;
    lhs(2,2)=0;
    lhs(2,3)=0;
    lhs(2,4)=0;
    lhs(2,5)=0;
    lhs(2,6)=0;
    lhs(2,7)=0;
    lhs(2,8)=clhs2;
    lhs(2,9)=0;
    lhs(2,10)=clhs3;
    lhs(2,11)=0;
    lhs(3,0)=0;
    lhs(3,1)=0;
    lhs(3,2)=0;
    lhs(3,3)=0;
    lhs(3,4)=0;
    lhs(3,5)=0;
    lhs(3,6)=0;
    lhs(3,7)=0;
    lhs(3,8)=0;
    lhs(3,9)=clhs2;
    lhs(3,10)=0;
    lhs(3,11)=clhs3;
    lhs(4,0)=0;
    lhs(4,1)=0;
    lhs(4,2)=0;
    lhs(4,3)=0;
    lhs(4,4)=0;
    lhs(4,5)=0;
    lhs(4,6)=0;
    lhs(4,7)=0;
    lhs(4,8)=DOperator(0,0);
    lhs(4,9)=0;
    lhs(4,10)=DOperator(1,0);
    lhs(4,11)=0;
    lhs(5,0)=0;
    lhs(5,1)=0;
    lhs(5,2)=0;
    lhs(5,3)=0;
    lhs(5,4)=0;
    lhs(5,5)=0;
    lhs(5,6)=0;
    lhs(5,7)=0;
    lhs(5,8)=0;
    lhs(5,9)=DOperator(0,0);
    lhs(5,10)=0;
    lhs(5,11)=DOperator(1,0);
    lhs(6,0)=0;
    lhs(6,1)=0;
    lhs(6,2)=0;
    lhs(6,3)=0;
    lhs(6,4)=0;
    lhs(6,5)=0;
    lhs(6,6)=0;
    lhs(6,7)=0;
    lhs(6,8)=DOperator(0,1);
    lhs(6,9)=0;
    lhs(6,10)=DOperator(1,1);
    lhs(6,11)=0;
    lhs(7,0)=0;
    lhs(7,1)=0;
    lhs(7,2)=0;
    lhs(7,3)=0;
    lhs(7,4)=0;
    lhs(7,5)=0;
    lhs(7,6)=0;
    lhs(7,7)=0;
    lhs(7,8)=0;
    lhs(7,9)=DOperator(0,1);
    lhs(7,10)=0;
    lhs(7,11)=DOperator(1,1);
    lhs(8,0)=clhs0;
    lhs(8,1)=0;
    lhs(8,2)=clhs2;
    lhs(8,3)=0;
    lhs(8,4)=DOperator(0,0);
    lhs(8,5)=0;
    lhs(8,6)=DOperator(0,1);
    lhs(8,7)=0;
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=0;
    lhs(9,1)=clhs0;
    lhs(9,2)=0;
    lhs(9,3)=clhs2;
    lhs(9,4)=0;
    lhs(9,5)=DOperator(0,0);
    lhs(9,6)=0;
    lhs(9,7)=DOperator(0,1);
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=clhs1;
    lhs(10,1)=0;
    lhs(10,2)=clhs3;
    lhs(10,3)=0;
    lhs(10,4)=DOperator(1,0);
    lhs(10,5)=0;
    lhs(10,6)=DOperator(1,1);
    lhs(10,7)=0;
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=0;
    lhs(11,1)=clhs1;
    lhs(11,2)=0;
    lhs(11,3)=clhs3;
    lhs(11,4)=0;
    lhs(11,5)=DOperator(1,0);
    lhs(11,6)=0;
    lhs(11,7)=DOperator(1,1);
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;


    return lhs;
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
boost::numeric::ublas::bounded_matrix<double, 9, 9> MeshTyingMortarCondition<3,4,ScalarValue>::CalculateLocalLHS<9>(
    const MortarConditionMatrices& rMortarConditionMatrices,
    DofData& rDofData,
    const unsigned int rMasterElementIndex,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    boost::numeric::ublas::bounded_matrix<double, 9, 9> lhs;

    // We get the mortar operators
    const boost::numeric::ublas::bounded_matrix<double, 3, 3>& MOperator = rMortarConditionMatrices.MOperator;
    const boost::numeric::ublas::bounded_matrix<double, 3, 3>& DOperator = rMortarConditionMatrices.DOperator;

    const double clhs0 =     -MOperator(0,0);
    const double clhs1 =     -MOperator(1,0);
    const double clhs2 =     -MOperator(2,0);
    const double clhs3 =     -MOperator(0,1);
    const double clhs4 =     -MOperator(1,1);
    const double clhs5 =     -MOperator(2,1);
    const double clhs6 =     -MOperator(0,2);
    const double clhs7 =     -MOperator(1,2);
    const double clhs8 =     -MOperator(2,2);

    lhs(0,0)=0;
    lhs(0,1)=0;
    lhs(0,2)=0;
    lhs(0,3)=0;
    lhs(0,4)=0;
    lhs(0,5)=0;
    lhs(0,6)=clhs0;
    lhs(0,7)=clhs1;
    lhs(0,8)=clhs2;
    lhs(1,0)=0;
    lhs(1,1)=0;
    lhs(1,2)=0;
    lhs(1,3)=0;
    lhs(1,4)=0;
    lhs(1,5)=0;
    lhs(1,6)=clhs3;
    lhs(1,7)=clhs4;
    lhs(1,8)=clhs5;
    lhs(2,0)=0;
    lhs(2,1)=0;
    lhs(2,2)=0;
    lhs(2,3)=0;
    lhs(2,4)=0;
    lhs(2,5)=0;
    lhs(2,6)=clhs6;
    lhs(2,7)=clhs7;
    lhs(2,8)=clhs8;
    lhs(3,0)=0;
    lhs(3,1)=0;
    lhs(3,2)=0;
    lhs(3,3)=0;
    lhs(3,4)=0;
    lhs(3,5)=0;
    lhs(3,6)=DOperator(0,0);
    lhs(3,7)=DOperator(1,0);
    lhs(3,8)=DOperator(2,0);
    lhs(4,0)=0;
    lhs(4,1)=0;
    lhs(4,2)=0;
    lhs(4,3)=0;
    lhs(4,4)=0;
    lhs(4,5)=0;
    lhs(4,6)=DOperator(0,1);
    lhs(4,7)=DOperator(1,1);
    lhs(4,8)=DOperator(2,1);
    lhs(5,0)=0;
    lhs(5,1)=0;
    lhs(5,2)=0;
    lhs(5,3)=0;
    lhs(5,4)=0;
    lhs(5,5)=0;
    lhs(5,6)=DOperator(0,2);
    lhs(5,7)=DOperator(1,2);
    lhs(5,8)=DOperator(2,2);
    lhs(6,0)=clhs0;
    lhs(6,1)=clhs3;
    lhs(6,2)=clhs6;
    lhs(6,3)=DOperator(0,0);
    lhs(6,4)=DOperator(0,1);
    lhs(6,5)=DOperator(0,2);
    lhs(6,6)=0;
    lhs(6,7)=0;
    lhs(6,8)=0;
    lhs(7,0)=clhs1;
    lhs(7,1)=clhs4;
    lhs(7,2)=clhs7;
    lhs(7,3)=DOperator(1,0);
    lhs(7,4)=DOperator(1,1);
    lhs(7,5)=DOperator(1,2);
    lhs(7,6)=0;
    lhs(7,7)=0;
    lhs(7,8)=0;
    lhs(8,0)=clhs2;
    lhs(8,1)=clhs5;
    lhs(8,2)=clhs8;
    lhs(8,3)=DOperator(2,0);
    lhs(8,4)=DOperator(2,1);
    lhs(8,5)=DOperator(2,2);
    lhs(8,6)=0;
    lhs(8,7)=0;
    lhs(8,8)=0;


    return lhs;
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
boost::numeric::ublas::bounded_matrix<double, 27, 27> MeshTyingMortarCondition<3,4,Vector3DValue>::CalculateLocalLHS<27>(
    const MortarConditionMatrices& rMortarConditionMatrices,
    DofData& rDofData,
    const unsigned int rMasterElementIndex,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    boost::numeric::ublas::bounded_matrix<double, 27, 27> lhs;

    // We get the mortar operators
    const boost::numeric::ublas::bounded_matrix<double, 3, 3>& MOperator = rMortarConditionMatrices.MOperator;
    const boost::numeric::ublas::bounded_matrix<double, 3, 3>& DOperator = rMortarConditionMatrices.DOperator;

    const double clhs0 =     -MOperator(0,0);
    const double clhs1 =     -MOperator(1,0);
    const double clhs2 =     -MOperator(2,0);
    const double clhs3 =     -MOperator(0,1);
    const double clhs4 =     -MOperator(1,1);
    const double clhs5 =     -MOperator(2,1);
    const double clhs6 =     -MOperator(0,2);
    const double clhs7 =     -MOperator(1,2);
    const double clhs8 =     -MOperator(2,2);

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
    lhs(0,10)=0;
    lhs(0,11)=0;
    lhs(0,12)=0;
    lhs(0,13)=0;
    lhs(0,14)=0;
    lhs(0,15)=0;
    lhs(0,16)=0;
    lhs(0,17)=0;
    lhs(0,18)=clhs0;
    lhs(0,19)=0;
    lhs(0,20)=0;
    lhs(0,21)=clhs1;
    lhs(0,22)=0;
    lhs(0,23)=0;
    lhs(0,24)=clhs2;
    lhs(0,25)=0;
    lhs(0,26)=0;
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
    lhs(1,10)=0;
    lhs(1,11)=0;
    lhs(1,12)=0;
    lhs(1,13)=0;
    lhs(1,14)=0;
    lhs(1,15)=0;
    lhs(1,16)=0;
    lhs(1,17)=0;
    lhs(1,18)=0;
    lhs(1,19)=clhs0;
    lhs(1,20)=0;
    lhs(1,21)=0;
    lhs(1,22)=clhs1;
    lhs(1,23)=0;
    lhs(1,24)=0;
    lhs(1,25)=clhs2;
    lhs(1,26)=0;
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
    lhs(2,10)=0;
    lhs(2,11)=0;
    lhs(2,12)=0;
    lhs(2,13)=0;
    lhs(2,14)=0;
    lhs(2,15)=0;
    lhs(2,16)=0;
    lhs(2,17)=0;
    lhs(2,18)=0;
    lhs(2,19)=0;
    lhs(2,20)=clhs0;
    lhs(2,21)=0;
    lhs(2,22)=0;
    lhs(2,23)=clhs1;
    lhs(2,24)=0;
    lhs(2,25)=0;
    lhs(2,26)=clhs2;
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
    lhs(3,10)=0;
    lhs(3,11)=0;
    lhs(3,12)=0;
    lhs(3,13)=0;
    lhs(3,14)=0;
    lhs(3,15)=0;
    lhs(3,16)=0;
    lhs(3,17)=0;
    lhs(3,18)=clhs3;
    lhs(3,19)=0;
    lhs(3,20)=0;
    lhs(3,21)=clhs4;
    lhs(3,22)=0;
    lhs(3,23)=0;
    lhs(3,24)=clhs5;
    lhs(3,25)=0;
    lhs(3,26)=0;
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
    lhs(4,10)=0;
    lhs(4,11)=0;
    lhs(4,12)=0;
    lhs(4,13)=0;
    lhs(4,14)=0;
    lhs(4,15)=0;
    lhs(4,16)=0;
    lhs(4,17)=0;
    lhs(4,18)=0;
    lhs(4,19)=clhs3;
    lhs(4,20)=0;
    lhs(4,21)=0;
    lhs(4,22)=clhs4;
    lhs(4,23)=0;
    lhs(4,24)=0;
    lhs(4,25)=clhs5;
    lhs(4,26)=0;
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
    lhs(5,10)=0;
    lhs(5,11)=0;
    lhs(5,12)=0;
    lhs(5,13)=0;
    lhs(5,14)=0;
    lhs(5,15)=0;
    lhs(5,16)=0;
    lhs(5,17)=0;
    lhs(5,18)=0;
    lhs(5,19)=0;
    lhs(5,20)=clhs3;
    lhs(5,21)=0;
    lhs(5,22)=0;
    lhs(5,23)=clhs4;
    lhs(5,24)=0;
    lhs(5,25)=0;
    lhs(5,26)=clhs5;
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
    lhs(6,10)=0;
    lhs(6,11)=0;
    lhs(6,12)=0;
    lhs(6,13)=0;
    lhs(6,14)=0;
    lhs(6,15)=0;
    lhs(6,16)=0;
    lhs(6,17)=0;
    lhs(6,18)=clhs6;
    lhs(6,19)=0;
    lhs(6,20)=0;
    lhs(6,21)=clhs7;
    lhs(6,22)=0;
    lhs(6,23)=0;
    lhs(6,24)=clhs8;
    lhs(6,25)=0;
    lhs(6,26)=0;
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
    lhs(7,10)=0;
    lhs(7,11)=0;
    lhs(7,12)=0;
    lhs(7,13)=0;
    lhs(7,14)=0;
    lhs(7,15)=0;
    lhs(7,16)=0;
    lhs(7,17)=0;
    lhs(7,18)=0;
    lhs(7,19)=clhs6;
    lhs(7,20)=0;
    lhs(7,21)=0;
    lhs(7,22)=clhs7;
    lhs(7,23)=0;
    lhs(7,24)=0;
    lhs(7,25)=clhs8;
    lhs(7,26)=0;
    lhs(8,0)=0;
    lhs(8,1)=0;
    lhs(8,2)=0;
    lhs(8,3)=0;
    lhs(8,4)=0;
    lhs(8,5)=0;
    lhs(8,6)=0;
    lhs(8,7)=0;
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(8,12)=0;
    lhs(8,13)=0;
    lhs(8,14)=0;
    lhs(8,15)=0;
    lhs(8,16)=0;
    lhs(8,17)=0;
    lhs(8,18)=0;
    lhs(8,19)=0;
    lhs(8,20)=clhs6;
    lhs(8,21)=0;
    lhs(8,22)=0;
    lhs(8,23)=clhs7;
    lhs(8,24)=0;
    lhs(8,25)=0;
    lhs(8,26)=clhs8;
    lhs(9,0)=0;
    lhs(9,1)=0;
    lhs(9,2)=0;
    lhs(9,3)=0;
    lhs(9,4)=0;
    lhs(9,5)=0;
    lhs(9,6)=0;
    lhs(9,7)=0;
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(9,12)=0;
    lhs(9,13)=0;
    lhs(9,14)=0;
    lhs(9,15)=0;
    lhs(9,16)=0;
    lhs(9,17)=0;
    lhs(9,18)=DOperator(0,0);
    lhs(9,19)=0;
    lhs(9,20)=0;
    lhs(9,21)=DOperator(1,0);
    lhs(9,22)=0;
    lhs(9,23)=0;
    lhs(9,24)=DOperator(2,0);
    lhs(9,25)=0;
    lhs(9,26)=0;
    lhs(10,0)=0;
    lhs(10,1)=0;
    lhs(10,2)=0;
    lhs(10,3)=0;
    lhs(10,4)=0;
    lhs(10,5)=0;
    lhs(10,6)=0;
    lhs(10,7)=0;
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(10,12)=0;
    lhs(10,13)=0;
    lhs(10,14)=0;
    lhs(10,15)=0;
    lhs(10,16)=0;
    lhs(10,17)=0;
    lhs(10,18)=0;
    lhs(10,19)=DOperator(0,0);
    lhs(10,20)=0;
    lhs(10,21)=0;
    lhs(10,22)=DOperator(1,0);
    lhs(10,23)=0;
    lhs(10,24)=0;
    lhs(10,25)=DOperator(2,0);
    lhs(10,26)=0;
    lhs(11,0)=0;
    lhs(11,1)=0;
    lhs(11,2)=0;
    lhs(11,3)=0;
    lhs(11,4)=0;
    lhs(11,5)=0;
    lhs(11,6)=0;
    lhs(11,7)=0;
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;
    lhs(11,12)=0;
    lhs(11,13)=0;
    lhs(11,14)=0;
    lhs(11,15)=0;
    lhs(11,16)=0;
    lhs(11,17)=0;
    lhs(11,18)=0;
    lhs(11,19)=0;
    lhs(11,20)=DOperator(0,0);
    lhs(11,21)=0;
    lhs(11,22)=0;
    lhs(11,23)=DOperator(1,0);
    lhs(11,24)=0;
    lhs(11,25)=0;
    lhs(11,26)=DOperator(2,0);
    lhs(12,0)=0;
    lhs(12,1)=0;
    lhs(12,2)=0;
    lhs(12,3)=0;
    lhs(12,4)=0;
    lhs(12,5)=0;
    lhs(12,6)=0;
    lhs(12,7)=0;
    lhs(12,8)=0;
    lhs(12,9)=0;
    lhs(12,10)=0;
    lhs(12,11)=0;
    lhs(12,12)=0;
    lhs(12,13)=0;
    lhs(12,14)=0;
    lhs(12,15)=0;
    lhs(12,16)=0;
    lhs(12,17)=0;
    lhs(12,18)=DOperator(0,1);
    lhs(12,19)=0;
    lhs(12,20)=0;
    lhs(12,21)=DOperator(1,1);
    lhs(12,22)=0;
    lhs(12,23)=0;
    lhs(12,24)=DOperator(2,1);
    lhs(12,25)=0;
    lhs(12,26)=0;
    lhs(13,0)=0;
    lhs(13,1)=0;
    lhs(13,2)=0;
    lhs(13,3)=0;
    lhs(13,4)=0;
    lhs(13,5)=0;
    lhs(13,6)=0;
    lhs(13,7)=0;
    lhs(13,8)=0;
    lhs(13,9)=0;
    lhs(13,10)=0;
    lhs(13,11)=0;
    lhs(13,12)=0;
    lhs(13,13)=0;
    lhs(13,14)=0;
    lhs(13,15)=0;
    lhs(13,16)=0;
    lhs(13,17)=0;
    lhs(13,18)=0;
    lhs(13,19)=DOperator(0,1);
    lhs(13,20)=0;
    lhs(13,21)=0;
    lhs(13,22)=DOperator(1,1);
    lhs(13,23)=0;
    lhs(13,24)=0;
    lhs(13,25)=DOperator(2,1);
    lhs(13,26)=0;
    lhs(14,0)=0;
    lhs(14,1)=0;
    lhs(14,2)=0;
    lhs(14,3)=0;
    lhs(14,4)=0;
    lhs(14,5)=0;
    lhs(14,6)=0;
    lhs(14,7)=0;
    lhs(14,8)=0;
    lhs(14,9)=0;
    lhs(14,10)=0;
    lhs(14,11)=0;
    lhs(14,12)=0;
    lhs(14,13)=0;
    lhs(14,14)=0;
    lhs(14,15)=0;
    lhs(14,16)=0;
    lhs(14,17)=0;
    lhs(14,18)=0;
    lhs(14,19)=0;
    lhs(14,20)=DOperator(0,1);
    lhs(14,21)=0;
    lhs(14,22)=0;
    lhs(14,23)=DOperator(1,1);
    lhs(14,24)=0;
    lhs(14,25)=0;
    lhs(14,26)=DOperator(2,1);
    lhs(15,0)=0;
    lhs(15,1)=0;
    lhs(15,2)=0;
    lhs(15,3)=0;
    lhs(15,4)=0;
    lhs(15,5)=0;
    lhs(15,6)=0;
    lhs(15,7)=0;
    lhs(15,8)=0;
    lhs(15,9)=0;
    lhs(15,10)=0;
    lhs(15,11)=0;
    lhs(15,12)=0;
    lhs(15,13)=0;
    lhs(15,14)=0;
    lhs(15,15)=0;
    lhs(15,16)=0;
    lhs(15,17)=0;
    lhs(15,18)=DOperator(0,2);
    lhs(15,19)=0;
    lhs(15,20)=0;
    lhs(15,21)=DOperator(1,2);
    lhs(15,22)=0;
    lhs(15,23)=0;
    lhs(15,24)=DOperator(2,2);
    lhs(15,25)=0;
    lhs(15,26)=0;
    lhs(16,0)=0;
    lhs(16,1)=0;
    lhs(16,2)=0;
    lhs(16,3)=0;
    lhs(16,4)=0;
    lhs(16,5)=0;
    lhs(16,6)=0;
    lhs(16,7)=0;
    lhs(16,8)=0;
    lhs(16,9)=0;
    lhs(16,10)=0;
    lhs(16,11)=0;
    lhs(16,12)=0;
    lhs(16,13)=0;
    lhs(16,14)=0;
    lhs(16,15)=0;
    lhs(16,16)=0;
    lhs(16,17)=0;
    lhs(16,18)=0;
    lhs(16,19)=DOperator(0,2);
    lhs(16,20)=0;
    lhs(16,21)=0;
    lhs(16,22)=DOperator(1,2);
    lhs(16,23)=0;
    lhs(16,24)=0;
    lhs(16,25)=DOperator(2,2);
    lhs(16,26)=0;
    lhs(17,0)=0;
    lhs(17,1)=0;
    lhs(17,2)=0;
    lhs(17,3)=0;
    lhs(17,4)=0;
    lhs(17,5)=0;
    lhs(17,6)=0;
    lhs(17,7)=0;
    lhs(17,8)=0;
    lhs(17,9)=0;
    lhs(17,10)=0;
    lhs(17,11)=0;
    lhs(17,12)=0;
    lhs(17,13)=0;
    lhs(17,14)=0;
    lhs(17,15)=0;
    lhs(17,16)=0;
    lhs(17,17)=0;
    lhs(17,18)=0;
    lhs(17,19)=0;
    lhs(17,20)=DOperator(0,2);
    lhs(17,21)=0;
    lhs(17,22)=0;
    lhs(17,23)=DOperator(1,2);
    lhs(17,24)=0;
    lhs(17,25)=0;
    lhs(17,26)=DOperator(2,2);
    lhs(18,0)=clhs0;
    lhs(18,1)=0;
    lhs(18,2)=0;
    lhs(18,3)=clhs3;
    lhs(18,4)=0;
    lhs(18,5)=0;
    lhs(18,6)=clhs6;
    lhs(18,7)=0;
    lhs(18,8)=0;
    lhs(18,9)=DOperator(0,0);
    lhs(18,10)=0;
    lhs(18,11)=0;
    lhs(18,12)=DOperator(0,1);
    lhs(18,13)=0;
    lhs(18,14)=0;
    lhs(18,15)=DOperator(0,2);
    lhs(18,16)=0;
    lhs(18,17)=0;
    lhs(18,18)=0;
    lhs(18,19)=0;
    lhs(18,20)=0;
    lhs(18,21)=0;
    lhs(18,22)=0;
    lhs(18,23)=0;
    lhs(18,24)=0;
    lhs(18,25)=0;
    lhs(18,26)=0;
    lhs(19,0)=0;
    lhs(19,1)=clhs0;
    lhs(19,2)=0;
    lhs(19,3)=0;
    lhs(19,4)=clhs3;
    lhs(19,5)=0;
    lhs(19,6)=0;
    lhs(19,7)=clhs6;
    lhs(19,8)=0;
    lhs(19,9)=0;
    lhs(19,10)=DOperator(0,0);
    lhs(19,11)=0;
    lhs(19,12)=0;
    lhs(19,13)=DOperator(0,1);
    lhs(19,14)=0;
    lhs(19,15)=0;
    lhs(19,16)=DOperator(0,2);
    lhs(19,17)=0;
    lhs(19,18)=0;
    lhs(19,19)=0;
    lhs(19,20)=0;
    lhs(19,21)=0;
    lhs(19,22)=0;
    lhs(19,23)=0;
    lhs(19,24)=0;
    lhs(19,25)=0;
    lhs(19,26)=0;
    lhs(20,0)=0;
    lhs(20,1)=0;
    lhs(20,2)=clhs0;
    lhs(20,3)=0;
    lhs(20,4)=0;
    lhs(20,5)=clhs3;
    lhs(20,6)=0;
    lhs(20,7)=0;
    lhs(20,8)=clhs6;
    lhs(20,9)=0;
    lhs(20,10)=0;
    lhs(20,11)=DOperator(0,0);
    lhs(20,12)=0;
    lhs(20,13)=0;
    lhs(20,14)=DOperator(0,1);
    lhs(20,15)=0;
    lhs(20,16)=0;
    lhs(20,17)=DOperator(0,2);
    lhs(20,18)=0;
    lhs(20,19)=0;
    lhs(20,20)=0;
    lhs(20,21)=0;
    lhs(20,22)=0;
    lhs(20,23)=0;
    lhs(20,24)=0;
    lhs(20,25)=0;
    lhs(20,26)=0;
    lhs(21,0)=clhs1;
    lhs(21,1)=0;
    lhs(21,2)=0;
    lhs(21,3)=clhs4;
    lhs(21,4)=0;
    lhs(21,5)=0;
    lhs(21,6)=clhs7;
    lhs(21,7)=0;
    lhs(21,8)=0;
    lhs(21,9)=DOperator(1,0);
    lhs(21,10)=0;
    lhs(21,11)=0;
    lhs(21,12)=DOperator(1,1);
    lhs(21,13)=0;
    lhs(21,14)=0;
    lhs(21,15)=DOperator(1,2);
    lhs(21,16)=0;
    lhs(21,17)=0;
    lhs(21,18)=0;
    lhs(21,19)=0;
    lhs(21,20)=0;
    lhs(21,21)=0;
    lhs(21,22)=0;
    lhs(21,23)=0;
    lhs(21,24)=0;
    lhs(21,25)=0;
    lhs(21,26)=0;
    lhs(22,0)=0;
    lhs(22,1)=clhs1;
    lhs(22,2)=0;
    lhs(22,3)=0;
    lhs(22,4)=clhs4;
    lhs(22,5)=0;
    lhs(22,6)=0;
    lhs(22,7)=clhs7;
    lhs(22,8)=0;
    lhs(22,9)=0;
    lhs(22,10)=DOperator(1,0);
    lhs(22,11)=0;
    lhs(22,12)=0;
    lhs(22,13)=DOperator(1,1);
    lhs(22,14)=0;
    lhs(22,15)=0;
    lhs(22,16)=DOperator(1,2);
    lhs(22,17)=0;
    lhs(22,18)=0;
    lhs(22,19)=0;
    lhs(22,20)=0;
    lhs(22,21)=0;
    lhs(22,22)=0;
    lhs(22,23)=0;
    lhs(22,24)=0;
    lhs(22,25)=0;
    lhs(22,26)=0;
    lhs(23,0)=0;
    lhs(23,1)=0;
    lhs(23,2)=clhs1;
    lhs(23,3)=0;
    lhs(23,4)=0;
    lhs(23,5)=clhs4;
    lhs(23,6)=0;
    lhs(23,7)=0;
    lhs(23,8)=clhs7;
    lhs(23,9)=0;
    lhs(23,10)=0;
    lhs(23,11)=DOperator(1,0);
    lhs(23,12)=0;
    lhs(23,13)=0;
    lhs(23,14)=DOperator(1,1);
    lhs(23,15)=0;
    lhs(23,16)=0;
    lhs(23,17)=DOperator(1,2);
    lhs(23,18)=0;
    lhs(23,19)=0;
    lhs(23,20)=0;
    lhs(23,21)=0;
    lhs(23,22)=0;
    lhs(23,23)=0;
    lhs(23,24)=0;
    lhs(23,25)=0;
    lhs(23,26)=0;
    lhs(24,0)=clhs2;
    lhs(24,1)=0;
    lhs(24,2)=0;
    lhs(24,3)=clhs5;
    lhs(24,4)=0;
    lhs(24,5)=0;
    lhs(24,6)=clhs8;
    lhs(24,7)=0;
    lhs(24,8)=0;
    lhs(24,9)=DOperator(2,0);
    lhs(24,10)=0;
    lhs(24,11)=0;
    lhs(24,12)=DOperator(2,1);
    lhs(24,13)=0;
    lhs(24,14)=0;
    lhs(24,15)=DOperator(2,2);
    lhs(24,16)=0;
    lhs(24,17)=0;
    lhs(24,18)=0;
    lhs(24,19)=0;
    lhs(24,20)=0;
    lhs(24,21)=0;
    lhs(24,22)=0;
    lhs(24,23)=0;
    lhs(24,24)=0;
    lhs(24,25)=0;
    lhs(24,26)=0;
    lhs(25,0)=0;
    lhs(25,1)=clhs2;
    lhs(25,2)=0;
    lhs(25,3)=0;
    lhs(25,4)=clhs5;
    lhs(25,5)=0;
    lhs(25,6)=0;
    lhs(25,7)=clhs8;
    lhs(25,8)=0;
    lhs(25,9)=0;
    lhs(25,10)=DOperator(2,0);
    lhs(25,11)=0;
    lhs(25,12)=0;
    lhs(25,13)=DOperator(2,1);
    lhs(25,14)=0;
    lhs(25,15)=0;
    lhs(25,16)=DOperator(2,2);
    lhs(25,17)=0;
    lhs(25,18)=0;
    lhs(25,19)=0;
    lhs(25,20)=0;
    lhs(25,21)=0;
    lhs(25,22)=0;
    lhs(25,23)=0;
    lhs(25,24)=0;
    lhs(25,25)=0;
    lhs(25,26)=0;
    lhs(26,0)=0;
    lhs(26,1)=0;
    lhs(26,2)=clhs2;
    lhs(26,3)=0;
    lhs(26,4)=0;
    lhs(26,5)=clhs5;
    lhs(26,6)=0;
    lhs(26,7)=0;
    lhs(26,8)=clhs8;
    lhs(26,9)=0;
    lhs(26,10)=0;
    lhs(26,11)=DOperator(2,0);
    lhs(26,12)=0;
    lhs(26,13)=0;
    lhs(26,14)=DOperator(2,1);
    lhs(26,15)=0;
    lhs(26,16)=0;
    lhs(26,17)=DOperator(2,2);
    lhs(26,18)=0;
    lhs(26,19)=0;
    lhs(26,20)=0;
    lhs(26,21)=0;
    lhs(26,22)=0;
    lhs(26,23)=0;
    lhs(26,24)=0;
    lhs(26,25)=0;
    lhs(26,26)=0;


    return lhs;
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
boost::numeric::ublas::bounded_matrix<double, 12, 12> MeshTyingMortarCondition<3,8,ScalarValue>::CalculateLocalLHS<12>(
    const MortarConditionMatrices& rMortarConditionMatrices,
    DofData& rDofData,
    const unsigned int rMasterElementIndex,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    boost::numeric::ublas::bounded_matrix<double, 12, 12> lhs;

    // We get the mortar operators
    const boost::numeric::ublas::bounded_matrix<double, 4, 4>& MOperator = rMortarConditionMatrices.MOperator;
    const boost::numeric::ublas::bounded_matrix<double, 4, 4>& DOperator = rMortarConditionMatrices.DOperator;

    const double clhs0 =     -MOperator(0,0);
    const double clhs1 =     -MOperator(1,0);
    const double clhs2 =     -MOperator(2,0);
    const double clhs3 =     -MOperator(3,0);
    const double clhs4 =     -MOperator(0,1);
    const double clhs5 =     -MOperator(1,1);
    const double clhs6 =     -MOperator(2,1);
    const double clhs7 =     -MOperator(3,1);
    const double clhs8 =     -MOperator(0,2);
    const double clhs9 =     -MOperator(1,2);
    const double clhs10 =     -MOperator(2,2);
    const double clhs11 =     -MOperator(3,2);
    const double clhs12 =     -MOperator(0,3);
    const double clhs13 =     -MOperator(1,3);
    const double clhs14 =     -MOperator(2,3);
    const double clhs15 =     -MOperator(3,3);

    lhs(0,0)=0;
    lhs(0,1)=0;
    lhs(0,2)=0;
    lhs(0,3)=0;
    lhs(0,4)=0;
    lhs(0,5)=0;
    lhs(0,6)=0;
    lhs(0,7)=0;
    lhs(0,8)=clhs0;
    lhs(0,9)=clhs1;
    lhs(0,10)=clhs2;
    lhs(0,11)=clhs3;
    lhs(1,0)=0;
    lhs(1,1)=0;
    lhs(1,2)=0;
    lhs(1,3)=0;
    lhs(1,4)=0;
    lhs(1,5)=0;
    lhs(1,6)=0;
    lhs(1,7)=0;
    lhs(1,8)=clhs4;
    lhs(1,9)=clhs5;
    lhs(1,10)=clhs6;
    lhs(1,11)=clhs7;
    lhs(2,0)=0;
    lhs(2,1)=0;
    lhs(2,2)=0;
    lhs(2,3)=0;
    lhs(2,4)=0;
    lhs(2,5)=0;
    lhs(2,6)=0;
    lhs(2,7)=0;
    lhs(2,8)=clhs8;
    lhs(2,9)=clhs9;
    lhs(2,10)=clhs10;
    lhs(2,11)=clhs11;
    lhs(3,0)=0;
    lhs(3,1)=0;
    lhs(3,2)=0;
    lhs(3,3)=0;
    lhs(3,4)=0;
    lhs(3,5)=0;
    lhs(3,6)=0;
    lhs(3,7)=0;
    lhs(3,8)=clhs12;
    lhs(3,9)=clhs13;
    lhs(3,10)=clhs14;
    lhs(3,11)=clhs15;
    lhs(4,0)=0;
    lhs(4,1)=0;
    lhs(4,2)=0;
    lhs(4,3)=0;
    lhs(4,4)=0;
    lhs(4,5)=0;
    lhs(4,6)=0;
    lhs(4,7)=0;
    lhs(4,8)=DOperator(0,0);
    lhs(4,9)=DOperator(1,0);
    lhs(4,10)=DOperator(2,0);
    lhs(4,11)=DOperator(3,0);
    lhs(5,0)=0;
    lhs(5,1)=0;
    lhs(5,2)=0;
    lhs(5,3)=0;
    lhs(5,4)=0;
    lhs(5,5)=0;
    lhs(5,6)=0;
    lhs(5,7)=0;
    lhs(5,8)=DOperator(0,1);
    lhs(5,9)=DOperator(1,1);
    lhs(5,10)=DOperator(2,1);
    lhs(5,11)=DOperator(3,1);
    lhs(6,0)=0;
    lhs(6,1)=0;
    lhs(6,2)=0;
    lhs(6,3)=0;
    lhs(6,4)=0;
    lhs(6,5)=0;
    lhs(6,6)=0;
    lhs(6,7)=0;
    lhs(6,8)=DOperator(0,2);
    lhs(6,9)=DOperator(1,2);
    lhs(6,10)=DOperator(2,2);
    lhs(6,11)=DOperator(3,2);
    lhs(7,0)=0;
    lhs(7,1)=0;
    lhs(7,2)=0;
    lhs(7,3)=0;
    lhs(7,4)=0;
    lhs(7,5)=0;
    lhs(7,6)=0;
    lhs(7,7)=0;
    lhs(7,8)=DOperator(0,3);
    lhs(7,9)=DOperator(1,3);
    lhs(7,10)=DOperator(2,3);
    lhs(7,11)=DOperator(3,3);
    lhs(8,0)=clhs0;
    lhs(8,1)=clhs4;
    lhs(8,2)=clhs8;
    lhs(8,3)=clhs12;
    lhs(8,4)=DOperator(0,0);
    lhs(8,5)=DOperator(0,1);
    lhs(8,6)=DOperator(0,2);
    lhs(8,7)=DOperator(0,3);
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=clhs1;
    lhs(9,1)=clhs5;
    lhs(9,2)=clhs9;
    lhs(9,3)=clhs13;
    lhs(9,4)=DOperator(1,0);
    lhs(9,5)=DOperator(1,1);
    lhs(9,6)=DOperator(1,2);
    lhs(9,7)=DOperator(1,3);
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=clhs2;
    lhs(10,1)=clhs6;
    lhs(10,2)=clhs10;
    lhs(10,3)=clhs14;
    lhs(10,4)=DOperator(2,0);
    lhs(10,5)=DOperator(2,1);
    lhs(10,6)=DOperator(2,2);
    lhs(10,7)=DOperator(2,3);
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=clhs3;
    lhs(11,1)=clhs7;
    lhs(11,2)=clhs11;
    lhs(11,3)=clhs15;
    lhs(11,4)=DOperator(3,0);
    lhs(11,5)=DOperator(3,1);
    lhs(11,6)=DOperator(3,2);
    lhs(11,7)=DOperator(3,3);
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;


    return lhs;
}

/***********************************************************************************/
/***********************************************************************************/

template< >
template< >
boost::numeric::ublas::bounded_matrix<double, 36, 36> MeshTyingMortarCondition<3,8,Vector3DValue>::CalculateLocalLHS<36>(
    const MortarConditionMatrices& rMortarConditionMatrices,
    DofData& rDofData,
    const unsigned int rMasterElementIndex,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    boost::numeric::ublas::bounded_matrix<double, 36, 36> lhs;

    // We get the mortar operators
    const boost::numeric::ublas::bounded_matrix<double, 4, 4>& MOperator = rMortarConditionMatrices.MOperator;
    const boost::numeric::ublas::bounded_matrix<double, 4, 4>& DOperator = rMortarConditionMatrices.DOperator;

    const double clhs0 =     -MOperator(0,0);
    const double clhs1 =     -MOperator(1,0);
    const double clhs2 =     -MOperator(2,0);
    const double clhs3 =     -MOperator(3,0);
    const double clhs4 =     -MOperator(0,1);
    const double clhs5 =     -MOperator(1,1);
    const double clhs6 =     -MOperator(2,1);
    const double clhs7 =     -MOperator(3,1);
    const double clhs8 =     -MOperator(0,2);
    const double clhs9 =     -MOperator(1,2);
    const double clhs10 =     -MOperator(2,2);
    const double clhs11 =     -MOperator(3,2);
    const double clhs12 =     -MOperator(0,3);
    const double clhs13 =     -MOperator(1,3);
    const double clhs14 =     -MOperator(2,3);
    const double clhs15 =     -MOperator(3,3);

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
    lhs(0,10)=0;
    lhs(0,11)=0;
    lhs(0,12)=0;
    lhs(0,13)=0;
    lhs(0,14)=0;
    lhs(0,15)=0;
    lhs(0,16)=0;
    lhs(0,17)=0;
    lhs(0,18)=0;
    lhs(0,19)=0;
    lhs(0,20)=0;
    lhs(0,21)=0;
    lhs(0,22)=0;
    lhs(0,23)=0;
    lhs(0,24)=clhs0;
    lhs(0,25)=0;
    lhs(0,26)=0;
    lhs(0,27)=clhs1;
    lhs(0,28)=0;
    lhs(0,29)=0;
    lhs(0,30)=clhs2;
    lhs(0,31)=0;
    lhs(0,32)=0;
    lhs(0,33)=clhs3;
    lhs(0,34)=0;
    lhs(0,35)=0;
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
    lhs(1,10)=0;
    lhs(1,11)=0;
    lhs(1,12)=0;
    lhs(1,13)=0;
    lhs(1,14)=0;
    lhs(1,15)=0;
    lhs(1,16)=0;
    lhs(1,17)=0;
    lhs(1,18)=0;
    lhs(1,19)=0;
    lhs(1,20)=0;
    lhs(1,21)=0;
    lhs(1,22)=0;
    lhs(1,23)=0;
    lhs(1,24)=0;
    lhs(1,25)=clhs0;
    lhs(1,26)=0;
    lhs(1,27)=0;
    lhs(1,28)=clhs1;
    lhs(1,29)=0;
    lhs(1,30)=0;
    lhs(1,31)=clhs2;
    lhs(1,32)=0;
    lhs(1,33)=0;
    lhs(1,34)=clhs3;
    lhs(1,35)=0;
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
    lhs(2,10)=0;
    lhs(2,11)=0;
    lhs(2,12)=0;
    lhs(2,13)=0;
    lhs(2,14)=0;
    lhs(2,15)=0;
    lhs(2,16)=0;
    lhs(2,17)=0;
    lhs(2,18)=0;
    lhs(2,19)=0;
    lhs(2,20)=0;
    lhs(2,21)=0;
    lhs(2,22)=0;
    lhs(2,23)=0;
    lhs(2,24)=0;
    lhs(2,25)=0;
    lhs(2,26)=clhs0;
    lhs(2,27)=0;
    lhs(2,28)=0;
    lhs(2,29)=clhs1;
    lhs(2,30)=0;
    lhs(2,31)=0;
    lhs(2,32)=clhs2;
    lhs(2,33)=0;
    lhs(2,34)=0;
    lhs(2,35)=clhs3;
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
    lhs(3,10)=0;
    lhs(3,11)=0;
    lhs(3,12)=0;
    lhs(3,13)=0;
    lhs(3,14)=0;
    lhs(3,15)=0;
    lhs(3,16)=0;
    lhs(3,17)=0;
    lhs(3,18)=0;
    lhs(3,19)=0;
    lhs(3,20)=0;
    lhs(3,21)=0;
    lhs(3,22)=0;
    lhs(3,23)=0;
    lhs(3,24)=clhs4;
    lhs(3,25)=0;
    lhs(3,26)=0;
    lhs(3,27)=clhs5;
    lhs(3,28)=0;
    lhs(3,29)=0;
    lhs(3,30)=clhs6;
    lhs(3,31)=0;
    lhs(3,32)=0;
    lhs(3,33)=clhs7;
    lhs(3,34)=0;
    lhs(3,35)=0;
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
    lhs(4,10)=0;
    lhs(4,11)=0;
    lhs(4,12)=0;
    lhs(4,13)=0;
    lhs(4,14)=0;
    lhs(4,15)=0;
    lhs(4,16)=0;
    lhs(4,17)=0;
    lhs(4,18)=0;
    lhs(4,19)=0;
    lhs(4,20)=0;
    lhs(4,21)=0;
    lhs(4,22)=0;
    lhs(4,23)=0;
    lhs(4,24)=0;
    lhs(4,25)=clhs4;
    lhs(4,26)=0;
    lhs(4,27)=0;
    lhs(4,28)=clhs5;
    lhs(4,29)=0;
    lhs(4,30)=0;
    lhs(4,31)=clhs6;
    lhs(4,32)=0;
    lhs(4,33)=0;
    lhs(4,34)=clhs7;
    lhs(4,35)=0;
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
    lhs(5,10)=0;
    lhs(5,11)=0;
    lhs(5,12)=0;
    lhs(5,13)=0;
    lhs(5,14)=0;
    lhs(5,15)=0;
    lhs(5,16)=0;
    lhs(5,17)=0;
    lhs(5,18)=0;
    lhs(5,19)=0;
    lhs(5,20)=0;
    lhs(5,21)=0;
    lhs(5,22)=0;
    lhs(5,23)=0;
    lhs(5,24)=0;
    lhs(5,25)=0;
    lhs(5,26)=clhs4;
    lhs(5,27)=0;
    lhs(5,28)=0;
    lhs(5,29)=clhs5;
    lhs(5,30)=0;
    lhs(5,31)=0;
    lhs(5,32)=clhs6;
    lhs(5,33)=0;
    lhs(5,34)=0;
    lhs(5,35)=clhs7;
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
    lhs(6,10)=0;
    lhs(6,11)=0;
    lhs(6,12)=0;
    lhs(6,13)=0;
    lhs(6,14)=0;
    lhs(6,15)=0;
    lhs(6,16)=0;
    lhs(6,17)=0;
    lhs(6,18)=0;
    lhs(6,19)=0;
    lhs(6,20)=0;
    lhs(6,21)=0;
    lhs(6,22)=0;
    lhs(6,23)=0;
    lhs(6,24)=clhs8;
    lhs(6,25)=0;
    lhs(6,26)=0;
    lhs(6,27)=clhs9;
    lhs(6,28)=0;
    lhs(6,29)=0;
    lhs(6,30)=clhs10;
    lhs(6,31)=0;
    lhs(6,32)=0;
    lhs(6,33)=clhs11;
    lhs(6,34)=0;
    lhs(6,35)=0;
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
    lhs(7,10)=0;
    lhs(7,11)=0;
    lhs(7,12)=0;
    lhs(7,13)=0;
    lhs(7,14)=0;
    lhs(7,15)=0;
    lhs(7,16)=0;
    lhs(7,17)=0;
    lhs(7,18)=0;
    lhs(7,19)=0;
    lhs(7,20)=0;
    lhs(7,21)=0;
    lhs(7,22)=0;
    lhs(7,23)=0;
    lhs(7,24)=0;
    lhs(7,25)=clhs8;
    lhs(7,26)=0;
    lhs(7,27)=0;
    lhs(7,28)=clhs9;
    lhs(7,29)=0;
    lhs(7,30)=0;
    lhs(7,31)=clhs10;
    lhs(7,32)=0;
    lhs(7,33)=0;
    lhs(7,34)=clhs11;
    lhs(7,35)=0;
    lhs(8,0)=0;
    lhs(8,1)=0;
    lhs(8,2)=0;
    lhs(8,3)=0;
    lhs(8,4)=0;
    lhs(8,5)=0;
    lhs(8,6)=0;
    lhs(8,7)=0;
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(8,12)=0;
    lhs(8,13)=0;
    lhs(8,14)=0;
    lhs(8,15)=0;
    lhs(8,16)=0;
    lhs(8,17)=0;
    lhs(8,18)=0;
    lhs(8,19)=0;
    lhs(8,20)=0;
    lhs(8,21)=0;
    lhs(8,22)=0;
    lhs(8,23)=0;
    lhs(8,24)=0;
    lhs(8,25)=0;
    lhs(8,26)=clhs8;
    lhs(8,27)=0;
    lhs(8,28)=0;
    lhs(8,29)=clhs9;
    lhs(8,30)=0;
    lhs(8,31)=0;
    lhs(8,32)=clhs10;
    lhs(8,33)=0;
    lhs(8,34)=0;
    lhs(8,35)=clhs11;
    lhs(9,0)=0;
    lhs(9,1)=0;
    lhs(9,2)=0;
    lhs(9,3)=0;
    lhs(9,4)=0;
    lhs(9,5)=0;
    lhs(9,6)=0;
    lhs(9,7)=0;
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(9,12)=0;
    lhs(9,13)=0;
    lhs(9,14)=0;
    lhs(9,15)=0;
    lhs(9,16)=0;
    lhs(9,17)=0;
    lhs(9,18)=0;
    lhs(9,19)=0;
    lhs(9,20)=0;
    lhs(9,21)=0;
    lhs(9,22)=0;
    lhs(9,23)=0;
    lhs(9,24)=clhs12;
    lhs(9,25)=0;
    lhs(9,26)=0;
    lhs(9,27)=clhs13;
    lhs(9,28)=0;
    lhs(9,29)=0;
    lhs(9,30)=clhs14;
    lhs(9,31)=0;
    lhs(9,32)=0;
    lhs(9,33)=clhs15;
    lhs(9,34)=0;
    lhs(9,35)=0;
    lhs(10,0)=0;
    lhs(10,1)=0;
    lhs(10,2)=0;
    lhs(10,3)=0;
    lhs(10,4)=0;
    lhs(10,5)=0;
    lhs(10,6)=0;
    lhs(10,7)=0;
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(10,12)=0;
    lhs(10,13)=0;
    lhs(10,14)=0;
    lhs(10,15)=0;
    lhs(10,16)=0;
    lhs(10,17)=0;
    lhs(10,18)=0;
    lhs(10,19)=0;
    lhs(10,20)=0;
    lhs(10,21)=0;
    lhs(10,22)=0;
    lhs(10,23)=0;
    lhs(10,24)=0;
    lhs(10,25)=clhs12;
    lhs(10,26)=0;
    lhs(10,27)=0;
    lhs(10,28)=clhs13;
    lhs(10,29)=0;
    lhs(10,30)=0;
    lhs(10,31)=clhs14;
    lhs(10,32)=0;
    lhs(10,33)=0;
    lhs(10,34)=clhs15;
    lhs(10,35)=0;
    lhs(11,0)=0;
    lhs(11,1)=0;
    lhs(11,2)=0;
    lhs(11,3)=0;
    lhs(11,4)=0;
    lhs(11,5)=0;
    lhs(11,6)=0;
    lhs(11,7)=0;
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;
    lhs(11,12)=0;
    lhs(11,13)=0;
    lhs(11,14)=0;
    lhs(11,15)=0;
    lhs(11,16)=0;
    lhs(11,17)=0;
    lhs(11,18)=0;
    lhs(11,19)=0;
    lhs(11,20)=0;
    lhs(11,21)=0;
    lhs(11,22)=0;
    lhs(11,23)=0;
    lhs(11,24)=0;
    lhs(11,25)=0;
    lhs(11,26)=clhs12;
    lhs(11,27)=0;
    lhs(11,28)=0;
    lhs(11,29)=clhs13;
    lhs(11,30)=0;
    lhs(11,31)=0;
    lhs(11,32)=clhs14;
    lhs(11,33)=0;
    lhs(11,34)=0;
    lhs(11,35)=clhs15;
    lhs(12,0)=0;
    lhs(12,1)=0;
    lhs(12,2)=0;
    lhs(12,3)=0;
    lhs(12,4)=0;
    lhs(12,5)=0;
    lhs(12,6)=0;
    lhs(12,7)=0;
    lhs(12,8)=0;
    lhs(12,9)=0;
    lhs(12,10)=0;
    lhs(12,11)=0;
    lhs(12,12)=0;
    lhs(12,13)=0;
    lhs(12,14)=0;
    lhs(12,15)=0;
    lhs(12,16)=0;
    lhs(12,17)=0;
    lhs(12,18)=0;
    lhs(12,19)=0;
    lhs(12,20)=0;
    lhs(12,21)=0;
    lhs(12,22)=0;
    lhs(12,23)=0;
    lhs(12,24)=DOperator(0,0);
    lhs(12,25)=0;
    lhs(12,26)=0;
    lhs(12,27)=DOperator(1,0);
    lhs(12,28)=0;
    lhs(12,29)=0;
    lhs(12,30)=DOperator(2,0);
    lhs(12,31)=0;
    lhs(12,32)=0;
    lhs(12,33)=DOperator(3,0);
    lhs(12,34)=0;
    lhs(12,35)=0;
    lhs(13,0)=0;
    lhs(13,1)=0;
    lhs(13,2)=0;
    lhs(13,3)=0;
    lhs(13,4)=0;
    lhs(13,5)=0;
    lhs(13,6)=0;
    lhs(13,7)=0;
    lhs(13,8)=0;
    lhs(13,9)=0;
    lhs(13,10)=0;
    lhs(13,11)=0;
    lhs(13,12)=0;
    lhs(13,13)=0;
    lhs(13,14)=0;
    lhs(13,15)=0;
    lhs(13,16)=0;
    lhs(13,17)=0;
    lhs(13,18)=0;
    lhs(13,19)=0;
    lhs(13,20)=0;
    lhs(13,21)=0;
    lhs(13,22)=0;
    lhs(13,23)=0;
    lhs(13,24)=0;
    lhs(13,25)=DOperator(0,0);
    lhs(13,26)=0;
    lhs(13,27)=0;
    lhs(13,28)=DOperator(1,0);
    lhs(13,29)=0;
    lhs(13,30)=0;
    lhs(13,31)=DOperator(2,0);
    lhs(13,32)=0;
    lhs(13,33)=0;
    lhs(13,34)=DOperator(3,0);
    lhs(13,35)=0;
    lhs(14,0)=0;
    lhs(14,1)=0;
    lhs(14,2)=0;
    lhs(14,3)=0;
    lhs(14,4)=0;
    lhs(14,5)=0;
    lhs(14,6)=0;
    lhs(14,7)=0;
    lhs(14,8)=0;
    lhs(14,9)=0;
    lhs(14,10)=0;
    lhs(14,11)=0;
    lhs(14,12)=0;
    lhs(14,13)=0;
    lhs(14,14)=0;
    lhs(14,15)=0;
    lhs(14,16)=0;
    lhs(14,17)=0;
    lhs(14,18)=0;
    lhs(14,19)=0;
    lhs(14,20)=0;
    lhs(14,21)=0;
    lhs(14,22)=0;
    lhs(14,23)=0;
    lhs(14,24)=0;
    lhs(14,25)=0;
    lhs(14,26)=DOperator(0,0);
    lhs(14,27)=0;
    lhs(14,28)=0;
    lhs(14,29)=DOperator(1,0);
    lhs(14,30)=0;
    lhs(14,31)=0;
    lhs(14,32)=DOperator(2,0);
    lhs(14,33)=0;
    lhs(14,34)=0;
    lhs(14,35)=DOperator(3,0);
    lhs(15,0)=0;
    lhs(15,1)=0;
    lhs(15,2)=0;
    lhs(15,3)=0;
    lhs(15,4)=0;
    lhs(15,5)=0;
    lhs(15,6)=0;
    lhs(15,7)=0;
    lhs(15,8)=0;
    lhs(15,9)=0;
    lhs(15,10)=0;
    lhs(15,11)=0;
    lhs(15,12)=0;
    lhs(15,13)=0;
    lhs(15,14)=0;
    lhs(15,15)=0;
    lhs(15,16)=0;
    lhs(15,17)=0;
    lhs(15,18)=0;
    lhs(15,19)=0;
    lhs(15,20)=0;
    lhs(15,21)=0;
    lhs(15,22)=0;
    lhs(15,23)=0;
    lhs(15,24)=DOperator(0,1);
    lhs(15,25)=0;
    lhs(15,26)=0;
    lhs(15,27)=DOperator(1,1);
    lhs(15,28)=0;
    lhs(15,29)=0;
    lhs(15,30)=DOperator(2,1);
    lhs(15,31)=0;
    lhs(15,32)=0;
    lhs(15,33)=DOperator(3,1);
    lhs(15,34)=0;
    lhs(15,35)=0;
    lhs(16,0)=0;
    lhs(16,1)=0;
    lhs(16,2)=0;
    lhs(16,3)=0;
    lhs(16,4)=0;
    lhs(16,5)=0;
    lhs(16,6)=0;
    lhs(16,7)=0;
    lhs(16,8)=0;
    lhs(16,9)=0;
    lhs(16,10)=0;
    lhs(16,11)=0;
    lhs(16,12)=0;
    lhs(16,13)=0;
    lhs(16,14)=0;
    lhs(16,15)=0;
    lhs(16,16)=0;
    lhs(16,17)=0;
    lhs(16,18)=0;
    lhs(16,19)=0;
    lhs(16,20)=0;
    lhs(16,21)=0;
    lhs(16,22)=0;
    lhs(16,23)=0;
    lhs(16,24)=0;
    lhs(16,25)=DOperator(0,1);
    lhs(16,26)=0;
    lhs(16,27)=0;
    lhs(16,28)=DOperator(1,1);
    lhs(16,29)=0;
    lhs(16,30)=0;
    lhs(16,31)=DOperator(2,1);
    lhs(16,32)=0;
    lhs(16,33)=0;
    lhs(16,34)=DOperator(3,1);
    lhs(16,35)=0;
    lhs(17,0)=0;
    lhs(17,1)=0;
    lhs(17,2)=0;
    lhs(17,3)=0;
    lhs(17,4)=0;
    lhs(17,5)=0;
    lhs(17,6)=0;
    lhs(17,7)=0;
    lhs(17,8)=0;
    lhs(17,9)=0;
    lhs(17,10)=0;
    lhs(17,11)=0;
    lhs(17,12)=0;
    lhs(17,13)=0;
    lhs(17,14)=0;
    lhs(17,15)=0;
    lhs(17,16)=0;
    lhs(17,17)=0;
    lhs(17,18)=0;
    lhs(17,19)=0;
    lhs(17,20)=0;
    lhs(17,21)=0;
    lhs(17,22)=0;
    lhs(17,23)=0;
    lhs(17,24)=0;
    lhs(17,25)=0;
    lhs(17,26)=DOperator(0,1);
    lhs(17,27)=0;
    lhs(17,28)=0;
    lhs(17,29)=DOperator(1,1);
    lhs(17,30)=0;
    lhs(17,31)=0;
    lhs(17,32)=DOperator(2,1);
    lhs(17,33)=0;
    lhs(17,34)=0;
    lhs(17,35)=DOperator(3,1);
    lhs(18,0)=0;
    lhs(18,1)=0;
    lhs(18,2)=0;
    lhs(18,3)=0;
    lhs(18,4)=0;
    lhs(18,5)=0;
    lhs(18,6)=0;
    lhs(18,7)=0;
    lhs(18,8)=0;
    lhs(18,9)=0;
    lhs(18,10)=0;
    lhs(18,11)=0;
    lhs(18,12)=0;
    lhs(18,13)=0;
    lhs(18,14)=0;
    lhs(18,15)=0;
    lhs(18,16)=0;
    lhs(18,17)=0;
    lhs(18,18)=0;
    lhs(18,19)=0;
    lhs(18,20)=0;
    lhs(18,21)=0;
    lhs(18,22)=0;
    lhs(18,23)=0;
    lhs(18,24)=DOperator(0,2);
    lhs(18,25)=0;
    lhs(18,26)=0;
    lhs(18,27)=DOperator(1,2);
    lhs(18,28)=0;
    lhs(18,29)=0;
    lhs(18,30)=DOperator(2,2);
    lhs(18,31)=0;
    lhs(18,32)=0;
    lhs(18,33)=DOperator(3,2);
    lhs(18,34)=0;
    lhs(18,35)=0;
    lhs(19,0)=0;
    lhs(19,1)=0;
    lhs(19,2)=0;
    lhs(19,3)=0;
    lhs(19,4)=0;
    lhs(19,5)=0;
    lhs(19,6)=0;
    lhs(19,7)=0;
    lhs(19,8)=0;
    lhs(19,9)=0;
    lhs(19,10)=0;
    lhs(19,11)=0;
    lhs(19,12)=0;
    lhs(19,13)=0;
    lhs(19,14)=0;
    lhs(19,15)=0;
    lhs(19,16)=0;
    lhs(19,17)=0;
    lhs(19,18)=0;
    lhs(19,19)=0;
    lhs(19,20)=0;
    lhs(19,21)=0;
    lhs(19,22)=0;
    lhs(19,23)=0;
    lhs(19,24)=0;
    lhs(19,25)=DOperator(0,2);
    lhs(19,26)=0;
    lhs(19,27)=0;
    lhs(19,28)=DOperator(1,2);
    lhs(19,29)=0;
    lhs(19,30)=0;
    lhs(19,31)=DOperator(2,2);
    lhs(19,32)=0;
    lhs(19,33)=0;
    lhs(19,34)=DOperator(3,2);
    lhs(19,35)=0;
    lhs(20,0)=0;
    lhs(20,1)=0;
    lhs(20,2)=0;
    lhs(20,3)=0;
    lhs(20,4)=0;
    lhs(20,5)=0;
    lhs(20,6)=0;
    lhs(20,7)=0;
    lhs(20,8)=0;
    lhs(20,9)=0;
    lhs(20,10)=0;
    lhs(20,11)=0;
    lhs(20,12)=0;
    lhs(20,13)=0;
    lhs(20,14)=0;
    lhs(20,15)=0;
    lhs(20,16)=0;
    lhs(20,17)=0;
    lhs(20,18)=0;
    lhs(20,19)=0;
    lhs(20,20)=0;
    lhs(20,21)=0;
    lhs(20,22)=0;
    lhs(20,23)=0;
    lhs(20,24)=0;
    lhs(20,25)=0;
    lhs(20,26)=DOperator(0,2);
    lhs(20,27)=0;
    lhs(20,28)=0;
    lhs(20,29)=DOperator(1,2);
    lhs(20,30)=0;
    lhs(20,31)=0;
    lhs(20,32)=DOperator(2,2);
    lhs(20,33)=0;
    lhs(20,34)=0;
    lhs(20,35)=DOperator(3,2);
    lhs(21,0)=0;
    lhs(21,1)=0;
    lhs(21,2)=0;
    lhs(21,3)=0;
    lhs(21,4)=0;
    lhs(21,5)=0;
    lhs(21,6)=0;
    lhs(21,7)=0;
    lhs(21,8)=0;
    lhs(21,9)=0;
    lhs(21,10)=0;
    lhs(21,11)=0;
    lhs(21,12)=0;
    lhs(21,13)=0;
    lhs(21,14)=0;
    lhs(21,15)=0;
    lhs(21,16)=0;
    lhs(21,17)=0;
    lhs(21,18)=0;
    lhs(21,19)=0;
    lhs(21,20)=0;
    lhs(21,21)=0;
    lhs(21,22)=0;
    lhs(21,23)=0;
    lhs(21,24)=DOperator(0,3);
    lhs(21,25)=0;
    lhs(21,26)=0;
    lhs(21,27)=DOperator(1,3);
    lhs(21,28)=0;
    lhs(21,29)=0;
    lhs(21,30)=DOperator(2,3);
    lhs(21,31)=0;
    lhs(21,32)=0;
    lhs(21,33)=DOperator(3,3);
    lhs(21,34)=0;
    lhs(21,35)=0;
    lhs(22,0)=0;
    lhs(22,1)=0;
    lhs(22,2)=0;
    lhs(22,3)=0;
    lhs(22,4)=0;
    lhs(22,5)=0;
    lhs(22,6)=0;
    lhs(22,7)=0;
    lhs(22,8)=0;
    lhs(22,9)=0;
    lhs(22,10)=0;
    lhs(22,11)=0;
    lhs(22,12)=0;
    lhs(22,13)=0;
    lhs(22,14)=0;
    lhs(22,15)=0;
    lhs(22,16)=0;
    lhs(22,17)=0;
    lhs(22,18)=0;
    lhs(22,19)=0;
    lhs(22,20)=0;
    lhs(22,21)=0;
    lhs(22,22)=0;
    lhs(22,23)=0;
    lhs(22,24)=0;
    lhs(22,25)=DOperator(0,3);
    lhs(22,26)=0;
    lhs(22,27)=0;
    lhs(22,28)=DOperator(1,3);
    lhs(22,29)=0;
    lhs(22,30)=0;
    lhs(22,31)=DOperator(2,3);
    lhs(22,32)=0;
    lhs(22,33)=0;
    lhs(22,34)=DOperator(3,3);
    lhs(22,35)=0;
    lhs(23,0)=0;
    lhs(23,1)=0;
    lhs(23,2)=0;
    lhs(23,3)=0;
    lhs(23,4)=0;
    lhs(23,5)=0;
    lhs(23,6)=0;
    lhs(23,7)=0;
    lhs(23,8)=0;
    lhs(23,9)=0;
    lhs(23,10)=0;
    lhs(23,11)=0;
    lhs(23,12)=0;
    lhs(23,13)=0;
    lhs(23,14)=0;
    lhs(23,15)=0;
    lhs(23,16)=0;
    lhs(23,17)=0;
    lhs(23,18)=0;
    lhs(23,19)=0;
    lhs(23,20)=0;
    lhs(23,21)=0;
    lhs(23,22)=0;
    lhs(23,23)=0;
    lhs(23,24)=0;
    lhs(23,25)=0;
    lhs(23,26)=DOperator(0,3);
    lhs(23,27)=0;
    lhs(23,28)=0;
    lhs(23,29)=DOperator(1,3);
    lhs(23,30)=0;
    lhs(23,31)=0;
    lhs(23,32)=DOperator(2,3);
    lhs(23,33)=0;
    lhs(23,34)=0;
    lhs(23,35)=DOperator(3,3);
    lhs(24,0)=clhs0;
    lhs(24,1)=0;
    lhs(24,2)=0;
    lhs(24,3)=clhs4;
    lhs(24,4)=0;
    lhs(24,5)=0;
    lhs(24,6)=clhs8;
    lhs(24,7)=0;
    lhs(24,8)=0;
    lhs(24,9)=clhs12;
    lhs(24,10)=0;
    lhs(24,11)=0;
    lhs(24,12)=DOperator(0,0);
    lhs(24,13)=0;
    lhs(24,14)=0;
    lhs(24,15)=DOperator(0,1);
    lhs(24,16)=0;
    lhs(24,17)=0;
    lhs(24,18)=DOperator(0,2);
    lhs(24,19)=0;
    lhs(24,20)=0;
    lhs(24,21)=DOperator(0,3);
    lhs(24,22)=0;
    lhs(24,23)=0;
    lhs(24,24)=0;
    lhs(24,25)=0;
    lhs(24,26)=0;
    lhs(24,27)=0;
    lhs(24,28)=0;
    lhs(24,29)=0;
    lhs(24,30)=0;
    lhs(24,31)=0;
    lhs(24,32)=0;
    lhs(24,33)=0;
    lhs(24,34)=0;
    lhs(24,35)=0;
    lhs(25,0)=0;
    lhs(25,1)=clhs0;
    lhs(25,2)=0;
    lhs(25,3)=0;
    lhs(25,4)=clhs4;
    lhs(25,5)=0;
    lhs(25,6)=0;
    lhs(25,7)=clhs8;
    lhs(25,8)=0;
    lhs(25,9)=0;
    lhs(25,10)=clhs12;
    lhs(25,11)=0;
    lhs(25,12)=0;
    lhs(25,13)=DOperator(0,0);
    lhs(25,14)=0;
    lhs(25,15)=0;
    lhs(25,16)=DOperator(0,1);
    lhs(25,17)=0;
    lhs(25,18)=0;
    lhs(25,19)=DOperator(0,2);
    lhs(25,20)=0;
    lhs(25,21)=0;
    lhs(25,22)=DOperator(0,3);
    lhs(25,23)=0;
    lhs(25,24)=0;
    lhs(25,25)=0;
    lhs(25,26)=0;
    lhs(25,27)=0;
    lhs(25,28)=0;
    lhs(25,29)=0;
    lhs(25,30)=0;
    lhs(25,31)=0;
    lhs(25,32)=0;
    lhs(25,33)=0;
    lhs(25,34)=0;
    lhs(25,35)=0;
    lhs(26,0)=0;
    lhs(26,1)=0;
    lhs(26,2)=clhs0;
    lhs(26,3)=0;
    lhs(26,4)=0;
    lhs(26,5)=clhs4;
    lhs(26,6)=0;
    lhs(26,7)=0;
    lhs(26,8)=clhs8;
    lhs(26,9)=0;
    lhs(26,10)=0;
    lhs(26,11)=clhs12;
    lhs(26,12)=0;
    lhs(26,13)=0;
    lhs(26,14)=DOperator(0,0);
    lhs(26,15)=0;
    lhs(26,16)=0;
    lhs(26,17)=DOperator(0,1);
    lhs(26,18)=0;
    lhs(26,19)=0;
    lhs(26,20)=DOperator(0,2);
    lhs(26,21)=0;
    lhs(26,22)=0;
    lhs(26,23)=DOperator(0,3);
    lhs(26,24)=0;
    lhs(26,25)=0;
    lhs(26,26)=0;
    lhs(26,27)=0;
    lhs(26,28)=0;
    lhs(26,29)=0;
    lhs(26,30)=0;
    lhs(26,31)=0;
    lhs(26,32)=0;
    lhs(26,33)=0;
    lhs(26,34)=0;
    lhs(26,35)=0;
    lhs(27,0)=clhs1;
    lhs(27,1)=0;
    lhs(27,2)=0;
    lhs(27,3)=clhs5;
    lhs(27,4)=0;
    lhs(27,5)=0;
    lhs(27,6)=clhs9;
    lhs(27,7)=0;
    lhs(27,8)=0;
    lhs(27,9)=clhs13;
    lhs(27,10)=0;
    lhs(27,11)=0;
    lhs(27,12)=DOperator(1,0);
    lhs(27,13)=0;
    lhs(27,14)=0;
    lhs(27,15)=DOperator(1,1);
    lhs(27,16)=0;
    lhs(27,17)=0;
    lhs(27,18)=DOperator(1,2);
    lhs(27,19)=0;
    lhs(27,20)=0;
    lhs(27,21)=DOperator(1,3);
    lhs(27,22)=0;
    lhs(27,23)=0;
    lhs(27,24)=0;
    lhs(27,25)=0;
    lhs(27,26)=0;
    lhs(27,27)=0;
    lhs(27,28)=0;
    lhs(27,29)=0;
    lhs(27,30)=0;
    lhs(27,31)=0;
    lhs(27,32)=0;
    lhs(27,33)=0;
    lhs(27,34)=0;
    lhs(27,35)=0;
    lhs(28,0)=0;
    lhs(28,1)=clhs1;
    lhs(28,2)=0;
    lhs(28,3)=0;
    lhs(28,4)=clhs5;
    lhs(28,5)=0;
    lhs(28,6)=0;
    lhs(28,7)=clhs9;
    lhs(28,8)=0;
    lhs(28,9)=0;
    lhs(28,10)=clhs13;
    lhs(28,11)=0;
    lhs(28,12)=0;
    lhs(28,13)=DOperator(1,0);
    lhs(28,14)=0;
    lhs(28,15)=0;
    lhs(28,16)=DOperator(1,1);
    lhs(28,17)=0;
    lhs(28,18)=0;
    lhs(28,19)=DOperator(1,2);
    lhs(28,20)=0;
    lhs(28,21)=0;
    lhs(28,22)=DOperator(1,3);
    lhs(28,23)=0;
    lhs(28,24)=0;
    lhs(28,25)=0;
    lhs(28,26)=0;
    lhs(28,27)=0;
    lhs(28,28)=0;
    lhs(28,29)=0;
    lhs(28,30)=0;
    lhs(28,31)=0;
    lhs(28,32)=0;
    lhs(28,33)=0;
    lhs(28,34)=0;
    lhs(28,35)=0;
    lhs(29,0)=0;
    lhs(29,1)=0;
    lhs(29,2)=clhs1;
    lhs(29,3)=0;
    lhs(29,4)=0;
    lhs(29,5)=clhs5;
    lhs(29,6)=0;
    lhs(29,7)=0;
    lhs(29,8)=clhs9;
    lhs(29,9)=0;
    lhs(29,10)=0;
    lhs(29,11)=clhs13;
    lhs(29,12)=0;
    lhs(29,13)=0;
    lhs(29,14)=DOperator(1,0);
    lhs(29,15)=0;
    lhs(29,16)=0;
    lhs(29,17)=DOperator(1,1);
    lhs(29,18)=0;
    lhs(29,19)=0;
    lhs(29,20)=DOperator(1,2);
    lhs(29,21)=0;
    lhs(29,22)=0;
    lhs(29,23)=DOperator(1,3);
    lhs(29,24)=0;
    lhs(29,25)=0;
    lhs(29,26)=0;
    lhs(29,27)=0;
    lhs(29,28)=0;
    lhs(29,29)=0;
    lhs(29,30)=0;
    lhs(29,31)=0;
    lhs(29,32)=0;
    lhs(29,33)=0;
    lhs(29,34)=0;
    lhs(29,35)=0;
    lhs(30,0)=clhs2;
    lhs(30,1)=0;
    lhs(30,2)=0;
    lhs(30,3)=clhs6;
    lhs(30,4)=0;
    lhs(30,5)=0;
    lhs(30,6)=clhs10;
    lhs(30,7)=0;
    lhs(30,8)=0;
    lhs(30,9)=clhs14;
    lhs(30,10)=0;
    lhs(30,11)=0;
    lhs(30,12)=DOperator(2,0);
    lhs(30,13)=0;
    lhs(30,14)=0;
    lhs(30,15)=DOperator(2,1);
    lhs(30,16)=0;
    lhs(30,17)=0;
    lhs(30,18)=DOperator(2,2);
    lhs(30,19)=0;
    lhs(30,20)=0;
    lhs(30,21)=DOperator(2,3);
    lhs(30,22)=0;
    lhs(30,23)=0;
    lhs(30,24)=0;
    lhs(30,25)=0;
    lhs(30,26)=0;
    lhs(30,27)=0;
    lhs(30,28)=0;
    lhs(30,29)=0;
    lhs(30,30)=0;
    lhs(30,31)=0;
    lhs(30,32)=0;
    lhs(30,33)=0;
    lhs(30,34)=0;
    lhs(30,35)=0;
    lhs(31,0)=0;
    lhs(31,1)=clhs2;
    lhs(31,2)=0;
    lhs(31,3)=0;
    lhs(31,4)=clhs6;
    lhs(31,5)=0;
    lhs(31,6)=0;
    lhs(31,7)=clhs10;
    lhs(31,8)=0;
    lhs(31,9)=0;
    lhs(31,10)=clhs14;
    lhs(31,11)=0;
    lhs(31,12)=0;
    lhs(31,13)=DOperator(2,0);
    lhs(31,14)=0;
    lhs(31,15)=0;
    lhs(31,16)=DOperator(2,1);
    lhs(31,17)=0;
    lhs(31,18)=0;
    lhs(31,19)=DOperator(2,2);
    lhs(31,20)=0;
    lhs(31,21)=0;
    lhs(31,22)=DOperator(2,3);
    lhs(31,23)=0;
    lhs(31,24)=0;
    lhs(31,25)=0;
    lhs(31,26)=0;
    lhs(31,27)=0;
    lhs(31,28)=0;
    lhs(31,29)=0;
    lhs(31,30)=0;
    lhs(31,31)=0;
    lhs(31,32)=0;
    lhs(31,33)=0;
    lhs(31,34)=0;
    lhs(31,35)=0;
    lhs(32,0)=0;
    lhs(32,1)=0;
    lhs(32,2)=clhs2;
    lhs(32,3)=0;
    lhs(32,4)=0;
    lhs(32,5)=clhs6;
    lhs(32,6)=0;
    lhs(32,7)=0;
    lhs(32,8)=clhs10;
    lhs(32,9)=0;
    lhs(32,10)=0;
    lhs(32,11)=clhs14;
    lhs(32,12)=0;
    lhs(32,13)=0;
    lhs(32,14)=DOperator(2,0);
    lhs(32,15)=0;
    lhs(32,16)=0;
    lhs(32,17)=DOperator(2,1);
    lhs(32,18)=0;
    lhs(32,19)=0;
    lhs(32,20)=DOperator(2,2);
    lhs(32,21)=0;
    lhs(32,22)=0;
    lhs(32,23)=DOperator(2,3);
    lhs(32,24)=0;
    lhs(32,25)=0;
    lhs(32,26)=0;
    lhs(32,27)=0;
    lhs(32,28)=0;
    lhs(32,29)=0;
    lhs(32,30)=0;
    lhs(32,31)=0;
    lhs(32,32)=0;
    lhs(32,33)=0;
    lhs(32,34)=0;
    lhs(32,35)=0;
    lhs(33,0)=clhs3;
    lhs(33,1)=0;
    lhs(33,2)=0;
    lhs(33,3)=clhs7;
    lhs(33,4)=0;
    lhs(33,5)=0;
    lhs(33,6)=clhs11;
    lhs(33,7)=0;
    lhs(33,8)=0;
    lhs(33,9)=clhs15;
    lhs(33,10)=0;
    lhs(33,11)=0;
    lhs(33,12)=DOperator(3,0);
    lhs(33,13)=0;
    lhs(33,14)=0;
    lhs(33,15)=DOperator(3,1);
    lhs(33,16)=0;
    lhs(33,17)=0;
    lhs(33,18)=DOperator(3,2);
    lhs(33,19)=0;
    lhs(33,20)=0;
    lhs(33,21)=DOperator(3,3);
    lhs(33,22)=0;
    lhs(33,23)=0;
    lhs(33,24)=0;
    lhs(33,25)=0;
    lhs(33,26)=0;
    lhs(33,27)=0;
    lhs(33,28)=0;
    lhs(33,29)=0;
    lhs(33,30)=0;
    lhs(33,31)=0;
    lhs(33,32)=0;
    lhs(33,33)=0;
    lhs(33,34)=0;
    lhs(33,35)=0;
    lhs(34,0)=0;
    lhs(34,1)=clhs3;
    lhs(34,2)=0;
    lhs(34,3)=0;
    lhs(34,4)=clhs7;
    lhs(34,5)=0;
    lhs(34,6)=0;
    lhs(34,7)=clhs11;
    lhs(34,8)=0;
    lhs(34,9)=0;
    lhs(34,10)=clhs15;
    lhs(34,11)=0;
    lhs(34,12)=0;
    lhs(34,13)=DOperator(3,0);
    lhs(34,14)=0;
    lhs(34,15)=0;
    lhs(34,16)=DOperator(3,1);
    lhs(34,17)=0;
    lhs(34,18)=0;
    lhs(34,19)=DOperator(3,2);
    lhs(34,20)=0;
    lhs(34,21)=0;
    lhs(34,22)=DOperator(3,3);
    lhs(34,23)=0;
    lhs(34,24)=0;
    lhs(34,25)=0;
    lhs(34,26)=0;
    lhs(34,27)=0;
    lhs(34,28)=0;
    lhs(34,29)=0;
    lhs(34,30)=0;
    lhs(34,31)=0;
    lhs(34,32)=0;
    lhs(34,33)=0;
    lhs(34,34)=0;
    lhs(34,35)=0;
    lhs(35,0)=0;
    lhs(35,1)=0;
    lhs(35,2)=clhs3;
    lhs(35,3)=0;
    lhs(35,4)=0;
    lhs(35,5)=clhs7;
    lhs(35,6)=0;
    lhs(35,7)=0;
    lhs(35,8)=clhs11;
    lhs(35,9)=0;
    lhs(35,10)=0;
    lhs(35,11)=clhs15;
    lhs(35,12)=0;
    lhs(35,13)=0;
    lhs(35,14)=DOperator(3,0);
    lhs(35,15)=0;
    lhs(35,16)=0;
    lhs(35,17)=DOperator(3,1);
    lhs(35,18)=0;
    lhs(35,19)=0;
    lhs(35,20)=DOperator(3,2);
    lhs(35,21)=0;
    lhs(35,22)=0;
    lhs(35,23)=DOperator(3,3);
    lhs(35,24)=0;
    lhs(35,25)=0;
    lhs(35,26)=0;
    lhs(35,27)=0;
    lhs(35,28)=0;
    lhs(35,29)=0;
    lhs(35,30)=0;
    lhs(35,31)=0;
    lhs(35,32)=0;
    lhs(35,33)=0;
    lhs(35,34)=0;
    lhs(35,35)=0;


    return lhs;
}


/****************************** END AD REPLACEMENT *********************************/
/***********************************************************************************/

/***************************** BEGIN AD REPLACEMENT ********************************/
/***********************************************************************************/

template<>
template<>
array_1d<double, 6> MeshTyingMortarCondition<2,3,ScalarValue>::CalculateLocalRHS<6>(
    const MortarConditionMatrices& rMortarConditionMatrices,
    DofData& rDofData,
    const unsigned int rMasterElementIndex,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    array_1d<double,6> rhs;

    // Initialize values
    const bounded_matrix<double, 2, ScalarValue> u1 = rDofData.u1;
    const bounded_matrix<double, 2, ScalarValue> u2 = rDofData.u2;

    const bounded_matrix<double, 2, ScalarValue> lm = rDofData.LagrangeMultipliers; 

    // Mortar operators
    const bounded_matrix<double, 2, 2>& MOperator = rMortarConditionMatrices.MOperator;
    const bounded_matrix<double, 2, 2>& DOperator = rMortarConditionMatrices.DOperator;


    rhs[0]=MOperator(0,0)*lm(0,0) + MOperator(1,0)*lm(1,0);
    rhs[1]=MOperator(0,1)*lm(0,0) + MOperator(1,1)*lm(1,0);
    rhs[2]=-(DOperator(0,0)*lm(0,0) + DOperator(1,0)*lm(1,0));
    rhs[3]=-(DOperator(0,1)*lm(0,0) + DOperator(1,1)*lm(1,0));
    rhs[4]=-DOperator(0,0)*u1(0,0) - DOperator(0,1)*u1(1,0) + MOperator(0,0)*u2(0,0) + MOperator(0,1)*u2(1,0);
    rhs[5]=-DOperator(1,0)*u1(0,0) - DOperator(1,1)*u1(1,0) + MOperator(1,0)*u2(0,0) + MOperator(1,1)*u2(1,0);


    return rhs;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
template<>
array_1d<double, 12> MeshTyingMortarCondition<2,3,Vector2DValue>::CalculateLocalRHS<12>(
    const MortarConditionMatrices& rMortarConditionMatrices,
    DofData& rDofData,
    const unsigned int rMasterElementIndex,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    array_1d<double,12> rhs;

    // Initialize values
    const bounded_matrix<double, 2, Vector2DValue> u1 = rDofData.u1;
    const bounded_matrix<double, 2, Vector2DValue> u2 = rDofData.u2;

    const bounded_matrix<double, 2, Vector2DValue> lm = rDofData.LagrangeMultipliers; 

    // Mortar operators
    const bounded_matrix<double, 2, 2>& MOperator = rMortarConditionMatrices.MOperator;
    const bounded_matrix<double, 2, 2>& DOperator = rMortarConditionMatrices.DOperator;


    rhs[0]=MOperator(0,0)*lm(0,0) + MOperator(1,0)*lm(1,0);
    rhs[1]=MOperator(0,0)*lm(0,1) + MOperator(1,0)*lm(1,1);
    rhs[2]=MOperator(0,1)*lm(0,0) + MOperator(1,1)*lm(1,0);
    rhs[3]=MOperator(0,1)*lm(0,1) + MOperator(1,1)*lm(1,1);
    rhs[4]=-(DOperator(0,0)*lm(0,0) + DOperator(1,0)*lm(1,0));
    rhs[5]=-(DOperator(0,0)*lm(0,1) + DOperator(1,0)*lm(1,1));
    rhs[6]=-(DOperator(0,1)*lm(0,0) + DOperator(1,1)*lm(1,0));
    rhs[7]=-(DOperator(0,1)*lm(0,1) + DOperator(1,1)*lm(1,1));
    rhs[8]=-DOperator(0,0)*u1(0,0) - DOperator(0,1)*u1(1,0) + MOperator(0,0)*u2(0,0) + MOperator(0,1)*u2(1,0);
    rhs[9]=-DOperator(0,0)*u1(0,1) - DOperator(0,1)*u1(1,1) + MOperator(0,0)*u2(0,1) + MOperator(0,1)*u2(1,1);
    rhs[10]=-DOperator(1,0)*u1(0,0) - DOperator(1,1)*u1(1,0) + MOperator(1,0)*u2(0,0) + MOperator(1,1)*u2(1,0);
    rhs[11]=-DOperator(1,0)*u1(0,1) - DOperator(1,1)*u1(1,1) + MOperator(1,0)*u2(0,1) + MOperator(1,1)*u2(1,1);


    return rhs;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
template<>
array_1d<double, 6> MeshTyingMortarCondition<2,4,ScalarValue>::CalculateLocalRHS<6>(
    const MortarConditionMatrices& rMortarConditionMatrices,
    DofData& rDofData,
    const unsigned int rMasterElementIndex,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    array_1d<double,6> rhs;

    // Initialize values
    const bounded_matrix<double, 2, ScalarValue> u1 = rDofData.u1;
    const bounded_matrix<double, 2, ScalarValue> u2 = rDofData.u2;

    const bounded_matrix<double, 2, ScalarValue> lm = rDofData.LagrangeMultipliers; 

    // Mortar operators
    const bounded_matrix<double, 2, 2>& MOperator = rMortarConditionMatrices.MOperator;
    const bounded_matrix<double, 2, 2>& DOperator = rMortarConditionMatrices.DOperator;


    rhs[0]=MOperator(0,0)*lm(0,0) + MOperator(1,0)*lm(1,0);
    rhs[1]=MOperator(0,1)*lm(0,0) + MOperator(1,1)*lm(1,0);
    rhs[2]=-(DOperator(0,0)*lm(0,0) + DOperator(1,0)*lm(1,0));
    rhs[3]=-(DOperator(0,1)*lm(0,0) + DOperator(1,1)*lm(1,0));
    rhs[4]=-DOperator(0,0)*u1(0,0) - DOperator(0,1)*u1(1,0) + MOperator(0,0)*u2(0,0) + MOperator(0,1)*u2(1,0);
    rhs[5]=-DOperator(1,0)*u1(0,0) - DOperator(1,1)*u1(1,0) + MOperator(1,0)*u2(0,0) + MOperator(1,1)*u2(1,0);


    return rhs;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
template<>
array_1d<double, 12> MeshTyingMortarCondition<2,4,Vector2DValue>::CalculateLocalRHS<12>(
    const MortarConditionMatrices& rMortarConditionMatrices,
    DofData& rDofData,
    const unsigned int rMasterElementIndex,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    array_1d<double,12> rhs;

    // Initialize values
    const bounded_matrix<double, 2, Vector2DValue> u1 = rDofData.u1;
    const bounded_matrix<double, 2, Vector2DValue> u2 = rDofData.u2;

    const bounded_matrix<double, 2, Vector2DValue> lm = rDofData.LagrangeMultipliers; 

    // Mortar operators
    const bounded_matrix<double, 2, 2>& MOperator = rMortarConditionMatrices.MOperator;
    const bounded_matrix<double, 2, 2>& DOperator = rMortarConditionMatrices.DOperator;


    rhs[0]=MOperator(0,0)*lm(0,0) + MOperator(1,0)*lm(1,0);
    rhs[1]=MOperator(0,0)*lm(0,1) + MOperator(1,0)*lm(1,1);
    rhs[2]=MOperator(0,1)*lm(0,0) + MOperator(1,1)*lm(1,0);
    rhs[3]=MOperator(0,1)*lm(0,1) + MOperator(1,1)*lm(1,1);
    rhs[4]=-(DOperator(0,0)*lm(0,0) + DOperator(1,0)*lm(1,0));
    rhs[5]=-(DOperator(0,0)*lm(0,1) + DOperator(1,0)*lm(1,1));
    rhs[6]=-(DOperator(0,1)*lm(0,0) + DOperator(1,1)*lm(1,0));
    rhs[7]=-(DOperator(0,1)*lm(0,1) + DOperator(1,1)*lm(1,1));
    rhs[8]=-DOperator(0,0)*u1(0,0) - DOperator(0,1)*u1(1,0) + MOperator(0,0)*u2(0,0) + MOperator(0,1)*u2(1,0);
    rhs[9]=-DOperator(0,0)*u1(0,1) - DOperator(0,1)*u1(1,1) + MOperator(0,0)*u2(0,1) + MOperator(0,1)*u2(1,1);
    rhs[10]=-DOperator(1,0)*u1(0,0) - DOperator(1,1)*u1(1,0) + MOperator(1,0)*u2(0,0) + MOperator(1,1)*u2(1,0);
    rhs[11]=-DOperator(1,0)*u1(0,1) - DOperator(1,1)*u1(1,1) + MOperator(1,0)*u2(0,1) + MOperator(1,1)*u2(1,1);


    return rhs;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
template<>
array_1d<double, 9> MeshTyingMortarCondition<3,4,ScalarValue>::CalculateLocalRHS<9>(
    const MortarConditionMatrices& rMortarConditionMatrices,
    DofData& rDofData,
    const unsigned int rMasterElementIndex,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    array_1d<double,9> rhs;

    // Initialize values
    const bounded_matrix<double, 3, ScalarValue> u1 = rDofData.u1;
    const bounded_matrix<double, 3, ScalarValue> u2 = rDofData.u2;

    const bounded_matrix<double, 3, ScalarValue> lm = rDofData.LagrangeMultipliers; 

    // Mortar operators
    const bounded_matrix<double, 3, 3>& MOperator = rMortarConditionMatrices.MOperator;
    const bounded_matrix<double, 3, 3>& DOperator = rMortarConditionMatrices.DOperator;


    rhs[0]=MOperator(0,0)*lm(0,0) + MOperator(1,0)*lm(1,0) + MOperator(2,0)*lm(2,0);
    rhs[1]=MOperator(0,1)*lm(0,0) + MOperator(1,1)*lm(1,0) + MOperator(2,1)*lm(2,0);
    rhs[2]=MOperator(0,2)*lm(0,0) + MOperator(1,2)*lm(1,0) + MOperator(2,2)*lm(2,0);
    rhs[3]=-(DOperator(0,0)*lm(0,0) + DOperator(1,0)*lm(1,0) + DOperator(2,0)*lm(2,0));
    rhs[4]=-(DOperator(0,1)*lm(0,0) + DOperator(1,1)*lm(1,0) + DOperator(2,1)*lm(2,0));
    rhs[5]=-(DOperator(0,2)*lm(0,0) + DOperator(1,2)*lm(1,0) + DOperator(2,2)*lm(2,0));
    rhs[6]=-DOperator(0,0)*u1(0,0) - DOperator(0,1)*u1(1,0) - DOperator(0,2)*u1(2,0) + MOperator(0,0)*u2(0,0) + MOperator(0,1)*u2(1,0) + MOperator(0,2)*u2(2,0);
    rhs[7]=-DOperator(1,0)*u1(0,0) - DOperator(1,1)*u1(1,0) - DOperator(1,2)*u1(2,0) + MOperator(1,0)*u2(0,0) + MOperator(1,1)*u2(1,0) + MOperator(1,2)*u2(2,0);
    rhs[8]=-DOperator(2,0)*u1(0,0) - DOperator(2,1)*u1(1,0) - DOperator(2,2)*u1(2,0) + MOperator(2,0)*u2(0,0) + MOperator(2,1)*u2(1,0) + MOperator(2,2)*u2(2,0);


    return rhs;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
template<>
array_1d<double, 27> MeshTyingMortarCondition<3,4,Vector3DValue>::CalculateLocalRHS<27>(
    const MortarConditionMatrices& rMortarConditionMatrices,
    DofData& rDofData,
    const unsigned int rMasterElementIndex,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    array_1d<double,27> rhs;

    // Initialize values
    const bounded_matrix<double, 3, Vector3DValue> u1 = rDofData.u1;
    const bounded_matrix<double, 3, Vector3DValue> u2 = rDofData.u2;

    const bounded_matrix<double, 3, Vector3DValue> lm = rDofData.LagrangeMultipliers; 

    // Mortar operators
    const bounded_matrix<double, 3, 3>& MOperator = rMortarConditionMatrices.MOperator;
    const bounded_matrix<double, 3, 3>& DOperator = rMortarConditionMatrices.DOperator;


    rhs[0]=MOperator(0,0)*lm(0,0) + MOperator(1,0)*lm(1,0) + MOperator(2,0)*lm(2,0);
    rhs[1]=MOperator(0,0)*lm(0,1) + MOperator(1,0)*lm(1,1) + MOperator(2,0)*lm(2,1);
    rhs[2]=MOperator(0,0)*lm(0,2) + MOperator(1,0)*lm(1,2) + MOperator(2,0)*lm(2,2);
    rhs[3]=MOperator(0,1)*lm(0,0) + MOperator(1,1)*lm(1,0) + MOperator(2,1)*lm(2,0);
    rhs[4]=MOperator(0,1)*lm(0,1) + MOperator(1,1)*lm(1,1) + MOperator(2,1)*lm(2,1);
    rhs[5]=MOperator(0,1)*lm(0,2) + MOperator(1,1)*lm(1,2) + MOperator(2,1)*lm(2,2);
    rhs[6]=MOperator(0,2)*lm(0,0) + MOperator(1,2)*lm(1,0) + MOperator(2,2)*lm(2,0);
    rhs[7]=MOperator(0,2)*lm(0,1) + MOperator(1,2)*lm(1,1) + MOperator(2,2)*lm(2,1);
    rhs[8]=MOperator(0,2)*lm(0,2) + MOperator(1,2)*lm(1,2) + MOperator(2,2)*lm(2,2);
    rhs[9]=-(DOperator(0,0)*lm(0,0) + DOperator(1,0)*lm(1,0) + DOperator(2,0)*lm(2,0));
    rhs[10]=-(DOperator(0,0)*lm(0,1) + DOperator(1,0)*lm(1,1) + DOperator(2,0)*lm(2,1));
    rhs[11]=-(DOperator(0,0)*lm(0,2) + DOperator(1,0)*lm(1,2) + DOperator(2,0)*lm(2,2));
    rhs[12]=-(DOperator(0,1)*lm(0,0) + DOperator(1,1)*lm(1,0) + DOperator(2,1)*lm(2,0));
    rhs[13]=-(DOperator(0,1)*lm(0,1) + DOperator(1,1)*lm(1,1) + DOperator(2,1)*lm(2,1));
    rhs[14]=-(DOperator(0,1)*lm(0,2) + DOperator(1,1)*lm(1,2) + DOperator(2,1)*lm(2,2));
    rhs[15]=-(DOperator(0,2)*lm(0,0) + DOperator(1,2)*lm(1,0) + DOperator(2,2)*lm(2,0));
    rhs[16]=-(DOperator(0,2)*lm(0,1) + DOperator(1,2)*lm(1,1) + DOperator(2,2)*lm(2,1));
    rhs[17]=-(DOperator(0,2)*lm(0,2) + DOperator(1,2)*lm(1,2) + DOperator(2,2)*lm(2,2));
    rhs[18]=-DOperator(0,0)*u1(0,0) - DOperator(0,1)*u1(1,0) - DOperator(0,2)*u1(2,0) + MOperator(0,0)*u2(0,0) + MOperator(0,1)*u2(1,0) + MOperator(0,2)*u2(2,0);
    rhs[19]=-DOperator(0,0)*u1(0,1) - DOperator(0,1)*u1(1,1) - DOperator(0,2)*u1(2,1) + MOperator(0,0)*u2(0,1) + MOperator(0,1)*u2(1,1) + MOperator(0,2)*u2(2,1);
    rhs[20]=-DOperator(0,0)*u1(0,2) - DOperator(0,1)*u1(1,2) - DOperator(0,2)*u1(2,2) + MOperator(0,0)*u2(0,2) + MOperator(0,1)*u2(1,2) + MOperator(0,2)*u2(2,2);
    rhs[21]=-DOperator(1,0)*u1(0,0) - DOperator(1,1)*u1(1,0) - DOperator(1,2)*u1(2,0) + MOperator(1,0)*u2(0,0) + MOperator(1,1)*u2(1,0) + MOperator(1,2)*u2(2,0);
    rhs[22]=-DOperator(1,0)*u1(0,1) - DOperator(1,1)*u1(1,1) - DOperator(1,2)*u1(2,1) + MOperator(1,0)*u2(0,1) + MOperator(1,1)*u2(1,1) + MOperator(1,2)*u2(2,1);
    rhs[23]=-DOperator(1,0)*u1(0,2) - DOperator(1,1)*u1(1,2) - DOperator(1,2)*u1(2,2) + MOperator(1,0)*u2(0,2) + MOperator(1,1)*u2(1,2) + MOperator(1,2)*u2(2,2);
    rhs[24]=-DOperator(2,0)*u1(0,0) - DOperator(2,1)*u1(1,0) - DOperator(2,2)*u1(2,0) + MOperator(2,0)*u2(0,0) + MOperator(2,1)*u2(1,0) + MOperator(2,2)*u2(2,0);
    rhs[25]=-DOperator(2,0)*u1(0,1) - DOperator(2,1)*u1(1,1) - DOperator(2,2)*u1(2,1) + MOperator(2,0)*u2(0,1) + MOperator(2,1)*u2(1,1) + MOperator(2,2)*u2(2,1);
    rhs[26]=-DOperator(2,0)*u1(0,2) - DOperator(2,1)*u1(1,2) - DOperator(2,2)*u1(2,2) + MOperator(2,0)*u2(0,2) + MOperator(2,1)*u2(1,2) + MOperator(2,2)*u2(2,2);


    return rhs;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
template<>
array_1d<double, 12> MeshTyingMortarCondition<3,8,ScalarValue>::CalculateLocalRHS<12>(
    const MortarConditionMatrices& rMortarConditionMatrices,
    DofData& rDofData,
    const unsigned int rMasterElementIndex,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    array_1d<double,12> rhs;

    // Initialize values
    const bounded_matrix<double, 4, ScalarValue> u1 = rDofData.u1;
    const bounded_matrix<double, 4, ScalarValue> u2 = rDofData.u2;

    const bounded_matrix<double, 4, ScalarValue> lm = rDofData.LagrangeMultipliers; 

    // Mortar operators
    const bounded_matrix<double, 4, 4>& MOperator = rMortarConditionMatrices.MOperator;
    const bounded_matrix<double, 4, 4>& DOperator = rMortarConditionMatrices.DOperator;


    rhs[0]=MOperator(0,0)*lm(0,0) + MOperator(1,0)*lm(1,0) + MOperator(2,0)*lm(2,0) + MOperator(3,0)*lm(3,0);
    rhs[1]=MOperator(0,1)*lm(0,0) + MOperator(1,1)*lm(1,0) + MOperator(2,1)*lm(2,0) + MOperator(3,1)*lm(3,0);
    rhs[2]=MOperator(0,2)*lm(0,0) + MOperator(1,2)*lm(1,0) + MOperator(2,2)*lm(2,0) + MOperator(3,2)*lm(3,0);
    rhs[3]=MOperator(0,3)*lm(0,0) + MOperator(1,3)*lm(1,0) + MOperator(2,3)*lm(2,0) + MOperator(3,3)*lm(3,0);
    rhs[4]=-(DOperator(0,0)*lm(0,0) + DOperator(1,0)*lm(1,0) + DOperator(2,0)*lm(2,0) + DOperator(3,0)*lm(3,0));
    rhs[5]=-(DOperator(0,1)*lm(0,0) + DOperator(1,1)*lm(1,0) + DOperator(2,1)*lm(2,0) + DOperator(3,1)*lm(3,0));
    rhs[6]=-(DOperator(0,2)*lm(0,0) + DOperator(1,2)*lm(1,0) + DOperator(2,2)*lm(2,0) + DOperator(3,2)*lm(3,0));
    rhs[7]=-(DOperator(0,3)*lm(0,0) + DOperator(1,3)*lm(1,0) + DOperator(2,3)*lm(2,0) + DOperator(3,3)*lm(3,0));
    rhs[8]=-DOperator(0,0)*u1(0,0) - DOperator(0,1)*u1(1,0) - DOperator(0,2)*u1(2,0) - DOperator(0,3)*u1(3,0) + MOperator(0,0)*u2(0,0) + MOperator(0,1)*u2(1,0) + MOperator(0,2)*u2(2,0) + MOperator(0,3)*u2(3,0);
    rhs[9]=-DOperator(1,0)*u1(0,0) - DOperator(1,1)*u1(1,0) - DOperator(1,2)*u1(2,0) - DOperator(1,3)*u1(3,0) + MOperator(1,0)*u2(0,0) + MOperator(1,1)*u2(1,0) + MOperator(1,2)*u2(2,0) + MOperator(1,3)*u2(3,0);
    rhs[10]=-DOperator(2,0)*u1(0,0) - DOperator(2,1)*u1(1,0) - DOperator(2,2)*u1(2,0) - DOperator(2,3)*u1(3,0) + MOperator(2,0)*u2(0,0) + MOperator(2,1)*u2(1,0) + MOperator(2,2)*u2(2,0) + MOperator(2,3)*u2(3,0);
    rhs[11]=-DOperator(3,0)*u1(0,0) - DOperator(3,1)*u1(1,0) - DOperator(3,2)*u1(2,0) - DOperator(3,3)*u1(3,0) + MOperator(3,0)*u2(0,0) + MOperator(3,1)*u2(1,0) + MOperator(3,2)*u2(2,0) + MOperator(3,3)*u2(3,0);


    return rhs;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
template<>
array_1d<double, 36> MeshTyingMortarCondition<3,8,Vector3DValue>::CalculateLocalRHS<36>(
    const MortarConditionMatrices& rMortarConditionMatrices,
    DofData& rDofData,
    const unsigned int rMasterElementIndex,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    array_1d<double,36> rhs;

    // Initialize values
    const bounded_matrix<double, 4, Vector3DValue> u1 = rDofData.u1;
    const bounded_matrix<double, 4, Vector3DValue> u2 = rDofData.u2;

    const bounded_matrix<double, 4, Vector3DValue> lm = rDofData.LagrangeMultipliers; 

    // Mortar operators
    const bounded_matrix<double, 4, 4>& MOperator = rMortarConditionMatrices.MOperator;
    const bounded_matrix<double, 4, 4>& DOperator = rMortarConditionMatrices.DOperator;


    rhs[0]=MOperator(0,0)*lm(0,0) + MOperator(1,0)*lm(1,0) + MOperator(2,0)*lm(2,0) + MOperator(3,0)*lm(3,0);
    rhs[1]=MOperator(0,0)*lm(0,1) + MOperator(1,0)*lm(1,1) + MOperator(2,0)*lm(2,1) + MOperator(3,0)*lm(3,1);
    rhs[2]=MOperator(0,0)*lm(0,2) + MOperator(1,0)*lm(1,2) + MOperator(2,0)*lm(2,2) + MOperator(3,0)*lm(3,2);
    rhs[3]=MOperator(0,1)*lm(0,0) + MOperator(1,1)*lm(1,0) + MOperator(2,1)*lm(2,0) + MOperator(3,1)*lm(3,0);
    rhs[4]=MOperator(0,1)*lm(0,1) + MOperator(1,1)*lm(1,1) + MOperator(2,1)*lm(2,1) + MOperator(3,1)*lm(3,1);
    rhs[5]=MOperator(0,1)*lm(0,2) + MOperator(1,1)*lm(1,2) + MOperator(2,1)*lm(2,2) + MOperator(3,1)*lm(3,2);
    rhs[6]=MOperator(0,2)*lm(0,0) + MOperator(1,2)*lm(1,0) + MOperator(2,2)*lm(2,0) + MOperator(3,2)*lm(3,0);
    rhs[7]=MOperator(0,2)*lm(0,1) + MOperator(1,2)*lm(1,1) + MOperator(2,2)*lm(2,1) + MOperator(3,2)*lm(3,1);
    rhs[8]=MOperator(0,2)*lm(0,2) + MOperator(1,2)*lm(1,2) + MOperator(2,2)*lm(2,2) + MOperator(3,2)*lm(3,2);
    rhs[9]=MOperator(0,3)*lm(0,0) + MOperator(1,3)*lm(1,0) + MOperator(2,3)*lm(2,0) + MOperator(3,3)*lm(3,0);
    rhs[10]=MOperator(0,3)*lm(0,1) + MOperator(1,3)*lm(1,1) + MOperator(2,3)*lm(2,1) + MOperator(3,3)*lm(3,1);
    rhs[11]=MOperator(0,3)*lm(0,2) + MOperator(1,3)*lm(1,2) + MOperator(2,3)*lm(2,2) + MOperator(3,3)*lm(3,2);
    rhs[12]=-(DOperator(0,0)*lm(0,0) + DOperator(1,0)*lm(1,0) + DOperator(2,0)*lm(2,0) + DOperator(3,0)*lm(3,0));
    rhs[13]=-(DOperator(0,0)*lm(0,1) + DOperator(1,0)*lm(1,1) + DOperator(2,0)*lm(2,1) + DOperator(3,0)*lm(3,1));
    rhs[14]=-(DOperator(0,0)*lm(0,2) + DOperator(1,0)*lm(1,2) + DOperator(2,0)*lm(2,2) + DOperator(3,0)*lm(3,2));
    rhs[15]=-(DOperator(0,1)*lm(0,0) + DOperator(1,1)*lm(1,0) + DOperator(2,1)*lm(2,0) + DOperator(3,1)*lm(3,0));
    rhs[16]=-(DOperator(0,1)*lm(0,1) + DOperator(1,1)*lm(1,1) + DOperator(2,1)*lm(2,1) + DOperator(3,1)*lm(3,1));
    rhs[17]=-(DOperator(0,1)*lm(0,2) + DOperator(1,1)*lm(1,2) + DOperator(2,1)*lm(2,2) + DOperator(3,1)*lm(3,2));
    rhs[18]=-(DOperator(0,2)*lm(0,0) + DOperator(1,2)*lm(1,0) + DOperator(2,2)*lm(2,0) + DOperator(3,2)*lm(3,0));
    rhs[19]=-(DOperator(0,2)*lm(0,1) + DOperator(1,2)*lm(1,1) + DOperator(2,2)*lm(2,1) + DOperator(3,2)*lm(3,1));
    rhs[20]=-(DOperator(0,2)*lm(0,2) + DOperator(1,2)*lm(1,2) + DOperator(2,2)*lm(2,2) + DOperator(3,2)*lm(3,2));
    rhs[21]=-(DOperator(0,3)*lm(0,0) + DOperator(1,3)*lm(1,0) + DOperator(2,3)*lm(2,0) + DOperator(3,3)*lm(3,0));
    rhs[22]=-(DOperator(0,3)*lm(0,1) + DOperator(1,3)*lm(1,1) + DOperator(2,3)*lm(2,1) + DOperator(3,3)*lm(3,1));
    rhs[23]=-(DOperator(0,3)*lm(0,2) + DOperator(1,3)*lm(1,2) + DOperator(2,3)*lm(2,2) + DOperator(3,3)*lm(3,2));
    rhs[24]=-DOperator(0,0)*u1(0,0) - DOperator(0,1)*u1(1,0) - DOperator(0,2)*u1(2,0) - DOperator(0,3)*u1(3,0) + MOperator(0,0)*u2(0,0) + MOperator(0,1)*u2(1,0) + MOperator(0,2)*u2(2,0) + MOperator(0,3)*u2(3,0);
    rhs[25]=-DOperator(0,0)*u1(0,1) - DOperator(0,1)*u1(1,1) - DOperator(0,2)*u1(2,1) - DOperator(0,3)*u1(3,1) + MOperator(0,0)*u2(0,1) + MOperator(0,1)*u2(1,1) + MOperator(0,2)*u2(2,1) + MOperator(0,3)*u2(3,1);
    rhs[26]=-DOperator(0,0)*u1(0,2) - DOperator(0,1)*u1(1,2) - DOperator(0,2)*u1(2,2) - DOperator(0,3)*u1(3,2) + MOperator(0,0)*u2(0,2) + MOperator(0,1)*u2(1,2) + MOperator(0,2)*u2(2,2) + MOperator(0,3)*u2(3,2);
    rhs[27]=-DOperator(1,0)*u1(0,0) - DOperator(1,1)*u1(1,0) - DOperator(1,2)*u1(2,0) - DOperator(1,3)*u1(3,0) + MOperator(1,0)*u2(0,0) + MOperator(1,1)*u2(1,0) + MOperator(1,2)*u2(2,0) + MOperator(1,3)*u2(3,0);
    rhs[28]=-DOperator(1,0)*u1(0,1) - DOperator(1,1)*u1(1,1) - DOperator(1,2)*u1(2,1) - DOperator(1,3)*u1(3,1) + MOperator(1,0)*u2(0,1) + MOperator(1,1)*u2(1,1) + MOperator(1,2)*u2(2,1) + MOperator(1,3)*u2(3,1);
    rhs[29]=-DOperator(1,0)*u1(0,2) - DOperator(1,1)*u1(1,2) - DOperator(1,2)*u1(2,2) - DOperator(1,3)*u1(3,2) + MOperator(1,0)*u2(0,2) + MOperator(1,1)*u2(1,2) + MOperator(1,2)*u2(2,2) + MOperator(1,3)*u2(3,2);
    rhs[30]=-DOperator(2,0)*u1(0,0) - DOperator(2,1)*u1(1,0) - DOperator(2,2)*u1(2,0) - DOperator(2,3)*u1(3,0) + MOperator(2,0)*u2(0,0) + MOperator(2,1)*u2(1,0) + MOperator(2,2)*u2(2,0) + MOperator(2,3)*u2(3,0);
    rhs[31]=-DOperator(2,0)*u1(0,1) - DOperator(2,1)*u1(1,1) - DOperator(2,2)*u1(2,1) - DOperator(2,3)*u1(3,1) + MOperator(2,0)*u2(0,1) + MOperator(2,1)*u2(1,1) + MOperator(2,2)*u2(2,1) + MOperator(2,3)*u2(3,1);
    rhs[32]=-DOperator(2,0)*u1(0,2) - DOperator(2,1)*u1(1,2) - DOperator(2,2)*u1(2,2) - DOperator(2,3)*u1(3,2) + MOperator(2,0)*u2(0,2) + MOperator(2,1)*u2(1,2) + MOperator(2,2)*u2(2,2) + MOperator(2,3)*u2(3,2);
    rhs[33]=-DOperator(3,0)*u1(0,0) - DOperator(3,1)*u1(1,0) - DOperator(3,2)*u1(2,0) - DOperator(3,3)*u1(3,0) + MOperator(3,0)*u2(0,0) + MOperator(3,1)*u2(1,0) + MOperator(3,2)*u2(2,0) + MOperator(3,3)*u2(3,0);
    rhs[34]=-DOperator(3,0)*u1(0,1) - DOperator(3,1)*u1(1,1) - DOperator(3,2)*u1(2,1) - DOperator(3,3)*u1(3,1) + MOperator(3,0)*u2(0,1) + MOperator(3,1)*u2(1,1) + MOperator(3,2)*u2(2,1) + MOperator(3,3)*u2(3,1);
    rhs[35]=-DOperator(3,0)*u1(0,2) - DOperator(3,1)*u1(1,2) - DOperator(3,2)*u1(2,2) - DOperator(3,3)*u1(3,2) + MOperator(3,0)*u2(0,2) + MOperator(3,1)*u2(1,2) + MOperator(3,2)*u2(2,2) + MOperator(3,3)*u2(3,2);


    return rhs;
}


/****************************** END AD REPLACEMENT *********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& CurrentProcessInfo 
    )
{
    KRATOS_TRY;  
    
    // Calculates the size of the system
    const unsigned int ConditionSize = MatrixSize * mPairSize; 
    
    if (rResult.size() != ConditionSize)
    {
        rResult.resize( ConditionSize, false );
    }
    
    unsigned int index = 0;
    
    /* ORDER - [ MASTER, SLAVE, LM ] */
    for ( unsigned int i_cond = 0;  i_cond < mPairSize; ++i_cond )
    {   
        // Master Nodes DoF Equation IDs
        GeometryType& current_master = mThisMasterConditions[i_cond]->GetGeometry( );
        
        if (TTensor == ScalarValue)
        {
            for ( unsigned int i_master = 0; i_master < NumNodes; i_master++ ) 
            {
                NodeType& master_node = current_master[i_master];
                rResult[index++] = master_node.GetDof( TEMPERATURE ).EquationId( );
//                 rResult[index++] = master_node.GetDof( mTyingVarScalar ).EquationId( );
            }
        }
        else
        {
            for ( unsigned int i_master = 0; i_master < NumNodes; i_master++ ) 
            {
                NodeType& master_node = current_master[i_master];
                rResult[index++] = master_node.GetDof( DISPLACEMENT_X ).EquationId( );
                rResult[index++] = master_node.GetDof( DISPLACEMENT_Y ).EquationId( );
                if (TDim == 3)
                {
                    rResult[index++] = master_node.GetDof( DISPLACEMENT_Z ).EquationId( );
                }
//                 for (unsigned int i_dof = 0; i_dof < TDim; i_dof++)
//                 {
//                     rResult[index++] = master_node.GetDof( mTyingVarVector[i_dof] ).EquationId( );
//                 }
            }
        }
        
        // Slave Nodes DoF Equation IDs
        if (TTensor == ScalarValue)
        {
            for ( unsigned int i_slave = 0; i_slave < NumNodes; i_slave++ ) 
            {
                NodeType& slave_node = this->GetGeometry()[i_slave];
                rResult[index++] = slave_node.GetDof( TEMPERATURE ).EquationId( );
//                 rResult[index++] = slave_node.GetDof( mTyingVarScalar ).EquationId( );
            }
        }
        else
        {
            for ( unsigned int i_slave = 0; i_slave < NumNodes; i_slave++ ) 
            {
                NodeType& slave_node = this->GetGeometry()[i_slave];
                rResult[index++] = slave_node.GetDof( DISPLACEMENT_X ).EquationId( );
                rResult[index++] = slave_node.GetDof( DISPLACEMENT_Y ).EquationId( );
                if (TDim == 3)
                {
                    rResult[index++] = slave_node.GetDof( DISPLACEMENT_Z ).EquationId( );
                }
//                 for (unsigned int i_dof = 0; i_dof < TDim; i_dof++)
//                 {
//                     rResult[index++] = slave_node.GetDof( mTyingVarVector[i_dof] ).EquationId( );
//                 }
            }
        }
        
        // Slave Nodes LM Equation IDs
        if (TTensor == ScalarValue)
        {
            for ( unsigned int i_slave = 0; i_slave < NumNodes; i_slave++ ) 
            {
                NodeType& slave_node = this->GetGeometry()[i_slave];
                rResult[index++] = slave_node.GetDof( SCALAR_LAGRANGE_MULTIPLIER ).EquationId( );
            }
        }
        else
        {
            for ( unsigned int i_slave = 0; i_slave < NumNodes; i_slave++ ) 
            {
                NodeType& slave_node = this->GetGeometry()[i_slave];
                rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_X ).EquationId( );
                rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Y ).EquationId( );
                if (TDim == 3)
                {
                    rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Z ).EquationId( );
                }
            }
        }
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim, TNumNodesElem, TTensor>::GetDofList(
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo 
)
{
    KRATOS_TRY;
    
    // TODO: You need the utility to get the dof
    
    // Calculates the size of the system
    const unsigned int ConditionSize = MatrixSize * mPairSize; 
    
    if (rConditionalDofList.size() != ConditionSize)
    {
        rConditionalDofList.resize( ConditionSize );
    }
    
    unsigned int index = 0;
    
    /* ORDER - [ MASTER, SLAVE, LM ] */
    for ( unsigned int i_cond = 0;  i_cond < mPairSize; ++i_cond )
    {   
        // Master Nodes DoF Equation IDs
        GeometryType& current_master = mThisMasterConditions[i_cond]->GetGeometry( );
        
        if (TTensor == ScalarValue)
        {
            for ( unsigned int i_master = 0; i_master < NumNodes; i_master++ ) 
            {
                NodeType& master_node = current_master[i_master];
                rConditionalDofList[index++] = master_node.pGetDof( TEMPERATURE );
//                 rConditionalDofList[index++] = master_node.pGetDof( mTyingVarScalar );
            }
        }
        else
        {
            for ( unsigned int i_master = 0; i_master < NumNodes; i_master++ ) 
            {
                NodeType& master_node = current_master[i_master];
                rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_X );
                rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_Y );
                if (TDim == 3)
                {
                    rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_Z );
                }
//                 for (unsigned int i_dof = 0; i_dof < TDim; i_dof++)
//                 {
//                     rConditionalDofList[index++] = master_node.pGetDof( mTyingVarVector[i_dof] );
//                 }
            }
        }
        
        // Slave Nodes DoF Equation IDs
        if (TTensor == ScalarValue)
        {
            for ( unsigned int i_slave = 0; i_slave < NumNodes; i_slave++ ) 
            {
                NodeType& slave_node = this->GetGeometry()[i_slave];
                rConditionalDofList[index++] = slave_node.pGetDof( TEMPERATURE );
//                 rConditionalDofList[index++] = slave_node.pGetDof( mTyingVarScalar );
            }
        }
        else
        {
            for ( unsigned int i_slave = 0; i_slave < NumNodes; i_slave++ ) 
            {
                NodeType& slave_node = this->GetGeometry()[i_slave];
                rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_X );
                rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_Y );
                if (TDim == 3)
                {
                    rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_Z );
                }
//                 for (unsigned int i_dof = 0; i_dof < TDim; i_dof++)
//                 {
//                     rConditionalDofList[index++] = slave_node.pGetDof( mTyingVarVector[i_dof] );
//                 }
            }
        }
        
        // Slave Nodes LM Equation IDs
        if (TTensor == ScalarValue)
        {
            for ( unsigned int i_slave = 0; i_slave < NumNodes; i_slave++ ) 
            {
                NodeType& slave_node = this->GetGeometry()[i_slave];
                rConditionalDofList[index++] = slave_node.pGetDof( SCALAR_LAGRANGE_MULTIPLIER );
            }
        }
        else
        {
            for ( unsigned int i_slave = 0; i_slave < NumNodes; i_slave++ ) 
            {
                NodeType& slave_node = this->GetGeometry()[i_slave];
                rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_X );
                rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Y );
                if (TDim == 3)
                {
                    rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Z );
                }
            }
        }
    }
    
    KRATOS_CATCH( "" );
}


//******************************* GET DOUBLE VALUE *********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::GetValueOnIntegrationPoints( 
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//******************************* GET ARRAY_1D VALUE *******************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::GetValueOnIntegrationPoints( 
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//******************************* GET VECTOR VALUE *********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::GetValueOnIntegrationPoints( 
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateOnIntegrationPoints( 
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;

//     // Create and initialize condition variables:
//     GeneralVariables rVariables;
//     
//     // Initialize the current DoF data
//     DofData rDofData;
//                                                
//     const unsigned int number_of_integration_pts =IntegrationPointsSlave.size();
//     if ( rOutput.size( ) != number_of_integration_pts )
//     {
//         rOutput.resize( number_of_integration_pts, false );
//     }
//     
//     const std::vector<double> zero_vector (number_of_integration_pts, 0.0);
//     rOutput = zero_vector;
// 
    // TODO: Add eventually
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateOnIntegrationPoints( 
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector< array_1d<double, 3 > >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
    
//     // Create and initialize condition variables:
//     GeneralVariables rVariables;
//     
//     // Initialize the current contact data
//     DofData rDofData;
//                                                                                                                         
//     const unsigned int number_of_integration_pts = IntegrationPointsSlave.size();
//     if ( rOutput.size() != number_of_integration_pts )
//     {
//         rOutput.resize( number_of_integration_pts );
//     }
//     
//     const array_1d<double, 3> zero_vector = ZeroVector(3);
//     for (unsigned int PointNumber = 0; PointNumber < number_of_integration_pts; PointNumber++)
//     {
//         rOutput[PointNumber] = zero_vector;
//     }
    
    // TODO: Add eventually
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodesElem, TensorValue TTensor>
void MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateOnIntegrationPoints( 
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

template class MeshTyingMortarCondition<2, 3, ScalarValue>;   // 2DLine/Triangle for scalar variables
template class MeshTyingMortarCondition<2, 3, Vector2DValue>; // 2DLine/Triangle for components variables
template class MeshTyingMortarCondition<2, 4, ScalarValue>;   // 2DLine/Quadrilateral for scalar variables
template class MeshTyingMortarCondition<2, 4, Vector2DValue>; // 2DLine/Quadrilateral for scalar variables
template class MeshTyingMortarCondition<3, 4, ScalarValue>;   // 3D Triangle/Tetrahedron for scalar variables
template class MeshTyingMortarCondition<3, 4, Vector3DValue>; // 3D Triangle/Tetrahedron for components variables
template class MeshTyingMortarCondition<3, 8, ScalarValue>;   // 3D Quadrilateral/Hexahedra for scalar variables
template class MeshTyingMortarCondition<3, 8, Vector3DValue>; // 3D Quadrilateral/Hexahedra for components variables

} // Namespace Kratos
