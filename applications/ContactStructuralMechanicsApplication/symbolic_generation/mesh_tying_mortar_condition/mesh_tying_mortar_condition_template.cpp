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

// replace_lhs

/****************************** END AD REPLACEMENT *********************************/
/***********************************************************************************/

/***************************** BEGIN AD REPLACEMENT ********************************/
/***********************************************************************************/

// replace_rhs

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
