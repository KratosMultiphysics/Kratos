// KRATOS  ___|  |       |       |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//           | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License: BSD License
//   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrándiz
//

// System includes

// External includes

// Project includes
/* Mortar includes */
#include "custom_conditions/ALM_mortar_contact_condition.h"

/* Additional includes */
#include <algorithm>

/* Utilities */
#include "utilities/math_utils.h"
#include "custom_utilities/exact_mortar_segmentation_utility.h"

namespace Kratos 
{
/**
 * Flags related to the condition computation 
 */
// Avoiding using the macro since this has a template parameter. If there was no template plase use the KRATOS_CREATE_LOCAL_FLAG macro
template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
const Kratos::Flags AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_RHS_VECTOR(Kratos::Flags::Create(0));
template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
const Kratos::Flags AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX(Kratos::Flags::Create(1));
template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
const Kratos::Flags AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS(Kratos::Flags::Create(2));
template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
const Kratos::Flags AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS(Kratos::Flags::Create(3));

/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
Condition::Pointer AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::Create( 
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties ) const
{
    KRATOS_ERROR << "You are calling to the base class method Create, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
    
    return boost::make_shared< AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional> >( NewId, this->GetGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
Condition::Pointer AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_ERROR << "You are calling to the base class method Create, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
    
    return boost::make_shared< AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional> >( NewId, pGeom, pProperties );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::~AugmentedLagrangianMethodMortarContactCondition( )
{
}

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::Initialize( ) 
{
    KRATOS_TRY;
    
    mIntegrationOrder = GetProperties().GetValue(INTEGRATION_ORDER_CONTACT);
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // First populate of the vector of master elements
    boost::shared_ptr<ConditionMap>& AllConditionSets = this->GetValue( CONTACT_SETS );
    mPairSize = AllConditionSets->size();
    mThisMasterElements.resize( mPairSize );
    
    unsigned int iCond = 0;
    for (auto iPair = AllConditionSets->begin(); iPair != AllConditionSets->end(); ++iPair )
    {
        mThisMasterElements[iCond] = iPair->first;
        iCond += 1;
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // NOTE: Add things if necessary
        
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // NOTE: Add things if necessary
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // NOTE: Add things if necessary
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateLocalSystem( 
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
    LocalSystem.CalculationFlags.Set(AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS, true);
    LocalSystem.CalculationFlags.Set(AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS, true);

    //Initialize sizes for the system components
    if ( rLHSVariables.size( ) != rLeftHandSideMatrices.size( ) )
    {
        rLeftHandSideMatrices.resize( rLHSVariables.size( ) );
    }

    if ( rRHSVariables.size( ) != rRightHandSideVectors.size( ) )
    {
        rRightHandSideVectors.resize( rRHSVariables.size( ) );
    }

    LocalSystem.CalculationFlags.Set(AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX, true);
    for ( unsigned int i = 0; i < rLeftHandSideMatrices.size( ); i++ )
    {
        // Note: rRightHandSideVectors.size() > 0
        this->InitializeSystemMatrices( rLeftHandSideMatrices[i], rRightHandSideVectors[0],LocalSystem.CalculationFlags );
    }

    LocalSystem.CalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_RHS_VECTOR, true );
    LocalSystem.CalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX, false ); // Temporarily only
    for ( unsigned int i = 0; i < rRightHandSideVectors.size( ); i++ )
    {
        // Note: rLeftHandSideMatrices.size() > 0
        this->InitializeSystemMatrices( rLeftHandSideMatrices[0], rRightHandSideVectors[i], LocalSystem.CalculationFlags  );
    }
    LocalSystem.CalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX, true ); // Reactivated again

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

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;

    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX, true );
    LocalSystem.CalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_RHS_VECTOR, true );

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

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateLeftHandSide( 
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX, true );

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

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateLeftHandSide( 
    std::vector< MatrixType >& rLeftHandSideMatrices,
    const std::vector< Variable< MatrixType > >& rLHSVariables,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX, true );
    LocalSystem.CalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS, true );

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

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateRightHandSide( 
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_RHS_VECTOR, true);

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

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateRightHandSide( 
    std::vector< VectorType >& rRightHandSideVectors,
    const std::vector< Variable< VectorType > >& rRHSVariables,
    ProcessInfo& rCurrentProcessInfo )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_RHS_VECTOR, true );
    LocalSystem.CalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS, true );

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

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>

void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::InitializeSystemMatrices( 
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    Flags& rCalculationFlags 
    )
{
    const unsigned int ConditionSize = this->CalculateConditionSize( );
    
    // Resizing as needed the LHS
    if ( rCalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX ) ) // Calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != ConditionSize )
        {
            rLeftHandSideMatrix.resize( ConditionSize, ConditionSize, false );
        }
        noalias( rLeftHandSideMatrix ) = ZeroMatrix( ConditionSize, ConditionSize ); // Resetting LHS
    }

    // Resizing as needed the RHS
    if ( rCalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_RHS_VECTOR ) ) // Calculation of the matrix is required
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

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateMassMatrix( 
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

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateDampingMatrix( 
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

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>

const unsigned int AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateConditionSize( )
{
    const unsigned int ConditionSize = mPairSize * MatrixSize;
    
    return ConditionSize;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim, TNumNodes, TFrictional>::CalculateConditionSystem( 
    LocalSystemComponents& rLocalSystem,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
    
    // We get the ALM variables
    const double PenaltyFactor = rCurrentProcessInfo[PENALTY_FACTOR];
    const double ScaleFactor = rCurrentProcessInfo[SCALE_FACTOR];
    
    // Create and initialize condition variables:#pragma omp critical
    GeneralVariables rVariables;
    
    // Create the current contact data
    DerivativeDataType rDerivativeData;
    rDerivativeData.Initialize(this->GetGeometry());
    
    // Create the mortar operators
    MortarConditionMatrices rThisMortarConditionMatrices;
    
    // Compute Ae and its derivative
    this->CalculateAeAndDeltaAe(rDerivativeData, rVariables, rCurrentProcessInfo); 
    
    // We call the exact integration utility
    ExactMortarIntegrationUtility<TDim, TNumNodes> IntUtil = ExactMortarIntegrationUtility<TDim, TNumNodes>(mIntegrationOrder);
    
    // Iterate over the master segments
    for (unsigned int PairIndex = 0; PairIndex < mPairSize; ++PairIndex)
    {   
        // The normal of the master condition
        const array_1d<double, 3> MasterNormal = mThisMasterElements[PairIndex]->GetValue(NORMAL);
        
        // Reading integration points
        ConditionArrayListType ConditionsPointsSlave;
        const bool IsInside = IntUtil.GetExactIntegration(this->GetGeometry(), this->GetValue(NORMAL), mThisMasterElements[PairIndex]->GetGeometry(), mThisMasterElements[PairIndex]->GetValue(NORMAL), ConditionsPointsSlave);
        
        IntegrationMethod ThisIntegrationMethod = GetIntegrationMethod();
        
        if (IsInside == true)
        {            
            // Initialize general variables for the current master element
            this->InitializeGeneralVariables( rVariables, rCurrentProcessInfo, PairIndex );
            
            // Update slave element info
            rDerivativeData.UpdateMasterPair(mThisMasterElements[PairIndex] );
            
            // Initialize the mortar operators
            rThisMortarConditionMatrices.Initialize();
            
//             // Compute the normal derivatives of the master
//             this->CalculateDeltaNormalMaster(rDerivativeData);
            
            for (unsigned int iGeom = 0; iGeom < ConditionsPointsSlave.size(); iGeom++)
            {
                std::vector<PointType::Pointer> PointsArray (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
                for (unsigned int iNode = 0; iNode < TDim; iNode++)
                {
                    PointType GlobalPoint;
                    GetGeometry().GlobalCoordinates(GlobalPoint, ConditionsPointsSlave[iGeom][iNode]);
                    PointsArray[iNode] = boost::make_shared<PointType>(GlobalPoint);
                }
                
                DecompositionType DecompGeom( PointsArray );
                
                const bool BadShape = (TDim == 2) ? false : ContactUtilities::HeronCheck(DecompGeom);
                
                if (BadShape == false)
                {
                    const GeometryType::IntegrationPointsArrayType& IntegrationPointsSlave = DecompGeom.IntegrationPoints( ThisIntegrationMethod );
                    
                    // Integrating the mortar operators
                    for ( unsigned int PointNumber = 0; PointNumber < IntegrationPointsSlave.size(); PointNumber++ )
                    {
                        const PointType LocalPointDecomp = IntegrationPointsSlave[PointNumber].Coordinates();
                        PointType LocalPointParent;
                        PointType GPGlobal;
                        DecompGeom.GlobalCoordinates(GPGlobal, LocalPointDecomp);
                        GetGeometry().PointLocalCoordinates(LocalPointParent, GPGlobal);
                        
                        // Calculate the kinematic variables
                        this->CalculateKinematics( rVariables, rDerivativeData, MasterNormal, LocalPointDecomp, LocalPointParent, DecompGeom);
                        
                        const double IntegrationWeight = IntegrationPointsSlave[PointNumber].Weight();
                    
                        if ( rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX ) ||
                                rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
                        {
                            /* Update the derivatives */
                            // Update the derivative of DetJ
                            this->CalculateDeltaDetJSlave(rVariables, rDerivativeData);
                            // Update the derivatives of the shape functions and the gap
//                             this->CalculateDeltaN(rVariables, rDerivativeData); // FIXME: This is the old version!!!!
                            // The derivatives of the dual shape function 
                            this->CalculateDeltaPhi(rVariables, rDerivativeData);
                            
                            this->CalculateMortarOperators(rThisMortarConditionMatrices, rVariables, rDerivativeData, IntegrationWeight);    
                        }
                        else // In case we are computing RHS we don't compute derivatives (not necessary)
                        {
                            this->CalculateMortarOperators(rThisMortarConditionMatrices, rVariables, IntegrationWeight);   
                        }
                    }
                }
            }
                    
//             // Debug
//             std::cout << "--------------------------------------------------" << std::endl;
//             KRATOS_WATCH(this->Id());
//             KRATOS_WATCH(PairIndex);
//             rThisMortarConditionMatrices.print();
            
            // Calculates the active/inactive combination pair
            const unsigned int ActiveInactive = GetActiveInactiveValue(this->GetGeometry());
            
            // Assemble of the matrix is required
            if ( rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX ) ||
                    rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
            {
                // Calculate the local contribution
                const bounded_matrix<double, MatrixSize, MatrixSize> LHS_contact_pair = this->CalculateLocalLHS( rThisMortarConditionMatrices, PairIndex, PenaltyFactor, ScaleFactor, ActiveInactive);
                
                // Contributions to stiffness matrix calculated on the reference config
                this->CalculateAndAddLHS( rLocalSystem, LHS_contact_pair, PairIndex );
                
//                 // Debug
//     //                 KRATOS_WATCH(LHS_contact_pair);
//                 LOG_MATRIX_PRETTY( LHS_contact_pair );
            }
            
            // Assemble of the vector is required
            if ( rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_RHS_VECTOR ) ||
                    rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
            {
                // Calculate the local contribution
                const array_1d<double, MatrixSize> RHS_contact_pair = this->CalculateLocalRHS( rThisMortarConditionMatrices, PairIndex, PenaltyFactor, ScaleFactor, ActiveInactive);
                
                // Contribution to previous step contact force and residuals vector
                this->CalculateAndAddRHS( rLocalSystem, RHS_contact_pair, PairIndex );
                
//                 // Debug
//     //                 KRATOS_WATCH(RHS_contact_pair);
//                 LOG_VECTOR_PRETTY( RHS_contact_pair );
            }
        }
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::InitializeGeneralVariables(
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

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateAeAndDeltaAe(
    DerivativeDataType& rDerivativeData,
    GeneralVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    bool IsIntegrated = false; 
    
    // We initilize the Ae components
    AeData rAeData;
    rAeData.Initialize();
    
    rDerivativeData.InitializeDeltaAeComponents();
    
    // We call the exact integration utility
    ExactMortarIntegrationUtility<TDim, TNumNodes> IntUtil = ExactMortarIntegrationUtility<TDim, TNumNodes>(mIntegrationOrder);
    
    for (unsigned int PairIndex = 0; PairIndex < mPairSize; ++PairIndex)
    {   
        // The normal of the master condition
        const array_1d<double, 3> MasterNormal = mThisMasterElements[PairIndex]->GetValue(NORMAL);

        // Reading integration points
        ConditionArrayListType ConditionsPointsSlave;
        const bool IsInside = IntUtil.GetExactIntegration(this->GetGeometry(), this->GetValue(NORMAL), mThisMasterElements[PairIndex]->GetGeometry(), mThisMasterElements[PairIndex]->GetValue(NORMAL), ConditionsPointsSlave);
        
        IntegrationMethod ThisIntegrationMethod = GetIntegrationMethod();
        
        if (IsInside == true)
        {
            IsIntegrated = true;
            
            // Initialize general variables for the current master element
            this->InitializeGeneralVariables( rVariables, rCurrentProcessInfo, PairIndex );
            
            // Update slave element info
            rDerivativeData.UpdateMasterPair(mThisMasterElements[PairIndex]);
            
            for (unsigned int iGeom = 0; iGeom < ConditionsPointsSlave.size(); iGeom++)
            {
                std::vector<PointType::Pointer> PointsArray (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
                for (unsigned int iNode = 0; iNode < TDim; iNode++)
                {
                    PointType GlobalPoint;
                    GetGeometry().GlobalCoordinates(GlobalPoint, ConditionsPointsSlave[iGeom][iNode]);
                    PointsArray[iNode] = boost::make_shared<PointType>(GlobalPoint);
                }
                
                DecompositionType DecompGeom( PointsArray );
                
                const bool BadShape = (TDim == 2) ? false : ContactUtilities::HeronCheck(DecompGeom);
                
                if (BadShape == false)
                {
                    const GeometryType::IntegrationPointsArrayType& IntegrationPointsSlave = DecompGeom.IntegrationPoints( ThisIntegrationMethod );
                    
                    // Integrating the mortar operators
                    for ( unsigned int PointNumber = 0; PointNumber < IntegrationPointsSlave.size(); PointNumber++ )
                    {
                        const PointType LocalPointDecomp = IntegrationPointsSlave[PointNumber].Coordinates();
                        PointType LocalPointParent;
                        PointType GPGlobal;
                        DecompGeom.GlobalCoordinates(GPGlobal, LocalPointDecomp);
                        GetGeometry().PointLocalCoordinates(LocalPointParent, GPGlobal);
                        
                        // Calculate the kinematic variables
                        this->CalculateKinematics( rVariables, rDerivativeData, MasterNormal, LocalPointDecomp, LocalPointParent, DecompGeom);
                        
                        // Update the derivative of DetJ
                        this->CalculateDeltaDetJSlave(rVariables, rDerivativeData); 
                        
                        // Integrate
                        const double IntegrationWeight = IntegrationPointsSlave[PointNumber].Weight();
                        this->CalculateDeltaAeComponents(rVariables, rDerivativeData, rAeData, IntegrationWeight);
                    }
                }
            }
        }
    }
    
    // We can consider the pair if at least one integration point is considered
    if (IsIntegrated == true)
    {
        this->CalculateDeltaAe(rDerivativeData, rAeData);
    }
}

/*********************************COMPUTE KINEMATICS*********************************/
/************************************************************************************/
 
template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateKinematics( 
    GeneralVariables& rVariables,
    const DerivativeDataType rDerivativeData,
    const array_1d<double, 3> MasterNormal,
    const PointType& LocalPointDecomp,
    const PointType& LocalPointParent,
    GeometryPointType& GeometryDecomp
    )
{   
    /* CALCULATE JACOBIAN AND JACOBIAN DETERMINANT */

    rVariables.jSlave = GeometryDecomp.Jacobian( rVariables.jSlave, LocalPointDecomp.Coordinates() );
    rVariables.DetJSlave = GeometryDecomp.DeterminantOfJacobian( LocalPointDecomp );
    
    /// SLAVE CONDITION ///
    
    // SHAPE FUNCTIONS 
    GetGeometry().ShapeFunctionsValues( rVariables.NSlave, LocalPointParent.Coordinates() );
    rVariables.PhiLagrangeMultipliers = prod(rDerivativeData.Ae, rVariables.NSlave);
//     rVariables.PhiLagrangeMultipliers = rVariables.NSlave; // TODO: This could be needed in the future to be different than the standart shape functions 
    
    // SHAPE FUNCTION DERIVATIVES
    GetGeometry().ShapeFunctionsLocalGradients( rVariables.DNDeSlave, LocalPointParent );
    
    // MASTER CONDITION
    this->MasterShapeFunctionValue( rVariables, MasterNormal, LocalPointParent);
}
 
/***********************************************************************************/
/*************** METHODS TO CALCULATE THE CONTACT CONDITION MATRICES ***************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::MasterShapeFunctionValue(
    GeneralVariables& rVariables,
    const array_1d<double, 3> MasterNormal,
    const PointType& LocalPoint 
    )
{    
    GeometryType& MasterSegment = rVariables.GetMasterElement( );

    PointType ProjectedGPGlobal;
    const array_1d<double,3> GPNormal = ContactUtilities::GaussPointNormal(rVariables.NSlave, GetGeometry());
    
    GeometryType::CoordinatesArrayType SlaveGPGlobal;
    this->GetGeometry( ).GlobalCoordinates( SlaveGPGlobal, LocalPoint );
    ContactUtilities::FastProjectDirection( MasterSegment, SlaveGPGlobal, ProjectedGPGlobal, MasterNormal, -GPNormal ); // The opposite direction
    
    GeometryType::CoordinatesArrayType ProjectedGPLocal;
    
    MasterSegment.PointLocalCoordinates( ProjectedGPLocal, ProjectedGPGlobal.Coordinates( ) ) ;

    // SHAPE FUNCTIONS 
    MasterSegment.ShapeFunctionsValues(         rVariables.NMaster,     ProjectedGPLocal );         
    MasterSegment.ShapeFunctionsLocalGradients( rVariables.DNDeMaster, ProjectedGPLocal );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateMortarOperators(
    MortarConditionMatrices& rThisMortarConditionMatrices,
    GeneralVariables& rVariables,
    DerivativeDataType& rDerivativeData,
    const double& rIntegrationWeight
    )
{
    /* DEFINITIONS */
    const double DetJSlave = rVariables.DetJSlave; 
    const VectorType Phi = rVariables.PhiLagrangeMultipliers;
    const VectorType N1  = rVariables.NSlave;
    const VectorType N2  = rVariables.NMaster;
    
    // Derivatives
    constexpr unsigned int Size1 =     (TNumNodes * TDim);
    constexpr unsigned int Size2 = 2 * (TNumNodes * TDim);

    const array_1d<double, Size1> DeltaJSlave  = rDerivativeData.DeltaJSlave;
    const array_1d<array_1d<double, TNumNodes >, Size1> DeltaPhi = rDerivativeData.DeltaPhi;
    const array_1d<array_1d<double, TNumNodes >, Size2> DeltaN1  = rDerivativeData.DeltaN1;
    const array_1d<array_1d<double, TNumNodes >, Size2> DeltaN2  = rDerivativeData.DeltaN2;
    
    // Mortar condition matrices - DOperator and MOperator
    bounded_matrix<double, TNumNodes, TNumNodes>& DOperator = rThisMortarConditionMatrices.DOperator;
    bounded_matrix<double, TNumNodes, TNumNodes>& MOperator = rThisMortarConditionMatrices.MOperator;
    
    // D and M directional derivatives
    array_1d<bounded_matrix<double, TNumNodes, TNumNodes>, Size2>& DeltaDOperator = rThisMortarConditionMatrices.DeltaDOperator;
    array_1d<bounded_matrix<double, TNumNodes, TNumNodes>, Size2>& DeltaMOperator = rThisMortarConditionMatrices.DeltaMOperator;
    
    for (unsigned int iSlave = 0; iSlave < TNumNodes; iSlave++)
    {
        const double phi = Phi[iSlave];
        
        for (unsigned int jSlave = 0; jSlave < TNumNodes; jSlave++)
        {
            const double n1  = N1[jSlave];
            const double n2  = N2[jSlave];
            
            DOperator(iSlave, jSlave) += DetJSlave * rIntegrationWeight * phi * n1;
            MOperator(iSlave, jSlave) += DetJSlave * rIntegrationWeight * phi * n2;
            
            for (unsigned int i = 0; i < TDim * TNumNodes; i++)
            {
                DeltaDOperator[i](iSlave, jSlave) += DeltaJSlave[i] * rIntegrationWeight * phi* n1        
                                                   + DetJSlave * rIntegrationWeight * DeltaPhi[i][iSlave] * n1
                                                   + DetJSlave * rIntegrationWeight * phi* DeltaN1[i][jSlave];
                                                                            
                DeltaMOperator[i](iSlave, jSlave) += DeltaJSlave[i] * rIntegrationWeight * phi* n2        
                                                   + DetJSlave * rIntegrationWeight * DeltaPhi[i][iSlave] * n2
                                                   + DetJSlave * rIntegrationWeight * phi* DeltaN2[i][jSlave];
            }
            for (unsigned int i = TDim * TNumNodes; i < 2 * TDim * TNumNodes; i++)
            {
                DeltaDOperator[i](iSlave, jSlave) += DetJSlave * rIntegrationWeight * phi * DeltaN1[i][jSlave];
                                                                            
                DeltaMOperator[i](iSlave, jSlave) += DetJSlave * rIntegrationWeight * phi * DeltaN2[i][jSlave];
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateMortarOperators(
    MortarConditionMatrices& rThisMortarConditionMatrices,
    GeneralVariables& rVariables,
    const double& rIntegrationWeight
    )
{
    /* DEFINITIONS */
    const double DetJSlave = rVariables.DetJSlave; 
    const VectorType Phi = rVariables.PhiLagrangeMultipliers;
    const VectorType N1  = rVariables.NSlave;
    const VectorType N2  = rVariables.NMaster;

    // Mortar condition matrices - DOperator and MOperator
    bounded_matrix<double, TNumNodes, TNumNodes>& DOperator = rThisMortarConditionMatrices.DOperator;
    bounded_matrix<double, TNumNodes, TNumNodes>& MOperator = rThisMortarConditionMatrices.MOperator;
    
    for (unsigned int iSlave = 0; iSlave < TNumNodes; iSlave++)
    {
        const double phi = Phi[iSlave];
        
        for (unsigned int jSlave = 0; jSlave < TNumNodes; jSlave++)
        {
            DOperator(iSlave, jSlave) += DetJSlave * rIntegrationWeight * phi * N1[jSlave];
            MOperator(iSlave, jSlave) += DetJSlave * rIntegrationWeight * phi * N2[jSlave];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateDeltaAeComponents(
    GeneralVariables& rVariables,
    DerivativeDataType& rDerivativeData,
    AeData& rAeData,
    const double& rIntegrationWeight
    )
{
    /* DEFINITIONS */
    const VectorType N1 = rVariables.NSlave;
    const double detJ   = rVariables.DetJSlave; 
     
    rAeData.De += rIntegrationWeight * this->ComputeDe( N1, detJ);
    rAeData.Me += rIntegrationWeight * this->ComputeMe( N1, detJ);
    
    for (unsigned int i = 0; i < TDim * TNumNodes; i++)
    {
        const double DeltaDetJ = rDerivativeData.DeltaJSlave[i];
        
        const bounded_matrix<double, TNumNodes, TNumNodes> DeltaMe = DeltaDetJ * outer_prod(N1, N1);
        bounded_matrix<double, TNumNodes, TNumNodes> DeltaDe;
        
        for (unsigned int iSlave = 0; iSlave < TNumNodes; iSlave++)
        {
            for (unsigned int jSlave = iSlave; jSlave < TNumNodes; jSlave++)
            {
                if (iSlave == jSlave)
                {
                    DeltaDe(iSlave, iSlave) = DeltaDetJ * N1[iSlave];
                }
                else
                {
                    DeltaDe(iSlave, jSlave) = 0.0;
                    DeltaDe(jSlave, iSlave) = 0.0;
                }
            }
        }
        
        rAeData.DeltaDe[i] += rIntegrationWeight * DeltaDe;
        rAeData.DeltaMe[i] += rIntegrationWeight * DeltaMe;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
bounded_matrix<double, TNumNodes, TNumNodes> AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::ComputeDe(
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
            else
            {
                De(i,j) = 0.0;
            }
        }
    }
    
    return De;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
bounded_matrix<double, TNumNodes, TNumNodes> AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::ComputeMe(
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

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateDeltaAe(
    DerivativeDataType& rDerivativeData,
    AeData& rAeData
    )
{        
    double auxdet = 0.0;
    const double Tolerance = std::numeric_limits<double>::epsilon();
    
    // We compute the norm
    const double NormMe = norm_frobenius(rAeData.Me);
    
    // Now we normalize the matrix
    const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes> NormalizedMe = rAeData.Me/NormMe;
    
    // We compute the normalized inverse
    const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes> NormalizedInvMe = MathUtils<double>::InvertMatrix<TNumNodes>(NormalizedMe, auxdet, Tolerance); 
    
    // Now we compute the inverse
    const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes> InvMe = NormalizedInvMe/NormMe;
    
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

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateAndAddLHS(
    LocalSystemComponents& rLocalSystem,
    const bounded_matrix<double, MatrixSize, MatrixSize>& LHS_contact_pair, 
    const unsigned int rPairIndex
    )
{
    /* DEFINITIONS */

    if ( rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
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
                KRATOS_ERROR <<  " CONDITION can not supply the required local system variable: " << rLeftHandSideVariables[i] << std::endl;
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

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::AssembleContactPairLHSToConditionSystem(
    const bounded_matrix<double, MatrixSize, MatrixSize>& rPairLHS,
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

template<>
bounded_matrix<double, 10, 10> AugmentedLagrangianMethodMortarContactCondition<2,2, false>::CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const unsigned int& rMasterElementIndex,
        const double& rPenaltyFactor,
        const double& rScaleFactor,
        const unsigned int& rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
    
    return ZeroMatrix(10, 10);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
bounded_matrix<double, 21, 21> AugmentedLagrangianMethodMortarContactCondition<3,3, false>::CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const unsigned int& rMasterElementIndex,
        const double& rPenaltyFactor,
        const double& rScaleFactor,
        const unsigned int& rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
    
    return ZeroMatrix(21, 21);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
bounded_matrix<double, 28, 28> AugmentedLagrangianMethodMortarContactCondition<3,4, false>::CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const unsigned int& rMasterElementIndex,
        const double& rPenaltyFactor,
        const double& rScaleFactor,
        const unsigned int& rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
    
    return ZeroMatrix(28, 28);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
bounded_matrix<double, 12, 12> AugmentedLagrangianMethodMortarContactCondition<2,2, true>::CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const unsigned int& rMasterElementIndex,
        const double& rPenaltyFactor,
        const double& rScaleFactor,
        const unsigned int& rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
    
    return ZeroMatrix(12, 12);
}

// /***********************************************************************************/
// /***********************************************************************************/
// 
// template<>
// bounded_matrix<double, 27, 27> AugmentedLagrangianMethodMortarContactCondition<3,3, true>::CalculateLocalLHS(
//         const MortarConditionMatrices& rMortarConditionMatrices,
//         const unsigned int& rMasterElementIndex,
//         const unsigned int& rActiveInactive
//         )
// {
//     KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
//     
//     return ZeroMatrix(27, 27);
// }
// 
// /***********************************************************************************/
// /***********************************************************************************/
// 
// template<>
// bounded_matrix<double, 36, 36> AugmentedLagrangianMethodMortarContactCondition<3,4, true>::CalculateLocalLHS(
//         const MortarConditionMatrices& rMortarConditionMatrices,
//         const unsigned int& rMasterElementIndex,
//         const unsigned int& rActiveInactive
//         )
// {
//     KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
//     
//     return ZeroMatrix(36, 36);
// }

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateAndAddRHS(
    LocalSystemComponents& rLocalSystem,
    const array_1d<double, MatrixSize>& RHS_contact_pair,
    const unsigned int rPairIndex
    )
{   
    if ( rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
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
                KRATOS_ERROR << " CONDITION can not supply the required local system variable: " << rRightHandSideVariables[i] << std::endl;
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

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::AssembleContactPairRHSToConditionSystem(
    const array_1d<double, MatrixSize>& rPairRHS,
    VectorType& rConditionRHS,
    const unsigned int rPairIndex
    )
{
    // Find location of the pair's master DOFs in ConditionRHS
    const unsigned int index_begin = rPairIndex * MatrixSize;
    const unsigned int index_end  = index_begin + MatrixSize;
    
    subrange( rConditionRHS, index_begin, index_end) += rPairRHS;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double,10> AugmentedLagrangianMethodMortarContactCondition<2,2, false>::CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const unsigned int& rMasterElementIndex,
        const double& rPenaltyFactor,
        const double& rScaleFactor,
        const unsigned int& rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
    
    return ZeroVector(10);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double,21> AugmentedLagrangianMethodMortarContactCondition<3,3, false>::CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const unsigned int& rMasterElementIndex,
        const double& rPenaltyFactor,
        const double& rScaleFactor,
        const unsigned int& rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
    
    return ZeroVector(21);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double,28> AugmentedLagrangianMethodMortarContactCondition<3,4, false>::CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const unsigned int& rMasterElementIndex,
        const double& rPenaltyFactor,
        const double& rScaleFactor,
        const unsigned int& rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
    
    return ZeroVector(28);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double,12> AugmentedLagrangianMethodMortarContactCondition<2,2, true>::CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const unsigned int& rMasterElementIndex,
        const double& rPenaltyFactor,
        const double& rScaleFactor,
        const unsigned int& rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
    
    return ZeroVector(12);
}

// /***********************************************************************************/
// /***********************************************************************************/
// 
// template<>
// array_1d<double,27> AugmentedLagrangianMethodMortarContactCondition<3,3, true>::CalculateLocalRHS(
//         const MortarConditionMatrices& rMortarConditionMatrices,
//         const unsigned int& rMasterElementIndex,
//         const unsigned int& rActiveInactive
//         )
// {
//     KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
//     
//     return ZeroVector(27);
// }
// 
// /***********************************************************************************/
// /***********************************************************************************/
// 
// template<>
// array_1d<double,36> AugmentedLagrangianMethodMortarContactCondition<3,4, true>::CalculateLocalRHS(
//         const MortarConditionMatrices& rMortarConditionMatrices,
//         const unsigned int& rMasterElementIndex,
//         const unsigned int& rActiveInactive
//         )
// {
//     KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
//     
//     return ZeroVector(36);
// }

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateDeltaDetJSlave(
   GeneralVariables& rVariables,
   DerivativeDataType& rDerivativeData
   ) // TODO: Do an explicit specialization!!!!
{
    if (TDim == 2)
    {
        // Fill up the elements corresponding to the slave DOFs - the rest remains zero
        for ( unsigned int iSlave = 0, i = 0; iSlave < TNumNodes; ++iSlave, i += TDim )
        {
            rDerivativeData.DeltaJSlave[i    ] = rVariables.jSlave( 0, 0 ) * rVariables.DNDeSlave( iSlave, 0) / rVariables.DetJSlave;
            rDerivativeData.DeltaJSlave[i + 1] = rVariables.jSlave( 1, 0 ) * rVariables.DNDeSlave( iSlave, 0) / rVariables.DetJSlave;
        }
    }
    else
    {
        const array_1d<double,TNumNodes>& DNDxi  = column( rVariables.DNDeSlave, 0 );
        const array_1d<double,TNumNodes>& DNDeta = column( rVariables.DNDeSlave, 1 );
        
        const array_1d<double,TDim>& Jxi  = column( rVariables.jSlave, 0 );
        const array_1d<double,TDim>& Jeta = column( rVariables.jSlave, 1 );
        
        const array_1d<double,TDim>& Normal = prod(trans(rDerivativeData.NormalMaster), rVariables.NSlave); // FIXME: Check this!!!!
        
        bounded_matrix<double, TDim, TDim> DeltaJxixJeta;
        
        for ( unsigned int iSlave = 0, i = 0; iSlave < TNumNodes; ++iSlave, i += TDim )
        {
            DeltaJxixJeta(0,0) = 0.0;
            DeltaJxixJeta(0,1) =  Jeta(2) * DNDxi(iSlave) - Jxi(2) * DNDeta(iSlave); 
            DeltaJxixJeta(0,2) = -Jeta(1) * DNDxi(iSlave) + Jxi(1) * DNDeta(iSlave); 
            DeltaJxixJeta(1,0) = -Jeta(2) * DNDxi(iSlave) + Jxi(2) * DNDeta(iSlave); 
            DeltaJxixJeta(1,1) = 0.0;
            DeltaJxixJeta(1,2) =  Jeta(0) * DNDxi(iSlave) - Jxi(0) * DNDeta(iSlave);
            DeltaJxixJeta(2,0) =  Jeta(1) * DNDxi(iSlave) - Jxi(1) * DNDeta(iSlave); 
            DeltaJxixJeta(2,1) = -Jeta(0) * DNDxi(iSlave) + Jxi(0) * DNDeta(iSlave); 
            DeltaJxixJeta(2,2) = 0.0;
            
            rDerivativeData.DeltaJSlave[i    ] = inner_prod( Normal, column( DeltaJxixJeta, 0 ) );
            rDerivativeData.DeltaJSlave[i + 1] = inner_prod( Normal, column( DeltaJxixJeta, 1 ) );
            rDerivativeData.DeltaJSlave[i + 2] = inner_prod( Normal, column( DeltaJxixJeta, 2 ) );
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional> // NOTE: Formulation taken from Mohamed Khalil work
bounded_matrix<double, TDim, TDim> AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::LocalDeltaNormal( // NOTE: Not the mean, look in the contact utilities 
    const GeometryType& CondGeometry,
    const unsigned int node_index
    )
{
    // Tolerance
    const double Tolerance = 1.0e-14;
        
    bounded_matrix<double, TDim, TDim> DeltaNeAdj;
    bounded_matrix<double, TDim, TDim> Ce;
    
    const bounded_matrix<double, TDim, TDim> I = IdentityMatrix(TDim, TDim);
    
    bounded_matrix<double, TDim, TDim> DeltaNormal = ZeroMatrix(TDim,TDim);
    
    const array_1d<double, 3> Ne = this->GetValue(NORMAL); // Normalized condition normal
    bounded_matrix<double, TDim, TDim> NeoNe = subrange( outer_prod( Ne, Ne ), 0, TDim, 0, TDim );
    
    if (TDim == 2)
    {
        const double NeNorm = this->GetGeometry( ).Length( ); // The norm of a geometry's normal is its characteristic dimension - length for 2D and area for 3D 
                        
        DeltaNeAdj( 0, 0 ) =  0.0;
        DeltaNeAdj( 0, 1 ) = -1.0;
        DeltaNeAdj( 1, 0 ) =  1.0;
        DeltaNeAdj( 1, 1 ) =  0.0;
        
        double DNDej = 0.0;
        if( node_index == 0 )
        {
            DNDej = - 0.5;
        }
        else if( node_index == 1 )
        {
            DNDej =   0.5;
        }
        
        Ce = prod( I - NeoNe, DeltaNeAdj ) / NeNorm; // In 2D, DeltaNeAdj is node-independent => evaluated outside the nodes loop
        
        DeltaNormal = - 2.0 * Ce * DNDej; // NOTE: Check why - 2???!!!, it was the only wayto ensure the same value as the symbolic. You will need to repeat this in 3D            
//         DeltaNormal = Ce * DNDej;             
    }
    else
    {
        const double NeNorm = this->GetGeometry( ).Area( ); // The norm of a geometry's normal is its characteristic dimension - length for 2D and area for 3D 
        
        MatrixType J = ZeroMatrix( 3, 2 ); // Jacobian [ 3D global x 2D local ]
        array_1d<double, 2> DNDej;
        array_1d<double, 3> LocalCoordsj;
        
        if( TNumNodes == 3 )    // linear triangle element
        {
            if( node_index == 0 )
            {
                LocalCoordsj[0] = 0.0;
                LocalCoordsj[1] = 0.0;
                DNDej[0] = - 1.0;
                DNDej[1] = - 1.0;
            }
            else if( node_index == 1 )
            {
                LocalCoordsj[0] = 1.0;
                LocalCoordsj[1] = 0.0;
                DNDej[0] = 1.0;
                DNDej[1] = 0.0;
            }
            else // node_index == 2
            {
                LocalCoordsj[0] = 0.0;
                LocalCoordsj[1] = 1.0;
                DNDej[0] = 0.0;
                DNDej[1] = 1.0;
            }
        }
        else if( TNumNodes == 4 )    // linear quad element 
        {
            if( node_index == 0 )
            {
                LocalCoordsj[0] = - 1.0;
                LocalCoordsj[1] = - 1.0;
                DNDej[0] = - 0.5;
                DNDej[1] = - 0.5;
            }
            else if( node_index == 1 )
            {
                LocalCoordsj[0] =   1.0;
                LocalCoordsj[1] = - 1.0;
                DNDej[0] =   0.5;
                DNDej[1] = - 0.5;
            }
            else if( node_index == 2 )
            {
                LocalCoordsj[0] =  1.0;
                LocalCoordsj[1] =  1.0;
                DNDej[0] =  0.5;
                DNDej[1] =  0.5;
            }
            else // node_index == 3
            {
                LocalCoordsj[0] = - 1.0;
                LocalCoordsj[1] =   1.0;
                DNDej[0] = - 0.5;
                DNDej[1] =   0.5;
            }
        }
        
        this->GetGeometry( ).Jacobian( J, LocalCoordsj );
        
        DeltaNeAdj(0,0) = 0.0;
        DeltaNeAdj(0,1) = +J(2,1) * DNDej[0] - J(2,0) * DNDej[1]; 
        DeltaNeAdj(0,2) = -J(1,1) * DNDej[0] + J(1,0) * DNDej[1]; 
        DeltaNeAdj(1,0) = -J(2,1) * DNDej[0] + J(2,0) * DNDej[1]; 
        DeltaNeAdj(1,1) = 0.0;                   
        DeltaNeAdj(1,2) = +J(0,1) * DNDej[0] - J(0,0) * DNDej[1]; 
        DeltaNeAdj(2,0) = +J(1,1) * DNDej[0] - J(1,0) * DNDej[1]; 
        DeltaNeAdj(2,1) = -J(0,1) * DNDej[0] + J(0,0) * DNDej[1]; 
        DeltaNeAdj(2,2) = 0.0;
        
        Ce = prod( I - NeoNe, DeltaNeAdj ) / NeNorm;
        DeltaNormal = Ce;
    }
    
    const double NeNorm = norm_2( Ne );
    const double NeNorm3 = NeNorm * NeNorm * NeNorm;
    
    if ( NeNorm3 > Tolerance )
    {
        const bounded_matrix<double, TDim, TDim> Cj = I / NeNorm - NeoNe / NeNorm3;
        DeltaNormal = prod( Cj, DeltaNormal );
    }
        
    return DeltaNormal; 
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional> // NOTE: Formulation taken from Mohamed Khalil work
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateDeltaNormalSlave(DerivativeDataType& rDerivativeData)
{
    if (TDim == 2) // TODO: Use explicit 
    {
        for ( unsigned int iSlave = 0, i = 0; iSlave < TNumNodes; ++iSlave, i += TDim )
        {
//             bounded_matrix<double, 2, 2> DeltaNormal = GetGeometry()[iSlave].GetValue(DELTA_NORMAL);
            bounded_matrix<double, 2, 2> DeltaNormal = this->LocalDeltaNormal(GetGeometry(), iSlave);
            for (unsigned iDoF = 0; iDoF < TDim; iDoF++) 
            {
                for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
                {
                    row(rDerivativeData.DeltaNormalSlave[iSlave * TDim + iDoF], iNode) = trans(column(DeltaNormal, iDoF)); 
                }
            }
        }
    }
    else
    {
        for ( unsigned int iSlave = 0, i = 0; iSlave < TNumNodes; ++iSlave, i += TDim )
        {
//             bounded_matrix<double, 3, 3> DeltaNormal = GetGeometry()[iSlave]->GetValue(DELTA_NORMAL);
            const bounded_matrix<double, 3, 3> DeltaNormal = this->LocalDeltaNormal(this->GetGeometry(), iSlave);
            
            // Calculate nodal tangents
            
            const MatrixType I = IdentityMatrix( TDim, TDim );
            
            array_1d<double, 2> DNDej;
            if( TNumNodes == 3 )    // linear triangle element // TODO: Use an enum
            {
                if( iSlave == 0 )
                {
                    DNDej[0] = -1.0;
                    DNDej[1] = -1.0;
                }
                else if( iSlave == 1 )
                {
                    DNDej[0] = 1.0;
                    DNDej[1] = 0.0;
                }
                else // iSlave == 2
                {
                    DNDej[0] = 0.0;
                    DNDej[1] = 1.0;
                }
            }
            else if( TNumNodes == 4 )    // linear quad element 
            {
                if( iSlave == 0 )
                {
                    DNDej[0] = -0.5;
                    DNDej[1] = -0.5;
                }
                else if( iSlave == 1 )
                {
                    DNDej[0] =  0.5;
                    DNDej[1] = -0.5;
                }
                else if( iSlave == 2 )
                {
                    DNDej[0] =  0.5;
                    DNDej[1] =  0.5;
                }
                else // iSlave == 3
                {
                    DNDej[0] = -0.5;
                    DNDej[1] =  0.5;
                }
            }
            
            for (unsigned iDoF = 0; iDoF < TDim; iDoF++) 
            {
                for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
                {
                    row(rDerivativeData.DeltaNormalSlave[iSlave * TDim + iDoF], iNode) = trans(column(DeltaNormal, iDoF)); 
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional> // NOTE: Formulation taken from Mohamed Khalil work
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateDeltaNormalMaster(DerivativeDataType& rDerivativeData)
{
    if (TDim == 2)
    {
        for ( unsigned int i_master = 0, i = 0; i_master < TNumNodes; ++i_master, i += TDim )
        {
//             bounded_matrix<double, 2, 2> DeltaNormal = GetGeometry[i_master].GetValue(DELTA_NORMAL);
            bounded_matrix<double, 2, 2> DeltaNormal = this->LocalDeltaNormal(rDerivativeData.MasterGeometry, i_master);
            for (unsigned iDoF = 0; iDoF < TDim; iDoF++) 
            {
                for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
                {
                    row(rDerivativeData.DeltaNormalMaster[i_master * TDim + iDoF], iNode) = trans(column(DeltaNormal, iDoF)); 
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
            for (unsigned iDoF = 0; iDoF < TDim; iDoF++) 
            {
                for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
                {
                    row(rDerivativeData.DeltaNormalMaster[i_master * TDim + iDoF], iNode)  = trans(column(DeltaNormal, iDoF)); 
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateDeltaN(
   GeneralVariables& rVariables,
   DerivativeDataType& rDerivativeData
   )
{
    static const unsigned int Size1 = (TNumNodes * TDim);
    
    // Shape functions
    const VectorType  N1 = rVariables.NSlave;
    const VectorType  N2 = rVariables.NMaster;
    
    // Coordinates
    const bounded_matrix<double, TNumNodes, TDim> u1 = rDerivativeData.u1;
    const bounded_matrix<double, TNumNodes, TDim> X1 = rDerivativeData.X1;
    const bounded_matrix<double, TNumNodes, TDim> u2 = rDerivativeData.u2;
    const bounded_matrix<double, TNumNodes, TDim> X2 = rDerivativeData.X2;
    
    // Normals
    const array_1d<double, TDim >  NormalSlaveg = prod(trans(rDerivativeData.NormalSlave), N1);
    const array_1d<double, TDim >  NormalMasterg = prod(trans(rDerivativeData.NormalMaster), N2);
    
    const array_1d<bounded_matrix<double, TNumNodes, TDim>, Size1> DNormalSlave = rDerivativeData.DeltaNormalSlave;
    const array_1d<bounded_matrix<double, TNumNodes, TDim>, Size1> DNormalMaster = rDerivativeData.DeltaNormalMaster;
    
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
      const double Dist = inner_prod(vector_nodes, NormalMasterg)/(inner_prod(NormalSlaveg, NormalMasterg) + tol);
      
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
      for (unsigned int iSlave = 0; iSlave < TNumNodes; iSlave++)
      {
         for (unsigned int i_dim = 0; i_dim < TDim; i_dim++)
         {
            const unsigned int iDoF = iSlave * TDim + i_dim;
            
            const array_1d<double, TNumNodes > DNormalSlaveg = prod(trans(DNormalSlave[iDoF]), N1);
            array_1d<double, TNumNodes > aux_vector = ZeroVector(TNumNodes);
            aux_vector[i_dim] = 1.0;
            const array_1d<double, TNumNodes > Deltavector_nodes = - N1[iSlave] * aux_vector;
            const double DeltaDist = ((inner_prod(Deltavector_nodes, NormalMasterg))* inner_prod(NormalSlaveg, NormalMasterg) - inner_prod(vector_nodes, NormalMasterg) * (inner_prod(DNormalSlaveg, NormalMasterg)))/std::pow(inner_prod(NormalSlaveg, NormalMasterg) + tol, 2);
            const array_1d<double, TNumNodes > Deltax1g = N1[iSlave] * aux_vector;
            const array_1d<double, TNumNodes > Deltax2g = Deltax1g + DeltaDist * NormalSlaveg + Dist * DNormalSlaveg;
            
            if (TDim == 2)
            {
               if (TNumNodes == 2)
               {
                  rDerivativeData.DeltaN2[iDoF][0] =  (Deltax2g[0] + Deltax2g[1])/div2;
                  rDerivativeData.DeltaN2[iDoF][1] =  - rDerivativeData.DeltaN2[iDoF][0];
               }
            }
            else
            {
               if (TNumNodes == 3)
               {
                  rDerivativeData.DeltaN2[iDoF][1] = - (mult1 * (Deltax2g[0] + Deltax2g[1]) + mult2 * (Deltax2g[0] + Deltax2g[2]))/div1;
                  rDerivativeData.DeltaN2[iDoF][2] =   (mult3 *  Deltax2g[0] + mult4 * Deltax2g[1] + mult5 * Deltax2g[2])/div2;
                  rDerivativeData.DeltaN2[iDoF][0] =  - rDerivativeData.DeltaN2[iDoF][1] - rDerivativeData.DeltaN2[iDoF][2];
               }
            }
         }
      }
      
      // Derivatives master
      for (unsigned int i_master = 0; i_master < TNumNodes; i_master++)
      {
         for (unsigned int i_dim = 0; i_dim < TDim; i_dim++)
         {
            const unsigned int iDoF = (TNumNodes + i_master) * TDim + i_dim;
            
            const array_1d<double, TNumNodes > DNormalMasterg = prod(trans(DNormalMaster[iDoF - TNumNodes * TDim]), N2);
            array_1d<double, TNumNodes > aux_vector = ZeroVector(TNumNodes);
//             if (i_master == 0) // NOTE: This is the way I considered in the symbolic
//             {
                aux_vector[i_dim] = 1.0;
//             }
            const array_1d<double, TNumNodes > Deltavector_nodes = N2[i_master] * aux_vector;
            const double DeltaDist = ((inner_prod(Deltavector_nodes, NormalMasterg) + inner_prod(vector_nodes, DNormalMasterg))* inner_prod(NormalSlaveg, NormalMasterg) - inner_prod(vector_nodes, NormalMasterg) * (inner_prod(NormalSlaveg, DNormalMasterg)))/std::pow(inner_prod(NormalSlaveg, NormalMasterg) + tol, 2);;
            const array_1d<double, TNumNodes > x2g = prod(trans(X2 + u2), N2);
            const array_1d<double, TNumNodes > Deltax2g = DeltaDist * NormalSlaveg;
            
            if (TDim == 2)
            {
               if (TNumNodes == 2)
               {
                   if (i_master == 0)
                   {
                       rDerivativeData.DeltaN2[iDoF][0] =  ( ((u2(1,0) + X2(1,0) + u2(1,1) + X2(1,1)) - x2g[0] - x2g[1]) + div2 * (Deltax2g[0] + Deltax2g[1]))/std::pow(div2, 2);
                   }
                   else 
                   {
                       rDerivativeData.DeltaN2[iDoF][0] =  ((-(u2(0,0) + X2(0,0) + u2(0,1) + X2(0,1)) + x2g[0] + x2g[1]) + div2 * (Deltax2g[0] + Deltax2g[1]))/std::pow(div2, 2);
                   }
                  
                  rDerivativeData.DeltaN2[iDoF][1] =  - rDerivativeData.DeltaN2[iDoF][0];
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
                   
                   if ((iDoF - TDim * TNumNodes) == 0)
                   {
                       rDerivativeData.DeltaN2[iDoF][1] =  ((u2(1, 1) - u2(1, 2) - u2(2, 1) + u2(2, 2) + X2(1, 1) - 
                                                           X2(1, 2) - X2(2, 1) + X2(2, 2)) * (multmaster0*multmaster2 + 
                                                           multmaster1*multmaster3) - (div1) * (u2(0, 1) - u2(0, 2) + 
                                                           X2(0, 1) - X2(0, 2) - x2g[1] + x2g[2] + 
                                                           mult1 * (-1 + Deltax2g[0] + Deltax2g[1]) + 
                                                           mult2 * (-1 + Deltax2g[0] + Deltax2g[2])) )/std::pow(div1, 2);
                                                           
                       rDerivativeData.DeltaN2[iDoF][2] =  ( -(u2(1, 1) - u2(1, 2) - u2(2, 1) + u2(2, 2) + X2(1, 1) - 
                                                            X2(1, 2) - X2(2, 1) + X2(2, 2)) * (coefmaster2 + (u2(0, 2) + X2(0, 2)) * 
                                                            multmaster5 + coefmaster1 - coefmaster5 - 
                                                            multmaster4 * multmaster6) + (div2) * (u2(1, 1) - u2(1, 2) + 
                                                            X2(1, 1) - X2(1, 2) - x2g[1] + x2g[2] + mult3 * Deltax2g[0] + 
                                                            mult4 * Deltax2g[1] + mult5 * Deltax2g[2]))/std::pow(div2, 2);
                   }
                   else if ((iDoF - TDim * TNumNodes) == 1)
                   {
                       rDerivativeData.DeltaN2[iDoF][1] =  ((-u2(1, 0) - u2(1, 2) + u2(2, 0) + u2(2, 2) - X2(1, 0) - 
                                                            X2(1, 2) + X2(2, 0) + X2(2, 2)) * (multmaster0*multmaster2 + 
                                                            multmaster1*multmaster3) - (div1) * (-u2(0, 0) - u2(0, 2) - 
                                                            X2(0, 0) - X2(0, 2) + x2g[0] + x2g[2] + 
                                                            mult1 * (-1 + Deltax2g[0] + Deltax2g[1]) + 
                                                            mult2 * (Deltax2g[0] + Deltax2g[2])) )/std::pow(div1, 2);
                                                            
                       rDerivativeData.DeltaN2[iDoF][2] =  (-(-u2(1, 0) - u2(1, 2) + u2(2, 0) + u2(2, 2) - X2(1, 0) - 
                                                             X2(1, 2) + X2(2, 0) + X2(2, 2)) * (coefmaster2 + (u2(0, 2) + X2(0, 2)) * 
                                                             multmaster5 + coefmaster1 - coefmaster5 - 
                                                             multmaster4 * multmaster6) + (div2) * (-u2(1, 0) - u2(1, 2) - 
                                                             X2(1, 0) - X2(1, 2) + x2g[0] + x2g[2] + mult3 * Deltax2g[0] + 
                                                             mult4 * Deltax2g[1] + mult5 * Deltax2g[2]) )/std::pow(div2, 2);
                   }
                   else if ((iDoF - TDim * TNumNodes) == 2)
                   {
                       rDerivativeData.DeltaN2[iDoF][1] =  ( (u2(1, 0) + u2(1, 1) - u2(2, 0) - u2(2, 1) + X2(1, 0) + 
                                                            X2(1, 1) - X2(2, 0) - X2(2, 1)) * (multmaster0*multmaster2 + 
                                                            multmaster1*multmaster3) - (div1) * (u2(0, 0) + u2(0, 1) + 
                                                            X2(0, 0) + X2(0, 1) - x2g[0] - x2g[1] + 
                                                            mult1 * (Deltax2g[0] + Deltax2g[1]) + 
                                                            mult2 * (-1 + Deltax2g[0] + Deltax2g[2])))/std::pow(div1, 2);
                        
                       rDerivativeData.DeltaN2[iDoF][2] =  (-(u2(1, 0) + u2(1, 1) - u2(2, 0) - u2(2, 1) + X2(1, 0) + 
                                                            X2(1, 1) - X2(2, 0) - X2(2, 1)) * (coefmaster2 + (u2(0, 2) + X2(0, 2)) * 
                                                            multmaster5 + coefmaster1 - coefmaster5 - 
                                                            multmaster4 * multmaster6) + (div2) * (u2(1, 0) + u2(1, 1) + 
                                                            X2(1, 0) + X2(1, 1) - x2g[0] - x2g[1] + mult3 * Deltax2g[0] + 
                                                            mult4 * Deltax2g[1] + mult5 * Deltax2g[2]) )/std::pow(div2, 2);
                   }
                   else if ((iDoF - TDim * TNumNodes) == 3)
                   {
                       rDerivativeData.DeltaN2[iDoF][1] =  ( (-u2(0, 1) + u2(0, 2) + u2(2, 1) - u2(2, 2) - X2(0, 1) + 
                                                            X2(0, 2) + X2(2, 1) - X2(2, 2)) * (multmaster0*multmaster2 + 
                                                            multmaster1 * multmaster3) - (div1) * (+mult1 * (Deltax2g[0] + Deltax2g[1]) + 
                                                            mult2 * (Deltax2g[0] + Deltax2g[2])))/std::pow(div1, 2);
                                                            
                       rDerivativeData.DeltaN2[iDoF][2] =  ( -(-u2(0, 1) + u2(0, 2) + u2(2, 1) - u2(2, 2) - X2(0, 1) + 
                                                              X2(0, 2) + X2(2, 1) - X2(2, 2)) * (coefmaster2 + (u2(0, 2) + X2(0, 2)) * 
                                                              multmaster5 + coefmaster1 - coefmaster5 - 
                                                              multmaster4 * multmaster6) + (div2) * (-u2(0, 1) + u2(0, 2) - 
                                                              X2(0, 1) + X2(0, 2) + x2g[1] - x2g[2] + mult3 * Deltax2g[0] + 
                                                              mult4 * Deltax2g[1] + mult5 * Deltax2g[2]))/std::pow(div2, 2);
                   }
                   else if ((iDoF - TDim * TNumNodes) == 4)
                   {
                       rDerivativeData.DeltaN2[iDoF][1] =  ( multmaster0*(multmaster0*multmaster2 + 
                                                           multmaster1 * multmaster3) - (div1) * (+mult1 * (Deltax2g[0] + Deltax2g[1]) + 
                                                           mult2 * (Deltax2g[0] + Deltax2g[2])))/std::pow(div1, 2);
                                                           
                       rDerivativeData.DeltaN2[iDoF][2] =  ( +mult1 * (coefmaster2 + (u2(0, 2) + X2(0, 2)) * multmaster5 + 
                                                            coefmaster1 - coefmaster5 - multmaster4 * multmaster6) + (div2) * (coefmaster3 + 
                                                            mult3 * Deltax2g[0] + mult4 * Deltax2g[1] + mult5 * Deltax2g[2]))/std::pow(div2, 2);
                   }
                   else if ((iDoF - TDim * TNumNodes) == 5)
                   {
                       rDerivativeData.DeltaN2[iDoF][1] =  ( multmaster1 *(multmaster0*multmaster2 + 
                                                           multmaster1 * multmaster3) - (div1) * (+mult1 * (Deltax2g[0] + Deltax2g[1]) + 
                                                           mult2 * (Deltax2g[0] + Deltax2g[2])))/std::pow(div1, 2);
                                                           
                       rDerivativeData.DeltaN2[iDoF][2] =  ( -multmaster1 *  (coefmaster2 + (u2(0, 2) + X2(0, 2)) * multmaster5 +
                                                            coefmaster1 - coefmaster5 -  multmaster4 * multmaster6) + (div2) * (coefmaster4 + 
                                                            mult3 * Deltax2g[0] + mult4 * Deltax2g[1] + mult5 * Deltax2g[2]))/std::pow(div2, 2);
                   }
                   else if ((iDoF - TDim * TNumNodes) == 6)
                   {
                       rDerivativeData.DeltaN2[iDoF][1] =  (mult3 *(multmaster0*multmaster2 + 
                                                          multmaster1*multmaster3) - (div1) * (-u2(0, 1) + u2(0, 2) - 
                                                          X2(0, 1) + X2(0, 2) + x2g[1] - x2g[2] + 
                                                          mult1 * (Deltax2g[0] + Deltax2g[1]) + 
                                                          mult2 * (Deltax2g[0] + Deltax2g[2])) )/std::pow(div1, 2);
                                                          
                       rDerivativeData.DeltaN2[iDoF][2] =  ( -mult3 * (coefmaster2 + (u2(0, 2) + X2(0, 2)) * multmaster5 + 
                                                            coefmaster1 - coefmaster5 - multmaster4 * multmaster6) + (div2) * (mult3 * Deltax2g[0] + 
                                                            mult4 * Deltax2g[1] + mult5 * Deltax2g[2]))/std::pow(div2, 2);
                   }
                   else if ((iDoF - TDim * TNumNodes) == 7)
                   {
                       rDerivativeData.DeltaN2[iDoF][1] =  (mult4 *(multmaster0*multmaster2 + 
                                                          multmaster1*multmaster3) - (div1) * (coefmaster3 + 
                                                          mult1 * (Deltax2g[0] + Deltax2g[1]) + 
                                                          mult2 * (Deltax2g[0] + Deltax2g[2])))/std::pow(div1, 2);
                                                          
                       rDerivativeData.DeltaN2[iDoF][2] =  ( -mult4 * (coefmaster2 + (u2(0, 2) + X2(0, 2)) * multmaster5 + 
                                                            coefmaster1 - coefmaster5 - multmaster4 * multmaster6) + (div2) * (mult3 * Deltax2g[0] + 
                                                            mult4 * Deltax2g[1] + mult5 * Deltax2g[2]))/std::pow(div2, 2);
                   }
                   else if ((iDoF - TDim * TNumNodes) == 8)
                   {
                       rDerivativeData.DeltaN2[iDoF][1] =  (mult5 *(multmaster0*multmaster2 + 
                                                          multmaster1*multmaster3) - (div1) * (coefmaster4 + 
                                                          mult1 * (Deltax2g[0] + Deltax2g[1]) + 
                                                          mult2 * (Deltax2g[0] + Deltax2g[2])))/std::pow(div1, 2);
                                                          
                       rDerivativeData.DeltaN2[iDoF][2] =  (-mult5 * (coefmaster2 + (u2(0, 2) + X2(0, 2)) * multmaster5 + 
                                                          coefmaster1 - coefmaster5 - multmaster4 * multmaster6) + (div2) * (mult3 * Deltax2g[0] + 
                                                          mult4 * Deltax2g[1] + mult5 * Deltax2g[2]) )/std::pow(div2, 2);
                   }

                  rDerivativeData.DeltaN2[iDoF][0] =  - rDerivativeData.DeltaN2[iDoF][1] - rDerivativeData.DeltaN2[iDoF][2];
               }
            }
         }
      }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateDeltaPhi(
   GeneralVariables& rVariables,
   DerivativeDataType& rDerivativeData
   )
{
    // Shape functions
    const VectorType N1 = rVariables.NSlave;
    
    for (unsigned int iSlave = 0; iSlave < TNumNodes; iSlave++)
    {
        for (unsigned int i_dim = 0; i_dim < TDim; i_dim++)
        {
            const unsigned int iDoF = iSlave * TDim + i_dim;
            
            rDerivativeData.DeltaPhi[iDoF] = prod(rDerivativeData.DeltaAe[iDoF], N1);;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& CurrentProcessInfo 
    )
{
    KRATOS_ERROR << "You are calling to the base class method EquationIdVector, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim, TNumNodes, TFrictional>::GetDofList(
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo 
)
{
    KRATOS_ERROR << "You are calling to the base class method GetDofList, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
}

//******************************* GET DOUBLE VALUE *********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::GetValueOnIntegrationPoints( 
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//******************************* GET ARRAY_1D VALUE *******************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::GetValueOnIntegrationPoints( 
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//******************************* GET VECTOR VALUE *********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::GetValueOnIntegrationPoints( 
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateOnIntegrationPoints( 
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;

    // TODO: Fill this!!!
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateOnIntegrationPoints( 
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector< array_1d<double, 3 > >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
    
    // TODO: Fill this!!!
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateOnIntegrationPoints( 
    const Variable<Vector>& rVariable, 
    std::vector<Vector>& rOutput, 
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
    
    // TODO: Fill this!!!
    
    KRATOS_CATCH( "" );
}

/**************************COMPUTE DELTA POSITION***********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
Matrix AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateDeltaPosition()
{
    KRATOS_TRY;

    Matrix DeltaPosition(TNumNodes, TDim);

    for ( unsigned int i = 0; i < TNumNodes; i++ )
    {
        const array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
        const array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);

        for ( unsigned int j = 0; j < TDim; j++ )
        {
            DeltaPosition(i,j) = (CurrentDisplacement[j] - PreviousDisplacement[j]);
        }
    }

    return DeltaPosition;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

// Frictionless cases
template class AugmentedLagrangianMethodMortarContactCondition<2, 2, false>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 3, false>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 4, false>;

// Frictional cases
template class AugmentedLagrangianMethodMortarContactCondition<2, 2, true>;
// template class AugmentedLagrangianMethodMortarContactCondition<3, 3, true>;
// template class AugmentedLagrangianMethodMortarContactCondition<3, 4, true>;

} // Namespace Kratos
