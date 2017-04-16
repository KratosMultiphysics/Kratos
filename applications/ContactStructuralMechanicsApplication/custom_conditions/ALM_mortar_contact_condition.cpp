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
#include "custom_utilities/contact_utilities.h"
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
    
    unsigned int i_cond = 0;
    for (auto ipair = AllConditionSets->begin(); ipair != AllConditionSets->end(); ++ipair )
    {
        mThisMasterElements[i_cond] = ipair->first;
        i_cond += 1;
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
    
    // Create and initialize condition variables:#pragma omp critical
    GeneralVariables rVariables;
    
    // Create the current contact data
    DerivativeData rDerivativeData;
    
    // Create the mortar operators
    MortarConditionMatrices rThisMortarConditionMatrices;
                     
    this->InitializeDerivativeData(rDerivativeData, rCurrentProcessInfo);
    
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
            
            // Update the contact data
            this->UpdateDerivativeData(rDerivativeData, PairIndex);
            
            // Initialize the mortar operators
            rThisMortarConditionMatrices.Initialize();
            
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
                        this->CalculateMortarOperators(rThisMortarConditionMatrices, rVariables, IntegrationWeight);
                    }
                }
            }
            
//             // Debug
//             std::cout << "--------------------------------------------------" << std::endl;
//             KRATOS_WATCH(this->Id());
//             KRATOS_WATCH(PairIndex);
//             rThisMortarConditionMatrices.print();
            
            // Setting the weighted gap
            // Mortar condition matrices - DOperator and MOperator
            const bounded_matrix<double, TNumNodes, TNumNodes>& DOperator = rThisMortarConditionMatrices.DOperator;
            const bounded_matrix<double, TNumNodes, TNumNodes>& MOperator = rThisMortarConditionMatrices.MOperator;
            
            // Current coordinates 
            const bounded_matrix<double, TNumNodes, TDim> x1 = ContactUtilities::GetCoordinates<TDim,TNumNodes>(this->GetGeometry());
            const bounded_matrix<double, TNumNodes, TDim> x2 = ContactUtilities::GetCoordinates<TDim,TNumNodes>(mThisMasterElements[PairIndex]->GetGeometry());
    
            const bounded_matrix<double, TNumNodes, TDim> Dx1Mx2 = prod(DOperator, x1) - prod(MOperator, x2); 
            
            for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
            {
                const array_1d<double, TDim> normal = subrange(GetGeometry()[iNode].GetValue(NORMAL), 0, TDim);
                const array_1d<double, TDim> aux_array = row(Dx1Mx2, iNode);
                                
                #pragma omp atomic 
                GetGeometry()[iNode].GetValue(WEIGHTED_GAP) += inner_prod(aux_array, - normal); 
            }
            
            if (TFrictional == true) // TODO: Check this!!!
            {
                // Current coordinates 
                const bounded_matrix<double, TNumNodes, TDim> x1Old = ContactUtilities::GetCoordinates<TDim,TNumNodes>(this->GetGeometry(), false, 1);
                const bounded_matrix<double, TNumNodes, TDim> x2Old = ContactUtilities::GetCoordinates<TDim,TNumNodes>(mThisMasterElements[PairIndex]->GetGeometry(), false, 1);
        
                const bounded_matrix<double, TNumNodes, TDim> Dx1OldMx2Old = prod(DOperator, x1Old) - prod(MOperator, x2Old); 
                
                for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
                {
                    // The tangents
                    const array_1d<double, TDim> TangentXi  = subrange(GetGeometry()[iNode].GetValue(TANGENT_XI),  0, TDim);
                    const array_1d<double, TDim> TangentEta = subrange(GetGeometry()[iNode].GetValue(TANGENT_ETA), 0, TDim);
                    
                    const array_1d<double, TDim> aux_array = row(Dx1OldMx2Old, iNode);
                                    
                    #pragma omp atomic 
                    GetGeometry()[iNode].GetValue(WEIGHTED_SLIP) += std::sqrt(std::pow(inner_prod(aux_array, TangentXi), 2) + std::pow(inner_prod(aux_array, TangentEta), 2)); 
                }
            }
        }
    }
    
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
    
    // Create and initialize condition variables:#pragma omp critical
    GeneralVariables rVariables;
    
    // Create the current contact data
    DerivativeData rDerivativeData;
    
    // Create the mortar operators
    MortarConditionMatrices rThisMortarConditionMatrices;
                                                                                  
    this->InitializeDerivativeData(rDerivativeData, rCurrentProcessInfo);
    
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
            
            // Update the contact data
            this->UpdateDerivativeData(rDerivativeData, PairIndex);
            
            // Initialize the mortar operators
            rThisMortarConditionMatrices.Initialize();
            
//             // Compute the normal derivatives of the master
//             this->CalculateDeltaNormalMaster(rDerivativeData);
            
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
                    
                        /* Update the derivatives */
                        // Update the derivative of DetJ
                        this->CalculateDeltaDetJSlave(rVariables, rDerivativeData);
                        // Update the derivatives of the shape functions and the gap
    //                     this->CalculateDeltaN(rVariables, rDerivativeData); // FIXME: This is the old version!!!!
                        // The derivatives of the dual shape function 
                        this->CalculateDeltaPhi(rVariables, rDerivativeData);
                        
                        const double IntegrationWeight = IntegrationPointsSlave[PointNumber].Weight();
                        
                        this->CalculateMortarOperators(rThisMortarConditionMatrices, rVariables, rDerivativeData, IntegrationWeight);       
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
                bounded_matrix<double, MatrixSize, MatrixSize> LHS_contact_pair = this->CalculateLocalLHS( rThisMortarConditionMatrices, PairIndex, ActiveInactive);
                
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
                const array_1d<double, MatrixSize> RHS_contact_pair = this->CalculateLocalRHS( rThisMortarConditionMatrices, PairIndex, ActiveInactive);
                
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
    DerivativeData& rDerivativeData,
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
            
            // Update the contact data
            this->UpdateDerivativeData(rDerivativeData, PairIndex);
            
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

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::InitializeDerivativeData(
    DerivativeData& rDerivativeData,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // Slave element info
    rDerivativeData.Initialize(GetGeometry());
    
    /* NORMALS */
    rDerivativeData.NormalSlave = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(GetGeometry(),  NORMAL); 
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::UpdateDerivativeData(
    DerivativeData& rDerivativeData,
    const unsigned int& rMasterElementIndex
    )
{    
    // Slave element info
    rDerivativeData.UpdateMasterPair(mThisMasterElements[rMasterElementIndex] );
    
    /* NORMALS */
    for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
    {
        array_1d<double,3> normal = GetGeometry()[iNode].GetValue(NORMAL);
    }
}

/*********************************COMPUTE KINEMATICS*********************************/
/************************************************************************************/
 
template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateKinematics( 
    GeneralVariables& rVariables,
    const DerivativeData rDerivativeData,
    const array_1d<double, 3> MasterNormal,
    const PointType& LocalPointDecomp,
    const PointType& LocalPointParent,
    GeometryPointType& GeometryDecomp
    )
{   
    /* CALCULATE JACOBIAN AND JACOBIAN DETERMINANT */

    rVariables.j_Slave = GeometryDecomp.Jacobian( rVariables.j_Slave, LocalPointDecomp.Coordinates() );
    rVariables.DetJSlave = GeometryDecomp.DeterminantOfJacobian( LocalPointDecomp );
    
    /// SLAVE CONDITION ///
    
    // SHAPE FUNCTIONS 
    GetGeometry().ShapeFunctionsValues( rVariables.N_Slave, LocalPointParent.Coordinates() );
    rVariables.Phi_LagrangeMultipliers = prod(rDerivativeData.Ae, rVariables.N_Slave);
//     rVariables.Phi_LagrangeMultipliers = rVariables.N_Slave; // TODO: This could be needed in the future to be different than the standart shape functions 
    
    // SHAPE FUNCTION DERIVATIVES
    GetGeometry().ShapeFunctionsLocalGradients( rVariables.DN_De_Slave, LocalPointParent );
    
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
    const array_1d<double,3> GPNormal = ContactUtilities::GaussPointNormal(rVariables.N_Slave, GetGeometry());
    
    GeometryType::CoordinatesArrayType SlaveGPGlobal;
    this->GetGeometry( ).GlobalCoordinates( SlaveGPGlobal, LocalPoint );
    ContactUtilities::FastProjectDirection( MasterSegment, SlaveGPGlobal, ProjectedGPGlobal, MasterNormal, -GPNormal ); // The opposite direction
    
    GeometryType::CoordinatesArrayType ProjectedGPLocal;
    
    MasterSegment.PointLocalCoordinates( ProjectedGPLocal, ProjectedGPGlobal.Coordinates( ) ) ;

    // SHAPE FUNCTIONS 
    MasterSegment.ShapeFunctionsValues(         rVariables.N_Master,     ProjectedGPLocal );         
    MasterSegment.ShapeFunctionsLocalGradients( rVariables.DN_De_Master, ProjectedGPLocal );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateMortarOperators(
    MortarConditionMatrices& rThisMortarConditionMatrices,
    GeneralVariables& rVariables,
    DerivativeData& rDerivativeData,
    const double& rIntegrationWeight
    )
{
    /* DEFINITIONS */
    const double DetJSlave = rVariables.DetJSlave; 
    const VectorType Phi = rVariables.Phi_LagrangeMultipliers;
    const VectorType N1  = rVariables.N_Slave;
    const VectorType N2  = rVariables.N_Master;
    
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
    
    for (unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++)
    {
        const double phi = Phi[i_slave];
        
        for (unsigned int j_slave = 0; j_slave < TNumNodes; j_slave++)
        {
            const double n1  = N1[j_slave];
            const double n2  = N2[j_slave];
            
            DOperator(i_slave, j_slave) += DetJSlave * rIntegrationWeight * phi * n1;
            MOperator(i_slave, j_slave) += DetJSlave * rIntegrationWeight * phi * n2;
            
            for (unsigned int i = 0; i < TDim * TNumNodes; i++)
            {
                DeltaDOperator[i](i_slave, j_slave) += DeltaJSlave[i] * rIntegrationWeight * phi* n1        
                                                     + DetJSlave * rIntegrationWeight * DeltaPhi[i][i_slave] * n1
                                                     + DetJSlave * rIntegrationWeight * phi* DeltaN1[i][j_slave];
                                                                            
                DeltaMOperator[i](i_slave, j_slave) += DeltaJSlave[i] * rIntegrationWeight * phi* n2        
                                                     + DetJSlave * rIntegrationWeight * DeltaPhi[i][i_slave] * n2
                                                     + DetJSlave * rIntegrationWeight * phi* DeltaN2[i][j_slave];
            }
            for (unsigned int i = TDim * TNumNodes; i < 2 * TDim * TNumNodes; i++)
            {
                DeltaDOperator[i](i_slave, j_slave) += DetJSlave * rIntegrationWeight * phi * DeltaN1[i][j_slave];
                                                                            
                DeltaMOperator[i](i_slave, j_slave) += DetJSlave * rIntegrationWeight * phi * DeltaN2[i][j_slave];
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
    const VectorType Phi = rVariables.Phi_LagrangeMultipliers;
    const VectorType N1  = rVariables.N_Slave;
    const VectorType N2  = rVariables.N_Master;

    // Mortar condition matrices - DOperator and MOperator
    bounded_matrix<double, TNumNodes, TNumNodes>& DOperator = rThisMortarConditionMatrices.DOperator;
    bounded_matrix<double, TNumNodes, TNumNodes>& MOperator = rThisMortarConditionMatrices.MOperator;
    
    for (unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++)
    {
        const double phi = Phi[i_slave];
        
        for (unsigned int j_slave = 0; j_slave < TNumNodes; j_slave++)
        {
            DOperator(i_slave, j_slave) += DetJSlave * rIntegrationWeight * phi * N1[j_slave];
            MOperator(i_slave, j_slave) += DetJSlave * rIntegrationWeight * phi * N2[j_slave];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateDeltaAeComponents(
    GeneralVariables& rVariables,
    DerivativeData& rDerivativeData,
    AeData& rAeData,
    const double& rIntegrationWeight
    )
{
    /* DEFINITIONS */
    const VectorType N1 = rVariables.N_Slave;
    const double detJ   = rVariables.DetJSlave; 
     
    rAeData.De += rIntegrationWeight * this->ComputeDe( N1, detJ);
    rAeData.Me += rIntegrationWeight * this->ComputeMe( N1, detJ);
    
    for (unsigned int i = 0; i < TDim * TNumNodes; i++)
    {
        const double DeltaDetJ = rDerivativeData.DeltaJSlave[i];
        
        const bounded_matrix<double, TNumNodes, TNumNodes> DeltaMe = DeltaDetJ * outer_prod(N1, N1);
        bounded_matrix<double, TNumNodes, TNumNodes> DeltaDe;
        
        for (unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++)
        {
            for (unsigned int j_slave = i_slave; j_slave < TNumNodes; j_slave++)
            {
                if (i_slave == j_slave)
                {
                    DeltaDe(i_slave, i_slave) = DeltaDetJ * N1[i_slave];
                }
                else
                {
                    DeltaDe(i_slave, j_slave) = 0.0;
                    DeltaDe(j_slave, i_slave) = 0.0;
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
    DerivativeData& rDerivativeData,
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
        const unsigned int& rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
    
    return ZeroMatrix(28, 28);
}

// /***********************************************************************************/
// /***********************************************************************************/
// 
// template<>
// bounded_matrix<double, 12, 12> AugmentedLagrangianMethodMortarContactCondition<2,2, true>::CalculateLocalLHS(
//         const MortarConditionMatrices& rMortarConditionMatrices,
//         const unsigned int& rMasterElementIndex,
//         const unsigned int& rActiveInactive
//         )
// {
//     KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
//     
//     return ZeroMatrix(12, 12);
// }
// 
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
        const unsigned int& rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
    
    return ZeroVector(28);
}

// /***********************************************************************************/
// /***********************************************************************************/
// 
// template<>
// array_1d<double,12> AugmentedLagrangianMethodMortarContactCondition<2,2, true>::CalculateLocalRHS(
//         const MortarConditionMatrices& rMortarConditionMatrices,
//         const unsigned int& rMasterElementIndex,
//         const unsigned int& rActiveInactive
//         )
// {
//     KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
//     
//     return ZeroVector(12);
// }
// 
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
   DerivativeData& rDerivativeData
   ) // TODO: Do an explicit specialization!!!!
{
    if (TDim == 2)
    {
        // Fill up the elements corresponding to the slave DOFs - the rest remains zero
        for ( unsigned int i_slave = 0, i = 0; i_slave < TNumNodes; ++i_slave, i += TDim )
        {
            rDerivativeData.DeltaJSlave[i    ] = rVariables.j_Slave( 0, 0 ) * rVariables.DN_De_Slave( i_slave, 0) / rVariables.DetJSlave;
            rDerivativeData.DeltaJSlave[i + 1] = rVariables.j_Slave( 1, 0 ) * rVariables.DN_De_Slave( i_slave, 0) / rVariables.DetJSlave;
        }
    }
    else
    {
        const array_1d<double,TNumNodes>& DN_Dxi  = column( rVariables.DN_De_Slave, 0 );
        const array_1d<double,TNumNodes>& DN_Deta = column( rVariables.DN_De_Slave, 1 );
        
        const array_1d<double,TDim>& J_xi    = column( rVariables.j_Slave, 0 );
        const array_1d<double,TDim>& J_eta   = column( rVariables.j_Slave, 1 );
        
        const array_1d<double,TDim>& n = prod(trans(rDerivativeData.NormalMaster), rVariables.N_Slave); // FIXME: Check this!!!!
        
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
            
            rDerivativeData.DeltaJSlave[i    ] = inner_prod( n, column( Delta_Jxi_x_Jeta, 0 ) );
            rDerivativeData.DeltaJSlave[i + 1] = inner_prod( n, column( Delta_Jxi_x_Jeta, 1 ) );
            rDerivativeData.DeltaJSlave[i + 2] = inner_prod( n, column( Delta_Jxi_x_Jeta, 2 ) );
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

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional> // NOTE: Formulation taken from Mohamed Khalil work
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateDeltaNormalSlave(DerivativeData& rDerivativeData)
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
                    row(rDerivativeData.DeltaNormalSlave[i_slave * TDim + i_dof], i_node)      =   trans(column(DeltaNormal, i_dof)); 
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
                    row(rDerivativeData.DeltaNormalSlave[i_slave * TDim + i_dof], i_node) = trans(column(DeltaNormal, i_dof)); 
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional> // NOTE: Formulation taken from Mohamed Khalil work
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateDeltaNormalMaster(DerivativeData& rDerivativeData)
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
                    row(rDerivativeData.DeltaNormalMaster[i_master * TDim + i_dof], i_node) = trans(column(DeltaNormal, i_dof)); 
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
                    row(rDerivativeData.DeltaNormalMaster[i_master * TDim + i_dof], i_node)  = trans(column(DeltaNormal, i_dof)); 
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
      for (unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++)
      {
         for (unsigned int i_dim = 0; i_dim < TDim; i_dim++)
         {
            const unsigned int i_dof = i_slave * TDim + i_dim;
            
            const array_1d<double, TNumNodes > DNormalSlaveg = prod(trans(DNormalSlave[i_dof]), N1);
            array_1d<double, TNumNodes > aux_vector = ZeroVector(TNumNodes);
            aux_vector[i_dim] = 1.0;
            const array_1d<double, TNumNodes > Deltavector_nodes = - N1[i_slave] * aux_vector;
            const double DeltaDist = ((inner_prod(Deltavector_nodes, NormalMasterg))* inner_prod(NormalSlaveg, NormalMasterg) - inner_prod(vector_nodes, NormalMasterg) * (inner_prod(DNormalSlaveg, NormalMasterg)))/std::pow(inner_prod(NormalSlaveg, NormalMasterg) + tol, 2);
            const array_1d<double, TNumNodes > Deltax1g = N1[i_slave] * aux_vector;
            const array_1d<double, TNumNodes > Deltax2g = Deltax1g + DeltaDist * NormalSlaveg + Dist * DNormalSlaveg;
            
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
            
            const array_1d<double, TNumNodes > DNormalMasterg = prod(trans(DNormalMaster[i_dof - TNumNodes * TDim]), N2);
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

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateDeltaPhi(
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
// template class AugmentedLagrangianMethodMortarContactCondition<2, 2, true>;
// template class AugmentedLagrangianMethodMortarContactCondition<3, 3, true>;
// template class AugmentedLagrangianMethodMortarContactCondition<3, 4, true>;

} // Namespace Kratos
