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
#include "custom_utilities/search_utilities.h"
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
= default;

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
    boost::shared_ptr<ConditionMap>& all_conditions_maps = this->GetValue( CONTACT_MAPS );
    mPairSize = all_conditions_maps->size();
    mThisMasterElements.resize( mPairSize );
    mThisMasterElementsActive.resize( mPairSize );
    
    unsigned int i_cond = 0;
    for (auto it_pair = all_conditions_maps->begin(); it_pair != all_conditions_maps->end(); ++it_pair )
    {
        mThisMasterElements[i_cond] = (it_pair->first);
        mThisMasterElementsActive[i_cond] = (it_pair->second);
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
    
    if (rCurrentProcessInfo[CONSIDER_PAIR_VARIATION] == true)
    {
        // We update the active/inactive pair
        boost::shared_ptr<ConditionMap>& all_conditions_maps = this->GetValue( CONTACT_MAPS );
        
        unsigned int i_cond = 0;
        for (auto it_pair = all_conditions_maps->begin(); it_pair != all_conditions_maps->end(); ++it_pair )
        {
            mThisMasterElementsActive[i_cond] = (it_pair->second);
            i_cond += 1;
        }
    }
        
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
    
    this->Set(VISITED, true);
    
    if (rCurrentProcessInfo[CONSIDER_PAIR_VARIATION] == true)
    {
        // Check pairs
        boost::shared_ptr<ConditionMap>& all_conditions_maps = this->GetValue( CONTACT_MAPS );
        GeometryType& this_geometry = GetGeometry();
        const double active_check_length = this_geometry.Length() * GetProperties().GetValue(ACTIVE_CHECK_FACTOR);
        SearchUtilities::ExactContactContainerChecker<TDim,TNumNodes>(all_conditions_maps, this_geometry, this->GetValue(NORMAL), active_check_length); 
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
    LocalSystemComponents local_system;

    // Calculation flags
    local_system.CalculationFlags.Set(AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS, true);
    local_system.CalculationFlags.Set(AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS, true);

    //Initialize sizes for the system components
    if ( rLHSVariables.size( ) != rLeftHandSideMatrices.size( ) )
    {
        rLeftHandSideMatrices.resize( rLHSVariables.size( ) );
    }

    if ( rRHSVariables.size( ) != rRightHandSideVectors.size( ) )
    {
        rRightHandSideVectors.resize( rRHSVariables.size( ) );
    }

    local_system.CalculationFlags.Set(AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX, true);
    for ( unsigned int i = 0; i < rLeftHandSideMatrices.size( ); i++ )
    {
        // Note: rRightHandSideVectors.size() > 0
        this->InitializeSystemMatrices( rLeftHandSideMatrices[i], rRightHandSideVectors[0],local_system.CalculationFlags );
    }

    local_system.CalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_RHS_VECTOR, true );
    local_system.CalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX, false ); // Temporarily only
    for ( unsigned int i = 0; i < rRightHandSideVectors.size( ); i++ )
    {
        // Note: rLeftHandSideMatrices.size() > 0
        this->InitializeSystemMatrices( rLeftHandSideMatrices[0], rRightHandSideVectors[i], local_system.CalculationFlags  );
    }
    local_system.CalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX, true ); // Reactivated again

    // Set Variables to Local system components
    local_system.SetLeftHandSideMatrices( rLeftHandSideMatrices );
    local_system.SetRightHandSideVectors( rRightHandSideVectors );

    local_system.SetLeftHandSideVariables( rLHSVariables );
    local_system.SetRightHandSideVariables( rRHSVariables );

    // Calculate condition system
    this->CalculateConditionSystem( local_system, rCurrentProcessInfo );
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
    LocalSystemComponents local_system;

    // Calculation flags
    local_system.CalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX, true );
    local_system.CalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_RHS_VECTOR, true );

    // Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, rRightHandSideVector, local_system.CalculationFlags );
    
    // Set Variables to Local system components
    local_system.SetLeftHandSideMatrix( rLeftHandSideMatrix );
    local_system.SetRightHandSideVector( rRightHandSideVector );

    // Calculate condition system
    this->CalculateConditionSystem( local_system, rCurrentProcessInfo );

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
    LocalSystemComponents local_system;

    // Calculation flags
    local_system.CalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX, true );

    VectorType right_hand_side_vector = Vector( );

    // Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, right_hand_side_vector, local_system.CalculationFlags );

    // Set Variables to Local system components
    local_system.SetLeftHandSideMatrix( rLeftHandSideMatrix );
    local_system.SetRightHandSideVector( right_hand_side_vector );

    // Calculate condition system
    this->CalculateConditionSystem( local_system, rCurrentProcessInfo );
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
    LocalSystemComponents local_system;

    // Calculation flags
    local_system.CalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX, true );
    local_system.CalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS, true );

    VectorType right_hand_side_vector = Vector( );

    // Initialize size for the system components
    for( unsigned int i = 0; i < rLeftHandSideMatrices.size( ); i++ )
    {
        this->InitializeSystemMatrices( rLeftHandSideMatrices[i], right_hand_side_vector, local_system.CalculationFlags );
    }

    local_system.SetLeftHandSideMatrices( rLeftHandSideMatrices );
    local_system.SetRightHandSideVector( right_hand_side_vector );

    // Calculate condition system
    this->CalculateConditionSystem( local_system, rCurrentProcessInfo );
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
    LocalSystemComponents local_system;

    // Calculation flags
    local_system.CalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_RHS_VECTOR, true);

    MatrixType left_hand_side_matrix = Matrix( );

    // Initialize size for the system components
    this->InitializeSystemMatrices( left_hand_side_matrix, rRightHandSideVector,local_system.CalculationFlags);

    //Set Variables to Local system components
    local_system.SetLeftHandSideMatrix( left_hand_side_matrix );
    local_system.SetRightHandSideVector( rRightHandSideVector );

    // Calculate condition system
    this->CalculateConditionSystem( local_system, rCurrentProcessInfo );
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
    LocalSystemComponents local_system;

    // Calculation flags
    local_system.CalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_RHS_VECTOR, true );
    local_system.CalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS, true );

    MatrixType left_hand_side_matrix = Matrix( );

    // Initialize size for the system components
    for( unsigned int i = 0; i < rRightHandSideVectors.size(); i++ )
    {
        this->InitializeSystemMatrices( left_hand_side_matrix, rRightHandSideVectors[i], local_system.CalculationFlags );
    }

    // Set Variables to Local system components
    local_system.SetLeftHandSideMatrix( left_hand_side_matrix );
    local_system.SetRightHandSideVectors( rRightHandSideVectors );

    // Calculate condition system
    this->CalculateConditionSystem( local_system, rCurrentProcessInfo );
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
    const unsigned int condition_size = this->CalculateConditionSize( );
    
    // Resizing as needed the LHS
    if ( rCalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX ) ) // Calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != condition_size )
        {
            rLeftHandSideMatrix.resize( condition_size, condition_size, false );
        }
        noalias( rLeftHandSideMatrix ) = ZeroMatrix( condition_size, condition_size ); // Resetting LHS
    }

    // Resizing as needed the RHS
    if ( rCalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_RHS_VECTOR ) ) // Calculation of the matrix is required
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
    const unsigned int condition_size = mPairSize * MatrixSize;
    
    return condition_size;
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
        
    // Create and initialize condition variables
    GeneralVariables rVariables;
    
    // Create the current contact data
    DerivativeDataType rDerivativeData;
    rDerivativeData.Initialize(this->GetGeometry(), rCurrentProcessInfo);
    
    // We compute the normal derivatives
    if (rCurrentProcessInfo[CONSIDER_NORMAL_VARIATION] == true)
    {
        // Compute the normal derivatives of the master
        this->CalculateDeltaNormalSlave(rDerivativeData);
    }
    
    // Create the mortar operators
    MortarConditionMatrices rThisMortarConditionMatrices;
    
    // We call the exact integration utility
    ExactMortarIntegrationUtility<TDim, TNumNodes, true>  integration_utility = ExactMortarIntegrationUtility<TDim, TNumNodes, true> (mIntegrationOrder);

//     // DEBUG
//     std::cout << "\nGraphics3D[{EdgeForm[{Thick,Dashed,Red}],FaceForm[],Triangle[{{" << GetGeometry()[0].X() << "," << GetGeometry()[0].Y() << "," << GetGeometry()[0].Z()  << "},{" << GetGeometry()[1].X() << "," << GetGeometry()[1].Y() << "," << GetGeometry()[1].Z()  << "},{" << GetGeometry()[2].X() << "," << GetGeometry()[2].Y() << "," << GetGeometry()[2].Z()  << "}}],Text[Style["<< this->Id() <<", Tiny],{"<< GetGeometry().Center().X() << "," << GetGeometry().Center().Y() << ","<< GetGeometry().Center().Z() << "}]}],";// << std::endl;
    
    // Iterate over the master segments
    for (unsigned int pair_index = 0; pair_index < mPairSize; ++pair_index)
    {   
        if (mThisMasterElementsActive[pair_index] == true)
        {
//             // DEBUG
//             if (mThisMasterElements[pair_index]->Is(VISITED) == false || mThisMasterElements[pair_index]->IsDefined(VISITED) == false)
//             {
//                 std::cout << "\nGraphics3D[{EdgeForm[{Thick,Dashed,Blue}],FaceForm[],Triangle[{{" << mThisMasterElements[pair_index]->GetGeometry()[0].X() << "," << mThisMasterElements[pair_index]->GetGeometry()[0].Y() << "," << mThisMasterElements[pair_index]->GetGeometry()[0].Z()  << "},{" << mThisMasterElements[pair_index]->GetGeometry()[1].X() << "," << mThisMasterElements[pair_index]->GetGeometry()[1].Y() << "," << mThisMasterElements[pair_index]->GetGeometry()[1].Z()  << "},{" << mThisMasterElements[pair_index]->GetGeometry()[2].X() << "," << mThisMasterElements[pair_index]->GetGeometry()[2].Y() << "," << mThisMasterElements[pair_index]->GetGeometry()[2].Z()  << "}}],Text[Style["<< mThisMasterElements[pair_index]->Id() <<", Tiny],{"<< mThisMasterElements[pair_index]->GetGeometry().Center().X() << "," << mThisMasterElements[pair_index]->GetGeometry().Center().Y() << ","<< mThisMasterElements[pair_index]->GetGeometry().Center().Z() << "}]}],";// << std::endl;
//                 
//                 mThisMasterElements[pair_index]->Set(VISITED, true);
//             }
            
            // The normal of the master condition
            const array_1d<double, 3>& master_normal = mThisMasterElements[pair_index]->GetValue(NORMAL);
            
            // Reading integration points
            ConditionArrayListType conditions_points_slave;
            const bool is_inside = integration_utility.GetExactIntegration(this->GetGeometry(), this->GetValue(NORMAL), mThisMasterElements[pair_index]->GetGeometry(), master_normal, conditions_points_slave);
            
            if (is_inside == true)
            {            
                IntegrationMethod this_integration_method = GetIntegrationMethod();
                
                // Initialize general variables for the current master element
                rVariables.Initialize();
                
                // Update slave element info
                rDerivativeData.UpdateMasterPair(mThisMasterElements[pair_index]);
                
                // Initialize the mortar operators
                rThisMortarConditionMatrices.Initialize();
                
//                 if (rCurrentProcessInfo[CONSIDER_NORMAL_VARIATION] == true && TDim == 2) // TODO: Can be needed for the shape function derivatives
//                 {
//                     // Compute the normal derivatives of the master
//                     this->CalculateDeltaNormalMaster(rDerivativeData, mThisMasterElements[pair_index]->GetGeometry());
//                 }
                
                const bool dual_LM = this->CalculateAeAndDeltaAe(rDerivativeData, rVariables, rCurrentProcessInfo, pair_index, conditions_points_slave, this_integration_method, master_normal);
                
                for (unsigned int i_geom = 0; i_geom < conditions_points_slave.size(); i_geom++)
                {
                    std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
                    array_1d<BelongType, TDim> belong_array;
                    for (unsigned int i_node = 0; i_node < TDim; i_node++)
                    {
                        PointType global_point;
                        GetGeometry().GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                        points_array[i_node] = boost::make_shared<PointType>(global_point);
                        belong_array[i_node] = conditions_points_slave[i_geom][i_node].GetBelong();
                    }
                    
                    DecompositionType decomp_geom( points_array );
                    
//                     // DEBUG
//                     std::cout << "\nGraphics3D[{Opacity[.3],Triangle[{{" << decomp_geom[0].X() << "," << decomp_geom[0].Y() << "," << decomp_geom[0].Z()  << "},{" << decomp_geom[1].X() << "," << decomp_geom[1].Y() << "," << decomp_geom[1].Z()  << "},{" << decomp_geom[2].X() << "," << decomp_geom[2].Y() << "," << decomp_geom[2].Z()  << "}}]}],";// << std::endl;
                    
                    const bool bad_shape = (TDim == 2) ? ContactUtilities::LengthCheck(decomp_geom, this->GetGeometry().Length() * 1.0e-6) : ContactUtilities::HeronCheck(decomp_geom);
                    
                    if (bad_shape == false)
                    {
                        Matrix delta_position;
                        delta_position = CalculateDeltaPosition(delta_position, conditions_points_slave[i_geom]);
                        
                        const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( this_integration_method );
                        
                        // Integrating the mortar operators
                        for ( unsigned int point_number = 0; point_number < integration_points_slave.size(); point_number++ )
                        {
                            const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
                            PointType local_point_parent;
                            PointType gp_global;
                            decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                            GetGeometry().PointLocalCoordinates(local_point_parent, gp_global);
                            
                            // Calculate the kinematic variables
                            this->CalculateKinematics( rVariables, rDerivativeData, master_normal, pair_index, local_point_decomp, local_point_parent, decomp_geom, dual_LM, delta_position);
                            
                            const double integration_weight = GetIntegrationWeight(rVariables, integration_points_slave, point_number);
                            
                            if ( rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX ) ||
                                    rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
                            {
                                /* Update the derivatives */
                                // Update the derivative of the integration vertex (just in 3D)
                                if (TDim == 3) this->CalculateDeltaCellVertex(rVariables, rDerivativeData, belong_array);
                                // Update the derivative of DetJ
                                this->CalculateDeltaDetjSlave(rVariables, rDerivativeData);
                                // Update the derivatives of the shape functions and the gap
//                                 this->CalculateDeltaN(rVariables, rDerivativeData); // FIXME: This is the old version!!!!
                                // The derivatives of the dual shape function 
                                this->CalculateDeltaPhi(rVariables, rDerivativeData);
                                
                                rThisMortarConditionMatrices.CalculateDeltaMortarOperators(rVariables, rDerivativeData, integration_weight);    
                            }
                            else // In case we are computing RHS we don't compute derivatives (not necessary)
                            {
                                rThisMortarConditionMatrices.CalculateMortarOperators(rVariables, integration_weight);   
                            }
                        }
                    }
                }
                        
//                 // Debug
//                 std::cout << "--------------------------------------------------" << std::endl;
//                 KRATOS_WATCH(this->Id());
//                 KRATOS_WATCH(pair_index);
//                 rThisMortarConditionMatrices.print();
                
                // Calculates the active/inactive combination pair
                const unsigned int active_inactive = GetActiveInactiveValue(this->GetGeometry());
                
                // Assemble of the matrix is required
                if ( rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX ) ||
                        rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
                {
                    // Calculate the local contribution
                    const bounded_matrix<double, MatrixSize, MatrixSize> LHS_contact_pair = this->CalculateLocalLHS( rThisMortarConditionMatrices, rDerivativeData, active_inactive);
                    
                    // Contributions to stiffness matrix calculated on the reference config
                    this->CalculateAndAddLHS( rLocalSystem, LHS_contact_pair, pair_index );
                    
//                     // Debug
//         //                 KRATOS_WATCH(LHS_contact_pair);
//                     LOG_MATRIX_PRETTY( LHS_contact_pair );
                }
                
                // Assemble of the vector is required
                if ( rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_RHS_VECTOR ) ||
                        rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
                {
                    // Calculate the local contribution
                    const array_1d<double, MatrixSize> RHS_contact_pair = this->CalculateLocalRHS( rThisMortarConditionMatrices, rDerivativeData, active_inactive);
                    
                    // Contribution to previous step contact force and residuals vector
                    this->CalculateAndAddRHS( rLocalSystem, RHS_contact_pair, pair_index );
                    
//                     // Debug
//         //                 KRATOS_WATCH(RHS_contact_pair);
//                     LOG_VECTOR_PRETTY( RHS_contact_pair );
                }
            }
        }
    }
    
    // Reseting flag
    if ((this)->Is(VISITED) == true)
    {
        (this)->Set(VISITED, false);
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
bool AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateAeAndDeltaAe(
    DerivativeDataType& rDerivativeData,
    GeneralVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo,
    const unsigned int PairIndex,
    ConditionArrayListType& ConditionsPointsSlave,
    IntegrationMethod ThisIntegrationMethod,
    const array_1d<double, 3>& MasterNormal
    )
{
    // We initilize the Ae components
    AeData rAeData;
    rAeData.Initialize();
    
    rDerivativeData.InitializeDeltaAeComponents();

    // Initialize general variables for the current master element
    rVariables.Initialize();
    
    // Update slave element info
    rDerivativeData.UpdateMasterPair(mThisMasterElements[PairIndex]);
    
    for (unsigned int i_geom = 0; i_geom < ConditionsPointsSlave.size(); i_geom++)
    {
        std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
        array_1d<BelongType, TDim> belong_array;
        for (unsigned int i_node = 0; i_node < TDim; i_node++)
        {
            PointType global_point;
            GetGeometry().GlobalCoordinates(global_point, ConditionsPointsSlave[i_geom][i_node]);
            points_array[i_node] = boost::make_shared<PointType>(global_point);
        }
        
        DecompositionType decomp_geom( points_array );
        
        const bool bad_shape = (TDim == 2) ? ContactUtilities::LengthCheck(decomp_geom, this->GetGeometry().Length() * 1.0e-6) : ContactUtilities::HeronCheck(decomp_geom);
        
        if (bad_shape == false)
        {
            Matrix delta_position;
            delta_position = CalculateDeltaPosition(delta_position, ConditionsPointsSlave[i_geom]);
            
            const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( ThisIntegrationMethod );
            
            // Integrating the mortar operators
            for ( unsigned int point_number = 0; point_number < integration_points_slave.size(); point_number++ )
            {
                const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
                PointType local_point_parent;
                PointType gp_global;
                decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                GetGeometry().PointLocalCoordinates(local_point_parent, gp_global);
                
                // Calculate the kinematic variables
                this->CalculateKinematics( rVariables, rDerivativeData, MasterNormal, PairIndex, local_point_decomp, local_point_parent, decomp_geom, false, delta_position);
                
                // Update the derivative of the integration vertex (just in 3D)
                if (TDim == 3) this->CalculateDeltaCellVertex(rVariables, rDerivativeData, belong_array);
                                
                // Update the derivative of DetJ
                this->CalculateDeltaDetjSlave(rVariables, rDerivativeData); 
                
                // Integrate
                const double integration_weight = GetIntegrationWeight(rVariables, integration_points_slave, point_number);
                
                rAeData.CalculateDeltaAeComponents(rVariables, rDerivativeData, integration_weight);
            }
        }
    }
    
    return rAeData.CalculateDeltaAe(rDerivativeData);
}

/*********************************COMPUTE KINEMATICS*********************************/
/************************************************************************************/
 
template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateKinematics( 
    GeneralVariables& rVariables,
    const DerivativeDataType rDerivativeData,
    const array_1d<double, 3> MasterNormal,
    const unsigned int PairIndex,
    const PointType& LocalPointDecomp,
    const PointType& LocalPointParent,
    GeometryPointType& GeometryDecomp,
    const bool DualLM,
    Matrix DeltaPosition
    )
{  
    /// SLAVE CONDITION ///
    /* CALCULATE JACOBIAN AND JACOBIAN DETERMINANT */
    rVariables.jSlave = GeometryDecomp.Jacobian( rVariables.jSlave, LocalPointDecomp.Coordinates(), DeltaPosition);
    rVariables.DetjSlave = MathUtils<double>::GeneralizedDet(rVariables.jSlave);
//     rVariables.DetjSlave = GeometryDecomp.DeterminantOfJacobian( LocalPointDecomp ); // TODO: Add this to the geometry.h
    
    /* SHAPE FUNCTIONS */
    GetGeometry().ShapeFunctionsValues( rVariables.NSlave, LocalPointParent.Coordinates() );
    if (DualLM == true)
    {
        rVariables.PhiLagrangeMultipliers = prod(rDerivativeData.Ae, rVariables.NSlave);
    }
    else
    {
        rVariables.PhiLagrangeMultipliers = rVariables.NSlave;
    }
    
    /* SHAPE FUNCTION DERIVATIVES */
    GetGeometry().ShapeFunctionsLocalGradients( rVariables.DNDeSlave, LocalPointParent );
    
    /// MASTER CONDITION ///
    this->MasterShapeFunctionValue( rVariables, MasterNormal, LocalPointParent, PairIndex);
}
 
/***********************************************************************************/
/*************** METHODS TO CALCULATE THE CONTACT CONDITION MATRICES ***************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::MasterShapeFunctionValue(
    GeneralVariables& rVariables,
    const array_1d<double, 3> MasterNormal,
    const PointType& LocalPoint,
    const unsigned int PairIndex
    )
{    
    GeometryType& master_geometry = mThisMasterElements[PairIndex]->GetGeometry();

    PointType projected_gp_global;
    const array_1d<double,3> gp_normal = ContactUtilities::GaussPointNormal(rVariables.NSlave, GetGeometry());
    
    GeometryType::CoordinatesArrayType slave_gp_global;
    this->GetGeometry( ).GlobalCoordinates( slave_gp_global, LocalPoint );
    ContactUtilities::FastProjectDirection( master_geometry, slave_gp_global, projected_gp_global, MasterNormal, -gp_normal ); // The opposite direction
    
    GeometryType::CoordinatesArrayType projected_gp_local;
    
    master_geometry.PointLocalCoordinates( projected_gp_local, projected_gp_global.Coordinates( ) ) ;

    // SHAPE FUNCTIONS 
    master_geometry.ShapeFunctionsValues(         rVariables.NMaster,    projected_gp_local );         
    master_geometry.ShapeFunctionsLocalGradients( rVariables.DNDeMaster, projected_gp_local );
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
        const DerivativeDataType& rDerivativeData,
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
        const DerivativeDataType& rDerivativeData,
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
        const DerivativeDataType& rDerivativeData,
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
        const DerivativeDataType& rDerivativeData,
        const unsigned int& rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
    
    return ZeroMatrix(12, 12);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
bounded_matrix<double, 27, 27> AugmentedLagrangianMethodMortarContactCondition<3, 3, true>::CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int& rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
    
    return ZeroMatrix(27, 27);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
bounded_matrix<double, 36, 36> AugmentedLagrangianMethodMortarContactCondition<3, 4, true>::CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int& rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
    
    return ZeroMatrix(36, 36);
}

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
array_1d<double,10> AugmentedLagrangianMethodMortarContactCondition<2, 2, false>::CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int& rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
    
    return ZeroVector(10);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double,21> AugmentedLagrangianMethodMortarContactCondition<3, 3, false>::CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int& rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
    
    return ZeroVector(21);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double,28> AugmentedLagrangianMethodMortarContactCondition<3, 4, false>::CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int& rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
    
    return ZeroVector(28);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double,12> AugmentedLagrangianMethodMortarContactCondition<2, 2, true>::CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int& rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
    
    return ZeroVector(12);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double,27> AugmentedLagrangianMethodMortarContactCondition<3, 3, true>::CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int& rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
    
    return ZeroVector(27);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double,36> AugmentedLagrangianMethodMortarContactCondition<3, 4, true>::CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int& rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
    
    return ZeroVector(36);
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateDeltaCellVertex(
   GeneralVariables& rVariables,
   DerivativeDataType& rDerivativeData,
   const array_1d<BelongType, TDim>& TheseBelongs
   ) 
{

}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateDeltaDetjSlave(
   GeneralVariables& rVariables,
   DerivativeDataType& rDerivativeData
   )
{
    if (TDim == 2)
    {
        // Fill up the elements corresponding to the slave DOFs - the rest remains zero
        for ( unsigned int i_slave = 0, i = 0; i_slave < TNumNodes; ++i_slave, i += TDim )
        {
            rDerivativeData.DeltaDetjSlave[i    ] = rVariables.jSlave( 0, 0 ) * rVariables.DNDeSlave( i_slave, 0) / rVariables.DetjSlave;
            rDerivativeData.DeltaDetjSlave[i + 1] = rVariables.jSlave( 1, 0 ) * rVariables.DNDeSlave( i_slave, 0) / rVariables.DetjSlave;
        }
    }
    else
    {
        const array_1d<double,TNumNodes>& DNDxi  = column( rVariables.DNDeSlave, 0 );
        const array_1d<double,TNumNodes>& DNDeta = column( rVariables.DNDeSlave, 1 );
        
        const array_1d<double,TDim>& Jxi  = column( rVariables.jSlave, 0 );
        const array_1d<double,TDim>& Jeta = column( rVariables.jSlave, 1 );
        
        const array_1d<double,TDim>& Normal = prod(trans(rDerivativeData.NormalMaster), rVariables.NSlave);
        
        bounded_matrix<double, TDim, TDim> DeltaJxixJeta;
        
        for ( unsigned int i_slave = 0, i = 0; i_slave < TNumNodes; ++i_slave, i += TDim )
        {
            DeltaJxixJeta(0,0) = 0.0;
            DeltaJxixJeta(0,1) =  Jeta(2) * DNDxi(i_slave) - Jxi(2) * DNDeta(i_slave); 
            DeltaJxixJeta(0,2) = -Jeta(1) * DNDxi(i_slave) + Jxi(1) * DNDeta(i_slave); 
            DeltaJxixJeta(1,0) = -Jeta(2) * DNDxi(i_slave) + Jxi(2) * DNDeta(i_slave); 
            DeltaJxixJeta(1,1) = 0.0;
            DeltaJxixJeta(1,2) =  Jeta(0) * DNDxi(i_slave) - Jxi(0) * DNDeta(i_slave);
            DeltaJxixJeta(2,0) =  Jeta(1) * DNDxi(i_slave) - Jxi(1) * DNDeta(i_slave); 
            DeltaJxixJeta(2,1) = -Jeta(0) * DNDxi(i_slave) + Jxi(0) * DNDeta(i_slave); 
            DeltaJxixJeta(2,2) = 0.0;
            
            rDerivativeData.DeltaDetjSlave[i    ] = inner_prod( Normal, column( DeltaJxixJeta, 0 ) );
            rDerivativeData.DeltaDetjSlave[i + 1] = inner_prod( Normal, column( DeltaJxixJeta, 1 ) );
            rDerivativeData.DeltaDetjSlave[i + 2] = inner_prod( Normal, column( DeltaJxixJeta, 2 ) );
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional> // NOTE: Formulation taken from Mohamed Khalil work
bounded_matrix<double, TDim, TDim> AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::LocalDeltaNormal( // NOTE: Not the mean, look in the contact utilities 
    const GeometryType& CondGeometry,
    const unsigned int NodeIndex
    )
{
    // Tolerance
    const double Tolerance = 1.0e-14;
        
    bounded_matrix<double, TDim, TDim> DeltaNeAdj;
    bounded_matrix<double, TDim, TDim> Ce;
    
    const bounded_matrix<double, TDim, TDim> I = IdentityMatrix(TDim, TDim);
    
    bounded_matrix<double, TDim, TDim> delta_normal = ZeroMatrix(TDim,TDim);
    
    const array_1d<double, 3>& Ne = this->GetValue(NORMAL); // Normalized condition normal
    bounded_matrix<double, TDim, TDim> NeoNe = subrange( outer_prod( Ne, Ne ), 0, TDim, 0, TDim );
    
    // Auxiliar value
    double NeNorm;
    
    if (TDim == 2)
    {
        NeNorm = CondGeometry.Length( ); // The norm of a geometry's normal is its characteristic dimension - length for 2D and area for 3D 
                        
        DeltaNeAdj( 0, 0 ) =  0.0;
        DeltaNeAdj( 0, 1 ) = -1.0;
        DeltaNeAdj( 1, 0 ) =  1.0;
        DeltaNeAdj( 1, 1 ) =  0.0;
        
        double DNDej = 0.0;
        if( NodeIndex == 0 )
        {
            DNDej = - 0.5;
        }
        else if( NodeIndex == 1 )
        {
            DNDej =   0.5;
        }
        
        Ce = prod( I - NeoNe, DeltaNeAdj ) / NeNorm; // In 2D, DeltaNeAdj is node-independent => evaluated outside the nodes loop
        
        delta_normal = - 2.0 * Ce * DNDej; // NOTE: Check why - 2???!!!, it was the only wayto ensure the same value as the symbolic. You will need to repeat this in 3D            
    //         delta_normal = Ce * DNDej;     
    }
    else
    {
        NeNorm = CondGeometry.Area( ); // The norm of a geometry's normal is its characteristic dimension - length for 2D and area for 3D 
        
        MatrixType J = ZeroMatrix( 3, 2 ); // Jacobian [ 3D global x 2D local ]
        array_1d<double, 2> DNDej;
        array_1d<double, 3> LocalCoordsj;
        
        if( TNumNodes == 3 )    // linear triangle element
        {
            if( NodeIndex == 0 )
            {
                LocalCoordsj[0] = 0.0;
                LocalCoordsj[1] = 0.0;
                DNDej[0] = - 1.0;
                DNDej[1] = - 1.0;
            }
            else if( NodeIndex == 1 )
            {
                LocalCoordsj[0] = 1.0;
                LocalCoordsj[1] = 0.0;
                DNDej[0] = 1.0;
                DNDej[1] = 0.0;
            }
            else // NodeIndex == 2
            {
                LocalCoordsj[0] = 0.0;
                LocalCoordsj[1] = 1.0;
                DNDej[0] = 0.0;
                DNDej[1] = 1.0;
            }
        }
        else if( TNumNodes == 4 )    // linear quad element 
        {
            if( NodeIndex == 0 )
            {
                LocalCoordsj[0] = - 1.0;
                LocalCoordsj[1] = - 1.0;
                DNDej[0] = - 0.5;
                DNDej[1] = - 0.5;
            }
            else if( NodeIndex == 1 )
            {
                LocalCoordsj[0] =   1.0;
                LocalCoordsj[1] = - 1.0;
                DNDej[0] =   0.5;
                DNDej[1] = - 0.5;
            }
            else if( NodeIndex == 2 )
            {
                LocalCoordsj[0] =  1.0;
                LocalCoordsj[1] =  1.0;
                DNDej[0] =  0.5;
                DNDej[1] =  0.5;
            }
            else // NodeIndex == 3
            {
                LocalCoordsj[0] = - 1.0;
                LocalCoordsj[1] =   1.0;
                DNDej[0] = - 0.5;
                DNDej[1] =   0.5;
            }
        }
        
        CondGeometry.Jacobian( J, LocalCoordsj );
        
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
        delta_normal = Ce;
    }
    
    NeNorm = norm_2( Ne );
    const double NeNorm3 = NeNorm * NeNorm * NeNorm;
    
    if ( NeNorm3 > Tolerance )
    {
        const bounded_matrix<double, TDim, TDim> Cj = I / NeNorm - NeoNe / NeNorm3;
        delta_normal = prod( Cj, delta_normal );
    }
        
    return delta_normal; 
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional> // NOTE: Formulation taken from Mohamed Khalil work
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateDeltaNormalSlave(DerivativeDataType& rDerivativeData)
{
    for ( unsigned int i_slave = 0, i = 0; i_slave < TNumNodes; ++i_slave, i += TDim )
    {
        const bounded_matrix<double, TDim, TDim>& delta_normal = GetGeometry()[i_slave].GetValue(DELTA_NORMAL);
//         const bounded_matrix<double, TDim, TDim> delta_normal = this->LocalDeltaNormal(GetGeometry(), i_slave);
        for (unsigned i_dof = 0; i_dof < TDim; i_dof++) 
        {
            for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
            {
                row(rDerivativeData.DeltaNormalSlave[i_slave * TDim + i_dof], i_node) = trans(column(delta_normal, i_dof)); 
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional> // NOTE: Formulation taken from Mohamed Khalil work
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateDeltaNormalMaster(
    DerivativeDataType& rDerivativeData,
    GeometryType& MasterGeometry
    )
{
    for ( unsigned int i_master = 0, i = 0; i_master < TNumNodes; ++i_master, i += TDim )
    {
//        const bounded_matrix<double, TDim, TDim>& delta_normal = GetGeometry[i_master].GetValue(DELTA_NORMAL);
        const bounded_matrix<double, TDim, TDim> delta_normal = this->LocalDeltaNormal(MasterGeometry, i_master);
        for (unsigned i_dof = 0; i_dof < TDim; i_dof++) 
        {
            for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
            {
                row(rDerivativeData.DeltaNormalMaster[i_master * TDim + i_dof], i_node) = trans(column(delta_normal, i_dof)); 
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
    const VectorType N1 = rVariables.NSlave;
    const VectorType N2 = rVariables.NMaster;
    
    // Coordinates
    const bounded_matrix<double, TNumNodes, TDim> u1 = rDerivativeData.u1;
    const bounded_matrix<double, TNumNodes, TDim> X1 = rDerivativeData.X1;
    const bounded_matrix<double, TNumNodes, TDim> u2 = rDerivativeData.u2;
    const bounded_matrix<double, TNumNodes, TDim> X2 = rDerivativeData.X2;
    
    // Normals
    const array_1d<double, TDim > NormalSlaveg = prod(trans(rDerivativeData.NormalSlave), N1);
    const array_1d<double, TDim > NormalMasterg = prod(trans(rDerivativeData.NormalMaster), N2);
    
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
            array_1d<double, TNumNodes > AuxVector = ZeroVector(TNumNodes);
            AuxVector[i_dim] = 1.0;
            const array_1d<double, TNumNodes > Deltavector_nodes = - N1[i_slave] * AuxVector;
            const double DeltaDist = ((inner_prod(Deltavector_nodes, NormalMasterg))* inner_prod(NormalSlaveg, NormalMasterg) - inner_prod(vector_nodes, NormalMasterg) * (inner_prod(DNormalSlaveg, NormalMasterg)))/std::pow(inner_prod(NormalSlaveg, NormalMasterg) + tol, 2);
            const array_1d<double, TNumNodes > Deltax1g = N1[i_slave] * AuxVector;
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
            array_1d<double, TNumNodes > AuxVector = ZeroVector(TNumNodes);
//             if (i_master == 0) // NOTE: This is the way I considered in the symbolic
//             {
                AuxVector[i_dim] = 1.0;
//             }
            const array_1d<double, TNumNodes > Deltavector_nodes = N2[i_master] * AuxVector;
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
   DerivativeDataType& rDerivativeData
   )
{
    // Shape functions
    const VectorType N1 = rVariables.NSlave;
    
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
Matrix AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateDeltaPosition(
        Matrix& DeltaPosition,
        const ConditionArrayType LocalCoordinates
        )
{
    KRATOS_TRY;

    DeltaPosition = ZeroMatrix(TDim, TDim);

    for ( unsigned int i_node = 0; i_node < TNumNodes; i_node++ )
    {
        const array_1d<double, 3 > delta_displacement = GetGeometry()[i_node].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[i_node].FastGetSolutionStepValue(DISPLACEMENT,1);
        
        for ( unsigned int j_node = 0; j_node < TDim; j_node++ )
        {
            Vector N;
            GetGeometry().ShapeFunctionsValues( N, LocalCoordinates[j_node].Coordinates() );

            for ( unsigned int j_dim = 0; j_dim < TDim; j_dim++ )
            {
                DeltaPosition(j_node, j_dim) += N[i_node] * delta_displacement[j_dim];
            }
        }
    }

    return DeltaPosition;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
double AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::GetIntegrationWeight(
    GeneralVariables& rVariables,
    const GeometryType::IntegrationPointsArrayType& ThisIntegrationMethod,
    const unsigned int PointNumber
    )
{
    return ThisIntegrationMethod[PointNumber].Weight();
}
    
/***********************************************************************************/
/***********************************************************************************/

// Frictionless cases
template class AugmentedLagrangianMethodMortarContactCondition<2, 2, false>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 3, false>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 4, false>;

// Frictional cases
template class AugmentedLagrangianMethodMortarContactCondition<2, 2, true>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 3, true>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 4, true>;

} // Namespace Kratos
