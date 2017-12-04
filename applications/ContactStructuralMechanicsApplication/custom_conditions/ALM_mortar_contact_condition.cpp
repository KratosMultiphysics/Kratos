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
#ifdef KRATOS_DEBUG
#include <iomanip>
#endif

// External includes

// Project includes
/* Mortar includes */
#include "custom_conditions/ALM_mortar_contact_condition.h"

/* Additional includes */
#include <algorithm>

/* Utilities */
#include "utilities/math_utils.h"
#include "custom_utilities/search_utilities.h"

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
    ConditionMap::Pointer& all_conditions_maps = this->GetValue( MAPPING_PAIRS );
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
        ConditionMap::Pointer& all_conditions_maps = this->GetValue( MAPPING_PAIRS );
        
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
        ConditionMap::Pointer& all_conditions_maps = this->GetValue( MAPPING_PAIRS );
        GeometryType& this_geometry = GetGeometry();
        const double active_check_length = this_geometry.Length() * rCurrentProcessInfo[ACTIVE_CHECK_FACTOR];
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
    for ( unsigned int i = 0; i < rLeftHandSideMatrices.size( ); ++i )
    {
        // Note: rRightHandSideVectors.size() > 0
        this->InitializeSystemMatrices( rLeftHandSideMatrices[i], rRightHandSideVectors[0],local_system.CalculationFlags );
    }

    local_system.CalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_RHS_VECTOR, true );
    local_system.CalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX, false ); // Temporarily only
    for ( unsigned int i = 0; i < rRightHandSideVectors.size( ); ++i )
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
    for( unsigned int i = 0; i < rLeftHandSideMatrices.size( ); ++i )
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
    for( unsigned int i = 0; i < rRightHandSideVectors.size(); ++i )
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
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::AddExplicitContribution(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // The slave geometry
    GeometryType& slave_geometry = GetGeometry();
    const array_1d<double, 3>& normal_slave = this->GetValue(NORMAL);
    
    // Create and initialize condition variables
    GeneralVariables rVariables;
    
    // Create the current contact data
    DerivativeDataType rDerivativeData;
    rDerivativeData.Initialize(slave_geometry, rCurrentProcessInfo);
    
    // Create the mortar operators
    MortarConditionMatrices rThisMortarConditionMatrices;
    
    // We call the exact integration utility
    IntegrationUtility integration_utility = IntegrationUtility (mIntegrationOrder);
    
    // If we consider the normal variation
    const bool consider_normal_variation = rCurrentProcessInfo[CONSIDER_NORMAL_VARIATION];
    
    // Iterate over the master segments
    for (unsigned int pair_index = 0; pair_index < mPairSize; ++pair_index)
    {   
        if (mThisMasterElementsActive[pair_index] == true)
        {
            GeometryType& master_geometry = mThisMasterElements[pair_index]->GetGeometry();
            
            // The normal of the master condition
            const array_1d<double, 3>& master_normal = mThisMasterElements[pair_index]->GetValue(NORMAL);
            
            // Reading integration points
            ConditionArrayListType conditions_points_slave;
            const bool is_inside = integration_utility.GetExactIntegration(slave_geometry, normal_slave, master_geometry, master_normal, conditions_points_slave);
            
            double integration_area;
            integration_utility.GetTotalArea(slave_geometry, conditions_points_slave, integration_area);
            
            if ((is_inside == true) && ((integration_area/slave_geometry.Area()) > 1.0e-3 * slave_geometry.Area()))
            {            
                IntegrationMethod this_integration_method = GetIntegrationMethod();
                
                // Initialize general variables for the current master element
                rVariables.Initialize();
                
                // Update slave element info
                rDerivativeData.UpdateMasterPair(mThisMasterElements[pair_index]);
                
                // Initialize the mortar operators
                rThisMortarConditionMatrices.Initialize();
                
                const bool dual_LM = DerivativesUtilitiesType::CalculateAeAndDeltaAe(slave_geometry, normal_slave, mThisMasterElements[pair_index], rDerivativeData, rVariables, consider_normal_variation, conditions_points_slave, this_integration_method, GetAxisymmetricCoefficient(rVariables));
                
                for (unsigned int i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom)
                {
                    std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
                    array_1d<BelongType, TDim> belong_array;
                    for (unsigned int i_node = 0; i_node < TDim; ++i_node)
                    {
                        PointType global_point;
                        slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                        points_array[i_node] = boost::make_shared<PointType>(global_point);
                        belong_array[i_node] = conditions_points_slave[i_geom][i_node].GetBelong();
                    }
                    
                    DecompositionType decomp_geom( points_array );
                    
                    const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, slave_geometry.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);
                    
                    if (bad_shape == false)
                    {
                        const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( this_integration_method );
                        
                        // Integrating the mortar operators
                        for ( unsigned int point_number = 0; point_number < integration_points_slave.size(); ++point_number )
                        {
                            // We compute the local coordinates 
                            const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
                            PointType local_point_parent;
                            PointType gp_global;
                            decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                            slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);
                            
                            // Calculate the kinematic variables
                            this->CalculateKinematics( rVariables, rDerivativeData, master_normal, pair_index, local_point_decomp, local_point_parent, decomp_geom, dual_LM);//, delta_position_slave);
                            
                            const double integration_weight = integration_points_slave[point_number].Weight() * GetAxisymmetricCoefficient(rVariables);
                            
                            rThisMortarConditionMatrices.CalculateMortarOperators(rVariables, integration_weight);   
                        }
                    }
                }
                
                // Setting the weighted gap
                // Mortar condition matrices - DOperator and MOperator
                const bounded_matrix<double, TNumNodes, TNumNodes>& DOperator = rThisMortarConditionMatrices.DOperator;
                const bounded_matrix<double, TNumNodes, TNumNodes>& MOperator = rThisMortarConditionMatrices.MOperator;
                
                // Current coordinates 
                const bounded_matrix<double, TNumNodes, TDim>& x1 = MortarUtilities::GetCoordinates<TDim,TNumNodes>(slave_geometry);
                const bounded_matrix<double, TNumNodes, TDim>& x2 = MortarUtilities::GetCoordinates<TDim,TNumNodes>(master_geometry);
        
                const bounded_matrix<double, TNumNodes, TDim> D_x1_M_x2 = prod(DOperator, x1) - prod(MOperator, x2); 
                
                for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
                {
                    const array_1d<double, 3>& normal = slave_geometry[i_node].GetValue(NORMAL);
                    const array_1d<double, TDim> aux_array = row(D_x1_M_x2, i_node);
                                    
                    double& weighted_gap = slave_geometry[i_node].FastGetSolutionStepValue(WEIGHTED_GAP);
                    
                    #pragma omp atomic
                    weighted_gap += inner_prod(aux_array, - subrange(normal, 0, TDim)); 
                }
                
                if (TFrictional == true) // TODO: Check this!!!
                {
                    // Old coordinates 
                    const bounded_matrix<double, TNumNodes, TDim>& x1_old = MortarUtilities::GetCoordinates<TDim,TNumNodes>(slave_geometry, false, 1);
                    const bounded_matrix<double, TNumNodes, TDim>& x2_old = MortarUtilities::GetCoordinates<TDim,TNumNodes>(master_geometry, false, 1);
            
                    const bounded_matrix<double, TNumNodes, TDim> D_x1_old_M_x2_old = prod(DOperator, x1_old) - prod(MOperator, x2_old); 
                    
                    for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
                    {
                        // We compute the tangent
                        const array_1d<double, 3>& normal = slave_geometry[i_node].GetValue(NORMAL);
                        const array_1d<double, 3>& lm = slave_geometry[i_node].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
                        const double lm_normal = inner_prod(normal, lm);
                        array_1d<double, 3> tangent_lm = lm - lm_normal * normal;
                        tangent_lm /= norm_2(tangent_lm); 
                        const array_1d<double, TDim>& tangent = subrange(tangent_lm, 0, TDim);
                        
                        const array_1d<double, TDim>& aux_array = row(D_x1_old_M_x2_old, i_node);
                                  
                        double& weighted_slip = slave_geometry[i_node].FastGetSolutionStepValue(WEIGHTED_SLIP);

                        #pragma omp atomic
                        weighted_slip += inner_prod(aux_array, tangent); 
                    }
                }
            }
        }
    }

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
        
    // The slave geometry
    GeometryType& slave_geometry = this->GetGeometry();
    const array_1d<double, 3>& normal_slave = this->GetValue(NORMAL);
    
    // Create and initialize condition variables
    GeneralVariables rVariables;
    
    // Create the current contact data
    DerivativeDataType rDerivativeData;
    rDerivativeData.Initialize(slave_geometry, rCurrentProcessInfo);
    
    const bool consider_normal_variation = rCurrentProcessInfo[CONSIDER_NORMAL_VARIATION];
    
    // We compute the normal derivatives
    if (consider_normal_variation == true)
    {
        // Compute the normal derivatives of the slave
        DerivativesUtilitiesType::CalculateDeltaNormalSlave(rDerivativeData, GetGeometry());
    }
    
    // Create the mortar operators
    MortarConditionMatrices rThisMortarConditionMatrices;
    
    // We call the exact integration utility
    IntegrationUtility  integration_utility = IntegrationUtility (mIntegrationOrder);
    
    // Iterate over the master segments
    for (unsigned int pair_index = 0; pair_index < mPairSize; ++pair_index)
    {   
        GeometryType& master_geometry = mThisMasterElements[pair_index]->GetGeometry();
        
        if (mThisMasterElementsActive[pair_index] == true)
        {
            // The normal of the master condition
            const array_1d<double, 3>& master_normal = mThisMasterElements[pair_index]->GetValue(NORMAL);
            
            // Reading integration points
            ConditionArrayListType conditions_points_slave;
            const bool is_inside = integration_utility.GetExactIntegration(slave_geometry, normal_slave, master_geometry, master_normal, conditions_points_slave);
            
            double integration_area;
            integration_utility.GetTotalArea(slave_geometry, conditions_points_slave, integration_area);
            
            if ((is_inside == true) && ((integration_area/slave_geometry.Area()) > 1.0e-3 * slave_geometry.Area()))
            {            
                IntegrationMethod this_integration_method = GetIntegrationMethod();
                
                // Initialize general variables for the current master element
                rVariables.Initialize();
                
                // Update slave element info
                rDerivativeData.UpdateMasterPair(mThisMasterElements[pair_index]);
                
                // Initialize the mortar operators
                rThisMortarConditionMatrices.Initialize();
                
                if (consider_normal_variation == true)
                {
                    // Compute the normal derivatives of the master
                    DerivativesUtilitiesType::CalculateDeltaNormalMaster(rDerivativeData, master_geometry);
                }
                
                const bool dual_LM =  DerivativesUtilitiesType::CalculateAeAndDeltaAe(slave_geometry, normal_slave, mThisMasterElements[pair_index], rDerivativeData, rVariables, consider_normal_variation, conditions_points_slave, this_integration_method, GetAxisymmetricCoefficient(rVariables));
                
            #ifdef KRATOS_DEBUG
                if (dual_LM == false)
                {
                    std::cout << "WARNING:: NOT USING DUAL LM. Integration area: " << integration_area << "\tOriginal area: " << slave_geometry.Area() << "\tRatio: " << integration_area/slave_geometry.Area() << std::endl;
//                     IntegrationUtility::MathematicaDebug(this->Id(), slave_geometry, mThisMasterElements[pair_index]->Id(), master_geometry, conditions_points_slave);
                }
            #endif
                
                for (unsigned int i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom)
                {
                    std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
                    array_1d<BelongType, TDim> belong_array;
                    for (unsigned int i_node = 0; i_node < TDim; ++i_node)
                    {
                        PointType global_point;
                        slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                        points_array[i_node] = boost::make_shared<PointType>(global_point);
                        belong_array[i_node] = conditions_points_slave[i_geom][i_node].GetBelong();
                    }
                    
                    DecompositionType decomp_geom( points_array );
                    
                    const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, slave_geometry.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);
                    
                    if (bad_shape == false)
                    {
//                         // Delta position
//                         Matrix delta_position_slave;
//                         delta_position_slave = DerivativesUtilitiesType::CalculateDeltaPosition(delta_position_slave, slave_geometry, conditions_points_slave[i_geom]);
                        
                        const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( this_integration_method );
                        
                        // Integrating the mortar operators
                        for ( unsigned int point_number = 0; point_number < integration_points_slave.size(); ++point_number )
                        {
                            // We reset the derivatives
                            rDerivativeData.ResetDerivatives();
                            
                            // We compute the local coordinates 
                            const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
                            PointType local_point_parent;
                            PointType gp_global;
                            decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                            slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);
                            
                            // Calculate the kinematic variables
                            this->CalculateKinematics( rVariables, rDerivativeData, master_normal, pair_index, local_point_decomp, local_point_parent, decomp_geom, dual_LM);//, delta_position_slave);
                            
                            const double integration_weight = integration_points_slave[point_number].Weight() * GetAxisymmetricCoefficient(rVariables);
                            
                            if ( rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX ) ||
                                    rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
                            {
                                /* Update the derivatives */
                                // Update the derivative of the integration vertex (just in 3D)
                                if (TDim == 3) DerivativesUtilitiesType::CalculateDeltaCellVertex(rVariables, rDerivativeData, belong_array, consider_normal_variation, slave_geometry, master_geometry, normal_slave);
                                // Update the derivative of DetJ
                                DerivativesUtilitiesType::CalculateDeltaDetjSlave(decomp_geom, rVariables, rDerivativeData);
                                // Update the derivatives of the shape functions and the gap
                                DerivativesUtilitiesType::CalculateDeltaN(rVariables, rDerivativeData, slave_geometry, master_geometry, normal_slave, master_normal, decomp_geom, local_point_decomp, local_point_parent, consider_normal_variation, dual_LM);
                                
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
                const unsigned int active_inactive = GetActiveInactiveValue(slave_geometry);
                
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
    /* SHAPE FUNCTIONS */
    GetGeometry().ShapeFunctionsValues( rVariables.NSlave, LocalPointParent.Coordinates() );
    rVariables.PhiLagrangeMultipliers = (DualLM == true) ? prod(rDerivativeData.Ae, rVariables.NSlave) : rVariables.NSlave;
    
    /* SHAPE FUNCTION DERIVATIVES */
    GetGeometry().ShapeFunctionsLocalGradients( rVariables.DNDeSlave, LocalPointParent );
    
    /* CALCULATE JACOBIAN AND JACOBIAN DETERMINANT */
    rVariables.jSlave = GeometryDecomp.Jacobian( rVariables.jSlave, LocalPointDecomp.Coordinates());//, DeltaPosition);
//     rVariables.DetjSlave = MathUtils<double>::GeneralizedDet(rVariables.jSlave);
    rVariables.DetjSlave = GeometryDecomp.DeterminantOfJacobian( LocalPointDecomp );
    
    if (rVariables.DetjSlave < 0.0)
    {
        KRATOS_ERROR << "WARNING:: CONDITION ID: " << this->Id() << " INVERTED. DETJ: " << rVariables.DetjSlave << std::endl;
    }
    
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
    const array_1d<double,3> gp_normal = MortarUtilities::GaussPointUnitNormal(rVariables.NSlave, GetGeometry());
    
    GeometryType::CoordinatesArrayType slave_gp_global;
    this->GetGeometry( ).GlobalCoordinates( slave_gp_global, LocalPoint );
    MortarUtilities::FastProjectDirection( master_geometry, slave_gp_global, projected_gp_global, MasterNormal, -gp_normal ); // The opposite direction
    
    GeometryType::CoordinatesArrayType projected_gp_local;
    
    master_geometry.PointLocalCoordinates( projected_gp_local, projected_gp_global.Coordinates( ) ) ;

    // SHAPE FUNCTIONS 
    master_geometry.ShapeFunctionsValues(         rVariables.NMaster,    projected_gp_local );         
    master_geometry.ShapeFunctionsLocalGradients( rVariables.DNDeMaster, projected_gp_local );
    
    // JACOBIAN
//     Matrix delta_position_master; // MASTER
//     delta_position_master = DerivativesUtilitiesType::CalculateDeltaPosition(delta_position_master, GetGeometry() master_geometry);
    rVariables.jMaster = master_geometry.Jacobian( rVariables.jMaster, projected_gp_local);//, delta_position_master); // Add delta Position
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
    if ( rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
    {
        /* COMPONENT-WISE LHS MATRIX */
        const std::vector<Variable<MatrixType> >& rLeftHandSideVariables = rLocalSystem.GetLeftHandSideVariables( );
        bool calculated;

        for ( unsigned int i = 0; i < rLeftHandSideVariables.size( ); ++i )
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
        const unsigned int rActiveInactive
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
        const unsigned int rActiveInactive
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
        const unsigned int rActiveInactive
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
        const unsigned int rActiveInactive
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
        const unsigned int rActiveInactive
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
        const unsigned int rActiveInactive
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

        for ( unsigned int i = 0; i < rRightHandSideVariables.size( ); ++i )
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
        const unsigned int rActiveInactive
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
        const unsigned int rActiveInactive
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
        const unsigned int rActiveInactive
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
        const unsigned int rActiveInactive
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
        const unsigned int rActiveInactive
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
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, you are evil, and your seed must be eradicated from the face of the earth" << std::endl;
    
    return ZeroVector(36);
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

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
double AugmentedLagrangianMethodMortarContactCondition< TDim, TNumNodes, TFrictional >::GetAxisymmetricCoefficient(const GeneralVariables& rVariables) const
{
    return 1.0;
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
