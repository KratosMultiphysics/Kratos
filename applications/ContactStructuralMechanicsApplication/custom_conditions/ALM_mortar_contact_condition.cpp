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
#include "utilities/exact_mortar_segmentation_utility.h"

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
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::AddExplicitContribution(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // The slave geometry
    GeometryType& slave_geometry = GetGeometry();
    
    // Create and initialize condition variables
    GeneralVariables rVariables;
    
    // Create the current contact data
    DerivativeDataType rDerivativeData;
    rDerivativeData.Initialize(slave_geometry, rCurrentProcessInfo);
    
    // Create the mortar operators
    MortarConditionMatrices rThisMortarConditionMatrices;
    
    // We call the exact integration utility
    ExactMortarIntegrationUtility<TDim, TNumNodes, true>  integration_utility = ExactMortarIntegrationUtility<TDim, TNumNodes, true> (mIntegrationOrder);
    
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
            const bool is_inside = integration_utility.GetExactIntegration(slave_geometry, this->GetValue(NORMAL), master_geometry, master_normal, conditions_points_slave);
            
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
                
                const bool dual_LM = this->CalculateAeAndDeltaAe(rDerivativeData, rVariables, rCurrentProcessInfo, pair_index, conditions_points_slave, this_integration_method, master_normal);
                
                for (unsigned int i_geom = 0; i_geom < conditions_points_slave.size(); i_geom++)
                {
                    std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
                    array_1d<BelongType, TDim> belong_array;
                    for (unsigned int i_node = 0; i_node < TDim; i_node++)
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
                        for ( unsigned int point_number = 0; point_number < integration_points_slave.size(); point_number++ )
                        {
                            const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
                            PointType local_point_parent;
                            PointType gp_global;
                            decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                            slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);
                            
                            // Calculate the kinematic variables
                            this->CalculateKinematics( rVariables, rDerivativeData, master_normal, pair_index, local_point_decomp, local_point_parent, decomp_geom, dual_LM);//, delta_position_slave);
                            
                            const double integration_weight = GetIntegrationWeight(rVariables, integration_points_slave, point_number);
                            
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
                
                for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
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
                    
                    for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
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
    GeometryType& slave_geometry = GetGeometry();
    
    // Create and initialize condition variables
    GeneralVariables rVariables;
    
    // Create the current contact data
    DerivativeDataType rDerivativeData;
    rDerivativeData.Initialize(slave_geometry, rCurrentProcessInfo);
    
    const bool& consider_normal_variation = rCurrentProcessInfo[CONSIDER_NORMAL_VARIATION];
    
    // We compute the normal derivatives
    if (consider_normal_variation == true)
    {
        // Compute the normal derivatives of the master
        this->CalculateDeltaNormalSlave(rDerivativeData);
    }
    
    // Create the mortar operators
    MortarConditionMatrices rThisMortarConditionMatrices;
    
    // We call the exact integration utility
    ExactMortarIntegrationUtility<TDim, TNumNodes, true>  integration_utility = ExactMortarIntegrationUtility<TDim, TNumNodes, true> (mIntegrationOrder);
    
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
            const bool is_inside = integration_utility.GetExactIntegration(slave_geometry, this->GetValue(NORMAL), master_geometry, master_normal, conditions_points_slave);
            
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
                    this->CalculateDeltaNormalMaster(rDerivativeData, master_geometry);
                }
                
                const bool dual_LM = this->CalculateAeAndDeltaAe(rDerivativeData, rVariables, rCurrentProcessInfo, pair_index, conditions_points_slave, this_integration_method, master_normal);
                
            #ifdef KRATOS_DEBUG
                if (dual_LM == false)
                {
                    std::cout << "WARNING:: NOT USING DUAL LM. Integration area: " << integration_area << "\tOriginal area: " << slave_geometry.Area() << "\tRatio: " << integration_area/slave_geometry.Area() << std::endl;
//                     std::cout << "Slave Condition ID: " << this->Id() << std::endl;
//                     for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
//                     {
//                         std::cout << "NODE ID: " << slave_geometry[i_node].Id() << "\tX: " << slave_geometry[i_node].X() << "\tY: " << slave_geometry[i_node].Y() << "\tZ: " << slave_geometry[i_node].Z() << std::endl;
//                     }
//                     std::cout << "Master Condition ID: " << mThisMasterElements[pair_index]->Id() << std::endl;
//                     for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
//                     {
//                         std::cout << "NODE ID: " << master_geometry[i_node].Id() << "\tX: " << master_geometry[i_node].X() << "\tY: " << master_geometry[i_node].Y() << "\tZ: " << master_geometry[i_node].Z() << std::endl;
//                     }
//                     std::cout << std::endl;
                    
// //                     // Mathematica debug
// //                     auto& slave_geometry = GetGeometry();
// //                     std::cout << "\nGraphics3D[{EdgeForm[{Thick,Dashed,Red}],FaceForm[],Polygon[{{";
// //             
// //                     for (unsigned int i = 0; i < TNumNodes; i++)
// //                     {
// //                         std::cout << slave_geometry[i].X() << "," << slave_geometry[i].Y() << "," << slave_geometry[i].Z();
// //                         
// //                         if (i < TNumNodes - 1) std::cout << "},{";
// //                     }
// //                     std::cout << "}}],Text[Style["<< this->Id() <<", Tiny],{"<< slave_geometry.Center().X() << "," << slave_geometry.Center().Y() << ","<< slave_geometry.Center().Z() << "}]}],";// << std::endl;
// //                     
// //                     std::cout << "\nGraphics3D[{EdgeForm[{Thick,Dashed,Blue}],FaceForm[],Polygon[{{";
// //             
// //                     for (unsigned int i = 0; i < TNumNodes; i++)
// //                     {
// //                         std::cout << master_geometry[i].X() << "," << master_geometry[i].Y() << "," << master_geometry[i].Z();
// //                         
// //                         if (i < TNumNodes - 1) std::cout << "},{";
// //                     }
// //                     std::cout << "}}],Text[Style["<< mThisMasterElements[pair_index]->Id() <<", Tiny],{"<< master_geometry.Center().X() << "," << master_geometry.Center().Y() << ","<< master_geometry.Center().Z() << "}]}],";// << std::endl;
// //                     
// //                     for (unsigned int i_geom = 0; i_geom < conditions_points_slave.size(); i_geom++)
// //                     {
// //                         std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
// //                         for (unsigned int i_node = 0; i_node < TDim; i_node++)
// //                         {
// //                             PointType global_point;
// //                             slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
// //                             points_array[i_node] = boost::make_shared<PointType>(global_point);
// //                         }
// //                         
// //                         DecompositionType decomp_geom( points_array );
// //                         
// //                         std::cout << "\nGraphics3D[{Opacity[.3],Triangle[{{"; 
// //                         for (unsigned int i = 0; i < 3; i++)
// //                         {
// //                             std::cout << std::setprecision(16) << decomp_geom[i].X() << "," << decomp_geom[i].Y() << "," << decomp_geom[i].Z();
// //                             
// //                             if (i < 2) std::cout << "},{";
// //                         }
// //                         std::cout << "}}]}],";// << std::endl;
// //                     }
// //                     
// //                     std::cout << std::endl;
                }
            #endif
                
                for (unsigned int i_geom = 0; i_geom < conditions_points_slave.size(); i_geom++)
                {
                    std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
                    array_1d<BelongType, TDim> belong_array;
                    for (unsigned int i_node = 0; i_node < TDim; i_node++)
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
//                         delta_position_slave = CalculateDeltaPosition(delta_position_slave, conditions_points_slave[i_geom]);
                        
                        const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( this_integration_method );
                        
                        // Integrating the mortar operators
                        for ( unsigned int point_number = 0; point_number < integration_points_slave.size(); point_number++ )
                        {
                            const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
                            PointType local_point_parent;
                            PointType gp_global;
                            decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                            slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);
                            
                            // Calculate the kinematic variables
                            this->CalculateKinematics( rVariables, rDerivativeData, master_normal, pair_index, local_point_decomp, local_point_parent, decomp_geom, dual_LM);//, delta_position_slave);
                            
                            const double integration_weight = GetIntegrationWeight(rVariables, integration_points_slave, point_number);
                            
                            if ( rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX ) ||
                                    rLocalSystem.CalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
                            {
                                /* Update the derivatives */
                                // Update the derivative of the integration vertex (just in 3D)
//                                 if (TDim == 3) this->CalculateDeltaCellVertex(rVariables, rDerivativeData, belong_array, consider_normal_variation, master_geometry);
                                // Update the derivative of DetJ
                                this->CalculateDeltaDetjSlave(rVariables, rDerivativeData);
                                // Update the derivatives of the shape functions and the gap
//                                 this->CalculateDeltaN(rVariables, rDerivativeData, master_geometry, master_normal, decomp_geom, local_point_decomp, local_point_parent);
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
//                 #ifdef KRATOS_DEBUG
//                     else
//                     {
//                         std::cout << "WARNING:: BAD SHAPE GEOMETRY" << std::endl;
//                         std::cout << "Slave Condition ID: " << this->Id() << std::endl;
//                         std::cout << "Master Condition ID: " << mThisMasterElements[pair_index]->Id() << std::endl;
//                         for (unsigned int i_node = 0; i_node < TDim; ++i_node)
//                         {
//                             std::cout << "X: " << decomp_geom[i_node].X() << "\tY: " << decomp_geom[i_node].Y() << "\tZ: " << decomp_geom[i_node].Z() << std::endl;
//                         }
//                         std::cout << std::endl;
//                         
// //                         // Mathematica debug
// //                         std::cout << "\nGraphics3D[{Opacity[.3],Triangle[{{"; 
// //                         for (unsigned int i = 0; i < 3; i++)
// //                         {
// //                             std::cout << std::setprecision(16) << decomp_geom[i_node].X() << "," << decomp_geom[i_node].Y() << "," << decomp_geom[i_node].Z();
// //                             
// //                             if (i < 2) std::cout << "},{";
// //                         }
// //                         std::cout << "}}]}],";// << std::endl;
//                     }
//                 #endif
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
    
    const bool& consider_normal_variation = rCurrentProcessInfo[CONSIDER_NORMAL_VARIATION];
    
    for (unsigned int i_geom = 0; i_geom < ConditionsPointsSlave.size(); i_geom++)
    {
        std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
        array_1d<BelongType, TDim> belong_array;
        for (unsigned int i_node = 0; i_node < TDim; i_node++)
        {
            PointType global_point;
            GetGeometry().GlobalCoordinates(global_point, ConditionsPointsSlave[i_geom][i_node]);
            points_array[i_node] = boost::make_shared<PointType>(global_point);
            belong_array[i_node] = ConditionsPointsSlave[i_geom][i_node].GetBelong();
        }
        
        DecompositionType decomp_geom( points_array );
        
        const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, this->GetGeometry().Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);
        
        if (bad_shape == false)
        {
//             /* Delta position */
//             Matrix delta_position_slave;
//             delta_position_slave = CalculateDeltaPosition(delta_position_slave, ConditionsPointsSlave[i_geom]);
            
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
                this->CalculateKinematics( rVariables, rDerivativeData, MasterNormal, PairIndex, local_point_decomp, local_point_parent, decomp_geom, false);//, delta_position_slave);
                
                // Update the derivative of the integration vertex (just in 3D)
                if (TDim == 3) this->CalculateDeltaCellVertex(rVariables, rDerivativeData, belong_array, consider_normal_variation, mThisMasterElements[PairIndex]->GetGeometry());
                                
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
    
    /* CALCULATE JACOBIAN AND JACOBIAN DETERMINANT */
    rVariables.jSlave = GeometryDecomp.Jacobian( rVariables.jSlave, LocalPointDecomp.Coordinates());//, DeltaPosition);
//     rVariables.DetjSlave = MathUtils<double>::GeneralizedDet(rVariables.jSlave);
    rVariables.DetjSlave = GeometryDecomp.DeterminantOfJacobian( LocalPointDecomp );
    
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
    Matrix delta_position_master; // MASTER
    delta_position_master = CalculateDeltaPosition(delta_position_master, master_geometry);
    rVariables.jMaster = master_geometry.Jacobian( rVariables.jMaster, projected_gp_local, delta_position_master); // Add delta Position
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
   const array_1d<BelongType, TDim>& TheseBelongs,
   const bool& ConsiderNormalVariation,
   GeometryType& MasterGeometry
   ) 
{
    // The normal and delta normal in the center of the element
    const array_1d<double, 3>& normal = this->GetValue(NORMAL);
    bounded_matrix<double, TDim, TDim> delta_normal = ZeroMatrix(TDim, TDim);
    
    const VectorType& N1 = rVariables.NSlave;
    const VectorType& N2 = rVariables.NMaster;
    
    if (ConsiderNormalVariation == true)
    {
        for ( unsigned int i_slave = 0, i = 0; i_slave < TNumNodes; ++i_slave, i += TDim )
        {
            delta_normal += this->LocalDeltaNormal(GetGeometry(), i_slave); // TODO: Check this!!!!
        }
        
        delta_normal /= static_cast<double>(TNumNodes);
    }
    
    for (unsigned i_belong = 0; i_belong < 3; i_belong++) 
    {
        if (TheseBelongs[i_belong] >= 2 * TNumNodes) // It belongs to an intersection
        {    
            // We compute the indexes
            unsigned int belong_index_slave_start, belong_index_slave_end, belong_index_master_start, belong_index_master_end;
            ConvertAuxHashIndex(static_cast<unsigned int>(TheseBelongs[i_belong]), belong_index_slave_start, belong_index_slave_end, belong_index_master_start, belong_index_master_end);
            
            // We add the initial value // NOTE: Added in the second part
//             for (unsigned i_dof = 0; i_dof < TDim; i_dof++) 
//             {
//                 LocalDeltaVertex(row(rDerivativeData.DeltaCellVertex[belong_index_slave_start * TDim + i_dof], i_belong), normal, delta_normal, N1, N2, i_dof, belong_index_slave_start, ConsiderNormalVariation, MasterGeometry);
//             }
            
            const array_1d<double, 3> xs1 = GetGeometry()[belong_index_slave_start].Coordinates();
            const array_1d<double, 3> xe1 = GetGeometry()[belong_index_slave_end].Coordinates();
            const array_1d<double, 3> xs2 = MasterGeometry[belong_index_master_start].Coordinates();
            const array_1d<double, 3> xe2 = MasterGeometry[belong_index_master_end].Coordinates();
            
            array_1d<double, 3> aux_num, aux_denom;
            MathUtils<double>::CrossProduct(aux_num, xs1 - xs2, xe2 - xs2);
            MathUtils<double>::CrossProduct(aux_denom, xe1 - xs1, xe2 - xs2);
            const double num = inner_prod(aux_num, normal);
            const double denom = inner_prod(aux_denom, normal);
            
            // We compute the first part
            array_1d<double, 3> aux_coords = (xe1 - xs1);
            array_1d<double, 3> aux_vertex_matrix;
            array_1d<double, 3> aux_cross_product;
            for (unsigned i_dof = 0; i_dof < TDim; i_dof++)
            {
                aux_vertex_matrix = ZeroVector(3);
                
                LocalDeltaVertex(aux_vertex_matrix, normal, delta_normal, N1, N2, i_dof, belong_index_slave_start, ConsiderNormalVariation, MasterGeometry, 1.0);
                LocalDeltaVertex(aux_vertex_matrix, normal, delta_normal, N1, N2, i_dof, belong_index_slave_end, ConsiderNormalVariation, MasterGeometry, - 1.0);   
            }
            
            MathUtils<double>::CrossProduct(aux_cross_product, aux_vertex_matrix, xe2 - xs2);
            const double coeff1 = inner_prod(aux_cross_product, normal);
            
            for (unsigned i_dof = 0; i_dof < TDim; i_dof++)
            {
                aux_vertex_matrix = ZeroVector(3);
                
                LocalDeltaVertex(aux_vertex_matrix, normal, delta_normal, N1, N2, i_dof, belong_index_master_end, ConsiderNormalVariation, MasterGeometry, 1.0);
                LocalDeltaVertex(aux_vertex_matrix, normal, delta_normal, N1, N2, i_dof, belong_index_master_start, ConsiderNormalVariation, MasterGeometry, - 1.0);   
            }

            MathUtils<double>::CrossProduct(aux_cross_product, xs1 - xs2, aux_vertex_matrix);
            const double coeff2 = inner_prod(aux_cross_product, normal);
            
            for (unsigned i_dof = 0; i_dof < TDim; i_dof++)
            {
                aux_vertex_matrix = ZeroVector(3);
                
                LocalDeltaVertex(aux_vertex_matrix, normal, delta_normal, N1, N2, i_dof, belong_index_slave_end, ConsiderNormalVariation, MasterGeometry, - 1.0);
                LocalDeltaVertex(aux_vertex_matrix, normal, delta_normal, N1, N2, i_dof, belong_index_slave_start, ConsiderNormalVariation, MasterGeometry, 1.0);   
            }

            MathUtils<double>::CrossProduct(aux_cross_product, aux_vertex_matrix, xe2 - xs2);
            const double coeff3 = inner_prod(aux_cross_product, normal);
            
            for (unsigned i_dof = 0; i_dof < TDim; i_dof++)
            {
                aux_vertex_matrix = ZeroVector(3);
                
                LocalDeltaVertex(aux_vertex_matrix, normal, delta_normal, N1, N2, i_dof, belong_index_master_end, ConsiderNormalVariation, MasterGeometry, - 1.0);
                LocalDeltaVertex(aux_vertex_matrix, normal, delta_normal, N1, N2, i_dof, belong_index_master_start, ConsiderNormalVariation, MasterGeometry, 1.0);   
            }

            MathUtils<double>::CrossProduct(aux_cross_product, xe1 - xs1, aux_vertex_matrix);
            const double coeff4 = inner_prod(aux_cross_product, normal);
            
            for ( unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave)
            {
                for (unsigned i_dof = 0; i_dof < TDim; i_dof++)
                {
                    row(rDerivativeData.DeltaCellVertex[i_slave * TDim + i_dof], i_belong) -= coeff1/denom * aux_coords; 
                    row(rDerivativeData.DeltaCellVertex[i_slave * TDim + i_dof], i_belong) -= coeff2/denom * aux_coords; 
                    row(rDerivativeData.DeltaCellVertex[i_slave * TDim + i_dof], i_belong) -= num * coeff3/std::pow(denom, 2.0) * aux_coords; 
                    row(rDerivativeData.DeltaCellVertex[i_slave * TDim + i_dof], i_belong) -= num * coeff4/std::pow(denom, 2.0) * aux_coords; 
                }
            }
            
            for ( unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave)
            {
                for (unsigned i_dof = 0; i_dof < TDim; i_dof++)
                {
                    row(rDerivativeData.DeltaCellVertex[i_slave * TDim + i_dof], i_belong) += inner_prod(aux_num,   trans(column(delta_normal, i_dof)))/std::pow(denom, 2.0) * aux_coords; 
                    row(rDerivativeData.DeltaCellVertex[i_slave * TDim + i_dof], i_belong) -= inner_prod(aux_denom, trans(column(delta_normal, i_dof)))/std::pow(denom, 2.0) * aux_coords; 
                }
            }
            
            // We compute the second part
            double coeff0 = num/denom;
            for (unsigned i_dof = 0; i_dof < TDim; i_dof++)
            {
                LocalDeltaVertex(row(rDerivativeData.DeltaCellVertex[belong_index_slave_end * TDim + i_dof], i_belong), normal, delta_normal, N1, N2, i_dof, belong_index_slave_end, ConsiderNormalVariation, MasterGeometry,     - coeff0);
                LocalDeltaVertex(row(rDerivativeData.DeltaCellVertex[belong_index_slave_start * TDim + i_dof], i_belong), normal, delta_normal, N1, N2, i_dof, belong_index_slave_start, ConsiderNormalVariation, MasterGeometry, 1.0 + coeff0);
            }
        }
        else // It belongs to a master/slave node
        {
            const unsigned int belong_index = static_cast<unsigned int>(TheseBelongs[i_belong]);
            
            for (unsigned i_dof = 0; i_dof < TDim; i_dof++)
            {
                LocalDeltaVertex(row(rDerivativeData.DeltaCellVertex[belong_index * TDim + i_dof], i_belong), normal, delta_normal, N1, N2, i_dof, belong_index, ConsiderNormalVariation, MasterGeometry);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::ConvertAuxHashIndex(
    const unsigned int& AuxIndex,
    unsigned int& BelongIndexSlaveStart, 
    unsigned int& BelongIndexSlaveEnd, 
    unsigned int& BelongIndexMasterStart, 
    unsigned int& BelongIndexMasterEnd
    )
{
    unsigned int index_to_decompose = AuxIndex - 2 * TNumNodes;
    
    BelongIndexMasterEnd = index_to_decompose/10000;
    index_to_decompose = std::fmod(index_to_decompose, 10000);
    BelongIndexMasterStart = index_to_decompose/1000;
    index_to_decompose = std::fmod(index_to_decompose, 1000);
    BelongIndexSlaveEnd = index_to_decompose/100;
    index_to_decompose = std::fmod(index_to_decompose, 100);
    BelongIndexSlaveStart = index_to_decompose/10;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::LocalDeltaVertex(
   array_1d<double, 3> DeltaVertexMatrix,
   const array_1d<double, 3>& Normal,
   const bounded_matrix<double, TDim, TDim>& DeltaNormal,
   const VectorType& N1,
   const VectorType& N2,
   const unsigned & iDoF,
   const unsigned & BelongIndex,
   const bool& ConsiderNormalVariation,
   GeometryType& MasterGeometry,
   const double Coeff 
   ) 
{
    const array_1d<double, 3> coords_center = GetGeometry().Center().Coordinates();
    array_1d<double, 3> coords_node;
    if (BelongIndex < TNumNodes)
    {
        coords_node = GetGeometry()[BelongIndex].Coordinates();
    }
    else
    {
        coords_node = MasterGeometry[BelongIndex - TNumNodes].Coordinates();
    }
    
    // The corresponding part to the nodal coordinates
    if (BelongIndex < TNumNodes) // SLAVE
    {
        array_1d<double, 3> aux_der(0.0);
        aux_der[iDoF] = 1.0;
        DeltaVertexMatrix += Coeff * aux_der * N1[BelongIndex]; 
    }
    
    // The corresponding part to the normal
    const double coordsxdeltanormal = inner_prod(coords_node - coords_center, column(DeltaNormal, iDoF));
    
    if (BelongIndex < TNumNodes) // SLAVE
    {
        const double factor_belong = - N1[BelongIndex]/static_cast<double>(TNumNodes) + N1[BelongIndex];   
        const double deltacoordsxnormal =  factor_belong * Normal[iDoF];
        
        DeltaVertexMatrix += - Coeff * Normal * (deltacoordsxnormal + coordsxdeltanormal);
    }
    else // MASTER
    { 
        const double deltacoordsxnormal = N2[BelongIndex - TNumNodes] * Normal[iDoF];
        
        DeltaVertexMatrix += - Coeff * Normal * deltacoordsxnormal;
    }
    
    // The corresponding part to delta normal
    const double coordsxnormal = - inner_prod(coords_node - coords_center, Normal);
    DeltaVertexMatrix += Coeff * coordsxnormal * trans(column(DeltaNormal, iDoF));
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
        
        const array_1d<double,TDim>& normal = prod(trans(rDerivativeData.NormalMaster), rVariables.NSlave);
        
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
            
            rDerivativeData.DeltaDetjSlave[i    ] = inner_prod( normal, column( DeltaJxixJeta, 0 ) );
            rDerivativeData.DeltaDetjSlave[i + 1] = inner_prod( normal, column( DeltaJxixJeta, 1 ) );
            rDerivativeData.DeltaDetjSlave[i + 2] = inner_prod( normal, column( DeltaJxixJeta, 2 ) );
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
    const double tolerance = 1.0e-14;
        
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
    
    if ( NeNorm3 > tolerance )
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
    for ( unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave)
    {
        const bounded_matrix<double, TDim, TDim>& delta_normal = GetGeometry()[i_slave].GetValue(DELTA_NORMAL);
//         const bounded_matrix<double, TDim, TDim> delta_normal = this->LocalDeltaNormal(GetGeometry(), i_slave);
        for (unsigned i_dof = 0; i_dof < TDim; i_dof++) 
        {
            row(rDerivativeData.DeltaNormalSlave[i_slave * TDim + i_dof], i_slave) = trans(column(delta_normal, i_dof)); 
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<> // NOTE: Formulation taken from Mohamed Khalil work
void AugmentedLagrangianMethodMortarContactCondition<2,2,false>::CalculateDeltaNormalMaster(
    DerivativeDataType& rDerivativeData,
    GeometryType& MasterGeometry
    )
{
    for ( unsigned int i_master = 0; i_master < 2; ++i_master )
    {
//        const bounded_matrix<double, 2, 2>& delta_normal = GetGeometry[i_master].GetValue(DELTA_NORMAL);
        const bounded_matrix<double, 2, 2> delta_normal = this->LocalDeltaNormal(MasterGeometry, i_master);
        for (unsigned i_dof = 0; i_dof < 2; i_dof++) 
        {
            row(rDerivativeData.DeltaNormalMaster[i_master * 2 + i_dof], i_master) = trans(column(delta_normal, i_dof));
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AugmentedLagrangianMethodMortarContactCondition<2,2,true>::CalculateDeltaNormalMaster(
    DerivativeDataType& rDerivativeData,
    GeometryType& MasterGeometry
    )
{
    for ( unsigned int i_master = 0; i_master < 2; ++i_master )
    {
//        const bounded_matrix<double, 2, 2>& delta_normal = GetGeometry[i_master].GetValue(DELTA_NORMAL);
        const bounded_matrix<double, 2, 2> delta_normal = this->LocalDeltaNormal(MasterGeometry, i_master);
        for (unsigned i_dof = 0; i_dof < 2; i_dof++) 
        {
            row(rDerivativeData.DeltaNormalMaster[i_master * 2 + i_dof], i_master) = trans(column(delta_normal, i_dof));
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AugmentedLagrangianMethodMortarContactCondition<3,3,false>::CalculateDeltaNormalMaster(
    DerivativeDataType& rDerivativeData,
    GeometryType& MasterGeometry
    ){}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AugmentedLagrangianMethodMortarContactCondition<3,3,true>::CalculateDeltaNormalMaster(
    DerivativeDataType& rDerivativeData,
    GeometryType& MasterGeometry
    ){}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AugmentedLagrangianMethodMortarContactCondition<3,4,false>::CalculateDeltaNormalMaster(
    DerivativeDataType& rDerivativeData,
    GeometryType& MasterGeometry
    ){}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AugmentedLagrangianMethodMortarContactCondition<3,4,true>::CalculateDeltaNormalMaster(
    DerivativeDataType& rDerivativeData,
    GeometryType& MasterGeometry
    ){}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateDeltaN(
   GeneralVariables& rVariables,
   DerivativeDataType& rDerivativeData,
   GeometryType& MasterGeometry,
   const array_1d<double, 3> MasterNormal,
   const DecompositionType& DecompGeom,
   const PointType& LocalPointDecomp,
   const PointType& LocalPointParent
   )
{
    /* Shape functions */
    const VectorType& N1 = rVariables.NSlave;
    const VectorType& N2 = rVariables.NMaster;
    /* Local gradients */
    const MatrixType& DNDe1 = rVariables.DNDeSlave;
    const MatrixType& DNDe2 = rVariables.DNDeMaster;
    /* Jacobians */
    const MatrixType& LHS1 = rVariables.jSlave;
    const MatrixType& LHS2 = rVariables.jMaster; 

    /* Shape function decomposition */
    VectorType N_decomp;
    DecompGeom.ShapeFunctionsValues( N_decomp, LocalPointDecomp.Coordinates() );
    
    /* LHSs */
    MatrixType inv_LHS1, inv_LHS2;
    double aux_det1, aux_det2;
    MathUtils<double>::GeneralizedInvertMatrix(LHS1, inv_LHS1, aux_det1);
    MathUtils<double>::GeneralizedInvertMatrix(LHS2, inv_LHS2, aux_det2);
    
    if (TDim == 2)
    {
        // TODO: Finish this!!!!
    }
    else
    {
        for ( unsigned int i_node = 0; i_node < 2 * TNumNodes; ++i_node)
        {
            for (unsigned i_dof = 0; i_dof < TDim; i_dof++) 
            {
                VectorType aux_RHS1 = ZeroVector(3);
                VectorType aux_RHS2 = ZeroVector(3);
            
                // Local contribution
                if (i_node < TNumNodes)
                {
                    aux_RHS1[i_dof] -= N1[i_node];
                }
                else
                {
                    aux_RHS2[i_dof] -= N2[i_node - TNumNodes];
                }
            
                // The vertex cell contribution
                if (i_node < TNumNodes)
                {
                    for(unsigned int i_belong = 0; i_belong < 3; i_belong++)
                    {
                        aux_RHS1 += N1[i_node] * N_decomp[i_belong] * row(rDerivativeData.DeltaCellVertex[i_node * TDim + i_dof], i_belong);
                    }
                }
                else
                {           
                    for(unsigned int i_belong = 0; i_belong < 3; i_belong++)
                    {
                        aux_RHS2 += N2[i_node - TNumNodes] * N_decomp[i_belong] * row(rDerivativeData.DeltaCellVertex[i_node * TDim + i_dof], i_belong);
                    }
                }
                
                // We compute the delta coordinates 
                const VectorType aux_delta_coords1 = prod(inv_LHS1, aux_RHS1);
                const VectorType aux_delta_coords2 = prod(inv_LHS2, aux_RHS2);
                
                // Now we can compute the delta shape functions // FIXME: Not improving converence (check this)
                double delta_position = 1.0; // TODO: Check the consistency of this!!!! (consider the delta or not)
                CalculateDeltaPosition(delta_position, MasterGeometry, i_node, i_dof);
                
                const double tolerance = std::numeric_limits<double>::epsilon();
                
                if (std::abs(aux_det1) > tolerance) rDerivativeData.DeltaN1[i_node * TDim + i_dof] = delta_position * (aux_delta_coords1[0] * column(DNDe1, 0) + aux_delta_coords1[1] * column(DNDe1, 1));
                if (std::abs(aux_det2) > tolerance) rDerivativeData.DeltaN2[i_node * TDim + i_dof] = delta_position * (aux_delta_coords2[0] * column(DNDe2, 0) + aux_delta_coords2[1] * column(DNDe2, 1));
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
    const VectorType& N1 = rVariables.NSlave;
    
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
Matrix& AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateDeltaPosition(
        Matrix& DeltaPosition,
        const ConditionArrayType& LocalCoordinates
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
Matrix& AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateDeltaPosition(
        Matrix& DeltaPosition,
        GeometryType& ThisGeometry
        )
{
    KRATOS_TRY;

    DeltaPosition = ZeroMatrix(TNumNodes, TDim);

    for ( unsigned int i_node = 0; i_node < TNumNodes; i_node++ )
    {
        const array_1d<double, 3 > delta_displacement = ThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT) - ThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT,1);
        
        for ( unsigned int i_dim = 0; i_dim < TDim; i_dim++ )
        {
            DeltaPosition(i_node, i_dim) += delta_displacement[i_dim];
        }
    }

    return DeltaPosition;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateDeltaPosition(
        VectorType& DeltaPosition,
        GeometryType& MasterGeometry,
        const unsigned int& IndexNode
        )
{
    KRATOS_TRY;

    if (IndexNode < TNumNodes)
    {
        DeltaPosition = GetGeometry()[IndexNode].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[IndexNode].FastGetSolutionStepValue(DISPLACEMENT,1);
    }
    else
    {
        DeltaPosition = MasterGeometry[IndexNode - TNumNodes].FastGetSolutionStepValue(DISPLACEMENT) - MasterGeometry[IndexNode - TNumNodes].FastGetSolutionStepValue(DISPLACEMENT,1);
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateDeltaPosition(
        VectorType& DeltaPosition,
        GeometryType& MasterGeometry,
        const unsigned int& IndexNode,
        const unsigned int& iDoF
        )
{
    KRATOS_TRY;

    DeltaPosition = ZeroVector(3);
    
    if (IndexNode < TNumNodes)
    {
        DeltaPosition[iDoF] = (GetGeometry()[IndexNode].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[IndexNode].FastGetSolutionStepValue(DISPLACEMENT,1))[iDoF];
    }
    else
    {
        DeltaPosition[iDoF] = (MasterGeometry[IndexNode - TNumNodes].FastGetSolutionStepValue(DISPLACEMENT) - MasterGeometry[IndexNode - TNumNodes].FastGetSolutionStepValue(DISPLACEMENT,1))[iDoF];
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::CalculateDeltaPosition(
        double& DeltaPosition,
        GeometryType& MasterGeometry,
        const unsigned int& IndexNode,
        const unsigned int& iDoF
        )
{
    KRATOS_TRY;
    
    if (IndexNode < TNumNodes)
    {
        DeltaPosition = (GetGeometry()[IndexNode].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[IndexNode].FastGetSolutionStepValue(DISPLACEMENT,1))[iDoF];
    }
    else
    {
        DeltaPosition = (MasterGeometry[IndexNode - TNumNodes].FastGetSolutionStepValue(DISPLACEMENT) - MasterGeometry[IndexNode - TNumNodes].FastGetSolutionStepValue(DISPLACEMENT,1))[iDoF];
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
double AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional>::GetIntegrationWeight(
    GeneralVariables& rVariables,
    const GeometryType::IntegrationPointsArrayType& ThisIntegrationMethod,
    const unsigned int& PointNumber
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
