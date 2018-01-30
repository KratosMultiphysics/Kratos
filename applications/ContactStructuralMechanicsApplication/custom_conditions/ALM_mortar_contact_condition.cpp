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

namespace Kratos 
{
/**
 * Flags related to the condition computation 
 */
// Avoiding using the macro since this has a template parameter. If there was no template plase use the KRATOS_CREATE_LOCAL_FLAG macro
template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
const Kratos::Flags AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::COMPUTE_RHS_VECTOR(Kratos::Flags::Create(0));
template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
const Kratos::Flags AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::COMPUTE_LHS_MATRIX(Kratos::Flags::Create(1));

/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
Condition::Pointer AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::Create( 
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties ) const
{
    KRATOS_ERROR << "You are calling to the base class method Create, check your condition declaration" << std::endl;
    
    return boost::make_shared< AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation> >( NewId, this->GetGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
Condition::Pointer AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_ERROR << "You are calling to the base class method Create, check your condition declaration" << std::endl;
    
    return boost::make_shared< AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation> >( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
Condition::Pointer AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties,
    GeometryType::Pointer pMasterGeom) const
{
    KRATOS_ERROR << "You are calling to the base class method Create, check your condition declaration" << std::endl;
    
    return boost::make_shared< AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation> >( NewId, pGeom, pProperties, pMasterGeom );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::~AugmentedLagrangianMethodMortarContactCondition( )
= default;

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::Initialize( ) 
{
    KRATOS_TRY;
    
    BaseType::Initialize();
    
    mIntegrationOrder = GetProperties().Has(INTEGRATION_ORDER_CONTACT) ? GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) : 2;
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // NOTE: Add things if necessary
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // We reset the flag
    this->Set(ISOLATED, false); // We set the corresponding flag
        
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    // NOTE: Add things if necessary
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;

    // Calculation flags
    mCalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::COMPUTE_LHS_MATRIX, true );
    mCalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::COMPUTE_RHS_VECTOR, true );
    
    // Resizing as needed the LHS
    if ( rLeftHandSideMatrix.size1() != MatrixSize || rLeftHandSideMatrix.size2() != MatrixSize )
            rLeftHandSideMatrix.resize( MatrixSize, MatrixSize, false );
    
    // Resizing as needed the RHS
    if ( rRightHandSideVector.size() != MatrixSize )
        rRightHandSideVector.resize( MatrixSize, false );
    
    // Calculate condition system
    CalculateConditionSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo );

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::CalculateLeftHandSide( 
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Calculation flags
    mCalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::COMPUTE_LHS_MATRIX, true );
    mCalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::COMPUTE_RHS_VECTOR, false);

    // Resizing as needed the LHS
    if ( rLeftHandSideMatrix.size1() != MatrixSize || rLeftHandSideMatrix.size2() != MatrixSize )
        rLeftHandSideMatrix.resize( MatrixSize, MatrixSize, false );
    
    // Creating an auxiliar vector
    VectorType aux_right_hand_side_vector = Vector();

    // Calculate condition system
    CalculateConditionSystem(rLeftHandSideMatrix, aux_right_hand_side_vector, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::CalculateRightHandSide( 
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    // Calculation flags
    mCalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::COMPUTE_LHS_MATRIX, false);
    mCalculationFlags.Set( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::COMPUTE_RHS_VECTOR, true);

    // Creating an auxiliar matrix
    MatrixType aux_left_hand_side_matrix = Matrix();
    
    // Resizing as needed the RHS
    if ( rRightHandSideVector.size() != MatrixSize )
        rRightHandSideVector.resize( MatrixSize, false );

    // Calculate condition system
    CalculateConditionSystem(aux_left_hand_side_matrix, rRightHandSideVector, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::CalculateMassMatrix( 
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

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::CalculateDampingMatrix( 
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

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::AddExplicitContribution(ProcessInfo& rCurrentProcessInfo)
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
    // TODO: Think about the limit
    IntegrationUtility integration_utility = IntegrationUtility (mIntegrationOrder, this->GetGeometry().Length());
    
    // If we consider the normal variation
    const NormalDerivativesComputation consider_normal_variation = static_cast<NormalDerivativesComputation>(rCurrentProcessInfo[CONSIDER_NORMAL_VARIATION]);
    
    // The master geometry
    GeometryType& master_geometry = this->GetPairedGeometry();
    
    // The normal of the master condition
    const array_1d<double, 3>& normal_master = this->GetValue(PAIRED_NORMAL);
    
    // Reading integration points
    ConditionArrayListType conditions_points_slave;
    const bool is_inside = (this->Is(ISOLATED)) ? false : integration_utility.GetExactIntegration(slave_geometry, normal_slave, master_geometry, normal_master, conditions_points_slave);
    
    double integration_area;
    integration_utility.GetTotalArea(slave_geometry, conditions_points_slave, integration_area);
    
    if ((is_inside == true) && ((integration_area/slave_geometry.Area()) > 1.0e-3 * slave_geometry.Area()))
    {            
        IntegrationMethod this_integration_method = GetIntegrationMethod();
        
        // Initialize general variables for the current master element
        rVariables.Initialize();
        
        // Update slave element info
        rDerivativeData.UpdateMasterPair(master_geometry);
        
        // Initialize the mortar operators
        rThisMortarConditionMatrices.Initialize();
        
        const bool dual_LM = DerivativesUtilitiesType::CalculateAeAndDeltaAe(slave_geometry, normal_slave, master_geometry, rDerivativeData, rVariables, consider_normal_variation, conditions_points_slave, this_integration_method, GetAxisymmetricCoefficient(rVariables));
        
        for (unsigned int i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom)
        {
            std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
            array_1d<BelongType, TDim> belong_array;
            for (unsigned int i_node = 0; i_node < TDim; ++i_node)
            {
                PointType global_point;
                slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                points_array[i_node] = PointType::Pointer( new PointType(global_point) );
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
                    this->CalculateKinematics( rVariables, rDerivativeData, normal_master, local_point_decomp, local_point_parent, decomp_geom, dual_LM);
                    
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
            const array_1d<double, 3>& normal = slave_geometry[i_node].FastGetSolutionStepValue(NORMAL);
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
                const array_1d<double, 3>& normal = slave_geometry[i_node].FastGetSolutionStepValue(NORMAL);
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

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
void AugmentedLagrangianMethodMortarContactCondition<TDim, TNumNodes, TFrictional, TNormalVariation>::CalculateConditionSystem( 
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
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
    
    const NormalDerivativesComputation consider_normal_variation = static_cast<NormalDerivativesComputation>(rCurrentProcessInfo[CONSIDER_NORMAL_VARIATION]);
    
    // We compute the normal derivatives
    if (TNormalVariation) DerivativesUtilitiesType::CalculateDeltaNormal(rDerivativeData.DeltaNormalSlave, GetGeometry());
    
    // Create the mortar operators
    MortarConditionMatrices rThisMortarConditionMatrices;
    
    // We call the exact integration utility
    // TODO: Think about the limit
    IntegrationUtility integration_utility = IntegrationUtility (mIntegrationOrder, this->GetGeometry().Length());
    
    // The master geometry
    GeometryType& master_geometry = this->GetPairedGeometry();
    const array_1d<double, 3>& normal_master = this->GetValue(PAIRED_NORMAL);
        
    // Reading integration points
    ConditionArrayListType conditions_points_slave;
    const bool is_inside = (this->Is(ISOLATED)) ? false : integration_utility.GetExactIntegration(slave_geometry, normal_slave, master_geometry, normal_master, conditions_points_slave);
    
    double integration_area;
    integration_utility.GetTotalArea(slave_geometry, conditions_points_slave, integration_area);
    
    if ((is_inside == true) && ((integration_area/slave_geometry.Area()) > 1.0e-3 * slave_geometry.Area()))
    {            
        IntegrationMethod this_integration_method = GetIntegrationMethod();
        
        // Initialize general variables for the current master element
        rVariables.Initialize();
        
        // Update slave element info
        rDerivativeData.UpdateMasterPair(master_geometry);
        
        // Initialize the mortar operators
        rThisMortarConditionMatrices.Initialize();
        
        if (TNormalVariation) DerivativesUtilitiesType::CalculateDeltaNormal(rDerivativeData.DeltaNormalMaster, master_geometry);
        
        const bool dual_LM =  DerivativesUtilitiesType::CalculateAeAndDeltaAe(slave_geometry, normal_slave, master_geometry, rDerivativeData, rVariables, consider_normal_variation, conditions_points_slave, this_integration_method, GetAxisymmetricCoefficient(rVariables));
        
    #ifdef KRATOS_DEBUG
        if (dual_LM == false)
            std::cout << "WARNING:: NOT USING DUAL LM. Integration area: " << integration_area << "\tOriginal area: " << slave_geometry.Area() << "\tRatio: " << integration_area/slave_geometry.Area() << std::endl;
    #endif
        
        for (unsigned int i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom)
        {
            std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
            array_1d<BelongType, TDim> belong_array;
            for (unsigned int i_node = 0; i_node < TDim; ++i_node)
            {
                PointType global_point;
                slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                points_array[i_node] = PointType::Pointer( new PointType(global_point));
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
                    // We reset the derivatives
                    rDerivativeData.ResetDerivatives();
                    
                    // We compute the local coordinates 
                    const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
                    PointType local_point_parent;
                    PointType gp_global;
                    decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                    slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);
                    
                    // Calculate the kinematic variables
                    this->CalculateKinematics( rVariables, rDerivativeData, normal_master, local_point_decomp, local_point_parent, decomp_geom, dual_LM);
                    
                    const double integration_weight = integration_points_slave[point_number].Weight() * GetAxisymmetricCoefficient(rVariables);
                    
                    if ( mCalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::COMPUTE_LHS_MATRIX ))
                    {
                        /* Update the derivatives */
                        // Update the derivative of the integration vertex (just in 3D)
                        if (TDim == 3) DerivativesUtilitiesType::CalculateDeltaCellVertex(rVariables, rDerivativeData, belong_array, consider_normal_variation, slave_geometry, master_geometry, normal_slave);
                        // Update the derivative of DetJ
                        DerivativesUtilitiesType::CalculateDeltaDetjSlave(decomp_geom, rVariables, rDerivativeData);
                        // Update the derivatives of the shape functions and the gap
                        DerivativesUtilitiesType::CalculateDeltaN(rVariables, rDerivativeData, slave_geometry, master_geometry, normal_slave, normal_master, decomp_geom, local_point_decomp, local_point_parent, consider_normal_variation, dual_LM);
                        
                        rThisMortarConditionMatrices.CalculateDeltaMortarOperators(rVariables, rDerivativeData, integration_weight);    
                    }
                    else // In case we are computing RHS we don't compute derivatives (not necessary)
                        rThisMortarConditionMatrices.CalculateMortarOperators(rVariables, integration_weight);
                }
            }
        }
            
        // Calculates the active/inactive combination pair
        const unsigned int active_inactive = GetActiveInactiveValue(slave_geometry);
        
        // Assemble of the matrix is required
        if ( mCalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::COMPUTE_LHS_MATRIX ) )
        {
            const bounded_matrix<double, MatrixSize, MatrixSize>& LHS_contact_pair = this->CalculateLocalLHS( rThisMortarConditionMatrices, rDerivativeData, active_inactive);
            rLeftHandSideMatrix = LHS_contact_pair;
        }
        
        // Assemble of the vector is required
        if ( mCalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::COMPUTE_RHS_VECTOR ))
        {
            const array_1d<double, MatrixSize>& RHS_contact_pair = this->CalculateLocalRHS( rThisMortarConditionMatrices, rDerivativeData, active_inactive);
            rRightHandSideVector = RHS_contact_pair;
        }
    }
    else //If not inside we fill we zero the local matrices
    {
        this->Set(ISOLATED, true); // We set the corresponding flag
        
        // Assemble of the matrix is required
        if ( mCalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::COMPUTE_LHS_MATRIX ) )
            rLeftHandSideMatrix = ZeroMatrix(MatrixSize, MatrixSize);
        
        // Assemble of the vector is required
        if ( mCalculationFlags.Is( AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::COMPUTE_RHS_VECTOR ))
            rRightHandSideVector = ZeroVector(MatrixSize);
    }
    
    KRATOS_CATCH( "" );
}

/*********************************COMPUTE KINEMATICS*********************************/
/************************************************************************************/
 
template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::CalculateKinematics( 
    GeneralVariables& rVariables,
    const DerivativeDataType& rDerivativeData,
    const array_1d<double, 3>& NormalMaster,
    const PointType& LocalPointDecomp,
    const PointType& LocalPointParent,
    GeometryPointType& GeometryDecomp,
    const bool DualLM
    )
{  
    /// SLAVE CONDITION ///
    /* SHAPE FUNCTIONS */
    GetGeometry().ShapeFunctionsValues( rVariables.NSlave, LocalPointParent.Coordinates() );
    rVariables.PhiLagrangeMultipliers = (DualLM == true) ? prod(rDerivativeData.Ae, rVariables.NSlave) : rVariables.NSlave;
    
    /* SHAPE FUNCTION DERIVATIVES */
    GetGeometry().ShapeFunctionsLocalGradients( rVariables.DNDeSlave, LocalPointParent );
    
    /* CALCULATE JACOBIAN AND JACOBIAN DETERMINANT */
    rVariables.jSlave = GeometryDecomp.Jacobian( rVariables.jSlave, LocalPointDecomp.Coordinates());
    rVariables.DetjSlave = GeometryDecomp.DeterminantOfJacobian( LocalPointDecomp );
    
    KRATOS_ERROR_IF(rVariables.DetjSlave < 0.0) << "WARNING:: CONDITION ID: " << this->Id() << " INVERTED. DETJ: " << rVariables.DetjSlave << std::endl;
    
    /// MASTER CONDITION ///
    this->MasterShapeFunctionValue( rVariables, NormalMaster, LocalPointParent);
}
 
/***********************************************************************************/
/*************** METHODS TO CALCULATE THE CONTACT CONDITION MATRICES ***************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::MasterShapeFunctionValue(
    GeneralVariables& rVariables,
    const array_1d<double, 3>& NormalMaster,
    const PointType& LocalPoint
    )
{    
    GeometryType& master_geometry = this->GetPairedGeometry();

    PointType projected_gp_global;
    const array_1d<double,3> gp_normal = MortarUtilities::GaussPointUnitNormal(rVariables.NSlave, GetGeometry());
    
    GeometryType::CoordinatesArrayType slave_gp_global;
    this->GetGeometry( ).GlobalCoordinates( slave_gp_global, LocalPoint );
    MortarUtilities::FastProjectDirection( master_geometry, slave_gp_global, projected_gp_global, NormalMaster, -gp_normal ); // The opposite direction
    
    GeometryType::CoordinatesArrayType projected_gp_local;
    
    master_geometry.PointLocalCoordinates( projected_gp_local, projected_gp_global.Coordinates( ) ) ;

    // SHAPE FUNCTIONS 
    master_geometry.ShapeFunctionsValues(         rVariables.NMaster,    projected_gp_local );         
    master_geometry.ShapeFunctionsLocalGradients( rVariables.DNDeMaster, projected_gp_local );
    
    // JACOBIAN
    rVariables.jMaster = master_geometry.Jacobian( rVariables.jMaster, projected_gp_local);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
bounded_matrix<double, 10, 10> AugmentedLagrangianMethodMortarContactCondition<2,2, false, false>::CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, check your condition definition" << std::endl;
    
    return ZeroMatrix(10, 10);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
bounded_matrix<double, 21, 21> AugmentedLagrangianMethodMortarContactCondition<3,3, false, false>::CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, check your condition definition" << std::endl;
    
    return ZeroMatrix(21, 21);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
bounded_matrix<double, 28, 28> AugmentedLagrangianMethodMortarContactCondition<3,4, false, false>::CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, check your condition definition" << std::endl;
    
    return ZeroMatrix(28, 28);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
bounded_matrix<double, 12, 12> AugmentedLagrangianMethodMortarContactCondition<2,2, true, false>::CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, check your condition definition" << std::endl;
    
    return ZeroMatrix(12, 12);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
bounded_matrix<double, 27, 27> AugmentedLagrangianMethodMortarContactCondition<3, 3, true, false>::CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, check your condition definition" << std::endl;
    
    return ZeroMatrix(27, 27);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
bounded_matrix<double, 36, 36> AugmentedLagrangianMethodMortarContactCondition<3, 4, true, false>::CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, check your condition definition" << std::endl;
    
    return ZeroMatrix(36, 36);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
bounded_matrix<double, 10, 10> AugmentedLagrangianMethodMortarContactCondition<2,2, false, true>::CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, check your condition definition" << std::endl;
    
    return ZeroMatrix(10, 10);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
bounded_matrix<double, 21, 21> AugmentedLagrangianMethodMortarContactCondition<3,3, false, true>::CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, check your condition definition" << std::endl;
    
    return ZeroMatrix(21, 21);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
bounded_matrix<double, 28, 28> AugmentedLagrangianMethodMortarContactCondition<3,4, false, true>::CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, check your condition definition" << std::endl;
    
    return ZeroMatrix(28, 28);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
bounded_matrix<double, 12, 12> AugmentedLagrangianMethodMortarContactCondition<2,2, true, true>::CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, check your condition definition" << std::endl;
    
    return ZeroMatrix(12, 12);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
bounded_matrix<double, 27, 27> AugmentedLagrangianMethodMortarContactCondition<3, 3, true, true>::CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, check your condition definition" << std::endl;
    
    return ZeroMatrix(27, 27);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
bounded_matrix<double, 36, 36> AugmentedLagrangianMethodMortarContactCondition<3, 4, true, true>::CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, check your condition definition" << std::endl;
    
    return ZeroMatrix(36, 36);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double,10> AugmentedLagrangianMethodMortarContactCondition<2, 2, false, false>::CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, check your condition definition" << std::endl;
    
    return ZeroVector(10);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double,21> AugmentedLagrangianMethodMortarContactCondition<3, 3, false, false>::CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, check your condition definition" << std::endl;
    
    return ZeroVector(21);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double,28> AugmentedLagrangianMethodMortarContactCondition<3, 4, false, false>::CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, check your condition definition" << std::endl;
    
    return ZeroVector(28);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double,12> AugmentedLagrangianMethodMortarContactCondition<2, 2, true, false>::CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, check your condition definition" << std::endl;
    
    return ZeroVector(12);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double,27> AugmentedLagrangianMethodMortarContactCondition<3, 3, true, false>::CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, check your condition definition" << std::endl;
    
    return ZeroVector(27);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double,36> AugmentedLagrangianMethodMortarContactCondition<3, 4, true, false>::CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, check your condition definition" << std::endl;
    
    return ZeroVector(36);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double,10> AugmentedLagrangianMethodMortarContactCondition<2, 2, false, true>::CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, check your condition definition" << std::endl;
    
    return ZeroVector(10);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double,21> AugmentedLagrangianMethodMortarContactCondition<3, 3, false, true>::CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, check your condition definition" << std::endl;
    
    return ZeroVector(21);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double,28> AugmentedLagrangianMethodMortarContactCondition<3, 4, false, true>::CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, check your condition definition" << std::endl;
    
    return ZeroVector(28);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double,12> AugmentedLagrangianMethodMortarContactCondition<2, 2, true, true>::CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, check your condition definition" << std::endl;
    
    return ZeroVector(12);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double,27> AugmentedLagrangianMethodMortarContactCondition<3, 3, true, true>::CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, check your condition definition" << std::endl;
    
    return ZeroVector(27);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double,36> AugmentedLagrangianMethodMortarContactCondition<3, 4, true, true>::CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const unsigned int rActiveInactive
        )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, check your condition definition" << std::endl;
    
    return ZeroVector(36);
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& CurrentProcessInfo 
    )
{
    KRATOS_ERROR << "You are calling to the base class method EquationIdVector, check your condition definition" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
void AugmentedLagrangianMethodMortarContactCondition<TDim, TNumNodes, TFrictional, TNormalVariation>::GetDofList(
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo 
)
{
    KRATOS_ERROR << "You are calling to the base class method GetDofList, check your condition definition" << std::endl;
}

//******************************* GET DOUBLE VALUE *********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::GetValueOnIntegrationPoints( 
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//******************************* GET ARRAY_1D VALUE *******************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::GetValueOnIntegrationPoints( 
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//******************************* GET VECTOR VALUE *********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::GetValueOnIntegrationPoints( 
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::CalculateOnIntegrationPoints( 
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints();
        
    if ( rOutput.size() != integration_points.size() )
        rOutput.resize( integration_points.size() );
    
    for (unsigned int point_number = 0; point_number < integration_points.size(); ++point_number)
    {
        rOutput[point_number] = 0.0;
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::CalculateOnIntegrationPoints( 
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector< array_1d<double, 3 > >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
    
    const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints();
        
    if ( rOutput.size() != integration_points.size() )
        rOutput.resize( integration_points.size() );
    
    for (unsigned int point_number = 0; point_number < integration_points.size(); ++point_number)
    {
        rOutput[point_number] = ZeroVector(3);
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::CalculateOnIntegrationPoints( 
    const Variable<Vector>& rVariable, 
    std::vector<Vector>& rOutput, 
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
    
    const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints();
        
    if ( rOutput.size() != integration_points.size() )
        rOutput.resize( integration_points.size() );
    
    for (unsigned int point_number = 0; point_number < integration_points.size(); ++point_number)
    {
        rOutput[point_number] = ZeroVector(3);
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
int AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation>::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int ierr = Condition::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    KRATOS_ERROR_IF(BaseType::mpPairedGeometry == nullptr) << "YOU HAVE NOT INITIALIZED THE PAIR GEOMETRY IN THE AugmentedLagrangianMethodMortarContactCondition" << std::endl;
    
    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
    KRATOS_CHECK_VARIABLE_KEY(WEIGHTED_GAP)
    KRATOS_CHECK_VARIABLE_KEY(NORMAL)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( unsigned int i = 0; i < TNumNodes; i++ )
    {
        Node<3> &rnode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(WEIGHTED_GAP,rnode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NORMAL,rnode)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode)
    }

    return ierr;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
double AugmentedLagrangianMethodMortarContactCondition< TDim, TNumNodes, TFrictional, TNormalVariation>::GetAxisymmetricCoefficient(const GeneralVariables& rVariables) const
{
    return 1.0;
}
    
/***********************************************************************************/
/***********************************************************************************/

// Frictionless cases
template class AugmentedLagrangianMethodMortarContactCondition<2, 2, false, false>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 3, false, false>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 4, false, false>;
template class AugmentedLagrangianMethodMortarContactCondition<2, 2, false, true>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 3, false, true>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 4, false, true>;

// Frictional cases
template class AugmentedLagrangianMethodMortarContactCondition<2, 2, true, false>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 3, true, false>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 4, true, false>;
template class AugmentedLagrangianMethodMortarContactCondition<2, 2, true, true>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 3, true, true>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 4, true, true>;

} // Namespace Kratos
