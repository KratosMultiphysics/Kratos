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
// #include <algorithm>
#ifdef KRATOS_DEBUG
#include <iomanip>
#endif

// External includes

// Project includes
#include "includes/global_variables.h"
#include "custom_conditions/mortar_contact_condition.h"

/* Utilities */
#include "utilities/geometrical_projection_utilities.h"
#include "utilities/math_utils.h"
#include "custom_utilities/mortar_explicit_contribution_utilities.h"

namespace Kratos
{

/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
Condition::Pointer MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties ) const
{
    KRATOS_ERROR << "You are calling to the base class method Create, check your condition declaration" << std::endl;

    return Kratos::make_intrusive< MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster> >( NewId, this->GetParentGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
Condition::Pointer MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_ERROR << "You are calling to the base class method Create, check your condition declaration" << std::endl;

    return Kratos::make_intrusive< MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster> >( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
Condition::Pointer MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties,
    GeometryType::Pointer pMasterGeom) const
{
    KRATOS_ERROR << "You are calling to the base class method Create, check your condition declaration" << std::endl;

    return Kratos::make_intrusive< MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster> >( NewId, pGeom, pProperties, pMasterGeom );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::~MortarContactCondition( )
= default;

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::Initialize(rCurrentProcessInfo);

    // We reset the ISOLATED flag
    this->Set(ISOLATED, false);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::InitializeSolutionStep(rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::InitializeNonLinearIteration(rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::FinalizeSolutionStep(rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::FinalizeNonLinearIteration(rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // Resizing as needed
    ResizeLHS(rLeftHandSideMatrix);
    ResizeRHS(rRightHandSideVector);

    // Calculate condition system
    CalculateConditionSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo );

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // Resizing as needed
    ResizeLHS(rLeftHandSideMatrix);

    // Creating an auxiliar vector
    VectorType aux_right_hand_side_vector = Vector();

    // Calculate condition system
    CalculateConditionSystem(rLeftHandSideMatrix, aux_right_hand_side_vector, rCurrentProcessInfo, true, false );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // Creating an auxiliar matrix
    MatrixType aux_left_hand_side_matrix = Matrix();

    // Resizing as needed
    ResizeRHS(rRightHandSideVector);

    // Calculate condition system
    CalculateConditionSystem(aux_left_hand_side_matrix, rRightHandSideVector, rCurrentProcessInfo, false );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    rMassMatrix.resize(0, 0, false);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    rDampingMatrix.resize(0, 0, false);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::AddExplicitContribution(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const IndexType integration_order = GetProperties().Has(INTEGRATION_ORDER_CONTACT) ? GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) : 2;
    MortarExplicitContributionUtilities<TDim, TNumNodes, FrictionalCase::FRICTIONLESS_PENALTY, TNormalVariation, TNumNodesMaster>::AddExplicitContributionOfMortarCondition(this, rCurrentProcessInfo, integration_order, IsAxisymmetric(), false);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    const Variable<double >& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const IndexType integration_order = this->GetProperties().Has(INTEGRATION_ORDER_CONTACT) ? this->GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) : 2;
    MortarExplicitContributionUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::ComputeNodalArea(this, rCurrentProcessInfo, rDestinationVariable, integration_order, false);
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    const Variable<array_1d<double, 3> >& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_ERROR << "You are calling to the base class method AddExplicitContribution, check your condition definition" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::CalculateConditionSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool ComputeLHS,
    const bool ComputeRHS
    )
{
    KRATOS_TRY;

    // The slave geometry
    const GeometryType& r_slave_geometry = this->GetParentGeometry();
    const array_1d<double, 3>& r_normal_slave = this->GetValue(NORMAL);

    // Create and initialize condition variables
    GeneralVariables general_variables;

    // Create the current contact data
    DerivativeDataType derivative_data;
    derivative_data.Initialize(r_slave_geometry, rCurrentProcessInfo);

    const NormalDerivativesComputation consider_normal_variation = static_cast<NormalDerivativesComputation>(rCurrentProcessInfo[CONSIDER_NORMAL_VARIATION]);

    // We compute the normal derivatives
    if (TNormalVariation) DerivativesUtilitiesType::CalculateDeltaNormalSlave(derivative_data.DeltaNormalSlave, GetParentGeometry());

    // Create the mortar operators
    MortarConditionMatrices mortar_operators;

    // We call the exact integration utility
    const auto& r_properties = this->GetProperties();
    const IndexType integration_order = r_properties.Has(INTEGRATION_ORDER_CONTACT) ? r_properties.GetValue(INTEGRATION_ORDER_CONTACT) : 2;
    const double distance_threshold = rCurrentProcessInfo.Has(DISTANCE_THRESHOLD) ? rCurrentProcessInfo[DISTANCE_THRESHOLD] : 1.0e24;
    const double zero_tolerance_factor = rCurrentProcessInfo.Has(ZERO_TOLERANCE_FACTOR) ? rCurrentProcessInfo[ZERO_TOLERANCE_FACTOR] : 1.0e0;
    const bool consider_tessellation = r_properties.Has(CONSIDER_TESSELLATION) ? r_properties[CONSIDER_TESSELLATION] : false;
    IntegrationUtility integration_utility = IntegrationUtility (integration_order, distance_threshold, 0, zero_tolerance_factor, consider_tessellation);

    // The master geometry
    const GeometryType& r_master_geometry = this->GetPairedGeometry();
    const array_1d<double, 3>& r_normal_master = this->GetPairedNormal();

    // Reading integration points
    ConditionArrayListType conditions_points_slave;
    const bool is_inside = CheckIsolatedElement(rCurrentProcessInfo[DELTA_TIME]) ? false : integration_utility.GetExactIntegration(r_slave_geometry, r_normal_slave, r_master_geometry, r_normal_master, conditions_points_slave);

    double integration_area;
    integration_utility.GetTotalArea(r_slave_geometry, conditions_points_slave, integration_area);

    const double geometry_area = r_slave_geometry.Area();
    if (is_inside && ((integration_area/geometry_area) > 1.0e-5)) {
        IntegrationMethod this_integration_method = GetIntegrationMethod();

        // Initialize general variables for the current master element
        general_variables.Initialize();

        // Update slave element info
        derivative_data.UpdateMasterPair(r_master_geometry, rCurrentProcessInfo);

        // Initialize the mortar operators
        mortar_operators.Initialize();

        if (TNormalVariation) DerivativesUtilitiesType::CalculateDeltaNormalMaster(derivative_data.DeltaNormalMaster, r_master_geometry);

        const bool dual_LM =  DerivativesUtilitiesType::CalculateAeAndDeltaAe(r_slave_geometry, r_normal_slave, r_master_geometry, derivative_data, general_variables, consider_normal_variation, conditions_points_slave, this_integration_method, GetAxisymmetricCoefficient(general_variables));

        for (IndexType i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom) {
            PointerVector< PointType > points_array(TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
            array_1d<BelongType, TDim> belong_array;
            PointType global_point;
            for (IndexType i_node = 0; i_node < TDim; ++i_node) {
                r_slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                points_array(i_node) = Kratos::make_shared<PointType>(PointType(global_point));
                belong_array[i_node] = conditions_points_slave[i_geom][i_node].GetBelong();
            }

            DecompositionType decomp_geom( points_array );

            const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, r_slave_geometry.Length() * 1.0e-12) : MortarUtilities::HeronCheck(decomp_geom);

            if (bad_shape == false) {
                const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( this_integration_method );

                // Integrating the mortar operators
                PointType local_point_parent, gp_global;
                for ( IndexType point_number = 0; point_number < integration_points_slave.size(); ++point_number ) {
                    // We reset the derivatives
                    derivative_data.ResetDerivatives();

                    // We compute the local coordinates
                    const PointType local_point_decomp = PointType(integration_points_slave[point_number].Coordinates());
                    decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                    r_slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);

                    // Calculate the kinematic variables
                    MortarExplicitContributionUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::CalculateKinematics(this, general_variables, derivative_data, r_normal_master, local_point_decomp, local_point_parent, decomp_geom, dual_LM);

                    const double integration_weight = integration_points_slave[point_number].Weight() * GetAxisymmetricCoefficient(general_variables);

                    if ( ComputeLHS) {
                        /* Update the derivatives */
                        // Update the derivative of the integration vertex (just in 3D)
                        if (TDim == 3) DerivativesUtilitiesType::CalculateDeltaCellVertex(general_variables, derivative_data, belong_array, consider_normal_variation, r_slave_geometry, r_master_geometry, r_normal_slave);
                        // Update the derivative of DetJ
                        DerivativesUtilitiesType::CalculateDeltaDetjSlave(decomp_geom, general_variables, derivative_data);
                        // Update the derivatives of the shape functions and the gap
                        DerivativesUtilitiesType::CalculateDeltaN(general_variables, derivative_data, r_slave_geometry, r_master_geometry, r_normal_slave, r_normal_master, decomp_geom, local_point_decomp, local_point_parent, consider_normal_variation, dual_LM);

                        mortar_operators.CalculateDeltaMortarOperators(general_variables, derivative_data, integration_weight);
                    } else // In case we are computing RHS we don't compute derivatives (not necessary)
                        mortar_operators.CalculateMortarOperators(general_variables, integration_weight);
                }
            }
        }

        // Calculates the active/inactive combination pair
        const IndexType active_inactive = GetActiveInactiveValue(r_slave_geometry);

        // Assemble of the matrix is required
        if ( ComputeLHS )
            this->CalculateLocalLHS( rLeftHandSideMatrix, mortar_operators, derivative_data, active_inactive, rCurrentProcessInfo);

        // Assemble of the vector is required
        if ( ComputeRHS)
            this->CalculateLocalRHS(rRightHandSideVector, mortar_operators, derivative_data, active_inactive, rCurrentProcessInfo);

    } else { //If not inside we fill we zero the local matrices
        this->Set(ISOLATED, true); // We set the corresponding flag

        // Assemble of the matrix is required
        if ( ComputeLHS )
            ZeroLHS(rLeftHandSideMatrix);

        // Assemble of the vector is required
        if ( ComputeRHS)
            ZeroRHS(rRightHandSideVector);
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
bool MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::CheckIsolatedElement(
    const double DeltaTime,
    const bool HalfJump
    )
{
    if (this->Is(ISOLATED))
        return true;

//     // We define the tolerance
//     const double tolerance = std::numeric_limits<double>::epsilon();
//
//     // Geometries
//     GeometryType& slave_geometry = this->GetParentGeometry();
//     GeometryType& master_geometry = this->GetPairedGeometry();
//
//     // Dynamic/static
//     const bool dynamic = slave_geometry[0].SolutionStepsDataHas(VELOCITY_X);
//     const double velocity_constant = HalfJump ? 0.25 : 0.5;
//     const double acceleration_constant = HalfJump ? 0.125 : 0.5;
//
//     // Some auxiliar values
//     PointType center_local_coords;
//     Vector N_slave, N_master;
//
//     // Slave geometry
//     slave_geometry.PointLocalCoordinates(center_local_coords, slave_geometry.Center());
//     slave_geometry.ShapeFunctionsValues( N_slave, center_local_coords.Coordinates() );
//
//     MatrixType delta_disp_mat_slave;
//     if (dynamic) {
//         delta_disp_mat_slave = DeltaTime * velocity_constant * (MortarUtilities::GetVariableMatrix<TDim, TNumNodes>(slave_geometry, VELOCITY, 0) + MortarUtilities::GetVariableMatrix<TDim, TNumNodes>(slave_geometry, VELOCITY, 1)) + std::pow(DeltaTime, 2) * acceleration_constant * MortarUtilities::GetVariableMatrix<TDim, TNumNodes>(slave_geometry, ACCELERATION, 1);
//     } else {
//         delta_disp_mat_slave = MortarUtilities::GetVariableMatrix<TDim, TNumNodes>(slave_geometry, DISPLACEMENT, 0) - MortarUtilities::GetVariableMatrix<TDim, TNumNodes>(slave_geometry, DISPLACEMENT, 1);
//     }
//
//     Vector delta_disp_vect_slave = prod(trans(delta_disp_mat_slave), N_slave);
//     if (TDim == 2){
//         delta_disp_vect_slave.resize(3, true);
//         delta_disp_vect_slave[2] = 0.0;
//     }
//
//     // Master geometry
//     master_geometry.PointLocalCoordinates(center_local_coords, master_geometry.Center());
//     master_geometry.ShapeFunctionsValues( N_master, center_local_coords.Coordinates() );
//
//     MatrixType delta_disp_mat_master;
//     if (dynamic) {
//         delta_disp_mat_master = DeltaTime * velocity_constant * (MortarUtilities::GetVariableMatrix<TDim, TNumNodes>(slave_geometry, VELOCITY, 0) + MortarUtilities::GetVariableMatrix<TDim, TNumNodes>(slave_geometry, VELOCITY, 1)) + std::pow(DeltaTime, 2) * acceleration_constant * MortarUtilities::GetVariableMatrix<TDim, TNumNodes>(slave_geometry, ACCELERATION, 1);
//     } else {
//         delta_disp_mat_master = MortarUtilities::GetVariableMatrix<TDim, TNumNodes>(slave_geometry, DISPLACEMENT, 0) - MortarUtilities::GetVariableMatrix<TDim, TNumNodes>(slave_geometry, DISPLACEMENT, 1);
//     }
//
//     Vector delta_disp_vect_master = prod(trans(delta_disp_mat_master), N_master); // TODO: Check multiplciation is consistent
//     if (TDim == 2){
//         delta_disp_vect_master.resize(3, true);
//         delta_disp_vect_slave[2] = 0.0;
//     }
//
//     // We will assume that if nothing it moves can be paired
//     const double norm_slave = norm_2(delta_disp_vect_slave);
//     const double norm_master = norm_2(delta_disp_vect_master);
//     if (norm_slave < tolerance && norm_master < tolerance)
//         return false;
//
//     // We define the normals
//     const array_1d<double, 3> normal_slave = this->GetValue(NORMAL);
//     const array_1d<double, 3>& r_normal_master = this->GetPairedNormal();
//
//     const double angle_slave = MathUtils<double>::VectorsAngle(delta_disp_vect_slave, normal_slave);
//     const double angle_master = MathUtils<double>::VectorsAngle(delta_disp_vect_master, normal_master);
//
//     // In case the both angles are in absolute value minor to angle threshold is active
//     const double angle_threshold = 1.025 * Globals::Pi/2; // We add some tolerance to the angle
//     const double ratio_master_slave = 5.0e-2;
//     if (std::abs(angle_slave) <= angle_threshold || std::abs(angle_master) <= angle_threshold )
//         return false;
//     else if (((std::abs(angle_slave) > angle_threshold) && (norm_slave < ratio_master_slave * norm_master))
//           || ((std::abs(angle_master) > angle_threshold) && (norm_master < ratio_master_slave * norm_slave))) { // In case the angle is greater Pi/2 and the other domain is not moving
// #ifdef KRATOS_DEBUG
//     KRATOS_INFO("ANGLE SLAVE") << "ABS(ANGLE SLAVE) : " << std::abs(angle_slave) << "NORM DELTA DISP: " << norm_slave << std::endl;
//     KRATOS_INFO("ANGLE MASTER") << "ABS(ANGLE MASTER) : " << std::abs(angle_master) << "NORM DELTA DISP: " << norm_master << std::endl;
// #endif
//         return true;
//     } else { // In case both are moving in oposition to the normal
// #ifdef KRATOS_DEBUG
//     KRATOS_INFO("ANGLE SLAVE") << "ABS(ANGLE SLAVE) : " << std::abs(angle_slave) << "NORM DELTA DISP: " << norm_slave << std::endl;
//     KRATOS_INFO("ANGLE MASTER") << "ABS(ANGLE MASTER) : " << std::abs(angle_master) << "NORM DELTA DISP: " << norm_master << std::endl;
// #endif
//         return true;
//     }
//
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::CalculateLocalLHS(
    Matrix& rLocalLHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DerivativeDataType& rDerivativeData,
    const IndexType rActiveInactive,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, check your condition definition" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::CalculateLocalRHS(
    Vector& rLocalRHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DerivativeDataType& rDerivativeData,
    const IndexType rActiveInactive,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, check your condition definition" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& CurrentProcessInfo
    ) const
{
    KRATOS_ERROR << "You are calling to the base class method EquationIdVector, check your condition definition" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::GetDofList(
    DofsVectorType& rConditionalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_ERROR << "You are calling to the base class method GetDofList, check your condition definition" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType &integration_points = GetParentGeometry().IntegrationPoints();

    if ( rOutput.size() != integration_points.size() )
        rOutput.resize( integration_points.size() );

    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
        rOutput[point_number] = 0.0;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector< array_1d<double, 3 > >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType &integration_points = GetParentGeometry().IntegrationPoints();

    if ( rOutput.size() != integration_points.size() )
        rOutput.resize( integration_points.size() );

    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
        rOutput[point_number] = ZeroVector(3);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType &integration_points = GetParentGeometry().IntegrationPoints();

    if ( rOutput.size() != integration_points.size() )
        rOutput.resize( integration_points.size() );

    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
        rOutput[point_number] = ZeroVector(3);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
int MortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int ierr = Condition::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    const GeometryType& r_current_slave = this->GetParentGeometry();
    KRATOS_ERROR_IF(r_current_slave.NumberOfGeometryParts() == 0) << "YOU HAVE NOT INITIALIZED THE PAIR GEOMETRY IN THE MortarContactCondition" << std::endl;

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
    KRATOS_CHECK_VARIABLE_KEY(WEIGHTED_GAP)
    KRATOS_CHECK_VARIABLE_KEY(NORMAL)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( IndexType i = 0; i < TNumNodes; ++i ) {
        const auto& r_node = r_current_slave[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(WEIGHTED_GAP,r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NORMAL,r_node)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, r_node)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, r_node)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, r_node)
    }

    return ierr;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
bool MortarContactCondition< TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::IsAxisymmetric() const
{
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
double MortarContactCondition< TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::GetAxisymmetricCoefficient(const GeneralVariables& rVariables) const
{
    return 1.0;
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition< TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::ResizeLHS(MatrixType& rLeftHandSideMatrix)
{
    // Resizing as needed the LHS
    if ( rLeftHandSideMatrix.size1() != MatrixSize || rLeftHandSideMatrix.size2() != MatrixSize )
            rLeftHandSideMatrix.resize( MatrixSize, MatrixSize, false );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition< TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::ResizeRHS(VectorType& rRightHandSideVector)
{
    // Resizing as needed the RHS
    if ( rRightHandSideVector.size() != MatrixSize )
        rRightHandSideVector.resize( MatrixSize, false );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition< TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::ZeroLHS(MatrixType& rLeftHandSideMatrix)
{
    rLeftHandSideMatrix = ZeroMatrix( MatrixSize, MatrixSize );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarContactCondition< TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::ZeroRHS(VectorType& rRightHandSideVector)
{
    rRightHandSideVector = ZeroVector( MatrixSize );
}

/***********************************************************************************/
/***********************************************************************************/

// Frictionless cases
template class MortarContactCondition<2, 2, FrictionalCase::FRICTIONLESS, false, 2>;
template class MortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS, false, 3>;
template class MortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS, false, 4>;
template class MortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS, false, 4>;
template class MortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS, false, 3>;
template class MortarContactCondition<2, 2, FrictionalCase::FRICTIONLESS, true, 2>;
template class MortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS, true, 3>;
template class MortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS, true, 4>;
template class MortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS, true, 4>;
template class MortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS, true, 3>;

// Frictionless components cases
template class MortarContactCondition<2, 2, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 2>;
template class MortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 3>;
template class MortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 4>;
template class MortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 4>;
template class MortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 3>;
template class MortarContactCondition<2, 2, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 2>;
template class MortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 3>;
template class MortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 4>;
template class MortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 4>;
template class MortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 3>;

// Frictional cases
template class MortarContactCondition<2, 2, FrictionalCase::FRICTIONAL, false, 2>;
template class MortarContactCondition<3, 3, FrictionalCase::FRICTIONAL, false, 3>;
template class MortarContactCondition<3, 4, FrictionalCase::FRICTIONAL, false, 4>;
template class MortarContactCondition<3, 3, FrictionalCase::FRICTIONAL, false, 4>;
template class MortarContactCondition<3, 4, FrictionalCase::FRICTIONAL, false, 3>;
template class MortarContactCondition<2, 2, FrictionalCase::FRICTIONAL, true, 2>;
template class MortarContactCondition<3, 3, FrictionalCase::FRICTIONAL, true, 3>;
template class MortarContactCondition<3, 4, FrictionalCase::FRICTIONAL, true, 4>;
template class MortarContactCondition<3, 3, FrictionalCase::FRICTIONAL, true, 4>;
template class MortarContactCondition<3, 4, FrictionalCase::FRICTIONAL, true, 3>;

// Frictionless penalty cases
template class MortarContactCondition<2, 2, FrictionalCase::FRICTIONLESS_PENALTY, false, 2>;
template class MortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS_PENALTY, false, 3>;
template class MortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS_PENALTY, false, 4>;
template class MortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS_PENALTY, false, 4>;
template class MortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS_PENALTY, false, 3>;
template class MortarContactCondition<2, 2, FrictionalCase::FRICTIONLESS_PENALTY, true, 2>;
template class MortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS_PENALTY, true, 3>;
template class MortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS_PENALTY, true, 4>;
template class MortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS_PENALTY, true, 4>;
template class MortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS_PENALTY, true, 3>;

// Frictional penalty cases
template class MortarContactCondition<2, 2, FrictionalCase::FRICTIONAL_PENALTY, false, 2>;
template class MortarContactCondition<3, 3, FrictionalCase::FRICTIONAL_PENALTY, false, 3>;
template class MortarContactCondition<3, 4, FrictionalCase::FRICTIONAL_PENALTY, false, 4>;
template class MortarContactCondition<3, 3, FrictionalCase::FRICTIONAL_PENALTY, false, 4>;
template class MortarContactCondition<3, 4, FrictionalCase::FRICTIONAL_PENALTY, false, 3>;
template class MortarContactCondition<2, 2, FrictionalCase::FRICTIONAL_PENALTY, true, 2>;
template class MortarContactCondition<3, 3, FrictionalCase::FRICTIONAL_PENALTY, true, 3>;
template class MortarContactCondition<3, 4, FrictionalCase::FRICTIONAL_PENALTY, true, 4>;
template class MortarContactCondition<3, 3, FrictionalCase::FRICTIONAL_PENALTY, true, 4>;
template class MortarContactCondition<3, 4, FrictionalCase::FRICTIONAL_PENALTY, true, 3>;

} // Namespace Kratos
