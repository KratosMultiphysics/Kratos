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
#include "custom_conditions/ALM_mortar_contact_condition.h"

/* Utilities */
#include "utilities/geometrical_projection_utilities.h"
#include "utilities/math_utils.h"

namespace Kratos
{

/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
Condition::Pointer AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties ) const
{
    KRATOS_ERROR << "You are calling to the base class method Create, check your condition declaration" << std::endl;

    return Kratos::make_shared< AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster> >( NewId, this->GetGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
Condition::Pointer AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_ERROR << "You are calling to the base class method Create, check your condition declaration" << std::endl;

    return Kratos::make_shared< AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster> >( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
Condition::Pointer AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties,
    GeometryType::Pointer pMasterGeom) const
{
    KRATOS_ERROR << "You are calling to the base class method Create, check your condition declaration" << std::endl;

    return Kratos::make_shared< AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster> >( NewId, pGeom, pProperties, pMasterGeom );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::~AugmentedLagrangianMethodMortarContactCondition( )
= default;

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::Initialize( )
{
    KRATOS_TRY;

    BaseType::Initialize();

    mIntegrationOrder = GetProperties().Has(INTEGRATION_ORDER_CONTACT) ? GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) : 2;

    // We reset the ISOLATED flag
    this->Set(ISOLATED, false);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    // NOTE: Add things if necessary

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    // NOTE: Add things if necessary

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
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
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo
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
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
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
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::CalculateMassMatrix(
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

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::CalculateDampingMatrix(
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

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::AddExplicitContribution(ProcessInfo& rCurrentProcessInfo)
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
    const double distance_threshold = rCurrentProcessInfo[DISTANCE_THRESHOLD];
    IntegrationUtility integration_utility = IntegrationUtility (mIntegrationOrder, distance_threshold);

    // If we consider the normal variation
    const NormalDerivativesComputation consider_normal_variation = static_cast<NormalDerivativesComputation>(rCurrentProcessInfo[CONSIDER_NORMAL_VARIATION]);

    // The master geometry
    GeometryType& master_geometry = this->GetPairedGeometry();

    // The normal of the master condition
    const array_1d<double, 3>& normal_master = this->GetValue(PAIRED_NORMAL);

    // Reading integration points
    ConditionArrayListType conditions_points_slave;
    const bool is_inside = integration_utility.GetExactIntegration(slave_geometry, normal_slave, master_geometry, normal_master, conditions_points_slave);

    double integration_area;
    integration_utility.GetTotalArea(slave_geometry, conditions_points_slave, integration_area);

    const double geometry_area = slave_geometry.Area();
    if (is_inside && ((integration_area/geometry_area) > 1.0e-3 * geometry_area)) {
        IntegrationMethod this_integration_method = GetIntegrationMethod();

        // Initialize general variables for the current master element
        rVariables.Initialize();

        // Update slave element info
        rDerivativeData.UpdateMasterPair(master_geometry, rCurrentProcessInfo);

        // Initialize the mortar operators
        rThisMortarConditionMatrices.Initialize();

        const bool dual_LM = DerivativesUtilitiesType::CalculateAeAndDeltaAe(slave_geometry, normal_slave, master_geometry, rDerivativeData, rVariables, consider_normal_variation, conditions_points_slave, this_integration_method, GetAxisymmetricCoefficient(rVariables));

        for (IndexType i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom) {
            PointerVector< PointType > points_array(TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
            array_1d<BelongType, TDim> belong_array;
            for (IndexType i_node = 0; i_node < TDim; ++i_node) {
                PointType global_point;
                slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                points_array(i_node) = Kratos::make_shared<PointType>(PointType(global_point));
                belong_array[i_node] = conditions_points_slave[i_geom][i_node].GetBelong();
            }

            DecompositionType decomp_geom( points_array );

            const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, slave_geometry.Length() * 1.0e-12) : MortarUtilities::HeronCheck(decomp_geom);

            if (bad_shape == false) {
                const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( this_integration_method );

                // Integrating the mortar operators
                for ( IndexType point_number = 0; point_number < integration_points_slave.size(); ++point_number ) {
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
        const BoundedMatrix<double, TNumNodes, TNumNodes>& DOperator = rThisMortarConditionMatrices.DOperator;
        const BoundedMatrix<double, TNumNodes, TNumNodesMaster>& MOperator = rThisMortarConditionMatrices.MOperator;

        // Current coordinates
        const BoundedMatrix<double, TNumNodes, TDim> x1 = MortarUtilities::GetCoordinates<TDim,TNumNodes>(slave_geometry);
        const BoundedMatrix<double, TNumNodesMaster, TDim> x2 = MortarUtilities::GetCoordinates<TDim,TNumNodesMaster>(master_geometry);

        const BoundedMatrix<double, TNumNodes, TDim> D_x1_M_x2 = prod(DOperator, x1) - prod(MOperator, x2);

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            const array_1d<double, 3>& normal = slave_geometry[i_node].FastGetSolutionStepValue(NORMAL);
            array_1d<double, TDim> aux_normal;
            for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
                aux_normal[i_dim] = normal[i_dim];
            }
            const array_1d<double, TDim> aux_array = row(D_x1_M_x2, i_node);

            double& weighted_gap = slave_geometry[i_node].FastGetSolutionStepValue(WEIGHTED_GAP);

            #pragma omp atomic
            weighted_gap += inner_prod(aux_array, - aux_normal);
        }

        // We reset the flag
        this->Set(ISOLATED, false);
    } else {
        // We set the flag
        this->Set(ISOLATED, true);
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void AugmentedLagrangianMethodMortarContactCondition<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::CalculateConditionSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool ComputeLHS,
    const bool ComputeRHS
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
    if (TNormalVariation) DerivativesUtilitiesType::CalculateDeltaNormalSlave(rDerivativeData.DeltaNormalSlave, GetGeometry());

    // Create the mortar operators
    MortarConditionMatrices rThisMortarConditionMatrices;

    // We call the exact integration utility
    const double distance_threshold = rCurrentProcessInfo[DISTANCE_THRESHOLD];
    IntegrationUtility integration_utility = IntegrationUtility (mIntegrationOrder, distance_threshold);

    // The master geometry
    GeometryType& master_geometry = this->GetPairedGeometry();
    const array_1d<double, 3>& normal_master = this->GetValue(PAIRED_NORMAL);

    // Reading integration points
    ConditionArrayListType conditions_points_slave;
    const bool is_inside = CheckIsolatedElement(rCurrentProcessInfo[DELTA_TIME]) ? false : integration_utility.GetExactIntegration(slave_geometry, normal_slave, master_geometry, normal_master, conditions_points_slave);

    double integration_area;
    integration_utility.GetTotalArea(slave_geometry, conditions_points_slave, integration_area);

    const double geometry_area = slave_geometry.Area();
    if (is_inside && ((integration_area/geometry_area) > 1.0e-3 * geometry_area)) {
        IntegrationMethod this_integration_method = GetIntegrationMethod();

        // Initialize general variables for the current master element
        rVariables.Initialize();

        // Update slave element info
        rDerivativeData.UpdateMasterPair(master_geometry, rCurrentProcessInfo);

        // Initialize the mortar operators
        rThisMortarConditionMatrices.Initialize();

        if (TNormalVariation) DerivativesUtilitiesType::CalculateDeltaNormalMaster(rDerivativeData.DeltaNormalMaster, master_geometry);

        const bool dual_LM =  DerivativesUtilitiesType::CalculateAeAndDeltaAe(slave_geometry, normal_slave, master_geometry, rDerivativeData, rVariables, consider_normal_variation, conditions_points_slave, this_integration_method, GetAxisymmetricCoefficient(rVariables));

    #ifdef KRATOS_DEBUG
        if (dual_LM == false)
            KRATOS_WARNING("No dual LM") << "NOT USING DUAL LM. Integration area: " << integration_area << "\tOriginal area: " << geometry_area << "\tRatio: " << integration_area/geometry_area << std::endl;
    #endif

        for (IndexType i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom) {
            PointerVector< PointType > points_array(TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
            array_1d<BelongType, TDim> belong_array;
            for (IndexType i_node = 0; i_node < TDim; ++i_node) {
                PointType global_point;
                slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                points_array(i_node) = Kratos::make_shared<PointType>(PointType(global_point));
                belong_array[i_node] = conditions_points_slave[i_geom][i_node].GetBelong();
            }

            DecompositionType decomp_geom( points_array );

            const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, slave_geometry.Length() * 1.0e-12) : MortarUtilities::HeronCheck(decomp_geom);

            if (bad_shape == false) {
                const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( this_integration_method );

                // Integrating the mortar operators
                for ( IndexType point_number = 0; point_number < integration_points_slave.size(); ++point_number ) {
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

                    if ( ComputeLHS) {
                        /* Update the derivatives */
                        // Update the derivative of the integration vertex (just in 3D)
                        if (TDim == 3) DerivativesUtilitiesType::CalculateDeltaCellVertex(rVariables, rDerivativeData, belong_array, consider_normal_variation, slave_geometry, master_geometry, normal_slave);
                        // Update the derivative of DetJ
                        DerivativesUtilitiesType::CalculateDeltaDetjSlave(decomp_geom, rVariables, rDerivativeData);
                        // Update the derivatives of the shape functions and the gap
                        DerivativesUtilitiesType::CalculateDeltaN(rVariables, rDerivativeData, slave_geometry, master_geometry, normal_slave, normal_master, decomp_geom, local_point_decomp, local_point_parent, consider_normal_variation, dual_LM);

                        rThisMortarConditionMatrices.CalculateDeltaMortarOperators(rVariables, rDerivativeData, integration_weight);
                    } else // In case we are computing RHS we don't compute derivatives (not necessary)
                        rThisMortarConditionMatrices.CalculateMortarOperators(rVariables, integration_weight);
                }
            }
        }

        // Calculates the active/inactive combination pair
        const IndexType active_inactive = GetActiveInactiveValue(slave_geometry);

        // Assemble of the matrix is required
        if ( ComputeLHS )
            this->CalculateLocalLHS( rLeftHandSideMatrix, rThisMortarConditionMatrices, rDerivativeData, active_inactive, rCurrentProcessInfo);

        // Assemble of the vector is required
        if ( ComputeRHS)
            this->CalculateLocalRHS(rRightHandSideVector, rThisMortarConditionMatrices, rDerivativeData, active_inactive, rCurrentProcessInfo);

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

/*********************************COMPUTE KINEMATICS*********************************/
/************************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::CalculateKinematics(
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

    KRATOS_ERROR_IF(rVariables.DetjSlave < 0.0) << "ERROR:: CONDITION ID: " << this->Id() << " INVERTED. DETJ: " << rVariables.DetjSlave << std::endl;

    /// MASTER CONDITION ///
    this->MasterShapeFunctionValue( rVariables, NormalMaster, LocalPointParent);
}

/***********************************************************************************/
/*************** METHODS TO CALCULATE THE CONTACT CONDITION MATRICES ***************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::MasterShapeFunctionValue(
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
    GeometricalProjectionUtilities::FastProjectDirection( master_geometry, slave_gp_global, projected_gp_global, NormalMaster, -gp_normal ); // The opposite direction

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

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
bool AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::CheckIsolatedElement(
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
//     GeometryType& slave_geometry = this->GetGeometry();
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
//     const array_1d<double, 3>& normal_master = this->GetValue(PAIRED_NORMAL);
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<2,2, FrictionalCase::FRICTIONLESS, false, 2>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3,3, FrictionalCase::FRICTIONLESS, false, 3>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3,4, FrictionalCase::FRICTIONLESS, false, 4>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3,3, FrictionalCase::FRICTIONLESS, false, 4>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3,4, FrictionalCase::FRICTIONLESS, false, 3>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<2,2, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 2>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3,3, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 3>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3,4, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 4>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3,3, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 4>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3,4, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 3>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<2,2, FrictionalCase::FRICTIONAL, false, 2>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONAL, false, 3>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONAL, false, 4>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONAL, false, 4>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONAL, false, 3>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<2,2, FrictionalCase::FRICTIONLESS, true, 2>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3,3, FrictionalCase::FRICTIONLESS, true, 3>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3,4, FrictionalCase::FRICTIONLESS, true, 4>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3,3, FrictionalCase::FRICTIONLESS, true, 4>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3,4, FrictionalCase::FRICTIONLESS, true, 3>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<2,2, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 2>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3,3, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 3>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3,4, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 4>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3,3, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 4>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3,4, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 3>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<2,2, FrictionalCase::FRICTIONAL, true, 2>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONAL, true, 3>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONAL, true, 4>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONAL, true, 4>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONAL, true, 3>::CalculateLocalLHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<2, 2, FrictionalCase::FRICTIONLESS, false, 2>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS, false, 3>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS, false, 4>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS, false, 4>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS, false, 3>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<2, 2, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 2>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 3>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 4>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 4>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 3>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<2, 2, FrictionalCase::FRICTIONAL, false, 2>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONAL, false, 3>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONAL, false, 4>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONAL, false, 4>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONAL, false, 3>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<2, 2, FrictionalCase::FRICTIONLESS, true, 2>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS, true, 3>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS, true, 4>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS, true, 4>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS, true, 3>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<2, 2, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 2>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 3>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 4>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 4>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 3>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<2, 2, FrictionalCase::FRICTIONAL, true, 2>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONAL, true, 3>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONAL, true, 4>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONAL, true, 4>::CalculateLocalRHS(
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

template<>
void AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONAL, true, 3>::CalculateLocalRHS(
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
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& CurrentProcessInfo
    )
{
    KRATOS_ERROR << "You are calling to the base class method EquationIdVector, check your condition definition" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void AugmentedLagrangianMethodMortarContactCondition<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::GetDofList(
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_ERROR << "You are calling to the base class method GetDofList, check your condition definition" << std::endl;
}

//******************************* GET DOUBLE VALUE *********************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::GetValueOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//******************************* GET ARRAY_1D VALUE *******************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::GetValueOnIntegrationPoints(
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//******************************* GET VECTOR VALUE *********************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::GetValueOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints();

    if ( rOutput.size() != integration_points.size() )
        rOutput.resize( integration_points.size() );

    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
        rOutput[point_number] = 0.0;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector< array_1d<double, 3 > >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints();

    if ( rOutput.size() != integration_points.size() )
        rOutput.resize( integration_points.size() );

    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
        rOutput[point_number] = ZeroVector(3);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints();

    if ( rOutput.size() != integration_points.size() )
        rOutput.resize( integration_points.size() );

    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
        rOutput[point_number] = ZeroVector(3);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
int AugmentedLagrangianMethodMortarContactCondition<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::Check( const ProcessInfo& rCurrentProcessInfo )
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
    for ( IndexType i = 0; i < TNumNodes; i++ ) {
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

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
double AugmentedLagrangianMethodMortarContactCondition< TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::GetAxisymmetricCoefficient(const GeneralVariables& rVariables) const
{
    return 1.0;
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void AugmentedLagrangianMethodMortarContactCondition< TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::ResizeLHS(MatrixType& rLeftHandSideMatrix)
{
    // Resizing as needed the LHS
    if ( rLeftHandSideMatrix.size1() != MatrixSize || rLeftHandSideMatrix.size2() != MatrixSize )
            rLeftHandSideMatrix.resize( MatrixSize, MatrixSize, false );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void AugmentedLagrangianMethodMortarContactCondition< TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::ResizeRHS(VectorType& rRightHandSideVector)
{
    // Resizing as needed the RHS
    if ( rRightHandSideVector.size() != MatrixSize )
        rRightHandSideVector.resize( MatrixSize, false );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void AugmentedLagrangianMethodMortarContactCondition< TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::ZeroLHS(MatrixType& rLeftHandSideMatrix)
{
    rLeftHandSideMatrix = ZeroMatrix( MatrixSize, MatrixSize );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void AugmentedLagrangianMethodMortarContactCondition< TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::ZeroRHS(VectorType& rRightHandSideVector)
{
    rRightHandSideVector = ZeroVector( MatrixSize );
}

/***********************************************************************************/
/***********************************************************************************/

// Frictionless cases
template class AugmentedLagrangianMethodMortarContactCondition<2, 2, FrictionalCase::FRICTIONLESS, false, 2>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS, false, 3>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS, false, 4>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS, false, 4>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS, false, 3>;
template class AugmentedLagrangianMethodMortarContactCondition<2, 2, FrictionalCase::FRICTIONLESS, true, 2>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS, true, 3>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS, true, 4>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS, true, 4>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS, true, 3>;

// Frictionless components cases
template class AugmentedLagrangianMethodMortarContactCondition<2, 2, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 2>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 3>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 4>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 4>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 3>;
template class AugmentedLagrangianMethodMortarContactCondition<2, 2, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 2>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 3>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 4>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 4>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 3>;

// Frictional cases
template class AugmentedLagrangianMethodMortarContactCondition<2, 2, FrictionalCase::FRICTIONAL, false, 2>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONAL, false, 3>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONAL, false, 4>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONAL, false, 4>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONAL, false, 3>;
template class AugmentedLagrangianMethodMortarContactCondition<2, 2, FrictionalCase::FRICTIONAL, true, 2>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONAL, true, 3>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONAL, true, 4>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 3, FrictionalCase::FRICTIONAL, true, 4>;
template class AugmentedLagrangianMethodMortarContactCondition<3, 4, FrictionalCase::FRICTIONAL, true, 3>;

} // Namespace Kratos
