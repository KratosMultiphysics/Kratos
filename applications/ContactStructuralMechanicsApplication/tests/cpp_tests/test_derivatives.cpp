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
#include "testing/testing.h"
#include "spaces/ublas_space.h"
#include "includes/properties.h"
#include "includes/model_part.h"

/* Utilities */
#include "utilities/geometrical_projection_utilities.h"
#include "utilities/mortar_utilities.h"
#include "utilities/exact_mortar_segmentation_utility.h"
#include "custom_utilities/derivatives_utilities.h"

namespace Kratos
{
    namespace Testing
    {
        typedef Point                                                    PointType;
        typedef Node<3>                                                   NodeType;
        typedef Geometry<NodeType>                                GeometryNodeType;
        typedef Geometry<PointType>                              GeometryPointType;
        typedef Vector                                                  VectorType;
        typedef Matrix                                                  MatrixType;
        typedef std::size_t                                              IndexType;

        ///Type definition for integration methods
        typedef GeometryData::IntegrationMethod                   IntegrationMethod;
        typedef IntegrationPoint<2>                            IntegrationPointType;
        typedef GeometryNodeType::IntegrationPointsArrayType integration_pointsType;

        enum CheckLevel {LEVEL_EXACT = 0, LEVEL_QUADRATIC_CONVERGENCE = 1, LEVEL_DEBUG = 2, LEVEL_FULL_DEBUG = 3};
        enum DerivateToCheck {CHECK_SHAPE_FUNCTION = 0, CHECK_JACOBIAN = 1, CHECK_PHI = 2, CHECK_NORMAL = 3};

        /**
         * This method is used to check the quadratic convergence of the derivatives
         * @param ThisModelPart The model part considered
         * @param SlaveCondition0 The slave condition in the reference configuration
         * @param MasterCondition0 The master condition in the reference configuration
         * @param SlaveCondition1 The slave condition in the current configuration
         * @param MasterCondition1 The master condition in the current configuration
         * @param NodesPerturbation Index of the node to pertubate
         * @param IndexPerturbation Index of the DoF to pertubate
         * @param Coeff Coefficient of perturbation
         * @param Derivative Derivative to check
         * @param Check Check type to consider
         */
        template<std::size_t TDim, std::size_t TNumNodes>
        static inline void TestDerivatives(
            ModelPart& ThisModelPart,
            Condition::Pointer SlaveCondition0,
            Condition::Pointer MasterCondition0,
            Condition::Pointer SlaveCondition1,
            Condition::Pointer MasterCondition1,
            const std::vector<IndexType> NodesPerturbation,
            IndexType IndexPerturbation,
            const std::vector<double> Coeff,
            const IndexType NumberIterations,
            const DerivateToCheck Derivative = CHECK_SHAPE_FUNCTION,
            const CheckLevel Check = LEVEL_QUADRATIC_CONVERGENCE
            )
        {
            // Type definitions
            typedef PointBelong<TNumNodes> PointBelongType;
            typedef array_1d<PointBelongType, TDim> ConditionArrayType;
            typedef typename std::vector<ConditionArrayType> ConditionArrayListType;
            typedef Line2D2<PointType> LineType;
            typedef Triangle3D3<PointType> TriangleType;
            typedef typename std::conditional<TDim == 2, LineType, TriangleType >::type DecompositionType;
            typedef typename std::conditional<TNumNodes == 2, PointBelongsLine2D2N, typename std::conditional<TNumNodes == 3, PointBelongsTriangle3D3N, PointBelongsQuadrilateral3D4N>::type>::type BelongType;
            typedef DerivativesUtilities<TDim, TNumNodes, false, true> DerivativesUtilitiesType;
            typedef ExactMortarIntegrationUtility<TDim, TNumNodes, true> IntegrationUtility;

            // Some definitions
            const NormalDerivativesComputation consider_normal_variation = static_cast<NormalDerivativesComputation>(ThisModelPart.GetProcessInfo()[CONSIDER_NORMAL_VARIATION]);
            const double quadratic_threshold = 1.8;
            const double tolerance = 1.0e-6;

            GeometryType& slave_geometry_0 = SlaveCondition0->GetGeometry();
            GeometryType& master_geometry_0 = MasterCondition0->GetGeometry();
            GeometryType& slave_geometry_1 = SlaveCondition1->GetGeometry();
            GeometryType& master_geometry_1 = MasterCondition1->GetGeometry();

            const array_1d<double, 3>& normal_slave_0 = SlaveCondition0->GetValue(NORMAL);
            const array_1d<double, 3>& normal_master_0 = MasterCondition0->GetValue(NORMAL);

            // Create and initialize condition variables
            MortarKinematicVariablesWithDerivatives<TDim, TNumNodes> rVariables0; // These are the kinematic variables for the initial configuration
            MortarKinematicVariablesWithDerivatives<TDim, TNumNodes> rVariables; // These are the kinematic variables for the current configuration

            // Create the initial contact data
            DerivativeData<TDim, TNumNodes, true> rDerivativeData0;
            rDerivativeData0.Initialize(slave_geometry_0, ThisModelPart.GetProcessInfo());

            // We call the exact integration utility
            IntegrationUtility integration_utility = IntegrationUtility (2);

            Vector error_vector_slave(NumberIterations, 0.0);
            Vector error_vector_master(NumberIterations, 0.0);
            for (IndexType iter = 0; iter < NumberIterations; ++iter) {
                ThisModelPart.GetProcessInfo()[STEP] = iter + 1;
                for (IndexType i_per = 0; i_per < NodesPerturbation.size(); ++i_per) {
                    // We add displacement to the corresponding node
                    array_1d<double, 3> aux_delta_disp = ZeroVector(3);
                    aux_delta_disp[IndexPerturbation] = static_cast<double>(iter + 1) * Coeff[i_per];
                    NodeType::Pointer node_to_move = ThisModelPart.pGetNode(NodesPerturbation[i_per]);
                    array_1d<double, 3>& current_disp = node_to_move->FastGetSolutionStepValue(DISPLACEMENT);
                    const array_1d<double, 3>& aux_previous_disp = node_to_move->FastGetSolutionStepValue(DISPLACEMENT);
                    array_1d<double, 3>& previous_disp = node_to_move->FastGetSolutionStepValue(DISPLACEMENT, 1);
                    previous_disp = aux_previous_disp;
                    current_disp = aux_delta_disp;
                    // Finally we move the mesh
                    noalias(node_to_move->Coordinates()) = node_to_move->GetInitialPosition().Coordinates() + node_to_move->FastGetSolutionStepValue(DISPLACEMENT);
                }

                if (consider_normal_variation != NO_DERIVATIVES_COMPUTATION) {
                    PointType aux_point;
                    aux_point.Coordinates() = ZeroVector(3);
                    const array_1d<double, 3>& normal_slave = slave_geometry_1.UnitNormal(aux_point);
                    SlaveCondition1->SetValue(NORMAL, normal_slave);
                    for (IndexType i_node = 0; i_node < slave_geometry_1.size(); ++i_node) {
                        GeometryType::CoordinatesArrayType point_local;
                        slave_geometry_1.PointLocalCoordinates( point_local, slave_geometry_1[i_node].Coordinates( ) ) ;
                        const array_1d<double, 3>& node_normal_slave = slave_geometry_1.UnitNormal(point_local);
                        array_1d<double, 3>& current_normal = slave_geometry_1[i_node].FastGetSolutionStepValue(NORMAL);
                        array_1d<double, 3>& previous_normal = slave_geometry_1[i_node].FastGetSolutionStepValue(NORMAL, 1);
                        previous_normal = current_normal;
                        current_normal = node_normal_slave;
                    }
                    const array_1d<double, 3>& normal_master = master_geometry_1.UnitNormal(aux_point);
                    MasterCondition1->SetValue(NORMAL, normal_master);
                    for (IndexType i_node = 0; i_node < master_geometry_1.size(); ++i_node) {
                        GeometryType::CoordinatesArrayType point_local;
                        master_geometry_1.PointLocalCoordinates( point_local, master_geometry_1[i_node].Coordinates( ) ) ;
                        const array_1d<double, 3>& node_normal_master = master_geometry_1.UnitNormal(point_local);
                        array_1d<double, 3>& current_normal = master_geometry_1[i_node].FastGetSolutionStepValue(NORMAL);
                        array_1d<double, 3>& previous_normal = master_geometry_1[i_node].FastGetSolutionStepValue(NORMAL, 1);
                        previous_normal = current_normal;
                        current_normal = node_normal_master;
                    }
                }

                // Create the current contact data
                DerivativeData<TDim, TNumNodes, true> rDerivativeData;
                rDerivativeData.Initialize(slave_geometry_1, ThisModelPart.GetProcessInfo());

                // We compute the normal derivatives
                if (consider_normal_variation == NODAL_ELEMENTAL_DERIVATIVES) {
                    // Compute the normal derivatives of the slave
                    DerivativesUtilitiesType::CalculateDeltaNormalSlave(rDerivativeData.DeltaNormalSlave, slave_geometry_1);
                    // Compute the normal derivatives of the master
                    DerivativesUtilitiesType::CalculateDeltaNormalMaster(rDerivativeData.DeltaNormalMaster, master_geometry_1);
                }

                const array_1d<double, 3>& normal_slave_1 = SlaveCondition1->GetValue(NORMAL);
                const array_1d<double, 3>& normal_master_1 = MasterCondition1->GetValue(NORMAL);

                // Reading integration points
                ConditionArrayListType conditions_points_slave0, conditions_points_slave;
                const bool is_inside0 = integration_utility.GetExactIntegration(slave_geometry_0, normal_slave_0, master_geometry_0, normal_master_0, conditions_points_slave0);
                const bool is_inside = integration_utility.GetExactIntegration(slave_geometry_1, normal_slave_1, master_geometry_1, normal_master_1, conditions_points_slave);

                IntegrationMethod this_integration_method = GeometryData::GI_GAUSS_2;

                if (is_inside && is_inside0) {
                    // if (Check == LEVEL_FULL_DEBUG) IntegrationUtility::MathematicaDebug(SlaveCondition1->Id(), slave_geometry_1, MasterCondition1->Id(), master_geometry_1, conditions_points_slave);

                    // Initialize general variables for the current master element
                    rVariables0.Initialize();
                    rVariables.Initialize();

                    // Update slave element info
                    rDerivativeData.UpdateMasterPair(MasterCondition1->GetGeometry(), ThisModelPart.GetProcessInfo());
                    rDerivativeData0.UpdateMasterPair(MasterCondition0->GetGeometry(), ThisModelPart.GetProcessInfo());

                    if (conditions_points_slave.size() == conditions_points_slave0.size()) {// Just in case we have the "same configuration"
                        DerivativesUtilitiesType::CalculateAeAndDeltaAe(slave_geometry_1, normal_slave_1, MasterCondition1->GetGeometry(), rDerivativeData, rVariables, consider_normal_variation, conditions_points_slave, this_integration_method);
                        DerivativesUtilitiesType::CalculateAeAndDeltaAe(slave_geometry_0, normal_slave_0, MasterCondition0->GetGeometry(), rDerivativeData0, rVariables0, consider_normal_variation, conditions_points_slave0, this_integration_method);

                        for (IndexType i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom) {
                            PointerVector< PointType > points_array(TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
                            PointerVector< PointType > points_array0(TDim);
                            array_1d<BelongType, TDim> belong_array;
                            for (IndexType i_node = 0; i_node < TDim; ++i_node) {
                                PointType global_point;
                                slave_geometry_1.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                                points_array(i_node) = Kratos::make_shared<PointType>(PointType(global_point));
                                belong_array[i_node] = conditions_points_slave[i_geom][i_node].GetBelong();
                                slave_geometry_0.GlobalCoordinates(global_point, conditions_points_slave0[i_geom][i_node]);
                                points_array0(i_node) = Kratos::make_shared<PointType>(PointType(global_point));
                            }

                            if (Check == LEVEL_DEBUG || Check == LEVEL_FULL_DEBUG) for (std::size_t i = 0; i < TDim; ++i) KRATOS_WATCH(static_cast<std::size_t>(belong_array[i]));

                            DecompositionType decomp_geom( points_array );
                            DecompositionType decomp_geom0( points_array0 );

                            const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, slave_geometry_0.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);
                            const bool bad_shape0 = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom0, slave_geometry_0.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom0);

                            if ((bad_shape == false) && (bad_shape0 == false)) {
                                const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( this_integration_method );

                                // Integrating the mortar operators
                                for ( IndexType point_number = 0; point_number < integration_points_slave.size(); ++point_number ) {
                                    // We reset the derivatives
                                    rDerivativeData.ResetDerivatives();

                                    // We compute the local coordinates
                                    const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
                                    PointType local_point_parent;
                                    PointType gp_global;

                                    // We compute the initial configuration
                                    decomp_geom0.GlobalCoordinates(gp_global, local_point_decomp);
                                    slave_geometry_0.PointLocalCoordinates(local_point_parent, gp_global);

                                    // SLAVE KINEMATIC COMPUTATIONS
                                    slave_geometry_0.ShapeFunctionsValues( rVariables0.NSlave, local_point_parent.Coordinates() );
                                    slave_geometry_0.ShapeFunctionsLocalGradients( rVariables0.DNDeSlave, local_point_parent );

                                    rVariables0.jSlave = decomp_geom0.Jacobian( rVariables0.jSlave, local_point_decomp.Coordinates());
                                    rVariables0.DetjSlave = decomp_geom0.DeterminantOfJacobian( local_point_decomp );

                                    // MASTER KINEMATIC COMPUTATIONS
                                    PointType projected_gp_global;
                                    array_1d<double,3> gp_normal = MortarUtilities::GaussPointUnitNormal(rVariables0.NSlave, slave_geometry_0);

                                    GeometryType::CoordinatesArrayType slave_gp_global;
                                    slave_geometry_0.GlobalCoordinates( slave_gp_global, local_point_parent );
                                    GeometricalProjectionUtilities::FastProjectDirection( master_geometry_0, slave_gp_global, projected_gp_global, normal_master_0, -gp_normal ); // The opposite direction

                                    GeometryType::CoordinatesArrayType projected_gp_local;

                                    master_geometry_0.PointLocalCoordinates( projected_gp_local, projected_gp_global.Coordinates( ) ) ;

                                    master_geometry_0.ShapeFunctionsValues( rVariables0.NMaster,    projected_gp_local );
                                    master_geometry_0.ShapeFunctionsLocalGradients( rVariables0.DNDeMaster, projected_gp_local );
                                    rVariables0.PhiLagrangeMultipliers = prod(rDerivativeData0.Ae, rVariables0.NSlave);
                                    rVariables0.jMaster = master_geometry_0.Jacobian( rVariables0.jMaster, projected_gp_local);

                                    // We compute the current configuration
                                    decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                                    slave_geometry_1.PointLocalCoordinates(local_point_parent, gp_global);

                                    // SLAVE KINEMATIC COMPUTATIONS
                                    slave_geometry_1.ShapeFunctionsValues( rVariables.NSlave, local_point_parent.Coordinates() );
                                    slave_geometry_1.ShapeFunctionsLocalGradients( rVariables.DNDeSlave, local_point_parent );

                                    rVariables.jSlave = decomp_geom.Jacobian( rVariables.jSlave, local_point_decomp.Coordinates());
                                    rVariables.DetjSlave = decomp_geom.DeterminantOfJacobian( local_point_decomp );

                                    // MASTER KINEMATIC COMPUTATIONS
                                    gp_normal = MortarUtilities::GaussPointUnitNormal(rVariables.NSlave, slave_geometry_1);

                                    slave_geometry_1.GlobalCoordinates( slave_gp_global, local_point_parent );
                                    GeometricalProjectionUtilities::FastProjectDirection( master_geometry_1, slave_gp_global, projected_gp_global, normal_master_1, -gp_normal ); // The opposite direction

                                    master_geometry_1.PointLocalCoordinates( projected_gp_local, projected_gp_global.Coordinates( ) ) ;

                                    master_geometry_1.ShapeFunctionsValues( rVariables.NMaster,    projected_gp_local );
                                    master_geometry_1.ShapeFunctionsLocalGradients( rVariables.DNDeMaster, projected_gp_local );
                                    rVariables.PhiLagrangeMultipliers = prod(rDerivativeData.Ae, rVariables.NSlave);
                                    rVariables.jMaster = master_geometry_1.Jacobian( rVariables.jMaster, projected_gp_local);

                                    // Now we compute the derivatives
                                    if (TDim == 3) DerivativesUtilitiesType::CalculateDeltaCellVertex(rVariables, rDerivativeData, belong_array, consider_normal_variation, slave_geometry_1, master_geometry_1, normal_slave_1);

                                    // Update the derivative of DetJ
                                    DerivativesUtilitiesType::CalculateDeltaDetjSlave(decomp_geom, rVariables, rDerivativeData);

                                    // Update the derivatives of the shape functions and the gap
                                    DerivativesUtilitiesType::CalculateDeltaN(rVariables, rDerivativeData, slave_geometry_1, master_geometry_1, normal_slave_1, normal_master_1, decomp_geom, local_point_decomp, local_point_parent, consider_normal_variation, true);

                                    if (Derivative == CHECK_SHAPE_FUNCTION) {
                                        // Now we compute the error of the delta N
                                        Vector aux_N_dx_slave  = rVariables0.NSlave;
                                        Vector aux_N_dx_master = rVariables0.NMaster;
                                        for (IndexType i_node = 0; i_node < (TDim - 1) * TNumNodes; ++i_node) {
                                            array_1d<double, 3> delta_disp;
                                            if (i_node < TNumNodes) delta_disp = slave_geometry_1[i_node].FastGetSolutionStepValue(DISPLACEMENT);
                                            else delta_disp = master_geometry_1[i_node - TNumNodes].FastGetSolutionStepValue(DISPLACEMENT);
                                            for (IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
                                                const auto& delta_n1 = rDerivativeData.DeltaN1[i_node * TDim + i_dof];
                                                const auto& delta_n2 = rDerivativeData.DeltaN2[i_node * TDim + i_dof];
                                                aux_N_dx_slave  += delta_n1 * delta_disp[i_dof];
                                                aux_N_dx_master += delta_n2 * delta_disp[i_dof];
                                            }
                                        }

                                        if (Check == LEVEL_DEBUG || Check == LEVEL_FULL_DEBUG) KRATOS_WATCH(rVariables0.NSlave)
                                        if (Check == LEVEL_DEBUG || Check == LEVEL_FULL_DEBUG) KRATOS_WATCH(rVariables.NSlave)
                                        if (Check == LEVEL_DEBUG || Check == LEVEL_FULL_DEBUG) KRATOS_WATCH(aux_N_dx_slave)
                                        if (Check == LEVEL_DEBUG || Check == LEVEL_FULL_DEBUG) KRATOS_WATCH(rVariables0.NMaster)
                                        if (Check == LEVEL_DEBUG || Check == LEVEL_FULL_DEBUG) KRATOS_WATCH(rVariables.NMaster)
                                        if (Check == LEVEL_DEBUG || Check == LEVEL_FULL_DEBUG) KRATOS_WATCH(aux_N_dx_master)

                                        error_vector_slave[iter] += norm_2(rVariables.NSlave  - aux_N_dx_slave);
                                        error_vector_master[iter] += norm_2(rVariables.NMaster - aux_N_dx_master);
                                    } else if (Derivative == CHECK_PHI) {
                                        // Now we compute the error of the delta Phi
                                        Matrix aux_Ae_dx_slave = rDerivativeData0.Ae;
                                        Vector aux_Phi_dx_slave  = rVariables0.PhiLagrangeMultipliers;
                                        Vector aux_N_dx_slave  = rVariables0.NSlave;

                                        for (IndexType i_node = 0; i_node < (TDim - 1) * TNumNodes; ++i_node) {
                                            array_1d<double, 3> delta_disp;
                                            if (i_node < TNumNodes) delta_disp = slave_geometry_1[i_node].FastGetSolutionStepValue(DISPLACEMENT);
                                            else delta_disp = master_geometry_1[i_node - TNumNodes].FastGetSolutionStepValue(DISPLACEMENT);
                                            for (IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
                                                const auto& delta_phi = rDerivativeData.DeltaPhi[i_node * TDim + i_dof];
                                                aux_Phi_dx_slave  += delta_phi * delta_disp[i_dof];
                                                const auto& delta_ae = rDerivativeData.DeltaAe[i_node * TDim + i_dof];
                                                aux_Ae_dx_slave += delta_ae * delta_disp[i_dof];
                                                const auto& delta_n1 = rDerivativeData.DeltaN1[i_node * TDim + i_dof];
                                                aux_N_dx_slave += delta_n1 * delta_disp[i_dof];
                                            }
                                        }

                                        if (Check == LEVEL_DEBUG || Check == LEVEL_FULL_DEBUG) KRATOS_WATCH(rDerivativeData.Ae - aux_Ae_dx_slave)
                                        if (Check == LEVEL_DEBUG || Check == LEVEL_FULL_DEBUG) KRATOS_WATCH(rVariables0.NSlave)
                                        if (Check == LEVEL_DEBUG || Check == LEVEL_FULL_DEBUG) KRATOS_WATCH(aux_N_dx_slave)
                                        if (Check == LEVEL_DEBUG || Check == LEVEL_FULL_DEBUG) KRATOS_WATCH(rVariables.NSlave)
                                        if (Check == LEVEL_DEBUG || Check == LEVEL_FULL_DEBUG) KRATOS_WATCH(rVariables0.PhiLagrangeMultipliers)
                                        if (Check == LEVEL_DEBUG || Check == LEVEL_FULL_DEBUG) KRATOS_WATCH(rVariables.PhiLagrangeMultipliers)
                                        if (Check == LEVEL_DEBUG || Check == LEVEL_FULL_DEBUG) KRATOS_WATCH(aux_Phi_dx_slave)

                                        error_vector_slave[iter] += norm_2(rVariables.PhiLagrangeMultipliers - aux_Phi_dx_slave);
                                    } else if (Derivative == CHECK_JACOBIAN) {
                                        // Now we compute the error of the delta Jacobian
                                        double aux_Detj_dx_slave = rVariables0.DetjSlave;
                                        for (IndexType i_node = 0; i_node < (TDim - 1) * TNumNodes; ++i_node) {
                                            array_1d<double, 3> delta_disp;
                                            if (i_node < TNumNodes) delta_disp = slave_geometry_1[i_node].FastGetSolutionStepValue(DISPLACEMENT);
                                            else delta_disp = master_geometry_1[i_node - TNumNodes].FastGetSolutionStepValue(DISPLACEMENT);
                                            for (IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
                                                const auto& delta_detj = rDerivativeData.DeltaDetjSlave[i_node * TDim + i_dof];
                                                aux_Detj_dx_slave += delta_detj * delta_disp[i_dof];
                                            }
                                        }

                                        if (Check == LEVEL_DEBUG || Check == LEVEL_FULL_DEBUG) KRATOS_WATCH(rVariables.DetjSlave)
                                        if (Check == LEVEL_DEBUG || Check == LEVEL_FULL_DEBUG) KRATOS_WATCH(aux_Detj_dx_slave)

                                        error_vector_slave[iter] += std::abs(rVariables.DetjSlave - aux_Detj_dx_slave);
                                    }
                                }
                            }
                        }

                        if (Derivative == CHECK_NORMAL) {
                            // Now we compute the error of the delta Normal
                            auto aux_Normal_dx_slave = rDerivativeData0.NormalSlave;
                            auto aux_Normal_dx_master = rDerivativeData0.NormalMaster;
                            array_1d<double, 2> delta_normal_0_error(2, 0.0);
                            array_1d<double, 3> aux_normal_0 = normal_slave_0;
                            const array_1d<array_1d<double, 3>, TDim * TNumNodes>& delta_normal_0 = DerivativesUtilitiesType::DeltaNormalCenter(slave_geometry_1);

                            for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
                                const array_1d<double, 3>& delta_disp_slave = slave_geometry_1[i_node].FastGetSolutionStepValue(DISPLACEMENT);
                                const array_1d<double, 3>& delta_disp_master = master_geometry_1[i_node].FastGetSolutionStepValue(DISPLACEMENT);

                                for (IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
                                    const auto& delta_normal_slave = rDerivativeData.DeltaNormalSlave[i_node * TDim + i_dof];
                                    const auto& delta_normal_master = rDerivativeData.DeltaNormalMaster[i_node * TDim + i_dof];
                                    aux_Normal_dx_slave  += delta_normal_slave * delta_disp_slave[i_dof];
                                    aux_normal_0 += delta_normal_0[i_node * TDim + i_dof] * delta_disp_slave[i_dof];
                                    aux_Normal_dx_master += delta_normal_master * delta_disp_master[i_dof];
                                }
                            }

                            const double error_delta_n0 = norm_2(aux_normal_0 - SlaveCondition1->GetValue(NORMAL));
                            if (Check == LEVEL_DEBUG || Check == LEVEL_FULL_DEBUG) KRATOS_WATCH(rDerivativeData0.NormalSlave)
                            if (Check == LEVEL_DEBUG || Check == LEVEL_FULL_DEBUG) KRATOS_WATCH(rDerivativeData.NormalSlave)
                            if (Check == LEVEL_DEBUG || Check == LEVEL_FULL_DEBUG) KRATOS_WATCH(aux_Normal_dx_slave)
                            if (Check == LEVEL_DEBUG || Check == LEVEL_FULL_DEBUG) KRATOS_WATCH(rDerivativeData0.NormalMaster)
                            if (Check == LEVEL_DEBUG || Check == LEVEL_FULL_DEBUG) KRATOS_WATCH(rDerivativeData.NormalMaster)
                            if (Check == LEVEL_DEBUG || Check == LEVEL_FULL_DEBUG) KRATOS_WATCH(aux_Normal_dx_master)

                            error_vector_slave[iter] += norm_frobenius(rDerivativeData.NormalSlave - aux_Normal_dx_slave);
                            error_vector_master[iter] += norm_frobenius(rDerivativeData.NormalMaster - aux_Normal_dx_master);

                            if (error_delta_n0 > 0.0) {
                                if (Check == LEVEL_DEBUG || Check == LEVEL_FULL_DEBUG) KRATOS_WATCH(error_delta_n0)
                                if (iter > 0) {
                                    delta_normal_0_error[0] = delta_normal_0_error[1];
                                    delta_normal_0_error[1] = error_delta_n0;
                                    const double log_coeff = std::log((static_cast<double>(iter) + 2.0)/(static_cast<double>(iter) + 1.0));
                                    const double slope = std::log(delta_normal_0_error[1]/delta_normal_0_error[0])/log_coeff;
                                    if (slope >= quadratic_threshold) {
                                        KRATOS_CHECK_GREATER_EQUAL(slope, quadratic_threshold);
                                    } else {
                                        KRATOS_INFO("SLOPE NORMAL 0") << "SLOPE NORMAL 0: " << slope << "\tCurrent:" << delta_normal_0_error[1] << "\tPrevious" << delta_normal_0_error[0] << std::endl;
                                    }
                                } else {
                                    delta_normal_0_error[1] = error_delta_n0;
                                }
                            }
                        }
                    } else
                        KRATOS_ERROR << "YOUR INITIAL SPLITTING DOES NOT COINCIDE WITH THE CURRENT ONE" << std::endl;
                } else
                    KRATOS_ERROR << "WRONG, YOU ARE SUPPOSED TO HAVE AN INTERSECTION" << std::endl;
            }

            if (Check == LEVEL_EXACT) { // LEVEL_EXACT SOLUTION
                for (IndexType iter = 0; iter < NumberIterations; ++iter) {
                    KRATOS_CHECK_LESS_EQUAL(error_vector_slave[iter], tolerance);
                    KRATOS_CHECK_LESS_EQUAL(error_vector_master[iter], tolerance);
                }
            } else if (Check == LEVEL_QUADRATIC_CONVERGENCE) {
                 for (IndexType iter = 0; iter < NumberIterations - 1; ++iter) {
                     const double log_coeff = std::log((static_cast<double>(iter) + 2.0)/(static_cast<double>(iter) + 1.0));

                     // First "exact" check
                     if (error_vector_slave[iter + 1] < tolerance) {
                         KRATOS_CHECK_LESS_EQUAL(error_vector_slave[iter + 1], tolerance);
                     } else {
                         const double slope_slave = std::log(error_vector_slave[iter + 1]/error_vector_slave[iter])/log_coeff;
                         if (slope_slave >= quadratic_threshold) {
                             KRATOS_CHECK_GREATER_EQUAL(slope_slave, quadratic_threshold);
                         } else {
                             KRATOS_INFO("SLOPE SLAVE") << "SLOPE SLAVE: " << slope_slave << "\tCurrent: " << error_vector_slave[iter + 1] << "\tPrevious: " << error_vector_slave[iter] << std::endl;
                         }
                     }

                     // First "exact" check
                     if (error_vector_master[iter + 1] < tolerance) {
                         KRATOS_CHECK_LESS_EQUAL(error_vector_master[iter + 1], tolerance);
                     } else {
                         const double slope_master = std::log(error_vector_master[iter + 1]/error_vector_master[iter])/log_coeff;
                         if (slope_master >= quadratic_threshold) {
                             KRATOS_CHECK_GREATER_EQUAL(slope_master, quadratic_threshold);
                         } else {
                             KRATOS_INFO("SLOPE MASTER") << "SLOPE MASTER: " << slope_master << "\tCurrent: " << error_vector_master[iter + 1] << "\tPrevious: " << error_vector_master[iter] << std::endl;
                         }
                     }
                 }
            } else { // LEVEL_DEBUG
                KRATOS_WATCH(error_vector_slave);
                KRATOS_WATCH(error_vector_master);
            }
        }


        /**
         * Checks if the derivatives of the dual shape functions work as expected
         * Case 1 of the Line2D2
         */

        KRATOS_TEST_CASE_IN_SUITE(DualShapeFunctionDerivativesLine1, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NO_DERIVATIVES_COMPUTATION);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, -1.0,0.0,0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2,  1.0,0.0,0.0);

            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3,  1.2,0.0,0.0);
            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, -0.8,0.0,0.0);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(5, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(6, p_node_2->X(), p_node_2->Y(), p_node_2->Z());

            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(7, p_node_3->X(), p_node_3->Y(), p_node_3->Z());
            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(8, p_node_4->X(), p_node_4->Y(), p_node_4->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (2);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            Line2D2 <NodeType> line_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (2);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            Line2D2 <NodeType> line0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = line_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 1, line_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 3, line0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (2);
            condition_nodes_1[0] = p_node_3;
            condition_nodes_1[1] = p_node_4;
            Line2D2 <NodeType> line_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (2);
            condition_nodes0_1[0] = p_node0_3;
            condition_nodes0_1[1] = p_node0_4;
            Line2D2 <NodeType> line0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = line_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 2, line_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 4, line0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 1);
            std::vector<double> coeff_perturbation(1, -5.0e-2);

            TestDerivatives<2, 2>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 1, coeff_perturbation, 6, CHECK_PHI, LEVEL_QUADRATIC_CONVERGENCE);
        }


        /**
         * Checks if the derivatives of the dual shape functions work as expected
         * Case 2 of the Line2D2
         */

        KRATOS_TEST_CASE_IN_SUITE(DualShapeFunctionDerivativesLine2, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NO_DERIVATIVES_COMPUTATION);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, -1.0,0.0,0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2,  1.0,0.0,0.0);

            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3,  1.2,0.0,0.0);
            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, -0.8,0.0,0.0);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(5, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(6, p_node_2->X(), p_node_2->Y(), p_node_2->Z());

            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(7, p_node_3->X(), p_node_3->Y(), p_node_3->Z());
            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(8, p_node_4->X(), p_node_4->Y(), p_node_4->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (2);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            Line2D2 <NodeType> line_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (2);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            Line2D2 <NodeType> line0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = line_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 1, line_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 3, line0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (2);
            condition_nodes_1[0] = p_node_3;
            condition_nodes_1[1] = p_node_4;
            Line2D2 <NodeType> line_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (2);
            condition_nodes0_1[0] = p_node0_3;
            condition_nodes0_1[1] = p_node0_4;
            Line2D2 <NodeType> line0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = line_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 2, line_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 4, line0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 1);
            std::vector<double> coeff_perturbation(1, -5.0e-2);

            TestDerivatives<2, 2>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 1, coeff_perturbation, 6, CHECK_PHI, LEVEL_QUADRATIC_CONVERGENCE);
        }


        /**
         * Checks if the derivatives of the jacobian work as expected
         * Case 1 of the Line2D2
         */

        KRATOS_TEST_CASE_IN_SUITE(JacobianDerivativesLine1, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NO_DERIVATIVES_COMPUTATION);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, -1.0,0.0,0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2,  1.0,0.0,0.0);

            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3,  1.2,1.0e-3,0.0);
            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, -0.8,1.0e-3,0.0);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(5, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(6, p_node_2->X(), p_node_2->Y(), p_node_2->Z());

            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(7, p_node_3->X(), p_node_3->Y(), p_node_3->Z());
            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(8, p_node_4->X(), p_node_4->Y(), p_node_4->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (2);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            Line2D2 <NodeType> line_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (2);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            Line2D2 <NodeType> line0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = line_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 1, line_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 3, line0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (2);
            condition_nodes_1[0] = p_node_3;
            condition_nodes_1[1] = p_node_4;
            Line2D2 <NodeType> line_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (2);
            condition_nodes0_1[0] = p_node0_3;
            condition_nodes0_1[1] = p_node0_4;
            Line2D2 <NodeType> line0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = line_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 2, line_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 4, line0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 1);
            std::vector<double> coeff_perturbation(1, -5.0e-1);
            TestDerivatives<2, 2>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 1, coeff_perturbation, 1, CHECK_JACOBIAN, LEVEL_QUADRATIC_CONVERGENCE);
        }


        /**
         * Checks if the derivatives of the jacobian work as expected
         * Case 2 of the Line2D2
         */

        KRATOS_TEST_CASE_IN_SUITE(JacobianDerivativesLine2, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NO_DERIVATIVES_COMPUTATION);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, -1.0,0.0,0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2,  1.0,0.0,0.0);

            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3,  1.2,1.0e-3,0.0);
            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, -0.8,1.0e-3,0.0);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(5, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(6, p_node_2->X(), p_node_2->Y(), p_node_2->Z());

            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(7, p_node_3->X(), p_node_3->Y(), p_node_3->Z());
            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(8, p_node_4->X(), p_node_4->Y(), p_node_4->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (2);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            Line2D2 <NodeType> line_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (2);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            Line2D2 <NodeType> line0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = line_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 1, line_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 3, line0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (2);
            condition_nodes_1[0] = p_node_3;
            condition_nodes_1[1] = p_node_4;
            Line2D2 <NodeType> line_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (2);
            condition_nodes0_1[0] = p_node0_3;
            condition_nodes0_1[1] = p_node0_4;
            Line2D2 <NodeType> line0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = line_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 2, line_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 4, line0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 1);
            std::vector<double> coeff_perturbation(1, -5.0e-2);
            TestDerivatives<2, 2>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 1, coeff_perturbation, 6, CHECK_JACOBIAN, LEVEL_QUADRATIC_CONVERGENCE);
        }


        /**
         * Checks if the derivatives of the normal work as expected
         * Case 1 of the Line2D2
         */

        KRATOS_TEST_CASE_IN_SUITE(NormalDerivativesLine1, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NODAL_ELEMENTAL_DERIVATIVES);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, -1.0,0.0,0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2,  1.0,0.0,0.0);

            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3,  1.2,1.0e-3,0.0);
            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, -0.8,1.0e-3,0.0);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(5, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(6, p_node_2->X(), p_node_2->Y(), p_node_2->Z());

            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(7, p_node_3->X(), p_node_3->Y(), p_node_3->Z());
            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(8, p_node_4->X(), p_node_4->Y(), p_node_4->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (2);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            Line2D2 <NodeType> line_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (2);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            Line2D2 <NodeType> line0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = line_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 1, line_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 3, line0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (2);
            condition_nodes_1[0] = p_node_3;
            condition_nodes_1[1] = p_node_4;
            Line2D2 <NodeType> line_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (2);
            condition_nodes0_1[0] = p_node0_3;
            condition_nodes0_1[1] = p_node0_4;
            Line2D2 <NodeType> line0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = line_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 2, line_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 4, line0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1,1);
            std::vector<double> coeff_perturbation(1,1.0e-1);
            //nodes_perturbed[0] = 3;
            TestDerivatives<2, 2>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 1, coeff_perturbation, 6, CHECK_NORMAL, LEVEL_QUADRATIC_CONVERGENCE);
        }


        /**
         * Checks if the derivatives of the normal work as expected
         * Case 2 of the Line2D2
         */

        KRATOS_TEST_CASE_IN_SUITE(NormalDerivativesLine2, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NODAL_ELEMENTAL_DERIVATIVES);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, -1.0,0.0,0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2,  1.0,5.0e-4,0.0);

            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3,  1.2,1.0e-3,0.0);
            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, -0.8,7.0e-4,0.0);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(5, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(6, p_node_2->X(), p_node_2->Y(), p_node_2->Z());

            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(7, p_node_3->X(), p_node_3->Y(), p_node_3->Z());
            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(8, p_node_4->X(), p_node_4->Y(), p_node_4->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (2);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            Line2D2 <NodeType> line_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (2);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            Line2D2 <NodeType> line0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = line_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 1, line_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 3, line0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (2);
            condition_nodes_1[0] = p_node_3;
            condition_nodes_1[1] = p_node_4;
            Line2D2 <NodeType> line_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (2);
            condition_nodes0_1[0] = p_node0_3;
            condition_nodes0_1[1] = p_node0_4;
            Line2D2 <NodeType> line0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = line_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 2, line_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 4, line0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1,1);
            std::vector<double> coeff_perturbation(1,1.0e-1);
            //nodes_perturbed[0] = 3;
            TestDerivatives<2, 2>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 1, coeff_perturbation, 6, CHECK_NORMAL, LEVEL_QUADRATIC_CONVERGENCE);
        }


        /**
         * Checks if the derivatives of the shape functions work as expected
         * Case 1 of the Line2D2
         */

        KRATOS_TEST_CASE_IN_SUITE(ShapeFunctionDerivativesLine1, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NO_DERIVATIVES_COMPUTATION);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, -1.0,0.0,0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2,  1.0,0.0,0.0);

            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3,  1.2,1.0e-3,0.0);
            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, -0.8,1.0e-3,0.0);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(5, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(6, p_node_2->X(), p_node_2->Y(), p_node_2->Z());

            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(7, p_node_3->X(), p_node_3->Y(), p_node_3->Z());
            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(8, p_node_4->X(), p_node_4->Y(), p_node_4->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (2);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            Line2D2 <NodeType> line_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (2);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            Line2D2 <NodeType> line0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = line_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 1, line_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 3, line0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (2);
            condition_nodes_1[0] = p_node_3;
            condition_nodes_1[1] = p_node_4;
            Line2D2 <NodeType> line_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (2);
            condition_nodes0_1[0] = p_node0_3;
            condition_nodes0_1[1] = p_node0_4;
            Line2D2 <NodeType> line0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = line_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 2, line_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 4, line0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 1);
            std::vector<double> coeff_perturbation(1, 5.0e-2);
            TestDerivatives<2, 2>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 1, coeff_perturbation, 6, CHECK_SHAPE_FUNCTION, LEVEL_QUADRATIC_CONVERGENCE);
        }


        /**
         * Checks if the derivatives of the shape functions work as expected
         * Case 2 of the Line2D2
         */

        KRATOS_TEST_CASE_IN_SUITE(ShapeFunctionDerivativesLine2, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NO_DERIVATIVES_COMPUTATION);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, -1.0,0.0,0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2,  1.0,0.0,0.0);

            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3,  1.2,1.0e-3,0.0);
            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, -0.8,1.0e-3,0.0);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(5, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(6, p_node_2->X(), p_node_2->Y(), p_node_2->Z());

            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(7, p_node_3->X(), p_node_3->Y(), p_node_3->Z());
            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(8, p_node_4->X(), p_node_4->Y(), p_node_4->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (2);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            Line2D2 <NodeType> line_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (2);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            Line2D2 <NodeType> line0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = line_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 1, line_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 3, line0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (2);
            condition_nodes_1[0] = p_node_3;
            condition_nodes_1[1] = p_node_4;
            Line2D2 <NodeType> line_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (2);
            condition_nodes0_1[0] = p_node0_3;
            condition_nodes0_1[1] = p_node0_4;
            Line2D2 <NodeType> line0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = line_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 2, line_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 4, line0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 4);
            std::vector<double> coeff_perturbation(1, 5.0e-2);
            TestDerivatives<2, 2>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 1, coeff_perturbation, 6, CHECK_SHAPE_FUNCTION, LEVEL_QUADRATIC_CONVERGENCE);
        }


        /**
         * Checks if the derivatives of the shape functions work as expected
         * Case 3 of the Line2D2
         */

        KRATOS_TEST_CASE_IN_SUITE(ShapeFunctionDerivativesLine3, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NO_DERIVATIVES_COMPUTATION);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, -1.0,0.0,0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2,  1.0,5.0e-4,0.0);

            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3,  1.2,1.0e-3,0.0);
            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, -0.8,7.0e-4,0.0);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(5, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(6, p_node_2->X(), p_node_2->Y(), p_node_2->Z());

            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(7, p_node_3->X(), p_node_3->Y(), p_node_3->Z());
            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(8, p_node_4->X(), p_node_4->Y(), p_node_4->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (2);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            Line2D2 <NodeType> line_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (2);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            Line2D2 <NodeType> line0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = line_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 1, line_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 3, line0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (2);
            condition_nodes_1[0] = p_node_3;
            condition_nodes_1[1] = p_node_4;
            Line2D2 <NodeType> line_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (2);
            condition_nodes0_1[0] = p_node0_3;
            condition_nodes0_1[1] = p_node0_4;
            Line2D2 <NodeType> line0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = line_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 2, line_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 4, line0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 1);
            std::vector<double> coeff_perturbation(1, 5.0e-2);
            TestDerivatives<2, 2>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 1, coeff_perturbation, 6, CHECK_SHAPE_FUNCTION, LEVEL_QUADRATIC_CONVERGENCE);
        }


        /**
         * Checks if the derivatives of the shape functions work as expected
         * Case 4 of the Line2D2
         */

        KRATOS_TEST_CASE_IN_SUITE(ShapeFunctionDerivativesLine4, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NO_DERIVATIVES_COMPUTATION);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, -1.0,0.0,0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2,  1.0,0.5,0.0);

            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3,  1.2,0.501,0.0);
            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, -0.8,0.001,0.0);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(5, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(6, p_node_2->X(), p_node_2->Y(), p_node_2->Z());

            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(7, p_node_3->X(), p_node_3->Y(), p_node_3->Z());
            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(8, p_node_4->X(), p_node_4->Y(), p_node_4->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (2);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            Line2D2 <NodeType> line_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (2);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            Line2D2 <NodeType> line0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = line_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 1, line_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 3, line0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (2);
            condition_nodes_1[0] = p_node_3;
            condition_nodes_1[1] = p_node_4;
            Line2D2 <NodeType> line_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (2);
            condition_nodes0_1[0] = p_node0_3;
            condition_nodes0_1[1] = p_node0_4;
            Line2D2 <NodeType> line0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = line_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 2, line_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition2D2N", 4, line0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 1);
            std::vector<double> coeff_perturbation(1, 5.0e-2);
//             TestDerivatives<2, 2>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 1, coeff_perturbation, 6, CHECK_SHAPE_FUNCTION, LEVEL_QUADRATIC_CONVERGENCE);
        }


        /**
         * Checks if the derivatives of the jacobian work as expected
         * Case 1 of the Triangle3D3
         */

        KRATOS_TEST_CASE_IN_SUITE(JacobianDerivativesTriangle1, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NO_DERIVATIVES_COMPUTATION);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, 0.0,0.0,0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.0,0.0,0.0);
            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 0.0,1.0,0.0);

            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, 0.0,1.0,1.0e-3);
            NodeType::Pointer p_node_5 = model_part.CreateNewNode(5, 0.0,0.0,1.0e-3);
            NodeType::Pointer p_node_6 = model_part.CreateNewNode(6, 1.0,0.0,1.0e-3);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(7, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(8, p_node_2->X(), p_node_2->Y(), p_node_2->Z());
            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(9, p_node_3->X(), p_node_3->Y(), p_node_3->Z());

            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(10, p_node_4->X(), p_node_4->Y(), p_node_4->Z());
            NodeType::Pointer p_node0_5 = model_part.CreateNewNode(11, p_node_5->X(), p_node_5->Y(), p_node_5->Z());
            NodeType::Pointer p_node0_6 = model_part.CreateNewNode(12, p_node_6->X(), p_node_6->Y(), p_node_6->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            Triangle3D3 <NodeType> triangle_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (3);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            condition_nodes0_0[2] = p_node0_3;
            Triangle3D3 <NodeType> triangle0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 1, triangle_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 3, triangle0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            condition_nodes_1[0] = p_node_4;
            condition_nodes_1[1] = p_node_5;
            condition_nodes_1[2] = p_node_6;
            Triangle3D3 <NodeType> triangle_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (3);
            condition_nodes0_1[0] = p_node0_4;
            condition_nodes0_1[1] = p_node0_5;
            condition_nodes0_1[2] = p_node0_6;
            Triangle3D3 <NodeType> triangle0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 2, triangle_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 4, triangle0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 4);
            std::vector<double> coeff_perturbation(1, -5.0e-1);
            TestDerivatives<3, 3>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 1, coeff_perturbation, 1, CHECK_JACOBIAN, LEVEL_EXACT);
        }

        /**
         * Checks if the derivatives of the jacobian work as expected
         * Case 2 of the Triangle3D3
         */

        KRATOS_TEST_CASE_IN_SUITE(JacobianDerivativesTriangle2, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NO_DERIVATIVES_COMPUTATION);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1,-0.1,0.1,1.0e-3);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.1,0.2,0.0);
            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 0.1,1.0,0.0);

            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4,-0.1,1.3,1.0e-3);
            NodeType::Pointer p_node_5 = model_part.CreateNewNode(5, 0.1,0.2,1.0e-3);
            NodeType::Pointer p_node_6 = model_part.CreateNewNode(6, 1.2,0.2,2.0e-3);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(7, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(8, p_node_2->X(), p_node_2->Y(), p_node_2->Z());
            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(9, p_node_3->X(), p_node_3->Y(), p_node_3->Z());

            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(10, p_node_4->X(), p_node_4->Y(), p_node_4->Z());
            NodeType::Pointer p_node0_5 = model_part.CreateNewNode(11, p_node_5->X(), p_node_5->Y(), p_node_5->Z());
            NodeType::Pointer p_node0_6 = model_part.CreateNewNode(12, p_node_6->X(), p_node_6->Y(), p_node_6->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            Triangle3D3 <NodeType> triangle_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (3);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            condition_nodes0_0[2] = p_node0_3;
            Triangle3D3 <NodeType> triangle0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 1, triangle_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 3, triangle0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            condition_nodes_1[0] = p_node_4;
            condition_nodes_1[1] = p_node_5;
            condition_nodes_1[2] = p_node_6;
            Triangle3D3 <NodeType> triangle_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (3);
            condition_nodes0_1[0] = p_node0_4;
            condition_nodes0_1[1] = p_node0_5;
            condition_nodes0_1[2] = p_node0_6;
            Triangle3D3 <NodeType> triangle0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 2, triangle_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 4, triangle0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 4);
            std::vector<double> coeff_perturbation(1, -5.0e-3);
            TestDerivatives<3, 3>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 0, coeff_perturbation, 6, CHECK_JACOBIAN, LEVEL_QUADRATIC_CONVERGENCE);
        }

        /**
         * Checks if the derivatives of the jacobian work as expected
         * Case 1 of the Quadrilateral3D4
         */

        KRATOS_TEST_CASE_IN_SUITE(JacobianDerivativesQuadrilateral1, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NO_DERIVATIVES_COMPUTATION);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, 0.0,0.2,1.0e-3);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.0,0.2,1.0e-3);
            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 1.1,1.1,0.0);
            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, 0.2,1.0,0.0);

            NodeType::Pointer p_node_5 = model_part.CreateNewNode(5,-0.1,1.0,1.0e-3);
            NodeType::Pointer p_node_6 = model_part.CreateNewNode(6, 1.0,1.1,1.0e-3);
            NodeType::Pointer p_node_7 = model_part.CreateNewNode(7, 1.0,0.1,2.0e-3);
            NodeType::Pointer p_node_8 = model_part.CreateNewNode(8, 0.0,0.1,2.0e-3);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(9, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(10, p_node_2->X(), p_node_2->Y(), p_node_2->Z());
            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(11, p_node_3->X(), p_node_3->Y(), p_node_3->Z());
            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(12, p_node_4->X(), p_node_4->Y(), p_node_4->Z());

            NodeType::Pointer p_node0_5 = model_part.CreateNewNode(13, p_node_5->X(), p_node_5->Y(), p_node_5->Z());
            NodeType::Pointer p_node0_6 = model_part.CreateNewNode(14, p_node_6->X(), p_node_6->Y(), p_node_6->Z());
            NodeType::Pointer p_node0_7 = model_part.CreateNewNode(15, p_node_7->X(), p_node_7->Y(), p_node_7->Z());
            NodeType::Pointer p_node0_8 = model_part.CreateNewNode(16, p_node_8->X(), p_node_8->Y(), p_node_8->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (4);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            condition_nodes_0[3] = p_node_4;
            Quadrilateral3D4 <NodeType> quadrilateral_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (4);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            condition_nodes0_0[2] = p_node0_3;
            condition_nodes0_0[3] = p_node0_4;
            Quadrilateral3D4 <NodeType> quadrilateral0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = quadrilateral_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", 1, quadrilateral_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", 3, quadrilateral0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (4);
            condition_nodes_1[0] = p_node_5;
            condition_nodes_1[1] = p_node_6;
            condition_nodes_1[2] = p_node_7;
            condition_nodes_1[3] = p_node_8;
            Quadrilateral3D4 <NodeType> quadrilateral_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (4);
            condition_nodes0_1[0] = p_node0_5;
            condition_nodes0_1[1] = p_node0_6;
            condition_nodes0_1[2] = p_node0_7;
            condition_nodes0_1[3] = p_node0_8;
            Quadrilateral3D4 <NodeType> quadrilateral0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = quadrilateral_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", 2, quadrilateral_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", 4, quadrilateral0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 5);
            std::vector<double> coeff_perturbation(1, -5.0e-3);
            TestDerivatives<3, 4>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 1, coeff_perturbation, 6, CHECK_JACOBIAN, LEVEL_QUADRATIC_CONVERGENCE);
        }

        /**
         * Checks if the derivatives of the shape functions work as expected
         * Case 1 of the Triangle3D3
         */

        KRATOS_TEST_CASE_IN_SUITE(ShapeFunctionDerivativesTriangle1, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NO_DERIVATIVES_COMPUTATION);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, 0.0,0.0,0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.0,0.0,0.0);
            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 0.0,1.0,0.0);

            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, 0.0,1.0,1.0e-3);
            NodeType::Pointer p_node_5 = model_part.CreateNewNode(5, 0.0,0.0,1.0e-3);
            NodeType::Pointer p_node_6 = model_part.CreateNewNode(6, 1.0,0.0,1.0e-3);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(7, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(8, p_node_2->X(), p_node_2->Y(), p_node_2->Z());
            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(9, p_node_3->X(), p_node_3->Y(), p_node_3->Z());

            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(10, p_node_4->X(), p_node_4->Y(), p_node_4->Z());
            NodeType::Pointer p_node0_5 = model_part.CreateNewNode(11, p_node_5->X(), p_node_5->Y(), p_node_5->Z());
            NodeType::Pointer p_node0_6 = model_part.CreateNewNode(12, p_node_6->X(), p_node_6->Y(), p_node_6->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            Triangle3D3 <NodeType> triangle_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (3);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            condition_nodes0_0[2] = p_node0_3;
            Triangle3D3 <NodeType> triangle0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 1, triangle_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 3, triangle0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            condition_nodes_1[0] = p_node_4;
            condition_nodes_1[1] = p_node_5;
            condition_nodes_1[2] = p_node_6;
            Triangle3D3 <NodeType> triangle_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (3);
            condition_nodes0_1[0] = p_node0_4;
            condition_nodes0_1[1] = p_node0_5;
            condition_nodes0_1[2] = p_node0_6;
            Triangle3D3 <NodeType> triangle0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 2, triangle_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 4, triangle0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 4);
            std::vector<double> coeff_perturbation(1, -5.0e-2);
            TestDerivatives<3, 3>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 1, coeff_perturbation, 6, CHECK_SHAPE_FUNCTION, LEVEL_EXACT);
        }

        /**
         * Checks if the derivatives of the shape functions work as expected
         * Case 2 of the Triangle3D3
         */

        KRATOS_TEST_CASE_IN_SUITE(ShapeFunctionDerivativesTriangle2, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NO_DERIVATIVES_COMPUTATION);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, 0.0,0.0,0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.0,0.0,0.0);
            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 0.0,1.0,0.0);

            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4,-0.1,1.0,1.0e-3);
            NodeType::Pointer p_node_5 = model_part.CreateNewNode(5, 0.0,0.0,1.0e-3);
            NodeType::Pointer p_node_6 = model_part.CreateNewNode(6, 1.0,0.0,1.0e-3);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(7, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(8, p_node_2->X(), p_node_2->Y(), p_node_2->Z());
            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(9, p_node_3->X(), p_node_3->Y(), p_node_3->Z());

            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(10, p_node_4->X(), p_node_4->Y(), p_node_4->Z());
            NodeType::Pointer p_node0_5 = model_part.CreateNewNode(11, p_node_5->X(), p_node_5->Y(), p_node_5->Z());
            NodeType::Pointer p_node0_6 = model_part.CreateNewNode(12, p_node_6->X(), p_node_6->Y(), p_node_6->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            Triangle3D3 <NodeType> triangle_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (3);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            condition_nodes0_0[2] = p_node0_3;
            Triangle3D3 <NodeType> triangle0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 1, triangle_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 3, triangle0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            condition_nodes_1[0] = p_node_4;
            condition_nodes_1[1] = p_node_5;
            condition_nodes_1[2] = p_node_6;
            Triangle3D3 <NodeType> triangle_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (3);
            condition_nodes0_1[0] = p_node0_4;
            condition_nodes0_1[1] = p_node0_5;
            condition_nodes0_1[2] = p_node0_6;
            Triangle3D3 <NodeType> triangle0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 2, triangle_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 4, triangle0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 4);
            std::vector<double> coeff_perturbation(1, -5.0e-2);
            TestDerivatives<3, 3>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 1, coeff_perturbation, 6, CHECK_SHAPE_FUNCTION, LEVEL_EXACT);
        }

        /**
         * Checks if the derivatives of the shape functions work as expected
         * Case 3 of the Triangle3D3
         */

        KRATOS_TEST_CASE_IN_SUITE(ShapeFunctionDerivativesTriangle3, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NO_DERIVATIVES_COMPUTATION);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1,-0.1,0.1,1.0e-3);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.1,0.2,0.0);
            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 0.1,1.0,0.0);

            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4,-0.1,1.3,1.0e-3);
            NodeType::Pointer p_node_5 = model_part.CreateNewNode(5, 0.1,0.2,1.0e-3);
            NodeType::Pointer p_node_6 = model_part.CreateNewNode(6, 1.2,0.2,2.0e-3);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(7, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(8, p_node_2->X(), p_node_2->Y(), p_node_2->Z());
            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(9, p_node_3->X(), p_node_3->Y(), p_node_3->Z());

            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(10, p_node_4->X(), p_node_4->Y(), p_node_4->Z());
            NodeType::Pointer p_node0_5 = model_part.CreateNewNode(11, p_node_5->X(), p_node_5->Y(), p_node_5->Z());
            NodeType::Pointer p_node0_6 = model_part.CreateNewNode(12, p_node_6->X(), p_node_6->Y(), p_node_6->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            Triangle3D3 <NodeType> triangle_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (3);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            condition_nodes0_0[2] = p_node0_3;
            Triangle3D3 <NodeType> triangle0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 1, triangle_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 3, triangle0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            condition_nodes_1[0] = p_node_4;
            condition_nodes_1[1] = p_node_5;
            condition_nodes_1[2] = p_node_6;
            Triangle3D3 <NodeType> triangle_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (3);
            condition_nodes0_1[0] = p_node0_4;
            condition_nodes0_1[1] = p_node0_5;
            condition_nodes0_1[2] = p_node0_6;
            Triangle3D3 <NodeType> triangle0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 2, triangle_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 4, triangle0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 4);
            std::vector<double> coeff_perturbation(1, -5.0e-3);
            TestDerivatives<3, 3>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 0, coeff_perturbation, 6, CHECK_SHAPE_FUNCTION, LEVEL_QUADRATIC_CONVERGENCE);
        }

        /**
         * Checks if the derivatives of the shape functions work as expected
         * Case 4 of the Triangle3D3
         */

        KRATOS_TEST_CASE_IN_SUITE(ShapeFunctionDerivativesTriangle4, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(ELEMENTAL_DERIVATIVES);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, 0.0,0.0,0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.0,0.0,0.0);
            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 0.0,1.0,0.0);

            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, 0.0,1.0,1.0e-3);
            NodeType::Pointer p_node_5 = model_part.CreateNewNode(5, 0.0,0.0,1.0e-3);
            NodeType::Pointer p_node_6 = model_part.CreateNewNode(6, 1.0,0.0,1.0e-3);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(7, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(8, p_node_2->X(), p_node_2->Y(), p_node_2->Z());
            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(9, p_node_3->X(), p_node_3->Y(), p_node_3->Z());

            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(10, p_node_4->X(), p_node_4->Y(), p_node_4->Z());
            NodeType::Pointer p_node0_5 = model_part.CreateNewNode(11, p_node_5->X(), p_node_5->Y(), p_node_5->Z());
            NodeType::Pointer p_node0_6 = model_part.CreateNewNode(12, p_node_6->X(), p_node_6->Y(), p_node_6->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            Triangle3D3 <NodeType> triangle_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (3);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            condition_nodes0_0[2] = p_node0_3;
            Triangle3D3 <NodeType> triangle0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 1, triangle_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 3, triangle0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            condition_nodes_1[0] = p_node_4;
            condition_nodes_1[1] = p_node_5;
            condition_nodes_1[2] = p_node_6;
            Triangle3D3 <NodeType> triangle_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (3);
            condition_nodes0_1[0] = p_node0_4;
            condition_nodes0_1[1] = p_node0_5;
            condition_nodes0_1[2] = p_node0_6;
            Triangle3D3 <NodeType> triangle0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 2, triangle_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 4, triangle0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 3);
            std::vector<double> coeff_perturbation(1, -1.0e-3);
//             TestDerivatives<3, 3>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 2, coeff_perturbation, 6, CHECK_SHAPE_FUNCTION, LEVEL_QUADRATIC_CONVERGENCE); // FIXME: Not quadratic convergence when normal derivative is taken into account
        }

        /**
         * Checks if the derivatives of the shape functions work as expected
         * Case 1 of the Quadrilateral3D4
         */

        KRATOS_TEST_CASE_IN_SUITE(ShapeFunctionDerivativesQuadrilateral1, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NO_DERIVATIVES_COMPUTATION);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, 0.0,0.0,0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.0,0.0,0.0);
            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 1.0,1.0,0.0);
            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, 0.0,1.0,0.0);

            NodeType::Pointer p_node_5 = model_part.CreateNewNode(5, 0.0,1.0,1.0e-3);
            NodeType::Pointer p_node_6 = model_part.CreateNewNode(6, 1.0,1.0,1.0e-3);
            NodeType::Pointer p_node_7 = model_part.CreateNewNode(7, 1.0,0.0,1.0e-3);
            NodeType::Pointer p_node_8 = model_part.CreateNewNode(8, 0.0,0.0,1.0e-3);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(9, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(10, p_node_2->X(), p_node_2->Y(), p_node_2->Z());
            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(11, p_node_3->X(), p_node_3->Y(), p_node_3->Z());
            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(12, p_node_4->X(), p_node_4->Y(), p_node_4->Z());

            NodeType::Pointer p_node0_5 = model_part.CreateNewNode(13, p_node_5->X(), p_node_5->Y(), p_node_5->Z());
            NodeType::Pointer p_node0_6 = model_part.CreateNewNode(14, p_node_6->X(), p_node_6->Y(), p_node_6->Z());
            NodeType::Pointer p_node0_7 = model_part.CreateNewNode(15, p_node_7->X(), p_node_7->Y(), p_node_7->Z());
            NodeType::Pointer p_node0_8 = model_part.CreateNewNode(16, p_node_8->X(), p_node_8->Y(), p_node_8->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (4);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            condition_nodes_0[3] = p_node_4;
            Quadrilateral3D4 <NodeType> quadrilateral_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (4);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            condition_nodes0_0[2] = p_node0_3;
            condition_nodes0_0[3] = p_node0_4;
            Quadrilateral3D4 <NodeType> quadrilateral0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = quadrilateral_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", 1, quadrilateral_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", 3, quadrilateral0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (4);
            condition_nodes_1[0] = p_node_5;
            condition_nodes_1[1] = p_node_6;
            condition_nodes_1[2] = p_node_7;
            condition_nodes_1[3] = p_node_8;
            Quadrilateral3D4 <NodeType> quadrilateral_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (4);
            condition_nodes0_1[0] = p_node0_5;
            condition_nodes0_1[1] = p_node0_6;
            condition_nodes0_1[2] = p_node0_7;
            condition_nodes0_1[3] = p_node0_8;
            Quadrilateral3D4 <NodeType> quadrilateral0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = quadrilateral_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", 2, quadrilateral_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", 4, quadrilateral0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 5);
            std::vector<double> coeff_perturbation(1, -5.0e-3);
            TestDerivatives<3, 4>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 1, coeff_perturbation, 6, CHECK_SHAPE_FUNCTION, LEVEL_QUADRATIC_CONVERGENCE);
        }

        /**
         * Checks if the derivatives of the shape functions work as expected
         * Case 2 of the Quadrilateral3D4
         */

        KRATOS_TEST_CASE_IN_SUITE(ShapeFunctionDerivativesQuadrilateral2, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NO_DERIVATIVES_COMPUTATION);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, 0.0,0.0,0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.0,0.0,0.0);
            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 1.0,1.0,0.0);
            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, 0.0,1.0,0.0);

            NodeType::Pointer p_node_5 = model_part.CreateNewNode(5,-0.1,1.0,1.0e-3);
            NodeType::Pointer p_node_6 = model_part.CreateNewNode(6, 1.0,1.0,1.0e-3);
            NodeType::Pointer p_node_7 = model_part.CreateNewNode(7, 1.0,0.0,1.0e-3);
            NodeType::Pointer p_node_8 = model_part.CreateNewNode(8, 0.0,0.0,1.0e-3);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(9, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(10, p_node_2->X(), p_node_2->Y(), p_node_2->Z());
            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(11, p_node_3->X(), p_node_3->Y(), p_node_3->Z());
            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(12, p_node_4->X(), p_node_4->Y(), p_node_4->Z());

            NodeType::Pointer p_node0_5 = model_part.CreateNewNode(13, p_node_5->X(), p_node_5->Y(), p_node_5->Z());
            NodeType::Pointer p_node0_6 = model_part.CreateNewNode(14, p_node_6->X(), p_node_6->Y(), p_node_6->Z());
            NodeType::Pointer p_node0_7 = model_part.CreateNewNode(15, p_node_7->X(), p_node_7->Y(), p_node_7->Z());
            NodeType::Pointer p_node0_8 = model_part.CreateNewNode(16, p_node_8->X(), p_node_8->Y(), p_node_8->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (4);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            condition_nodes_0[3] = p_node_4;
            Quadrilateral3D4 <NodeType> quadrilateral_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (4);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            condition_nodes0_0[2] = p_node0_3;
            condition_nodes0_0[3] = p_node0_4;
            Quadrilateral3D4 <NodeType> quadrilateral0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = quadrilateral_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", 1, quadrilateral_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", 3, quadrilateral0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (4);
            condition_nodes_1[0] = p_node_5;
            condition_nodes_1[1] = p_node_6;
            condition_nodes_1[2] = p_node_7;
            condition_nodes_1[3] = p_node_8;
            Quadrilateral3D4 <NodeType> quadrilateral_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (4);
            condition_nodes0_1[0] = p_node0_5;
            condition_nodes0_1[1] = p_node0_6;
            condition_nodes0_1[2] = p_node0_7;
            condition_nodes0_1[3] = p_node0_8;
            Quadrilateral3D4 <NodeType> quadrilateral0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = quadrilateral_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", 2, quadrilateral_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", 4, quadrilateral0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 5);
            std::vector<double> coeff_perturbation(1, -5.0e-3);
            TestDerivatives<3, 4>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 1, coeff_perturbation, 6, CHECK_SHAPE_FUNCTION, LEVEL_QUADRATIC_CONVERGENCE);
        }

        /**
         * Checks if the derivatives of the shape functions work as expected
         * Case 3 of the Quadrilateral3D4
         */

        KRATOS_TEST_CASE_IN_SUITE(ShapeFunctionDerivativesQuadrilateral3, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NO_DERIVATIVES_COMPUTATION);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, 0.0,0.2,1.0e-3);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.0,0.2,1.0e-3);
            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 1.1,1.1,0.0);
            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, 0.2,1.0,0.0);

            NodeType::Pointer p_node_5 = model_part.CreateNewNode(5,-0.1,1.0,1.0e-3);
            NodeType::Pointer p_node_6 = model_part.CreateNewNode(6, 1.0,1.1,1.0e-3);
            NodeType::Pointer p_node_7 = model_part.CreateNewNode(7, 1.0,0.1,2.0e-3);
            NodeType::Pointer p_node_8 = model_part.CreateNewNode(8, 0.0,0.1,2.0e-3);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(9, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(10, p_node_2->X(), p_node_2->Y(), p_node_2->Z());
            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(11, p_node_3->X(), p_node_3->Y(), p_node_3->Z());
            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(12, p_node_4->X(), p_node_4->Y(), p_node_4->Z());

            NodeType::Pointer p_node0_5 = model_part.CreateNewNode(13, p_node_5->X(), p_node_5->Y(), p_node_5->Z());
            NodeType::Pointer p_node0_6 = model_part.CreateNewNode(14, p_node_6->X(), p_node_6->Y(), p_node_6->Z());
            NodeType::Pointer p_node0_7 = model_part.CreateNewNode(15, p_node_7->X(), p_node_7->Y(), p_node_7->Z());
            NodeType::Pointer p_node0_8 = model_part.CreateNewNode(16, p_node_8->X(), p_node_8->Y(), p_node_8->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (4);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            condition_nodes_0[3] = p_node_4;
            Quadrilateral3D4 <NodeType> quadrilateral_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (4);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            condition_nodes0_0[2] = p_node0_3;
            condition_nodes0_0[3] = p_node0_4;
            Quadrilateral3D4 <NodeType> quadrilateral0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = quadrilateral_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", 1, quadrilateral_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", 3, quadrilateral0_0, p_cond_prop);
            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (4);
            condition_nodes_1[0] = p_node_5;
            condition_nodes_1[1] = p_node_6;
            condition_nodes_1[2] = p_node_7;
            condition_nodes_1[3] = p_node_8;
            Quadrilateral3D4 <NodeType> quadrilateral_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (4);
            condition_nodes0_1[0] = p_node0_5;
            condition_nodes0_1[1] = p_node0_6;
            condition_nodes0_1[2] = p_node0_7;
            condition_nodes0_1[3] = p_node0_8;
            Quadrilateral3D4 <NodeType> quadrilateral0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = quadrilateral_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", 2, quadrilateral_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", 4, quadrilateral0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 5);
            std::vector<double> coeff_perturbation(1, -5.0e-3);
            TestDerivatives<3, 4>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 1, coeff_perturbation, 6, CHECK_SHAPE_FUNCTION, LEVEL_QUADRATIC_CONVERGENCE);
        }

        /**
         * Checks if the derivatives of the dual shape functions work as expected
         * Case 1 of the Triangle3D3
         */

        KRATOS_TEST_CASE_IN_SUITE(DualShapeFunctionDerivativesTriangle1, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NO_DERIVATIVES_COMPUTATION);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, 0.0,0.0,0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.0,0.0,0.0);
            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 0.0,1.0,0.0);

            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, 0.0,1.0,1.0e-3);
            NodeType::Pointer p_node_5 = model_part.CreateNewNode(5, 0.0,0.0,1.0e-3);
            NodeType::Pointer p_node_6 = model_part.CreateNewNode(6, 1.0,0.0,1.0e-3);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(7, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(8, p_node_2->X(), p_node_2->Y(), p_node_2->Z());
            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(9, p_node_3->X(), p_node_3->Y(), p_node_3->Z());

            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(10, p_node_4->X(), p_node_4->Y(), p_node_4->Z());
            NodeType::Pointer p_node0_5 = model_part.CreateNewNode(11, p_node_5->X(), p_node_5->Y(), p_node_5->Z());
            NodeType::Pointer p_node0_6 = model_part.CreateNewNode(12, p_node_6->X(), p_node_6->Y(), p_node_6->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            Triangle3D3 <NodeType> triangle_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (3);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            condition_nodes0_0[2] = p_node0_3;
            Triangle3D3 <NodeType> triangle0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 1, triangle_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 3, triangle0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            condition_nodes_1[0] = p_node_4;
            condition_nodes_1[1] = p_node_5;
            condition_nodes_1[2] = p_node_6;
            Triangle3D3 <NodeType> triangle_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (3);
            condition_nodes0_1[0] = p_node0_4;
            condition_nodes0_1[1] = p_node0_5;
            condition_nodes0_1[2] = p_node0_6;
            Triangle3D3 <NodeType> triangle0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 2, triangle_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 4, triangle0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 4);
            std::vector<double> coeff_perturbation(1, -5.0e-2);
            TestDerivatives<3, 3>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 1, coeff_perturbation, 6, CHECK_PHI, LEVEL_EXACT);
        }

        /**
         * Checks if the derivatives of the dual shape functions work as expected
         * Case 2 of the Triangle3D3
         */

        KRATOS_TEST_CASE_IN_SUITE(DualShapeFunctionDerivativesTriangle2, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NO_DERIVATIVES_COMPUTATION);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1,-0.1,0.1,1.0e-3);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.1,0.2,0.0);
            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 0.1,1.0,0.0);

            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4,-0.1,1.3,1.0e-3);
            NodeType::Pointer p_node_5 = model_part.CreateNewNode(5, 0.1,0.2,1.0e-3);
            NodeType::Pointer p_node_6 = model_part.CreateNewNode(6, 1.2,0.2,2.0e-3);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(7, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(8, p_node_2->X(), p_node_2->Y(), p_node_2->Z());
            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(9, p_node_3->X(), p_node_3->Y(), p_node_3->Z());

            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(10, p_node_4->X(), p_node_4->Y(), p_node_4->Z());
            NodeType::Pointer p_node0_5 = model_part.CreateNewNode(11, p_node_5->X(), p_node_5->Y(), p_node_5->Z());
            NodeType::Pointer p_node0_6 = model_part.CreateNewNode(12, p_node_6->X(), p_node_6->Y(), p_node_6->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            Triangle3D3 <NodeType> triangle_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (3);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            condition_nodes0_0[2] = p_node0_3;
            Triangle3D3 <NodeType> triangle0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 1, triangle_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 3, triangle0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            condition_nodes_1[0] = p_node_4;
            condition_nodes_1[1] = p_node_5;
            condition_nodes_1[2] = p_node_6;
            Triangle3D3 <NodeType> triangle_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (3);
            condition_nodes0_1[0] = p_node0_4;
            condition_nodes0_1[1] = p_node0_5;
            condition_nodes0_1[2] = p_node0_6;
            Triangle3D3 <NodeType> triangle0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 2, triangle_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 4, triangle0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 4);
            std::vector<double> coeff_perturbation(1, -5.0e-3);
            TestDerivatives<3, 3>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 0, coeff_perturbation, 6, CHECK_PHI, LEVEL_QUADRATIC_CONVERGENCE);
        }

        /**
         * Checks if the derivatives of the dual shape functions work as expected
         * Case 1 of the Quadrilateral3D4
         */

        KRATOS_TEST_CASE_IN_SUITE(DualShapeFunctionDerivativesQuadrilateral1, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NO_DERIVATIVES_COMPUTATION);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, 0.0,0.2,1.0e-3);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.0,0.2,1.0e-3);
            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 1.1,1.1,0.0);
            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, 0.2,1.0,0.0);

            NodeType::Pointer p_node_5 = model_part.CreateNewNode(5,-0.1,1.0,1.0e-3);
            NodeType::Pointer p_node_6 = model_part.CreateNewNode(6, 1.0,1.1,1.0e-3);
            NodeType::Pointer p_node_7 = model_part.CreateNewNode(7, 1.0,0.1,2.0e-3);
            NodeType::Pointer p_node_8 = model_part.CreateNewNode(8, 0.0,0.1,2.0e-3);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(9, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(10, p_node_2->X(), p_node_2->Y(), p_node_2->Z());
            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(11, p_node_3->X(), p_node_3->Y(), p_node_3->Z());
            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(12, p_node_4->X(), p_node_4->Y(), p_node_4->Z());

            NodeType::Pointer p_node0_5 = model_part.CreateNewNode(13, p_node_5->X(), p_node_5->Y(), p_node_5->Z());
            NodeType::Pointer p_node0_6 = model_part.CreateNewNode(14, p_node_6->X(), p_node_6->Y(), p_node_6->Z());
            NodeType::Pointer p_node0_7 = model_part.CreateNewNode(15, p_node_7->X(), p_node_7->Y(), p_node_7->Z());
            NodeType::Pointer p_node0_8 = model_part.CreateNewNode(16, p_node_8->X(), p_node_8->Y(), p_node_8->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (4);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            condition_nodes_0[3] = p_node_4;
            Quadrilateral3D4 <NodeType> quadrilateral_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (4);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            condition_nodes0_0[2] = p_node0_3;
            condition_nodes0_0[3] = p_node0_4;
            Quadrilateral3D4 <NodeType> quadrilateral0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = quadrilateral_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", 1, quadrilateral_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", 3, quadrilateral0_0, p_cond_prop);
            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (4);
            condition_nodes_1[0] = p_node_5;
            condition_nodes_1[1] = p_node_6;
            condition_nodes_1[2] = p_node_7;
            condition_nodes_1[3] = p_node_8;
            Quadrilateral3D4 <NodeType> quadrilateral_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (4);
            condition_nodes0_1[0] = p_node0_5;
            condition_nodes0_1[1] = p_node0_6;
            condition_nodes0_1[2] = p_node0_7;
            condition_nodes0_1[3] = p_node0_8;
            Quadrilateral3D4 <NodeType> quadrilateral0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = quadrilateral_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", 2, quadrilateral_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", 4, quadrilateral0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 5);
            std::vector<double> coeff_perturbation(1, -5.0e-3);
            TestDerivatives<3, 4>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 1, coeff_perturbation, 6, CHECK_PHI, LEVEL_QUADRATIC_CONVERGENCE);
        }

        /**
         * Checks if the derivatives of the normal work as expected
         * Case 1 of the Triangle3D3
         */

        KRATOS_TEST_CASE_IN_SUITE(NormalDerivativesTriangle1, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NODAL_ELEMENTAL_DERIVATIVES);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, 0.0,0.0,0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.0,0.0,0.0);
            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 0.0,1.0,0.0);

            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, 0.0,1.0,1.0e-3);
            NodeType::Pointer p_node_5 = model_part.CreateNewNode(5, 0.0,0.0,1.0e-3);
            NodeType::Pointer p_node_6 = model_part.CreateNewNode(6, 1.0,0.0,1.0e-3);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(7, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(8, p_node_2->X(), p_node_2->Y(), p_node_2->Z());
            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(9, p_node_3->X(), p_node_3->Y(), p_node_3->Z());

            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(10, p_node_4->X(), p_node_4->Y(), p_node_4->Z());
            NodeType::Pointer p_node0_5 = model_part.CreateNewNode(11, p_node_5->X(), p_node_5->Y(), p_node_5->Z());
            NodeType::Pointer p_node0_6 = model_part.CreateNewNode(12, p_node_6->X(), p_node_6->Y(), p_node_6->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            Triangle3D3 <NodeType> triangle_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (3);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            condition_nodes0_0[2] = p_node0_3;
            Triangle3D3 <NodeType> triangle0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 1, triangle_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 3, triangle0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            condition_nodes_1[0] = p_node_4;
            condition_nodes_1[1] = p_node_5;
            condition_nodes_1[2] = p_node_6;
            Triangle3D3 <NodeType> triangle_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (3);
            condition_nodes0_1[0] = p_node0_4;
            condition_nodes0_1[1] = p_node0_5;
            condition_nodes0_1[2] = p_node0_6;
            Triangle3D3 <NodeType> triangle0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 2, triangle_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 4, triangle0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 4);
            std::vector<double> coeff_perturbation(1, 5.0e-3);
            TestDerivatives<3, 3>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 2, coeff_perturbation, 6, CHECK_NORMAL, LEVEL_QUADRATIC_CONVERGENCE);
        }

        /**
         * Checks if the derivatives of the normal work as expected
         * Case 2 of the Triangle3D3
         */

        KRATOS_TEST_CASE_IN_SUITE(NormalDerivativesTriangle2, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NODAL_ELEMENTAL_DERIVATIVES);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, 0.0,0.0,0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.0,0.0,0.0);
            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 0.0,1.0,0.0);

            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, 0.0,1.0,1.0e-3);
            NodeType::Pointer p_node_5 = model_part.CreateNewNode(5, 0.0,0.0,1.0e-3);
            NodeType::Pointer p_node_6 = model_part.CreateNewNode(6, 1.0,0.0,1.0e-3);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(7, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(8, p_node_2->X(), p_node_2->Y(), p_node_2->Z());
            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(9, p_node_3->X(), p_node_3->Y(), p_node_3->Z());

            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(10, p_node_4->X(), p_node_4->Y(), p_node_4->Z());
            NodeType::Pointer p_node0_5 = model_part.CreateNewNode(11, p_node_5->X(), p_node_5->Y(), p_node_5->Z());
            NodeType::Pointer p_node0_6 = model_part.CreateNewNode(12, p_node_6->X(), p_node_6->Y(), p_node_6->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            Triangle3D3 <NodeType> triangle_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (3);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            condition_nodes0_0[2] = p_node0_3;
            Triangle3D3 <NodeType> triangle0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 1, triangle_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 3, triangle0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            condition_nodes_1[0] = p_node_4;
            condition_nodes_1[1] = p_node_5;
            condition_nodes_1[2] = p_node_6;
            Triangle3D3 <NodeType> triangle_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (3);
            condition_nodes0_1[0] = p_node0_4;
            condition_nodes0_1[1] = p_node0_5;
            condition_nodes0_1[2] = p_node0_6;
            Triangle3D3 <NodeType> triangle0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 2, triangle_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 4, triangle0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 3);
            std::vector<double> coeff_perturbation(1, -5.0e-2);
            TestDerivatives<3, 3>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 2, coeff_perturbation, 1, CHECK_NORMAL, LEVEL_QUADRATIC_CONVERGENCE);
        }

        /**
         * Checks if the derivatives of the normal work as expected
         * Case 3 of the Triangle3D3
         */

        KRATOS_TEST_CASE_IN_SUITE(NormalDerivativesTriangle3, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NODAL_ELEMENTAL_DERIVATIVES);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, 0.0,0.0,0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.0,0.0,0.0);
            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 0.0,1.0,0.0);

            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, 0.0,1.0,1.0e-3);
            NodeType::Pointer p_node_5 = model_part.CreateNewNode(5, 0.0,0.0,1.0e-3);
            NodeType::Pointer p_node_6 = model_part.CreateNewNode(6, 1.0,0.0,1.0e-3);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(7, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(8, p_node_2->X(), p_node_2->Y(), p_node_2->Z());
            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(9, p_node_3->X(), p_node_3->Y(), p_node_3->Z());

            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(10, p_node_4->X(), p_node_4->Y(), p_node_4->Z());
            NodeType::Pointer p_node0_5 = model_part.CreateNewNode(11, p_node_5->X(), p_node_5->Y(), p_node_5->Z());
            NodeType::Pointer p_node0_6 = model_part.CreateNewNode(12, p_node_6->X(), p_node_6->Y(), p_node_6->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            Triangle3D3 <NodeType> triangle_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (3);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            condition_nodes0_0[2] = p_node0_3;
            Triangle3D3 <NodeType> triangle0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 1, triangle_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 3, triangle0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            condition_nodes_1[0] = p_node_4;
            condition_nodes_1[1] = p_node_5;
            condition_nodes_1[2] = p_node_6;
            Triangle3D3 <NodeType> triangle_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (3);
            condition_nodes0_1[0] = p_node0_4;
            condition_nodes0_1[1] = p_node0_5;
            condition_nodes0_1[2] = p_node0_6;
            Triangle3D3 <NodeType> triangle0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 2, triangle_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 4, triangle0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(2);
            std::vector<double> coeff_perturbation(2, 1.0e-1);
            nodes_perturbed[0] = 4;
            nodes_perturbed[1] = 5;
            TestDerivatives<3, 3>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 2, coeff_perturbation, 6, CHECK_NORMAL, LEVEL_QUADRATIC_CONVERGENCE);
        }

        /**
         * Checks if the derivatives of the normal work as expected
         * Case 4 of the Triangle3D3
         */

        KRATOS_TEST_CASE_IN_SUITE(NormalDerivativesTriangle4, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NODAL_ELEMENTAL_DERIVATIVES);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, 0.0,0.0,0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.0,0.0,0.0);
            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 0.0,1.0,0.0);

            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, 0.0,1.0,1.0e-3);
            NodeType::Pointer p_node_5 = model_part.CreateNewNode(5, 0.0,0.0,1.0e-3);
            NodeType::Pointer p_node_6 = model_part.CreateNewNode(6, 1.0,0.0,1.0e-3);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(7, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(8, p_node_2->X(), p_node_2->Y(), p_node_2->Z());
            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(9, p_node_3->X(), p_node_3->Y(), p_node_3->Z());

            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(10, p_node_4->X(), p_node_4->Y(), p_node_4->Z());
            NodeType::Pointer p_node0_5 = model_part.CreateNewNode(11, p_node_5->X(), p_node_5->Y(), p_node_5->Z());
            NodeType::Pointer p_node0_6 = model_part.CreateNewNode(12, p_node_6->X(), p_node_6->Y(), p_node_6->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            Triangle3D3 <NodeType> triangle_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (3);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            condition_nodes0_0[2] = p_node0_3;
            Triangle3D3 <NodeType> triangle0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 1, triangle_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 3, triangle0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            condition_nodes_1[0] = p_node_4;
            condition_nodes_1[1] = p_node_5;
            condition_nodes_1[2] = p_node_6;
            Triangle3D3 <NodeType> triangle_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (3);
            condition_nodes0_1[0] = p_node0_4;
            condition_nodes0_1[1] = p_node0_5;
            condition_nodes0_1[2] = p_node0_6;
            Triangle3D3 <NodeType> triangle0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 2, triangle_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 4, triangle0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(2);
            std::vector<double> coeff_perturbation(2);
            nodes_perturbed[0] = 4;
            nodes_perturbed[1] = 5;
            coeff_perturbation[0] = 1.0e-1;
            coeff_perturbation[1] = 5.0e-2;
            TestDerivatives<3, 3>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 2, coeff_perturbation, 6, CHECK_NORMAL, LEVEL_QUADRATIC_CONVERGENCE);
        }

        /**
         * Checks if the derivatives of the normal work as expected
         * Case 5 of the Triangle3D3
         */

        KRATOS_TEST_CASE_IN_SUITE(NormalDerivativesTriangle5, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NODAL_ELEMENTAL_DERIVATIVES);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1,-0.1,0.1,1.0e-3);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.1,0.2,0.0);
            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 0.1,1.0,0.0);

            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4,-0.1,1.3,1.0e-3);
            NodeType::Pointer p_node_5 = model_part.CreateNewNode(5, 0.1,0.2,1.0e-3);
            NodeType::Pointer p_node_6 = model_part.CreateNewNode(6, 1.2,0.2,2.0e-3);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(7, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(8, p_node_2->X(), p_node_2->Y(), p_node_2->Z());
            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(9, p_node_3->X(), p_node_3->Y(), p_node_3->Z());

            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(10, p_node_4->X(), p_node_4->Y(), p_node_4->Z());
            NodeType::Pointer p_node0_5 = model_part.CreateNewNode(11, p_node_5->X(), p_node_5->Y(), p_node_5->Z());
            NodeType::Pointer p_node0_6 = model_part.CreateNewNode(12, p_node_6->X(), p_node_6->Y(), p_node_6->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            Triangle3D3 <NodeType> triangle_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (3);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            condition_nodes0_0[2] = p_node0_3;
            Triangle3D3 <NodeType> triangle0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 1, triangle_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 3, triangle0_0, p_cond_prop);

            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            condition_nodes_1[0] = p_node_4;
            condition_nodes_1[1] = p_node_5;
            condition_nodes_1[2] = p_node_6;
            Triangle3D3 <NodeType> triangle_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (3);
            condition_nodes0_1[0] = p_node0_4;
            condition_nodes0_1[1] = p_node0_5;
            condition_nodes0_1[2] = p_node0_6;
            Triangle3D3 <NodeType> triangle0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 2, triangle_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 4, triangle0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 1);
            std::vector<double> coeff_perturbation(1, 5.0e-2);
            TestDerivatives<3, 3>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 2, coeff_perturbation, 6, CHECK_NORMAL, LEVEL_QUADRATIC_CONVERGENCE);
        }

         /**
         * Checks if the derivatives of the normal functions work as expected
         * Case 1 of the Quadrilateral3D4
         */

        KRATOS_TEST_CASE_IN_SUITE(NormalDerivativesQuadrilateral1, KratosContactStructuralMechanicsFastSuite)
        {
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            model_part.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] = static_cast<int>(NODAL_ELEMENTAL_DERIVATIVES);

            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);

            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(NORMAL);

            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);

            // First we create the nodes
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, 0.0,0.2,1.0e-3);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.0,0.2,1.0e-3);
            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 1.1,1.1,0.0);
            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, 0.2,1.0,0.0);

            NodeType::Pointer p_node_5 = model_part.CreateNewNode(5,-0.1,1.0,1.0e-3);
            NodeType::Pointer p_node_6 = model_part.CreateNewNode(6, 1.0,1.1,1.0e-3);
            NodeType::Pointer p_node_7 = model_part.CreateNewNode(7, 1.0,0.1,2.0e-3);
            NodeType::Pointer p_node_8 = model_part.CreateNewNode(8, 0.0,0.1,2.0e-3);

            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(9, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(10, p_node_2->X(), p_node_2->Y(), p_node_2->Z());
            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(11, p_node_3->X(), p_node_3->Y(), p_node_3->Z());
            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(12, p_node_4->X(), p_node_4->Y(), p_node_4->Z());

            NodeType::Pointer p_node0_5 = model_part.CreateNewNode(13, p_node_5->X(), p_node_5->Y(), p_node_5->Z());
            NodeType::Pointer p_node0_6 = model_part.CreateNewNode(14, p_node_6->X(), p_node_6->Y(), p_node_6->Z());
            NodeType::Pointer p_node0_7 = model_part.CreateNewNode(15, p_node_7->X(), p_node_7->Y(), p_node_7->Z());
            NodeType::Pointer p_node0_8 = model_part.CreateNewNode(16, p_node_8->X(), p_node_8->Y(), p_node_8->Z());

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (4);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            condition_nodes_0[3] = p_node_4;
            Quadrilateral3D4 <NodeType> quadrilateral_0( PointerVector<NodeType>{condition_nodes_0} );

            std::vector<NodeType::Pointer> condition_nodes0_0 (4);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            condition_nodes0_0[2] = p_node0_3;
            condition_nodes0_0[3] = p_node0_4;
            Quadrilateral3D4 <NodeType> quadrilateral0_0( PointerVector<NodeType>{condition_nodes0_0} );

            const array_1d<double, 3>& normal_0 = quadrilateral_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", 1, quadrilateral_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", 3, quadrilateral0_0, p_cond_prop);
            p_cond0_0->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_0->GetGeometry().PointLocalCoordinates( point_local, p_cond0_0->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_0->GetGeometry().UnitNormal(point_local);
                p_cond0_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_0->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<NodeType::Pointer> condition_nodes_1 (4);
            condition_nodes_1[0] = p_node_5;
            condition_nodes_1[1] = p_node_6;
            condition_nodes_1[2] = p_node_7;
            condition_nodes_1[3] = p_node_8;
            Quadrilateral3D4 <NodeType> quadrilateral_1( PointerVector<NodeType>{condition_nodes_1} );

            std::vector<NodeType::Pointer> condition_nodes0_1 (4);
            condition_nodes0_1[0] = p_node0_5;
            condition_nodes0_1[1] = p_node0_6;
            condition_nodes0_1[2] = p_node0_7;
            condition_nodes0_1[3] = p_node0_8;
            Quadrilateral3D4 <NodeType> quadrilateral0_1( PointerVector<NodeType>{condition_nodes0_1} );

            const array_1d<double, 3>& normal_1 = quadrilateral_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", 2, quadrilateral_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D4N", 4, quadrilateral0_1, p_cond_prop);
            p_cond0_1->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            for (IndexType i_node = 0; i_node < p_cond0_0->GetGeometry().size(); ++i_node) {
                GeometryType::CoordinatesArrayType point_local;
                p_cond0_1->GetGeometry().PointLocalCoordinates( point_local, p_cond0_1->GetGeometry()[i_node].Coordinates( ) ) ;
                const array_1d<double, 3>& node_normal = p_cond0_1->GetGeometry().UnitNormal(point_local);
                p_cond0_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
                p_cond_1->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL) = node_normal;
            }

            std::vector<IndexType> nodes_perturbed(1, 5);
            std::vector<double> coeff_perturbation(1, -5.0e-2);
            TestDerivatives<3, 4>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, nodes_perturbed, 2, coeff_perturbation, 6, CHECK_NORMAL, LEVEL_QUADRATIC_CONVERGENCE);
        }
    } // namespace Testing
}  // namespace Kratos.
