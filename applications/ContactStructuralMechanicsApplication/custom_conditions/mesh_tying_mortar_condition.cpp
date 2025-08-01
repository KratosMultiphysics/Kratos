// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/geometrical_projection_utilities.h"
#include "custom_conditions/mesh_tying_mortar_condition.h"

namespace Kratos
{
/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
Condition::Pointer MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive< MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster> >( NewId, this->GetParentGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
Condition::Pointer MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive< MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster> >( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
Condition::Pointer MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties,
    GeometryType::Pointer pMasterGeom) const
{
    return Kratos::make_intrusive< MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster> >( NewId, pGeom, pProperties, pMasterGeom);
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::~MeshTyingMortarCondition( )
= default;

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Initialize BaseType
    BaseType::Initialize(rCurrentProcessInfo);

    // We get the unkown variable
    const std::string r_variable_name = GetProperties().Has(TYING_VARIABLE) ? GetProperties().GetValue(TYING_VARIABLE) : "DISPLACEMENT";
    mpDoFVariables.clear();
    mpLMVariables.clear();
    if (KratosComponents<Variable<double>>::Has(r_variable_name)) {
        mpDoFVariables.push_back(&KratosComponents<Variable<double>>::Get(r_variable_name));
        mpLMVariables.push_back(&SCALAR_LAGRANGE_MULTIPLIER);
    } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(r_variable_name)) {
        mpDoFVariables.reserve(TDim);
        mpLMVariables.reserve(TDim);
        mpDoFVariables.push_back(&KratosComponents<Variable<double>>::Get(r_variable_name + "_X"));
        mpLMVariables.push_back(&VECTOR_LAGRANGE_MULTIPLIER_X);
        mpDoFVariables.push_back(&KratosComponents<Variable<double>>::Get(r_variable_name + "_Y"));
        mpLMVariables.push_back(&VECTOR_LAGRANGE_MULTIPLIER_Y);
        // In case of 3D, we add the Z variable
        if constexpr (TDim == 3) {
            mpDoFVariables.push_back(&KratosComponents<Variable<double>>::Get(r_variable_name + "_Z"));
            mpLMVariables.push_back(&VECTOR_LAGRANGE_MULTIPLIER_Z);
        }
    } else {
        KRATOS_ERROR << "Compatible variables are: double or array_1d<double, 3> " << std::endl;
    }

    // We define the integration method
    const auto& r_properties = this->GetProperties();
    const IndexType integration_order = r_properties.Has(INTEGRATION_ORDER_CONTACT) ? r_properties.GetValue(INTEGRATION_ORDER_CONTACT) : 2;

    // The slave geometry
    GeometryType& r_slave_geometry = this->GetParentGeometry();
    const array_1d<double, 3>& r_normal_slave = this->GetValue(NORMAL);

    // Create and initialize condition variables:
    GeneralVariables rVariables;

    // Create Ae matrix
    MatrixDualLM Ae;

    // The master geometry
    GeometryType& r_master_geometry = this->GetPairedGeometry();
    const array_1d<double, 3>& r_normal_master = this->GetPairedNormal();

    // Initialize general variables for the current master element
    rVariables.Initialize();

    // Initialize the mortar operators
    mMortarConditionMatrices.Initialize();

    // We call the exact integration utility
    const double distance_threshold = 1.0e24;
    const double zero_tolerance_factor = 1.0e0;
    const bool consider_tessellation = r_properties.Has(CONSIDER_TESSELLATION) ? r_properties[CONSIDER_TESSELLATION] : false;
    IntegrationUtility integration_utility = IntegrationUtility (integration_order, distance_threshold, 0, zero_tolerance_factor, consider_tessellation);

    // Reading integration points
    ConditionArrayListType conditions_points_slave;
    const bool is_inside = integration_utility.GetExactIntegration(r_slave_geometry, r_normal_slave, r_master_geometry, r_normal_master, conditions_points_slave);

    double integration_area;
    integration_utility.GetTotalArea(r_slave_geometry, conditions_points_slave, integration_area);

    if ((is_inside == true) && ((integration_area/r_slave_geometry.Area()) > 1.0e-3 * r_slave_geometry.Area())) {
        IntegrationMethod this_integration_method = GetIntegrationMethod();

        // Initialize general variables for the current master element
        rVariables.Initialize();

        // Initialize the mortar operators
        mMortarConditionMatrices.Initialize();

        const bool dual_LM = CalculateAe(r_normal_master, Ae, rVariables, conditions_points_slave, this_integration_method);

        for (IndexType i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom) {
            PointerVector< PointType > points_array(TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
            for (IndexType i_node = 0; i_node < TDim; ++i_node) {
                PointType global_point;
                r_slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                points_array(i_node) = Kratos::make_shared<PointType>(PointType(global_point));
            }

            DecompositionType decomp_geom( points_array );

            bool bad_shape;
            if constexpr (TDim == 2) {
                bad_shape = MortarUtilities::LengthCheck(decomp_geom, r_slave_geometry.Length() * CheckThresholdCoefficient);
            } else { 
                bad_shape = MortarUtilities::HeronCheck(decomp_geom);
            }

            if (!bad_shape) {
                const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( this_integration_method );

                // Integrating the mortar operators
                for ( IndexType point_number = 0; point_number < integration_points_slave.size(); ++point_number ) {
                    // We compute the local coordinates
                    const auto local_point_decomp = PointType{integration_points_slave[point_number].Coordinates()};
                    PointType local_point_parent;
                    PointType gp_global;
                    decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                    r_slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);

                    // Calculate the kinematic variables
                    this->CalculateKinematics( rVariables, Ae, r_normal_master, local_point_decomp, local_point_parent, decomp_geom, dual_LM);

                    const double integration_weight = integration_points_slave[point_number].Weight();

                    mMortarConditionMatrices.CalculateMortarOperators(rVariables, integration_weight);
                }
            }
        }
    } else { // We deactivate
        this->Set(ACTIVE, false);
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::InitializeSolutionStep(rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::InitializeNonLinearIteration(rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::FinalizeSolutionStep(rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::FinalizeNonLinearIteration(rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // Compute the matrix size
    const SizeType matrix_size = mpDoFVariables.size() * (2 * TNumNodes + TNumNodesMaster);

    // Resizing as needed the LHS
    if ( rLeftHandSideMatrix.size1() != matrix_size || rLeftHandSideMatrix.size2() != matrix_size ) {
        rLeftHandSideMatrix.resize( matrix_size, matrix_size, false );
    }

    // Resizing as needed the RHS
    if ( rRightHandSideVector.size() != matrix_size ) {
        rRightHandSideVector.resize( matrix_size, false );
    }

    // Calculate condition system
    CalculateConditionSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo );

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // Compute the matrix size
    const SizeType matrix_size = mpDoFVariables.size() * (2 * TNumNodes + TNumNodesMaster);

    // Resizing as needed the LHS
    if ( rLeftHandSideMatrix.size1() != matrix_size || rLeftHandSideMatrix.size2() != matrix_size ) {
        rLeftHandSideMatrix.resize( matrix_size, matrix_size, false );
    }

    // Creating an auxiliary vector
    VectorType aux_right_hand_side_vector = Vector();

    // Calculate condition system
    CalculateConditionSystem(rLeftHandSideMatrix, aux_right_hand_side_vector, rCurrentProcessInfo, true, false );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // Creating an auxiliary matrix
    MatrixType aux_left_hand_side_matrix = Matrix();

    // Compute the matrix size
    const SizeType matrix_size = mpDoFVariables.size() * (2 * TNumNodes + TNumNodesMaster);

    // Resizing as needed the RHS
    if ( rRightHandSideVector.size() != matrix_size ) {
        rRightHandSideVector.resize( matrix_size, false );
    }

    // Calculate condition system
    CalculateConditionSystem(aux_left_hand_side_matrix, rRightHandSideVector, rCurrentProcessInfo, false );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::CalculateMassMatrix(
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

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::CalculateDampingMatrix(
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

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MeshTyingMortarCondition<TDim, TNumNodes, TNumNodesMaster>::CalculateConditionSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool ComputeLHS,
    const bool ComputeRHS
    )
{
    KRATOS_TRY;

    // Create the current DoF data
    DofData dof_data(mpDoFVariables.size());

    // Initialize the DoF data
    this->InitializeDofData(dof_data);

    // Update slave element info
    dof_data.UpdateMasterPair(this->GetPairedGeometry(), mpDoFVariables);

    // Assemble of the matrix is required
    if ( ComputeLHS ) {
        // Calculate the local contribution
        this->CalculateLocalLHS(rLeftHandSideMatrix, mMortarConditionMatrices, dof_data, rCurrentProcessInfo);
    }

    // Assemble of the vector is required
    if ( ComputeRHS) {
        // Calculate the local contribution
        this->CalculateLocalRHS( rRightHandSideVector, mMortarConditionMatrices, dof_data, rCurrentProcessInfo);
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
bool MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::CalculateAe(
    const array_1d<double, 3>& rNormalMaster,
    MatrixDualLM& rAe,
    GeneralVariables& rVariables,
    const ConditionArrayListType& rConditionsPointsSlave,
    const IntegrationMethod ThisIntegrationMethod
    )
{
    // We initialize the Ae components
    AeData rAeData;
    rAeData.Initialize();

    rAe = ZeroMatrix(TNumNodes, TNumNodes);

    // The slave geometry
    GeometryType& r_slave_geometry = this->GetParentGeometry();

    // Initialize general variables for the current master element
    rVariables.Initialize();

    // Calculating the proportion between the integrated area and segment area
    for (IndexType i_geom = 0; i_geom < rConditionsPointsSlave.size(); ++i_geom) {
        PointerVector< PointType > points_array(TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
        for (IndexType i_node = 0; i_node < TDim; ++i_node) {
            PointType global_point;
            r_slave_geometry.GlobalCoordinates(global_point, rConditionsPointsSlave[i_geom][i_node]);
            points_array(i_node) = Kratos::make_shared<PointType>(PointType(global_point));
        }

        DecompositionType decomp_geom( points_array );

        bool bad_shape;
        if constexpr (TDim == 2) {
            bad_shape = MortarUtilities::LengthCheck(decomp_geom, r_slave_geometry.Length() * CheckThresholdCoefficient);
        } else { 
            bad_shape = MortarUtilities::HeronCheck(decomp_geom);
        }

        if (!bad_shape) {
            const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( ThisIntegrationMethod );

            // Integrating the mortar operators
            for ( IndexType point_number = 0; point_number < integration_points_slave.size(); ++point_number ) {
                // We compute the local coordinates
                const auto local_point_decomp = PointType{integration_points_slave[point_number].Coordinates()};
                PointType local_point_parent;
                PointType gp_global;
                decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                r_slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);

                // Calculate the kinematic variables
                // We compute the current configuration
                decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                r_slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);

                this->CalculateKinematics( rVariables, rAe, rNormalMaster, local_point_decomp, local_point_parent, decomp_geom, false);

                // Integrate
                const double integration_weight = integration_points_slave[point_number].Weight();

                rAeData.CalculateAeComponents(rVariables, integration_weight);
            }
        }
    }

    return rAeData.CalculateAe(rAe);
}

/*********************************COMPUTE KINEMATICS*********************************/
/************************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::CalculateKinematics(
    GeneralVariables& rVariables,
    const MatrixDualLM& rAe,
    const array_1d<double, 3>& rNormalMaster,
    const PointType& rLocalPointDecomp,
    const PointType& rLocalPointParent,
    const GeometryPointType& rGeometryDecomp,
    const bool DualLM
    )
{
    /// SLAVE CONDITION ///
    /* SHAPE FUNCTIONS */
    GetParentGeometry().ShapeFunctionsValues( rVariables.NSlave, rLocalPointParent.Coordinates() );
    rVariables.PhiLagrangeMultipliers = (DualLM == true) ? prod(rAe, rVariables.NSlave) : rVariables.NSlave;

    /* CALCULATE JACOBIAN AND JACOBIAN DETERMINANT */
    rVariables.DetjSlave = rGeometryDecomp.DeterminantOfJacobian( rLocalPointDecomp );

    KRATOS_ERROR_IF(rVariables.DetjSlave < 0.0) << "WARNING:: CONDITION ID: " << this->Id() << " INVERTED. DETJ: " << rVariables.DetjSlave << std::endl;

    /// MASTER CONDITION ///
    this->MasterShapeFunctionValue( rVariables, rNormalMaster, rLocalPointParent);
}

/***********************************************************************************/
/*************** METHODS TO CALCULATE THE CONTACT CONDITION MATRICES ***************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::MasterShapeFunctionValue(
    GeneralVariables& rVariables,
    const array_1d<double, 3>& rNormalMaster,
    const PointType& rLocalPoint
    )
{
    GeometryType& r_master_geometry = this->GetPairedGeometry();

    PointType projected_gp_global;
    const array_1d<double,3> gp_normal = MortarUtilities::GaussPointUnitNormal(rVariables.NSlave, GetParentGeometry());

    GeometryType::CoordinatesArrayType slave_gp_global;
    this->GetParentGeometry().GlobalCoordinates( slave_gp_global, rLocalPoint );
    GeometricalProjectionUtilities::FastProjectDirection( r_master_geometry, PointType{slave_gp_global}, projected_gp_global, rNormalMaster, -gp_normal ); // The opposite direction

    GeometryType::CoordinatesArrayType projected_gp_local;

    r_master_geometry.PointLocalCoordinates(projected_gp_local, projected_gp_global.Coordinates( ) ) ;

    // SHAPE FUNCTIONS
    r_master_geometry.ShapeFunctionsValues( rVariables.NMaster, projected_gp_local );
}

/***************************** BEGIN AD REPLACEMENT ********************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::CalculateLocalLHS(
    Matrix& rLocalLHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData& rDofData,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We get the mortar operators
    const auto& r_MOperator = rMortarConditionMatrices.MOperator;
    const auto& r_DOperator = rMortarConditionMatrices.DOperator;

    // Clear matrix
    rLocalLHS.clear();

    // Get the DoF size
    const SizeType dof_size = mpDoFVariables.size();

    // Get the scale factor
    const double scale_factor = rCurrentProcessInfo.Has(SCALE_FACTOR) ? rCurrentProcessInfo[SCALE_FACTOR] : rCurrentProcessInfo.Has(BUILD_SCALE_FACTOR) ? rCurrentProcessInfo[BUILD_SCALE_FACTOR] : 1.0;

    // Initial index 
    IndexType initial_row_index = 0;
    const IndexType initial_column_index = dof_size * (TNumNodes + TNumNodesMaster);

    // Iterate over the number of dofs on master side
    for (IndexType j = 0; j < TNumNodesMaster; ++j) {
        for (IndexType k = 0; k < TNumNodes; ++k) {
            const double value = - scale_factor * r_MOperator(k, j);
            for (IndexType i = 0; i < dof_size; ++i) {
                rLocalLHS(initial_row_index + j * dof_size + i, initial_column_index + k * dof_size + i) = value;
                rLocalLHS(initial_column_index + k * dof_size + i, initial_row_index + j * dof_size + i) = value;
            }
        }
    }
    
    // Update initial index
    initial_row_index = dof_size * TNumNodesMaster;

    // Iterate over the number of dofs on slave side
    for (IndexType j = 0; j < TNumNodes; ++j) {
        for (IndexType k = 0; k < TNumNodes; ++k) {
            const double value = scale_factor * r_DOperator(k, j);
            for (IndexType i = 0; i < dof_size; ++i) {
                rLocalLHS(initial_row_index + j * dof_size + i, initial_column_index + k * dof_size + i) = value;
                rLocalLHS(initial_column_index + k * dof_size + i, initial_row_index + j * dof_size + i) = value;
            }
        }
    }

    // // Debugging
    // LOG_MATRIX_PRETTY(rLocalLHS);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::CalculateLocalRHS(
    Vector& rLocalRHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DofData& rDofData,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // Initialize values
    const auto& r_u1 = rDofData.u1;
    const auto& r_u2 = rDofData.u2;
    const auto& r_lm = rDofData.LagrangeMultipliers;

    // Mortar operators
    const auto& r_MOperator = rMortarConditionMatrices.MOperator;
    const auto& r_DOperator = rMortarConditionMatrices.DOperator;

    // Get the DoF size
    const SizeType dof_size = mpDoFVariables.size();

    // Get the scale factor
    const double scale_factor = rCurrentProcessInfo.Has(SCALE_FACTOR) ? rCurrentProcessInfo[SCALE_FACTOR] : rCurrentProcessInfo.Has(BUILD_SCALE_FACTOR) ? rCurrentProcessInfo[BUILD_SCALE_FACTOR] : 1.0;

    // Initial index 
    IndexType initial_index = 0;

    // Master side
    const Matrix Mlm = scale_factor * prod(trans(r_MOperator), r_lm);
    for (IndexType i = 0; i < TNumNodesMaster; ++i) {
        for (IndexType j = 0; j < dof_size; ++j) {
            rLocalRHS[initial_index + i * dof_size + j] = Mlm(i, j);
        }
    }

    // Slave side
    initial_index = TNumNodesMaster * dof_size;
    const Matrix Dlm = scale_factor * prod(trans(r_DOperator), r_lm);
    for (IndexType i = 0; i < TNumNodes; ++i) {
        for (IndexType j = 0; j < dof_size; ++j) {
            rLocalRHS[initial_index + i * dof_size + j] = - Dlm(i, j);
        }
    }

    // LM slave side
    initial_index = (TNumNodes + TNumNodesMaster) * dof_size;
    const Matrix Du1Mu2 = scale_factor * (prod(r_DOperator, r_u1) - prod(r_MOperator, r_u2));
    for (IndexType i = 0; i < TNumNodes; ++i) {
        for (IndexType j = 0; j < dof_size; ++j) {
            rLocalRHS[initial_index + i * dof_size + j] = - Du1Mu2(i, j);
        }
    }

    // // Debugging
    // LOG_VECTOR_PRETTY(rLocalRHS);
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    // Compute the matrix size
    const SizeType matrix_size = mpDoFVariables.size() * (2 * TNumNodes + TNumNodesMaster);

    if (rResult.size() != matrix_size) {
        rResult.resize( matrix_size, false );
    }

    IndexType index = 0;

    /* ORDER - [ MASTER, SLAVE, LM ] */
    // Master Nodes DoF Equation IDs
    const GeometryType& r_current_master = this->GetPairedGeometry();
    for ( IndexType i_master = 0; i_master < TNumNodesMaster; ++i_master ) {
        const Node& r_master_node = r_current_master[i_master];
        for (auto& p_var : mpDoFVariables) {
            rResult[index++] = r_master_node.GetDof(*p_var).EquationId( );
        }
    }


    // Slave Nodes DoF Equation IDs
    const GeometryType& r_current_slave = this->GetParentGeometry();
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        const Node& r_slave_node = r_current_slave[i_slave];
        for (auto& p_var : mpDoFVariables) {
            rResult[index++] = r_slave_node.GetDof(*p_var).EquationId( );
        }
    }

    // Slave Nodes LM Equation IDs
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        const Node& r_slave_node = r_current_slave[i_slave];
        for (auto& p_var : mpLMVariables) {
            rResult[index++] = r_slave_node.GetDof(*p_var).EquationId( );
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MeshTyingMortarCondition<TDim, TNumNodes, TNumNodesMaster>::GetDofList(
    DofsVectorType& rConditionalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    // Compute the matrix size
    const SizeType matrix_size = mpDoFVariables.size() * (2 * TNumNodes + TNumNodesMaster);

    if (rConditionalDofList.size() != matrix_size) {
        rConditionalDofList.resize( matrix_size );
    }

    IndexType index = 0;

    /* ORDER - [ MASTER, SLAVE, LM ] */
    // Master Nodes DoF Equation IDs
    const GeometryType& r_current_master = this->GetPairedGeometry();
    for ( IndexType i_master = 0; i_master < TNumNodesMaster; ++i_master ) {
        const Node& r_master_node = r_current_master[i_master];
        for (auto& p_var : mpDoFVariables) {
            rConditionalDofList[index++] = r_master_node.pGetDof(*p_var);
        }
    }


    // Slave Nodes DoF Equation IDs
    const GeometryType& r_current_slave = this->GetParentGeometry();
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        const Node& r_slave_node = r_current_slave[i_slave];
        for (auto& p_var : mpDoFVariables) {
            rConditionalDofList[index++] = r_slave_node.pGetDof(*p_var);
        }
    }

    // Slave Nodes LM Equation IDs
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        const Node& r_slave_node = r_current_slave[i_slave];
        for (auto& p_var : mpLMVariables) {
            rConditionalDofList[index++] = r_slave_node.pGetDof(*p_var);
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType& r_integration_points = GetParentGeometry().IntegrationPoints();

    if ( rOutput.size() != r_integration_points.size() ) {
        rOutput.resize( r_integration_points.size() );
    }

    for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
        rOutput[point_number] = 0.0;
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector< array_1d<double, 3 > >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType& r_integration_points = GetParentGeometry().IntegrationPoints();

    if ( rOutput.size() != r_integration_points.size() ) {
        rOutput.resize( r_integration_points.size() );
    }

    for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
        rOutput[point_number] = ZeroVector(3);
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType& r_integration_points = GetParentGeometry().IntegrationPoints();

    if ( rOutput.size() != r_integration_points.size() ) {
        rOutput.resize( r_integration_points.size() );
    }

    for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
        rOutput[point_number] = ZeroVector(3);
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
int MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int ierr = Condition::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    KRATOS_ERROR_IF(this->GetParentGeometry().NumberOfGeometryParts() == 0) << "YOU HAVE NOT INITIALIZED THE PAIR GEOMETRY IN THE MeshTyingMortarCondition" << std::endl;

    return ierr;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, PairedCondition );
    rSerializer.save("MortarConditionMatrices", mMortarConditionMatrices);
    rSerializer.save("DoFVariables", mpDoFVariables);
    rSerializer.save("LMVariables", mpLMVariables);
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MeshTyingMortarCondition<TDim,TNumNodes, TNumNodesMaster>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, PairedCondition );
    rSerializer.load("MortarConditionMatrices", mMortarConditionMatrices);
    rSerializer.load("DoFVariables", mpDoFVariables);
    rSerializer.load("LMVariables", mpLMVariables);
}

/***********************************************************************************/
/***********************************************************************************/

template class MeshTyingMortarCondition<2, 2, 2>; // 2D Line/Line
template class MeshTyingMortarCondition<3, 3, 3>; // 3D Triangle/Triangle
template class MeshTyingMortarCondition<3, 4, 4>; // 3D Quadrilateral/Quadrilateral
template class MeshTyingMortarCondition<3, 3, 4>; // 3D Triangle/Quadrilateral
template class MeshTyingMortarCondition<3, 4, 3>; // 3D Quadrilateral/Triangle

} // Namespace Kratos
