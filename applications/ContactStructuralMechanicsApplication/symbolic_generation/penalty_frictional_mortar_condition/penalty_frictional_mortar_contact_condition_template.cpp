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
/* Mortar includes */
#include "custom_conditions/penalty_frictional_mortar_contact_condition.h"

namespace Kratos
{
/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
Condition::Pointer PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesPointerType pProperties ) const
{
    return Kratos::make_shared< PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster > >( NewId, this->GetGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
Condition::Pointer PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::Create(
    IndexType NewId,
    GeometryPointerType pGeom,
    PropertiesPointerType pProperties) const
{
    return Kratos::make_shared< PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster> >( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
Condition::Pointer PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::Create(
    IndexType NewId,
    GeometryPointerType pGeom,
    PropertiesType::Pointer pProperties,
    GeometryType::Pointer pMasterGeom ) const
{
    return Kratos::make_shared< PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster> >( NewId, pGeom, pProperties, pMasterGeom );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::~PenaltyMethodFrictionalMortarContactCondition( )
= default;

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
void PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::Initialize( )
{
    KRATOS_TRY;

    BaseType::Initialize();

    // We initailize the previous mortar operators
    mPreviousMortarOperators.Initialize();

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
void PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    BaseType::FinalizeSolutionStep(rCurrentProcessInfo);

    // We "save" the mortar operator for the next step
    // The slave geometry
    GeometryType& slave_geometry = this->GetGeometry();
    const array_1d<double, 3>& normal_slave = this->GetValue(NORMAL);

    // Create and initialize condition variables
    GeneralVariables rVariables;

    // Create the current contact data
    DerivativeDataType rDerivativeData;
    rDerivativeData.Initialize(slave_geometry, rCurrentProcessInfo);

    // We call the exact integration utility
    const double distance_threshold = rCurrentProcessInfo[DISTANCE_THRESHOLD];
    IntegrationUtility integration_utility = IntegrationUtility (BaseType::mIntegrationOrder, distance_threshold);

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
        IntegrationMethod this_integration_method = this->GetIntegrationMethod();

        // Initialize general variables for the current master element
        rVariables.Initialize();

        // Update slave element info
        rDerivativeData.UpdateMasterPair(master_geometry, rCurrentProcessInfo);

        // Initialize the mortar operators
        mPreviousMortarOperators.Initialize();

        const bool dual_LM = DerivativesUtilitiesType::CalculateAeAndDeltaAe(slave_geometry, normal_slave, master_geometry, rDerivativeData, rVariables, consider_normal_variation, conditions_points_slave, this_integration_method, this->GetAxisymmetricCoefficient(rVariables));

        for (IndexType i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom) {
            std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
            array_1d<BelongType, TDim> belong_array;
            for (IndexType i_node = 0; i_node < TDim; ++i_node) {
                PointType global_point;
                slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                points_array[i_node] = Kratos::make_shared<PointType>(PointType(global_point));
                belong_array[i_node] = conditions_points_slave[i_geom][i_node].GetBelong();
            }

            DecompositionType decomp_geom( points_array );

            const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, slave_geometry.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);

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

                    const double integration_weight = integration_points_slave[point_number].Weight() * this->GetAxisymmetricCoefficient(rVariables);

                    mPreviousMortarOperators.CalculateMortarOperators(rVariables, integration_weight);
                }
            }
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
void PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::AddExplicitContribution(ProcessInfo& rCurrentProcessInfo)
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

    // Create the mortar operators
    MortarConditionMatrices rThisMortarConditionMatrices;

    // We call the exact integration utility
    const double distance_threshold = rCurrentProcessInfo[DISTANCE_THRESHOLD];
    IntegrationUtility integration_utility = IntegrationUtility (BaseType::mIntegrationOrder, distance_threshold);

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
        IntegrationMethod this_integration_method = this->GetIntegrationMethod();

        // Initialize general variables for the current master element
        rVariables.Initialize();

        // Update slave element info
        rDerivativeData.UpdateMasterPair(master_geometry,rCurrentProcessInfo);

        // Initialize the mortar operators
        rThisMortarConditionMatrices.Initialize();

        const bool dual_LM = DerivativesUtilitiesType::CalculateAeAndDeltaAe(slave_geometry, normal_slave, master_geometry, rDerivativeData, rVariables, consider_normal_variation, conditions_points_slave, this_integration_method, this->GetAxisymmetricCoefficient(rVariables));

        for (IndexType i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom) {
            std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
            array_1d<BelongType, TDim> belong_array;
            for (IndexType i_node = 0; i_node < TDim; ++i_node) {
                PointType global_point;
                slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                points_array[i_node] = Kratos::make_shared<PointType>(PointType(global_point));
                belong_array[i_node] = conditions_points_slave[i_geom][i_node].GetBelong();
            }

            DecompositionType decomp_geom( points_array );

            const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, slave_geometry.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);

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

                    const double integration_weight = integration_points_slave[point_number].Weight() * this->GetAxisymmetricCoefficient(rVariables);

                    rThisMortarConditionMatrices.CalculateMortarOperators(rVariables, integration_weight);
                }
            }
        }

        // Setting the weighted gap
        // Mortar condition matrices - DOperator and MOperator
        const BoundedMatrix<double, TNumNodes, TNumNodes>& DOperator = rThisMortarConditionMatrices.DOperator;
        const BoundedMatrix<double, TNumNodes, TNumNodes>& MOperator = rThisMortarConditionMatrices.MOperator;

        // Current coordinates
        const BoundedMatrix<double, TNumNodes, TDim> x1 = MortarUtilities::GetCoordinates<TDim,TNumNodes>(slave_geometry);
        const BoundedMatrix<double, TNumNodes, TDim> x2 = MortarUtilities::GetCoordinates<TDim,TNumNodes>(master_geometry);

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

        // Setting the weighted slip
        // The increment of time
        const double delta_time = rCurrentProcessInfo[DELTA_TIME];

        // Delta mortar condition matrices - DOperator and MOperator
        const BoundedMatrix<double, TNumNodes, TNumNodes> DeltaDOperator = DOperator - mPreviousMortarOperators.DOperator;
        const BoundedMatrix<double, TNumNodes, TNumNodes> DeltaMOperator = MOperator - mPreviousMortarOperators.MOperator;

        // Old coordinates
        const BoundedMatrix<double, TNumNodes, TDim> x1_old = MortarUtilities::GetCoordinates<TDim,TNumNodes>(slave_geometry, false, 1);
        const BoundedMatrix<double, TNumNodes, TDim> x2_old = MortarUtilities::GetCoordinates<TDim,TNumNodes>(master_geometry, false, 1);

        const BoundedMatrix<double, TNumNodes, TDim> D_x1_old_M_x2_old = prod(DOperator, x1_old) - prod(MOperator, x2_old);

        const BoundedMatrix<double, TNumNodes, TDim> delta_D_x1_M_x2 = prod(DeltaDOperator, x1) - prod(DeltaMOperator, x2);

        // The tangent matrix
        const BoundedMatrix<double, TNumNodes, TDim> tangent_slave = MortarUtilities::ComputeTangentMatrix<TNumNodes, TDim>(slave_geometry);

        // The estimation of the slip time derivative
        const BoundedMatrix<double, TNumNodes, TDim> slip_time_derivative = (D_x1_old_M_x2_old - D_x1_M_x2)/delta_time - delta_D_x1_M_x2/delta_time;

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            // We compute the tangent
            const array_1d<double, TDim>& tangent_node = row(tangent_slave, i_node);
            const array_1d<double, TDim>& slip_time_derivative_node = row(slip_time_derivative, i_node);
            const double slip_node = delta_time * inner_prod(tangent_node, slip_time_derivative_node);

            // The weighted slip
            array_1d<double, 3>& weighted_slip = slave_geometry[i_node].FastGetSolutionStepValue(WEIGHTED_SLIP);

            for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
                #pragma omp atomic
                weighted_slip[i_dim] += slip_node * tangent_node[i_dim];
            }
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

template< SizeType TDim, SizeType TNumNodes, bool TNormalVariation, SizeType TNumNodesMaster >
void PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation,TNumNodesMaster>::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    Variable<double>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // TODO: Add something if necessary

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, bool TNormalVariation, SizeType TNumNodesMaster >
void PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation,TNumNodesMaster>::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    Variable<array_1d<double, 3>>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // Computing the force residual
    if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL) {
        if (this->Is(ACTIVE)) {
            // Getting geometries
            GeometryType& r_slave_geometry = this->GetGeometry();
            GeometryType& r_master_geometry = this->GetPairedGeometry();

            for ( IndexType i_master = 0; i_master < TNumNodesMaster; ++i_master ) {
                NodeType& r_master_node = r_master_geometry[i_master];
                const IndexType index = TDim * i_master;

                array_1d<double, 3>& r_force_residual = r_master_node.FastGetSolutionStepValue(FORCE_RESIDUAL);

                for (IndexType j = 0; j < TDim; ++j) {
                    #pragma omp atomic
                    r_force_residual[j] += rRHSVector[index + j];
                }
            }
            for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
                NodeType& r_slave_node = r_slave_geometry[i_slave];
                const IndexType index = TDim * (TNumNodesMaster + i_slave);

                array_1d<double, 3>& r_force_residual = r_slave_node.FastGetSolutionStepValue(FORCE_RESIDUAL);

                for (IndexType j = 0; j < TDim; ++j) {
                    #pragma omp atomic
                    r_force_residual[j] += rRHSVector[index + j];
                }
            }
        }
    }

    KRATOS_CATCH("")
}

/***************************** BEGIN AD REPLACEMENT ********************************/
/***********************************************************************************/

// replace_lhs

/****************************** END AD REPLACEMENT *********************************/
/***********************************************************************************/

/***************************** BEGIN AD REPLACEMENT ********************************/
/***********************************************************************************/

// replace_rhs

/****************************** END AD REPLACEMENT *********************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
void PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& CurrentProcessInfo
    )
{
    KRATOS_TRY;

    if (rResult.size() != MatrixSize)
        rResult.resize( MatrixSize, false );

    IndexType index = 0;

    /* ORDER - [ MASTER, SLAVE ] */
    GeometryType& r_slave_geometry = this->GetGeometry();
    GeometryType& r_master_geometry = this->GetPairedGeometry();

    // Master Nodes Displacement Equation IDs
    for ( IndexType i_master = 0; i_master < TNumNodes; ++i_master ) { // NOTE: Assuming same number of nodes for master and slave
        NodeType& r_master_node = r_master_geometry[i_master];
        rResult[index++] = r_master_node.GetDof( DISPLACEMENT_X ).EquationId( );
        rResult[index++] = r_master_node.GetDof( DISPLACEMENT_Y ).EquationId( );
        if (TDim == 3) rResult[index++] = r_master_node.GetDof( DISPLACEMENT_Z ).EquationId( );
    }

    // Slave Nodes Displacement Equation IDs
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        NodeType& r_slave_node = r_slave_geometry[ i_slave ];
        rResult[index++] = r_slave_node.GetDof( DISPLACEMENT_X ).EquationId( );
        rResult[index++] = r_slave_node.GetDof( DISPLACEMENT_Y ).EquationId( );
        if (TDim == 3) rResult[index++] = r_slave_node.GetDof( DISPLACEMENT_Z ).EquationId( );
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
void PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::GetDofList(
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY;

    if (rConditionalDofList.size() != MatrixSize)
        rConditionalDofList.resize( MatrixSize );

    IndexType index = 0;

    /* ORDER - [ MASTER, SLAVE ] */
    GeometryType& r_slave_geometry = this->GetGeometry();
    GeometryType& r_master_geometry = this->GetPairedGeometry();

    // Master Nodes Displacement Equation IDs
    for ( IndexType i_master = 0; i_master < TNumNodes; ++i_master ){ // NOTE: Assuming same number of nodes for master and slave
        NodeType& r_master_node = r_master_geometry[i_master];
        rConditionalDofList[index++] = r_master_node.pGetDof( DISPLACEMENT_X );
        rConditionalDofList[index++] = r_master_node.pGetDof( DISPLACEMENT_Y );
        if (TDim == 3) rConditionalDofList[index++] = r_master_node.pGetDof( DISPLACEMENT_Z );
    }

    // Slave Nodes Displacement Equation IDs
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        NodeType& r_slave_node = r_slave_geometry[ i_slave ];
        rConditionalDofList[index++] = r_slave_node.pGetDof( DISPLACEMENT_X );
        rConditionalDofList[index++] = r_slave_node.pGetDof( DISPLACEMENT_Y );
        if (TDim == 3) rConditionalDofList[index++] = r_slave_node.pGetDof( DISPLACEMENT_Z );
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
int PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int ierr = BaseType::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(NORMAL)
    KRATOS_CHECK_VARIABLE_KEY(WEIGHTED_SLIP)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    GeometryType& r_slave_geometry = this->GetGeometry();
    for ( IndexType i = 0; i < TNumNodes; ++i ) {
        NodeType& r_node = r_slave_geometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NORMAL, r_node)
    }

    return ierr;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template class PenaltyMethodFrictionalMortarContactCondition<2, 2, false, 2>;
template class PenaltyMethodFrictionalMortarContactCondition<3, 3, false, 3>;
template class PenaltyMethodFrictionalMortarContactCondition<3, 4, false, 4>;
template class PenaltyMethodFrictionalMortarContactCondition<3, 3, false, 4>;
template class PenaltyMethodFrictionalMortarContactCondition<3, 4, false, 3>;
template class PenaltyMethodFrictionalMortarContactCondition<2, 2, true,  2>;
template class PenaltyMethodFrictionalMortarContactCondition<3, 3, true,  3>;
template class PenaltyMethodFrictionalMortarContactCondition<3, 4, true,  4>;
template class PenaltyMethodFrictionalMortarContactCondition<3, 3, true,  4>;
template class PenaltyMethodFrictionalMortarContactCondition<3, 4, true,  3>;

} // Namespace Kratos
