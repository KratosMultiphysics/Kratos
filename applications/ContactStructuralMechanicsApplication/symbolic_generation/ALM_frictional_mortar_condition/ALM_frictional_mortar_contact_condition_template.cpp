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
#include "custom_conditions/ALM_frictional_mortar_contact_condition.h"

namespace Kratos 
{
/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
Condition::Pointer AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::Create( 
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesPointerType pProperties ) const
{
    return Kratos::make_shared< AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation > >( NewId, this->GetGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
Condition::Pointer AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::Create(
    IndexType NewId,
    GeometryPointerType pGeom,
    PropertiesPointerType pProperties) const
{
    return Kratos::make_shared< AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation> >( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
Condition::Pointer AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::Create(
    IndexType NewId,
    GeometryPointerType pGeom,
    PropertiesType::Pointer pProperties,
    GeometryType::Pointer pMasterGeom ) const
{
    return Kratos::make_shared< AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation> >( NewId, pGeom, pProperties, pMasterGeom );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::~AugmentedLagrangianMethodFrictionalMortarContactCondition( )
= default;

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::Initialize( )
{
    KRATOS_TRY;

    BaseType::Initialize();

    // We initailize the previous mortar operators
    mPreviousMortarOperators.Initialize();

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
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

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::AddExplicitContribution(ProcessInfo& rCurrentProcessInfo)
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
            const array_1d<double, TDim> aux_array = row(D_x1_M_x2, i_node);

            double& weighted_gap = slave_geometry[i_node].FastGetSolutionStepValue(WEIGHTED_GAP);

            #pragma omp atomic
            weighted_gap += inner_prod(aux_array, - subrange(normal, 0, TDim));
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

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& CurrentProcessInfo 
    )
{
    KRATOS_TRY;   
    
    if (rResult.size() != MatrixSize)
        rResult.resize( MatrixSize, false );
    
    IndexType index = 0;
    
    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    GeometryType& current_master = this->GetPairedGeometry();;
    
    // Master Nodes Displacement Equation IDs
    for ( IndexType i_master = 0; i_master < TNumNodes; ++i_master ) { // NOTE: Assuming same number of nodes for master and slave
        NodeType& master_node = current_master[i_master];
        rResult[index++] = master_node.GetDof( DISPLACEMENT_X ).EquationId( );
        rResult[index++] = master_node.GetDof( DISPLACEMENT_Y ).EquationId( );
        if (TDim == 3) rResult[index++] = master_node.GetDof( DISPLACEMENT_Z ).EquationId( );
    }

    // Slave Nodes Displacement Equation IDs
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        NodeType& slave_node = this->GetGeometry()[ i_slave ];
        rResult[index++] = slave_node.GetDof( DISPLACEMENT_X ).EquationId( );
        rResult[index++] = slave_node.GetDof( DISPLACEMENT_Y ).EquationId( );
        if (TDim == 3) rResult[index++] = slave_node.GetDof( DISPLACEMENT_Z ).EquationId( );
    }

    // Slave Nodes  Lambda Equation IDs
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        NodeType& slave_node = this->GetGeometry()[ i_slave ];
        rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_X ).EquationId( );
        rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Y ).EquationId( );
        if (TDim == 3) rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Z ).EquationId( );
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::GetDofList(
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo 
)
{
    KRATOS_TRY;
    
    if (rConditionalDofList.size() != MatrixSize)
        rConditionalDofList.resize( MatrixSize );
    
    IndexType index = 0;
    
    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    GeometryType& current_master = this->GetPairedGeometry();;

    // Master Nodes Displacement Equation IDs
    for ( IndexType i_master = 0; i_master < TNumNodes; ++i_master ){ // NOTE: Assuming same number of nodes for master and slave
        NodeType& master_node = current_master[i_master];
        rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_X );
        rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_Y );
        if (TDim == 3) rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_Z );
    }

    // Slave Nodes Displacement Equation IDs
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        NodeType& slave_node = this->GetGeometry()[ i_slave ];
        rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_X );
        rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_Y );
        if (TDim == 3) rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_Z );
    }

    // Slave Nodes Lambda Equation IDs
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        NodeType& slave_node = this->GetGeometry()[ i_slave ];
        rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_X );
        rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Y );
        if (TDim == 3) rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Z );
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
int AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int ierr = BaseType::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(NORMAL)
    KRATOS_CHECK_VARIABLE_KEY(VECTOR_LAGRANGE_MULTIPLIER)
    KRATOS_CHECK_VARIABLE_KEY(WEIGHTED_SLIP)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( IndexType i = 0; i < TNumNodes; ++i ) {
        Node<3> &rnode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NORMAL,rnode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VECTOR_LAGRANGE_MULTIPLIER,rnode)

        KRATOS_CHECK_DOF_IN_NODE(VECTOR_LAGRANGE_MULTIPLIER_X, rnode)
        KRATOS_CHECK_DOF_IN_NODE(VECTOR_LAGRANGE_MULTIPLIER_Y, rnode)
        KRATOS_CHECK_DOF_IN_NODE(VECTOR_LAGRANGE_MULTIPLIER_Z, rnode)
    }

    return ierr;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template class AugmentedLagrangianMethodFrictionalMortarContactCondition<2, 2, false>;
template class AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 3, false>;
template class AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 4, false>;
template class AugmentedLagrangianMethodFrictionalMortarContactCondition<2, 2, true>;
template class AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 3, true>;
template class AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 4, true>;

} // Namespace Kratos
