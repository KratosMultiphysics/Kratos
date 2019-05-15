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
#include "custom_utilities/mortar_explicit_contribution_utilities.h"
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
void PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    BaseType::InitializeSolutionStep(rCurrentProcessInfo);

    // Compute the previous mortar operators
    if (!mPreviousMortarOperatorsInitialized) {
        ComputePreviousMortarOperators(rCurrentProcessInfo);
        mPreviousMortarOperatorsInitialized = true;
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
void PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    BaseType::FinalizeSolutionStep(rCurrentProcessInfo);

    // Compute the previous mortar operators
    ComputePreviousMortarOperators(rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
void PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::AddExplicitContribution(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    MortarExplicitContributionUtilities<TDim, TNumNodes, FrictionalCase::FRICTIONAL_PENALTY, TNormalVariation, TNumNodesMaster>::AddExplicitContributionOfMortarFrictionalCondition(this, rCurrentProcessInfo, mPreviousMortarOperators, BaseType::mIntegrationOrder, false, false);

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
            GeometryType& r_master_geometry = this->GetGeometry();

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


/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, bool TNormalVariation, SizeType TNumNodesMaster >
void PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation,TNumNodesMaster>::ComputePreviousMortarOperators( ProcessInfo& rCurrentProcessInfo)
{
    // We "save" the mortar operator for the next step
    // The slave geometry
    GeometryType& r_slave_geometry = this->GetGeometry();
    const array_1d<double, 3>& r_normal_slave = this->GetValue(NORMAL);

    // Create and initialize condition variables
    GeneralVariables rVariables;

    // Create the current contact data
    DerivativeDataType rDerivativeData;
    rDerivativeData.Initialize(r_slave_geometry, rCurrentProcessInfo);

    // We call the exact integration utility
    const double distance_threshold = rCurrentProcessInfo[DISTANCE_THRESHOLD];
    IntegrationUtility integration_utility = IntegrationUtility (BaseType::mIntegrationOrder, distance_threshold);

    // If we consider the normal variation
    const NormalDerivativesComputation consider_normal_variation = static_cast<NormalDerivativesComputation>(rCurrentProcessInfo[CONSIDER_NORMAL_VARIATION]);

    // The master geometry
    GeometryType& r_master_geometry = this->GetPairedGeometry();

    // The normal of the master condition
    const array_1d<double, 3>& r_normal_master = this->GetValue(PAIRED_NORMAL);

    // Reading integration points
    ConditionArrayListType conditions_points_slave;
    const bool is_inside = integration_utility.GetExactIntegration(r_slave_geometry, r_normal_slave, r_master_geometry, r_normal_master, conditions_points_slave);

    double integration_area;
    integration_utility.GetTotalArea(r_slave_geometry, conditions_points_slave, integration_area);

    const double geometry_area = r_slave_geometry.Area();
    if (is_inside && ((integration_area/geometry_area) > 1.0e-5)) {
        IntegrationMethod this_integration_method = this->GetIntegrationMethod();

        // Initialize general variables for the current master element
        rVariables.Initialize();

        // Update slave element info
        rDerivativeData.UpdateMasterPair(r_master_geometry, rCurrentProcessInfo);

        // Initialize the mortar operators
        mPreviousMortarOperators.Initialize();

        const bool dual_LM = DerivativesUtilitiesType::CalculateAeAndDeltaAe(r_slave_geometry, r_normal_slave, r_master_geometry, rDerivativeData, rVariables, consider_normal_variation, conditions_points_slave, this_integration_method, this->GetAxisymmetricCoefficient(rVariables));

        for (IndexType i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom) {
            PointerVector<PointType> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
            array_1d<BelongType, TDim> belong_array;
            for (IndexType i_node = 0; i_node < TDim; ++i_node) {
                PointType global_point;
                r_slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                points_array(i_node) = Kratos::make_shared<PointType>(PointType(global_point));
                belong_array[i_node] = conditions_points_slave[i_geom][i_node].GetBelong();
            }

            DecompositionType decomp_geom( points_array );

            const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, r_slave_geometry.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);

            if (bad_shape == false) {
                const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( this_integration_method );

                // Integrating the mortar operators
                for ( IndexType point_number = 0; point_number < integration_points_slave.size(); ++point_number ) {
                    // We compute the local coordinates

                    const PointType local_point_decomp = PointType(integration_points_slave[point_number].Coordinates());
                    PointType local_point_parent;
                    PointType gp_global;
                    decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                    r_slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);

                    // Calculate the kinematic variables
                    MortarExplicitContributionUtilities<TDim, TNumNodes, FrictionalCase::FRICTIONAL_PENALTY, TNormalVariation, TNumNodesMaster>::CalculateKinematics(this, rVariables, rDerivativeData, r_normal_master, local_point_decomp, local_point_parent, decomp_geom, dual_LM);

                    const double integration_weight = integration_points_slave[point_number].Weight() * this->GetAxisymmetricCoefficient(rVariables);

                    mPreviousMortarOperators.CalculateMortarOperators(rVariables, integration_weight);
                }
            }
        }
    }
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
