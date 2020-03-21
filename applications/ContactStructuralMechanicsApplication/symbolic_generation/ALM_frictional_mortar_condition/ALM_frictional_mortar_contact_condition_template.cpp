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
#include "custom_conditions/ALM_frictional_mortar_contact_condition.h"

namespace Kratos
{
/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
Condition::Pointer AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesPointerType pProperties ) const
{
    return Kratos::make_intrusive< AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster> >( NewId, this->GetParentGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
Condition::Pointer AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::Create(
    IndexType NewId,
    GeometryPointerType pGeom,
    PropertiesPointerType pProperties) const
{
    return Kratos::make_intrusive< AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster> >( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
Condition::Pointer AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::Create(
    IndexType NewId,
    GeometryPointerType pGeom,
    PropertiesType::Pointer pProperties,
    GeometryType::Pointer pMasterGeom ) const
{
    return Kratos::make_intrusive< AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster> >( NewId, pGeom, pProperties, pMasterGeom );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::~AugmentedLagrangianMethodFrictionalMortarContactCondition( )
= default;

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::Initialize( )
{
    KRATOS_TRY;

    BaseType::Initialize();

    // We initailize the previous mortar operators
    mPreviousMortarOperators.Initialize();
    mPreviousMortarOperatorsInitialized = false;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
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
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    BaseType::FinalizeSolutionStep(rCurrentProcessInfo);

    // Setting previous mortar operators flag
    mPreviousMortarOperatorsInitialized = false;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::AddExplicitContribution(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const IndexType integration_order = this->GetProperties().Has(INTEGRATION_ORDER_CONTACT) ? this->GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) : 2;
    MortarExplicitContributionUtilities<TDim, TNumNodes, FrictionalCase::FRICTIONAL_PENALTY, TNormalVariation, TNumNodesMaster>::AddExplicitContributionOfMortarFrictionalCondition(this, rCurrentProcessInfo, mPreviousMortarOperators, integration_order, false, false);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, bool TNormalVariation, SizeType TNumNodesMaster >
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation,TNumNodesMaster>::ComputePreviousMortarOperators( ProcessInfo& rCurrentProcessInfo)
{
    const IndexType integration_order = this->GetProperties().Has(INTEGRATION_ORDER_CONTACT) ? this->GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) : 2;
    MortarExplicitContributionUtilities<TDim, TNumNodes, FrictionalCase::FRICTIONAL_PENALTY, TNormalVariation, TNumNodesMaster>::ComputePreviousMortarOperators(this, rCurrentProcessInfo, mPreviousMortarOperators, integration_order, false);
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
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::EquationIdVector(
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
    for ( IndexType i_master = 0; i_master < TNumNodesMaster; ++i_master ) { // NOTE: Assuming same number of nodes for master and slave
        NodeType& master_node = current_master[i_master];
        rResult[index++] = master_node.GetDof( DISPLACEMENT_X ).EquationId( );
        rResult[index++] = master_node.GetDof( DISPLACEMENT_Y ).EquationId( );
        if (TDim == 3) rResult[index++] = master_node.GetDof( DISPLACEMENT_Z ).EquationId( );
    }

    // Slave Nodes Displacement Equation IDs
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        NodeType& slave_node = this->GetParentGeometry()[ i_slave ];
        rResult[index++] = slave_node.GetDof( DISPLACEMENT_X ).EquationId( );
        rResult[index++] = slave_node.GetDof( DISPLACEMENT_Y ).EquationId( );
        if (TDim == 3) rResult[index++] = slave_node.GetDof( DISPLACEMENT_Z ).EquationId( );
    }

    // Slave Nodes  Lambda Equation IDs
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        NodeType& slave_node = this->GetParentGeometry()[ i_slave ];
        rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_X ).EquationId( );
        rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Y ).EquationId( );
        if (TDim == 3) rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Z ).EquationId( );
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::GetDofList(
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
    for ( IndexType i_master = 0; i_master < TNumNodesMaster; ++i_master ){ // NOTE: Assuming same number of nodes for master and slave
        NodeType& master_node = current_master[i_master];
        rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_X );
        rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_Y );
        if (TDim == 3) rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_Z );
    }

    // Slave Nodes Displacement Equation IDs
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        NodeType& slave_node = this->GetParentGeometry()[ i_slave ];
        rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_X );
        rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_Y );
        if (TDim == 3) rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_Z );
    }

    // Slave Nodes Lambda Equation IDs
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        NodeType& slave_node = this->GetParentGeometry()[ i_slave ];
        rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_X );
        rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Y );
        if (TDim == 3) rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Z );
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
int AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::Check( const ProcessInfo& rCurrentProcessInfo )
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
        Node<3> &rnode = this->GetParentGeometry()[i];
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

template class AugmentedLagrangianMethodFrictionalMortarContactCondition<2, 2, false, 2>;
template class AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 3, false, 3>;
template class AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 4, false, 4>;
template class AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 3, false, 4>;
template class AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 4, false, 3>;
template class AugmentedLagrangianMethodFrictionalMortarContactCondition<2, 2, true,  2>;
template class AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 3, true,  3>;
template class AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 4, true,  4>;
template class AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 3, true,  4>;
template class AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 4, true,  3>;

} // Namespace Kratos
