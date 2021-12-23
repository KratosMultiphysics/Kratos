// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:		 BSD License
//					 license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
/* Mortar includes */
#include "custom_conditions/ALM_frictionless_mortar_contact_condition.h"

namespace Kratos
{
/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, bool TNormalVariation, SizeType TNumNodesMaster >
Condition::Pointer AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation,TNumNodesMaster>::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesPointerType pProperties ) const
{
    return Kratos::make_intrusive< AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes, TNormalVariation, TNumNodesMaster> >( NewId, this->GetParentGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, bool TNormalVariation, SizeType TNumNodesMaster >
Condition::Pointer AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation,TNumNodesMaster>::Create(
    IndexType NewId,
    GeometryPointerType pGeom,
    PropertiesPointerType pProperties) const
{
    return Kratos::make_intrusive<  AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes, TNormalVariation, TNumNodesMaster > >( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, bool TNormalVariation, SizeType TNumNodesMaster >
Condition::Pointer AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation,TNumNodesMaster>::Create(
    IndexType NewId,
    GeometryPointerType pGeom,
    PropertiesType::Pointer pProperties,
    GeometryType::Pointer pMasterGeom) const
{
    return Kratos::make_intrusive<  AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes, TNormalVariation, TNumNodesMaster > >( NewId, pGeom, pProperties, pMasterGeom );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, bool TNormalVariation, SizeType TNumNodesMaster >
AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes, TNormalVariation, TNumNodesMaster>::~AugmentedLagrangianMethodFrictionlessMortarContactCondition( )
= default;

/***************************** BEGIN AD REPLACEMENT ********************************/
/***********************************************************************************/

// replace_lhs

/****************************** END AD REPLACEMENT *********************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::CalculateLocalRHS(
    Vector& rLocalRHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DerivativeDataType& rDerivativeData,
    const IndexType rActiveInactive,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    StaticCalculateLocalRHS(this, rLocalRHS, rMortarConditionMatrices, rDerivativeData, rActiveInactive, rCurrentProcessInfo);
}

/***************************** BEGIN AD REPLACEMENT ********************************/
/***********************************************************************************/

// replace_rhs

/****************************** END AD REPLACEMENT *********************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, bool TNormalVariation, SizeType TNumNodesMaster >
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation,TNumNodesMaster>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& CurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    if (rResult.size() != MatrixSize)
        rResult.resize( MatrixSize, false );

    IndexType index = 0;

    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    const GeometryType& r_current_master = this->GetPairedGeometry();
    const GeometryType& r_current_slave = this->GetParentGeometry();

    // Master Nodes Displacement Equation IDs
    for ( IndexType i_master = 0; i_master < TNumNodesMaster; ++i_master ) { // NOTE: Assuming same number of nodes for master and slave
        const NodeType& r_master_node = r_current_master[i_master];
        rResult[index++] = r_master_node.GetDof( DISPLACEMENT_X ).EquationId( );
        rResult[index++] = r_master_node.GetDof( DISPLACEMENT_Y ).EquationId( );
        if (TDim == 3) rResult[index++] = r_master_node.GetDof( DISPLACEMENT_Z ).EquationId( );
    }

    // Slave Nodes Displacement Equation IDs
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        const NodeType& r_slave_node = r_current_slave[ i_slave ];
        rResult[index++] = r_slave_node.GetDof( DISPLACEMENT_X ).EquationId( );
        rResult[index++] = r_slave_node.GetDof( DISPLACEMENT_Y ).EquationId( );
        if (TDim == 3) rResult[index++] = r_slave_node.GetDof( DISPLACEMENT_Z ).EquationId( );
    }

    // Slave Nodes  Lambda Equation IDs
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        const NodeType& r_slave_node = r_current_slave[ i_slave ];
        rResult[index++] = r_slave_node.GetDof( LAGRANGE_MULTIPLIER_CONTACT_PRESSURE ).EquationId( );
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, bool TNormalVariation, SizeType TNumNodesMaster >
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation,TNumNodesMaster>::GetDofList(
    DofsVectorType& rConditionalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    if (rConditionalDofList.size() != MatrixSize)
        rConditionalDofList.resize( MatrixSize );

    IndexType index = 0;

    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    const GeometryType& r_current_master = this->GetPairedGeometry();
    const GeometryType& r_current_slave = this->GetParentGeometry();

    // Master Nodes Displacement Equation IDs
    for ( IndexType i_master = 0; i_master < TNumNodesMaster; ++i_master ) { // NOTE: Assuming same number of nodes for master and slave
        const NodeType& r_master_node = r_current_master[i_master];
        rConditionalDofList[index++] = r_master_node.pGetDof( DISPLACEMENT_X );
        rConditionalDofList[index++] = r_master_node.pGetDof( DISPLACEMENT_Y );
        if (TDim == 3) rConditionalDofList[index++] = r_master_node.pGetDof( DISPLACEMENT_Z );
    }

    // Slave Nodes Displacement Equation IDs
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        const NodeType& r_slave_node = r_current_slave[ i_slave ];
        rConditionalDofList[index++] = r_slave_node.pGetDof( DISPLACEMENT_X );
        rConditionalDofList[index++] = r_slave_node.pGetDof( DISPLACEMENT_Y );
        if (TDim == 3) rConditionalDofList[index++] = r_slave_node.pGetDof( DISPLACEMENT_Z );
    }

    // Slave Nodes Lambda Equation IDs
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave )  {
        const NodeType& r_slave_node = r_current_slave[ i_slave ];
        rConditionalDofList[index++] = r_slave_node.pGetDof( LAGRANGE_MULTIPLIER_CONTACT_PRESSURE );
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, bool TNormalVariation, SizeType TNumNodesMaster >
int AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation,TNumNodesMaster>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int ierr = BaseType::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    const GeometryType& r_current_slave = this->GetParentGeometry();
    for ( IndexType i = 0; i < TNumNodes; ++i ) {
        const NodeType &r_node = r_current_slave[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE,r_node)

        KRATOS_CHECK_DOF_IN_NODE(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE, r_node)
    }

    return ierr;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template class AugmentedLagrangianMethodFrictionlessMortarContactCondition<2, 2, false>;
template class AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 3, false, 3>;
template class AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 4, false, 4>;
template class AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 3, false, 4>;
template class AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 4, false, 3>;
template class AugmentedLagrangianMethodFrictionlessMortarContactCondition<2, 2, true>;
template class AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 3, true, 3>;
template class AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 4, true, 4>;
template class AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 3, true, 4>;
template class AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 4, true, 3>;

} // Namespace Kratos
