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
#include "custom_conditions/ALM_frictionless_mortar_contact_condition.h"

namespace Kratos 
{
/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
Condition::Pointer AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation>::Create( 
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesPointerType pProperties ) const
{
    return boost::make_shared< AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes, TNormalVariation > >( NewId, this->GetGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
Condition::Pointer AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation>::Create(
    IndexType NewId,
    GeometryPointerType pGeom,
    PropertiesPointerType pProperties) const
{
    return boost::make_shared<  AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes, TNormalVariation > >( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
Condition::Pointer AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation>::Create(
    IndexType NewId,
    GeometryPointerType pGeom,
    PropertiesType::Pointer pProperties,
    GeometryType::Pointer pMasterGeom) const
{
    return boost::make_shared<  AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes, TNormalVariation > >( NewId, pGeom, pProperties, pMasterGeom );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes, TNormalVariation>::~AugmentedLagrangianMethodFrictionlessMortarContactCondition( )
= default;

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

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& CurrentProcessInfo 
    )
{
    KRATOS_TRY;   
    
    if (rResult.size() != MatrixSize)
    {
        rResult.resize( MatrixSize, false );
    }
    
    unsigned int index = 0;
    
    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    GeometryType& current_master = this->GetPairedGeometry();;
    
    // Master Nodes Displacement Equation IDs
    for ( unsigned int i_master = 0; i_master < TNumNodes; i_master++ ) // NOTE: Assuming same number of nodes for master and slave
    {
        NodeType& master_node = current_master[i_master];
        rResult[index++] = master_node.GetDof( DISPLACEMENT_X ).EquationId( );
        rResult[index++] = master_node.GetDof( DISPLACEMENT_Y ).EquationId( );
        if (TDim == 3) rResult[index++] = master_node.GetDof( DISPLACEMENT_Z ).EquationId( );
    }

    // Slave Nodes Displacement Equation IDs
    for ( unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++ )
    {
        NodeType& slave_node = this->GetGeometry()[ i_slave ];
        rResult[index++] = slave_node.GetDof( DISPLACEMENT_X ).EquationId( );
        rResult[index++] = slave_node.GetDof( DISPLACEMENT_Y ).EquationId( );
        if (TDim == 3) rResult[index++] = slave_node.GetDof( DISPLACEMENT_Z ).EquationId( );
    }

    // Slave Nodes  Lambda Equation IDs
    for ( unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++ )
    {
        NodeType& slave_node = this->GetGeometry()[ i_slave ];
        rResult[index++] = slave_node.GetDof( NORMAL_CONTACT_STRESS ).EquationId( );
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation>::GetDofList(
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;
    
    if (rConditionalDofList.size() != MatrixSize)
    {
        rConditionalDofList.resize( MatrixSize );
    }
    
    unsigned int index = 0;
    
    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    GeometryType& current_master = this->GetPairedGeometry();;

    // Master Nodes Displacement Equation IDs
    for ( unsigned int i_master = 0; i_master < TNumNodes; ++i_master ) // NOTE: Assuming same number of nodes for master and slave
    {
        NodeType& master_node = current_master[i_master];
        rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_X );
        rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_Y );
        if (TDim == 3) rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_Z );
    }

    // Slave Nodes Displacement Equation IDs
    for ( unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave )
    {
        NodeType& slave_node = this->GetGeometry()[ i_slave ];
        rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_X );
        rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_Y );
        if (TDim == 3) rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_Z );
    }

    // Slave Nodes Lambda Equation IDs
    for ( unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave )
    {
        NodeType& slave_node = this->GetGeometry()[ i_slave ];
        rConditionalDofList[index++] = slave_node.pGetDof( NORMAL_CONTACT_STRESS );
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
int AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation>::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int ierr = BaseType::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(NORMAL)
    KRATOS_CHECK_VARIABLE_KEY(NORMAL_CONTACT_STRESS)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( unsigned int i = 0; i < TNumNodes; ++i )
    {
        Node<3> &rnode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NORMAL,rnode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NORMAL_CONTACT_STRESS,rnode)

        KRATOS_CHECK_DOF_IN_NODE(NORMAL_CONTACT_STRESS, rnode)
    }

    return ierr;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template class AugmentedLagrangianMethodFrictionlessMortarContactCondition<2, 2, false>;
template class AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 3, false>;
template class AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 4, false>;
template class AugmentedLagrangianMethodFrictionlessMortarContactCondition<2, 2, true>;
template class AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 3, true>;
template class AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 4, true>;

} // Namespace Kratos
