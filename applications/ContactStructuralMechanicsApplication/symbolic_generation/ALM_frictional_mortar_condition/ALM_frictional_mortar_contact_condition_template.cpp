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

/* Utilities */
#include "custom_utilities/contact_utilities.h"

namespace Kratos 
{
/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
Condition::Pointer AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::Create( 
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesPointerType pProperties ) const
{
    return boost::make_shared< AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation> >( NewId, this->GetGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
Condition::Pointer AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::Create(
    IndexType NewId,
    GeometryPointerType pGeom,
    PropertiesPointerType pProperties) const
{
    return boost::make_shared< AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation> >( NewId, pGeom, pProperties );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::~AugmentedLagrangianMethodFrictionalMortarContactCondition( )
{
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

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& CurrentProcessInfo 
    )
{
    KRATOS_TRY;   
    
    ConditionMap::Pointer& all_conditions_maps = this->GetValue( MAPPING_PAIRS );
    
    // Calculates the size of the system
    const unsigned int condition_size = (TDim * ( TNumNodes + TNumNodes) + TNumNodes)* all_conditions_maps->size(); 
    
    if (rResult.size() != condition_size)
    {
        rResult.resize( condition_size, false );
    }
    
    unsigned int index = 0;
    
    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    for (auto it_pair = all_conditions_maps->begin(); it_pair != all_conditions_maps->end(); ++it_pair )
    {
        GeometryType& current_master = (it_pair->first)->GetGeometry( );
        
        // Master Nodes Displacement Equation IDs
        for ( unsigned int i_master = 0; i_master < TNumNodes; i_master++ ) // NOTE: Assuming same number of nodes for master and slave
        {
            NodeType& master_node = current_master[i_master];
            rResult[index++] = master_node.GetDof( DISPLACEMENT_X ).EquationId( );
            rResult[index++] = master_node.GetDof( DISPLACEMENT_Y ).EquationId( );
            if (TDim == 3)
            {
                rResult[index++] = master_node.GetDof( DISPLACEMENT_Z ).EquationId( );
            }
        }

        // Slave Nodes Displacement Equation IDs
        for ( unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++ )
        {
            NodeType& slave_node = this->GetGeometry()[ i_slave ];
            rResult[index++] = slave_node.GetDof( DISPLACEMENT_X ).EquationId( );
            rResult[index++] = slave_node.GetDof( DISPLACEMENT_Y ).EquationId( );
            if (TDim == 3)
            {
                rResult[index++] = slave_node.GetDof( DISPLACEMENT_Z ).EquationId( );
            }
        }

        // Slave Nodes  Lambda Equation IDs
        for ( unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++ )
        {
            NodeType& slave_node = this->GetGeometry()[ i_slave ];
            rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_X ).EquationId( );
            rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Y ).EquationId( );
            if (TDim == 3)
            {
                rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Z ).EquationId( );
            }
        }
        
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::GetDofList(
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo 
)
{
    KRATOS_TRY;
    
    ConditionMap::Pointer& all_conditions_maps = this->GetValue( MAPPING_PAIRS );
    
    // Calculates the size of the system
    const unsigned int condition_size = (TDim * ( TNumNodes + TNumNodes) + TNumNodes)* all_conditions_maps->size(); 
    
    if (rConditionalDofList.size() != condition_size)
    {
        rConditionalDofList.resize( condition_size );
    }
    
    unsigned int index = 0;
    
    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    for (auto it_pair = all_conditions_maps->begin(); it_pair != all_conditions_maps->end(); ++it_pair )
    {
        GeometryType& current_master = (it_pair->first)->GetGeometry( );

        // Master Nodes Displacement Equation IDs
        for ( unsigned int i_master = 0; i_master < TNumNodes; i_master++ ) // NOTE: Assuming same number of nodes for master and slave
        {
            NodeType& master_node = current_master[i_master];
            rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_X );
            rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_Y );
            if (TDim == 3)
            {
                rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_Z );
            }
        }

        // Slave Nodes Displacement Equation IDs
        for ( unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++ )
        {
            NodeType& slave_node = this->GetGeometry()[ i_slave ];
            rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_X );
            rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_Y );
            if (TDim == 3)
            {
                rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_Z );
            }
        }

        // Slave Nodes Lambda Equation IDs
        for ( unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++ )
        {
            NodeType& slave_node = this->GetGeometry()[ i_slave ];
            rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_X );
            rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Y );
            if (TDim == 3)
            {
                rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Z );
            }
        }
    }
    
    KRATOS_CATCH( "" );
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
