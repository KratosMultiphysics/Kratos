// KRATOS  ___|  |       |       |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//           | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License: BSD License
//   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrándiz
//

// System includes

// External includes

// Project includes
/* Mortar includes */
#include "custom_conditions/ALM_frictionless_mortar_contact_condition.h"

/* Utilities */
#include "custom_utilities/contact_utilities.h"

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
    return boost::make_shared< AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes, TNormalVariation> >( NewId, this->GetGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
Condition::Pointer AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation>::Create(
    IndexType NewId,
    GeometryPointerType pGeom,
    PropertiesPointerType pProperties) const
{
    return boost::make_shared< AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes, TNormalVariation> >( NewId, pGeom, pProperties );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes, TNormalVariation>::~AugmentedLagrangianMethodFrictionlessMortarContactCondition( )
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
void AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& CurrentProcessInfo 
    )
{
    KRATOS_TRY;   
    
    boost::shared_ptr<ConditionSet>& AllConditionSets = this->GetValue( CONTACT_SETS );
    
    // Calculates the size of the system
    const unsigned int ConditionSize = (TDim * ( TNumNodes + TNumNodes) + TNumNodes)* AllConditionSets->size(); 
    
    if (rResult.size() != ConditionSize)
    {
        rResult.resize( ConditionSize, false );
    }
    
    unsigned int index = 0;
    
    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    for (auto itPair = AllConditionSets->begin(); itPair != AllConditionSets->end(); ++itPair )
    {
        GeometryType& CurrentMaster = (*itPair)->GetGeometry( );
        
        // Master Nodes Displacement Equation IDs
        for ( unsigned int iMaster = 0; iMaster < TNumNodes; iMaster++ ) // NOTE: Assuming same number of nodes for master and slave
        {
            NodeType& MasterNode = CurrentMaster[iMaster];
            rResult[index++] = MasterNode.GetDof( DISPLACEMENT_X ).EquationId( );
            rResult[index++] = MasterNode.GetDof( DISPLACEMENT_Y ).EquationId( );
            if (TDim == 3)
            {
                rResult[index++] = MasterNode.GetDof( DISPLACEMENT_Z ).EquationId( );
            }
        }

        // Slave Nodes Displacement Equation IDs
        for ( unsigned int iSlave = 0; iSlave < TNumNodes; iSlave++ )
        {
            NodeType& SlaveNode = this->GetGeometry()[ iSlave ];
            rResult[index++] = SlaveNode.GetDof( DISPLACEMENT_X ).EquationId( );
            rResult[index++] = SlaveNode.GetDof( DISPLACEMENT_Y ).EquationId( );
            if (TDim == 3)
            {
                rResult[index++] = SlaveNode.GetDof( DISPLACEMENT_Z ).EquationId( );
            }
        }

        // Slave Nodes  Lambda Equation IDs
        for ( unsigned int iSlave = 0; iSlave < TNumNodes; iSlave++ )
        {
            NodeType& SlaveNode = this->GetGeometry()[ iSlave ];
            rResult[index++] = SlaveNode.GetDof( NORMAL_CONTACT_STRESS ).EquationId( );
        }
        
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
    
    boost::shared_ptr<ConditionSet>& AllConditionSets = this->GetValue( CONTACT_SETS );
    
    // Calculates the size of the system
    const unsigned int ConditionSize = (TDim * ( TNumNodes + TNumNodes) + TNumNodes)* AllConditionSets->size(); 
    
    if (rConditionalDofList.size() != ConditionSize)
    {
        rConditionalDofList.resize( ConditionSize );
    }
    
    unsigned int index = 0;
    
    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    for (auto itPair = AllConditionSets->begin(); itPair != AllConditionSets->end(); ++itPair )
    {
        GeometryType& CurrentMaster = (*itPair)->GetGeometry( );

        // Master Nodes Displacement Equation IDs
        for ( unsigned int iMaster = 0; iMaster < TNumNodes; iMaster++ ) // NOTE: Assuming same number of nodes for master and slave
        {
            NodeType& MasterNode = CurrentMaster[iMaster];
            rConditionalDofList[index++] = MasterNode.pGetDof( DISPLACEMENT_X );
            rConditionalDofList[index++] = MasterNode.pGetDof( DISPLACEMENT_Y );
            if (TDim == 3)
            {
                rConditionalDofList[index++] = MasterNode.pGetDof( DISPLACEMENT_Z );
            }
        }

        // Slave Nodes Displacement Equation IDs
        for ( unsigned int iSlave = 0; iSlave < TNumNodes; iSlave++ )
        {
            NodeType& SlaveNode = this->GetGeometry()[ iSlave ];
            rConditionalDofList[index++] = SlaveNode.pGetDof( DISPLACEMENT_X );
            rConditionalDofList[index++] = SlaveNode.pGetDof( DISPLACEMENT_Y );
            if (TDim == 3)
            {
                rConditionalDofList[index++] = SlaveNode.pGetDof( DISPLACEMENT_Z );
            }
        }

        // Slave Nodes Lambda Equation IDs
        for ( unsigned int iSlave = 0; iSlave < TNumNodes; iSlave++ )
        {
            NodeType& SlaveNode = this->GetGeometry()[ iSlave ];
            rConditionalDofList[index++] = SlaveNode.pGetDof( NORMAL_CONTACT_STRESS );
        }
    }
    
    KRATOS_CATCH( "" );
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
