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
#include "custom_conditions/ALM_frictional_mortar_contact_condition.h"

/* Utilities */
#include "custom_utilities/contact_utilities.h"

namespace Kratos 
{
/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes>::Create( 
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesPointerType pProperties ) const
{
    return boost::make_shared< AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes> >( NewId, this->GetGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes>::Create(
    IndexType NewId,
    GeometryPointerType pGeom,
    PropertiesPointerType pProperties) const
{
    return boost::make_shared< AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes> >( NewId, pGeom, pProperties );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes>::~AugmentedLagrangianMethodFrictionalMortarContactCondition( )
{
}

/***************************** BEGIN AD REPLACEMENT ********************************/
/***********************************************************************************/

template<>
bounded_matrix<double, 12, 12> AugmentedLagrangianMethodFrictionalMortarContactCondition<2,2>::CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const double& rPenaltyFactor,
        const double& rScaleFactor,
        const unsigned int& rActiveInactive
        )
{
    bounded_matrix<double,12,12> lhs = ZeroMatrix(12, 12);

    return lhs;
}


/****************************** END AD REPLACEMENT *********************************/
/***********************************************************************************/

/***************************** BEGIN AD REPLACEMENT ********************************/
/***********************************************************************************/

template<>
array_1d<double, 12> AugmentedLagrangianMethodFrictionalMortarContactCondition<2,2>::CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const double& rPenaltyFactor,
        const double& rScaleFactor,
        const unsigned int& rActiveInactive
        )
{
    array_1d<double,12> rhs(0.0,12);

    return rhs;
}


/****************************** END AD REPLACEMENT *********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& CurrentProcessInfo 
    )
{
    KRATOS_TRY;   
    
    boost::shared_ptr<ConditionMap>& AllConditionSets = this->GetValue( CONTACT_SETS );
    
    // Calculates the size of the system
    const unsigned int ConditionSize = (TDim * ( TNumNodes + TNumNodes + TNumNodes) ) * AllConditionSets->size(); 
    
    if (rResult.size() != ConditionSize)
    {
        rResult.resize( ConditionSize, false );
    }
    
    unsigned int index = 0;
    
    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    for (auto ipair = AllConditionSets->begin(); ipair != AllConditionSets->end(); ++ipair )
    {
        GeometryType& current_master = (ipair->first)->GetGeometry( );
        
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

template< unsigned int TDim, unsigned int TNumNodes>
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim, TNumNodes>::GetDofList(
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo 
)
{
    KRATOS_TRY;
    
    boost::shared_ptr<ConditionMap>& AllConditionSets = this->GetValue( CONTACT_SETS );
    
    // Calculates the size of the system
    const unsigned int ConditionSize = (TDim * ( TNumNodes + TNumNodes + TNumNodes) ) * AllConditionSets->size(); 
    
    if (rConditionalDofList.size() != ConditionSize)
    {
        rConditionalDofList.resize( ConditionSize );
    }
    
    unsigned int index = 0;
    
    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    for (auto ipair = AllConditionSets->begin(); ipair != AllConditionSets->end(); ++ipair )
    {
        GeometryType& current_master = (ipair->first)->GetGeometry( );

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

template class AugmentedLagrangianMethodFrictionalMortarContactCondition<2, 2>;
// template class AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 3>;
// template class AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 4>;

} // Namespace Kratos
