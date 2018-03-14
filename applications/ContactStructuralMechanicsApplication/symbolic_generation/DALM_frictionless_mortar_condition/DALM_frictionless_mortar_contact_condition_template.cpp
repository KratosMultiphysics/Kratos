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
#include "custom_conditions/DALM_frictionless_mortar_contact_condition.h"

namespace Kratos 
{
/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
Condition::Pointer DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation>::Create( 
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesPointerType pProperties ) const
{
    return Kratos::make_shared< DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes, TNormalVariation > >( NewId, this->GetGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
Condition::Pointer DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation>::Create(
    IndexType NewId,
    GeometryPointerType pGeom,
    PropertiesPointerType pProperties) const
{
    return Kratos::make_shared<  DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes, TNormalVariation > >( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
Condition::Pointer DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation>::Create(
    IndexType NewId,
    GeometryPointerType pGeom,
    PropertiesType::Pointer pProperties,
    GeometryType::Pointer pMasterGeom) const
{
    return Kratos::make_shared<  DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes, TNormalVariation > >( NewId, pGeom, pProperties, pMasterGeom );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes, TNormalVariation>::~DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition( )
= default;

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
void DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation>::CalculateLocalLHS(
    Matrix& rLocalLHS, 
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DerivativeDataType& rDerivativeData,
    const unsigned int rActiveInactive
    ) 
{
    BaseType::CalculateLocalLHS(rLocalLHS, rMortarConditionMatrices, rDerivativeData, rActiveInactive);
    
    constexpr std::size_t initial_index = TDim * (TNumNodes + TNumNodes);
    
    // Initialize the zero values
    for (std::size_t i = 0; i < TDim * (TNumNodes + TNumNodes) + TNumNodes; ++i) {
        for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node) {
            rLocalLHS(i, initial_index + TNumNodes + i_node) = 0.0;
            rLocalLHS(initial_index + TNumNodes + i_node, i) = 0.0;
        }
    }
    for (std::size_t i = TDim * (TNumNodes + TNumNodes) + TNumNodes; i < TDim * (TNumNodes + TNumNodes) + 2 * TNumNodes; ++i)
        for (std::size_t j = TDim * (TNumNodes + TNumNodes) + TNumNodes; j < TDim * (TNumNodes + TNumNodes) + 2 * TNumNodes; ++j)
            rLocalLHS(j, j) = 0.0;
        
    // Double ALM contribution
    CalculateLocalLHSDALM(rLocalLHS, rMortarConditionMatrices, rDerivativeData, rActiveInactive);
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
void DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation>::CalculateLocalRHS(
    Vector& rLocalRHS, 
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DerivativeDataType& rDerivativeData,
    const unsigned int rActiveInactive
    ) 
{
    BaseType::CalculateLocalRHS(rLocalRHS, rMortarConditionMatrices, rDerivativeData, rActiveInactive);
    
    constexpr std::size_t initial_index = TDim * (TNumNodes + TNumNodes);
    
    // We iterate over the nodes
    for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node) {
        rLocalRHS[initial_index + TNumNodes + i_node] = 0.0;
    }
    
    // Double ALM contribution
    CalculateLocalRHSDALM(rLocalRHS, rMortarConditionMatrices, rDerivativeData, rActiveInactive);
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
void DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;   
    
    if (rResult.size() != MatrixSize)
        rResult.resize( MatrixSize, false );
    
    unsigned int index = 0;
    
    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    GeometryType& current_master = this->GetPairedGeometry();
    
    // Master Nodes Displacement Equation IDs
    for ( unsigned int i_master = 0; i_master < TNumNodes; ++i_master ) { // NOTE: Assuming same number of nodes for master and slave
        NodeType& master_node = current_master[i_master];
        rResult[index++] = master_node.GetDof( DISPLACEMENT_X ).EquationId( );
        rResult[index++] = master_node.GetDof( DISPLACEMENT_Y ).EquationId( );
        if (TDim == 3) rResult[index++] = master_node.GetDof( DISPLACEMENT_Z ).EquationId( );
    }

    // Slave Nodes Displacement Equation IDs
    for ( unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        NodeType& slave_node = this->GetGeometry()[ i_slave ];
        rResult[index++] = slave_node.GetDof( DISPLACEMENT_X ).EquationId( );
        rResult[index++] = slave_node.GetDof( DISPLACEMENT_Y ).EquationId( );
        if (TDim == 3) rResult[index++] = slave_node.GetDof( DISPLACEMENT_Z ).EquationId( );
    }

    // Slave Nodes  Lambda Equation IDs
    for ( unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        NodeType& slave_node = this->GetGeometry()[ i_slave ];
        rResult[index++] = slave_node.GetDof( NORMAL_CONTACT_STRESS ).EquationId( );
    }
    
    // Master Nodes  Lambda Equation IDs
    for ( unsigned int i_master = 0; i_master < TNumNodes; ++i_master ) {
        NodeType& master_node = current_master[i_master];
        rResult[index++] = master_node.GetDof( NORMAL_CONTACT_STRESS ).EquationId( );
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
void DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation>::GetDofList(
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;
    
    if (rConditionalDofList.size() != MatrixSize)
        rConditionalDofList.resize( MatrixSize );
    
    unsigned int index = 0;
    
    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    GeometryType& current_master = this->GetPairedGeometry();;

    // Master Nodes Displacement Equation IDs
    for ( unsigned int i_master = 0; i_master < TNumNodes; ++i_master ) { // NOTE: Assuming same number of nodes for master and slave
        NodeType& master_node = current_master[i_master];
        rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_X );
        rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_Y );
        if (TDim == 3) rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_Z );
    }

    // Slave Nodes Displacement Equation IDs
    for ( unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        NodeType& slave_node = this->GetGeometry()[ i_slave ];
        rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_X );
        rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_Y );
        if (TDim == 3) rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_Z );
    }

    // Slave Nodes Lambda Equation IDs
    for ( unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave )  {
        NodeType& slave_node = this->GetGeometry()[ i_slave ];
        rConditionalDofList[index++] = slave_node.pGetDof( NORMAL_CONTACT_STRESS );
    }
    
    // Master Nodes Lambda Equation IDs
    for ( unsigned int i_master = 0; i_master < TNumNodes; ++i_master ) { 
        NodeType& master_node = current_master[i_master];
        rConditionalDofList[index++] = master_node.pGetDof( NORMAL_CONTACT_STRESS );
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
int DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation>::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int ierr = BaseType::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(NORMAL)
    KRATOS_CHECK_VARIABLE_KEY(NORMAL_CONTACT_STRESS)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( unsigned int i = 0; i < TNumNodes; ++i ) {
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

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
void DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation>::ResizeLHS(MatrixType& rLeftHandSideMatrix)
{
    // Resizing as needed the LHS
    if ( rLeftHandSideMatrix.size1() != MatrixSize || rLeftHandSideMatrix.size2() != MatrixSize )
            rLeftHandSideMatrix.resize( MatrixSize, MatrixSize, false );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
void DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation>::ResizeRHS(VectorType& rRightHandSideVector)
{
    // Resizing as needed the RHS
    if ( rRightHandSideVector.size() != MatrixSize )
        rRightHandSideVector.resize( MatrixSize, false );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
void DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation>::ZeroLHS(MatrixType& rLeftHandSideMatrix)
{
    rLeftHandSideMatrix = ZeroMatrix( MatrixSize, MatrixSize );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation >
void DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation>::ZeroRHS(VectorType& rRightHandSideVector)
{
    rRightHandSideVector = ZeroVector( MatrixSize );
}
    
/***********************************************************************************/
/***********************************************************************************/

template class DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition<2, 2, false>;
template class DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 3, false>;
template class DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 4, false>;
template class DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition<2, 2, true>;
template class DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 3, true>;
template class DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 4, true>;

} // Namespace Kratos
