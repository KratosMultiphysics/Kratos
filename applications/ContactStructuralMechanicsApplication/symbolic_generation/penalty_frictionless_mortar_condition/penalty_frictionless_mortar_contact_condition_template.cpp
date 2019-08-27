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
#include "custom_conditions/penalty_frictionless_mortar_contact_condition.h"

namespace Kratos
{
/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, bool TNormalVariation, SizeType TNumNodesMaster >
Condition::Pointer PenaltyMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation,TNumNodesMaster>::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesPointerType pProperties ) const
{
    return Kratos::make_intrusive< PenaltyMethodFrictionlessMortarContactCondition<TDim,TNumNodes, TNormalVariation, TNumNodesMaster > >( NewId, this->GetParentGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, bool TNormalVariation, SizeType TNumNodesMaster >
Condition::Pointer PenaltyMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation,TNumNodesMaster>::Create(
    IndexType NewId,
    GeometryPointerType pGeom,
    PropertiesPointerType pProperties) const
{
    return Kratos::make_intrusive<  PenaltyMethodFrictionlessMortarContactCondition<TDim,TNumNodes, TNormalVariation, TNumNodesMaster > >( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, bool TNormalVariation, SizeType TNumNodesMaster >
Condition::Pointer PenaltyMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation,TNumNodesMaster>::Create(
    IndexType NewId,
    GeometryPointerType pGeom,
    PropertiesType::Pointer pProperties,
    GeometryType::Pointer pMasterGeom) const
{
    return Kratos::make_intrusive<  PenaltyMethodFrictionlessMortarContactCondition<TDim,TNumNodes, TNormalVariation, TNumNodesMaster > >( NewId, pGeom, pProperties, pMasterGeom );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, bool TNormalVariation, SizeType TNumNodesMaster >
PenaltyMethodFrictionlessMortarContactCondition<TDim,TNumNodes, TNormalVariation, TNumNodesMaster>::~PenaltyMethodFrictionlessMortarContactCondition( )
= default;

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, bool TNormalVariation, SizeType TNumNodesMaster >
void PenaltyMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation,TNumNodesMaster>::AddExplicitContribution(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const IndexType integration_order = this->GetProperties().Has(INTEGRATION_ORDER_CONTACT) ? this->GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) : 2;
    MortarExplicitContributionUtilities<TDim, TNumNodes, FrictionalCase::FRICTIONLESS_PENALTY, TNormalVariation, TNumNodesMaster>::AddExplicitContributionOfMortarCondition(this, rCurrentProcessInfo, integration_order, false, true);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, bool TNormalVariation, SizeType TNumNodesMaster >
void PenaltyMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation,TNumNodesMaster>::AddExplicitContribution(
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
void PenaltyMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation,TNumNodesMaster>::AddExplicitContribution(
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
            GeometryType& r_slave_geometry = this->GetParentGeometry();
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

template< SizeType TDim, SizeType TNumNodes, bool TNormalVariation, SizeType TNumNodesMaster >
void PenaltyMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation,TNumNodesMaster>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& CurrentProcessInfo
    )
{
    KRATOS_TRY;

    if (rResult.size() != MatrixSize)
        rResult.resize( MatrixSize, false );

    IndexType index = 0;

    /* ORDER - [ MASTER, SLAVE ] */
    GeometryType& r_slave_geometry = this->GetParentGeometry();
    GeometryType& r_master_geometry = this->GetPairedGeometry();

    // Master Nodes Displacement Equation IDs
    for ( IndexType i_master = 0; i_master < TNumNodesMaster; ++i_master ) { // NOTE: Assuming same number of nodes for master and slave
        NodeType& r_master_node = r_master_geometry[i_master];
        rResult[index++] = r_master_node.GetDof( DISPLACEMENT_X ).EquationId( );
        rResult[index++] = r_master_node.GetDof( DISPLACEMENT_Y ).EquationId( );
        if (TDim == 3) rResult[index++] = r_master_node.GetDof( DISPLACEMENT_Z ).EquationId( );
    }

    // Slave Nodes Displacement Equation IDs
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        NodeType& r_slave_node = r_slave_geometry[i_slave];
        rResult[index++] = r_slave_node.GetDof( DISPLACEMENT_X ).EquationId( );
        rResult[index++] = r_slave_node.GetDof( DISPLACEMENT_Y ).EquationId( );
        if (TDim == 3) rResult[index++] = r_slave_node.GetDof( DISPLACEMENT_Z ).EquationId( );
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, bool TNormalVariation, SizeType TNumNodesMaster >
void PenaltyMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation,TNumNodesMaster>::GetDofList(
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    if (rConditionalDofList.size() != MatrixSize)
        rConditionalDofList.resize( MatrixSize );

    IndexType index = 0;

    /* ORDER - [ MASTER, SLAVE ] */
    GeometryType& r_slave_geometry = this->GetParentGeometry();
    GeometryType& r_master_geometry = this->GetPairedGeometry();

    // Master Nodes Displacement Equation IDs
    for ( IndexType i_master = 0; i_master < TNumNodesMaster; ++i_master ) { // NOTE: Assuming same number of nodes for master and slave
        NodeType& r_master_node = r_master_geometry[i_master];
        rConditionalDofList[index++] = r_master_node.pGetDof( DISPLACEMENT_X );
        rConditionalDofList[index++] = r_master_node.pGetDof( DISPLACEMENT_Y );
        if (TDim == 3) rConditionalDofList[index++] = r_master_node.pGetDof( DISPLACEMENT_Z );
    }

    // Slave Nodes Displacement Equation IDs
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        NodeType& r_slave_node = r_slave_geometry[i_slave];
        rConditionalDofList[index++] = r_slave_node.pGetDof( DISPLACEMENT_X );
        rConditionalDofList[index++] = r_slave_node.pGetDof( DISPLACEMENT_Y );
        if (TDim == 3) rConditionalDofList[index++] = r_slave_node.pGetDof( DISPLACEMENT_Z );
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, bool TNormalVariation, SizeType TNumNodesMaster >
int PenaltyMethodFrictionlessMortarContactCondition<TDim,TNumNodes,TNormalVariation,TNumNodesMaster>::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int ierr = BaseType::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(NORMAL)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    GeometryType& r_slave_geometry = this->GetParentGeometry();
    for ( IndexType i = 0; i < TNumNodes; ++i ) {
        NodeType& r_node = r_slave_geometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NORMAL,r_node)
    }

    return ierr;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template class PenaltyMethodFrictionlessMortarContactCondition<2, 2, false>;
template class PenaltyMethodFrictionlessMortarContactCondition<3, 3, false, 3>;
template class PenaltyMethodFrictionlessMortarContactCondition<3, 4, false, 4>;
template class PenaltyMethodFrictionlessMortarContactCondition<3, 3, false, 4>;
template class PenaltyMethodFrictionlessMortarContactCondition<3, 4, false, 3>;
template class PenaltyMethodFrictionlessMortarContactCondition<2, 2, true>;
template class PenaltyMethodFrictionlessMortarContactCondition<3, 3, true, 3>;
template class PenaltyMethodFrictionlessMortarContactCondition<3, 4, true, 4>;
template class PenaltyMethodFrictionlessMortarContactCondition<3, 3, true, 4>;
template class PenaltyMethodFrictionlessMortarContactCondition<3, 4, true, 3>;

} // Namespace Kratos
