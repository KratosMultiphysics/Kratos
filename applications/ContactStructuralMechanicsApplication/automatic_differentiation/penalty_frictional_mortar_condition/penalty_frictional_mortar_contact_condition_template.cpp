// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
/* Mortar includes */
#include "custom_utilities/mortar_explicit_contribution_utilities.h"
#include "custom_conditions/penalty_frictional_mortar_contact_condition.h"
#include "utilities/atomic_utilities.h"

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
    return Kratos::make_intrusive< PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster > >( NewId, this->GetParentGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
Condition::Pointer PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::Create(
    IndexType NewId,
    GeometryPointerType pGeom,
    PropertiesPointerType pProperties) const
{
    return Kratos::make_intrusive< PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster> >( NewId, pGeom, pProperties );
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
    return Kratos::make_intrusive< PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster> >( NewId, pGeom, pProperties, pMasterGeom );
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
void PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::Initialize(rCurrentProcessInfo);

    // We initailize the previous mortar operators
    mPreviousMortarOperators.Initialize();

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
void PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
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
void PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
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
void PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::AddExplicitContribution(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    MortarExplicitContributionUtilities<TDim, TNumNodes, FrictionalCase::FRICTIONAL_PENALTY, TNormalVariation, TNumNodesMaster>::AddExplicitContributionOfMortarFrictionalCondition(this, rCurrentProcessInfo, mPreviousMortarOperators, false, false);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, bool TNormalVariation, SizeType TNumNodesMaster >
void PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation,TNumNodesMaster>::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    const Variable<double>& rDestinationVariable,
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
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // Computing the force residual
    if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL) {
        if (this->Is(ACTIVE)) {
            // Getting geometries
            GeometryType& r_slave_geometry = this->GetParentGeometry();
            GeometryType& r_master_geometry = this->GetParentGeometry();

            for ( IndexType i_master = 0; i_master < TNumNodesMaster; ++i_master ) {
                Node& r_master_node = r_master_geometry[i_master];
                const IndexType index = TDim * i_master;

                array_1d<double, 3>& r_force_residual = r_master_node.FastGetSolutionStepValue(FORCE_RESIDUAL);

                for (IndexType j = 0; j < TDim; ++j) {
                    AtomicAdd(r_force_residual[j], rRHSVector[index + j]);
                }
            }
            for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
                Node& r_slave_node = r_slave_geometry[i_slave];
                const IndexType index = TDim * (TNumNodesMaster + i_slave);

                array_1d<double, 3>& r_force_residual = r_slave_node.FastGetSolutionStepValue(FORCE_RESIDUAL);

                for (IndexType j = 0; j < TDim; ++j) {
                    AtomicAdd(r_force_residual[j], rRHSVector[index + j]);
                }
            }
        }
    }

    KRATOS_CATCH("")
}


/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, bool TNormalVariation, SizeType TNumNodesMaster >
void PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation,TNumNodesMaster>::ComputePreviousMortarOperators(const ProcessInfo& rCurrentProcessInfo)
{
    const IndexType integration_order = this->GetProperties().Has(INTEGRATION_ORDER_CONTACT) ? this->GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) : 2;
    MortarExplicitContributionUtilities<TDim, TNumNodes, FrictionalCase::FRICTIONAL_PENALTY, TNormalVariation, TNumNodesMaster>::ComputePreviousMortarOperators(this, rCurrentProcessInfo, mPreviousMortarOperators, integration_order, false);
}

/***************************** BEGIN AD REPLACEMENT ********************************/
/***********************************************************************************/

// replace_lhs

/****************************** END AD REPLACEMENT *********************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
void PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::CalculateLocalRHS(
    Vector& rLocalRHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DerivativeDataType& rDerivativeData,
    const IndexType rActiveInactive,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    StaticCalculateLocalRHS(this, mPreviousMortarOperators, GetFrictionCoefficient(), rLocalRHS, rMortarConditionMatrices, rDerivativeData, rActiveInactive, rCurrentProcessInfo);
}

/***************************** BEGIN AD REPLACEMENT ********************************/
/***********************************************************************************/

// replace_rhs

/****************************** END AD REPLACEMENT *********************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
void PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& CurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    if (rResult.size() != MatrixSize)
        rResult.resize( MatrixSize, false );

    IndexType index = 0;

    /* ORDER - [ MASTER, SLAVE ] */
    const GeometryType& r_slave_geometry = this->GetParentGeometry();
    const GeometryType& r_master_geometry = this->GetPairedGeometry();

    // Master Nodes Displacement Equation IDs
    for ( IndexType i_master = 0; i_master < TNumNodesMaster; ++i_master ) { // NOTE: Assuming same number of nodes for master and slave
        const Node& r_master_node = r_master_geometry[i_master];
        rResult[index++] = r_master_node.GetDof( DISPLACEMENT_X ).EquationId( );
        rResult[index++] = r_master_node.GetDof( DISPLACEMENT_Y ).EquationId( );
        if constexpr (TDim == 3) rResult[index++] = r_master_node.GetDof( DISPLACEMENT_Z ).EquationId( );
    }

    // Slave Nodes Displacement Equation IDs
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        const Node& r_slave_node = r_slave_geometry[ i_slave ];
        rResult[index++] = r_slave_node.GetDof( DISPLACEMENT_X ).EquationId( );
        rResult[index++] = r_slave_node.GetDof( DISPLACEMENT_Y ).EquationId( );
        if constexpr (TDim == 3) rResult[index++] = r_slave_node.GetDof( DISPLACEMENT_Z ).EquationId( );
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
void PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::GetDofList(
    DofsVectorType& rConditionalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    if (rConditionalDofList.size() != MatrixSize)
        rConditionalDofList.resize( MatrixSize );

    IndexType index = 0;

    /* ORDER - [ MASTER, SLAVE ] */
    const GeometryType& r_slave_geometry = this->GetParentGeometry();
    const GeometryType& r_master_geometry = this->GetPairedGeometry();

    // Master Nodes Displacement Equation IDs
    for ( IndexType i_master = 0; i_master < TNumNodesMaster; ++i_master ){ // NOTE: Assuming same number of nodes for master and slave
        const Node& r_master_node = r_master_geometry[i_master];
        rConditionalDofList[index++] = r_master_node.pGetDof( DISPLACEMENT_X );
        rConditionalDofList[index++] = r_master_node.pGetDof( DISPLACEMENT_Y );
        if constexpr (TDim == 3) rConditionalDofList[index++] = r_master_node.pGetDof( DISPLACEMENT_Z );
    }

    // Slave Nodes Displacement Equation IDs
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        const Node& r_slave_node = r_slave_geometry[ i_slave ];
        rConditionalDofList[index++] = r_slave_node.pGetDof( DISPLACEMENT_X );
        rConditionalDofList[index++] = r_slave_node.pGetDof( DISPLACEMENT_Y );
        if constexpr (TDim == 3) rConditionalDofList[index++] = r_slave_node.pGetDof( DISPLACEMENT_Z );
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
int PenaltyMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation, TNumNodesMaster>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int ierr = BaseType::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    const GeometryType& r_slave_geometry = this->GetParentGeometry();
    for ( IndexType i = 0; i < TNumNodes; ++i ) {
        const Node& r_node = r_slave_geometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(WEIGHTED_SLIP, r_node)
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
