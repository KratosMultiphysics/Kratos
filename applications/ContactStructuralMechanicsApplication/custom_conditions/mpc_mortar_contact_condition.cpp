//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/global_variables.h"
#include "custom_conditions/mpc_mortar_contact_condition.h"

/* Utilities */
#include "utilities/geometrical_projection_utilities.h"
#include "utilities/math_utils.h"
#include "custom_utilities/mortar_explicit_contribution_utilities.h"

namespace Kratos
{

/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
Condition::Pointer MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive< MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster> >( NewId, this->GetParentGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
Condition::Pointer MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive< MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster> >( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
Condition::Pointer MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties,
    GeometryType::Pointer pMasterGeom) const
{
    return Kratos::make_intrusive< MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster> >( NewId, pGeom, pProperties, pMasterGeom );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::~MPCMortarContactCondition( )
= default;

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::Initialize( )
{
    KRATOS_TRY;

    BaseType::Initialize();

    // We initailize the previous mortar operators
    if (this->Is(SLIP)) {
        mPreviousMortarOperators.Initialize();
    }

    // Getting constraint (if possible)
    ProcessInfo aux_process_info;
    if (this->Has(CONSTRAINT_POINTER)) {
        auto p_const = this->GetValue(CONSTRAINT_POINTER);

        const IndexType integration_order = this->GetProperties().Has(INTEGRATION_ORDER_CONTACT) ? this->GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) : 2;
        MortarConditionMatrices mortar_operators;
        const bool dual_LM = MortarExplicitContributionUtilities<TDim, TNumNodes, FrictionalCase::FRICTIONLESS_PENALTY, false, TNumNodesMaster>::ComputePreviousMortarOperators(this, aux_process_info, mortar_operators, integration_order, false);

        // The relation matrix and constant vector
        Matrix relation_matrix = ZeroMatrix(TDim * TNumNodes, TDim * TNumNodesMaster);// + TDim * TNumNodes);
        Vector constant_vector = ZeroVector(TDim * TNumNodes);

        // Initialize matrices
        if (this->Is(SLIP)) {
            UpdateConstraintFrictional(mortar_operators, relation_matrix, constant_vector, aux_process_info, dual_LM);
        } else if (this->Is(RIGID)) {
            UpdateConstraintTying(mortar_operators, relation_matrix, constant_vector, aux_process_info, dual_LM);
        } else {
            UpdateConstraintFrictionless(mortar_operators, relation_matrix, constant_vector, aux_process_info, dual_LM);
        }

        // Updating DoF database
        ConstraintDofDatabaseUpdate(relation_matrix, constant_vector, aux_process_info);

        // Update local system
        p_const->SetLocalSystem(relation_matrix, constant_vector, aux_process_info);
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    BaseType::InitializeSolutionStep(rCurrentProcessInfo);

    // Only update if needed
    const bool is_blocked = this->IsDefined(BLOCKED) ? this->Is(BLOCKED) : false;
    if (!is_blocked) {
        // Compute the previous mortar operators
        if (!mPreviousMortarOperatorsInitialized && this->Is(SLIP)) {
            ComputePreviousMortarOperators(rCurrentProcessInfo);
            mPreviousMortarOperatorsInitialized = true;
        }

        // Getting constraint (if possible)
        if (this->Has(CONSTRAINT_POINTER)) {
            auto p_const = this->GetValue(CONSTRAINT_POINTER);

            if (p_const != nullptr) {
                const IndexType integration_order = this->GetProperties().Has(INTEGRATION_ORDER_CONTACT) ? this->GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) : 2;
                MortarConditionMatrices mortar_operators;
                const bool dual_LM = MortarExplicitContributionUtilities<TDim, TNumNodes, FrictionalCase::FRICTIONLESS_PENALTY, false, TNumNodesMaster>::ComputePreviousMortarOperators(this, rCurrentProcessInfo, mortar_operators, integration_order, false);

                // The relation matrix and constant vector
                Matrix relation_matrix = ZeroMatrix(TDim * TNumNodes, TDim * TNumNodesMaster);// + TDim * TNumNodes);
                Vector constant_vector = ZeroVector(TDim * TNumNodes);

                // Initialize matrices
                if (this->Is(SLIP)) {
                    UpdateConstraintFrictional(mortar_operators, relation_matrix, constant_vector, rCurrentProcessInfo, dual_LM);
                } else if (this->Is(RIGID)) {
                    UpdateConstraintTying(mortar_operators, relation_matrix, constant_vector, rCurrentProcessInfo, dual_LM);
                } else {
                    UpdateConstraintFrictionless(mortar_operators, relation_matrix, constant_vector, rCurrentProcessInfo, dual_LM);
                }

                // Updating DoF database
                ConstraintDofDatabaseUpdate(relation_matrix, constant_vector, rCurrentProcessInfo);

                // Update local system
                p_const->SetLocalSystem(relation_matrix, constant_vector, rCurrentProcessInfo);
            }
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Only update if needed
    const bool is_interactive = this->IsDefined(INTERACTION) ? this->Is(INTERACTION) : false;
    if (is_interactive) {
        // Getting constraint (if possible)
        if (this->Has(CONSTRAINT_POINTER)) {
            auto p_const = this->GetValue(CONSTRAINT_POINTER);

            if (p_const != nullptr) {
                const IndexType integration_order = this->GetProperties().Has(INTEGRATION_ORDER_CONTACT) ? this->GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) : 2;
                MortarConditionMatrices mortar_operators;
                const bool dual_LM = MortarExplicitContributionUtilities<TDim, TNumNodes, FrictionalCase::FRICTIONLESS_PENALTY, false, TNumNodesMaster>::ComputePreviousMortarOperators(this, rCurrentProcessInfo, mortar_operators, integration_order, false);

                // The relation matrix and constant vector
                Matrix relation_matrix = ZeroMatrix(TDim * TNumNodes, TDim * TNumNodesMaster);// + TDim * TNumNodes);
                Vector constant_vector = ZeroVector(TDim * TNumNodes);

                // Initialize matrices
                if (this->Is(SLIP)) {
                    UpdateConstraintFrictional(mortar_operators, relation_matrix, constant_vector, rCurrentProcessInfo, dual_LM);
                } else if (this->Is(RIGID)) {
                    UpdateConstraintTying(mortar_operators, relation_matrix, constant_vector, rCurrentProcessInfo, dual_LM);
                } else {
                    UpdateConstraintFrictionless(mortar_operators, relation_matrix, constant_vector, rCurrentProcessInfo, dual_LM);
                }

                // Updating DoF database
                ConstraintDofDatabaseUpdate(relation_matrix, constant_vector, rCurrentProcessInfo);

                // Update local system
                p_const->SetLocalSystem(relation_matrix, constant_vector, rCurrentProcessInfo);
            }
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    BaseType::FinalizeSolutionStep(rCurrentProcessInfo);

    // Compute the previous mortar operators
    if (this->Is(SLIP)) {
        ComputePreviousMortarOperators(rCurrentProcessInfo);
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    // TODO: Add somethig if necessary

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

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
        NodeType& r_slave_node = r_slave_geometry[ i_slave ];
        rResult[index++] = r_slave_node.GetDof( DISPLACEMENT_X ).EquationId( );
        rResult[index++] = r_slave_node.GetDof( DISPLACEMENT_Y ).EquationId( );
        if (TDim == 3) rResult[index++] = r_slave_node.GetDof( DISPLACEMENT_Z ).EquationId( );
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::GetDofList(
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    if (rConditionalDofList.size() != MatrixSize)
        rConditionalDofList.resize( MatrixSize );

    IndexType index = 0;

    /* ORDER - [ MASTER, SLAVE ] */
    GeometryType& r_slave_geometry = this->GetParentGeometry();
    GeometryType& r_master_geometry = this->GetPairedGeometry();

    // Master Nodes Displacement Equation IDs
    for ( IndexType i_master = 0; i_master < TNumNodesMaster; ++i_master ){ // NOTE: Assuming same number of nodes for master and slave
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

    KRATOS_CATCH("")
}


/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    // Resizing as needed the LHS
    if ( rLeftHandSideMatrix.size1() != MatrixSize || rLeftHandSideMatrix.size2() != MatrixSize )
            rLeftHandSideMatrix.resize( MatrixSize, MatrixSize, false );

    noalias(rLeftHandSideMatrix) = ZeroMatrix( MatrixSize, MatrixSize );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    // Resizing as needed the RHS
    if ( rRightHandSideVector.size() != MatrixSize )
        rRightHandSideVector.resize( MatrixSize, false );

    noalias(rRightHandSideVector) = ZeroVector( MatrixSize );

    // Adding slip frictional force
    if (this->Is(SLIP)) {
        // We got some auxiliar values
        const GeometryType& r_slave_geometry = this->GetParentGeometry();

        for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
            const NodeType& r_slave_node = r_slave_geometry[ i_slave ];
            if (r_slave_node.Is(SLIP)) {
                const array_1d<double, 3>& r_total_force = r_slave_node.FastGetSolutionStepValue(REACTION);
                const array_1d<double, 3>& r_normal = r_slave_node.FastGetSolutionStepValue(NORMAL);
                const double contact_force = inner_prod(r_total_force, r_normal);
                const double mu = r_slave_node.GetValue(FRICTION_COEFFICIENT);
                const array_1d<double, 3> tangent_force = r_total_force - contact_force * r_normal;
                const double norm_tangent_force = norm_2(tangent_force);
                const array_1d<double, 3>& r_slip = r_slave_node.FastGetSolutionStepValue(WEIGHTED_SLIP);
                const double norm_slip = norm_2(r_slip);
                array_1d<double, 3> tangent_direction;
                if (norm_tangent_force > std::numeric_limits<double>::epsilon()) {
                    noalias(tangent_direction) = tangent_force/norm_tangent_force;
                } else if (norm_slip > std::numeric_limits<double>::epsilon()) {
                    noalias(tangent_direction) = r_slip/norm_slip; // TODO: Check sign
                } else {
                    noalias(tangent_direction) = ZeroVector(3);
                }
                const array_1d<double, 3> friction_tangent_force = - mu * contact_force * tangent_direction;
                for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
                    rRightHandSideVector[TNumNodesMaster * TDim + i_slave * TDim + i_dim] += friction_tangent_force[i_dim]; // TODO: Check sign
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    rMassMatrix.resize(0, 0, false);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    rDampingMatrix.resize(0, 0, false);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::AddExplicitContribution(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const IndexType integration_order = this->GetProperties().Has(INTEGRATION_ORDER_CONTACT) ? this->GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) : 2;
    Vector slave_mortar_weigth(TNumNodes);
    if (this->Is(SLIP)) {
        const auto mortar_operators = MortarExplicitContributionUtilities<TDim, TNumNodes, FrictionalCase::FRICTIONAL_PENALTY, false, TNumNodesMaster>::AddExplicitContributionOfMortarFrictionalCondition(this, rCurrentProcessInfo, mPreviousMortarOperators, integration_order, false, false);
    } else {
        const auto mortar_operators = MortarExplicitContributionUtilities<TDim, TNumNodes, FrictionalCase::FRICTIONLESS_PENALTY, false, TNumNodesMaster>::AddExplicitContributionOfMortarCondition(this, rCurrentProcessInfo, integration_order, false, false);
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::ConstraintDofDatabaseUpdate(
    Matrix& rRelationMatrix,
    Vector& rConstantVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    auto p_const = this->GetValue(CONSTRAINT_POINTER);

    if (p_const != nullptr) {
        // Checking relation matrix
        std::vector<IndexType> slave_dofs_OK, master_dofs_OK;
        const SizeType potential_slave_size = rRelationMatrix.size1();
        const SizeType potential_master_size = rRelationMatrix.size2();
        for (IndexType i = 0; i < potential_slave_size; ++i) {
            double aux_value = 0.0;
            for (IndexType j = 0; j < rRelationMatrix.size2(); ++j) {
                aux_value += std::abs(rRelationMatrix(i, j));
            }
            if (aux_value > 1.0e-4) slave_dofs_OK.push_back(i);
        }
        for (IndexType j = 0; j < potential_master_size; ++j) {
            double aux_value = 0.0;
            for (IndexType i = 0; i < rRelationMatrix.size1(); ++i) {
                aux_value += std::abs(rRelationMatrix(i, j));
            }
            if (aux_value > 1.0e-4) master_dofs_OK.push_back(j);
        }

        // Rebuild the relation matrix and the constant vector
        const SizeType slave_size = slave_dofs_OK.size();
        const SizeType master_size = master_dofs_OK.size();
        if (slave_size < potential_slave_size || master_size < potential_master_size) {
            Matrix auxiliar_matrix(slave_size, master_size);
            Vector auxiliar_vector(slave_size);
            IndexType slave_counter = 0;
            for (IndexType i : slave_dofs_OK) {
                IndexType master_counter = 0;
                for (IndexType j : master_dofs_OK) {
                    auxiliar_matrix(slave_counter, master_counter) = rRelationMatrix(i, j);
                    ++master_counter;
                }
                auxiliar_vector[slave_counter] = rConstantVector[i];
                ++slave_counter;
            }
            rRelationMatrix.resize(slave_size, master_size, false);
            rConstantVector.resize(slave_size, false);

            noalias(rRelationMatrix) = auxiliar_matrix;
            noalias(rConstantVector) = auxiliar_vector;
        }

        // We got some auxiliar values
        const GeometryType& r_slave_geometry = this->GetParentGeometry();
        const GeometryType& r_master_geometry = this->GetPairedGeometry();

        // Create dof list
        /* MASTER */
        BaseType::DofsVectorType auxiliar_master_dof_vector, master_dof_vector;

        master_dof_vector.resize(0);
        master_dof_vector.reserve(TDim * TNumNodesMaster);// + TDim * TNumNodes);

        if (master_size < potential_master_size) {
            auxiliar_master_dof_vector.resize(0);
            auxiliar_master_dof_vector.reserve(TDim * TNumNodes);

            for (auto& r_node_master : r_master_geometry) {
                auxiliar_master_dof_vector.push_back(r_node_master.pGetDof(DISPLACEMENT_X));
                auxiliar_master_dof_vector.push_back(r_node_master.pGetDof(DISPLACEMENT_Y));
                if (TDim == 3) auxiliar_master_dof_vector.push_back(r_node_master.pGetDof(DISPLACEMENT_Z));
            }

            master_dof_vector.resize(0);
            master_dof_vector.reserve(master_size);
            for (IndexType j : master_dofs_OK) {
                master_dof_vector.push_back(auxiliar_master_dof_vector[j]);
            }
        } else {
            for (auto& r_node_master : r_master_geometry) {
                master_dof_vector.push_back(r_node_master.pGetDof(DISPLACEMENT_X));
                master_dof_vector.push_back(r_node_master.pGetDof(DISPLACEMENT_Y));
                if (TDim == 3) master_dof_vector.push_back(r_node_master.pGetDof(DISPLACEMENT_Z));
            }
//             for (auto& r_node_slave : r_slave_geometry) {
//                 master_dof_vector.push_back(r_node_slave.pGetDof(DISPLACEMENT_X));
//                 master_dof_vector.push_back(r_node_slave.pGetDof(DISPLACEMENT_Y));
//                 if (TDim == 3) master_dof_vector.push_back(r_node_slave.pGetDof(DISPLACEMENT_Z));
//             }
        }

        /* SLAVE */
        BaseType::DofsVectorType auxiliar_slave_dof_vector, slave_dof_vector;

        if (slave_size < potential_slave_size) {
            auxiliar_slave_dof_vector.resize(0);
            auxiliar_slave_dof_vector.reserve(TDim * TNumNodes);

            for (auto& r_node_slave : r_slave_geometry) {
                auxiliar_slave_dof_vector.push_back(r_node_slave.pGetDof(DISPLACEMENT_X));
                auxiliar_slave_dof_vector.push_back(r_node_slave.pGetDof(DISPLACEMENT_Y));
                if (TDim == 3) auxiliar_slave_dof_vector.push_back(r_node_slave.pGetDof(DISPLACEMENT_Z));
            }

            slave_dof_vector.resize(0);
            slave_dof_vector.reserve(slave_size);
            for (IndexType i : slave_dofs_OK) {
                slave_dof_vector.push_back(auxiliar_slave_dof_vector[i]);
            }
        } else {
            slave_dof_vector.resize(0);
            slave_dof_vector.reserve(TDim * TNumNodes);

            for (auto& r_node_slave : r_slave_geometry) {
                slave_dof_vector.push_back(r_node_slave.pGetDof(DISPLACEMENT_X));
                slave_dof_vector.push_back(r_node_slave.pGetDof(DISPLACEMENT_Y));
                if (TDim == 3) slave_dof_vector.push_back(r_node_slave.pGetDof(DISPLACEMENT_Z));
            }
        }

        p_const->SetDofList(slave_dof_vector, master_dof_vector, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::UpdateConstraintFrictionless(
    MortarConditionMatrices& rMortarConditionMatrices,
    Matrix& rRelationMatrix,
    Vector& rConstantVector,
    ProcessInfo& rCurrentProcessInfo,
    const bool DualLM
    )
{
    // We got some auxiliar values
    const GeometryType& r_slave_geometry = this->GetParentGeometry();
    const GeometryType& r_master_geometry = this->GetPairedGeometry();

    // Current coordinates
    const BoundedMatrix<double, TNumNodes, TDim> x1 = MortarUtilities::GetCoordinates<TDim,TNumNodes>(r_slave_geometry);
//     const BoundedMatrix<double, TNumNodes, TDim> x1_0 = MortarUtilities::GetCoordinates<TDim,TNumNodes>(r_slave_geometry, false);
    const BoundedMatrix<double, TNumNodes, TDim> u1_0 = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(r_slave_geometry, DISPLACEMENT, 1);
    const BoundedMatrix<double, TNumNodesMaster, TDim> x2 = MortarUtilities::GetCoordinates<TDim,TNumNodesMaster>(r_master_geometry);
//     const BoundedMatrix<double, TNumNodesMaster, TDim> x2_0 = MortarUtilities::GetCoordinates<TDim,TNumNodes>(r_master_geometry, false);
    const BoundedMatrix<double, TNumNodesMaster, TDim> u2_0 = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(r_master_geometry, DISPLACEMENT, 1);

    // Mortar condition matrices - DOperator and MOperator
    const BoundedMatrix<double, TNumNodes, TNumNodes>& r_DOperator = rMortarConditionMatrices.DOperator;
    const BoundedMatrix<double, TNumNodes, TNumNodesMaster>& r_MOperator = rMortarConditionMatrices.MOperator;

    // Inverse DOperator
    BoundedMatrix<double, TNumNodes, TNumNodes> inverse_DOperator;
    if (DualLM) { // If dual LM inverse is trivial
        for (IndexType i = 0; i < TNumNodes; ++i) {
            for (IndexType j = 0; j < TNumNodes; ++j) {
                if (i == j) inverse_DOperator(i, i) = 1.0/r_DOperator(i, i);
                else inverse_DOperator(i, j) = 0.0;
            }
        }
    } else { // Otherwise regular inverse
        double det;
        MathUtils<double>::InvertMatrix(r_DOperator, inverse_DOperator, det, -1.0);
    }

    // Multiplying in order to obtain the proper relation matrix
//     const BoundedMatrix<double, TNumNodes, TNumNodesMaster> D_inv_M = prod(inverse_DOperator, r_MOperator);
    BoundedMatrix<double, TNumNodes, TNumNodesMaster> D_inv_M = prod(inverse_DOperator, r_MOperator);
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        for (IndexType j_node = 0; j_node < TNumNodesMaster; ++j_node) {
            if (D_inv_M(i_node, j_node) < -1.0e-8) { // If there is a negative value we deactivate all the row (not in the same section)
                for (IndexType j_node = 0; j_node < TNumNodesMaster; ++j_node) {
                    D_inv_M(i_node, j_node) = 0.0;
                }
                break;
            }
        }
    }

    // We create the relation matrix
    // TODO: Check this!!!
//     Vector lumping_factors;
//     lumping_factors = r_slave_geometry.LumpingFactors(lumping_factors);
//     const double domain_size = r_slave_geometry.DomainSize();
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        const bool node_is_active = r_slave_geometry[i_node].IsDefined(ACTIVE) ? r_slave_geometry[i_node].Is(ACTIVE) : true;
        if (node_is_active) {
//             const double weight_coeff = (lumping_factors[i_node] * domain_size)/r_slave_geometry[i_node].GetValue(NODAL_MAUX);
            double weight_coeff = 1.0/r_slave_geometry[i_node].GetValue(NODAL_PAUX);
            const array_1d<double, 3>& r_normal = r_slave_geometry[i_node].FastGetSolutionStepValue(NORMAL);
            for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
                for (IndexType j_node = 0; j_node < TNumNodesMaster; ++j_node) {
                    for (IndexType j_dim = 0; j_dim < TDim; ++j_dim) {
//                         const double coeff = r_normal[i_dim] * r_normal[j_dim] * aux_D_inv_M(i_node, j_node);
                        const double coeff = r_normal[i_dim] * r_normal[j_dim] * D_inv_M(i_node, j_node);
                        if (std::abs(coeff) > 1.0e-12) rRelationMatrix(i_node * TDim + i_dim, j_node * TDim + j_dim) = weight_coeff * coeff;
                    }
                }
            }
        }
//         else {
//             for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
//                 rRelationMatrix(i_node * TDim + i_dim, (TNumNodesMaster + i_node) * TDim + i_dim) = 1.0;
//             }
//         }
    }

    // We create the constact vector
    // The gap matrix and the normal matrix
    array_1d<double, TNumNodes> aux_gap_array;
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        const double nodal_area = r_slave_geometry[i_node].Has(NODAL_AREA) ? r_slave_geometry[i_node].GetValue(NODAL_AREA) : 1.0;
        aux_gap_array[i_node] = r_slave_geometry[i_node].FastGetSolutionStepValue(WEIGHTED_GAP)/nodal_area;
    }

    // Manual "Initial" gap by previous displacement
    const BoundedMatrix<double, TNumNodes, TDim> u10_D_inv_M_u20 = u1_0 - prod(D_inv_M, u2_0);
    array_1d<double, TDim> aux_normal, aux_array;
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        const array_1d<double, 3>& r_normal = r_slave_geometry[i_node].FastGetSolutionStepValue(NORMAL);
        for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
            aux_normal[i_dim] = r_normal[i_dim];
        }
        noalias(aux_array) = row(u10_D_inv_M_u20, i_node);
        aux_gap_array[i_node] += inner_prod(aux_array, aux_normal);
    }

//     // Manual "Initial" gap by position
//     const BoundedMatrix<double, TNumNodes, TDim> x10_D_inv_M_x20 = x1_0 - prod(D_inv_M, x2_0);
//     array_1d<double, TDim> aux_normal, aux_array;
//     for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
//         const array_1d<double, 3>& r_normal = r_slave_geometry[i_node].FastGetSolutionStepValue(NORMAL);
//         for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
//             aux_normal[i_dim] = r_normal[i_dim];
//         }
//         noalias(aux_array) = row(x10_D_inv_M_x20, i_node);
//         aux_gap_array[i_node] += inner_prod(aux_array, aux_normal);
//     }

    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
//         const double weight_coeff = (lumping_factors[i_node] * domain_size)/r_slave_geometry[i_node].GetValue(NODAL_MAUX);
        const double weight_coeff = 1.0/r_slave_geometry[i_node].GetValue(NODAL_PAUX);
        const array_1d<double, 3>& r_normal_slave_node = r_slave_geometry[i_node].FastGetSolutionStepValue(NORMAL);
        for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
            rConstantVector[i_node * TDim + i_dim] = weight_coeff * aux_gap_array[i_node] * r_normal_slave_node[i_dim];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::UpdateConstraintFrictional(
    MortarConditionMatrices& rMortarConditionMatrices,
    Matrix& rRelationMatrix,
    Vector& rConstantVector,
    ProcessInfo& rCurrentProcessInfo,
    const bool DualLM
    )
{
    // We got some auxiliar values
    const GeometryType& r_slave_geometry = this->GetParentGeometry();
    const GeometryType& r_master_geometry = this->GetPairedGeometry();

    // Current coordinates
    const BoundedMatrix<double, TNumNodes, TDim> x1 = MortarUtilities::GetCoordinates<TDim,TNumNodes>(r_slave_geometry);
//     const BoundedMatrix<double, TNumNodes, TDim> x1_0 = MortarUtilities::GetCoordinates<TDim,TNumNodes>(r_slave_geometry, false);
    const BoundedMatrix<double, TNumNodes, TDim> u1_0 = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(r_slave_geometry, DISPLACEMENT, 1);
    const BoundedMatrix<double, TNumNodesMaster, TDim> x2 = MortarUtilities::GetCoordinates<TDim,TNumNodesMaster>(r_master_geometry);
//     const BoundedMatrix<double, TNumNodesMaster, TDim> x2_0 = MortarUtilities::GetCoordinates<TDim,TNumNodes>(r_master_geometry, false);
    const BoundedMatrix<double, TNumNodesMaster, TDim> u2_0 = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(r_master_geometry, DISPLACEMENT, 1);

    // Mortar condition matrices - DOperator and MOperator
    const BoundedMatrix<double, TNumNodes, TNumNodes>& r_DOperator = rMortarConditionMatrices.DOperator;
    const BoundedMatrix<double, TNumNodes, TNumNodesMaster>& r_MOperator = rMortarConditionMatrices.MOperator;

    // Inverse DOperator
    BoundedMatrix<double, TNumNodes, TNumNodes> inverse_DOperator;
    if (DualLM) { // If dual LM inverse is trivial
        for (IndexType i = 0; i < TNumNodes; ++i) {
            for (IndexType j = 0; j < TNumNodes; ++j) {
                if (i == j) inverse_DOperator(i, i) = 1.0/r_DOperator(i, i);
                else inverse_DOperator(i, j) = 0.0;
            }
        }
    } else { // Otherwise regular inverse
        double det;
        MathUtils<double>::InvertMatrix(r_DOperator, inverse_DOperator, det, -1.0);
    }

    // Multiplying in order to obtain the proper relation matrix
//     const BoundedMatrix<double, TNumNodes, TNumNodesMaster> D_inv_M = prod(inverse_DOperator, r_MOperator);
    BoundedMatrix<double, TNumNodes, TNumNodesMaster> D_inv_M = prod(inverse_DOperator, r_MOperator);
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        for (IndexType j_node = 0; j_node < TNumNodesMaster; ++j_node) {
            if (D_inv_M(i_node, j_node) < -1.0e-8) { // If there is a negative value we deactivate all the row (not in the same section)
                for (IndexType j_node = 0; j_node < TNumNodesMaster; ++j_node) {
                    D_inv_M(i_node, j_node) = 0.0;
                }
                break;
            }
        }
    }

    // We create the relation matrix
//     Vector lumping_factors;
//     lumping_factors = r_slave_geometry.LumpingFactors(lumping_factors);
//     const double domain_size = r_slave_geometry.DomainSize();
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        const bool node_is_active = r_slave_geometry[i_node].IsDefined(ACTIVE) ? r_slave_geometry[i_node].Is(ACTIVE) : true;
        if (node_is_active) {
//             const double weight_coeff = (lumping_factors[i_node] * domain_size)/r_slave_geometry[i_node].GetValue(NODAL_MAUX);
            double weight_coeff = 1.0/r_slave_geometry[i_node].GetValue(NODAL_PAUX);
            const bool is_slip = r_slave_geometry[i_node].IsDefined(SLIP) ? r_slave_geometry[i_node].Is(SLIP) : false;
            if (is_slip) { // SLIP // TODO: Add nodal forces
                const array_1d<double, 3>& r_normal = r_slave_geometry[i_node].FastGetSolutionStepValue(NORMAL);
                for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
                    for (IndexType j_node = 0; j_node < TNumNodesMaster; ++j_node) {
                        for (IndexType j_dim = 0; j_dim < TDim; ++j_dim) {
                            const double coeff = r_normal[i_dim] * r_normal[j_dim] * D_inv_M(i_node, j_node);
                            if (std::abs(coeff) > 1.0e-12) rRelationMatrix(i_node * TDim + i_dim, j_node * TDim + j_dim) = weight_coeff * coeff;
                        }
                    }
                }
            } else { // STICK // TODO: ADD the contribution of slip to constant vector
                for (IndexType j_node = 0; j_node < TNumNodesMaster; ++j_node) {
                    const double coeff = D_inv_M(i_node, j_node);
                    if (std::abs(coeff) > 1.0e-12) {
                        for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
                            rRelationMatrix(i_node * TDim + i_dim, j_node * TDim + i_dim) = weight_coeff * coeff;
                        }
                    }
                }
            }
        }
    }

    // We create the constact vector
    // The gap matrix and the normal matrix
    array_1d<double, TNumNodes> aux_gap_array;
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        const double nodal_area = r_slave_geometry[i_node].Has(NODAL_AREA) ? r_slave_geometry[i_node].GetValue(NODAL_AREA) : 1.0;
        aux_gap_array[i_node] = r_slave_geometry[i_node].FastGetSolutionStepValue(WEIGHTED_GAP)/nodal_area;
    }

    // Manual "Initial" gap by previous displacement
    const BoundedMatrix<double, TNumNodes, TDim> u10_D_inv_M_u20 = u1_0 - prod(D_inv_M, u2_0);
    array_1d<double, TDim> aux_normal, aux_array;
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        const array_1d<double, 3>& r_normal = r_slave_geometry[i_node].FastGetSolutionStepValue(NORMAL);
        for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
            aux_normal[i_dim] = r_normal[i_dim];
        }
        noalias(aux_array) = row(u10_D_inv_M_u20, i_node);
        aux_gap_array[i_node] += inner_prod(aux_array, aux_normal);
    }

    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
//         const double weight_coeff = (lumping_factors[i_node] * domain_size)/r_slave_geometry[i_node].GetValue(NODAL_MAUX);
        const double weight_coeff = 1.0/r_slave_geometry[i_node].GetValue(NODAL_PAUX);
        const array_1d<double, 3>& r_normal_slave_node = r_slave_geometry[i_node].FastGetSolutionStepValue(NORMAL);
        for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
            rConstantVector[i_node * TDim + i_dim] = weight_coeff * aux_gap_array[i_node] * r_normal_slave_node[i_dim];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::UpdateConstraintTying(
    MortarConditionMatrices& rMortarConditionMatrices,
    Matrix& rRelationMatrix,
    Vector& rConstantVector,
    ProcessInfo& rCurrentProcessInfo,
    const bool DualLM
    )
{
    // We got some auxiliar values
    const GeometryType& r_slave_geometry = this->GetParentGeometry();
    const GeometryType& r_master_geometry = this->GetPairedGeometry();

    // Current coordinates
    const BoundedMatrix<double, TNumNodes, TDim> x1 = MortarUtilities::GetCoordinates<TDim,TNumNodes>(r_slave_geometry);
//     const BoundedMatrix<double, TNumNodes, TDim> x1_0 = MortarUtilities::GetCoordinates<TDim,TNumNodes>(r_slave_geometry, false);
    const BoundedMatrix<double, TNumNodes, TDim> u1_0 = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(r_slave_geometry, DISPLACEMENT, 1);
    const BoundedMatrix<double, TNumNodesMaster, TDim> x2 = MortarUtilities::GetCoordinates<TDim,TNumNodesMaster>(r_master_geometry);
//     const BoundedMatrix<double, TNumNodesMaster, TDim> x2_0 = MortarUtilities::GetCoordinates<TDim,TNumNodes>(r_master_geometry, false);
    const BoundedMatrix<double, TNumNodesMaster, TDim> u2_0 = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(r_master_geometry, DISPLACEMENT, 1);

    // Mortar condition matrices - DOperator and MOperator
    const BoundedMatrix<double, TNumNodes, TNumNodes>& r_DOperator = rMortarConditionMatrices.DOperator;
    const BoundedMatrix<double, TNumNodes, TNumNodesMaster>& r_MOperator = rMortarConditionMatrices.MOperator;

    // Inverse DOperator
    BoundedMatrix<double, TNumNodes, TNumNodes> inverse_DOperator;
    if (DualLM) { // If dual LM inverse is trivial
        for (IndexType i = 0; i < TNumNodes; ++i) {
            for (IndexType j = 0; j < TNumNodes; ++j) {
                if (i == j) inverse_DOperator(i, i) = 1.0/r_DOperator(i, i);
                else inverse_DOperator(i, j) = 0.0;
            }
        }
    } else { // Otherwise regular inverse
        double det;
        MathUtils<double>::InvertMatrix(r_DOperator, inverse_DOperator, det, -1.0);
    }

    // Multiplying in order to obtain the proper relation matrix
//     const BoundedMatrix<double, TNumNodes, TNumNodesMaster> D_inv_M = prod(inverse_DOperator, r_MOperator);
    BoundedMatrix<double, TNumNodes, TNumNodesMaster> D_inv_M = prod(inverse_DOperator, r_MOperator);
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        for (IndexType j_node = 0; j_node < TNumNodesMaster; ++j_node) {
            if (D_inv_M(i_node, j_node) < -1.0e-8) { // If there is a negative value we deactivate all the row (not in the same section)
                for (IndexType j_node = 0; j_node < TNumNodesMaster; ++j_node) {
                    D_inv_M(i_node, j_node) = 0.0;
                }
                break;
            }
        }
    }

    // We create the relation matrix
//     Vector lumping_factors;
//     lumping_factors = r_slave_geometry.LumpingFactors(lumping_factors);
//     const double domain_size = r_slave_geometry.DomainSize();
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        const bool node_is_active = r_slave_geometry[i_node].IsDefined(ACTIVE) ? r_slave_geometry[i_node].Is(ACTIVE) : true;
        if (node_is_active) {
//             const double weight_coeff = (lumping_factors[i_node] * domain_size)/r_slave_geometry[i_node].GetValue(NODAL_MAUX);
            double weight_coeff = 1.0/r_slave_geometry[i_node].GetValue(NODAL_PAUX);
            for (IndexType j_node = 0; j_node < TNumNodesMaster; ++j_node) {
                const double coeff = D_inv_M(i_node, j_node);
                if (std::abs(coeff) > 1.0e-12) {
                    for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
                        rRelationMatrix(i_node * TDim + i_dim, j_node * TDim + i_dim) = weight_coeff * coeff;
                    }
                }
            }
        }
    }

    // We create the constact vector
    // The gap matrix and the normal matrix
    array_1d<double, TNumNodes> aux_gap_array;
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        const double nodal_area = r_slave_geometry[i_node].Has(NODAL_AREA) ? r_slave_geometry[i_node].GetValue(NODAL_AREA) : 1.0;
        aux_gap_array[i_node] = r_slave_geometry[i_node].FastGetSolutionStepValue(WEIGHTED_GAP)/nodal_area;
    }

    // Manual "Initial" gap by previous displacement
    const BoundedMatrix<double, TNumNodes, TDim> u10_D_inv_M_u20 = u1_0 - prod(D_inv_M, u2_0);
    array_1d<double, TDim> aux_normal, aux_array;
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        const array_1d<double, 3>& r_normal = r_slave_geometry[i_node].FastGetSolutionStepValue(NORMAL);
        for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
            aux_normal[i_dim] = r_normal[i_dim];
        }
        noalias(aux_array) = row(u10_D_inv_M_u20, i_node);
        aux_gap_array[i_node] += inner_prod(aux_array, aux_normal);
    }

    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
//         const double weight_coeff = (lumping_factors[i_node] * domain_size)/r_slave_geometry[i_node].GetValue(NODAL_MAUX);
        const double weight_coeff = 1.0/r_slave_geometry[i_node].GetValue(NODAL_PAUX);
        const array_1d<double, 3>& r_normal_slave_node = r_slave_geometry[i_node].FastGetSolutionStepValue(NORMAL);
        for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
            rConstantVector[i_node * TDim + i_dim] = weight_coeff * aux_gap_array[i_node] * r_normal_slave_node[i_dim];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::ComputePreviousMortarOperators( ProcessInfo& rCurrentProcessInfo)
{
    const IndexType integration_order = this->GetProperties().Has(INTEGRATION_ORDER_CONTACT) ? this->GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) : 2;
    MortarExplicitContributionUtilities<TDim, TNumNodes, FrictionalCase::FRICTIONAL_PENALTY, false, TNumNodesMaster>::ComputePreviousMortarOperators(this, rCurrentProcessInfo, mPreviousMortarOperators, integration_order, false);
}

//******************************* GET DOUBLE VALUE *********************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::GetValueOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//******************************* GET ARRAY_1D VALUE *******************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::GetValueOnIntegrationPoints(
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//******************************* GET VECTOR VALUE *********************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::GetValueOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType &integration_points = GetParentGeometry().IntegrationPoints();

    if ( rOutput.size() != integration_points.size() )
        rOutput.resize( integration_points.size() );

    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
        rOutput[point_number] = 0.0;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector< array_1d<double, 3 > >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType &integration_points = GetParentGeometry().IntegrationPoints();

    if ( rOutput.size() != integration_points.size() )
        rOutput.resize( integration_points.size() );

    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
        rOutput[point_number] = ZeroVector(3);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType &integration_points = GetParentGeometry().IntegrationPoints();

    if ( rOutput.size() != integration_points.size() )
        rOutput.resize( integration_points.size() );

    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
        rOutput[point_number] = ZeroVector(3);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
int MPCMortarContactCondition<TDim,TNumNodes,TNumNodesMaster>::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int ierr = Condition::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    KRATOS_ERROR_IF(this->GetParentGeometry().NumberOfGeometryParts() == 0) << "YOU HAVE NOT INITIALIZED THE PAIR GEOMETRY IN THE MeshTyingMortarCondition" << std::endl;

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
    KRATOS_CHECK_VARIABLE_KEY(WEIGHTED_GAP)
    KRATOS_CHECK_VARIABLE_KEY(NORMAL)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( IndexType i = 0; i < TNumNodes; i++ ) {
        Node<3> &rnode = this->GetParentGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(WEIGHTED_GAP,rnode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NORMAL,rnode)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode)
    }

    return ierr;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

// Cases
template class MPCMortarContactCondition<2, 2, 2>;
template class MPCMortarContactCondition<3, 3, 3>;
template class MPCMortarContactCondition<3, 4, 4>;
template class MPCMortarContactCondition<3, 3, 4>;
template class MPCMortarContactCondition<3, 4, 3>;

} // Namespace Kratos
