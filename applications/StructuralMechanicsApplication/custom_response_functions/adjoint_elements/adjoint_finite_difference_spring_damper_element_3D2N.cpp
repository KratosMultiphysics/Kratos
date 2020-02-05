// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//


#include "adjoint_finite_difference_spring_damper_element_3D2N.h"
#include "custom_elements/spring_damper_element_3D2N.hpp"
#include "structural_mechanics_application_variables.h"
#include "includes/checks.h"

namespace Kratos
{

template <class TPrimalElement>
void AdjointFiniteDifferenceSpringDamperElement<TPrimalElement>::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(this->Has(NODAL_DISPLACEMENT_STIFFNESS) || this->Has(NODAL_ROTATIONAL_STIFFNESS)) <<
        "Neither NODAL_DISPLACEMENT_STIFFNESS nor NODAL_ROTATIONAL_STIFFNESS available!" << std::endl;

    // As the stiffness parameters are saved in the non-historical database of the element these parameters
    // have to be explicitly transfered from the adjoint to the primal element. Please note: if the stiffness parameters
    // would be part of the element properties this transferring would be not necessary. The stiffness parameters
    // are needed by the primal element in order to compute later on the element stiffness matrix for the
    // adjoint problem and the sensitivity matrix as the element contribution to the pseudo-load.
    const auto& r_const_this = *this;
    this->pGetPrimalElement()->SetValue(NODAL_DISPLACEMENT_STIFFNESS, r_const_this.GetValue(NODAL_DISPLACEMENT_STIFFNESS));
    this->pGetPrimalElement()->SetValue(NODAL_ROTATIONAL_STIFFNESS, r_const_this.GetValue(NODAL_ROTATIONAL_STIFFNESS));

    BaseType::InitializeSolutionStep(rCurrentProcessInfo);

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferenceSpringDamperElement<TPrimalElement>::CalculateSensitivityMatrix(
                                            const Variable<double>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const SizeType number_of_nodes = this->GetGeometry().PointsNumber();
    const SizeType dimension = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);
    const SizeType num_dofs_per_node = dimension * 2;
    const SizeType local_size = number_of_nodes * num_dofs_per_node;
    if ((rOutput.size1() != 0) || (rOutput.size2() != local_size)) {
            rOutput.resize(0, local_size, false);
    }
    noalias(rOutput) = ZeroMatrix(0, local_size);

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferenceSpringDamperElement<TPrimalElement>::CalculateSensitivityMatrix(
                                            const Variable<array_1d<double,3>>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const SizeType number_of_nodes = this->GetGeometry().PointsNumber();
    const SizeType dimension = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);
    const SizeType num_dofs_per_node = dimension * 2;
    const SizeType local_size = number_of_nodes * num_dofs_per_node;

    if(rDesignVariable == SHAPE_SENSITIVITY) {
        if ((rOutput.size1() != dimension*number_of_nodes) || (rOutput.size2() != local_size)) {
                rOutput.resize(dimension*number_of_nodes, local_size, false);
        }
        noalias(rOutput) = ZeroMatrix(dimension*number_of_nodes, local_size);
    }
    else if (rDesignVariable == NODAL_ROTATIONAL_STIFFNESS || rDesignVariable == NODAL_DISPLACEMENT_STIFFNESS ) {
        if ((rOutput.size1() != dimension) || (rOutput.size2() != local_size)) {
                rOutput.resize(dimension, local_size, false);
        }

        // reset original stiffness parameters before computing the derivatives
        this->pGetPrimalElement()->SetValue(NODAL_ROTATIONAL_STIFFNESS, rDesignVariable.Zero());
        this->pGetPrimalElement()->SetValue(NODAL_DISPLACEMENT_STIFFNESS, rDesignVariable.Zero());

        ProcessInfo process_info = rCurrentProcessInfo;
        Vector RHS;
        for(IndexType dir_i = 0; dir_i < dimension; ++dir_i) {
            // The following approach assumes a linear dependency between RHS and spring stiffness
            array_1d<double, 3> perturbed_nodal_stiffness = ZeroVector(3);
            perturbed_nodal_stiffness[dir_i] = 1.0;
            this->pGetPrimalElement()->SetValue(rDesignVariable, perturbed_nodal_stiffness);
            this->pGetPrimalElement()->CalculateRightHandSide(RHS, process_info);

            KRATOS_ERROR_IF_NOT(RHS.size() == local_size) << "Size of the pseudo-load does not fit!" << std::endl;
            for(IndexType i = 0; i < RHS.size(); ++i) {
                rOutput(dir_i, i) = RHS[i];
            }
        }

        // give original stiffness parameters back.
        const auto& r_const_this = *this;
        this->pGetPrimalElement()->SetValue(NODAL_DISPLACEMENT_STIFFNESS, r_const_this.GetValue(NODAL_DISPLACEMENT_STIFFNESS));
        this->pGetPrimalElement()->SetValue(NODAL_ROTATIONAL_STIFFNESS, r_const_this.GetValue(NODAL_ROTATIONAL_STIFFNESS));
    }
    else {
        if ((rOutput.size1() != 0) || (rOutput.size2() != local_size)) {
            rOutput.resize(0, local_size, false);
        }
        noalias(rOutput) = ZeroMatrix(0, local_size);
    }

    KRATOS_CATCH("")
}

template <class TPrimalElement>
int AdjointFiniteDifferenceSpringDamperElement<TPrimalElement>::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(this->pGetPrimalElement()) << "Primal element pointer is nullptr!" << std::endl;

    // Check that all required variables have been registered

    // The displacement terms
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(NODAL_MASS)
    KRATOS_CHECK_VARIABLE_KEY(NODAL_DISPLACEMENT_STIFFNESS)
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_DISPLACEMENT)

    // The rotational terms
    KRATOS_CHECK_VARIABLE_KEY(ROTATION)
    KRATOS_CHECK_VARIABLE_KEY(ANGULAR_VELOCITY)
    KRATOS_CHECK_VARIABLE_KEY(ANGULAR_ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(NODAL_INERTIA)
    KRATOS_CHECK_VARIABLE_KEY(NODAL_ROTATIONAL_STIFFNESS)

    KRATOS_CHECK_VARIABLE_KEY(NODAL_DAMPING_RATIO)
    KRATOS_CHECK_VARIABLE_KEY(NODAL_ROTATIONAL_DAMPING_RATIO)
    KRATOS_CHECK_VARIABLE_KEY(VOLUME_ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_ROTATION)

    // Verify that the dofs exist
    for ( std::size_t i = 0; i < this->GetGeometry().size(); i++ ) {
        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        NodeType& rnode = this->GetGeometry()[i];

        // The displacement terms
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_DISPLACEMENT,rnode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VOLUME_ACCELERATION,rnode)

        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_X,rnode)
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Y,rnode)
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Z,rnode)

        // The rotational terms
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_ROTATION,rnode)

        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_X,rnode)
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_Y,rnode)
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_Z,rnode)
    }

    return 0;

    KRATOS_CATCH( "Problem in the Check in the AdjointFiniteDifferenceSpringDamperElement!" )
}

template <class TPrimalElement>
void AdjointFiniteDifferenceSpringDamperElement<TPrimalElement>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

template <class TPrimalElement>
void AdjointFiniteDifferenceSpringDamperElement<TPrimalElement>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

template class AdjointFiniteDifferenceSpringDamperElement<SpringDamperElement3D2N>;

} // namespace Kratos.


