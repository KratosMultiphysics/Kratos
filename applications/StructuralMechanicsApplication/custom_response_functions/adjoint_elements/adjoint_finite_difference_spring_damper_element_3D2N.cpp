// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//


#include "adjoint_finite_difference_spring_damper_element_3D2N.h"
#include "custom_elements/spring_damper_element.hpp"
#include "structural_mechanics_application_variables.h"
#include "includes/checks.h"

namespace Kratos
{

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
    else if (this->Has(rDesignVariable) && (rDesignVariable == NODAL_ROTATIONAL_STIFFNESS || rDesignVariable == NODAL_DISPLACEMENT_STIFFNESS)) {
        if ((rOutput.size1() != dimension) || (rOutput.size2() != local_size)) {
                rOutput.resize(dimension, local_size, false);
        }

        // save original stiffness parameters
        const auto variable_value = this->pGetPrimalElement()->GetValue(rDesignVariable);

        // reset original stiffness parameters before computing the derivatives
        this->pGetPrimalElement()->SetValue(rDesignVariable, rDesignVariable.Zero());

        Vector RHS;
        for(IndexType dir_i = 0; dir_i < dimension; ++dir_i) {
            // The following approach assumes a linear dependency between RHS and spring stiffness
            array_1d<double, 3> perturbed_nodal_stiffness = ZeroVector(3);
            perturbed_nodal_stiffness[dir_i] = 1.0;
            this->pGetPrimalElement()->SetValue(rDesignVariable, perturbed_nodal_stiffness);
            this->pGetPrimalElement()->CalculateRightHandSide(RHS, rCurrentProcessInfo);

            KRATOS_ERROR_IF_NOT(RHS.size() == local_size) << "Size of the pseudo-load does not fit!" << std::endl;
            for(IndexType i = 0; i < RHS.size(); ++i) {
                rOutput(dir_i, i) = RHS[i];
            }
        }

        // give original stiffness parameters back
        this->pGetPrimalElement()->SetValue(rDesignVariable, variable_value);
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
int AdjointFiniteDifferenceSpringDamperElement<TPrimalElement>::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(this->pGetPrimalElement()) << "Primal element pointer is nullptr!" << std::endl;

    // Verify that the dofs exist
    for ( std::size_t i = 0; i < this->GetGeometry().size(); i++ ) {
        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        const NodeType& rnode = this->GetGeometry()[i];

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

template class AdjointFiniteDifferenceSpringDamperElement<SpringDamperElement<3>>;

} // namespace Kratos.


