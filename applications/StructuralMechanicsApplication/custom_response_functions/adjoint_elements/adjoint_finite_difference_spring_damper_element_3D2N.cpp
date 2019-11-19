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
    this->pGetPrimalElement()->SetValue(NODAL_DISPLACEMENT_STIFFNESS, this->GetValue(NODAL_DISPLACEMENT_STIFFNESS));
    this->pGetPrimalElement()->SetValue(NODAL_ROTATIONAL_STIFFNESS, this->GetValue(NODAL_ROTATIONAL_STIFFNESS));
    BaseType::InitializeSolutionStep(rCurrentProcessInfo);
}

template <class TPrimalElement>
void AdjointFiniteDifferenceSpringDamperElement<TPrimalElement>::CalculateSensitivityMatrix(
                                            const Variable<double>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
{


    if(rDesignVariable == SHAPE_SENSITIVITY) {
        BaseType::CalculateSensitivityMatrix(rDesignVariable, rOutput, rCurrentProcessInfo);
    }
    else {
        std::cout << "CalculateSensitivityMatrix of AdjointFiniteDifferenceSpringDamperElement is called!" << std::endl;
    }

}

template <class TPrimalElement>
void AdjointFiniteDifferenceSpringDamperElement<TPrimalElement>::CalculateSensitivityMatrix(
                                            const Variable<array_1d<double,3>>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    std::cout << "CalculateSensitivityMatrix of AdjointFiniteDifferenceSpringDamperElement is called!" << std::endl;
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

    KRATOS_CATCH( "Problem in the Check in the SpringDamperElement3D2N" )
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


