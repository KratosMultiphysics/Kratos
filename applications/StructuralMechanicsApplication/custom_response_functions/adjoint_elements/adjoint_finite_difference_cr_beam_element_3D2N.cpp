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


#include "adjoint_finite_difference_cr_beam_element_3D2N.h"
#include "structural_mechanics_application_variables.h"
#include "custom_response_functions/response_utilities/stress_response_definitions.h"
#include "includes/checks.h"
#include "custom_elements/cr_beam_element_linear_3D2N.hpp"


namespace Kratos
{

template <class TPrimalElement>
void AdjointFiniteDifferenceCrBeamElement<TPrimalElement>::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
					      std::vector< array_1d<double, 3 > >& rOutput,
					      const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == ADJOINT_CURVATURE || rVariable == ADJOINT_STRAIN) {
        const double E = this->GetProperties()[YOUNG_MODULUS];
        const double nu = this->GetProperties()[POISSON_RATIO];
        const double G = E / (2.0 * (1.0 + nu));
        const double A = this->GetProperties()[CROSS_AREA];
        const double J = this->GetProperties()[TORSIONAL_INERTIA];
        const double Iy = this->GetProperties()[I22];
        const double Iz = this->GetProperties()[I33];

        if (rVariable == ADJOINT_CURVATURE) {
            this->CalculateAdjointFieldOnIntegrationPoints(MOMENT, rOutput, rCurrentProcessInfo);

            for (IndexType i = 0; i < rOutput.size(); ++i) {
                rOutput[i][0] *=  1.0 / (G * J);
                rOutput[i][1] *= -1.0 / (E * Iy);
                rOutput[i][2] *= -1.0 / (E * Iz);
            }
        } else if (rVariable == ADJOINT_STRAIN) {
            this->CalculateAdjointFieldOnIntegrationPoints(FORCE, rOutput, rCurrentProcessInfo);

            KRATOS_WARNING_IF("ADJOINT_STRAIN", (this->GetProperties().Has(AREA_EFFECTIVE_Y) || this->GetProperties().Has(AREA_EFFECTIVE_Z)))
                        << "Not available for Timoschenko beam!" << std::endl;

            for (IndexType i = 0; i < rOutput.size(); ++i) {
                rOutput[i][0] *= 1.0 / (E * A);
                rOutput[i][1] *= 0.0;
                rOutput[i][2] *= 0.0;
            }
        }
    } else {
        this->CalculateAdjointFieldOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }

    KRATOS_CATCH("")
}

template <class TPrimalElement>
int AdjointFiniteDifferenceCrBeamElement<TPrimalElement>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    int return_value = BaseType::Check(rCurrentProcessInfo);

    KRATOS_ERROR_IF_NOT(this->mHasRotationDofs) << "Adjoint beam element does not have rotation dofs!" << std::endl;
    KRATOS_ERROR_IF_NOT(this->mpPrimalElement) << "Primal element pointer is nullptr!" << std::endl;

    //TODO: Check() of primal element should be called, but is not possible because of DOF check!

    KRATOS_ERROR_IF(this->GetGeometry().WorkingSpaceDimension() != 3 || this->GetGeometry().size() != 2)
    << "The beam element works only in 3D and with 2 noded elements" << "" << std::endl;

    KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_CHECK_VARIABLE_KEY(CROSS_AREA);
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_ROTATION);

    // check dofs
    const GeometryType& r_geom = this->GetGeometry();
    for (IndexType i = 0; i < r_geom.size(); i++)
    {
        const auto& r_node = r_geom[i];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ROTATION, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_DISPLACEMENT, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_ROTATION, r_node);

        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Y, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Z, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_Y, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_Z, r_node);
    }

    const double numerical_limit = std::numeric_limits<double>::epsilon();

    KRATOS_ERROR_IF(this->GetProperties().Has(CROSS_AREA) == false || this->GetProperties()[CROSS_AREA] <= numerical_limit)
    << "CROSS_AREA not provided for this element" << this->Id() << std::endl;
    KRATOS_ERROR_IF(this->GetProperties().Has(YOUNG_MODULUS) == false || this->GetProperties()[YOUNG_MODULUS] <= numerical_limit)
    << "YOUNG_MODULUS not provided for this element" << this->Id() << std::endl;
    KRATOS_ERROR_IF_NOT( this->GetProperties().Has(DENSITY) )
    << "DENSITY not provided for this element" << this->Id() << std::endl;
    KRATOS_ERROR_IF_NOT( this->GetProperties().Has(POISSON_RATIO) )
    << "POISSON_RATIO not provided for this element" << this->Id() << std::endl;
    KRATOS_ERROR_IF_NOT( this->GetProperties().Has(TORSIONAL_INERTIA) )
    << "TORSIONAL_INERTIA not provided for this element" << this->Id() << std::endl;

    KRATOS_ERROR_IF_NOT( this->GetProperties().Has(I22) )
    << "I22 not provided for this element" << this->Id() << std::endl;
    KRATOS_ERROR_IF_NOT( this->GetProperties().Has(I33) )
    << "I33 not provided for this element" << this->Id() << std::endl;

    return return_value;

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferenceCrBeamElement<TPrimalElement>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

template <class TPrimalElement>
void AdjointFiniteDifferenceCrBeamElement<TPrimalElement>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

template class AdjointFiniteDifferenceCrBeamElement<CrBeamElementLinear3D2N>;

} // namespace Kratos.


