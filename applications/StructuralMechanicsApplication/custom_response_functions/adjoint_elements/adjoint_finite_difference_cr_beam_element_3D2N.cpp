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


namespace Kratos
{

AdjointFiniteDifferenceCrBeamElement::AdjointFiniteDifferenceCrBeamElement(Element::Pointer pPrimalElement)
    : AdjointFiniteDifferencingBaseElement(pPrimalElement, true)
{
}

AdjointFiniteDifferenceCrBeamElement::~AdjointFiniteDifferenceCrBeamElement()
{
}

void AdjointFiniteDifferenceCrBeamElement::Calculate(const Variable<Vector >& rVariable,
                        Vector& rOutput,
                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if(rVariable == STRESS_ON_GP || rVariable == STRESS_ON_NODE)
    {
        TracedStressType traced_stress_type = static_cast<TracedStressType>(this->GetValue(TRACED_STRESS_TYPE));

        std::vector< array_1d<double, 3 > > stress_vector;
        int direction_1 = 0;
        bool stress_is_moment = true;

        switch (traced_stress_type)
        {
            case TracedStressType::MX:
            {
                direction_1 = 0;
                break;
            }
            case TracedStressType::MY:
            {
                direction_1 = 1;
                break;
            }
            case TracedStressType::MZ:
            {
                direction_1 = 2;
                break;
            }
            case TracedStressType::FX:
            {
                direction_1 = 0;
                stress_is_moment = false;
                break;
            }
            case TracedStressType::FY:
            {
                direction_1 = 1;
                stress_is_moment = false;
                break;
            }
            case TracedStressType::FZ:
            {
                direction_1 = 2;
                stress_is_moment = false;
                break;
            }
            default:
                KRATOS_ERROR << "Invalid stress type! Stress type not supported for this element!" << std::endl;
        }

        if(stress_is_moment)
            mpPrimalElement->CalculateOnIntegrationPoints(MOMENT, stress_vector, rCurrentProcessInfo);
        else
            mpPrimalElement->CalculateOnIntegrationPoints(FORCE, stress_vector, rCurrentProcessInfo);

        if(rVariable == STRESS_ON_GP)
        {
            const SizeType  GP_num = GetGeometry().IntegrationPointsNumber(Kratos::GeometryData::GI_GAUSS_3);

            rOutput.resize(GP_num);
            for(IndexType i = 0; i < GP_num ; i++)
            {
                rOutput(i) = stress_vector[i][direction_1];
            }
        }
        else if(rVariable == STRESS_ON_NODE)
        {
            rOutput.resize(2);
            rOutput(0) = 2 * stress_vector[0][direction_1] - stress_vector[1][direction_1];
            rOutput(1) = 2 * stress_vector[2][direction_1] - stress_vector[1][direction_1];
        }
    }
    else
    {
        rOutput.resize(3);
        rOutput.clear();
    }

    KRATOS_CATCH("")
}

int AdjointFiniteDifferenceCrBeamElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    int return_value = AdjointFiniteDifferencingBaseElement::Check(rCurrentProcessInfo);

    KRATOS_ERROR_IF_NOT(mpPrimalElement) << "Primal element pointer is nullptr!" << std::endl;

    //TODO: Check() of primal element should be called, but is not possible because of DOF check!

    KRATOS_ERROR_IF(GetGeometry().WorkingSpaceDimension() != 3 || GetGeometry().size() != 2)
    << "The beam element works only in 3D and with 2 noded elements" << "" << std::endl;

    KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_CHECK_VARIABLE_KEY(CROSS_AREA);
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_ROTATION);

    // check dofs
    GeometryType& r_geom = GetGeometry();
    for (IndexType i = 0; i < r_geom.size(); i++)
    {
        auto& r_node = r_geom[i];

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

void AdjointFiniteDifferenceCrBeamElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, AdjointFiniteDifferencingBaseElement);
}

void AdjointFiniteDifferenceCrBeamElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, AdjointFiniteDifferencingBaseElement);
}

} // namespace Kratos.


