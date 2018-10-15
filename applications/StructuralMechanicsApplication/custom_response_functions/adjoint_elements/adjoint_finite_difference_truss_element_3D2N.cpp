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


#include "adjoint_finite_difference_truss_element_3D2N.h"
#include "structural_mechanics_application_variables.h"
#include "custom_response_functions/response_utilities/stress_response_definitions.h"
#include "includes/checks.h"


namespace Kratos
{

AdjointFiniteDifferenceTrussElement::AdjointFiniteDifferenceTrussElement(Element::Pointer pPrimalElement)
    : AdjointFiniteDifferencingBaseElement(pPrimalElement)
{
}

AdjointFiniteDifferenceTrussElement::~AdjointFiniteDifferenceTrussElement()
{
}

void AdjointFiniteDifferenceTrussElement::CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
                                    Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if(rStressVariable == STRESS_ON_GP)
    {
        const SizeType num_nodes = mpPrimalElement->GetGeometry().PointsNumber();
        const SizeType dimension = mpPrimalElement->GetGeometry().WorkingSpaceDimension();
        const SizeType num_dofs = num_nodes * dimension;
        const SizeType num_GP = (mpPrimalElement->GetGeometry().IntegrationPoints()).size();
        if ( (rOutput.size1() != num_dofs) || (rOutput.size2() != num_GP ) )
            rOutput.resize(num_dofs, num_GP);

        double derivative_pre_factor;
        this->GetDerivativePreFactor(derivative_pre_factor, rCurrentProcessInfo);

        Vector length_derivative_vector;
        this->CalculateCurrentLengthDisplacementDerivative(length_derivative_vector);

        for(IndexType i = 0; i < num_dofs; ++i)
        {
             for(IndexType j = 0; j < num_GP; ++j)
                rOutput(i, j) = length_derivative_vector[i] * derivative_pre_factor;
        }
    }
    else
        KRATOS_ERROR << "Stress displacement derivative only available for Gauss-points quantities!" << std::endl;

    KRATOS_CATCH("")
}

int AdjointFiniteDifferenceTrussElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    int return_value = AdjointFiniteDifferencingBaseElement::Check(rCurrentProcessInfo);

    KRATOS_ERROR_IF_NOT(mpPrimalElement) << "Primal element pointer is nullptr!" << std::endl;

    //TODO: Check() of primal element should be called, but is not possible because of DOF check!

    KRATOS_ERROR_IF(this->GetGeometry().WorkingSpaceDimension() != 3 || this->GetGeometry().size() != 2)
        << "The truss element works only in 3D and with 2 noded elements" << "" << std::endl;

    CheckVariables();
    CheckDofs();
    CheckProperties(rCurrentProcessInfo);

    KRATOS_ERROR_IF(this->GetGeometry().Length() < std::numeric_limits<double>::epsilon())
        << "Element #" << this->Id() << " has a length of zero!" << std::endl;

    return return_value;

    KRATOS_CATCH("")
}


void AdjointFiniteDifferenceTrussElement::CheckVariables()
{
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_CHECK_VARIABLE_KEY(CROSS_AREA);
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_DISPLACEMENT);
}

void AdjointFiniteDifferenceTrussElement::CheckDofs()
{
    const SizeType number_of_nodes = this->GetGeometry().size();

    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        const NodeType &r_node = this->GetGeometry()[i];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_DISPLACEMENT, r_node);

        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Y, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Z, r_node);
    }
}

void AdjointFiniteDifferenceTrussElement::CheckProperties(const ProcessInfo& rCurrentProcessInfo)
{
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    const PropertiesType & r_properties = GetProperties();

    KRATOS_ERROR_IF(r_properties.Has(CROSS_AREA) == false || r_properties[CROSS_AREA] <= numerical_limit)
    << "CROSS_AREA not provided for this element" << this->Id() << std::endl;
    KRATOS_ERROR_IF(r_properties.Has(YOUNG_MODULUS) == false || r_properties[YOUNG_MODULUS] <= numerical_limit)
    << "YOUNG_MODULUS not provided for this element" << this->Id() << std::endl;
    KRATOS_ERROR_IF_NOT( r_properties.Has(DENSITY) )
    << "DENSITY not provided for this element" << this->Id() << std::endl;

    KRATOS_ERROR_IF_NOT(r_properties.Has(CONSTITUTIVE_LAW))
    << "CONSTITUTIVE_LAW not provided for element " << this->Id() << std::endl;
    const ConstitutiveLaw::Pointer& cl = r_properties[CONSTITUTIVE_LAW];
    KRATOS_ERROR_IF(cl == nullptr)
    << "CONSTITUTIVE_LAW not provided for element " << this->Id() << std::endl;
    cl->Check(r_properties ,this->GetGeometry(),rCurrentProcessInfo);
}

double AdjointFiniteDifferenceTrussElement::CalculateReferenceLength()
{
  KRATOS_TRY;
  const double numerical_limit = std::numeric_limits<double>::epsilon();
  const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
  const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
  const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
  const double L = std::sqrt(dx * dx + dy * dy + dz * dz);

  KRATOS_DEBUG_ERROR_IF(L <= numerical_limit)
   << "Reference Length of element" << this->Id() << "~ 0" << std::endl;
  return L;
  KRATOS_CATCH("")
}

double AdjointFiniteDifferenceTrussElement::CalculateCurrentLength()
{
    KRATOS_TRY;
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    const double du =
        this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X) -
        this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
    const double dv =
        this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y) -
        this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
    const double dw =
        this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Z) -
        this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z);
    const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
    const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
    const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
    const double l = std::sqrt((du + dx) * (du + dx) + (dv + dy) * (dv + dy) +
                               (dw + dz) * (dw + dz));

    KRATOS_DEBUG_ERROR_IF(l <= numerical_limit)
        << "Current Length of element" << this->Id() << "~ 0" << std::endl;
    return l;
    KRATOS_CATCH("")
}

void AdjointFiniteDifferenceTrussElement::CalculateCurrentLengthDisplacementDerivative(Vector& rDerivativeVector)
{
    KRATOS_TRY;

    const SizeType number_of_nodes = this->GetGeometry().size();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();
    const SizeType num_dofs = number_of_nodes * dimension;

    if (rDerivativeVector.size() != num_dofs)
        rDerivativeVector.resize(num_dofs, false);

    const double l = CalculateCurrentLength();
    const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
    const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
    const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
    const double u1 = this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
    const double u2 = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X);
    const double v1 = this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
    const double v2 = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y);
    const double w1 = this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z);
    const double w2 = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Z);

    rDerivativeVector[0] = (u1 - u2 - dx) / l;
    rDerivativeVector[1] = (v1 - v2 - dy) / l;
    rDerivativeVector[2] = (w1 - w2 - dz) / l;
    rDerivativeVector[3] = -1.0 * rDerivativeVector[0];
    rDerivativeVector[4] = -1.0 * rDerivativeVector[1];
    rDerivativeVector[5] = -1.0 * rDerivativeVector[2];

    KRATOS_CATCH("")
}

void AdjointFiniteDifferenceTrussElement::GetDerivativePreFactor(double& rDerivativePreFactor, const ProcessInfo& rCurrentProcessInfo)
{
    TracedStressType traced_stress_type = static_cast<TracedStressType>(this->GetValue(TRACED_STRESS_TYPE));

    switch (traced_stress_type)
    {
        case TracedStressType::FX:
        {
            rDerivativePreFactor = this->CalculateDerivativePreFactorFX(rCurrentProcessInfo);
            break;
        }
        case TracedStressType::PK2:
        {
            rDerivativePreFactor = this->CalculateDerivativePreFactorPK2(rCurrentProcessInfo);
            break;
        }
        default:
            KRATOS_ERROR << "Invalid stress type! Stress type not supported for this element!" << std::endl;
    }
}

double AdjointFiniteDifferenceTrussElement::CalculateDerivativePreFactorFX(const ProcessInfo& rCurrentProcessInfo)
{
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    const double E = mpPrimalElement->GetProperties()[YOUNG_MODULUS];
    const double A = mpPrimalElement->GetProperties()[CROSS_AREA];
    const double l_0 = CalculateReferenceLength();
    const double l = CalculateCurrentLength();
    double prestress = 0.0;
    if (mpPrimalElement->GetProperties().Has(TRUSS_PRESTRESS_PK2))
        prestress = mpPrimalElement->GetProperties()[TRUSS_PRESTRESS_PK2];
    std::vector<Vector> GL_strain;
    mpPrimalElement->CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, GL_strain, rCurrentProcessInfo);
    const double GL_strain_X = GL_strain[0][0]; //one Gauss-Point result is enough due to constant strains.

    double derivative_pre_factor = A / l_0 * (E * GL_strain_X + prestress + E * l * l / (l_0 * l_0));

    KRATOS_DEBUG_ERROR_IF(std::abs(derivative_pre_factor) <= numerical_limit)
        << "Derivative pre-factor of " << this->Id() << "~ 0" << std::endl;

    return derivative_pre_factor;
}

double AdjointFiniteDifferenceTrussElement::CalculateDerivativePreFactorPK2(const ProcessInfo& rCurrentProcessInfo)
{
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    const double E = mpPrimalElement->GetProperties()[YOUNG_MODULUS];
    const double l = CalculateCurrentLength();
    const double l_0 = CalculateReferenceLength();
    double derivative_pre_factor = E * l / (l_0 * l_0);

    KRATOS_DEBUG_ERROR_IF(derivative_pre_factor<=numerical_limit)
        << "Derivative pre-factor of " << this->Id() << "~ 0" << std::endl;

    return derivative_pre_factor;
}

void AdjointFiniteDifferenceTrussElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, AdjointFiniteDifferencingBaseElement);
}

void AdjointFiniteDifferenceTrussElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, AdjointFiniteDifferencingBaseElement);
}

} // namespace Kratos.


