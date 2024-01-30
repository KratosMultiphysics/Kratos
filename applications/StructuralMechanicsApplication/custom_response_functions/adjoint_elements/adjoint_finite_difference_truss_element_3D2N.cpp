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


#include "adjoint_finite_difference_truss_element_3D2N.h"
#include "structural_mechanics_application_variables.h"
#include "custom_response_functions/response_utilities/stress_response_definitions.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "includes/checks.h"
#include "custom_elements/truss_element_3D2N.hpp"
#include "custom_elements/truss_element_linear_3D2N.hpp"




namespace Kratos
{

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElement<TPrimalElement>::CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
                                    Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if(rStressVariable == STRESS_ON_GP)
    {
        const SizeType num_nodes = this->mpPrimalElement->GetGeometry().PointsNumber();
        const SizeType dimension = this->mpPrimalElement->GetGeometry().WorkingSpaceDimension();
        const SizeType num_dofs = num_nodes * dimension;
        const SizeType num_GP = (this->mpPrimalElement->GetGeometry().IntegrationPoints()).size();
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

template <class TPrimalElement>
int AdjointFiniteDifferenceTrussElement<TPrimalElement>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    int return_value = BaseType::Check(rCurrentProcessInfo);

    KRATOS_ERROR_IF_NOT(this->mpPrimalElement) << "Primal element pointer is nullptr!" << std::endl;

    //TODO: Check() of primal element should be called, but is not possible because of DOF check!

    KRATOS_ERROR_IF(this->GetGeometry().WorkingSpaceDimension() != 3 || this->GetGeometry().size() != 2)
        << "The truss element works only in 3D and with 2 noded elements" << "" << std::endl;

    this->CheckDofs();
    this->CheckProperties(rCurrentProcessInfo);

    KRATOS_ERROR_IF(StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this)
         < std::numeric_limits<double>::epsilon())
        << "Element #" << this->Id() << " has a length of zero!" << std::endl;

    return return_value;

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElement<TPrimalElement>::CheckDofs() const
{
    const SizeType number_of_nodes = this->GetGeometry().size();

    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        const auto& r_node = this->GetGeometry()[i];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_DISPLACEMENT, r_node);

        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Y, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Z, r_node);
    }
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElement<TPrimalElement>::CheckProperties(const ProcessInfo& rCurrentProcessInfo) const
{
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    const PropertiesType & r_properties = this->GetProperties();

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

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElement<TPrimalElement>::CalculateCurrentLengthDisplacementDerivative(Vector& rDerivativeVector)
{
    KRATOS_TRY;

    const SizeType number_of_nodes = this->GetGeometry().size();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();
    const SizeType num_dofs = number_of_nodes * dimension;

    if (rDerivativeVector.size() != num_dofs)
        rDerivativeVector.resize(num_dofs, false);

    const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
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

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElement<TPrimalElement>::GetDerivativePreFactor(double& rDerivativePreFactor, const ProcessInfo& rCurrentProcessInfo)
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

template <class TPrimalElement>
double AdjointFiniteDifferenceTrussElement<TPrimalElement>::CalculateDerivativePreFactorFX(const ProcessInfo& rCurrentProcessInfo)
{
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    const double E = this->mpPrimalElement->GetProperties()[YOUNG_MODULUS];
    const double A = this->mpPrimalElement->GetProperties()[CROSS_AREA];
    const double l_0 = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
    double prestress = 0.0;
    if (this->mpPrimalElement->GetProperties().Has(TRUSS_PRESTRESS_PK2))
        prestress = this->mpPrimalElement->GetProperties()[TRUSS_PRESTRESS_PK2];
    std::vector<Vector> GL_strain;
    this->mpPrimalElement->CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, GL_strain, rCurrentProcessInfo);
    const double GL_strain_X = GL_strain[0][0]; //one Gauss-Point result is enough due to constant strains.

    double derivative_pre_factor = A / l_0 * (E * GL_strain_X + prestress + E * l * l / (l_0 * l_0));

    KRATOS_DEBUG_ERROR_IF(std::abs(derivative_pre_factor) <= numerical_limit)
        << "Derivative pre-factor of " << this->Id() << "~ 0" << std::endl;

    return derivative_pre_factor;
}

template <class TPrimalElement>
double AdjointFiniteDifferenceTrussElement<TPrimalElement>::CalculateDerivativePreFactorPK2(const ProcessInfo& rCurrentProcessInfo)
{
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    const double E = this->mpPrimalElement->GetProperties()[YOUNG_MODULUS];
    const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
    const double l_0 = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    double derivative_pre_factor = E * l / (l_0 * l_0);

    KRATOS_DEBUG_ERROR_IF(derivative_pre_factor<=numerical_limit)
        << "Derivative pre-factor of " << this->Id() << "~ 0" << std::endl;

    return derivative_pre_factor;
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElement<TPrimalElement>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElement<TPrimalElement>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

template class AdjointFiniteDifferenceTrussElement<TrussElement3D2N>;
template class AdjointFiniteDifferenceTrussElement<TrussElementLinear3D2N>;

} // namespace Kratos.


