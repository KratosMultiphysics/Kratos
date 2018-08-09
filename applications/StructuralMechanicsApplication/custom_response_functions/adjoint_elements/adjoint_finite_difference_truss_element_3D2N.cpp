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


namespace Kratos
{

AdjointFiniteDifferenceTrussElement::AdjointFiniteDifferenceTrussElement(Element::Pointer pPrimalElement)
    : AdjointFiniteDifferenceTrussElementLinear(pPrimalElement)
{
}

AdjointFiniteDifferenceTrussElement::~AdjointFiniteDifferenceTrussElement()
{
}

void AdjointFiniteDifferenceTrussElement::Calculate(const Variable<Vector >& rVariable,
                        Vector& rOutput,
                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    
    if(rVariable == STRESS_ON_GP)
    {
        TracedStressType traced_stress_type = static_cast<TracedStressType>(this->GetValue(TRACED_STRESS_TYPE));

        const SizeType  GP_num = (mpPrimalElement->GetGeometry().IntegrationPoints()).size();
        if (rOutput.size() != GP_num) 
            rOutput.resize(GP_num, false);

        switch (traced_stress_type)
        {
            case TracedStressType::FX:
            {
                std::vector< array_1d<double, 3 > > force_vector;
                mpPrimalElement->GetValueOnIntegrationPoints(FORCE, force_vector, rCurrentProcessInfo);
                for(IndexType i = 0; i < GP_num ; ++i)
                    rOutput(i) = force_vector[i][0];
                break;
            }
            case TracedStressType::PK2X:
            {
                std::vector<Vector> stress_vector;
                mpPrimalElement->GetValueOnIntegrationPoints(PK2_STRESS_VECTOR, stress_vector, rCurrentProcessInfo);
                for(IndexType i = 0; i < GP_num ; ++i)
                    rOutput(i) = stress_vector[i][0];
                break;
            }
            default:
                KRATOS_ERROR << "Invalid stress type! Stress type not supported for this element!" << std::endl;
        }
    }
    else
    {
        rOutput.resize(1);
        rOutput.clear();
    }

    KRATOS_CATCH("")
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
        rOutput.resize(num_dofs, num_GP);
        rOutput.clear();

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

double AdjointFiniteDifferenceTrussElement::GetPerturbationSizeModificationFactor(const Variable<array_1d<double,3>>& rDesignVariable)
{
    KRATOS_TRY;

    if(rDesignVariable == SHAPE)
        return this->CalculateReferenceLength();
    else
        return 1.0;

    KRATOS_CATCH("")
}

double AdjointFiniteDifferenceTrussElement::CalculateReferenceLength() 
{
    KRATOS_TRY;
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
    const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
    const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
    const double L = std::sqrt(dx * dx + dy * dy + dz * dz);

    KRATOS_ERROR_IF(L<=numerical_limit)
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

    KRATOS_ERROR_IF(l<=numerical_limit)
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

    KRATOS_ERROR_IF_NOT(num_dofs == 6) << "This implementation expects 6 dofs!" << std::endl;

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

    rDerivativeVector[0] = -1.0 * (u2 - u1 + dx) / l;
    rDerivativeVector[1] = -1.0 * (v2 - v1 + dy) / l;
    rDerivativeVector[2] = -1.0 * (w2 - w1 + dz) / l;
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
        case TracedStressType::PK2X:
        {
            rDerivativePreFactor = this->CalculateDerivativePreFactorPK2X(rCurrentProcessInfo);
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
    const double L0 = CalculateReferenceLength();
    const double l = CalculateCurrentLength();
    double prestress = 0.00;
    if (mpPrimalElement->GetProperties().Has(TRUSS_PRESTRESS_PK2))
        prestress = mpPrimalElement->GetProperties()[TRUSS_PRESTRESS_PK2];
    std::vector<Vector> GL_strain;  
    mpPrimalElement->CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, GL_strain, rCurrentProcessInfo);
    const double GL_strain_X = GL_strain[0][0]; //one Gauss-Point result is enough due to constant strains.
        
    double derivative_pre_factor = A / L0 * (E * GL_strain_X + prestress + E * l * l / (L0 * L0));

    KRATOS_ERROR_IF(derivative_pre_factor<=numerical_limit)
        << "Derivative pre-factor of " << this->Id() << "~ 0" << std::endl;

    return derivative_pre_factor;
}

double AdjointFiniteDifferenceTrussElement::CalculateDerivativePreFactorPK2X(const ProcessInfo& rCurrentProcessInfo)
{
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    const double E = mpPrimalElement->GetProperties()[YOUNG_MODULUS];
    const double l = CalculateCurrentLength();
    const double L0 = CalculateReferenceLength();
    double derivative_pre_factor = E * l / (L0 * L0);

    KRATOS_ERROR_IF(derivative_pre_factor<=numerical_limit)
        << "Derivative pre-factor of " << this->Id() << "~ 0" << std::endl;

    return derivative_pre_factor;
}

void AdjointFiniteDifferenceTrussElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, AdjointFiniteDifferenceTrussElementLinear);
}

void AdjointFiniteDifferenceTrussElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, AdjointFiniteDifferenceTrussElementLinear);
}

} // namespace Kratos.


