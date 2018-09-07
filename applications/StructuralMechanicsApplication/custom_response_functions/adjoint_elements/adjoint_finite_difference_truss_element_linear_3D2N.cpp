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


#include "adjoint_finite_difference_truss_element_linear_3D2N.h"
#include "structural_mechanics_application_variables.h"
#include "custom_response_functions/response_utilities/stress_response_definitions.h"


namespace Kratos
{

AdjointFiniteDifferenceTrussElementLinear::AdjointFiniteDifferenceTrussElementLinear(Element::Pointer pPrimalElement)
    : AdjointFiniteDifferenceTrussElement(pPrimalElement)
{
}

AdjointFiniteDifferenceTrussElementLinear::~AdjointFiniteDifferenceTrussElementLinear()
{
}

void AdjointFiniteDifferenceTrussElementLinear::Calculate(const Variable<Vector >& rVariable,
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
                mpPrimalElement->CalculateOnIntegrationPoints(FORCE, force_vector, rCurrentProcessInfo);
                for(IndexType i = 0; i < GP_num ; ++i)
                    rOutput(i) = force_vector[i][0];
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

void AdjointFiniteDifferenceTrussElementLinear::CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
                                    Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    TracedStressType traced_stress_type = static_cast<TracedStressType>(this->GetValue(TRACED_STRESS_TYPE));

    if(traced_stress_type  == TracedStressType::FX)
    {
        // ensure that adjoint load is determined without influence of pre-stress
        // pre-stress does not cancel out when computing this derivative with unit-displacements!
        Properties::Pointer p_global_properties = mpPrimalElement->pGetProperties();

        Properties::Pointer p_local_property(Kratos::make_shared<Properties>(*p_global_properties));
        mpPrimalElement->SetProperties(p_local_property);

        p_local_property->SetValue(TRUSS_PRESTRESS_PK2, 0.0);

        AdjointFiniteDifferencingBaseElement::CalculateStressDisplacementDerivative(rStressVariable,
                                           rOutput, rCurrentProcessInfo);

        mpPrimalElement->SetProperties(p_global_properties);
    }
    else
        AdjointFiniteDifferencingBaseElement::CalculateStressDisplacementDerivative(rStressVariable,
                                   rOutput, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

void AdjointFiniteDifferenceTrussElementLinear::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, AdjointFiniteDifferenceTrussElement);
}

void AdjointFiniteDifferenceTrussElementLinear::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, AdjointFiniteDifferenceTrussElement);
}

} // namespace Kratos.


