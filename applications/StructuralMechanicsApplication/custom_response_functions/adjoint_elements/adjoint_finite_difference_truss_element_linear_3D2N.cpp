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
#include "custom_elements/truss_element_linear_3D2N.hpp"


namespace Kratos
{

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
                                    Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    TracedStressType traced_stress_type = static_cast<TracedStressType>(this->GetValue(TRACED_STRESS_TYPE));

    if(traced_stress_type  == TracedStressType::FX)
    {
        // ensure that adjoint load is determined without influence of pre-stress
        // pre-stress does not cancel out when computing this derivative with unit-displacements!
        Properties::Pointer p_global_properties = this->mpPrimalElement->pGetProperties();

        Properties::Pointer p_local_property(Kratos::make_shared<Properties>(*p_global_properties));
        this->mpPrimalElement->SetProperties(p_local_property);

        p_local_property->SetValue(TRUSS_PRESTRESS_PK2, 0.0);

        AdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateStressDisplacementDerivative(rStressVariable,
                                           rOutput, rCurrentProcessInfo);

        this->mpPrimalElement->SetProperties(p_global_properties);
    }
    else
        AdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateStressDisplacementDerivative(rStressVariable,
                                   rOutput, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

template class AdjointFiniteDifferenceTrussElementLinear<TrussElementLinear3D2N>;

} // namespace Kratos.


