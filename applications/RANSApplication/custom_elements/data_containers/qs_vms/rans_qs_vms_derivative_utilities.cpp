//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "utilities/element_size_calculator.h"

// Application includes
// k-epsilon
#include "custom_elements/data_containers/k_epsilon/k_element_data.h"
#include "custom_elements/data_containers/k_epsilon/epsilon_element_data.h"
// k-omega
#include "custom_elements/data_containers/k_omega/k_element_data.h"
#include "custom_elements/data_containers/k_omega/omega_element_data.h"
// k-omega-sst
#include "custom_elements/data_containers/k_omega_sst/k_element_data.h"
#include "custom_elements/data_containers/k_omega_sst/omega_element_data.h"


// Include base h
#include "rans_qs_vms_derivative_utilities.h"

namespace Kratos
{
namespace RANSQSVMSDerivativeHelperUtilities
{
template <unsigned int TDim>
using AdjointVariableInformationType = AdjointVariableInformation::VariableInformation<TDim>;

template <unsigned int TDim>
std::vector<AdjointVariableInformationType<TDim>> GetEffectiveViscosityDependentVariables();

template <>
std::vector<AdjointVariableInformationType<2>> GetEffectiveViscosityDependentVariables<2>()
{
    return std::vector<AdjointVariableInformationType<2>>{
        AdjointVariableInformation::VelocityInformation<2>(0),
        AdjointVariableInformation::VelocityInformation<2>(1),
    };
}

template <>
std::vector<AdjointVariableInformationType<3>> GetEffectiveViscosityDependentVariables<3>()
{
    return std::vector<AdjointVariableInformationType<3>>{
        AdjointVariableInformation::VelocityInformation<3>(0),
        AdjointVariableInformation::VelocityInformation<3>(1),
        AdjointVariableInformation::VelocityInformation<3>(2)};
}
} // namespace RANSQSVMSDerivativeHelperUtilities

template <unsigned int TDim>
template <unsigned int TNumNodes, class TElementData>
array_1d<double, TDim> RansQSVMSDerivativeUtilities<TDim>::TurbulenceVariableDerivative<TNumNodes, TElementData>::CalculateEffectiveVelocityDerivative(
    const array_1d<double, TDim>& rVelocity) const
{
    array_1d<double, TDim> result = ZeroVector(TDim);
    return result;
}

template <unsigned int TDim>
template <unsigned int TNumNodes, class TElementData>
const Variable<double>& RansQSVMSDerivativeUtilities<TDim>::TurbulenceVariableDerivative<TNumNodes, TElementData>::GetDerivativeVariable() const
{
    return TElementData::GetScalarVariable();
}

template <unsigned int TDim>
template <unsigned int TNumNodes, class TElementData>
double RansQSVMSDerivativeUtilities<TDim>::TurbulenceVariableDerivative<TNumNodes, TElementData>::CalculateElementLengthDerivative(
    const double ElementLength) const
{
    return 0.0;
}

template <unsigned int TDim>
template <unsigned int TNumNodes, class TElementData>
void RansQSVMSDerivativeUtilities<TDim>::TurbulenceVariableDerivative<TNumNodes, TElementData>::CalculateStrainRateDerivative(
    Vector& rOutput,
    const Matrix& rNodalVelocity) const
{
    rOutput.clear();
}

template <unsigned int TDim>
template <unsigned int TNumNodes>
std::vector<AdjointVariableInformation::VariableInformation<TDim>> RansQSVMSDerivativeUtilities<TDim>::KOmegaSSTVelocityDerivative<TNumNodes>::GetEffectiveViscosityDependentVariables() const
{
    return RANSQSVMSDerivativeHelperUtilities::GetEffectiveViscosityDependentVariables<TDim>();
}

template <unsigned int TDim>
template <unsigned int TNumNodes>
std::vector<AdjointVariableInformation::VariableInformation<TDim>> RansQSVMSDerivativeUtilities<TDim>::KOmegaSSTShapeDerivative<TNumNodes>::GetEffectiveViscosityDependentVariables() const
{
    return RANSQSVMSDerivativeHelperUtilities::GetEffectiveViscosityDependentVariables<TDim>();
}

// template instantiations

// k-epsilon
template class RansQSVMSDerivativeUtilities<2>::TurbulenceVariableDerivative<3, KEpsilonElementData::KElementData<2>>;
template class RansQSVMSDerivativeUtilities<2>::TurbulenceVariableDerivative<3, KEpsilonElementData::EpsilonElementData<2>>;

template class RansQSVMSDerivativeUtilities<3>::TurbulenceVariableDerivative<4, KEpsilonElementData::KElementData<3>>;
template class RansQSVMSDerivativeUtilities<3>::TurbulenceVariableDerivative<4, KEpsilonElementData::EpsilonElementData<3>>;

// k-omega
template class RansQSVMSDerivativeUtilities<2>::TurbulenceVariableDerivative<3, KOmegaElementData::KElementData<2>>;
template class RansQSVMSDerivativeUtilities<2>::TurbulenceVariableDerivative<3, KOmegaElementData::OmegaElementData<2>>;

template class RansQSVMSDerivativeUtilities<3>::TurbulenceVariableDerivative<4, KOmegaElementData::KElementData<3>>;
template class RansQSVMSDerivativeUtilities<3>::TurbulenceVariableDerivative<4, KOmegaElementData::OmegaElementData<3>>;

// k-omega-sst
template class RansQSVMSDerivativeUtilities<2>::TurbulenceVariableDerivative<3, KOmegaSSTElementData::KElementData<2>>;
template class RansQSVMSDerivativeUtilities<2>::TurbulenceVariableDerivative<3, KOmegaSSTElementData::OmegaElementData<2>>;

template class RansQSVMSDerivativeUtilities<3>::TurbulenceVariableDerivative<4, KOmegaSSTElementData::KElementData<3>>;
template class RansQSVMSDerivativeUtilities<3>::TurbulenceVariableDerivative<4, KOmegaSSTElementData::OmegaElementData<3>>;

// k-omega-sst specific velocity derivatives
template class RansQSVMSDerivativeUtilities<2>::KOmegaSSTVelocityDerivative<3>;
template class RansQSVMSDerivativeUtilities<3>::KOmegaSSTVelocityDerivative<4>;

// k-omega-sst specific shape derivatives
template class RansQSVMSDerivativeUtilities<2>::KOmegaSSTShapeDerivative<3>;
template class RansQSVMSDerivativeUtilities<3>::KOmegaSSTShapeDerivative<4>;

template class RansQSVMSDerivativeUtilities<2>;
template class RansQSVMSDerivativeUtilities<3>;

} // namespace Kratos
