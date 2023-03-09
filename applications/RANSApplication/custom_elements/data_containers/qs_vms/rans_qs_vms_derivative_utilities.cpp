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
#include "custom_elements/data_containers/k_epsilon/k_element_data.h"
#include "custom_elements/data_containers/k_epsilon/epsilon_element_data.h"

// Include base h
#include "rans_qs_vms_derivative_utilities.h"

namespace Kratos
{
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

// template instantiations

template class RansQSVMSDerivativeUtilities<2>::TurbulenceVariableDerivative<3, KEpsilonElementData::KElementData<2>>;
template class RansQSVMSDerivativeUtilities<2>::TurbulenceVariableDerivative<3, KEpsilonElementData::EpsilonElementData<2>>;

template class RansQSVMSDerivativeUtilities<3>::TurbulenceVariableDerivative<4, KEpsilonElementData::KElementData<3>>;
template class RansQSVMSDerivativeUtilities<3>::TurbulenceVariableDerivative<4, KEpsilonElementData::EpsilonElementData<3>>;

template class RansQSVMSDerivativeUtilities<2>;
template class RansQSVMSDerivativeUtilities<3>;

} // namespace Kratos
