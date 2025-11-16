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

using DependentVariablesListType = std::vector<
                                        std::tuple<
                                            const Variable<double>&,
                                            std::vector<const Variable<double>*>
                                        >
                                    >;

template <unsigned int TComponentIndex>
const DependentVariablesListType GetEffectiveViscosityDependentVariablesForVelocity();

template <>
const DependentVariablesListType GetEffectiveViscosityDependentVariablesForVelocity<0>()
{
    return DependentVariablesListType{
        std::tuple_cat(
            std::tie(VELOCITY_X),
            std::make_tuple(std::vector<const Variable<double>*>{
                &VELOCITY_GRADIENT_TENSOR_XX,
                &VELOCITY_GRADIENT_TENSOR_XY,
                &VELOCITY_GRADIENT_TENSOR_XZ
            })
        )
    };
}

template <>
const DependentVariablesListType GetEffectiveViscosityDependentVariablesForVelocity<1>()
{
    return DependentVariablesListType{
        std::tuple_cat(
            std::tie(VELOCITY_Y),
            std::make_tuple(std::vector<const Variable<double>*>{
                &VELOCITY_GRADIENT_TENSOR_YX,
                &VELOCITY_GRADIENT_TENSOR_YY,
                &VELOCITY_GRADIENT_TENSOR_YZ
            })
        )
    };
}

template <>
const DependentVariablesListType GetEffectiveViscosityDependentVariablesForVelocity<2>()
{
    return DependentVariablesListType{
        std::tuple_cat(
            std::tie(VELOCITY_Z),
            std::make_tuple(std::vector<const Variable<double>*>{
                &VELOCITY_GRADIENT_TENSOR_ZX,
                &VELOCITY_GRADIENT_TENSOR_ZY,
                &VELOCITY_GRADIENT_TENSOR_ZZ
            })
        )
    };
}

template <unsigned int TDim>
const DependentVariablesListType GetEffectiveViscosityDependentVariablesForShape();

template <>
const DependentVariablesListType GetEffectiveViscosityDependentVariablesForShape<2>()
{
    return DependentVariablesListType{
        std::tuple_cat(
            std::tie(VELOCITY_X),
            std::make_tuple(std::vector<const Variable<double>*>{
                &VELOCITY_GRADIENT_TENSOR_XX,
                &VELOCITY_GRADIENT_TENSOR_XY
            })
        ),
        std::tuple_cat(
            std::tie(VELOCITY_Y),
            std::make_tuple(std::vector<const Variable<double>*>{
                &VELOCITY_GRADIENT_TENSOR_YX,
                &VELOCITY_GRADIENT_TENSOR_YY
            })
        )
    };
}

template <>
const DependentVariablesListType GetEffectiveViscosityDependentVariablesForShape<3>()
{
    return DependentVariablesListType{
        std::tuple_cat(
            std::tie(VELOCITY_X),
            std::make_tuple(std::vector<const Variable<double>*>{
                &VELOCITY_GRADIENT_TENSOR_XX,
                &VELOCITY_GRADIENT_TENSOR_XY,
                &VELOCITY_GRADIENT_TENSOR_XZ
            })
        ),
        std::tuple_cat(
            std::tie(VELOCITY_Y),
            std::make_tuple(std::vector<const Variable<double>*>{
                &VELOCITY_GRADIENT_TENSOR_YX,
                &VELOCITY_GRADIENT_TENSOR_YY,
                &VELOCITY_GRADIENT_TENSOR_YZ,
            })
        ),
        std::tuple_cat(
            std::tie(VELOCITY_Z),
            std::make_tuple(std::vector<const Variable<double>*>{
                &VELOCITY_GRADIENT_TENSOR_ZX,
                &VELOCITY_GRADIENT_TENSOR_ZY,
                &VELOCITY_GRADIENT_TENSOR_ZZ,
            })
        )
    };
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
template<unsigned int TNumNodes, unsigned int TComponentIndex>
typename RansQSVMSDerivativeUtilities<TDim>::DependentVariablesListType RansQSVMSDerivativeUtilities<TDim>::KOmegaSSTVelocityDerivative<TNumNodes, TComponentIndex>::GetEffectiveViscosityDependentVariables() const
{
    return RANSQSVMSDerivativeHelperUtilities::GetEffectiveViscosityDependentVariablesForVelocity<TComponentIndex>();
}

template <unsigned int TDim>
template<unsigned int TNumNodes, unsigned int TComponentIndex>
typename RansQSVMSDerivativeUtilities<TDim>::DependentVariablesListType RansQSVMSDerivativeUtilities<TDim>::KOmegaSSTShapeDerivative<TNumNodes, TComponentIndex>::GetEffectiveViscosityDependentVariables() const
{
    return RANSQSVMSDerivativeHelperUtilities::GetEffectiveViscosityDependentVariablesForShape<TDim>();
}

// template instantiations

// k-epsilon
template class RansQSVMSDerivativeUtilities<2>::TurbulenceVariableDerivative<3, KEpsilonElementData::KElementData<2>>;
template class RansQSVMSDerivativeUtilities<2>::TurbulenceVariableDerivative<3, KEpsilonElementData::EpsilonElementData<2>>;
template class RansQSVMSDerivativeUtilities<2>::TurbulenceVariableDerivative<4, KEpsilonElementData::KElementData<2>>;
template class RansQSVMSDerivativeUtilities<2>::TurbulenceVariableDerivative<4, KEpsilonElementData::EpsilonElementData<2>>;

template class RansQSVMSDerivativeUtilities<3>::TurbulenceVariableDerivative<4, KEpsilonElementData::KElementData<3>>;
template class RansQSVMSDerivativeUtilities<3>::TurbulenceVariableDerivative<4, KEpsilonElementData::EpsilonElementData<3>>;
template class RansQSVMSDerivativeUtilities<3>::TurbulenceVariableDerivative<8, KEpsilonElementData::KElementData<3>>;
template class RansQSVMSDerivativeUtilities<3>::TurbulenceVariableDerivative<8, KEpsilonElementData::EpsilonElementData<3>>;

// k-omega
template class RansQSVMSDerivativeUtilities<2>::TurbulenceVariableDerivative<3, KOmegaElementData::KElementData<2>>;
template class RansQSVMSDerivativeUtilities<2>::TurbulenceVariableDerivative<3, KOmegaElementData::OmegaElementData<2>>;
template class RansQSVMSDerivativeUtilities<2>::TurbulenceVariableDerivative<4, KOmegaElementData::KElementData<2>>;
template class RansQSVMSDerivativeUtilities<2>::TurbulenceVariableDerivative<4, KOmegaElementData::OmegaElementData<2>>;

template class RansQSVMSDerivativeUtilities<3>::TurbulenceVariableDerivative<4, KOmegaElementData::KElementData<3>>;
template class RansQSVMSDerivativeUtilities<3>::TurbulenceVariableDerivative<4, KOmegaElementData::OmegaElementData<3>>;
template class RansQSVMSDerivativeUtilities<3>::TurbulenceVariableDerivative<8, KOmegaElementData::KElementData<3>>;
template class RansQSVMSDerivativeUtilities<3>::TurbulenceVariableDerivative<8, KOmegaElementData::OmegaElementData<3>>;

// k-omega-sst
template class RansQSVMSDerivativeUtilities<2>::TurbulenceVariableDerivative<3, KOmegaSSTElementData::KElementData<2>>;
template class RansQSVMSDerivativeUtilities<2>::TurbulenceVariableDerivative<3, KOmegaSSTElementData::OmegaElementData<2>>;
template class RansQSVMSDerivativeUtilities<2>::TurbulenceVariableDerivative<4, KOmegaSSTElementData::KElementData<2>>;
template class RansQSVMSDerivativeUtilities<2>::TurbulenceVariableDerivative<4, KOmegaSSTElementData::OmegaElementData<2>>;

template class RansQSVMSDerivativeUtilities<3>::TurbulenceVariableDerivative<4, KOmegaSSTElementData::KElementData<3>>;
template class RansQSVMSDerivativeUtilities<3>::TurbulenceVariableDerivative<4, KOmegaSSTElementData::OmegaElementData<3>>;
template class RansQSVMSDerivativeUtilities<3>::TurbulenceVariableDerivative<8, KOmegaSSTElementData::KElementData<3>>;
template class RansQSVMSDerivativeUtilities<3>::TurbulenceVariableDerivative<8, KOmegaSSTElementData::OmegaElementData<3>>;

// k-omega-sst specific velocity derivatives
template class RansQSVMSDerivativeUtilities<2>::KOmegaSSTVelocityDerivative<3, 0>;
template class RansQSVMSDerivativeUtilities<2>::KOmegaSSTVelocityDerivative<3, 1>;
template class RansQSVMSDerivativeUtilities<2>::KOmegaSSTVelocityDerivative<4, 0>;
template class RansQSVMSDerivativeUtilities<2>::KOmegaSSTVelocityDerivative<4, 1>;

template class RansQSVMSDerivativeUtilities<3>::KOmegaSSTVelocityDerivative<4, 0>;
template class RansQSVMSDerivativeUtilities<3>::KOmegaSSTVelocityDerivative<4, 1>;
template class RansQSVMSDerivativeUtilities<3>::KOmegaSSTVelocityDerivative<4, 2>;
template class RansQSVMSDerivativeUtilities<3>::KOmegaSSTVelocityDerivative<8, 0>;
template class RansQSVMSDerivativeUtilities<3>::KOmegaSSTVelocityDerivative<8, 1>;
template class RansQSVMSDerivativeUtilities<3>::KOmegaSSTVelocityDerivative<8, 2>;

// k-omega-sst specific shape derivatives
template class RansQSVMSDerivativeUtilities<2>::KOmegaSSTShapeDerivative<3, 0>;
template class RansQSVMSDerivativeUtilities<2>::KOmegaSSTShapeDerivative<3, 1>;
template class RansQSVMSDerivativeUtilities<2>::KOmegaSSTShapeDerivative<4, 0>;
template class RansQSVMSDerivativeUtilities<2>::KOmegaSSTShapeDerivative<4, 1>;

template class RansQSVMSDerivativeUtilities<3>::KOmegaSSTShapeDerivative<4, 0>;
template class RansQSVMSDerivativeUtilities<3>::KOmegaSSTShapeDerivative<4, 1>;
template class RansQSVMSDerivativeUtilities<3>::KOmegaSSTShapeDerivative<4, 2>;
template class RansQSVMSDerivativeUtilities<3>::KOmegaSSTShapeDerivative<8, 0>;
template class RansQSVMSDerivativeUtilities<3>::KOmegaSSTShapeDerivative<8, 1>;
template class RansQSVMSDerivativeUtilities<3>::KOmegaSSTShapeDerivative<8, 2>;

template class RansQSVMSDerivativeUtilities<2>;
template class RansQSVMSDerivativeUtilities<3>;

} // namespace Kratos
