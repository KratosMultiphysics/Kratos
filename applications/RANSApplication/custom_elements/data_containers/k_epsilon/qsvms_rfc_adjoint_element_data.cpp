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
#include "includes/checks.h"
#include "geometries/geometry_data.h"
#include "includes/process_info.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "rans_application_variables.h"

// Include base h
#include "qsvms_rfc_adjoint_element_data.h"

namespace Kratos
{
namespace KEpsilonElementData
{

namespace {
template <unsigned int TDim>
std::vector<const Variable<double>*> GetDofsList();

template <unsigned int TDim>
std::vector<const Variable<double>*> GetSecondDerivativeVariablesList();

template <>
std::vector<const Variable<double>*> GetDofsList<2>()
{
    return {
        &ADJOINT_FLUID_VECTOR_1_X,
        &ADJOINT_FLUID_VECTOR_1_Y,
        &ADJOINT_FLUID_SCALAR_1,
        &RANS_SCALAR_1_ADJOINT_1,
        &RANS_SCALAR_2_ADJOINT_1};
}

template <>
std::vector<const Variable<double>*> GetSecondDerivativeVariablesList<2>()
{
    return {
        &RELAXED_ACCELERATION_X,
        &RELAXED_ACCELERATION_Y,
        &PRESSURE,
        &RANS_AUXILIARY_VARIABLE_1,
        &RANS_AUXILIARY_VARIABLE_2};
}

template <>
std::vector<const Variable<double>*> GetDofsList<3>()
{
    return {
        &ADJOINT_FLUID_VECTOR_1_X,
        &ADJOINT_FLUID_VECTOR_1_Y,
        &ADJOINT_FLUID_VECTOR_1_Z,
        &ADJOINT_FLUID_SCALAR_1,
        &RANS_SCALAR_1_ADJOINT_1,
        &RANS_SCALAR_2_ADJOINT_1};
}

template <>
std::vector<const Variable<double>*> GetSecondDerivativeVariablesList<3>()
{
    return {
        &RELAXED_ACCELERATION_X,
        &RELAXED_ACCELERATION_Y,
        &RELAXED_ACCELERATION_Z,
        &PRESSURE,
        &RANS_AUXILIARY_VARIABLE_1,
        &RANS_AUXILIARY_VARIABLE_2};
}
} // namespace

template <unsigned int TDim, unsigned int TNumNodes>
void QSVMSRFCAdjointElementData<TDim, TNumNodes>::Check(
    const Element& rElement,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    Fluid::Check(rElement, rProcessInfo);
    TurbulenceModelEquation1::TResidualsDerivatives::Check(rElement, rProcessInfo);
    TurbulenceModelEquation2::TResidualsDerivatives::Check(rElement, rProcessInfo);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
std::vector<const Variable<double>*> QSVMSRFCAdjointElementData<TDim, TNumNodes>::GetDofVariablesList()
{
    return GetDofsList<TDim>();
}

template <unsigned int TDim, unsigned int TNumNodes>
std::vector<const Variable<double>*> QSVMSRFCAdjointElementData<TDim, TNumNodes>::GetPrimalSecondDerivativeVariablesList()
{
    return GetSecondDerivativeVariablesList<TDim>();
}

// template instantiations

template class QSVMSRFCAdjointElementData<2, 3>;
template class QSVMSRFCAdjointElementData<3, 4>;

} // namespace KEpsilonElementData

} // namespace Kratos
