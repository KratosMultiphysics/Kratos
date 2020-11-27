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
#include "geometries/geometry_data.h"
#include "includes/process_info.h"

// Application includes
#include "fluid_dynamics_application_variables.h"

// Include base h
#include "qs_vms_adjoint_element_data.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <unsigned int TDim, unsigned int TNumNodes>
int QSVMSAdjointElementData<TDim, TNumNodes>::Check(
    const GeometryType& rGeometry,
    const ProcessInfo& rProcessInfo)
{
    // return TResidualDerivatives::Check(rGeometry, rProcessInfo);
    return 0;
}

template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod QSVMSAdjointElementData<TDim, TNumNodes>::GetIntegrationMethod()
{
    return TResidualDerivatives::GetIntegrationMethod();
}

template <>
std::vector<const Variable<double>*> QSVMSAdjointElementData<2, 3>::GetDofVariablesList()
{
    return {&ADJOINT_FLUID_VECTOR_1_X, &ADJOINT_FLUID_VECTOR_1_Y, &ADJOINT_FLUID_SCALAR_1};
}

template <>
std::vector<const Variable<double>*> QSVMSAdjointElementData<2, 4>::GetDofVariablesList()
{
    return {&ADJOINT_FLUID_VECTOR_1_X, &ADJOINT_FLUID_VECTOR_1_Y, &ADJOINT_FLUID_SCALAR_1};
}

template <>
std::vector<const Variable<double>*> QSVMSAdjointElementData<3, 4>::GetDofVariablesList()
{
    return {&ADJOINT_FLUID_VECTOR_1_X, &ADJOINT_FLUID_VECTOR_1_Y, &ADJOINT_FLUID_VECTOR_1_Z, &ADJOINT_FLUID_SCALAR_1};
}

template <>
std::vector<const Variable<double>*> QSVMSAdjointElementData<3, 8>::GetDofVariablesList()
{
    return {&ADJOINT_FLUID_VECTOR_1_X, &ADJOINT_FLUID_VECTOR_1_Y, &ADJOINT_FLUID_VECTOR_1_Z, &ADJOINT_FLUID_SCALAR_1};
}

// template instantiations

template class QSVMSAdjointElementData<2, 3>;
template class QSVMSAdjointElementData<2, 4>;

template class QSVMSAdjointElementData<3, 4>;
template class QSVMSAdjointElementData<3, 8>;

} // namespace Kratos
