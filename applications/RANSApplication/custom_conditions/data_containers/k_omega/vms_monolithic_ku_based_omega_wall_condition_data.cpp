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
#include "vms_monolithic_ku_based_omega_wall_condition_data.h"

namespace Kratos
{
namespace KOmegaWallConditionData
{

namespace {
template <unsigned int TDim>
std::vector<const Variable<double>*> GetDofsList();

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
} // namespace

template <unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedOmegaUBasedWallConditionData<TDim, TNumNodes>::Check(
    const Condition& rCondition,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    Fluid::TResidualsDerivatives::Check(rCondition, rProcessInfo);
    TurbulenceModelEquation2::TResidualsDerivatives::Check(rCondition, rProcessInfo);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void VMSMonolithicKBasedOmegaUBasedWallConditionData<TDim, TNumNodes>::InitializeCondition(
    Condition& rCondition,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    Fluid::TResidualsDerivatives::InitializeCondition(rCondition, rProcessInfo);
    TurbulenceModelEquation2::TResidualsDerivatives::InitializeCondition(rCondition, rProcessInfo);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
std::vector<const Variable<double>*> VMSMonolithicKBasedOmegaUBasedWallConditionData<TDim, TNumNodes>::GetDofVariablesList()
{
    return GetDofsList<TDim>();
}

// template instantiations

template class VMSMonolithicKBasedOmegaUBasedWallConditionData<2, 2>;
template class VMSMonolithicKBasedOmegaUBasedWallConditionData<3, 3>;

} // namespace KOmegaWallConditionData

} // namespace Kratos
