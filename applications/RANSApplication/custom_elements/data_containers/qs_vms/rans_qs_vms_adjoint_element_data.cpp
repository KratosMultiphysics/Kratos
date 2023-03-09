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
#include "custom_elements/data_containers/k_epsilon/k_element_data.h"
#include "custom_elements/data_containers/k_epsilon/epsilon_element_data.h"

// Include base h
#include "rans_qs_vms_adjoint_element_data.h"

namespace Kratos
{
namespace {
template <unsigned int TDim>
std::vector<const Variable<double>*> GetQSVMSAdjointDofsList();

template <>
std::vector<const Variable<double>*> GetQSVMSAdjointDofsList<2>()
{
    return {&ADJOINT_FLUID_VECTOR_1_X, &ADJOINT_FLUID_VECTOR_1_Y, &ADJOINT_FLUID_SCALAR_1};
}

template <>
std::vector<const Variable<double>*> GetQSVMSAdjointDofsList<3>()
{
    return {&ADJOINT_FLUID_VECTOR_1_X, &ADJOINT_FLUID_VECTOR_1_Y,
            &ADJOINT_FLUID_VECTOR_1_Z, &ADJOINT_FLUID_SCALAR_1};
}
} // namespace

template <unsigned int TDim, unsigned int TNumNodes, class TTurbulenceModelElementData1, class TTurbulenceModelElementData2>
void RansQSVMSAdjointElementData<TDim, TNumNodes, TTurbulenceModelElementData1, TTurbulenceModelElementData2>::Check(
    const Element& rElement,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    TResidualsDerivatives::Check(rElement, rProcessInfo);

    const auto& dofs_list = GetDofVariablesList();

    const auto& r_geometry = rElement.GetGeometry();
    for (IndexType c = 0; c < TNumNodes; ++c) {
        const auto& r_node = r_geometry[c];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_FLUID_VECTOR_1, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_FLUID_VECTOR_2, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_FLUID_VECTOR_3, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_FLUID_SCALAR_1, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(AUX_ADJOINT_FLUID_VECTOR_1, r_node);

        for (const auto p_variable : dofs_list) {
            KRATOS_CHECK_DOF_IN_NODE((*p_variable), r_node);
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TTurbulenceModelElementData1, class TTurbulenceModelElementData2>
std::vector<const Variable<double>*> RansQSVMSAdjointElementData<TDim, TNumNodes, TTurbulenceModelElementData1, TTurbulenceModelElementData2>::GetDofVariablesList()
{
    return GetQSVMSAdjointDofsList<TDim>();
}

template <unsigned int TDim, unsigned int TNumNodes, class TTurbulenceModelElementData1, class TTurbulenceModelElementData2>
GeometryData::IntegrationMethod RansQSVMSAdjointElementData<TDim, TNumNodes, TTurbulenceModelElementData1, TTurbulenceModelElementData2>::GetIntegrationMethod()
{
    return TResidualsDerivatives::GetIntegrationMethod();
}

// template instantiations

template class RansQSVMSAdjointElementData<2, 3, KEpsilonElementData::KElementData<2>, KEpsilonElementData::EpsilonElementData<2>>;
template class RansQSVMSAdjointElementData<2, 4, KEpsilonElementData::KElementData<2>, KEpsilonElementData::EpsilonElementData<2>>;

template class RansQSVMSAdjointElementData<3, 4, KEpsilonElementData::KElementData<3>, KEpsilonElementData::EpsilonElementData<3>>;
template class RansQSVMSAdjointElementData<3, 8, KEpsilonElementData::KElementData<3>, KEpsilonElementData::EpsilonElementData<3>>;

} // namespace Kratos
