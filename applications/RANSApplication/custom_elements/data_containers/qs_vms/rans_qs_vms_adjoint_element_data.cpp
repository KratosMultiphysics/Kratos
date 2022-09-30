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

// Element data containers
#include "custom_elements/data_containers/k_epsilon/k_element_data.h"
#include "custom_elements/data_containers/k_epsilon/epsilon_element_data.h"
#include "custom_elements/data_containers/k_omega/k_element_data.h"
#include "custom_elements/data_containers/k_omega/omega_element_data.h"
#include "custom_elements/data_containers/k_omega_sst/k_element_data.h"
#include "custom_elements/data_containers/k_omega_sst/omega_element_data.h"
#include "custom_elements/data_containers/qs_vms/qs_vms_derivative_utilities.h"
#include "custom_elements/data_containers/qs_vms/rans_qs_vms_derivative_utilities.h"

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

template <unsigned int TDim, unsigned int TNumNodes, template<unsigned int, unsigned int> class TVelocityDerivative, template<unsigned int, unsigned int> class TShapeDerivative, class TTurbulenceModelElementData1, class TTurbulenceModelElementData2>
void RansQSVMSAdjointElementData<TDim, TNumNodes, TVelocityDerivative, TShapeDerivative, TTurbulenceModelElementData1, TTurbulenceModelElementData2>::Check(
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

template <unsigned int TDim, unsigned int TNumNodes, template<unsigned int, unsigned int> class TVelocityDerivative, template<unsigned int, unsigned int> class TShapeDerivative, class TTurbulenceModelElementData1, class TTurbulenceModelElementData2>
std::vector<const Variable<double>*> RansQSVMSAdjointElementData<TDim, TNumNodes, TVelocityDerivative, TShapeDerivative, TTurbulenceModelElementData1, TTurbulenceModelElementData2>::GetDofVariablesList()
{
    return GetQSVMSAdjointDofsList<TDim>();
}

template <unsigned int TDim, unsigned int TNumNodes, template<unsigned int, unsigned int> class TVelocityDerivative, template<unsigned int, unsigned int> class TShapeDerivative, class TTurbulenceModelElementData1, class TTurbulenceModelElementData2>
GeometryData::IntegrationMethod RansQSVMSAdjointElementData<TDim, TNumNodes, TVelocityDerivative, TShapeDerivative, TTurbulenceModelElementData1, TTurbulenceModelElementData2>::GetIntegrationMethod()
{
    return TResidualsDerivatives::GetIntegrationMethod();
}

// template instantiations
// k-epsilon
template class RansQSVMSAdjointElementData<2, 3, QSVMSDerivativeUtilities<2>::VelocityDerivative, QSVMSDerivativeUtilities<2>::ShapeDerivative, KEpsilonElementData::KElementData<2>, KEpsilonElementData::EpsilonElementData<2>>;
template class RansQSVMSAdjointElementData<2, 4, QSVMSDerivativeUtilities<2>::VelocityDerivative, QSVMSDerivativeUtilities<2>::ShapeDerivative, KEpsilonElementData::KElementData<2>, KEpsilonElementData::EpsilonElementData<2>>;

template class RansQSVMSAdjointElementData<3, 4, QSVMSDerivativeUtilities<3>::VelocityDerivative, QSVMSDerivativeUtilities<3>::ShapeDerivative, KEpsilonElementData::KElementData<3>, KEpsilonElementData::EpsilonElementData<3>>;
template class RansQSVMSAdjointElementData<3, 8, QSVMSDerivativeUtilities<3>::VelocityDerivative, QSVMSDerivativeUtilities<3>::ShapeDerivative, KEpsilonElementData::KElementData<3>, KEpsilonElementData::EpsilonElementData<3>>;

// k-omega
template class RansQSVMSAdjointElementData<2, 3, QSVMSDerivativeUtilities<2>::VelocityDerivative, QSVMSDerivativeUtilities<2>::ShapeDerivative, KOmegaElementData::KElementData<2>, KOmegaElementData::OmegaElementData<2>>;
template class RansQSVMSAdjointElementData<2, 4, QSVMSDerivativeUtilities<2>::VelocityDerivative, QSVMSDerivativeUtilities<2>::ShapeDerivative, KOmegaElementData::KElementData<2>, KOmegaElementData::OmegaElementData<2>>;

template class RansQSVMSAdjointElementData<3, 4, QSVMSDerivativeUtilities<3>::VelocityDerivative, QSVMSDerivativeUtilities<3>::ShapeDerivative, KOmegaElementData::KElementData<3>, KOmegaElementData::OmegaElementData<3>>;
template class RansQSVMSAdjointElementData<3, 8, QSVMSDerivativeUtilities<3>::VelocityDerivative, QSVMSDerivativeUtilities<3>::ShapeDerivative, KOmegaElementData::KElementData<3>, KOmegaElementData::OmegaElementData<3>>;

// k-omega-sst
template class RansQSVMSAdjointElementData<2, 3, RansQSVMSDerivativeUtilities<2>::KOmegaSSTVelocityDerivative, RansQSVMSDerivativeUtilities<2>::KOmegaSSTShapeDerivative, KOmegaSSTElementData::KElementData<2>, KOmegaSSTElementData::OmegaElementData<2>>;
template class RansQSVMSAdjointElementData<2, 4, RansQSVMSDerivativeUtilities<2>::KOmegaSSTVelocityDerivative, RansQSVMSDerivativeUtilities<2>::KOmegaSSTShapeDerivative, KOmegaSSTElementData::KElementData<2>, KOmegaSSTElementData::OmegaElementData<2>>;

template class RansQSVMSAdjointElementData<3, 4, RansQSVMSDerivativeUtilities<3>::KOmegaSSTVelocityDerivative, RansQSVMSDerivativeUtilities<3>::KOmegaSSTShapeDerivative, KOmegaSSTElementData::KElementData<3>, KOmegaSSTElementData::OmegaElementData<3>>;
template class RansQSVMSAdjointElementData<3, 8, RansQSVMSDerivativeUtilities<3>::KOmegaSSTVelocityDerivative, RansQSVMSDerivativeUtilities<3>::KOmegaSSTShapeDerivative, KOmegaSSTElementData::KElementData<3>, KOmegaSSTElementData::OmegaElementData<3>>;


} // namespace Kratos
