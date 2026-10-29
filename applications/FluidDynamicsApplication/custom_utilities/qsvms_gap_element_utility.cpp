//     |  /           |
//     ' /   __| _` | __|  _ \   __|
//     . \  |   (   | |   (   |\__ `
//    _|\_\_|  \__,_|\__|\___/ ____/
//                    Multi-Physics
//
//   Fluid Dynamics Application
//
//   License:         BSD License
//                    Kratos default license: kratos/license.txt
//
//   Main author:     Nicolò Antonelli
//

#include "custom_utilities/qsvms_gap_element_utility.h"
#include "custom_elements/qs_vms.h"
#include "custom_elements/data_containers/time_integrated_qs_vms/time_integrated_qs_vms_data.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
Element::Pointer QSVMSGapElementUtility::Create(
    Element::IndexType NewId,
    Element& rParent,
    Geometry<Node>::Pointer pQuadratureGeometry)
{
    KRATOS_TRY
    using QSVMSQ1 = QSVMS<TimeIntegratedQSVMSData<2, 4>>;
    KRATOS_ERROR_IF(dynamic_cast<QSVMSQ1*>(&rParent) == nullptr)
        << "QSVMS gap integration requires TimeIntegratedQSVMS2D4N." << std::endl;
    KRATOS_ERROR_IF_NOT(pQuadratureGeometry) << "Missing gap quadrature geometry." << std::endl;
    const auto& r_parent_geometry = rParent.GetGeometry();
    KRATOS_ERROR_IF(r_parent_geometry.PointsNumber() != 4 || pQuadratureGeometry->PointsNumber() != 4)
        << "QSVMS gap integration requires four parent Q1 nodes." << std::endl;
    for (std::size_t i = 0; i < 4; ++i) {
        KRATOS_ERROR_IF(pQuadratureGeometry->pGetPoint(i) != r_parent_geometry.pGetPoint(i))
            << "Gap quadrature must use the parent nodes in parent order." << std::endl;
    }
    auto p_element = rParent.Create(NewId, pQuadratureGeometry, rParent.pGetProperties());
    // QSVMS reads this coefficient from the element rather than Properties.
    if (rParent.Has(C_SMAGORINSKY)) {
        p_element->SetValue(C_SMAGORINSKY, rParent.GetValue(C_SMAGORINSKY));
    }
    p_element->Set(ACTIVE, true);
    return p_element;
    KRATOS_CATCH("")
}
} // namespace Kratos
