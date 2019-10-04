//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:             BSD License
//                               Kratos default license:
//kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#include "custom_processes/update_pressure_value_pfem_conditions_process.h"

namespace Kratos {

template <SizeType TDim>
UpdatePressureValuePfemConditionsProcess<TDim>::UpdatePressureValuePfemConditionsProcess(
    ModelPart& rModelPart)
    : mrModelPart(rModelPart)
{
}

/***********************************************************************************/
/***********************************************************************************/
template <SizeType TDim>
void UpdatePressureValuePfemConditionsProcess<TDim>::Execute()
{
    auto& pfem_cond_sub_model = mrModelPart.GetSubModelPart("PFEMPressureConditions");
    const auto it_cond_begin = pfem_cond_sub_model.ConditionsBegin();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(pfem_cond_sub_model.Conditions().size()); i++) {
        auto it_cond = it_cond_begin + i;
        auto& r_geometry = it_cond->GetGeometry();
        double average_pressure = 0.0;

        for (unsigned int i = 0; i < r_geometry.PointsNumber(); i++) {
            const auto &r_node = r_geometry[i];
            const double nodal_pressure = r_node.FastGetSolutionStepValue(PRESSURE);
            average_pressure += nodal_pressure;
        }
        average_pressure /= r_geometry.PointsNumber();
        it_cond->SetValue(POSITIVE_FACE_PRESSURE, average_pressure);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class UpdatePressureValuePfemConditionsProcess<2>;
template class UpdatePressureValuePfemConditionsProcess<3>;

}  // namespace Kratos