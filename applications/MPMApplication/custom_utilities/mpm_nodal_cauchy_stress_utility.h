//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Laura Moreno Martínez

#pragma once

// System includes
#include <limits>
#include <vector>

// Project includes
#include "includes/model_part.h"
#include "includes/variables.h"
#include "utilities/atomic_utilities.h"
#include "utilities/parallel_utilities.h"
#include "mpm_application_variables.h"

namespace Kratos
{

class MPMNodalCauchyStressUtility
{
public:
    static void CalculateNodalCauchyStress(
        ModelPart& rMaterialPointModelPart,
        ModelPart& rGridModelPart)
    {
        const ProcessInfo& r_current_process_info = rMaterialPointModelPart.GetProcessInfo();

        ResetNodalCauchyStress(rGridModelPart, r_current_process_info);
        AddMaterialPointStressContributions(rMaterialPointModelPart, r_current_process_info);
        NormalizeNodalCauchyStress(rGridModelPart);
    }

private:
    static std::size_t GetStressSize(const ProcessInfo& rCurrentProcessInfo)
    {
        const bool is_axisymmetric = (rCurrentProcessInfo.Has(IS_AXISYMMETRIC))
            ? rCurrentProcessInfo.GetValue(IS_AXISYMMETRIC)
            : false;

        return (rCurrentProcessInfo[DOMAIN_SIZE] == 3) ? 6 : (is_axisymmetric ? 4 : 3);
    }

    static void ResetNodalCauchyStress(
        ModelPart& rGridModelPart,
        const ProcessInfo& rCurrentProcessInfo)
    {
        const std::size_t stress_size = GetStressSize(rCurrentProcessInfo);

        block_for_each(rGridModelPart.Nodes(), [&](Node& rNode) {
            auto& r_nodal_cauchy_stress_vector = rNode.FastGetSolutionStepValue(CAUCHY_STRESS_VECTOR);
            if (r_nodal_cauchy_stress_vector.size() != stress_size) {
                r_nodal_cauchy_stress_vector.resize(stress_size, false);
            }
            noalias(r_nodal_cauchy_stress_vector) = ZeroVector(stress_size);
        });
    }

    static void AddMaterialPointStressContributions(
        ModelPart& rMaterialPointModelPart,
        const ProcessInfo& rCurrentProcessInfo)
    {
        auto& r_elements = rMaterialPointModelPart.Elements();

        IndexPartition<std::size_t>(r_elements.size()).for_each([&](std::size_t i_elem) {
            auto it_elem = r_elements.begin() + i_elem;
            auto& r_geometry = it_elem->GetGeometry();
            const std::size_t number_of_nodes = r_geometry.PointsNumber();
            const Matrix& r_N = r_geometry.ShapeFunctionsValues();

            std::vector<Vector> cauchy_stress_values;
            it_elem->CalculateOnIntegrationPoints(MP_CAUCHY_STRESS_VECTOR, cauchy_stress_values, rCurrentProcessInfo);

            std::vector<double> mass_values;
            it_elem->CalculateOnIntegrationPoints(MP_MASS, mass_values, rCurrentProcessInfo);

            for (std::size_t i = 0; i < number_of_nodes; ++i) {

                Vector nodal_cauchy_stress_vector = cauchy_stress_values[0];
                nodal_cauchy_stress_vector *= r_N(0,i) * mass_values[0];

                AtomicAddVector(r_geometry[i].FastGetSolutionStepValue(CAUCHY_STRESS_VECTOR), nodal_cauchy_stress_vector);
            }
        });
    }

    static void NormalizeNodalCauchyStress(
        ModelPart& rGridModelPart)
    {
        block_for_each(rGridModelPart.Nodes(), [&](Node& rNode) {
            const double nodal_mass = rNode.FastGetSolutionStepValue(NODAL_MASS);
            if (nodal_mass > std::numeric_limits<double>::epsilon()) {
                rNode.FastGetSolutionStepValue(CAUCHY_STRESS_VECTOR) /= nodal_mass;
            }
        });
    }
};

} // namespace Kratos
