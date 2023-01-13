//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:     Reza Najian Asl, https://github.com/RezaNajian
//                   Suneth Warnakulasuriya
//


// System includes
#include <tuple>
#include <sstream>

// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/variable_utils.h"
#include "utilities/atomic_utilities.h"
#include "utilities/openmp_utils.h"

// Application includes

// Include base h
#include "linear_strain_energy_response_utilities.h"

namespace Kratos
{

double LinearStrainEnergyResponseUtilities::CalculateStrainEnergy(ModelPart& rModelPart)
{
    KRATOS_TRY

    using tls_type = std::tuple<Matrix, Vector, Vector>;

    const double local_value = block_for_each<SumReduction<double>>(rModelPart.Elements(), tls_type(), [&](auto& rElement, tls_type& rTLS) {
        if (!rElement.IsDefined(ACTIVE) || (rElement.Is(ACTIVE))) {
            Matrix& r_lhs = std::get<0>(rTLS);
            Vector& r_rhs = std::get<1>(rTLS);
            Vector& r_u = std::get<2>(rTLS);

            rElement.GetValuesVector(r_u);
            rElement.CalculateLocalSystem(r_lhs, r_rhs, rModelPart.GetProcessInfo());

            // Compute strain energy
            return 0.5 * inner_prod(r_u, prod(r_lhs, r_u));
        } else {
            return 0.0;
        }
    });

    return rModelPart.GetCommunicator().GetDataCommunicator().SumAll(local_value);

    KRATOS_CATCH("");
}

void LinearStrainEnergyResponseUtilities::CalculateStrainEnergyShapeSensitivity(
    ModelPart& rModelPart,
    const double Delta,
    const Variable<array_1d<double, 3>>& rOutputSensitivityVariable)
{
    KRATOS_TRY

    VariableUtils().SetNonHistoricalVariableToZero(rOutputSensitivityVariable, rModelPart.Nodes());

    using tls_type = std::tuple<Vector, Vector, Vector, Element::Pointer>;

    const auto& r_process_info = rModelPart.GetProcessInfo();

    const int max_id = block_for_each<MaxReduction<int>>(rModelPart.GetRootModelPart().Nodes(), [](const auto& rNode) {
        return rNode.Id();
    });

    std::vector<std::string> model_part_names;

    block_for_each(rModelPart.Elements(), tls_type(), [&](auto& rElement, tls_type& rTLS) {
        if (!rElement.IsDefined(ACTIVE) || (rElement.Is(ACTIVE))) {
            Vector& r_u = std::get<0>(rTLS);
            Vector& r_ref_rhs = std::get<1>(rTLS);
            Vector& r_perturbed_rhs = std::get<2>(rTLS);
            Element::Pointer& p_element = std::get<3>(rTLS);

            auto& r_geometry = rElement.GetGeometry();
            const auto domain_size = r_geometry.WorkingSpaceDimension();

            rElement.GetValuesVector(r_u);

            // calculate the reference value
            rElement.CalculateRightHandSide(r_ref_rhs, r_process_info);

            if (!p_element) {

                std::stringstream name;
                const int thread_id = OpenMPUtils::ThisThread();
                name << rModelPart.Name() << "_temp_" << thread_id;
                GeometryType::PointsArrayType nodes;
                #pragma omp critical
                {
                    model_part_names.push_back(name.str());
                    auto& tls_model_part = rModelPart.CreateSubModelPart(name.str());

                    for (IndexType i = 0; i < r_geometry.size(); ++i) {
                        nodes.push_back(tls_model_part.CreateNewNode((thread_id + 1) * max_id + i + 1, r_geometry[i]));
                    }
                }

                p_element = rElement.Create(1, nodes, rElement.pGetProperties());
                p_element->Initialize(r_process_info);
                p_element->InitializeSolutionStep(r_process_info);
            }

            *p_element = static_cast<const Element&>(rElement);

            // TODO: For some reason, assignment operator of the element
            //       assigns current position to initial position of the lhs element
            //       hence these are corrected manually. But need to be checked
            //       in the element/node/geometry assignment operators.
            p_element->GetGeometry().SetData(r_geometry.GetData());
            p_element->SetFlags(rElement.GetFlags());
            for (IndexType i = 0; i < r_geometry.size(); ++i) {
                auto& r_node = p_element->GetGeometry()[i];
                const auto& r_orig_node = r_geometry[i];
                r_node = r_orig_node;
            }

            // now calculate perturbed
            for (IndexType i = 0; i < r_geometry.size(); ++i) {
                auto& r_orig_node_sensitivity = r_geometry[i].GetValue(rOutputSensitivityVariable);

                auto& r_node = p_element->GetGeometry()[i];
                auto& r_coordinates = r_node.Coordinates();
                auto& r_initial_coordintes = r_node.GetInitialPosition();

                r_initial_coordintes[0] += Delta;
                r_coordinates[0] += Delta;
                p_element->CalculateRightHandSide(r_perturbed_rhs, r_process_info);
                r_initial_coordintes[0] -= Delta;
                r_coordinates[0] -= Delta;
                AtomicAdd<double>(r_orig_node_sensitivity[0], 0.5 * inner_prod(r_u, r_perturbed_rhs - r_ref_rhs) / Delta);

                r_initial_coordintes[1] += Delta;
                r_coordinates[1] += Delta;
                p_element->CalculateRightHandSide(r_perturbed_rhs, r_process_info);
                r_initial_coordintes[1] -= Delta;
                r_coordinates[1] -= Delta;
                AtomicAdd<double>(r_orig_node_sensitivity[1], 0.5 * inner_prod(r_u, r_perturbed_rhs - r_ref_rhs) / Delta);

                if (domain_size == 3) {
                    r_initial_coordintes[2] += Delta;
                    r_coordinates[2] += Delta;
                    p_element->CalculateRightHandSide(r_perturbed_rhs, r_process_info);
                    r_initial_coordintes[2] -= Delta;
                    r_coordinates[2] -= Delta;
                    AtomicAdd<double>(r_orig_node_sensitivity[2], 0.5 * inner_prod(r_u, r_perturbed_rhs - r_ref_rhs) / Delta);
                }
            }
        }
    });

    // now clear the temp model part
    VariableUtils().SetFlag(TO_ERASE, false, rModelPart.Nodes());
    for (const auto& r_model_part_name : model_part_names) {
        auto& tls_model_part = rModelPart.GetSubModelPart(r_model_part_name);
        for (auto& r_node : tls_model_part.Nodes()) {
            r_node.Set(TO_ERASE, true);
        }
    }
    rModelPart.RemoveNodesFromAllLevels(TO_ERASE);

    // Assemble nodal result
    rModelPart.GetCommunicator().AssembleNonHistoricalData(rOutputSensitivityVariable);

    KRATOS_CATCH("");
}

void LinearStrainEnergyResponseUtilities::CalculateStrainEnergyElementPropertiesSensitivity(
    ModelPart& rModelPart,
    const double Delta,
    const Variable<double>& rPrimalVariable,
    const Variable<double>& rOutputSensitivityVariable)
{
    KRATOS_TRY

    VariableUtils().SetNonHistoricalVariableToZero(rOutputSensitivityVariable, rModelPart.Elements());

    using tls_type = std::tuple<Vector, Vector, Vector>;

    const auto& r_process_info = rModelPart.GetProcessInfo();

    block_for_each(rModelPart.Elements(), tls_type(), [&](auto& rElement, tls_type& rTLS) {
        if (!rElement.IsDefined(ACTIVE) || (rElement.Is(ACTIVE))) {
            Vector& r_u = std::get<0>(rTLS);
            Vector& r_ref_rhs = std::get<1>(rTLS);
            Vector& r_perturbed_rhs = std::get<2>(rTLS);

            rElement.GetValuesVector(r_u);

            // calculate the reference value
            rElement.CalculateRightHandSide(r_ref_rhs, r_process_info);

            // now calculate perturbed
            auto& r_properties = rElement.GetProperties();
            r_properties[rPrimalVariable] += Delta;
            rElement.CalculateRightHandSide(r_perturbed_rhs, r_process_info);
            r_properties[rPrimalVariable] -= Delta;

            // now calculate the sensitivity
            rElement.GetValue(rOutputSensitivityVariable) += 0.5 * inner_prod(r_u, r_perturbed_rhs - r_ref_rhs) / Delta;
        }
    });

    KRATOS_CATCH("");
}

}