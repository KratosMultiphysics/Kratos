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

// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/variable_utils.h"

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

            auto& r_geometry = rElement.GetGeometry();
            const auto domain_size = r_geometry.WorkingSpaceDimension();

            // now calculate perturbed
            for (auto& r_node : r_geometry) {
                r_node.SetLock();

                auto& r_output_value = r_node.GetValue(rOutputSensitivityVariable);
                auto& r_coordinates = r_node.Coordinates();
                auto& r_initial_coordintes = r_node.GetInitialPosition();

                r_initial_coordintes[0] += Delta;
                r_coordinates[0] += Delta;
                rElement.CalculateRightHandSide(r_perturbed_rhs, r_process_info);
                r_initial_coordintes[0] -= Delta;
                r_coordinates[0] -= Delta;
                r_output_value[0] += 0.5 * inner_prod(r_u, r_perturbed_rhs - r_ref_rhs) / Delta;

                r_initial_coordintes[1] += Delta;
                r_coordinates[1] += Delta;
                rElement.CalculateRightHandSide(r_perturbed_rhs, r_process_info);
                r_initial_coordintes[1] -= Delta;
                r_coordinates[1] -= Delta;
                r_output_value[1] += 0.5 * inner_prod(r_u, r_perturbed_rhs - r_ref_rhs) / Delta;

                if (domain_size == 3) {
                    r_initial_coordintes[2] += Delta;
                    r_coordinates[2] += Delta;
                    rElement.CalculateRightHandSide(r_perturbed_rhs, r_process_info);
                    r_initial_coordintes[2] -= Delta;
                    r_coordinates[2] -= Delta;
                    r_output_value[2] += 0.5 * inner_prod(r_u, r_perturbed_rhs - r_ref_rhs) / Delta;
                }

                r_node.UnSetLock();
            }
        }
    });

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