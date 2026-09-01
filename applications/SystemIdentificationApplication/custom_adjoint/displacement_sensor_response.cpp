//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
//
//  Main authors:    Máté Kelemen
//

// --- Kratos Core Includes ---
#include "custom_adjoint/displacement_sensor_response.hpp"
#include "utilities/brute_force_point_locator.h"


namespace Kratos {


SensorResponse::Pointer DisplacementSensorResponse::Create(
    const ModelPart& rDomainModelPart,
    const ModelPart& rSensorModelPart,
    IndexType Id,
    Parameters SensorParameters) const {
        KRATOS_TRY
            DisplacementSensorResponse prototype;
            SensorParameters.ValidateAndAssignDefaults(prototype.GetDefaultParameters());

            // Parse sensor name.
            const std::string name = SensorParameters["name"].GetString();
            KRATOS_ERROR_IF(name.empty()) << "displacement sensor response must have a name";

            // Parse sensor direction.
            const auto& direction = SensorParameters["direction"].GetVector();
            KRATOS_ERROR_IF_NOT(direction.size() == 3) << std::format(
                "expecting \"{}\" displacement sensor to have a \"direction\" of exactly 3 components, but got {}",
                name,
                direction.size());
            const double direction_norm = norm_2(direction);
            KRATOS_ERROR_IF_NOT(direction_norm) << std::format(
                "degenerate \"direction\" for displacement sensor {}",
                name);

            // Parse sensor position.
            const auto& position = SensorParameters["position"].GetVector();
            KRATOS_ERROR_IF_NOT(position.size() == 3) << std::format(
                "expecting \"{}\" displacement sensor to have a \"position\" of exactly 3 components, but got {}",
                name,
                position.size());

            Point pos(position[0], position[1], position[2]);
            Vector dummy_shape_functions;

            // Parse design variables.
            const auto design_variable_names = SensorParameters["design_variables"].GetStringArray();
            std::vector<const Variable<double>*> design_variables;
            design_variables.reserve(design_variable_names.size());
            std::transform(
                design_variable_names.begin(),
                design_variable_names.end(),
                std::back_inserter(design_variables),
                [] (const std::string& r_design_variable_name) -> const Variable<double>* {
                    const bool has_variable = KratosComponents<Variable<double>>::Has(r_design_variable_name);
                    KRATOS_ERROR_IF_NOT(has_variable) << std::format("design variable \"{}\" is not registered", r_design_variable_name);
                    return &KratosComponents<Variable<double>>::Get(r_design_variable_name);
                });

            // Find which element this sensor is located in.
            const auto element_id = BruteForcePointLocator(
                const_cast<ModelPart&>(rDomainModelPart)
            ).FindElement(pos, dummy_shape_functions);
            const auto& rp_element = rDomainModelPart.pGetElement(element_id);

            // Construct a node for this sensor.
            ModelPart& r_mutable_sensor_model_part = const_cast<ModelPart&>(rSensorModelPart);
            auto p_node = r_mutable_sensor_model_part.CreateNewNode(Id, position[0], position[1], position[2]);

            // Construct the sensor.
            const auto p_sensor = Kratos::make_shared<DisplacementSensorResponse>(
                design_variables,
                name,
                p_node,
                rp_element);

            for (std::size_t i_component=0ul; i_component<mDirection.size(); ++i_component)
                p_sensor->mDirection[i_component] = direction[i_component] / direction_norm;

            return p_sensor;
        KRATOS_CATCH("")
}


double DisplacementSensorResponse::ComputeValue(
    const Element& rElement,
    const ProcessInfo& rProcessInfo,
    int iBuffer) const {
        KRATOS_ERROR_IF_NOT(iBuffer == 0)
            << "requested buffer " << iBuffer << ", which is not supported";

        KRATOS_TRY
            if (rElement.Id() == this->GetContainingElement().Id()) {
                double directional_displacement = 0.0;
                const auto& r_geometry = this->GetContainingElement().GetGeometry();

                // Map the sensor's location from physical to the element's parametric space.
                const array_1d<double,3>& r_physical_coordinates = *this->GetNode();
                array_1d<double,3> parametric_coordinates {0.0, 0.0, 0.0};
                r_geometry.PointLocalCoordinates(
                    parametric_coordinates,
                    r_physical_coordinates);

                // Compute shape function values at the sensor's parametric coordinates.
                Vector shape_function_values;
                r_geometry.ShapeFunctionsValues(
                    shape_function_values,
                    parametric_coordinates);

                array_1d<double, 3> displacement = ZeroVector(3);
                for (IndexType i = 0; i < r_geometry.size(); ++i) {
                    displacement += r_geometry[i].GetSolutionStepValue(DISPLACEMENT) * shape_function_values[i];
                }

                directional_displacement = inner_prod(displacement, mDirection);
                return directional_displacement;
            } else {
                return 0.0;
            }
        KRATOS_CATCH("")
}


void DisplacementSensorResponse::GetStateVariables(
    std::vector<IAdjoint::DynamicVariable>& rVariables,
    const Element& rElement,
    const ProcessInfo& rProcessInfo) const {
        KRATOS_TRY
            rElement.GetInfluencingVariables(
                rVariables,
                rProcessInfo);
            rVariables.erase(
                std::remove_if(
                    rVariables.begin(),
                    rVariables.end(),
                    [] (const IAdjoint::DynamicVariable& r_variable) -> bool {
                        return r_variable.SourceKey() != DISPLACEMENT.SourceKey();
                    }),
                rVariables.end());
        KRATOS_CATCH("")
}


void DisplacementSensorResponse::ComputeDerivative(
    Vector& rOutput,
    const Element& rElement,
    std::span<const IAdjoint::DynamicVariable> Variables,
    const ProcessInfo& rProcessInfo,
    int iBuffer) const {
        KRATOS_ERROR_IF_NOT(iBuffer == 0)
            << "requested buffer " << iBuffer << ", which is not supported";

        KRATOS_TRY
            rOutput = ZeroVector(Variables.size());
            if (rElement.Id() == this->GetContainingElement().Id()) {
                const Geometry<Node>& r_geometry = rElement.GetGeometry();

                // Map the sensor's location from physical to the element's parametric space.
                const array_1d<double,3>& r_physical_coordinates = *this->GetNode();
                array_1d<double,3> parametric_coordinates {0.0, 0.0, 0.0};
                rElement.GetGeometry().PointLocalCoordinates(
                    parametric_coordinates,
                    r_physical_coordinates);

                // Compute shape function values at the sensor's parametric coordinates.
                Vector shape_function_values;
                rElement.GetGeometry().ShapeFunctionsValues(
                    shape_function_values,
                    parametric_coordinates);

                // Build the output.
                for (std::size_t i_variable=0ul; i_variable<Variables.size(); ++i_variable) {
                    const auto& r_variable = Variables[i_variable];
                    const auto i_node = r_variable.GetDynamicIndex();
                    if (i_node != -1 && r_variable.SourceKey() == DISPLACEMENT.SourceKey()) {
                        KRATOS_ERROR_IF_NOT(i_node < r_geometry.size())
                            << "variable " << r_variable.Name() << " "
                            << "with index " << i_node << " "
                            << "in geometry of " << r_geometry.size() << " nodes";
                        rOutput[i_variable] = shape_function_values[i_node] * mDirection[r_variable.GetComponentIndex()];
                    } /*if i_node != -1 && DISPLACEMENT_* */ else {
                        // Check whether the requested variable is a design variable.
                        const auto& r_design_variable_types = this->GetDesignVariableTypes();
                        const auto it_design_variable_type = std::find_if(
                            r_design_variable_types.begin(),
                            r_design_variable_types.end(),
                            [&r_variable] (const Variable<double>* p_variable_type) -> bool {
                                return p_variable_type->Key() == r_variable;
                            });
                        const bool is_design_variable = it_design_variable_type != r_design_variable_types.end();

                        // Error if the requested variable is neither a state nor a design variable.
                        KRATOS_ERROR_IF_NOT(is_design_variable) << "unsupported variable " << r_variable.Name();
                    }
                } // for variable
            }
        KRATOS_CATCH("")
}


Parameters DisplacementSensorResponse::GetDefaultParameters() const {
    return Parameters(R"({
        "position" : [0, 0, 0],
        "direction" : [0, 0, 0],
        "name" : "",
        "design_variables" : []
    })");
}


} // namespace Kratos
