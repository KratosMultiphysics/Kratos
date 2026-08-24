//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

// --- Kratos Core Includes ---
#include "custom_adjoint/displacement_sensor_response.hpp"


namespace Kratos {


double DisplacementSensorResponse::ComputeValue(
    const ModelPart& rModelPart,
    int iBuffer) const {
        KRATOS_ERROR_IF_NOT(iBuffer == 0)
            << "requested buffer " << iBuffer << ", which is not supported";

        KRATOS_TRY
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
            return rModelPart.GetCommunicator().GetDataCommunicator().SumAll(directional_displacement);
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


void DisplacementSensorResponse::GetStateVariables(
    std::vector<IAdjoint::DynamicVariable>& rVariables,
    const Condition& rCondition,
    const ProcessInfo& rProcessInfo) const {
        rVariables.clear();
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
                    KRATOS_ERROR_IF_NOT(i_node < r_geometry.size())
                        << "variable " << r_variable.Name() << " "
                        << "with index " << i_node << " "
                        << "in geometry of " << r_geometry.size() << " nodes";

                    if (r_variable.SourceKey() == DISPLACEMENT.SourceKey()) {
                        // Requested variable is a state variable.
                        rOutput[i_variable] = shape_function_values[i_node] * mDirection[r_variable.GetComponentIndex()];
                    } else {
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


void DisplacementSensorResponse::ComputeDerivative(
    Vector& rOutput,
    const Condition& rCondition,
    std::span<const IAdjoint::DynamicVariable> Variables,
    const ProcessInfo& rProcessInfo,
    int iBuffer) const {
        KRATOS_ERROR_IF_NOT(iBuffer == 0)
            << "requested buffer " << iBuffer << ", which is not supported";
        for (const auto& r_variable : Variables)
            KRATOS_ERROR_IF_NOT(r_variable.SourceKey() == DISPLACEMENT.SourceKey())
                << "unsupported variable " << r_variable.Name();
        rOutput = ZeroVector(Variables.size());
}


} // namespace Kratos
