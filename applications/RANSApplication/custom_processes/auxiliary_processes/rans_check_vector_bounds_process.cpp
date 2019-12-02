//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_utilities/rans_check_utilities.h"

// Include base h
#include "rans_check_vector_bounds_process.h"

namespace Kratos
{
RansCheckVectorBoundsProcess::RansCheckVectorBoundsProcess(Model& rModel, Parameters rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "variable_name"   : "PLEASE_SPECIFY_SCALAR_VARIABLE",
            "component_type"  : "magnitude"
        })");

    mrParameters.ValidateAndAssignDefaults(default_parameters);

    mVariableName = mrParameters["variable_name"].GetString();
    mModelPartName = mrParameters["model_part_name"].GetString();

    if (mrParameters["component_type"].GetString() == "magnitude")
        mVectorComponent = VectorComponent::Magnitude;
    else if (mrParameters["component_type"].GetString() == "x")
        mVectorComponent = VectorComponent::X;
    else if (mrParameters["component_type"].GetString() == "y")
        mVectorComponent = VectorComponent::Y;
    else if (mrParameters["component_type"].GetString() == "z")
        mVectorComponent = VectorComponent::Z;
    else
        KRATOS_ERROR
            << "Vector component_type type not found. [ component_type = "
            << mrParameters["component_type"].GetString() << " ]\n.";

    KRATOS_CATCH("");
}

RansCheckVectorBoundsProcess::~RansCheckVectorBoundsProcess()
{
}

int RansCheckVectorBoundsProcess::Check()
{
    KRATOS_TRY

    const Variable<array_1d<double, 3>> vector_variable =
        KratosComponents<Variable<array_1d<double, 3>>>::Get(mVariableName);

        RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(
        mrModel.GetModelPart(mModelPartName), vector_variable);

    return 0;

    KRATOS_CATCH("");
}

void RansCheckVectorBoundsProcess::Execute()
{
    KRATOS_TRY
    const Variable<array_1d<double, 3>> vector_variable =
        KratosComponents<Variable<array_1d<double, 3>>>::Get(mVariableName);

    const ModelPart::NodesContainerType& r_nodes =
        mrModel.GetModelPart(mModelPartName).Nodes();

    array_1d<double, 3> vector_weights;

    switch (mVectorComponent)
    {
    case VectorComponent::Magnitude:
        vector_weights[0] = 1.0;
        vector_weights[1] = 1.0;
        vector_weights[2] = 1.0;
        break;
    case VectorComponent::X:
        vector_weights[0] = 1.0;
        vector_weights[1] = 0.0;
        vector_weights[2] = 0.0;
        break;
    case VectorComponent::Y:
        vector_weights[0] = 0.0;
        vector_weights[1] = 1.0;
        vector_weights[2] = 0.0;
        break;
    case VectorComponent::Z:
        vector_weights[0] = 0.0;
        vector_weights[1] = 0.0;
        vector_weights[2] = 1.0;
        break;
    }

    double min_value = 0.0;
    double max_value = 0.0;
    bool is_initialized = false;

    const int number_of_nodes = r_nodes.size();
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        const ModelPart::NodeType& r_node = *(r_nodes.begin() + i_node);
        const array_1d<double, 3>& vector_value =
            r_node.FastGetSolutionStepValue(vector_variable);

        double current_value = 0.0;
        for (int dim = 0; dim < 3; ++dim)
            current_value += std::pow(vector_value[dim] * vector_weights[dim], 2);

        current_value = std::pow(current_value, 0.5);

        if (!is_initialized)
        {
            min_value = current_value;
            max_value = current_value;
            is_initialized = true;
        }

        if (min_value > current_value)
            min_value = current_value;
        if (max_value < current_value)
            max_value = current_value;
    }

    KRATOS_INFO(this->Info())
        << vector_variable.Name() << " is bounded between [ " << min_value
        << ", " << max_value << " ] in " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

std::string RansCheckVectorBoundsProcess::Info() const
{
    return std::string("RansCheckVectorBoundsProcess");
}

void RansCheckVectorBoundsProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansCheckVectorBoundsProcess::PrintData(std::ostream& rOStream) const
{
}

} // namespace Kratos.
