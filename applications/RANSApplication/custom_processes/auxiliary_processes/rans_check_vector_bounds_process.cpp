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
#include <limits>

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_utilities/rans_check_utilities.h"
#include "utilities/openmp_utils.h"

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

int RansCheckVectorBoundsProcess::Check()
{
    KRATOS_TRY

    const Variable<array_1d<double, 3>>& vector_variable =
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
    const Variable<array_1d<double, 3>>& vector_variable =
        KratosComponents<Variable<array_1d<double, 3>>>::Get(mVariableName);

    const ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);
    const ModelPart::NodesContainerType& r_nodes = r_model_part.Nodes();

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

    double min_value = std::numeric_limits<double>::max();
    double max_value = std::numeric_limits<double>::lowest();

    const int number_of_nodes = r_nodes.size();
    const int number_of_threads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::PartitionVector node_partition;
    OpenMPUtils::DivideInPartitions(number_of_nodes, number_of_threads, node_partition);

    Vector max_values(number_of_threads, max_value);
    Vector min_values(number_of_threads, min_value);

#pragma omp parallel
    {
        const int k = OpenMPUtils::ThisThread();

        auto nodes_begin = r_nodes.begin() + node_partition[k];
        auto nodes_end = r_nodes.begin() + node_partition[k + 1];

        for (auto itNode = nodes_begin; itNode != nodes_end; ++itNode)
        {
            const array_1d<double, 3>& vector_value =
                itNode->FastGetSolutionStepValue(vector_variable);

            double value = 0.0;
            for (int dim = 0; dim < 3; ++dim)
                value += std::pow(vector_value[dim] * vector_weights[dim], 2);
            value = std::pow(value, 0.5);

            min_values[k] = std::min(min_values[k], value);
            max_values[k] = std::max(max_values[k], value);
        }
    }

    for (int i = 0; i < number_of_threads; ++i)
    {
        min_value = std::min(min_value, min_values[i]);
        max_value = std::max(max_value, max_values[i]);
    }

    const Communicator& r_communicator = r_model_part.GetCommunicator();
    min_value = r_communicator.GetDataCommunicator().MinAll(min_value);
    max_value = r_communicator.GetDataCommunicator().MaxAll(max_value);

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
