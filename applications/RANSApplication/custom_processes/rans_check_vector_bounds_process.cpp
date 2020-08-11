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
#include <limits>

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "utilities/parallel_utilities.h"

// Include base h
#include "rans_check_vector_bounds_process.h"

namespace Kratos
{
RansCheckVectorBoundsProcess::RansCheckVectorBoundsProcess(
    Model& rModel,
    Parameters rParameters)
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

    if (mrParameters["component_type"].GetString() == "magnitude") {
        mVectorComponent = VectorComponent::Magnitude;
    } else if (mrParameters["component_type"].GetString() == "x") {
        mVectorComponent = VectorComponent::X;
    } else if (mrParameters["component_type"].GetString() == "y") {
        mVectorComponent = VectorComponent::Y;
    } else if (mrParameters["component_type"].GetString() == "z") {
        mVectorComponent = VectorComponent::Z;
    } else {
        KRATOS_ERROR
            << "Vector component_type type not found. [ component_type = "
            << mrParameters["component_type"].GetString() << " ]\n.";
    }

    KRATOS_CATCH("");
}

int RansCheckVectorBoundsProcess::Check()
{
    KRATOS_TRY

    const Variable<array_1d<double, 3>>& vector_variable =
        KratosComponents<Variable<array_1d<double, 3>>>::Get(mVariableName);

    const auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(vector_variable))
        << mVariableName << " is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    return 0;

    KRATOS_CATCH("");
}

void RansCheckVectorBoundsProcess::Execute()
{
    KRATOS_TRY
    const Variable<array_1d<double, 3>>& r_vector_variable =
        KratosComponents<Variable<array_1d<double, 3>>>::Get(mVariableName);

    const auto& r_model_part = mrModel.GetModelPart(mModelPartName);
    const auto& r_nodes = r_model_part.Nodes();

    array_1d<double, 3> vector_weights;

    switch (mVectorComponent) {
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

    class MinMaxReducer{
        public:
            typedef std::tuple<double,double> value_type;
            double min_value = std::numeric_limits<double>::max();
            double max_value = std::numeric_limits<double>::lowest();

            value_type GetValue()
            {
                value_type values;
                std::get<0>(values) = min_value;
                std::get<1>(values) = max_value;
                return values;
            }

            void LocalReduce(const double Value){
                min_value = std::min(min_value, Value);
                max_value = std::max(max_value, Value);
            }
            void ThreadSafeReduce(MinMaxReducer& rOther){
                #pragma omp critical
                {
                    min_value = std::min(min_value, rOther.min_value);
                    max_value = std::max(max_value, rOther.max_value);
                }
            }
    };

    const int number_of_nodes = r_nodes.size();

    double min_value, max_value;
    std::tie(min_value, max_value) =
        IndexPartition<int>(number_of_nodes).for_each<MinMaxReducer>([&](const int i_node) -> double {
            const array_1d<double, 3>& r_vector_value =
                (r_nodes.begin() + i_node)->FastGetSolutionStepValue(r_vector_variable);

            double value = 0.0;
            for (int dim = 0; dim < 3; ++dim) {
                value += std::pow(r_vector_value[dim] * vector_weights[dim], 2);
            }
            value = std::pow(value, 0.5);

            return value;
        });

    const auto& r_communicator = r_model_part.GetCommunicator();
    min_value = r_communicator.GetDataCommunicator().MinAll(min_value);
    max_value = r_communicator.GetDataCommunicator().MaxAll(max_value);

    KRATOS_INFO(this->Info())
        << r_vector_variable.Name() << " is bounded between [ " << min_value
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
