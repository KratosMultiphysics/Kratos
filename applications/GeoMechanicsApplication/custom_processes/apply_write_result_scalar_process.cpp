// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#include "apply_write_result_scalar_process.h"
#include "includes/model_part.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

ApplyWriteScalarProcess::ApplyWriteScalarProcess(ModelPart& model_part, Parameters rParameters)
    : Process(Flags()), mrModelPart(model_part)
{
    KRATOS_TRY

    // only include validation with c++11 since raw_literals do not exist in c++03
    Parameters default_parameters(R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "append_file" : false
            }  )");

    // Some values need to be mandatorily prescribed since no meaningful default value exist.
    // For this reason try accessing to them So that an error is thrown if they don't exist
    rParameters["variable_name"];
    rParameters["model_part_name"];
    rParameters["append_file"];

    // Now validate agains defaults -- this also ensures no type mismatch
    rParameters.ValidateAndAssignDefaults(default_parameters);

    mModelPartName     = rParameters["model_part_name"].GetString();
    mVariableName      = rParameters["variable_name"].GetString();
    mAppendFile        = rParameters["append_file"].GetBool();
    mTimeUnitConverter = model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];

    KRATOS_CATCH("")
}

void ApplyWriteScalarProcess::Execute() {}

/// this function is designed for being called at the beginning of the computations
/// right after reading the model and the groups
void ApplyWriteScalarProcess::ExecuteInitialize()
{
    KRATOS_TRY

    const SizeType nNodes = mrModelPart.NumberOfNodes();

    if (nNodes > 0) {
        const Variable<double>& var  = KratosComponents<Variable<double>>::Get(mVariableName);
        const double            Time = mrModelPart.GetProcessInfo()[TIME] / mTimeUnitConverter;

        ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();
        mOutFile.resize(nNodes);

        for (IndexType i = 0; i < nNodes; ++i) {
            ModelPart::NodesContainerType::iterator it = it_begin + i;

            const IndexType nodeId = it->Id();
            std::string fileName = mModelPartName + "_" + std::to_string(nodeId) + "_" + mVariableName + ".res";

            if (mAppendFile) {
                // append instead of overwrite
                mOutFile[i].open(fileName, std::ios::app);
            } else {
                // open a new file and overwrite
                mOutFile[i].open(fileName, std::ios::trunc); // overwrite
                mOutFile[i] << "Time"
                            << "   " << mVariableName << "\n";
                const double value = it->FastGetSolutionStepValue(var);
                mOutFile[i] << Time << "   " << value << "\n";
            }
        }
    }

    KRATOS_CATCH("")
}

void ApplyWriteScalarProcess::ExecuteFinalizeSolutionStep()
{
    KRATOS_TRY

    const IndexType nNodes = mrModelPart.NumberOfNodes();

    if (nNodes > 0) {
        const Variable<double>& var  = KratosComponents<Variable<double>>::Get(mVariableName);
        const double            Time = mrModelPart.GetProcessInfo()[TIME] / mTimeUnitConverter;
        ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

        for (IndexType i = 0; i < nNodes; ++i) {
            ModelPart::NodesContainerType::iterator it = it_begin + i;

            const double value = it->FastGetSolutionStepValue(var);
            mOutFile[i] << Time << "   " << value << "\n";
        }
    }

    KRATOS_CATCH("")
}

void ApplyWriteScalarProcess::ExecuteFinalize()
{
    KRATOS_TRY

    for (auto& outFile : mOutFile) {
        outFile.close();
    }

    KRATOS_CATCH("")
}

std::string ApplyWriteScalarProcess::Info() const { return "ApplyWriteScalarProcess"; }

} // namespace Kratos