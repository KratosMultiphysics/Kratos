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

#pragma once

#include "includes/table.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyWriteScalarProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyWriteScalarProcess);

    /// Defining a table with double argument and result type as table type.
    using TableType = Table<double,double>;
    using IndexType = std::size_t;

    ApplyWriteScalarProcess(ModelPart& model_part,
                            Parameters rParameters
                            ) : Process(Flags()) , mrModelPart(model_part)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "append_file" : false
            }  )" );

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["variable_name"];
        rParameters["model_part_name"];
        rParameters["append_file"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mModelPartName = rParameters["model_part_name"].GetString();
        mVariableName = rParameters["variable_name"].GetString();
        mAppendFile   = rParameters["append_file"].GetBool();
        mTimeUnitConverter = model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];

        KRATOS_CATCH("")
    }

    ~ApplyWriteScalarProcess() override = default;
    ApplyWriteScalarProcess(const ApplyWriteScalarProcess&) = delete;
    ApplyWriteScalarProcess& operator=(const ApplyWriteScalarProcess&) = delete;
    ApplyWriteScalarProcess(ApplyWriteScalarProcess&&) = delete;
    ApplyWriteScalarProcess& operator=(ApplyWriteScalarProcess&&) = delete;

    /// Execute method is used to execute the ApplyWriteScalarProcess algorithms.
    void Execute() override
    {
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY

        const SizeType nNodes = mrModelPart.NumberOfNodes();

        if (nNodes > 0) {
            const Variable<double> &var = KratosComponents< Variable<double> >::Get(mVariableName);
            const double Time = mrModelPart.GetProcessInfo()[TIME]/mTimeUnitConverter;

            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();
            mOutFile.resize(nNodes);

            for (IndexType i = 0; i<nNodes; ++i) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                const IndexType nodeId = it->Id();
                std::string fileName = mModelPartName + "_" + std::to_string(nodeId) + "_" + mVariableName + ".res";

                if (mAppendFile) {
                    // append instead of overwrite
                    mOutFile[i].open(fileName, std::ios::app);
                } else {
                    // open a new file and overwrite
                    mOutFile[i].open(fileName, std::ios::trunc); // overwrite
                    mOutFile[i] << "Time" << "   " << mVariableName << "\n";
                    const double value = it->FastGetSolutionStepValue(var);
                    mOutFile[i] << Time << "   " << value << "\n";
                }
            }
        }
        
        KRATOS_CATCH("")
    }

    /**
     * @brief This function will be executed at every time step AFTER performing the solve phase
     */
    void ExecuteFinalizeSolutionStep() override
    {
        KRATOS_TRY

        const IndexType nNodes = mrModelPart.NumberOfNodes();

        if (nNodes > 0) {
            const Variable<double> &var = KratosComponents< Variable<double> >::Get(mVariableName);
            const double Time = mrModelPart.GetProcessInfo()[TIME]/mTimeUnitConverter;
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

            for (IndexType i = 0; i<nNodes; ++i) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                const double value = it->FastGetSolutionStepValue(var);
                mOutFile[i] << Time << "   " << value << "\n";
            }
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief This function is designed for being called at the end of the computations
     */
    void ExecuteFinalize() override
    {
        KRATOS_TRY

        for (auto& outFile : mOutFile) {
            outFile.close();
        }

        KRATOS_CATCH("")
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyWriteScalarProcess";
    }

private:
    /// Member Variables
    ModelPart& mrModelPart;
    std::string mVariableName;
    std::string mModelPartName;
    bool mAppendFile;
    double mTimeUnitConverter;
    std::vector<std::ofstream> mOutFile;

}; // Class ApplyWriteScalarProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyWriteScalarProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyWriteScalarProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}