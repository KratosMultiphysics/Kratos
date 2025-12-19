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
#include "processes/process.h"

namespace Kratos
{

class ModelPart;
class Parameters;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyWriteScalarProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyWriteScalarProcess);

    /// Defining a table with double argument and result type as table type.
    using TableType = Table<double, double>;
    using IndexType = std::size_t;

    ApplyWriteScalarProcess(ModelPart& model_part, Parameters rParameters);
    ~ApplyWriteScalarProcess() override                                = default;
    ApplyWriteScalarProcess(const ApplyWriteScalarProcess&)            = delete;
    ApplyWriteScalarProcess& operator=(const ApplyWriteScalarProcess&) = delete;
    ApplyWriteScalarProcess(ApplyWriteScalarProcess&&)                 = delete;
    ApplyWriteScalarProcess& operator=(ApplyWriteScalarProcess&&)      = delete;

    /// Execute method is used to execute the ApplyWriteScalarProcess algorithms.
    void Execute() override;

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override;

    /**
     * @brief This function will be executed at every time step AFTER performing the solve phase
     */
    void ExecuteFinalizeSolutionStep() override;

    /**
     * @brief This function is designed for being called at the end of the computations
     */
    void ExecuteFinalize() override;

    /// Turn back information as a string.
    std::string Info() const override;

private:
    /// Member Variables
    ModelPart&                 mrModelPart;
    std::string                mVariableName;
    std::string                mModelPartName;
    bool                       mAppendFile;
    double                     mTimeUnitConverter;
    std::vector<std::ofstream> mOutFile;

}; // Class ApplyWriteScalarProcess

/// input stream function
inline std::istream& operator>>(std::istream& rIStream, ApplyWriteScalarProcess& rThis);

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream, const ApplyWriteScalarProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos