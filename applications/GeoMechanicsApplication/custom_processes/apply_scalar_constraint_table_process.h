// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf,
//                   Marjan Fathian
//
#pragma once

#include "processes/process.h"

#include <string>

namespace Kratos
{
class Model;
class ModelPart;
class Parameters;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyScalarConstraintTableProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyScalarConstraintTableProcess);

    ApplyScalarConstraintTableProcess(Model& rModel, const Parameters& rProcessSettings);

    ~ApplyScalarConstraintTableProcess() override = default;

    ApplyScalarConstraintTableProcess(const ApplyScalarConstraintTableProcess&)            = delete;
    ApplyScalarConstraintTableProcess& operator=(const ApplyScalarConstraintTableProcess&) = delete;

    using ProcessUniquePointer = std::unique_ptr<Process>;

    void        ExecuteInitialize() override;
    void        ExecuteInitializeSolutionStep() override;
    void        ExecuteFinalize() override;
    std::string Info() const override;

private:
    void MakeInternalProcess(ModelPart& rModelPart, const Parameters& rProcessSettings);
    void MakeProcessForFluidPressureType(ModelPart&               rModelPart,
                                         const Parameters&        rProcessSettings,
                                         std::vector<std::string> NamesOfSettingsToCopy);
    void MakeScalarConstraintProcess(ModelPart&               rModelPart,
                                     const Parameters&        rProcessSettings,
                                     std::vector<std::string> NamesOfSettingsToCopy);
    void MakeProcessForHydrostaticFluidPressure(ModelPart&               rModelPart,
                                                const Parameters&        rProcessSettings,
                                                std::vector<std::string> NamesOfSettingsToCopy);
    void MakeProcessForPhreaticLine(ModelPart&               rModelPart,
                                    const Parameters&        rProcessSettings,
                                    std::vector<std::string> NamesOfSettingsToCopy);
    void MakeProcessForPhreaticMultiLine(ModelPart&               rModelPart,
                                         const Parameters&        rProcessSettings,
                                         std::vector<std::string> NamesOfSettingsToCopy);
    void MakeProcessForPhreaticSurface(ModelPart&               rModelPart,
                                       const Parameters&        rProcessSettings,
                                       std::vector<std::string> NamesOfSettingsToCopy);
    void MakeProcessForInterpolatedLine(ModelPart&               rModelPart,
                                        const Parameters&        rProcessSettings,
                                        std::vector<std::string> NamesOfSettingsToCopy);
    void AppendOptionalFluidParameters(const Parameters&         rProcessSettings,
                                       std::vector<std::string>& rNamesOfSettingsToCopy) const;

    template <typename TableProcessType, typename ConstantProcessType>
    void       InstantiateProcessByTablePresence(ModelPart&                 rModelPart,
                                                 const Parameters&          rProcessSettings,
                                                 std::vector<std::string>&& rNamesOfSettingsToCopy);
    Parameters PrepareProcessParameters(const ModelPart&           rModelPart,
                                        const Parameters&          rProcessSettings,
                                        std::vector<std::string>&& rNamesOfSettingsToCopy);

    std::vector<std::unique_ptr<Process>> mProcesses;
};

} // namespace Kratos
