// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse,
//                   Anne van de Graaf,
//                   Gennady Markelov
//

#pragma once

#include "processes/process.h"

#include <memory>
#include <string>
#include <vector>

namespace Kratos
{

class ModelPart;
class Parameters;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyNormalLoadTableProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyNormalLoadTableProcess);

    ApplyNormalLoadTableProcess(ModelPart& rModelPart, const Parameters& rProcessSettings);

    ~ApplyNormalLoadTableProcess() override = default;

    ApplyNormalLoadTableProcess(const ApplyNormalLoadTableProcess&)            = delete;
    ApplyNormalLoadTableProcess& operator=(const ApplyNormalLoadTableProcess&) = delete;

    void        ExecuteInitialize() override;
    void        ExecuteInitializeSolutionStep() override;
    void        ExecuteFinalize() override;
    std::string Info() const override;

private:
    void        MakeInternalProcesses(const Parameters& rProcessSettings);
    void        MakeProcessForNormalComponent(const Parameters& rProcessSettings);
    void        MakeProcessForTangentialComponent(const Parameters& rProcessSettings);
    void        MakeProcessForUniformFluidPressureType(const Parameters& rProcessSettings,
                                                       const std::vector<std::string>& NamesOfSettingsToCopy);
    void        MakeProcessForHydrostaticFluidPressureType(const Parameters&        rProcessSettings,
                                                           std::vector<std::string> NamesOfSettingsToCopy);
    void        MakeProcessForPhreaticLineFluidPressureType(const Parameters& rProcessSettings,
                                                            std::vector<std::string> NamesOfSettingsToCopy);
    void        MakeProcessForPhreaticSurfaceFluidPressureType(const Parameters& rProcessSettings,
                                                               std::vector<std::string> NamesOfSettingsToCopy);
    static bool IsNormalComponentActive(const Parameters& rProcessSettings);
    static bool IsTangentialComponentActive(const Parameters& rProcessSettings);
    static bool IsComponentActive(const Parameters& rProcessSettings, int componentNumber);

    ModelPart&                            mrModelPart;
    std::vector<std::unique_ptr<Process>> mProcesses;
};

} // namespace Kratos
