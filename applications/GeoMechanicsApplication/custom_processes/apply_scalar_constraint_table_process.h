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

namespace Kratos
{

class ModelPart;
class Parameters;


class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyScalarConstraintTableProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyScalarConstraintTableProcess);

    ApplyScalarConstraintTableProcess(ModelPart&        rModelPart,
                                      const Parameters& rProcessSettings);

    ~ApplyScalarConstraintTableProcess() override;

    ApplyScalarConstraintTableProcess(const ApplyScalarConstraintTableProcess&) = delete;
    ApplyScalarConstraintTableProcess& operator=(const ApplyScalarConstraintTableProcess&) = delete;

    using ProcessUniquePointer = std::unique_ptr<Process>;

    void ExecuteInitialize() override;
    void ExecuteInitializeSolutionStep() override;

    std::string Info() const override;

private:
    void MakeInternalProcess(const Parameters& rProcessSettings);
    void MakeProcessForFluidPressureType(const Parameters&        rProcessSettings,
                                         std::vector<std::string> NamesOfSettingsToCopy);
    void MakeScalarConstraintProcess(const Parameters&        rProcessSettings,
                                     std::vector<std::string> NamesOfSettingsToCopy);
    void MakeProcessForHydrostaticFluidPressure(const Parameters&        rProcessSettings,
                                                std::vector<std::string> NamesOfSettingsToCopy);
    void MakeProcessForPhreaticLine(const Parameters&        rProcessSettings,
                                    std::vector<std::string> NamesOfSettingsToCopy);
    void MakeProcessForPhreaticMultiLine(const Parameters&        rProcessSettings,
                                         std::vector<std::string> NamesOfSettingsToCopy);
    void MakeProcessForPhreaticSurface(const Parameters&        rProcessSettings,
                                       std::vector<std::string> NamesOfSettingsToCopy);
    void MakeProcessForInterpolatedLine(const Parameters&        rProcessSettings,
                                        std::vector<std::string> NamesOfSettingsToCopy);

    ModelPart& mrModelPart;
    ProcessUniquePointer mProcess;
};

}
