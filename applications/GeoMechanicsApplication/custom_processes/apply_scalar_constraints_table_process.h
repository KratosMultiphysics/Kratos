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


class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyScalarConstraintsTableProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyScalarConstraintsTableProcess);

    ApplyScalarConstraintsTableProcess(ModelPart&        rModelPart,
                                       const Parameters& rProcessSettings);

    ~ApplyScalarConstraintsTableProcess() override;

    ApplyScalarConstraintsTableProcess(const ApplyScalarConstraintsTableProcess&) = delete;
    ApplyScalarConstraintsTableProcess& operator=(const ApplyScalarConstraintsTableProcess&) = delete;

    using ProcessUniquePointer = std::unique_ptr<Process>;

    void ExecuteInitialize() override;
    void ExecuteInitializeSolutionStep() override;

private:
    void MakeInternalProcess(const Parameters& rProcessSettings);
    void MakeProcessForFluidPressureType(const Parameters&        rProcessSettings,
                                         std::vector<std::string> NamesOfSettingsToCopy);
    void MakeScalarConstraintsProcess(const Parameters&        rProcessSettings,
                                      std::vector<std::string> NamesOfSettingsToCopy);
    void MakeProcessForHydrostaticFluidPressure(const Parameters&        rProcessSettings,
                                                std::vector<std::string> NamesOfSettingsToCopy);
    void MakeProcessForPhreaticLine(const Parameters&        rProcessSettings,
                                    std::vector<std::string> NamesOfSettingsToCopy);
    void MakeProcessForPhreaticSurface(const Parameters&        rProcessSettings,
                                       std::vector<std::string> NamesOfSettingsToCopy);
    void MakeProcessForInterpolatedLine(const Parameters&        rProcessSettings,
                                        std::vector<std::string> NamesOfSettingsToCopy);

    ModelPart& mrModelPart;
    ProcessUniquePointer mProcess;
};

}
