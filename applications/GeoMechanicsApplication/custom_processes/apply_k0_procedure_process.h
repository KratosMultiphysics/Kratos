// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//

#pragma once

#include "processes/process.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{

class Element;
class ModelPart;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyK0ProcedureProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyK0ProcedureProcess);

    ApplyK0ProcedureProcess(ModelPart& model_part, const Parameters& rK0Settings);
    ~ApplyK0ProcedureProcess() override = default;

    void ExecuteInitialize() override;
    void ExecuteFinalize() override;

    void ExecuteFinalizeSolutionStep() override;
    std::string Info() const override;

private:
    [[nodiscard]] bool UseStandardProcedure() const;
    void CalculateK0Stresses(Element& rElement);

    ModelPart& mrModelPart;
    const Parameters mSettings;
};

} // namespace Kratos
