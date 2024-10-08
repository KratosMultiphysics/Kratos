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

#include "containers/array_1d.h"
#include "includes/element.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

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
    int  Check() override;

    void                      ExecuteFinalizeSolutionStep() override;
    [[nodiscard]] std::string Info() const override;

private:
    [[nodiscard]] bool  UseStandardProcedure() const;
    array_1d<double, 3> CreateK0Vector(const Element::PropertiesType& rProp) const;
    void                CalculateK0Stresses(Element& rElement);

    ModelPart&       mrModelPart;
    const Parameters mSettings;
};

} // namespace Kratos
