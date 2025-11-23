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

    ApplyK0ProcedureProcess(Model& rModel, Parameters K0Settings);
    ~ApplyK0ProcedureProcess() override = default;

    void ExecuteInitialize() override;
    void ExecuteFinalize() override;
    int  Check() override;

    void                      ExecuteFinalizeSolutionStep() override;
    [[nodiscard]] std::string Info() const override;

private:
    [[nodiscard]] bool                       UseStandardProcedure() const;
    [[nodiscard]] static array_1d<double, 3> CreateK0Vector(const Element::PropertiesType& rProp);
    static void CalculateK0Stresses(Element& rElement, const ProcessInfo& rProcessInfo);
    static void CheckK0(const Properties& rProperties, IndexType ElementId);
    static void CheckK0MainDirection(const Properties& rProperties, IndexType ElementId);
    static void CheckOCRorPOP(const Properties& rProperties, IndexType ElementId);
    static void CheckPhi(const Properties& rProperties, IndexType ElementId);
    static void CheckPoissonUnloadingReloading(const Properties& rProperties, IndexType ElementId);
    static void CheckSufficientMaterialParameters(const Properties& rProperties, IndexType ElementId);

    std::vector<std::reference_wrapper<ModelPart>> mrModelParts;
    const Parameters                               mSettings;
};

} // namespace Kratos
