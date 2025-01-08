// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors: Marjan Fathian
//                Wijtze Pieter Kikstra
//                Anne van de Graaf
//

#pragma once

#include "geo_mechanics_application_variables.h"
#include "includes/element.h"
#include "includes/kratos_export_api.h"
#include "processes/process.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyCPhiReductionProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyCPhiReductionProcess);

    ApplyCPhiReductionProcess(ModelPart& rModelPart, const Parameters&)
        : Process(Flags()), mrModelPart(rModelPart)
    {
    }

    void ExecuteInitializeSolutionStep() override;
    void ExecuteFinalizeSolutionStep() override;
    void ExecuteFinalize() override;
    int  Check() override;

private:
    ModelPart& mrModelPart;
    double     mReductionFactor         = 1.0;
    double     mPreviousReductionFactor = 1.0;
    double     mReductionIncrement      = 0.1;

    [[nodiscard]] double GetAndCheckPhi(const Element::PropertiesType& rProp) const;

    [[nodiscard]] double ComputeReducedPhi(double Phi) const;

    [[nodiscard]] double GetAndCheckC(const Element::PropertiesType& rProp) const;

    void SetCPhiAtElement(Element& rElement, double ReducedPhi, double ReducedC) const;

    void SetValueAtElement(Element& rElement, const Variable<Vector>& rVariable, const Vector& rValue) const;

    [[nodiscard]] bool IsStepRestarted() const;
};

} // namespace Kratos