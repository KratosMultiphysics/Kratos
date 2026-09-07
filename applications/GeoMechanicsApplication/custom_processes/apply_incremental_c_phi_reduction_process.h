// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//

#pragma once

#include "includes/element.h"
#include "includes/kratos_export_api.h"
#include "processes/process.h"

#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace Kratos
{


class Model;
class Parameters;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyIncrementalCPhiReductionProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyIncrementalCPhiReductionProcess);

    ApplyIncrementalCPhiReductionProcess(Model& rModel, const Parameters& rProcessSettings);
    void                      ExecuteInitializeSolutionStep() override;
    void                      ExecuteFinalizeSolutionStep() override;
    void                      ExecuteFinalize() override;
    int                       Check() override;
    [[nodiscard]] std::string Info() const override;
    [[nodiscard]] double      GetTrialFactor() const;
    [[nodiscard]] bool        HasLastConvergedFactor() const;
    [[nodiscard]] double      GetLastConvergedFactor() const;
    [[nodiscard]] std::size_t GetTrialNumber() const;
    [[nodiscard]] bool        IsFinished() const;

private:
    std::vector<std::reference_wrapper<ModelPart>> mrModelParts;

    std::string mIncrementStrategy;
    double mTrialFactor;
    std::optional<double> mLastConvergedFactor;
    double mFactorIncrement;

    std::size_t mTrialNumber;
    std::size_t mMaxTrials;

    [[nodiscard]] double CalculateNewFactor() const;

    void ApplyReduction();

    void AdvanceFactor();

    [[nodiscard]] static double GetAndCheckCohesion(const Properties& rProperties);

    [[nodiscard]] static double GetAndCheckFrictionAngle(const Properties& rProperties,
                                                         IndexType         PropertyId);

    [[nodiscard]] static double CalculateReducedCohesion(
        double Cohesion,
        double Factor);

    [[nodiscard]] static double CalculateReducedFrictionAngle(
        double FrictionAngle,
        double Factor);

    static void SetCPhiAtElement(Element& rElement, double ReducedFrictionAngle, double ReducedCohesion);

    static void InitializeParametersForInternalMohrCoulombModel(Element& rElement);
};

} // namespace Kratos
