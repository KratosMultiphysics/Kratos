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

#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{
class ModelPart;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyConstantBoundaryHydrostaticPressureProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyConstantBoundaryHydrostaticPressureProcess);

    ApplyConstantBoundaryHydrostaticPressureProcess(ModelPart& model_part, Parameters rParameters);
    ApplyConstantBoundaryHydrostaticPressureProcess(const ApplyConstantBoundaryHydrostaticPressureProcess&) = delete;
    ApplyConstantBoundaryHydrostaticPressureProcess& operator=(const ApplyConstantBoundaryHydrostaticPressureProcess&) = delete;
    ~ApplyConstantBoundaryHydrostaticPressureProcess() override = default;

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override;

    /// Turn back information as a string.
    std::string Info() const override;

    ModelPart& GetModelPart();

    [[nodiscard]] const std::string& GetVariableName() const;

    [[nodiscard]] unsigned int GetGravityDirection() const;

    [[nodiscard]] double GetReferenceCoordinate() const;

    [[nodiscard]] double GetSpecificWeight() const;

private:
    ModelPart&   mrModelPart;
    std::string  mVariableName;
    bool         mIsFixed;
    bool         mIsFixedProvided;
    unsigned int mGravityDirection;
    double       mReferenceCoordinate;
    double       mSpecificWeight;
};

} // namespace Kratos