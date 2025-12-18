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

#include "geo_mechanics_application_variables.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "processes/process.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyConstantBoundaryPhreaticLinePressureProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyConstantBoundaryPhreaticLinePressureProcess);

    /// Constructor
    ApplyConstantBoundaryPhreaticLinePressureProcess(ModelPart& model_part, Parameters rParameters);
    ApplyConstantBoundaryPhreaticLinePressureProcess(const ApplyConstantBoundaryPhreaticLinePressureProcess&) = delete;
    ApplyConstantBoundaryPhreaticLinePressureProcess& operator=(const ApplyConstantBoundaryPhreaticLinePressureProcess&) = delete;
    ~ApplyConstantBoundaryPhreaticLinePressureProcess() override = default;

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override;

    /// Turn back information as a string.
    std::string Info() const override;

protected:
    [[nodiscard]] const ModelPart& GetModelPart() const;

    [[nodiscard]] ModelPart& GetModelPart();

    [[nodiscard]] const std::string& GetVariableName() const;

    [[nodiscard]] unsigned int GetGravityDirection() const;

    [[nodiscard]] unsigned int GetHorizontalDirection() const;

    [[nodiscard]] double GetSpecificWeight() const;

    [[nodiscard]] const Vector3& GetFirstReferenceCoordinate() const;

    [[nodiscard]] double GetSlope() const;

    [[nodiscard]] double GetMinHorizontalCoordinate() const;

    [[nodiscard]] double GetMaxHorizontalCoordinate() const;

private:
    ModelPart&   mrModelPart;
    std::string  mVariableName;
    bool         mIsFixed;
    bool         mIsFixedProvided;
    unsigned int mGravityDirection;
    unsigned int mHorizontalDirection;
    double       mSpecificWeight;
    unsigned int mOutOfPlaneDirection;
    Vector3      mFirstReferenceCoordinate;
    Vector3      mSecondReferenceCoordinate;
    double       mSlope;
    double       mMinHorizontalCoordinate;
    double       mMaxHorizontalCoordinate;
};

} // namespace Kratos