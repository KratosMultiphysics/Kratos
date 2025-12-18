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
#include "processes/process.h"

namespace Kratos
{
class ModelPart;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyConstantBoundaryPhreaticSurfacePressureProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyConstantBoundaryPhreaticSurfacePressureProcess);

    ApplyConstantBoundaryPhreaticSurfacePressureProcess(ModelPart& model_part, Parameters rParameters);
    ApplyConstantBoundaryPhreaticSurfacePressureProcess(const ApplyConstantBoundaryPhreaticSurfacePressureProcess&) = delete;
    ApplyConstantBoundaryPhreaticSurfacePressureProcess& operator=(
        const ApplyConstantBoundaryPhreaticSurfacePressureProcess&) = delete;
    ~ApplyConstantBoundaryPhreaticSurfacePressureProcess() override = default;

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override;

    /// Turn back information as a string.
    std::string Info() const override;

protected:
    /// Member Variables
    ModelPart&   mrModelPart;
    std::string  mVariableName;
    bool         mIsFixed;
    bool         mIsFixedProvided;
    unsigned int mGravityDirection;
    double       mSpecificWeight;
    Vector3      mFirstReferenceCoordinate;
    Vector3      mSecondReferenceCoordinate;
    Vector3      mThirdReferenceCoordinate;
    Vector3      mNormalVector;
    double       mEqRHS;

private:
    void calculateEquationParameters();
};

} // namespace Kratos