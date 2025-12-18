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
#include "includes/node.h"
#include "processes/process.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyConstantPhreaticLinePressureProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyConstantPhreaticLinePressureProcess);

    ApplyConstantPhreaticLinePressureProcess(ModelPart& model_part, Parameters rParameters);
    ApplyConstantPhreaticLinePressureProcess(const ApplyConstantPhreaticLinePressureProcess&) = delete;
    ApplyConstantPhreaticLinePressureProcess& operator=(const ApplyConstantPhreaticLinePressureProcess&) = delete;
    ~ApplyConstantPhreaticLinePressureProcess() override = default;

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
    bool         mIsSeepage;
    unsigned int mGravityDirection;
    double       mSpecificWeight;
    unsigned int mOutOfPlaneDirection;
    unsigned int mHorizontalDirection;
    Vector3      mFirstReferenceCoordinate;
    Vector3      mSecondReferenceCoordinate;
    double       mSlope;
    double       mMinHorizontalCoordinate;
    double       mMaxHorizontalCoordinate;
    double       mPressureTensionCutOff;

    double CalculatePressure(const Node& rNode) const;
};

} // namespace Kratos