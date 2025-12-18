// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#pragma once

#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{
class ModelPart;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyConstantHydrostaticPressureProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyConstantHydrostaticPressureProcess);

    ApplyConstantHydrostaticPressureProcess(ModelPart& model_part, Parameters rParameters);

    ApplyConstantHydrostaticPressureProcess& operator=(const ApplyConstantHydrostaticPressureProcess&) = delete;
    ~ApplyConstantHydrostaticPressureProcess() override = default;

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override;

    const std::string& GetName() const;

    /// Turn back information as a string.
    std::string Info() const override;

protected:
    /// Member Variables
    ModelPart&   mrModelPart;
    std::string  mModelPartName;
    std::string  mVariableName;
    bool         mIsFixed;
    bool         mIsFixedProvided;
    bool         mIsSeepage;
    unsigned int mGravityDirection;
    double       mReferenceCoordinate;
    double       mSpecificWeight;
    double       mPressureTensionCutOff;
};

} // namespace Kratos