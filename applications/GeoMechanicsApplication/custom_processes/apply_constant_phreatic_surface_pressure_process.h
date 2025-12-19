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

#include "geo_mechanics_application_variables.h"
#include "processes/process.h"

namespace Kratos
{

class ModelPart;
class Node;
class Parameters;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyConstantPhreaticSurfacePressureProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyConstantPhreaticSurfacePressureProcess);

    ApplyConstantPhreaticSurfacePressureProcess(ModelPart& model_part, Parameters rParameters);
    ApplyConstantPhreaticSurfacePressureProcess(const ApplyConstantPhreaticSurfacePressureProcess&) = delete;
    ApplyConstantPhreaticSurfacePressureProcess& operator=(const ApplyConstantPhreaticSurfacePressureProcess&) = delete;
    ~ApplyConstantPhreaticSurfacePressureProcess() override = default;

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
    Vector3      mFirstReferenceCoordinate;
    Vector3      mSecondReferenceCoordinate;
    Vector3      mThirdReferenceCoordinate;
    Vector3      mNormalVector;
    double       mEqRHS;
    double       mPressureTensionCutOff;

private:
    void calculateEquationParameters();
};

} // namespace Kratos