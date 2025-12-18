// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Richard Faasse
//

#pragma once

#include "processes/process.h"

namespace Kratos
{

class ModelPart;

///@name Kratos Classes
///@{

/**
 * @class GeoApplyConstantScalarValueProcess
 * @brief A class to apply a constant scalar value to nodes in a model part for a given variable.
 * @details This process applies a single scalar value to all nodes in a model part for a specified
 * variable. When the variable is fixed, it will remain constant throughout the stage. If the
 * variable is not fixed, this process will set the specified value as an initial condition.
 * @author Riccardo Rossi
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoApplyConstantScalarValueProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(GeoApplyConstantScalarValueProcess);

    GeoApplyConstantScalarValueProcess(ModelPart& rModelPart, const Parameters& rParameters);
    ~GeoApplyConstantScalarValueProcess() override = default;

    void                      ExecuteInitialize() override;
    void                      ExecuteInitializeSolutionStep() override;
    void                      ExecuteFinalize() override;
    [[nodiscard]] std::string Info() const override;

protected:
    ModelPart&  mrModelPart;
    std::string mVariableName;

private:
    double mDoubleValue   = 0.0;
    int    mIntValue      = 0;
    bool   mBoolValue     = false;
    bool   mIsFixed       = false;
    bool   mIsInitialized = false;
};

} // namespace Kratos.