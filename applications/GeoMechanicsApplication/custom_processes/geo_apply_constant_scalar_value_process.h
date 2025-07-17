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
//

#pragma once

#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @class GeoApplyConstantScalarValueProcess
 * @brief A class to apply a constant scalar value to nodes in a model part for a given variable.
 * @details This process applies a single scalar value to nodes in a model part for a specified
 * variable. When the variable is fixed, it will remain constant throughout the stage. If the
 * variable is not fixed, this process will set the specified value as an initial condition.
 * @author Riccardo Rossi
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoApplyConstantScalarValueProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(GeoApplyConstantScalarValueProcess);

    GeoApplyConstantScalarValueProcess(Model& rModel, Parameters ThisParameters);
    GeoApplyConstantScalarValueProcess(ModelPart& rModelPart, Parameters ThisParameters);
    GeoApplyConstantScalarValueProcess(ModelPart&              rModelPart,
                                       const Variable<double>& rVariable,
                                       double                  DoubleValue,
                                       bool                    IsFixed);
    GeoApplyConstantScalarValueProcess(ModelPart& rModelPart, const Variable<int>& rVariable, int IntValue, bool IsFixed);

    GeoApplyConstantScalarValueProcess(ModelPart& rModelPart, const Variable<bool>& rVariable, bool BoolValue, bool IsFixed);

    ~GeoApplyConstantScalarValueProcess() override = default;

    void ExecuteInitialize() override;
    void ExecuteFinalize() override;

    [[nodiscard]] std::string Info() const override { return "GeoApplyConstantScalarValueProcess"; }

protected:
    ModelPart&  mrModelPart;
    std::string mVariableName;

private:
    double mDoubleValue = 0.0;
    int    mIntValue    = 0;
    bool   mBoolValue   = false;
    bool   mIsFixed     = false;

    template <class TVarType>
    void InternalApplyValue(const TVarType& rVariable, bool ToBeFixed, typename TVarType::Type Value);

    template <class TVarType>
    void InternalApplyValueWithoutFixing(const TVarType& rVariable, typename TVarType::Type Value);
};

inline std::ostream& operator<<(std::ostream& rOStream, const GeoApplyConstantScalarValueProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.