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
 * @details This process applies a single scalar value to all nodes in a model part for a specified
 * variable. When the variable is fixed, it will remain constant throughout the stage. If the
 * variable is not fixed, this process will set the specified value as an initial condition.
 * @author Riccardo Rossi
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoApplyConstantScalarValueProcess : public Process
{
public:
    KRATOS_DEFINE_LOCAL_FLAG(VARIABLE_IS_FIXED);

    KRATOS_CLASS_POINTER_DEFINITION(GeoApplyConstantScalarValueProcess);

    GeoApplyConstantScalarValueProcess(Model& rModel, Parameters ThisParameters);
    GeoApplyConstantScalarValueProcess(ModelPart& rModelPart, Parameters ThisParameters);
    GeoApplyConstantScalarValueProcess(ModelPart&              rModelPart,
                                       const Variable<double>& rVariable,
                                       const double            DoubleValue,
                                       const Flags             Options);
    GeoApplyConstantScalarValueProcess(ModelPart&           rModelPart,
                                       const Variable<int>& rVariable,
                                       const int            IntValue,
                                       const Flags          options);

    GeoApplyConstantScalarValueProcess(ModelPart&            rModelPart,
                                       const Variable<bool>& rVariable,
                                       const bool            BoolValue,
                                       const Flags           options);

    ~GeoApplyConstantScalarValueProcess() override = default;

    void             ExecuteInitializeSolutionStep() override;
    void             ExecuteFinalize() override;
    const Parameters GetDefaultParameters() const override;

    [[nodiscard]] std::string Info() const override { return "GeoApplyConstantScalarValueProcess"; }

    void PrintInfo(std::ostream& rOStream) const override { rOStream << Info(); }

    void PrintData(std::ostream& rOStream) const override {}

protected:
    ModelPart&  mrModelPart;          /// Reference to the model part.
    std::string mVariableName;        /// Name of the variable.
    double      mDoubleValue = 0.0;   /// Double value.
    int         mIntValue    = 0;     /// Integer value.
    bool        mBoolValue   = false; /// Boolean value.
    bool        mInitialized = false;

private:
    template <class TVarType>
    void InternalApplyValue(const TVarType& rVariable, const bool ToBeFixed, const typename TVarType::Type Value);

    template <class TVarType>
    void InternalApplyValueWithoutFixing(const TVarType& rVariable, const typename TVarType::Type Value);
};

inline std::ostream& operator<<(std::ostream& rOStream, const GeoApplyConstantScalarValueProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.