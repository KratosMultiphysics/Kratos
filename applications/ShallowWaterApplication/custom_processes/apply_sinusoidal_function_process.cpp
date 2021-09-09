//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "includes/checks.h"
#include "utilities/parallel_utilities.h"
#include "apply_sinusoidal_function_process.h"


namespace Kratos
{

template< class TVarType >
ApplySinusoidalFunctionProcess<TVarType>::ApplySinusoidalFunctionProcess(
    ModelPart& rThisModelPart,
    TVarType& rThisVariable,
    Parameters& rThisParameters)
     : mrModelPart(rThisModelPart)
     , mrVariable(rThisVariable)
{
    rThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    const double pi = std::acos(-1);

    // Initialization of member variables
    mAmplitude = rThisParameters["amplitude"].GetDouble();
    mPeriod = rThisParameters["period"].GetDouble();
    mAngularFrequency = 2 * pi / mPeriod;
    mPhase = rThisParameters["phase_shift"].GetDouble() * mAngularFrequency;
    mVerticalShift = rThisParameters["vertical_shift"].GetDouble();
    mSmoothTime = std::max(rThisParameters["smooth_time"].GetDouble(), std::numeric_limits<double>::epsilon());
}


template< class TVarType >
int ApplySinusoidalFunctionProcess<TVarType>::Check()
{
    if (mrModelPart.Nodes().size() != 0) {
        const auto& r_node = *mrModelPart.NodesBegin();
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(mrVariable, r_node);
    }
    KRATOS_CHECK(mPeriod >= std::numeric_limits<double>::epsilon());
    return 0;
}


template<>
int ApplySinusoidalFunctionProcess< Variable< array_1d<double, 3> > >::Check()
{
    if (mrModelPart.Nodes().size() != 0) {
        const auto& r_node = *mrModelPart.NodesBegin();
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(mrVariable, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NORMAL, r_node);
    }
    KRATOS_CHECK(mPeriod >= std::numeric_limits<double>::epsilon());
    return 0;
}


template< class TVarType >
void ApplySinusoidalFunctionProcess<TVarType>::ExecuteInitializeSolutionStep()
{
    const double pi = std::acos(-1);
    double time = mrModelPart.GetProcessInfo().GetValue(TIME);
    double smooth = 2 * std::atan(time / mSmoothTime) / pi;
    double value = smooth * Function(time);
    block_for_each(mrModelPart.Nodes(), [&](NodeType& rNode){
        rNode.FastGetSolutionStepValue(mrVariable) = value;
    });
}


template<>
void ApplySinusoidalFunctionProcess<Variable<array_1d<double, 3>>>::ExecuteInitializeSolutionStep()
{
    const double pi = std::acos(-1);
    double time = mrModelPart.GetProcessInfo().GetValue(TIME);
    double smooth = 2 * std::atan(time / mSmoothTime) / pi;
    double modulus = smooth * Function(time);
    block_for_each(mrModelPart.Nodes(), [&](NodeType& rNode){
        array_1d<double, 3> direction = rNode.FastGetSolutionStepValue(NORMAL);
        noalias(rNode.FastGetSolutionStepValue(mrVariable)) = -modulus * direction;
    });
}


template< class TVarType >
double ApplySinusoidalFunctionProcess<TVarType>::Function(double& rTime)
{
    return mAmplitude * std::sin(mAngularFrequency * rTime - mPhase) + mVerticalShift;
}


template<class TVarType>
const Parameters ApplySinusoidalFunctionProcess<TVarType>::GetDefaultParameters() const
{
    Parameters default_parameters = Parameters(R"(
    {
        "amplitude"       : 1.0,
        "period"          : 1.0,
        "phase_shift"     : 0.0,
        "vertical_shift"  : 0.0,
        "smooth_time"     : 0.0
    })");
    return default_parameters;
}


template class ApplySinusoidalFunctionProcess< Variable<double> >;
template class ApplySinusoidalFunctionProcess< Variable< array_1d<double, 3> > >;

}  // namespace Kratos.
