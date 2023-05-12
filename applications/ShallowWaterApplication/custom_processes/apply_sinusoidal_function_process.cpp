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
    mDirection = rThisParameters["direction"].GetVector();
    mDirection /= norm_2(mDirection);
    mAmplitude = rThisParameters["amplitude"].GetDouble();
    const double period = rThisParameters["period"].GetDouble();
    const double wavelength = rThisParameters["wavelength"].GetDouble();
    mFrequency = 2 * pi / period;
    mWavenumber = 2 * pi / wavelength;
    mPhase = rThisParameters["phase"].GetDouble();
    mShift = rThisParameters["shift"].GetDouble();
    mSmoothTime = std::max(rThisParameters["smooth_time"].GetDouble(), std::numeric_limits<double>::epsilon());
    mSmoothTimes = rThisParameters["smooth_time_centers"].GetVector();
}


template< class TVarType >
int ApplySinusoidalFunctionProcess<TVarType>::Check()
{
    if (mrModelPart.Nodes().size() != 0) {
        const auto& r_node = *mrModelPart.NodesBegin();
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(mrVariable, r_node);
    }
    KRATOS_CHECK(mFrequency < std::numeric_limits<double>::max());
    KRATOS_CHECK(mWavenumber < std::numeric_limits<double>::max());
    KRATOS_CHECK(mFrequency > 0.0);
    KRATOS_CHECK(mWavenumber > 0.0);
    KRATOS_CHECK(mDirection.size() == 3);
    KRATOS_CHECK(norm_2(mDirection) > 0.0);
    return 0;
}


template< class TVarType >
void ApplySinusoidalFunctionProcess<TVarType>::ExecuteInitializeSolutionStep()
{
    double time = mrModelPart.GetProcessInfo().GetValue(TIME);
    const double smooth = SmoothFactor(time);
    block_for_each(mrModelPart.Nodes(), [&](NodeType& rNode){
        double value = smooth * Function(rNode, time);
        rNode.FastGetSolutionStepValue(mrVariable) = value;
    });
}


template<>
void ApplySinusoidalFunctionProcess<Variable<array_1d<double,3>>>::ExecuteInitializeSolutionStep()
{
    const double time = mrModelPart.GetProcessInfo().GetValue(TIME);
    const double smooth = SmoothFactor(time);
    block_for_each(mrModelPart.Nodes(), [&](NodeType& rNode){
        double modulus = smooth * Function(rNode, time);
        noalias(rNode.FastGetSolutionStepValue(mrVariable)) = modulus * mDirection;
    });
}


template<class TVarType>
double ApplySinusoidalFunctionProcess<TVarType>::Function(const array_1d<double,3>& rCoordinates, const double& rTime)
{
    const double x = inner_prod(mDirection, rCoordinates);
    return mAmplitude * std::sin(mFrequency * rTime - mWavenumber * x + mPhase) + mShift;
}


template<class TVarType>
double ApplySinusoidalFunctionProcess<TVarType>::SmoothFactor(const double& rTime)
{
    const double pi = std::acos(-1);
    double smooth = 1.0;
    for (auto origin : mSmoothTimes)
    {
        const double distance = std::abs(rTime - origin);
        smooth = std::min(smooth, 2 * std::atan(distance / mSmoothTime) / pi);
    }
    return smooth;
}


template<class TVarType>
const Parameters ApplySinusoidalFunctionProcess<TVarType>::GetDefaultParameters() const
{
    Parameters default_parameters = Parameters(R"(
    {
        "direction"           : [1.0, 0.0, 0.0],
        "amplitude"           : 1.0,
        "period"              : 1.0,
        "wavelength"          : 1.0,
        "phase"               : 0.0,
        "shift"               : 0.0,
        "smooth_time"         : 0.0,
        "smooth_time_centers" : [0.0]
    })");
    return default_parameters;
}


template class ApplySinusoidalFunctionProcess< Variable<double> >;
template class ApplySinusoidalFunctionProcess< Variable< array_1d<double, 3> > >;

}  // namespace Kratos.
