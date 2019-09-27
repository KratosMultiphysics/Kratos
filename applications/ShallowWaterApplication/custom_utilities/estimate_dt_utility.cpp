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
#include "estimate_dt_utility.h"
#include "shallow_water_application_variables.h"


namespace Kratos
{

EstimateDtShallow::EstimateDtShallow(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
) : mrModelPart(rThisModelPart)
{

    Parameters default_parameters(R"({
        "automatic_time_step"   : true,
        "time_step"             : 1.0,
        "courant_number"        : 1.0,
        "consider_froude"       : true,
        "minimum_delta_time"    : 1e-4,
        "maximum_delta_time"    : 1e+6
    })");

    ThisParameters.ValidateAndAssignDefaults(default_parameters);

    mEstimateDt = ThisParameters["automatic_time_step"].GetBool();
    mConstantDt = ThisParameters["time_step"].GetDouble();
    mCourant = ThisParameters["courant_number"].GetDouble();
    mConsiderFroude = ThisParameters["consider_froude"].GetBool();
    mMinDt = ThisParameters["minimum_delta_time"].GetDouble();
    mMaxDt = ThisParameters["maximum_delta_time"].GetDouble();

}


double EstimateDtShallow::EstimateDt() const
{
    if (mEstimateDt) {
        return EstimateTimeStep();
    } else {
        return mConstantDt;
    }
}


double EstimateDtShallow::EstimateTimeStep() const
{
    const double gravity = mrModelPart.GetProcessInfo().GetValue(GRAVITY_Z);
    const auto nodes_begin = mrModelPart.NodesBegin();
    double min_characteristic_time = std::numeric_limits<double>::max();

    #pragma omp parallel for shared(min_characteristic_time)
    for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfNodes()); ++i)
    {
        const double local_thread_min_characteristic_time = NodalCharacteristicTime(*(nodes_begin+i), gravity);
        #pragma omp critical
        {
            min_characteristic_time = std::min(min_characteristic_time, local_thread_min_characteristic_time);
        }
    }

    double current_time = min_characteristic_time * mCourant;

    if (current_time < mMinDt) {current_time = mMinDt;}
    else if (current_time > mMaxDt) {current_time = mMaxDt;}

    return current_time;
}


double EstimateDtShallow::NodalCharacteristicTime(const Node<3>& rNode, double gravity) const
{
    double velocity = norm_2(rNode.FastGetSolutionStepValue(VELOCITY));
    double wave_vel = 0.0;
    const double epsilon = std::numeric_limits<double>::epsilon();
    if (mConsiderFroude)
    {
        wave_vel = std::sqrt(gravity * rNode.FastGetSolutionStepValue(HEIGHT));
    }
    return rNode.FastGetSolutionStepValue(NODAL_H) / (velocity + wave_vel + epsilon);
}

}  // namespace Kratos.
