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
#include "includes/model_part.h"
#include "shallow_water_application_variables.h"


namespace Kratos
{

EstimateTimeStepUtility::EstimateTimeStepUtility(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
) : mrModelPart(rThisModelPart)
{

    Parameters default_parameters(R"({
        "automatic_time_step"   : true,
        "time_step"             : 1.0,
        "courant_number"        : 1.0,
        "minimum_delta_time"    : 1e-4,
        "maximum_delta_time"    : 1e+6
    })");

    ThisParameters.ValidateAndAssignDefaults(default_parameters);

    mEstimateDt = ThisParameters["automatic_time_step"].GetBool();
    mConstantDt = ThisParameters["time_step"].GetDouble();
    mCourant = ThisParameters["courant_number"].GetDouble();
    mMinDt = ThisParameters["minimum_delta_time"].GetDouble();
    mMaxDt = ThisParameters["maximum_delta_time"].GetDouble();

}

double EstimateTimeStepUtility::Execute() const
{
    if (mEstimateDt) {
        return EstimateTimeStep();
    } else {
        return mConstantDt;
    }
}

double EstimateTimeStepUtility::EstimateTimeStep() const
{
    const double gravity = mrModelPart.GetProcessInfo().GetValue(GRAVITY_Z);
    const auto elements_begin = mrModelPart.ElementsBegin();
    double min_characteristic_time = std::numeric_limits<double>::max();

    #pragma omp parallel for shared(min_characteristic_time)
    for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfElements()); ++i)
    {
        const auto& r_geometry = (elements_begin + i)->GetGeometry();
        const double local_thread_min_characteristic_time = ElementCharacteristicTime(r_geometry, gravity);
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

double EstimateTimeStepUtility::ElementCharacteristicTime(const GeometryType& rGeometry, double Gravity) const
{
    array_1d<double,3> velocity = ZeroVector(3);
    double height = 0.0;
    for (auto& r_node : rGeometry)
    {
        velocity += r_node.FastGetSolutionStepValue(VELOCITY);
        height += r_node.FastGetSolutionStepValue(HEIGHT);
    }
    const double lambda = norm_2(velocity) + std::sqrt(Gravity * height);
    const double epsilon = std::numeric_limits<double>::epsilon();
    return rGeometry.Length() / (lambda + epsilon);
}

}  // namespace Kratos.
