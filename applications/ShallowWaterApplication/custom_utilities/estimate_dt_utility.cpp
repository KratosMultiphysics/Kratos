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
#include "utilities/openmp_utils.h"
#include "processes/find_nodal_h_process.h"
#include "shallow_water_application_variables.h"


namespace Kratos
{

EstimateDtShallow::EstimateDtShallow(
    ModelPart& rThisModelPart,
    Parameters& rThisParameters
) : mrModelPart(rThisModelPart)
{

    Parameters default_parameters(R"({
        "automatic_time_step"   : true,
        "time_step"             : 1.0,
        "courant_number"        : 1.0,
        "consider_froude"       : true,
        "compute_nodal_h"       : false,
        "minimum_delta_time"    : 1e-4,
        "maximum_delta_time"    : 1e+6
    })");

    rThisParameters.ValidateAndAssignDefaults(default_parameters);

    mEstimateDt = rThisParameters["automatic_time_step"].GetBool();
    mConstantDt = rThisParameters["time_step"].GetDouble();
    mCourant = rThisParameters["courant_number"].GetDouble();
    mConsiderFroude = rThisParameters["consider_froude"].GetBool();
    mComputeNodalH = rThisParameters["compute_nodal_h"].GetBool();
    mMinDt = rThisParameters["minimum_delta_time"].GetDouble();
    mMaxDt = rThisParameters["maximum_delta_time"].GetDouble();
    mEpsilon = std::numeric_limits<double>::epsilon();

    if (mComputeNodalH)
    {
        FindNodalHProcess<true>(mrModelPart)();
    }
}


double EstimateDtShallow::EstimateDt()
{
    double time_step = mConstantDt;
    if (mEstimateDt)
    {
        time_step = EstimateTimeStep();
    }
    return time_step;
}


double EstimateDtShallow::EstimateTimeStep()
{
    mGravity = mrModelPart.GetProcessInfo().GetValue(GRAVITY_Z);

    unsigned int num_threads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::PartitionVector NodesPartition;
    OpenMPUtils::DivideInPartitions(mrModelPart.NumberOfNodes(),num_threads, NodesPartition);

    std::vector<double> min_characteristic_time(num_threads, 0.0);

    #pragma omp parallel shared(min_characteristic_time)
    {
        int k = OpenMPUtils::ThisThread();
        ModelPart::NodeIterator NodeBegin = mrModelPart.NodesBegin() + NodesPartition[k];
        ModelPart::NodeIterator NodeEnd = mrModelPart.NodesBegin() + NodesPartition[k+1];

        double min_thread_time = std::numeric_limits<double>::infinity();

        for(ModelPart::NodeIterator itNode = NodeBegin; itNode != NodeEnd; ++itNode)
        {
            double characteristic_time = NodalCharacteristicTime(*itNode);
            if (characteristic_time < min_thread_time)
            {
                min_thread_time = characteristic_time;
            }
        }

        min_characteristic_time[k] = min_thread_time;
    }

    // Reduce to maximum the thread results
    // Note that MSVC14 does not support max reductions, which are part of OpenMP 3.1
    double current_time = min_characteristic_time[0];
    for (unsigned int k = 1; k < num_threads; k++)
    {
        if (min_characteristic_time[k] < current_time)
        {
            current_time = min_characteristic_time[k];
        }
    }
    current_time *= mCourant;

    if (current_time < mMinDt) {current_time = mMinDt;}
    else if (current_time > mMaxDt) {current_time = mMaxDt;}

    return current_time;
}


double EstimateDtShallow::NodalCharacteristicTime(Node<3>& rNode)
{
    double velocity = norm_2(rNode.FastGetSolutionStepValue(VELOCITY));
    double wave_vel = 0.0;
    if (mConsiderFroude)
    {
        wave_vel = std::sqrt(mGravity * rNode.FastGetSolutionStepValue(HEIGHT));
    }
    return rNode.FastGetSolutionStepValue(NODAL_H) / (velocity + wave_vel + mEpsilon);
}

}  // namespace Kratos.
