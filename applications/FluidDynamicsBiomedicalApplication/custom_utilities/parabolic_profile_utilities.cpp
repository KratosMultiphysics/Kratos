//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Eduardo Soudah
//                   Ruben Zorilla
//

// System includes

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "fluid_dynamics_biomedical_application_variables.h"
#include "parabolic_profile_utilities.h"

namespace Kratos {

void ParabolicProfileUtilities::ImposeParabolicInlet(
    ModelPart& rModelPart,
    const double MaxParabolaValue)
{
    // Impose the parabolic inlet profile
    ImposeParabolicProfile(rModelPart, MaxParabolaValue);
}

void ParabolicProfileUtilities::ImposeParabolicInlet(
    ModelPart& rModelPart,
    const GenericFunctionUtility::Pointer MaxParabolaValue)
{
    // Impose the parabolic inlet profile
    ImposeParabolicProfile(rModelPart, MaxParabolaValue);
}

template<class TInputType>
void ParabolicProfileUtilities::ImposeParabolicProfile(
    ModelPart &rModelPart,
    const TInputType& rMaxParabolaValue)
{
    // Get time value from model part ProcessInfo
    const auto& r_process_info = rModelPart.GetProcessInfo();
    KRATOS_ERROR_IF_NOT(r_process_info.Has(TIME)) << "TIME is not present in '" << rModelPart.FullName() << "' ProcessInfo." << std::endl;
    const double time = r_process_info[TIME];

    // Get the maximum distance to the wall
    const double max_dist = block_for_each<MaxReduction<double>>(rModelPart.Nodes(), [](NodeType& rNode){
        return rNode.GetValue(WALL_DISTANCE);
    });
    KRATOS_ERROR_IF(std::abs(max_dist) < 1.0e-12) << "WALL_DISTANCE is close to zero." << std::endl;

    // Impose the parabolic profile values
    block_for_each(rModelPart.Nodes(), array_1d<double,3>(), [time, max_dist, &rMaxParabolaValue](NodeType& rNode, array_1d<double,3>& rTLSValue){
        //Calculate distance for each node
        const double wall_dist = rNode.GetValue(WALL_DISTANCE) < 0.0 ? 0.0 : rNode.GetValue(WALL_DISTANCE);

        // Calculate the inlet direction
        const auto& r_n = rNode.GetValue(INLET_NORMAL);
        const double n_norm = norm_2(r_n);
        if (n_norm > 1.0e-12) {
            rTLSValue = -r_n / n_norm;
        } else {
            KRATOS_WARNING("ImposeParabolicProfile") << "Node " << rNode.Id() << " INLET_NORMAL is close to zero." << std::endl;
            rTLSValue = -r_n;
        }

        // Calculate the inlet value module
        const double value_in = GetMaxParabolaValue(time, rNode, rMaxParabolaValue) * (1.0-(std::pow(max_dist-wall_dist,2)/std::pow(max_dist,2)));
        rTLSValue *= value_in;

        // Set and fix the VELOCITY DOFs
        rNode.FastGetSolutionStepValue(VELOCITY) = rTLSValue;
        rNode.Fix(VELOCITY_X);
        rNode.Fix(VELOCITY_Y);
        rNode.Fix(VELOCITY_Z);
    });
}

void ParabolicProfileUtilities::FreeParabolicInlet(ModelPart& rModelPart)
{
    block_for_each(rModelPart.Nodes(), [](NodeType& rNode){
        rNode.Free(VELOCITY_X);
        rNode.Free(VELOCITY_Y);
        rNode.Free(VELOCITY_Z);
    });
}

template<>
double ParabolicProfileUtilities::GetMaxParabolaValue<double>(
    const double Time,
    const NodeType& rNode,
    const double& rMaxParabolaValue)
{
    return rMaxParabolaValue;
}

template<>
double ParabolicProfileUtilities::GetMaxParabolaValue<GenericFunctionUtility::Pointer>(
    const double Time,
    const NodeType& rNode,
    const GenericFunctionUtility::Pointer& rMaxParabolaValue)
{
    if(!rMaxParabolaValue->UseLocalSystem()) {
        return rMaxParabolaValue->CallFunction(rNode.X(), rNode.Y(), rNode.Z(), Time, rNode.X0(), rNode.Y0(), rNode.Z0());
    } else {
        return rMaxParabolaValue->RotateAndCallFunction(rNode.X(), rNode.Y(), rNode.Z(), Time, rNode.X0(), rNode.Y0(), rNode.Z0());
    }
}

template KRATOS_API(FLUID_DYNAMICS_BIOMEDICAL_APPLICATION) void ParabolicProfileUtilities::ImposeParabolicProfile<double>(ModelPart&, const double&);
template KRATOS_API(FLUID_DYNAMICS_BIOMEDICAL_APPLICATION) void ParabolicProfileUtilities::ImposeParabolicProfile<GenericFunctionUtility::Pointer>(ModelPart&, const GenericFunctionUtility::Pointer&);

} //namespace Kratos
