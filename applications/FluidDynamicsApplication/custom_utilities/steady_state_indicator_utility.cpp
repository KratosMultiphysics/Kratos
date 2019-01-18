//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//                   Jordi Cotela
//

#include "includes/variables.h"
#include "includes/kratos_parameters.h"
#include "steady_state_indicator_utility.h"


namespace Kratos {
    SteadyStateIndicatorUtility::SteadyStateIndicatorUtility(ModelPart &ModelPart):
      mrModelPart(ModelPart){
        mChangeInVelocity = 0.0;
        mChangeInPressure = 0.0;
        mIsSteady = false;
    }

    void SteadyStateIndicatorUtility::EstimateQuantityChangesInTime(){
        int NumNodes = mrModelPart.NumberOfNodes();
        double DomainArea = 0.0;
        for (int i = 0; i < NumNodes; ++i)
        {
            auto inode = mrModelPart.NodesBegin() + i;
            double NodalArea = inode->FastGetSolutionStepValue(NODAL_AREA);
            mChangeInVelocity += this->NodalVelocityChange((*inode)) * NodalArea;
            mChangeInPressure += this->NodalPressureChange((*inode)) * NodalArea;
            DomainArea += NodalArea;
        }
        //std::cout << "mChangeInVelocity" << mChangeInVelocity << std::endl;
        mChangeInVelocity = mChangeInVelocity / DomainArea;
        mChangeInPressure = mChangeInPressure / DomainArea;
    }
    
    double SteadyStateIndicatorUtility::NodalVelocityChange(ModelPart::NodeType& inode){
        array_1d<double, 3 > prev_ave_vel = inode.FastGetSolutionStepValue(TIME_AVERAGED_VELOCITY,2);
        array_1d<double, 3 > averaged_vel = inode.FastGetSolutionStepValue(TIME_AVERAGED_VELOCITY);
        double NodalChangeInVelocity = 0.0;
        for (int i = 0; i<3; ++i){
            NodalChangeInVelocity += (averaged_vel[i]-prev_ave_vel[i])*(averaged_vel[i]-prev_ave_vel[i]);}
        NodalChangeInVelocity = sqrt(NodalChangeInVelocity);
        return NodalChangeInVelocity;
    }

    double SteadyStateIndicatorUtility::NodalPressureChange(ModelPart::NodeType& inode){
        double prev_ave_p = inode.FastGetSolutionStepValue(TIME_AVERAGED_PRESSURE,2);
        double averaged_p = inode.FastGetSolutionStepValue(TIME_AVERAGED_PRESSURE);
        return sqrt((averaged_p - prev_ave_p)*(averaged_p - prev_ave_p));
    }

}