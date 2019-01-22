//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Mengjie Zhao
//
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "utilities/openmp_utils.h"
#include "processes/find_nodal_h_process.h"

// Application includes
#include "time_averaging_process.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

/* Public functions *******************************************************/


TimeAveragingProcess::TimeAveragingProcess(
    ModelPart& rModelPart)
    : Process(), mrModelPart(rModelPart) {
        std::cout << "TimeAveragingProcess: Initialized." << std::endl;
}

void TimeAveragingProcess::ExecuteInitialize() {
    int NumNodes = mrModelPart.NumberOfNodes();
    for (int i = 0; i < NumNodes; ++i)
    {
        auto inode = mrModelPart.NodesBegin() + i;
        inode->SetValue(TIME_AVERAGED_VELOCITY,inode->FastGetSolutionStepValue(VELOCITY));
        inode->SetValue(TIME_AVERAGED_PRESSURE,inode->FastGetSolutionStepValue(PRESSURE));
    }

}

void TimeAveragingProcess::ExecuteFinalizeSolutionStep() {
    int NumNodes = mrModelPart.NumberOfNodes();
    for (int i = 0; i < NumNodes; ++i)
    {
        ModelPart::NodeIterator inode = mrModelPart.NodesBegin() + i;
        this->AverageVelocity(inode);
        this->AveragePressure(inode);
    }
    std::cout << "TimeAveragingProcess: Nodal Quantities Averaged." << std::endl;
}

/* Protected functions ****************************************************/

void TimeAveragingProcess::AverageVelocity(ModelPart::NodeIterator& inode){
    array_1d<double,3> current_vel = inode->FastGetSolutionStepValue(VELOCITY);
    array_1d<double,3> &prev_ave_vel = inode->FastGetSolutionStepValue(TIME_AVERAGED_VELOCITY,1);
    array_1d<double, 3 > &averaged_vel = inode->FastGetSolutionStepValue(TIME_AVERAGED_VELOCITY);
    double t =  mrModelPart.GetProcessInfo()[TIME];
    double dt =  mrModelPart.GetProcessInfo()[DELTA_TIME];
	int time_step = mrModelPart.GetProcessInfo()[STEP];
    if (time_step == 1){
		noalias(averaged_vel) = current_vel;
    } else {
        noalias(averaged_vel) = ( prev_ave_vel * (t - dt) + current_vel * dt ) / t;
    }
}

void TimeAveragingProcess::AveragePressure(ModelPart::NodeIterator& inode){
    double& current_pre = inode->FastGetSolutionStepValue(PRESSURE);
    double& prev_ave_pre = inode->FastGetSolutionStepValue(TIME_AVERAGED_PRESSURE,1);
    double& averaged_pre = inode->FastGetSolutionStepValue(TIME_AVERAGED_PRESSURE);
    double t =  mrModelPart.GetProcessInfo()[TIME];
    double dt =  mrModelPart.GetProcessInfo()[DELTA_TIME];
    int time_step = mrModelPart.GetProcessInfo()[STEP];
    // when current pressure does not converge at current step, ignore current value for averaging
    if (time_step == 1){
		averaged_pre = current_pre;
    }
    else {
        averaged_pre = ( prev_ave_pre * (t - dt) + current_pre * dt ) / t;
    }
    inode->SetValue(TIME_AVERAGED_PRESSURE,averaged_pre);
}


/* Private functions ****************************************************/

};  // namespace Kratos.
