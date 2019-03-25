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
    ModelPart& rModelPart, 
    const bool is_velocity_averaged,  
    const bool is_pressure_averaged,  
    const bool is_reaction_averaged)
    : Process(), mrModelPart(rModelPart) {
        mDoAveragedVelocity = is_velocity_averaged;
        mDoAveragedPressure = is_pressure_averaged;
        mDoAveragedReaction = is_reaction_averaged;
        
        if (mDoAveragedVelocity == 1){
            std::cout << "TimeAveragingProcess: VELOCITY will be averaged." << std::endl;
            if (TIME_AVERAGED_VELOCITY.Key() == 0){
                rModelPart.AddNodalSolutionStepVariable(TIME_AVERAGED_VELOCITY);
                std::cout << "TimeAveragingProcess: TIME_AVERAGED_VELOCITY added." << std::endl;
            }
        }

        if (mDoAveragedPressure == 1){
            std::cout << "TimeAveragingProcess: PRESSURE will be averaged." << std::endl;
            if (TIME_AVERAGED_PRESSURE.Key() == 0){
                rModelPart.AddNodalSolutionStepVariable(TIME_AVERAGED_PRESSURE);
                std::cout << "TimeAveragingProcess: TIME_AVERAGED_PRESSURE added." << std::endl;
            }
        }

        if (mDoAveragedReaction == 1){
            std::cout << "TimeAveragingProcess: REACTION will be averaged." << std::endl;
            if (TIME_AVERAGED_REACTION.Key() == 0){
                rModelPart.AddNodalSolutionStepVariable(TIME_AVERAGED_REACTION);
                std::cout << "TimeAveragingProcess: TIME_AVERAGED_REACTION added." << std::endl;
            }
        }

        std::cout << "TimeAveragingProcess: Initialized." << std::endl;
}

void TimeAveragingProcess::ExecuteInitialize() {

    int NumNodes = mrModelPart.NumberOfNodes();
    for (int i = 0; i < NumNodes; ++i)
    {
        auto inode = mrModelPart.NodesBegin() + i;
        if (mDoAveragedVelocity  == true )
            inode->SetValue(TIME_AVERAGED_VELOCITY,inode->FastGetSolutionStepValue(VELOCITY));
        if (mDoAveragedPressure  == true )
            inode->SetValue(TIME_AVERAGED_PRESSURE,inode->FastGetSolutionStepValue(PRESSURE));
        if (mDoAveragedReaction  == true )
            inode->SetValue(TIME_AVERAGED_REACTION,inode->FastGetSolutionStepValue(REACTION));
    }

}


void TimeAveragingProcess::ExecuteFinalizeSolutionStep() {
    int NumNodes = mrModelPart.NumberOfNodes();
    for (int i = 0; i < NumNodes; ++i)
    {
        ModelPart::NodeIterator inode = mrModelPart.NodesBegin() + i;
        if (mDoAveragedVelocity == true ){
            this->AverageVelocity(inode);
        }
        if (mDoAveragedPressure == true ){
            this->AveragePressure(inode);
        }
        if (mDoAveragedReaction == true ){
            this->AverageReaction(inode);
        }
    }
    std::cout << "TimeAveragingProcess: Nodal Quantities Averaged." << std::endl;
}

/* Protected functions ****************************************************/

void TimeAveragingProcess::AverageVelocity(ModelPart::NodeIterator& inode){
    double T =  mrModelPart.GetProcessInfo()[TIME];
    double Dt =  mrModelPart.GetProcessInfo()[DELTA_TIME];
    double TimeDtep = mrModelPart.GetProcessInfo()[STEP];

    array_1d<double,3> current_vel = inode->FastGetSolutionStepValue(VELOCITY);
    array_1d<double,3> &prev_ave_vel = inode->FastGetSolutionStepValue(TIME_AVERAGED_VELOCITY,1);
    array_1d<double,3 > &averaged_vel = inode->FastGetSolutionStepValue(TIME_AVERAGED_VELOCITY);

    if (T == 0.0){
		noalias(averaged_vel) = current_vel;
    } else {
        noalias(averaged_vel) = ( prev_ave_vel * (T - Dt) + current_vel * Dt ) / T;
    }
}

void TimeAveragingProcess::AverageReaction(ModelPart::NodeIterator& inode){

    double T =  mrModelPart.GetProcessInfo()[TIME];
    double Dt =  mrModelPart.GetProcessInfo()[DELTA_TIME];
    double TimeDtep = mrModelPart.GetProcessInfo()[STEP];
    
    array_1d<double,3> current_reaction = inode->FastGetSolutionStepValue(REACTION);
    array_1d<double,3> &prev_ave_reaction = inode->FastGetSolutionStepValue(TIME_AVERAGED_REACTION,1);
    array_1d<double,3> &averaged_reaction = inode->FastGetSolutionStepValue(TIME_AVERAGED_REACTION);

    if (T == 0.0){
		noalias(averaged_reaction) = averaged_reaction;
    } else {
        noalias(averaged_reaction) = ( prev_ave_reaction * (T - Dt) + current_reaction * Dt ) / T;
    }
}

void TimeAveragingProcess::AveragePressure(ModelPart::NodeIterator& inode){
    double T =  mrModelPart.GetProcessInfo()[TIME];
    double Dt =  mrModelPart.GetProcessInfo()[DELTA_TIME];
    double TimeDtep = mrModelPart.GetProcessInfo()[STEP];
    
    double current_pre = inode->FastGetSolutionStepValue(PRESSURE);
    double &prev_ave_pre = inode->FastGetSolutionStepValue(TIME_AVERAGED_PRESSURE,1);
    double &averaged_pre = inode->FastGetSolutionStepValue(TIME_AVERAGED_PRESSURE);

    if (T == 0.0){
		averaged_pre = current_pre;
    }
    else {
        averaged_pre = ( prev_ave_pre * (T - Dt) + current_pre * Dt ) / T;
    }
    inode->SetValue(TIME_AVERAGED_PRESSURE,averaged_pre);
}


/* Private functions ****************************************************/

};  // namespace Kratos.
