#include "includes/variables.h"
#include "includes/kratos_parameters.h"
#include "steady_state_indicator_utility.h"


namespace Kratos {
    SteadyStateIndicatorUtility::SteadyStateIndicatorUtility(ModelPart &ModelPart):
      mrModelPart(ModelPart){
        mChangeInVelocity = 0.0;
        mChangeInPressure = 0.0;
        mChangeInTimeAveragedVelocity = 0.0;
        mChangeInTimeAveragedPressure = 0.0;
        mIsSteady = false;
    }

    void SteadyStateIndicatorUtility::EstimateTimeAveragedQuantityChangesInTime(){
        int NumNodes = mrModelPart.NumberOfNodes();
        double DomainArea = 0.0;
        for (int i = 0; i < NumNodes; ++i)
        {
            auto inode = mrModelPart.NodesBegin() + i;
            double NodalArea = inode->FastGetSolutionStepValue(NODAL_AREA);
            mChangeInTimeAveragedVelocity += this->NodalTimeAveragedVelocityChange((*inode)) * NodalArea;
            mChangeInTimeAveragedPressure += this->NodalTimeAveragedPressureChange((*inode)) * NodalArea;
            DomainArea += NodalArea;
        }
        mChangeInTimeAveragedVelocity = mChangeInTimeAveragedVelocity / DomainArea;
        mChangeInTimeAveragedPressure = mChangeInTimeAveragedPressure / DomainArea;
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
        mChangeInVelocity = mChangeInVelocity / DomainArea;
        mChangeInPressure = mChangeInPressure / DomainArea;
    }
    
    double SteadyStateIndicatorUtility::NodalVelocityChange(ModelPart::NodeType& inode){
        array_1d<double, 3 > prev_vel = inode.FastGetSolutionStepValue(VELOCITY,1);
        array_1d<double, 3 > curr_vel = inode.FastGetSolutionStepValue(VELOCITY);
        double change_in_percentage = 0.0;
        double prev_velocity_magnitude = 0.0;
        double curr_velocity_magnitude = 0.0;

        for (int i = 0; i<3; ++i){
            prev_velocity_magnitude += sqrt(prev_vel[i]*prev_vel[i]);
            curr_velocity_magnitude += sqrt(curr_vel[i]*curr_vel[i]);
        }

        if (prev_velocity_magnitude == 0.0)
            change_in_percentage = 0.0;
        else
            change_in_percentage = sqrt((curr_velocity_magnitude - prev_velocity_magnitude)*(curr_velocity_magnitude - prev_velocity_magnitude))/fabs(prev_velocity_magnitude);
        
        return change_in_percentage;
    }

    double SteadyStateIndicatorUtility::NodalTimeAveragedVelocityChange(ModelPart::NodeType& inode){
        array_1d<double, 3 > prev_ave_vel = inode.FastGetSolutionStepValue(TIME_AVERAGED_VELOCITY,1);
        array_1d<double, 3 > curr_ave_vel = inode.FastGetSolutionStepValue(TIME_AVERAGED_VELOCITY);
        double change_in_percentage = 0.0;
        double prev_velocity_magnitude = 0.0;
        double curr_velocity_magnitude = 0.0;

        for (int i = 0; i<3; ++i){
            prev_velocity_magnitude += sqrt(prev_ave_vel[i]*prev_ave_vel[i]);
            curr_velocity_magnitude += sqrt(curr_ave_vel[i]*curr_ave_vel[i]);
        }

        if (prev_velocity_magnitude == 0.0)
            change_in_percentage = 0.0;
        else
            change_in_percentage = sqrt((curr_velocity_magnitude - prev_velocity_magnitude)*(curr_velocity_magnitude - prev_velocity_magnitude))/fabs(prev_velocity_magnitude);
        
        return change_in_percentage;
    }

    double SteadyStateIndicatorUtility::NodalPressureChange(ModelPart::NodeType& inode){
        double prev_p = inode.FastGetSolutionStepValue(PRESSURE,1);
        double curr_p = inode.FastGetSolutionStepValue(PRESSURE);
        double change_in_percentage = 0.0;
        if (prev_p == 0.0)
            change_in_percentage = 0.0;
        else
            change_in_percentage = sqrt((curr_p - prev_p)*(curr_p - prev_p))/fabs(prev_p);
        return change_in_percentage;
    }

    double SteadyStateIndicatorUtility::NodalTimeAveragedPressureChange(ModelPart::NodeType& inode){
        double prev_ave_p = inode.FastGetSolutionStepValue(TIME_AVERAGED_PRESSURE,1);
        double curr_ave_p = inode.FastGetSolutionStepValue(TIME_AVERAGED_PRESSURE);
        double change_in_percentage = 0.0;
        if (prev_ave_p == 0.0)
            change_in_percentage = 0.0;
        else
            change_in_percentage = sqrt((curr_ave_p - prev_ave_p)*(curr_ave_p - prev_ave_p))/fabs(prev_ave_p);
        return change_in_percentage;
    }

}