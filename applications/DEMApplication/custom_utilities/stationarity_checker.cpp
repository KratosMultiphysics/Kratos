// Author: Salva Latorre

#include "stationarity_checker.h"

namespace Kratos {

    StationarityChecker::StationarityChecker() {
        mPreviousChangeTime = 0.0;
    }
    
    StationarityChecker::~StationarityChecker() {}
    
    bool StationarityChecker::CheckIfItsTimeToChangeGravity(ModelPart& rSpheresModelPart,
                                       const double velocity_threshold_for_gravity_change,
                                       const double min_time_between_changes,
                                       const double max_time_between_changes) {

        const double& current_time = rSpheresModelPart.GetProcessInfo()[TIME];
                        
        if (current_time < mPreviousChangeTime + min_time_between_changes) return false;
        if (current_time > mPreviousChangeTime + max_time_between_changes) {
            mPreviousChangeTime  = current_time;
            return true;
        }

        const size_t number_of_nodes = rSpheresModelPart.Nodes().size();
        double max_squared_velocity = 0.0;
        for (size_t i = 0; i < number_of_nodes; i++) {
            
            const auto node_it = rSpheresModelPart.Nodes().begin() + i;
            auto& vel = node_it->FastGetSolutionStepValue(VELOCITY);
            const double node_i_squared_velocity_module = vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2];
            if (node_i_squared_velocity_module > max_squared_velocity) max_squared_velocity = node_i_squared_velocity_module;
        }

        if (max_squared_velocity < velocity_threshold_for_gravity_change * velocity_threshold_for_gravity_change) {
            mPreviousChangeTime  = current_time;
            return true;
        } else return false;
    }

    std::string StationarityChecker::Info() const {
        std::stringstream buffer;
        buffer << "StationarityChecker" ;
        return buffer.str();
    }

    void StationarityChecker::PrintInfo(std::ostream& rOStream) const {rOStream << "StationarityChecker";}

    void StationarityChecker::PrintData(std::ostream& rOStream) const {}
    
} // namespace Kratos
