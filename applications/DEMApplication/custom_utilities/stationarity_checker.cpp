// Author: Salva Latorre

#include "custom_elements/spheric_continuum_particle.h"
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

    bool StationarityChecker::CheckIfVariableIsNullInModelPart(ModelPart& r_spheres_modelPart, const Variable<double>& r_var,
                                                               const double tolerance, const bool ignore_isolated_particles) {

        KRATOS_ERROR_IF(!r_spheres_modelPart.HasNodalSolutionStepVariable(r_var)) << "Variable " << r_var.Name() << " is not added to the nodes of the ModelPart. Steadiness cannot be assessed with this variable" << std::endl;

        typedef ModelPart::ElementsContainerType ElementsArrayType;
        ElementsArrayType& pElements = r_spheres_modelPart.Elements();

        for (int k = 0; k < (int)pElements.size(); k++) {

            ElementsArrayType::iterator it = pElements.ptr_begin() + k;
            Element* p_element = &(*it);
            SphericContinuumParticle* p_sphere = dynamic_cast<SphericContinuumParticle*>(p_element);
            const double value = p_sphere->GetGeometry()[0].FastGetSolutionStepValue(r_var);

            const unsigned int initial_cont_neigh_size = p_sphere->mContinuumInitialNeighborsSize;
            if (ignore_isolated_particles) {
                if (!initial_cont_neigh_size) {
                    continue;
                } else {
                    unsigned int damaged_bonds = 0;
                    for (int l = 0; l < (int)initial_cont_neigh_size; l++) {
                        if (p_sphere->mIniNeighbourFailureId[l]) {
                            damaged_bonds++;
                        }
                    }
                    if (damaged_bonds == initial_cont_neigh_size) {
                        continue;
                    }
                }
            } 
            
            if (std::abs(value) > tolerance) {
                return false;
            }
        }
        return true;
    }

    std::string StationarityChecker::Info() const {
        std::stringstream buffer;
        buffer << "StationarityChecker" ;
        return buffer.str();
    }

    void StationarityChecker::PrintInfo(std::ostream& rOStream) const {rOStream << "StationarityChecker";}

    void StationarityChecker::PrintData(std::ostream& rOStream) const {}

} // namespace Kratos
