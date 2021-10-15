// Author: Joaquin Gonzalez-Usua (jgonzalez@cimne.upc.edu)
// Date: April 2021
//This drag law is explained and implemented in detail in the following paper:
//Lattice Boltzmann investigation on fluid flows through packed beds: Interaction between fluid rheology and bed properties
//Qi et al. 2020
#include "swimming_DEM_application.h"
#include "dallavalle_drag_law.h"

namespace Kratos {

    DragLaw::Pointer DallavalleDragLaw::Clone() const {
        DallavalleDragLaw::Pointer p_clone(new DallavalleDragLaw(*this));
        return p_clone;
    }

    void DallavalleDragLaw::Initialize(const ProcessInfo& r_process_info) {}

    std::string DallavalleDragLaw::GetTypeOfLaw() {
        std::string type_of_law = "Dallavalle drag law";
        return type_of_law;
    }

    void DallavalleDragLaw::ComputeForce(SphericParticle* p_particle,
                                       const double reynolds_number,
                                       double particle_radius,
                                       double fluid_density,
                                       double fluid_kinematic_viscosity,
                                       array_1d<double, 3>& minus_slip_velocity,
                                       array_1d<double, 3>& drag_force,
                                       const ProcessInfo& r_current_process_info)
    {

        double drag_coeff;

        const double equivalent_diameter = this->CalculateEquivalentDiameter(p_particle);

        Geometry<Node<3> >& r_geometry = p_particle->GetGeometry();
        const double eps = r_geometry[0].FastGetSolutionStepValue(FLUID_FRACTION_PROJECTED);

        double weighting_sum = this->CalculateWeightingSum(p_particle, equivalent_diameter);

        if (reynolds_number < 0.01){
            return StokesDragLaw::ComputeForce(p_particle,
                                               reynolds_number,
                                               particle_radius,
                                               fluid_density,
                                               fluid_kinematic_viscosity,
                                               minus_slip_velocity,
                                               drag_force,
                                               r_current_process_info);
        }

        double y = (2 *  particle_radius) / equivalent_diameter;
        const double weighting_parameter = 0.5 * (eps) / weighting_sum + 0.5 * y + 0.5 * (1 - eps) * std::pow(y, 2);
        const double norm_minus_slip_velocity = SWIMMING_MODULUS_3(minus_slip_velocity);
        const double mean_reynolds_particle = eps * norm_minus_slip_velocity * equivalent_diameter / fluid_kinematic_viscosity;
        array_1d<double, 3>& slip_vel =  r_geometry[0].FastGetSolutionStepValue(SLIP_VELOCITY);
        slip_vel = minus_slip_velocity;
        const double beta = 2.65 * (eps + 1) - (5.3 - 3.5 * eps) * std::pow(eps, 2) * std::exp(-std::pow(1.5 - std::log(mean_reynolds_particle), 2) / 2);

        drag_coeff = std::pow((2.654 / (1.0 + 3.213) + 4.8 / (std::sqrt(mean_reynolds_particle))),2);
        noalias(drag_force) = 1.0 / 8.0 * drag_coeff * Globals::Pi * fluid_density * y * weighting_parameter * std::pow(equivalent_diameter, 2) * norm_minus_slip_velocity * minus_slip_velocity * std::pow(eps, 2.0 - beta);
    }

    double DallavalleDragLaw::CalculateEquivalentDiameter(SphericParticle* p_particle){

        const double particle_radius = p_particle->GetRadius();
        double particle_equivalent_diameter = this->GetParticleMassFraction(p_particle) / (2 * particle_radius);
        double particle_diameter = 2 * particle_radius;
        std::vector<double> vector_diameter;
        vector_diameter.push_back(particle_diameter);

        for (unsigned int i = 0; i < p_particle->mNeighbourElements.size(); ++i){
            const double neigh_diameter = 2 *  p_particle->mNeighbourElements[i]->GetRadius();
            std::vector<double>::iterator it = std::find(vector_diameter.begin(), vector_diameter.end(), neigh_diameter);
            if(it == vector_diameter.end()){
                vector_diameter.push_back(neigh_diameter);
                particle_equivalent_diameter += this->GetParticleMassFraction(p_particle->mNeighbourElements[i])
                                            / (neigh_diameter);

            }
        }

        return 1.0 / particle_equivalent_diameter;
    }

    double DallavalleDragLaw::GetParticleMassFraction(SphericParticle* p_particle){

        double particle_i_component_volume = p_particle->CalculateVolume();
        double sum_volume = particle_i_component_volume;
        for (unsigned int i = 0; i < p_particle->mNeighbourElements.size(); ++i){
            if (p_particle->mNeighbourElements[i]->GetRadius() == p_particle->GetRadius()){
                particle_i_component_volume += p_particle->mNeighbourElements[i]->CalculateVolume();
            }
            sum_volume += p_particle->mNeighbourElements[i]->CalculateVolume();
        }
        return particle_i_component_volume / sum_volume;
    }

    double DallavalleDragLaw::CalculateWeightingSum(SphericParticle* p_particle, const double& equivalent_diameter){

        const double mass_fraction = this->GetParticleMassFraction(p_particle);
        const double particle_diameter = 2 *  p_particle->GetRadius();

        std::vector<double> vector_diameter;
        vector_diameter.push_back(particle_diameter);

        double sum_parameter = mass_fraction * std::pow(equivalent_diameter / particle_diameter, 2);

        for (unsigned int i = 0; i < p_particle->mNeighbourElements.size(); ++i){
            const double neigh_diameter = 2 *  p_particle->mNeighbourElements[i]->GetRadius();
            std::vector<double>::iterator it = std::find(vector_diameter.begin(), vector_diameter.end(), neigh_diameter);
            if(it == vector_diameter.end()){
                vector_diameter.push_back(neigh_diameter);
                sum_parameter += this->GetParticleMassFraction(p_particle->mNeighbourElements[i]) * std::pow(equivalent_diameter / neigh_diameter, 2);
            }
        }

        return sum_parameter;
        }

} // namespace Kratos