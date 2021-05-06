// Author: Joaquin Gonzalez-Usua (jgonzalez@cimne.upc.edu)
// Date: April 2021
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
        double eps = r_geometry[0].FastGetSolutionStepValue(FLUID_FRACTION_PROJECTED);

        if (eps > 0.999){
            eps = 0.9;
        }

        double weighting_sum = this->CalculateWeightingSum(p_particle);

        const double y = (2 * particle_radius) / equivalent_diameter;
        const double weighting_parameter = 0.5 * (1 + eps) / weighting_sum + 0.5 * y + 0.5 * (1 - eps) * std::pow(y, 2);
        const double norm_minus_slip_velocity = SWIMMING_MODULUS_3(minus_slip_velocity);
        const double mean_reynolds_particle = eps * norm_minus_slip_velocity * equivalent_diameter / fluid_kinematic_viscosity;

        const double beta = 2.65 * (eps + 1) - (5.3 - 3.5 * eps) * std::pow(eps, 2) * std::exp(-std::pow(1.5 - std::log(mean_reynolds_particle), 2) / 2);

        drag_coeff = std::pow((2.654 / (1.0 + 3.213) + 4.8 / (std::sqrt(reynolds_number))),2);

        noalias(drag_force) = (1.0 / 8.0) * drag_coeff * Globals::Pi * fluid_density * weighting_parameter * std::pow(equivalent_diameter, 2) * norm_minus_slip_velocity * minus_slip_velocity * std::pow(eps, 2.0 - beta);
    }

    double DallavalleDragLaw::CalculateEquivalentDiameter(SphericParticle* p_particle){

        double particle_equivalent_diameter = this->GetParticleMassFraction(p_particle) / (2 * p_particle->GetRadius());

        for (unsigned int i = 0; i < p_particle->mNeighbourElements.size(); ++i){
            particle_equivalent_diameter += this->GetParticleMassFraction(p_particle->mNeighbourElements[i])
                                            / (2 * p_particle->mNeighbourElements[i]->GetRadius());
        }

        return 1.0 / particle_equivalent_diameter;
    }

    double DallavalleDragLaw::GetParticleMassFraction(SphericParticle* p_particle){

        double sum_masses = p_particle->GetMass();
        for (unsigned int i = 0; i < p_particle->mNeighbourElements.size(); ++i){
            sum_masses += p_particle->mNeighbourElements[i]->GetMass();
        }
        return p_particle->GetMass()/sum_masses;
    }

    double DallavalleDragLaw::CalculateWeightingSum(SphericParticle* p_particle){

    double sum_parameter = this->GetParticleMassFraction(p_particle) / std::pow( 2 * p_particle->GetRadius() / this->CalculateEquivalentDiameter(p_particle), 2);
    for (unsigned int i = 0; i < p_particle->mNeighbourElements.size(); ++i){
        sum_parameter += this->GetParticleMassFraction(p_particle->mNeighbourElements[i])
                        / std::pow(2 * p_particle->mNeighbourElements[i]->GetRadius() / this->CalculateEquivalentDiameter(p_particle->mNeighbourElements[i]), 2);
    }
    return sum_parameter;
    }

} // namespace Kratos