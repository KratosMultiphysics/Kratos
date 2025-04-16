//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Chengshun Shang (cshang@cimne.upc.edu)
//

// Description:
// This model can be found as Model Type C in [Jun Ai, 2011, Assessment of rolling resistance models in discrete element simulations]
// ATTENTION: Current implementation only works for spherical particles!

#include "DEM_rolling_friction_model_elastic_platic_resistance.h"
#include "custom_utilities/GeometryFunctions.h"

namespace Kratos{

    DEMRollingFrictionModel::Pointer DEMRollingFrictionModelElasticPlasticResistance::Clone() const{
        DEMRollingFrictionModel::Pointer p_clone(new DEMRollingFrictionModelElasticPlasticResistance(*this));
        return p_clone;
    }

    std::unique_ptr<DEMRollingFrictionModel> DEMRollingFrictionModelElasticPlasticResistance::CloneUnique() {
        return Kratos::make_unique<DEMRollingFrictionModelElasticPlasticResistance>();
    }

    void DEMRollingFrictionModelElasticPlasticResistance::Check(Properties::Pointer pProp) const {
        
        if(!pProp->Has(DEM_ROLLING_FRICTION_MODEL_NAME)) 
        {
            KRATOS_ERROR << "Variable DEM_ROLLING_FRICTION_MODEL_NAME should be presented."<<std::endl;
        }
        
    }

    void DEMRollingFrictionModelElasticPlasticResistance::InitializeContact(SphericParticle* const p_element, SphericParticle* const p_neighbor, const double indentation) {
        
        //Get equivalent Radius
        const double my_radius       = p_element->GetRadius();
        const double other_radius    = p_neighbor->GetRadius();
        const double radius_sum      = my_radius + other_radius;
        const double radius_sum_inv  = 1.0 / radius_sum;
        const double equiv_radius    = my_radius * other_radius * radius_sum_inv;

        //Get equivalent Young's Modulus
        const double my_young        = p_element->GetYoung();
        const double other_young     = p_neighbor->GetYoung();
        const double my_poisson      = p_element->GetPoisson();
        const double other_poisson   = p_neighbor->GetPoisson();
        const double equiv_young     = my_young * other_young / (other_young * (1.0 - my_poisson * my_poisson) + my_young * (1.0 - other_poisson * other_poisson));
        //Get equivalent Shear Modulus
        const double my_shear_modulus = 0.5 * my_young / (1.0 + my_poisson);
        const double other_shear_modulus = 0.5 * other_young / (1.0 + other_poisson);
        const double equiv_shear = 1.0 / ((2.0 - my_poisson)/my_shear_modulus + (2.0 - other_poisson)/other_shear_modulus);

        //Normal and Tangent elastic constants
        const double sqrt_equiv_radius_and_indentation = sqrt(equiv_radius * indentation);
        double Kn = 2.0 * equiv_young * sqrt_equiv_radius_and_indentation;
        mKt = 4.0 * equiv_shear * Kn / equiv_young;
    }

    void DEMRollingFrictionModelElasticPlasticResistance::ComputeRollingFriction(SphericParticle* p_element, SphericParticle* p_neighbor, const ProcessInfo& r_process_info, double LocalContactForce[3], double indentation, array_1d<double, 3>& mContactMoment)
    {
        array_1d<double, 3> elementRelAngularVelocity;
        noalias(elementRelAngularVelocity) = p_element->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY) - p_neighbor->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
        if (elementRelAngularVelocity[0] || elementRelAngularVelocity[1] || elementRelAngularVelocity[2]){

            // Normalize relative angular velocity
            array_1d<double, 3> elementRelAngularVelocity_normalise = elementRelAngularVelocity;
            GeometryFunctions::normalize(elementRelAngularVelocity_normalise);

            const array_1d<double, 3>& my_delta_rotation = p_element->GetGeometry()[0].FastGetSolutionStepValue(DELTA_ROTATION);
            const array_1d<double, 3>& other_delta_rotation = p_neighbor->GetGeometry()[0].FastGetSolutionStepValue(DELTA_ROTATION);
            array_1d<double, 3> delta_theta = ZeroVector(3);
            noalias(delta_theta) = my_delta_rotation - other_delta_rotation;

            // Calculate delta_theta for bending
            array_1d<double, 3> delta_theta_t;
            array_1d<double, 3> delta_theta_n;
            array_1d<double, 3> normal_contact_vector;
            noalias(normal_contact_vector) = p_element->GetGeometry()[0].Coordinates() - p_neighbor->GetGeometry()[0].Coordinates();
            GeometryFunctions::normalize(normal_contact_vector);
            GeometryFunctions::CrossProduct(delta_theta, normal_contact_vector, delta_theta_n);
            noalias(delta_theta_t) = delta_theta - delta_theta_n;

            // Get rolling friction coefficient
            Properties& r_properties = p_element->GetProperties().GetSubProperties(p_neighbor->GetProperties().Id());
            const double rolling_friction_coefficient = r_properties[ROLLING_FRICTION];

            const double equivalent_radius = p_element->GetRadius() * p_neighbor->GetRadius() / (p_element->GetRadius() + p_neighbor->GetRadius());
            const double arm_length  = equivalent_radius;

            // Get normal contact force
            const double force = std::abs(LocalContactForce[2]);
            double max_rolling_friction_moment = rolling_friction_coefficient * force * arm_length;

            const double Kr = mKt * equivalent_radius * equivalent_radius; //TODO: This must be improved
            m_rolling_friction_moment[0] -= Kr * delta_theta[0];
            m_rolling_friction_moment[1] -= Kr * delta_theta[1];
            m_rolling_friction_moment[2] -= Kr * delta_theta[2];

            // Check if the rolling friction moment exceeds the maximum value
            double m_rolling_friction_moment_modlule = GeometryFunctions::module(m_rolling_friction_moment);

            if (m_rolling_friction_moment_modlule > max_rolling_friction_moment) {
                // Normalize the rolling friction moment to the maximum value
                m_rolling_friction_moment[0] *= max_rolling_friction_moment / m_rolling_friction_moment_modlule;
                m_rolling_friction_moment[1] *= max_rolling_friction_moment / m_rolling_friction_moment_modlule;
                m_rolling_friction_moment[2] *= max_rolling_friction_moment / m_rolling_friction_moment_modlule;
            }

            // Discount rolling friction moment from contact moment
            mContactMoment[0] += m_rolling_friction_moment[0];
            mContactMoment[1] += m_rolling_friction_moment[1];
            mContactMoment[2] += m_rolling_friction_moment[2];

            // Compute energy dissipation from rolling friction
            double& inelastic_rollingresistance_energy = p_element->GetInelasticRollingResistanceEnergy();
            array_1d<double, 3> rolling_friction_moment = ZeroVector(3);
            rolling_friction_moment[0] = m_rolling_friction_moment[0];
            rolling_friction_moment[1] = m_rolling_friction_moment[1];
            rolling_friction_moment[2] = m_rolling_friction_moment[2];
            CalculateInelasticRollingResistanceEnergy(inelastic_rollingresistance_energy, rolling_friction_moment, elementRelAngularVelocity, r_process_info[DELTA_TIME]);
        }
    }

    void DEMRollingFrictionModelElasticPlasticResistance::ComputeRollingFrictionWithWall(SphericParticle* p_element, Condition* const wall, const ProcessInfo& r_process_info, double LocalContactForce[3], double indentation, array_1d<double, 3>& mContactMoment)
    {
        array_1d<double, 3> element1AngularVelocity;
        noalias(element1AngularVelocity) = p_element->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
        if (element1AngularVelocity[0] || element1AngularVelocity[1] || element1AngularVelocity[2]){
            
            // Normalize angular velocity
            array_1d<double, 3> element1AngularVelocity_normalise = element1AngularVelocity;
            GeometryFunctions::normalize(element1AngularVelocity_normalise);

            // Get rolling friction coefficient
            Properties& r_properties = p_element->GetProperties().GetSubProperties(wall->GetProperties().Id());
            const double rolling_friction_coefficient = r_properties[ROLLING_FRICTION];

            // Get normal contact force
            const double force = std::abs(LocalContactForce[2]);

            // Calculate arm length (particle radius discounted by identation)
            const double arm_length = p_element->GetRadius() - indentation;

            // Calculate rolling friction moment
            array_1d<double, 3> rolling_friction_moment;
            rolling_friction_moment[0] = -element1AngularVelocity_normalise[0] * rolling_friction_coefficient * force * arm_length;
            rolling_friction_moment[1] = -element1AngularVelocity_normalise[1] * rolling_friction_coefficient * force * arm_length;
            rolling_friction_moment[2] = -element1AngularVelocity_normalise[2] * rolling_friction_coefficient * force * arm_length;

            // Discount rolling friction moment from contact moment
            mContactMoment[0] += rolling_friction_moment[0];
            mContactMoment[1] += rolling_friction_moment[1];
            mContactMoment[2] += rolling_friction_moment[2];

            // Compute energy dissipation from rolling friction
            double& inelastic_rollingresistance_energy = p_element->GetInelasticRollingResistanceEnergy();
            CalculateInelasticRollingResistanceEnergyWithWall(inelastic_rollingresistance_energy, rolling_friction_moment, element1AngularVelocity, r_process_info[DELTA_TIME]);
        }
    }

    void DEMRollingFrictionModelElasticPlasticResistance::CalculateInelasticRollingResistanceEnergy(double& inelastic_rollingresistance_energy, const array_1d<double, 3>& rolling_friction_moment, const array_1d<double, 3>& relative_angular_velocity, double dt)
    {
      // Each particle in a contact with another particle receives half the contact energy
      const double rollingresistance_energy = std::abs(DEM_INNER_PRODUCT_3(rolling_friction_moment, relative_angular_velocity)) * dt;
      inelastic_rollingresistance_energy += 0.50 * rollingresistance_energy;
    }

    void DEMRollingFrictionModelElasticPlasticResistance::CalculateInelasticRollingResistanceEnergyWithWall(double& inelastic_rollingresistance_energy, const array_1d<double, 3>& rolling_friction_moment, const array_1d<double, 3>& relative_angular_velocity, double dt)
    {
      // Each particle in a contact with a wall receives all the contact energy
      const double rollingresistance_energy = std::abs(DEM_INNER_PRODUCT_3(rolling_friction_moment, relative_angular_velocity)) * dt;
      inelastic_rollingresistance_energy += rollingresistance_energy;
    }

}//namespace Kratos