//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Chengshun Shang (cshang@cimne.upc.edu)
//

// Description:
// This model can be found as Model Type A in [Jun Ai, 2011, Assessment of rolling resistance models in discrete element simulations]
// ATTENTION: Current implementation only works for spherical particles!

#include "DEM_rolling_friction_model_constant_torque.h"
#include "custom_utilities/GeometryFunctions.h"

namespace Kratos{

    DEMRollingFrictionModel::Pointer DEMRollingFrictionModelConstantTorque::Clone() const{
        DEMRollingFrictionModel::Pointer p_clone(new DEMRollingFrictionModelConstantTorque(*this));
        return p_clone;
    }

    std::unique_ptr<DEMRollingFrictionModel> DEMRollingFrictionModelConstantTorque::CloneUnique() {
        return Kratos::make_unique<DEMRollingFrictionModelConstantTorque>();
    }

    void DEMRollingFrictionModelConstantTorque::Check(Properties::Pointer pProp) const {
        
        if(!pProp->Has(DEM_ROLLING_FRICTION_MODEL_NAME)) 
        {
            KRATOS_ERROR << "Variable DEM_ROLLING_FRICTION_MODEL_NAME should be presented."<<std::endl;
        }
        
    }

    void DEMRollingFrictionModelConstantTorque::ComputeRollingFriction(SphericParticle* p_element, 
                                                                        SphericParticle* p_neighbor, 
                                                                        const ProcessInfo& r_process_info, 
                                                                        double LocalContactForce[3], 
                                                                        double indentation, 
                                                                        array_1d<double, 3>& mContactMoment, 
                                                                        double LocalCoordSystem2[3],
                                                                        array_1d<double, 3>& OldRollingFrictionMoment)
    {
        array_1d<double, 3> elementRelAngularVelocity;
        noalias(elementRelAngularVelocity) = p_element->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY) - p_neighbor->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
        if (elementRelAngularVelocity[0] || elementRelAngularVelocity[1] || elementRelAngularVelocity[2]){
          
            // Normalize relative angular velocity
            array_1d<double, 3> elementRelAngularVelocity_normalise = elementRelAngularVelocity;
            GeometryFunctions::normalize(elementRelAngularVelocity_normalise);

            // Get rolling friction coefficient
            Properties& r_properties = p_element->GetProperties().GetSubProperties(p_neighbor->GetProperties().Id());
            const double rolling_friction_coefficient = r_properties[ROLLING_FRICTION];

            // Get normal contact force
            const double force = std::abs(LocalContactForce[2]);

            //Calculate arm length (particle radius discounted by identation)
            //const double my_young    = p_element->GetYoung();
            //const double other_young = p_neighbor->GetYoung();
            //const double arm_length  = p_element->GetRadius() - indentation * other_young / (other_young + my_young);
            const double min_radius  = std::min(p_element->GetRadius(), p_neighbor->GetRadius());
            const double arm_length  = min_radius;

            // Calculate rolling friction moment
            array_1d<double, 3> rolling_friction_moment;
            rolling_friction_moment[0] = -elementRelAngularVelocity_normalise[0] * rolling_friction_coefficient * force * arm_length;
            rolling_friction_moment[1] = -elementRelAngularVelocity_normalise[1] * rolling_friction_coefficient * force * arm_length;
            rolling_friction_moment[2] = -elementRelAngularVelocity_normalise[2] * rolling_friction_coefficient * force * arm_length;

            // Discount rolling friction moment from contact moment
            mContactMoment[0] += rolling_friction_moment[0];
            mContactMoment[1] += rolling_friction_moment[1];
            mContactMoment[2] += rolling_friction_moment[2];

            // Compute energy dissipation from rolling friction
            double& inelastic_rollingresistance_energy = p_element->GetInelasticRollingResistanceEnergy();
            CalculateInelasticRollingResistanceEnergy(inelastic_rollingresistance_energy, rolling_friction_moment, elementRelAngularVelocity, r_process_info[DELTA_TIME]);
        }
    }

    void DEMRollingFrictionModelConstantTorque::ComputeRollingFrictionWithWall(SphericParticle* p_element, 
                                                                                Condition* const wall, 
                                                                                const ProcessInfo& r_process_info, 
                                                                                double LocalContactForce[3], 
                                                                                double indentation, 
                                                                                array_1d<double, 3>& mContactMoment, 
                                                                                double LocalCoordSystem2[3],
                                                                                array_1d<double, 3>& OldRollingFrictionMoment)
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

    void DEMRollingFrictionModelConstantTorque::CalculateInelasticRollingResistanceEnergy(double& inelastic_rollingresistance_energy, const array_1d<double, 3>& rolling_friction_moment, const array_1d<double, 3>& relative_angular_velocity, double dt)
    {
      // Each particle in a contact with another particle receives half the contact energy
      const double rollingresistance_energy = std::abs(DEM_INNER_PRODUCT_3(rolling_friction_moment, relative_angular_velocity)) * dt;
      inelastic_rollingresistance_energy += 0.50 * rollingresistance_energy;
    }

    void DEMRollingFrictionModelConstantTorque::CalculateInelasticRollingResistanceEnergyWithWall(double& inelastic_rollingresistance_energy, const array_1d<double, 3>& rolling_friction_moment, const array_1d<double, 3>& relative_angular_velocity, double dt)
    {
      // Each particle in a contact with a wall receives all the contact energy
      const double rollingresistance_energy = std::abs(DEM_INNER_PRODUCT_3(rolling_friction_moment, relative_angular_velocity)) * dt;
      inelastic_rollingresistance_energy += rollingresistance_energy;
    }

}//namespace Kratos