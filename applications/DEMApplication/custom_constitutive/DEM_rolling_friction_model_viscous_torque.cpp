//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

// Description:
// This model can be found as Model Type B in [Jun Ai, 2011, Assessment of rolling resistance models in discrete element simulations]
// ATTENTION: Current implementation only works for spherical particles!

// Project includes
#include "DEM_rolling_friction_model_viscous_torque.h"

namespace Kratos {

  DEMRollingFrictionModel::Pointer DEMRollingFrictionModelViscousTorque::Clone() const {
    DEMRollingFrictionModel::Pointer p_clone(new DEMRollingFrictionModelViscousTorque(*this));
    return p_clone;
  }

  std::unique_ptr<DEMRollingFrictionModel> DEMRollingFrictionModelViscousTorque::CloneUnique() {
    return Kratos::make_unique<DEMRollingFrictionModelViscousTorque>();
  }

  void DEMRollingFrictionModelViscousTorque::Check(Properties::Pointer pProp) const {
    if (!pProp->Has(DEM_ROLLING_FRICTION_MODEL_NAME)) {
      KRATOS_ERROR << "Variable DEM_ROLLING_FRICTION_MODEL_NAME should be presented." << std::endl;
    }
  }

  void DEMRollingFrictionModelViscousTorque::ComputeRollingFriction(SphericParticle* p_element, SphericParticle* p_neighbor, const ProcessInfo& r_process_info, double LocalContactForce[3], double indentation, array_1d<double, 3>& mContactMoment)
  {
    // Get rolling friction coefficient
    Properties& r_properties = p_element->GetProperties().GetSubProperties(p_neighbor->GetProperties().Id());
    const double rolling_friction_coefficient = r_properties[ROLLING_FRICTION];

    // Get normal contact force
    const double force = std::abs(LocalContactForce[2]);

    // Calculate arm length (particle radius discounted by identation)
    const double my_young         = p_element->GetYoung();
    const double other_young      = p_neighbor->GetYoung();
    const double identation_ratio = other_young / (other_young + my_young);
    const double my_arm_length    = p_element->GetRadius()  - indentation * identation_ratio;
    const double other_arm_length = p_neighbor->GetRadius() - indentation * (1.0 - identation_ratio);

    // Calculate relative angular velocity
    array_1d<double, 3> elementRelAngularVelocity;
    const array_1d<double, 3>& my_angular_vel    = p_element->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
    const array_1d<double, 3>& other_angular_vel = p_neighbor->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
    noalias(elementRelAngularVelocity) = my_angular_vel - other_angular_vel;

    // Calculate relative translational velocity at the contact between two particles due to relative rotation
    array_1d<double, 3> relative_translational_velocity;
    relative_translational_velocity[0] = my_arm_length * my_angular_vel[0] - other_arm_length * other_angular_vel[0];
    relative_translational_velocity[1] = my_arm_length * my_angular_vel[1] - other_arm_length * other_angular_vel[1];
    relative_translational_velocity[2] = my_arm_length * my_angular_vel[2] - other_arm_length * other_angular_vel[2];

    // Calculate rolling friction moment
    array_1d<double, 3> rolling_friction_moment;
    rolling_friction_moment[0] = -rolling_friction_coefficient * force * my_arm_length * relative_translational_velocity[0];
    rolling_friction_moment[1] = -rolling_friction_coefficient * force * my_arm_length * relative_translational_velocity[1];
    rolling_friction_moment[2] = -rolling_friction_coefficient * force * my_arm_length * relative_translational_velocity[2];

    // Discount rolling friction moment from contact moment
    mContactMoment[0] += rolling_friction_moment[0];
    mContactMoment[1] += rolling_friction_moment[1];
    mContactMoment[2] += rolling_friction_moment[2];

    // Compute energy dissipation from rolling friction
    double& inelastic_rollingresistance_energy = p_element->GetInelasticRollingResistanceEnergy();
    CalculateInelasticRollingResistanceEnergy(inelastic_rollingresistance_energy, rolling_friction_moment, elementRelAngularVelocity, r_process_info[DELTA_TIME]);
  }

  void DEMRollingFrictionModelViscousTorque::ComputeRollingFrictionWithWall(SphericParticle* p_element, Condition* const wall, const ProcessInfo& r_process_info, double LocalContactForce[3], double indentation, array_1d<double, 3>& mContactMoment)
  {
    // ATTENTION: It considers only the motion of the particle!

    // Get rolling friction coefficient
    Properties& r_properties = p_element->GetProperties().GetSubProperties(wall->GetProperties().Id());
    const double rolling_friction_coefficient = r_properties[ROLLING_FRICTION];

    // Get normal contact force
    const double force = std::abs(LocalContactForce[2]);

    // Calculate arm length (particle radius discounted by identation)
    const double arm_length = p_element->GetRadius() - indentation;

    // Calculate relative translational velocity at the contact due to relative rotation
    const array_1d<double, 3>& angular_vel = p_element->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
    array_1d<double, 3> relative_translational_velocity;
    relative_translational_velocity[0] = arm_length * angular_vel[0];
    relative_translational_velocity[1] = arm_length * angular_vel[1];
    relative_translational_velocity[2] = arm_length * angular_vel[2];

    // Calculate rolling friction moment
    array_1d<double, 3> rolling_friction_moment;
    rolling_friction_moment[0] = -rolling_friction_coefficient * force * arm_length * relative_translational_velocity[0];
    rolling_friction_moment[1] = -rolling_friction_coefficient * force * arm_length * relative_translational_velocity[1];
    rolling_friction_moment[2] = -rolling_friction_coefficient * force * arm_length * relative_translational_velocity[2];

    // Discount rolling friction moment from contact moment
    mContactMoment[0] += rolling_friction_moment[0];
    mContactMoment[1] += rolling_friction_moment[1];
    mContactMoment[2] += rolling_friction_moment[2];

    // Compute energy dissipation from rolling friction
    double& inelastic_rollingresistance_energy = p_element->GetInelasticRollingResistanceEnergy();
    CalculateInelasticRollingResistanceEnergyWithWall(inelastic_rollingresistance_energy, rolling_friction_moment, angular_vel, r_process_info[DELTA_TIME]);
  }

  void DEMRollingFrictionModelViscousTorque::CalculateInelasticRollingResistanceEnergy(double& inelastic_rollingresistance_energy, const array_1d<double, 3>& rolling_friction_moment, const array_1d<double, 3>& relative_angular_velocity, double dt)
  {
    // Each particle in a contact with another particle receives half the contact energy
    const double rollingresistance_energy = std::abs(DEM_INNER_PRODUCT_3(rolling_friction_moment, relative_angular_velocity)) * dt;
    inelastic_rollingresistance_energy += 0.50 * rollingresistance_energy;
  }

  void DEMRollingFrictionModelViscousTorque::CalculateInelasticRollingResistanceEnergyWithWall(double& inelastic_rollingresistance_energy, const array_1d<double, 3>& rolling_friction_moment, const array_1d<double, 3>& relative_angular_velocity, double dt)
  {
    // Each particle in a contact with a wall receives all the contact energy
    const double rollingresistance_energy = std::abs(DEM_INNER_PRODUCT_3(rolling_friction_moment, relative_angular_velocity)) * dt;
    inelastic_rollingresistance_energy += rollingresistance_energy;
  }

} // namespace Kratos