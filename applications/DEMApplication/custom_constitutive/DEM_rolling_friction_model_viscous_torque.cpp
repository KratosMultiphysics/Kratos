/////////////////////////////////////////////////
// Main author: Rafael Rangel (CIMNE)
// Email: rrangel@cimne.upc.edu
// Date: Oct 2022
//
// Description:
// This model can be found as Model Type B in [Jun Ai, 2011, Assessment of rolling resistance models in discrete element simulations]
/////////////////////////////////////////////////

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

  void DEMRollingFrictionModelViscousTorque::ComputeRollingFriction(SphericParticle* p_element, SphericParticle* p_neighbor, double LocalCoordSystem_2[3], double LocalContactForce[3], array_1d<double, 3>& mContactMoment)
  {
    // Get rolling friction coefficient
    Properties& properties_of_this_contact = p_element->GetProperties().GetSubProperties(p_neighbor->GetProperties().Id());
    const double rolling_friction_coefficient = properties_of_this_contact[ROLLING_FRICTION];

    // Calculate equivalent radius (this only works for sphere particles)
    const double my_radius    = p_element->GetRadius();
    const double other_radius = p_neighbor->GetRadius();
    const double equiv_radius = my_radius * other_radius / (my_radius + other_radius);

    // Get normal contact force
    const double force = fabs(LocalContactForce[2]);

    // Calculate relative translational velocity at the contact between two particles due to relative rotation
    const array_1d<double, 3>& my_angular_vel    = p_element->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
    const array_1d<double, 3>& other_angular_vel = p_neighbor->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);

    array_1d<double, 3> my_arm_vector;
    my_arm_vector[0] = -LocalCoordSystem_2[0] * my_radius;
    my_arm_vector[1] = -LocalCoordSystem_2[1] * my_radius;
    my_arm_vector[2] = -LocalCoordSystem_2[2] * my_radius;

    array_1d<double, 3> other_arm_vector;
    other_arm_vector[0] = LocalCoordSystem_2[0] * other_radius;
    other_arm_vector[1] = LocalCoordSystem_2[1] * other_radius;
    other_arm_vector[2] = LocalCoordSystem_2[2] * other_radius;

    array_1d<double, 3> my_cross_product;
    array_1d<double, 3> other_cross_product;
    GeometryFunctions::CrossProduct(my_angular_vel, my_arm_vector, my_cross_product);
    GeometryFunctions::CrossProduct(other_angular_vel, other_arm_vector, other_cross_product);

    array_1d<double, 3> relative_velocity;
    noalias(relative_velocity) = other_cross_product - my_cross_product;

    const double relative_velocity_modulus = DEM_MODULUS_3(relative_velocity);

    // Calculate rolling friction moment
    mContactMoment[0] -= rolling_friction_coefficient * equiv_radius * force * relative_velocity_modulus;
    mContactMoment[1] -= rolling_friction_coefficient * equiv_radius * force * relative_velocity_modulus;
    mContactMoment[2] -= rolling_friction_coefficient * equiv_radius * force * relative_velocity_modulus;
  }

  void DEMRollingFrictionModelViscousTorque::ComputeRollingFrictionWithWall(SphericParticle* p_element, Condition* const wall, double LocalCoordSystem_2[3], double LocalContactForce[3], array_1d<double, 3>& mContactMoment)
  {
    // Get rolling friction coefficient
    Properties& properties_of_this_contact = p_element->GetProperties().GetSubProperties(wall->GetProperties().Id());
    const double rolling_friction_coefficient = properties_of_this_contact[ROLLING_FRICTION];

    // Get particle radius
    const double radius = p_element->GetRadius();

    // Get normal contact force
    const double force = fabs(LocalContactForce[2]);

    // Calculate relative translational velocity at the contact due to relative rotation
    const array_1d<double, 3>& angular_vel = p_element->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);

    array_1d<double, 3> arm_vector;
    arm_vector[0] = -LocalCoordSystem_2[0] * radius;
    arm_vector[1] = -LocalCoordSystem_2[1] * radius;
    arm_vector[2] = -LocalCoordSystem_2[2] * radius;

    array_1d<double, 3> cross_product;
    GeometryFunctions::CrossProduct(angular_vel, arm_vector, cross_product);
    const double relative_velocity_modulus = DEM_MODULUS_3(cross_product);

    // Calculate rolling friction moment
    mContactMoment[0] -= rolling_friction_coefficient * radius * force * relative_velocity_modulus;
    mContactMoment[1] -= rolling_friction_coefficient * radius * force * relative_velocity_modulus;
    mContactMoment[2] -= rolling_friction_coefficient * radius * force * relative_velocity_modulus;
  }

} // namespace Kratos