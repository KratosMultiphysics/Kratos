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

  void DEMRollingFrictionModelViscousTorque::ComputeRollingFriction(SphericParticle* p_element, SphericParticle* p_neighbor, double LocalCoordSystem_2[3], double LocalContactForce[3], double indentation, array_1d<double, 3>& mContactMoment)
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

    // Calculate relative translational velocity at the contact between two particles due to relative rotation
    array_1d<double, 3> my_arm_vector;
    my_arm_vector[0] = -LocalCoordSystem_2[0] * my_arm_length;
    my_arm_vector[1] = -LocalCoordSystem_2[1] * my_arm_length;
    my_arm_vector[2] = -LocalCoordSystem_2[2] * my_arm_length;

    array_1d<double, 3> other_arm_vector;
    other_arm_vector[0] = LocalCoordSystem_2[0] * other_arm_length;
    other_arm_vector[1] = LocalCoordSystem_2[1] * other_arm_length;
    other_arm_vector[2] = LocalCoordSystem_2[2] * other_arm_length;

    const array_1d<double, 3>& my_angular_vel    = p_element->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
    const array_1d<double, 3>& other_angular_vel = p_neighbor->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);

    array_1d<double, 3> my_cross_product;
    array_1d<double, 3> other_cross_product;
    GeometryFunctions::CrossProduct(my_angular_vel, my_arm_vector, my_cross_product);
    GeometryFunctions::CrossProduct(other_angular_vel, other_arm_vector, other_cross_product);

    array_1d<double, 3> relative_velocity;
    noalias(relative_velocity) = other_cross_product - my_cross_product;

    const double relative_velocity_modulus = DEM_MODULUS_3(relative_velocity);

    // Calculate rolling friction moment
    mContactMoment[0] -= rolling_friction_coefficient * force * my_arm_length * relative_velocity_modulus;
    mContactMoment[1] -= rolling_friction_coefficient * force * my_arm_length * relative_velocity_modulus;
    mContactMoment[2] -= rolling_friction_coefficient * force * my_arm_length * relative_velocity_modulus;
  }

  void DEMRollingFrictionModelViscousTorque::ComputeRollingFrictionWithWall(SphericParticle* p_element, Condition* const wall, double LocalCoordSystem_2[3], double LocalContactForce[3], double indentation, array_1d<double, 3>& mContactMoment)
  {
    // Get rolling friction coefficient
    Properties& r_properties = p_element->GetProperties().GetSubProperties(wall->GetProperties().Id());
    const double rolling_friction_coefficient = r_properties[ROLLING_FRICTION];

    // Get normal contact force
    const double force = std::abs(LocalContactForce[2]);

    // Calculate arm length (particle radius discounted by identation)
    const double arm_length = p_element->GetRadius() - indentation;

    // Calculate relative translational velocity at the contact due to relative rotation
    // ATTENTION: It considers only the motion of the particle!
    array_1d<double, 3> arm_vector;
    arm_vector[0] = -LocalCoordSystem_2[0] * arm_length;
    arm_vector[1] = -LocalCoordSystem_2[1] * arm_length;
    arm_vector[2] = -LocalCoordSystem_2[2] * arm_length;
    
    const array_1d<double, 3>& angular_vel = p_element->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);

    array_1d<double, 3> cross_product;
    GeometryFunctions::CrossProduct(angular_vel, arm_vector, cross_product);
    const double relative_velocity_modulus = DEM_MODULUS_3(cross_product);

    // Calculate rolling friction moment
    mContactMoment[0] -= rolling_friction_coefficient * force * arm_length * relative_velocity_modulus;
    mContactMoment[1] -= rolling_friction_coefficient * force * arm_length * relative_velocity_modulus;
    mContactMoment[2] -= rolling_friction_coefficient * force * arm_length * relative_velocity_modulus;
  }

} // namespace Kratos