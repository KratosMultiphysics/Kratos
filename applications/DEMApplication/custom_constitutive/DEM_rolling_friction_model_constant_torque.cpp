/////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: chengshun.shang1996@gmail.com
// Date: Aug 2022
/////////////////////////////////////////////////

//This model can be found as Model Type A in [Jun Ai, 2011, Assessment of rolling resistance models in discrete element simulations]

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

    void DEMRollingFrictionModelConstantTorque::ComputeRollingFriction(SphericParticle* p_element, SphericParticle* p_neighbor, double LocalContactForce[3], array_1d<double, 3>& mContactMoment, double indentation)
    {
        array_1d<double, 3> elementRelAngularVelocity;
        noalias(elementRelAngularVelocity) = p_element->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY) - p_neighbor->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
        if (elementRelAngularVelocity[0] || elementRelAngularVelocity[1] || elementRelAngularVelocity[2]){
            const double my_young      = p_element->GetYoung(); 
            const double other_young   = p_neighbor->GetYoung();
            const double my_arm_length = p_element->GetInteractionRadius() - indentation * other_young / (other_young + my_young);
            //TODO: This should be updated for non-spherical particle
            double contact_center_point_to_element1_mass_center_distance = my_arm_length; //Here, this only works for sphere particles
                       
            GeometryFunctions::normalize(elementRelAngularVelocity);

            array_1d<double, 3> elementRelAngularVelocity_normalise = elementRelAngularVelocity;

            Properties& properties_of_this_contact = p_element->GetProperties().GetSubProperties(p_neighbor->GetProperties().Id());

            mContactMoment[0] -= elementRelAngularVelocity_normalise[0] * fabs(LocalContactForce[2]) * contact_center_point_to_element1_mass_center_distance * properties_of_this_contact[ROLLING_FRICTION]; 

            mContactMoment[1] -= elementRelAngularVelocity_normalise[1] * fabs(LocalContactForce[2]) * contact_center_point_to_element1_mass_center_distance * properties_of_this_contact[ROLLING_FRICTION]; 

            mContactMoment[2] -= elementRelAngularVelocity_normalise[2] * fabs(LocalContactForce[2]) * contact_center_point_to_element1_mass_center_distance * properties_of_this_contact[ROLLING_FRICTION]; 

        } 
    }

    void DEMRollingFrictionModelConstantTorque::ComputeRollingFrictionWithWall(double LocalContactForce[3], SphericParticle* p_element, Condition* const wall, double indentation, array_1d<double, 3>& mContactMoment)
    {
        array_1d<double, 3> element1AngularVelocity;
        noalias(element1AngularVelocity) = p_element->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
        if (element1AngularVelocity[0] || element1AngularVelocity[1] || element1AngularVelocity[2]){

            const double my_young      = p_element->GetYoung(); 
            const double walls_young   = wall->GetProperties()[YOUNG_MODULUS];
            
            double arm_length = p_element->GetInteractionRadius() - indentation* walls_young / (walls_young + my_young); 
            
            GeometryFunctions::normalize(element1AngularVelocity);

            array_1d<double, 3> element1AngularVelocity_normalise = element1AngularVelocity;

            Properties& properties_of_this_contact = p_element->GetProperties().GetSubProperties(wall->GetProperties().Id());

            mContactMoment[0] -= element1AngularVelocity_normalise[0] * fabs(LocalContactForce[2]) * arm_length * properties_of_this_contact[ROLLING_FRICTION]; 

            mContactMoment[1] -= element1AngularVelocity_normalise[1] * fabs(LocalContactForce[2]) * arm_length * properties_of_this_contact[ROLLING_FRICTION]; 

            mContactMoment[2] -= element1AngularVelocity_normalise[2] * fabs(LocalContactForce[2]) * arm_length * properties_of_this_contact[ROLLING_FRICTION]; 

        } 
    }


}//namespace Kratos