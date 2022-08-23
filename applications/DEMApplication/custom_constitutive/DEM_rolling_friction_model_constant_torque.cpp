/////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: chengshun.shang1996@gmail.com
// Date: Aug 2022
/////////////////////////////////////////////////

#include "DEM_rolling_friction_model_constant_torque.h"

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

    void DEMRollingFrictionModelConstantTorque::ComputeRollingFriction(SphericParticle* p_element, SphericParticle* p_neighbor, double LocalContactForce[3], array_1d<double, 3>& mContactMoment)
    {
        array_1d<double, 3> elementRelAngularVelocity;
        noalias(elementRelAngularVelocity) = p_element->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY) - p_neighbor->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
        if (elementRelAngularVelocity[0] || elementRelAngularVelocity[1] || elementRelAngularVelocity[2]){
            array_1d<double, 3> other_to_me_vect;
            noalias(other_to_me_vect) = p_element->GetGeometry()[0].Coordinates() - p_neighbor->GetGeometry()[0].Coordinates();
            double bond_center_point_to_element1_mass_center_distance = DEM_MODULUS_3(other_to_me_vect) / 2; //Here, this only works for sphere particles
            
            double elementRelAngularVelocity_modulus = sqrt(elementRelAngularVelocity[0] * elementRelAngularVelocity[0] + 
                                                    elementRelAngularVelocity[1] * elementRelAngularVelocity[1] +
                                                    elementRelAngularVelocity[2] * elementRelAngularVelocity[2]);

            array_1d<double, 3> elementRelAngularVelocity_normalise;
            elementRelAngularVelocity_normalise[0] = elementRelAngularVelocity[0] / elementRelAngularVelocity_modulus;
            elementRelAngularVelocity_normalise[1] = elementRelAngularVelocity[1] / elementRelAngularVelocity_modulus;
            elementRelAngularVelocity_normalise[2] = elementRelAngularVelocity[2] / elementRelAngularVelocity_modulus;

            Properties& properties_of_this_contact = p_element->GetProperties().GetSubProperties(p_neighbor->GetProperties().Id());

            mContactMoment[0] -= elementRelAngularVelocity_normalise[0] * fabs(LocalContactForce[2]) * bond_center_point_to_element1_mass_center_distance * properties_of_this_contact[ROLLING_FRICTION]; 

            mContactMoment[1] -= elementRelAngularVelocity_normalise[1] * fabs(LocalContactForce[2]) * bond_center_point_to_element1_mass_center_distance * properties_of_this_contact[ROLLING_FRICTION]; 

            mContactMoment[2] -= elementRelAngularVelocity_normalise[2] * fabs(LocalContactForce[2]) * bond_center_point_to_element1_mass_center_distance * properties_of_this_contact[ROLLING_FRICTION]; 

        } 
    }

    void DEMRollingFrictionModelConstantTorque::ComputeRollingFrictionWithWall(double LocalContactForce[3], SphericParticle* p_element, Condition* const wall, double indentation, array_1d<double, 3>& mContactMoment)
    {
        array_1d<double, 3> element1AngularVelocity;
        noalias(element1AngularVelocity) = p_element->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
        if (element1AngularVelocity[0] || element1AngularVelocity[1] || element1AngularVelocity[2]){

            double arm_length = p_element->GetInteractionRadius() - indentation; 
            
            double element1AngularVelocity_modulus = sqrt(element1AngularVelocity[0] * element1AngularVelocity[0] + 
                                                    element1AngularVelocity[1] * element1AngularVelocity[1] +
                                                    element1AngularVelocity[2] * element1AngularVelocity[2]);

            array_1d<double, 3> element1AngularVelocity_normalise;
            element1AngularVelocity_normalise[0] = element1AngularVelocity[0] / element1AngularVelocity_modulus;
            element1AngularVelocity_normalise[1] = element1AngularVelocity[1] / element1AngularVelocity_modulus;
            element1AngularVelocity_normalise[2] = element1AngularVelocity[2] / element1AngularVelocity_modulus;

            Properties& properties_of_this_contact = p_element->GetProperties().GetSubProperties(wall->GetProperties().Id());

            mContactMoment[0] -= element1AngularVelocity_normalise[0] * fabs(LocalContactForce[2]) * arm_length * properties_of_this_contact[ROLLING_FRICTION]; 

            mContactMoment[1] -= element1AngularVelocity_normalise[1] * fabs(LocalContactForce[2]) * arm_length * properties_of_this_contact[ROLLING_FRICTION]; 

            mContactMoment[2] -= element1AngularVelocity_normalise[2] * fabs(LocalContactForce[2]) * arm_length * properties_of_this_contact[ROLLING_FRICTION]; 

        } 
    }


}//namespace Kratos