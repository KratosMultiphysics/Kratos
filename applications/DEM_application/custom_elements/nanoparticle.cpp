//
// Author: Miguel Angel Celigueta    maceli@cimne.upc.edu
//

// System includes
#include <string>
#include <iostream>
#include <iomanip> // to improve std::cout precision

// Project includes
#include "nanoparticle.h"

namespace Kratos {
      
      /*void NanoParticle::ComputeAdditionalForces(array_1d<double, 3>& additionally_applied_force,
                                                      array_1d<double, 3>& additionally_applied_moment,
                                                      ProcessInfo& r_current_process_info,
                                                      const array_1d<double,3>& gravity)
        {
            KRATOS_TRY
            
            array_1d<double, 3> brownian_motion_force; brownian_motion_force.clear();
            array_1d<double, 3> van_der_waals_force; van_der_waals_force.clear();
            array_1d<double, 3> double_layer_force; double_layer_force.clear();

            ComputeBrownianMotionForce(brownian_motion_force, r_current_process_info);
            ComputeVanDerWaalsForce(van_der_waals_force, r_current_process_info);
            ComputeDoubleLayerForce(double_layer_force, r_current_process_info);
            
            additionally_applied_force += brownian_motion_force + van_der_waals_force + double_layer_force;


            //Now add the contribution of base class function (gravity or other forces added in upper levels):
            SphericParticle::ComputeAdditionalForces(additionally_applied_force, additionally_applied_moment, r_current_process_info, gravity);

            KRATOS_CATCH( "" )
        }
      
      
        void NanoParticle::ComputeBrownianMotionForce(array_1d<double, 3>& brownian_motion_force, ProcessInfo& r_current_process_info){
            KRATOS_TRY
            
            
            KRATOS_CATCH( "" )
        }
        void NanoParticle::ComputeVanDerWaalsForce(array_1d<double, 3>& van_der_waals_force, ProcessInfo& r_current_process_info){
            KRATOS_TRY
            
            
            KRATOS_CATCH( "" )            
        }
        void NanoParticle::ComputeDoubleLayerForce(array_1d<double, 3>& double_layer_force, ProcessInfo& r_current_process_info){
            KRATOS_TRY
            
            
            KRATOS_CATCH( "" )            
        }*/
      
    
} // namespace Kratos
