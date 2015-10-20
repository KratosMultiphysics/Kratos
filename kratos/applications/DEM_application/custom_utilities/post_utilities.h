#ifndef POST_UTILITIES_H
#define POST_UTILITIES_H

// System includes

// Project includes

#include "utilities/timer.h"
#include "custom_utilities/create_and_destroy.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_elements/Particle_Contact_Element.h"
#include "includes/variables.h"
#include "DEM_application.h"

/* System includes */

#include <limits>
#include <iostream>
#include <iomanip>
#include <iostream>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "utilities/openmp_utils.h"

namespace Kratos {
    
class PostUtilities {
    
public:
    
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::NodesContainerType NodesContainerType;

    /// Default constructor.       
    
    PostUtilities() {};

    /// Destructor.

    virtual ~PostUtilities() {};
    
    void AddModelPartToModelPart(ModelPart& rCompleteModelPart, ModelPart& rModelPartToAdd)
    {
        ////WATCH OUT! This function respects the existing Id's!
        KRATOS_TRY;      
      
        //preallocate the memory needed
        int tot_nodes = rCompleteModelPart.Nodes().size() + rModelPartToAdd.GetCommunicator().LocalMesh().Nodes().size();
        int tot_elements = rCompleteModelPart.Elements().size() + rModelPartToAdd.GetCommunicator().LocalMesh().Elements().size();
        rCompleteModelPart.Nodes().reserve( tot_nodes );
        rCompleteModelPart.Elements().reserve( tot_elements );
        for (ModelPart::NodesContainerType::ptr_iterator node_it = rModelPartToAdd.GetCommunicator().LocalMesh().Nodes().ptr_begin(); node_it != rModelPartToAdd.GetCommunicator().LocalMesh().Nodes().ptr_end(); node_it++)
	{			
            rCompleteModelPart.Nodes().push_back( *node_it );
	}
      
        for (ModelPart::ElementsContainerType::ptr_iterator elem_it = rModelPartToAdd.GetCommunicator().LocalMesh().Elements().ptr_begin(); elem_it != rModelPartToAdd.GetCommunicator().LocalMesh().Elements().ptr_end(); elem_it++)
	{			
	    rCompleteModelPart.Elements().push_back( *elem_it );          
	}
      
        KRATOS_CATCH("");
    
    }    
    
    
    array_1d<double,3> VelocityTrap(ModelPart& rModelPart, const array_1d<double,3>& low_point, const array_1d<double,3>& high_point){
        
        ElementsArrayType& pElements = rModelPart.GetCommunicator().LocalMesh().Elements();

        OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), pElements.size(), this->GetElementPartition());          

        double velocity_X = 0.0, velocity_Y = 0.0, velocity_Z = 0.0;
        
        int number_of_elements = 0;
          
        #pragma omp parallel for reduction(+: velocity_X, velocity_Y, velocity_Z, number_of_elements)
        
        for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){
            
            ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
            ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

            for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it) {
                 
                array_1d<double,3> coor = (it)->GetGeometry()[0].Coordinates();                  
                            
                if (coor[0] >= low_point[0] && coor[0] <= high_point[0] && 
                    coor[1] >= low_point[1] && coor[1] <= high_point[1] && 
                    coor[2] >= low_point[2] && coor[2] <= high_point[2]) {
                    
                    velocity_X += (it)->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X);
                    velocity_Y += (it)->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y);
                    velocity_Z += (it)->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Z);
                    
                    number_of_elements++;
                }
                  
            } //elements for
            
            for (int i = 0; i < 3; ++i) {

                if (high_point[i] < low_point[i]) {
                    KRATOS_THROW_ERROR(std::logic_error, "Check the limits of the Velocity Trap Box. Maximum coordinates smaller than minimum coordinates.", "");
                }
            }

        } //parallel for
              
        if (number_of_elements) {
            velocity_X /= number_of_elements;
            velocity_Y /= number_of_elements;
            velocity_Z /= number_of_elements;
        }
        
        array_1d<double,3> velocity; 
        
        velocity[0] = velocity_X;
        velocity[1] = velocity_Y;
        velocity[2] = velocity_Z;
        
        return velocity;
        
    }//VelocityTrap
    
    
    void IntegrationOfForces(ModelPart::NodesContainerType& mesh_nodes , array_1d<double, 3>& total_forces, array_1d<double, 3>& rotation_center,
                             array_1d<double, 3>& total_moment) {
                
        for (ModelPart::NodesContainerType::ptr_iterator node_pointer_it = mesh_nodes.ptr_begin();
            node_pointer_it != mesh_nodes.ptr_end(); ++node_pointer_it) {
                
            const array_1d<double, 3>& contact_forces_summed_at_structure_point = (*node_pointer_it)->FastGetSolutionStepValue(CONTACT_FORCES);
            noalias(total_forces) += contact_forces_summed_at_structure_point;
            
            array_1d<double, 3> vector_from_structure_center_to_structure_point;
            noalias(vector_from_structure_center_to_structure_point) = (*node_pointer_it)->Coordinates() - rotation_center;
            
            array_1d<double, 3> moment_to_add;
            GeometryFunctions::CrossProduct(vector_from_structure_center_to_structure_point, contact_forces_summed_at_structure_point, moment_to_add);
            
            noalias(total_moment) += moment_to_add;
        }
    }
    
    
    double QuasiStaticAdimensionalNumber(ModelPart& rParticlesModelPart, ModelPart& rContactModelPart, ProcessInfo& rCurrentProcessInfo ){
        
        double adimensional_value = 0.0;

        ElementsArrayType& pParticleElements = rParticlesModelPart.GetCommunicator().LocalMesh().Elements();
        
        OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), pParticleElements.size(), this->GetElementPartition());          
  
        array_1d<double,3> particle_forces;
        
        const array_1d<double,3>& gravity = rCurrentProcessInfo[GRAVITY];
        
        double total_force = 0.0;
        
        //#pragma omp parallel for
        #pragma omp parallel for reduction(+:total_force)
        for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){
            ElementsArrayType::iterator it_begin = pParticleElements.ptr_begin() + this->GetElementPartition()[k];
            ElementsArrayType::iterator it_end   = pParticleElements.ptr_begin() + this->GetElementPartition()[k + 1];

            for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it) {    
                             
                Element::GeometryType& geom = it->GetGeometry();
              
                if( geom[0].IsNot(DEMFlags::FIXED_VEL_X) && geom[0].IsNot(DEMFlags::FIXED_VEL_Y) && geom[0].IsNot(DEMFlags::FIXED_VEL_Z) )
                {
                    particle_forces  = geom[0].FastGetSolutionStepValue(TOTAL_FORCES);
                    double mass = geom[0].FastGetSolutionStepValue(NODAL_MASS);

                
                    particle_forces[0] += mass * gravity[0];
                    particle_forces[1] += mass * gravity[1];
                    particle_forces[2] += mass * gravity[2];

                    double module = 0.0;
                    GeometryFunctions::module(particle_forces, module);
                
                    total_force += module;
                
                } //if
              
            }//balls
            
        }//paralel 
          
        ElementsArrayType& pContactElements        = rContactModelPart.GetCommunicator().LocalMesh().Elements();

        OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), pContactElements.size(), this->GetElementPartition());          
  
        array_1d<double,3> contact_forces;
        double total_elastic_force = 0.0;
        
        #pragma omp parallel for reduction(+:total_elastic_force)
        for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){
            ElementsArrayType::iterator it_begin = pContactElements.ptr_begin() + this->GetElementPartition()[k];
            ElementsArrayType::iterator it_end   = pContactElements.ptr_begin() + this->GetElementPartition()[k + 1];

            for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                                                        
                Element::GeometryType& geom = it->GetGeometry();
              
                if (geom[0].IsNot(DEMFlags::FIXED_VEL_X) && geom[0].IsNot(DEMFlags::FIXED_VEL_Y) && geom[0].IsNot(DEMFlags::FIXED_VEL_Z) &&
                    geom[1].IsNot(DEMFlags::FIXED_VEL_X) && geom[1].IsNot(DEMFlags::FIXED_VEL_Y) && geom[1].IsNot(DEMFlags::FIXED_VEL_Z)) {                
   
                    contact_forces  = it->GetValue(LOCAL_CONTACT_FORCE);                
                    double module = 0.0;
                    GeometryFunctions::module(contact_forces, module);
                    total_elastic_force += module;                
                }              
            }  
        }
   
        if (total_elastic_force != 0.0)
        {
            adimensional_value =  total_force/total_elastic_force;   
        }
        else
        {
            KRATOS_THROW_ERROR(std::runtime_error,"There are no elastic forces= ", total_elastic_force)
        }
  
        return adimensional_value;
        
    }//QuasiStaticAdimensionalNumber
    
    vector<unsigned int>& GetElementPartition(){return (mElementPartition);};
    
protected:
    
    vector<unsigned int> mElementPartition;
    
}; // Class PostUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // POST_UTILITIES_H
