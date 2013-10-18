#ifndef POST_UTILITIES_H
#define POST_UTILITIES_H

// /* External includes */

// System includes

// Project includes
#include "utilities/timer.h"
#include "custom_utilities/create_and_destroy.h"

#include "custom_elements/spheric_swimming_particle.h"
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

namespace Kratos
{
class PostUtilities
{
public:
    
    typedef ModelPart::ElementsContainerType                          ElementsArrayType;

    /// Default constructor.       
    
    PostUtilities() {};

    /// Destructor.

    virtual ~PostUtilities() {};
    
    array_1d<double, 3> VelocityTrap(ModelPart& rModelPart, const array_1d<double, 3> low_point, const array_1d<double, 3> high_point){
        
          ElementsArrayType& pElements        = rModelPart.GetCommunicator().LocalMesh().Elements();

          OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), pElements.size(), this->GetElementPartition());          

          double velocity_X = 0.0, velocity_Y = 0.0, velocity_Z = 0.0;
          int number_of_elements = 0;
          
          //#pragma omp parallel for
          #pragma omp parallel for reduction(+:velocity_X,velocity_Y,velocity_Z,number_of_elements) 
          for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){
              ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
              ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

              for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                 
                 array_1d<double, 3 > coor = (it)->GetGeometry()(0)->Coordinates();                  
                            
                 if( coor[0] >= low_point[0] && coor[0] <= high_point[0]  && 
                     coor[1] >= low_point[1] && coor[1] <= high_point[1]  && 
                     coor[2] >= low_point[2] && coor[2] <= high_point[2] ) {
                     
                     velocity_X += (it)->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY_X);
                     velocity_Y += (it)->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY_Y);
                     velocity_Z += (it)->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY_Z);
                     number_of_elements ++;
                 }
                  
              } //elements for

          } //parallel for
              
          if (number_of_elements){
              velocity_X /= number_of_elements;
              velocity_Y /= number_of_elements;
              velocity_Z /= number_of_elements;
          }
          array_1d<double, 3> velocity; 
          velocity[0]= velocity_X;
          velocity[1]= velocity_Y;
          velocity[2]= velocity_Z;
          return velocity;
        
    }//VelocityTrap
    
    vector<unsigned int>&    GetElementPartition(){return (mElementPartition);};
    
    protected:

      
        vector<unsigned int>                mElementPartition;
    
    

}; // Class PostUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{



} // namespace Kratos.

#endif // POST_UTILITIES_H
