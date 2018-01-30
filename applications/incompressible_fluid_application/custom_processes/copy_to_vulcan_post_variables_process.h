// Kratos Multi-Physics
// 
// Copyright (c) 2015, Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
// 
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer 
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement: 
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
// 	
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY 
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#if !defined(KRATOS_COPY_TO_VULCAN_POST_VARIABLES_PROCESS_INCLUDED)
#define KRATOS_COPY_TO_VULCAN_POST_VARIABLES_PROCESS_INCLUDED


#include <string>
#include <iostream>
#include <algorithm>

#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "processes/process.h"
#include "utilities/geometry_utilities.h"
#include "geometries/triangle_2d_3.h"
#include "incompressible_fluid_application.h"
#include "includes/c2c_variables.h"

namespace Kratos
{ 
	class CopyToVulcanPostVariablesProcess 
		: public Process 
	 {
	   public:

	      CopyToVulcanPostVariablesProcess(ModelPart& rThisModelPart)
			:Process(), mrModelPart(rThisModelPart)
		{
		}

	      /// Destructor.
	      virtual ~CopyToVulcanPostVariablesProcess()
		{
		}
      

	      ///@}
	      ///@name Operators 
	      ///@{

	      void operator()()
		{
		  Execute();
		}

		int Check()
		{
            if (mrModelPart.NodesBegin()->SolutionStepsDataHas(MATERIAL) == false)
                KRATOS_THROW_ERROR(std::logic_error, "ERROR! Add MATERIAL variable!!!!!! ", "");
            if (mrModelPart.NodesBegin()->SolutionStepsDataHas(DISTANCE) == false)
                KRATOS_THROW_ERROR(std::logic_error, "ERROR! Add DISTANCE variable!!!!!! ", "");
            if (mrModelPart.NodesBegin()->SolutionStepsDataHas(PRESSURE) == false)
                KRATOS_THROW_ERROR(std::logic_error, "ERROR! Add PRESSURE variable!!!!!! ", "");
            if (mrModelPart.NodesBegin()->SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_THROW_ERROR(std::logic_error, "ERROR! Add VELOCITY variable!!!!!!", "");
            if (mrModelPart.NodesBegin()->SolutionStepsDataHas(VELOCITIES) == false)
                KRATOS_THROW_ERROR(std::logic_error, "ERROR! Add VELOCITIES variable!!!!!!", "");
            if (mrModelPart.NodesBegin()->SolutionStepsDataHas(TEMPERATURE) == false)
                KRATOS_THROW_ERROR(std::logic_error, "ERROR! Add VELOCITIES variable!!!!!!", "");

			return 0;
		}
		
	   virtual void Execute()
		 {
			  Check();
			  double max_distance = 0.00;
			  double min_distance = -1.00e-6;
			  //	  double max_temp = mrModelPart.GetProcessInfo()[AUX_INDEX];
			  //double min_temp = mrModelPart.GetProcessInfo()[AMBIENT_TEMPERATURE];

			 for(ModelPart::NodeIterator i_node = mrModelPart.NodesBegin() ; i_node != mrModelPart.NodesEnd() ; i_node++)
			 {
				 double distance = i_node->GetSolutionStepValue(DISTANCE);
				 if(distance > max_distance)
					 max_distance = distance;
				 else if(distance < min_distance)
					 min_distance = distance;

			 }

			 double distance_norm_inverse = 1.00 / (max_distance - min_distance);
			 //find mean interface temperature
			 //double mean_interface_T = 0.0;
			 //double cnt=0.0;
			 //for(ModelPart::NodeIterator i_node = mrModelPart.NodesBegin() ; i_node != mrModelPart.NodesEnd() ; i_node++)
			 //{
				// double distance = i_node->GetSolutionStepValue(DISTANCE);
				// if(distance <= 0.0 && distance > 2.0*min_distance )
				// {
				//	 cnt += 1.0;
				//	 mean_interface_T += i_node->FastGetSolutionStepValue(TEMPERATURE);
				// }
			 //}

			 //if(cnt != 0.0)
				// mean_interface_T/=cnt;
			 //else{
				// KRATOS_WATCH("INSIDE COPY TO VULCAN >>>> NO interface node is found");
				// mean_interface_T = mrModelPart.GetProcessInfo()[SOLID_TEMPERATURE];
			 //    }



			 for(ModelPart::NodeIterator i_node = mrModelPart.NodesBegin() ; i_node != mrModelPart.NodesEnd() ; i_node++)
			 {
				 double distance = i_node->GetSolutionStepValue(DISTANCE) * distance_norm_inverse;
				 // the distance is between -1 and 1 where < 0 is material while
				 // the material is between 0 and 1 with >=0.5 is material!
				 i_node->FastGetSolutionStepValue(MATERIAL) = (1.00 - distance) * 0.5;
				 
				 if(distance <= 0)
				 {
					 i_node->FastGetSolutionStepValue(VELOCITIES) = i_node->FastGetSolutionStepValue(VELOCITY);
					 i_node->FastGetSolutionStepValue(PRESSURES) = i_node->FastGetSolutionStepValue(PRESSURE);

					 double temp = i_node->FastGetSolutionStepValue(TEMPERATURE);
					/* if(temp > max_temp)
							temp = max_temp;
					 if (temp < min_temp)
							temp = min_temp;*/
					 i_node->FastGetSolutionStepValue(TEMPERATURES) = temp;

				 }
				 else
				 {
					 i_node->FastGetSolutionStepValue(VELOCITIES) = ZeroVector(3);
					 i_node->FastGetSolutionStepValue(PRESSURES) = 0.00;
                                         i_node->FastGetSolutionStepValue(LAST_AIR) = 42.0f; // charlie Borrar esto que lo he puesto para debugear
                                         i_node->FastGetSolutionStepValue(MOULD_VFACT) = 42.0f; // charlie Borrar esto que lo he puesto para debugear
					 //i_node->FastGetSolutionStepValue(TEMPERATURES) = i_node->FastGetSolutionStepValue(TEMPERATURE);
				 }

			 }
		 }

		private:
			ModelPart& mrModelPart;
	      ///@}

	};

}//namespace kratos

#endif // KRATOS_COPY_TO_VULCAN_POST_VARIABLES_PROCESS_INCLUDED
