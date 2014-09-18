/*
==============================================================================
KratosPFEMApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-11-19 15:38:01 $
//   Revision:            $Revision: 1.1 $
//
//
#if !defined(KRATOS_COPY_TO_VULCAN_POST_VARIABLES_PROCESS_INCLUDED)
#define KRATOS_COPY_TO_VULCAN_POST_VARIABLES_PROCESS_INCLUDED


#include <string>
#include <iostream>
#include <algorithm>

#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
#include "geometries/triangle_2d_3.h"
//#include "kratos/applications/MeshingApplication/meshing_application.h"

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
                KRATOS_ERROR(std::logic_error, "ERROR! Add MATERIAL variable!!!!!! ", "");
            if (mrModelPart.NodesBegin()->SolutionStepsDataHas(DISTANCE) == false)
                KRATOS_ERROR(std::logic_error, "ERROR! Add DISTANCE variable!!!!!! ", "");
            if (mrModelPart.NodesBegin()->SolutionStepsDataHas(PRESSURE) == false)
                KRATOS_ERROR(std::logic_error, "ERROR! Add PRESSURE variable!!!!!! ", "");
            if (mrModelPart.NodesBegin()->SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_ERROR(std::logic_error, "ERROR! Add VELOCITY variable!!!!!!", "");
            if (mrModelPart.NodesBegin()->SolutionStepsDataHas(VELOCITIES) == false)
                KRATOS_ERROR(std::logic_error, "ERROR! Add VELOCITIES variable!!!!!!", "");
            if (mrModelPart.NodesBegin()->SolutionStepsDataHas(TEMPERATURE) == false)
                KRATOS_ERROR(std::logic_error, "ERROR! Add VELOCITIES variable!!!!!!", "");

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
