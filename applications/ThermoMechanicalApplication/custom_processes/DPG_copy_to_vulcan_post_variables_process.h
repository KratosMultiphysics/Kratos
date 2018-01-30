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
#if !defined(KRATOS_DPG_COPY_TO_VULCAN_POST_VARIABLES_PROCESS_INCLUDED)
#define KRATOS_DPG_COPY_TO_VULCAN_POST_VARIABLES_PROCESS_INCLUDED


#include <string>
#include <iostream>
#include <algorithm>

#include "includes/define.h"
#include "includes/c2c_variables.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
#include "geometries/triangle_2d_3.h"
//#include "kratos/applications/MeshingApplication/meshing_application.h"

namespace Kratos
{ 
	class DPGCopyToVulcanPostVariablesProcess 
		: public Process 
	 {
	   public:

	      DPGCopyToVulcanPostVariablesProcess(ModelPart& rThisModelPart, const int filling, const int solidification, const int SI )
			:Process(), mrModelPart(rThisModelPart), mrfilling(filling), mrsolidificaion(solidification),mrUnits(SI)
		{
		}

	      /// Destructor.
	      virtual ~DPGCopyToVulcanPostVariablesProcess()
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
            if (mrModelPart.NodesBegin()->SolutionStepsDataHas(SOLIDFRACTION) == false)
                KRATOS_THROW_ERROR(std::logic_error, "ERROR! Add SOLIDFRACTION variable!!!!!!", "");

			return 0;
		}
		
	   virtual void Execute()
		 {
		   if(mrfilling == 1){
			  Check();
			  double max_distance = 0.00;
			  double min_distance = -1.00e-6;
			  //double max_temp = mrModelPart.GetProcessInfo()[AUX_INDEX];
			 // double min_temp = mrModelPart.GetProcessInfo()[AMBIENT_TEMPERATURE];
			 // double min_mat_temp = 100000.0;

		//	 for(ModelPart::NodeIterator i_node = mrModelPart.NodesBegin() ; i_node != mrModelPart.NodesEnd() ; i_node++)
			 #pragma omp parallel for
		        for (int k = 0; k< static_cast<int> (mrModelPart.Nodes().size()); k++)
		        {
				ModelPart::NodesContainerType::iterator i_node = mrModelPart.NodesBegin() + k;
				 double distance = i_node->GetSolutionStepValue(DISTANCE);
				// double slip_flag = i_node->FastGetSolutionStepValue(IS_SLIP);
				 if(distance > max_distance)
					 max_distance = distance;
				 else if(distance < min_distance)
					 min_distance = distance;

				 //save max_vel in each node
				 if(distance<0.0){
				 const  array_1d<double, 3 > node_vel = i_node->FastGetSolutionStepValue(VELOCITY);
				 double& max_vel = i_node->FastGetSolutionStepValue(MAX_VEL);
				 double nomr_vel = norm_2(node_vel);
				 if( nomr_vel> max_vel )
					 max_vel = nomr_vel;
				 }

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
			 //Find Meean_interface_temp
			 int elem_size = mrModelPart.Elements().size();
			 double mean_interface_temp = 0.0;
			 double cnt_intr = 1.0;
			 for (int i = 0; i < elem_size; i++)
				{
					PointerVector< Element>::iterator it_elem = mrModelPart.ElementsBegin() + i;
					//KRATOS_WATCH(it_elem->GetValue(AUX_INDEX));
					if(it_elem->GetValue(AUX_INDEX) == 1.0){
						Geometry<Node < 3 > >& element_geometry = it_elem->GetGeometry();
						for( int ii = 0; ii < 4; ii++){
							double nd_dist = element_geometry[ii].FastGetSolutionStepValue(DISTANCE);
							//KRATOS_WATCH(nd_dist);
							if (nd_dist < 0.0){
								mean_interface_temp += element_geometry[ii].FastGetSolutionStepValue(TEMPERATURE);
								cnt_intr += 1.0;
							}
						}

					}

			    }	
			 mean_interface_temp /= cnt_intr;

			 double min_wet_temp = mrModelPart.GetProcessInfo()[FLUID_TEMPERATURE]; // TO-DO someday change the temperature in the inlet
			 min_wet_temp=1000000.0; // this a temporal fixing
			 double max_wet_temp = mrModelPart.GetProcessInfo()[SOLID_TEMPERATURE];

			 #pragma omp parallel for
		        for (int k = 0; k< static_cast<int> (mrModelPart.Nodes().size()); k++)
		        {
				ModelPart::NodesContainerType::iterator i_node = mrModelPart.NodesBegin() + k;
			 	 double distance = i_node->GetSolutionStepValue(DISTANCE);
				 double slip_falg = i_node->GetSolutionStepValue(IS_SLIP);
				 if(distance <=0.0 && slip_falg!=20.0 && slip_falg!=30.0)
				 {
					 double temp_node = i_node->FastGetSolutionStepValue(TEMPERATURE);
					 if (temp_node < min_wet_temp)
						 min_wet_temp = temp_node;
					 else if(temp_node > max_wet_temp)
						max_wet_temp = temp_node;
				 }
				}
				//min_wet_temp = 0.998*min_wet_temp + 0.00*max_wet_temp;
//			 for(ModelPart::NodeIterator i_node = mrModelPart.NodesBegin() ; i_node != mrModelPart.NodesEnd() ; i_node++)



			(mrUnits == 1) ? SIUnitsOutput(distance_norm_inverse,min_wet_temp, min_wet_temp-(max_wet_temp-min_wet_temp)*0.10) : USUnitsOutput(distance_norm_inverse, min_wet_temp, min_wet_temp-(max_wet_temp-min_wet_temp)*0.10);



		   }

		   if(mrsolidificaion == 1)
			   CopySolidificationVariables(mrModelPart);
		 }

		private:
			ModelPart& mrModelPart;
			const int mrfilling;
			const int mrsolidificaion;
			const int mrUnits;

			void CopySolidificationVariables(ModelPart& rThisModelPart)
			{

			if(mrUnits == 1){
			 #pragma omp parallel for
		     for (int k = 0; k< static_cast<int> (mrModelPart.Nodes().size()); k++)
		      {
			    ModelPart::NodesContainerType::iterator i_node = mrModelPart.NodesBegin() + k;

				i_node->FastGetSolutionStepValue(TEMPERATURES) = i_node->FastGetSolutionStepValue(TEMPERATURE);
				i_node->FastGetSolutionStepValue(SOLIDFRACTION) = i_node->FastGetSolutionStepValue(SOLIDFRACTION);

			  }
			}
			else
			{
			 #pragma omp parallel for
		     for (int k = 0; k< static_cast<int> (mrModelPart.Nodes().size()); k++)
		      {
			    ModelPart::NodesContainerType::iterator i_node = mrModelPart.NodesBegin() + k;

				i_node->FastGetSolutionStepValue(TEMPERATURES_US) = (1.8 * i_node->FastGetSolutionStepValue(TEMPERATURE) + 32.0);
				i_node->FastGetSolutionStepValue(SOLIDFRACTION) = i_node->FastGetSolutionStepValue(SOLIDFRACTION);

				 // SOLIDIF_MODULUS is already written en CM take care to convert it to inch
				 double si_sdf_mod = i_node->FastGetSolutionStepValue(SOLIDIF_MODULUS); 
				 i_node->FastGetSolutionStepValue(SOLIDIF_MODULUS_US) = si_sdf_mod / 2.54 ; 

				 //shrinkage_porosity
				 double si_shrinkage = i_node->FastGetSolutionStepValue(SHRINKAGE_POROSITY); 
				 i_node->FastGetSolutionStepValue(SHRINKAGE_POROSITY_US) = 16387064.0  * si_shrinkage; //m^3 -> in^3

 
			  }

			}

			}
           /////////////////////////////////////////
			void SIUnitsOutput(double distance_norm_inverse, double min_wet_temp, double const print_temperature)
			{

			 #pragma omp parallel for
		      for (int k = 0; k< static_cast<int> (mrModelPart.Nodes().size()); k++)
		       {
				ModelPart::NodesContainerType::iterator i_node = mrModelPart.NodesBegin() + k;
			 	 double distance = i_node->GetSolutionStepValue(DISTANCE) * distance_norm_inverse;
				 //double slip_falg = i_node->GetSolutionStepValue(IS_SLIP);
				 // the distance is between -1 and 1 where < 0 is material while
				 // the material is between 0 and 1 with >=0.5 is material!
				 i_node->FastGetSolutionStepValue(MATERIAL) = (1.00 - distance) * 0.5;
				 //i_node->FastGetSolutionStepValue(TEMPERATURES) = i_node->FastGetSolutionStepValue(TEMPERATURE);
				 
				 if(distance <= 0)
				 {
					 i_node->FastGetSolutionStepValue(VELOCITIES) = i_node->FastGetSolutionStepValue(VELOCITY);
					 i_node->FastGetSolutionStepValue(PRESSURES) = i_node->FastGetSolutionStepValue(PRESSURE);

					 double& temp = i_node->FastGetSolutionStepValue(TEMPERATURES);
					 temp = i_node->FastGetSolutionStepValue(TEMPERATURE);
// 					 i_node->FastGetSolutionStepValue(SOLIDFRACTION) = i_node->FastGetSolutionStepValue(DP_ALPHA1);

					 if(temp<min_wet_temp)
							 temp = min_wet_temp;
					 /*if(slip_falg==30.0 || slip_falg==20.0)
					 {
						 double temp = i_node->FastGetSolutionStepValue(TEMPERATURE);
						 if(temp<min_wet_temp)
							 i_node->FastGetSolutionStepValue(TEMPERATURES) = min_wet_temp;
					 }*/


				 }
				 else
				 {
					 i_node->FastGetSolutionStepValue(VELOCITIES) = ZeroVector(3);
					 i_node->FastGetSolutionStepValue(PRESSURES) = 0.00;
					 //i_node->FastGetSolutionStepValue(TEMPERATURES) = mean_interface_temp-10.0;
					// i_node->FastGetSolutionStepValue(TEMPERATURES) = 0.95*min_wet_temp;
					 i_node->FastGetSolutionStepValue(TEMPERATURES) = print_temperature; //0.998*min_wet_temp;
					 //i_node->FastGetSolutionStepValue(TEMPERATURES) = i_node->FastGetSolutionStepValue(TEMPERATURE);
// 					 i_node->FastGetSolutionStepValue(SOLIDFRACTION) = 0.00;
				 }

			  }


			}
         ////////////////////////////////////////////////////
			void USUnitsOutput(double distance_norm_inverse, double min_wet_temp, double const print_temperature)
			{
				min_wet_temp = min_wet_temp* 1.8 + 32.0;

			 #pragma omp parallel for
		      for (int k = 0; k< static_cast<int> (mrModelPart.Nodes().size()); k++)
		       {
				ModelPart::NodesContainerType::iterator i_node = mrModelPart.NodesBegin() + k;
			 	 double distance = i_node->GetSolutionStepValue(DISTANCE) * distance_norm_inverse;
				 //double slip_falg = i_node->GetSolutionStepValue(IS_SLIP);
				 i_node->FastGetSolutionStepValue(MATERIAL) = (1.00 - distance) * 0.5;
				 

				 if(distance <= 0)
				 {
					 i_node->FastGetSolutionStepValue(VELOCITIES) = i_node->FastGetSolutionStepValue(VELOCITY);
					 i_node->FastGetSolutionStepValue(PRESSURES) = i_node->FastGetSolutionStepValue(PRESSURE);

					 double& temp = i_node->FastGetSolutionStepValue(TEMPERATURES_US);
					 temp = (1.8 * i_node->FastGetSolutionStepValue(TEMPERATURE) + 32.0);
					 //i_node->FastGetSolutionStepValue(SOLIDFRACTION) = i_node->FastGetSolutionStepValue(DP_ALPHA1);

					 if(temp<min_wet_temp)
							 temp = min_wet_temp;


				 }
				 else
				 {
					 i_node->FastGetSolutionStepValue(VELOCITIES) = ZeroVector(3);
					 i_node->FastGetSolutionStepValue(PRESSURES) = 0.00;
					 i_node->FastGetSolutionStepValue(TEMPERATURES_US) = (1.8 * print_temperature + 32.0) ;
					 //i_node->FastGetSolutionStepValue(TEMPERATURES) = (1.8 * print_temperature + 32.0) ;
					 //i_node->FastGetSolutionStepValue(SOLIDFRACTION) = 0.00;
				 }

			  }


			}



	      ///@}

	};

}//namespace kratos

#endif // KRATOS_DPG__COPY_TO_VULCAN_POST_VARIABLES_PROCESS_INCLUDED
