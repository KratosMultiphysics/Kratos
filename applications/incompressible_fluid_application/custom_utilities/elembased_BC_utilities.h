/*
==============================================================================
KratosIncompressibleFluidApplicationApplication 
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
//   Last Modified by:    $Author: antonia $
//   Date:                $Date: 2009-01-14 16:24:38 $
//   Revision:            $Revision: 1.11 $
//
//


#if !defined(KRATOS_ELEMBASED_BC_UTILITIES_H_INCLUDED)
#define  KRATOS_ELEMBASED_BC_UTILITIES_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// #include <omp.h>


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
//#include "geometries/geometry.h"
#include "utilities/geometry_utilities.h"
#include "incompressible_fluid_application.h"


namespace Kratos
{
					
//REMEMBER to clear the AUX index loops that are now usless********************************			
	class ElemBasedBCUtilities
	{
		public:
			//constructor and destructor
			ElemBasedBCUtilities(    ModelPart& mr_model_part
				    )
			:mr_model_part(mr_model_part){};
			~ElemBasedBCUtilities(){};


			//********************************
			// Definition of the fluid domain.and of the element intersected by the free surface
			 void SetDividedElem_2D()
			 {
			 KRATOS_TRY
				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();
								inode++)	
				{	
					 inode->FastGetSolutionStepValue(AUX_INDEX) = 0;
				}
				
				for( ModelPart::ElementsContainerType::iterator iel = mr_model_part.ElementsBegin();
								iel != mr_model_part.ElementsEnd();
								iel++)
				{
					 Geometry<Node<3> >& geom = iel->GetGeometry();
					 double Area = GeometryUtils::CalculateVolume2D(geom);
					 double toll = 0.15*sqrt(Area * 2.30940108);//0.15*(h in a equailateral triangle of given area)

					 array_1d<double,3> dist = ZeroVector(3);
					 for (unsigned int i = 0; i < geom.size(); i++)
						dist[i] = geom[i].GetSolutionStepValue(DISTANCE);
						
					 //no fluid element
					 if( dist[0] > -toll && dist[1] > -toll && dist[2] > -toll  ) 	
						{iel->GetValue(IS_DIVIDED) = 0.0;
						
						}
					 //fluid element
					 //else if( dist[0] < 0.0 && dist[1] < 0.0 && dist[2] < 0.0  ) 	
					 else if( dist[0] < toll && dist[1] < toll && dist[2] < toll  ) 
						{iel->GetValue(IS_DIVIDED) = -1.0;

						}
					 //half fluid half no fluid element
					 else
						{
						      iel->GetValue(IS_DIVIDED) = 1.0;
// // 						      KRATOS_WATCH("IS DIVIDED ELEMENTS*******************")
// // 						      KRATOS_WATCH(iel->Id())
						      for(unsigned int i=0; i< geom.size() ; i++)
						      {
							     geom[i].FastGetSolutionStepValue(AUX_INDEX) = 1;
						      }
//  						     // KRATOS_WATCH(iel->Id());

						}
				
				 }
			 KRATOS_CATCH("")
			 }



			void SetPressureAndVelocityFixities()
			{
			KRATOS_TRY	
				
				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();
								inode++)	
				{	
					 if(inode->FastGetSolutionStepValue(DISTANCE) > 0.0) //
					 {
					       inode->Fix(VELOCITY_X);
					       inode->Fix(VELOCITY_Y);
					       inode->Fix(VELOCITY_Z);
					       inode->Fix(PRESSURE);
					 }
				}

				for( ModelPart::ElementsContainerType::iterator iel = mr_model_part.ElementsBegin();
								iel != mr_model_part.ElementsEnd();
								iel++)	
				{
					 Geometry<Node<3> >& geom = iel->GetGeometry();
					 if(iel->GetValue(IS_DIVIDED) == 1.0)
					 {
						for (unsigned int i = 0; i < geom.size(); i++)
						{
							 if( geom[i].FastGetSolutionStepValue(IS_STRUCTURE) != 1.0)
							 {
								geom[i].Free(VELOCITY_X);
								geom[i].Free(VELOCITY_Y);
								geom[i].Free(VELOCITY_Z);
								geom[i].Free(PRESSURE);
							 }
						}
					 }
				}
				
				
				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();
								inode++)	
				{	
					 if(inode->IsFixed(PRESSURE) == true && inode->FastGetSolutionStepValue(IS_STRUCTURE) != 1.0) 
					 {
					       inode->FastGetSolutionStepValue(VELOCITY_X) = 0.0;
					       inode->FastGetSolutionStepValue(VELOCITY_X) = 0.0;
					       inode->FastGetSolutionStepValue(VELOCITY_X) = 0.0;
					       inode->FastGetSolutionStepValue(PRESSURE) = 0.0;
					       KRATOS_WATCH(inode->Id());
					 }
				}
			 KRATOS_CATCH("")
			 }


// 			void FixPressure()
// 			{
// 			KRATOS_TRY	
// 				
// 				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
// 								inode != mr_model_part.NodesEnd();
// 								inode++)	
// 				{
// 					if( inode->FastGetSolutionStepValue(DISTANCE) > 0.0) //all the NON fluid domain
// 					{
// 						inode->FastGetSolutionStepValue(PRESSURE) = 0.0;//metti alla fine
// 						inode->Fix(PRESSURE);
// // 						KRATOS_WATCH(inode->Id());
// // 						KRATOS_WATCH(inode->FastGetSolutionStepValue(PRESSURE));
// 					}
// 				}
// 
// 				//free pressure of the nodes that are part of elements cut by the free surface (we want to include them in the calculation).
// 				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
// 								inode != mr_model_part.NodesEnd();
// 								inode++)	
// 				{
// 				    	if (inode->FastGetSolutionStepValue(DISTANCE) < 0.0)
// 					{
// 						WeakPointerVector< Node<3> >& neighb_nodes = inode->GetValue(NEIGHBOUR_NODES); 
// 						for( WeakPointerVector< Node<3> >::iterator j =	neighb_nodes.begin(); j != neighb_nodes.end(); j++) 
// 						{ 
// 							if(j->FastGetSolutionStepValue(DISTANCE) > 0.0)
// 							{
// 							    j->Free(PRESSURE);
// 							}
// 						}  
// 						
// 					}
// 				}
// 
// 
// 
// 			KRATOS_CATCH("")
// 			}
// 
// 			void FixVelocity(double extrapolation_distance)
// 			{
// 			KRATOS_TRY
// 
// 
// 				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
// 								inode != mr_model_part.NodesEnd();
// 								inode++)	
// 				{
// 					if( inode->FastGetSolutionStepValue(DISTANCE) > extrapolation_distance) //OUT of the extrapolation domain
// 					{
// 						inode->FastGetSolutionStepValue(VELOCITY_X) = 0.0;
// 						inode->FastGetSolutionStepValue(VELOCITY_Y) = 0.0;
// 						inode->FastGetSolutionStepValue(VELOCITY_Z) = 0.0;
// 						inode->Fix(VELOCITY_X);
// 						inode->Fix(VELOCITY_Y);
// 						inode->Fix(VELOCITY_Z);
// 
// 					}
// 				}
// 
// 			KRATOS_CATCH("")
// 			}
// 
// 			void FreePressure()
// 			{
// 			KRATOS_TRY	
// 
// 				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
// 								inode != mr_model_part.NodesEnd();
// 								inode++)	
// 				{
// 
// 					inode->Free(PRESSURE);
// 				}
// 
// 			KRATOS_CATCH("")
// 			}
// 
// 			void FreeVelocity(double extrapolation_distance)
// 			{
// 			KRATOS_TRY	
// 
// 				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
// 								inode != mr_model_part.NodesEnd();
// 								inode++)	
// 				{
// 					inode->Free(VELOCITY_X);
// 					inode->Free(VELOCITY_Y);
// 					inode->Free(VELOCITY_Z);				
// 				}
// 			KRATOS_CATCH("")
// 			}

			void FreePressureAndVelocity()
			{
			KRATOS_TRY	

				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();
								inode++)	
				{	
					if(inode->GetSolutionStepValue(IS_STRUCTURE) != 1.0 && inode->GetSolutionStepValue(DISTANCE) > 0.0)
					{
						inode->Free(PRESSURE);
						inode->Free(VELOCITY_X);
						inode->Free(VELOCITY_Y);
						inode->Free(VELOCITY_Z);	
/*KRATOS_WATCH("freed nodes")
KRATOS_WATCH(inode->Id())	*/		
					}
				}
			KRATOS_CATCH("")
			}

			void SetToZeroPressureAndVelocity(double extrapolation_distance)
			{
			KRATOS_TRY	
				for( ModelPart::ElementsContainerType::iterator iel = mr_model_part.ElementsBegin();
								iel != mr_model_part.ElementsEnd();
								iel++)	
				{
					 Geometry<Node<3> >& geom = iel->GetGeometry();
					 //if it is a NO fluid element
					 if(iel->GetValue(IS_DIVIDED) == 0.0)
					 {
						for (unsigned int i = 0; i < geom.size(); i++)
						{	 //if the node is out of the fluid domain and its nodes are not part of a divided element
							 //control if it has some fixities otherwise v = 0 and p = 0.
// 							 if( geom[i].FastGetSolutionStepValue(AUX_INDEX) != 1)
							 if(geom[i].FastGetSolutionStepValue(DISTANCE) > extrapolation_distance)
							 {
								  geom[i].GetSolutionStepValue(PRESSURE) = 0.0;
								  if((geom[i].GetDof(VELOCITY_X)).IsFixed() == false )
									   geom[i].GetSolutionStepValue(VELOCITY_X) = 0.0;
// KRATOS_WATCH(geom[i].Id());}
// 								  else
// 								  {  KRATOS_WATCH("ELEMENTI CON BC IN VEL_X");
// 								      KRATOS_WATCH(iel->Id())}
								  if((geom[i].GetDof(VELOCITY_Y)).IsFixed() == false  )
								  	   geom[i].GetSolutionStepValue(VELOCITY_Y) = 0.0;
// 								  else
// 								  {  KRATOS_WATCH("ELEMENTI CON BC IN VEL_Y");
// 								      KRATOS_WATCH(iel->Id())}
								  if((geom[i].GetDof(VELOCITY_Z)).IsFixed() == false )
								  	   geom[i].GetSolutionStepValue(VELOCITY_Z) = 0.0;	

							 }
						}
					}
				}
			KRATOS_CATCH("")
			}

		private:
			ModelPart& mr_model_part;
				
	};
} //namespace Kratos

#endif //KRATOS_ELEMBASED_BC_UTILITIES_H_INCLUDED defined


