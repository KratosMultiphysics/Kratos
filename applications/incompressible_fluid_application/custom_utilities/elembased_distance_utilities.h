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


#if !defined(KRATOS_ELEMBASED_DISTANCE_UTILITIES_H_INCLUDED)
#define  KRATOS_ELEMBASED_DISTANCE_UTILITIES_H_INCLUDED


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
					
			
	class ElemBasedDistanceUtilities
	{
		public:
			//constructor and destructor
			ElemBasedDistanceUtilities( ModelPart& mr_model_part
						   )
			: mr_model_part(mr_model_part){};
			~ElemBasedDistanceUtilities(){};
			

			void IdentifyFreeSurface()
			{
				KRATOS_TRY
				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();inode++)	
				{
					inode->GetValue(IS_VISITED) = 0;
				}
				
				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();
								inode++)	
				{
					
					if (inode->FastGetSolutionStepValue(DISTANCE) > 0.0)
					{	
						WeakPointerVector< Node<3> >& neighb_nodes = inode->GetValue(NEIGHBOUR_NODES); 
						for( WeakPointerVector< Node<3> >::iterator j =	neighb_nodes.begin(); j != neighb_nodes.end(); j++) 
						{ 
							if(j->FastGetSolutionStepValue(DISTANCE) < 0.0)
							{

								inode->GetValue(IS_VISITED) = 1;
// 								KRATOS_WATCH("ENTRAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
// 								KRATOS_WATCH(inode->Id())
							}
						} 
					}
					if (inode->FastGetSolutionStepValue(DISTANCE) == 0.0)
					{	inode->GetValue(IS_VISITED) = 1; 
// 						KRATOS_WATCH(inode->Id())
					}
				}
				KRATOS_CATCH("")	
			}


			void ChangeSignToDistance()
			{
				KRATOS_TRY	

				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();
								inode++)	
				{
					double dist = inode->FastGetSolutionStepValue(DISTANCE);
					inode->FastGetSolutionStepValue(DISTANCE) = -dist;
// 					KRATOS_WATCH(dist)
				}

				KRATOS_CATCH("")
			}

			void MarkNodesByDistance(double min, double max )
			{
				KRATOS_TRY	

				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();
								inode++)	
				{
					double dist = inode->FastGetSolutionStepValue(DISTANCE);
					if(dist > min && dist < max)
						inode->GetValue(IS_VISITED) = 1;
					else
						inode->GetValue(IS_VISITED) = 0;
				}

				KRATOS_CATCH("")
			}

			void SaveScalarVariableToOldStep(Variable<double>& rVar)
			{
				KRATOS_TRY	

				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();
								inode++)	
				{
					inode->FastGetSolutionStepValue(rVar,1) = inode->FastGetSolutionStepValue(rVar);
				}

				KRATOS_CATCH("")
			}

			void MarkExternalAndMixedNodes( )
			{
				KRATOS_TRY	

				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();
								inode++)	
				{
					inode->GetValue(IS_VISITED) = 0;
				}

				//detect the nodes outside the fluid surface
				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();
								inode++)	
				{
					if( inode->FastGetSolutionStepValue(DISTANCE) >= 0.0) //candidates are only the ones outside the fluid domain and the boundary
					{
						inode->GetValue(IS_VISITED) = 1;
					}
				}
				KRATOS_CATCH("")
			}

			void MarkInternalAndMixedNodes( )
			{
				KRATOS_TRY	

				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();
								inode++)	
				{
					inode->GetValue(IS_VISITED) = 0;
				}

				//detect the nodes inside the fluid surface
				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();
								inode++)	
				{
					if( inode->FastGetSolutionStepValue(DISTANCE) <= 0.0) //candidates are the ones inside the fluid domain and the boundary
					{
						inode->GetValue(IS_VISITED) = 1; 
					}
				}
				KRATOS_CATCH("")
			}
			
		
		


		private:
			ModelPart& mr_model_part;
			
			
	};


		// 			//*******************************
// 			//function to free dynamic memory
// 			void Clear()
// 			{
// 			KRATOS_TRY
// 			KRATOS_CATCH("")
// 			}
			//generate a model part with all of the elements and nodes between min_dist and max_dist
} //namespace Kratos

#endif //KRATOS_ELEMBASED_DISTANCE_UTILITIES_H_INCLUDED defined


