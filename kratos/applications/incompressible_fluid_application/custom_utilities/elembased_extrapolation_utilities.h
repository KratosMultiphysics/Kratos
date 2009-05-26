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


#if !defined(KRATOS_ELEMBASED_EXTRAPOLATION_UTILITIES_H_INCLUDED)
#define  KRATOS_ELEMBASED_EXTRAPOLATION_UTILITIES_H_INCLUDED

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
					
			
	class ElemBasedExtrapolationUtilities
	{
		public:
			//constructor and destructor
			ElemBasedExtrapolationUtilities(    ModelPart& mr_model_part
				    )
			:mr_model_part(mr_model_part){};
			~ElemBasedExtrapolationUtilities(){};


			//********************************
			//function to compute coefficients
			void ExtrapolateVelocities(unsigned int extrapolation_layers)
			{
			KRATOS_TRY	

				typedef Node<3> 				PointType;
				typedef PointerVector<PointType >       	PointVector;
				typedef PointVector::iterator 		PointIterator;

				//reset is visited flag
				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();
								inode++)	
				{
					inode->GetValue(IS_VISITED) = 0;
				}

//Fill layers begin			
				//generate a container with the layers to be extrapolated
				std::vector< PointVector > layers(extrapolation_layers);

				//detect the nodes inside the fluid surface
				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();
								inode++)	
				{
					if( inode->FastGetSolutionStepValue(DISTANCE) <= 0.0) //candidates are only the ones inside the fluid domain
					{
						WeakPointerVector< Node<3> >& neighb_nodes = inode->GetValue(NEIGHBOUR_NODES); 
						for( WeakPointerVector< Node<3> >::iterator i =	neighb_nodes.begin(); i != neighb_nodes.end(); i++) 
						{ 
							if(i->FastGetSolutionStepValue(DISTANCE) > 0) //add the node as free surface if one of its neighb is outside
							{
								if( inode->GetValue(IS_VISITED) == 0)
								{
									layers[0].push_back( *(inode.base() ) );	
									inode->GetValue(IS_VISITED) = 1;
								}
							}
						} 
					}
				}


				//fill the following layers by neighbour relationships
				//each layer fills the following
				for(unsigned int il = 0; il<extrapolation_layers-1; il++)
				{
					for( PointIterator iii=(layers[il]).begin(); iii!=(layers[il]).end(); iii++)
					{
						WeakPointerVector< Node<3> >& neighb_nodes = iii->GetValue(NEIGHBOUR_NODES); 
						for(WeakPointerVector< Node<3> >::iterator jjj=neighb_nodes.begin(); jjj !=neighb_nodes.end(); jjj++) //destination = origin1 + value * Minv*origin
						{ 
							if( jjj->FastGetSolutionStepValue(DISTANCE) > 0 &&
							jjj->GetValue(IS_VISITED) == 0.0 )
							{
								layers[il+1].push_back( Node<3>::Pointer( *(jjj.base() ) ) );
								jjj->GetValue(IS_VISITED) = double(il+2.0);
							}
						}
					}
				}
//Fill layers end			

				//perform extrapolation layer by layer by making an average 
				//of the neighbours of lower order
				array_1d<double,3> aux;	
				for(unsigned int il = 1; il<extrapolation_layers; il++)
				{
					for( PointIterator iii=layers[il].begin(); iii!=layers[il].end(); iii++)
					{
						//extrapolate the average velocity
						noalias(aux) = ZeroVector(3);
						double avg_number = 0.0;
						
						WeakPointerVector< Node<3> >& neighb_nodes = iii->GetValue(NEIGHBOUR_NODES); 
						for(WeakPointerVector< Node<3> >::iterator i=neighb_nodes.begin(); 	i !=neighb_nodes.end(); i++) 
						{ 
							if(i->GetValue(IS_VISITED) < il+1 && i->GetValue(IS_VISITED) > 0)
							{
								noalias(aux) += i->FastGetSolutionStepValue(VELOCITY);
								avg_number += 1.0;
							}
						} 
						
						array_1d<double,3>& Vel = iii->FastGetSolutionStepValue(VELOCITY);

						//we have to respect bc on velocity that otherwise will be deleted
						if(iii->IsFixed(VELOCITY_X) == false && avg_number != 0.0)
							 Vel[0] += aux[0]/avg_number;
						if(iii->IsFixed(VELOCITY_Y) == false && avg_number != 0.0)
							 Vel[1] += aux[1]/avg_number;
						if(iii->IsFixed(VELOCITY_Z) == false && avg_number != 0.0)
							 Vel[2] += aux[2]/avg_number;

// 						noalias( iii->FastGetSolutionStepValue(VELOCITY) ) = aux;
					}
				}

				

			KRATOS_CATCH("")
			}


		private:
			ModelPart& mr_model_part;
				
	};
} //namespace Kratos

#endif //KRATOS_ELEMBASED_EXTRAPOLATION_UTILITIES_H_INCLUDED defined


