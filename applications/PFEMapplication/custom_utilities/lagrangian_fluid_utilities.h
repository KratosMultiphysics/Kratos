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
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:31 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_LARANGIAN_FLUID_UTILITIES_INCLUDED )
#define  KRATOS_LARANGIAN_FLUID_UTILITIES_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// External includes 


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "custom_utilities/geometry_utilities.h"
#include "geometries/tetrahedra_3d_4.h"
#include "incompressible_fluid_application.h"

namespace Kratos
{
	class PfemUtils
	{
	public:


		//**********************************************************************************************
		//**********************************************************************************************
		double EstimateDeltaTime(double dt_min, double dt_max, ModelPart::NodesContainerType& rNodes)
		{
			KRATOS_TRY

				array_1d<double,3> xc, dist;
			double deltatime = 1000.0;
			array_1d<double,3> vaux;

			for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
			{
				const array_1d<double,3>& vc = in->FastGetSolutionStepValue(VELOCITY);

				if((in->GetValue(NEIGHBOUR_NODES)).size() != 0)
				{
					xc[0] = in->X();
					xc[1] = in->Y();
					xc[2] = in->Z();

					double dt_node = 0.00;
					for( WeakPointerVector< Node<3> >::iterator i = in->GetValue(NEIGHBOUR_NODES).begin();
						i != in->GetValue(NEIGHBOUR_NODES).end(); i++)
					{
						dist[0] = i->X();
						dist[1] = i->Y();
						dist[2] = i->Z();
						dist -= xc;

						noalias(vaux) =  i->FastGetSolutionStepValue(VELOCITY);
						noalias(vaux) -=  vc;
						double vnorm2 = inner_prod(vaux,vaux);

						if(vnorm2 != 0.0)
						{
							double nom = inner_prod(dist,vaux);
							double dt_side = fabs(nom/vnorm2);
							if(dt_side > dt_node) 
								dt_node = dt_side;
							else if(dt_node < deltatime)
								deltatime = dt_node;
						} 
					}
				}
			}	

			if(deltatime > dt_max) deltatime = dt_max;
			else if (deltatime < dt_min) deltatime = dt_min; 

			return deltatime;

			KRATOS_CATCH("")
		}

		//**********************************************************************************************
		//**********************************************************************************************
		void MarkOuterNodes(const array_1d<double,3>& corner1, const array_1d<double,3>& corner2,
			ModelPart::NodesContainerType& rNodes)
		{			
			KRATOS_TRY;

			//add a big number to the id of the nodes to be erased
			int n_erased = 0;
			int added_id = 100000000;
			double xmax, xmin, ymax,ymin, zmax, zmin;



			if(corner1[0] > corner2[0])
			{	xmax = corner1[0]; xmin = corner2[0]; 	}
			else
			{	xmax = corner2[0]; xmin = corner1[0]; 	}

			if(corner1[1] > corner2[1])
			{	ymax = corner1[1]; ymin = corner2[1]; 	}
			else
			{	ymax = corner2[1]; ymin = corner1[1]; 	}

			if(corner1[2] > corner2[2])
			{	zmax = corner1[2]; zmin = corner2[2]; 	}
			else
			{	zmax = corner2[2]; zmin = corner1[2]; 	}



			for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
			{
				bool erase = false;
				double& x = in->X(); double& y = in->Y(); double& z = in->Z();

				if(x<xmin || x>xmax)		erase = true;
				else if(y<ymin || y>ymax)	erase = true;
				else if(z<zmin || z>zmax)	erase = true;

				if(erase == true)
				{
					n_erased += 1;
					in->GetValue(ERASE_FLAG) = true;
				}
			}

			KRATOS_CATCH("")
		}


		//**********************************************************************************************
		//**********************************************************************************************
		//calculates the ale spatial acceleration, taking in account both the spatial and the time dependent
		//componenents acc = (v-vm)*dv/dx
		template< int TDim > void CalculateSpatialALEAcceleration(ModelPart& ThisModelPart)
		{			
			KRATOS_TRY;


			double factor = 1.00/(TDim + 1);
			const int nnodes = TDim+1;

			//reset the acceleration on the nodes 
			ModelPart::NodesContainerType& rNodes = ThisModelPart.Nodes();
			for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
			{
				noalias(in->FastGetSolutionStepValue(PRESSURE)) = 0.00;
			}

			//compute the projection (first step
			for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin(); 
				i!=ThisModelPart.ElementsEnd(); i++)
			{	
				//calculating shape functions values
				Geometry< Node<3> >& geom = i->GetGeometry();
				double Volume;
				
				GeometryUtils::CalculateGeometryData(geom,DN_DX,N,Volume);

				double el_pressure = i->GetValue(PRESSURE);

				for(int I = 0; I<nnodes; I++)
				{
					geom[I].FastGetSolutionStepValue(PRESSURE) += factor * el_pressure;

				}	
			}

			//finalize the projection of the pressure
			for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
			{
				double area = in->FastGetSolutionStepValue(NODAL_AREA);
				in->FastGetSolutionStepValue(PRESSURE) /= area;
			} 

			KRATOS_CATCH("")
		}

		//**********************************************************************************************
		//*		//**********************************************************************************************
		void CalculateNodalArea(ModelPart& ThisModelPart, int domain_size)
		{
			//set to zero the nodal area
			for(ModelPart::NodesContainerType::iterator in = ThisModelPart.NodesBegin(); 
				in!=ThisModelPart.NodesEnd(); in++)
			{
				in->FastGetSolutionStepValue(NODAL_AREA) = 0.00;
			}

			if(domain_size == 2)
			{
				double area = 0.0;
				for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin(); 
					i!=ThisModelPart.ElementsEnd(); i++)
				{	
					//calculating shape functions values
					Geometry< Node<3> >& geom = i->GetGeometry();
					
					area = GeometryUtils::CalculateVolume2D(geom);
					area *= 0.333333333333333333333333333;
					
					geom[0].FastGetSolutionStepValue(NODAL_AREA) += area;
					geom[1].FastGetSolutionStepValue(NODAL_AREA) += area;
					geom[2].FastGetSolutionStepValue(NODAL_AREA) += area;
				}
			}
			else if(domain_size == 3)
			{
				for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin(); 
					i!=ThisModelPart.ElementsEnd(); i++)
				{	
					double vol;
					//calculating shape functions values
					Geometry< Node<3> >& geom = i->GetGeometry();
					
					vol = GeometryUtils::CalculateVolume3D(geom);
					vol *= 0.25;
					
					geom[0].FastGetSolutionStepValue(NODAL_AREA) += vol;
					geom[1].FastGetSolutionStepValue(NODAL_AREA) += vol;
					geom[2].FastGetSolutionStepValue(NODAL_AREA) += vol;
					geom[3].FastGetSolutionStepValue(NODAL_AREA) += vol;
				}
			}
			//r


			
		}


		//**********************************************************************************************
		//**********************************************************************************************
		//imposes the velocity that corresponds to a 
		void MoveLonelyNodes(ModelPart& ThisModelPart)
		{			
			KRATOS_TRY;

			double Dt = ThisModelPart.GetProcessInfo()[DELTA_TIME];
			double accfactor = 0.5*Dt*Dt;

			array_1d<double,3> DeltaDisp;

			const array_1d<double,3> body_force = ThisModelPart.ElementsBegin()->GetProperties()[BODY_FORCE];
			for(ModelPart::NodeIterator i = ThisModelPart.NodesBegin() ; 
				i != ThisModelPart.NodesEnd() ; ++i)
			{
				if(
					(i)->FastGetSolutionStepValue(IS_STRUCTURE) == 0 && //if it is not a wall node
					(i)->GetValue(NEIGHBOUR_ELEMENTS).size() == 0 //and it is lonely
					)
				{
					//set to zero the pressure
					(i)->FastGetSolutionStepValue(PRESSURE) = 0;

					const array_1d<double,3>& old_vel = (i)->FastGetSolutionStepValue(VELOCITY,1);
					array_1d<double,3>& vel = (i)->FastGetSolutionStepValue(VELOCITY);
					array_1d<double,3>& acc = (i)->FastGetSolutionStepValue(ACCELERATION);

					noalias(acc) = body_force;

					noalias(vel) = old_vel;
					noalias(vel) += Dt * acc ;

					//calculate displacements
					noalias(DeltaDisp) = Dt * vel;

					array_1d<double,3>& disp = i->FastGetSolutionStepValue(DISPLACEMENT);
					noalias(disp) = i->FastGetSolutionStepValue(DISPLACEMENT,1);
					noalias(disp) += DeltaDisp;

				}

			}

			KRATOS_CATCH("")
		}

		//**********************************************************************************************
		//**********************************************************************************************
		double CalculateVolume(ModelPart& ThisModelPart, int domain_size)
		{			
			KRATOS_TRY;

			//auxiliary vectors
			double Atot = 0.00;

			if(domain_size == 2)
			{
				for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin(); 
					i!=ThisModelPart.ElementsEnd(); i++)
				{	
					//calculating shape functions values
					Geometry< Node<3> >& geom = i->GetGeometry();
					
					Atot += GeometryUtils::CalculateVolume2D(geom);
				}
			}
			else if(domain_size == 3)
			{
				for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin(); 
					i!=ThisModelPart.ElementsEnd(); i++)
				{	
					//calculating shape functions values
					Geometry< Node<3> >& geom = i->GetGeometry();
					
					Atot += GeometryUtils::CalculateVolume3D(geom);
				}
			}
			//return the total area
			return Atot;
			KRATOS_CATCH("");
		}



	private:

		//aux vars
		static boost::numeric::ublas::bounded_matrix<double,3,3> msJ; //local jacobian
		static boost::numeric::ublas::bounded_matrix<double,3,3> msJinv; //inverse jacobian
		static array_1d<double,3> msc; //center pos
		static array_1d<double,3> ms_rhs; //center pos



	};

}  // namespace Kratos.

#endif // KRATOS_LARANGIAN_FLUID_UTILITIES_INCLUDED  defined 


