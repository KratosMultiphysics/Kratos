/*
==============================================================================
KratosULFApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Pawel Ryzhakov
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
//   Date:                $Date: 2009-01-15 14:50:34 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_EMBEDDED_LOCATOR_PROCESS_INCLUDED )
#define  KRATOS_EMBEDDED_LOCATOR_PROCESS_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <algorithm>
#include <boost/timer.hpp>

// External includes 


// Project includes
#include "meshing_application.h"
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "external_includes/gid_mesh_library.h"


namespace Kratos
{

	///@name Kratos Globals
	///@{ 

	///@} 
	///@name Type Definitions
	///@{ 


	///@} 
	///@name  Enum's
	///@{

	///@}
	///@name  Functions 
	///@{

	///@}
	///@name Kratos Classes
	///@{

	/// Short class definition.
	/** Detail class definition.
		calculate the nodal H for all the nodes depending on the min distance
		of the neighbouring nodes.

		lonely nodes are given the average value of the H
	*/

	class EmbeddedMeshLocatorProcess 
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of EmbeddedMeshLocatorProcess
		KRATOS_CLASS_POINTER_DEFINITION(EmbeddedMeshLocatorProcess);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		EmbeddedMeshLocatorProcess(ModelPart& model_part)
			: mr_model_part(model_part)
		{
		}

		/// Destructor.
		virtual ~EmbeddedMeshLocatorProcess()
		{
		}
  

		///@}
		///@name Operators 
		///@{

		void operator()()
		{
			Execute();
		}
  

		///@}
		///@name Operations
		///@{

		//this function calculates distances of the nodes of the given model part mr_model_part from the surface mesh (embedded_model_part) and stores it nodally
		//this permits to identify the intersected elements
		void Locate(ModelPart& embedded_model_part)
		{
			KRATOS_TRY

			ModelPart::NodesContainerType& embeddedNodes = embedded_model_part.Nodes();
			ModelPart::NodesContainerType& volumeNodes = mr_model_part.Nodes();
			//reset the COUNTER 
			for(ModelPart::NodesContainerType::iterator in = embeddedNodes.begin(); in!=embeddedNodes.end(); in++)
			{
				in->FastGetSolutionStepValue(COUNTER)=0;
			}
			
			for(ModelPart::NodesContainerType::iterator in = volumeNodes.begin(); in!=volumeNodes.end(); in++)
			{
				in->FastGetSolutionStepValue(COUNTER)=0;
			}
			
			
			//the variable COUNTER stores for the skin nodes adn volume nodes store their positions in the Gid lists of surface and volume nodes respectively.
			
			///////////////////////////////////////////////////////////////////////////////////////////////////////////
			//count the nummber of the "skin" nodes
			unsigned int n_skin_nodes=embedded_model_part.Nodes().size();
			
			//each node is represented by its three coordinates
			double* skin_nodes = NULL;
			skin_nodes = new double[n_skin_nodes*3];				
			
			unsigned int count=0;
			unsigned int base=0;
			
			for(ModelPart::NodesContainerType::iterator in = embeddedNodes.begin(); in!=embeddedNodes.end(); in++)
			{
			//each node has three coords
			base=count*3;
			skin_nodes[base]=in->X();
			skin_nodes[base+1]=in->Y();
			skin_nodes[base+2]=in->Z();
			in->FastGetSolutionStepValue(COUNTER)=count;
			count++;
			}


			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			//count the number of the "skin" faces (must be triangles) 
			//check if these are called CONDITIONS in kratos!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			unsigned int n_skin_faces=embedded_model_part.Conditions().size();
			//KRATOS_WATCH(embedded_model_part.Conditions().size())
			//each face is represented by "numbers" of its three nodes (vertices)
			int* skin_faces = NULL;
			skin_faces = new int[n_skin_faces*3];
			//we reuse the counter			
			count=0;
			ModelPart::ConditionsContainerType& embeddedConds = embedded_model_part.Conditions();
			
			
			for(ModelPart::ConditionsContainerType::iterator ic = embeddedConds.begin() ; 
					ic != embeddedConds.end() ; ++ic)
			{
			//each face has three nodes
			base=count*3;
			//connectivities - according to their position in the skin_nodes list 
			//(three consecutive coords are corresponding to one node, i.e. one entry in the skin_faces list)
			unsigned int position_0=ic->GetGeometry()[0].FastGetSolutionStepValue(COUNTER);
			unsigned int position_1=ic->GetGeometry()[1].FastGetSolutionStepValue(COUNTER);
			unsigned int position_2=ic->GetGeometry()[2].FastGetSolutionStepValue(COUNTER);
			skin_faces[base]=position_0;
			skin_faces[base+1]=position_1;
			skin_faces[base+2]=position_2;
			count++;
			}
			

			
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			//count the nummber of the volume nodes
			unsigned int n_vol_nodes=mr_model_part.Nodes().size();
			
			double* vol_nodes = NULL;
			vol_nodes = new double[n_vol_nodes*3];
			  
			count=0;

			double search_radius=0.0;
			double temp=0.0;
			//and fill in the list of the volume nodes
			for(ModelPart::NodesContainerType::iterator in = volumeNodes.begin(); in!=volumeNodes.end(); in++)
			{
			//each node has three coords
			base=count*3;
			vol_nodes[base]=in->X();
			vol_nodes[base+1]=in->Y(); 
			vol_nodes[base+2]=in->Z();
			in->FastGetSolutionStepValue(COUNTER)=count;
			count++;			
			} 
			//compute the largest edge of an element in the volume mesh
			double h_max=0.0;		
			double h_real;
			for(ModelPart::NodesContainerType::iterator in = volumeNodes.begin(); in!=volumeNodes.end(); in++)
			{
				if((in->GetValue(NEIGHBOUR_NODES)).size() != 0)
				{
					double xc = in->X();
					double yc = in->Y();
					double zc = in->Z();

					double h = 10000000.0;
					for( WeakPointerVector< Node<3> >::iterator i = in->GetValue(NEIGHBOUR_NODES).begin();
									i != in->GetValue(NEIGHBOUR_NODES).end(); i++)
					{
						double x = i->X();
						double y = i->Y();
						double z = i->Z();
						double l = (x-xc)*(x-xc);
						l += (y-yc)*(y-yc);
						l += (z-zc)*(z-zc);

						//if(l>h) h = l;
						if(l<h) h = l;
					}

					h = sqrt(h);	

				if (h>h_max)
					h_max=h;			

				} 
				
			}
			
			KRATOS_WATCH(h_max)
			//the search radius should be at least 4 times the element size
			search_radius=4.0*h_max;
			KRATOS_WATCH("=========================EMBEDDED=INPUT=INFO================================================================")	
			KRATOS_WATCH(n_skin_faces)
			KRATOS_WATCH(n_skin_nodes)
			KRATOS_WATCH(n_vol_nodes)
			KRATOS_WATCH(search_radius)
			
			if (search_radius<=0.0)
				KRATOS_ERROR(std::logic_error,"error: YOUR SEARCH RADIUS FOR FINDING INTERSECTIONS IS ZERO OR NEGATIVE!!!!!! CHECK IF NODAL H WAS COMPUTED","");
			unsigned int TDim=3;
			void* hinput=GetGidInputHandleWithBoundaryMesh(TDim, n_skin_nodes, skin_nodes,n_skin_faces,skin_faces, 3);
						 //void* hinput=GetGidInputHandleWithBoundaryMesh(3,8,coord,12,conec,3);
			AddNodesInMesh(hinput, n_vol_nodes, vol_nodes);

			//////////////////TIMER///////////////////////////////
			boost::timer distance_calc_time;

			void* hout=GetGidOutputHandle();

			
			 
			//GiDML_GetNodesDistanceInTetrahedraMesh(hinput, hout); 1 - cuadrado, 0 - con signo
			GiDML_GetNodesDistanceRadiusInTetrahedraMesh(hinput, hout, search_radius, 0);
			std::cout << "Distance from embedded object: calculation time = " << distance_calc_time.elapsed() << std::endl;
			KRATOS_WATCH("============================================================================================================")
			
			double* distances=NULL;
			//compute the distances using GiD library
			distances=new double[n_vol_nodes];
			distances=GetAttributesOfNodes(hout);
			
			//store the distances in the nodes of KRATOS
			//count=0;
			for(ModelPart::NodesContainerType::iterator in = volumeNodes.begin(); in!=volumeNodes.end(); in++)
			{
			//each node has three coords
			count=in->FastGetSolutionStepValue(COUNTER);
			in->FastGetSolutionStepValue(DISTANCE)=distances[count]; 
			}
			
			//clearing memory
			delete [] vol_nodes;  // When done, free memory pointed to by a.
			delete [] skin_nodes;
			delete [] skin_faces;
			delete [] distances;
			vol_nodes=NULL;
			skin_nodes=NULL;
			skin_faces=NULL;
			distances=NULL;

			//now we assign to the nodes that lie outside of the search radius the correct sign.
			temp=0;
			double max_dist=0;
			for(ModelPart::NodesContainerType::iterator in = volumeNodes.begin(); in!=volumeNodes.end(); in++)
			{
				temp=in->FastGetSolutionStepValue(DISTANCE);
				if (temp>max_dist)
					max_dist=temp;
			}
			KRATOS_WATCH(max_dist)
			const double default_distance=max_dist;//search_radius+1.0;
			

			for(ModelPart::NodesContainerType::iterator in = volumeNodes.begin(); in!=volumeNodes.end(); in++)
			{
			double distance;
			double distance_temp;
			distance=in->FastGetSolutionStepValue(DISTANCE);
			if (distance==default_distance)
				{				
				//KRATOS_WATCH(distance)	
				WeakPointerVector< Node<3> >& neighb = in->GetValue(NEIGHBOUR_NODES);
				bool inner_elem=false;
				//int i=0;
				//while (neighb[i]>0 && i<neighb.size())
				//{
				//i++;
				//}
				//if (i<neighb.size()) //if there was a node with a negative distance in the patch, set this distance to all the nodes (except for it
				for (unsigned int i=0;i<neighb.size();i++)
					{
					distance_temp=neighb[i].FastGetSolutionStepValue(DISTANCE);
					if (distance_temp<0.0)
						{
						inner_elem=true;						
						}
					}
				if (inner_elem==true)
					{
						for (unsigned int i=0;i<neighb.size();i++)
						{
						if (neighb[i].FastGetSolutionStepValue(DISTANCE)==default_distance)
							neighb[i].FastGetSolutionStepValue(DISTANCE)=-default_distance;
						}
					}
					

				}
			}
			bool uncertain=true;
			while (uncertain==true)
			{
			uncertain=false;
			for(ModelPart::NodesContainerType::iterator in = volumeNodes.begin(); in!=volumeNodes.end(); in++)
			{
			double distance;
			double distance_temp;
			distance=in->FastGetSolutionStepValue(DISTANCE);
			if (distance==-default_distance)	
				{
				WeakPointerVector< Node<3> >& neighb= in->GetValue(NEIGHBOUR_NODES);				
				for (unsigned int i=0;i<neighb.size();i++)
					{
					distance_temp=neighb[i].FastGetSolutionStepValue(DISTANCE);
					if (distance_temp==default_distance)
						{
						neighb[i].FastGetSolutionStepValue(DISTANCE)=-default_distance;
						KRATOS_WATCH("Second round of negative nodes")
						uncertain=true;
						}

					}

				}

			}
			}
			
					
			KRATOS_CATCH("")
		}


		///@}
		///@name Access
		///@{ 


		///@}
		///@name Inquiry
		///@{


		///@}      
		///@name Input and output
		///@{

		/// Turn back information as a string.
		virtual std::string Info() const
		{
			return "EmbeddedMeshLocatorProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "EmbeddedMeshLocatorProcess";
		}

		/// Print object's data.
		virtual void PrintData(std::ostream& rOStream) const
		{
		}


		///@}      
		///@name Friends
		///@{


		///@}

	protected:
		///@name Protected static Member Variables 
		///@{ 


		///@} 
		///@name Protected member Variables 
		///@{ 


		///@} 
		///@name Protected Operators
		///@{ 


		///@} 
		///@name Protected Operations
		///@{ 


		///@} 
		///@name Protected  Access 
		///@{ 


		///@}      
		///@name Protected Inquiry 
		///@{ 


		///@}    
		///@name Protected LifeCycle 
		///@{ 


		///@}

	private:
		///@name Static Member Variables 
		///@{ 


		///@} 
		///@name Member Variables 
		///@{ 
		ModelPart& mr_model_part;
		


		///@} 
		///@name Private Operators
		///@{ 

		///@} 
		///@name Private Operations
		///@{ 


		///@} 
		///@name Private  Access 
		///@{ 


		///@}    
		///@name Private Inquiry 
		///@{ 


		///@}    
		///@name Un accessible methods 
		///@{ 

		/// Assignment operator.
		EmbeddedMeshLocatorProcess& operator=(EmbeddedMeshLocatorProcess const& rOther){};

		/// Copy constructor.
		//EmbeddedMeshLocatorProcess(EmbeddedMeshLocatorProcess const& rOther);


		///@}    

	}; // Class EmbeddedMeshLocatorProcess 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		EmbeddedMeshLocatorProcess& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const EmbeddedMeshLocatorProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_EMBEDDED_LOCATOR_PROCESS_INCLUDED  defined 


