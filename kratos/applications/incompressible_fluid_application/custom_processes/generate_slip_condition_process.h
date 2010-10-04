/*
==============================================================================
KratosIncompressibleFluidApplication 
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
//   Date:                $Date: 2007-03-06 10:30:32 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_GENERATE_SLIP_CONDITION_PROCESS_H_INCLUDED )
#define  KRATOS_GENERATE_SLIP_CONDITION_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// External includes 


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/model_part.h"
#include "geometries/point_3d.h"
#include "processes/find_conditions_neighbours_process.h"
#include "incompressible_fluid_application.h"

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
	//assign a size to each node depending on the distance fof the neighbouring nodes
	*/

	class GenerateSlipConditionProcess 
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of SubscaleEstimatorProcess
		KRATOS_CLASS_POINTER_DEFINITION(GenerateSlipConditionProcess);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		GenerateSlipConditionProcess(ModelPart& model_part,int domain_size)
			: mr_model_part(model_part),mdomain_size(domain_size)
		{
		}

		/// Destructor.
		virtual ~GenerateSlipConditionProcess()
		{
		}


		///@}
		///@name Operators 
		///@{



		///@}
		///@name Operations
		///@{

		//generate a list of new nodes refining accordingly to the REFINE_FLAG
		void Execute()
		{
			KRATOS_TRY

			mslip_nodes.clear();
			
			for(ModelPart::NodeIterator in = mr_model_part.NodesBegin() ; 
				in != mr_model_part.NodesEnd() ; ++in)
			{
				in->GetValue(IS_STRUCTURE) = 0;
				array_1d<double,3>& n = in->GetValue(NORMAL);
				noalias(n) = ZeroVector(3);
				
			}

			 //calculate area normals face-by-face
			array_1d<double, 3 > area_normal;
			//2D case
			if (mdomain_size == 2)
			{
			    for (ModelPart::ConditionsContainerType::iterator cond_it = mr_model_part.ConditionsBegin(); cond_it != mr_model_part.ConditionsEnd(); cond_it++)
				CalculateNormal2D(cond_it, area_normal);
			}//3D case
			else if (mdomain_size == 3)
			{
			    //help vectors for cross product
			    array_1d<double, 3 > v1;
			    array_1d<double, 3 > v2;
			    for (ModelPart::ConditionsContainerType::iterator cond_it = mr_model_part.ConditionsBegin(); cond_it != mr_model_part.ConditionsEnd(); cond_it++)
				CalculateNormal3D(cond_it, area_normal, v1, v2);
			}

			//loop over all faces
			const double node_factor = 1.0 / double(mdomain_size);
			for (ModelPart::ConditionsContainerType::iterator cond_it = mr_model_part.ConditionsBegin(); cond_it != mr_model_part.ConditionsEnd(); cond_it++)
			{
			    //get geometry data of the face
			    Geometry<Node < 3 > >& face_geometry = cond_it->GetGeometry();

			    //reference for area normal of the face
			    array_1d<double, 3 > & face_normal = cond_it->GetValue(NORMAL);

			    //slip condition
			    if (cond_it->GetValue(IS_STRUCTURE) == true)
				for (int i = 0; i < mdomain_size; i++)
				{
				    face_geometry[i].GetValue(IS_STRUCTURE) = 1;
				    face_geometry[i].GetValue(NORMAL) += node_factor*face_normal;
				    mslip_nodes.push_back(Node<3>::Pointer(face_geometry(i)));
				}
			}
			
 			mslip_nodes.Unique();
			
			
			
			
			
			
			//****************************************************************
			/*			//now generate the list of external edges of the model. On such edges we will
			//impose that the velocity is parallel to the edges
			for(ModelPart::NodeIterator in = mslip_nodes.begin(); in != mslip_nodes.end() ; ++in)
			{
				in->GetValue(IS_VISITED)=0.0;
			}
			
			//find the conditions neighbours
			(FindConditionsNeighboursProcess(mr_model_part, mdomain_size, 10)).Execute();
			
			if(mdomain_size==3)
			{
			      //detect edges whose faces make an angle larger than 45
			      for (ModelPart::ConditionsContainerType::iterator cond_it = mr_model_part.ConditionsBegin(); cond_it != mr_model_part.ConditionsEnd(); cond_it++)
			      {
				  //get geometry data of the face
				  Geometry<Node < 3 > >& face_geometry = cond_it->GetGeometry();

				  //reference for area normal of the face
				  const array_1d<double, 3 > & face_normal = cond_it->GetValue(NORMAL);
				  double An = norm_2(face_normal);

				  unsigned int current_id = cond_it->Id();

				  //slip condition
				  if (cond_it->GetValue(IS_STRUCTURE) == 1.0) //this is a slip face --> now look for its neighbours
				  {
				      const WeakPointerVector<Condition>& neighb = cond_it->GetValue(NEIGHBOUR_CONDITIONS);

				      //check for neighbour zero
				      if(neighb[0].Id() != current_id) //check if the neighbour exists
					  CornerDectectionHelper(face_geometry,face_normal,An,neighb,1,2, 0);

				      //check for neighbour one
				      if(neighb[0].Id() != current_id) //check if the neighbour exists
					  CornerDectectionHelper(face_geometry,face_normal,An,neighb,2,0, 1 );

				      //check for neighbour two
				      if(neighb[0].Id() != current_id) //check if the neighbour exists
					  CornerDectectionHelper(face_geometry,face_normal,An,neighb,0,1, 2 );

				  }
			      }
			      
			      mConflictiveNodes.clear();
// 			      mEdgeDirection.clear();
			      
			      //loop on nodes and create the list of conflicting ones
			      for(ModelPart::NodeIterator in = mslip_nodes.begin(); in != mslip_nodes.end() ; ++in)
			      {
				      if(in->GetValue(IS_VISITED)==1.0 && AllVelocitiesFixed(in)==false &&in->IsFixed(PRESSURE)==false) //node is potentially conflictive
				      {
					  WeakPointerVector<Node<3> > ConflictiveNeigh;
					  WeakPointerVector<Node<3> >& neighb = in->GetValue(NEIGHBOUR_NODES);
					  
					  unsigned int n_conflictive_neighbours = 0;
					  for(WeakPointerVector<Node<3> >::iterator nit = neighb.begin(); nit!=neighb.end(); nit++)
					  {
					      if(nit->GetValue(IS_VISITED) == 1.0)
					      {
						n_conflictive_neighbours++;
						ConflictiveNeigh.push_back( (*nit.base()).lock() );
					      }
					  }
					  
					  //now decide what to do
					  if(n_conflictive_neighbours == 0)
					  {
					      KRATOS_ERROR(std::logic_error, "it is impossible that a conflictive node does not have any conflictive neighbour","");
					  }
					  else if(n_conflictive_neighbours == 1)
					  {
					      mConflictiveNodes.push_back( *(in.base()) );
					      array_1d<double,3> dir = ConflictiveNeigh[0].Coordinates() - in->Coordinates();
					      dir /= norm_2(dir);
					
					      //make it orthogonal to the normal
					      const array_1d<double,3> n = in->GetValue(NORMAL);
					      double A2 = inner_prod(n,n);
					      double scalar_prod = inner_prod(n,dir);
					      noalias(dir) -= scalar_prod/A2 * n;
					      dir /= norm_2(dir);
					      
					      in->GetValue(AUX_VECTOR) = dir;

					      
// 					      noalias(in->GetValue(AUX_VECTOR)) = ZeroVector(3);	
// 					      in->Fix(VELOCITY_X);
// 					      in->Fix(VELOCITY_Y);
// 					      in->Fix(VELOCITY_Z);
// 					      in->Fix(FRACT_VEL_X);
// 					      in->Fix(FRACT_VEL_Y);
// 					      in->Fix(FRACT_VEL_Z);
					  }
					  else if(n_conflictive_neighbours == 2)
					  {
// 					    noalias(in->GetValue(AUX_VECTOR)) = ZeroVector(3);	
// 					    in->Fix(VELOCITY_X);
// 					      in->Fix(VELOCITY_Y);
// 					      in->Fix(VELOCITY_Z);
// 					      in->Fix(FRACT_VEL_X);
// 					      in->Fix(FRACT_VEL_Y);
// 					      in->Fix(FRACT_VEL_Z);
					      
					      mConflictiveNodes.push_back( *(in.base()) );
 					      array_1d<double,3> dir = ConflictiveNeigh[1].Coordinates() - ConflictiveNeigh[0].Coordinates();
					      dir /= norm_2(dir);

					      
					      
					      //make it orthogonal to the normal
					      const array_1d<double,3> n = in->GetValue(NORMAL);
					      double A2 = inner_prod(n,n);
					      double scalar_prod = inner_prod(n,dir);
					      noalias(dir) -= scalar_prod/A2 * n;
					      dir /= norm_2(dir);
					      
					      in->GetValue(AUX_VECTOR) = dir;
					  }
					  else if(n_conflictive_neighbours > 2)
					  {
					      noalias(in->GetValue(AUX_VECTOR)) = ZeroVector(3);	
					      in->Fix(VELOCITY_X);
					      in->Fix(VELOCITY_Y);
					      in->Fix(VELOCITY_Z);
					      in->Fix(FRACT_VEL_X);
					      in->Fix(FRACT_VEL_Y);
					      in->Fix(FRACT_VEL_Z);
					  }
				      }
			      }
			}
			
			//here we generate the slip BCs
			//generate new point conditions on all of the nodes
			Properties::Pointer properties = mr_model_part.GetMesh().pGetProperties(0);
			Condition const& rReferenceBoundaryCondition =  KratosComponents<Condition>::Get("NoSlipFractStep");
			for(ModelPart::NodesContainerType::iterator it = mslip_nodes.begin(); it != mslip_nodes.end(); it++)
			{
			  
			  //mark the node
			  it->GetValue(IS_STRUCTURE) = 1.0;
			  
// 			  unsigned int id = (mr_model_part.Conditions().end()-1)->Id() + 1;
			  
// 			  Point3D<Node < 3 > > geom(  *(it.base())  );
			    
// 			  Condition::Pointer p_cond = rReferenceBoundaryCondition.Create(id, geom, properties);
// 			  (mr_model_part.Conditions()).push_back(p_cond);
			}
			
	*/		
			

			KRATOS_CATCH("")
		}

		void SetNormalVelocityToZero(Variable< array_1d<double,3> >& rVar)
		{
			KRATOS_TRY
		
			for(ModelPart::NodeIterator in = mslip_nodes.begin() ; 
				in != mslip_nodes.end() ; ++in)
			{
// 				const double press = in->FastGetSolutionStepValue(PRESSURE);
				const array_1d<double,3>& n = in->GetValue(NORMAL);
				double A2 = inner_prod(n,n);
				
				array_1d<double,3>& vel = in->FastGetSolutionStepValue(rVar);
				double scalar_prod = inner_prod(n,vel);
				scalar_prod/=A2;		
				noalias(vel) -= scalar_prod*n;
				
// 				array_1d<double,3>& fv = in->FastGetSolutionStepValue(FRACT_VEL);
// 				scalar_prod = inner_prod(n,fv);
// 				scalar_prod/=A2;		
// 				noalias(fv) -= scalar_prod*n;
				
			}

			KRATOS_CATCH("")
		}
		
		void ApplyEdgeConstraints(Variable< array_1d<double,3> >& rVar)
		{
			KRATOS_TRY
			
// 			unsigned int counter = 0;
// 		
// 			for(ModelPart::NodeIterator in = mConflictiveNodes.begin() ; 
// 				in != mConflictiveNodes.end() ; ++in)
// 			{
// 				const array_1d<double,3>& d = in->GetValue(AUX_VECTOR);	
// 				
// 				array_1d<double,3>& vel = in->FastGetSolutionStepValue(rVar);					
// 				double beta = inner_prod(d,vel);				
// 				noalias(vel) = beta*d;							
// 				
// 				counter++;
// 			}

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
			return "GenerateSlipConditionProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "GenerateSlipConditionProcess";
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
		int mdomain_size;
		ModelPart::NodesContainerType mslip_nodes;
		ModelPart::NodesContainerType mConflictiveNodes;
// 		std::vector< array_1d<double,3> > mEdgeDirection;

		///@} 
		///@name Private Operators
		///@{ 
		//***********************************************************
		//functions to calculate area normals for boundary conditions
		void CalculateNormal2D(ModelPart::ConditionsContainerType::iterator cond_it, array_1d<double, 3 > & area_normal)
		{
		    Geometry<Node < 3 > >& face_geometry = (cond_it)->GetGeometry();

		    area_normal[0] = face_geometry[1].Y() - face_geometry[0].Y();
		    area_normal[1] = -(face_geometry[1].X() - face_geometry[0].X());
		    area_normal[2] = 0.00;

		    noalias((cond_it)->GetValue(NORMAL)) = area_normal;
		}

		void CalculateNormal3D(ModelPart::ConditionsContainerType::iterator cond_it, array_1d<double, 3 > & area_normal, array_1d<double, 3 > & v1, array_1d<double, 3 > & v2)
		{
		    Geometry<Node < 3 > >& face_geometry = (cond_it)->GetGeometry();

		    v1[0] = face_geometry[1].X() - face_geometry[0].X();
		    v1[1] = face_geometry[1].Y() - face_geometry[0].Y();
		    v1[2] = face_geometry[1].Z() - face_geometry[0].Z();

		    v2[0] = face_geometry[2].X() - face_geometry[0].X();
		    v2[1] = face_geometry[2].Y() - face_geometry[0].Y();
		    v2[2] = face_geometry[2].Z() - face_geometry[0].Z();

		    MathUtils<double>::CrossProduct(area_normal, v1, v2);
		    area_normal *= -0.5;

		    noalias((cond_it)->GetValue(NORMAL)) = area_normal;
		}
		
		
		//**************************************
		void CornerDectectionHelper(Geometry< Node<3> >& face_geometry,
					    const array_1d<double,3>& face_normal,
					    const double An,
					    const WeakPointerVector<Condition>& neighb,
					    const unsigned int i1,
					    const unsigned int i2,
					    const unsigned int neighb_index
					    )
		{
		    //double acceptable_angle = 45.0/180.0*3.1415926; //angles of less than 45 deg will be accepted
		    double acceptable_cos = 0.707;  // cos(acceptable_angle);

// 		    if(face_geometry[i1].Id() < face_geometry[i2].Id()) //we do this to add the face ones
// 		    {
			const array_1d<double, 3 > & neighb_normal = neighb[neighb_index].GetValue(NORMAL);
			double neighb_An = norm_2(neighb_normal);

			double cos_normal = 1.0/(An*neighb_An) * inner_prod(face_normal,neighb_normal);

			//if the angle is too big between the two normals then the edge in the middle is a corner
			if(cos_normal < acceptable_cos)
			{
			  face_geometry[i1].GetValue(IS_VISITED) = 1.0;
			  face_geometry[i2].GetValue(IS_VISITED) = 1.0;
			}
// 		    }


		}
		
		bool AllVelocitiesFixed(ModelPart::NodesContainerType::iterator& it)
		{
		    if(it->IsFixed(VELOCITY_X)==true && it->IsFixed(VELOCITY_Y)==true && it->IsFixed(VELOCITY_Z)==true)
		      return true;
		    else
		      return false;
		}


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



		///@}    

	}; // Class SubscaleEstimatorProcess 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		GenerateSlipConditionProcess& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const GenerateSlipConditionProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_GENERATE_SLIP_CONDITION_PROCESS_H_INCLUDED  defined 


