//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2008-02-25 19:05:54 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_TETGEN_REFINER_H_INCLUDED )
#define  KRATOS_TETGEN_REFINER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <stdlib.h>
#include <timer.hpp>



#include "tetgen.h" // Defined tetgenio, tetrahedralize().

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "geometries/tetrahedra_3d_4.h"
#include "meshing_application.h"

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
	*/
	class TetGenRefine  
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of TetGenRefine
		KRATOS_CLASS_POINTER_DEFINITION(TetGenRefine);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		TetGenRefine():
		mJ(ZeroMatrix(3,3)), //local jacobian
		mJinv(ZeroMatrix(3,3)), //inverse jacobian
		mC(ZeroVector(3)), //dimension = number of nodes
		mRhs(ZeroVector(3)){} //dimension = number of nodes

		/// Destructor.
		virtual ~TetGenRefine(){}


		///@}
		///@name Operators 
		///@{


		///@}
		///@name Operations
		///@{


		//*******************************************************************************************
		//*******************************************************************************************
		void RefineMesh(
			ModelPart& ThisModelPart , 
			Element const& rReferenceElement, 
			Condition const& rReferenceBoundaryCondition,
			double alpha_param = 1.4)
		{

			KRATOS_TRY

			boost::timer auxiliary;

			//clearing elements
			tetgenio in, out;

			// All indices start from 1.
			in.firstnumber = 1;

			//reserving memory for the node list
			ModelPart::SizeType number_of_points = ThisModelPart.Nodes().size();
			in.numberofpoints = number_of_points;
			in.pointlist = new REAL[in.numberofpoints * 3];

			//reserving memory for attribute list
			in.numberofpointattributes = ThisModelPart.GetNodalSolutionStepTotalDataSize();
			typedef double* pointer_to_double_type;
			in.pointattributelist = new pointer_to_double_type[number_of_points] ;
			//reserving memory for the element list
			in.numberoftetrahedra =  ThisModelPart.Elements().size();
			in.tetrahedronlist = new int[in.numberoftetrahedra * 4];

			//reserving area for volume constraints
			in.tetrahedronvolumelist = new REAL[in.numberoftetrahedra];

			//writing the point coordinates in a vector
			ModelPart::NodesContainerType::iterator nodes_begin = ThisModelPart.NodesBegin();

			//reorder node Ids
			for(unsigned int i = 0; i<ThisModelPart.Nodes().size(); i++)
			{
				(nodes_begin + i)->Id() = i+1;
			}

			//give the coordinates to the mesher
			for(unsigned int i = 0; i<ThisModelPart.Nodes().size(); i++)
			{
				//set node coordinates
				int base = i*3;
				in.pointlist[base] = (nodes_begin + i)->X();
				in.pointlist[base+1] = (nodes_begin + i)->Y();
				in.pointlist[base+2] = (nodes_begin + i)->Z();

				//set node attributes
				in.pointattributelist[i] = (nodes_begin + i)->SolutionStepData().Data();
				//...
				//...
			}

			//provide the list of tethraedra to the mesher
			unsigned int i = 0;
			for(ModelPart::ElementsContainerType iel =  ThisModelPart.ElementsBegin(); iel!=  ThisModelPart.ElementsEnd(); iel++)
			{
				Geometry< Node<3> >& geom = iel->GetGeometry(); 
				int base = i*4;
				in.tetrahedronlist[base]   = geom[0]->Id();
				in.tetrahedronlist[base+1] = geom[1]->Id();
				in.tetrahedronlist[base+2] = geom[2]->Id();
				in.tetrahedronlist[base+3] = geom[3]->Id();

				double h =  geom[0].FastGetSolutionStepValue(NODAL_H);
				h +=  geom[1].FastGetSolutionStepValue(NODAL_H);
				h +=  geom[2].FastGetSolutionStepValue(NODAL_H);
				h +=  geom[3].FastGetSolutionStepValue(NODAL_H);
				h *= 0.25;

				//a given tethraedra of edge h has volume V = sqrt(2)/12 * h*h*h = 0.11785113 * h*h*h
				double desired_volume = 0.11785113 * h*h*h;

				//set the desired area in the list
				in.tetrahedronvolumelist[i] = desired_volume;
				
				//update counter
				i = i + 1;
			}	

			//erase elements and conditions
			ThisModelPart.Elements().clear();
			ThisModelPart.Conditions().clear();

			//regenerate the mesh (r)
			//refine the mesh using a quality criteria (Q) without splitting the boundaries (Y) and without jettisoning nodes
			tetrahedralize("rqQYJ", &in, &out); 

			double first_part_time = auxiliary.elapsed();
			std::cout << "mesh generation time = " << first_part_time << std::endl;

			//freeing unnecessary memory
			//in.deinitialize();
			//in.initialize();

			//generate the new node list
			unsigned int original_size = ThisModelPart.Nodes().size();
			if( out.numberofpoints > original_size;  )
			{
				for(unsigned int i = original_size; i < out.numberofpoints; i++)
				{
					int id = i
					int base = i*3;
					double& x = out.pointlist[base];		
					double& y = out.pointlist[base+1];
					double& z = out.pointlist[base+2];
					
					//create new node with the x,y,z corresponding and the attribute list given
					ThisModelPart.AddNode(id,x,y,z, attributelist);
				}
			}

			//generate Kratos Tetrahedra3D4
			boost::timer adding_elems;
			nodes_begin = ThisModelPart.NodesBegin();
			(ThisModelPart.Elements()).reserve(out.numberoftetrahedra);
			for(int iii = 0; iii< out.numberoftetrahedra; iii++)
			{
				int id = iii + 1;
				int base = iii * 4;
				Tetrahedra3D4<Node<3> > geom(
					*( (nodes_begin +  out.tetrahedronlist[base]-1).base() 		), 
					*( (nodes_begin +  out.tetrahedronlist[base+1]-1).base() 	), 
					*( (nodes_begin +  out.tetrahedronlist[base+2]-1).base() 	), 
					*( (nodes_begin +  out.tetrahedronlist[base+3]-1).base() 	) 
					);	

				 ThisModelPart.GetMesh().pGetProperties(1)
				Element::Pointer p_element = rReferenceElement.Create(id, geom, ThisModelPart.GetMesh().pGetProperties(prop_id));
				(ThisModelPart.Elements()).push_back(p_element);

			}
			ThisModelPart.Elements().Sort();
			std::cout << "time for adding elems" << adding_elems.elapsed() << std::endl;;

			//***********************************************************************************
			//***********************************************************************************
			//mark boundary nodes 
			//reset the boundary flag
			for(ModelPart::NodesContainerType::const_iterator in = ThisModelPart.NodesBegin(); in!=ThisModelPart.NodesEnd(); in++)
			{
				in->FastGetSolutionStepValue(IS_BOUNDARY) = 0;
			}

			//***********************************************************************************
			//***********************************************************************************
			boost::timer adding_faces;

			(ThisModelPart.Conditions()).reserve(outnew.numberoftrifaces   );

			//creating the faces
			for(ModelPart::ElementsContainerType::const_iterator iii = ThisModelPart.ElementsBegin();
				iii != ThisModelPart.ElementsEnd(); iii++)
			{
				
				int base = ( iii->Id() - 1 )*4;
				
				//create boundary faces and mark the boundary nodes 
				//each face is opposite to the corresponding node number so
				// 0 ----- 1 2 3
				// 1 ----- 0 3 2
				// 2 ----- 0 1 3
				// 3 ----- 0 2 1

				//node 1

				if( outnew.neighborlist[base] == -1)
				{
					CreateBoundaryFace(1, 2, 3, ThisModelPart,   0, *(iii.base()), properties );
				}
				if( outnew.neighborlist[base+1] == -1)
				{
					CreateBoundaryFace(0,3,2, ThisModelPart,   1, *(iii.base()), properties );
				}
				if( outnew.neighborlist[base+2] == -1)
				{
					CreateBoundaryFace(0,1,3, ThisModelPart,   2, *(iii.base()), properties );
				}
				if( outnew.neighborlist[base+3] == -1)
				{
					CreateBoundaryFace(0,2,1, ThisModelPart,   3, *(iii.base()), properties );
				}

			}
			std::cout << "time for adding faces" << adding_faces.elapsed() << std::endl;;
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
		virtual std::string Info() const{return "";}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const{}

		/// Print object's data.
		virtual void PrintData(std::ostream& rOStream) const{}


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
		 boost::numeric::ublas::bounded_matrix<double,3,3> mJ; //local jacobian
		 boost::numeric::ublas::bounded_matrix<double,3,3> mJinv; //inverse jacobian
		 array_1d<double,3> mC; //center pos
		 array_1d<double,3> mRhs; //center pos


		void CreateBoundaryFace(const int& i1, const int& i2, const int& i3, ModelPart& ThisModelPart, const int& outer_node_id, Element::Pointer origin_element, Properties::Pointer properties)
		{
			KRATOS_TRY

			Geometry<Node<3> >& geom = origin_element->GetGeometry();
			//mark the nodes as free surface
			geom[i1].FastGetSolutionStepValue(IS_BOUNDARY) = 1;
			geom[i2].FastGetSolutionStepValue(IS_BOUNDARY) = 1;
			geom[i3].FastGetSolutionStepValue(IS_BOUNDARY) = 1;

			//generate a face condition
			Condition::NodesArrayType temp;
			temp.reserve(3);
			temp.push_back(geom(i1)); 
			temp.push_back(geom(i2));
			temp.push_back(geom(i3));
			Geometry< Node<3> >::Pointer cond = Geometry< Node<3> >::Pointer(new Geometry< Node<3> >(temp) );
			int id = (origin_element->Id()-1)*4;
			Condition::Pointer p_cond = Condition::Pointer(new Condition(id, cond, properties) );

			//assigning the neighbour node
			(p_cond->GetValue(NEIGHBOUR_NODES)).clear();
			(p_cond->GetValue(NEIGHBOUR_NODES)).push_back( Node<3>::WeakPointer( geom(outer_node_id) ) );
			(p_cond->GetValue(NEIGHBOUR_ELEMENTS)).clear();
			(p_cond->GetValue(NEIGHBOUR_ELEMENTS)).push_back( Element::WeakPointer( origin_element ) );
			ThisModelPart.Conditions().push_back(p_cond);
			KRATOS_CATCH("")

		}


		//*****************************************************************************************
		//*****************************************************************************************
		void CalculateElementData(const array_1d<double,3>& c1,
			const array_1d<double,3>& c2,
			const array_1d<double,3>& c3,
			const array_1d<double,3>& c4,
			double& volume,
			array_1d<double,3>& xc,
			double& radius,
			double& hmin,
			double& hmax
			)
		{
			KRATOS_TRY
			const double& x0 = c1[0]; const double& y0 = c1[1]; const double& z0 = c1[2];
			const double& x1 = c2[0]; const double& y1 = c2[1]; const double& z1 = c2[2];
			const double& x2 = c3[0]; const double& y2 = c3[1]; const double& z2 = c3[2];
			const double& x3 = c4[0]; const double& y3 = c4[1]; const double& z3 = c4[2];

			//calculate min side lenght 
			//(use xc as a auxiliary vector) !!!!!!!!!!!!
			double aux;
			noalias(xc) = c4; noalias(xc)-=c1; aux = inner_prod(xc,xc); hmin = aux; hmax = aux;
			noalias(xc) = c3; noalias(xc)-=c1; aux = inner_prod(xc,xc); if(aux<hmin) hmin = aux; else if(aux>hmax) hmax = aux;
			noalias(xc) = c2; noalias(xc)-=c1; aux = inner_prod(xc,xc); if(aux<hmin) hmin = aux;
			noalias(xc) = c4; noalias(xc)-=c2; aux = inner_prod(xc,xc); if(aux<hmin) hmin = aux;
			noalias(xc) = c3; noalias(xc)-=c2; aux = inner_prod(xc,xc); if(aux<hmin) hmin = aux;
			noalias(xc) = c4; noalias(xc)-=c3; aux = inner_prod(xc,xc); if(aux<hmin) hmin = aux;
			hmin = sqrt(hmin);
			hmax = sqrt(hmax);

			//calculation of the jacobian
			mJ(0,0) = x1-x0; mJ(0,1) = y1-y0; mJ(0,2) = z1-z0;
			mJ(1,0) = x2-x0; mJ(1,1) = y2-y0; mJ(1,2) = z2-z0;
			mJ(2,0) = x3-x0; mJ(2,1) = y3-y0; mJ(2,2) = z3-z0;

			//inverse of the jacobian
			//first column
			mJinv(0,0) = mJ(1,1)*mJ(2,2) - mJ(1,2)*mJ(2,1);
			mJinv(1,0) = -mJ(1,0)*mJ(2,2) + mJ(1,2)*mJ(2,0);
			mJinv(2,0) = mJ(1,0)*mJ(2,1) - mJ(1,1)*mJ(2,0);		
			//second column
			mJinv(0,1) = -mJ(0,1)*mJ(2,2) + mJ(0,2)*mJ(2,1);
			mJinv(1,1) = mJ(0,0)*mJ(2,2) - mJ(0,2)*mJ(2,0);
			mJinv(2,1) = -mJ(0,0)*mJ(2,1) + mJ(0,1)*mJ(2,0);
			//third column
			mJinv(0,2) = mJ(0,1)*mJ(1,2) - mJ(0,2)*mJ(1,1);
			mJinv(1,2) = -mJ(0,0)*mJ(1,2) + mJ(0,2)*mJ(1,0);
			mJinv(2,2) = mJ(0,0)*mJ(1,1) - mJ(0,1)*mJ(1,0);
			//calculation of determinant (of the input matrix)

			double detJ = mJ(0,0)*mJinv(0,0) 
				+ mJ(0,1)*mJinv(1,0)
				+ mJ(0,2)*mJinv(2,0);

			volume = detJ * 0.16666666667;

			if(volume < 1e-3 * hmax*hmax*hmax)  //this is a sliver and we will remove it
			{
				//very bad element // ser a very high center for it to be eliminated
				xc[0]=0.0; xc[1]=0.0; xc[2]=0.0;
				radius = 10000000;
			}
			else
			{

				double x0_2 = x0*x0 + y0*y0 + z0*z0;
				double x1_2 = x1*x1 + y1*y1 + z1*z1; 
				double x2_2 = x2*x2 + y2*y2 + z2*z2; 
				double x3_2 = x3*x3 + y3*y3 + z3*z3; 

				//finalizing the calculation of the inverted matrix
				mJinv /= detJ;

				//calculating the RHS
				//calculating the RHS
				mRhs[0] = 0.5*(x1_2 - x0_2);
				mRhs[1] = 0.5*(x2_2 - x0_2);
				mRhs[2] = 0.5*(x3_2 - x0_2);

				//calculate position of the center
				noalias(xc) = prod(mJinv,mRhs);

				//calculate radius
				radius = pow(xc[0] - x0,2);
				radius		  += pow(xc[1] - y0,2);
				radius		  += pow(xc[2] - z0,2);
				radius = sqrt(radius);
			}
			KRATOS_CATCH("")
		}

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
		TetGenRefine& operator=(TetGenRefine const& rOther);


		///@}    

	}; // Class TetGenRefine 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		TetGenRefine& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const TetGenRefine& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_TETGEN_REFINER_H_INCLUDED  defined 


