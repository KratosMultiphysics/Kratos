//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: kazem $
//   Date:                $Date: 2009-01-15 14:50:34 $
//   Revision:            $Revision: 1.8 $
//




#if !defined(KRATOS_TETGEN_PFEM_CONTACT_H_INCLUDED )
#define  KRATOS_TETGEN_PFEM_CONTACT_H_INCLUDED
  


// System includes
#include <string>
#include <iostream> 
#include <stdlib.h>
#include <boost/timer.hpp>



#include "tetgen.h" // Defined tetgenio, tetrahedralize().

// Project includes
#include "includes/define.h"
#include "utilities/geometry_utilities.h"
#include "includes/model_part.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "meshing_application.h"
#include "processes/node_erase_process.h"

#include "spatial_containers/spatial_containers.h"
//#include "containers/bucket.h"
//#include "containers/kd_tree.h"
//#include "external_includes/trigen_refine.h"
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
	class TetGenPfemContact  
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of TetGenPfemRefineFace
		KRATOS_CLASS_POINTER_DEFINITION(TetGenPfemContact);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
	    TetGenPfemContact() {}

		/// Destructor.
		virtual ~TetGenPfemContact(){}


		///@}
		///@name Operators 
		///@{


		///@}
		///@name Operations
		///@{


		//*******************************************************************************************
		//*******************************************************************************************
		void ReGenerateMesh(
			ModelPart& ThisModelPart ,
			Element const& rReferenceElement,
			Condition const& rReferenceBoundaryCondition 
			)
		{

			KRATOS_TRY
KRATOS_WATCH(">>>>>>> INSIDE MESHER <<<<<<<<<");
		     std::vector <int> shell_list;
		     shell_list.clear();
		     int shell_num = 0;
		     ModelPart::NodesContainerType shell_nodes;

		    for(ModelPart::ElementsContainerType::iterator elem = ThisModelPart.ElementsBegin(); 
					elem!=ThisModelPart.ElementsEnd(); elem++)
			  { 
				Geometry< Node<3> >& geom = elem->GetGeometry();
//KRATOS_WATCH(geom);

				shell_nodes.push_back(geom(0));
				shell_nodes.push_back(geom(1));
				shell_nodes.push_back(geom(2));

			  }

			  shell_nodes.Unique();

		    ModelPart::NodesContainerType::iterator nodes_begin = shell_nodes.begin();
		    for(unsigned int ii=0; ii < shell_nodes.size(); ii++)
			    ( nodes_begin + ii ) ->SetId( ii + 1 );

		    for(ModelPart::ElementsContainerType::iterator elem = ThisModelPart.ElementsBegin(); 
					elem!=ThisModelPart.ElementsEnd(); elem++)
			  { 
				++shell_num;
				Geometry< Node<3> >& geom = elem->GetGeometry();
				shell_list.push_back(geom[0].Id());
				shell_list.push_back(geom[1].Id());
				shell_list.push_back(geom[2].Id());
			  }

		    //mesh generation
		     tetgenio in_shell, out_shell;

 		     in_shell.firstnumber = 1;
		     in_shell.numberofpoints = shell_nodes.size();
		     in_shell.pointlist = new REAL[in_shell.numberofpoints * 3];


                     tetgenio::facet *f;
                     tetgenio::polygon *p;
			//give the corrdinates to the mesher
			for(unsigned int i = 0; i<shell_nodes.size(); i++)
			{
				int base = i*3;
				in_shell.pointlist[base] = (nodes_begin + i)->X();
				in_shell.pointlist[base+1] = (nodes_begin + i)->Y();
				in_shell.pointlist[base+2] = (nodes_begin + i)->Z();
			}


                       if(shell_num != 0){
                            int cnt = 0;
                            in_shell.numberoffacets = shell_num;
                            in_shell.facetmarkerlist = new int[in_shell.numberoffacets];
                            in_shell.facetlist = new tetgenio::facet[in_shell.numberoffacets];
                            for(int ii=0; ii<shell_num;++ii)
                            {
                                f = &in_shell.facetlist[ii];
                                f->numberofpolygons = 1;
                                f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
                                f->numberofholes = 0;
                                f->holelist = NULL;


                                p = &f->polygonlist[0];
                                p->numberofvertices = 3;
                                p->vertexlist = new int[p->numberofvertices];
                                p->vertexlist[0] = shell_list[cnt];
                                p->vertexlist[1] = shell_list[cnt + 1];
                                p->vertexlist[2] = shell_list[cnt + 2];
                                cnt +=3;
			
				in_shell.facetmarkerlist[ii] = 5;
                            }

                        //holes
                            in_shell.numberofholes = 0;
			    in_shell.holelist = NULL;

                        }
//KRATOS_WATCH("@@@@@@@@@@@@@ before meshing @@@@@@@@@@@@@@@@");
		//	in_shell.save_nodes("shell_mesh_in");
		//	in_shell.save_poly("shell_mesh_in");
			//char tetgen_options[] = "VMYYJ";pA
			char tetgen_options[] = "pYYJnQ";

                   //KRATOS_WATCH(in_shell.numberofpoints);

			tetrahedralize(tetgen_options, &in_shell, &out_shell); //with option to remove slivers

		  //     out_shell.save_nodes("shell_mesh_out");
		  //      out_shell.save_elements("shell_mesh_out");
		        //out_shell.save_faces("shell_mesh_out");
		//	out_shell.save_neighbors("shell_mesh_out");

			//deinitialize
			in_shell.deinitialize();
			in_shell.initialize();

			//clear elements
			ThisModelPart.Elements().clear();
			//ThisModelPart.Conditions().clear();
			

			//add contact elements
			boost::timer adding_elems;
			Properties::Pointer properties = ThisModelPart.GetMesh().pGetProperties(1);
			(ThisModelPart.Elements()).reserve(out_shell.numberoftetrahedra);
	KRATOS_WATCH(out_shell.numberoftetrahedra);
			for(int iii = 0; iii< out_shell.numberoftetrahedra; iii++)
			{
				int id = iii + 1;
				int base = iii * 4;
	//KRATOS_WATCH("inside 1");
				Tetrahedra3D4<Node<3> > geom(
					*( (nodes_begin +  out_shell.tetrahedronlist[base]-1).base() 		), 
					*( (nodes_begin +  out_shell.tetrahedronlist[base+1]-1).base() 	), 
					*( (nodes_begin +  out_shell.tetrahedronlist[base+2]-1).base() 	), 
					*( (nodes_begin +  out_shell.tetrahedronlist[base+3]-1).base() 	) 
					);
//KRATOS_WATCH(geom);
//KRATOS_WATCH("AFTER GEOM");
				Element::Pointer p_element = rReferenceElement.Create(id, geom, properties);
				p_element->GetValue(IS_CONTACT_MASTER) = 10;
				p_element->GetValue(IS_WATER_ELEMENT) = -10.0;
					//KRATOS_WATCH("inside 12");
							      
				(ThisModelPart.Elements()).push_back(p_element);
//KRATOS_WATCH(p_element->GetGeometry());
//KRATOS_WATCH("After get geometry");

			}
KRATOS_WATCH("After adding elements");
			ThisModelPart.Elements().Sort();
			
			boost::timer adding_neighb;
//			//filling the neighbour list 
			ModelPart::ElementsContainerType::const_iterator el_begin = ThisModelPart.ElementsBegin();
			for(ModelPart::ElementsContainerType::const_iterator iii = ThisModelPart.ElementsBegin();
				iii != ThisModelPart.ElementsEnd(); iii++)
			{
				//Geometry< Node<3> >& geom = iii->GetGeometry();
				int base = ( iii->Id() - 1 )*4;

				(iii->GetValue(NEIGHBOUR_ELEMENTS)).resize(4);
				WeakPointerVector< Element >& neighb = iii->GetValue(NEIGHBOUR_ELEMENTS);

				for(int i = 0; i<4; i++)
				{
					int index = out_shell.neighborlist[base+i];
					if(index > 0)
						neighb(i) = *((el_begin + index-1).base());
					else
						neighb(i) = *(iii.base()); 
				}
			}
			std::cout << "time for adding neigbours" << adding_neighb.elapsed() << std::endl;;

			//***********************************************************************************
			//***********************************************************************************
			//mark boundary nodes 
			//reset the boundary flag
/*
			for(ModelPart::NodesContainerType::const_iterator in = ThisModelPart.NodesBegin(); in!=ThisModelPart.NodesEnd(); in++)
			{
				in->FastGetSolutionStepValue(IS_BOUNDARY) = 0;
			}
*/

			//***********************************************************************************
			//***********************************************************************************
/*			boost::timer adding_faces;
KRATOS_WATCH(out_shell.numberoftrifaces);
			(ThisModelPart.Conditions()).reserve(out_shell.numberoftrifaces   );

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
KRATOS_WATCH(rReferenceBoundaryCondition);
KRATOS_WATCH(properties);
				//if( neighb(0).expired()  );
				if( out_shell.neighborlist[base] == -1)
				{
					CreateBoundaryFace(1, 2, 3, ThisModelPart,   0, *(iii.base()), properties,rReferenceBoundaryCondition );
				}
				//if(neighb(1).expired() );
				if( out_shell.neighborlist[base+1] == -1)
				{
					CreateBoundaryFace(0,3,2, ThisModelPart,   1, *(iii.base()), properties, rReferenceBoundaryCondition );
				}
				if( out_shell.neighborlist[base+2] == -1)
				//if(neighb(2).expired() );
				{
					CreateBoundaryFace(0,1,3, ThisModelPart,   2, *(iii.base()), properties, rReferenceBoundaryCondition );
				}
				if( out_shell.neighborlist[base+3] == -1)
				//if(neighb(3).expired() );
				{
					CreateBoundaryFace(0,2,1, ThisModelPart,   3, *(iii.base()), properties, rReferenceBoundaryCondition );
				}

			}
						
			std::cout << "time for adding faces" << adding_faces.elapsed() << std::endl;;
*/
			out_shell.deinitialize();
			out_shell.initialize();

			shell_nodes.clear();
			shell_list.clear();

KRATOS_WATCH(">>>>>>>>>>>>>>>>>>>>< BYE BYE MESHER <<<<<<<<<<<<<<<<<<<<<<<<");
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
/*
	      void CreateBoundaryFace(const int& i1, const int& i2, const int& i3, ModelPart& ThisModelPart, const int& outer_node_id, Element::Pointer origin_element, Properties::Pointer properties,Condition const& rReferenceBoundaryCondition)
		{
			KRATOS_TRY
KRATOS_WATCH("inside create boundary");
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
KRATOS_WATCH("face is created 1");
KRATOS_WATCH(geom);
KRATOS_WATCH(properties);
KRATOS_WATCH(rReferenceBoundaryCondition);
			Geometry< Node<3> >::Pointer cond = Geometry< Node<3> >::Pointer(new Triangle3D3< Node<3> >(temp) );
			//Geometry< Node<3> >::Pointer cond = Geometry< Node<3> >::Pointer(new Triangle3D< Node<3> >(temp) );
			int id = (origin_element->Id()-1)*4;
			//Condition::Pointer p_cond = Condition::Pointer(new Condition(id, cond, properties) );
			//Condition::Pointer p_cond = rReferenceBoundaryCondition::Pointer(new Condition(id, cond, properties) );
                        Condition::Pointer p_cond = rReferenceBoundaryCondition.Create(id, temp, properties);
			//assigning the neighbour node
KRATOS_WATCH("face is created 2");
			(p_cond->GetValue(NEIGHBOUR_NODES)).clear();
			(p_cond->GetValue(NEIGHBOUR_NODES)).push_back( Node<3>::WeakPointer( geom(outer_node_id) ) );
			(p_cond->GetValue(NEIGHBOUR_ELEMENTS)).clear();
			(p_cond->GetValue(NEIGHBOUR_ELEMENTS)).push_back( Element::WeakPointer( origin_element ) );
			ThisModelPart.Conditions().push_back(p_cond);
			KRATOS_CATCH("")
KRATOS_WATCH("face is created");
		}
*/
		
		//////////////////////////////////////////////////////////////////////////////////////

		///@} 
		///@name Member Variables 
		///@{ 


		///@} 
		///@name Private Operators
		///@{ 

		///@} 
		///@name Private Operations
		///@{ 
		

  
		
		//**********************************************************************************************
		//**********************************************************************************************

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
		TetGenPfemContact& operator=(TetGenPfemContact const& rOther);


		///@}    

	}; // Class TetGenPfemModeler 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		TetGenPfemContact& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const TetGenPfemContact& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_TETGEN_PFEM_CONTACT_H_INCLUDED 


