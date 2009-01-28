//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: antonia $
//   Date:                $Date: 2008-10-10 15:50:38 $
//   Revision:            $Revision: 1.5 $
//
//


#if !defined(KRATOS_TETGEN_MODELER_H_INCLUDED )
#define  KRATOS_TETGEN_MODELER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <stdlib.h>
#include <boost/timer.hpp>



#include "tetgen.h" // Defined tetgenio, tetrahedralize().

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "geometries/tetrahedra_3d_4.h"
#include "meshing_application.h"

#include "spatial_containers/spatial_containers.h"

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
	class TetGenModeler  
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of TetGenModeler
		KRATOS_CLASS_POINTER_DEFINITION(TetGenModeler);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
	    TetGenModeler() :
		mJ(ZeroMatrix(3,3)), //local jacobian
		mJinv(ZeroMatrix(3,3)), //inverse jacobian
		mc(ZeroVector(3)), //dimension = number of nodes
		mRhs(ZeroVector(3)){} //dimension = number of nodes

		/// Destructor.
		virtual ~TetGenModeler(){}


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
			Condition const& rReferenceBoundaryCondition,
			double alpha_param = 1.4)
		{

			KRATOS_TRY

			boost::timer auxiliary;

			//clearing elements
			ThisModelPart.Elements().clear();
			ThisModelPart.Conditions().clear();


			tetgenio in, out, in2, outnew;

			// All indices start from 1.
 			in.firstnumber = 1;

			in.numberofpoints = ThisModelPart.Nodes().size();
			in.pointlist = new REAL[in.numberofpoints * 3];

			//writing the point coordinates in a vector
			ModelPart::NodesContainerType::iterator nodes_begin = ThisModelPart.NodesBegin();

			//reorder node Ids
			for(unsigned int i = 0; i<ThisModelPart.Nodes().size(); i++)
			{
				(nodes_begin + i)->Id() = i+1;
			}

			//give the corrdinates to the mesher
			for(unsigned int i = 0; i<ThisModelPart.Nodes().size(); i++)
			{
				int base = i*3;
				in.pointlist[base] = (nodes_begin + i)->X();
				in.pointlist[base+1] = (nodes_begin + i)->Y();
				in.pointlist[base+2] = (nodes_begin + i)->Z();
			}
			
			char tetgen_options[] = "SJ";
			tetrahedralize(tetgen_options, &in, &out); //with option to remove slivers

			double first_part_time = auxiliary.elapsed();
			std::cout << "mesh generation time = " << first_part_time << std::endl;

			//generate Kratos Tetrahedra3D4
			int el_number = out.numberoftetrahedra;

			boost::timer alpha_shape_time;
			//calculate r , center for all of the elements
			std::vector<int> preserved_list(el_number);
			array_1d<double,3> x1,x2,x3,x4, xc;
			int point_base;
			int number_of_preserved_elems = 0;
			for(int el = 0; el< el_number; el++)
			{
				int base = el * 4;

				//coordinates
				point_base = (out.tetrahedronlist[base] - 1)*3;
				x1[0] = out.pointlist[point_base]; 
				x1[1] = out.pointlist[point_base+1]; 
				x1[2] = out.pointlist[point_base+2];

				point_base = (out.tetrahedronlist[base+1] - 1)*3;
				x2[0] = out.pointlist[point_base]; 
				x2[1] = out.pointlist[point_base+1]; 
				x2[2] = out.pointlist[point_base+2];

				point_base = (out.tetrahedronlist[base+2] - 1)*3;
				x3[0] = out.pointlist[point_base]; 
				x3[1] = out.pointlist[point_base+1]; 
				x3[2] = out.pointlist[point_base+2];

				point_base = (out.tetrahedronlist[base+3] - 1)*3;
				x4[0] = out.pointlist[point_base]; 
				x4[1] = out.pointlist[point_base+1]; 
				x4[2] = out.pointlist[point_base+2];

				//calculate the geometrical data 
				//degenerate elements are given a very high radius
				double geometrical_hmin, geometrical_hmax;
				double radius;
				double vol;
				CalculateElementData(x1,x2,x3,x4,vol,xc,radius,geometrical_hmin,geometrical_hmax);

				//calculate the prescribed h
				double prescribed_h = (nodes_begin + out.tetrahedronlist[base]-1)->FastGetSolutionStepValue(NODAL_H);
				prescribed_h += (nodes_begin + out.tetrahedronlist[base+1]-1)->FastGetSolutionStepValue(NODAL_H);
				prescribed_h += (nodes_begin + out.tetrahedronlist[base+2]-1)->FastGetSolutionStepValue(NODAL_H);
				prescribed_h += (nodes_begin + out.tetrahedronlist[base+3]-1)->FastGetSolutionStepValue(NODAL_H);
				prescribed_h *= 0.25;

				//check the number of nodes on the wall
				int nb = int( (nodes_begin + out.tetrahedronlist[base]-1)->FastGetSolutionStepValue(IS_STRUCTURE) );
				nb += int( (nodes_begin + out.tetrahedronlist[base+1]-1)->FastGetSolutionStepValue(IS_STRUCTURE) );
				nb += int( (nodes_begin + out.tetrahedronlist[base+2]-1)->FastGetSolutionStepValue(IS_STRUCTURE) );
				nb += int((nodes_begin + out.tetrahedronlist[base+3]-1)->FastGetSolutionStepValue(IS_STRUCTURE) );

				//check the number of nodes of boundary
				int nboundary = int( (nodes_begin + out.tetrahedronlist[base]-1)->FastGetSolutionStepValue(IS_FREE_SURFACE) );
				nboundary += int( (nodes_begin + out.tetrahedronlist[base+1]-1)->FastGetSolutionStepValue(IS_FREE_SURFACE) );
				nboundary += int( (nodes_begin + out.tetrahedronlist[base+2]-1)->FastGetSolutionStepValue(IS_FREE_SURFACE) );
				nboundary += int((nodes_begin + out.tetrahedronlist[base+3]-1)->FastGetSolutionStepValue(IS_FREE_SURFACE) );
				
				//check the number of nodes of boundary
				int nfluid = int( (nodes_begin + out.tetrahedronlist[base]-1)->FastGetSolutionStepValue(IS_FLUID) );
				nfluid += int( (nodes_begin + out.tetrahedronlist[base+1]-1)->FastGetSolutionStepValue(IS_FLUID) );
				nfluid += int( (nodes_begin + out.tetrahedronlist[base+2]-1)->FastGetSolutionStepValue(IS_FLUID) );
				nfluid += int((nodes_begin + out.tetrahedronlist[base+3]-1)->FastGetSolutionStepValue(IS_FLUID) );

				//cases:
				//4 nodes on the wall - elminate
				// at least one node of boundary OR at least one node NOT of fluid --> pass alpha shape
				
				if(nb == 4) // 4 nodes on the wall
					preserved_list[el] = false;
				else if (nboundary != 0 || nfluid != 4) //close to the free surface or external
				{
					if( radius  < prescribed_h * alpha_param && //alpha shape says to preserve
						nb!=4) //the nodes are not all on the boundary
					{
						preserved_list[el] = true; //preserve!!
						number_of_preserved_elems += 1;
					}					
				}
				else
				{
					preserved_list[el] = true;
					number_of_preserved_elems += 1;
				}

			}
			std::cout << "time for passing alpha shape" << alpha_shape_time.elapsed() << std::endl;

			//freeing unnecessary memory
			in.deinitialize();
			in.initialize();

			//setting the new new ids in a aux vector			
//			tetgenio in2;

  			in2.firstnumber = 1;
			in2.numberofpoints = ThisModelPart.Nodes().size();
			in2.pointlist = new REAL[in2.numberofpoints * 3];

			for(unsigned int i = 0; i<ThisModelPart.Nodes().size(); i++)
			{
				int base = i*3;
				in2.pointlist[base] = (nodes_begin + i)->X();
				in2.pointlist[base+1] = (nodes_begin + i)->Y();
				in2.pointlist[base+2] = (nodes_begin + i)->Z();
			}
std::cout << "qui" << std::endl;
			in2.numberoftetrahedra = number_of_preserved_elems;
			in2.tetrahedronlist = new int[in2.numberoftetrahedra * 4];

			int counter = 0;
			for(int el = 0; el< el_number; el++)
			{
				if( preserved_list[el] == true ) 
				{
					//saving the compact element list
					int new_base = counter*4;
					int old_base = el*4;
					in2.tetrahedronlist[new_base] = out.tetrahedronlist[old_base];
					in2.tetrahedronlist[new_base+1] = out.tetrahedronlist[old_base+1];
					in2.tetrahedronlist[new_base+2] = out.tetrahedronlist[old_base+2];
					in2.tetrahedronlist[new_base+3] = out.tetrahedronlist[old_base+3];

					counter += 1;
				}
			}

			//creating a new mesh
			boost::timer mesh_recreation_time;

			//freeing unnecessary memory
			out.deinitialize();
			out.initialize();
			

// 			std::cout << "Starting with nodal data copy..." << std::endl;
// 			unsigned int data_size = ThisModelPart.GetNodalSolutionStepTotalDataSize();
// 			in2.numberofpointattributes = data_size;
// 			in2.pointattributelist = new REAL[in2.numberofpoints*data_size];

// 			REAL* i_attribute = in2.pointattributelist;
			
// 			for(ModelPart::NodeIterator i_node = ThisModelPart.NodesBegin() ; i_node != ThisModelPart.NodesEnd() ; i_node++)
// 			{
// 			    memcpy(i_attribute, i_node->SolutionStepData().Data(), data_size*sizeof(double));
// 			    i_attribute += data_size;
// 			}

// 			std::cout << "Nodal data has been copied!!!" << std::endl;
			
			

//			tetgenio outnew;

// 			tetrahedralize("rn", &in2, &outnew); 
			char regeneration_options[] = "rqYVnJ";
			tetrahedralize(regeneration_options, &in2, &outnew);

 			std::cout << "mesh recreation time" << mesh_recreation_time.elapsed() << std::endl;;

			for(int i = in2.numberofpoints ; i < outnew.numberofpoints ; i++)
			{
				int base = i*3;
				ThisModelPart.CreateNewNode(i+1, outnew.pointlist[base], outnew.pointlist[base + 1], outnew.pointlist[base + 2]);
 			}
			
// 			i_attribute = outnew.pointattributelist + in2.numberofpoints * data_size;
//  			for(unsigned int i = in2.numberofpoints ; i < outnew.numberofpoints ; i++)
//  			{
// // // 			    std::cout << "node #" << i + 1 << " : ";
// // // 			    for(int j = 0 ; j < data_size ; j++)
// // // 			    {
// // // 				*(i_attribute + j) = 1000;
// // // 				std::cout << *(i_attribute + j) << ", ";
// // // 			    }
// // // 			    std::cout << std::endl;
// // 			    int base = i*3;
//  			    ThisModelPart.CreateNewNode(i+1, outnew.pointlist[base], outnew.pointlist[base + 1], outnew.pointlist[base + 2], i_attribute , 0);
// // 			    i_attribute += data_size;
// // // 			    KRATOS_WATCH(ThisModelPart.GetNode(i+1).GetSolutionStepValue(TEMPERATURE));
// // 			}
// // 			std::cout << "new nodes created!!!" << std::endl;

			//putting the new nodes in a spatial container
			typedef Node<3> PointType;
			typedef Node<3>::Pointer PointPointerType;
			typedef std::vector<PointType::Pointer>           PointVector;
			typedef std::vector<PointType::Pointer>::iterator PointIterator;
			typedef std::vector<double>               DistanceVector;
			typedef std::vector<double>::iterator     DistanceIterator;


			typedef Bucket<3, PointType, PointVector, PointPointerType, PointIterator, DistanceIterator > BucketType;

			typedef Tree< KDTreePartition<BucketType> > kd_tree; //Kdtree;

		
			//creating an auxiliary vector with pointers to the new nodes
			PointVector new_nodes;
			for(int i = in2.numberofpoints ; i < outnew.numberofpoints ; i++)
			{
				new_nodes.push_back( *(ThisModelPart.NodesBegin()+i).base() );
			}
						
			//filling the spatial container
			unsigned int bucket_size = 2;
			kd_tree  nodes_tree2(new_nodes.begin(),new_nodes.end(),bucket_size);
			
			//performing the interpolation - all of the nodes in this list will be preserved
			unsigned int max_results = 20;
			PointVector results(max_results);
			
			int data_size = ThisModelPart.GetNodalSolutionStepTotalDataSize();
			double* work_array;
			
			for(int el = 0; el< in2.numberoftetrahedra; el++)
			{
				int base = el * 4;

				//coordinates
				point_base = (in2.tetrahedronlist[base] - 1)*3;
				x1[0] = in2.pointlist[point_base]; 
				x1[1] = in2.pointlist[point_base+1]; 
				x1[2] = in2.pointlist[point_base+2];

				point_base = (in2.tetrahedronlist[base+1] - 1)*3;
				x2[0] = in2.pointlist[point_base]; 
				x2[1] = in2.pointlist[point_base+1]; 
				x2[2] = in2.pointlist[point_base+2];

				point_base = (in2.tetrahedronlist[base+2] - 1)*3;
				x3[0] = in2.pointlist[point_base]; 
				x3[1] = in2.pointlist[point_base+1]; 
				x3[2] = in2.pointlist[point_base+2];

				point_base = (in2.tetrahedronlist[base+3] - 1)*3;
				x4[0] = in2.pointlist[point_base]; 
				x4[1] = in2.pointlist[point_base+1]; 
				x4[2] = in2.pointlist[point_base+2];

				//calculate the geometrical data 
				//degenerate elements are given a very high radius
				double geometrical_hmin, geometrical_hmax;
				double radius;
				double vol;
				CalculateElementData(x1,x2,x3,x4,vol,xc,radius,geometrical_hmin,geometrical_hmax);
				
				//find all of the nodes in a radius centered in xc
				Node<3> xxx(12345,xc);
				nodes_tree2.SearchInRadius(xxx,radius,results.begin(),max_results);
				
				//for each of the nodes in radius find if it is inside the element or not
				//if inside interpolate
// 				for(PointIterator it=results.begin(); it!=results.end(); it++)
// 				{
// 					CheckIfInside(test_flag,N,xcheck);
// 					
// 					if(test_flag == true) //node is inside
// 					{
// 						
// 					}
// 				}
				
				



			}
			
			//cleaning unnecessary data
			in2.deinitialize();
			in2.initialize();

			//***********************************************************************************
			//***********************************************************************************
			boost::timer adding_elems;
			//add preserved elements to the kratos
			Properties::Pointer properties = ThisModelPart.GetMesh().pGetProperties(1);
			nodes_begin = ThisModelPart.NodesBegin();
			(ThisModelPart.Elements()).reserve(outnew.numberoftetrahedra);
			
			KRATOS_WATCH(properties);
			
std::cout << "ln334" << std::endl;
			
			for(int iii = 0; iii< outnew.numberoftetrahedra; iii++)
			{
				int id = iii + 1;
				int base = iii * 4;
				Tetrahedra3D4<Node<3> > geom(
					*( (nodes_begin +  outnew.tetrahedronlist[base]-1).base() 		), 
					*( (nodes_begin +  outnew.tetrahedronlist[base+1]-1).base() 	), 
					*( (nodes_begin +  outnew.tetrahedronlist[base+2]-1).base() 	), 
					*( (nodes_begin +  outnew.tetrahedronlist[base+3]-1).base() 	) 
					);


#ifdef _DEBUG
ModelPart::NodesContainerType& ModelNodes = ThisModelPart.Nodes();
				if( *(ModelNodes).find( outnew.tetrahedronlist[base]).base() == *(ThisModelPart.Nodes().end()).base() ) 
					KRATOS_ERROR(std::logic_error,"trying to use an inexisting node","");
				if( *(ModelNodes).find( outnew.tetrahedronlist[base+1]).base() == *(ThisModelPart.Nodes().end()).base() ) 
					KRATOS_ERROR(std::logic_error,"trying to use an inexisting node","");
				if( *(ModelNodes).find( outnew.tetrahedronlist[base+2]).base() == *(ThisModelPart.Nodes().end()).base() ) 
					KRATOS_ERROR(std::logic_error,"trying to use an inexisting node","");
				if( *(ModelNodes).find( outnew.tetrahedronlist[base+3]).base() == *(ThisModelPart.Nodes().end()).base() ) 
					KRATOS_ERROR(std::logic_error,"trying to use an inexisting node","");
#endif

				Element::Pointer p_element = rReferenceElement.Create(id, geom, properties);
				(ThisModelPart.Elements()).push_back(p_element);

			}
			std::cout << "time for adding elems" << adding_elems.elapsed() << std::endl;;
			ThisModelPart.Elements().Sort();
	
std::cout << "ln368" << std::endl;

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
					int index = outnew.neighborlist[base+i];
					if(index > 0)
						neighb(i) = *((el_begin + index-1).base());
					else
						neighb(i) = Element::WeakPointer();
				}
			}
			std::cout << "time for adding neigbours" << adding_neighb.elapsed() << std::endl;;


			//***********************************************************************************
			//***********************************************************************************
			//mark boundary nodes 
			//reset the boundary flag
			for(ModelPart::NodesContainerType::const_iterator in = ThisModelPart.NodesBegin(); in!=ThisModelPart.NodesEnd(); in++)
			{
				in->FastGetSolutionStepValue(IS_BOUNDARY) = 0;
			}

// 			//reset is fluid flag for the structure nodes
// 			for(ModelPart::NodesContainerType::const_iterator in = ThisModelPart.NodesBegin(); in!=ThisModelPart.NodesEnd(); in++)
// 			{
// 				if(in->FastGetSolutionStepValue(IS_STRUCTURE) )
// 					in->FastGetSolutionStepValue(IS_FLUID) = 0.0;
// 			}

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

				//if( neighb(0).expired()  );
				if( outnew.neighborlist[base] == -1)
				{
					CreateBoundaryFace(1, 2, 3, ThisModelPart,   0, *(iii.base()), properties );
				}
				//if(neighb(1).expired() );
				if( outnew.neighborlist[base+1] == -1)
				{
					CreateBoundaryFace(0,3,2, ThisModelPart,   1, *(iii.base()), properties );
				}
				if( outnew.neighborlist[base+2] == -1)
				//if(neighb(2).expired() );
				{
					CreateBoundaryFace(0,1,3, ThisModelPart,   2, *(iii.base()), properties );
				}
				if( outnew.neighborlist[base+3] == -1)
				//if(neighb(3).expired() );
				{
					CreateBoundaryFace(0,2,1, ThisModelPart,   3, *(iii.base()), properties );
				}

			}
			std::cout << "time for adding faces" << adding_faces.elapsed() << std::endl;;

			//setting free the memory
//			in.deinitialize();
//			outnew.deinitialize();

			//mark the nodes of structure which become "wet"
// 			for(ModelPart::NodesContainerType::const_iterator in = ThisModelPart.NodesBegin(); in!=ThisModelPart.NodesEnd(); in++)
// 			{
// 				//marking wet nodes
// 				if(in->FastGetSolutionStepValue(IS_STRUCTURE) )
// 					if( (in->GetValue(NEIGHBOUR_ELEMENTS)).size() != 0)
// 						in->FastGetSolutionStepValue(IS_FLUID) = 1.0;
// 				//marking as free surface the lonely nodes
// 				else
// 					if( (in->GetValue(NEIGHBOUR_ELEMENTS)).size() == 0)
// 						in->FastGetSolutionStepValue(IS_BOUNDARY) = 1.0;
// 					
// 			}


			double second_part_time = auxiliary.elapsed();
			std::cout << "second part time = " << second_part_time - first_part_time << std::endl;


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
		array_1d<double,3> mc; //center pos
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
			const double x0 = c1[0]; const double y0 = c1[1]; const double z0 = c1[2];
			const double x1 = c2[0]; const double y1 = c2[1]; const double z1 = c2[2];
			const double x2 = c3[0]; const double y2 = c3[1]; const double z2 = c3[2];
			const double x3 = c4[0]; const double y3 = c4[1]; const double z3 = c4[2];

// 			KRATOS_WATCH("111111111111111111111");
			//calculate min side lenght 
			//(use xc as a auxiliary vector) !!!!!!!!!!!!
			double aux;
			noalias(xc) = c4; noalias(xc)-=c1; aux = inner_prod(xc,xc); hmin = aux; hmax = aux;
			noalias(xc) = c3; noalias(xc)-=c1; aux = inner_prod(xc,xc); if(aux<hmin) hmin = aux; else if(aux>hmax) hmax = aux;
			noalias(xc) = c2; noalias(xc)-=c1; aux = inner_prod(xc,xc); if(aux<hmin) hmin = aux; else if(aux>hmax) hmax = aux;
			noalias(xc) = c4; noalias(xc)-=c2; aux = inner_prod(xc,xc); if(aux<hmin) hmin = aux; else if(aux>hmax) hmax = aux;
			noalias(xc) = c3; noalias(xc)-=c2; aux = inner_prod(xc,xc); if(aux<hmin) hmin = aux; else if(aux>hmax) hmax = aux;
			noalias(xc) = c4; noalias(xc)-=c3; aux = inner_prod(xc,xc); if(aux<hmin) hmin = aux; else if(aux>hmax) hmax = aux;
			hmin = sqrt(hmin);
			hmax = sqrt(hmax);

			//calculation of the jacobian
			mJ(0,0) = x1-x0; mJ(0,1) = y1-y0; mJ(0,2) = z1-z0;
			mJ(1,0) = x2-x0; mJ(1,1) = y2-y0; mJ(1,2) = z2-z0;
			mJ(2,0) = x3-x0; mJ(2,1) = y3-y0; mJ(2,2) = z3-z0;
// 			KRATOS_WATCH("33333333333333333333333");

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

// 			KRATOS_WATCH("44444444444444444444444444");
			double detJ = mJ(0,0)*mJinv(0,0) 
				+ mJ(0,1)*mJinv(1,0)
				+ mJ(0,2)*mJinv(2,0);

			volume = detJ * 0.16666666667;
// 			KRATOS_WATCH("55555555555555555555555");

			if(volume < 1e-3 * hmax*hmax*hmax)  //this is a sliver and we will remove it
			{
// 			KRATOS_WATCH("666666666666666666666666666");
				//very bad element // ser a very high center for it to be eliminated
				xc[0]=0.0; xc[1]=0.0; xc[2]=0.0;
				radius = 10000000;
// 			KRATOS_WATCH("777777777777777777777777");
			}
			else
			{

// 			KRATOS_WATCH("888888888888888888888888");
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
// 			KRATOS_WATCH("999999999999999999999999999");
			}
			KRATOS_CATCH("")
		}

		template< class T, std::size_t dim >
				class PointDistance{
					public:
						double operator()( T const& p1, T const& p2 ){
							double dist = 0.0;
							for( std::size_t i = 0 ; i < dim ; i++){
								double tmp = p1[i] - p2[i];
								dist += tmp*tmp;
							}
							return sqrt(dist);
						}
				};		

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
		TetGenModeler& operator=(TetGenModeler const& rOther);


		///@}    

	}; // Class TetGenModeler 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		TetGenModeler& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const TetGenModeler& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_TETGEN_MODELER_H_INCLUDED  defined 


