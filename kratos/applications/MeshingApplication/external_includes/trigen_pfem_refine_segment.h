//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2009-01-22 17:13:57 $
//   Revision:            $Revision: 1.5 $
//
//


#if !defined(KRATOS_TRIGEN_PFEM_REFINE_SEGMENT_H_INCLUDED )
#define  KRATOS_TRIGEN_PFEM_REFINE_SEGMENT_H_INCLUDED

 

// System includes
#include <string>
#include <iostream> 
#include <stdlib.h>

#if !defined(KRATOS_TRIANGLE_EXTERNAL_H_INCLUDED)
#define  KRATOS_TRIANGLE_EXTERNAL_H_INCLUDED
#include "triangle.h" 
#endif

#include <boost/timer.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "geometries/triangle_2d_3.h"
#include "meshing_application.h"
#include "processes/node_erase_process.h" 
#include "processes/find_nodal_neighbours_process.h" 
#include "spatial_containers/spatial_containers.h"



namespace Kratos
{
	extern "C" {
		void triangulate(char *, struct triangulateio *, struct triangulateio *,struct triangulateio *);
		//void trifree();
	}


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
	class TriGenPFEMRefineSegment  

	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of TriGenModeler
		KRATOS_CLASS_POINTER_DEFINITION(TriGenPFEMRefineSegment);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		TriGenPFEMRefineSegment() :
		    mJ(ZeroMatrix(2,2)), //local jacobian
		    mJinv(ZeroMatrix(2,2)), //inverse jacobian
		    mC(ZeroVector(2)), //dimension = number of nodes
		    mRhs(ZeroVector(2)){}
		    //mpNodeEraseProcess(NULL){} //dimension = number of nodes

		/// Destructor.
		virtual ~TriGenPFEMRefineSegment(){}


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
			NodeEraseProcess& node_erase, bool rem_nodes = true, bool add_nodes=true,
			double my_alpha = 1.4, double h_factor=0.5)
		{

			KRATOS_TRY
			if (ThisModelPart.NodesBegin()->SolutionStepsDataHas(IS_FREE_SURFACE)==false )
				KRATOS_ERROR(std::logic_error,"Add  ----IS_FREE_SURFACE---- variable!!!!!! ERROR","");
			if (ThisModelPart.NodesBegin()->SolutionStepsDataHas(IS_STRUCTURE)==false )
				KRATOS_ERROR(std::logic_error,"Add  ----IS_STRUCTURE---- variable!!!!!! ERROR","");
			if (ThisModelPart.NodesBegin()->SolutionStepsDataHas(IS_BOUNDARY)==false )
				KRATOS_ERROR(std::logic_error,"Add  ----IS_BOUNDARY---- variable!!!!!! ERROR","");
			if (ThisModelPart.NodesBegin()->SolutionStepsDataHas(IS_FLUID)==false )
				KRATOS_ERROR(std::logic_error,"Add  ----IS_FLUID---- variable!!!!!! ERROR","");
			if (ThisModelPart.NodesBegin()->SolutionStepsDataHas(IS_WATER)==false )
				KRATOS_ERROR(std::logic_error,"Add  ----IS_WATER---- variable!!!!!! ERROR","");
			if (ThisModelPart.NodesBegin()->SolutionStepsDataHas(IS_INTERFACE)==false )
				KRATOS_ERROR(std::logic_error,"Add  ----IS_INTERFACE---- variable!!!!!! ERROR","");

			KRATOS_WATCH("Trigen PFEM Refining Segment Mesher")
			boost::timer auxiliary;
			

//clearing elements
			KRATOS_WATCH(ThisModelPart.Elements().size());
//KRATOS_WATCH("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
			//****************************************************
			//filling the interface list before removing the elements and nodes
			//****************************************************
			std::vector <int> mid_seg_list;
		        int seg_num = 0;
		        int num_interface = 0;
			SegmentDetecting(ThisModelPart,  mid_seg_list, seg_num, num_interface, h_factor);


			ThisModelPart.Elements().clear();
			ThisModelPart.Conditions().clear();

			////////////////////////////////////////////////////////////
			typedef Node<3> PointType;
			typedef Node<3>::Pointer PointPointerType;
			typedef std::vector<PointType::Pointer>           PointVector;
			typedef PointVector::iterator PointIterator;
			typedef std::vector<double>               DistanceVector;
			typedef std::vector<double>::iterator     DistanceIterator;

			int step_data_size = ThisModelPart.GetNodalSolutionStepDataSize();
			KRATOS_WATCH(step_data_size);
			// bucket types
			//typedef Bucket<3, PointType, ModelPart::NodesContainerType, PointPointerType, PointIterator, DistanceIterator > BucketType;
			//typedef Bins< 3, PointType, PointVector, PointPointerType, PointIterator, DistanceIterator > StaticBins;
			// bucket types
			typedef Bucket<3, PointType, PointVector, PointPointerType, PointIterator, DistanceIterator > BucketType;
		
					
			//*************
			// DynamicBins;	
			typedef Tree< KDTreePartition<BucketType> > kd_tree; //Kdtree;
			//typedef Tree< StaticBins > Bin; 			     //Binstree;
			unsigned int bucket_size = 20;

			//performing the interpolation - all of the nodes in this list will be preserved
			unsigned int max_results = 100;
			//PointerVector<PointType> res(max_results);
			//NodeIterator res(max_results);
			PointVector res(max_results);
			DistanceVector res_distances(max_results);
			Node<3> work_point(0,0.0,0.0,0.0);

			//****************************************************
			//filling the interface list before removing the nodes
			//****************************************************
		/*	for(ModelPart::NodesContainerType::const_iterator in = ThisModelPart.NodesBegin();
			     in != ThisModelPart.NodesEnd(); in++)
				{
					in->FastGetSolutionStepValue(IS_VISITED) = 0.0;
					in->FastGetSolutionStepValue(AUX_INDEX) = 0.0;
				}
			std::vector <int> mid_seg_list;
			//list_of_nodes.reserve(ThisModelPart.Nodes().size());
		        int seg_num = 0;


			for(ModelPart::NodesContainerType::const_iterator in = ThisModelPart.NodesBegin();
			     in != ThisModelPart.NodesEnd(); in++)
				{
			if(in->FastGetSolutionStepValue(IS_INTERFACE) == 1.0)
				  {

				 WeakPointerVector< Node<3> >& neighb = in->GetValue(NEIGHBOUR_NODES);
				 int num_of_intr_neigh = 0;
					//KRATOS_WATCH(neighb.size());		
			         for( WeakPointerVector< Node<3> >::iterator ngh_ind = neighb.begin(); ngh_ind!=neighb.end(); ngh_ind++)
				    {
					if(ngh_ind->FastGetSolutionStepValue(IS_INTERFACE) == 1.0) num_of_intr_neigh++;
				    }	

			  if(num_of_intr_neigh <= 2)
				{   
			         for( WeakPointerVector< Node<3> >::iterator ngh_ind = neighb.begin(); ngh_ind!=neighb.end(); ngh_ind++)
				    {
				
				if(ngh_ind->FastGetSolutionStepValue(IS_INTERFACE) == 1.0)   
				      {
					 if(ngh_ind->FastGetSolutionStepValue(IS_VISITED) != 10.0)
					{
//KRATOS_WATCH(ngh_ind->Id());
					  mid_seg_list.push_back(in->Id());
					  mid_seg_list.push_back(ngh_ind->Id());
					  seg_num++;
					  //row += 2;
					}
					  ngh_ind->FastGetSolutionStepValue(AUX_INDEX) += 1.0;
				      }
				    }
				  in->FastGetSolutionStepValue(IS_VISITED) = 10.0;
				 }
			 else
				{
				  in->FastGetSolutionStepValue(IS_VISITED) = 5.0;
				}
//KRATOS_WATCH("(((((((((after neighbour loop))))))))))))))))");

				   }
				}

			//conecting remaining interface nodes( those which have more than two interface neighbors)
			for(ModelPart::NodesContainerType::const_iterator in = ThisModelPart.NodesBegin();
			     in != ThisModelPart.NodesEnd(); in++)
				{
			         if(in->FastGetSolutionStepValue(IS_INTERFACE) == 1.0 && in->FastGetSolutionStepValue(IS_VISITED) == 5.0)
				   {
				     if(in->FastGetSolutionStepValue(AUX_INDEX) < 2.0)
				       {
				        WeakPointerVector< Node<3> >& neighb = in->GetValue(NEIGHBOUR_NODES);
			                  for( WeakPointerVector< Node<3> >::iterator ngh_ind = neighb.begin(); ngh_ind!=neighb.end(); ngh_ind++)
				             {
 					       if(ngh_ind->FastGetSolutionStepValue(IS_INTERFACE) == 1.0 && ngh_ind->FastGetSolutionStepValue(IS_VISITED) != 10.0 && ngh_ind->FastGetSolutionStepValue(AUX_INDEX) < 2.0)
					         {
						 	 mid_seg_list.push_back(in->Id());
					 		 mid_seg_list.push_back(ngh_ind->Id());
					 		 seg_num++;
							KRATOS_WATCH("**** A COMPLEX BOUNDARY IS DETECTED *******");
						 }
					  ngh_ind->FastGetSolutionStepValue(AUX_INDEX) += 1.0;
					      }

					}
				           in->FastGetSolutionStepValue(IS_VISITED) = 10.0;
				    }
				 }*/

 			//if the remove_node switch is activated, we check if the nodes got too close
			if (rem_nodes==true)
			{
	
				PointVector list_of_nodes;
				list_of_nodes.reserve(ThisModelPart.Nodes().size());
				for(ModelPart::NodesContainerType::iterator i_node = ThisModelPart.NodesBegin() ; i_node != ThisModelPart.NodesEnd() ; i_node++)
				{
						(list_of_nodes).push_back(*(i_node.base()));
				}

				kd_tree  nodes_tree1(list_of_nodes.begin(),list_of_nodes.end(), bucket_size);
					
				unsigned int n_points_in_radius;			
				//radius means the distance, closer than which no node shall be allowd. if closer -> mark for erasing
				double radius;
				
				for(ModelPart::NodesContainerType::const_iterator in = ThisModelPart.NodesBegin();
					in != ThisModelPart.NodesEnd(); in++)
					{
					radius=h_factor*in->FastGetSolutionStepValue(NODAL_H);
				
					work_point[0]=in->X();
					work_point[1]=in->Y();
					work_point[2]=in->Z();
							
					n_points_in_radius = nodes_tree1.SearchInRadius(work_point, radius, res.begin(),res_distances.begin(), max_results);
						if (n_points_in_radius>1)
						{
					    if (in->FastGetSolutionStepValue(IS_BOUNDARY)==0.0 && in->FastGetSolutionStepValue(IS_STRUCTURE)==0.0 && in->FastGetSolutionStepValue(IS_INTERFACE)==0.0)
							{
								//look if we are already erasing any of the other nodes 
								double erased_nodes = 0;
								for(PointIterator i=res.begin(); i!=res.begin() + n_points_in_radius ; i++)
									erased_nodes += (*i)->GetValue(ERASE_FLAG);

							/*								
								//to avoid remove of near boundary nodes
								double center_flag=in->FastGetSolutionStepValue(IS_WATER);
								for(PointIterator i=res.begin(); i!=res.begin() + n_points_in_radius ; i++)								
								{
								double ngh_flag = (*i)->FastGetSolutionStepValue(IS_WATER);
									if(center_flag != ngh_flag)
										erased_nodes++;	
								}
							*/
								if( erased_nodes < 1.0) //we cancel the node if no other nodes are being erased
									in->GetValue(ERASE_FLAG)=1;

							
							}
							else if ( (in)->FastGetSolutionStepValue(IS_STRUCTURE)!=1.0) //boundary nodes will be removed if they get REALLY close to another boundary node (0.2 * h_factor)
							{
								//here we loop over the neighbouring nodes and if there are nodes
								//with IS_BOUNDARY=1 which are closer than 0.2*nodal_h from our we remove the node we are considering
								unsigned int k = 0;
								unsigned int counter = 0;
								for(PointIterator i=res.begin(); i!=res.begin() + n_points_in_radius ; i++)
									{
										if ( (*i)->FastGetSolutionStepValue(IS_BOUNDARY,1)==1.0 && res_distances[k] < 0.2*radius && res_distances[k] > 0.0 )
										{
	// 										KRATOS_WATCH( res_distances[k] );
											counter += 1;
										}
										k++;
									}
								if(counter > 0)
									in->GetValue(ERASE_FLAG)=1;
							}
						}
				
					}
				//not erase INTERFACE node
				for(ModelPart::NodesContainerType::const_iterator in = ThisModelPart.NodesBegin();
					in != ThisModelPart.NodesEnd(); in++)
					{
					if((in)->FastGetSolutionStepValue(IS_INTERFACE) == 1.0)
						in->GetValue(ERASE_FLAG) = 0;

					}
						
				node_erase.Execute();				

				KRATOS_WATCH("Number of nodes after erasing")
				KRATOS_WATCH(ThisModelPart.Nodes().size())	
			}
			
			//////////////////////////////////////////////////////////////////////////////////////////////////////////
			//												  	//
			//		Now we shall pass the Alpha Shape for the second time, having the "bad nodes" removed	//
			//////////////////////////////////////////////////////////////////////////////////////////////////////////
			//creating the containers for the input and output
			struct triangulateio in_mid;
			struct triangulateio out_mid;
			struct triangulateio vorout_mid;

			initialize_triangulateio(in_mid);
			initialize_triangulateio(out_mid);
			initialize_triangulateio(vorout_mid);
			
			//assigning the size of the input containers
						
			in_mid.numberofpoints = ThisModelPart.Nodes().size();
			in_mid.pointlist = (REAL *) malloc(in_mid.numberofpoints * 2 * sizeof(REAL));

			//reorder the id's and give the coordinates to the mesher
			ModelPart::NodesContainerType::iterator nodes_begin = ThisModelPart.NodesBegin();
			std::vector <int> reorder_mid_seg_list;
			 if(seg_num != 0)
				reorder_mid_seg_list.resize(seg_num*2,false);

			for(unsigned int i = 0; i<ThisModelPart.Nodes().size(); i++)
			{ 
				int base = i*2;
				//int base = ((nodes_begin + i)->Id()   -  1 ) * 2;

				//from now on it is consecutive
				int pr_id = (nodes_begin + i)->Id();

				(nodes_begin + i)->SetId(i+1);
//				(nodes_begin + i)->Id() = i+1;

				in_mid.pointlist[base] = (nodes_begin + i)->X();
				in_mid.pointlist[base+1] = (nodes_begin + i)->Y();

                                Node<3>::DofsContainerType& node_dofs = (nodes_begin + i)->GetDofs();

                                for(Node<3>::DofsContainerType::iterator iii = node_dofs.begin();    iii != node_dofs.end(); iii++)
                                {
                                    iii->SetId(i+1);
//                                    iii->Id() = i+1;
                                }
				//reordering segment list
				if(seg_num != 0)
				 {
				   if( (nodes_begin + i)->FastGetSolutionStepValue(IS_INTERFACE) == 1.0 )
				      {
					for(int jj = 0; jj < seg_num*2; ++jj)
					   if( mid_seg_list[jj] == pr_id)
					  	reorder_mid_seg_list[jj] = (nodes_begin + i)->Id();

				      }	
				 }
			}
//(for(int ii = 0; ii<seg_num*2; ii++){
//KRATOS_WATCH(mid_seg_list[ii])
//KRATOS_WATCH(reorder_mid_seg_list[ii])}
			//in_mid.numberoftriangles = ThisModelPart.Elements().size();
			//in_mid.trianglelist = (int*) malloc(in_mid.numberoftriangles * 3 * sizeof(int));

			//***********preserving the list of interface segments


		        if(seg_num != 0){
				in_mid.numberofsegments = seg_num;
				in_mid.segmentlist = (int*) malloc(in_mid.numberofsegments * 2 * sizeof(int));
				in_mid.segmentmarkerlist = (int*) malloc(in_mid.numberofsegments * sizeof(int));
				//in2.segmentlist = seg_list;
					  for( int ii = 0; ii < seg_num*2; ++ii)
						in_mid.segmentlist[ii] = reorder_mid_seg_list[ii];

					  for( int jj = 0; jj < seg_num ; ++jj)
						in_mid.segmentmarkerlist[jj] = 5;

					}
			//***********end ofsegments
			KRATOS_WATCH(seg_num);
			//********************************REGIONAL ATTRIBUTE
			/*int num_water_node= 0;
			for(unsigned int i = 0; i<ThisModelPart.Nodes().size(); i++)
			{ 
				if((nodes_begin + i)->FastGetSolutionStepValue(IS_WATER) == 1.0)
					num_water_node++;

			}
			in_mid.numberofregions = num_water_node;
			in_mid.regionlist = (REAL *) malloc(in_mid.numberofregions * 3 * sizeof(REAL));

			int base = 0;
			for(unsigned int i = 0; i<ThisModelPart.Nodes().size(); i++)
			{ 

				if((nodes_begin + i)->FastGetSolutionStepValue(IS_WATER) == 1.0)
				{
				  in_mid.regionlist[base] = (nodes_begin + i)->X();
				  in_mid.regionlist[base + 1] = (nodes_begin + i)->Y();
				  in_mid.regionlist[base + 2] = 14.0;
				  base+=3;
				}			
			}*/

			/*for(unsigned int i = 0; i<ThisModelPart.Nodes().size() ; i++)
			 { 
			   if( (nodes_begin + i)->FastGetSolutionStepValue(IS_INTERFACE) != 1.0)
			     {
				if((nodes_begin + i)->FastGetSolutionStepValue(IS_WATER) == 0.0)
					in_mid.regionlist[i] = 7.0;
				else
					in_mid.regionlist[i] = 14.0;					
			     }
			   else
				in_mid.regionlist[i] = 0.0;		

			}*/
			//********************************end of REGIONAL ATTRIBUTE
			//for(unsigned int ii=0;ii<mid_seg_list.size(); ++ii)
			//	KRATOS_WATCH(in_mid.segmentlist[ii]);	

			// "P" suppresses the output .poly file. Saves disk space, but you 				
			//lose the ability to maintain constraining segments  on later refinements of the mesh. 
			// "B" Suppresses boundary markers in the output .node, .poly, and .edge output files
			// "n" outputs a list of triangles neighboring each triangle.
			// "e" outputs edge list (i.e. all the "connectivities")
			//char options1[] = "Pne";
			char options1[] = "pcne";
			KRATOS_WATCH(ThisModelPart.Nodes().size());
			
			triangulate(options1, &in_mid, &out_mid, &vorout_mid);
			//print out the mesh generation time
			std::cout<<"mesh generation time = "<<auxiliary.elapsed();
			//number of newly generated triangles
			unsigned int el_number=out_mid.numberoftriangles;
			KRATOS_WATCH("*********NUMBER OF ELEMENTS***********");
			KRATOS_WATCH(el_number);

			//prepairing for alpha shape passing : creating necessary arrays
			//list of preserved elements is created: at max el_number can be preserved (all elements)
			std::vector<int> preserved_list1(el_number);
			preserved_list1.resize(el_number);

			array_1d<double,3> x1,x2,x3,xc;
			///*******************************************************************************************
			///*******************************************************************************************
			//region flag check
		/*	for(unsigned int i = 0; i<num_water_node*3 ; i++)
			 { 
				KRATOS_WATCH(in_mid.regionlist[i]);
			 }
			for(unsigned int el = 0; el< el_number; el++)
			{
				KRATOS_WATCH("*""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""*");

				KRATOS_WATCH(out_mid.triangleattributelist[el]);
			}		*/
			//end of region flag check
			///*******************************************************************************************
			///*******************************************************************************************
			//int number_of_preserved_elems=0;
			int number_of_preserved_elems=0;
			int point_base;
			//loop for passing alpha shape			
			for(unsigned int el = 0; el< el_number; el++)
			{
				int base = el * 3;

				//coordinates
				point_base = (out_mid.trianglelist[base] - 1)*2;
				x1[0] = out_mid.pointlist[point_base]; 
				x1[1] = out_mid.pointlist[point_base+1]; 

				point_base = (out_mid.trianglelist[base+1] - 1)*2;
				x2[0] = out_mid.pointlist[point_base]; 
				x2[1] = out_mid.pointlist[point_base+1]; 

				point_base = (out_mid.trianglelist[base+2] - 1)*2;
				x3[0] = out_mid.pointlist[point_base]; 
				x3[1] = out_mid.pointlist[point_base+1]; 

				//here we shall temporarily save the elements and pass them afterwards to the alpha shape
				Geometry<Node<3> > temp;

				temp.push_back( *((ThisModelPart.Nodes()).find( out_mid.trianglelist[base]).base()	) );
				temp.push_back( *((ThisModelPart.Nodes()).find( out_mid.trianglelist[base+1]).base()	) );
				temp.push_back( *((ThisModelPart.Nodes()).find( out_mid.trianglelist[base+2]).base()	) );

				int number_of_structure_nodes = int( temp[0].FastGetSolutionStepValue(IS_STRUCTURE) );
				number_of_structure_nodes += int( temp[1].FastGetSolutionStepValue(IS_STRUCTURE) );
				number_of_structure_nodes += int( temp[2].FastGetSolutionStepValue(IS_STRUCTURE) );

				//check the number of nodes of boundary
				int nfs = int( temp[0].FastGetSolutionStepValue(IS_FREE_SURFACE) );
				nfs += int( temp[1].FastGetSolutionStepValue(IS_FREE_SURFACE) );
				nfs += int( temp[2].FastGetSolutionStepValue(IS_FREE_SURFACE) );
				
				//check the number of nodes of boundary
				int nfluid = int( temp[0].FastGetSolutionStepValue(IS_FLUID) );
				nfluid += int( temp[1].FastGetSolutionStepValue(IS_FLUID) );
				nfluid += int( temp[2].FastGetSolutionStepValue(IS_FLUID) );

				//check a three nodede interface element
				int num_interface = 0;
				for( int ii = 0; ii<3; ++ii)
					if(temp[ii].FastGetSolutionStepValue(IS_INTERFACE) == 1.0)
						num_interface++;

				/*if(num_interface == 3)
					{
						KRATOS_WATCH("OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO inside mesher w 3 interface");
				               // KRATOS_WATCH(out_mid.triangleattributelist[el]);						
					}
				*/

				//first check that we are working with fluid elements, otherwise throw an error
				//if (nfluid!=3)
				//	KRATOS_ERROR(std::logic_error,"THATS NOT FLUID or NOT TRIANGLE!!!!!! ERROR","");
				//otherwise perform alpha shape check

				
				if(number_of_structure_nodes!=3) //if it is = 3 it is a completely fixed element -> do not add it
				{
					if (nfs != 0 || nfluid != 3)  //in this case it is close to the surface so i should use alpha shape 
					{
						
						if( AlphaShape(my_alpha,temp) && number_of_structure_nodes!=3) //if alpha shape says to preserve
						{
							preserved_list1[el] = true;
							number_of_preserved_elems += 1;
														
						}
					}
					else //internal triangle --- should be ALWAYS preserved
					{							
						double bigger_alpha = my_alpha*5.0;
						if( AlphaShape(bigger_alpha,temp) && number_of_structure_nodes!=3) 
							{
							preserved_list1[el] = true;
							number_of_preserved_elems += 1;							
							}
					}				
				}
				else 
					preserved_list1[el] = false;
			}
			//freeing memory 


			//NOW WE SHALL PERFORM ADAPTIVE REMESHING, i.e. insert and remove nodes based upon mesh quality
			// and prescribed sizes
			struct triangulateio in2;
			struct triangulateio out2;
			struct triangulateio vorout2;
			
			initialize_triangulateio(in2);
			initialize_triangulateio(out2);
			initialize_triangulateio(vorout2);

//			in2.firstnumber = 1;
			in2.numberofpoints = ThisModelPart.Nodes().size();
			in2.pointlist = (REAL *) malloc(in2.numberofpoints * 2 * sizeof(REAL));

			//writing the input point list
			for(unsigned int i = 0; i<ThisModelPart.Nodes().size(); i++)
			{
				int base = i*2;
				in2.pointlist[base] = (nodes_begin + i)->X();
				in2.pointlist[base+1] = (nodes_begin + i)->Y();
			}
			in2.numberoftriangles=number_of_preserved_elems;

			in2.trianglelist = (int *) malloc(in2.numberoftriangles * 3 * sizeof(int));
			in2.trianglearealist = (REAL *) malloc(in2.numberoftriangles * sizeof(REAL));
//			in2.trianglelist = new int[in2.numberoftriangles * 3];
//			in2.trianglearealist = new double[in2.numberoftriangles];


			int counter = 0;
			for (unsigned int el=0; el<el_number; el++)
			{
				if( preserved_list1[el] == true )
				{
					//saving the list of ONLY preserved triangles, the ones that passed alpha-shape check
					int new_base = counter*3;
					int old_base = el*3;
					//copying in case it is preserved					
					in2.trianglelist[new_base] = out_mid.trianglelist[old_base];
					in2.trianglelist[new_base+1] = out_mid.trianglelist[old_base+1];
					in2.trianglelist[new_base+2] = out_mid.trianglelist[old_base+2];

					//calculate the prescribed h
					double prescribed_h = (nodes_begin + out_mid.trianglelist[old_base]-1)->FastGetSolutionStepValue(NODAL_H);
					prescribed_h += (nodes_begin + out_mid.trianglelist[old_base+1]-1)->FastGetSolutionStepValue(NODAL_H);
					prescribed_h += (nodes_begin + out_mid.trianglelist[old_base+2]-1)->FastGetSolutionStepValue(NODAL_H);
					prescribed_h *= 0.3333;
					//if h is the height of a equilateral triangle, the area is 1/2*h*h					
					in2.trianglearealist[counter] = 0.5*(1.5*prescribed_h*1.5*prescribed_h);
					counter += 1;
				}
			
			}
			//***********preserving the list of interface segments
			/*std::vector <int> seg_list;
			//list_of_nodes.reserve(ThisModelPart.Nodes().size());
		         seg_num = 0;
			//int row = 0;
	//	FindNodalNeighboursProcess(ThisModelPart,10,10).Execute();
			for(ModelPart::NodesContainerType::const_iterator in = ThisModelPart.NodesBegin();
			     in != ThisModelPart.NodesEnd(); in++)
				{
					in->FastGetSolutionStepValue(IS_VISITED) = 0.0;
				}

			for(ModelPart::NodesContainerType::const_iterator in = ThisModelPart.NodesBegin();
			     in != ThisModelPart.NodesEnd(); in++)
				{
					//KRATOS_WATCH("&&&&&&&&&&&&&&&&& inside node loops&&&&&&&&&&&&&&&&&&&&&&&&&");
				if(in->FastGetSolutionStepValue(IS_INTERFACE) == 1.0)
				  {
					//KRATOS_WATCH("@@@@@@@@@@@@@@@@2 an interface is detected@@@@@@@@@@@@@@@");

				 WeakPointerVector< Node<3> >& neighb = in->GetValue(NEIGHBOUR_NODES);
										     
			         for( WeakPointerVector< Node<3> >::iterator ngh_ind = neighb.begin(); ngh_ind!=neighb.end(); ngh_ind++)
				    {

				    if(ngh_ind->FastGetSolutionStepValue(IS_INTERFACE) == 1.0 && ngh_ind->FastGetSolutionStepValue(IS_VISITED) != 10.0)
					{
					  seg_list.push_back(in->Id());
					  seg_list.push_back(ngh_ind->Id());
					  seg_num++;
					  //row += 2;
					}
				    }
				  in->FastGetSolutionStepValue(IS_VISITED) = 10.0;
				   }
				}

*/
		        if(seg_num){
				in2.numberofsegments = seg_num;
				in2.segmentlist = (int*) malloc(in2.numberofsegments * 2 * sizeof(int));
				in2.segmentmarkerlist = (int*) malloc(in2.numberofsegments * sizeof(int));
				//in2.segmentlist = seg_list;
					  for( int ii = 0; ii < seg_num*2; ++ii)
						//in2.segmentlist[ii] = seg_list[ii];
						in2.segmentlist[ii] = reorder_mid_seg_list[ii];

					  for( int jj = 0; jj < seg_num ; ++jj)
						in2.segmentmarkerlist[jj] = 5;

					}
			//***********end ofsegments
			KRATOS_WATCH(seg_num);
			//for(unsigned int ii=0;ii<reorder_mid_seg_list.size(); ++ii)
			//	KRATOS_WATCH(in2.segmentlist[ii]);
                        //here we generate a new mesh adding/removing nodes, by initializing "q"-quality mesh and "a"-area constraint switches
			//
			// MOST IMPORTANT IS "r" switch, that refines previously generated mesh!!!!!!!!!!(that is the one given inside in2)
			//char mesh_regen_opts[] = "YYJaqrn";
			//char mesh_regen_opts[] = "YJq1.4arn";Yjaqrpcn
			//Y means dont insert node in boundary
			//YY meand no node on edges (inside)
		
			if (add_nodes==true)
				{
				char mesh_regen_opts[] = "Yjaqrpcn";
				triangulate(mesh_regen_opts, &in2, &out2, &vorout2);
				KRATOS_WATCH("Adaptive remeshing with segment executed")
				}
			else 
				{
				char mesh_regen_opts[] = "YJrn";
				triangulate(mesh_regen_opts, &in2, &out2, &vorout2);
				KRATOS_WATCH("Non-Adaptive remeshing executed")
				}

			//segment check
			/*unsigned int input_num_of_marker = in2.numberofsegments;
			unsigned int output_num_of_marker = out2.numberofsegments;
			unsigned int in_point_attr = in2.numberofpointattributes;
			unsigned int out_point_attr = out2.numberofpointattributes;*/

			/*for( int ii = 0; ii<num_of_marker; ii++)
				KRATOS_WATCH(out_mid.segmentmarkerlist[ii]);*/
			/*KRATOS_WATCH(input_num_of_marker);
			KRATOS_WATCH(output_num_of_marker);
			KRATOS_WATCH(in_point_attr);
			KRATOS_WATCH(out_point_attr);
			

			for(int ii = 0; ii<out2.numberofpoints; ii++)
				{
				KRATOS_WATCH("********* new pointr ***********")
				KRATOS_WATCH(out2.pointmarkerlist[ii])
				KRATOS_WATCH(out2.pointlist[ii*2])
				KRATOS_WATCH(out2.pointlist[ii*2 + 1])

				}

*/

			//end of segment check





			//and now we shall find out where the new nodes belong to
			//defintions for spatial search
			typedef Node<3> PointType;
			typedef Node<3>::Pointer PointPointerType;
			typedef std::vector<PointType::Pointer>           PointVector;
			//typedef std::vector<PointType::Pointer>::iterator PointIterator;
			typedef PointVector::iterator PointIterator;
			typedef std::vector<double>               DistanceVector;
			typedef std::vector<double>::iterator     DistanceIterator;


			typedef Bucket<3, PointType, PointVector, PointPointerType, PointIterator, DistanceIterator > BucketType;

			typedef Tree< KDTreePartition<BucketType> > kd_tree; //Kdtree;

			//int step_data_size = ThisModelPart.GetNodalSolutionStepDataSize();

			//creating an auxiliary list for the new nodes 
			PointVector list_of_new_nodes;
			
			//node to get the DOFs from					
			Node<3>::DofsContainerType& reference_dofs = (ThisModelPart.NodesBegin())->GetDofs();

			double z = 0.0;
			int n_points_before_refinement = in2.numberofpoints;
			//if points were added, we add them as nodes to the ModelPart
			if (out2.numberofpoints > n_points_before_refinement )
			{
			for(int i = n_points_before_refinement; i<out2.numberofpoints;i++)
				{
					int id=i+1;
					int base = i*2;
					double& x= out2.pointlist[base];
					double& y= out2.pointlist[base+1];

					Node<3>::Pointer pnode = ThisModelPart.CreateNewNode(id,x,y,z);

                                        pnode->SetBufferSize(ThisModelPart.NodesBegin()->GetBufferSize() );
							
					list_of_new_nodes.push_back( pnode );

					//assigning interface flag of segment
					//segment are assigned to 5 so every new node on them have also flag 5
					if(out2.pointmarkerlist[i] == 5)
					     {
						pnode->FastGetSolutionStepValue(IS_INTERFACE) = 1.0;
						KRATOS_WATCH("NEW INTERFACE ON SEGMEMNT is DETECTED");
						//KRATOS_WATCH(pnode->Id());
					     }
					else
						pnode->FastGetSolutionStepValue(IS_INTERFACE) = 0.0;
											
					//KRATOS_WATCH(pnode->FastGetSolutionStepValue(IS_INTERFACE));
					//std::cout << "new node id = " << pnode->Id() << std::endl;
					//generating the dofs
					for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin();    iii != reference_dofs.end(); iii++)
					{
						Node<3>::DofType& rDof = *iii;
						Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );
						
						(p_new_dof)->FreeDof();
//                                                (p_new_dof)->EquationId() = -1;

					}
					
				}
			}	
			
			std::cout << "During refinement we added " << out2.numberofpoints-n_points_before_refinement<< " nodes " <<std::endl;
			//unsigned int bucket_size = 20;
			//performing the interpolation - all of the nodes in this list will be preserved
			//unsigned int max_results = 50;
			//PointVector results(max_results);
			//DistanceVector results_distances(max_results);
			array_1d<double,3> N;
				
			//WHAT ARE THOSE????
// 			Node<3> work_point(0,0.0,0.0,0.0);
 			unsigned int MaximumNumberOfResults = 500;
			PointVector Results(MaximumNumberOfResults);
			DistanceVector ResultsDistances(MaximumNumberOfResults);

			step_data_size = ThisModelPart.GetNodalSolutionStepDataSize();


			if(out2.numberofpoints-n_points_before_refinement > 0) //if we added points
			{
						
				kd_tree  nodes_tree2(list_of_new_nodes.begin(),list_of_new_nodes.end(),bucket_size);
				nodes_begin = ThisModelPart.NodesBegin();
								
				for(int el = 0; el< in2.numberoftriangles; el++)
				{	
					int base = el * 3;
					//coordinates
					point_base = (in2.trianglelist[base] - 1)*2;
					x1[0] = in2.pointlist[point_base]; 
					x1[1] = in2.pointlist[point_base+1]; 

					point_base = (in2.trianglelist[base+1] - 1)*2;
					x2[0] = in2.pointlist[point_base]; 
					x2[1] = in2.pointlist[point_base+1]; 
					
					point_base = (in2.trianglelist[base+2] - 1)*2;
					x3[0] = in2.pointlist[point_base]; 
					x3[1] = in2.pointlist[point_base+1]; 
		
				
					//find the center and "radius" of the element
					double xc,  yc, radius;	
					CalculateCenterAndSearchRadius( x1[0], x1[1],
									x2[0], x2[1],
									x3[0], x3[1],
									xc,yc,radius);
									
					//find all of the new nodes within the radius
					int number_of_points_in_radius;
					work_point.X() = xc; work_point.Y() = yc; work_point.Z() = 0.0;
					
					number_of_points_in_radius = nodes_tree2.SearchInRadius(work_point, radius, Results.begin(),
							ResultsDistances.begin(),  MaximumNumberOfResults);

					Triangle2D3<Node<3> > geom(
						*( (nodes_begin +  in2.trianglelist[base]-1).base() 	), 
						*( (nodes_begin +  in2.trianglelist[base+1]-1).base() 	), 
						*( (nodes_begin +  in2.trianglelist[base+2]-1).base() 	)
						);	
	
					//check if inside and eventually interpolate
					for( PointIterator it_found = Results.begin(); it_found != Results.begin() + number_of_points_in_radius; it_found++)
					{
						bool is_inside = false;						
						is_inside = CalculatePosition(x1[0], x1[1],
										x2[0], x2[1],
										x3[0], x3[1],
										(*it_found)->X(),(*it_found)->Y(),N);
						
						
						if(is_inside == true)
						{
							double regionflag = 0.0;
							//regionflag = out_mid.triangleattributelist[el];
							Interpolate(  geom,  N, step_data_size, *(it_found ), regionflag);
							
						}
					}
				}		
			}

			ThisModelPart.Elements().clear();
			
			//set the coordinates to the original value
			for( PointVector::iterator it =  list_of_new_nodes.begin(); it!=list_of_new_nodes.end(); it++)
			{				
				const array_1d<double,3>& disp = (*it)->FastGetSolutionStepValue(DISPLACEMENT);
				(*it)->X0() = (*it)->X() - disp[0];
				(*it)->Y0() = (*it)->Y() - disp[1];
				(*it)->Z0() = 0.0;	
			}
			//cleaning unnecessary data
			//in2.deinitialize();
			//in2.initialize();
			//free( in2.pointlist );
			
			//add preserved elements to the kratos
			Properties::Pointer properties = ThisModelPart.GetMesh().pGetProperties(1);
			nodes_begin = ThisModelPart.NodesBegin();
			(ThisModelPart.Elements()).reserve(out2.numberoftriangles);
						
			for(int iii = 0; iii< out2.numberoftriangles; iii++)
			{
				int id = iii + 1;
				int base = iii * 3;
				Triangle2D3<Node<3> > geom(
					*( (nodes_begin +  out2.trianglelist[base]-1).base() 	), 
					*( (nodes_begin +  out2.trianglelist[base+1]-1).base() 	), 
					*( (nodes_begin +  out2.trianglelist[base+2]-1).base() 	)
					);


#ifdef _DEBUG
ModelPart::NodesContainerType& ModelNodes = ThisModelPart.Nodes();
				if( *(ModelNodes).find( out2.trianglelist[base]).base() == *(ThisModelPart.Nodes().end()).base() ) 
					KRATOS_ERROR(std::logic_error,"trying to use an inexisting node","");
				if( *(ModelNodes).find( out2.trianglelist[base+1]).base() == *(ThisModelPart.Nodes().end()).base() ) 
					KRATOS_ERROR(std::logic_error,"trying to use an inexisting node","");
				if( *(ModelNodes).find( out2.trianglelist[base+2]).base() == *(ThisModelPart.Nodes().end()).base() ) 
					KRATOS_ERROR(std::logic_error,"trying to use an inexisting node","");
				
#endif

				Element::Pointer p_element = rReferenceElement.Create(id, geom, properties);
				(ThisModelPart.Elements()).push_back(p_element);

			}
			ThisModelPart.Elements().Sort();

			//filling the neighbour list 
			ModelPart::ElementsContainerType::const_iterator el_begin = ThisModelPart.ElementsBegin();
			for(ModelPart::ElementsContainerType::const_iterator iii = ThisModelPart.ElementsBegin();
				iii != ThisModelPart.ElementsEnd(); iii++)
			{
				//Geometry< Node<3> >& geom = iii->GetGeometry();
				int base = ( iii->Id() - 1 )*3;

				(iii->GetValue(NEIGHBOUR_ELEMENTS)).resize(3);
				WeakPointerVector< Element >& neighb = iii->GetValue(NEIGHBOUR_ELEMENTS);
				for(int i = 0; i<3; i++)
				{
					int index = out2.neighborlist[base+i];
					if(index > 0)
						neighb(i) = *((el_begin + index-1).base());
					else
						neighb(i) = Element::WeakPointer();
				}
			}

			//reset the boundary flag
			for(ModelPart::NodesContainerType::const_iterator in = ThisModelPart.NodesBegin(); in!=ThisModelPart.NodesEnd(); in++)
			{
				in->FastGetSolutionStepValue(IS_BOUNDARY) = 0;
			}
			//filling the elemental neighbours list (from now on the elements list can not change)
			ModelPart::ElementsContainerType::iterator elements_end = ThisModelPart.Elements().end();

			ThisModelPart.Elements().Unique();
			
			//now the boundary faces
			for(ModelPart::ElementsContainerType::iterator iii = ThisModelPart.ElementsBegin();	iii != ThisModelPart.ElementsEnd(); iii++)
			{
				int base = ( iii->Id() - 1 )*3;
				
				ModelPart::ElementsContainerType::iterator el_neighb;
				/*each face is opposite to the corresponding node number so
				 0 ----- 2 1 
				 1 ----- 0 2
				 2 ----- 1 0 
				*/

				////finding boundaries and creating the "skin"
				//
				//********************************************************************
				//first face
				el_neighb = (ThisModelPart.Elements()).find( out2.neighborlist[base] );
				if( el_neighb == elements_end )
				{
					//std::cout << "node0" << std::endl;
					//if no neighnour is present => the face is free surface
					iii->GetGeometry()[1].FastGetSolutionStepValue(IS_BOUNDARY) = 1;
					iii->GetGeometry()[2].FastGetSolutionStepValue(IS_BOUNDARY) = 1;
				
					//Generate condition
					Condition::NodesArrayType temp;
					temp.reserve(2);
					temp.push_back(iii->GetGeometry()(2)); 
					temp.push_back(iii->GetGeometry()(1));
				
					Geometry< Node<3> >::Pointer cond = 
						Geometry< Node<3> >::Pointer(new Geometry< Node<3> >(temp) );
					int id = (iii->Id()-1)*3;
					
					//Condition::Pointer p_cond = 
					//	Condition::Pointer(new Condition(id, cond, properties) );	

					Condition::Pointer p_cond = rReferenceBoundaryCondition.Create(id, temp, properties);
					(ThisModelPart.Conditions()).push_back(p_cond);


				}
				
				//********************************************************************
				//second face
				el_neighb = (ThisModelPart.Elements()).find( out2.neighborlist[base+1] );
				//if( el != ThisModelPart.Elements().end() )
				if( el_neighb == elements_end )
				{
					//if no neighnour is present => the face is free surface
					iii->GetGeometry()[2].FastGetSolutionStepValue(IS_BOUNDARY) = 1;
					iii->GetGeometry()[0].FastGetSolutionStepValue(IS_BOUNDARY) = 1;

					//Generate condition
					Condition::NodesArrayType temp;
					temp.reserve(2);
					temp.push_back(iii->GetGeometry()(0)); 
					temp.push_back(iii->GetGeometry()(2));
					
					Geometry< Node<3> >::Pointer cond = 
						Geometry< Node<3> >::Pointer(new Geometry< Node<3> >(temp) );
					int id = (iii->Id()-1)*3+1;
					//
					//Condition::Pointer p_cond = 
					//	Condition::Pointer(new Condition(id, cond, properties) );

					Condition::Pointer p_cond = rReferenceBoundaryCondition.Create(id, temp, properties);
					(ThisModelPart.Conditions()).push_back(p_cond);
	

				}

				//********************************************************************
				//third face
				el_neighb = (ThisModelPart.Elements()).find( out2.neighborlist[base+2] );
				if( el_neighb == elements_end )
				{
					//if no neighnour is present => the face is free surface
					iii->GetGeometry()[0].FastGetSolutionStepValue(IS_BOUNDARY) = 1;
					iii->GetGeometry()[1].FastGetSolutionStepValue(IS_BOUNDARY) = 1;
					
//					Generate condition
					Condition::NodesArrayType temp;
					temp.reserve(2);
					temp.push_back(iii->GetGeometry()(1)); 
					temp.push_back(iii->GetGeometry()(0));
					Geometry< Node<3> >::Pointer cond = 
						Geometry< Node<3> >::Pointer(new Geometry< Node<3> >(temp) );
					int id = (iii->Id()-1)*3+2;
					
					//Condition::Pointer p_cond = 
					//	Condition::Pointer(new Condition(id, cond, properties) );

					Condition::Pointer p_cond = rReferenceBoundaryCondition.Create(id, temp, properties);
					(ThisModelPart.Conditions()).push_back(p_cond);

					

				}

			
			}



			clean_triangulateio(in2);

			clean_triangulateio(out2);

			clean_triangulateio(vorout2);
			KRATOS_WATCH("Finished remeshing with Trigen_PFEM_Refine_segment")

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
		 boost::numeric::ublas::bounded_matrix<double,2,2> mJ; //local jacobian
		 boost::numeric::ublas::bounded_matrix<double,2,2> mJinv; //inverse jacobian
		 array_1d<double,2> mC; //center pos
		 array_1d<double,2> mRhs; //center pos
		 //NodeEraseProcess* mpNodeEraseProcess;


		///@} 
		///@name Private Operators
		///@{ 
		//returns false if it should be removed
		bool AlphaShape(double alpha_param, Geometry<Node<3> >& pgeom)
			//bool AlphaShape(double alpha_param, Triangle2D<Node<3> >& pgeom)
		{
			KRATOS_TRY
			
			
			double x0 = pgeom[0].X();
			double x1 = pgeom[1].X();
			double x2 = pgeom[2].X();
			
			double y0 = pgeom[0].Y();
			double y1 = pgeom[1].Y();
			double y2 = pgeom[2].Y();
						
			mJ(0,0)=2.0*(x1-x0);	mJ(0,1)=2.0*(y1-y0);
			mJ(1,0)=2.0*(x2-x0);	mJ(1,1)=2.0*(y2-y0);
			
			
			double detJ = mJ(0,0)*mJ(1,1)-mJ(0,1)*mJ(1,0);
						
			mJinv(0,0) =  mJ(1,1); mJinv(0,1) = -mJ(0,1);
			mJinv(1,0) = -mJ(1,0); mJinv(1,1) =  mJ(0,0);
		
			bounded_matrix<double,2,2> check;

//calculate average h
			double h;
			h =  pgeom[0].FastGetSolutionStepValue(NODAL_H);
			h += pgeom[1].FastGetSolutionStepValue(NODAL_H);
			h += pgeom[2].FastGetSolutionStepValue(NODAL_H);
			h *= 0.333333333;
		
			
			if(detJ < 5e-3*h*h) 
			{
				//std::cout << "detJ = " << detJ << std::endl;
				////mark as boundary
				pgeom[0].GetSolutionStepValue(IS_BOUNDARY) = 1;
				pgeom[1].GetSolutionStepValue(IS_BOUNDARY) = 1;
				pgeom[2].GetSolutionStepValue(IS_BOUNDARY) = 1;
				return false;
			}
			
			else
			{

				double x0_2 = x0*x0 + y0*y0;
				double x1_2 = x1*x1 + y1*y1; 
				double x2_2 = x2*x2 + y2*y2; 

				//finalizing the calculation of the inverted matrix
				//std::cout<<"MATR INV"<<MatrixInverse(mJ)<<std::endl;
				mJinv /= detJ;
				//calculating the RHS
				mRhs[0] = (x1_2 - x0_2);
				mRhs[1] = (x2_2 - x0_2);

				//calculate position of the center
				noalias(mC) = prod(mJinv,mRhs);

				double radius = sqrt(pow(mC[0]-x0,2)+pow(mC[1]-y0,2));

				
				if (radius < h*alpha_param)
				{
					return true;
				}
				else
				{
					return false;
				}
			}
			

			KRATOS_CATCH("")
		}
		//AUXILLIARY FCTNS
		inline void CalculateCenterAndSearchRadius(const double x0, const double y0,
						const double x1, const double y1,
    						const double x2, const double y2,
	  					double& xc, double& yc, double& R
					       )
		{
			xc = 0.3333333333333333333*(x0+x1+x2);
			yc = 0.3333333333333333333*(y0+y1+y2);
			 
			double R1 = (xc-x0)*(xc-x0) + (yc-y0)*(yc-y0);
			double R2 = (xc-x1)*(xc-x1) + (yc-y1)*(yc-y1);
			double R3 = (xc-x2)*(xc-x2) + (yc-y2)*(yc-y2);
			
			R = R1;
			if(R2 > R) R = R2;
			if(R3 > R) R = R3;
			
			R = sqrt(R);
		}
		
		inline double CalculateVol(	const double x0, const double y0,
						const double x1, const double y1,
    						const double x2, const double y2
					  )
		{
			return 0.5*( (x1-x0)*(y2-y0)- (y1-y0)*(x2-x0) );
		}
		
		inline bool CalculatePosition(	const double x0, const double y0,
						const double x1, const double y1,
   						const double x2, const double y2,
						const double xc, const double yc,
						array_1d<double,3>& N		
					  )
		{
			double area = CalculateVol(x0,y0,x1,y1,x2,y2);
			double inv_area = 0.0;
			if(area < 0.000000000001)
			  {
				KRATOS_ERROR(std::logic_error,"element with zero area found","");
			  }
			else
			  {
				inv_area = 1.0 / area;
			  }
			
			  
			  N[0] = CalculateVol(x1,y1,x2,y2,xc,yc) * inv_area;
			  N[1] = CalculateVol(x2,y2,x0,y0,xc,yc) * inv_area;
			  N[2] = CalculateVol(x0,y0,x1,y1,xc,yc) * inv_area;
			  
/*			  N[0] = CalculateVol(x0,y0,x1,y1,xc,yc) * inv_area;
			N[1] = CalculateVol(x1,y1,x2,y2,xc,yc) * inv_area;
			N[2] = CalculateVol(x2,y2,x0,y0,xc,yc) * inv_area;*/
			
			if(N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0) //if the xc yc is inside the triangle
			//if(N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0) //if the xc yc is inside the triangle return true
				return true;
			
			return false;
		}	
		//**********************************************************************************************
		//**********************************************************************************************			
		void Interpolate( Triangle2D3<Node<3> >& geom, const array_1d<double,3>& N, 
				  unsigned int step_data_size,
      				Node<3>::Pointer pnode, double region_flag)
		{
			unsigned int buffer_size = pnode->GetBufferSize();
			//KRATOS_WATCH("Buffer size")
			//KRATOS_WATCH(buffer_size)
//KRATOS_WATCH(pnode->FastGetSolutionStepValue(IS_INTERFACE));
			//here we want to keep the flag of aaded segment node and not interpolate it
			double aux_interface = pnode->FastGetSolutionStepValue(IS_INTERFACE);

			for(unsigned int step = 0; step<buffer_size; step++)
			{	
				
						//getting the data of the solution step
				double* step_data = (pnode)->SolutionStepData().Data(step);
				
											
				double* node0_data = geom[0].SolutionStepData().Data(step);
				double* node1_data = geom[1].SolutionStepData().Data(step);
				double* node2_data = geom[2].SolutionStepData().Data(step);
								
				//copying this data in the position of the vector we are interested in
				for(unsigned int j= 0; j<step_data_size; j++)
				{ 

					step_data[j] = N[0]*node0_data[j] + N[1]*node1_data[j] + N[2]*node2_data[j];
								

				}						
			}
			if (N[0]==0.0 && N[1]==0.0 && N[2]==0.0) 
					KRATOS_ERROR(std::logic_error,"SOMETHING's wrong with the added nodes!!!!!! ERROR","");

			//if ( pnode->FastGetSolutionStepValue(BULK_MODULUS)==0.0) 
			//		KRATOS_ERROR(std::logic_error,"SOMETHING's wrong with the added nodes!!!!!! ERROR","");
			
			//now we assure that the flag variables are set coorect!! since we add nodes inside of the fluid volume only
			//we manually reset the IS_BOUNDARY, IS_FLUID, IS_STRUCTURE, IS_FREE_SURFACE values in a right way
			//not to have values, like 0.33 0.66 resulting if we would have been interpolating them in the same way 		
			//as the normal variables, like Velocity etc		
			
			
			pnode->FastGetSolutionStepValue(IS_BOUNDARY)=0.0;
			pnode->FastGetSolutionStepValue(IS_STRUCTURE)=0.0;
			pnode->GetValue(ERASE_FLAG)=0.0;
			pnode->FastGetSolutionStepValue(IS_FREE_SURFACE)=0.0;
			pnode->FastGetSolutionStepValue(IS_FLUID)=1.0;
			//pnode->FastGetSolutionStepValue(IS_VISITED)=0.0;

			//pnode->FastGetSolutionStepValue(IS_INTERFACE)=1.0;
			pnode->FastGetSolutionStepValue(IS_INTERFACE) = aux_interface;

			double same_colour = 0.0;
			for(int ii= 0; ii<= 2; ++ii)
				if(geom[ii].FastGetSolutionStepValue(IS_WATER) == 0.0)
							same_colour++;
			if(same_colour == 3.0)
				pnode->FastGetSolutionStepValue(IS_WATER) = 0.0;
			else
				pnode->FastGetSolutionStepValue(IS_WATER) = 1.0;

			/*if(region_flag == 14.0)
				{
				pnode->FastGetSolutionStepValue(IS_WATER) = 1.0;
				KRATOS_WATCH("RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR");
				KRATOS_WATCH(region_flag);
				}*/
			/*if(pnode->FastGetSolutionStepValue(IS_WATER)!= 1.0 && pnode->FastGetSolutionStepValue(IS_WATER)!= 0.0 && 
pnode->FastGetSolutionStepValue(IS_WATER)!=-1.0)	   
				{
				pnode->FastGetSolutionStepValue(IS_WATER) = 1.0;
				KRATOS_WATCH("$$$$$$$$$$$$$$$  A NEAR BOUNDARY NODE IS INSERTED $$$$$$$$$$$$$$$$$$$$$$");
				}*/
			/*if(pnode->FastGetSolutionStepValue(DISTANCE) <= 0.0)
				{pnode->FastGetSolutionStepValue(IS_WATER) = 0.0;}
			else
				{pnode->FastGetSolutionStepValue(IS_WATER) = 1.0;}*/

			if(aux_interface) 
					pnode->FastGetSolutionStepValue(IS_WATER) = 0.0;   
				

//KRATOS_WATCH(pnode->FastGetSolutionStepValue(IS_INTERFACE));


		}

		//**********************************************************************************************
		//**********************************************************************************************
		void SegmentDetecting(ModelPart& ThisModelPart, std::vector <int>& seg_list, int& seg_num, int& num_interface, double h_factor)
		{
		KRATOS_TRY;
			int Tdim = 2;
                        seg_list.clear();
			num_interface = 0;
			std::vector <int> nodes_of_bad_segments;
			std::vector <int> raw_seg_list;

		if(h_factor > 0.1)
			h_factor = 0.5;

		//delete interface flag
		   for(ModelPart::NodeIterator ind = ThisModelPart.NodesBegin(); ind != ThisModelPart.NodesEnd(); ++ind)
			ind->FastGetSolutionStepValue(IS_INTERFACE) = 0.0;



			for(ModelPart::ElementsContainerType::iterator elem = ThisModelPart.ElementsBegin(); 
					elem!=ThisModelPart.ElementsEnd(); elem++)
			  {
	
			    if(elem->GetValue(IS_WATER_ELEMENT) == 0)
				{

				 WeakPointerVector< Element >& neighbor_els = elem->GetValue(NEIGHBOUR_ELEMENTS);
				Geometry< Node<3> >& geom = elem->GetGeometry();

				 for(int ii=0; ii<(Tdim+1); ++ii)
					{

					 if(neighbor_els[ii].GetValue(IS_WATER_ELEMENT) == 1 && neighbor_els[ii].Id() != elem->Id())
					  {

					   if(ii == 0) // 2,1
						{
				if(geom[1].FastGetSolutionStepValue(IS_STRUCTURE) == 0.0 && geom[2].FastGetSolutionStepValue(IS_STRUCTURE) == 0.0)
						      {  

							//check for bad segments
							 double length = 0.0;
					 length = pow(geom[1].X()-geom[2].X(),2) + pow(geom[1].Y()-geom[2].Y(),2) +pow(geom[1].Z()-geom[2].Z(),2);
							 length = sqrt(length);

							 if(length < h_factor*geom[2].FastGetSolutionStepValue(NODAL_H))
							    {
								//detect bad segment and mark to earase first node
								nodes_of_bad_segments.push_back(geom[2].Id());
								nodes_of_bad_segments.push_back(geom[1].Id());

								
								//earse bad node
								geom[2].FastGetSolutionStepValue(IS_INTERFACE) = 0.0;
							        geom[2].GetValue(ERASE_FLAG) = 1.0;
							    }
							 else
							    {
								//one segment is detected
								++seg_num;
								raw_seg_list.push_back(geom[2].Id());
								raw_seg_list.push_back(geom[1].Id());
							        geom[2].FastGetSolutionStepValue(IS_INTERFACE) = 1.0;
							        geom[1].FastGetSolutionStepValue(IS_INTERFACE) = 1.0;
							    }
						      }
						}
					   if(ii == 1) // 0,2
						{
				if(geom[0].FastGetSolutionStepValue(IS_STRUCTURE) == 0.0 && geom[2].FastGetSolutionStepValue(IS_STRUCTURE) == 0.0)
						      {  
							//check for bad segments
							 double length = 0.0;
					 length = pow(geom[0].X()-geom[2].X(),2) + pow(geom[0].Y()-geom[2].Y(),2) +pow(geom[0].Z()-geom[2].Z(),2);
							 length = sqrt(length);

							 if(length < h_factor*geom[0].FastGetSolutionStepValue(NODAL_H))
							    {
								//detect bad segment and mark to earase first node
								nodes_of_bad_segments.push_back(geom[0].Id());
								nodes_of_bad_segments.push_back(geom[2].Id());

								//earse bad node
								geom[0].FastGetSolutionStepValue(IS_INTERFACE) = 0.0;
							        geom[0].GetValue(ERASE_FLAG) = 1.0;
							    }
							 else
							    {
								//one segment is detected
								++seg_num;
								raw_seg_list.push_back(geom[0].Id());
								raw_seg_list.push_back(geom[2].Id());
								geom[0].FastGetSolutionStepValue(IS_INTERFACE) = 1.0;
								geom[2].FastGetSolutionStepValue(IS_INTERFACE) = 1.0;
							    }
						      }
						}
					   if(ii == 2) // 1,0
						{
				if(geom[0].FastGetSolutionStepValue(IS_STRUCTURE) == 0.0 && geom[1].FastGetSolutionStepValue(IS_STRUCTURE) == 0.0)
						      {  
							//check for bad segments
							 double length = 0.0;
					 length = pow(geom[0].X()-geom[1].X(),2) + pow(geom[0].Y()-geom[1].Y(),2) +pow(geom[0].Z()-geom[1].Z(),2);
							 length = sqrt(length);

							 if(length < h_factor*geom[1].FastGetSolutionStepValue(NODAL_H))
							    {
								//detect bad segment and mark to earase first node
								nodes_of_bad_segments.push_back(geom[1].Id());
								nodes_of_bad_segments.push_back(geom[0].Id());

								//earse bad node
								geom[1].FastGetSolutionStepValue(IS_INTERFACE) = 0.0;
							        geom[1].GetValue(ERASE_FLAG) = 1.0;
							    }
							  else
							    {
								//one segment is detected
								++seg_num;
								raw_seg_list.push_back(geom[1].Id());
								raw_seg_list.push_back(geom[0].Id());
								geom[1].FastGetSolutionStepValue(IS_INTERFACE) = 1.0;
								geom[0].FastGetSolutionStepValue(IS_INTERFACE) = 1.0;
							    }
						      }
						}

					  }
					}

			  	}
			  }

		         //count interface nodes
			 for(ModelPart::NodeIterator ind = ThisModelPart.NodesBegin(); ind != ThisModelPart.NodesEnd(); ++ind)
				if(ind->FastGetSolutionStepValue(IS_INTERFACE) == 1.0)
					num_interface++;

/*			for(unsigned int seg=0; seg < raw_seg_list.size(); seg+=2)
			   {
					seg_list.push_back(raw_seg_list[seg]);
					seg_list.push_back(raw_seg_list[seg + 1]);

			   }


*/


			//merge for mark to erase node (the node marked for erase is relpaced byanother node of deleted segment on remaining segment)
			for(int ii=0; ii<nodes_of_bad_segments.size(); ii+=2)
			   {
				KRATOS_WATCH("!!!!!!!!!!!!!!!!!!! BAD NODES !!!!!!!!!!!!!!!!!!!!!!");

			    int bad_node = nodes_of_bad_segments[ii];	

//COORDINTAES of bad node
KRATOS_WATCH(ThisModelPart.Nodes()[bad_node].X());	
KRATOS_WATCH(ThisModelPart.Nodes()[bad_node].Y());
//ens of COORDINTAES
	  
			    for(int jj=0; jj<raw_seg_list.size(); jj++)
			      {
			//merge in segment list
			       if(raw_seg_list[jj] == bad_node )
					{
					raw_seg_list[jj] = nodes_of_bad_segments[ii + 1];
					KRATOS_WATCH("OOOOOOOOOOOOOOOOOO MERGE IS DONE OOOOOOOOOOOOOOOOOOOOOO");
					}
			      }
			//possible replace in nodes_of_bad_segments for adjacent bad segments
			    for(int kk=0; kk<nodes_of_bad_segments.size(); kk+=2)
			      {					
			       if(nodes_of_bad_segments[kk + 1] == bad_node )
					{
					nodes_of_bad_segments[kk + 1] = nodes_of_bad_segments[ii + 1];
					KRATOS_WATCH("HHHHHHHHHHHHHHHHHH BAD_SEGMENT LIST UPDATED HHHHHHHHHHHHHHHHHHHHHHHH");
					}
			      }


			    }


			//fill de final list
			seg_num = 0 ;
			for(int seg=0; seg < raw_seg_list.size(); seg+=2)
			   {
				if(raw_seg_list[seg] != raw_seg_list[seg + 1])
				    {
					seg_list.push_back(raw_seg_list[seg]);
					seg_list.push_back(raw_seg_list[seg + 1]);
					seg_num++;
				    }
			   }

		KRATOS_CATCH("");
		}
		//**********************************************************************************************
		//**********************************************************************************************
		void initialize_triangulateio( triangulateio& tr )
		{
			tr.pointlist                  = (REAL*) NULL;
			tr.pointattributelist         = (REAL*) NULL;
			tr.pointmarkerlist            = (int*) NULL;
			tr.numberofpoints             = 0;
			tr.numberofpointattributes    = 0;
			tr.trianglelist               = (int*) NULL;
			tr.triangleattributelist      = (REAL*) NULL;
			tr.trianglearealist           = (REAL*) NULL;
			tr.neighborlist               = (int*) NULL;
			tr.numberoftriangles          = 0;
			tr.numberofcorners            = 3;
			tr.numberoftriangleattributes = 0;
			tr.segmentlist                = (int*) NULL;
			tr.segmentmarkerlist          = (int*) NULL;
			tr.numberofsegments           = 0;
			tr.holelist                   = (REAL*) NULL;
			tr.numberofholes              = 0;
			tr.regionlist                 = (REAL*) NULL;
			tr.numberofregions            = 0;
			tr.edgelist                   = (int*) NULL;
			tr.edgemarkerlist             = (int*) NULL;
			tr.normlist                   = (REAL*) NULL;
			tr.numberofedges              = 0;
		};

		void clean_triangulateio( triangulateio& tr )
		{
                        if(tr.pointlist != NULL) free(tr.pointlist );
			if(tr.pointattributelist != NULL) free(tr.pointattributelist );
			if(tr.pointmarkerlist != NULL) free(tr.pointmarkerlist   );
			if(tr.trianglelist != NULL) free(tr.trianglelist  );
			if(tr.triangleattributelist != NULL) free(tr.triangleattributelist );
			if(tr.trianglearealist != NULL) free(tr.trianglearealist );
			if(tr.neighborlist != NULL) free(tr.neighborlist   );
			if(tr.segmentlist != NULL) free(tr.segmentlist    );
			if(tr.segmentmarkerlist != NULL) free(tr.segmentmarkerlist   );
			if(tr.holelist != NULL) free(tr.holelist      );
			if(tr.regionlist != NULL) free(tr.regionlist  );
			if(tr.edgelist != NULL) free(tr.edgelist   );
			if(tr.edgemarkerlist != NULL) free(tr.edgemarkerlist   );
			if(tr.normlist != NULL) free(tr.normlist  );
		};
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
		TriGenPFEMRefineSegment& operator=(TriGenPFEMRefineSegment const& rOther);


		///@}    

	}; // Class TriGenPFEMRefineSegment 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		TriGenPFEMRefineSegment& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const TriGenPFEMRefineSegment& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_TRIGEN_PFEM_MODELER_H_INCLUDED  defined 


				
