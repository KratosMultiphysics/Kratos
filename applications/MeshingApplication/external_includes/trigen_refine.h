//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2009-01-15 14:50:34 $
//   Revision:            $Revision: 1.11 $
//
//


#if !defined(KRATOS_TRIGEN_REFINE_H_INCLUDED )
#define  KRATOS_TRIGEN_REFINE_H_INCLUDED

 

// System includes
#include <string>
#include <iostream> 
#include <stdlib.h>

#if !defined(KRATOS_TRIANGLE_EXTERNAL_H_INCLUDED)
#define  KRATOS_TRIANGLE_EXTERNAL_H_INCLUDED
#include "triangle.h" 
#endif

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "geometries/triangle_2d_3.h"
#include "meshing_application.h"

#include "spatial_containers/spatial_containers.h"



namespace Kratos
{
	extern "C" {
		void triangulate(char *, struct triangulateio *, struct triangulateio *,struct triangulateio *);
		//void trifree();
	}

	template< class T, std::size_t dim >
			class DistanceCalculator{
				public:
					double operator()( T const& p1, T const& p2 ){
						double dist = 0.0;
						for( std::size_t i = 0 ; i < dim ; i++){
							double tmp = p1[i] - p2[i];
							dist += tmp*tmp;
						}
						return dist;
					}
			};
			
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
	class TriGenCDTrefine  
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of TriGenCDTrefine
		KRATOS_CLASS_POINTER_DEFINITION(TriGenCDTrefine);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		TriGenCDTrefine() {} //

		/// Destructor.
		virtual ~TriGenCDTrefine(){}


		///@}
		///@name Operators 
		///@{


		///@}
		///@name Operations
		///@{


		//*******************************************************************************************
		//*******************************************************************************************
		void RefineCDT(
			ModelPart& ThisModelPart , 
   			bool prescribe_element_size,
			Element const& rReferenceElement
			)
		{

			KRATOS_TRY


			struct triangulateio in;
			struct triangulateio mid;
			struct triangulateio out;
			struct triangulateio vorout;
			
			initialize_triangulateio(in);
			initialize_triangulateio(mid);
			initialize_triangulateio(out);
			initialize_triangulateio(vorout);
			
			//*********************************************************************
			//input mesh			
			in.numberofpoints = ThisModelPart.Nodes().size();
			in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
			
			//writing the points coordinates in a vector
			ModelPart::NodesContainerType::iterator nodes_begin = ThisModelPart.NodesBegin();
			for(unsigned int i = 0; i<ThisModelPart.Nodes().size(); i++)
			{
				int base = i*2;

				//from now on it is consecutive
				(nodes_begin + i)->Id() = i+1;

				in.pointlist[base] = (nodes_begin + i)->X();
				in.pointlist[base+1] = (nodes_begin + i)->Y();
			}			
			
			in.numberoftriangles = ThisModelPart.Elements().size();
			in.trianglelist = (int*) malloc(in.numberoftriangles * 3 * sizeof(int));
			
			//copying the elements in the input file
			ModelPart::ElementsContainerType::iterator elem_begin = ThisModelPart.ElementsBegin();
			for(unsigned int el = 0; el<ThisModelPart.Elements().size(); el++)
			{
				int base = el * 3;
				
				Geometry<Node<3> >& geom = (elem_begin+el)->GetGeometry();
				
				//saving the connectivity
				in.trianglelist[base] = geom[0].Id();
				in.trianglelist[base+1] = geom[1].Id();
				in.trianglelist[base+2] = geom[2].Id();
				
			}
			
			if (prescribe_element_size == true)
			{
				in.trianglearealist = (REAL *) malloc(in.numberoftriangles * sizeof(REAL));
				for(unsigned int el = 0; el<ThisModelPart.Elements().size(); el++)
				{
					Geometry<Node<3> >& geom = (elem_begin+el)->GetGeometry();
					
					//assigning desired element size
					double h_wish = geom[0].FastGetSolutionStepValue(NODAL_H)
								+ geom[1].FastGetSolutionStepValue(NODAL_H)
								+ geom[2].FastGetSolutionStepValue(NODAL_H);
					h_wish *= 0.33333333333333333;
				
					double wish_area = h_wish*h_wish*0.5;
				
					if(wish_area == 0.0)
						KRATOS_ERROR(std::logic_error,"asking for a 0 area ... impossible to obtain","");
					in.trianglearealist[el] = wish_area;
				}
			}
								
			std::vector< double> preserved(in.numberofpoints);
			
			//read and regenerate the existing mesh ... to generate the boundaries
			if (prescribe_element_size == true) //remesh according to the size function
			{
				char options[] = "rcEe";
				triangulate(options, &in, &mid, &vorout);
				KRATOS_WATCH("sono aaaaa");
				
/*				mid.numberoftriangles = 0;
				free(mid.trianglelist );
				mid.trianglelist  = (int*) NULL;*/
				
				KRATOS_WATCH("sono qui");
				KRATOS_WATCH(mid.numberoftriangles);
				KRATOS_WATCH(mid.numberofpoints);
				KRATOS_WATCH(mid.numberofsegments);
				KRATOS_WATCH(mid.numberofholes);
				
				//save the wish size for all the nodes
				std::vector< double> h_sizes(in.numberofpoints);
				for(int i = 0; i<in.numberofpoints; i++)
				{
					preserved[i] = true;
					h_sizes[i] = (nodes_begin+i)->FastGetSolutionStepValue(NODAL_H);
				}
				KRATOS_WATCH( "ffff ");
				
				
// preserved[ 14  ] = false;
				//coarsen the mesh
				//step2 - coarsen the mesh
// 				for(int edge_it = 0; edge_it < mid.numberofedges; edge_it++) //looping over all edges 
// 				{
// 					int base = edge_it * 2;
// 					
// 					int index1 = mid.edgelist[base]  - 1;
// 					int index2 = mid.edgelist[base+1] -1;
// 					
// 					//if one of the two vertexes is marked for removal don't do anything more
// 					if( preserved[ index1  ] == true && preserved[ index2  ] == true ) 
// 					{
// 						double dx = mid.pointlist[index1*2]-mid.pointlist[index2*2];
// 						double dy = mid.pointlist[index1*2+1]-mid.pointlist[index2*2+1];
// 						double edge_lenght2 = dx*dx + dy*dy;
// 						
// 						double& h1 = h_sizes[index1];
// 						double& h2 = h_sizes[index2];
// 						double wish_size = 0.5*( h1 + h2 );
// 						wish_size *= wish_size;
// 						
// 						if( edge_lenght2  < 0.25*wish_size) //a node should be removed
// 						{
// 							if( h1 < h2) //we should remove the first
// 							{
// 								//if the first is not boundary remove it else try to remove the second
// 								if( mid.pointmarkerlist[index1] != 1)
// 									preserved[ index1  ] = false;
// 								else if( mid.pointmarkerlist[index2] != 1)
// 									preserved[ index2  ] = false;
// 							}
// 							else //trying to remove the second
// 							{
// 								if( mid.pointmarkerlist[index2] != 1)
// 									preserved[ index2  ] = false;
// 								else if( mid.pointmarkerlist[index1] != 1)
// 									preserved[ index1  ] = false;
// 							}
// 							KRATOS_WATCH("removing node");
// 									
// 						}
// 						
// 					}
// 				}
				
				
				KRATOS_WATCH( "eeeee ");
				
				//prepare a compacted list
				std::vector< int > new_ids(in.numberofpoints);
				int counter = 1;
				for(int i=0; i<in.numberofpoints; i++)
					if(preserved[ i ] == true)
						new_ids[i] = counter++;
					else
						new_ids[i] = -1000; //we want to break the memory to get a segfault
				int preserved_nodes = counter-1;
				KRATOS_WATCH( in.numberofpoints );
				KRATOS_WATCH( mid.numberofpoints );
				KRATOS_WATCH( preserved_nodes );
				KRATOS_WATCH( "ddddd ");
				//changing ids in the segment list
				for( int i=0; i<mid.numberofsegments; i++)
				{
					int base = i*2;
/*					KRATOS_WATCH(base);
					KRATOS_WATCH(mid.segmentlist[base]);
					KRATOS_WATCH(mid.segmentlist[base+1]);
					KRATOS_WATCH(new_ids[ mid.segmentlist[base]-1 ]);
					KRATOS_WATCH(new_ids[ mid.segmentlist[base+1]-1 ]);*/
					
					mid.segmentlist[base] = new_ids[ mid.segmentlist[base]-1 ]; 
					mid.segmentlist[base+1] = new_ids[ mid.segmentlist[base+1]-1 ]; 
				}				
				KRATOS_WATCH( "ooooo ");
				//step3 - prepare input for final remeshing
				mid.numberofpoints = preserved_nodes;
				free( mid.pointlist );
				mid.pointlist = (REAL *) malloc(mid.numberofpoints * 2 * sizeof(REAL));
				int destination_counter = 0;
				int origin_counter = 0;
				for( int i=0; i<in.numberofpoints; i++)
				{
					if(preserved[ i ] == true)
					{
						mid.pointlist[ destination_counter++ ] = in.pointlist[ origin_counter++ ]; //X
						mid.pointlist[ destination_counter++ ] = in.pointlist[ origin_counter++ ]; //Y
					}
					else
					{
						origin_counter+=2;
					}
				}
				mid.numberofpointattributes = 0	;
				mid.numberofedges = 0;
				free(mid.edgelist);
				mid.edgelist = (int*) NULL;
				KRATOS_WATCH( "aoaoaoaoa ");		
				
				
				//free the memory used in the first step
				clean_triangulateio(in);
				clean_triangulateio(vorout);
				initialize_triangulateio(vorout);
				
				//char regen_options[] = "pqY";
				char regen_options[] = "pY"; //not introducing any new point
				triangulate(regen_options, &mid, &out, &vorout);
				
				KRATOS_WATCH("sono bbb");
				KRATOS_WATCH(out.numberoftriangles);
				KRATOS_WATCH(out.numberofpoints);
				KRATOS_WATCH(out.numberofsegments);
				KRATOS_WATCH(out.numberofholes);
				
				
			 	
			}
			else //mesh just to ensure quality
			{
// 				KRATOS_ERROR(std::logic_error,"we should rework this","");
				char options[] = "rqYB";
				triangulate(options, &in, &out, &vorout);
			}
						
			
			//******************* removing the nodes *************************
			unsigned last_id_counter = (ThisModelPart.Nodes().end()-1)->Id();
			last_id_counter += 10000000;
			int removed = 0;
			nodes_begin = ThisModelPart.NodesBegin();
			for( unsigned int i = 0; i<ThisModelPart.Nodes().size(); i++)
			{
				if(preserved[i] == false)
				{
					(nodes_begin+i)->Id() += last_id_counter;
					removed++;
				}
			}
			
			(ThisModelPart.Nodes()).Sort();
			
			KRATOS_WATCH(removed);
			
			//int final_size = ThisModelPart.Nodes().size() - removed;
			(ThisModelPart.Nodes()).erase( ThisModelPart.NodesEnd()-removed,ThisModelPart.NodesEnd());
			
			KRATOS_WATCH( ThisModelPart.Nodes().size() );
			
			
			//
			
			//free the memory used in the first step
// 			clean_triangulateio(in);
			
											
			//*******************************************************************
			//properties to be used in the generation
			Properties::Pointer properties = ThisModelPart.GetMesh().pGetProperties(1);

			//generate Kratos triangle
			int el_number = out.numberoftriangles;
			
			//defintions for spatial search
			typedef Node<3> PointType;
			typedef Node<3>::Pointer PointPointerType;
			typedef std::vector<PointType::Pointer>           PointVector;
			typedef std::vector<PointType::Pointer>::iterator PointIterator;
			typedef std::vector<double>               DistanceVector;
			typedef std::vector<double>::iterator     DistanceIterator;

			typedef Bucket<3, PointType, PointVector, PointPointerType, PointIterator, DistanceIterator > BucketType;

			typedef Tree< KDTreePartition<BucketType> > kd_tree; //Kdtree;

			//creating an auxiliary list for the new nodes 
			PointVector list_of_new_nodes;
			
			//node to get the DOFs from					
			Node<3>::DofsContainerType& reference_dofs = (ThisModelPart.NodesBegin())->GetDofs();

// nodes_begin = ThisModelPart.NodesBegin();
// for(int i = 0; i<mid.numberofpoints;i++)
// {
// 	int base = i*2;
// 	out.pointlist[base];
// 	KRATOS_WATCH((nodes_begin+i)->Id()	);
// 	KRATOS_WATCH((nodes_begin+i)->X()	);
// 	KRATOS_WATCH((nodes_begin+i)->Y()	);
// 	KRATOS_WATCH(out.pointlist[base]	);
// 	KRATOS_WATCH(out.pointlist[base+1]	);
// }
				

			//generate list of nodes
			double z = 0.0;
			int original_size = mid.numberofpoints;
			if( out.numberofpoints > original_size  )
			{
				for(int i = original_size; i < out.numberofpoints; i++)
				{
					int id = i+1;
					int base = i*2;
					double& x = out.pointlist[base];		
					double& y = out.pointlist[base+1];
					
					//create new node with the x,y,z corresponding and the attribute list given
					Node<3>::Pointer pnode = ThisModelPart.CreateNewNode(id,x,y,z);

					//putting the new node also in an auxiliary list
					list_of_new_nodes.push_back( pnode );
					
					std::cout << "new node id = " << pnode->Id() << std::endl;
					//generating the dofs
					for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin();    iii != reference_dofs.end(); iii++)
					{
						Node<3>::DofType& rDof = *iii;
						Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );
						
						(p_new_dof)->FreeDof();
					}
					
				}
			}	
			std::cout << out.numberofpoints - original_size << " nodes added" << std::endl;
			
			if(out.numberofpoints - original_size > 0) //if we added points
			{
				//create a spatial database with the list of new nodes
				unsigned int bucket_size = 20;
				kd_tree  nodes_tree2(list_of_new_nodes.begin(),list_of_new_nodes.end(),bucket_size);
							
				//work arrays
				Node<3> work_point(0,0.0,0.0,0.0);
				unsigned int MaximumNumberOfResults = 100;
				PointVector Results(MaximumNumberOfResults);
				DistanceVector ResultsDistances(MaximumNumberOfResults);
				array_1d<double,3> N;
				int step_data_size = ThisModelPart.GetNodalSolutionStepDataSize();
				
				//loop over all of the elements in the "old" list to perform the interpolation
				for( ModelPart::ElementsContainerType::iterator el_it = ThisModelPart.ElementsBegin();
							el_it != ThisModelPart.ElementsEnd(); el_it++)
				{
					Geometry<Node<3> >&geom = el_it->GetGeometry();
					
					//find the center and "radius" of the element
					double xc,  yc, radius;	
					CalculateCenterAndSearchRadius( geom[0].X(), geom[0].Y(),
									geom[1].X(), geom[1].Y(),
									geom[2].X(), geom[2].Y(),
									xc,yc,radius);
									
					//find all of the new nodes within the radius
					int number_of_points_in_radius;
					work_point.X() = xc; work_point.Y() = yc;
	/*				KRATOS_WATCH(xc);				
					KRATOS_WATCH(yc);				
					KRATOS_WATCH(radius);				
					KRATOS_WATCH(work_point);*/
					number_of_points_in_radius = nodes_tree2.SearchInRadius(work_point, radius, Results.begin(),
							ResultsDistances.begin(),  MaximumNumberOfResults);
					
	
					//check if inside and eventually interpolate
					for( PointIterator it_found = Results.begin(); it_found != Results.begin() + number_of_points_in_radius; it_found++)
					{
						bool is_inside = false;
						is_inside = CalculatePosition(geom[0].X(), geom[0].Y(),
										geom[1].X(), geom[1].Y(),
										geom[2].X(), geom[2].Y(),
										(*it_found)->X(),(*it_found)->Y(),N);
						
						//KRATOS_WATCH(N);
						
						//if the node falls inside the element interpolate
						if(is_inside == true)
						{
							Interpolate(  el_it,  N, step_data_size, *(it_found.base() ) );
// 							KRATOS_WATCH("dddddddddddddddddd");	
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
			
			
			
					
// 			std::cout << "calcolare la pos originale!!!!" << std::endl;
			
			//note that the node list can not be changed starting from here
			ModelPart::NodesContainerType& ModelNodes = ThisModelPart.Nodes();
			//generate kratos elements (conditions are not touched)
			for(int el = 0; el< el_number; el++)
			{
				int id = el + 1;
				int base = el * 3;

				Geometry<Node<3> > temp;
				temp.push_back( *((ModelNodes).find( out.trianglelist[base]).base()	) );
				temp.push_back( *((ModelNodes).find( out.trianglelist[base+1]).base()	) );
				temp.push_back( *((ModelNodes).find( out.trianglelist[base+2]).base()	) );
				
				Element::Pointer p_element = rReferenceElement.Create(id, temp, properties);
				ThisModelPart.Elements().push_back(p_element);
				
				double area = CalculateVol(temp[0].X(),temp[0].Y(),temp[1].X(),temp[1].Y(),temp[2].X(),temp[2].Y() );
				if( area < 0.0 )
				{
					KRATOS_WATCH( p_element->Id() );
					KRATOS_WATCH( area  );
				}
			}
			
			

			
			
			



/*			free( in.pointmarkerlist);
			free( in.pointmarkerlist );
			free( in.regionlist );*/
	


//			free( out.pointattributelist );
//			free( out.trianglelist );
// 			free( out.triangleattributelist );
//			free( out.neighborlist );

 			clean_triangulateio(mid);
			clean_triangulateio(out);
			clean_triangulateio(vorout);

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
			free(tr.pointlist );
			free(tr.pointattributelist );
			free(tr.pointmarkerlist   );
			free(tr.trianglelist  );
			free(tr.triangleattributelist );
			free(tr.trianglearealist );
			free(tr.neighborlist   );
			free(tr.segmentlist    );
			free(tr.segmentmarkerlist   );
			free(tr.holelist      );
			free(tr.regionlist  );
			free(tr.edgelist   );
			free(tr.edgemarkerlist   );
			free(tr.normlist  );
		};
		
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
			if(area == 0.0)
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
			
			if(N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0) //if the xc yc is inside the triangle return true
				return true;
			
			return false;
		}	
			
		void Interpolate( ModelPart::ElementsContainerType::iterator el_it, const array_1d<double,3>& N, 
				  int step_data_size,
      				Node<3>::Pointer pnode)
		{
			Geometry< Node<3> >& geom = el_it->GetGeometry();
			
			unsigned int buffer_size = pnode->GetBufferSize();
			
			for(unsigned int step = 0; step<buffer_size; step++)
			{
						//getting the data of the solution step
				double* step_data = (pnode)->SolutionStepData().Data(step);
				
				double* node0_data = geom[0].SolutionStepData().Data(step);
				double* node1_data = geom[1].SolutionStepData().Data(step);
				double* node2_data = geom[2].SolutionStepData().Data(step);
					
				//copying this data in the position of the vector we are interested in
				for(int j= 0; j< step_data_size; j++)
				{
					step_data[j] = N[0]*node0_data[j] + N[1]*node1_data[j] + N[2]*node2_data[j];
				}						
			}				
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
		TriGenCDTrefine& operator=(TriGenCDTrefine const& rOther);


		///@}    

	}; // Class TriGenCDTrefine 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		TriGenCDTrefine& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const TriGenCDTrefine& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_TRIGEN_REFINE_H_INCLUDED  defined 


