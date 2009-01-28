//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: antonia $
//   Date:                $Date: 2008-10-13 08:56:42 $
//   Revision:            $Revision: 1.5 $
//
//
//README::::look to the key word "VERSION" if you want to find all the points where you have to change something so that you can pass from a kdtree to a bin data search structure;

#if !defined(KRATOS_PROJECTION )
#define  KRATOS_PROJECTION

// /* External includes */
// #include "boost/smart_ptr.hpp"

// System includes
#include <string>
#include <iostream> 
#include <stdlib.h>


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "geometries/triangle_2d_3.h"
#include "utilities/timer.h"

// #include "geometries/tetrahedra_3d_4.h"

#include "meshing_application.h"

// #include "containers/kratos_spacial_search.h"


//Database includes
#include "spatial_containers/spatial_containers.h"


namespace Kratos
{

	template< class T, std::size_t dim >
			class DistanceCalculator
			{
				public:
					double operator()( T const& p1, T const& p2 ){
						double dist = 0.0;
						for( std::size_t i = 0 ; i < dim ; i++){
							double tmp = p1[i] - p2[i];
							dist += tmp*tmp;
						}
						return dist; //square distance because it is easier to work without the square root//
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
	//class MeshTransfer
	template<std::size_t TDim>
	class MeshTransfer  
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of MeshTransfer
		KRATOS_CLASS_POINTER_DEFINITION(MeshTransfer<TDim>);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		MeshTransfer() {} //

		/// Destructor.
		virtual ~MeshTransfer(){}


		///@}
		///@name Operators 
		///@{


		///@}
		///@name Operations
		///@{

		//If you want to pass the whole model part
		//**********************************************************************
		//**********************************************************************
		void Direct_Interpolation(
			ModelPart& rOrigin_ModelPart , 
			ModelPart& rDestination_ModelPart 
			)
 		{
			KRATOS_TRY

 			//**********************************************************************
// 			numofpts_origin = rOrigin_ModelPart.Nodes().size();
// 			numofpts_destination = rDestination_ModelPart.Nodes().size();
// 					
 			//*******************************************************************
			//properties to be used in the generation
			Properties::Pointer properties = rDestination_ModelPart.GetMesh().pGetProperties(1);

			//defintions for spatial search
			typedef Node<3> PointType;
			typedef Node<3>::Pointer PointTypePointer;
			typedef std::vector<PointType::Pointer>           PointVector;
			typedef std::vector<PointType::Pointer>::iterator PointIterator;
			typedef std::vector<double>               DistanceVector;
			typedef std::vector<double>::iterator     DistanceIterator;


			//creating an auxiliary list for the new nodes 
			PointVector list_of_new_nodes;

// 			KRATOS_WATCH("STARTING KDTREE CONSTRUCTION");
			//starting calculating time of construction of the kdtree
			boost::timer kdtree_construction;

			//*************
			// Bucket types
			   typedef Bucket< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
// 			   typedef Bins< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > StaticBins;
// 			   typedef BinsDynamic< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > DynamicBins;
			//*************
			// DynamicBins;
			
			   typedef Tree< KDTreePartition<BucketType> > tree; 		//Kdtree;
// 			   typedef Tree< OCTreePartition<BucketType> > tree; 		//Octree;
// 			   typedef Tree< StaticBins > tree; 		     		//Binstree;
// 			   typedef Tree< KDTreePartition<StaticBins> > tree; 		//KdtreeBins;
// 			   typedef typename KdtreeBins::Partitions SubPartitions;
// 			   typedef Tree< OCTreePartition<StaticBins> > tree; 		//OctreeBins;
/*			
			   typedef Bins< TDim, PointType, stdPointVector> stdBins;
			   typedef Tree< Bins<TDim,PointType,stdPointVector> > tree; 	//stdStaticBins;*/


			for(ModelPart::NodesContainerType::iterator node_it = rDestination_ModelPart.NodesBegin();
						node_it != rDestination_ModelPart.NodesEnd(); ++node_it)
			{
				//PointType::Pointer pnode(new PointType(*node_it));
 				Node<3>::Pointer pnode = *(node_it.base());

				//putting the nodes of the destination_model part in an auxiliary list
				list_of_new_nodes.push_back( pnode );
			}

			std::cout << "kdt constructin time " << kdtree_construction.elapsed() << std::endl;
			//finishing calculating time of construction of the kdtree	
// 			KRATOS_WATCH("FINISHING KDTREE CONSTRUCTION");

			//create a spatial database with the list of new nodes
			unsigned int bucket_size = 20;
			tree nodes_tree(list_of_new_nodes.begin(),list_of_new_nodes.end(),bucket_size);

		
			//work arrays
			Node<3> work_point(0,0.0,0.0,0.0);
			unsigned int MaximumNumberOfResults = 10000;
			PointVector Results(MaximumNumberOfResults);
			DistanceVector ResultsDistances(MaximumNumberOfResults);
			array_1d<double,3> N; //Shape functions vector//
			int step_data_size = rDestination_ModelPart.GetNodalSolutionStepDataSize();

			for(ModelPart::NodesContainerType::iterator node_it = rDestination_ModelPart.NodesBegin();
						node_it != rDestination_ModelPart.NodesEnd(); ++node_it)
			{
				//Setting to zero the whole model part
				Clear(node_it,  step_data_size );
			}
			
			//loop over all of the elements in the "old" list to perform the interpolation
			for( ModelPart::ElementsContainerType::iterator el_it = rOrigin_ModelPart.ElementsBegin();
						el_it != rOrigin_ModelPart.ElementsEnd(); el_it++)
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

				//look between the new nodes which of them is inside the radius of the circumscribed cyrcle
				number_of_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, Results.begin(),
 						ResultsDistances.begin(),  MaximumNumberOfResults);
// 				if((el_it) -> Id()==245)
// 				{
// 				KRATOS_WATCH("Entra en el elemento 245");	
// 				}
// 				if((el_it) -> Id()==307)
// 				{
// 				KRATOS_WATCH("Entra en el elemento 307!!!!!!!!!!!!!");	
// 						KRATOS_WATCH(geom[0].X());
// 						KRATOS_WATCH(geom[0].Y());
// 						KRATOS_WATCH(geom[1].X());
// 						KRATOS_WATCH(geom[1].Y());
// 					std::cout << "velocity correction time " << vel_time.elapsed() << std::endl;	KRATOS_WATCH(geom[2].X());
// 						KRATOS_WATCH(geom[2].Y());
// 				}

				//check if inside 
				for( PointIterator it_found = Results.begin(); it_found != Results.begin() + number_of_points_in_radius; it_found++)
				{	
 				
					bool is_inside = false;
					//once we are sure the node in inside the circle we have to see if it is inside the triangle i.e. if all the Origin element shape functions are >1
					is_inside = CalculatePosition(geom[0].X(), geom[0].Y(),
									geom[1].X(), geom[1].Y(),
									geom[2].X(), geom[2].Y(),
									(*it_found)->X(),(*it_found)->Y(),N);
				
// 					
// 					if((*it_found)->Id() == 20163)
// 					{
// 						KRATOS_WATCH("*********************************  GOOOOOD");
// 						KRATOS_WATCH("NODE 20163");
// 						KRATOS_WATCH(is_inside);
//  						KRATOS_WATCH(geom[0].X());
// 						KRATOS_WATCH(geom[0].Y());
// 						KRATOS_WATCH(geom[1].X());
// 						KRATOS_WATCH(geom[1].Y());
// 						KRATOS_WATCH(geom[2].X());
// 						KRATOS_WATCH(geom[2].Y());
// 						KRATOS_WATCH(N);
// 					}
					//if the node falls inside the element interpolate
					if(is_inside == true)
					{
						//Interpolating all the rVariables of the rOrigin_ModelPart to get their nodal value in the rDestination_ModelPart
						Interpolate(  el_it,  N, step_data_size, *(it_found.base() ) );
						
					}
				}
 			}		
			
			KRATOS_CATCH("")
		}


		//If you want to pass only one variable
		//**********************************************************************
		//**********************************************************************
		template<class TDataType>
		void Direct_Variable_Interpolation(
			ModelPart& rOrigin_ModelPart , 
			ModelPart& rDestination_ModelPart,
			Variable<TDataType>& rVariable 
			)
 		{

			KRATOS_TRY

 			//*******************************************************************
			//properties to be used in the generation
			Properties::Pointer properties = rDestination_ModelPart.GetMesh().pGetProperties(1);

			//defintions for spatial search
			typedef Node<3> PointType;
			typedef Node<3>::Pointer PointTypePointer;
			typedef std::vector<PointType::Pointer>           PointVector;
			typedef std::vector<PointType::Pointer>::iterator PointIterator;
			typedef std::vector<double>               DistanceVector;
			typedef std::vector<double>::iterator     DistanceIterator;


			//creating an auxiliary list for the new nodes 
			PointVector list_of_new_nodes;

// 			KRATOS_WATCH("STARTING KDTREE CONSTRUCTION");
			//starting calculating time of construction of the kdtree
			boost::timer kdtree_construction;

			//*************
			// Bucket types
			   typedef Bucket< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
// 			   typedef Bins< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > StaticBins;
// 			   typedef BinsDynamic< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > DynamicBins;
			//*************
			// DynamicBins;
			
			   typedef Tree< KDTreePartition<BucketType> > tree; 		//Kdtree;
// 			   typedef Tree< OCTreePartition<BucketType> > tree; 		//Octree;
// 			   typedef Tree< StaticBins > tree; 		     		//Binstree;
// 			   typedef Tree< KDTreePartition<StaticBins> > tree; 		//KdtreeBins;
// 			   typedef typename KdtreeBins::Partitions SubPartitions;
// 			   typedef Tree< OCTreePartition<StaticBins> > tree; 		//OctreeBins;
/*			
			   typedef Bins< TDim, PointType, stdPointVector> stdBins;
			   typedef Tree< Bins<TDim,PointType,stdPointVector> > tree; 	//stdStaticBins;*/

			

			for(ModelPart::NodesContainerType::iterator node_it = rDestination_ModelPart.NodesBegin();
						node_it != rDestination_ModelPart.NodesEnd(); ++node_it)
			{

				ClearVariables(node_it, rVariable);			


				//PointType::Pointer pnode(new PointType(*node_it));
				Node<3>::Pointer pnode = *(node_it.base());

				//putting the nodes of the destination_model part in an auxiliary list
				list_of_new_nodes.push_back( pnode );
			}

			std::cout << "kdt constructin time " << kdtree_construction.elapsed() << std::endl;
			//finishing calculating time of construction of the kdtree	
// 			KRATOS_WATCH("FINISHING KDTREE CONSTRUCTION");
			
			//create a spatial database with the list of new nodes
			unsigned int bucket_size = 20;
			tree nodes_tree(list_of_new_nodes.begin(),list_of_new_nodes.end(),bucket_size);
	
			//work arrays
			Node<3> work_point(0,0.0,0.0,0.0);
			unsigned int MaximumNumberOfResults = 10000;
			PointVector Results(MaximumNumberOfResults);
			DistanceVector ResultsDistances(MaximumNumberOfResults);
			array_1d<double,3> N; //Shape functions vector//
			//int step_data_size = rDestination_ModelPart.GetNodalSolutionStepDataSize();
			//unsigned int TDim = 3;

			//loop over all of the elements in the "old" list to perform the interpolation
			for( ModelPart::ElementsContainerType::iterator el_it = rOrigin_ModelPart.ElementsBegin();
						el_it != rOrigin_ModelPart.ElementsEnd(); el_it++)
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

				//look between the new nodes which of them is inside the radius of the circumscribed cyrcle
				number_of_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, Results.begin(),
 						ResultsDistances.begin(),  MaximumNumberOfResults);

				//check if inside 
				for( PointIterator it_found = Results.begin(); it_found != Results.begin() + number_of_points_in_radius; it_found++)
				{	
 				
					bool is_inside = false;
					//once we are sure the node in inside the circle we have to see if it is inside the triangle i.e. if all the Origin element shape functions are >1
					is_inside = CalculatePosition(geom[0].X(), geom[0].Y(),
									geom[1].X(), geom[1].Y(),
									geom[2].X(), geom[2].Y(),
									(*it_found)->X(),(*it_found)->Y(),N);
				
					//if the node falls inside the element interpolate
					if(is_inside == true)
					{//CANCELLA insert the variable TDim
						//Interpolating all the rVariables of the rOrigin_ModelPart to get their nodal value in the rDestination_ModelPart
						Interpolate(  el_it,  N, *(it_found.base() ) , rVariable );
						
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

		/// Turn back information as a stemplate<class T, std::size_t dim> tring.
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
		///@name Protected static Member rVariables 
		///@{ 


		///@} 
		///@name Protected member rVariables 
		///@{ template<class T, std::size_t dim> 


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
		///@name Static Member rVariables 
		///@{ 


		///@} 
		///@name Member rVariables 
		///@{ 
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
			
			R = 1.01 * sqrt(R);
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
		
			
			if(N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <=1.0 && N[1]<= 1.0 && N[2] <= 1.0) //if the xc yc is inside the triangle return true
				return true;
			
			return false;
		}	
		

//el_it		     	Element iterator
//N			Shape functions
//step_data_size	
//pnode			pointer to the node
		//total model part
		void Interpolate( 
				ModelPart::ElementsContainerType::iterator el_it, 
				const array_1d<double,3>& N, 
				int step_data_size,
      				Node<3>::Pointer pnode)
		{
			//Geometry element of the rOrigin_ModelPart
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
		


		//array1D
		 void Interpolate( 
				ModelPart::ElementsContainerType::iterator el_it, 
				const array_1d<double,3>& N, 
      				Node<3>::Pointer pnode,
				Variable<array_1d<double,3> >& rVariable)
		{
			//Geometry element of the rOrigin_ModelPart
			Geometry< Node<3> >& geom = el_it->GetGeometry();
			
			unsigned int buffer_size = pnode->GetBufferSize();

			for(unsigned int step = 0; step<buffer_size; step++)
			{
				//getting the data of the solution step
				array_1d<double,3>& step_data = (pnode)->FastGetSolutionStepValue(rVariable , step);
				//Reference or no reference???//CANCELLA
				array_1d<double,3>& node0_data = geom[0].FastGetSolutionStepValue(rVariable , step);
				array_1d<double,3>& node1_data = geom[1].FastGetSolutionStepValue(rVariable , step);
				array_1d<double,3>& node2_data = geom[2].FastGetSolutionStepValue(rVariable , step);
					
				//copying this data in the position of the vector we are interested in
				for(unsigned int j= 0; j< TDim; j++)
				{
					step_data[j] = N[0]*node0_data[j] + N[1]*node1_data[j] + N[2]*node2_data[j];
				}						
			}				
		}

		//scalar
		void Interpolate( 
				ModelPart::ElementsContainerType::iterator el_it, 
				const array_1d<double,3>& N, 
      				Node<3>::Pointer pnode,
				Variable<double>& rVariable)
		{
			//Geometry element of the rOrigin_ModelPart
			Geometry< Node<3> >& geom = el_it->GetGeometry();
			
			unsigned int buffer_size = pnode->GetBufferSize();
		//facendo un loop sugli step temporali step_data come salva i dati al passo anteriore? Cioś dove passiamo l'informazione ai nodi???
			for(unsigned int step = 0; step<buffer_size; step++)
			{
				//getting the data of the solution step
				double& step_data = (pnode)->FastGetSolutionStepValue(rVariable , step);
				//Reference or no reference???//CANCELLA
				double& node0_data = geom[0].FastGetSolutionStepValue(rVariable , step);
				double& node1_data = geom[1].FastGetSolutionStepValue(rVariable , step);
				double& node2_data = geom[2].FastGetSolutionStepValue(rVariable , step);
					
				//copying this data in the position of the vector we are interested in
				
				step_data = N[0]*node0_data + N[1]*node1_data + N[2]*node2_data;
										
			}				
		}

		inline void Clear(ModelPart::NodesContainerType::iterator node_it,  int step_data_size )	
		{
			unsigned int buffer_size = node_it->GetBufferSize();
			
			for(unsigned int step = 0; step<buffer_size; step++)
			{
				//getting the data of the solution step
				double* step_data = (node_it)->SolutionStepData().Data(step);
					
				//copying this data in the position of the vector we are interested in
				for(int j= 0; j< step_data_size; j++)
				{
					step_data[j] = 0.0;
				}						
			}	

		}

		inline void ClearVariables(ModelPart::NodesContainerType::iterator node_it , Variable<array_1d<double,3> >& rVariable)
		{
			array_1d<double, 3>& Aux_var = node_it->FastGetSolutionStepValue(rVariable, 0);
					
			noalias(Aux_var) = ZeroVector(3);
			
		}	


		inline void ClearVariables(ModelPart::NodesContainerType::iterator node_it,  Variable<double>& rVariable)	
		{
			double& Aux_var = node_it->FastGetSolutionStepValue(rVariable, 0);
					
			Aux_var = 0.0;

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
		MeshTransfer& operator=(MeshTransfer const& rOther);


		///@}    

	}; // Class MeshTransfer 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 




	/// output stream function
	template<std::size_t TDim>
	inline std::ostream& operator << (std::ostream& rOStream, 
		const MeshTransfer<TDim>& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_PROJECTION  defined 


