//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: antonia $
//   Date:                $Date: 2012-03-08 08:56:42 $
//
//

#if !defined(KRATOS_BINBASED_SEEPAGE_COUPLED_MAPPING )
#define  KRATOS_BINBASED_SEEPAGE_COUPLED_MAPPING

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

#include "PFEM_application.h"

//Database includes
#include "spatial_containers/spatial_containers.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/binbased_nodes_in_element_locator.h"

//iNCLUDE THE DRAG UTILITIED TO CALCULATE THE SEEPAGE DRAG	
// #include "custom_utilities/drag_utilities.h"

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

	/// This class allows the interpolation between non-matching meshes in 2D and 3D in the case of the specific application of coupling a
	///seepage problem with a granular structural response.
	/** @author  Antonia Larese De Tetto <antoldt@cimne.upc.edu>
	* 
	* For every node of the destination model part it is checked in which element of the origin model part it is
	* contained and a linear interpolation is performed
	*
	* The data structure used by default is a bin, 
	*
	* For a more general tool that allows the mapping between 2 and 3D non-matching meshes, please see /kratos/applications/MeshingApplication/custom_utilities/projection.h
	*/

	//class BinBasedSeepageCoupledMapping
	template<std::size_t TDim >
	class BinBasedSeepageCoupledMapping
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of BinBasedSeepageCoupledMapping
		KRATOS_CLASS_POINTER_DEFINITION(BinBasedSeepageCoupledMapping<TDim >);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		BinBasedSeepageCoupledMapping() {} //

		/// Destructor.
		virtual ~BinBasedSeepageCoupledMapping(){}


		///@}
		///@name Operators 
		///@{


		///@}
		///@name Operations
		///@{


		//If you want to pass only one variable
		//**********************************************************************
		//**********************************************************************
		/// Interpolate from the fixed to the moving mesh
		/**
		  * @param rOrigin_ModelPart: the model part  all the variable should be taken from
		  * @param rDestination_ModelPart: the destination model part where we want to know the values of the variables
		  * @param bin_of_objects_fixed: precomputed bin of objects (elelments of the fixed mesh). It is to be constructed separately @see binbased_nodes_in_element_locator 

		  */
		// Form fixed to moving model part
// 		template<class TDataType>
		void InterpolationFromFixedMesh(
			ModelPart& rFixed_ModelPart , 
			ModelPart& rMoving_ModelPart,
// 			Variable<TDataType>& rOriginVariable ,
// 			Variable<TDataType>& rDestinationVariable,
			BinBasedFastPointLocator<TDim>& bin_of_objects_fixed
		)
		{

		    KRATOS_TRY
 		    KRATOS_WATCH("Interpolate From Fixed Mesh*************************************")
		    //Clear all the variables to be mapped 
		    for(ModelPart::NodesContainerType::iterator node_it = rMoving_ModelPart.NodesBegin();
					    node_it != rMoving_ModelPart.NodesEnd(); ++node_it)
		    {

			    ClearVariables(node_it, WATER_PRESSURE);	
    			    ClearVariables(node_it, SEEPAGE_DRAG);			
			    ClearVariables(node_it, PRESS_PROJ);			
		    }		    
		    


		    array_1d<double, TDim + 1 > N;
		    const int max_results = 10000;
		    typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);			
		    const int nparticles = rMoving_ModelPart.Nodes().size();
			
		    #pragma omp parallel for firstprivate(results,N)
		    for (int i = 0; i < nparticles; i++)
		    {
			ModelPart::NodesContainerType::iterator iparticle = rMoving_ModelPart.NodesBegin() + i;
			Node < 3 > ::Pointer pparticle = *(iparticle.base());
			typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
			Element::Pointer pelement;
			// look for the fixed element in which the moving node falls into
			bool is_found = bin_of_objects_fixed.FindPointOnMesh(pparticle->Coordinates(), N, pelement, result_begin, max_results);
			//interpolate the variables	
			
			if (is_found == true)
			{
				//Interpolate(  el_it,  N, *it_found , rOriginVariable , rDestinationVariable  );
				Interpolate(  pelement,  N, pparticle, PRESSURE , WATER_PRESSURE  );
				Interpolate(  pelement,  N, pparticle, SEEPAGE_DRAG , SEEPAGE_DRAG  );
				Interpolate(  pelement,  N, pparticle, PRESS_PROJ , PRESS_PROJ  );
			}
		    }
					
   
		    
		    
		    KRATOS_CATCH("")
		}
				
				
	
		/// Interpolate form the moving  to the fixed mesh
		/**
		  * @param rOrigin_ModelPart: the model part  all the variable should be taken from
		  * @param rDestination_ModelPart: the destination model part where we want to know the values of the variables
		  * @param bin_of_nodes_fixed: precomputed bin of nodes of the fixed mesh. It is to be constructed separately @see binbased_nodes_in_element_locator 
		  */
		// From moving to fixed model part
// 		template<class TDataType>
		void InterpolationFromMovingMesh(
			ModelPart& rMoving_ModelPart , 
			ModelPart& rFixed_ModelPart,
// 			Variable<TDataType>& rMovingDomainVariable ,
// 			Variable<TDataType>& rFixedDomainVariable,
			BinBasedNodesInElementLocator<TDim>& bin_of_nodes_fixed //this is a bin of objects which contains the FIXED model part
			
			)
		{
		    KRATOS_TRY

//  		    KRATOS_WATCH("Transfer From Moving Mesh*************************************")
// 		    if (rMoving_ModelPart.NodesBegin()->SolutionStepsDataHas(rMovingDomainVariable) == false)
// 			KRATOS_ERROR(std::logic_error, "Add  MovingDomain VARIABLE!!!!!! ERROR", "");
// 		    if (rFixed_ModelPart.NodesBegin()->SolutionStepsDataHas(rFixedDomainVariable) == false)
// 			KRATOS_ERROR(std::logic_error, "Add  FixedDomain VARIABLE!!!!!! ERROR", "");

		    
		    //clearing all the mapped variables
		    for(ModelPart::NodesContainerType::iterator node_it = rFixed_ModelPart.NodesBegin();
					    node_it != rFixed_ModelPart.NodesEnd(); ++node_it)
		    {
			    ClearVariables(node_it, POROSITY);					    
			    ClearVariables(node_it, DIAMETER);
			    ClearVariables(node_it, STRUCTURE_VELOCITY);
		    }		    
		    
		    //defintions for spatial search
		    typedef typename BinBasedNodesInElementLocator<TDim>::PointType PointType;
		    typedef typename PointType::Pointer PointTypePointer;
		    typedef typename BinBasedNodesInElementLocator<TDim>::PointIterator PointIterator;
		    typedef typename BinBasedNodesInElementLocator<TDim>::DistanceIterator DistanceIterator;
		    typedef typename BinBasedNodesInElementLocator<TDim>::PointVector PointVector;
		    typedef typename BinBasedNodesInElementLocator<TDim>::DistanceVector DistanceVector;
		    
		    const unsigned int max_results = 10000;
		    Matrix Nmat(max_results,TDim+1);
		    boost::numeric::ublas::vector<int> positions(max_results);
		    PointVector work_results(max_results);
		    DistanceVector work_distances(max_results);
		    Node<3> work_point(0,0.0,0.0,0.0);
		    
		    for(ModelPart::ElementsContainerType::iterator elem_it = rMoving_ModelPart.ElementsBegin(); elem_it != rMoving_ModelPart.ElementsEnd(); ++elem_it)
		    {
		      //look for all the fixed nodes in a moving element
		      unsigned int nfound = bin_of_nodes_fixed.FindNodesInElement(*(elem_it.base()), positions, Nmat, max_results, work_results.begin(), work_distances.begin(), work_point);        			
			for(unsigned int k=0; k<nfound; k++)
			{
			    PointIterator it = work_results.begin() + positions[k];
			    
			    
			    array_1d<double,TDim+1> N = row(Nmat,k);
			    Interpolate(  *(elem_it.base()),  N, *(it.base()), POROSITY , POROSITY);
			    Interpolate(  *(elem_it.base()),  N, *(it.base()), DIAMETER , DIAMETER);
			    Interpolate(  *(elem_it.base()),  N, *(it.base()), VELOCITY , STRUCTURE_VELOCITY);
			}

		     }

		    for(ModelPart::NodesContainerType::iterator node_it = rFixed_ModelPart.NodesBegin();
					    node_it != rFixed_ModelPart.NodesEnd(); ++node_it)
		    {
		      Node < 3 > ::Pointer pnode = *(node_it.base());
		      double& porosity = pnode->FastGetSolutionStepValue(POROSITY);
		      double& diameter = pnode->FastGetSolutionStepValue(DIAMETER);
		      if(porosity == 0.0){porosity = 1.0;}
		      if(diameter == 0.0){diameter = 1.0;}
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



		inline void CalculateCenterAndSearchRadius(Geometry<Node<3> >&geom,
						double& xc, double& yc, double& zc, double& R, array_1d<double,3>& N		
					      )
		{
			double x0 = geom[0].X();double  y0 = geom[0].Y(); 
			double x1 = geom[1].X();double  y1 = geom[1].Y(); 
			double x2 = geom[2].X();double  y2 = geom[2].Y(); 
			

			xc = 0.3333333333333333333*(x0+x1+x2);
			yc = 0.3333333333333333333*(y0+y1+y2);
			zc = 0.0;

			double R1 = (xc-x0)*(xc-x0) + (yc-y0)*(yc-y0);
			double R2 = (xc-x1)*(xc-x1) + (yc-y1)*(yc-y1);
			double R3 = (xc-x2)*(xc-x2) + (yc-y2)*(yc-y2);
			
			R = R1;
			if(R2 > R) R = R2;
			if(R3 > R) R = R3;
			
			R = 1.01 * sqrt(R);
		}
		//***************************************
		//***************************************
		inline void CalculateCenterAndSearchRadius(Geometry<Node<3> >&geom,
						double& xc, double& yc, double& zc, double& R, array_1d<double,4>& N		

					      )
		{
			double x0 = geom[0].X();double  y0 = geom[0].Y();double  z0 = geom[0].Z();
			double x1 = geom[1].X();double  y1 = geom[1].Y();double  z1 = geom[1].Z();
			double x2 = geom[2].X();double  y2 = geom[2].Y();double  z2 = geom[2].Z();
			double x3 = geom[3].X();double  y3 = geom[3].Y();double  z3 = geom[3].Z();	


			xc = 0.25*(x0+x1+x2+x3);
			yc = 0.25*(y0+y1+y2+y3);
			zc = 0.25*(z0+z1+z2+z3);			 

			double R1 = (xc-x0)*(xc-x0) + (yc-y0)*(yc-y0) + (zc-z0)*(zc-z0);
			double R2 = (xc-x1)*(xc-x1) + (yc-y1)*(yc-y1) + (zc-z1)*(zc-z1);
			double R3 = (xc-x2)*(xc-x2) + (yc-y2)*(yc-y2) + (zc-z2)*(zc-z2);
			double R4 = (xc-x3)*(xc-x3) + (yc-y3)*(yc-y3) + (zc-z3)*(zc-z3);
			
			R = R1;
			if(R2 > R) R = R2;
			if(R3 > R) R = R3;
			if(R4 > R) R = R4;
			
			R = sqrt(R);
		}
		//***************************************
		//***************************************
		inline double CalculateVol(	const double x0, const double y0,
						const double x1, const double y1,
						const double x2, const double y2
					  )
		{
			return 0.5*( (x1-x0)*(y2-y0)- (y1-y0)*(x2-x0) );
		}
		//***************************************
		//***************************************
		inline double CalculateVol(	const double x0, const double y0, const double z0,
						const double x1, const double y1, const double z1,
						const double x2, const double y2, const double z2,
						const double x3, const double y3, const double z3
					  )
		{
			double x10 = x1 - x0;
			double y10 = y1 - y0;
			double z10 = z1 - z0;

			double x20 = x2 - x0;
			double y20 = y2 - y0;
			double z20 = z2 - z0;

			double x30 = x3 - x0;
			double y30 = y3 - y0;
			double z30 = z3 - z0;

			double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;
			return  detJ*0.1666666666666666666667;
			
			//return 0.5*( (x1-x0)*(y2-y0)- (y1-y0)*(x2-x0) );
		}
		//***************************************
		//***************************************
		inline bool CalculatePosition(	Geometry<Node<3> >&geom,
						const double xc, const double yc, const double zc,
						array_1d<double,3>& N		
					  )
		{
			double x0 = geom[0].X();double  y0 = geom[0].Y(); 
			double x1 = geom[1].X();double  y1 = geom[1].Y(); 
			double x2 = geom[2].X();double  y2 = geom[2].Y(); 
			
			double area = CalculateVol(x0,y0,x1,y1,x2,y2);
			double inv_area = 0.0;
			if(area == 0.0)
			  {

// 				KRATOS_ERROR(std::logic_error,"element with zero area found","");
				//The interpolated node will not be inside an elemente with zero area
				return false;
				
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
		
		//***************************************
		//***************************************

		inline bool CalculatePosition(	Geometry<Node<3> >&geom,
						const double xc, const double yc, const double zc,
						array_1d<double,4>& N		
					  )
		{ 
			
			double x0 = geom[0].X();double  y0 = geom[0].Y();double  z0 = geom[0].Z();
			double x1 = geom[1].X();double  y1 = geom[1].Y();double  z1 = geom[1].Z();
			double x2 = geom[2].X();double  y2 = geom[2].Y();double  z2 = geom[2].Z();
			double x3 = geom[3].X();double  y3 = geom[3].Y();double  z3 = geom[3].Z();	
	
			double vol = CalculateVol(x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3);

			double inv_vol = 0.0;
			if(vol < 0.0000000000001)
			  {

// 				KRATOS_ERROR(std::logic_error,"element with zero vol found","");
				//The interpolated node will not be inside an elemente with zero volume
				return false;
// 				KRATOS_WATCH("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
			  }
			else
			  {
				inv_vol = 1.0 / vol;
			  }

			  N[0] = CalculateVol(x1,y1,z1,x3,y3,z3,x2,y2,z2,xc,yc,zc) * inv_vol;
			  N[1] = CalculateVol(x3,y3,z3,x0,y0,z0,x2,y2,z2,xc,yc,zc) * inv_vol;
			  N[2] = CalculateVol(x3,y3,z3,x1,y1,z1,x0,y0,z0,xc,yc,zc) * inv_vol;
			  N[3] = CalculateVol(x0,y0,z0,x1,y1,z1,x2,y2,z2,xc,yc,zc) * inv_vol;

						
			if(N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[3] >=0.0 &&
			  N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0 && N[3] <=1.0)
			//if the xc yc zc is inside the tetrahedron return true
				return true;
			
			return false;
		}	
//el_it		     	Element iterator
//N			Shape functions
//step_data_size	
//pnode			pointer to the node


		//projecting an array1D 2Dversion
		void Interpolate( 
				Element::Pointer el_it, 
				const array_1d<double,3>& N, 
				Node<3>::Pointer pnode,
				Variable<array_1d<double,3> >& rOriginVariable,
				Variable<array_1d<double,3> >& rDestinationVariable)
		{
// 		  		  KRATOS_ERROR(std::logic_error,"INTERPOLATE ARRAY 2D","")

			//Geometry element of the rOrigin_ModelPart
			Geometry< Node<3> >& geom = el_it->GetGeometry();
			
			unsigned int buffer_size = pnode->GetBufferSize();

			for(unsigned int step = 0; step<buffer_size; step++)
			{
				//getting the data of the solution step
				array_1d<double,3>& step_data = (pnode)->FastGetSolutionStepValue(rDestinationVariable , step);
				//Reference or no reference???//CANCELLA
				const array_1d<double,3>& node0_data = geom[0].FastGetSolutionStepValue(rOriginVariable , step);
				const array_1d<double,3>& node1_data = geom[1].FastGetSolutionStepValue(rOriginVariable , step);
				const array_1d<double,3>& node2_data = geom[2].FastGetSolutionStepValue(rOriginVariable , step);
					
				//copying this data in the position of the vector we are interested in
				for(unsigned int j= 0; j< TDim; j++)
				{
					step_data[j] = N[0]*node0_data[j] + N[1]*node1_data[j] + N[2]*node2_data[j];
				}						
			}	
// 			pnode->GetValue(IS_VISITED) = 1.0;
			
		}

		//projecting an array1D 3Dversion
		void Interpolate( 
				Element::Pointer el_it, 
				const array_1d<double,4>& N, 
				Node<3>::Pointer pnode,
				Variable<array_1d<double,3> >& rOriginVariable,
				Variable<array_1d<double,3> >& rDestinationVariable)

		{
// 		  	KRATOS_ERROR(std::logic_error,"INTERPOLATE ARRAY 3D","")

			//Geometry element of the rOrigin_ModelPart
			Geometry< Node<3> >& geom = el_it->GetGeometry();
			
			unsigned int buffer_size = pnode->GetBufferSize();

			for(unsigned int step = 0; step<buffer_size; step++)
			{
				//getting the data of the solution step
				array_1d<double,3>& step_data = (pnode)->FastGetSolutionStepValue(rDestinationVariable , step);
				//Reference or no reference???//CANCELLA
				const array_1d<double,3>& node0_data = geom[0].FastGetSolutionStepValue(rOriginVariable , step);
				const array_1d<double,3>& node1_data = geom[1].FastGetSolutionStepValue(rOriginVariable , step);
				const array_1d<double,3>& node2_data = geom[2].FastGetSolutionStepValue(rOriginVariable , step);
				const array_1d<double,3>& node3_data = geom[3].FastGetSolutionStepValue(rOriginVariable , step);

				//copying this data in the position of the vector we are interested in
				for(unsigned int j= 0; j< TDim; j++)
				{
					step_data[j] = N[0]*node0_data[j] + N[1]*node1_data[j] + N[2]*node2_data[j] + N[3]*node3_data[j];
				}						
			}
// 			pnode->GetValue(IS_VISITED) = 1.0;
				
		}
		//projecting a scalar 2Dversion
		void Interpolate( 
				Element::Pointer el_it, 
				const array_1d<double,3>& N, 
				Node<3>::Pointer pnode,
				Variable<double>& rOriginVariable,
				Variable<double>& rDestinationVariable)
		{
// 		  	  KRATOS_ERROR(std::logic_error,"INTERPOLATE SCALAR 2D","")

			//Geometry element of the rOrigin_ModelPart
			Geometry< Node<3> >& geom = el_it->GetGeometry();
			
			unsigned int buffer_size = pnode->GetBufferSize();
		//facendo un loop sugli step temporali step_data come salva i dati al passo anteriore? Cioś dove passiamo l'informazione ai nodi???
			for(unsigned int step = 0; step<buffer_size; step++)
			{
				//getting the data of the solution step
				double& step_data = (pnode)->FastGetSolutionStepValue(rDestinationVariable , step);
				//Reference or no reference???//CANCELLA
				const double node0_data = geom[0].FastGetSolutionStepValue(rOriginVariable , step);
				const double node1_data = geom[1].FastGetSolutionStepValue(rOriginVariable , step);
				const double node2_data = geom[2].FastGetSolutionStepValue(rOriginVariable , step);
					
				//copying this data in the position of the vector we are interested in
				
				step_data = N[0]*node0_data + N[1]*node1_data + N[2]*node2_data;
										
			}	
// // 			pnode->GetValue(IS_VISITED) = 1.0;
			
		}
		//projecting a scalar 3Dversion
		void Interpolate( 
				Element::Pointer el_it, 
				const array_1d<double,4>& N, 
				Node<3>::Pointer pnode,
				Variable<double>& rOriginVariable,
				Variable<double>& rDestinationVariable)
		{
			//Geometry element of the rOrigin_ModelPart
			Geometry< Node<3> >& geom = el_it->GetGeometry();
			
			unsigned int buffer_size = pnode->GetBufferSize();
		//facendo un loop sugli step temporali step_data come salva i dati al passo anteriore? Cioś dove passiamo l'informazione ai nodi???
			for(unsigned int step = 0; step<buffer_size; step++)
			{
				//getting the data of the solution step
				double& step_data = (pnode)->FastGetSolutionStepValue(rDestinationVariable , step);
				//Reference or no reference???//CANCELLA
				const double node0_data = geom[0].FastGetSolutionStepValue(rOriginVariable , step);
				const double node1_data = geom[1].FastGetSolutionStepValue(rOriginVariable , step);
				const double node2_data = geom[2].FastGetSolutionStepValue(rOriginVariable , step);
				const double node3_data = geom[3].FastGetSolutionStepValue(rOriginVariable , step);

				//copying this data in the position of the vector we are interested in
				
				step_data = N[0]*node0_data + N[1]*node1_data + N[2]*node2_data + N[3]*node3_data; 
// 				KRATOS_WATCH(step_data)					
			}	
// 			pnode->GetValue(IS_VISITED) = 1.0;
			
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
		BinBasedSeepageCoupledMapping& operator=(BinBasedSeepageCoupledMapping const& rOther);


		///@}    

	}; // Class BinBasedSeepageCoupledMapping 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 




	/// output stream function
	template<std::size_t TDim>
	inline std::ostream& operator << (std::ostream& rOStream, 
		const BinBasedSeepageCoupledMapping<TDim>& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_BINBASED_SEEPAGE_COUPLED_MAPPING  defined 


