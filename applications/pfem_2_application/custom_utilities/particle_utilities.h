/*
==============================================================================
KratosTestApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2010
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

The  above  copyright  notice  and  this permission  notice  sKRATOS_WATCH(disp);hall  be
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
//   Last Modified by:    $Author: julio.marti $
//   Date:                $Date: 2007-03-06 10:30:31 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_PARTICLES_UTILITIES_INCLUDED )
#define  KRATOS_PARTICLES_UTILITIES_INCLUDED

#define PRESSURE_ON_EULERIAN_MESH
#define USE_FEW_PARTICLES

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
#include "geometries/tetrahedra_3d_4.h"
#include "pfem_2_application.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/timer.h"
#include "processes/node_erase_process.h"
#include "utilities/binbased_fast_point_locator.h"
//#include "utilities/enrichment_utilities.h"

#include <boost/timer.hpp>
#include "utilities/timer.h"

#ifdef _OPENMP
#include "omp.h"
#endif

namespace Kratos
{
  
  template< class T, std::size_t dim >
    class DistanceCalculator1
    {
    public:
      
      double operator()(T const& p1, T const& p2)
      {
	double dist = 0.0;
        for (std::size_t i = 0; i < dim; i++)
	  {
            double tmp = p1[i] - p2[i];
            dist += tmp*tmp;
	  }
        return dist; //square distance because it is easier to work without the square root//
      }
      
    };
  
  template<std::size_t TDim> class ParticleUtils
    {
    public:
      KRATOS_CLASS_POINTER_DEFINITION(ParticleUtils<TDim>);
      
      
      void EstimateTime(ModelPart& rEulerianModelPart,const double max_dt)
      {
        KRATOS_TRY
	  // KRATOS_THROW_ERROR(std::logic_error,  "NEGATIVE VALUE OF Time step estimated" , "");
	  //initializee dt with max dt
	  //initialize dt with incredible value
	  double /*dt, glob_min_dt,*/ dummy;
	// 	    double h, nu;
        array_1d<double,3> N = ZeroVector(3);
        array_1d<double,3> aux = ZeroVector(3); //dimension = number of nodes
        array_1d<double,3> vel = ZeroVector(3); //dimension = number of nodes
        boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX = ZeroMatrix(3,2);
        array_1d<double,2> ms_vel_gauss = ZeroVector(2); //dimesion coincides with space dimension
	
        //initialize it with given value
	//        glob_min_dt=max_dt;
	
	
	//        dt=0.0;
        for(ModelPart::ElementsContainerType::iterator im = rEulerianModelPart.ElementsBegin() ; im !=rEulerianModelPart.ElementsEnd() ; ++im)
	  {
            GeometryUtils::CalculateGeometryData(im->GetGeometry(),DN_DX,N,dummy);
	    
            double h = sqrt(2.00*dummy);
	    
            array_1d<double,3> const& v = im->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
            ms_vel_gauss[0] = v[0];
            ms_vel_gauss[1] = v[1];
	    
            //direction of the height is stored in the auxilliary vector
            for (unsigned int i=1; i<3; i++)
	      {
                array_1d<double,3> const& vi = im->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
                ms_vel_gauss[0] += vi[0];
                ms_vel_gauss[1] += vi[1];
	      }
            ms_vel_gauss *=0.3333;

            double norm_u = ms_vel_gauss[0]*ms_vel_gauss[0] + ms_vel_gauss[1]*ms_vel_gauss[1];
            norm_u = sqrt(norm_u);
	    
            double courant= norm_u * max_dt / h;
	    
	    
            double& counter = im->GetValue(POISSON_RATIO);
            counter = courant;
	    
	  }
	

        KRATOS_CATCH("");
      }
      
      void VisualizationModelPart(ModelPart& rCompleteModelPart, ModelPart& rEulerianModelPart, ModelPart & rLagrangianModelPart)
      {
        KRATOS_TRY;
	
        rCompleteModelPart.Elements() = rEulerianModelPart.Elements();
        rCompleteModelPart.Nodes() = rEulerianModelPart.Nodes();
	
	
        unsigned int id;
        if(rEulerianModelPart.Nodes().size()!= 0)
	  id = (rEulerianModelPart.Nodes().end() - 1)->Id() + 1;
        else
	  id = 1;
	
        //preallocate the memory needed
        int tot_nodes = rEulerianModelPart.Nodes().size() + rLagrangianModelPart.Nodes().size();
        rCompleteModelPart.Nodes().reserve( tot_nodes );
	
        //note that here we renumber the nodes
        for (ModelPart::NodesContainerType::iterator node_it = rLagrangianModelPart.NodesBegin();
	     node_it != rLagrangianModelPart.NodesEnd(); node_it++)
	  {
	    node_it->SetId(id++);
            rCompleteModelPart.AddNode(*(node_it.base()));
	  }
	
        KRATOS_CATCH("");
      }
      
      
      void TransferToEulerianMesh_Face_Heat_Flux(ModelPart& rEulerianModelPart, ModelPart & rLagrangianModelPart)
      {
        KRATOS_TRY
	  
	  //defintions for spatial search
	  typedef Node < 3 > PointType;
        typedef Node < 3 > ::Pointer PointTypePointer;
        typedef std::vector<PointType::Pointer> PointVector;
        typedef std::vector<PointType::Pointer>::iterator PointIterator;
        typedef std::vector<double> DistanceVector;
        typedef std::vector<double>::iterator DistanceIterator;
	
        //creating an auxiliary list for the new nodes
        PointVector list_of_nodes;
	
        //*************
        // Bucket types
        typedef Bucket< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
	
        typedef Tree< KDTreePartition<BucketType> > tree; //Kdtree;
	
	
        //starting calculating time of construction of the kdtree
        boost::timer kdtree_construction;
	
        for (ModelPart::NodesContainerType::iterator node_it = rLagrangianModelPart.NodesBegin();
	     node_it != rLagrangianModelPart.NodesEnd(); ++node_it)
	  {
	    PointTypePointer pnode = *(node_it.base());
	  
	    //putting the nodes of the destination_model part in an auxiliary list
	    list_of_nodes.push_back(pnode);
	  }
	
        std::cout << "kdt constructin time " << kdtree_construction.elapsed() << std::endl;
	
        //create a spatial database with the list of new nodes
        unsigned int bucket_size = 20;
        tree nodes_tree(list_of_nodes.begin(), list_of_nodes.end(), bucket_size);
	
        //work arrays
        Node < 3 > work_point(0, 0.0, 0.0, 0.0);
        unsigned int MaximumNumberOfResults = 10000;
        PointVector Results(MaximumNumberOfResults);
        DistanceVector SquaredResultsDistances(MaximumNumberOfResults);
	
	
        if (rEulerianModelPart.NodesBegin()->SolutionStepsDataHas(NODAL_H) == false)
	  KRATOS_THROW_ERROR(std::logic_error, "Add  ----NODAL_H---- variable!!!!!! ERROR", "");
	
        double sigma = 0.0;
        if (TDim == 2)
	  sigma = 10.0 / (7.0 * 3.1415926);
        else
	  sigma = 1.0 / 3.1415926;
	
	for (ModelPart::NodesContainerType::iterator node_it = rEulerianModelPart.NodesBegin(); node_it != rEulerianModelPart.NodesEnd(); node_it++)
	  {
	    if((node_it)->FastGetSolutionStepValue(IS_FREE_SURFACE)==1 and (node_it)->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET)==0)
	      { //IS_FREE_SURFACE
		work_point.X() = node_it->X();
		work_point.Y() = node_it->Y();
		work_point.Z() = node_it->Z();
		
		double radius = 1.5 * node_it->FastGetSolutionStepValue(NODAL_H);
		
		//find all of the new nodes within the radius
		int number_of_points_in_radius;
		
		//look between the new nodes which of them is inside the radius of the circumscribed cyrcle
		number_of_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, Results.begin(), SquaredResultsDistances.begin(), MaximumNumberOfResults);
		
		if (number_of_points_in_radius > 0)
		  {
		    
		    double& temperature = (node_it)->FastGetSolutionStepValue(FACE_HEAT_FLUX);
		    //double temperature=0.0;
		    
		    double temperature_aux = 0.0;
		    
		    double tot_weight = 0.0;
		    
		    for (int k = 0; k < number_of_points_in_radius; k++)
		      {
			double distance = sqrt(*(SquaredResultsDistances.begin() + k));
			double weight = SPHCubicKernel(sigma, distance, radius);
			
			PointIterator it_found = Results.begin() + k;
			
			if((*it_found)->FastGetSolutionStepValue(IS_INTERFACE)==1) //MATERIAL_VARIABLE
			  {
			    temperature_aux += weight * (*it_found)->FastGetSolutionStepValue(INCIDENT_RADIATION_FUNCTION);//);//FACE_HEAT_FLUX
			    tot_weight += weight;
			  }
		      }
		    if(tot_weight>0.0)
		      {
			temperature_aux /= tot_weight;
			
			temperature +=(0.5 * temperature_aux * 1.00); //1.5 //1.25
		      }
		  }
	      }
            else
	      {
	      }
	  }
        KRATOS_CATCH("")
	  }
      
      void TransferToEulerianMesh(ModelPart& rEulerianModelPart, ModelPart & rLagrangianModelPart)
      {
        KRATOS_TRY
	  
	  //defintions for spatial search
	  typedef Node < 3 > PointType;
        typedef Node < 3 > ::Pointer PointTypePointer;
        typedef std::vector<PointType::Pointer> PointVector;
        typedef std::vector<PointType::Pointer>::iterator PointIterator;
        typedef std::vector<double> DistanceVector;
        typedef std::vector<double>::iterator DistanceIterator;
	
        //creating an auxiliary list for the new nodes
        PointVector list_of_nodes;

        //*************
        // Bucket types
        typedef Bucket< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
	
        typedef Tree< KDTreePartition<BucketType> > tree; //Kdtree;
	
	
        //starting calculating time of construction of the kdtree
        boost::timer kdtree_construction;
	
        for (ModelPart::NodesContainerType::iterator node_it = rLagrangianModelPart.NodesBegin();
	     node_it != rLagrangianModelPart.NodesEnd(); ++node_it)
	  {
            PointTypePointer pnode = *(node_it.base());
	    
            //putting the nodes of the destination_model part in an auxiliary list
            list_of_nodes.push_back(pnode);
	  }
	
        std::cout << "kdt constructin time " << kdtree_construction.elapsed() << std::endl;
	
        //create a spatial database with the list of new nodes
        unsigned int bucket_size = 20;
        tree nodes_tree(list_of_nodes.begin(), list_of_nodes.end(), bucket_size);
	
        //work arrays
        Node < 3 > work_point(0, 0.0, 0.0, 0.0);
        unsigned int MaximumNumberOfResults = 10000;
        PointVector Results(MaximumNumberOfResults);
        DistanceVector SquaredResultsDistances(MaximumNumberOfResults);
	
	
        if (rEulerianModelPart.NodesBegin()->SolutionStepsDataHas(NODAL_H) == false)
	  KRATOS_THROW_ERROR(std::logic_error, "Add  ----NODAL_H---- variable!!!!!! ERROR", "");
        
        double sigma = 0.0;
        if (TDim == 2)
	  sigma = 10.0 / (7.0 * 3.1415926);
        else
	  sigma = 1.0 / 3.1415926;
	
        for (ModelPart::NodesContainerType::iterator node_it = rEulerianModelPart.NodesBegin(); node_it != rEulerianModelPart.NodesEnd(); node_it++)
	  {
	    if((node_it)->FastGetSolutionStepValue(IS_INTERFACE)==1)
	      { //IS_FREE_SURFACE
		work_point.X() = node_it->X();
		work_point.Y() = node_it->Y();
		work_point.Z() = node_it->Z();
		//KRATOS_THROW_ERROR(std::logic_error, "Add  ----NODAL_H---- variable!!!!!! ERROR", "");
		double radius = 1.5 * node_it->FastGetSolutionStepValue(NODAL_H);
		
            	//find all of the new nodes within the radius
            	int number_of_points_in_radius;
		
            	//look between the new nodes which of them is inside the radius of the circumscribed cyrcle
            	number_of_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, Results.begin(), SquaredResultsDistances.begin(), MaximumNumberOfResults);
		
            	if (number_of_points_in_radius > 0)
		  {
		    
		    //double& temperature = (node_it)->FastGetSolutionStepValue(TEMPERATURE);
		    
		    double temperature_aux = 0.0;
		    
		    double tot_weight = 0.0;
		    
		    for (int k = 0; k < number_of_points_in_radius; k++)
		      {
			
			double distance = sqrt(*(SquaredResultsDistances.begin() + k));
			double weight = SPHCubicKernel(sigma, distance, radius);
			
			PointIterator it_found = Results.begin() + k;
			
			if((*it_found)->FastGetSolutionStepValue(IS_FREE_SURFACE)==1) //MATERIAL_VARIABLE
			  {
			    double tempp=0.0;
			    tempp=(*it_found)->FastGetSolutionStepValue(YCH4);
			    //KRATOS_THROW_ERROR(std::logic_error, "nodo without temperature", "");
			    if(tempp<298.0) tempp=298.0;
			    //else tempp=(*it_found)->FastGetSolutionStepValue(YCH4);
			    temperature_aux += weight * tempp;//temperature
			    tot_weight += weight;
			    //KRATOS_THROW_ERROR(std::logic_error, "Add  ----NODAL_H---- variable!!!!!! ERROR", "");
			  }
		      }
		    if(tot_weight>0.0)
		      {
			temperature_aux /= tot_weight;
			(node_it)->FastGetSolutionStepValue(FUEL)=temperature_aux;
		      }
		    else 
		      {
			KRATOS_WATCH(tot_weight);
			KRATOS_WATCH((node_it)->X());
			KRATOS_WATCH((node_it)->Y());
			if((node_it)->FastGetSolutionStepValue(TEMPERATURE)<298.0) (node_it)->FastGetSolutionStepValue(FUEL)=298.0;
			else (node_it)->FastGetSolutionStepValue(FUEL)=(node_it)->FastGetSolutionStepValue(TEMPERATURE);
		      }
		  }
	      }
            else
	      {
		(node_it)->FastGetSolutionStepValue(FUEL)=(node_it)->FastGetSolutionStepValue(TEMPERATURE);
	      }
	  }
        KRATOS_CATCH("")
	  }
      
      void TransferToEulerianMeshShapeBased_aux_3D(ModelPart& rEulerianModelPart, ModelPart & rLagrangianModelPart, BinBasedFastPointLocator<TDim>& node_locator)
      {
        KRATOS_TRY
	  
	  //defintions for spatial search
	  //typedef Node < 3 > PointType;
	  //typedef Node < 3 > ::Pointer PointTypePointer;
	  
	  array_1d<double, TDim + 1 > N;
        const int max_results = 1000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
        const int nparticles = rLagrangianModelPart.Nodes().size();
	
#pragma omp parallel for firstprivate(results,N)
        for (int i = 0; i < nparticles; i++)
	  {
            ModelPart::NodesContainerType::iterator iparticle = rLagrangianModelPart.NodesBegin() + i;
	    
            Node < 3 > ::Pointer pparticle = *(iparticle.base());
            typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
	    
            Element::Pointer pelement;
	    
            bool is_found = node_locator.FindPointOnMesh(pparticle->Coordinates(), N, pelement, result_begin, max_results);
	    
            if (is_found == true)
	      {
                Geometry<Node<3> >& geom = pelement->GetGeometry();
		
                boost::numeric::ublas::bounded_matrix<double, 4, 3 > msDN_DX;
		
	    	array_1d<double, 4 > msN;
		
	    	double Area=0.0;
	    	GeometryUtils::CalculateGeometryData(geom, msDN_DX, msN, Area);
		
		
		array_1d<double, 3 > qrad=ZeroVector(3);
		
		double temmp=0.0;
		for (unsigned int jj = 0; jj < 3; jj++)
		  {
		    for (unsigned int kk = 0; kk < 4; kk++)
		      {
			temmp=geom[kk].FastGetSolutionStepValue(TEMPERATURE);
			if(temmp<298.0) temmp=298.0;
			
			qrad[jj] += msDN_DX(kk, jj) * temmp;//geom[kk].FastGetSolutionStepValue(TEMPERATURE);
		      }
		  }
		
		
        	//double faceheatflux=0.0;
		
		
		(iparticle)->FastGetSolutionStepValue(NORMAL) *=(-1.0);
		(iparticle)->FastGetSolutionStepValue(FACE_HEAT_FLUX) += fabs( (iparticle)->FastGetSolutionStepValue(NORMAL_X) * qrad[0] + (iparticle)->FastGetSolutionStepValue(NORMAL_Y) * qrad[1] + (iparticle)->FastGetSolutionStepValue(NORMAL_Z) * qrad[2]) *0.0131;
	      }	     
	  }	
        KRATOS_CATCH("")
	  }
      
      
      void TransferToEulerianMeshShapeBased_aux(ModelPart& rEulerianModelPart, ModelPart & rLagrangianModelPart, BinBasedFastPointLocator<TDim>& node_locator)
      {
        KRATOS_TRY
	  
	  array_1d<double, TDim + 1 > N;
        const int max_results = 1000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
        const int nparticles = rLagrangianModelPart.Nodes().size();
	
#pragma omp parallel for firstprivate(results,N)
        for (int i = 0; i < nparticles; i++)
	  {
            ModelPart::NodesContainerType::iterator iparticle = rLagrangianModelPart.NodesBegin() + i;
	    
            Node < 3 > ::Pointer pparticle = *(iparticle.base());
            typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
	    
            Element::Pointer pelement;
	    
            bool is_found = node_locator.FindPointOnMesh(pparticle->Coordinates(), N, pelement, result_begin, max_results);
	    
            if (is_found == true)
	      {
                Geometry<Node<3> >& geom = pelement->GetGeometry();
                
                boost::numeric::ublas::bounded_matrix<double, 3, 2 > msDN_DX;
		
	    	array_1d<double, 3 > msN;
	    	array_1d<double, 3 > N;
		
	    	double Area=0.0;
	    	GeometryUtils::CalculateGeometryData(geom, msDN_DX, msN, Area);
		
		int s0=0;
		int s1=0;
		int s2=0;
		int sum=0;
		if(geom[0].FastGetSolutionStepValue(IS_INTERFACE)>0.5) s0=1; 	//IS_INTERFACE
		
		if(geom[1].FastGetSolutionStepValue(IS_INTERFACE)>0.5) s1=1;
		
		if(geom[2].FastGetSolutionStepValue(IS_INTERFACE)>0.5) s2=1;
		
		sum=s0 + s1 + s2; 
		
		array_1d<double, 2 > qrad=ZeroVector(2);
		array_1d<double, 2 > qrad_P1=ZeroVector(2);
        	array_1d<double,2>  interface_segment=ZeroVector(2);
        	array_1d<double,2>  normaledge1=ZeroVector(2);
		
		for (unsigned int jj = 0; jj < 2; jj++)
		  {
		    for (unsigned int kk = 0; kk < 3; kk++)
		      {
			qrad[jj] += msDN_DX(kk, jj) * geom[kk].FastGetSolutionStepValue(TEMPERATURE);
			qrad_P1[jj] += msDN_DX(kk, jj) * geom[kk].FastGetSolutionStepValue(INCIDENT_RADIATION_FUNCTION);			
		      }
		  }
		
        	double faceheatflux=0.0;
        	//double faceheatflux_P1=0.0;
        	if(sum==2)
		  {
		    if((geom[1].FastGetSolutionStepValue(IS_INTERFACE)>0.5 && geom[0].FastGetSolutionStepValue(IS_INTERFACE)>0.5))  //IS_INTERFACE
		      {
            		double norm=0.0;
	    		//KRATOS_THROW_ERROR(std::logic_error, "element with zero vol found", ""); 
            		interface_segment[0] = (geom[0].X()-geom[1].X());
            		interface_segment[1] = (geom[0].Y()-geom[1].Y());
            		norm = sqrt(  pow((interface_segment[0]),2) + pow((interface_segment[1]),2));
            		//double area1=norm;
            		normaledge1(0)= -interface_segment[1]/norm;
            		normaledge1(1)= interface_segment[0]/norm;
            		faceheatflux += fabs(1.0*(qrad[0]*normaledge1(0)+qrad[1]*normaledge1(1))*0.0131);
		      }
		    
		    if((geom[1].FastGetSolutionStepValue(IS_INTERFACE)>0.5 && geom[2].FastGetSolutionStepValue(IS_INTERFACE)>0.5))
		      {
			double norm=0.0;
			interface_segment[0] = (geom[1].X()-geom[2].X());
			interface_segment[1] = (geom[1].Y()-geom[2].Y());
			norm = sqrt(  pow((interface_segment[0]),2) + pow((interface_segment[1]),2));
			//double area1=norm;
			normaledge1(0)= -interface_segment[1]/norm;
			normaledge1(1)= interface_segment[0]/norm;
			faceheatflux += fabs(1.0*(qrad[0]*normaledge1(0)+qrad[1]*normaledge1(1))*0.0131);
		      }
		    if((geom[2].FastGetSolutionStepValue(IS_INTERFACE)>0.5 && geom[0].FastGetSolutionStepValue(IS_INTERFACE)>0.5))
		      {
            		double norm=0.0;
			interface_segment[0] = (geom[2].X()-geom[0].X());
			interface_segment[1] = (geom[2].Y()-geom[0].Y());
			norm = sqrt(  pow((interface_segment[0]),2) + pow((interface_segment[1]),2));
			normaledge1(0)= -interface_segment[1]/norm;
			normaledge1(1)= interface_segment[0]/norm;
			faceheatflux += fabs(1.0*(qrad[0]*normaledge1(0)+qrad[1]*normaledge1(1))*0.0131);
		      }
		  }
		if(sum==1)
		  {
		    if((geom[1].FastGetSolutionStepValue(IS_INTERFACE)<0.5 && geom[0].FastGetSolutionStepValue(IS_INTERFACE)<0.5))  //IS_INTERFACE
		      {
            		double norm=0.0;
			interface_segment[0] = (geom[0].X()-geom[1].X());
            		interface_segment[1] = (geom[0].Y()-geom[1].Y());
            		norm = sqrt(  pow((interface_segment[0]),2) + pow((interface_segment[1]),2));
			normaledge1(0)= -interface_segment[1]/norm;
            		normaledge1(1)= interface_segment[0]/norm;
			faceheatflux += fabs(1.0*(qrad[0]*normaledge1(0)+qrad[1]*normaledge1(1))*0.0131);
		      }		    
		    if((geom[1].FastGetSolutionStepValue(IS_INTERFACE)<0.5 && geom[2].FastGetSolutionStepValue(IS_INTERFACE)<0.5))
		      {
			double norm=0.0;
			interface_segment[0] = (geom[1].X()-geom[2].X());
			interface_segment[1] = (geom[1].Y()-geom[2].Y());
			norm = sqrt(  pow((interface_segment[0]),2) + pow((interface_segment[1]),2));
			normaledge1(0)= -interface_segment[1]/norm;
			normaledge1(1)= interface_segment[0]/norm;
			faceheatflux += fabs(1.0*(qrad[0]*normaledge1(0)+qrad[1]*normaledge1(1))*0.0131);
		      }		    
		    if((geom[2].FastGetSolutionStepValue(IS_INTERFACE)<0.5 && geom[0].FastGetSolutionStepValue(IS_INTERFACE)<0.5))
		      {
            		double norm=0.0;
			interface_segment[0] = (geom[2].X()-geom[0].X());
			interface_segment[1] = (geom[2].Y()-geom[0].Y());
			norm = sqrt(  pow((interface_segment[0]),2) + pow((interface_segment[1]),2));
			normaledge1(0)= -interface_segment[1]/norm;
			normaledge1(1)= interface_segment[0]/norm;
			faceheatflux += fabs(1.0*(qrad[0]*normaledge1(0)+qrad[1]*normaledge1(1))*0.0131);
		      }
		  }
		
		(iparticle)->FastGetSolutionStepValue(FACE_HEAT_FLUX)+=(faceheatflux /*+ fhf+ faceheatflux_P1*/);
	      }
	  }
	KRATOS_CATCH("")
	  }
      ///3D
      void CalculateNormal(ModelPart& full_model_part){
	KRATOS_TRY
	  //resetting the normals
	  array_1d<double,3> zero;
	noalias(zero) = ZeroVector(3);
	
	for(ModelPart::NodesContainerType::const_iterator in = full_model_part.NodesBegin(); in!=full_model_part.NodesEnd(); in++)
	  {
	    in->FastGetSolutionStepValue(NORMAL) = zero;
	  }
	
	array_1d<double,3> v1;
	array_1d<double,3> v2;
	//array_1d<double,3>& An =zero;
	//double area_normal=0.0;
	array_1d<double,3> area_normal;
	
	for(ModelPart::ElementsContainerType::iterator iii = full_model_part.ElementsBegin(); iii != full_model_part.ElementsEnd(); iii++)
	  {
	    if( iii->GetGeometry()[1].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1.0 && iii->GetGeometry()[2].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1.0 && iii->GetGeometry()[3].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1.0)
	      {		
		v1[0] =  iii->GetGeometry()[1].X() -iii->GetGeometry()[3].X();
		v1[1] = iii->GetGeometry()[1].Y() - iii->GetGeometry()[3].Y();
		v1[2] = iii->GetGeometry()[1].Z() - iii->GetGeometry()[3].Z();
		
		v2[0] = iii->GetGeometry()[2].X() - iii->GetGeometry()[3].X();
		v2[1] = iii->GetGeometry()[2].Y() - iii->GetGeometry()[3].Y();
		v2[2] = iii->GetGeometry()[2].Z() - iii->GetGeometry()[3].Z();
		
		MathUtils<double>::CrossProduct(area_normal,v1,v2);
		//area_normal *= -0.5;
		
		array_1d<double,3> msAuxVec = ZeroVector(3);
		double c0 = fabs(area_normal[0]);
		double c1 = fabs(area_normal[1]);
		double c2 = fabs(area_normal[2]);
		msAuxVec[0]=c0;
		msAuxVec[1]=c1;
		msAuxVec[2]=c2;
		//					double norm_c =norm_2(msAuxVec);
		
		double norm_u = msAuxVec[0]*msAuxVec[0] + msAuxVec[1]*msAuxVec[1] + msAuxVec[2]*msAuxVec[2];
		double norm_c =sqrt(norm_u);
		
		iii->GetGeometry()[1].FastGetSolutionStepValue(NORMAL) += area_normal/ norm_c;
		iii->GetGeometry()[2].FastGetSolutionStepValue(NORMAL) += area_normal/ norm_c;
		iii->GetGeometry()[3].FastGetSolutionStepValue(NORMAL) += area_normal/ norm_c;
		
	      }
	    if( iii->GetGeometry()[0].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1.0 && iii->GetGeometry()[3].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1.0 && iii->GetGeometry()[2].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1.0)
	      {
		v1[0] =  iii->GetGeometry()[0].X() -iii->GetGeometry()[2].X();
		v1[1] = iii->GetGeometry()[0].Y() - iii->GetGeometry()[2].Y();
		v1[2] = iii->GetGeometry()[0].Z() - iii->GetGeometry()[2].Z();
		
		v2[0] = iii->GetGeometry()[3].X() - iii->GetGeometry()[2].X();
		v2[1] = iii->GetGeometry()[3].Y() - iii->GetGeometry()[2].Y();
		v2[2] = iii->GetGeometry()[3].Z() - iii->GetGeometry()[2].Z();
		MathUtils<double>::CrossProduct(area_normal,v1,v2);
		//area_normal *= -0.5;
		array_1d<double,3> msAuxVec = ZeroVector(3);
		double c0 = fabs(area_normal[0]);
		double c1 = fabs(area_normal[1]);
		double c2 = fabs(area_normal[2]);
		msAuxVec[0]=c0;
		msAuxVec[1]=c1;
		msAuxVec[2]=c2;
		//double norm_c =norm_2(msAuxVec);
		
		double norm_u = msAuxVec[0]*msAuxVec[0] + msAuxVec[1]*msAuxVec[1] + msAuxVec[2]*msAuxVec[2];
		double norm_c =sqrt(norm_u);
		
		iii->GetGeometry()[0].FastGetSolutionStepValue(NORMAL) += area_normal/ norm_c;
		iii->GetGeometry()[3].FastGetSolutionStepValue(NORMAL) += area_normal/ norm_c;
		iii->GetGeometry()[2].FastGetSolutionStepValue(NORMAL) += area_normal/ norm_c;
		
	      }
	    
	    if( iii->GetGeometry()[0].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1.0 && iii->GetGeometry()[1].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1.0 && iii->GetGeometry()[3].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1.0)
	      {
		
		v1[0] =  iii->GetGeometry()[0].X() -iii->GetGeometry()[3].X();
		v1[1] = iii->GetGeometry()[0].Y() - iii->GetGeometry()[3].Y();
		v1[2] = iii->GetGeometry()[0].Z() - iii->GetGeometry()[3].Z();
		
		v2[0] = iii->GetGeometry()[1].X() - iii->GetGeometry()[3].X();
		v2[1] = iii->GetGeometry()[1].Y() - iii->GetGeometry()[3].Y();
		v2[2] = iii->GetGeometry()[1].Z() - iii->GetGeometry()[3].Z();
		
		MathUtils<double>::CrossProduct(area_normal,v1,v2);
		//area_normal *= -0.5;
		array_1d<double,3> msAuxVec = ZeroVector(3);
		double c0 = fabs(area_normal[0]);
		double c1 = fabs(area_normal[1]);
		double c2 = fabs(area_normal[2]);
		msAuxVec[0]=c0;
		msAuxVec[1]=c1;
		msAuxVec[2]=c2;
		
		double norm_u = msAuxVec[0]*msAuxVec[0] + msAuxVec[1]*msAuxVec[1] + msAuxVec[2]*msAuxVec[2];
		double norm_c =sqrt(norm_u);
		
		iii->GetGeometry()[0].FastGetSolutionStepValue(NORMAL) += area_normal/ norm_c;
		iii->GetGeometry()[1].FastGetSolutionStepValue(NORMAL) += area_normal/ norm_c;
		iii->GetGeometry()[3].FastGetSolutionStepValue(NORMAL) += area_normal/ norm_c;				
		
	      }
	    if( iii->GetGeometry()[0].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1.0 && iii->GetGeometry()[2].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1.0 && iii->GetGeometry()[1].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1.0){
	      
	      v1[0] =  iii->GetGeometry()[0].X() -iii->GetGeometry()[1].X();
	      v1[1] = iii->GetGeometry()[0].Y() - iii->GetGeometry()[1].Y();
	      v1[2] = iii->GetGeometry()[0].Z() - iii->GetGeometry()[1].Z();
	      
	      v2[0] = iii->GetGeometry()[2].X() - iii->GetGeometry()[1].X();
	      v2[1] = iii->GetGeometry()[2].Y() - iii->GetGeometry()[1].Y();
	      v2[2] = iii->GetGeometry()[2].Z() - iii->GetGeometry()[1].Z();
	      
	      MathUtils<double>::CrossProduct(area_normal,v1,v2);
	      //area_normal *= -0.5;	
	      array_1d<double,3> msAuxVec = ZeroVector(3);
	      double c0 = fabs(area_normal[0]);
	      double c1 = fabs(area_normal[1]);
	      double c2 = fabs(area_normal[2]);
	      msAuxVec[0]=c0;
	      msAuxVec[1]=c1;
	      msAuxVec[2]=c2;
	      //					double norm_c =norm_2(msAuxVec);
	      
	      double norm_u = msAuxVec[0]*msAuxVec[0] + msAuxVec[1]*msAuxVec[1] + msAuxVec[2]*msAuxVec[2];
	      double norm_c =sqrt(norm_u);
	      
	      iii->GetGeometry()[0].FastGetSolutionStepValue(NORMAL) += area_normal/ norm_c;
	      iii->GetGeometry()[2].FastGetSolutionStepValue(NORMAL) += area_normal/ norm_c;
	      iii->GetGeometry()[1].FastGetSolutionStepValue(NORMAL) += area_normal/ norm_c;	
	    }
	    
	  }
	for(ModelPart::NodesContainerType::iterator iii = full_model_part.NodesBegin(); iii != full_model_part.NodesEnd(); iii++)
	  {
	    if(iii->FastGetSolutionStepValue(IS_FREE_SURFACE)==1.0){
	      array_1d<double,3>& value_y1 = iii->FastGetSolutionStepValue(NORMAL);
	      double norm_y1 =norm_2(value_y1);
	      value_y1 /=norm_y1;
	    }	
	  }
	
	KRATOS_CATCH("")	
	  }
      
      //restarting the step from the beginning
      void RestartStep(ModelPart & rModelPart)
      {
        KRATOS_TRY;
	
        //setting the variables to their value at the beginning of the time step
        rModelPart.OverwriteSolutionStepData(1, 0);
	
        //setting the coordinates to their value at the beginning of the step
        for (ModelPart::NodesContainerType::iterator node_it = rModelPart.NodesBegin();node_it != rModelPart.NodesEnd(); node_it++)
	  {
            array_1d<double, 3 > & coords = node_it->Coordinates();
            const array_1d<double, 3 > & old_disp = node_it->FastGetSolutionStepValue(DISPLACEMENT, 1);
	    
            coords[0] = node_it->X0() + old_disp[0];
            coords[1] = node_it->Y0() + old_disp[1];
            coords[2] = node_it->Z0() + old_disp[2];
	  }
	
	
        KRATOS_CATCH("");
      }
      
      
      void MoveMesh_Streamlines_freesurfaceflows(ModelPart& rModelPart, unsigned int substeps)
      {      
	const double dt = rModelPart.GetProcessInfo()[DELTA_TIME];
	//KRATOS_THROW_ERROR(std::logic_error, "element with zero vol found", "");
	BinBasedFastPointLocator<TDim> SearchStructure(rModelPart);
	SearchStructure.UpdateSearchDatabase();
	
	//do movement
	array_1d<double, 3 > veulerian;
	double temperature=0.0;
	array_1d<double, 3 > acc_particle;
	array_1d<double, TDim + 1 > N;
	const int max_results = 10000;
	typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
	
	const int nparticles = rModelPart.Nodes().size();
	
#pragma omp parallel for firstprivate(results,N,veulerian,acc_particle)
        for (int i = 0; i < nparticles; i++)
	  {
            //int substep = 0;
            int subdivisions = 5;
	    
	    ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;
	    Node < 3 > ::Pointer pparticle = *(iparticle.base());
	    //small_dt = dt / subdivisions;
	    
	    bool do_move = true;
	    bool first_time=false;
	    iparticle->FastGetSolutionStepValue(DISTANCE)=0.0;
	    
	    iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY) = iparticle->FastGetSolutionStepValue(VELOCITY,1);  //AUX_VEL
	    
	    if(iparticle->FastGetSolutionStepValue(IS_STRUCTURE) == 1.0) do_move = false;	
	    
	    //iparticle->FastGetSolutionStepValue(TEMPERATURE) = 298.0;
	    if( do_move == true  ) //note that we suppose the velocity components to be all fixed
	      {
		array_1d<double,3> old_position = pparticle->Coordinates();
		array_1d<double,3> current_position = pparticle->Coordinates();
		noalias(iparticle->GetInitialPosition()) = old_position;
		iparticle->FastGetSolutionStepValue(DISPLACEMENT,1) = ZeroVector(3);
		//array_1d<double, 3 > & vel_particle = iparticle->FastGetSolutionStepValue(VELOCITY); 
		//subdivisions=10;
		const double small_dt = dt / subdivisions;
		//
		for (int substep = 0; substep < subdivisions; substep++)
		  {
		    typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
		    Element::Pointer pelement;
		    bool is_found = SearchStructure.FindPointOnMesh(current_position, N, pelement, result_begin, max_results);
		    iparticle->Set(TO_ERASE, true);
		    //(iparticle)->GetValue(ERASE_FLAG) = true;
		    //KRATOS_WATCH(is_found);
		    if (is_found == true)
		      {
			Geometry< Node < 3 > >& geom = pelement->GetGeometry();
			//int nn=0;
			noalias(veulerian) = ZeroVector(3); //0.0;//N[0] * geom[0].FastGetSolutionStepValue(VELOCITY,1);  
			temperature=0.0;//N[0] * geom[0].FastGetSolutionStepValue(TEMPERATURE);  
			for (unsigned int k = 0; k < geom.size(); k++)
			  {
			    noalias(veulerian) += N[k] * geom[k].FastGetSolutionStepValue(VELOCITY,1);
			  }
			first_time=true;
			noalias(current_position) += small_dt*veulerian;
			pparticle->Set(TO_ERASE, false);
			iparticle->FastGetSolutionStepValue(DISTANCE) += small_dt;
			iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY)=veulerian;
		      }
		    else
		      {
			double time1=iparticle->FastGetSolutionStepValue(DISTANCE);
			array_1d<double,3> acc;
			acc[0] =  0.0;
			acc[1] = -10.0;
			acc[2] =  0.0;
			if( first_time == false && iparticle->FastGetSolutionStepValue(IS_STRUCTURE) == 0.0 ) 
			  {
			    noalias(current_position) += small_dt *iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY);
			    //noalias(current_position) += small_dt * small_dt * acc;	
			    pparticle->Set(TO_ERASE, false);
			  }	
			else
			  {
			    time1 -=small_dt;
			    //double tiempo_restante=dt-time1;
			    noalias(current_position) += small_dt *iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY);
			    //noalias(current_position) += small_dt * small_dt * acc;	
			    pparticle->Set(TO_ERASE, false);
			  }
			
		       }
		  }//for
		
		//update the displacement BUT DO NOT OVERWRITE THE POSITION!!                
		iparticle->FastGetSolutionStepValue(DISPLACEMENT) = current_position - iparticle->GetInitialPosition(); 
		//KRATOS_WATCH(iparticle->FastGetSolutionStepValue(DISPLACEMENT));
	      }//move
	  }
	//compute mesh velocity
        for(ModelPart::NodesContainerType::iterator it = rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); it++)
	  {
	    //array_1d<double,3>& dn = it->FastGetSolutionStepValue(DISPLACEMENT,1);
	    array_1d<double,3>& dn1 = it->FastGetSolutionStepValue(DISPLACEMENT);
	    noalias(it->Coordinates()) = it->GetInitialPosition();
	    noalias(it->Coordinates()) += dn1; 
	  }
      }
      
      void MoveLonelyNodes(ModelPart& ThisModelPart)
      {
        KRATOS_TRY;
        double Dt = ThisModelPart.GetProcessInfo()[DELTA_TIME];
	
        array_1d<double,3> DeltaDisp, acc;
	
        for(ModelPart::NodeIterator i = ThisModelPart.NodesBegin() ;
	    i != ThisModelPart.NodesEnd() ; ++i)
	  {
            if(
	       (i)->FastGetSolutionStepValue(IS_STRUCTURE) == 0 && 
	       (i)->GetValue(NEIGHBOUR_ELEMENTS).size() == 0 &&
	       ((i)->GetDof(VELOCITY_X).IsFixed() == false || (i)->GetDof(VELOCITY_Y).IsFixed() == false || (i)->GetDof(VELOCITY_Z).IsFixed() == false) 
	       )
	      {
                //i->Set(TO_ERASE,true);
                //set to zero the pressure
                (i)->FastGetSolutionStepValue(PRESSURE) = 0;
		const array_1d<double,3>& old_vel = (i)->FastGetSolutionStepValue(VELOCITY,1);
                array_1d<double,3>& vel = (i)->FastGetSolutionStepValue(VELOCITY);
                //array_1d<double,3>& acc = (i)->FastGetSolutionStepValue(ACCELERATION);
		
                noalias(acc) =  (i)->FastGetSolutionStepValue(BODY_FORCE);
                acc[0]= 0.0;
	        acc[1]= -10.0;
                acc[2]= 0.0; 
                noalias(vel) = old_vel;
                noalias(vel) += Dt * acc ;
		
                //calculate displacements
                //noalias(DeltaDisp) = Dt * vel;
		
                //array_1d<double,3>& disp = i->FastGetSolutionStepValue(DISPLACEMENT);
                //noalias(disp) = i->FastGetSolutionStepValue(DISPLACEMENT,1);
                //noalias(disp) += DeltaDisp;
		noalias(i->Coordinates()) += Dt * Dt * acc;
		
	      }
	    
	  }
	
        KRATOS_CATCH("")
	  }
      
      

    private:
      
      inline double SPHCubicKernel(const double sigma, const double r, const double hmax)
      {
        double h_half = 0.5 * hmax;
        const double s = r / h_half;
        const double coeff = sigma / pow(h_half, static_cast<int>(TDim));
	
        if (s <= 1.0)
	  return coeff * (1.0 - 1.5 * s * s + 0.75 * s * s * s);
        else if (s <= 2.0)
	  return 0.25 * coeff * pow(2.0 - s, 3);
        else
            return 0.0;
      }
      
      inline void CalculateCenterAndSearchRadius(Geometry<Node < 3 > >&geom, double& xc, double& yc, double& zc, double& R, array_1d<double, 3 > & N )
      {
        double x0 = geom[0].X();
        double y0 = geom[0].Y();
        double x1 = geom[1].X();
        double y1 = geom[1].Y();
        double x2 = geom[2].X();
        double y2 = geom[2].Y();
	
	xc = 0.3333333333333333333 * (x0 + x1 + x2);
        yc = 0.3333333333333333333 * (y0 + y1 + y2);
        zc = 0.0;
	
        double R1 = (xc - x0)*(xc - x0) + (yc - y0)*(yc - y0);
        double R2 = (xc - x1)*(xc - x1) + (yc - y1)*(yc - y1);
        double R3 = (xc - x2)*(xc - x2) + (yc - y2)*(yc - y2);
	
        R = R1;
        if (R2 > R) R = R2;
        if (R3 > R) R = R3;
	
        R = 1.01 * sqrt(R);
      }
      
      inline void CalculateCenterAndSearchRadius(Geometry<Node < 3 > >&geom, double& xc, double& yc, double& zc, double& R, array_1d<double, 4 > & N )
      {
        double x0 = geom[0].X();
        double y0 = geom[0].Y();
        double z0 = geom[0].Z();
        double x1 = geom[1].X();
        double y1 = geom[1].Y();
        double z1 = geom[1].Z();
        double x2 = geom[2].X();
        double y2 = geom[2].Y();
        double z2 = geom[2].Z();
        double x3 = geom[3].X();
        double y3 = geom[3].Y();
        double z3 = geom[3].Z();

	
        xc = 0.25 * (x0 + x1 + x2 + x3);
        yc = 0.25 * (y0 + y1 + y2 + y3);
        zc = 0.25 * (z0 + z1 + z2 + z3);
	
        double R1 = (xc - x0)*(xc - x0) + (yc - y0)*(yc - y0) + (zc - z0)*(zc - z0);
        double R2 = (xc - x1)*(xc - x1) + (yc - y1)*(yc - y1) + (zc - z1)*(zc - z1);
        double R3 = (xc - x2)*(xc - x2) + (yc - y2)*(yc - y2) + (zc - z2)*(zc - z2);
        double R4 = (xc - x3)*(xc - x3) + (yc - y3)*(yc - y3) + (zc - z3)*(zc - z3);
	
        R = R1;
        if (R2 > R) R = R2;
        if (R3 > R) R = R3;
        if (R4 > R) R = R4;

        R = sqrt(R);
      }

      
      inline bool CalculatePosition(Geometry<Node < 3 > >&geom,const double xc, const double yc, const double zc,  array_1d<double, 4 > & N )
      {
    
	double x0 = geom[0].X();
	double y0 = geom[0].Y();
	double z0 = geom[0].Z();
	double x1 = geom[1].X();
	double y1 = geom[1].Y();
	double z1 = geom[1].Z();
	double x2 = geom[2].X();
	double y2 = geom[2].Y();
	double z2 = geom[2].Z();
	double x3 = geom[3].X();
	double y3 = geom[3].Y();
	double z3 = geom[3].Z();
    
	double vol = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);
    
	double inv_vol = 0.0;
	if (vol < 0.0000000000001)
	  {
	    KRATOS_THROW_ERROR(std::logic_error, "element with zero vol found", "");
	  }
	else
	  {
	    inv_vol = 1.0 / vol;
	  }
	
	N[0] = CalculateVol(x1, y1, z1, x3, y3, z3, x2, y2, z2, xc, yc, zc) * inv_vol;
	N[1] = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, xc, yc, zc) * inv_vol;
	N[2] = CalculateVol(x3, y3, z3, x1, y1, z1, x0, y0, z0, xc, yc, zc) * inv_vol;
	N[3] = CalculateVol(x3, y3, z3, x0, y0, z0, x2, y2, z2, xc, yc, zc) * inv_vol;
    
    
	if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[3] >= 0.0 &&
	    N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0 && N[3] <= 1.0)
	  //if the xc yc zc is inside the tetrahedron return true
	  return true;
    
	return false;
      }
  
      inline double CalculateVol(const double x0, const double y0,
				 const double x1, const double y1,
				 const double x2, const double y2
				 )
      {
	return 0.5 * ((x1 - x0)*(y2 - y0)- (y1 - y0)*(x2 - x0));
      }
  
  
      inline double CalculateVol(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1, const double x2, const double y2, const double z2,const double x3, const double y3, const double z3 )
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
	return detJ * 0.1666666666666666666667;
      }
  
      void ComputeGaussPointPositions(Geometry< Node < 3 > >& geom, boost::numeric::ublas::bounded_matrix<double, 4, 3 > & pos, boost::numeric::ublas::bounded_matrix<double, 4, 3 > & N)
      {
	double one_third = 1.0 / 3.0;
	double one_sixt = 1.0 / 6.0;
	double two_third = 2.0 * one_third;
      
	N(0, 0) = one_sixt;
	N(0, 1) = one_sixt;
	N(0, 2) = two_third;
	N(1, 0) = two_third;
	N(1, 1) = one_sixt;
	N(1, 2) = one_sixt;
	N(2, 0) = one_sixt;
	N(2, 1) = two_third;
	N(2, 2) = one_sixt;
	N(3, 0) = one_third;
	N(3, 1) = one_third;
	N(3, 2) = one_third;
      
      
	//first
	pos(0, 0) = one_sixt * geom[0].X() + one_sixt * geom[1].X() + two_third * geom[2].X();
	pos(0, 1) = one_sixt * geom[0].Y() + one_sixt * geom[1].Y() + two_third * geom[2].Y();
	pos(0, 2) = one_sixt * geom[0].Z() + one_sixt * geom[1].Z() + two_third * geom[2].Z();
      
	//second
	pos(1, 0) = two_third * geom[0].X() + one_sixt * geom[1].X() + one_sixt * geom[2].X();
	pos(1, 1) = two_third * geom[0].Y() + one_sixt * geom[1].Y() + one_sixt * geom[2].Y();
	pos(1, 2) = two_third * geom[0].Z() + one_sixt * geom[1].Z() + one_sixt * geom[2].Z();
      
	//third
	pos(2, 0) = one_sixt * geom[0].X() + two_third * geom[1].X() + one_sixt * geom[2].X();
	pos(2, 1) = one_sixt * geom[0].Y() + two_third * geom[1].Y() + one_sixt * geom[2].Y();
	pos(2, 2) = one_sixt * geom[0].Z() + two_third * geom[1].Z() + one_sixt * geom[2].Z();
      
	//fourth
	pos(3, 0) = one_third * geom[0].X() + one_third * geom[1].X() + one_third * geom[2].X();
	pos(3, 1) = one_third * geom[0].Y() + one_third * geom[1].Y() + one_third * geom[2].Y();
	pos(3, 2) = one_third * geom[0].Z() + one_third * geom[1].Z() + one_third * geom[2].Z();

      }
  
      void ComputeGaussPointPositions(Geometry< Node < 3 > >& geom, boost::numeric::ublas::bounded_matrix<double, 16, 3 > & pos, boost::numeric::ublas::bounded_matrix<double, 16, 3 > & N)
      {
	//lower diagonal terms
	double ypos = 1.0 / 12.0;
	int pos_counter = 0;
	for (unsigned int i = 0; i < 4; i++)
	  {
	    double xpos = 1.0 / 12.0;
	    for (unsigned int j = 0; j < 4 - i; j++)
	      {
		double N1 = xpos;
		double N2 = ypos;
		double N3 = 1.0 - xpos - ypos;

		pos(pos_counter, 0) = N1 * geom[0].X() + N2 * geom[1].X() + N3 * geom[2].X();
		pos(pos_counter, 1) = N1 * geom[0].Y() + N2 * geom[1].Y() + N3 * geom[2].Y();
		pos(pos_counter, 2) = N1 * geom[0].Z() + N2 * geom[1].Z() + N3 * geom[2].Z();
	    
		N(pos_counter, 0) = N1;
		N(pos_counter, 1) = N2;
		N(pos_counter, 2) = N3;
	    
		xpos += 1.0 / 4.0;
		pos_counter += 1;
	    
	      }
	    ypos += 1.0 / 4.0;
	  }
    
	//lower diagonal terms
	ypos = 2.0 / 12.0;
	// pos_counter = 8;
	for (unsigned int i = 0; i < 3; i++)
	  {
	    double xpos = 2.0 / 12.0;
	    for (unsigned int j = 0; j < 4 - i; j++)
	      {
		double N1 = xpos;
		double N2 = ypos;
		double N3 = 1.0 - xpos - ypos;
	    
		pos(pos_counter, 0) = N1 * geom[0].X() + N2 * geom[1].X() + N3 * geom[2].X();
		pos(pos_counter, 1) = N1 * geom[0].Y() + N2 * geom[1].Y() + N3 * geom[2].Y();
		pos(pos_counter, 2) = N1 * geom[0].Z() + N2 * geom[1].Z() + N3 * geom[2].Z();
	    
		N(pos_counter, 0) = N1;
		N(pos_counter, 1) = N2;
		N(pos_counter, 2) = N3;
	    
		xpos += 1.0 / 4.0;
		pos_counter += 1;
	    
	      }
	    ypos += 1.0 / 4.0;
	  }
      }
  
      void ConsistentMassMatrix(const double A, boost::numeric::ublas::bounded_matrix<double, 3, 3 > & M)
      {
	double c1 = A / 12.0;
	double c2 = 2.0 * c1;
	M(0, 0) = c2;
	M(0, 1) = c1;
	M(0, 2) = c1;
	M(1, 0) = c1;
	M(1, 1) = c2;
	M(1, 2) = c1;
	M(2, 0) = c1;
	M(2, 1) = c1;
	M(2, 2) = c2;
      }
      //<-------><-------><-------><-------><-------><-------><-------><-------><-------><-------><------->
      //<-------><-------><-------><-------><-------><-------><-------><-------><-------><-------><------->
      void CalculateInterfaceNormal(boost::numeric::ublas::bounded_matrix<double, 3, 2 >& rPoints, array_1d<double,3>&  rDistances, array_1d<double,2>&  normal, double & interface_area, array_1d<double,3>&  Ninterface, boost::numeric::ublas::bounded_matrix<double, 2, 2 >& rInterfacePoints)
      {
	double sign_correction=1.0;
      
      
      
	boost::numeric::ublas::bounded_matrix<double, 2, 2 > InterfacePoints;
      
	array_1d<bool,3>  cut_edges;
      
	array_1d<double,2>  interface_segment=ZeroVector(2);
      
	if ((rDistances(0)*rDistances(1))<0.0) cut_edges[0]=true;//edge 12 is cut        
      
	else         cut_edges[0]=false;
      
      
      
	if ((rDistances(1)*rDistances(2))<0.0) cut_edges[1]=true;//edge 23 is cut. 
      
	else         cut_edges[1]=false;
      
      
      
	if ((rDistances(2)*rDistances(0))<0.0) cut_edges[2]=true;//edge 13 is cut. 
      
	else         cut_edges[2]=false;
      
      
      
      
      
	if (cut_edges[0])
	
	  {
	  
	    if (rDistances(0)>0.0) sign_correction=1.0;
	  
	    else sign_correction=-1.0;
	  
	  
      
	    const double relative_position = fabs(rDistances(1)/(rDistances(1)-rDistances(0) ) ); 
	  
	    InterfacePoints(0,0) = relative_position*rPoints(0,0) +  (1.0-relative_position)*rPoints(1,0);
	  
	    InterfacePoints(0,1) = relative_position*rPoints(0,1) +  (1.0-relative_position)*rPoints(1,1);
	  
	  
	  
	    if (cut_edges[1])
	    
	      {
	      
		const double relative_position2 = fabs(rDistances(2)/(rDistances(1)-rDistances(2) ) ); 
	      
		InterfacePoints(1,0) = relative_position2*rPoints(1,0) +  (1.0-relative_position2)*rPoints(2,0);
	      
		InterfacePoints(1,1) = relative_position2*rPoints(1,1) +  (1.0-relative_position2)*rPoints(2,1);
	      
	      }
      
	    else
	    
	      {
	      
		const double relative_position2 = fabs(rDistances(0)/(rDistances(2)-rDistances(0) ) ); 
	      
		InterfacePoints(1,0) = relative_position2*rPoints(2,0) +  (1.0-relative_position2)*rPoints(0,0);
	      
		InterfacePoints(1,1) = relative_position2*rPoints(2,1) +  (1.0-relative_position2)*rPoints(0,1);
	      
	      }
	  
	  }
      
	else
	
	  {
	  
	    if (rDistances(1)>0.0) sign_correction=1.0;
	  
	    else sign_correction=-1.0;
	  
	  
	  
	    const double relative_position = fabs(rDistances(2)/(rDistances(2)-rDistances(1) ) ); 
	  
	    InterfacePoints(0,0) = relative_position*rPoints(1,0) +  (1.0-relative_position)*rPoints(2,0);
	  
	    InterfacePoints(0,1) = relative_position*rPoints(1,1) +  (1.0-relative_position)*rPoints(2,1);
	  
	  
	  
	    const double relative_position2 = fabs(rDistances(0)/(rDistances(2)-rDistances(0) ) ); 
	  
	    InterfacePoints(1,0) = relative_position2*rPoints(2,0) +  (1.0-relative_position2)*rPoints(0,0);
	  
	    InterfacePoints(1,1) = relative_position2*rPoints(2,1) +  (1.0-relative_position2)*rPoints(0,1);
	  
	  }        
      
	interface_segment[0] = (InterfacePoints(1,0)-InterfacePoints(0,0));
      
	interface_segment[1] = (InterfacePoints(1,1)-InterfacePoints(0,1));        
      
      
	const double norm = sqrt(  pow((interface_segment[0]),2) + pow((interface_segment[1]),2));
      
  
      
	normal(0)= -interface_segment[1]*sign_correction/norm;
      
	normal(1)= interface_segment[0]*sign_correction/norm;
      
	//KRATOS_WATCH(interface_segment)
      
	//KRATOS_WATCH(InterfacePoints)
      
	interface_area=norm;
      
	rInterfacePoints(0,0)=InterfacePoints(0,0);
	rInterfacePoints(0,1)=InterfacePoints(0,1);
	rInterfacePoints(1,0)=InterfacePoints(1,0);
	rInterfacePoints(1,1)=InterfacePoints(1,1);
      
	const double x_interface = 0.5*(InterfacePoints(0,0)+InterfacePoints(1,0));
      
	const double y_interface = 0.5*(InterfacePoints(0,1)+InterfacePoints(1,1));
      
	// CalculatePosition(rPoints, x_interface, y_interface, 0.0, Ninterface );
      
	///meto aqui el CalculatePosition(rPoints, x_interface, y_interface, 0.0, Ninterface );
	double x0 = rPoints(0,0);
	double y0 = rPoints(0,1);
	double x1 = rPoints(1,0);
	double y1 = rPoints(1,1);
	double x2 = rPoints(2,0);
	double y2 = rPoints(2,1);
	double area = CalculateVol(x0, y0, x1, y1, x2, y2);
	double inv_area = 0.0;
	if (area == 0.0)
	  {
	    KRATOS_THROW_ERROR(std::logic_error, "element with zero area found", "");
	  }
	else
	  {
	    inv_area = 1.0 / area;
	  }
      
	Ninterface[0]= CalculateVol(x1, y1, x2, y2, x_interface, y_interface) * inv_area;
	Ninterface[1] = CalculateVol(x2, y2, x0, y0, x_interface, y_interface) * inv_area;
	Ninterface[2] = CalculateVol(x0, y0, x1, y1, x_interface, y_interface) * inv_area;
      
      }
    
      bool CalculatePosition(Geometry<Node < 3 > >&geom,const double xc, const double yc, const double zc,array_1d<double, 3 > & N	)
      {
	double x0 = geom[0].X();
	double y0 = geom[0].Y();
	double x1 = geom[1].X();
	double y1 = geom[1].Y();
	double x2 = geom[2].X();
	double y2 = geom[2].Y();
      
	double area = CalculateVol(x0, y0, x1, y1, x2, y2);
	double inv_area = 0.0;
	if (area == 0.0)
	  {
	    KRATOS_THROW_ERROR(std::logic_error, "element with zero area found", "");
	  }
	else
	  {
	    inv_area = 1.0 / area;
	  }
    
    
	N[0] = CalculateVol(x1, y1, x2, y2, xc, yc) * inv_area;
	N[1] = CalculateVol(x2, y2, x0, y0, xc, yc) * inv_area;
	N[2] = CalculateVol(x0, y0, x1, y1, xc, yc) * inv_area;
    
    
	if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0) //if the xc yc is inside the triangle return true
	  return true;
    
	return false;
      }

      template<class T>
	bool InvertMatrix(const T& input, T& inverse)

	{
	  typedef permutation_matrix<std::size_t> pmatrix;
  
	  // create a working copy of the input
  
	  T A(input);
  
	  // create a permutation matrix for the LU-factorization
  
	  pmatrix pm(A.size1());
  
	  // perform LU-factorization
  
	  int res = lu_factorize(A, pm);
  
	  if (res != 0)
    
	    return false;
  
	  // create identity matrix of "inverse"
  
	  inverse.assign(identity_matrix<double> (A.size1()));
  
	  // backsubstitute to get the inverse
  
	  lu_substitute(A, pm, inverse);
	  
	  return true;
    

	}

    };

} // namespace Kratos.

#endif // KRATOS_LAGRANGIAN_PARTICLES_UTILITIES_INCLUDED  defined


