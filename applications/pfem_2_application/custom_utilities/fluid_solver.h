/*
* File:   edgebased_levelset.h
* Author: rrossi
*
* Created on July 31, 2009, 10:51 AM
*/
 
/*
==============================================================================
KratosPFEMApplication
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
//   Last Modified by:    $Author: antonia $
//   Date:                $Date: 2009-01-14 16:24:38 $
//   Revision:            $Revision: 1.11 $
//
//


#if !defined(KRATOS_EDGEBASED_DIFUSSION_SOLVER_H_INCLUDED)
#define  KRATOS_EDGEBASED_DIFUSSION_SOLVER_H_INCLUDED

//#define SPLIT_OSS
#define SYMM_PRESS


// System includes 
#include <string>
#include <iostream>
#include <algorithm>

// #include <omp.h>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
//#include "geometries/geometry.h"
#include "utilities/geometry_utilities.h"
#include "streamlinesdiffusion_application.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/timer.h"

namespace Kratos {

    template<unsigned int TDim, class MatrixContainer, class TSparseSpace, class TLinearSolver>
    class DifussionSolver {
    public:
	//name for the self defined structure
	typedef EdgesStructureType<TDim> CSR_Tuple;
	typedef std::vector<CSR_Tuple> EdgesVectorType;

	//name for row start and column index vectors
	typedef std::vector<unsigned int> IndicesVectorType;
	//defining matrix type for test calculations
	typedef std::vector< array_1d<double, TDim> > CalcVectorType;
	//defining type for local storage of nodal values
	typedef std::vector<double> ValuesVectorType;

	//defining types for matrix operations
	typedef typename TSparseSpace::MatrixType TSystemMatrixType;
	typedef typename TSparseSpace::VectorType TSystemVectorType;

	//constructor and destructor

	DifussionSolver(MatrixContainer& mr_matrix_container,
		ModelPart& mr_model_part
		)
	: mr_matrix_container(mr_matrix_container),
	  mr_model_part(mr_model_part)


	{
	};

	~DifussionSolver() {
	};

	//***********************************
	//function to initialize fluid solver

	void Initialize(double h_input, double r_input
		) {
	    KRATOS_TRY
		merror=0.0;
		mstep=0;
	    //get number of nodes
	    unsigned int n_nodes = mr_model_part.Nodes().size(); 
	    unsigned int n_edges = mr_matrix_container.GetNumberEdges();
	    //size data vectors
	    mWork.resize(n_nodes);
	    mDUMMY_UNKNOWN.resize(n_nodes);                            //debemos crear entonces mUNKNOWN
	    mToBeIntegratedNodeList.resize(n_nodes);
	    mDUMMY_UNKNOWN_OLD.resize(n_nodes);
	    mPOINT_SOURCE.resize(n_nodes); 
	    mCONDUCTIVITY.resize(n_nodes); 
	    mSMOOTHED_CONDUCTIVITY.resize(n_nodes);
	    mRequiredNodesForIntegration.resize(n_nodes);
	    //read  from Kratos
	    mr_matrix_container.FillScalarFromDatabase(TEMPERATURE,mDUMMY_UNKNOWN, mr_model_part.Nodes());              //cá ta!    ESTA LINEA ES NECESARIA!!!!!!!!!!!!!!
	    mr_matrix_container.FillScalarFromDatabase(POINT_HEAT_SOURCE,mPOINT_SOURCE, mr_model_part.Nodes());              //cá ta!    ESTA LINEA ES NECESARIA!!!!!!!!!!!!!!
	    mr_matrix_container.FillScalarFromDatabase(CONDUCTIVITY,mCONDUCTIVITY, mr_model_part.Nodes());       

	    mr_matrix_container.FillOldScalarFromDatabase(TEMPERATURE,mDUMMY_UNKNOWN_OLD, mr_model_part.Nodes());   
	    mFirstStep = true;

	    //loop to categorize boundary nodes
	    m_vector.resize(n_nodes);
	    for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
		    inode != mr_model_part.NodesEnd();
		    inode++) {
			inode->FastGetSolutionStepValue(CONDUCTIVITY)=1.0;
			int index = inode->FastGetSolutionStepValue(AUX_INDEX);	
			if (inode->IsFixed(TEMPERATURE)) {
			    mFixedTemperature.push_back(index);                                                             //fixed (dirichlet) BC
				}
	    }
/*
	    for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
		    inode != mr_model_part.NodesEnd();
		    inode++) {
			int index = inode->FastGetSolutionStepValue(AUX_INDEX);	
			if (inode->IsFixed(DUMMY_POINT_SOURCE)) {
			    mFixedSource.push_back(index);                                                             //fixed (dirichlet) BC
				}
	    }
*/
		m_vector.resize(n_nodes); //mass vector

	    for (int i_node = 0; i_node < n_nodes; i_node++) 
		{     //erase this later
			mRequiredNodesForIntegration[i_node]=1; } //meaning does not require storing g_i and m_i of this node,


      h =  h_input; 
	  reference_distance = r_input;
	  inv_int_dist = 1.0 / (reference_distance+h*0.5) ;
	  
	  /* 
      double average_nodal_area =  pow(h,3); 		//these three lines were used to estimate the number of nodes contained in a sphere of r_input. Now unused (it has the exact size)
      double average_nodes_in_radius = (4.0/3.0)  * 3.1416  * pow( reference_distance,3)  / average_nodal_area;
      full_vector_lenght = int ( 1.5  * n_nodes *average_nodes_in_radius ) + n_nodes;
      */


	  BuildPointClouds(); //(found below). it calculates the size of the phi vector and also fills it.
	  
	  //mlinearized_mass=0.0;
	  //now we calculate the new mass of each node using the Fi matrix we've just obtained
      int absolute_counter = 0;
	  for (int i_node=0; i_node != n_nodes ; ++i_node) {
		
		//double linearized_m_integral=0.0;
        double m_integral= 0.0; //g_vector[i_node];
        //KRATOS_WATCH(i_node);
        for (int counter=0; counter !=mToBeIntegratedNodeList[i_node]; ++counter) {
			m_integral +=  (mr_matrix_container.GetLumpedMass()[(mg_int_vector[absolute_counter])]) * (mg_double_vector[absolute_counter]) ;
			//KRATOS_WATCH(mg_double_vector[absolute_counter])
		    ++absolute_counter;
		}
        m_vector[i_node]=m_integral;
	  }
	   
      //we can further reduce the time consumed in each timestep by modifing the fi matrix , divining each row by the integrated mass we've just obtained:
      //actually it's a long vector
      //this way we avoid diving by the mass in each timestep.
      absolute_counter = 0;
      for (int i_node=0; i_node != n_nodes ; ++i_node) {
		
        double  mi =  m_vector[i_node]; //g_vector[i_node];
        for (int counter=0; counter !=mToBeIntegratedNodeList[i_node]; ++counter) {
			//if (i_node==10) std::cout << mg_int_vector[absolute_counter] << mg_double_vector[absolute_counter] << std::endl; 
			
			mg_double_vector[absolute_counter] *= (1.0 /  mi) ;
			//KRATOS_WATCH(mg_double_vector[absolute_counter])
		    ++absolute_counter;
		}
	  } 
       
        
	  KRATOS_CATCH("")
	}
	



	//*************************************************************************
	//function to solve fluid equations - fractional step 2: calculate pressure

	void SolveStep2(typename TLinearSolver::Pointer pLinearSolver) {               //used in each timestep
	    KRATOS_TRY

	    //PREREQUISITES

	    //allocate memory for variables
	    ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
	    int n_nodes = rNodes.size();
	    //unknown and right-hand side vector
	    TSystemVectorType g_vector, rhs;
	    g_vector.resize(n_nodes);
	    m_vector.resize(n_nodes);
	    rhs.resize(n_nodes);
	    array_1d<double, TDim> work_array;           //revisar si sirve dps
	    //read time step size from Kratos
	    ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
	    double delta_t = CurrentProcessInfo[DELTA_TIME];
	    double tiempo = CurrentProcessInfo[TIME];
	    ++mstep;
		const double sq_reference_distance= reference_distance*reference_distance;
//#ifdef _OPENMP

	    double time_inv = 1.0/delta_t;

	   mr_matrix_container.FillScalarFromDatabase(TEMPERATURE,mDUMMY_UNKNOWN, mr_model_part.Nodes());         
	   mr_matrix_container.FillOldScalarFromDatabase(CONDUCTIVITY,mCONDUCTIVITY, mr_model_part.Nodes());       
	   mr_matrix_container.FillScalarFromDatabase(POINT_HEAT_SOURCE,mPOINT_SOURCE, mr_model_part.Nodes());              
       mr_matrix_container.FillOldScalarFromDatabase(TEMPERATURE,mDUMMY_UNKNOWN_OLD, mr_model_part.Nodes());   
	   
	   ModelPart::NodesContainerType::iterator it_begin = rNodes.begin();
	   
	   //GUARDA CON ESTA LíNEA!
	   //for (int i_node = 0; i_node < n_nodes; i_node++) { mPOINT_SOURCE[i_node]= mPOINT_SOURCE[i_node] * ( 0.6 + 0.5 * cos(tiempo*20.0)); }

	   calculate_smoothed_conductivity();
	   
		//first we calculate the nodal gs in the current time step
	   for (int i_node = 0; i_node < n_nodes; i_node++) {
			if (mRequiredNodesForIntegration[i_node]==1)   {
				g_vector[i_node]  = 0.0;
				//ModelPart::NodesContainerType::iterator node_it = it_begin + i_node  ;  
				double node_conductivity = mSMOOTHED_CONDUCTIVITY[i_node];
				double l_ii = 0.0;
				//loop over all neighbours
				for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++) {
					unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
					CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];

					//compute laplacian operator
					double sum_l_ikjk;
					edge_ij.CalculateScalarLaplacian(sum_l_ikjk);
					
					//we modify sum_l_ikjk to take into account the mean diffusivity;
					sum_l_ikjk *= (mSMOOTHED_CONDUCTIVITY[j_neighbour] + node_conductivity)*0.5;
					
					g_vector[i_node]  -= sum_l_ikjk * mDUMMY_UNKNOWN_OLD[j_neighbour] ;   //
					l_ii -= sum_l_ikjk;
					//KRATOS_WATCH(mCONDUCTIVITY[j_neighbour]*mCONDUCTIVITY[i_node]);
				}

 				
				g_vector[i_node] -= l_ii *  mDUMMY_UNKNOWN_OLD[i_node]  ;
				
				

				
				
				g_vector[i_node]  += mPOINT_SOURCE[i_node];
				
				//LET's try deviding (later, when assembling the whole full vector of g) the conductivity of the neighbour node and this one.
				//g_vector[i_node] *=  mCONDUCTIVITY[i_node];
				
				//m_vector[i_node] = mr_matrix_container.GetLumpedMass()[i_node];
			}
        }
       //listo, ya tengo todos los valores de g
       
       

              ///************************************
           for (unsigned int i_pressure = 0; i_pressure < mFixedTemperature.size(); i_pressure++) {           //FIJAMOS CONDICIONES DE CONTORNO
               unsigned int i_node = mFixedTemperature[i_pressure];
               //KRATOS_WATCH(mFixedTemperature[i_pressure]);
               //mL(i_node, i_node) = huge;    
               //Geometry<Node < 3 > >& face_geometry = (1)->GetGeometry();                                                                
               //g_vector[i_node] =0.0; // READD
           }
       ///*******************************************************

       
      
      //ahora tengo q sumar las contribuciones de los gj dentro del radio de integración
      int absolute_counter= 0;
      for (int i_node = 0; i_node < n_nodes; i_node++) {     
		ModelPart::NodesContainerType::iterator node_it = it_begin + i_node  ;  
		double node_conductivity = mCONDUCTIVITY[i_node];
        double g_integral= 0.0; //g_vector[i_node];
        for (int counter=0; counter !=mToBeIntegratedNodeList[i_node]; ++counter) {
			g_integral +=  (g_vector[(mg_int_vector[absolute_counter])]) * (mg_double_vector[absolute_counter])  ;
		    ++absolute_counter;
		}
        node_it->FastGetSolutionStepValue(G_VALUE) =   g_integral  * delta_t; // / mCONDUCTIVITY[i_node]; // + mDUMMY_UNKNOWN_OLD[i_node]; //i removed the 1/mass since it is already inside the modified fi matrix

      }

	  //now we reset the values to the imposed ones
      for (unsigned int i_pressure = 0; i_pressure < mFixedTemperature.size(); i_pressure++) {           //FIJAMOS CONDICIONES DE CONTORNO
               unsigned int i_node = mFixedTemperature[i_pressure];
               //KRATOS_WATCH(mFixedTemperature[i_pressure]);
               //mL(i_node, i_node) = huge;    
               //Geometry<Node < 3 > >& face_geometry = (1)->GetGeometry();                                                                
              //mDUMMY_UNKNOWN[i_node]=mDUMMY_UNKNOWN_OLD[i_node]; // the values are already loaded in mDUMMY_UNKNOWN
              ModelPart::NodesContainerType::iterator node_it = it_begin + i_node  ; 
              node_it->FastGetSolutionStepValue(G_VALUE) =  0.0; //node_it->FastGetSolutionStepValue(G_VALUE)/100.0; //(node_it->FastGetSolutionStepValue(G_VALUE))/10.0;
      }


	  //mr_matrix_container.WriteScalarToDatabase(TEMPERATURE, mDUMMY_UNKNOWN, rNodes);                           //cambiar
 /*
	 double smallest_mass=100000000000000000.0;
	 double biggest_mass=0.0;
     double m_total=0.0;
     int node_of_smallest_mass=0;
	 for (int i_node = 0; i_node < n_nodes; i_node++) {     
		//ModelPart::NodesContainerType::iterator node_it = it_begin + i_node  ;  
        //double m_total= 0.0; //g_vector[i_node];
        m_total+=m_vector[i_node];
        if (m_vector[i_node]<smallest_mass) {smallest_mass=m_vector[i_node]; node_of_smallest_mass=i_node+1;}
        if (m_vector[i_node]>biggest_mass) biggest_mass=m_vector[i_node]; 

      }
      m_total *= 1.0/n_nodes;
	  //double predicted= m_total*m_total;
	  */
	  
	  /*
	  //KRATOS_WATCH(m_total);
	  KRATOS_WATCH(mlinearized_mass);
	  KRATOS_WATCH(mcentral_distance);
	  KRATOS_WATCH(msmallest_mass);
	  KRATOS_WATCH(mbiggest_mass);
	  //KRATOS_WATCH(node_of_smallest_mass);

	
		//calculemos el error (sólo para este ejemplo)
		double zero_time=5.0;
		//double zero_time=0.05;
		//double perm=0.05;
		double perm=0.95; 
		double mod_time=zero_time+tiempo; //-2.0*delta_t;
		double max_temp = 1000.0 / (4.0 * 3.14159 * perm * mod_time);
		//double max_temp = 1 / (4.0 * 3.14159 * perm * mod_time);
		double sum_errores=0.0;
		int tot_elem_errores=0;
		
		for (int i_node = 0; i_node < n_nodes; i_node++) {
			ModelPart::NodesContainerType::iterator node_it = it_begin + i_node  ;             
			array_1d<double, 3>  coordLocal;
			coordLocal[0] = node_it->X();
			coordLocal[1] = node_it->Y();
			if (mDUMMY_UNKNOWN[i_node]<(-0.01)){
				KRATOS_WATCH(i_node);
				KRATOS_ERROR(std::logic_error, "i stopped here!", "");
			}
			
			//double radius_2= pow((coordLocal[0]- 25.0),2)+pow((coordLocal[1]-25.0),2);
			double radius_2= pow((coordLocal[0]),2)+pow((coordLocal[1]),2);
	
			double temperatura=1000.0 * pow(2.7182818,(-radius_2 / (mod_time * 4.0 * perm ) ) ) / (4.0 * 3.14159 * perm * mod_time);
			//double temperatura= pow(2.7182818,(-radius_2 / (mod_time * 4.0 * perm ) ) ) / (4.0 * 3.14159 * perm * mod_time);
			double errorr = 0.0;
			if (radius_2<80.0)
				errorr = (mDUMMY_UNKNOWN[i_node] - temperatura) / max_temp;
				//errorr = ((mass_integral-exact_mass_integral)/exact_mass_integral) * 100.0;
				//(it_begin+i_node)->GetSolutionStepValue(Q_FLUX)[1] = errorr;
				//sum_errores+=fabs(errorr); ++tot_elem_errores;
				sum_errores+=pow(errorr,2); ++tot_elem_errores;
		}
		double errorrr = sqrt(sum_errores/tot_elem_errores); KRATOS_WATCH(errorrr);
		double sq_errorrr=sum_errores/tot_elem_errores;

		if (merror<errorrr){
			//KRATOS_WATCH(merror);
			//KRATOS_ERROR(std::logic_error, "i stopped here!", "");
			merror=errorrr;
			merrortime=tiempo;
		}
		sq_errorrr=pow(merror,2);
		KRATOS_WATCH(merror);
		KRATOS_WATCH(sq_errorrr)
		KRATOS_WATCH(merrortime);	
		*/
		
		//merror=errorrr; 
		//KRATOS_WATCH(mIntegredNodes);
	
	
    KRATOS_CATCH("")
    
	}











		void BuildPointClouds()
 		{
			KRATOS_TRY
			KRATOS_WATCH(TDim)
 			//**********************************************************************
// 			numofpts_origin = rOrigin_ModelPart.Nodes().size();
// 			numofpts_destination = rDestination_ModelPart.Nodes().size();
// 					
 			//*******************************************************************
			//properties to be used in the generation
			Properties::Pointer properties = mr_model_part.GetMesh().pGetProperties(1);

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

			mlinearized_mass=0.0;
			mcentral_distance=0.0;
			msmallest_mass=10000000000000000000.0;
			mbiggest_mass=0.0;
			int nodes_to_calculate_mass=0;
			
			for(ModelPart::NodesContainerType::iterator node_it = mr_model_part.NodesBegin();
						node_it != mr_model_part.NodesEnd(); ++node_it)
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
            //int hola = list_of_new_nodes.begin() - list_of_new_nodes.end();
		
			//work arrays
			Node<3> work_point(0,0.0,0.0,0.0);
			unsigned int MaximumNumberOfResults = 10000;
			PointVector Results(MaximumNumberOfResults);
			//PointIterator Results(MaximumNumberOfResults);
			//PointVector global_point_vector(full_vector_lenght);
			DistanceVector ResultsDistances(MaximumNumberOfResults);
			array_1d<double,TDim+1> N; //Shape functions vector//
			boost::timer search_and_interpolation_time;
			//loop over all of the elements in the "old" list to perform the interpolation
			
			int full_vector_lenght = 0;
			//first we go through all the nodes to determine the size of the vector containing the phi values and Ids
			for( ModelPart::NodesContainerType::iterator node_it = mr_model_part.NodesBegin();
						node_it != mr_model_part.NodesEnd(); node_it++)
			{
				work_point.X() = node_it->X(); work_point.Y() = node_it->Y(); work_point.Z() = node_it->Z();

				//find all of the new nodes within the radius
				int number_of_points_in_radius;
				
				//look between the new nodes which of them is inside the radius of the circumscribed cyrcle
				number_of_points_in_radius = nodes_tree.SearchInRadius(work_point, reference_distance, Results.begin(),
 						ResultsDistances.begin(),  MaximumNumberOfResults);
				
				full_vector_lenght += number_of_points_in_radius;
			}
			
			//now we resize the vector
			KRATOS_WATCH(full_vector_lenght)
			mg_int_vector.resize(full_vector_lenght); 
			mg_double_vector.resize(full_vector_lenght);
			double temp_mass=0.0;
			double central_temp_mass;
			//done. now we can make a second loop but this time saving the Ids and phi functions that we'll need later.
			int absolute_counter=0;
			
			
			
			for( ModelPart::NodesContainerType::iterator node_it = mr_model_part.NodesBegin();
						node_it != mr_model_part.NodesEnd(); node_it++)
			{

				work_point.X() = node_it->X(); work_point.Y() = node_it->Y(); work_point.Z() = node_it->Z();
				//KRATOS_WATCH(work_point.X());
				//find all of the new nodes within the radius
				int number_of_points_in_radius;
				
				//look between the new nodes which of them is inside the radius of the circumscribed cyrcle
				number_of_points_in_radius = nodes_tree.SearchInRadius(work_point, reference_distance, Results.begin(),
 						ResultsDistances.begin(),  MaximumNumberOfResults);

				int node_number = node_it -  mr_model_part.NodesBegin();                  //saving the number of neighbours required in this node.
				mToBeIntegratedNodeList[node_number] = number_of_points_in_radius;
                temp_mass=0.0;
                //central_temp_mass;
                int index=0;
                //KRATOS_WATCH("another node")
                vector<bool> neighbourhood(6, false);
                //bool useful_node_to_calculate_mass=true;
                
				for( PointIterator it_found = Results.begin(); it_found != Results.begin() + number_of_points_in_radius; it_found++)
				{	
					mg_int_vector[absolute_counter] =  (*it_found)->Id() - 1; //Results[index] - list_of_new_nodes.begin();
					
					
					
					double x_point=(mr_model_part.NodesBegin()+(*it_found)->Id() - 1)->X();
					double y_point=(mr_model_part.NodesBegin()+(*it_found)->Id() - 1)->Y();
					double z_point=(mr_model_part.NodesBegin()+(*it_found)->Id() - 1)->Z();
					double distance_to_compare=reference_distance*0.1;
					if (( x_point-work_point.X())>distance_to_compare) neighbourhood[0]=true;
					if (( y_point-work_point.Y())>distance_to_compare) neighbourhood[1]=true;
					if ((-x_point+work_point.X())>distance_to_compare) neighbourhood[2]=true;
					if ((-y_point+work_point.Y())>distance_to_compare) neighbourhood[3]=true;
					if (( z_point-work_point.Z())>distance_to_compare) neighbourhood[4]=true;
					if ((-z_point+work_point.Z())>distance_to_compare) neighbourhood[5]=true;
					

					
					//KRATOS_WATCH(x_point);
					double distance = sqrt(ResultsDistances[index])*inv_int_dist;
					
					double correction_factor;
					double rd_over_h = reference_distance/h;
					if (rd_over_h>2.0){ 
						if (rd_over_h<6.0) correction_factor = 1.0 -  ( (-0.5 + 0.25 * reference_distance/h ) *  distance) ; //when the  reference_distance (integration radius R )  is larger than 2h the phi function must be lowered.
						else  correction_factor = 1.0 - distance;
					}
					else {correction_factor=1.0;}
					correction_factor=1.0-distance;
					 
					//this factor could be set to (1- distance), avoiding the mesh dependancy, but it would not be as high as possible for low R/h factors.
					
					//KRATOS_WATCH(ResultsDistances[index]);
					//double weight = (1.0  - 9.876 * pow(distance,2) + 22.915  * pow(distance,3) - 21.042 * pow(distance,4 ) +  7.048 * pow(distance,5 )) * correction_factor; 
					//double weight= 1.0- 1.0 * distance - 1.0 * pow(distance,2) + 1.0 * pow(distance,3);
					//double weight = - 0.25 * pow(distance,3) + 1.25 * pow(distance,2) - 2.0 *distance + 1.0 ;
					
					double weight= 1.0- 2.0*distance + pow(distance,2);
					mg_double_vector[absolute_counter] = weight;
					
					if (TDim==2) 
					{
						if (ResultsDistances[index]<0.0000000001) {
							central_temp_mass =  pow(mr_matrix_container.GetLumpedMass()[(mg_int_vector[absolute_counter])],(1.0/2.0));
							temp_mass += pow(mr_matrix_container.GetLumpedMass()[(mg_int_vector[absolute_counter])],(1.0/2.0));
							//mlinearized_mass += pow(mr_matrix_container.GetLumpedMass()[(mg_int_vector[absolute_counter])],(1.0/2.0));
							mcentral_distance += pow(mr_matrix_container.GetLumpedMass()[(mg_int_vector[absolute_counter])],(1.0/2.0)); }
						else 
						{	
							temp_mass += 2.0 * (mr_matrix_container.GetLumpedMass()[(mg_int_vector[absolute_counter])]) * weight / (6.28318531*sqrt(ResultsDistances[index])) ;
							//mlinearized_mass += 2.0 * (mr_matrix_container.GetLumpedMass()[(mg_int_vector[absolute_counter])]) * weight / (6.28318531*sqrt(ResultsDistances[index])) ;
						}
					}
					else //TDim==3
					{
						if (ResultsDistances[index]<0.0000000001) {
							central_temp_mass = pow(mr_matrix_container.GetLumpedMass()[(mg_int_vector[absolute_counter])],(1.0/3.0));
							temp_mass += pow(mr_matrix_container.GetLumpedMass()[(mg_int_vector[absolute_counter])],(1.0/3.0));
							//mlinearized_mass += pow(mr_matrix_container.GetLumpedMass()[(mg_int_vector[absolute_counter])],(1.0/3.0));
							mcentral_distance += pow(mr_matrix_container.GetLumpedMass()[(mg_int_vector[absolute_counter])],(1.0/3.0)); }
						else 
						{
							temp_mass += 2.0 * (mr_matrix_container.GetLumpedMass()[(mg_int_vector[absolute_counter])]) * weight / (12.5663706*ResultsDistances[index]) ;
							//mlinearized_mass += 2.0 * (mr_matrix_container.GetLumpedMass()[(mg_int_vector[absolute_counter])]) * weight / (12.5663706*ResultsDistances[index]) ;
						}
					}

					//KRATOS_WATCH(
					//KRATOS_WATCH(weight);
					++absolute_counter;
					++index;
				 }
				 bool useful_node_to_calculate_mass=true;
				 int positive_points=6;
				 for(int i=0; i!=6; i++)
				 {
					if (neighbourhood[i]==false) {useful_node_to_calculate_mass=false; positive_points += -1;}
				 }
				 
				 if(useful_node_to_calculate_mass==true)
				 {
					 mlinearized_mass +=temp_mass;
					 nodes_to_calculate_mass++;
				 }
				 else
				 {
					 temp_mass-= central_temp_mass;
					 if (positive_points==5) temp_mass *= 2.0;
					 else if (positive_points==4) temp_mass *=4.0;
					 else temp_mass *=8.0;
					 temp_mass +=central_temp_mass;
					 mlinearized_mass +=temp_mass;
					 nodes_to_calculate_mass++;
					 useful_node_to_calculate_mass==true;
				 }
				 
				 
				    
				if ((absolute_counter-1) >= full_vector_lenght) KRATOS_ERROR(std::logic_error, "oops, more neighbours than expected!", "");
				if (temp_mass>mbiggest_mass) mbiggest_mass=temp_mass;
				if (temp_mass<msmallest_mass ) msmallest_mass=temp_mass;
				
				
				
			}
			
			KRATOS_WATCH(absolute_counter);	
			mlinearized_mass*=1.0/	nodes_to_calculate_mass;
			mcentral_distance*=1.0/	mr_model_part.Nodes().size();
			std::cout << "search and interpolation time " << search_and_interpolation_time.elapsed() << std::endl;
			
			
			
			
			
			ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes(); 
			//mIntegredNodes=0;
			ModelPart::NodesContainerType::iterator it_begin = rNodes.begin();
			unsigned int n_nodes = mr_model_part.Nodes().size(); 
			
			for (int i_node = 0; i_node < n_nodes; i_node++) {     
		//mToBeIntegratedNodeList[i_node]=0;
				ModelPart::NodesContainerType::iterator node_it = it_begin + i_node  ; 

				//work_point.X() = node_it->X(); work_point.Y() = node_it->Y(); work_point.Z() = node_it->Z();        
				array_1d<double, 3>  coordLocal;
				coordLocal[0] = node_it->X();
				coordLocal[1] = node_it->Y();
				coordLocal[2] = node_it->Z();
				
				
				//m_vector[i_node] = mr_matrix_container.GetLumpedMass()[i_node];
				
				//double weight;
				double local_sq_distance=0.0;
				double sq_distance=0.0;  
				unsigned int conectividades [10000];    //guarda las conectividades para no repetir aportes
				//unsigned int aportantes [10000];    //guarda las conectividades para no repetir aportes
				//for ( int contador = 0; contador != 10000; contador++) { conectividades[contador]=0; aportantes[contador]=0; }
				//unsigned int cant_neigh = 0;
				//conectividades[0] = i_node ; //la posición cero es el nodo local
				//++cant_neigh; // el nodo local ya está ocupando la primer posición
				
				//loop over all neighbours
				for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++) {
					unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
					
					ModelPart::NodesContainerType::iterator node_it1 = it_begin + j_neighbour  ;             
					array_1d<double, 3>  coordNeigh1;
					coordNeigh1[0] = node_it1->X();
					coordNeigh1[1] = node_it1->Y();
					coordNeigh1[2] = node_it1->Z();
					
					//conectividades[cant_neigh] = j_neighbour ; //we save this node
					//++cant_neigh;
					
					local_sq_distance =  pow((coordNeigh1[0] - coordLocal[0]),2) +  pow((coordNeigh1[1] - coordLocal[1]),2) +  pow((coordNeigh1[2] - coordLocal[2]),2)  ;
					if (local_sq_distance>sq_distance)	{ //ok, esto significa que este nodo es útil y que además sus vecinos puede que lo sean también
						sq_distance=local_sq_distance;
						node_it->GetValue(FLAG_VARIABLE)=sqrt(sq_distance);
					}

				}
			}
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			KRATOS_CATCH("")
		}





    private:
    

    
    void calculate_smoothed_conductivity()
    {
		//ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
	    int n_nodes =  mr_model_part.Nodes().size();
	  int absolute_counter= 0;
      for (int i_node = 0; i_node < n_nodes; i_node++) {     
		ModelPart::NodesContainerType::iterator node_it =  mr_model_part.NodesBegin() + i_node  ;  
		//double node_conductivity = mCONDUCTIVITY[i_node];
		double extra_weight = double(mToBeIntegratedNodeList[i_node])*2.0;//1.0;
        double smoothed_conductivity= 0.0; //g_vector[i_node];
        for (int counter=0; counter !=mToBeIntegratedNodeList[i_node]; ++counter) {
			smoothed_conductivity +=  (mCONDUCTIVITY[(mg_int_vector[absolute_counter])])  ;
			if 	(mg_int_vector[absolute_counter]==i_node) smoothed_conductivity +=  extra_weight*(mCONDUCTIVITY[(mg_int_vector[absolute_counter])])  ;
		    ++absolute_counter;
		}
        node_it->FastGetSolutionStepValue(PRESSURE) =   smoothed_conductivity/(3.0*double(mToBeIntegratedNodeList[i_node])); // / mCONDUCTIVITY[i_node]; // + mDUMMY_UNKNOWN_OLD[i_node]; //i removed the 1/mass since it is already inside the modified fi matrix
		mSMOOTHED_CONDUCTIVITY[i_node]=smoothed_conductivity/(3.0*double(mToBeIntegratedNodeList[i_node]));
      }
		
		
	}
		
    
    
    
    
    
    
    
	MatrixContainer& mr_matrix_container;
	ModelPart& mr_model_part;

	CalcVectorType mWork; //, mDUMMY_UNKNOWN, mvel_n1 ; //, mvel_n, mvel_n1, mx;
	//pressure vector p at time steps n and n+1
	ValuesVectorType mDUMMY_UNKNOWN, mPOINT_SOURCE, mDUMMY_UNKNOWN_OLD, mCONDUCTIVITY, mSMOOTHED_CONDUCTIVITY;
	TSystemVectorType m_vector; //mPn, mPn1;                                                  //poner algo acá !//////////////////////      PONER ALGO ACA
	
	//minimum length of the edges surrounding edges surrounding each nodal point
	ValuesVectorType mHmin;
	

	//flag for first time step
	bool mFirstStep;
	int mstep;
	double reference_distance, h, inv_int_dist, merror, merrortime, mlinearized_mass, mcentral_distance, msmallest_mass, mbiggest_mass; //o tiene que ser const?
	int full_vector_lenght;
	int mIntegredNodes;
	//flag to differentiate interior and boundary nodes
	ValuesVectorType mNodalFlag;
	//lists of nodes with different types of boundary conditions
	//IndicesVectorType // mFixedTemperature; //mSlipBoundaryList, mPressureOutletList, mFixedVelocities;             //FIXED TEMPERATURES !!!!! FIJARSE COMO ERA ORIGINALMENTE
	//CalcVectorType mFixedVelocitiesValues;
	ValuesVectorType mFixedTemperature;
	ValuesVectorType mFixedSource;

	vector<int> mg_int_vector;
	vector<double> mg_double_vector;


	//variables for edge BCs
	IndicesVectorType medge_nodes;
	CalcVectorType medge_nodes_direction;
	IndicesVectorType mcorner_nodes, mToBeIntegratedNodeList, mRequiredNodesForIntegration;






    };
} //namespace Kratos

#endif //KRATOS_EDGEBASED_DIFUSSION_SOLVER_H_INCLUDED defined


