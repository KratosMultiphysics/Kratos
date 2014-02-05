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
//   Date:                $Date: 2008-10-13 06:58:23 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_ESTIMATE_TIME_STEP )
#define  KRATOS_ESTIMATE_TIME_STEP



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
//#include "geometries/tetrahedra_3d_4.h"
#include "geometries/point.h"
#include "thermo_mechanical_application.h"
// #include "custom_conditions/environment_contact.h"
//#include "includes/variables.h"
#include "utilities/openmp_utils.h"
#include "includes/convection_diffusion_settings.h"




namespace Kratos
{
 	
	template<unsigned int TDim>
	class EstimateTimeStep
	{
	public:

		//**********************************************************************************************
		//**********************************************************************************************
		//

	  double ComputeDt(ModelPart& ThisModelPart, const double dist_max, const double CFL, const double dt_min ,const double dt_max  )
	  {			
	    KRATOS_TRY

	      const unsigned int NumNodes = TDim +1;

	      int NumThreads = OpenMPUtils::GetNumThreads();
	      OpenMPUtils::PartitionVector ElementPartition;
	      OpenMPUtils::DivideInPartitions(ThisModelPart.NumberOfElements(),NumThreads,ElementPartition);

	      std::vector<double> MaxProj(NumThreads,0.0);
      // NumThreads = 1;
	      #pragma omp parallel shared(MaxProj)
	      {
		  int k = OpenMPUtils::ThisThread();
		  ModelPart::ElementIterator ElemBegin = ThisModelPart.ElementsBegin() + ElementPartition[k];
		  ModelPart::ElementIterator ElemEnd = ThisModelPart.ElementsBegin() + ElementPartition[k+1];

		  double& rMaxProj = MaxProj[k];

		  double Area;
		  array_1d<double, NumNodes> N;
		  array_1d<double, NumNodes> dist_vec;
		  boost::numeric::ublas::bounded_matrix<double, NumNodes, TDim> DN_DX;

		  for( ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
		  {
		      // Get the element's geometric parameters
		      Geometry< Node<3> >& rGeom = itElem->GetGeometry();
				  double ele_dist = 0.0;
				  for (unsigned int kk = 0; kk < rGeom.size(); kk++){
					  double dist = rGeom[kk].FastGetSolutionStepValue(DISTANCE);
					  dist_vec[kk] = dist;
					ele_dist +=  fabs(dist);
				  }
		    ele_dist /= NumNodes;
		    if(ele_dist <= dist_max){
		      GeometryUtils::CalculateGeometryData(rGeom, DN_DX, N, Area);

		      // Elemental Velocity
		      array_1d<double,TDim> ElementVel = N[0]*itElem->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
		      for (unsigned int i = 1; i < NumNodes; ++i)
			  noalias(ElementVel) += N[i]*rGeom[i].FastGetSolutionStepValue(VELOCITY);

		      //compute normal velocity
		      array_1d<double, TDim>  grad_dist = prod(trans(DN_DX),dist_vec);
		      double norm_grad_dist = norm_2(grad_dist); //grad_dist[0]*grad_dist[0] + grad_dist[1]*grad_dist[1] + grad_dist[2]*grad_dist[2];
      // 		norm_grad_dist = sqrt(norm_grad_dist);
		      
		      double normal_speed = inner_prod(ElementVel,grad_dist);
		      if(norm_grad_dist > 0.0){
			normal_speed /= norm_grad_dist;
			noalias(ElementVel) = (normal_speed/norm_grad_dist)*grad_dist;
		      }
		      else 
			noalias(ElementVel) = ZeroVector(3);
		      
		      
      // 		// Velocity norm
      //                 double VelNorm = ElementVel[0]*ElementVel[0];
      //                 for (unsigned int d = 1; d < TDim; ++d)
      //                     VelNorm += ElementVel[d]*ElementVel[d];
      //                 VelNorm = sqrt(VelNorm);


		      // Maximum element size along the direction of velocity
		      for (unsigned int i = 0; i < NumNodes; ++i)
		      {
			  double Proj = 0.0;
			  for (unsigned int d = 0; d < TDim; ++d)
			      Proj += ElementVel[d]*DN_DX(i,d);
			  Proj = fabs(Proj);
			  if (Proj > rMaxProj) rMaxProj = Proj;
		      }
		    }
		  }
	      }

	      // Obtain the maximum projected element size (compare thread results)
	      double Max = 0.0;
	      for (int k = 0; k < NumThreads; ++k)
		  if (Max < MaxProj[k]) Max = MaxProj[k];

	      // Dt to obtain desired CFL
	      double dt = CFL / Max;
	      if(dt > dt_max)
		  dt = dt_max;
	      else if(dt < dt_min)
		  dt = dt_min;

	      //perform mpi sync if needed
	      double global_dt = dt;
	      ThisModelPart.GetCommunicator().MinAll(global_dt);
	      dt = global_dt;

	      return dt;

	    KRATOS_CATCH("")	
	  }
	  
         /* Compute solidificatio nand cooling DT */
	  double ComputeSolidificationCoolingDt(ModelPart& ThisModelPart, 
						const double solidification_percent, 
						const double max_cooling_delta_temp, 
						const double dt_min,
						const double dt_max)
	  {			
	    KRATOS_TRY
	    
//	   const int NumNodes = TDim +1;
// 	   is_cold = 0;
	   
	   const double current_dt = ThisModelPart.GetProcessInfo()[DELTA_TIME];
	   const int is_solidified = ThisModelPart.GetProcessInfo()[IS_SOLIDIFIED];
	   int global_is_solidified = is_solidified;
// 	   ThisModelPart.GetCommunicator().MinAll(global_is_solidified);	   
	   
	   if(global_is_solidified == 0)
	   {
         //pre solidification Dt
        int is_hot = CheckMaxTemperature(ThisModelPart);
		if( is_hot == 1 ){
			KRATOS_WATCH("<<<<<<<<<<<<<<<<<<<<pre solidification Dt>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
			double max_presolodification_delta_tem = 10.0;
			int node_size = ThisModelPart.Nodes().size();
			double max_delta_temp = 0.0;
			for (int ii = 0; ii < node_size; ii++)
			   {
					 ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + ii;

			 double current_temp = it->FastGetSolutionStepValue(TEMPERATURE);
			 double old_temp = it->FastGetSolutionStepValue(TEMPERATURE,1);	
			 current_temp -= old_temp;
		 
			 if( current_temp > max_delta_temp)
			   max_delta_temp = current_temp;
			 
			   }
	       
			ThisModelPart.GetCommunicator().MaxAll(max_delta_temp);
	   
			if( max_delta_temp > 0.0 ){
			  double new_delta_time = max_presolodification_delta_tem / max_delta_temp;
			  new_delta_time *= current_dt; 

			  if( new_delta_time > dt_max)
			new_delta_time = dt_max;
			  else if( new_delta_time < dt_min)
			new_delta_time = dt_min;

		      
			  return new_delta_time;
			}
			else 
			{
			  return current_dt;
			}

		}
		else//solidification Dt
		{

	    double current_solidified_volume = 0.0;
	    double old_solidified_volume = 0.0;
	    double tot_vol = 0.0;
	    
	    int node_size = ThisModelPart.Nodes().size();	    
	    for (int ii = 0; ii < node_size; ii++)
	       {
                 ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + ii;
		 double vol = it->FastGetSolutionStepValue(NODAL_VOLUME);
		 double current_S = it->FastGetSolutionStepValue(SOLID_FRACTION);
		 double old_S = it->FastGetSolutionStepValue(SOLID_FRACTION,1);	
		 
		 current_solidified_volume += vol*current_S;
		 old_solidified_volume += vol*old_S;
		 tot_vol += vol;

		 //filling solidifiacation time
		 double is_visited = it->FastGetSolutionStepValue(IS_VISITED);
		 if(is_visited == 0.0 && current_S == 1.0){
			 it->FastGetSolutionStepValue(IS_VISITED) = 1.0;

			 double solid_time =  ThisModelPart.GetProcessInfo()[TIME];
			 double modulus = ThisModelPart.GetProcessInfo()[K0];

			 it->FastGetSolutionStepValue(SOLIDIF_TIME) = solid_time;
			 it->FastGetSolutionStepValue(SOLIDIF_MODULUS) = modulus*sqrt(solid_time);
		 }
		 
	       }
	       
	    ThisModelPart.GetCommunicator().SumAll(current_solidified_volume);	 
	    ThisModelPart.GetCommunicator().SumAll(old_solidified_volume);
	    ThisModelPart.GetCommunicator().SumAll(tot_vol);

	    if(tot_vol == 0.0) 
                KRATOS_ERROR(std::logic_error, "inside ComputeSolidificationCoolingDt: total volume is zero!", "")
	      
	    if(current_solidified_volume == tot_vol){
	      ThisModelPart.GetProcessInfo()[IS_SOLIDIFIED] = 1;

	      return current_dt;}
	    else
	    {
	      double  delta_solid = current_solidified_volume - old_solidified_volume;
			      
	      if( delta_solid > 0.0 )
	      {
		delta_solid /= tot_vol;
			      KRATOS_WATCH(delta_solid);
		double new_dt = (solidification_percent /  delta_solid) * current_dt;
		if( new_dt > dt_max)
		  new_dt = 1.5 * current_dt;//dt_max;
		else if( new_dt < dt_min)
		  new_dt = dt_min;

		
		return new_dt;	      
	      }
	      else{
		return current_dt;}
	    }
	   }
	    
	   }
	   //coling delta_t
	   else
	   {
        double cooling_dt_max = 30.0*dt_max;
	    int node_size = ThisModelPart.Nodes().size();
	    double max_delta_temp = 0.0;
	    for (int ii = 0; ii < node_size; ii++)
	       {
                 ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + ii;

		 double current_temp = it->FastGetSolutionStepValue(TEMPERATURE);
		 double old_temp = it->FastGetSolutionStepValue(TEMPERATURE,1);	
		 current_temp -= old_temp;
		 
		 if( current_temp > max_delta_temp)
		   max_delta_temp = current_temp;
			 
	       }
	       
	    ThisModelPart.GetCommunicator().MaxAll(max_delta_temp);
	   
	    if( max_delta_temp > 0.0 ){
	      double new_delta_time = max_cooling_delta_temp / max_delta_temp;
	      new_delta_time *= current_dt; 

	      if( new_delta_time > cooling_dt_max)
		new_delta_time = cooling_dt_max;
	      else if( new_delta_time < dt_min)
		new_delta_time = dt_min;

		      
	      return new_delta_time;
	    }
	    else 
	    {
// 	      is_cold = 1;
	
	      return current_dt;
	    }
	   }
	    KRATOS_CATCH("")	
	  }
	  /////FIND AN ESTIMATION FOR SOLIDIFICATION TIME 
	  double EstimateSolidificationTime(ModelPart& ThisModelPart)
	  {			
	    KRATOS_TRY	  
	    
	    double solidification_time = 0.0;
	    
//             ConvectionDiffusionSettings::Pointer my_settings = ThisModelPart.GetProcessInfo().GetValue(CONVECTION_DIFFUSION_SETTINGS);
// 	
//             const Variable<double>& rDensityVar = my_settings->GetDensityVariable();	    
//             const Variable<double>& rTransferCoefficientVar = my_settings->GetTransferCoefficientVariable();
	    
	    
// 	    ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + 100;
	    const double density = ThisModelPart.GetProcessInfo()[DENSITY];
	    const double cc= ThisModelPart.GetProcessInfo()[SPECIFIC_HEAT];
	    const double htc= ThisModelPart.GetProcessInfo()[HTC];	    
    	    
	    const double TT_solid = ThisModelPart.GetProcessInfo()[SOLID_TEMPERATURE];
	    const double TT_liquid = ThisModelPart.GetProcessInfo()[FLUID_TEMPERATURE];		  
	    
	    const double LL = ThisModelPart.GetProcessInfo()[LATENT_HEAT];	

	    double tot_vol = 0.0;
            double tot_area = 0.0;	    
	    int node_size = ThisModelPart.Nodes().size();	    
	    for (int ii = 0; ii < node_size; ii++)
	       {
                 ModelPart::NodesContainerType::iterator it_nd = ThisModelPart.NodesBegin() + ii;
		 double vol = it_nd->FastGetSolutionStepValue(NODAL_VOLUME);
		 tot_vol += vol;
		 double area =  it_nd->FastGetSolutionStepValue(NODAL_PAUX);
		 tot_area += area;
	       }
	    
	    if ( tot_area == 0.0 || tot_vol == 0.0)
	      KRATOS_ERROR(std::invalid_argument,"AREA or VOLUME is Zero", "");
	    
	    
	    solidification_time = 2.0 * density * ( cc * ( TT_liquid - TT_solid) + LL) / (htc * TT_solid);
	    solidification_time *= pow(tot_vol/tot_area , 0.8);
	    
	    return solidification_time;
	    
	    KRATOS_CATCH("")	
	  }	

	 double CheckStopTemperature(ModelPart& ThisModelPart,const double stop_temperature)
	  {			
	    KRATOS_TRY	
		const double TT_liquid = ThisModelPart.GetProcessInfo()[FLUID_TEMPERATURE];	
		double delta_max_temp = TT_liquid - stop_temperature;
		double sum_temp = 0.0;

	    int node_size = ThisModelPart.Nodes().size();	    
	    for (int ii = 0; ii < node_size; ii++)
	       {
            ModelPart::NodesContainerType::iterator it_nd = ThisModelPart.NodesBegin() + ii;
		    double temp = it_nd->FastGetSolutionStepValue(TEMPERATURE);
			sum_temp += temp;
			//if( temp > stop_temperature)
			//	return 0.0;
		   }
		sum_temp /= double(node_size);
		sum_temp -= stop_temperature;

		double cooled_percent = 100.0*(1.0 - sum_temp/delta_max_temp);

		return cooled_percent;
	    KRATOS_CATCH("")	
	  }	

	//**********************************************************************************************
	//**********************************************************************************************
	//
	 double ComputeSurfaceWaveDt(ModelPart& ThisModelPart, const double total_volume, const double edge_size,const double max_Dt)
	  {			
	    KRATOS_TRY
        const double cutted_area = ThisModelPart.GetProcessInfo()[CUTTED_AREA];
		const double wet_volume = ThisModelPart.GetProcessInfo()[WET_VOLUME];
		double compare_percente = 0.01;

		if( wet_volume < compare_percente * total_volume || cutted_area <= 0.0)
			return max_Dt;

		double empty_volume = total_volume - wet_volume;

		double dist_1 = empty_volume/cutted_area;
		double dist_2 = cutted_area/dist_1;
		double dist_3 = sqrt(cutted_area);
		
		double max_dist = dist_1;
		max_dist = (dist_1 > dist_2) ? dist_1 : dist_2;
		max_dist = (max_dist > dist_3) ? max_dist : dist_3;
		
		//max_dist = 5.0*dist_3;

		double inv_sqrt_gravity = 0.319275428;
		double wave_dt =inv_sqrt_gravity * edge_size/sqrt(max_dist);

		return ((wave_dt>max_Dt) ? max_Dt : wave_dt);
	    KRATOS_CATCH("")	
	  }	

	private:

	 int CheckMaxTemperature(ModelPart& ThisModelPart)
	 {
	    double last_temp = ThisModelPart.GetTable(3).Data().back().first;
		double is_hot_point = 1.0;
		int node_size = ThisModelPart.Nodes().size();

		for (int ii = 0; ii < node_size; ii++)
		    {
		      ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + ii;

		      double current_temp = it->FastGetSolutionStepValue(TEMPERATURE);
	      
		      if( current_temp < last_temp){
			is_hot_point = 0.0;
			break;
		      }		      
		    }
	       
		ThisModelPart.GetCommunicator().MinAll(is_hot_point);

		return (is_hot_point==1.0)? 1 : 0;


	 }

	};

}  // namespace Kratos.

#endif // ASSIGN_NO_SLIP_CONDITION  defined 


