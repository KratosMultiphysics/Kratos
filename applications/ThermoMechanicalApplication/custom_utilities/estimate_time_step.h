//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//
//  Main authors:    Kazem Kamran
//                   Jordi Rubio
//


#if !defined(KRATOS_ESTIMATE_TIME_STEP )
#define  KRATOS_ESTIMATE_TIME_STEP



// System includes
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath> // Added by Jordi Rubio
// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
//#include "geometries/tetrahedra_3d_4.h"
#include "geometries/point.h"
#include "thermo_mechanical_application.h"
#include "includes/c2c_variables.h"

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

	KRATOS_CLASS_POINTER_DEFINITION(EstimateTimeStep<TDim>);

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
		  BoundedMatrix<double, NumNodes, TDim> DN_DX;

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
		    if(ele_dist <= dist_max)
			{
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
			  if(norm_grad_dist > 0.0)
			  {
				normal_speed /= norm_grad_dist;
				noalias(ElementVel) = (normal_speed/norm_grad_dist)*grad_dist;
			  }
			  else	noalias(ElementVel) = ZeroVector(3);


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
	  //*******************************************************************
	  //************  COMPUTE SOLIDIFICATION DT         *******************
	  //*******************************************************************
	  double ComputeSolidificationCoolingDt(ModelPart& ThisModelPart,
		  const double solidification_percent,
		  const double max_cooling_delta_temp,
		  const double change_in_shrinkage,
		  const double limit_of_mushy_zone,
		  const bool improve_solidification_tracking,
		  const double dt_min,
		  const double dt_max)
	  {
		  KRATOS_TRY

			  //	   const int NumNodes = TDim +1;
			  // 	   is_cold = 0;

		  const double current_dt = ThisModelPart.GetProcessInfo()[DELTA_TIME];
		  const int is_solidified = ThisModelPart.GetProcessInfo()[IS_SOLIDIFIED];
		  //const bool improve_solidification_tracking = true;
		  const double low_boundary_over_mushy_zone = 5.0;
		  int global_is_solidified = is_solidified;
		  // 	   ThisModelPart.GetCommunicator().MinAll(global_is_solidified);

		  //**********************************
		  //**** NOT ALL NODES ARE LIQUID ****
		  //**********************************
		  if(global_is_solidified == 0)
		  {
				//**********************************************************************
				//**** CHECK  IF ALL NODES ARE OVER LIQUIDUS TEMPERATURE            ****
				//**********************************************************************

			  //pre solidification Dt
			  int is_hot = CheckMaxTemperature(ThisModelPart);
			  if( is_hot == 1 )
			  {
				  	//**************************************
					//**** WHEN EVERYTHING IS LIQUID    ****
					//**************************************
				  // if so, then use a maximum cooling objective of 10.0
				  double max_presolodification_delta_tem = std::min(10.0,max_cooling_delta_temp);
				  int node_size = ThisModelPart.Nodes().size();
				  double max_delta_temp = 0.0;
				  std::vector<double> mdelta(OpenMPUtils::GetNumThreads(),0.0);
#pragma omp parallel for shared(mdelta)
				  for (int ii = 0; ii < node_size; ii++)
				  {
					  ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + ii;

					  double current_temp = it->FastGetSolutionStepValue(TEMPERATURE);
					  double old_temp = it->FastGetSolutionStepValue(TEMPERATURE,1);
					  // Get the Maximum on each thread
					  //max_delta_temp=std::max(-current_temp + old_temp,max_delta_temp);
					  double& md = mdelta[OpenMPUtils::ThisThread()];
					  md= std::max(-current_temp + old_temp, md);
				  }
				  //workaround because VS does not support omp 4.0
				  for (int i = 0; i < OpenMPUtils::GetNumThreads(); i++)
				  {
					  max_delta_temp = std::max(max_delta_temp, mdelta[i]);
				  }
				  ThisModelPart.GetCommunicator().MaxAll(max_delta_temp);
				  // Now we can keep on
				  if( max_delta_temp > 0.0 )
				  {
					  double new_delta_time = std::min(1.5, max_presolodification_delta_tem / max_delta_temp); //
					  new_delta_time *= current_dt;
					  // new_delta_time = 1.5*current_dt; // Previous transformationn

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
				//**************************************
				//**** WHEN NOT EVERYTHING IS LIQUID    ****
				//**************************************
				  double current_solidified_volume = 0.0;
				  double old_solidified_volume = 0.0;
				  double current_over_mushy_zone=0.0 ;
				  double old_over_mushy_zone=0.0;
				  double tot_vol = 0.0;

				  std::vector<double> mdelta(OpenMPUtils::GetNumThreads(),0.0);
				  //double max_delta_temp=0.0;
				  int node_size = ThisModelPart.Nodes().size();
#pragma omp parallel for reduction(+:current_solidified_volume,old_solidified_volume, current_over_mushy_zone, old_over_mushy_zone,tot_vol)
				  for (int ii = 0; ii < node_size; ii++)
				  {
					  // Now we look for the Solidifcation Volume
					  ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + ii;
					  double vol = it->GetValue(NODAL_VOLUME);
					  double& md= mdelta[OpenMPUtils::ThisThread()];
					  if (vol > md) { md = vol; }
					  double current_S = it->FastGetSolutionStepValue(SOLIDFRACTION);
					  double old_S = it->FastGetSolutionStepValue(SOLIDFRACTION,1);
					  if(current_S>=limit_of_mushy_zone) {	current_over_mushy_zone+= vol;  }
					  if(old_S>=limit_of_mushy_zone) {	old_over_mushy_zone+= vol;  }
					  current_solidified_volume += vol*current_S;
					  old_solidified_volume += vol*old_S;
					  tot_vol += vol;

					  //filling solidifiacation time
					  //double is_visited = it->FastGetSolutionStepValue(IS_VISITED);
					  //if(is_visited == 0.0 && current_S == 1.0){
						 // it->FastGetSolutionStepValue(IS_VISITED) = 1.0;

						  //double solid_time =  ThisModelPart.GetProcessInfo()[TIME];
						  //double modulus = ThisModelPart.GetProcessInfo()[K0];

						  //it->FastGetSolutionStepValue(SOLIDIF_TIME) = solid_time;
						  //it->FastGetSolutionStepValue(SOLIDIF_MODULUS) = modulus*sqrt(solid_time);
					  //}
					  // Now for the maximum change in temperature

					  //double current_temp = it->FastGetSolutionStepValue(TEMPERATURE);
					  //double old_temp = it->FastGetSolutionStepValue(TEMPERATURE,1);
					//  max_delta_temp=std::max(-current_temp + old_temp,max_delta_temp);

				  }
				  //workaround because VS does not support omp 4.0
				  double max_nodal_volume;
				  for (int i = 0; i < OpenMPUtils::GetNumThreads(); i++)
				  {
					  max_nodal_volume = std::max(max_nodal_volume, mdelta[i]);
				  }

				  ThisModelPart.GetCommunicator().SumAll(current_solidified_volume);
				  ThisModelPart.GetCommunicator().SumAll(old_solidified_volume);
				  ThisModelPart.GetCommunicator().SumAll(tot_vol);
				  ThisModelPart.GetCommunicator().SumAll(current_over_mushy_zone);
				  ThisModelPart.GetCommunicator().SumAll(old_over_mushy_zone);
				  ThisModelPart.GetCommunicator().MaxAll(max_nodal_volume);

				  if(tot_vol == 0.0)   KRATOS_THROW_ERROR(std::logic_error, "inside ComputeSolidificationCoolingDt: total volume is zero!", "")

					if(current_solidified_volume >= tot_vol)
					{
						ThisModelPart.GetProcessInfo()[IS_SOLIDIFIED] = 1;
						return current_dt;
					}
					else
					{

						double delta_solid = current_solidified_volume - old_solidified_volume;
						double delta_over_mushy_zone=current_over_mushy_zone-old_over_mushy_zone;
						if( delta_solid > 0.0 || delta_over_mushy_zone > 0.0 )
						{
							delta_solid /= tot_vol;
							delta_over_mushy_zone /= tot_vol;
							if(delta_solid<=0.0){delta_solid=delta_over_mushy_zone*solidification_percent/change_in_shrinkage;}  //It sets the value so that in next step we get exactly the same value
							double solidification_ratio = solidification_percent / delta_solid;
							double limited_change_in_shrinkage_ratio = 0.0;
							if (delta_over_mushy_zone <= 0.0)
							{
								delta_over_mushy_zone = delta_solid*change_in_shrinkage / solidification_percent; // It sets the value so that in next step we get exactly the same value
								limited_change_in_shrinkage_ratio = change_in_shrinkage;
							}
							else
							{
								if (improve_solidification_tracking == true)
								{
									double lower_limit = (max_nodal_volume/ tot_vol)*low_boundary_over_mushy_zone;
									double liquid_percent = 1.00 - (current_over_mushy_zone / tot_vol);
									limited_change_in_shrinkage_ratio = std::min(change_in_shrinkage, liquid_percent*0.5);//delta_over_mushy_zone*0.8);
									limited_change_in_shrinkage_ratio = std::max(limited_change_in_shrinkage_ratio, lower_limit);
									KRATOS_WATCH(lower_limit)
									KRATOS_WATCH(change_in_shrinkage)
									KRATOS_WATCH(limited_change_in_shrinkage_ratio)

								}
								else
								{
									limited_change_in_shrinkage_ratio = change_in_shrinkage;
								}
							}
							double change_in_shrinkage_ratio= limited_change_in_shrinkage_ratio / delta_over_mushy_zone;
							double K = std::min(solidification_ratio, change_in_shrinkage_ratio);
							double new_dt = std::min(1.5, K ) * current_dt;
							if( new_dt > dt_max) new_dt = dt_max;
							else if( new_dt < dt_min) new_dt = dt_min;
							return new_dt;
						}
						else
						{
						return current_dt;
						}
					}
			  }
		  }

		  else //coling delta_t
		  {
			  //double cooling_dt_max = dt_max;//30.0*dt_max;
			  int node_size = ThisModelPart.Nodes().size();
			  double max_delta_temp = 0.0;
			  for (int ii = 0; ii < node_size; ii++)
			  {
				  ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + ii;

				  double current_temp = it->FastGetSolutionStepValue(TEMPERATURE);
				  double old_temp = it->FastGetSolutionStepValue(TEMPERATURE,1);
				  max_delta_temp=std::max(-current_temp+old_temp,max_delta_temp);
			  }

			  ThisModelPart.GetCommunicator().MaxAll(max_delta_temp);

			  if( max_delta_temp > 0.0 ){
				  double new_delta_time = max_cooling_delta_temp / max_delta_temp;
				  new_delta_time *= current_dt;

				  if( new_delta_time > dt_max)
					  new_delta_time = dt_max;
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
	  //************************************************************************************************
	  //*************** FIND AN ESTIMATION FOR SOLIDIFICATION TIME. NO VIRTUAL MOLD *******************
	  //************************************************************************************************
	   double EstimateSolidificationTimeNoVirtualMould(ModelPart& ThisModelPart)
	  {
	    KRATOS_TRY

	    double solidification_time = 0.0;
	    const double density = ThisModelPart.GetProcessInfo()[DENSITY];
	    const double cc= ThisModelPart.GetProcessInfo()[SPECIFIC_HEAT];
	    const double htc= ThisModelPart.GetProcessInfo()[HTC];
		const double mould_temperature=ThisModelPart.GetProcessInfo()[MOLD_AVERAGE_TEMPERATURE];

	    const double TT_solid = ThisModelPart.GetProcessInfo()[SOLID_TEMPERATURE];
	    const double TT_liquid = ThisModelPart.GetProcessInfo()[FLUID_TEMPERATURE];

	    const double LL = ThisModelPart.GetProcessInfo()[LATENT_HEAT];
	    double tot_vol = 0.0;
        double tot_area = 0.0;
	    int node_size = ThisModelPart.Nodes().size();
	    for (int ii = 0; ii < node_size; ii++)
		// Compute Part Volume and Area
	    {
			ModelPart::NodesContainerType::iterator it_nd = ThisModelPart.NodesBegin() + ii;
			double vol = it_nd->GetValue(NODAL_VOLUME);
			tot_vol += vol;
			double area =  it_nd->FastGetSolutionStepValue(NODAL_PAUX);
			tot_area += area;
	    }

	    // Check Area is not 0
		if ( tot_area == 0.0 || tot_vol == 0.0)   KRATOS_THROW_ERROR(std::invalid_argument,"AREA or VOLUME is Zero", "");
	    // Formula for stimating solidification time
	    solidification_time = density * ( cc * ( TT_liquid - TT_solid) + LL) / (htc * 0.5*(TT_solid-mould_temperature));
		solidification_time *= pow(tot_vol/tot_area , 0.8);
		return solidification_time;
		KRATOS_CATCH("")
	   }
	  //************************************************************************************************
	  //**************** FIND AN ESTIMATION FOR SOLIDIFICATION TIME. VIRTUAL MOLD *********************
	  //************************************************************************************************
		/* For solving this we are going to suppose that we dissipate all the energy through the mould outer surface.
		 We Estimate the inner Energy of the System as the SUM of 3 contributions. The energy loss needed to cool the
		 mart, the energy loss needed to make the part chage its phase and the energy needed to cool the mould. The
		 contribution of the mould only is considered if positive.
		 All terms are linearized with respect to the temperature, so that.
			E_1=V_{part}*C_{part}*\rho_{part}*T
			E_2=V_{mould}*V_{fact}*\rho_{mould}*C_{mould}*max(0,(T_{mould}-T_{end})/(T_{ini}-T_{end}))*T
			E_3=LH*\rho*V_{part}*(T-T_{end})/(T_{ini}-T_{end})
		Now we have that, the q (Energy time derivative is)
			dE/dt=HTC_{env}*Sfact*A_{part}*(T-T_{env})
		We can set it as ODE, and solve analytacally (recall this is a linealization, but will be enough for our purpose)
			dE_1/dT=V_{part}*C_{part}*\rho_{part}
			dE_2/dT=V_{mould}*V_{fact}*\rho_{mould}*C_{mould}*max(0,(T_{mould}-T_{end})/(T_{ini}-T_{end}))
			dE_3/dT=LH*\rho*V_{part}/(T_{ini}-T_{end})
		Solving the ODE,we have that:
			\Delta t= =( (dE_1/dT+dE_2/dT+dE_2/dT)/ HTC_{env}*Sfact*A_{part} )*ln( ((T_{ini}- T_{env})/(T_{end}- T_{env}) )
		As starting temperature we suppose the Average temperature, and as initial mould temperature, we suppose the initial mouls temperature
		 */
	  double EstimateSolidificationTime(ModelPart& ThisModelPart)
	  {
	    KRATOS_TRY

	    double solidification_time = 0.0;

		 // Auxiliaty variables
		double dE1=0.0;
		double dE2=0.0;
		double dE3=0.0;
		double DENOM=0.0;

		// Environment and part variables
		const double ambient_temperature=ThisModelPart.GetProcessInfo()[AMBIENT_TEMPERATURE];
	    const double LL = ThisModelPart.GetProcessInfo()[LATENT_HEAT];
		const double density = ThisModelPart.GetProcessInfo()[DENSITY];
	    const double cc= ThisModelPart.GetProcessInfo()[SPECIFIC_HEAT];
		const double initial_temperature= ThisModelPart.GetProcessInfo()[AVERAGE_TEMPERATURE];
		const double stop_temperature= ThisModelPart.GetProcessInfo()[SOLID_TEMPERATURE];

		// Loop Over the nodes - Compute E1 term
	    double tot_vol = 0.0;
	    int node_size = ThisModelPart.Nodes().size();
	    for (int ii = 0; ii < node_size; ii++)
	    {
			ModelPart::NodesContainerType::iterator it_nd = ThisModelPart.NodesBegin() + ii;
			double vol = it_nd->GetValue(NODAL_VOLUME);
			tot_vol += vol;
			vol=pow(vol,1.0);
			// dE1 - First Term
			//dE_1/dT=V_{part}*C_{part}*\rho_{part}
			dE1+= vol*density*cc;
			// dE3 - Third Term
			//dE_3/dT=LH*\rho*V_{part}/(T_{ini}-T_{end})
			dE3+=LL*density*vol/(initial_temperature-stop_temperature);
	    }
		double tot_area=0.0;
		double avg_conductivity=0.0;
		double avg_density=0.0;
		double avg_sheat=0.0;
		 // Loop over the conditions Compute E2, E3 and Denom term
		 for (ModelPart::ConditionIterator itCond = ThisModelPart.ConditionsBegin(); itCond != ThisModelPart.ConditionsEnd(); itCond++ )
        {
			// Generate the Geometry of the condition
			Condition::GeometryType& rGeom = itCond->GetGeometry();
			const double mould_density= itCond->GetProperties()[MOLD_DENSITY];
			const double mould_specific_heat= itCond->GetProperties()[MOLD_SPECIFIC_HEAT];
			const double mould_thickness = itCond->GetProperties()[MOLD_THICKNESS];
			const double mould_vfact= itCond->GetProperties()[MOLD_VFACT];
			const double mould_sfact= itCond->GetProperties()[MOLD_SFACT];
			const double mould_htc_env= itCond->GetProperties()[MOLD_HTC_ENVIRONMENT];
			const double mould_cond=itCond->GetProperties()[MOLD_CONDUCTIVITY];
			//const double mould_conductivity = itCond->GetProperties()[MOLD_CONDUCTIVITY];
			const double mould_temperature = itCond->GetProperties()[MOLD_TEMPERATURE];
			double tarea=rGeom.DomainSize();
			const double condition_area=pow(tarea,1.0);
			tot_area+=tarea;
			// dE2 - Second Term
			//dE_2/dT=V_{mould}*V_{fact}*\rho_{mould}*C_{mould}*max(0,(T_{mould}-T_{end})/(T_{ini}-T_{end}))
			double aux =condition_area*mould_thickness*mould_vfact*mould_density*mould_specific_heat;
			double aux2=(mould_temperature- stop_temperature)/(initial_temperature-stop_temperature ) ;
			dE2+=std::max(aux*aux2,0.0);
			// Denom.
			// HTC_{env}*Sfact*A_{part} )
			DENOM+= mould_htc_env*condition_area*mould_sfact;
			// To be used by Chorinov Formulation
			avg_conductivity+=mould_cond*tarea;
			avg_density+=mould_density*tarea;
			avg_sheat+=mould_specific_heat*tarea;

		 }

		 solidification_time = ((dE1+dE2+dE3)/DENOM)*log( (initial_temperature-ambient_temperature)/(stop_temperature-ambient_temperature) );

		// Now we will compute Chorinov's Rule
		avg_conductivity/=tot_area;
		avg_density/=tot_area;
		avg_sheat/=tot_area;

		double solidification_time_chorinov=0.0;
		//const double htc= ThisModelPart.GetProcessInfo()[HTC];
		//const double mould_temperature=ThisModelPart.GetProcessInfo()[MOLD_AVERAGE_TEMPERATURE];
		solidification_time_chorinov=pow(density*LL/fabs(initial_temperature-stop_temperature),2)*(3.1416/(4*avg_conductivity*avg_density*avg_sheat));
		solidification_time_chorinov*=1+(cc*pow((initial_temperature-stop_temperature)/LL,2));
		solidification_time_chorinov*=pow(tot_vol/tot_area,1.5);
		return std::max(solidification_time,solidification_time_chorinov);
	    KRATOS_CATCH("")
	  }

	 double CheckStopTemperature(ModelPart& ThisModelPart,const double stop_temperature)
	  {
	    KRATOS_TRY
		const double avg_temp = ThisModelPart.GetProcessInfo()[AVERAGE_TEMPERATURE];
		double delta_max_temp = avg_temp - stop_temperature;
		double sum_temp = 0.0;

	    int node_size = ThisModelPart.Nodes().size();
#pragma omp parallel for reduction(+:sum_temp)
	    for (int ii = 0; ii < node_size; ii++)
	       {
            ModelPart::NodesContainerType::iterator it_nd = ThisModelPart.NodesBegin() + ii;
		    double temp = it_nd->FastGetSolutionStepValue(TEMPERATURE);
			sum_temp += std::max(temp,0.99999*stop_temperature); // before temp

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

	 int CheckIsInTransition(ModelPart& ThisModelPart)
	 {
		double const sol_temp = ThisModelPart.GetProcessInfo()[SOLID_TEMPERATURE] ;
		double const liq_temp = ThisModelPart.GetProcessInfo()[FLUID_TEMPERATURE] ;
		int is_in_range_point = 0;
		int node_size = ThisModelPart.Nodes().size();

		for (int ii = 0; ii < node_size; ii++)
		    {
		      ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + ii;

		      double current_temp = it->FastGetSolutionStepValue(TEMPERATURE);

			  if( current_temp >sol_temp && current_temp<liq_temp){
				is_in_range_point = 1;
				break;
		      }
		    }

		ThisModelPart.GetCommunicator().MaxAll(is_in_range_point);

		return (is_in_range_point==1)? 1 : 0;
	 }

	//**********************************************************************************************
	//**********************************************************************************************
	//
	 double EstimateCoolingTime(ModelPart& ThisModelPart,const double stop_temperature)
	 {
		 /* For solving this we are going to suppose that we dissipate all the energy through the mould outer surface.
		 We Estimate the inner Energy of the System as the SUM of 3 contributions. The energy loss needed to cool the
		 mart, the energy loss needed to make the part chage its phase and the energy needed to cool the mould. The
		 contribution of the mould only is considered if positive.
		 All terms are linearized with respect to the temperature, so that.
			E_1=V_{part}*C_{part}*\rho_{part}*T
			E_2=V_{mould}*V_{fact}*\rho_{mould}*C_{mould}*max(0,(T_{mould}-T_{end})/(T_{ini}-T_{end}))*T
			E_3=LH*\rho*V_{part}*(T-T_{end})/(T_{ini}-T_{end})
		Now we have that, the q (Energy time derivative is)
			dE/dt=HTC_{env}*Sfact*A_{part}*(T-T_{env})
		We can set it as ODE, and solve analytacally (recall this is a linealization, but will be enough for our purpose)
			dE_1/dT=V_{part}*C_{part}*\rho_{part}
			dE_2/dT=V_{mould}*V_{fact}*\rho_{mould}*C_{mould}*max(0,(T_{mould}-T_{end})/(T_{ini}-T_{end}))
			dE_3/dT=LH*\rho*V_{part}/(T_{ini}-T_{end})
		Solving the ODE,we have that:
			\Delta t= =( (dE_1/dT+dE_2/dT+dE_2/dT)/ HTC_{env}*Sfact*A_{part} )*ln( ((T_{ini}- T_{env})/(T_{end}- T_{env}) )
		As starting temperature we suppose the Average temperature, and as initial mould temperature, we suppose the initial mouls temperature
		 */
		 // Auxiliaty variables
		double CoolingTime;
		double dE1=0.0;
		double dE2=0.0;
		double dE3=0.0;
		double DENOM=0.0;

		// Environment and part variables
		const double ambient_temperature=ThisModelPart.GetProcessInfo()[AMBIENT_TEMPERATURE];
	    const double LL = ThisModelPart.GetProcessInfo()[LATENT_HEAT];
		const double density = ThisModelPart.GetProcessInfo()[DENSITY];
	    const double cc= ThisModelPart.GetProcessInfo()[SPECIFIC_HEAT];
		const double initial_temperature= ThisModelPart.GetProcessInfo()[AVERAGE_TEMPERATURE];

		// Loop Over the nodes - Compute E1 term
	    double tot_vol = 0.0;
	    int node_size = ThisModelPart.Nodes().size();
	    for (int ii = 0; ii < node_size; ii++)
	    {
			ModelPart::NodesContainerType::iterator it_nd = ThisModelPart.NodesBegin() + ii;
			double vol = it_nd->GetValue(NODAL_VOLUME);
			tot_vol += vol;
			vol=pow(vol,0.8);
			// dE1 - First Term
			//dE_1/dT=V_{part}*C_{part}*\rho_{part}
			dE1+= vol*density*cc;
			// dE3 - Third Term
			//dE_3/dT=LH*\rho*V_{part}/(T_{ini}-T_{end})
			dE3+=LL*density*vol/(initial_temperature-stop_temperature);
	    }
		double tot_area=0.0;
		double avg_conductivity=0.0;
		double avg_density=0.0;
		double avg_sheat=0.0;
		double avg_env_htc=0.0;
		 // Loop over the conditions Compute E2, E3 and Denom term
		 for (ModelPart::ConditionIterator itCond = ThisModelPart.ConditionsBegin(); itCond != ThisModelPart.ConditionsEnd(); itCond++ )
        {
			// Generate the Geometry of the condition
			Condition::GeometryType& rGeom = itCond->GetGeometry();
			const double mould_density= itCond->GetProperties()[MOLD_DENSITY];
			const double mould_specific_heat= itCond->GetProperties()[MOLD_SPECIFIC_HEAT];
			const double mould_thickness = itCond->GetProperties()[MOLD_THICKNESS];
			const double mould_vfact= itCond->GetProperties()[MOLD_VFACT];
			const double mould_sfact= itCond->GetProperties()[MOLD_SFACT];
			const double mould_htc_env= itCond->GetProperties()[MOLD_HTC_ENVIRONMENT];
			const double mould_conductivity = itCond->GetProperties()[MOLD_CONDUCTIVITY];
			const double mould_temperature = itCond->GetProperties()[MOLD_TEMPERATURE];

			double tarea=rGeom.DomainSize();
			tot_area+=tarea;

			const double condition_area=pow(fabs(rGeom.DomainSize()),1.0);
			// dE2 - Second Term
			//dE_2/dT=V_{mould}*V_{fact}*\rho_{mould}*C_{mould}*max(0,(T_{mould}-T_{end})/(T_{ini}-T_{end}))
			double aux =condition_area*mould_thickness*mould_vfact*mould_density*mould_specific_heat;
			double aux2=(mould_temperature- stop_temperature)/(initial_temperature-stop_temperature ) ;
			dE2+=std::max(aux*aux2,0.0);
			// Denom.
			// HTC_{env}*Sfact*A_{part} )
			DENOM+= mould_htc_env*condition_area*mould_sfact;

			// To be used by Chorinov Formulation
			avg_conductivity+=mould_conductivity*tarea;
			avg_density+=mould_density*tarea;
			avg_sheat+=mould_specific_heat*tarea;
			avg_env_htc+=tarea*mould_htc_env*mould_sfact;
		 }

		 CoolingTime = ((dE1+dE2+dE3)/DENOM)*log( (initial_temperature-ambient_temperature)/(stop_temperature-ambient_temperature) );

		// Now we will compute Chorinov's Rule
		avg_conductivity/=tot_area;
		avg_density/=tot_area;
		avg_sheat/=tot_area;
		avg_env_htc/=tot_area;
		double cooling_time_chorinov=0.0;
		//const double htc= ThisModelPart.GetProcessInfo()[HTC];
		const double solid_temp= ThisModelPart.GetProcessInfo()[SOLID_TEMPERATURE];
		//const double mould_temperature=ThisModelPart.GetProcessInfo()[MOLD_AVERAGE_TEMPERATURE];
		cooling_time_chorinov=pow(density*LL/fabs(initial_temperature-solid_temp),2)*(3.1416/(4*avg_conductivity*avg_density*avg_sheat));
		cooling_time_chorinov*=1+(cc*pow((initial_temperature-solid_temp)/LL,2));
		cooling_time_chorinov*=pow(tot_vol/tot_area,2);
		// Now we compute the time from solidification to cooling
		double time_to_cool=(tot_vol*cc*density*fabs(solid_temp-stop_temperature))/(avg_env_htc*tot_area*(0.5*initial_temperature+0.5*stop_temperature-ambient_temperature));
		cooling_time_chorinov+=time_to_cool;

		return std::max(CoolingTime,cooling_time_chorinov);
	 }

	 /////////////////////////////////////////////////////////////////////////
	 int CheckMinTemperature(ModelPart& ThisModelPart)
	 {
	    double last_temp = 1e9;
		//double is_hot_point = 1.0;
		int node_size = ThisModelPart.Nodes().size();

		for (int ii = 0; ii < node_size; ii++)
		{
			ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + ii;
			double current_temp = it->FastGetSolutionStepValue(TEMPERATURE);
			if( current_temp < last_temp){last_temp = current_temp;   }
		}
		ThisModelPart.GetCommunicator().MinAll(last_temp);
		return last_temp;
	 }
	 //private:
	//	 /////////////////////////////////////////////////////////////////////////////
	 int CheckMaxTemperature(ModelPart& ThisModelPart)
	 {
		 double last_temp = ThisModelPart.GetProcessInfo()[FLUID_TEMPERATURE]; //GetTable(3).Data().back().first;
		 double is_hot_point = 1.0;
		 int node_size = ThisModelPart.Nodes().size();
		 //KRATOS_WATCH(omp_get_max_threads())
		 std::vector<double> local_is_hot_point(OpenMPUtils::GetNumThreads(),1.0);
// #pragma omp parallel for shared(local_is_hot_point)
		 for (int ii = 0; ii < node_size; ii++)
		 {
			 ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + ii;

			 double current_temp = it->FastGetSolutionStepValue(TEMPERATURE);

			 if (current_temp < last_temp) {
				 //is_hot_point = 0.0;
				 local_is_hot_point[OpenMPUtils::ThisThread()]= 0.0;
				 break;
			 }
		 }
		 //Now we finf the minimum is_hot_point among threads
		 for (int ii = 0; ii < OpenMPUtils::GetNumThreads(); ii++)
		 {
			 is_hot_point = std::min(is_hot_point, local_is_hot_point[ii]);
		 }


		 ThisModelPart.GetCommunicator().MinAll(is_hot_point);

		 return (is_hot_point == 1.0) ? 1 : 0;
	 }

	};

}  // namespace Kratos.

#endif // ASSIGN_NO_SLIP_CONDITION  defined


