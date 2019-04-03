//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Kazem Kamran
//                   Riccardo Rossi
//                   Pooyan Dadvand
//                   Jordi Rubio
//


#if !defined(KRATOS_BIPHASIC_FILLING_UTILITIES_INCLUDED )
#define  KRATOS_BIPHASIC_FILLING_UTILITIES_INCLUDED



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
#include "utilities/enrichment_utilities.h"
#include "utilities/timer.h"
//#include "geometries/tetrahedra_3d_4.h"
#include "geometries/point.h"
#include "thermo_mechanical_application.h"
// #include "custom_conditions/environment_contact.h"
//#include "includes/variables.h"
#include "../incompressible_fluid_application/custom_utilities/parallel_extrapolation_utilities.h"
#include "includes/kratos_flags.h"
#include "includes/c2c_variables.h"
#include "includes/cfd_variables.h"


namespace Kratos
{

class BiphasicFillingUtilities
{
public:

	KRATOS_CLASS_POINTER_DEFINITION(BiphasicFillingUtilities);

    //**********************************************************************************************
    //**********************************************************************************************
    double CreateAutoExitAssignAirSmagorinsky(ModelPart& ThisModelPart, double y_wall, double C_Smagorinsky)
    {
        KRATOS_TRY;
        /*			AirSmagorinskey(ThisModelPart, C_Smagorinsky);*/
        int node_size = ThisModelPart.Nodes().size();
        double is_exit = 0.0;
        for (int ii = 0; ii < node_size; ii++)
        {
            ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + ii;
            double str_flag = it->GetValue(IS_STRUCTURE);
            double slip_flag = it->GetSolutionStepValue(IS_SLIP);

            if (str_flag == 0.0 && slip_flag>=10.0)
            {
                is_exit = 1.0;
// 					 return 1.0;
            }
        }
        //syncronoze
        ThisModelPart.GetCommunicator().MaxAll(is_exit);
        if(is_exit == 1.0)
            return 1.0;


        // if there is no dry node
        double is_dry_node = 0.0;
        #pragma omp parallel for firstprivate(node_size)
        for (int ii = 0; ii < node_size; ii++)
        {
            ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + ii;
            double dist = it->FastGetSolutionStepValue(DISTANCE);
            double slip_flag = it->GetSolutionStepValue(IS_SLIP);
            double str_flag = it->GetValue(IS_STRUCTURE);

            if(dist > 0.0)
            {
                is_dry_node = 1.0;
                //slip_flag=10.0 refers to the boundary nodes with well-defined normal
                if(slip_flag==10.0 && str_flag!=0.0)
                {
                    it->SetValue(IS_STRUCTURE,0.0);
                    it->SetValue(Y_WALL,y_wall);
                }
            }
        }
        //syncronoze
        ThisModelPart.GetCommunicator().MaxAll(is_dry_node);

        //assign smagorinsky at air element
        AirSmagorinskey(ThisModelPart, C_Smagorinsky);


        return is_dry_node;

        KRATOS_CATCH("")
    }
    //**********************************************************************************************
    /*for node in fluid_model_part.Nodes:
    	slip_flag = node.GetSolutionStepValue(IS_SLIP)
    	nd_dist = node.GetSolutionStepValue(DISTANCE)
            if((slip_flag == 20.0 or slip_flag == 30.0 )):#
    	  if(nd_dist< 0.0):
    	    node.SetValue(IS_STRUCTURE,1.0)
    	    node.SetValue(Y_WALL,y_wall_val)
    	  else:
    	    node.SetValue(IS_STRUCTURE,0.0)
    	    node.SetValue(Y_WALL,y_wall_val*y_wall_fac)

    	if(slip_flag == 10.0):
    	  if(nd_dist< 0.0):
    	    node.SetValue(Y_WALL,y_wall_val)
    	    node.SetValue(IS_STRUCTURE,1.0)
    	  else:
    	    node.SetValue(Y_WALL,y_wall_val*y_wall_fac)
    	    node.SetValue(IS_STRUCTURE,0.0)*/
    //**********************************************************************************************
    double AssignSmoothBoundaryAirExit(ModelPart& ThisModelPart, bool air_exit_flag, const double y_wall_val, const double  y_wall_fac)
    {
        KRATOS_TRY;
// 			int node_size = ThisModelPart.Nodes().size();
        int node_size = ThisModelPart.GetCommunicator().LocalMesh().Nodes().size();
        double is_str = 1.0;
        if(air_exit_flag) is_str = 0.0;

        int wet_nodes = 0;

        #pragma omp parallel for firstprivate(node_size)
        for (int ii = 0; ii < node_size; ii++)
        {

// 			    ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + ii;
            ModelPart::NodesContainerType::iterator it = ThisModelPart.GetCommunicator().LocalMesh().NodesBegin() + ii;

            double dist = it->FastGetSolutionStepValue(DISTANCE);
            double slip_flag = it->GetSolutionStepValue(IS_SLIP);
            if(dist<0.0)
            {
                #pragma omp atomic
                wet_nodes++;
            }

            if(slip_flag == 20.0 || slip_flag == 30.0 )//edges(20) and corners(30) are automatic air exits till they are wetten
                if(dist<0.0)
                {
                    it->SetValue(IS_STRUCTURE,1.0);
                    it->SetValue(Y_WALL,y_wall_val);
                }
                else
                {
                    it->SetValue(IS_STRUCTURE,0.0);
                    it->SetValue(Y_WALL,y_wall_val*y_wall_fac);
                }
            else if(slip_flag == 10.0)//smooth boundaries(10), if dry, can be air exit or not
            {
                if(dist<0.0)
                {
                    it->SetValue(IS_STRUCTURE,1.0);
                    it->SetValue(Y_WALL,y_wall_val);
                }
                else
                {
                    it->SetValue(IS_STRUCTURE,is_str);
                    it->SetValue(Y_WALL,y_wall_val*y_wall_fac);
                }//y_wall_val*y_wall_fac
            }

            //filling time
            //double is_visited = it->FastGetSolutionStepValue(IS_VISITED);
            //if(is_visited == 0.0 && dist<=0.0)
            //{
            //    it->FastGetSolutionStepValue(IS_VISITED) = 1.0;

            //    double filling_time =  ThisModelPart.GetProcessInfo()[TIME];
            //    it->FastGetSolutionStepValue(FILLTIME) = filling_time;
            //}


        }
        //syncronoze
        ThisModelPart.GetCommunicator().SumAll(wet_nodes);
        ThisModelPart.GetCommunicator().SumAll(node_size);
        double filling_percent = 0.0;

        if(wet_nodes != 0) filling_percent = 100.0*double(wet_nodes)/double(node_size);

        return filling_percent;
        KRATOS_CATCH("")
    }

    //**********************************************************************************************
    //**********************************************************************************************
    double ComputeFillPercentage(ModelPart& ThisModelPart, const double corrected_time )
    {
        KRATOS_TRY

        int node_size = ThisModelPart.GetCommunicator().LocalMesh().Nodes().size();

        int wet_nodes = 0;

        #pragma omp parallel for firstprivate(node_size) reduction(+:wet_nodes)
        for (int ii = 0; ii < node_size; ii++)
        {
            ModelPart::NodesContainerType::iterator it = ThisModelPart.GetCommunicator().LocalMesh().NodesBegin() + ii;

            double dist = it->FastGetSolutionStepValue(DISTANCE);

            if(dist<0.0)
                wet_nodes++;


            //filling time
            double is_visited = it->FastGetSolutionStepValue(IS_VISITED);
            if(is_visited == 0.0 && dist<=0.0)
            {
                it->FastGetSolutionStepValue(IS_VISITED) = 1.0;

                //double filling_time =  ThisModelPart.GetProcessInfo()[TIME];
                //it->FastGetSolutionStepValue(FILLTIME) = filling_time*time_correction_factor;
				it->FastGetSolutionStepValue(FILLTIME) =corrected_time;
            }
        }

        //syncronoze
        ThisModelPart.GetCommunicator().SumAll(wet_nodes);
        ThisModelPart.GetCommunicator().SumAll(node_size);
        double filling_percent = 0.0;

        if(wet_nodes != 0) filling_percent = 100.0*double(wet_nodes)/double(node_size);

        return filling_percent;
        KRATOS_CATCH("")
    }

    //**********************************************************************************************
    //**********************************************************************************************
    void ApplyFluidProperties(ModelPart& ThisModelPart, const double water_mu, const double water_density ,const double air_mu,const double air_density)
    {
        KRATOS_TRY;
        int node_size = ThisModelPart.Nodes().size();

        #pragma omp parallel for firstprivate(node_size)
        for (int ii = 0; ii < node_size; ii++)
        {
            ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + ii;
            double dist = it->FastGetSolutionStepValue(DISTANCE);
            if(dist<=0.0)
            {
                it->FastGetSolutionStepValue(DENSITY) = water_density;
                it->FastGetSolutionStepValue(VISCOSITY) = water_mu;
            }
            else
            {
                it->FastGetSolutionStepValue(DENSITY) = air_density;
                it->FastGetSolutionStepValue(VISCOSITY) = air_mu;
            }
        }
        KRATOS_CATCH("")
    }
    //**********************************************************************************************
    //**********************************************************************************************
    void DistanceFarRegionCorrection(ModelPart& ThisModelPart, const double max_dist)
    {
        KRATOS_TRY;
        int node_size = ThisModelPart.Nodes().size();
       // double max_cutted_elem_size = 0.0;
        //max_cutted_elem_size = ComputeCharactristicCuttedLength(ThisModelPart);

        #pragma omp parallel for firstprivate(node_size)
        for (int ii = 0; ii < node_size; ii++)
        {
            ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + ii;
            double& current_dist = it->FastGetSolutionStepValue(DISTANCE);
            const double old_dist = it->FastGetSolutionStepValue(DISTANCE,1);

           // if( fabs(old_dist) >= CFL*max_cutted_elem_size && current_dist*old_dist <= 0.0)
            if( fabs(old_dist) >= max_dist && current_dist*old_dist <= 0.0)
                current_dist = old_dist;
        }
        KRATOS_CATCH("")
    }
	//**********************************************************************************************
	//**********************************************************************************************
	double sign(const double& a)
	{
		if(a < 0) return -1.0;
		else return 1.0;
	}
    //**********************************************************************************************
    //**********************************************************************************************
    void VolumeCorrection(ModelPart& ThisModelPart, const double Net_volume, const double max_correction, const bool CorrectNegativeVolume=false)
    {
        KRATOS_TRY

        double wet_volume = 0.0;
		double wet_volume_old=wet_volume;
        double cutted_area = 0.0;
		double wet_volume_left=0.0;
		double wet_volume_right=0.0;
		double tol=1e-5;
		double tolv=5e-3;
		double lower_correction;
		double upper_correction;
        int node_size = ThisModelPart.GetCommunicator().LocalMesh().Nodes().size();

		// First we compute the Total Volume of the Fluid
		// First we compute the Total Volume of the Fluid
        ComputeWetVolumeAndCuttedArea(ThisModelPart, wet_volume, cutted_area);

		// Now we compute the difference between the Total Volume and the volume that has enetered through the inlet
        double volume_difference = (fabs(Net_volume) - wet_volume);
		// First guess in correction
		double correction = volume_difference/cutted_area;
		double correction_old=correction;
		//double signcorrection=sign(correction);
		//Maximum signed correction
		double maximum_signed_correction=fabs(max_correction)*(sign(correction)); //Way of obtaining sign(x)
		// If correction is greater than signed correction, we keep the maximum correction.If the wet volume exceeds the correction, then we start from maximum_signed_correction
		bool exit_loop=false; // Just to skip the iterations
		if((correction>maximum_signed_correction)&&(correction>0)){
			correction=maximum_signed_correction;
			ComputeVolumeAndCuttedAreaInDistance(ThisModelPart,wet_volume,cutted_area,correction);
			lower_correction=0.0;
			ComputeVolumeAndCuttedAreaInDistance(ThisModelPart,wet_volume_left,cutted_area,lower_correction);
			upper_correction=correction;
			wet_volume_right=wet_volume;
			if(fabs(Net_volume)>wet_volume_right){exit_loop=true;}
			}
		if((correction<maximum_signed_correction)&&(correction<0)){
			correction=maximum_signed_correction;
			ComputeVolumeAndCuttedAreaInDistance(ThisModelPart,wet_volume,cutted_area,correction);
			upper_correction=0.0;
			ComputeVolumeAndCuttedAreaInDistance(ThisModelPart,wet_volume_right,cutted_area,upper_correction);
			lower_correction=correction;
			wet_volume_left=wet_volume;
			if(fabs(Net_volume)<wet_volume_left){exit_loop=true;}
			}
		// Now we find the left and right limits
		if(exit_loop==false)
		{
			double extreme_correction;
			ComputeVolumeAndCuttedAreaInDistance(ThisModelPart,wet_volume,cutted_area,correction);
			if(correction>0)
			{
				extreme_correction=fabs(max_correction);
				if(wet_volume<=fabs(Net_volume)){
					lower_correction=correction;
					wet_volume_left=wet_volume;
					upper_correction=extreme_correction;
					ComputeVolumeAndCuttedAreaInDistance(ThisModelPart,wet_volume_right,cutted_area,extreme_correction);
				}
				else
				{
					lower_correction=0.0;
					ComputeVolumeAndCuttedAreaInDistance(ThisModelPart,wet_volume_left,cutted_area,0.0);
					wet_volume_right=wet_volume;
					upper_correction=correction;
				}
			}
			else
			{
				extreme_correction=-fabs(max_correction);
				if(wet_volume>=fabs(Net_volume)){
					upper_correction=correction;
					wet_volume_right=wet_volume;
					lower_correction=extreme_correction;
					ComputeVolumeAndCuttedAreaInDistance(ThisModelPart,wet_volume_left,cutted_area,extreme_correction);
				}
				else
				{
					upper_correction=0.0;
					ComputeVolumeAndCuttedAreaInDistance(ThisModelPart,wet_volume_right,cutted_area,0.0);
					wet_volume_left=wet_volume;
					lower_correction=correction;
				}
			}
			ComputeVolumeAndCuttedAreaInDistance(ThisModelPart,wet_volume,cutted_area,correction);
			ComputeVolumeAndCuttedAreaInDistance(ThisModelPart,wet_volume_right,cutted_area,upper_correction);
			ComputeVolumeAndCuttedAreaInDistance(ThisModelPart,wet_volume_left,cutted_area,lower_correction);
		}

		// Now we loop until convergence
		unsigned int iteration=0;
		//double inc_correction=1000.0;
		while((iteration<10)&&(exit_loop==false))
			{
			correction_old=correction;
			wet_volume_old=wet_volume;
			double aux_vol_r=wet_volume_right-fabs(Net_volume);
			double aux_vol_l=wet_volume_left-fabs(Net_volume);
			correction=(aux_vol_r*lower_correction-aux_vol_l*upper_correction)/(wet_volume_right-wet_volume_left);
			if((correction<lower_correction)||(correction>upper_correction)){KRATOS_WATCH("ERROR CORRECTING VOLUME IN VOLUME_CORRECTION");}
			ComputeVolumeAndCuttedAreaInDistance(ThisModelPart,wet_volume,cutted_area,correction);
			volume_difference = fabs(Net_volume) - wet_volume;
			if(fabs(Net_volume)>wet_volume){
				lower_correction=correction;
				wet_volume_left=wet_volume;
				}
			else{
				upper_correction=correction;
				wet_volume_right=wet_volume;
				}
			//Now the middle point just in case
			double middle_point=(lower_correction+upper_correction)/2.0;
			double wet_volume_middle=0.0;
			ComputeVolumeAndCuttedAreaInDistance(ThisModelPart,wet_volume_middle,cutted_area,middle_point);
			if(wet_volume_middle<fabs(Net_volume)){
				lower_correction=middle_point;
				wet_volume_left=wet_volume_middle;
			}
			else{
				upper_correction=middle_point;
				wet_volume_right=wet_volume_middle;
			}

			//inc_correction=upper_correction-lower_correction;

			iteration++;
			if((fabs(correction_old-correction)<tol)&&((fabs(wet_volume-wet_volume_old)/wet_volume_old)<tolv)){
				exit_loop=true;
				//std::cout << "Volume Correction performed: it= "<< iteration <<" Correction =" << correction << " Wet_volume =" << wet_volume << " Net Volume =" << fabs(Net_volume) << std::endl;
				}
			}

		//BLOCK TO BE MODIFIED, JUST COMPUTE CORRECTION IF IT SUPPOSED TO DO IT
		//Now we set the correction to be 0 if it is negative -> if it is positive, the distance of point near 0 positive becomes negotive, so that the front advances
		if(CorrectNegativeVolume==false && correction<0)
		{
			correction=0.0;
			ComputeVolumeAndCuttedAreaInDistance(ThisModelPart,wet_volume,cutted_area,correction);
		}
		// END OF BLOCK TO BE MODIFIED
		// Now we correct the distances
        #pragma omp parallel for firstprivate(node_size)
        for (int ii = 0; ii < node_size; ii++)
        {
			ModelPart::NodesContainerType::iterator it = ThisModelPart.GetCommunicator().LocalMesh().NodesBegin() + ii;
			it->FastGetSolutionStepValue(DISTANCE) -= correction;
        }

        ThisModelPart.GetCommunicator().SynchronizeVariable(DISTANCE);

        ThisModelPart.GetProcessInfo()[CUTTED_AREA] =cutted_area ;
        ThisModelPart.GetProcessInfo()[WET_VOLUME] = wet_volume;
        if (ThisModelPart.GetCommunicator().MyPID() == 0)
            std::cout << "Volume correction : " << correction << " in " << iteration<< std::endl;

        //std::cout << "Volume Correction " << " Net volume: "<< fabs(Net_volume) << " wet volume: " << wet_volume << " percent: "<< wet_volume/fabs(Net_volume)<< " Area: "<< cutted_area << std::endl;
        KRATOS_CATCH("")
    }
    //**********************************************************************************************
    //**********************************************************************************************
    void PosetiveVolumeCorrection(ModelPart& ThisModelPart, const double Net_volume, const double max_correction)
    {
        KRATOS_TRY;

        double wet_volume = 0.0;
        double cutted_area = 0.0;
// 	      int node_size = ThisModelPart.Nodes().size();
        int node_size = ThisModelPart.GetCommunicator().LocalMesh().Nodes().size();

// #pragma omp parallel for firstprivate(node_size) reduction(+:wet_volume,cutted_area )
// 		for (int ii = 0; ii < node_size; ii++)
// 		{
// 		  ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + ii;
//
// 		  wet_volume += it->FastGetSolutionStepValue(WET_VOLUME);
// 		  cutted_area += it->FastGetSolutionStepValue(CUTTED_AREA);
// 		}

        ComputePosVolumeAndCuttedArea(ThisModelPart, wet_volume, cutted_area);


        double volume_difference = fabs(Net_volume) - wet_volume;
        double correction = volume_difference/cutted_area;
        if(correction > max_correction)
            correction = max_correction;
        if(correction < -max_correction)
            correction = -max_correction;

        ThisModelPart.GetProcessInfo()[CUTTED_AREA] =cutted_area ;
        ThisModelPart.GetProcessInfo()[WET_VOLUME] = wet_volume;

        const double liquidus_temp = ThisModelPart.GetProcessInfo()[FLUID_TEMPERATURE];

        //volume loss is just corrected
         if(volume_difference > 0.0)
         {
//             TODO: this is not correct in MPI parallel
            #pragma omp parallel for firstprivate(node_size)
            for (int ii = 0; ii < node_size; ii++)
            {
// 		  ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + ii;
                ModelPart::NodesContainerType::iterator it = ThisModelPart.GetCommunicator().LocalMesh().NodesBegin() + ii;


		//  double alpha = it->FastGetSolutionStepValue(DP_ALPHA1);
		double dist = it->FastGetSolutionStepValue(DISTANCE);
                if(dist < 0 && (dist+correction)>0 )
                    it->FastGetSolutionStepValue(TEMPERATURE) = liquidus_temp;

                dist += correction;

		}
	}

        std::cout << "Volume Correction " << " Net volume: "<< fabs(Net_volume) << " wet volume: " << wet_volume << " percent: "<< wet_volume/fabs(Net_volume)<< " Area: "<< cutted_area << std::endl;
        KRATOS_CATCH("")
    }
    //**********************************************************************************************
    //**********************************************************************************************
    void ComputeNetInletVolume(ModelPart& ThisModelPart)
    {
        KRATOS_TRY;
        double net_input = 0.0;

/*        int node_size = ThisModelPart.GetCommunicator().LocalMesh().Nodes().size();
        #pragma omp parallel for firstprivate(node_size) reduction(+:net_input)
        for (int ii = 0; ii < node_size; ii++)
        {
            ModelPart::NodesContainerType::iterator it = ThisModelPart.GetCommunicator().LocalMesh().NodesBegin() + ii;
            //double str_flag = it->GetValue(IS_STRUCTURE);
            //double slip_flag = it->GetSolutionStepValue(IS_SLIP);
            //double distance = it->GetSolutionStepValue(DISTANCE);

           // if ( (str_flag != 0.0 || slip_flag == 0.0) && distance < 0.0 )
            if ( it->Is(INLET) )
            {
                const array_1d<double, 3> vel = it->FastGetSolutionStepValue(VELOCITY);
                const array_1d<double, 3> normal = it->FastGetSolutionStepValue(NORMAL);

                net_input += inner_prod(vel,normal);
            }
        }
        //syncronoze
        ThisModelPart.GetCommunicator().SumAll(net_input);*/
        for (ModelPart::ConditionIterator iCond = ThisModelPart.ConditionsBegin(); iCond != ThisModelPart.ConditionsEnd(); iCond++)
        {
            if (iCond->GetValue(IS_INLET) != 0.0)
            {
                Geometry< Node<3> >& rGeometry = iCond->GetGeometry();
                array_1d<double,3> v1, v2, AreaNormal;
                v1[0] = rGeometry[1].X() - rGeometry[0].X();
                v1[1] = rGeometry[1].Y() - rGeometry[0].Y();
                v1[2] = rGeometry[1].Z() - rGeometry[0].Z();

                v2[0] = rGeometry[2].X() - rGeometry[0].X();
                v2[1] = rGeometry[2].Y() - rGeometry[0].Y();
                v2[2] = rGeometry[2].Z() - rGeometry[0].Z();

                MathUtils<double>::CrossProduct(AreaNormal,v1,v2);
                AreaNormal *= 0.5;

                array_1d<double,3> Velocity(3,0.0);
                for (unsigned int i = 0; i < rGeometry.PointsNumber(); i++)
                    Velocity += rGeometry[i].FastGetSolutionStepValue(VELOCITY);
                Velocity /= 3.0;

                net_input -= Velocity[0]*AreaNormal[0] + Velocity[1]*AreaNormal[1] + Velocity[2]*AreaNormal[2];
            }
        }
        ThisModelPart.GetCommunicator().SumAll(net_input);

        ProcessInfo& CurrentProcessInfo = ThisModelPart.GetProcessInfo();
        const double delta_t = CurrentProcessInfo[DELTA_TIME];
        double& net_volume = CurrentProcessInfo[NET_INPUT_MATERIAL];
        net_volume += (net_input*delta_t);


        KRATOS_CATCH("")
    }


        //**********************************************************************************************
    //**********************************************************************************************
    /**This function applies a velocity reduction. Velocity is not allowed to be greater in modulus
     * than          old_vel_norm + max_acc_modulus * dt
     * the function is designed to palliate the effect of errors in the solution of the linear system
     * which result in unphysical velocity peaks, typically located at edges.
     */
    void ApplyVelocityLimitation(ModelPart& ThisModelPart, const double max_acc_modulus)
    {
        KRATOS_TRY;
//          double net_input = 0.0;
        int node_size = ThisModelPart.Nodes().size();
        const double dt = ThisModelPart.GetProcessInfo()[DELTA_TIME];

        #pragma omp parallel for firstprivate(node_size)
        //reduction(+:net_input)
        for (int ii = 0; ii < node_size; ii++)
        {
            ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + ii;
            if(!it->IsFixed(VELOCITY_X) && !it->IsFixed(VELOCITY_Y) &&  !it->IsFixed(VELOCITY_Z))
            {
                array_1d<double, 3>& vel = it->FastGetSolutionStepValue(VELOCITY);
                const array_1d<double, 3>& old_vel = it->FastGetSolutionStepValue(VELOCITY,1);
                const double current_vel_norm = norm_2(vel);
                const double old_vel_norm = norm_2(old_vel);

                const double slip_flag = it->FastGetSolutionStepValue(IS_SLIP);

                if(slip_flag > 11.0) //edge or corners -- here we reduce by a factor of 6 the max acceleration
                {
                    const double acceptable_vel_norm = old_vel_norm + 0.1666667*max_acc_modulus*dt;

                    const double ratio = current_vel_norm/acceptable_vel_norm;

                    //velocity is reduced if too high
                    if(ratio > 1.0) vel /= ratio;
                }
                else
                {
                    const double acceptable_vel_norm = old_vel_norm + max_acc_modulus*dt;

                    const double ratio = current_vel_norm/acceptable_vel_norm;

                    //velocity is reduced if too high
                    if(ratio > 1.0) vel /= ratio;
                }
            }

        }

        KRATOS_CATCH("")
    }


    //**********************************************************************************************
    //**********************************************************************************************
    void ComputeNodalVolume(ModelPart& ThisModelPart)
    {
        KRATOS_TRY;
        //first of all set to zero the nodal variables to be updated nodally
        for (ModelPart::NodeIterator i = ThisModelPart.NodesBegin();
                i != ThisModelPart.NodesEnd(); ++i)
        {
            (i)->GetValue(NODAL_VOLUME) = 0.00;
        }

        for (ModelPart::ElementIterator i = ThisModelPart.ElementsBegin();
                i != ThisModelPart.ElementsEnd(); ++i)
        {
            Geometry< Node<3> >& rGeometry = i->GetGeometry();
            double volume = 0.25 * rGeometry.DomainSize()/3.0;//Attention DomainSize() Returns JAcobian/2.0, Volume is Jacobian/6.0

            for (int jj =0; jj<4; ++jj)
                rGeometry[jj].GetValue(NODAL_VOLUME) += volume;
        }

        //ThisModelPart.GetCommunicator().AssembleNonHistoricalData(NODAL_VOLUME);


        KRATOS_CATCH("")
    }
    //**********************************************************************************************
    //**********************************************************************************************
    int SolidificationDuringFilling(ModelPart& ThisModelPart, double BandWidth)
	{
	  KRATOS_TRY;

        //Check for stop criteria
        return StopSolidifCriteria(ThisModelPart,BandWidth);


	  KRATOS_CATCH("")
	}
	//**********************************************************************************************
	//**********************************************************************************************
	void LastStepExtrapolations(ModelPart& ThisModelPart, const double corrected_time )
	{
	  KRATOS_TRY;
      //Check if there is any dry node
	  /* Defining for Solid_fraction_extrapolation*/

	  int node_size = ThisModelPart.GetCommunicator().LocalMesh().Nodes().size();
	  double is_wet = 1.0;
	  #pragma omp parallel for firstprivate(node_size)
	  for (int ii = 0; ii < node_size; ii++)
		{
			ModelPart::NodesContainerType::iterator it = ThisModelPart.GetCommunicator().LocalMesh().NodesBegin() + ii;
			double dist = it->FastGetSolutionStepValue(DISTANCE);
			if(dist>=0.0)
			{
				is_wet = 0.0;
				it->GetValue(FILLTIME)=corrected_time;
				//dist = -1.0;
			}
	  }

        //syncronoze
        ThisModelPart.GetCommunicator().MinAll(is_wet);

		//If there is a dry node then do extrapolation for velocity
		if(is_wet == 0.0)
		{
			ParallelExtrapolationUtilities<3>::ExtrapolateTemperature(ThisModelPart, DISTANCE, TEMPERATURE, NODAL_AREA,10);
			//ParallelExtrapolationUtilities<3>::ParallelExtrapolationUtilities().ExtrapolateVelocity(ThisModelPart, DISTANCE, VELOCITY, NODAL_AREA,10);
			#pragma omp parallel for firstprivate(node_size)
			for (int ii = 0; ii < node_size; ii++)
			{
				ModelPart::NodesContainerType::iterator it = ThisModelPart.GetCommunicator().LocalMesh().NodesBegin() + ii;
				double& dist = it->FastGetSolutionStepValue(DISTANCE);
				if (dist >= 0.0)
				{
					dist = -1.0;


				}

    //        //filling time
    //        double is_visited = it->FastGetSolutionStepValue(IS_VISITED);
    //        if(is_visited == 0.0 && dist<=0.0)
    //        {
    //            it->FastGetSolutionStepValue(IS_VISITED) = 1.0;

    //            //double filling_time =  ThisModelPart.GetProcessInfo()[TIME];
    //            //it->FastGetSolutionStepValue(FILLTIME) = filling_time * time_correction_factor;
				//it->GetValue(FILLTIME) =corrected_time;
    //        }
			}
		}

	  KRATOS_CATCH("")
	}
	//**********************************************************************************************
	//**********************************************************************************************
	void ViscosityBasedSolidification(ModelPart& ThisModelPart, double ViscosityFactor)
	{
	   KRATOS_TRY;

	  int node_size = ThisModelPart.GetCommunicator().LocalMesh().Nodes().size();

	  #pragma omp parallel for firstprivate(node_size)
	  for (int ii = 0; ii < node_size; ii++)
	  {
		ModelPart::NodesContainerType::iterator it = ThisModelPart.GetCommunicator().LocalMesh().NodesBegin() + ii;
		//double temperature = it->FastGetSolutionStepValue(TEMPERATURE);
		const double dist = it->FastGetSolutionStepValue(DISTANCE);
		if(dist<=0.0)
		{
			const double alpha = 1.0 + pow(it->FastGetSolutionStepValue(SOLIDFRACTION),2);
			double& visc = it->FastGetSolutionStepValue(VISCOSITY);
                        visc *= alpha*ViscosityFactor;
			//visc *= (ViscosityFactor - (ViscosityFactor -1.0)*(1.0 - alpha));
		}

	  }
	  KRATOS_CATCH("")
	}
	//**********************************************************************************************
	//**********************************************************************************************
	//void MacroPorosityToShrinkageComputation(ModelPart& ThisModelPart, ModelPart::NodesContainerType& visited_nodes, unsigned int division_number)
	//{
	//  KRATOS_TRY;

	//  double max_porosity = -1000.0;
	//  double min_porosity = 1000000.0;
 //     ModelPart::NodesContainerType& r_nodes = ThisModelPart.Nodes();
 //	  for(ModelPart::NodesContainerType::iterator i_node = r_nodes.begin(); i_node!=r_nodes.end(); i_node++)
	//  {
	//	i_node->Set(NOT_VISITED);
	//	const double mcp = i_node->FastGetSolutionStepValue(MACRO_POROSITY);
	//	if(max_porosity <= mcp)
	//		max_porosity = mcp;
	//	if(min_porosity >= mcp && mcp != 0.0)
	//		min_porosity = mcp;
	//  }

	// // double min_porosity = 0.01*max_porosity;
	///*	  double step_length = (max_porosity - min_porosity)/double(division_number);

	//  double floor_mp =  min_porosity;// + step_length;

 // for(unsigned int cnt = 1; cnt <= division_number; cnt++)
	//  {
	//	  double cnt_val = min_porosity + double(cnt) * step_length;

 //
 //		  for(ModelPart::NodesContainerType::iterator i_node = r_nodes.begin(); i_node!=r_nodes.end(); i_node++)
	//	   {
 //            const double nd_mcp = i_node->FastGetSolutionStepValue(MACRO_POROSITY);
	//		 if( nd_mcp >= floor_mp && i_node->IsNot(VISITED) )
	//		 {
	//			if( nd_mcp <= cnt_val)
	//			{
	//				i_node->Set(VISITED);
	//				i_node->FastGetSolutionStepValue(SHRINKAGE_POROSITY) = cnt_val;
	//				visited_nodes.push_back(*(i_node.base()));
	//			}
	//		 }
	//	   }
	//
	//  }*/
 //
	//double floor_mp =  min_porosity;// + step_length;
 //	for(ModelPart::NodesContainerType::iterator i_node = r_nodes.begin(); i_node!=r_nodes.end(); i_node++)
 //     {
 //       const double nd_mcp = i_node->FastGetSolutionStepValue(MACRO_POROSITY);
	//	if( nd_mcp >= floor_mp )
	//		 {
	//				i_node->FastGetSolutionStepValue(SHRINKAGE_POROSITY) = nd_mcp;
	//				visited_nodes.push_back(*(i_node.base()));

	//		 }

	//   }

	///*for(ModelPart::NodesContainerType::iterator i_node = r_nodes.begin(); i_node!=r_nodes.end(); i_node++)
	//{
 //       const double nd_mcp = i_node->FastGetSolutionStepValue(MACRO_POROSITY);
	//	if(i_node->IsNot(VISITED) && nd_mcp <= floor_mp)
	//		i_node->FastGetSolutionStepValue(SHRINKAGE_POROSITY) = std::numeric_limits<double>::infinity();
	//}*/


	//  KRATOS_CATCH("")
	//}
	////**********************************************************************************************
    //**********************************************************************************************
    void ComputePosetiveVolume(ModelPart& ThisModelPart)
    {
        KRATOS_TRY;

  //     //set as active the internal nodes
  //     int node_size = ThisModelPart.Nodes().size();
  //      #pragma omp parallel for
  //      for (int i = 0; i < node_size; i++)
  //      {
  //          ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + i;
  //          it->FastGetSolutionStepValue(NODAL_MASS) = 0.0;
  //      }
  //      int elem_size = ThisModelPart.Elements().size();
  //      array_1d<double, 4 > N;
  //      BoundedMatrix <double, 4, 3> DN_DX;
		//#pragma omp parallel for private(DN_DX,N) firstprivate(elem_size)
  //      for (int i = 0; i < elem_size; i++)
  //        {
  //              PointerVector< Element>::iterator it = ThisModelPart.ElementsBegin() + i;

  //              Geometry<Node < 3 > >&geom = it->GetGeometry();
  //                  double Volume;
  //                  GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);

  //                  for (unsigned int k = 0; k < 4; k++)
  //                  {
  //                          geom[k].SetLock();
  //                          geom[k].FastGetSolutionStepValue(NODAL_MASS) += Volume*0.25;
  //                          geom[k].UnSetLock();
  //                  }
		// }

  //      double net_input = 0.0;
		//ModelPart::NodesContainerType& r_nodes = ThisModelPart.Nodes();
 	//	for(ModelPart::NodesContainerType::iterator i_node = r_nodes.begin(); i_node!=r_nodes.end(); i_node++)
		// {
		//	 double distance = i_node->GetSolutionStepValue(DISTANCE);
		//	 if(distance >= 0.0)
		//	 {
		//		 double nd_vol = i_node->GetSolutionStepValue(NODAL_MASS);
		//		 net_input += nd_vol;
		//	 }
		//}
		//ProcessInfo& CurrentProcessInfo = ThisModelPart.GetProcessInfo();
		//double& net_volume = CurrentProcessInfo[NET_INPUT_MATERIAL];
		//net_volume = net_input;
		//KRATOS_WATCH(net_volume);

  //      //int node_size = ThisModelPart.GetCommunicator().LocalMesh().Nodes().size();
  //      //#pragma omp parallel for firstprivate(node_size) reduction(+:net_input)
  //      //for (int ii = 0; ii < node_size; ii++)
  //      //{
  //      //    ModelPart::NodesContainerType::iterator it = ThisModelPart.GetCommunicator().LocalMesh().NodesBegin() + ii;
  //      //    double str_flag = it->GetValue(IS_STRUCTURE);
  //      //    double slip_flag = it->GetSolutionStepValue(IS_SLIP);
  //      //    double distance = it->GetSolutionStepValue(DISTANCE);

  //      //   // if ( (str_flag != 0.0 || slip_flag == 0.0) && distance < 0.0 )
  //      //    if ( it->Is(INLET) )
  //      //    {
  //      //        const array_1d<double, 3> vel = it->FastGetSolutionStepValue(VELOCITY);
  //      //        const array_1d<double, 3> normal = it->FastGetSolutionStepValue(NORMAL);

  //      //        net_input += inner_prod(vel,normal);
  //      //    }
  //      //}
  //      ////syncronoze
  //      //ThisModelPart.GetCommunicator().SumAll(net_input);

  //      //ProcessInfo& CurrentProcessInfo = ThisModelPart.GetProcessInfo();
  //      //const double delta_t = CurrentProcessInfo[DELTA_TIME];
  //      //double& net_volume = CurrentProcessInfo[NET_INPUT_MATERIAL];
  //      //net_volume += (net_input*delta_t);


        KRATOS_CATCH("")
    }
        //**********************************************************************************************
    //**********************************************************************************************
    /**This function correct temperature in case that temperature at a node is maximum than
     * FLUID_TEMPERATURE which is the inlet temperature. It simply replaces the temperature of the
     * previous step in this case.
     */
    void ApplyTemperatureLimitation(ModelPart& ThisModelPart, const double max_temperature, const double min_temperature)
    {
        KRATOS_TRY;
        int node_size = ThisModelPart.Nodes().size();


        #pragma omp parallel for firstprivate(node_size)
        for (int ii = 0; ii < node_size; ii++)
        {
            ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + ii;
			double& current_temp = it->FastGetSolutionStepValue(TEMPERATURE);

			if( current_temp > max_temperature || current_temp < min_temperature )//1.05*fluid_temp
			{
				double old_temp = it->FastGetSolutionStepValue(TEMPERATURE,1);
				current_temp = old_temp;
			}
        }

        KRATOS_CATCH("")
    }
    //**********************************************************************************************
    //**********************************************************************************************
    double CheckIfAllNodesAreWet(ModelPart& ThisModelPart)
    {
        KRATOS_TRY;
        int node_size = ThisModelPart.Nodes().size();

        // if there is no dry node
        double is_dry_node = 0.0;
        #pragma omp parallel for firstprivate(node_size)
        for (int ii = 0; ii < node_size; ii++)
        {
            ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + ii;
            double dist = it->FastGetSolutionStepValue(DISTANCE);

            if(dist > 0.0)
                is_dry_node = 1.0;

        }
        //syncronoze
        ThisModelPart.GetCommunicator().MaxAll(is_dry_node);

        return is_dry_node;

        KRATOS_CATCH("")
    }
	//**********************************************************************************************
	//**********************************************************************************************
	/* THIS FUNCTION COMPUTES THE WET VOLUME */
	double ComputeWetVolume(ModelPart& ThisModelPart)
	{
		double wet_volume = 0.0;
        double cutted_area = 0.0;
		// First we compute the Total Volume of the Fluid
		// First we compute the Total Volume of the Fluid
        ComputeWetVolumeAndCuttedArea(ThisModelPart, wet_volume, cutted_area);
		return wet_volume;
	}

	//**********************************************************************************************
    //**********************************************************************************************
    double ComputePartVolume(ModelPart& ThisModelPart)
    {
		double vol=0.0;
		ModelPart::ElementIterator ibegin = ThisModelPart.ElementsBegin();
		unsigned int size_to_loop=ThisModelPart.Elements().size();

        KRATOS_TRY;
		#pragma omp parallel for reduction(+: vol)
        for (int k = 0; k < static_cast<int>(size_to_loop); ++k)
        {
			ModelPart::ElementIterator i = ibegin+k;
            Geometry< Node<3> >& rGeometry = i->GetGeometry();
			double elem_volume = rGeometry.DomainSize(); // Looks like now it is repaired 3.0;//Attention DomainSize() Returns JAcobian/2.0, Volume is Jacobian/6.0
			vol+=elem_volume;
        }
        KRATOS_CATCH("")
		return vol;
    }
	//**********************************************************************************************
    //**********************************************************************************************
    double ComputePartArea(ModelPart& ThisModelPart)
    {
		double area=0.0;
		ModelPart::ConditionIterator ibegin = ThisModelPart.ConditionsBegin();
		unsigned int size_to_loop=ThisModelPart.Conditions().size();
        KRATOS_TRY;
		#pragma omp parallel for reduction(+: area)
        for (int k=0; k <static_cast<int>(size_to_loop); ++k)
        {
			ModelPart::ConditionIterator i=ibegin+k;
            Geometry< Node<3> >& rGeometry = i->GetGeometry();
			double condition_area = rGeometry.DomainSize(); // Looks like now it is repaired 3.0;//Attention DomainSize() Returns JAcobian/2.0, Volume is Jacobian/6.0
			area+=condition_area;
        }
        KRATOS_CATCH("")
		return area;
    }
	//**********************************************************************************************
    //**********************************************************************************************
    double ComputePartInletArea(ModelPart& ThisModelPart)
    {
		double area=0.0;
		ModelPart::ConditionIterator ibegin = ThisModelPart.ConditionsBegin();
		unsigned int size_to_loop=ThisModelPart.Conditions().size();
        KRATOS_TRY;
		#pragma omp parallel for reduction(+: area)
        for (int k=0; k <static_cast<int>(size_to_loop); ++k)
        {
			ModelPart::ConditionIterator i=ibegin+k;
            Geometry< Node<3> >& rGeometry = i->GetGeometry();
			double condition_area = rGeometry.DomainSize(); // Looks like now it is repaired 3.0;//Attention DomainSize() Returns JAcobian/2.0, Volume is Jacobian/6.0
			if(i->GetValue(IS_INLET)>0.0){ area+=condition_area;}
        }
        KRATOS_CATCH("")
		return area;
    }
	//**********************************************************************************************
	//**********************************************************************************************
	double ComputePartMaxh(ModelPart& ThisModelPart)
	{
		KRATOS_TRY
		double h_max = 0.0;
		for (ModelPart::ElementsContainerType::iterator it = ThisModelPart.ElementsBegin(); it != ThisModelPart.ElementsEnd(); it++)
		{
			Geometry<Node<3> >&geom = it->GetGeometry();
			double h = 0.0;

			for (unsigned int i = 0; i<4; i++)
			{
				double xc = geom[i].X();
				double yc = geom[i].Y();
				double zc = geom[i].Z();
				for (unsigned int j = i + 1; j<4; j++)
				{
					double x = geom[j].X();
					double y = geom[j].Y();
					double z = geom[j].Z();
					double l = (x - xc)*(x - xc);
					l += (y - yc)*(y - yc);
					l += (z - zc)*(z - zc);
					if (l > h) h = l;
				}
			}

			h = sqrt(h);

			if (h > h_max) h_max = h;

		}

		ThisModelPart.GetCommunicator().MaxAll(h_max);

		return h_max;

		KRATOS_CATCH("");
	}
	//**********************************************************************************************
	//**********************************************************************************************
	double ComputePartAvgh(ModelPart& ThisModelPart)
	{
		KRATOS_TRY
		double h_avg = 0.0;
		unsigned int n_edges = 0; // It will count held
		for (ModelPart::ElementsContainerType::iterator it = ThisModelPart.ElementsBegin(); it != ThisModelPart.ElementsEnd(); it++)
		{
			Geometry<Node<3> >&geom = it->GetGeometry();
			//double h = 0.0;
			for (unsigned int i = 0; i<4; i++)
			{
				double xc = geom[i].X();
				double yc = geom[i].Y();
				double zc = geom[i].Z();
				for (unsigned int j = i + 1; j<4 ; j++)
				{
					double x = geom[j].X();
					double y = geom[j].Y();
					double z = geom[j].Z();
					double l = (x - xc)*(x - xc);
					l += (y - yc)*(y - yc);
					l += (z - zc)*(z - zc);
					l = sqrt(l);
					h_avg += l;
					n_edges += 1;
				}
			}

			//h = sqrt(h);
			//h_avg += h;

		}
		double n = static_cast<double>(n_edges);

		ThisModelPart.GetCommunicator().SumAll(h_avg);
		ThisModelPart.GetCommunicator().SumAll(n);
		h_avg /= n;
		return h_avg;

		KRATOS_CATCH("");
	}

	//**********************************************************************************************
	//**********************************************************************************************
private:
// 	void AssignDecelerateFactor(ModelPart& ThisModelPart)
// 	{
// 	   KRATOS_TRY;
//
//       //double solidus_temp = ThisModelPart.GetTable(3).Data().front().first;
//       //double liquidus_temp = ThisModelPart.GetTable(3).Data().back().first;
// 	  double solidus_temp=ThisModelPart.GetProcessInfo().GetValue(FLUID_TEMPERATURE);
// 	  double liquidus_temp=ThisModelPart.GetProcessInfo().GetValue(SOLID_TEMPERATURE);
// 	  //double temp_t = liquidus_temp;
// 	  //if ( solidus_temp > liquidus_temp)
// 	  //{
// 		 // liquidus_temp = solidus_temp;
// 		 // solidus_temp = temp_t;
// 	  //}
//
//
// 	  int node_size = ThisModelPart.GetCommunicator().LocalMesh().Nodes().size();
//
// 	  #pragma omp parallel for firstprivate(node_size)
// 	  for (int ii = 0; ii < node_size; ii++)
// 		{
// 			double alpha = 0.0;
// 			ModelPart::NodesContainerType::iterator it = ThisModelPart.GetCommunicator().LocalMesh().NodesBegin() + ii;
//
// 			double dist = it->FastGetSolutionStepValue(DISTANCE);
// 			if( dist<=0.0)
// 			{
// 			  double temperature = it->FastGetSolutionStepValue(TEMPERATURE);
//
//
// 			  if(temperature < solidus_temp) alpha = 1.0;
// 			  else if( temperature >= solidus_temp && temperature < liquidus_temp) alpha = 1.0-(temperature - solidus_temp)/(liquidus_temp - solidus_temp);
//
// 			}
// 			it->FastGetSolutionStepValue(DP_ALPHA1) = alpha;
// 		}
//
// 	  KRATOS_CATCH("")
// 	}
	//**********************************************************************************************
	//**********************************************************************************************
// 	void VelocityReduction(ModelPart& ThisModelPart)
// 	{
// 	   KRATOS_TRY;
//
// 	  int node_size = ThisModelPart.GetCommunicator().LocalMesh().Nodes().size();
//
// 	  #pragma omp parallel for firstprivate(node_size)
// 	  for (int ii = 0; ii < node_size; ii++)
// 	  {
// 		ModelPart::NodesContainerType::iterator it = ThisModelPart.GetCommunicator().LocalMesh().NodesBegin() + ii;
// 		//double temperature = it->FastGetSolutionStepValue(TEMPERATURE);
// 		double alpha = it->FastGetSolutionStepValue(DP_ALPHA1);
//
// 		if(alpha >= 0.9){
// 			it->FastGetSolutionStepValue(VELOCITY_X) = 0.0;
// 			it->FastGetSolutionStepValue(VELOCITY_Y) = 0.0;
// 			it->FastGetSolutionStepValue(VELOCITY_Z) = 0.0;
// 			it->Fix(VELOCITY_X);
// 			it->Fix(VELOCITY_Y);
// 			it->Fix(VELOCITY_Z);
// 		}
// 		else if(alpha<= 0.9 && alpha > 0.1)
// 		{
// 		  double& Vx = it->FastGetSolutionStepValue(VELOCITY_X);
// 		  double& Vy = it->FastGetSolutionStepValue(VELOCITY_Y);
// 		  double& Vz = it->FastGetSolutionStepValue(VELOCITY_Z);
// 		  Vx *= (1.0-alpha);
// 		  Vy *= (1.0-alpha);
// 		  Vz *= (1.0-alpha);
//
// 		  it->Fix(VELOCITY_X);
// 		  it->Fix(VELOCITY_Y);
// 		  it->Fix(VELOCITY_Z);
// 		}
// 	  }
// 	  KRATOS_CATCH("")
// 	}

	//**********************************************************************************************
	//**********************************************************************************************
	int StopSolidifCriteria(ModelPart& ThisModelPart, double ref_dist)
	{
	   KRATOS_TRY;

	  int is_hot = 0;
	  int node_size = ThisModelPart.GetCommunicator().LocalMesh().Nodes().size();

      //CAN NOT DO THIS IN PARALLEL! it is just wrong!!!
	  for (int ii = 0; ii < node_size; ii++)
	  {
		ModelPart::NodesContainerType::iterator it = ThisModelPart.GetCommunicator().LocalMesh().NodesBegin() + ii;
		//double temperature = it->FastGetSolutionStepValue(TEMPERATURE);
		double distance = it->FastGetSolutionStepValue(DISTANCE);
		if(distance <= 0.0 && fabs(distance) <= ref_dist){
			double solid_fraction = it->FastGetSolutionStepValue(SOLIDFRACTION);
			if(solid_fraction < 0.9)
            {
				is_hot=1;
                break;
            }
		}
	  }


	  ThisModelPart.GetCommunicator().MaxAll(is_hot);
	  return is_hot;

	  KRATOS_CATCH("")
	}

	//**********************************************************************************************
	//**********************************************************************************************

    void AirSmagorinskey(ModelPart& ThisModelPart, double C_Smagorinsky)
    {
        int elem_size = ThisModelPart.Elements().size();

        #pragma omp parallel for firstprivate(elem_size)
        for(int ii = 0; ii<elem_size; ii++)
        {
            PointerVector< Element>::iterator iel=ThisModelPart.ElementsBegin()+ii;
            double dist_sign = 1.0;
            Geometry< Node<3> >& geom = iel->GetGeometry();
            for(unsigned int i =0; i<geom.size(); i++)
            {
                double dist = geom[i].FastGetSolutionStepValue(DISTANCE);
                if(dist_sign*dist < 0.0)
                {
                    dist_sign = -1.0;
                    break;
                }
            }
            // to be sure to not apply to the cutted elements and just to all air elements
            if(dist_sign == 1.0)
                iel->SetValue(C_SMAGORINSKY, C_Smagorinsky );

        }
    }
    //**********************************************************************************************
    //**********************************************************************************************
    void ComputeWetVolumeAndCuttedArea(ModelPart& ThisModelPart, double& wet_volume, double& cutted_area)
    {
        KRATOS_TRY;
        int elem_size = ThisModelPart.Elements().size();
        double wetvol = 0.0;
        double cutare = 0.0;

        #pragma omp parallel for firstprivate(elem_size) reduction(+:wetvol,cutare)
        for(int ii = 0; ii<elem_size; ii++)
        {
            PointerVector< Element>::iterator iel=ThisModelPart.ElementsBegin()+ii;

            // Calculate this element's geometric parameters
            double Area;
            array_1d<double, 4> N;
            BoundedMatrix<double, 4, 3> DN_DX;
            GeometryUtils::CalculateGeometryData(iel->GetGeometry(), DN_DX, N, Area);
            //get position of the cut surface
            Vector distances(4);
            Matrix Nenriched(6, 1);
            Vector volumes(6);
            Matrix coords(4, 3);
            Matrix Ngauss(6, 4);
            Vector signs(6);
            std::vector< Matrix > gauss_gradients(6);
            //fill coordinates
            for (unsigned int i = 0; i < 4; i++)
            {
                const array_1d<double, 3 > & xyz = iel->GetGeometry()[i].Coordinates();
                volumes[i] = 0.0;
                distances[i] = iel->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
                for (unsigned int j = 0; j < 3; j++)
                    coords(i, j) = xyz[j];
            }
            for (unsigned int i = 0; i < 6; i++)
                gauss_gradients[i].resize(1, 3, false);

            array_1d<double,6> edge_areas;
            unsigned int ndivisions = EnrichmentUtilities::CalculateTetrahedraEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched,edge_areas);

            if(ndivisions == 1)
            {
                if( signs[0] < 0.0)
                    wetvol += volumes[0];
            }
            else
            {
                double ele_wet_volume=0.0;
                for (unsigned int kk = 0; kk < ndivisions; kk++)
                {
                    if( signs[kk]<0.0 )
                        ele_wet_volume += volumes[kk];
                }
                wetvol += ele_wet_volume;

                for(unsigned int i=0; i<6; i++)
                    cutare += edge_areas[i];
                //cutare += 1.80140543 * pow(ele_wet_volume,0.666666666667); // equilateral tetrahedraon is considered
            }
        }
        //syncronoze
        ThisModelPart.GetCommunicator().SumAll(wetvol);
        ThisModelPart.GetCommunicator().SumAll(cutare);

        wet_volume = wetvol;
        cutted_area = cutare;

        KRATOS_CATCH("")
    }
	 //**********************************************************************************************
    //**********************************************************************************************
    void ComputeVolumeAndCuttedAreaInDistance(ModelPart& ThisModelPart, double& wet_volume, double& cutted_area, const double& reference_distance)
    {
        KRATOS_TRY;
        int elem_size = ThisModelPart.Elements().size();
        double wetvol = 0.0;
        double cutare = 0.0;

        #pragma omp parallel for firstprivate(elem_size) reduction(+:wetvol,cutare)
        for(int ii = 0; ii<elem_size; ii++)
        {
            PointerVector< Element>::iterator iel=ThisModelPart.ElementsBegin()+ii;

            // Calculate this element's geometric parameters
            double Area;
            array_1d<double, 4> N;
            BoundedMatrix<double, 4, 3> DN_DX;
            GeometryUtils::CalculateGeometryData(iel->GetGeometry(), DN_DX, N, Area);
            //get position of the cut surface
            Vector distances = ZeroVector(4);
            Matrix Nenriched = ZeroMatrix(6, 1);
            Vector volumes = ZeroVector(6);
            Matrix coords = ZeroMatrix(4, 3);
            Matrix Ngauss = ZeroMatrix(6, 4);
            Vector signs = ZeroVector(6);
            std::vector< Matrix > gauss_gradients(6);
            //fill coordinates
            for (unsigned int i = 0; i < 4; i++)
            {
                const array_1d<double, 3 > & xyz = iel->GetGeometry()[i].Coordinates();
                volumes[i] = 0.0;
                distances[i] = (iel->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE))-reference_distance;
                for (unsigned int j = 0; j < 3; j++)
                    coords(i, j) = xyz[j];
            }
            for (unsigned int i = 0; i < 6; i++)
            {
                gauss_gradients[i].resize(1, 3, false);
                gauss_gradients[i] = ZeroMatrix(1,3);
            }

            array_1d<double,6> edge_areas;
            unsigned int ndivisions = EnrichmentUtilities::CalculateTetrahedraEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched,edge_areas);

            if(ndivisions == 1)
            {
                if( signs[0] < 0.0)
                    wetvol += volumes[0];
            }
            else
            {
                double ele_wet_volume=0.0;
                for (unsigned int kk = 0; kk < ndivisions; kk++)
                {
                    if( signs[kk]<0.0 )
                        ele_wet_volume += volumes[kk];
                }
                wetvol += ele_wet_volume;

                for(unsigned int i=0; i<6; i++)
                    cutare += edge_areas[i];
                //cutare += 1.80140543 * pow(ele_wet_volume,0.666666666667); // equilateral tetrahedraon is considered
            }
        }
        //syncronoze
        ThisModelPart.GetCommunicator().SumAll(wetvol);
        ThisModelPart.GetCommunicator().SumAll(cutare);

        wet_volume = wetvol;
        cutted_area = cutare;

        KRATOS_CATCH("")
    }
    //**********************************************************************************************
    //**********************************************************************************************
    void ComputePosVolumeAndCuttedArea(ModelPart& ThisModelPart, double& pos_volume, double& cutted_area)
    {
        KRATOS_TRY;
        int elem_size = ThisModelPart.Elements().size();
        double wetvol = 0.0;
        double cutare = 0.0;

        #pragma omp parallel for firstprivate(elem_size) reduction(+:wetvol,cutare)
        for(int ii = 0; ii<elem_size; ii++)
        {
            PointerVector< Element>::iterator iel=ThisModelPart.ElementsBegin()+ii;

            // Calculate this element's geometric parameters
            double Area;
            array_1d<double, 4> N;
            BoundedMatrix<double, 4, 3> DN_DX;
            GeometryUtils::CalculateGeometryData(iel->GetGeometry(), DN_DX, N, Area);
            //get position of the cut surface
            Vector distances(4);
            Matrix Nenriched(6, 1);
            Vector volumes(6);
            Matrix coords(4, 3);
            Matrix Ngauss(6, 4);
            Vector signs(6);
            std::vector< Matrix > gauss_gradients(6);
            //fill coordinates
            for (unsigned int i = 0; i < 4; i++)
            {
                const array_1d<double, 3 > & xyz = iel->GetGeometry()[i].Coordinates();
                volumes[i] = 0.0;
                distances[i] = iel->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
                for (unsigned int j = 0; j < 3; j++)
                    coords(i, j) = xyz[j];
            }
            for (unsigned int i = 0; i < 6; i++)
                gauss_gradients[i].resize(1, 3, false);

            array_1d<double,6> edge_areas;
            unsigned int ndivisions = EnrichmentUtilities::CalculateTetrahedraEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched,edge_areas);

            if(ndivisions == 1)
            {
                if( signs[0] > 0.0)
                    wetvol += volumes[0];
            }
            else
            {
                double ele_wet_volume=0.0;
                for (unsigned int kk = 0; kk < ndivisions; kk++)
                {
                    if( signs[kk]>0.0 )
                        ele_wet_volume += volumes[kk];
                }
                wetvol += ele_wet_volume;

                for(unsigned int i=0; i<6; i++)
                    cutare += edge_areas[i];
                //cutare += 1.80140543 * pow(ele_wet_volume,0.666666666667); // equilateral tetrahedraon is considered
            }
        }
        //syncronoze
        ThisModelPart.GetCommunicator().SumAll(wetvol);
        ThisModelPart.GetCommunicator().SumAll(cutare);

        pos_volume = wetvol;
        cutted_area = cutare;

        KRATOS_CATCH("")
    }
    double ComputeCharactristicCuttedLength(ModelPart& ThisModelPart)
    {
        KRATOS_TRY;

        double max_cutted_len = 0.0;
        double cnt = 0.0;
        array_1d<double,4> dist;
        int elem_size = ThisModelPart.Elements().size();

        #pragma omp parallel for private(dist) firstprivate(elem_size)
        for (int i = 0; i < elem_size; i++)
        {
            PointerVector< Element>::iterator it = ThisModelPart.ElementsBegin() + i;

            Geometry<Node < 3 > >& element_geometry = it->GetGeometry();

            for (unsigned int i = 0; i < 4; i++)
                dist[i] = element_geometry[i].FastGetSolutionStepValue(DISTANCE);

            bool is_divided = IsDivided(dist);
            if (is_divided == true)
            {
                #pragma omp atomic
                cnt+= 1.0;

                double max_pos_dist = 0.0;
                double min_neg_dist = 0.0;
                for (unsigned int ii = 0; ii < 4; ii++)
                {
                    if ( dist[ii] > max_pos_dist)
                        max_pos_dist = dist[ii];
                    else if( dist[ii] < min_neg_dist)
                        min_neg_dist =  dist[ii];
                }

                double this_ele_dist = max_pos_dist -  min_neg_dist ;
                if(this_ele_dist > max_cutted_len)
                    max_cutted_len = this_ele_dist;

            }
        }
        return max_cutted_len;

        KRATOS_CATCH("")
    }

    bool IsDivided(array_1d<double,4>& dist)
    {
        unsigned int positive = 0;
        unsigned int negative = 0;

        for(unsigned int i=0; i<4; i++)
        {
            if(dist[i] >= 0)
                positive++;
            else
                negative++;
        }

        bool is_divided = false;
        if(positive > 0 && negative>0)
            is_divided = true;

        return is_divided;
    }

	void CorrectTemperatureInAddedVolume(ModelPart& ThisModelPart, const double correction)
	{

	  //compute min temp in wet part ( edge and corners are excluded)
	  double min_wet_temp = ThisModelPart.GetProcessInfo()[FLUID_TEMPERATURE];

	  #pragma omp parallel for
      for (int k = 0; k< static_cast<int> (ThisModelPart.Nodes().size()); k++)
		 {
			ModelPart::NodesContainerType::iterator i_node = ThisModelPart.NodesBegin() + k;
			double distance = i_node->GetSolutionStepValue(DISTANCE);
			double slip_falg = i_node->GetSolutionStepValue(IS_SLIP);
			if(distance <=0.0 && slip_falg!=20.0 && slip_falg!=30.0)
				 {
				  double temp_node = i_node->FastGetSolutionStepValue(TEMPERATURE);
				  if (temp_node < min_wet_temp)
						 min_wet_temp = temp_node;
				 }

		 }

	   //assign to TEMPERATURE in the added zone
       #pragma omp parallel for
	   for (int k = 0; k< static_cast<int> (ThisModelPart.Nodes().size()); k++)
		  {
			ModelPart::NodesContainerType::iterator i_node = ThisModelPart.NodesBegin() + k;
			double distance = i_node->GetSolutionStepValue(DISTANCE);

			if(distance >= 0.0 && distance <= 5.0 * correction)
			{
				i_node->FastGetSolutionStepValue(TEMPERATURE) = min_wet_temp;
				i_node->FastGetSolutionStepValue(TEMPERATURE,1) = min_wet_temp;

			}

	      }

	}
};

}  // namespace Kratos.

#endif // KRATOS_BIPHASIC_FILLING_UTILITIES_INCLUDED  defined


