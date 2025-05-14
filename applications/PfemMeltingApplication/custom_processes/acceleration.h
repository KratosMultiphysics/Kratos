//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Julio Marti
//


#if !defined(KRATOS_FLUID_PRESSURE_CALCULATE_PROCESS_INCLUDED )
#define  KRATOS_FLUID_PRESSURE_CALCULATE_PROCESS_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "utilities/geometry_utilities.h"
#include "pfem_melting_application.h"
#include "custom_elements/hypo.h"


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
	Update the PRESSURE_FORCE on the nodes


*/

class FluidPressureCalculateProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of HypoelasticStressCalculateProcess
    KRATOS_CLASS_POINTER_DEFINITION(FluidPressureCalculateProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FluidPressureCalculateProcess(ModelPart& model_part, unsigned int domain_size)
        : mr_model_part(model_part),mdomain_size(domain_size)
    {
    }

    /// Destructor.
    ~FluidPressureCalculateProcess() override
    {
    }


    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {
        KRATOS_TRY

        ProcessInfo& proc_info = mr_model_part.GetProcessInfo();


	if (mdomain_size==2)
	   {
           //Matrix dummy=ZeroMatrix(2,2);
           double Area;
           double dummy;
           array_1d<double, 3> N;

           BoundedMatrix<double, 3, 2> DN_DX;
           

           //THIS SHOULD BE EXECUTED ONLY FOR THOSE ELEMENTS OF THE MONOLITHIC MODEL THAT ARE IDENTIFIED TO BELONG TO THE SOLID domain
	   //before the first step we initialize Cauchy stress to zero
//	   if (proc_info[TIME]==0.0)
// 	      {
//                for(ModelPart::NodesContainerType::iterator i = mr_model_part.NodesBegin(); i!=mr_model_part.NodesEnd(); i++)
//            	{
//                i->FastGetSolutionStepValue(PRESSUREAUX)=0.0;
//                i->FastGetSolutionStepValue(NODAL_VOLUME) = 0.0;
//          	}
//	      }
	   //and now we actually compute it
//	   else
//	     {
	        for(ModelPart::NodesContainerType::iterator i = mr_model_part.NodesBegin(); i!=mr_model_part.NodesEnd(); i++)
            	{
                i->FastGetSolutionStepValue(PRESSUREAUX)=0.0;
                i->FastGetSolutionStepValue(NODAL_VOLUME) = 0.0;
          	}
	     //}



	     int solid_nodes = 0;
             for(ModelPart::ElementsContainerType::iterator im = mr_model_part.ElementsBegin() ; im != mr_model_part.ElementsEnd() ; ++im)
                 {
	         //IN A MONOLITHIC FLUID-SOLID MODEL WE WANT TO EXCEUTE THIS FUNCTION ONLY FOR THE SOLID ELEMENTS
                 //if(im->GetGeometry()[0].Is(STRUCTURE) && im->GetGeometry()[1].Is(STRUCTURE) && im->GetGeometry()[2].Is(STRUCTURE))
                 //if(im->GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE) == 1.0 && im->GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE) == 1.0 && im->GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE) == 1.0 )

                    solid_nodes = 0;
                 	
                 		Geometry< Node<3> >& geom = im->GetGeometry();

				for(unsigned int i2 = 0; i2 < geom.size(); i2++)
				{
					if(geom[i2].FastGetSolutionStepValue(IS_SOLID) == 1) solid_nodes += 1;
					//KRATOS_WATCH("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")

					//KRATOS_WATCH("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")										
					//KRATOS_WATCH(solid_nodes)
				}	    
                    
                	if (solid_nodes<3){
                               //KRATOS_ERROR << "AddTimeIntegratedLHS is not implemented." << std::endl; 
                    		im->Calculate(NODAL_VOLUME,dummy,proc_info);


		     		GeometryUtils::CalculateGeometryData(im->GetGeometry(), DN_DX, N, Area);
                    		Geometry< Node >& geom = im->GetGeometry();
                    
                    		const double dt = mr_model_part.GetProcessInfo()[DELTA_TIME];
                            

  		        	double K=-0.333333333333333333*(geom[0].FastGetSolutionStepValue(BULK_MODULUS)+geom[1].FastGetSolutionStepValue(BULK_MODULUS)+geom[2].FastGetSolutionStepValue(BULK_MODULUS));
                       	
                       	//KRATOS_WATCH(K)
                       	
                       	//KRATOS_ERROR << "AddTimeIntegratedLHS is not implemented." << std::endl;                	
                       	//K=-100000.0;
                       	K=-25000.0;
                       	K=-500.0;

				//calculate the divergence
				const array_1d<double,3>& v0 = geom[0].FastGetSolutionStepValue(VELOCITY);
				const array_1d<double,3>& v1 = geom[1].FastGetSolutionStepValue(VELOCITY);
				const array_1d<double,3>& v2 = geom[2].FastGetSolutionStepValue(VELOCITY);

				double div_v = DN_DX(0,0)*v0[0] + DN_DX(0,1)*v0[1] + DN_DX(1,0)*v1[0] + DN_DX(1,1)*v1[1] + DN_DX(2,0)*v2[0] + DN_DX(2,1)*v2[1];
				double dp_el = K * dt * div_v * Area;
		               //dp_el = K * div_v * Area;    
				array_1d<double,3> pn = ZeroVector(3); //dimension = number of nodes
				pn[0] = geom[0].FastGetSolutionStepValue(PRESSURE,1);
				pn[1] = geom[1].FastGetSolutionStepValue(PRESSURE,1);
				pn[2] = geom[2].FastGetSolutionStepValue(PRESSURE,1);
				double diag_term = Area/6.0;
				double out_term = Area/12.0;


				geom[0].FastGetSolutionStepValue(PRESSUREAUX) += dp_el *0.333333333333333333 + pn[0] * 0.333333333333333333 * Area;
				geom[1].FastGetSolutionStepValue(PRESSUREAUX) += dp_el *0.333333333333333333 + pn[1] * 0.333333333333333333 * Area;
				geom[2].FastGetSolutionStepValue(PRESSUREAUX) += dp_el *0.333333333333333333 + pn[2] * 0.333333333333333333 * Area;
				}
		}
                    
                    for(ModelPart::NodesContainerType::iterator i = mr_model_part.NodesBegin(); i!=mr_model_part.NodesEnd(); i++)
            	{
            	
            	 const double& ar= i->FastGetSolutionStepValue(NODAL_VOLUME);
                //KRATOS_WATCH(i->FastGetSolutionStepValue(NODAL_VOLUME))


            	 double aux=ar;
            	 if(aux<0.00000000000001) aux=1.0;
            	 //aux=1.0;
            	 i->FastGetSolutionStepValue(PRESSUREAUX)/=aux;
                //i->FastGetSolutionStepValue(PRESSURE,1)=i->FastGetSolutionStepValue(PRESSUREAUX);
                i->FastGetSolutionStepValue(PRESSURE)=i->FastGetSolutionStepValue(PRESSUREAUX);
                //if(i->FastGetSolutionStepValue(IS_FREE_SURFACE)==1) i->FastGetSolutionStepValue(PRESSURE)=0.0;    
                 }
                 
                 
            //KRATOS_ERROR << "AddTimeIntegratedLHS is not implemented." << std::endl; 




}



	ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
	ModelPart& model_part=BaseType::GetModelPart();
	double dt = model_part.GetProcessInfo()[DELTA_TIME];
	array_1d<double, 3 > zero = ZeroVector(3);

	int number_of_threads=0;
	#ifdef _OPENMP
		number_of_threads = omp_get_max_threads();
	#else
		number_of_threads = 1;
	#endif

	for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();i != BaseType::GetModelPart().NodesEnd(); ++i)
	 {
	   array_1d<double, 3 > zero = ZeroVector(3);
	   i->FastGetSolutionStepValue(FORCE)=ZeroVector(3);
	   (i)->FastGetSolutionStepValue(NODAL_MASS)=0.0;
		    i->FastGetSolutionStepValue(ANGULAR_ACCELERATION)=ZeroVector(3);
	   //nodal_mass = 0.0;
	 }

	rCurrentProcessInfo[FRACTIONAL_STEP] = 6;
	for (ModelPart::ElementIterator i = BaseType::GetModelPart().ElementsBegin(); i != BaseType::GetModelPart().ElementsEnd(); ++i)
	 {
		    (i)->InitializeSolutionStep(BaseType::GetModelPart().GetProcessInfo());
	 }


		    //boost::numeric::ublas::bounded_matrix<double,4,3> DN_DX;
		    //array_1d<double,4> N;
		    //if nd=2
		    boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX;
		    array_1d<double,3> N;


	double volume;
	//array_1d<double,3> aux0, aux1, aux2, aux3; //this are sized to 3 even in 2D!!
	for (typename ModelPart::ElementsContainerType::iterator it=model_part.ElementsBegin(); it!=model_part.ElementsEnd(); ++it)
	{
		Geometry< Node<3> >& geom = it->GetGeometry();


		const array_1d<double, 3 > & fv0 = geom[0].FastGetSolutionStepValue(VELOCITY);
		const array_1d<double, 3 > & fv1 = geom[1].FastGetSolutionStepValue(VELOCITY);
		const array_1d<double, 3 > & fv2 = geom[2].FastGetSolutionStepValue(VELOCITY);


		const double rho0 = geom[0].FastGetSolutionStepValue(DENSITY);
		const double rho1 = geom[1].FastGetSolutionStepValue(DENSITY);
		const double rho2 = geom[2].FastGetSolutionStepValue(DENSITY);
		//const double rho3 = geom[3].FastGetSolutionStepValue(DENSITY);

		//double density = 0.25 *(rho0 + rho1 + rho2 + rho3);
		double density = 0.3333333 *(rho0 + rho1 + rho2 );
		       
		GeometryUtils::CalculateGeometryData(geom, DN_DX, N, volume);


	       array_1d<double, 3 > aux;
	       array_1d<double, 3 > vel_gauss;
	       noalias(aux) = fv0;
	       noalias(vel_gauss) = N[0] * aux;

	       noalias(aux) = fv1;
	       noalias(vel_gauss) += N[1] * aux;

	       noalias(aux) = fv2;
	       noalias(vel_gauss) += N[2] * aux;

	       //noalias(aux) = fv3;
	       //noalias(vel_gauss) += N[3] * aux;

	       //CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
	       //array_1d<double, 4 > u_DN;
	       //boost::numeric::ublas::bounded_matrix<double, 4, 4 > LeftHandSideMatrix;

	       array_1d<double, 3 > u_DN;
	       boost::numeric::ublas::bounded_matrix<double, 3, 3 > LeftHandSideMatrix;
		           
		       

	       noalias(u_DN) = prod(DN_DX, vel_gauss);
	       noalias(LeftHandSideMatrix) = outer_prod(N, u_DN);

	       //multiplication by the Volume
	       LeftHandSideMatrix *= (volume * density);



	   //====================================================================
	   //====================================================================
	   //====================================================================
	   //con la estabilizacion

		    const double nu0 = geom[0].FastGetSolutionStepValue(VISCOSITY);
		    const double nu1 = geom[1].FastGetSolutionStepValue(VISCOSITY);
		    const double nu2 = geom[2].FastGetSolutionStepValue(VISCOSITY);



		    boost::numeric::ublas::bounded_matrix<double, 3, 3 > WorkMatrix;
		    LeftHandSideMatrix *= 0.0;


	   //calculating average density and viscosity
	   double nu = 0.33333333333333*(nu0 + nu1 + nu2 );
	   
	   //ms_vel_gauss[i] =  msN[0]*(fv0[i]) + msN[1]*(fv1[i]) +  msN[2]*(fv2[i]);
	   //but with one integration N=0.333333333
	   vel_gauss[0] =  0.33333333333333*(fv0[0]+fv1[0]+fv2[0]);
	   vel_gauss[1] =  0.33333333333333*(fv0[1]+fv1[1]+fv2[1]);

	   //calculating parameter tau (saved internally to each element)
	   double h = sqrt(2.00*volume);
	   double norm_u = vel_gauss[0]*vel_gauss[0] + vel_gauss[1]*vel_gauss[1];
	   norm_u = sqrt(norm_u);
	   double tau = 1.00 / ( 1.0 * 4.00*nu/(h*h) + 2.00*norm_u/h + 1.0 * 1.0/dt);


		    //tau *=0.0;
	   //====================================================================
	   //calculation of convective and stabilizing terms ... using 3 gauss points (on the sides)
	   // double c1 = 4.00; double c2 = 2.00;

	   //double norm_u = 0.0;
	   //double tau = 0.0;
	   double area_density_third = volume * density * 0.33333333333333333333;

	   // ******************* GAUSS1 ***************************
	   N[0] = 0.5;
	   N[1] = 0.5;
	   N[2] = 0.0;

	   // vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points)
	   //(note that the fractional step vel is used)
	   vel_gauss[0] = N[0] * fv0[0] + N[1] * fv1[0] + N[2] * fv2[0];
	   vel_gauss[1] = N[0] * fv0[1] + N[1] * fv1[1] + N[2] * fv2[1];
	   // KRATOS_WATCH(vel_gauss);

	   //calculating parameter tau
	   //norm_u = vel_gauss[0] * vel_gauss[0] + vel_gauss[1] * vel_gauss[1];
	   //norm_u = sqrt(norm_u);

	//         tau = CalculateTau(h, nu, norm_u, rCurrentProcessInfo);
	   //tau = CalculateTau(DN_DX, vel_gauss, h, nu, norm_u, rCurrentProcessInfo);

	   //CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
	   noalias(u_DN) = prod(DN_DX, vel_gauss);
	   noalias(WorkMatrix) = outer_prod(N, u_DN);

	   //CONVECTION STABILIZING CONTRIBUTION (Suu)
	   noalias(WorkMatrix) += tau * outer_prod(u_DN, u_DN);

	   //adding gauss point contribution
	   noalias(LeftHandSideMatrix) += area_density_third * WorkMatrix;

	   
	   // ******************* GAUSS2 ***************************
	   N[0] = 0.0;
	   N[1] = 0.5;
	   N[2] = 0.5;

	   // vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points)
	   //(note that the fractional step vel is used)
	   vel_gauss[0] = N[0] * fv0[0] + N[1] * fv1[0] + N[2] * fv2[0];
	   vel_gauss[1] = N[0] * fv0[1] + N[1] * fv1[1] + N[2] * fv2[1];
	   // KRATOS_WATCH(vel_gauss);

	   //calculating parameter tau
	   //norm_u = vel_gauss[0] * vel_gauss[0] + vel_gauss[1] * vel_gauss[1];
	   //norm_u = sqrt(norm_u);

	//         tau = CalculateTau(h, nu, norm_u, rCurrentProcessInfo);
	   //tau = CalculateTau(DN_DX, vel_gauss, h, nu, norm_u, rCurrentProcessInfo);

	   //CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
	   noalias(u_DN) = prod(DN_DX, vel_gauss);
	   noalias(WorkMatrix) = outer_prod(N, u_DN);

	   //CONVECTION STABILIZING CONTRIBUTION (Suu)
	   noalias(WorkMatrix) += tau * outer_prod(u_DN, u_DN);

	   //adding gauss point contribution
	   noalias(LeftHandSideMatrix) += area_density_third * WorkMatrix;

	   // ******************* GAUSS3 ***************************
	   N[0] = 0.5;
	   N[1] = 0.0;
	   N[2] = 0.5;

	   // vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points)
	   //(note that the fractional step vel is used)
	   vel_gauss[0] = N[0] * fv0[0] + N[1] * fv1[0] + N[2] * fv2[0];
	   vel_gauss[1] = N[0] * fv0[1] + N[1] * fv1[1] + N[2] * fv2[1];
	   // KRATOS_WATCH(vel_gauss);

	   //calculating parameter tau
	   //norm_u = vel_gauss[0] * vel_gauss[0] + vel_gauss[1] * vel_gauss[1];
	   //norm_u = sqrt(norm_u);


	//         tau = CalculateTau(h, nu, norm_u, rCurrentProcessInfo);
	   //tau = CalculateTau(DN_DX, vel_gauss, h, nu, norm_u, rCurrentProcessInfo);

	   // KRATOS_WATCH(tau);

	   //CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
	   noalias(u_DN) = prod(DN_DX, vel_gauss);
	   noalias(WorkMatrix) = outer_prod(N, u_DN);

	   //CONVECTION STABILIZING CONTRIBUTION (Suu)
	   noalias(WorkMatrix) += tau * outer_prod(u_DN, u_DN);

	   //adding gauss point contribution
	   noalias(LeftHandSideMatrix) += area_density_third * WorkMatrix;


	//
	//
	//

	       Vector TempVecX(3);
	       Vector TempVecY(3);
	       Vector TempVecZ(3);
	       //Vector LumpedM(3);
		        Vector temp_vec_np(3);
	       temp_vec_np[0]=geom[0].FastGetSolutionStepValue(VELOCITY_X);
	       temp_vec_np[1]=geom[1].FastGetSolutionStepValue(VELOCITY_X);
	       temp_vec_np[2]=geom[2].FastGetSolutionStepValue(VELOCITY_X);
	       //temp_vec_np[3]=geom[3].FastGetSolutionStepValue(VELOCITY_X,1);

		        TempVecX=prod(LeftHandSideMatrix, temp_vec_np);


	       geom[0].FastGetSolutionStepValue(FORCE_X) += TempVecX[0];
	       geom[1].FastGetSolutionStepValue(FORCE_X) += TempVecX[1];
	       geom[2].FastGetSolutionStepValue(FORCE_X) += TempVecX[2];
	       //geom[3].FastGetSolutionStepValue(FORCE_X) += TempVecX[3];



	       //rhs term
	       
	       temp_vec_np[0]=geom[0].FastGetSolutionStepValue(VELOCITY_Y);
	       temp_vec_np[1]=geom[1].FastGetSolutionStepValue(VELOCITY_Y);
	       temp_vec_np[2]=geom[2].FastGetSolutionStepValue(VELOCITY_Y);
	       //temp_vec_np[3]=geom[3].FastGetSolutionStepValue(VELOCITY_Y,1);

		        TempVecY=prod(LeftHandSideMatrix, temp_vec_np);

	       geom[0].FastGetSolutionStepValue(FORCE_Y) += TempVecY[0];
	       geom[1].FastGetSolutionStepValue(FORCE_Y) += TempVecY[1];
	       geom[2].FastGetSolutionStepValue(FORCE_Y) += TempVecY[2];
	       //geom[3].FastGetSolutionStepValue(FORCE_Y) += TempVecY[3];


	       
	/*        temp_vec_np[0]=geom[0].FastGetSolutionStepValue(VELOCITY_Z);
	       temp_vec_np[1]=geom[1].FastGetSolutionStepValue(VELOCITY_Z);
	       temp_vec_np[2]=geom[2].FastGetSolutionStepValue(VELOCITY_Z);
	       //temp_vec_np[3]=geom[3].FastGetSolutionStepValue(VELOCITY_Z,1);
	       TempVecZ=prod(LeftHandSideMatrix, temp_vec_np);

	       geom[0].FastGetSolutionStepValue(FORCE_Z) += TempVecZ[0];
	       geom[1].FastGetSolutionStepValue(FORCE_Z) += TempVecZ[1];
	       geom[2].FastGetSolutionStepValue(FORCE_Z) += TempVecZ[2];
	*/
	       //geom[3].FastGetSolutionStepValue(FORCE_Z) += TempVecZ[3];


	/*
	double rho = 0.25*(geom[0].FastGetSolutionStepValue(DENSITY) + geom[1].FastGetSolutionStepValue(DENSITY)+ geom[2].FastGetSolutionStepValue(DENSITY)+ geom[3].FastGetSolutionStepValue(DENSITY));

		       double nu = 0.25*(geom[0].FastGetSolutionStepValue(VISCOSITY) + geom[1].FastGetSolutionStepValue(VISCOSITY)+ geom[2].FastGetSolutionStepValue(VISCOSITY)+ geom[3].FastGetSolutionStepValue(VISCOSITY));

	      double p_avg = 0.25*(geom[0].FastGetSolutionStepValue(PRESSURE) + geom[1].FastGetSolutionStepValue(PRESSURE)+ geom[2].FastGetSolutionStepValue(PRESSURE)+ geom[3].FastGetSolutionStepValue(PRESSURE));


		        Vector TempVecX(4);
	       Vector TempVecY(4);
	       Vector TempVecZ(4);
	       //Vector LumpedM(3);
		        Vector temp_vec_np(4);
		        double fcomp=0.0;
	       temp_vec_np[0]=geom[0].FastGetSolutionStepValue(VELOCITY_X,1);
	       temp_vec_np[1]=geom[1].FastGetSolutionStepValue(VELOCITY_X,1);
	       temp_vec_np[2]=geom[2].FastGetSolutionStepValue(VELOCITY_X,1);
	       temp_vec_np[3]=geom[3].FastGetSolutionStepValue(VELOCITY_X,1);

	fcomp = 0.25*(geom[0].FastGetSolutionStepValue(BODY_FORCE)[0] + geom[1].FastGetSolutionStepValue(BODY_FORCE)[0]+ geom[2].FastGetSolutionStepValue(BODY_FORCE)[0]+ geom[3].FastGetSolutionStepValue(BODY_FORCE)[0]);



		 noalias(LeftHandSideMatrix) = nu * prod(DN_DX, trans(DN_DX));
		        TempVecX=(-1.0) * prod(LeftHandSideMatrix, temp_vec_np);
		        //KRATOS_WATCH(TempVecX)
		        TempVecX[0] += DN_DX(0, 0) * p_avg;
		        TempVecX[1] += DN_DX(1, 0) * p_avg;
		        TempVecX[2] += DN_DX(2, 0) * p_avg;
		        TempVecX[3] += DN_DX(3, 0) * p_avg;

	TempVecX += rho * fcomp * N;
	TempVecX *=volume;

		        geom[0].FastGetSolutionStepValue(FORCE_X) += TempVecX[0];
	       geom[1].FastGetSolutionStepValue(FORCE_X) += TempVecX[1];
	       geom[2].FastGetSolutionStepValue(FORCE_X) += TempVecX[2];
	       geom[3].FastGetSolutionStepValue(FORCE_X) += TempVecX[3];


	fcomp = 0.25*(geom[0].FastGetSolutionStepValue(BODY_FORCE)[1] + geom[1].FastGetSolutionStepValue(BODY_FORCE)[1]+ geom[2].FastGetSolutionStepValue(BODY_FORCE)[1]+ geom[3].FastGetSolutionStepValue(BODY_FORCE)[1]);



	       temp_vec_np[0]=geom[0].FastGetSolutionStepValue(VELOCITY_Y,1);
	       temp_vec_np[1]=geom[1].FastGetSolutionStepValue(VELOCITY_Y,1);
	       temp_vec_np[2]=geom[2].FastGetSolutionStepValue(VELOCITY_Y,1);
	       temp_vec_np[3]=geom[3].FastGetSolutionStepValue(VELOCITY_Y,1);

	TempVecY=(-1.0) * prod(LeftHandSideMatrix, temp_vec_np);
		        TempVecY[0] += DN_DX(0, 1) * p_avg;
		        TempVecY[1] += DN_DX(1, 1) * p_avg;
		        TempVecY[2] += DN_DX(2, 1) * p_avg;
		        TempVecY[3] += DN_DX(3, 1) * p_avg;

		        TempVecY += density * fcomp * N;
		        TempVecY *=volume;


	geom[0].FastGetSolutionStepValue(FORCE_Y) += TempVecY[0];
	       geom[1].FastGetSolutionStepValue(FORCE_Y) += TempVecY[1];
	       geom[2].FastGetSolutionStepValue(FORCE_Y) += TempVecY[2];
	       geom[3].FastGetSolutionStepValue(FORCE_Y) += TempVecY[3];

	fcomp = 0.25*(geom[0].FastGetSolutionStepValue(BODY_FORCE)[2] + geom[1].FastGetSolutionStepValue(BODY_FORCE)[2]+ geom[2].FastGetSolutionStepValue(BODY_FORCE)[2]+ geom[3].FastGetSolutionStepValue(BODY_FORCE)[2]);

	       temp_vec_np[0]=geom[0].FastGetSolutionStepValue(VELOCITY_Z,1);
	       temp_vec_np[1]=geom[1].FastGetSolutionStepValue(VELOCITY_Z,1);
	       temp_vec_np[2]=geom[2].FastGetSolutionStepValue(VELOCITY_Z,1);
	       temp_vec_np[3]=geom[3].FastGetSolutionStepValue(VELOCITY_Z,1);

	TempVecZ=(-1.0) * prod(LeftHandSideMatrix, temp_vec_np);
	TempVecZ[0] += DN_DX(0, 2) * p_avg;
		        TempVecZ[1] += DN_DX(1, 2) * p_avg;
		        TempVecZ[2] += DN_DX(2, 2) * p_avg;
		        TempVecZ[3] += DN_DX(3, 2) * p_avg;

	TempVecZ += density * fcomp * N;
	TempVecZ *=volume;

	geom[0].FastGetSolutionStepValue(FORCE_Z) += TempVecZ[0];
	       geom[1].FastGetSolutionStepValue(FORCE_Z) += TempVecZ[1];
	       geom[2].FastGetSolutionStepValue(FORCE_Z) += TempVecZ[2];
	       geom[3].FastGetSolutionStepValue(FORCE_Z) += TempVecZ[3];*/
	}



	for (ModelPart::NodeIterator it = BaseType::GetModelPart().NodesBegin();it != BaseType::GetModelPart().NodesEnd(); ++it)
		{
		    //VELOCITY = VELOCITY + dt * Minv * AUX_VECTOR

	   double A = (it)->FastGetSolutionStepValue(NODAL_MASS);
	   if(A<0.0000000000000001){
	     A=1.0;
	   }
		    double dt_Minv = 1.0 / A;



		    array_1d<double,3>& temp = it->FastGetSolutionStepValue(FORCE);
		    it->FastGetSolutionStepValue(ANGULAR_ACCELERATION_X)= dt_Minv*temp[0];
		    it->FastGetSolutionStepValue(ANGULAR_ACCELERATION_Y)= dt_Minv*temp[1];
		    it->FastGetSolutionStepValue(ANGULAR_ACCELERATION_Z)= 0.0; //dt_Minv*temp[2];
	   
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

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "FluidPressureCalculateProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FluidPressureCalculateProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }


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
    ModelPart& mr_model_part;
    unsigned int mdomain_size;


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
//		FluidPressureCalculateProcess& operator=(HypoelasticStressCalculateProcess const& rOther);

    /// Copy constructor.
//		FluidPressureCalculateProcess(HypoelasticStressCalculateProcess const& rOther);


    ///@}

}; // Class FluidPressureCalculateProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  FluidPressureCalculateProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FluidPressureCalculateProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_FLUID_CALCULATE_PROCESS_INCLUDED  defined


