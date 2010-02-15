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
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2008-02-14 09:41:09 $
//   Revision:            $Revision: 1.2 $
//
// 


// System includes 


// External includes 


// Project includes 																							////perche' qui non ho external includes e i Project includes son diversi?
#include "includes/define.h"
#include "custom_conditions/no_slip_condition_2d.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	NoSlipCondition2D::NoSlipCondition2D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)																																	///
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	NoSlipCondition2D::NoSlipCondition2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{

	}

	Condition::Pointer NoSlipCondition2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new NoSlipCondition2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
KRATOS_WATCH("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
	}

	NoSlipCondition2D::~NoSlipCondition2D()													////e' il destructor questo?
	{
	}


	//************************************************************************************
	//************************************************************************************
	void NoSlipCondition2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		//calculation flags																								
		MatrixType temp = Matrix();
		CalculateLocalSystem(temp, rRightHandSideVector,  rCurrentProcessInfo);

	}

	//************************************************************************************
	//************************************************************************************
	void NoSlipCondition2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		//calculation flags
	  //		bool CalculateStiffnessMatrixFlag = true;
	  //bool CalculateResidualVectorFlag = true;
		//KRATOS_WATCH("@@@@@@@@ INSIDE CONDITION@@@@@@@@@@@@@@@@@@@");

        int nodes_number = 1;
        int dim = 2;
        unsigned int matsize = nodes_number * dim ;

        if (rLeftHandSideMatrix.size1() != matsize)
            rLeftHandSideMatrix.resize(matsize, matsize);

        if (rRightHandSideVector.size() != matsize)
            rRightHandSideVector.resize(matsize);


        noalias(rLeftHandSideMatrix) = ZeroMatrix(matsize, matsize);
        noalias(rRightHandSideVector) = ZeroVector(matsize);



// 			KRATOS_WATCH("@@@@@@@@finished total@@@@@@@@");
//    			KRATOS_WATCH("rRightHandSideVector");	
//    			KRATOS_WATCH(rRightHandSideVector);
			
			//	if(mean_ex_p !=GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_PRESSURE))
			//		mean_ex_p = 0.0;
			//KRATOS_WATCH(mean_ex_p);
			
/*			for(unsigned int ii = 0; ii< 2; ++ii)
				{

					int id = (2 + 1)*(ii);
					rRightHandSideVector[id] = mean_ex_p * An[0]* 0.5;
					rRightHandSideVector[id + 1] = mean_ex_p * An[1]* 0.5;
					rRightHandSideVector[id + 2] = 0.0;
				//KRATOS_WATCH(An);
				}*/
// 			KRATOS_WATCH(An);
//KRATOS_WATCH(p0);
//KRATOS_WATCH(p1);
//KRATOS_WATCH("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
 		//	KRATOS_WATCH(rRightHandSideVector);
		


		
	}

        //************************************************************************************
	//************************************************************************************

	void NoSlipCondition2D::CalculateLocalVelocityContribution(MatrixType& rDampMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

	int nodes_number = 1;
	int dim = 2;
	unsigned int matsize = nodes_number*(dim);

	if(rDampMatrix.size1() != matsize)
			rDampMatrix.resize(matsize,matsize,false); //false says not to preserve existing storage!!


	noalias(rDampMatrix) = ZeroMatrix(matsize,matsize); 


			//calculate normal to element:
			array_1d<double,2> area_normal;
                        area_normal = GetGeometry()[0].FastGetSolutionStepValue(NORMAL);
			array_1d<double,2> vel,mesh_vel;
			array_1d<double,2> normal;
			double area = area_normal[0]*area_normal[0] + area_normal[1]*area_normal[1];
			area = pow(area,0.5);
			normal[0] = area_normal[0] / area;
			normal[1] = area_normal[1] / area;

			//effective velocity:
			vel = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY); 
			
			mesh_vel = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);

			vel -= mesh_vel;
			double mod_vel = vel[0]*vel[0] + vel[1]*vel[1];
			mod_vel = pow(mod_vel,0.5);
			
			//calcolo la K moltiplicativa:
			double delta_t = rCurrentProcessInfo[DELTA_TIME];
			double v = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
			double K = 10.0*(mod_vel/area + v/(area*area));
			
    				
// noalias(rDampMatrix) += (K*area)*outer_prod(normal,normal);
			
			for(unsigned int i=0; i<2 ; i++)
			{
				for(unsigned int s=0; s<2 ; s++)

				    {
					rDampMatrix(i,s) += K*normal[i]*normal[s]*area;
				    }

			}

			for(unsigned int i=0; i<2 ; i++)

			{
				for(unsigned int q=0; q<2 ; q++)
					{
 					rRightHandSideVector[i] +=  rDampMatrix(i,q)*mesh_vel[q];
					}
			}

 		CalculateForce(area, rRightHandSideVector);

		KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************


      	void NoSlipCondition2D::CalculateForce(const double area, VectorType& rRightHandSideVector)
	{ 
		array_1d<double,2> vel;
		vel = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY) - GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);		 

		//determination of velocity's module and versor:
		double mod_vel = vel[0]*vel[0] + vel[1]*vel[1];
		mod_vel = sqrt(mod_vel);
		array_1d<double,2> vers_vel;
  
		vers_vel[0] = vel[0] / mod_vel;
		vers_vel[1] = vel[1] / mod_vel;

 		if (mod_vel == 0.0)
 		{
			vers_vel[0] = 0.0;
			vers_vel[1] = 0.0;
 		}

		//parameters:
		double k = 0.4;
		double a = 1.0/k;
		double B = 5.0;
		double v =  GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
		double density =  GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
		double mu = v * density;
		double toll = 0.0000001;
		double ym = 0.0825877;    //0.0093823   //ym respectively for u=0,2920 and u=2.57
		double y_mas_incercept = 10.9931899;
		double y_mas, y=1;
			
		//get mod_uthaw from the laminar law in order to see where our ym is (logarithmic or laminar region):
		double mod_uthaw = sqrt(mod_vel*v/ym);
		y_mas = ym*mod_uthaw/v;

 		if (y_mas>y_mas_incercept)
 		{
			//begin cicle to calculate the real u_thaw's module:
			int i=0;
 			while(fabs(y)>toll && i<100)           				 
 			{
 				y = mod_uthaw*(a*log(ym*mod_uthaw/v)+B)-mod_vel;
 				double y1 = a*log(ym*mod_uthaw/v)+B+a;
 				double mod_uthaw_new = mod_uthaw-y/y1;
 				mod_uthaw = mod_uthaw_new;
				i++;
 			 }

			if (fabs(y)>toll)
				KRATOS_WATCH("NOT CONVERGED IN NEWTON RAPSON LOOP!"); 
// 			KRATOS_WATCH("finished loop"); 
// 			KRATOS_WATCH(y)
// 			KRATOS_WATCH(fabs(y))
// 			KRATOS_WATCH(mod_uthaw)
// 			KRATOS_WATCH(vers_vel[0])
// 			KRATOS_WATCH(vers_vel[1])
// 			KRATOS_WATCH(area)


			//creating the rRightHandSideVector:
 			rRightHandSideVector[0] -= vers_vel[0]*area*mod_uthaw*mod_uthaw;
 			rRightHandSideVector[1] -= vers_vel[1]*area*mod_uthaw*mod_uthaw;
 		}	
 	        else
 			{
			rRightHandSideVector[0] -= vers_vel[0]*area*mu*mod_vel/ym;
			rRightHandSideVector[1] -= vers_vel[1]*area*mu*mod_vel/ym;
 			};

		//calculate pressure contributions:
		array_1d<double,2> area_normal;
                area_normal = GetGeometry()[0].FastGetSolutionStepValue(NORMAL);

		double p, p_x, p_y;
		p = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
		p_x = -area_normal[0]*p;
		p_y = -area_normal[1]*p;

		GetGeometry()[0].FastGetSolutionStepValue(REACTION_X) = -rRightHandSideVector[0] + p_x;
		GetGeometry()[0].FastGetSolutionStepValue(REACTION_Y) = -rRightHandSideVector[1] + p_y;
}
	//************************************************************************************
	//************************************************************************************

	void NoSlipCondition2D::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
										ProcessInfo& rCurrentProcessInfo,
										bool CalculateStiffnessMatrixFlag,
										bool CalculateResidualVectorFlag)
	{
		KRATOS_TRY

/*		unsigned int number_of_nodes = GetGeometry().size();
		unsigned int MatSize=number_of_nodes;

		//calculate lenght
		double x21 = GetGeometry()[1].X() - GetGeometry()[0].X();
		double y21 = GetGeometry()[1].Y() - GetGeometry()[0].Y();
		double lenght = x21*x21 + y21*y21;
		lenght = sqrt(lenght);

		double dt = rCurrentProcessInfo[DELTA_TIME];

		//proposal by riccardo
		double density = 0.5*(GetGeometry()[0].FastGetSolutionStepValue(DENSITY) + GetGeometry()[1].FastGetSolutionStepValue(DENSITY));

		double sound_speed = 10.0;
//		double K = density * sound_speed; //K = density * sound_speed*sound_speed

		//chapuza to make this of the same order of the other terms in the system independently on the deltatime
//		double K = density * 0.0001 / (dt*dt);  //.this would give a c = 1 for dt = 0.01...which we tested to be good

//		K = 1000;

	


		if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
		{
			if(rLeftHandSideMatrix.size1() != MatSize )
				rLeftHandSideMatrix.resize(MatSize,MatSize,false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize,MatSize);	
//			if(GetGeometry()[0].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1)
//				rLeftHandSideMatrix(0,0)  = 0.5 * lenght*lenght / (  K *dt);
//			if(GetGeometry()[1].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1)
//				rLeftHandSideMatrix(1,1) = 0.5 * lenght*lenght / (  K * dt);
			if(GetGeometry()[0].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1)
				rLeftHandSideMatrix(0,0)  = 0.5 * lenght  / ( density * sound_speed * dt );
			if(GetGeometry()[1].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1)
				rLeftHandSideMatrix(1,1) = 0.5 * lenght / (  density * sound_speed *dt );
		}

		if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
		{
			if(rRightHandSideVector.size() != MatSize )
				rRightHandSideVector.resize(MatSize,false);

//			noalias(rRightHandSideVector) = ZeroVector(2);

			array_1d<double,2> temp;
//			temp[0] = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1) - GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
//			temp[1] = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1) - GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
			temp[0] = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_OLD_IT) - GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
			temp[1] = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT) - GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
			noalias(rRightHandSideVector) =  prod( rLeftHandSideMatrix , temp);
		}
//		KRATOS_WATCH(rLeftHandSideMatrix);
//		KRATOS_WATCH(rRightHandSideVector);
*/
		KRATOS_CATCH("")
	}



	//************************************************************************************
	//************************************************************************************
	void NoSlipCondition2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dim = 2;
		unsigned int node_size = dim;

		
			if(rResult.size() != number_of_nodes*node_size)
				rResult.resize(number_of_nodes*node_size,false);	

			for (unsigned int i=0;i<number_of_nodes;i++)
			{
				rResult[i*node_size] = GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
				rResult[i*node_size+1] = GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();

			}
		KRATOS_CATCH("")
		
	}

	//************************************************************************************
	//************************************************************************************
	  void NoSlipCondition2D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dim = 2;
		unsigned int node_size = dim;

			if(ElementalDofList.size() != number_of_nodes*node_size)
				ElementalDofList.resize(number_of_nodes*node_size);	

			for (unsigned int i=0;i<number_of_nodes;i++)
			{
				ElementalDofList[i*node_size] = GetGeometry()[i].pGetDof(VELOCITY_X);
				ElementalDofList[i*node_size+1] = GetGeometry()[i].pGetDof(VELOCITY_Y);
			}
		KRATOS_CATCH("");
	}

    //*************************************************************************************
    //*************************************************************************************

    void NoSlipCondition2D::GetFirstDerivativesVector(Vector& values, int Step) {

        const unsigned int number_of_nodes = GetGeometry().size();

        unsigned int MatSize = 2 ;
        if (values.size() != MatSize) values.resize(MatSize, false);
        for (unsigned int i = 0; i < number_of_nodes; i++) {
            unsigned int index = i * 2;
            values[index] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_X, Step);
            values[index + 1] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_Y, Step);
            //values[index + 2] = GetGeometry()[i].GetSolutionStepValue(PRESSURE, Step);

        }


    }
    //************************************************************************************
    //************************************************************************************

    void NoSlipCondition2D::GetSecondDerivativesVector(Vector& values, int Step) {
        const unsigned int number_of_nodes = GetGeometry().size();
        unsigned int MatSize =2;
        if (values.size() != MatSize) values.resize(MatSize, false);
        for (unsigned int i = 0; i < number_of_nodes; i++) {
            unsigned int index = i * 2;
            values[index] = 0.0;
            values[index + 1] = 0.0;
            //values[index + 2] = 0.0;
        }
    }
    //************************************************************************************
    //************************************************************************************	  
} // Namespace Kratos


