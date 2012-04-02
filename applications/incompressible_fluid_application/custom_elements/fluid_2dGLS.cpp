/*
==============================================================================
KratosR1IncompressibleFluidApplication 
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
//   Last modified by:    $Author: jmarti $
//   Date:                $Date: 2008-11-26 08:17:41 $
//   Revision:            $Revision: 1.0 $
//
//
 
//#define GRADPN_FORM

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/fluid_2dGLS.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{

    //************************************************************************************
    //************************************************************************************

    Fluid2DGLS::Fluid2DGLS(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    Fluid2DGLS::Fluid2DGLS(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {
    }

    Element::Pointer Fluid2DGLS::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        KRATOS_TRY

        return Element::Pointer(new Fluid2DGLS(NewId, GetGeometry().Create(ThisNodes), pProperties));
        KRATOS_CATCH("");
    }

    Fluid2DGLS::~Fluid2DGLS()
    {
    }

    //************************************************************************************
    //************************************************************************************

    void Fluid2DGLS::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

                int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];

        if (FractionalStepNumber < 3) //first step of the fractional step solution
        {
//            int ComponentIndex = FractionalStepNumber - 1;

        } else if (FractionalStepNumber == 4)//second step of the fractional step solution
        {
            Stage2(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
        }

        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    void Fluid2DGLS::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {

		KRATOS_TRY
		
		///////////////////////NECESSARY LOCALS///////////////////////////////////////////
		boost::numeric::ublas::bounded_matrix<double,3,3> msWorkMatrix = ZeroMatrix(3,3);
		boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX = ZeroMatrix(3,2);
		array_1d<double,3> msN = ZeroVector(3); //dimension = number of nodes
		boost::numeric::ublas::bounded_matrix<double,6,2> msShapeFunc = ZeroMatrix(6,2);
		boost::numeric::ublas::bounded_matrix<double,2,6> msConvOp = ZeroMatrix(2,6);
	       	boost::numeric::ublas::bounded_matrix<double,6,6> msAuxMat = ZeroMatrix(6,6);
		array_1d<double,6> msAuxVec = ZeroVector(6); //dimension = number of nodes
		array_1d<double,2> ms_adv_vel = ZeroVector(2); //dimesion coincides with space dimension
		array_1d<double,2> ms_vel_gauss = ZeroVector(2); //dimesion coincides with space dimension
		array_1d<double,3> ms_temp_vec_np = ZeroVector(3); //dimension = number of nodes
		array_1d<double,3> ms_aux0 = ZeroVector(3); //dimension = number of nodes
		array_1d<double,3> ms_aux1 = ZeroVector(3); //dimension = number of nodes
		boost::numeric::ublas::bounded_matrix<double,6,6> msAuxMat1 = ZeroMatrix(6,6);
    		boost::numeric::ublas::bounded_matrix<double,6,3> msAuxMat2 = ZeroMatrix(6,3);
		boost::numeric::ublas::bounded_matrix<double,2,2> msGrad_ug = ZeroMatrix(2,2);
		array_1d<double,6> msStabMomRes = ZeroVector(6); //dimension = number of nodes
		boost::numeric::ublas::bounded_matrix<double,6,3> msGradOp = ZeroMatrix(6,3);
		
		boost::numeric::ublas::bounded_matrix<double,3,3> msWorkMatrix1 = ZeroMatrix(3,3);
		///////////////////////////////////////////////////////////////////////////////////
		
		if(rRightHandSideVector.size() != 3)
		{
			rRightHandSideVector.resize(3,false);
		}

		double dt = rCurrentProcessInfo[DELTA_TIME];
		
		//fract. vel, that is calculated in the first Fractional Step.. but is saved inside the "VELOCITY" VARIABLE
		//so, u_n os VELOCITY, 1 and u_n-1 VELOCITY,2 
		const array_1d<double,3>& fv0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
		const array_1d<double,3>& fv0_old = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,1);
		//const array_1d<double,3>& fv0_n_1 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,2);
		const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
		const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
		double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
		//double p0_old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_OLD_IT);
		double p0_old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1); 	 	 	 	
		const array_1d<double,3>& ff0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);	

		const array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
		const array_1d<double,3>& fv1_old = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,1);
		//const array_1d<double,3>& fv1_n_1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,2);
		const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
		const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
		double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE); 
		//double p1_old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT); 	 	
		double p1_old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1); 	 	 	 	
		const array_1d<double,3>& ff1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);	

		const array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);	
		const array_1d<double,3>& fv2_old = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,1);
		//const array_1d<double,3>& fv2_n_1 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,2);
		const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
		const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
		double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE); 
		double p2_old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1); 	 	 	 	
		//old iteration can be used if we want to iterate between 1st and 2nd fractional steps
		//double p2_old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT); 	 	
		const array_1d<double,3>& ff2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);	


		//in msAuxVec we store the velocity, (not the frac step vel, but u_n, the one that enters the stabilization)
		msAuxVec[0]=fv0[0];
		msAuxVec[1]=fv0[1];
		msAuxVec[2]=fv1[0];
		msAuxVec[3]=fv1[1];
		msAuxVec[4]=fv2[0];
		msAuxVec[5]=fv2[1];

		//getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Area);

		//calculating average density and viscosity
		double nu = 0.33333333333333*(nu0 + nu1 + nu2 );
 		double density = 0.33333333333333*(rho0 + rho1 + rho2 );

		//ms_vel_gauss[i] =  msN[0]*(fv0[i]) + msN[1]*(fv1[i]) +  msN[2]*(fv2[i]);
		//but with one integration N=0.333333333
		ms_vel_gauss[0] =  0.33333333333333*(fv0[0]+fv1[0]+fv2[0]);
		ms_vel_gauss[1] =  0.33333333333333*(fv0[1]+fv1[1]+fv2[1]);

		//calculating parameter tau (saved internally to each element)
		double h = sqrt(2.00*Area);
		double norm_u = ms_vel_gauss[0]*ms_vel_gauss[0] + ms_vel_gauss[1]*ms_vel_gauss[1];
		norm_u = sqrt(norm_u);
		double tau = 1.00 / ( 4.00*nu/(h*h) +2.00*norm_u/h + 1.0/dt);
		//tau=0.0;		
		noalias(msWorkMatrix)=prod(msDN_DX,trans(msDN_DX));
		noalias(msWorkMatrix1) = (dt/2.0 + tau) * Area*msWorkMatrix;


		//////////////////////////////////////////////////////////
		////////////		AND NOW RHS	//////////////////
		//////////////////////////////////////////////////////////	
	
		//Dirichlet contribution  (that is: LHS*p_new)
		ms_temp_vec_np[0] = p0; 
		ms_temp_vec_np[1] = p1; 
		ms_temp_vec_np[2] = p2; 
		//LHS is already multiplied by AREA
		noalias(rRightHandSideVector) = -prod(msWorkMatrix1,ms_temp_vec_np);
		
		//NOW RHS-=dt L p_old		
		//changing the meaning of temp_vec_np
		ms_temp_vec_np[0] = p0_old; 
		ms_temp_vec_np[1] = p1_old; 
		ms_temp_vec_np[2] = p2_old; 

		noalias(rRightHandSideVector) += Area* dt/2.0 * (prod(msWorkMatrix,ms_temp_vec_np)) ;
		
		//***************************************************************************
		
		//here we have the Du_tila term
		double Gaux;
		Gaux =  msDN_DX(0,0)*fv0[0] + msDN_DX(0,1)*fv0[1];
		Gaux += msDN_DX(1,0)*fv1[0] + msDN_DX(1,1)*fv1[1];
		Gaux += msDN_DX(2,0)*fv2[0] + msDN_DX(2,1)*fv2[1];
		rRightHandSideVector[0] -= density*Area*Gaux * msN[0]; 
		rRightHandSideVector[1] -= density*Area*Gaux * msN[1]; 
		rRightHandSideVector[2] -= density*Area*Gaux * msN[2]; 
		
		ms_aux0=0.33333333333333333*(ff0+ff1+ff2);
		//ms_aux1 - is the product of: (nabla q, f)
		ms_aux1[0]=msDN_DX(0,0)*ms_aux0[0]+msDN_DX(0,1)*ms_aux0[1];
		ms_aux1[1]=msDN_DX(1,0)*ms_aux0[0]+msDN_DX(1,1)*ms_aux0[1];
		ms_aux1[2]=msDN_DX(2,0)*ms_aux0[0]+msDN_DX(2,1)*ms_aux0[1];
		//KRATOS_WATCH(temp)
		//rRightHandSideVector += tau*density*Area*ms_aux1;
		
		
		//RHS += -tau*nablaN*du_gausspoint/dt
			
		//we reuse ms_vel_gauss to store the accelerations( (u_n - u_n-1)/dt)		
		
		ms_vel_gauss[0]=0.33333333333*(fv0[0]+fv1[0]+fv2[0]-fv0_old[0]-fv1_old[0]-fv2_old[0])/dt;
		ms_vel_gauss[1]=0.33333333333*(fv0[1]+fv1[1]+fv2[1]-fv0_old[1]-fv1_old[1]-fv2_old[1])/dt;
		
		//and now we reuse ms_aux1

		ms_aux1=prod(msDN_DX,ms_vel_gauss);
		
		noalias(rRightHandSideVector) -= tau*density*Area*ms_aux1;
		
 
		//and now the stabilization term referring to the convective operator
		//RHS+=nablaq*convop (convetcion of the Gauss point vel)
		
		//contains d_ug/dx
		//x-component of u derived with resp to x, remember msAuxVec stores velocity
		msGrad_ug(0,0)=msDN_DX(0,0)*msAuxVec[0]+msDN_DX(1,0)*msAuxVec[2]+msDN_DX(2,0)*msAuxVec[4];
		//x-component of u derived with resp to y
		msGrad_ug(0,1)=msDN_DX(0,1)*msAuxVec[0]+msDN_DX(1,1)*msAuxVec[2]+msDN_DX(2,1)*msAuxVec[4];
		//y-component of u derived with resp to x
		msGrad_ug(1,0)=msDN_DX(0,0)*msAuxVec[1]+msDN_DX(1,0)*msAuxVec[3]+msDN_DX(2,0)*msAuxVec[5];
		//y-component of u derived with resp to y
		msGrad_ug(1,1)=msDN_DX(0,1)*msAuxVec[1]+msDN_DX(1,1)*msAuxVec[3]+msDN_DX(2,1)*msAuxVec[5];
		
		array_1d<double,2> a;
		a[0]=0.33333333333333*(msAuxVec[0]+msAuxVec[2]+msAuxVec[4])*msGrad_ug(0,0)+0.33333333333333*(msAuxVec[1]+msAuxVec[3]+msAuxVec[5])*msGrad_ug(0,1);
		a[1]=0.33333333333333*(msAuxVec[0]+msAuxVec[2]+msAuxVec[4])*msGrad_ug(1,0)+0.33333333333333*(msAuxVec[1]+msAuxVec[3]+msAuxVec[5])*msGrad_ug(1,1);
		
		//we again reuse ms_aux0
		//noalias(ms_aux0) = prod(msDN_DX,a);
		//noalias(rRightHandSideVector) -double rho_gauss=0.0;


		

		
		KRATOS_CATCH("")


    }

    //************************************************************************************
    //************************************************************************************
    //calculation by component of the fractional step velocity corresponding to the first stage

    void Fluid2DGLS::Stage2(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
		KRATOS_TRY
		
		///////////////////////NECESSARY LOCALS///////////////////////////////////////////
		boost::numeric::ublas::bounded_matrix<double,3,3> msWorkMatrix = ZeroMatrix(3,3);
		boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX = ZeroMatrix(3,2);
		array_1d<double,3> msN = ZeroVector(3); //dimension = number of nodes
		boost::numeric::ublas::bounded_matrix<double,6,2> msShapeFunc = ZeroMatrix(6,2);
		boost::numeric::ublas::bounded_matrix<double,2,6> msConvOp = ZeroMatrix(2,6);
	       	boost::numeric::ublas::bounded_matrix<double,6,6> msAuxMat = ZeroMatrix(6,6);
		array_1d<double,6> msAuxVec = ZeroVector(6); //dimension = number of nodes
		array_1d<double,2> ms_adv_vel = ZeroVector(2); //dimesion coincides with space dimension
		array_1d<double,2> ms_vel_gauss = ZeroVector(2); //dimesion coincides with space dimension
		array_1d<double,3> ms_temp_vec_np = ZeroVector(3); //dimension = number of nodes
		array_1d<double,3> ms_aux0 = ZeroVector(3); //dimension = number of nodes
		array_1d<double,3> ms_aux1 = ZeroVector(3); //dimension = number of nodes
		boost::numeric::ublas::bounded_matrix<double,6,6> msAuxMat1 = ZeroMatrix(6,6);
    		boost::numeric::ublas::bounded_matrix<double,6,3> msAuxMat2 = ZeroMatrix(6,3);
		boost::numeric::ublas::bounded_matrix<double,2,2> msGrad_ug = ZeroMatrix(2,2);
		array_1d<double,6> msStabMomRes = ZeroVector(6); //dimension = number of nodes
		boost::numeric::ublas::bounded_matrix<double,6,3> msGradOp = ZeroMatrix(6,3);
		array_1d<double, 2 > vel_gauss;
		boost::numeric::ublas::bounded_matrix<double,3,3> Tres = ZeroMatrix(3,3);
		boost::numeric::ublas::bounded_matrix<double,3,3> msResta = IdentityMatrix(3,3);

		unsigned int TDim = 2;

 		//const double Cs = this->GetValue(C_SMAGORINSKY);
	

		//KRATOS_WATCH(Cs);

		//KRATOS_ERROR(std::logic_error, "method not implemented", "");

		mThisIntegrationMethod= GeometryData::GI_GAUSS_1;



		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

		const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);//GetGeometry().CalculateShapeFunctionsIntegrationPointsValues(GeometryData::GI_GAUSS_1);
		///////////////////////////////////////////////////////////////////////////////////
		
		if(rRightHandSideVector.size() != 3)
		{
			rLeftHandSideMatrix.resize(3,3,false);
			rRightHandSideVector.resize(3,false);
		}

		double dt = rCurrentProcessInfo[DELTA_TIME];
		
		//fract. vel, that is calculated in the first Fractional Step.. but is saved inside the "VELOCITY" VARIABLE
		//so, u_n os VELOCITY, 1 and u_n-1 VELOCITY,2 
		const array_1d<double,3>& fv0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
		const array_1d<double,3>& fv0_old = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,1);
		//const array_1d<double,3>& fv0_n_1 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,2);
		const array_1d<double, 3 > & proj0 = GetGeometry()[0].FastGetSolutionStepValue(PRESS_PROJ);
		const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
		const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
		double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
		//double p0_old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_OLD_IT);
		double p0_old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1); 	 	 	 	
		const array_1d<double,3>& ff0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);	

		const array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
		const array_1d<double,3>& fv1_old = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,1);
		//const array_1d<double,3>& fv1_n_1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,2);
		const array_1d<double, 3 > & proj1 = GetGeometry()[1].FastGetSolutionStepValue(PRESS_PROJ);
		const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
		const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
		double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE); 
		//double p1_old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT); 	 	
		double p1_old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1); 	 	 	 	
		const array_1d<double,3>& ff1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);	

		const array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);	
		const array_1d<double,3>& fv2_old = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,1);
		//const array_1d<double,3>& fv2_n_1 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,2);
		const array_1d<double, 3 > & proj2 = GetGeometry()[2].FastGetSolutionStepValue(PRESS_PROJ);
		const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
		const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
		double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE); 
		double p2_old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1); 	 	 	 	
		//old iteration can be used if we want to iterate between 1st and 2nd fractional steps
		//double p2_old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT); 	 	
		const array_1d<double,3>& ff2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);	

		int counti=0.0;
//		double density = 0.0;

		const unsigned int number_of_points = GetGeometry().size();
		
		
		//KRATOS_WATCH(counti);

		//in msAuxVec we store the velocity, (not the frac step vel, but u_n, the one that enters the stabilization)
		msAuxVec[0]=fv0[0];
		msAuxVec[1]=fv0[1];
		msAuxVec[2]=fv1[0];
		msAuxVec[3]=fv1[1];
		msAuxVec[4]=fv2[0];
		msAuxVec[5]=fv2[1];

		//getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Area);

		//KRATOS_WATCH("AREA");
		//KRATOS_WATCH(Area);
		//calculating average density and viscosity
		double nu = 0.33333333333333*(nu0 + nu1 + nu2 );
 		double density = 0.33333333333333*(rho0 + rho1 + rho2 );

		

		//ms_vel_gauss[i] =  msN[0]*(fv0[i]) + msN[1]*(fv1[i]) +  msN[2]*(fv2[i]);
		//but with one integration N=0.333333333
		ms_vel_gauss[0] =  0.33333333333333*(fv0[0]+fv1[0]+fv2[0]);
		ms_vel_gauss[1] =  0.33333333333333*(fv0[1]+fv1[1]+fv2[1]);

		//calculating parameter tau (saved internally to each element)
		double h = sqrt(2.00*Area);
		double norm_u = ms_vel_gauss[0]*ms_vel_gauss[0] + ms_vel_gauss[1]*ms_vel_gauss[1];
		norm_u = sqrt(norm_u);


		double rho_gauss=0.0;


		int countip=0.0;
		int countipp=0.0;


		//for(unsigned int i = 0; i < 3.0 ; i++) 
		//{
		//	if(GetGeometry()[i].FastGetSolutionStepValue(FLAG_VARIABLE)==1.0) counti++;
		//}

		//for(unsigned int i = 0; i < 3.0 ; i++) 
		//{
		//	if(GetGeometry()[i].FastGetSolutionStepValue(FLAG_VARIABLE)==10.0) countip++;
		//}

		//if(countip==3.0)	{density=10.0;} 
		//else density=1.0;


		//for(unsigned int i = 0; i < 3.0 ; i++) 
		//{
		//	if(GetGeometry()[i].FastGetSolutionStepValue(FLAG_VARIABLE)==6.0) countipp++;
		//}

		//if(countipp==3.0) {density=6.0;}

		//if(countip==3.0)	{density=1000.0;} 
		//else density=1.0;


		
		//tau=0.0;
		
		//AND NOW WE ADD THE RESPECTIVE CONTRIBUTIONS TO THE RHS AND LHS of THE SECOND FRAC STEP
		//we use Backward Euler for this step, therefore stab. contribution no RHS +=Tau1*(gradQ, residual)
		//								   and LHS +=Tau1*(gradQ, gradP)
		//laplacian term	       L = Dt * gradN * trans(gradN);
		//stabilization term       Spp = tau * gradN * trans(gradN);
		//WATCH OUT for DIVISION with RHO - check if it changes or not in case of Momentum being the primary Variable
		//
		//	msWorkMatrix stores the element laplacian
		//
		//noalias(msWorkMatrix)=prod(msDN_DX,trans(msDN_DX));
		//noalias(rLeftHandSideMatrix) = (dt/(2.0*density) + tau/density) * Area*msWorkMatrix;
		//rhs consists of D*u_tilda (divergence of the Fractional velocity) and the term: Tau1*(nabla_q, residual_gausspoint)
		//fv is u_tilda
		
		//////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////	
	
		//Dirichlet contribution  (that is: LHS*p_new)
		//ms_temp_vec_np[0] = p0; 
		//ms_temp_vec_np[1] = p1; 
		//ms_temp_vec_np[2] = p2; 
		//LHS is already multiplied by AREA
		//noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix,ms_temp_vec_np);
		
		//NOW RHS-=dt L p_old		
		//changing the meaning of temp_vec_np
		//ms_temp_vec_np[0] = p0_old; 
		//ms_temp_vec_np[1] = p1_old; 
		//ms_temp_vec_np[2] = p2_old; 

		//noalias(rRightHandSideVector) += Area* dt/(2.0 * density ) * (prod(msWorkMatrix,ms_temp_vec_np)) ;
		
		//***************************************************************************
		




		/*vel_gauss[0] = msN[0] * proj0[0] + msN[1] * proj1[0] + msN[2] * proj2[0];
        	vel_gauss[1] = msN[0] * proj0[1] + msN[1] * proj1[1] + msN[2] * proj2[1];
        	vel_gauss *= tau;
        	noalias(rRightHandSideVector) += prod(msDN_DX, vel_gauss)*Area;*/// *density;

		
		/*if(countip>0.0) {density=10.0;}
		if(counti>0.0) {density=1.0;}*/

		//density = 0.33333333333333*(rho0 + rho1 + rho2 );


		//double Gaux;
		//Gaux =  msDN_DX(0,0)*fv0[0] + msDN_DX(0,1)*fv0[1];
		//Gaux += msDN_DX(1,0)*fv1[0] + msDN_DX(1,1)*fv1[1];
		//Gaux += msDN_DX(2,0)*fv2[0] + msDN_DX(2,1)*fv2[1];
		//rRightHandSideVector[0] -= Area*Gaux * msN[0]; 
		//rRightHandSideVector[1] -= Area*Gaux * msN[1]; 
		//rRightHandSideVector[2] -= Area*Gaux * msN[2]; 

		//vel_gauss[0] = msN[0] * proj0[0] + msN[1] * proj1[0] + msN[2] * proj2[0];
       		//vel_gauss[1] = msN[0] * proj0[1] + msN[1] * proj1[1] + msN[2] * proj2[1];
       		//vel_gauss *= tau;
       		//noalias(rRightHandSideVector) += prod(msDN_DX, vel_gauss) * Area ;



		///////////////////////////////////////
		///////////////////////////////////////
		//mInvJ0.resize(integration_points.size());
		//mDetJ0.resize(integration_points.size(),false);
		//KRATOS_WATCH("ACAAAAAAAAAAA");
		double total_density=0.0;
 		//const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

		//GeometryType::JacobiansType J0;
		//J0 = GetGeometry().Jacobian(J0, mThisIntegrationMethod);  


		//density = this->GetValue(DENSITY);
		//KRATOS_WATCH(density);

		density = this->GetValue(DENSITY);
		//density = 0.33333333333333*(rho0 + rho1 + rho2 );
		//if(density>1.01 and density<999.0) density=500.0;
		/*noalias(rLeftHandSideMatrix) = ZeroMatrix(TDim + 1, TDim + 1);
		noalias(rRightHandSideVector) =  ZeroVector(TDim + 1);*/
		double tau = 1.00 / ( 4.00* nu/(density * h*h) +  2.00*norm_u/h +  1.0/dt);


		noalias(msWorkMatrix)=prod(msDN_DX,trans(msDN_DX));
		noalias(rLeftHandSideMatrix) = (dt/(2.0*density) + tau/density) * msWorkMatrix * Area;

		ms_temp_vec_np[0] = p0; 
		ms_temp_vec_np[1] = p1; 
		ms_temp_vec_np[2] = p2; 
		//LHS is already multiplied by AREA

		noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix,ms_temp_vec_np);


		ms_temp_vec_np[0] = p0_old; 
		ms_temp_vec_np[1] = p1_old; 
		ms_temp_vec_np[2] = p2_old; 

		noalias(rRightHandSideVector) +=  dt/(2.0 * density ) * (prod(msWorkMatrix,ms_temp_vec_np)) * Area;


		double nodal_1 = 0.333333333333 * Area * density ;
		double nodal_2 = 0.333333333333 * Area * density ;
		double nodal_3 = 0.333333333333 * Area * density ;

        	GetGeometry()[0].SetLock();
		GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA) += nodal_1;
		GetGeometry()[0].UnSetLock();

		GetGeometry()[1].SetLock();
		GetGeometry()[1].FastGetSolutionStepValue(NODAL_AREA) += nodal_2;
		GetGeometry()[1].UnSetLock();

		GetGeometry()[2].SetLock();
		GetGeometry()[2].FastGetSolutionStepValue(NODAL_AREA) += nodal_3;
		GetGeometry()[2].UnSetLock();


		double Gaux;
		Gaux =  msDN_DX(0,0)*fv0[0] + msDN_DX(0,1)*fv0[1];
		Gaux += msDN_DX(1,0)*fv1[0] + msDN_DX(1,1)*fv1[1];
		Gaux += msDN_DX(2,0)*fv2[0] + msDN_DX(2,1)*fv2[1];
		rRightHandSideVector[0] -= Area * Gaux * msN[0]; 
		rRightHandSideVector[1] -= Area * Gaux * msN[1]; 
		rRightHandSideVector[2] -= Area * Gaux * msN[2];



		vel_gauss[0] = msN[0] * proj0[0] + msN[1] * proj1[0] + msN[2] * proj2[0];
        	vel_gauss[1] = msN[0] * proj0[1] + msN[1] * proj1[1] + msN[2] * proj2[1];
        	vel_gauss *= tau;
        	noalias(rRightHandSideVector) += prod(msDN_DX, vel_gauss) * Area;


				
			


		////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////

/*		for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
		{

			//KRATOS_WATCH(integration_points.size());
			//KRATOS_WATCH(PointNumber);
			unsigned int number_of_nodes = GetGeometry().PointsNumber();
			unsigned int dimension = GetGeometry().WorkingSpaceDimension();
			const Vector& N=row(Ncontainer,PointNumber);
			
			MathUtils<double>::InvertMatrix(J0[PointNumber],mInvJ0[PointNumber],mDetJ0[PointNumber]);
			double Weight = integration_points[PointNumber].Weight()* mDetJ0[PointNumber];
				
			noalias(msDN_DX) = prod(DN_De[PointNumber],mInvJ0[PointNumber]);
			//KRATOS_WATCH(msDN_DX);
			//double rho_gauss=0.0;
			//density =0.0;
			//for (unsigned int i=0;i<number_of_nodes;i++) density += N[i]*GetGeometry()[i].FastGetSolutionStepValue(DENSITY);
			density = N[0]*GetGeometry()[0].FastGetSolutionStepValue(DENSITY) + N[1]*GetGeometry()[1].FastGetSolutionStepValue(DENSITY) +N[2]*GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
			density = this->GetValue(DENSITY);	
		
			double tau = 1.00 / ( 4.00* nu/(density * h*h) +  2.00*norm_u/h +  1.0/dt);
			noalias(msWorkMatrix)=prod(msDN_DX,trans(msDN_DX));


			noalias(rLeftHandSideMatrix) += (dt/(2.0*density) + tau/density) * msWorkMatrix * Weight;
		
			//NOW RHS-=dt L p_old		
			//changing the meaning of temp_vec_np
			ms_temp_vec_np[0] = p0_old; 
			ms_temp_vec_np[1] = p1_old; 
			ms_temp_vec_np[2] = p2_old; 

			noalias(rRightHandSideVector) +=  dt/(2.0 * density ) * (prod(msWorkMatrix,ms_temp_vec_np)) * Weight;


			Tres(0,0)=N[0]*N[0];
			Tres(0,1)=N[0]*N[1];
			Tres(0,2)=N[0]*N[2];
			Tres(1,0)=N[1]*N[0];
			Tres(1,1)=N[1]*N[1];
			Tres(1,2)=N[1]*N[2];
			Tres(2,0)=N[2]*N[0];
			Tres(2,1)=N[2]*N[1];
			Tres(2,2)=N[2]*N[2];
			//KRATOS_WATCH(Tres);

			double nodal_1 = (N[0]*N[0] + N[0]*N[1] + N[0]*N[2]) * Weight * density ;
			double nodal_2 = (N[1]*N[0] + N[1]*N[1] + N[1]*N[2]) * Weight * density ;
			double nodal_3 = (N[2]*N[0] + N[2]*N[1] + N[2]*N[2]) * Weight * density ;

        		GetGeometry()[0].SetLock();
			GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA) += nodal_1;
			GetGeometry()[0].UnSetLock();

			GetGeometry()[1].SetLock();
			GetGeometry()[1].FastGetSolutionStepValue(NODAL_AREA) += nodal_2;
			GetGeometry()[1].UnSetLock();

			GetGeometry()[2].SetLock();
			GetGeometry()[2].FastGetSolutionStepValue(NODAL_AREA) += nodal_3;
			GetGeometry()[2].UnSetLock();

			
			//density=1.0;
			double Gaux;
			Gaux =  msDN_DX(0,0)*fv0[0] + msDN_DX(0,1)*fv0[1];
			Gaux += msDN_DX(1,0)*fv1[0] + msDN_DX(1,1)*fv1[1];
			Gaux += msDN_DX(2,0)*fv2[0] + msDN_DX(2,1)*fv2[1];
			rRightHandSideVector[0] -= Weight * Gaux * N[0]; 
			rRightHandSideVector[1] -= Weight * Gaux * N[1]; 
			rRightHandSideVector[2] -= Weight * Gaux * N[2];



			vel_gauss[0] = N[0] * proj0[0] + N[1] * proj1[0] + N[2] * proj2[0];
        		vel_gauss[1] = N[0] * proj0[1] + N[1] * proj1[1] + N[2] * proj2[1];
        		vel_gauss *= tau;
        		noalias(rRightHandSideVector) += prod(msDN_DX, vel_gauss) * Weight;


				
		}	

		//Dirichlet contribution  (that is: LHS*p_new)
		ms_temp_vec_np[0] = -p0; 
		ms_temp_vec_np[1] = -p1; 
		ms_temp_vec_np[2] = -p2; 
		//LHS is already multiplied by AREA
		noalias(rRightHandSideVector) += prod(rLeftHandSideMatrix,ms_temp_vec_np);

*/		
		KRATOS_CATCH("")


    }

    //************************************************************************************
    //************************************************************************************
    // this subroutine calculates the nodal contributions for the explicit steps of the
    // fractional step procedure

    void Fluid2DGLS::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];
        boost::numeric::ublas::bounded_matrix<double, 3, 3 > MassFactors = 1.0 / 3.0 * IdentityMatrix(3, 3);
        boost::numeric::ublas::bounded_matrix<double, 3, 3 > WorkMatrix;
        boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX;
        array_1d<double, 3 > N;
        array_1d<double, 2 > vel_gauss;
        array_1d<double, 3 > temp_vec_np;
        array_1d<double, 3 > u_DN;
        array_1d<double, 3 > aux0;
        array_1d<double, 3 > aux1;
        array_1d<double, 3 > aux2;

        //getting data for the given geometry
        double Area;
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);
	//KRATOS_ERROR(std::logic_error, "method not implemented", "");
        if (FractionalStepNumber == 5) //calculation of stabilization terms
        {
		//KRATOS_ERROR(std::logic_error, "method not implemented", "");
		///////////////////////NECESSARY LOCALS///////////////////////////////////////////
		boost::numeric::ublas::bounded_matrix<double,3,3> msWorkMatrix = ZeroMatrix(3,3);
		array_1d<double,6> GalerkinRHS = ZeroVector(6); //dimension = number of nodes
		boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX = ZeroMatrix(3,2);
		array_1d<double,3> msN = ZeroVector(3); //dimension = number of nodes
		boost::numeric::ublas::bounded_matrix<double,6,2> msShapeFunc = ZeroMatrix(6,2);
		boost::numeric::ublas::bounded_matrix<double,2,6> msConvOp = ZeroMatrix(2,6);
	       	boost::numeric::ublas::bounded_matrix<double,6,6> msAuxMat = ZeroMatrix(6,6);
		array_1d<double,6> msAuxVec = ZeroVector(6); //dimension = number of nodes
		array_1d<double,2> ms_adv_vel = ZeroVector(2); //dimesion coincides with space dimension
		array_1d<double,2> ms_vel_gauss = ZeroVector(2); //dimesion coincides with space dimension

		///////////////////////////////////////////////////////////////////////////////////

	
		//first we compute  the force term and pressure gradient terms:
		//getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Area);

		//getting the velocity on the nodes and other necessary variabless
		const array_1d<double,3> vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,1);
		double p_n0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1);
		const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
		const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
		
		const array_1d<double,3> vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,1);
		double p_n1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1);
		const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
		const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
		
		const array_1d<double,3>& vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,1);
		double p_n2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1);
		const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
		const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
			


		//====================================================================
		//calculating viscosity and density
		double nu = 0.333333333333333333333333*(nu0 + nu1 + nu2 );
		double density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );

		int counti=0.0;
		int countip=0.0;
		int countipp=0.0;



		//for(unsigned int i = 0; i < 3.0 ; i++) 
		//{
		//	if(GetGeometry()[i].FastGetSolutionStepValue(FLAG_VARIABLE)==1.0) counti++;
		//}

		for(unsigned int i = 0; i < 3.0 ; i++) 
		{
			if(GetGeometry()[i].FastGetSolutionStepValue(FLAG_VARIABLE)==10.0) countip++;
		}



		//for(unsigned int i = 0; i < 3.0 ; i++) 
		//{
		//	if(GetGeometry()[i].FastGetSolutionStepValue(FLAG_VARIABLE)==6.0) countipp++;
		//}

		//if(countipp==3.0) {density=6.0;}

		//if(countip==3.0)	{density=1000.0;} 
		//else density=1.0;

/*		if(countip>0.0) {density=10.0;}
		if(counti>0.0) {density=1.0;}*/


		//density = 0.33333333333333*(rho0 + rho1 + rho2 );
		//density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );
		
		//density = this->GetValue(DENSITY);
		//KRATOS_WATCH(density);
		
		density = this->GetValue(DENSITY);
		//density = 0.33333333333333*(rho0 + rho1 + rho2 );
		//if(density>1.01 and density<999.0) density=500.0;
		//density = 0.33333333333333*(rho0 + rho1 + rho2 );
		//if(density>1.01 and density<900.0) density=500.0;
		//KRATOS_WATCH(density);
		noalias(msWorkMatrix) = Area* /*density **/ nu * prod(msDN_DX,trans(msDN_DX));
				
		//x comp
		GalerkinRHS[0]=-1.0*(msWorkMatrix(0,0)*vel0[0]+msWorkMatrix(0,1)*vel1[0]+msWorkMatrix(0,2)*vel2[0]);
		//y comp
		GalerkinRHS[1]=-1.0*(msWorkMatrix(0,0)*vel0[1]+msWorkMatrix(0,1)*vel1[1]+msWorkMatrix(0,2)*vel2[1]);

		//x comp
		GalerkinRHS[2]=-1.0*(msWorkMatrix(1,0)*vel0[0]+msWorkMatrix(1,1)*vel1[0]+msWorkMatrix(1,2)*vel2[0]);
		//y comp
		GalerkinRHS[3]=-1.0*(msWorkMatrix(1,0)*vel0[1]+msWorkMatrix(1,1)*vel1[1]+msWorkMatrix(1,2)*vel2[1]);

		//x comp
		GalerkinRHS[4]=-1.0*(msWorkMatrix(2,0)*vel0[0]+msWorkMatrix(2,1)*vel1[0]+msWorkMatrix(2,2)*vel2[0]);
		//y comp
		GalerkinRHS[5]=-1.0*(msWorkMatrix(2,0)*vel0[1]+msWorkMatrix(2,1)*vel1[1]+msWorkMatrix(2,2)*vel2[1]);

		//////////////////////////////////////
		//////////////////////////////////////
	        /*GetGeometry()[0].SetLock();
		array_1d<double,3>& rmu0 = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
		rmu0[0] += GalerkinRHS[0];
		rmu0[1] += GalerkinRHS[1];
		GetGeometry()[0].UnSetLock();

		GetGeometry()[1].SetLock();
		array_1d<double,3>& rmu1 = GetGeometry()[1].FastGetSolutionStepValue(ANGULAR_VELOCITY);
		rmu1[0] += GalerkinRHS[2];
		rmu1[1] += GalerkinRHS[3];
		GetGeometry()[1].UnSetLock();
	
		GetGeometry()[2].SetLock();
		array_1d<double,3>& rmu2 = GetGeometry()[2].FastGetSolutionStepValue(ANGULAR_VELOCITY);
		rmu2[0] += GalerkinRHS[4];
		rmu2[1] += GalerkinRHS[5];
		GetGeometry()[2].UnSetLock();*/
		//////////////////////////////////////
		//////////////////////////////////////

		
		//density=1.0;


		ms_adv_vel[0] = msN[0]*(vel0[0])+msN[1]*(vel1[0])+msN[2]*(vel2[0]);
		ms_adv_vel[1] = msN[0]*(vel0[1])+msN[1]*(vel1[1])+msN[2]*(vel2[1]);

	
		
		const array_1d<double,3> body_force = 0.333333333333333*(GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE)+						GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE) +	GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE));
		unsigned int number_of_nodes=3;
		
		for(unsigned int i = 0; i<number_of_nodes; i++)
		{
			GalerkinRHS[i*2] += body_force[0]* density * Area * 0.3333333333333;
			GalerkinRHS[i*2+1] += body_force[1] * density * Area * 0.3333333333333;
		}
		
		
		double p_avg=0.333333333333*(p_n0+p_n1+p_n2)*Area;



		GalerkinRHS[0]+=msDN_DX(0,0)*p_avg;
		GalerkinRHS[1]+=msDN_DX(0,1)*p_avg; 

		GalerkinRHS[2]+=msDN_DX(1,0)*p_avg; 
		GalerkinRHS[3]+=msDN_DX(1,1)*p_avg; 

		GalerkinRHS[4]+=msDN_DX(2,0)*p_avg; 
		GalerkinRHS[5]+=msDN_DX(2,1)*p_avg; 
		


		//////////////////////////////////////
		//////////////////////////////////////
	        /*GetGeometry()[0].SetLock();
		array_1d<double,3>& rmp0 = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_ACCELERATION);
		rmp0[0] += msDN_DX(0,0)*p_avg;
		rmp0[1] += msDN_DX(0,1)*p_avg; 
		GetGeometry()[0].UnSetLock();

		GetGeometry()[1].SetLock();
		array_1d<double,3>& rmp1 = GetGeometry()[1].FastGetSolutionStepValue(ANGULAR_ACCELERATION);
		rmp1[0] += msDN_DX(1,0)*p_avg; 
		rmp1[1] += msDN_DX(1,1)*p_avg; 
		GetGeometry()[1].UnSetLock();
	
		GetGeometry()[2].SetLock();
		array_1d<double,3>& rmp2 = GetGeometry()[2].FastGetSolutionStepValue(ANGULAR_ACCELERATION);
		rmp2[0] += msDN_DX(2,0)*p_avg; 
		rmp2[1] += msDN_DX(2,1)*p_avg; 
		GetGeometry()[2].UnSetLock();*/
		//////////////////////////////////////
		//////////////////////////////////////


		GetGeometry()[0].SetLock();
		//array_1d<double,3>& rhsg0 = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_ACCELERATION);
		//rhsg0[0] += body_force[0]* Area * 0.3333333333333;// *density ;
		//rhsg0[1] += body_force[1]* Area * 0.3333333333333;// *density ;

		array_1d<double,3>& rhs0 = GetGeometry()[0].FastGetSolutionStepValue(FORCE);
		rhs0[0] += GalerkinRHS[0];
		rhs0[1] += GalerkinRHS[1];

		GetGeometry()[0].UnSetLock();

		GetGeometry()[1].SetLock();
		//array_1d<double,3>& rhsg1 = GetGeometry()[1].FastGetSolutionStepValue(ANGULAR_ACCELERATION);
		//rhsg1[0] += body_force[0]* Area * 0.3333333333333;// *density ;
		//rhsg1[1] += body_force[1]* Area * 0.3333333333333;// *density;

		array_1d<double,3>& rhs1 = GetGeometry()[1].FastGetSolutionStepValue(FORCE);
		rhs1[0] += GalerkinRHS[2];
		rhs1[1] += GalerkinRHS[3];
		GetGeometry()[1].UnSetLock();
	
		GetGeometry()[2].SetLock();
		//array_1d<double,3>& rhsg2 = GetGeometry()[2].FastGetSolutionStepValue(ANGULAR_ACCELERATION);
		//rhsg2[0] += body_force[0]* Area * 0.3333333333333;//*density;
		//rhsg2[1] += body_force[1]* Area * 0.3333333333333;//*density;

		array_1d<double,3>& rhs2 = GetGeometry()[2].FastGetSolutionStepValue(FORCE);
		rhs2[0] += GalerkinRHS[4];
		rhs2[1] += GalerkinRHS[5];
		GetGeometry()[2].UnSetLock();


		
		/*GetGeometry()[0].SetLock();
		array_1d<double,3>& rhs0 = GetGeometry()[0].FastGetSolutionStepValue(FORCE);
		rhs0[0] += GalerkinRHS[0];
		rhs0[1] += GalerkinRHS[1];
		GetGeometry()[0].UnSetLock();*/

/*		GetGeometry()[1].SetLock();
		array_1d<double,3>& rhs1 = GetGeometry()[1].FastGetSolutionStepValue(FORCE);
		rhs1[0] += GalerkinRHS[2];
		rhs1[1] += GalerkinRHS[3];
		GetGeometry()[1].UnSetLock();*/
	
		/*GetGeometry()[2].SetLock();
		array_1d<double,3>& rhs2 = GetGeometry()[2].FastGetSolutionStepValue(FORCE);
		rhs2[0] += GalerkinRHS[4];
		rhs2[1] += GalerkinRHS[5];
		GetGeometry()[2].UnSetLock();*/


		double nodal_contrib = 0.333333333333333333333333333 * Area;

        	GetGeometry()[0].SetLock();
        	GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib * density;
		GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MASS) += nodal_contrib ;
        	GetGeometry()[0].UnSetLock();

        	GetGeometry()[1].SetLock();
        	GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib * density;
		GetGeometry()[1].FastGetSolutionStepValue(PARTICLE_MASS) += nodal_contrib ;
        	GetGeometry()[1].UnSetLock();

        	GetGeometry()[2].SetLock();
        	GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib * density;
		GetGeometry()[2].FastGetSolutionStepValue(PARTICLE_MASS) += nodal_contrib;
        	GetGeometry()[2].UnSetLock();






        }
        else if (FractionalStepNumber == 6) //calculation of velocities
        {
	//KRATOS_ERROR(std::logic_error, "method not implemented", "");

		///////////////////////NECESSARY LOCALS///////////////////////////////////////////
		boost::numeric::ublas::bounded_matrix<double,3,3> msWorkMatrix = ZeroMatrix(3,3);
		array_1d<double,6> GalerkinRHS = ZeroVector(6); //dimension = number of nodes
		boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX = ZeroMatrix(3,2);
		array_1d<double,3> msN = ZeroVector(3); //dimension = number of nodes
		boost::numeric::ublas::bounded_matrix<double,6,2> msShapeFunc = ZeroMatrix(6,2);
		boost::numeric::ublas::bounded_matrix<double,2,6> msConvOp = ZeroMatrix(2,6);
	       	boost::numeric::ublas::bounded_matrix<double,6,6> msAuxMat = ZeroMatrix(6,6);
		array_1d<double,6> msAuxVec = ZeroVector(6); //dimension = number of nodes
		array_1d<double,2> ms_adv_vel = ZeroVector(2); //dimesion coincides with space dimension
		array_1d<double,2> ms_vel_gauss = ZeroVector(2); //dimesion coincides with space dimension

		///////////////////////////////////////////////////////////////////////////////////

	
		//first we compute  the force term and pressure gradient terms:
		//getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Area);

		//getting the velocity on the nodes and other necessary variabless
		const array_1d<double,3> vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
		double p_n0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
		const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
		const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
		
		const array_1d<double,3> vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
		double p_n1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
		const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
		const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
		
		const array_1d<double,3>& vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);
		double p_n2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
		const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
		const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
			
		//====================================================================
		//calculating viscosity and density
		double nu = 0.333333333333333333333333*(nu0 + nu1 + nu2 );
		double density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );
		
		noalias(msWorkMatrix) = Area*density*nu * prod(msDN_DX,trans(msDN_DX));
				
		//x comp
		GalerkinRHS[0]=-1.0*(msWorkMatrix(0,0)*vel0[0]+msWorkMatrix(0,1)*vel1[0]+msWorkMatrix(0,2)*vel2[0]);
		//y comp
		GalerkinRHS[1]=-1.0*(msWorkMatrix(0,0)*vel0[1]+msWorkMatrix(0,1)*vel1[1]+msWorkMatrix(0,2)*vel2[1]);

		//x comp
		GalerkinRHS[2]=-1.0*(msWorkMatrix(1,0)*vel0[0]+msWorkMatrix(1,1)*vel1[0]+msWorkMatrix(1,2)*vel2[0]);
		//y comp
		GalerkinRHS[3]=-1.0*(msWorkMatrix(1,0)*vel0[1]+msWorkMatrix(1,1)*vel1[1]+msWorkMatrix(1,2)*vel2[1]);

		//x comp
		GalerkinRHS[4]=-1.0*(msWorkMatrix(2,0)*vel0[0]+msWorkMatrix(2,1)*vel1[0]+msWorkMatrix(2,2)*vel2[0]);
		//y comp
		GalerkinRHS[5]=-1.0*(msWorkMatrix(2,0)*vel0[1]+msWorkMatrix(2,1)*vel1[1]+msWorkMatrix(2,2)*vel2[1]);

		ms_adv_vel[0] = msN[0]*(vel0[0])+msN[1]*(vel1[0])+msN[2]*(vel2[0]);
		ms_adv_vel[1] = msN[0]*(vel0[1])+msN[1]*(vel1[1])+msN[2]*(vel2[1]);

	
		
		const array_1d<double,3> body_force = 0.333333333333333*(GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE)+	GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE) +	GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE));
		unsigned int number_of_nodes=3;
		for(unsigned int i = 0; i<number_of_nodes; i++)
		{
			GalerkinRHS[i*2] += body_force[0]* density * Area * 0.3333333333333;
			GalerkinRHS[i*2+1] += body_force[1] * density * Area * 0.3333333333333;
		}
		
		
		double p_avg=0.333333333333*(p_n0+p_n1+p_n2)*Area;
		GalerkinRHS[0]+=msDN_DX(0,0)*p_avg;
		GalerkinRHS[1]+=msDN_DX(0,1)*p_avg; 

		GalerkinRHS[2]+=msDN_DX(1,0)*p_avg; 
		GalerkinRHS[3]+=msDN_DX(1,1)*p_avg; 

		GalerkinRHS[4]+=msDN_DX(2,0)*p_avg; 
		GalerkinRHS[5]+=msDN_DX(2,1)*p_avg; 
		
		
		GetGeometry()[0].SetLock();
		array_1d<double,3>& rhs0 = GetGeometry()[0].FastGetSolutionStepValue(FORCE);
		rhs0[0] += GalerkinRHS[0];
		rhs0[1] += GalerkinRHS[1];
		GetGeometry()[0].UnSetLock();

		GetGeometry()[1].SetLock();
		array_1d<double,3>& rhs1 = GetGeometry()[1].FastGetSolutionStepValue(FORCE);
		rhs1[0] += GalerkinRHS[2];
		rhs1[1] += GalerkinRHS[3];
		GetGeometry()[1].UnSetLock();
	
		GetGeometry()[2].SetLock();
		array_1d<double,3>& rhs2 = GetGeometry()[2].FastGetSolutionStepValue(FORCE);
		rhs2[0] += GalerkinRHS[4];
		rhs2[1] += GalerkinRHS[5];
		GetGeometry()[2].UnSetLock();


		double nodal_contrib = 0.333333333333333333333333333 * Area*density;

        	GetGeometry()[0].SetLock();
        	GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
        	GetGeometry()[0].UnSetLock();

        	GetGeometry()[1].SetLock();
        	GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
        	GetGeometry()[1].UnSetLock();

        	GetGeometry()[2].SetLock();
        	GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
        	GetGeometry()[2].UnSetLock();




        }

        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************

    void Fluid2DGLS::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
    {
        unsigned int number_of_nodes = GetGeometry().PointsNumber();
        if (rResult.size() != number_of_nodes)
            rResult.resize(number_of_nodes, false);

        unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

        if (FractionalStepNumber == 1) //step 1
            for (unsigned int i = 0; i < number_of_nodes; i++)
                rResult[i] = GetGeometry()[i].GetDof(FRACT_VEL_X).EquationId();
        else if (FractionalStepNumber == 2) //step 2
            for (unsigned int i = 0; i < number_of_nodes; i++)
                rResult[i] = GetGeometry()[i].GetDof(FRACT_VEL_Y).EquationId();

        else if (FractionalStepNumber == 4) // pressure correction step
            for (unsigned int i = 0; i < number_of_nodes; i++)
                rResult[i] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
    }

    //************************************************************************************
    //************************************************************************************

    void Fluid2DGLS::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
    {
        unsigned int number_of_nodes = GetGeometry().PointsNumber();
        if (ElementalDofList.size() != number_of_nodes)
            ElementalDofList.resize(number_of_nodes);

        unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

        if (FractionalStepNumber == 1) //step 1
            for (unsigned int i = 0; i < number_of_nodes; i++)
                ElementalDofList[i] = GetGeometry()[i].pGetDof(FRACT_VEL_X);
        else if (FractionalStepNumber == 2) //step 2
            for (unsigned int i = 0; i < number_of_nodes; i++)
                ElementalDofList[i] = GetGeometry()[i].pGetDof(FRACT_VEL_Y);
        else if (FractionalStepNumber == 4) // pressure correction step
            for (unsigned int i = 0; i < number_of_nodes; i++)
                ElementalDofList[i] = GetGeometry()[i].pGetDof(PRESSURE);

    }

    //************************************************************************************
    //************************************************************************************

    double Fluid2DGLS::ComputeSmagorinskyViscosity(const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX,
            const double& h,
            const double& C,
            const double nu
            )
    {
        boost::numeric::ublas::bounded_matrix<double, 2, 2 > dv_dx = ZeroMatrix(2, 2);

        // Compute Symmetric Grad(u). Note that only the lower half of the matrix is filled
        for (unsigned int k = 0; k < 3; ++k)
        {
            const array_1d< double, 3 > & rNodeVel = this->GetGeometry()[k].FastGetSolutionStepValue(FRACT_VEL);
            for (unsigned int i = 0; i < 2; ++i)
            {
                for (unsigned int j = 0; j < i; ++j) // Off-diagonal
                    dv_dx(i, j) += 0.5 * (DN_DX(k, j) * rNodeVel[i] + DN_DX(k, i) * rNodeVel[j]);
                dv_dx(i, i) += DN_DX(k, i) * rNodeVel[i]; // Diagonal
            }
        }

        // Norm[ Grad(u) ]
        double NormS(0.0);
        for (unsigned int i = 0; i < 2; ++i)
        {
            for (unsigned int j = 0; j < i; ++j)
                NormS += 2.0 * dv_dx(i, j) * dv_dx(i, j); // Using symmetry, lower half terms of the matrix are added twice
            NormS += dv_dx(i, i) * dv_dx(i, i); // Diagonal terms
        }

        NormS = sqrt(NormS);

        // Total Viscosity
        return 2.0 * C * C * h * h * NormS;
    }

    //************************************************************************************

    int Fluid2DGLS::Check(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        //check the area

        //check if if is in the XY plane

        //check that no variable has zero key


        return 0;


        KRATOS_CATCH("");
    }


    //#undef GRADPN_FORM
} // Namespace Kratos


