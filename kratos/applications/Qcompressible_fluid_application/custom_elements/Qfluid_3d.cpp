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
//   Date:                $Date: 2009-01-23 14:33:59 $
//   Revision:            $Revision: 1.1 $
//
//


// System includes 
//#define GRADPN_FORM //the grad(pn) is used instead of the G(pn) in doing the splitting
 
// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/Qfluid_3d.h"
#include "Qcompressible_fluid_application.h"
#include "utilities/math_utils.h" 
#include "utilities/geometry_utilities.h" 

namespace Kratos 
{
	//static variables
	boost::numeric::ublas::bounded_matrix<double,4,4> QFluid3D::msaux_matrix;
	boost::numeric::ublas::bounded_matrix<double,4,4> QFluid3D::msMassFactors;
	boost::numeric::ublas::bounded_matrix<double,4,4> QFluid3D::msMassC;
	boost::numeric::ublas::bounded_matrix<double,3,3> QFluid3D::msF;
	boost::numeric::ublas::bounded_matrix<double,4,3> QFluid3D::msDN_DX;
  //boost::numeric::ublas::bounded_matrix<double,3,6> UpdatedLagrangianFluid::ms_temp;
  	array_1d<double,4> QFluid3D::msN; //dimension = number of nodes
	array_1d<double,3> QFluid3D::ms_aux; //dimension coincides with space dimension
	array_1d<double,3> QFluid3D::ms_vel_gauss; //dimesion coincides with space dimension
  	array_1d<double,4> QFluid3D::ms_temp_vec_np; //dimension = number of nodes
  	array_1d<double,4> QFluid3D::ms_press; //dimension = number of nodes
  	array_1d<double,4> QFluid3D::ms_u_DN; //dimension = number of nodes
	
	
	/*std::vector< Matrix > mInvJ;
	Vector mDetJ;*/
	


	//************************************************************************************
	//************************************************************************************
	QFluid3D::QFluid3D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	QFluid3D::QFluid3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
	  
	  msMassFactors(0,0) = 0.25;msMassFactors(0,1) = 0.0;msMassFactors(0,2) = 0.0;msMassFactors(0,3) = 0.0;
	  msMassFactors(1,0) = 0.0;msMassFactors(1,1) = 0.25;msMassFactors(1,2) = 0.0;msMassFactors(1,3) = 0.0;
	  msMassFactors(2,0) = 0.0;msMassFactors(2,1) = 0.0;msMassFactors(2,2) = 0.25;msMassFactors(2,3) = 0.0;
	  msMassFactors(3,0) = 0.0;msMassFactors(3,1) = 0.0;msMassFactors(3,2) = 0.0;msMassFactors(3,3) = 0.25;
	  //InitializeAuxiliaries();

	  msMassC(0,0) = 1/10 ;
	  msMassC(0,1) = 1/20;
	  msMassC(0,2) = 1/20;
	  msMassC(0,3) = 1/20;
	  msMassC(1,0) = 1/20;
	  msMassC(1,1) = 1/10;
	  msMassC(1,2) = 1/20;
	  msMassC(1,3) = 1/20;
	  msMassC(2,0) = 1/20;
	  msMassC(2,1) = 1/20;
	  msMassC(2,2) = 1/10;
	  msMassC(2,3) = 1/20;
	  msMassC(3,0) = 1/20;
	  msMassC(3,1) = 1/20;
	  msMassC(3,2) = 1/20;
	  msMassC(3,3) = 1/10;
	}
  
	Element::Pointer QFluid3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		KRATOS_TRY
		return Element::Pointer(new QFluid3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
		KRATOS_CATCH("");
	}

	QFluid3D::~QFluid3D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void QFluid3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];

		if(FractionalStepNumber < 3) //first step of the fractional step solution
		{
			Stage1(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
		}
		else if (FractionalStepNumber == 4)//second step of the fractional step solution
		{
			Stage2(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
		}

		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void QFluid3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_ERROR(std::logic_error,  "method not implemented" , "");
	}

	//************************************************************************************
	//************************************************************************************
	//calculation by component of the fractional step velocity corresponding to the first stage
	void QFluid3D::Stage1(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
										   ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		const unsigned int number_of_points = 4;
		const unsigned int dim = 3;	
		const unsigned int number_of_nodes = GetGeometry().size();
		unsigned int matsize = number_of_points *dim;
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

		if(rLeftHandSideMatrix.size1() != matsize)
			rLeftHandSideMatrix.resize(matsize,matsize,false);

		if(rRightHandSideVector.size() != matsize)
			rRightHandSideVector.resize(matsize,false);
		
		//getting data for the given geometry
		double Volume;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Volume);
		//CalculateGeometryData(msDN_DX,msN,Volume);

 		int n=0;
		//getting the velocity on the nodes
		const array_1d<double,3>& fv0 = GetGeometry()[0].FastGetSolutionStepValue(FRACT_VEL);
		const array_1d<double,3>& v0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
  		double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1);//(PRESSURE_OLD_IT);
    		double p0oldaux = GetGeometry()[0].FastGetSolutionStepValue(PRESSUREAUX,1);//_OLD_IT);
		const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
		const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
                const double T0 = GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);
      		const int f0 = GetGeometry()[0].FastGetSolutionStepValue(IS_FLUID);    
    		const int f0i = GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE);
    		const double k0 = GetGeometry()[0].FastGetSolutionStepValue(BULK_MODULUS);
    		//const double& f0j = GetGeometry()[0].FastGetSolutionStepValue(FLAG_VARIABLE);
    		const int f0a = GetGeometry()[0].FastGetSolutionStepValue(FLAG_VARIABLE); 
	    	const int f_0a = GetGeometry()[0].FastGetSolutionStepValue(MATERIAL_VARIABLE);
     		const double A0 = GetGeometry()[0].FastGetSolutionStepValue(ARRHENIUS);
    		int j=0;
		if(f0a==1) j++;


		const array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL);
		const array_1d<double,3>& v1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
  		double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1);//(PRESSURE_OLD_IT);
   		double p1oldaux = GetGeometry()[1].FastGetSolutionStepValue(PRESSUREAUX,1);//_OLD_IT);
	        const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
		const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
                const double T1 = GetGeometry()[1].FastGetSolutionStepValue(TEMPERATURE);
      		const int f1 = GetGeometry()[1].FastGetSolutionStepValue(IS_FLUID);    
    		const int f1i = GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE);
    		const double k1 = GetGeometry()[1].FastGetSolutionStepValue(BULK_MODULUS);
    		//const double& f0j = GetGeometry()[0].FastGetSolutionStepValue(FLAG_VARIABLE);
    		const int f1a = GetGeometry()[1].FastGetSolutionStepValue(FLAG_VARIABLE); 
	    	const int f_1a = GetGeometry()[1].FastGetSolutionStepValue(MATERIAL_VARIABLE);
     		const double A1 = GetGeometry()[1].FastGetSolutionStepValue(ARRHENIUS);
    		if(f1a==1) j++;


		const array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL);
		const array_1d<double,3>& v2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);
  		double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1);//(PRESSURE_OLD_IT);
   		double p2oldaux = GetGeometry()[2].FastGetSolutionStepValue(PRESSUREAUX,1);//_OLD_IT);
	        const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
		const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
                const double T2 = GetGeometry()[2].FastGetSolutionStepValue(TEMPERATURE);
      		const int f2 = GetGeometry()[2].FastGetSolutionStepValue(IS_FLUID);    
    		const int f2i = GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE);
    		const double k2 = GetGeometry()[2].FastGetSolutionStepValue(BULK_MODULUS);
    		//const double& f0j = GetGeometry()[0].FastGetSolutionStepValue(FLAG_VARIABLE);
    		const int f2a = GetGeometry()[2].FastGetSolutionStepValue(FLAG_VARIABLE); 
	    	const int f_2a = GetGeometry()[2].FastGetSolutionStepValue(MATERIAL_VARIABLE);
     		const double A2 = GetGeometry()[2].FastGetSolutionStepValue(ARRHENIUS);
    		if(f2a==1) j++;

		const array_1d<double,3>& fv3 = GetGeometry()[3].FastGetSolutionStepValue(FRACT_VEL);
		const array_1d<double,3>& v3 = GetGeometry()[3].FastGetSolutionStepValue(VELOCITY);
  		double p3old = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE,1);//(PRESSURE_OLD_IT);
   		double p3oldaux = GetGeometry()[3].FastGetSolutionStepValue(PRESSUREAUX,1);//_OLD_IT);
	        const double nu3 = GetGeometry()[3].FastGetSolutionStepValue(VISCOSITY);
		const double rho3 = GetGeometry()[3].FastGetSolutionStepValue(DENSITY);
                const double T3 = GetGeometry()[3].FastGetSolutionStepValue(TEMPERATURE);
      		const int f3 = GetGeometry()[3].FastGetSolutionStepValue(IS_FLUID);    
    		const int f3i = GetGeometry()[3].FastGetSolutionStepValue(IS_INTERFACE);
    		const double k3 = GetGeometry()[3].FastGetSolutionStepValue(BULK_MODULUS);
    		//const double& f0j = GetGeometry()[0].FastGetSolutionStepValue(FLAG_VARIABLE);
    		const int f3a = GetGeometry()[3].FastGetSolutionStepValue(FLAG_VARIABLE); 
	    	const int f_3a = GetGeometry()[3].FastGetSolutionStepValue(MATERIAL_VARIABLE);
     		const double A3 = GetGeometry()[3].FastGetSolutionStepValue(ARRHENIUS);
    		if(f3a==1) j++;





		//deformation gradient
 		noalias(msF)= ZeroMatrix(3,3);
 
               	Matrix kronecker(3,3);
               	noalias(kronecker)=ZeroMatrix(3,3);
                                    
                for(int i=0; i<3;i++){ 
                     	kronecker(i,i)=1.0;
               		}
// 			//deformation gradient
 		for(int i=0; i<dim; i++)
 			for(int j=0; j<dim; j++)
 				for(int node=0; node<number_of_nodes; node++)
 				{
 					msF(i,j) += (GetGeometry()[node]).GetSolutionStepValue(DESP)(i)	*msDN_DX(node,j);//(DISPLACEMENT)(i)
					/*double s= (GetGeometry()[node]).GetSolutionStepValue(DESP)(i);
					KRATOS_WATCH(s);*/				
				}
// 
 		for(int i=0; i<dim; i++) msF(i,i)= kronecker(i,i);
/*********************/
/**************/
		MathUtils<double>::InvertMatrix(msF,mInvJ[0],mDetJ[0]);


		//calculating viscosity
		double nu = 0.25*(nu0 + nu1 + nu2 + nu3);
 		double density = 0.25*(rho0 + rho1 + rho2 + rho3);
                double temperature= 0.25*(T0 + T1 + T2 + T3);
		//getting the BDF2 coefficients (not fixed to allow variable time step)
		//the coefficients INCLUDE the time step
		const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

		double K=10000000;//1000000.0;//500000;//K=50000000
		double k=K/BDFcoeffs[0];

		/*KRATOS_WATCH(mInvJ[0]);	
		KRATOS_WATCH(msDN_DX);	*/
///////////////////
//		noalias(msDN_DX)=prod(msDN_DX,mInvJ[0]);

///////////////////no considero el jacobiano, lo supongo igual a 1.
		/*KRATOS_WATCH(msDN_DX);*/	
		CalculateViscousMatrix(rLeftHandSideMatrix, msDN_DX, nu, k,density);
		//KRATOS_WATCH(rLeftHandSideMatrix);	

			//INERTIA CONTRIBUTION
		noalias(msaux_matrix) = BDFcoeffs[0] * msMassFactors * density;
		//KRATOS_WATCH(msaux_matrix);	
		//adding all contributions to the stiffness matrix
		ExpandAndAddReducedMatrix(rLeftHandSideMatrix,msaux_matrix,dim);

		//multiplication by the area
		rLeftHandSideMatrix *= Volume;
		//KRATOS_WATCH(rLeftHandSideMatrix);	

		// *****************************************
		//CALCULATION OF THE RHS
		const array_1d<double,3>& force0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);
		const array_1d<double,3>& force1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);
		const array_1d<double,3>& force2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);		
		const array_1d<double,3>& force3 = GetGeometry()[3].FastGetSolutionStepValue(BODY_FORCE);		

		array_1d<double,4> rhs_aux;
		for( unsigned int component_index = 0; component_index < dim; component_index++)
		  {
			//external forces (component)
		    double force_component = 0.25*(force0[component_index] + force1[component_index] + force2[component_index] + force3[component_index]);
			noalias(rhs_aux) = (force_component )*msN *(1-0.01*(temperature-273))* density;//poner aqui BOUSINESQ

			//adding pressure gradient (integrated by parts)
			double p_avg = p0old + p1old + p2old + p3old;
			p_avg *= 0.25 ;
			rhs_aux[0] += msDN_DX(0,component_index)*p_avg; 
			rhs_aux[1] += msDN_DX(1,component_index)*p_avg; 
			rhs_aux[2] += msDN_DX(2,component_index)*p_avg;
			rhs_aux[3] += msDN_DX(3,component_index)*p_avg;
			
			//adding the inertia terms
			// RHS += M*vhistory 
			//calculating the historical velocity
			
			noalias(ms_temp_vec_np) = ZeroVector(4);
			//KRATOS_WATCH(ms_temp_vec_np);	
			for(unsigned int iii = 0; iii<number_of_points; iii++)
			  {
			    const array_1d<double,3>& v = (GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY,1) );
			    ms_temp_vec_np[iii] = BDFcoeffs[0]*v[component_index]*density;
			  }
			
			noalias(rhs_aux) += prod(msMassFactors,ms_temp_vec_np) ;
			//KRATOS_WATCH(rhs_aux);	
			//RHS += Suy * proj[component] 
			/*double proj_component = proj0[component_index] 
								  + proj1[component_index] 
								  + proj2[component_index]
								  + proj3[component_index];
								  proj_component *= 0.25;
								  noalias(rhs_aux) += (tau*proj_component)*ms_u_DN;*/   
			
			//writing the rhs_aux in its place
			for( unsigned int i = 0; i < number_of_points; i++)
			  {
			    rRightHandSideVector[i*dim + component_index] = rhs_aux[i];
			  }
		  }
		
		//multiplying by area
		rRightHandSideVector *= Volume;
		//KRATOS_WATCH(rRightHandSideVector);	
		//LHS dirichlet contribution
		Vector fvvect(12);
		for( unsigned int component_index = 0; component_index < dim; component_index++)
		  {
		    fvvect[0 + component_index] = fv0[component_index];
		    fvvect[3 + component_index] = fv1[component_index];
		    fvvect[6 + component_index] = fv2[component_index];
		    fvvect[9 + component_index] = fv3[component_index];
		  }
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,fvvect);
		//KRATOS_WATCH(rRightHandSideVector)
		//KRATOS_WATCH(rLeftHandSideMatrix)
		//KRATOS_WATCH(rRightHandSideVector);	
		KRATOS_CATCH("");
	}
  
  //************************************************************************************
  //************************************************************************************
  //calculation by component of the fractional step velocity corresponding to the first stage
  void QFluid3D::Stage2(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
			ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY;

		unsigned int number_of_points = 4;

		if(rLeftHandSideMatrix.size1() != number_of_points)
			rLeftHandSideMatrix.resize(number_of_points,number_of_points);

		if(rRightHandSideVector.size() != number_of_points)
			rRightHandSideVector.resize(number_of_points);

		//getting data for the given geometry
		double Volume;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Volume);
		//CalculateGeometryData(msDN_DX,msN,Volume);

		
		const array_1d<double,3>& fv0 = GetGeometry()[0].FastGetSolutionStepValue(FRACT_VEL);
		
		double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE); 
   		double p0n = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1); 
		double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_OLD_IT); 
		const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
		const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
		
		const array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL);
		
		double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
 		double p1n = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1);  
		double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT); 
		const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
		const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
		
		const array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL);
		
		double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE); 
		double p2n = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1);
		double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT); 
		const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
		const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);

		const array_1d<double,3>& fv3 = GetGeometry()[3].FastGetSolutionStepValue(FRACT_VEL);
		
		double p3 = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE); 
		double p3n = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE,1);
		double p3old = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE_OLD_IT); 
		const double nu3 = GetGeometry()[3].FastGetSolutionStepValue(VISCOSITY);
		const double rho3 = GetGeometry()[3].FastGetSolutionStepValue(DENSITY);

		
		//calculating avergage density and viscosity
		double nu = 0.25*(nu0 + nu1 + nu2 + nu3);
 		double density = 0.25*(rho0 + rho1 + rho2 + rho3);

			
		//getting the BDF2 coefficients (nhttp://www.google.com.ar/ot fixed to allow variable time step)
		//the coefficients INCLUDE the time step
		const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

		//CALCULATION OF THE LEFT HAND SIDE
		noalias(rLeftHandSideMatrix) = msMassFactors/*((1.00/BDFcoeffs[0] + tau)/density) * prod(msDN_DX,trans(msDN_DX))*/;
		//KRATOS_WATCH(rLeftHandSideMatrix);	
		//calculation of the RHS
		// RHS = -G*vfrac

 		double K=1000000;//K=50000000;
    		double k=K/BDFcoeffs[0];

		double Gaux;
		Gaux =  msDN_DX(0,0)*fv0[0] + msDN_DX(0,1)*fv0[1] + msDN_DX(0,2)*fv0[2];
		Gaux += msDN_DX(1,0)*fv1[0] + msDN_DX(1,1)*fv1[1] + msDN_DX(1,2)*fv1[2];
		Gaux += msDN_DX(2,0)*fv2[0] + msDN_DX(2,1)*fv2[1] + msDN_DX(2,2)*fv2[2];
		Gaux += msDN_DX(3,0)*fv3[0] + msDN_DX(3,1)*fv3[1] + msDN_DX(3,2)*fv3[2];

		rRightHandSideVector[0] = - Gaux * msN[0] * k; 
		rRightHandSideVector[1] = - Gaux * msN[1] * k; 
		rRightHandSideVector[2] = - Gaux * msN[2] * k; 
		rRightHandSideVector[3] = - Gaux * msN[3] * k; 
		
		//KRATOS_WATCH(rRightHandSideVector);	
  		ms_temp_vec_np[0] = p0n; /*no se si tiene sentido*/
    		ms_temp_vec_np[1] = p1n; 
    		ms_temp_vec_np[2] = p2n; 
		ms_temp_vec_np[3] = p3n;

		noalias(ms_press) = prod( /*msMassC*/msMassFactors,ms_temp_vec_np);
    		noalias(rRightHandSideVector) += ms_press; 
    		//KRATOS_WATCH(rRightHandSideVector);	

		// RHS += dt * L * pold 
		ms_temp_vec_np[0] = p0old; 
		ms_temp_vec_np[1] = p1old; 
		ms_temp_vec_np[2] = p2old; 
		ms_temp_vec_np[3] = p3old; 
		noalias(ms_press) = prod(/*msMassC*/msMassFactors,ms_temp_vec_np);
		noalias(rRightHandSideVector) -= ms_press; 
		//KRATOS_WATCH(rRightHandSideVector);	
		//multiplicating by the Volume
		rLeftHandSideMatrix *= 0/*Volume*/;
		rRightHandSideVector *= 0/*Volume*/;

		//adding contributions to nodal Volumes following the corresponding lumping term
		/*double nodal_contrib = 0.25 * Volume;
		GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA) += nodal_contrib;
		GetGeometry()[1].FastGetSolutionStepValue(NODAL_AREA) += nodal_contrib;
		GetGeometry()[2].FastGetSolutionStepValue(NODAL_AREA) += nodal_contrib;
		GetGeometry()[3].FastGetSolutionStepValue(NODAL_AREA) += nodal_contrib;*/

		KRATOS_CATCH("");
	}
	  
	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the 
	// fractional step procedure
	void QFluid3D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];
 
		//getting data for the given geometry
		double Volume;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Volume);
		//CalculateGeometryData(msDN_DX,msN,Volume);
		
		if(FractionalStepNumber  == 5) //calculation of stab	//adding contributions to nodal Volumes following the corresponding lumping term
		
		{

					}		
		else if(FractionalStepNumber == 6) //calculation of velocities
		{
		const unsigned int number_of_nodes = GetGeometry().size();
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
		const unsigned int dim = 3;	
		const Vector& BDFcoeffs = CurrentProcessInfo[BDF_COEFFICIENTS];

		double K=10000000;//1000000.0;//500000;//K=50000000
		double k=K/BDFcoeffs[0];
		int n=0;

		//calculating actual jacobian
/*		GeometryType::JacobiansType J;
		J = GetGeometry().Jacobian(J);
		KRATOS_WATCH(J);	*/
		mInvJ.resize(integration_points.size());	
		mDetJ.resize(integration_points.size());	

		/*MathUtils<double>::InvertMatrix(J[0],mInvJ[0],mDetJ[0]);

		KRATOS_WATCH(J[0]);	
		KRATOS_WATCH(mInvJ[0]);	
		KRATOS_WATCH(mDetJ[0]);	*/
		
/**************/
		//deformation gradient
 		noalias(msF)= ZeroMatrix(3,3);
 
               	Matrix kronecker(3,3);
               	noalias(kronecker)=ZeroMatrix(3,3);
                                    
                for(int i=0; i<3;i++){ 
                     	kronecker(i,i)=1.0;
               		}
// 			//deformation gradient
 		for(int i=0; i<dim; i++)
 			for(int j=0; j<dim; j++)
 				for(int node=0; node<number_of_nodes; node++)
 				{
 					msF(i,j) += (GetGeometry()[node]).GetSolutionStepValue(DESP)(i)	*msDN_DX(node,j);
					/*double s= (GetGeometry()[node]).GetSolutionStepValue(DESP)(i);
					KRATOS_WATCH(s);*/		
 				}
// 
 		for(int i=0; i<dim; i++) msF(i,i)= kronecker(i,i);
/*********************/

		MathUtils<double>::InvertMatrix(msF,mInvJ[0],mDetJ[0]);
		/*KRATOS_WATCH(msF);	
		KRATOS_WATCH(mInvJ[0]);	
		KRATOS_WATCH(mDetJ[0]);	

		KRATOS_WATCH(mInvJ[0]);	
		KRATOS_WATCH(msDN_DX);*/	
///////////////////
	//	noalias(msDN_DX)=prod(msDN_DX,mInvJ[0]);
///////////////////no considero el jacobiano, lo supongo igual a 1.
		/*KRATOS_WATCH(msDN_DX);*/


		const array_1d<double,3>& fv0 = GetGeometry()[0].FastGetSolutionStepValue(FRACT_VEL);
		const double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_OLD_IT);
                const int f0i = GetGeometry()[0].FastGetSolutionStepValue(IS_STRUCTURE); 
                if(f0i) n++;
		const array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL);
		const double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT);
		const int f1i = GetGeometry()[1].FastGetSolutionStepValue(IS_STRUCTURE);  
 		if(f1i) n++;
		const array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL);
		const double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT);
		const int f2i = GetGeometry()[2].FastGetSolutionStepValue(IS_STRUCTURE);    
		if(f2i) n++;
		const array_1d<double,3>& fv3 = GetGeometry()[3].FastGetSolutionStepValue(FRACT_VEL);
		const double p3old = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE_OLD_IT);
		const int f3i = GetGeometry()[2].FastGetSolutionStepValue(IS_STRUCTURE); 
                if(f3i) n++;
		//const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

		double p_avg = p0old + p1old + p2old + p3old;
		p_avg *= 0.25 ;
		
	
		if(n==4) 
			{

		GetGeometry()[0].FastGetSolutionStepValue(NODAL_PRESS) += 0;
		GetGeometry()[1].FastGetSolutionStepValue(NODAL_PRESS) += 0;
 	   	GetGeometry()[2].FastGetSolutionStepValue(NODAL_PRESS) += 0;			    	  			GetGeometry()[3].FastGetSolutionStepValue(NODAL_PRESS) += 0;	    	  
		
 	   	double nodal_contrib = 0.25 * Volume;
 	   	GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) += 0;
 	   	GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS) += 0;
 	   	GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS) += 0;
	   	GetGeometry()[3].FastGetSolutionStepValue(NODAL_MASS) += 0;
			}
		else
			{
		double aux0=0;
		double aux1=0;
		double aux2=0;
		double aux3=0;

	    	//double dts=0.01;
	    	//double K=1000000;//K=50000000;
		//double k=K*dts;///BDFcoeffs[0];
       

		double Gaux;
		Gaux =  msDN_DX(0,0)*fv0[0] + msDN_DX(0,1)*fv0[1] + msDN_DX(0,2)*fv0[2];
		Gaux += msDN_DX(1,0)*fv1[0] + msDN_DX(1,1)*fv1[1] + msDN_DX(1,2)*fv1[2];
		Gaux += msDN_DX(2,0)*fv2[0] + msDN_DX(2,1)*fv2[1] + msDN_DX(2,2)*fv2[2];
		Gaux += msDN_DX(3,0)*fv3[0] + msDN_DX(3,1)*fv3[1] + msDN_DX(3,2)*fv3[2];

		aux0 = - Gaux * msN[0] * k + /*p0old*/p_avg * msN[0]; 
		aux1 = - Gaux * msN[1] * k + /*p1old*/p_avg * msN[1]; 
		aux2 = - Gaux * msN[2] * k + /*p2old*/p_avg * msN[2]; 
		aux3 = - Gaux * msN[3] * k + /*p3old*/p_avg * msN[3]; 

		GetGeometry()[0].FastGetSolutionStepValue(NODAL_PRESS) += aux0 * Volume;
		GetGeometry()[1].FastGetSolutionStepValue(NODAL_PRESS) += aux1 * Volume;
 	   	GetGeometry()[2].FastGetSolutionStepValue(NODAL_PRESS) += aux2 * Volume;		    	  
		GetGeometry()[3].FastGetSolutionStepValue(NODAL_PRESS) += aux3 * Volume;	    	  
		

 

 	   	double nodal_contrib = 0.25 * Volume;
 	   	GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
 	   	GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
 	   	GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
	   	GetGeometry()[3].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;


			} 
		}


		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void QFluid3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
		const unsigned int dim = 3;

		unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

		if(FractionalStepNumber == 1) //step 1
		{
			if(rResult.size() != number_of_nodes*dim)
				rResult.resize(number_of_nodes*dim,false);	

			for (unsigned int i=0;i<number_of_nodes;i++)
			{
				rResult[i*dim] = GetGeometry()[i].GetDof(FRACT_VEL_X).EquationId();
				rResult[i*dim+1] = GetGeometry()[i].GetDof(FRACT_VEL_Y).EquationId();
				rResult[i*dim+2] = GetGeometry()[i].GetDof(FRACT_VEL_Z).EquationId();
			}
		}
		else if(FractionalStepNumber == 4) // pressure correction step
		{
			if(rResult.size() != number_of_nodes)
				rResult.resize(number_of_nodes,false);	
			for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
		}

		}

	//************************************************************************************
	//************************************************************************************
  void QFluid3D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
  {
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dim = 3;
    
    unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];
    
    if(FractionalStepNumber == 1) //step 1
      {
			if(ElementalDofList.size() != number_of_nodes*dim)
			  ElementalDofList.resize(number_of_nodes*dim);	
			
			for (unsigned int i=0;i<number_of_nodes;i++)
			  {
			    ElementalDofList[i*dim] = GetGeometry()[i].pGetDof(FRACT_VEL_X);
			    ElementalDofList[i*dim+1] = GetGeometry()[i].pGetDof(FRACT_VEL_Y);
			    ElementalDofList[i*dim+2] = GetGeometry()[i].pGetDof(FRACT_VEL_Z);
			  }
      }
    else if(FractionalStepNumber == 4) // pressure correction step
      {
	if(ElementalDofList.size() != number_of_nodes)
				ElementalDofList.resize(number_of_nodes);	
	
	for (unsigned int i=0;i<number_of_nodes;i++)
	  ElementalDofList[i] = GetGeometry()[i].pGetDof(PRESSURE);
      }
  }
  
  
  //************************************************************************************
  //************************************************************************************
  void QFluid3D::CalculateViscousMatrix(MatrixType& K, const boost::numeric::ublas::bounded_matrix<double,4,3>& DN_DX, const double& nu, const double& k, const double& density)
  {
    //horrorful ... just to make sure it works
    Matrix B(6,12);
    noalias(B) = ZeroMatrix(6,12);
    unsigned int start;
    double betta=0.666;
    
    double c;
    c=k; /*betta=0;*/
    double p;
    for(unsigned int i = 0; i<4; i++)
      {
      start=i*3;
      B(0,start) = DN_DX(i,0);
      B(1,start+1) = DN_DX(i,1);
      B(2,start+2) = DN_DX(i,2);
      B(3,start) = DN_DX(i,1);	B(3,start+1) = DN_DX(i,0);
      B(4,start) = DN_DX(i,2);	B(4,start+2) = DN_DX(i,0);
      B(5,start+1) = DN_DX(i,2);	B(5,start+2) = DN_DX(i,1);
      }
    
      Matrix ms_constitutive_matrix(6,6);
      noalias(ms_constitutive_matrix) = ZeroMatrix(6,6);
      ms_constitutive_matrix(0,0) = c;          ms_constitutive_matrix(0,1) = c ;	ms_constitutive_matrix(0,2) = c;
      ms_constitutive_matrix(0,3)= 0.0;	        ms_constitutive_matrix(0,4) = 0.0 ;	ms_constitutive_matrix(0,5) = 0.0;
				
      ms_constitutive_matrix(1,0) = c; 	        ms_constitutive_matrix(1,1) = c;	ms_constitutive_matrix(1,2) = c;
      ms_constitutive_matrix(1,3) = 0.0; 	ms_constitutive_matrix(1,4) = 0.0;	ms_constitutive_matrix(1,5) = 0.0;

      ms_constitutive_matrix(2,0) = c;	        ms_constitutive_matrix(2,1) = c;	ms_constitutive_matrix(2,2) = c;
      ms_constitutive_matrix(2,3) = 0.0;	ms_constitutive_matrix(2,4) = 0.0;	ms_constitutive_matrix(2,5) = 0.0;
		
      ms_constitutive_matrix(3,0) = 0.0;	ms_constitutive_matrix(3,1) = 0.0;	ms_constitutive_matrix(3,2) = 0.0;
      ms_constitutive_matrix(3,3) = 0.0;	ms_constitutive_matrix(3,4) = 0.0;	ms_constitutive_matrix(3,5) = 0.0;

      ms_constitutive_matrix(4,0) = 0.0;	ms_constitutive_matrix(4,1) = 0.0;	ms_constitutive_matrix(4,2) = 0.0;
      ms_constitutive_matrix(4,3) = 0.0;	ms_constitutive_matrix(4,4) = 0.0;	ms_constitutive_matrix(4,5) = 0.0;
		
      ms_constitutive_matrix(5,0) = 0.0;	ms_constitutive_matrix(5,1) = 0.0;	ms_constitutive_matrix(5,2) = 0.0;
      ms_constitutive_matrix(5,3) = 0.0;	ms_constitutive_matrix(5,4) = 0.0;	ms_constitutive_matrix(5,5) = 0.0;


 /*     Matrix msCapx(6,6);
      noalias(msCapx) = ZeroMatrix(6,6);
      msCapx(0,0)=2.0*nu; msCapx(1,1)=2.0*nu; msCapx(2,2)=2.0*nu;
      msCapx(3,3)=nu;     msCapx(4,4)=nu;     msCapx(5,5)=nu;*/
      
      Matrix temp(6,12);
      noalias(temp) = prod(ms_constitutive_matrix,B);
      
      //noalias(K) = prod(trans(B),temp);
	const double& a = nu;
      noalias(ms_constitutive_matrix) = ZeroMatrix(6,6);

      ms_constitutive_matrix(0,0) = (4.0/3.0)*a+c;	ms_constitutive_matrix(0,1) = -2.0/3.0*a+c;	ms_constitutive_matrix(0,2) = -2.0/3.0*a+c;
      ms_constitutive_matrix(0,3) = 0.0;			ms_constitutive_matrix(0,4) = 0.0;			ms_constitutive_matrix(0,5) = 0.0;
		
      ms_constitutive_matrix(1,0) = -2.0/3.0*a+c; 	ms_constitutive_matrix(1,1) = 4.0/3.0*a+c;	ms_constitutive_matrix(1,2) = -2.0/3.0*a+c;
      ms_constitutive_matrix(1,3) = 0.0;		 	ms_constitutive_matrix(1,4) = 0.0;			ms_constitutive_matrix(1,5) = 0.0;

      ms_constitutive_matrix(2,0) = -2.0/3.0*a+c;	ms_constitutive_matrix(2,1) = -2.0/3.0*a+c;	ms_constitutive_matrix(2,2) = 4.0/3.0*a+c;
      ms_constitutive_matrix(2,3) = 0.0;			ms_constitutive_matrix(2,4) = 0.0;			ms_constitutive_matrix(2,5) = 0.0;
		
      ms_constitutive_matrix(3,0) = 0.0;			ms_constitutive_matrix(3,1) = 0.0;			ms_constitutive_matrix(3,2) = 0.0;
      ms_constitutive_matrix(3,3) = a;			ms_constitutive_matrix(3,4) = 0.0;			ms_constitutive_matrix(3,5) = 0.0;
		
		ms_constitutive_matrix(4,0) = 0.0;			ms_constitutive_matrix(4,1) = 0.0;			ms_constitutive_matrix(4,2) = 0.0;
		ms_constitutive_matrix(4,3) = 0.0;			ms_constitutive_matrix(4,4) = a;			ms_constitutive_matrix(4,5) = 0.0;
		
		ms_constitutive_matrix(5,0) = 0.0;			ms_constitutive_matrix(5,1) = 0.0;			ms_constitutive_matrix(5,2) = 0.0;
		ms_constitutive_matrix(5,3) = 0.0;			ms_constitutive_matrix(5,4) = 0.0;			ms_constitutive_matrix(5,5) = a;

temp = prod( ms_constitutive_matrix , B);

noalias(K) = prod(trans(B),temp); 
   // KRATOS_WATCH(K);
    


/*    K(0,0) = 2.0 * nu * pow(DN_DX(0,0), 2) + nu * pow(DN_DX(0,1), 2)+ nu * pow(DN_DX(0,2), 2)-(betta*nu-c)*DN_DX(0,0)*DN_DX(0,0);
    p=K(0,0);
    K(0,1) = DN_DX(0,1) * nu * DN_DX(0,0)-(betta*nu-c)*DN_DX(0,0)*DN_DX(0,1);
    p=K(0,1);
    K(0,2) = DN_DX(0,2) * nu * DN_DX(0,0)-(betta*nu-c)*DN_DX(0,0)*DN_DX(0,2);
    p=K(0,2);
    
    K(0,3) = 2.00 * DN_DX(0,0) * nu * DN_DX(1,0) + DN_DX(0,1) * nu * DN_DX(1,1) + DN_DX(0,2) * nu * DN_DX(1,2)-(betta*nu-c)*DN_DX(0,0)*DN_DX(1,0);
    p=K(0,3);    
    K(0,4) = DN_DX(0,1) * nu * DN_DX(1,0)-(betta*nu-c)*DN_DX(0,0)*DN_DX(1,1);
    p=K(0,4);    
    K(0,5) = DN_DX(0,2) * nu * DN_DX(1,0)-(betta*nu-c)*DN_DX(0,0)*DN_DX(1,2);
    p=K(0,5);
    K(0,6) = 2.00 * DN_DX(0,0) * nu * DN_DX(2,0) + DN_DX(0,1) * nu * DN_DX(2,1) + DN_DX(0,2) * nu * DN_DX(2,2)-(betta*nu-c)*DN_DX(0,0)*DN_DX(2,0);
    p=K(0,6);    
    K(0,7) = DN_DX(0,1) * nu * DN_DX(2,0)-(betta*nu-c)*DN_DX(0,0)*DN_DX(2,1);
    p=K(0,7);
    K(0,8) = DN_DX(0,2) * nu * DN_DX(2,0)-(betta*nu-c)*DN_DX(0,0)*DN_DX(2,2);
    p=K(0,8);
    K(0,9) = 2.00 * DN_DX(0,0) * nu * DN_DX(3,0) + DN_DX(0,1) * nu * DN_DX(3,1) + DN_DX(0,2) * nu * DN_DX(3,2)-(betta*nu-c)*DN_DX(0,0)*DN_DX(3,0);
	   p=K(0,9);
    K(0,10) = DN_DX(0,1) * nu * DN_DX(3,0)-(betta*nu-c)*DN_DX(0,0)*DN_DX(3,1);
   p=K(0,10);    
	K(0,11) = DN_DX(0,2) * nu * DN_DX(3,0)-(betta*nu-c)*DN_DX(0,0)*DN_DX(3,2);  
       p=K(0,11);


////////////////////
    K(1,0) = DN_DX(0,0) * nu * DN_DX(0,1)-(betta*nu-c)*DN_DX(0,1)*DN_DX(0,0);
   p=K(1,0);    
K(1,1) = 2.00 * nu * pow(DN_DX(0,1), 2) + pow(DN_DX(0,0), 2) * nu + nu * pow(DN_DX(0,2), 2)-(betta*nu-c)*DN_DX(0,1)*DN_DX(0,1);
p=K(1,1);   
 K(1,2) = DN_DX(0,2) * nu * DN_DX(0,1)-(betta*nu-c)*DN_DX(0,1)*DN_DX(0,2);    
    p=K(1,2);


    K(1,3) = DN_DX(0,0) * nu * DN_DX(1,1)-(betta*nu-c)*DN_DX(0,1)*DN_DX(1,0);
p=K(1,3);    
K(1,4) = 2.00 * DN_DX(0,1) * nu * DN_DX(1,1) + DN_DX(0,0) * nu * DN_DX(1,0) + DN_DX(0,2) * nu * DN_DX(1,2)-(betta*nu-c)*DN_DX(0,1)*DN_DX(1,1);
p=K(1,4);    
K(1,5) = DN_DX(0,2) * nu * DN_DX(1,1)-(betta*nu-c)*DN_DX(0,1)*DN_DX(1,2);
    p=K(1,5);


    K(1,6) = DN_DX(0,0) * nu * DN_DX(2,1)-(betta*nu-c)*DN_DX(0,1)*DN_DX(2,1);
p=K(1,6);    
K(1,7) = 2.00 * DN_DX(0,1) * nu * DN_DX(2,1) + DN_DX(0,0) * nu * DN_DX(2,0) + DN_DX(0,2) * nu * DN_DX(2,2)-(betta*nu-c)*DN_DX(0,1)*DN_DX(2,1);
p=K(1,7);    
K(1,8) = DN_DX(0,2) * nu * DN_DX(2,1)-(betta*nu-c)*DN_DX(0,1)*DN_DX(2,2);
p=K(1,8);        


    K(1,9) = DN_DX(0,0) * nu * DN_DX(3,1)-(betta*nu-c)*DN_DX(0,1)*DN_DX(3,1);
p=K(1,9);    
K(1,10) = 2.00 * DN_DX(0,1) * nu * DN_DX(3,1) + DN_DX(0,0) * nu * DN_DX(3,0) + DN_DX(0,2) * nu * DN_DX(3,2)-(betta*nu-c)*DN_DX(0,1)*DN_DX(3,1);
p=K(1,10);    
K(1,11) = DN_DX(0,2) * nu * DN_DX(3,1)-(betta*nu-c)*DN_DX(0,1)*DN_DX(3,2);
p=K(1,11);



     /////////////////
    K(2,0) = DN_DX(0,2) * nu * DN_DX(0,0)-(betta*nu-c)*DN_DX(0,2)*DN_DX(0,0)-(betta*nu-c)*DN_DX(0,2)*DN_DX(0,0);
    p=K(2,0);    
    K(2,1) = DN_DX(0,2) * nu * DN_DX(0,1)-(betta*nu-c)*DN_DX(0,2)*DN_DX(0,1);
    p=K(2,1);
    K(2,2) = 2.00 * nu * pow(DN_DX(0,2), 2) + pow(DN_DX(0,0), 2) * nu + nu * pow(DN_DX(0,1), 2)-(betta*nu-c)*DN_DX(0,2)*DN_DX(0,2);
    p=K(2,2);   
    

    K(2,3) = DN_DX(0,0) * nu * DN_DX(1,2)-(betta*nu-c)*DN_DX(0,2)*DN_DX(1,0);
    p=K(2,3);
    K(2,4) = DN_DX(0,1) * nu * DN_DX(1,2)-(betta*nu-c)*DN_DX(0,2)*DN_DX(1,1);
    p=K(2,4);   
    K(2,5) = 2.00 * DN_DX(0,2) * nu * DN_DX(1,2) + DN_DX(0,0) * nu * DN_DX(1,0) + DN_DX(0,1) * nu * DN_DX(1,1)-(betta*nu-c)*DN_DX(0,2)*DN_DX(1,2);
    p=K(2,5);   
    

    K(2,6) = DN_DX(0,0) * nu * DN_DX(2,2)-(betta*nu-c)*DN_DX(0,2)*DN_DX(2,0);
    p=K(2,6);    
    K(2,7) = DN_DX(0,1) * nu * DN_DX(2,2)-(betta*nu-c)*DN_DX(0,2)*DN_DX(2,1);
    p=K(2,7);    
    K(2,8) = 2.00 * DN_DX(0,2) * nu * DN_DX(2,2) + DN_DX(0,0) * nu * DN_DX(2,0) + DN_DX(0,1) * nu * DN_DX(2,1)-(betta*nu-c)*DN_DX(0,2)*DN_DX(2,2);
    p=K(2,8); 
    

    K(2,9) = DN_DX(0,0) * nu * DN_DX(3,2)-(betta*nu-c)*DN_DX(0,2)*DN_DX(3,0);
    p=K(2,9);    
    K(2,10) = DN_DX(0,1) * nu * DN_DX(3,2)-(betta*nu-c)*DN_DX(0,2)*DN_DX(3,1);
    p=K(2,10);    
    K(2,11) = 2.00 * DN_DX(0,2) * nu * DN_DX(3,2) + DN_DX(0,0) * nu * DN_DX(3,0) + DN_DX(0,1) * nu * DN_DX(3,1)-(betta*nu-c)*DN_DX(0,2)*DN_DX(3,2);
    p=K(2,11);
    //////////////


    K(3,0) = 2.00 * DN_DX(0,0) * nu * DN_DX(1,0) + DN_DX(0,1) * nu * DN_DX(1,1) + DN_DX(0,2) * nu * DN_DX(1,2)-(betta*nu-c)*DN_DX(1,0)*DN_DX(0,0);
    p=K(3,0);
    K(3,1) = DN_DX(0,0) * nu * DN_DX(1,1)-(betta*nu-c)*DN_DX(1,0)*DN_DX(0,1);
    p=K(3,1);
    K(3,2) = DN_DX(0,0) * nu * DN_DX(1,2)-(betta*nu-c)*DN_DX(1,0)*DN_DX(0,2);
    p=K(3,2); 

    K(3,3) = 2.00 * pow(DN_DX(1,0), 2) * nu + nu * pow(DN_DX(1,1), 2) + nu * pow(DN_DX(1,2), 2)-(betta*nu-c)*DN_DX(1,0)*DN_DX(1,0);
    p=K(3,3);    
    K(3,4) = DN_DX(1,1) * nu * DN_DX(1,0)-(betta*nu-c)*DN_DX(1,0)*DN_DX(1,1);
    p=K(3,4);    
    K(3,5) = DN_DX(1,2) * nu * DN_DX(1,0)-(betta*nu-c)*DN_DX(1,0)*DN_DX(1,2);
    p=K(3,5);


    K(3,6) = 2.00 * DN_DX(1,0) * nu * DN_DX(2,0) + DN_DX(1,1) * nu * DN_DX(2,1) + DN_DX(1,2) * nu * DN_DX(2,2)-(betta*nu-c)*DN_DX(1,0)*DN_DX(2,0);
    p=K(3,6);    
    K(3,7) = DN_DX(1,1) * nu * DN_DX(2,0)-(betta*nu-c)*DN_DX(1,0)*DN_DX(2,1);
    p=K(3,7);    
    K(3,8) = DN_DX(1,2) * nu * DN_DX(2,0)-(betta*nu-c)*DN_DX(1,0)*DN_DX(2,2);
    p=K(3,8);



    K(3,9) = 2.00 * DN_DX(1,0) * nu * DN_DX(3,0) + DN_DX(1,1) * nu * DN_DX(3,1) + DN_DX(1,2) * nu * DN_DX(3,2)-(betta*nu-c)*DN_DX(1,0)*DN_DX(3,0);
    p=K(3,9);    
    K(3,10) = DN_DX(1,1) * nu * DN_DX(3,0)-(betta*nu-c)*DN_DX(1,0)*DN_DX(3,1);
    p=K(3,10);   
    K(3,11) = DN_DX(1,2) * nu * DN_DX(3,0)-(betta*nu-c)*DN_DX(1,0)*DN_DX(3,2);
    p=K(3,11);
////////

    K(4,0) = DN_DX(0,1) * nu * DN_DX(1,0)-(betta*nu-c)*DN_DX(1,1)*DN_DX(0,0);
    p=K(4,0);    
    K(4,1) = 2.00 * DN_DX(0,1) * nu * DN_DX(1,1) + DN_DX(0,0) * nu * DN_DX(1,0) + DN_DX(1,2) * nu * DN_DX(0,2)-(betta*nu-c)*DN_DX(1,1)*DN_DX(0,1);
    p=K(4,1);
    K(4,2) = DN_DX(0,1) * nu * DN_DX(1,2)-(betta*nu-c)*DN_DX(1,1)*DN_DX(0,2);
    p=K(4,2); 
    

    K(4,3) = DN_DX(1,1) * nu * DN_DX(1,0)-(betta*nu-c)*DN_DX(1,1)*DN_DX(1,0);
    p=K(4,3);    
    K(4,4) = 2.00 * nu * pow(DN_DX(1,1), 2) + pow(DN_DX(1,0), 2) * nu + nu * pow(DN_DX(1,2), 2)-(betta*nu-c)*DN_DX(1,1)*DN_DX(1,1);
    p=K(4,4);    
    K(4,5) = DN_DX(1,2) * nu * DN_DX(1,1)-(betta*nu-c)*DN_DX(1,1)*DN_DX(1,2);
    p=K(4,5);    
   

   K(4,6) = DN_DX(1,0) * nu * DN_DX(2,1)-(betta*nu-c)*DN_DX(1,1)*DN_DX(2,0);
   p=K(4,6);    
   K(4,7) = 2.00 * DN_DX(1,1) * nu * DN_DX(2,1) + DN_DX(1,0) * nu * DN_DX(2,0) + DN_DX(1,2) * nu * DN_DX(2,2)-(betta*nu-c)*DN_DX(1,1)*DN_DX(2,1);
   p=K(4,7);    
   K(4,8) = DN_DX(1,2) * nu * DN_DX(2,1)-(betta*nu-c)*DN_DX(1,1)*DN_DX(2,2);
   p=K(4,8);
   

   K(4,9) = DN_DX(1,0) * nu * DN_DX(3,1)-(betta*nu-c)*DN_DX(1,1)*DN_DX(3,0);
   p=K(4,9);    
   K(4,10) = 2.00 * DN_DX(1,1) * nu * DN_DX(3,1) + DN_DX(1,0) * nu * DN_DX(3,0) + DN_DX(1,2) * nu * DN_DX(3,2)-(betta*nu-c)*DN_DX(1,1)*DN_DX(3,1);
   p=K(4,10);    
   K(4,11) = DN_DX(1,2) * nu * DN_DX(3,1)-(betta*nu-c)*DN_DX(1,1)*DN_DX(3,2);
   p=K(4,11);       
   
////////////////

   K(5,0) = DN_DX(0,2) * nu * DN_DX(1,0)-(betta*nu-c)*DN_DX(1,2)*DN_DX(0,0);
   p=K(5,0);    
   K(5,1) = DN_DX(0,2) * nu * DN_DX(1,1)-(betta*nu-c)*DN_DX(1,2)*DN_DX(0,1);
   p=K(5,1);
   K(5,2) = 2.00 * DN_DX(0,2) * nu * DN_DX(1,2) + DN_DX(0,0) * nu * DN_DX(1,0) + DN_DX(0,1) * nu * DN_DX(1,1)-(betta*nu-c)*DN_DX(1,2)*DN_DX(0,2);
   p=K(5,2); 


   K(5,3) = DN_DX(1,2) * nu * DN_DX(1,0)-(betta*nu-c)*DN_DX(1,2)*DN_DX(1,0);
   p=K(5,3);    
   K(5,4) = DN_DX(1,2) * nu * DN_DX(1,1)-(betta*nu-c)*DN_DX(1,2)*DN_DX(1,1);
   p=K(5,4);    
   K(5,5) = 2.00 * nu * pow(DN_DX(1,2), 2) + pow(DN_DX(1,0), 2) * nu + nu * pow(DN_DX(1,1), 2)-(betta*nu-c)*DN_DX(1,2)*DN_DX(1,2);
   p=K(5,5);
   

   K(5,6) = DN_DX(1,0) * nu * DN_DX(2,2)-(betta*nu-c)*DN_DX(1,2)*DN_DX(2,0);
   p=K(5,6);    
   K(5,7) = DN_DX(1,1) * nu * DN_DX(2,2)-(betta*nu-c)*DN_DX(1,2)*DN_DX(2,1);
   p=K(5,7);    
   K(5,8) = 2.00 * DN_DX(1,2) * nu * DN_DX(2,2) + DN_DX(1,0) * nu * DN_DX(2,0) + DN_DX(1,1) * nu * DN_DX(2,1)-(betta*nu-c)*DN_DX(1,2)*DN_DX(2,2);
   p=K(5,8);

    
   K(5,9) = DN_DX(1,0) * nu * DN_DX(3,2)-(betta*nu-c)*DN_DX(1,2)*DN_DX(3,0);
   p=K(5,9);    
   K(5,10) = DN_DX(1,1) * nu * DN_DX(3,2)-(betta*nu-c)*DN_DX(1,2)*DN_DX(3,1);
   p=K(5,10);    
   K(5,11) = 2.00 * DN_DX(1,2) * nu * DN_DX(3,2) + DN_DX(1,0) * nu * DN_DX(3,0) + DN_DX(1,1) * nu * DN_DX(3,1)-(betta*nu-c)*DN_DX(1,2)*DN_DX(3,2);
   p=K(5,11);    
//////////////


    K(6,0) = 2.00 * DN_DX(0,0) * nu * DN_DX(2,0) + DN_DX(0,1) * nu * DN_DX(2,1) + DN_DX(0,2) * nu * DN_DX(2,2)-(betta*nu-c)*DN_DX(2,0)*DN_DX(0,0);
    p=K(6,0);    
    K(6,1) = DN_DX(0,0) * nu * DN_DX(2,1)-(betta*nu-c)*DN_DX(2,0)*DN_DX(0,1);
    p=K(6,1);
    K(6,2) = DN_DX(0,0) * nu * DN_DX(2,2)-(betta*nu-c)*DN_DX(2,0)*DN_DX(0,2);
    p=K(6,2);
    

    K(6,3) = 2.00 * DN_DX(1,0) * nu * DN_DX(2,0) + DN_DX(1,1) * nu * DN_DX(2,1) + DN_DX(1,2) * nu * DN_DX(2,2)-(betta*nu-c)*DN_DX(2,0)*DN_DX(1,0);
    p=K(6,3);
    K(6,4) = DN_DX(1,0) * nu * DN_DX(2,1)-(betta*nu-c)*DN_DX(2,0)*DN_DX(1,1);
    p=K(6,4);    
    K(6,5) = DN_DX(1,0) * nu * DN_DX(2,2)-(betta*nu-c)*DN_DX(2,0)*DN_DX(1,2);
    p=K(6,5);
    
    K(6,6) = 2.00 * pow(DN_DX(2,0), 2) * nu + nu * pow(DN_DX(2,1), 2) + nu * pow(DN_DX(2,2), 2)-(betta*nu-c)*DN_DX(2,0)*DN_DX(2,0);
    p=K(6,6);
    K(6,7) = DN_DX(2,1) * nu * DN_DX(2,0)-(betta*nu-c)*DN_DX(2,0)*DN_DX(2,1);
    p=K(6,7);    
    K(6,8) = DN_DX(2,2) * nu * DN_DX(2,0)-(betta*nu-c)*DN_DX(2,0)*DN_DX(2,2);
    p=K(6,8);    
    

    K(6,9) = 2.00 * DN_DX(2,0) * nu * DN_DX(3,0) + DN_DX(2,1) * nu * DN_DX(3,1) + DN_DX(2,2) * nu * DN_DX(3,2)-(betta*nu-c)*DN_DX(2,0)*DN_DX(3,0);
    p=K(6,9);    
    K(6,10) = DN_DX(2,1) * nu * DN_DX(3,0)-(betta*nu-c)*DN_DX(2,0)*DN_DX(3,1);
    p=K(6,10);    
    K(6,11) = DN_DX(2,2) * nu * DN_DX(3,0)-(betta*nu-c)*DN_DX(2,0)*DN_DX(3,2);
    p=K(6,11);    

////////////

    K(7,0) = DN_DX(0,1) * nu * DN_DX(2,0)-(betta*nu-c)*DN_DX(2,1)*DN_DX(0,0);
    p=K(7,0);    
    K(7,1) = 2.00 * DN_DX(0,1) * nu * DN_DX(2,1) + DN_DX(0,0) * nu * DN_DX(2,0) + DN_DX(0,2) * nu * DN_DX(2,2)-(betta*nu-c)*DN_DX(2,1)*DN_DX(0,1);
    p=K(7,1);
    K(7,2) = DN_DX(0,1) * nu * DN_DX(2,2)-(betta*nu-c)*DN_DX(2,1)*DN_DX(0,2);
    p=K(7,2);
    


    K(7,3) = DN_DX(1,1) * nu * DN_DX(2,0)-(betta*nu-c)*DN_DX(2,1)*DN_DX(1,0);
    p=K(7,3);    
    K(7,4) = 2.00 * DN_DX(1,1) * nu * DN_DX(2,1) + DN_DX(1,0) * nu * DN_DX(2,0) + DN_DX(1,2) * nu * DN_DX(2,2)-(betta*nu-c)*DN_DX(2,1)*DN_DX(1,1);
    p=K(7,4);
    K(7,5) = DN_DX(1,1) * nu * DN_DX(2,2)-(betta*nu-c)*DN_DX(2,1)*DN_DX(1,2);
    p=K(7,5);
    

    K(7,6) = DN_DX(2,1) * nu * DN_DX(2,0)-(betta*nu-c)*DN_DX(2,1)*DN_DX(2,0);
    p=K(7,6);    
    K(7,7) = 2.00 * nu * pow(DN_DX(2,1), 2) + pow(DN_DX(2,0), 2) * nu + nu * pow(DN_DX(2,2), 2)-(betta*nu-c)*DN_DX(2,1)*DN_DX(2,1);
    p=K(7,7);    
    K(7,8) = DN_DX(2,2) * nu * DN_DX(2,1)-(betta*nu-c)*DN_DX(2,1)*DN_DX(2,2);
    p=K(7,8);    
    

    K(7,9) = DN_DX(2,0) * nu * DN_DX(3,1)-(betta*nu-c)*DN_DX(2,1)*DN_DX(3,0);
    p=K(7,9);    
    K(7,10) = 2.00 * DN_DX(2,1) * nu * DN_DX(3,1) + DN_DX(2,0) * nu * DN_DX(3,0) + DN_DX(2,2) * nu * DN_DX(3,2)-(betta*nu-c)*DN_DX(2,0)*DN_DX(0,0)-(betta*nu-c)*DN_DX(2,1)  *DN_DX(3,1);
    p=K(7,10);
    K(7,11) = DN_DX(2,2) * nu * DN_DX(3,1)-(betta*nu-c)*DN_DX(2,0)*DN_DX(0,0)-(betta*nu-c)*DN_DX(2,1)*DN_DX(3,2);
    p=K(7,11);     
////////////////

/////////////
    K(8,0) = DN_DX(0,2) * nu * DN_DX(2,0)-(betta*nu-c)*DN_DX(2,2)*DN_DX(0,0);
    p=K(8,0);
    K(8,1) = DN_DX(0,2) * nu * DN_DX(2,1)-(betta*nu-c)*DN_DX(2,2)*DN_DX(0,1);
    p=K(8,1);
    K(8,2) = 2.00 * DN_DX(0,2) * nu * DN_DX(2,2) + DN_DX(0,0) * nu * DN_DX(2,0) + DN_DX(0,1) * nu * DN_DX(2,1)-(betta*nu-c)*DN_DX(2,2)*DN_DX(0,2);
    p=K(8,2);
    

    K(8,3) = DN_DX(1,2) * nu * DN_DX(2,0)-(betta*nu-c)*DN_DX(2,2)*DN_DX(1,0);
    p=K(8,3);    
    K(8,4) = DN_DX(1,2) * nu * DN_DX(2,1)-(betta*nu-c)*DN_DX(2,2)*DN_DX(1,1);
    p=K(8,4);    
    K(8,5) = 2.00 * DN_DX(1,2) * nu * DN_DX(2,2) + DN_DX(1,0) * nu * DN_DX(2,0) + DN_DX(1,1) * nu * DN_DX(2,1)-(betta*nu-c)*DN_DX(2,2)*DN_DX(1,2);
    p=K(8,5);
    


    K(8,6) = DN_DX(2,2) * nu * DN_DX(2,0)-(betta*nu-c)*DN_DX(2,2)*DN_DX(2,0);
    p=K(8,6);    
    K(8,7) = DN_DX(2,2) * nu * DN_DX(2,1)-(betta*nu-c)*DN_DX(2,2)*DN_DX(2,1);
    p=K(8,7);    
    K(8,8) = 2.00 * nu * pow(DN_DX(2,2), 2) + pow(DN_DX(2,0), 2) * nu + nu * pow(DN_DX(2,1), 2)-(betta*nu-c)*DN_DX(2,2)*DN_DX(2,2);
    p=K(8,8);   
    
    K(8,9) = DN_DX(2,0) * nu * DN_DX(3,2)-(betta*nu-c)*DN_DX(2,2)*DN_DX(3,0);
    p=K(8,9);    
    K(8,10) = DN_DX(2,1) * nu * DN_DX(3,2)-(betta*nu-c)*DN_DX(2,2)*DN_DX(3,1);
    p=K(8,10);    
    K(8,11) = 2.00 * DN_DX(2,2) * nu * DN_DX(3,2) + DN_DX(2,0) * nu * DN_DX(3,0) + DN_DX(2,1) * nu * DN_DX(3,1)-(betta*nu-c)*DN_DX(2,2)*DN_DX(3,2);
    p=K(8,11); 
///////////////
         
    K(9,0) = 2.00 * DN_DX(0,0) * nu * DN_DX(3,0) + DN_DX(0,1) * nu * DN_DX(3,1) + DN_DX(0,2) * nu * DN_DX(3,2)-(betta*nu-c)*DN_DX(3,0)*DN_DX(0,0);
    p=K(9,0);    
    K(9,1) = DN_DX(0,0) * nu * DN_DX(3,1)-(betta*nu-c)*DN_DX(3,0)*DN_DX(0,1);
    p=K(9,1);
    K(9,2) = DN_DX(0,0) * nu * DN_DX(3,2)-(betta*nu-c)*DN_DX(3,0)*DN_DX(0,2);
    p=K(9,2); 
    
    K(9,3) = 2.00 * DN_DX(1,0) * nu * DN_DX(3,0) + DN_DX(1,1) * nu * DN_DX(3,1) + DN_DX(1,2) * nu * DN_DX(3,2)-(betta*nu-c)*DN_DX(3,0)*DN_DX(1,0);
    p=K(9,3);    
    K(9,4) = DN_DX(1,0) * nu * DN_DX(3,1)-(betta*nu-c)*DN_DX(3,0)*DN_DX(1,1);
    p=K(9,4);    
    K(9,5) = DN_DX(1,0) * nu * DN_DX(3,2)-(betta*nu-c)*DN_DX(3,0)*DN_DX(1,2);
    p=K(9,5);
    

    K(9,6) = 2.00 * DN_DX(2,0) * nu * DN_DX(3,0) + DN_DX(2,1) * nu * DN_DX(3,1) + DN_DX(2,2) * nu * DN_DX(3,2)-(betta*nu-c)*DN_DX(3,0)*DN_DX(2,0);
    p=K(9,6);    
    K(9,7) = DN_DX(2,0) * nu * DN_DX(3,1)-(betta*nu-c)*DN_DX(3,0)*DN_DX(2,1);
    p=K(9,7);    
    K(9,8) = DN_DX(2,0) * nu * DN_DX(3,2)-(betta*nu-c)*DN_DX(3,0)*DN_DX(2,2);
    p=K(9,8);    
    

    K(9,9) = 2.00 * pow(DN_DX(3,0), 2) * nu + nu * pow(DN_DX(3,1), 2) + nu * pow(DN_DX(3,2), 2)-(betta*nu-c)*DN_DX(3,0)*DN_DX(3,0);
    p=K(9,9);    
    K(9,10) = DN_DX(3,1) * nu * DN_DX(3,0)-(betta*nu-c)*DN_DX(3,0)*DN_DX(3,1);
    p=K(9,10);    
    K(9,11) = DN_DX(3,2) * nu * DN_DX(3,0)-(betta*nu-c)*DN_DX(3,0)*DN_DX(3,2);
    p=K(9,11);       
    

    K(10,0) = DN_DX(0,1) * nu * DN_DX(3,0)-(betta*nu-c)*DN_DX(3,1)*DN_DX(0,0);
    p=K(10,0);    
    K(10,1) = 2.00 * DN_DX(0,1) * nu * DN_DX(3,1) + DN_DX(0,0) * nu * DN_DX(3,0) + DN_DX(0,2) * nu * DN_DX(3,2)-(betta*nu-c)*DN_DX(3,1)*DN_DX(0,1);
    p=K(10,1);
    K(10,2) = DN_DX(0,1) * nu * DN_DX(3,2)-(betta*nu-c)*DN_DX(3,1)*DN_DX(0,2);
    p=K(10,2);    
    

    K(10,3) = DN_DX(1,1) * nu * DN_DX(3,0)-(betta*nu-c)*DN_DX(3,1)*DN_DX(1,0);
    p=K(10,3);    
    K(10,4) = 2.00 * DN_DX(1,1) * nu * DN_DX(3,1) + DN_DX(1,0) * nu * DN_DX(3,0) + DN_DX(1,2) * nu * DN_DX(3,2)-(betta*nu-c)*DN_DX(3,1)*DN_DX(1,1);
    p=K(10,4);    
    K(10,5) = DN_DX(1,1) * nu * DN_DX(3,2)-(betta*nu-c)*DN_DX(3,1)*DN_DX(1,2);
    p=K(10,5);    
    
    K(10,6) = DN_DX(2,1) * nu * DN_DX(3,0)-(betta*nu-c)*DN_DX(3,1)*DN_DX(2,0);
    p=K(10,6);    
    K(10,7) = 2.00 * DN_DX(2,1) * nu * DN_DX(3,1) + DN_DX(2,0) * nu * DN_DX(3,0) + DN_DX(2,2) * nu * DN_DX(3,2)-(betta*nu-c)*DN_DX(3,1)*DN_DX(2,1);
    p=K(10,7);    
    K(10,8) = DN_DX(2,1) * nu * DN_DX(3,2)-(betta*nu-c)*DN_DX(3,1)*DN_DX(2,2);
    p=K(10,8);    
    
    K(10,9) = DN_DX(3,1) * nu * DN_DX(3,0)-(betta*nu-c)*DN_DX(3,1)*DN_DX(3,0);
    p=K(10,9);    
    K(10,10) = 2.00 * nu * pow(DN_DX(3,1), 2) + pow(DN_DX(3,0), 2) * nu + nu * pow(DN_DX(3,2), 2)-(betta*nu-c)*DN_DX(3,1)*DN_DX(3,1);
    p=K(10,10);    
    K(10,11) = DN_DX(3,2) * nu * DN_DX(3,1)-(betta*nu-c)*DN_DX(3,1)*DN_DX(3,2);  
    p=K(10,11); 
////////////

    K(11,0) = DN_DX(0,2) * nu * DN_DX(3,0)-(betta*nu-c)*DN_DX(3,2)*DN_DX(0,0);
    p=K(11,0);
    K(11,1) = DN_DX(0,2) * nu * DN_DX(3,1)-(betta*nu-c)*DN_DX(3,2)*DN_DX(0,1);
    p=K(11,1);
    K(11,2) = 2.00 * DN_DX(0,2) * nu * DN_DX(3,2) + DN_DX(0,0) * nu * DN_DX(3,0) + DN_DX(0,1) * nu * DN_DX(3,1)-(betta*nu-c)*DN_DX(3,2)*DN_DX(0,2);
    p=K(11,2);
    
    K(11,3) = DN_DX(1,2) * nu * DN_DX(3,0)-(betta*nu-c)*DN_DX(3,2)*DN_DX(1,0);
    p=K(11,3);    
    K(11,4) = DN_DX(1,2) * nu * DN_DX(3,1)-(betta*nu-c)*DN_DX(3,2)*DN_DX(1,1);
    p=K(11,4);    
    K(11,5) = 2.00 * DN_DX(1,2) * nu * DN_DX(3,2) + DN_DX(1,0) * nu * DN_DX(3,0) + DN_DX(1,1) * nu * DN_DX(3,1)-(betta*nu-c)*DN_DX(3,2)*DN_DX(1,2);
    p=K(11,5);
    
    K(11,6) = DN_DX(2,2) * nu * DN_DX(3,0)-(betta*nu-c)*DN_DX(3,2)*DN_DX(2,0);
    p=K(11,6);    
    K(11,7) = DN_DX(2,2) * nu * DN_DX(3,1)-(betta*nu-c)*DN_DX(3,2)*DN_DX(2,1);
    p=K(11,7);    
    K(11,8) = 2.00 * DN_DX(2,2) * nu * DN_DX(3,2) + DN_DX(2,0) * nu * DN_DX(3,0) + DN_DX(2,1) * nu * DN_DX(3,1)-(betta*nu-c)*DN_DX(3,2)*DN_DX(2,2);
    p=K(11,8);  
    
    K(11,9) = DN_DX(3,2) * nu * DN_DX(3,0)-(betta*nu-c)*DN_DX(3,2)*DN_DX(3,0);
    p=K(11,9);    
    K(11,10) = DN_DX(3,2) * nu * DN_DX(3,1)-(betta*nu-c)*DN_DX(3,2)*DN_DX(3,1);
    p=K(11,10);    
    K(11,11) = 2.00 * nu * pow(DN_DX(3,2), 2) + pow(DN_DX(3,0), 2) * nu + nu * pow(DN_DX(3,1), 2)-(betta*nu-c)*DN_DX(3,2)*DN_DX(3,2);
    p=K(11,11);*/

    //filling the symmetric part
    /*for(unsigned int i = 1; i<K.size1(); i++)
      for(unsigned int j = 0; j<i; j++)
      K(i,j) = K(j,i);*/
    
  }
  
  //***********************************************************************
  //***********************************************************************
  //performs the Kroneker product of the Reduced Matrix with the identity matrix of 
  //size "dimension" ADDING to the destination matrix
  inline void  QFluid3D::ExpandAndAddReducedMatrix(
			MatrixType& Destination,
			boost::numeric::ublas::bounded_matrix<double,4,4>& ReducedMatrix,
			const unsigned int dimension)
  {
    KRATOS_TRY
      for (unsigned int i=0;i<ReducedMatrix.size2();i++)
	{
	  int rowindex = i*dimension;
	  for (unsigned int j=0;j<ReducedMatrix.size2();j++)
				{
				  unsigned int colindex = j*dimension;
				  for(unsigned int ii=0;ii<dimension;ii++)
				    Destination(rowindex+ii,colindex+ii)+=ReducedMatrix(i,j);
				}
	}
    KRATOS_CATCH("")
      }
  
  //***********************************************************************
  //***********************************************************************
  
	inline double QFluid3D::CalculateH(double Volume)
	{
	  double h = pow(6.00*Volume,0.3333333);
	  
	  //double h = 0.0;
	  //for(int j=0; j<3; j++)
	  //{
	  //	for(int i = j+1; i<4; i++)
	  //	{ 
		//		double l;
		//		l =  pow(GetGeometry()[i].X() - GetGeometry()[j].X()   ,2);
		//		l += pow(GetGeometry()[i].Y() - GetGeometry()[j].Y()   ,2);
		//		l += pow(GetGeometry()[i].Z() - GetGeometry()[j].Z()   ,2);
	  
		//		if(l > h) h = l;
		//	}
		//}
		//h = sqrt(h);
		
	  return h;
	}
  
} // Namespace Kratos



