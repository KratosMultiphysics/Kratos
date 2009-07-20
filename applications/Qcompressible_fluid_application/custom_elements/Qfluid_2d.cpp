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
 
//#define GRADPN_FORM
//#define STOKES

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/Qfluid_2d.h"
#include "utilities/math_utils.h"
#include "Qcompressible_fluid_application.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{
  boost::numeric::ublas::bounded_matrix<double,3,3> QFluid2D::msaux_matrix;
  boost::numeric::ublas::bounded_matrix<double,3,3> QFluid2D::msMassFactors;
  boost::numeric::ublas::bounded_matrix<double,3,3> QFluid2D::msMassC;
  boost::numeric::ublas::bounded_matrix<double,2,2> QFluid2D::msF;
  boost::numeric::ublas::bounded_matrix<double,3,2> QFluid2D::msDN_DX;
  boost::numeric::ublas::bounded_matrix<double,3,2> QFluid2D::msDN_DXF;
  //boost::numeric::ublas::bounded_matrix<double,3,6> QFluid2D::ms_temp;
  array_1d<double,3> QFluid2D::msN; //dimension = number of nodes
  array_1d<double,2> QFluid2D::ms_vel_gauss; //dimesion coincides with space dimension
  array_1d<double,3> QFluid2D::ms_press; 
  array_1d<double,3> QFluid2D::ms_temp_vec_np; //dimension = number of nodes
  array_1d<double,3> QFluid2D::ms_u_DN; //dimension = number of nodes
  
  /*std::vector< Matrix > mInvJ;
  Vector mDetJ;*/
  //************************************************************************************
  //************************************************************************************
  QFluid2D::QFluid2D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
  {		
    //DO NOT ADD DOFS HERE!!!
  }
  
  //************************************************************************************
  //************************************************************************************
  QFluid2D::QFluid2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
  {
        
    //filling the mass factors
    msMassC(0,0) = 1.00/6.00;  msMassC(0,1) = 1.00/12.00; msMassC(0,2) = 1.00/12.00;
    msMassC(1,0) = 1.00/12.00; msMassC(1,1) = 1.00/6.00;  msMassC(1,2) = 1.00/12.00;
    msMassC(2,0) = 1.00/12.00; msMassC(2,1) = 1.00/12.00; msMassC(2,2) = 1.00/6.00;
    
    
    
    
    msMassFactors(0,0) = 1.00/3.00; msMassFactors(0,1) = 0.00;		msMassFactors(0,2) = 0.00;
    msMassFactors(1,0) = 0.00;		msMassFactors(1,1) = 1.00/3.00; msMassFactors(1,2) = 0.00;
    msMassFactors(2,0) = 0.00;		msMassFactors(2,1) = 0.00;		msMassFactors(2,2) = 1.00/3.00;		
   
  }
  
  Element::Pointer QFluid2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
  {
    KRATOS_TRY
      //ddddddddd      
      return Element::Pointer(new QFluid2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
  }
  
  QFluid2D::~QFluid2D()
  {
  }
  
  //************************************************************************************
  //************************************************************************************
  void QFluid2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY
      
      int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];
	  
    if(FractionalStepNumber < 3) //first step of the fractional step solution
      {
	
	Stage1(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
	
      }
    else if (FractionalStepNumber == 4)//second step of the fractional step solution
      {
	//Stage2(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
      }
    
    KRATOS_CATCH("")
      }
  
  //************************************************************************************
  //************************************************************************************
  void QFluid2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_ERROR(std::logic_error,  "method not implemented" , "");
  }
  
  //************************************************************************************
  //************************************************************************************
  //calculation by component of the fractional step velocity corresponding to the first stage
  
  void QFluid2D::Stage1(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY;
    
    unsigned int number_of_points = 3;
    unsigned int dim = 2;
    unsigned int matsize = number_of_points *dim;
    const unsigned int number_of_nodes = GetGeometry().size();    
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
    int j=0;
    if(rLeftHandSideMatrix.size1() != matsize)
      rLeftHandSideMatrix.resize(matsize,matsize,false);
    
    if(rRightHandSideVector.size() != matsize)
      rRightHandSideVector.resize(matsize,false);
    
    //getting data for the given geometry
    double Area, T;
    //CalculateGeometryData(msDN_DX,msN,Area);
    GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);
    
    
    /*mInvJ.resize(integration_points.size());	
      mDetJ.resize(integration_points.size());*/	
    
    //getting the velocity vector on the nodes
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
    const double A0 = GetGeometry()[0].FastGetSolutionStepValue(ARRHENIUS);
    if(f0a==1) j++;

    //if(f0j==1) n++;
    //if(f0i==2)	
    const array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL);
    const array_1d<double,3>& v1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
    double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1);//_OLD_IT);
    double p1oldaux = GetGeometry()[1].FastGetSolutionStepValue(PRESSUREAUX,1);//_OLD_IT);

    const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
    const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
    const double T1 = GetGeometry()[1].FastGetSolutionStepValue(TEMPERATURE);    
    const int f1 = GetGeometry()[1].FastGetSolutionStepValue(IS_FLUID); 
    const int f1i = GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE);
    const double k1 = GetGeometry()[1].FastGetSolutionStepValue(BULK_MODULUS);    
    //const double& f1j = GetGeometry()[1].FastGetSolutionStepValue(FLAG_VARIABLE);   
    //if(f1j==1) n++;
    //if(f1i==2)
    const int f1a = GetGeometry()[1].FastGetSolutionStepValue(FLAG_VARIABLE); 
    const double A1 = GetGeometry()[1].FastGetSolutionStepValue(ARRHENIUS);
    if(f1a==1) j++;

    const array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL);
    const array_1d<double,3>& v2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);
    double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1);//_OLD_IT);
    double p2oldaux = GetGeometry()[2].FastGetSolutionStepValue(PRESSUREAUX,1);//_OLD_IT);
    
    const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
    const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
    const double T2 = GetGeometry()[2].FastGetSolutionStepValue(TEMPERATURE);
    const int f2 = GetGeometry()[2].FastGetSolutionStepValue(IS_FLUID); 
    const int f2i = GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE);    
    const double k2 = GetGeometry()[2].FastGetSolutionStepValue(BULK_MODULUS);
    //const double& f2j = GetGeometry()[2].FastGetSolutionStepValue(FLAG_VARIABLE);  
    const int f2a = GetGeometry()[2].FastGetSolutionStepValue(FLAG_VARIABLE);
    const double A2 = GetGeometry()[2].FastGetSolutionStepValue(ARRHENIUS); 
    if(f2a==1) j++;
    

    //if(f2j==1) n++;
    
    
    //deformation gradient
    noalias(msF)= ZeroMatrix(2,2);
    Matrix kronecker(2,2);
    noalias(kronecker)=ZeroMatrix(2,2);
    
    for(int i=0; i<2;i++){ 
      kronecker(i,i)=1.0;
    }
    //deformation gradient
    for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
	for(int node=0; node<number_of_nodes; node++)
	  {
	    msF(i,j) += (GetGeometry()[node]).GetSolutionStepValue(DESP)(i)	*msDN_DX(node,j);//(DISPLACEMENT)(i)
	    /*double s= (GetGeometry()[node]).GetSolutionStepValue(DESP)(i);
	      KRATOS_WATCH(s);*/				
	  }
    // 
    for(int i=0; i<dim; i++) msF(i,i) = kronecker(i,i);
    /*********************/
    /**************/
    
    //calculating actual jacobian
    /*		GeometryType::JacobiansType J;
		GetGeometry().Jacobian(J); */
      
      double dd=0.0;
      Matrix Jinvv (2,2);
      MathUtils<double>::InvertMatrix(msF/*J[0]*/,Jinvv,dd);
      //MathUtils<double>::InvertMatrix(msF,mInvJ[0],mDetJ[0]);
      if (dd<=0.0)
	{
	  std::cout<<"ERROR........ negative det"<<std::endl;
	  std::cout<<"ERROR........ negative det"<<std::endl;
	  std::cout<<"ERROR........ negative det"<<std::endl;
	  std::cout<<"ERROR........ negative det"<<std::endl;
	  std::cout<<"ERROR........ negative det"<<std::endl;
	  std::cout<<"ERROR........ negative det"<<std::endl;
	  std::cout<<"ERROR........ negative det"<<std::endl;
	}
      
      
      //calculating viscosity
      double nu =0;
      double density =0;
      double K =0;
	    /*nu = 0.333333333333333333333333*(nu0 + nu1 + nu2 );
	    density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );*/
	    K = 0.3333333333333333333333*(k0 + k1 + k2 );

     if(j==3){
	    nu = 0.333333333333333333333333*(nu0 + nu1 + nu2 );
	    density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );
	    K = 0.3333333333333333333333*(k0 + k1 + k2 );
	   

	}
	 else{

		/*if(f0a == 1) {
			GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE)= 1;
 		}

		if(f1a == 1) {
			GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE)= 1;
 		}

		if(f2a == 1) {
			GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE)= 1;
 		}
		*/
		 	/*nu = 0.333333333333333333333333*(nu0 + nu1 + nu2 );
		    	density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );*/
			nu = 0.333333333333333333333333*(nu0 + nu1 + nu2 );
	    		density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );
	    		K = 0.3333333333333333333333*(k0 + k1 + k2 );
		    	/*K = 0.3333333333333333333333*(k0 + k1 + k2 );	*/		
			nu = 0.001;
			density = 1.3; //0.3333333333333333333333*(rho0 + rho1 + rho2 );
			/*K= 100000;*/	
		
		}
	/*nu = 0.333333333333333333333333*(nu0 + nu1 + nu2 );
	density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );
	K = 0.3333333333333333333333*(k0 + k1 + k2 );	*/
	/*nu = 0.333333333333333333333333*(nu0 + nu1 + nu2 );
	density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );
	K = 0.3333333333333333333333*(k0 + k1 + k2 );*/

     const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
    	
    //nu=0;
    //noalias(msDN_DX)=prod(msDN_DX,mInvJ[0]);
    //msDN_DX=prod(msDN_DX,Jinvv);

    double k=K/BDFcoeffs[0];

    CalculateViscousMatrix(rLeftHandSideMatrix, msDN_DX, nu, k,density);
	/*KRATOS_WATCH(rLeftHandSideMatrix);*/
       noalias(msaux_matrix) = BDFcoeffs[0] * msMassFactors * density;
    /*KRATOS_WATCH(msaux_matrix);*/
    //adding all contributions to the stiffness matrix
    ExpandAndAddReducedMatrix(rLeftHandSideMatrix,msaux_matrix,dim);

/*KRATOS_WATCH(rLeftHandSideMatrix);*/  
  /*KRATOS_WATCH(rLeftHandSideMatrix);*/
    //multiplication by the area
    //rLeftHandSideMatrix *= Area*density;
    /*KRATOS_WATCH(rLeftHandSideMatrix);*/
    rLeftHandSideMatrix *= Area /** density*/;
/*KRATOS_WATCH(rLeftHandSideMatrix);*/
    // *****************************************
    
    //CALCULATION OF THE RHS
    const array_1d<double,3>& force0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);
    const array_1d<double,3>& force1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);
    const array_1d<double,3>& force2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);		

    
    array_1d<double,3> rhs_aux;
    for( unsigned int component_index = 0; component_index < dim; component_index++)
      {
	//external forces (component)
	double force_component = 0.3333333333333333*(force0[component_index] + force1[component_index] + force2[component_index]);
	noalias(rhs_aux) = force_component * msN * density/**(1-0.002*(T-273))*/;
	
	//adding pressure gradient (integrated by parts)
	double p_avg=0; double p_arrhenius=0;
	if(j==3){
	p_avg = p0old + p1old + p2old;
	}
	else{
	p_avg = p0oldaux + p1oldaux + p2oldaux;
	}
	
	p_avg *= 0.3333333333333333/*/density*/;
	rhs_aux[0] += msDN_DX(0,component_index)*(p_avg ); 
	rhs_aux[1] += msDN_DX(1,component_index)*(p_avg ); 
	rhs_aux[2] += msDN_DX(2,component_index)*(p_avg );
	
	//KRATOS_WATCH(rhs_aux);
	
	//adding the inertia terms
	// RHS += M*vhistory 
	//calculating the historical velocity
	noalias(ms_temp_vec_np) = ZeroVector(3);		
	for(unsigned int iii = 0; iii<number_of_points; iii++)
	  {
	    const array_1d<double,3>& v = (GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY,1) );
	    ms_temp_vec_np[iii] = BDFcoeffs[0]*v[component_index]*density;
	  }
	
	noalias(rhs_aux) += prod(msMassFactors,ms_temp_vec_np) ;
	//KRATOS_WATCH(rhs_aux);

	//writing the rhs_aux in its place
	for( unsigned int i = 0; i < number_of_points; i++)
	  {
	    rRightHandSideVector[i*dim + component_index] = rhs_aux[i];
	  }
      }
		
    //multiplying by area
   
    rRightHandSideVector *= (Area/* * density*/);
    
    //LHS stabilization contribution //restando los terminos de la interfase
    Vector fvvect(6);
    //KRATOS_WATCH(rRightHandSideVector);

    /*    KRATOS_WATCH(VELOCITY);*/
    
    	/*KRATOS_WATCH(fv0);
	KRATOS_WATCH(fv1);
	KRATOS_WATCH(fv2);
	KRATOS_WATCH(v0);
	KRATOS_WATCH(v1);
	KRATOS_WATCH(v2);

    KRATOS_WATCH(rhs_aux);*/
    
    for( unsigned int component_index = 0; component_index < dim; component_index++)
      {
	fvvect[0 + component_index] = fv0[component_index];
	fvvect[2 + component_index] = fv1[component_index];
	fvvect[4 + component_index] = fv2[component_index];
      }
    /*KRATOS_WATCH(fvvect);
      KRATOS_WATCH(rRightHandSideVector);
    KRATOS_WATCH(rRightHandSideVector);
	KRATOS_WATCH(prod(rLeftHandSideMatrix,fvvect));*/
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,fvvect);
    /*KRATOS_WATCH(rRightHandSideVector);*/
    /*KRATOS_WATCH(this->Id());

/*    KRATOS_WATCH(rLeftHandSideMatrix);
    KRATOS_WATCH(rRightHandSideVector);*/

   /* double nodal_contrib = 0.333333333333333333333333333 * Area;
    GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
    GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
    GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib; */
    KRATOS_CATCH("");
  }
  
  
  
  
  //************************************************************************************
  //************************************************************************************
  //calculation by component of the fractional step velocity corresponding to the first stage
  void QFluid2D::Stage2(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
			ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY;
    
     
   KRATOS_CATCH("");
  }
  
  //************************************************************************************
  //************************************************************************************
  // this subroutine calculates the nodal contributions for the explicit steps of the 
  // fractional step procedure
  void QFluid2D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
  	{
	KRATOS_TRY
	int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];
	  
	//getting data for the given geometry
	double Area;
	GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Area);
	if(FractionalStepNumber == 6) //calculation of velocities
	{
	      
	      //Calculo de la presion
		const unsigned int number_of_nodes = GetGeometry().size();
	      	const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
	      	const unsigned int dim = 2;
	      	const Vector& BDFcoeffs = CurrentProcessInfo[BDF_COEFFICIENTS];
	      	int n=0; int j=0;
		//deformation gradient
 		noalias(msF)= ZeroMatrix(2,2);
 
               	Matrix kronecker(2,2);
               	noalias(kronecker)=ZeroMatrix(2,2);
                                    
                for(int i=0; i<2;i++){ 
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
 		for(int i=0; i<dim; i++) msF(i,i) = kronecker(i,i);
/*********************/

		//athUtils<double>::InvertMatrix(msF,mInvJ[0],mDetJ[0]);
	        double dd=0.0;
		Matrix Jinvv (2,2);
                MathUtils<double>::InvertMatrix(msF/*J[0]*/,Jinvv,dd);
     //MathUtils<double>::InvertMatrix(msF,mInvJ[0],mDetJ[0]);
		if (dd<=0.0)
		{
			std::cout<<"ERROR........ negative det"<<std::endl;
			std::cout<<"ERROR........ negative det"<<std::endl;
			std::cout<<"ERROR........ negative det"<<std::endl;
			std::cout<<"ERROR........ negative det"<<std::endl;
			std::cout<<"ERROR........ negative det"<<std::endl;
			std::cout<<"ERROR........ negative det"<<std::endl;
			std::cout<<"ERROR........ negative det"<<std::endl;
		}
		/*KRATOS_WATCH(msF);*/	
///////////////////
		//msDN_DX=prod(msDN_DX,Jinvv);
///////////////////no considero el jacobiano, lo supongo igual a 1.
	      
	      double& p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
	      double p0n = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1); 
	      double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1);//_OLD_IT);
 	      double p0oldaux = GetGeometry()[0].FastGetSolutionStepValue(PRESSUREAUX,1);//_OLD_IT);

	      array_1d<double,3>& fv0 = GetGeometry()[0].FastGetSolutionStepValue(FRACT_VEL);
	      const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
	      const int f0i = GetGeometry()[0].FastGetSolutionStepValue(IS_STRUCTURE); 
	      const int f0a = GetGeometry()[0].FastGetSolutionStepValue(FLAG_VARIABLE); 
	      if(f0i) n++;
	      if(f0a==1) j++;
	      const double A0 = GetGeometry()[0].FastGetSolutionStepValue(ARRHENIUS);
	      const double A00 = GetGeometry()[0].FastGetSolutionStepValue(ARRHENIUSAUX);
      	      const double k0 = GetGeometry()[0].FastGetSolutionStepValue(BULK_MODULUS);

	      double& p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
	      double p1n = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1); 
	      double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1);//_OLD_IT);
	      double p1oldaux = GetGeometry()[1].FastGetSolutionStepValue(PRESSUREAUX,1);//_OLD_IT);
	      array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL);
	      const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
	      const int f1i = GetGeometry()[1].FastGetSolutionStepValue(IS_STRUCTURE);  
	      const int f1a = GetGeometry()[1].FastGetSolutionStepValue(FLAG_VARIABLE); 
	      if(f1i) n++;
	      if(f1a==1) j++;	
	      const double A1 = GetGeometry()[1].FastGetSolutionStepValue(ARRHENIUS);
      	      const double A11 = GetGeometry()[1].FastGetSolutionStepValue(ARRHENIUSAUX);
       	      const double k1 = GetGeometry()[1].FastGetSolutionStepValue(BULK_MODULUS);

	      double& p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
	      double p2n = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1); 
	      double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1);//_OLD_IT);
	      double p2oldaux = GetGeometry()[2].FastGetSolutionStepValue(PRESSUREAUX,1);//_OLD_IT);
	      array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL);
	      const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
	      const int f2i = GetGeometry()[2].FastGetSolutionStepValue(IS_STRUCTURE);    
	      const int f2a = GetGeometry()[2].FastGetSolutionStepValue(FLAG_VARIABLE); 	     
	      if(f2i) n++;
	      if(f2a==1) j++;	


	      const double A2 = GetGeometry()[2].FastGetSolutionStepValue(ARRHENIUS);
	      const double A22 = GetGeometry()[2].FastGetSolutionStepValue(ARRHENIUSAUX);
      	      const double k2 = GetGeometry()[2].FastGetSolutionStepValue(BULK_MODULUS);
	      double density = 0.33333333333333333333333*(rho0 + rho1 + rho2);
	      
	      double Aver = 0.33333333333333*(A0 + A1 + A2 );
	      double K = 0.33333333333333*(k0 + k1 + k2 );
	      double k=K/BDFcoeffs[0];
	      if(n==3)
		{
			GetGeometry()[0].FastGetSolutionStepValue(NODAL_PRESS) += 0;
			GetGeometry()[1].FastGetSolutionStepValue(NODAL_PRESS) += 0;
         		GetGeometry()[2].FastGetSolutionStepValue(NODAL_PRESS) +=0; 
			GetGeometry()[0].FastGetSolutionStepValue(NODAL_PRESSAUX) +=0; 
			GetGeometry()[1].FastGetSolutionStepValue(NODAL_PRESSAUX) +=0; 
			GetGeometry()[2].FastGetSolutionStepValue(NODAL_PRESSAUX) +=0; 		


 	   		double nodal_contrib = 0.33333333333333 * Area;
 	   		GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) += 0;
 	   		GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS) += 0;
 	   		GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS) += 0;
			GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASSAUX) += 0;
			GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASSAUX) += 0;
			GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASSAUX) += 0;


	   	}
		else			
		{		  
		  	double aux0=0;
		  	double aux1=0;
		  	double aux2=0;
		  	if(j==3)
			{
		  		double Gaux;
		  		//k=100000/BDFcoeffs[0]; /*Aver=0;*/
		  		double p_avg = p0old + p1old + p2old;
	          		p_avg *= 0.33333333333333 ;
		  		Gaux =  msDN_DX(0,0)*fv0[0] + msDN_DX(0,1)*fv0[1]; //LO COMENTE YO
		  		Gaux += msDN_DX(1,0)*fv1[0] + msDN_DX(1,1)*fv1[1];
		  		Gaux += msDN_DX(2,0)*fv2[0] + msDN_DX(2,1)*fv2[1];
		  		aux0 = - Gaux * msN[0] * k + p_avg * msN[0] + Aver * msN[0] * k; 
		  		aux1 = - Gaux * msN[1] * k + p_avg * msN[1] + Aver * msN[1] * k; 
		  		aux2 = - Gaux * msN[2] * k + p_avg * msN[2] + Aver * msN[2] * k; 
		  
		 	 	GetGeometry()[0].FastGetSolutionStepValue(NODAL_PRESS) += aux0 * Area;
		  		GetGeometry()[1].FastGetSolutionStepValue(NODAL_PRESS) += aux1 * Area;
		  		GetGeometry()[2].FastGetSolutionStepValue(NODAL_PRESS) += aux2 * Area;		    	  
		  
		  		double nodal_contrib = 0.33333333333333 * Area;
		  	
		  		GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
		  		GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
		  		GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
  			}
			else
			{
		    	// k=100000/BDFcoeffs[0]; 
				//k=100000/BDFcoeffs[0];
	
				double Gaux;
				double p_avg = p0oldaux + p1oldaux + p2oldaux;
	        		p_avg *= 0.33333333333333 ;
		  		//aqui genera masa

//////////////////////////
 				double Aver0 = 0; double  Aver1 = 0; double Aver2 = 0;
			/*Aver=0;*/	int s=0;
						
				Aver = A00 + A11 + A22; //Aver0 + Aver1 + Aver2;	
			//Aver /= s;
				Aver *= 0.33333333333333 ;
				Aver *= -1;

//////////////////////////

		  		Gaux =  msDN_DX(0,0)*fv0[0] + msDN_DX(0,1)*fv0[1]; //LO COMENTE YO
		  		Gaux += msDN_DX(1,0)*fv1[0] + msDN_DX(1,1)*fv1[1];
		  		Gaux += msDN_DX(2,0)*fv2[0] + msDN_DX(2,1)*fv2[1];
		  		aux0 = - Gaux * msN[0] * k + p_avg * msN[0] + Aver * msN[0] * k; 
		  		aux1 = - Gaux * msN[1] * k + p_avg * msN[1] + Aver * msN[1] * k; 
		  		aux2 = - Gaux * msN[2] * k + p_avg * msN[2] + Aver * msN[2] * k; 
		  
		  		GetGeometry()[0].FastGetSolutionStepValue(NODAL_PRESSAUX) += aux0 * Area;
		  		GetGeometry()[1].FastGetSolutionStepValue(NODAL_PRESSAUX) += aux1 * Area;
		  		GetGeometry()[2].FastGetSolutionStepValue(NODAL_PRESSAUX) += aux2 * Area;		    	  
		  
		  		double nodal_contribp = 0.33333333333333 * Area;
		  
				double SS1 = GetGeometry()[0].FastGetSolutionStepValue(NODAL_PRESSAUX); 
				double SS2 = GetGeometry()[1].FastGetSolutionStepValue(NODAL_PRESSAUX); 
				double SS3 = GetGeometry()[2].FastGetSolutionStepValue(NODAL_PRESSAUX); 
		
				double ll1=GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASSAUX) += nodal_contribp;
				double ll2=GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASSAUX) += nodal_contribp;
				double ll3=GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASSAUX) += nodal_contribp;
			}
		}
		}
		else if(FractionalStepNumber == 5)
		{
	      		const Vector& BDFcoeffs = CurrentProcessInfo[BDF_COEFFICIENTS];
	      		int n=0; int j=0;
		
			double& p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSUREAUX);
			double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSUREAUX,1);//_OLD_IT);
	      		const double rho0old = GetGeometry()[0].FastGetSolutionStepValue(DENSITY,1);
	      		const int f0i = GetGeometry()[0].FastGetSolutionStepValue(IS_STRUCTURE); 
	      		const int f0a = GetGeometry()[0].FastGetSolutionStepValue(FLAG_VARIABLE); 
      	      		const double k0 = GetGeometry()[0].FastGetSolutionStepValue(BULK_MODULUS);
	      		if(f0i) n++;
	      		if(f0a==1) j++;


	      		double& p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSUREAUX);
	        	double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSUREAUX,1);//_OLD_IT);
	  		const double rho1old = GetGeometry()[1].FastGetSolutionStepValue(DENSITY,1);
	        	const double k1 = GetGeometry()[1].FastGetSolutionStepValue(BULK_MODULUS);
	      		const int f1i = GetGeometry()[1].FastGetSolutionStepValue(IS_STRUCTURE);  
	      		const int f1a = GetGeometry()[1].FastGetSolutionStepValue(FLAG_VARIABLE); 
	      		if(f1i) n++;
	      		if(f1a==1) j++;


			double& p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSUREAUX);
	      		double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSUREAUX,1);//_OLD_IT);
	      		const double rho2old = GetGeometry()[2].FastGetSolutionStepValue(DENSITY,1);
	      		const int f2i = GetGeometry()[2].FastGetSolutionStepValue(IS_STRUCTURE);    
	      		const int f2a = GetGeometry()[2].FastGetSolutionStepValue(FLAG_VARIABLE); 
			const double k2 = GetGeometry()[1].FastGetSolutionStepValue(BULK_MODULUS);	     
	      		if(f2i) n++;
	      		if(f2a==1) j++;	


	      		double densityold = 0.33333333333333333333333*(rho0old + rho1old + rho2old);
	      		double p_avg = 0.33333333333333333333333*(p0 + p1 + p2);
			double p_avg_1 = 0.33333333333333333333333*(p0old + p1old + p2old);
	      		double K = 0.33333333333333*(k0 + k1 + k2 );
	      		double k=K/BDFcoeffs[0];
	       		
	
	  		double aux0=0;
		  	double aux1=0;
		  	double aux2=0;
		  	if(j<3){
		  		double Gaux;
		  		k=100000/*/BDFcoeffs[0]*/; /*Aver=0;*/
				aux0 =  (densityold/k) * (p_avg-p_avg_1) * msN[0] + densityold * msN[0]; 
		  		aux1 =  (densityold/k) * (p_avg-p_avg_1) * msN[1] + densityold * msN[1]; 
		  		aux2 =  (densityold/k) * (p_avg-p_avg_1) * msN[2] + densityold * msN[2]; 
		  
		  		GetGeometry()[0].FastGetSolutionStepValue(NODAL_DENSITYAUX) += aux0 * Area;
		  		GetGeometry()[1].FastGetSolutionStepValue(NODAL_DENSITYAUX) += aux1 * Area;
		  		GetGeometry()[2].FastGetSolutionStepValue(NODAL_DENSITYAUX) += aux2 * Area;		    	  
				}
		}
	KRATOS_CATCH("");
}  
  //************************************************************************************
  //************************************************************************************
  void QFluid2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
  {
    
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dim = 2;
    unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];
    
    if(FractionalStepNumber == 1) //step 1
      {
	if(rResult.size() != number_of_nodes*dim)
	  rResult.resize(number_of_nodes*dim,false);	
	
	for (unsigned int i=0;i<number_of_nodes;i++)
	  {
	    rResult[i*dim] = GetGeometry()[i].GetDof(FRACT_VEL_X).EquationId();
	    rResult[i*dim+1] = GetGeometry()[i].GetDof(FRACT_VEL_Y).EquationId();
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
  void QFluid2D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
  {
    
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dim = 2;
    
    unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];
    
    if(FractionalStepNumber == 1) //step 1
      {
	if(ElementalDofList.size() != number_of_nodes*dim)
	  ElementalDofList.resize(number_of_nodes*dim);	
	
	for (unsigned int i=0;i<number_of_nodes;i++)
	  {
	    ElementalDofList[i*dim] = GetGeometry()[i].pGetDof(FRACT_VEL_X);
	    ElementalDofList[i*dim+1] = GetGeometry()[i].pGetDof(FRACT_VEL_Y);
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
/*	inline void Fluid2D::CalculateGeometryData(boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, array_1d<double,3>& N, double& Area)
	{
		double det_J = GetGeometry()[1].X() * GetGeometry()[2].Y()
		- GetGeometry()[1].X() * GetGeometry()[0].Y()
		- GetGeometry()[0].X() * GetGeometry()[2].Y()
		- GetGeometry()[1].Y() * GetGeomtanque_12.601.post.binetry()[2].X() 
		+ GetGeometry()[1].Y() * GetGeometry()[0].X() 
		+ GetGeometry()[0].Y() * GetGeometry()[2].X();

		double temp = 1.00/det_J;
		
		Area = fabs(det_J) * 0.5;
		
		// X derivatives
		DN_DX(0, 0) = (GetGeometry()[1].Y() - GetGeometry()[2].Y()) *temp;	
		DN_DX(1, 0) = (GetGeometry()[2].Y() - GetGeometry()[0].Y()) *temp;
		DN_DX(2, 0) = (GetGeometry()[0].Y() - GetGeometry()[1].Y()) *temp;
		// Y derivatives
		DN_DX(0, 1) = (GetGeometry()[2].X() - GetGeometry()[1].X()) *temp;
		DN_DX(1, 1) = (GetGeometry()[0].X() - GetGeometry()[2].X()) *temp;
		DN_DX(2, 1) = (GetGeometry()[1].X() - GetGeometry()[0].X()) *temp;

		//shape functions
		N[0] = 0.333333333333333333333333333;
		N[1] = 0.333333333333333333333333333;
		N[2] = 0.333333333333333333333333333;
	}
*/
//#undef GRADPN_FORM

  void QFluid2D::CalculateViscousMatrix(MatrixType& K, const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double& nu, const double& k, const double& density)
  {

	int dim=2;
	Matrix B(3,6);
	Matrix ms_temp(3,6);
	noalias(ms_temp) = ZeroMatrix(3,6);
	noalias(B) = ZeroMatrix(3,6);
	
		for (unsigned int i=0;i<3;i++)
		{
			unsigned int index = dim*i;
			B(0,index+0)=DN_DX(i,0);					B(0,index+1)= 0.0;
			B(1,index+0)=0.0;						B(1,index+1)= DN_DX(i,1);
			B(2,index+0)= DN_DX(i,1);					B(2,index+1)= DN_DX(i,0); 
		}			
				
		Matrix ms_constitutive_matrix(3,3);
		ms_constitutive_matrix = ZeroMatrix(3,3);
		//constitutive tensor
		ms_constitutive_matrix(0,0) = k/2+(4.0/3.0)*nu;	ms_constitutive_matrix(0,1) = k/2 -2.0/3.0*nu;	ms_constitutive_matrix(0,2) = 0.0;
		ms_constitutive_matrix(1,0) = k/2-2.0/3.0*nu; 	ms_constitutive_matrix(1,1) = k/2+4.0/3.0*nu;	ms_constitutive_matrix(1,2) = 0.0;
		ms_constitutive_matrix(2,0) = 0.0;	        ms_constitutive_matrix(2,1) = 0.0;	ms_constitutive_matrix(2,2) = 0.0+nu;
			
		//calculating viscous contributions
		ms_temp = prod( ms_constitutive_matrix , B);
		noalias(K) = prod( trans(B) , ms_temp);
		


				
		
		
    /*double betta=0.666;
    
    K(0,0) = 2.0 * pow(DN_DX(0,0), 2) * nu + nu * pow(DN_DX(0,1), 2)-(betta*nu-k)*DN_DX(0,0)*DN_DX(0,0); //k=K(0,0);
    K(0,1) = DN_DX(0,1) * nu * DN_DX(0,0)-(betta*nu-k)*DN_DX(0,0)*DN_DX(0,1);//k=K(0,1);
    
    K(0,2) = 2.0 * DN_DX(0,0) * nu * DN_DX(1,0) + DN_DX(0,1) * nu * DN_DX(1,1)-(betta*nu-k)*DN_DX(0,0)*DN_DX(1,0);//k=K(0,2);
    K(0,3) = DN_DX(0,1) * nu * DN_DX(1,0)-(betta*nu-k)*DN_DX(0,0)*DN_DX(1,1);//k=K(0,3);
    K(0,4) = 2.0 * DN_DX(0,0) * nu * DN_DX(2,0) + DN_DX(0,1) * nu * DN_DX(2,1)-(betta*nu-k)*DN_DX(0,0)*DN_DX(2,0);//k=K(0,4);
    K(0,5) = DN_DX(0,1) * nu * DN_DX(2,0)-(betta*nu-k)*DN_DX(0,0)*DN_DX(2,1);//k=K(0,5);
    K(1,0) = DN_DX(0,1) * nu * DN_DX(0,0)-(betta*nu-k)*DN_DX(0,1)*DN_DX(0,0);
    K(1,1) = 2.0 * nu * pow(DN_DX(0,1), 2) + pow(DN_DX(0,0), 2) * nu-(betta*nu-k)*DN_DX(0,1)*DN_DX(0,1);//k=K(1,1);
    K(1,2) = DN_DX(0,0) * nu * DN_DX(1,1)-(betta*nu-k)*DN_DX(0,1)*DN_DX(1,0);//k=K(1,2);
    K(1,3) = 2.0 * DN_DX(0,1) * nu * DN_DX(1,1) + DN_DX(0,0) * nu * DN_DX(1,0)-(betta*nu-k)*DN_DX(0,1)*DN_DX(1,1);//k=K(1,3);
    K(1,4) = DN_DX(0,0) * nu * DN_DX(2,1)-(betta*nu-k)*DN_DX(0,1)*DN_DX(2,0);//k=K(1,4);
    K(1,5) = 2.0 * DN_DX(0,1) * nu * DN_DX(2,1) + DN_DX(0,0) * nu * DN_DX(2,0)-(betta*nu-k)*DN_DX(0,1)*DN_DX(2,1);//k=K(1,5);
    K(2,0) = 2.0 * DN_DX(0,0) * nu * DN_DX(1,0) + DN_DX(0,1) * nu * DN_DX(1,1)-(betta*nu-k)*DN_DX(1,0)*DN_DX(0,0);
    K(2,1) = DN_DX(0,0) * nu * DN_DX(1,1)-(betta*nu-k)*DN_DX(1,0)*DN_DX(0,1);
    K(2,2) = 2.0 * pow(DN_DX(1,0), 2) * nu + nu * pow(DN_DX(1,1), 2)-(betta*nu-k)*DN_DX(1,0)*DN_DX(1,0);//k=K(2,2);
    K(2,3) = DN_DX(1,1) * nu * DN_DX(1,0)-(betta*nu-k)*DN_DX(1,0)*DN_DX(1,1);//k=K(2,3);
    K(2,4) = 2.0 * DN_DX(1,0) * nu * DN_DX(2,0) + DN_DX(1,1) * nu * DN_DX(2,1)-(betta*nu-k)*DN_DX(1,0)*DN_DX(2,0);//k=K(2,4);
    K(2,5) = DN_DX(1,1) * nu * DN_DX(2,0)-(betta*nu-k)*DN_DX(1,0)*DN_DX(2,1);//k=K(2,5);
    K(3,0) = DN_DX(0,1) * nu * DN_DX(1,0)-(betta*nu-k)*DN_DX(1,1)*DN_DX(0,0);
    K(3,1) = 2.0 * DN_DX(0,1) * nu * DN_DX(1,1) + DN_DX(0,0) * nu * DN_DX(1,0)-(betta*nu-k)*DN_DX(1,1)*DN_DX(0,1);
    K(3,2) = DN_DX(1,1) * nu * DN_DX(1,0)-(betta*nu-k)*DN_DX(1,1)*DN_DX(1,0);
    K(3,3) = 2.0 * nu * pow(DN_DX(1,1), 2) + pow(DN_DX(1,0), 2) * nu-(betta*nu-k)*DN_DX(1,1)*DN_DX(1,1);//k=K(3,3);
    K(3,4) = DN_DX(1,0) * nu * DN_DX(2,1)-(betta*nu-k)*DN_DX(1,1)*DN_DX(2,0);//k=K(3,4);
    K(3,5) = 2.0 * DN_DX(1,1) * nu * DN_DX(2,1) + DN_DX(1,0) * nu * DN_DX(2,0)-(betta*nu-k)*DN_DX(1,1)*DN_DX(2,1);//k=K(3,5);
    K(4,0) = 2.0 * DN_DX(0,0) * nu * DN_DX(2,0) + DN_DX(0,1) * nu * DN_DX(2,1)-(betta*nu-k)*DN_DX(2,0)*DN_DX(0,0);
    K(4,1) = DN_DX(0,0) * nu * DN_DX(2,1)-(betta*nu-k)*DN_DX(2,0)*DN_DX(0,1);
    K(4,2) = 2.0 * DN_DX(1,0) * nu * DN_DX(2,0) + DN_DX(1,1) * nu * DN_DX(2,1)-(betta*nu-k)*DN_DX(2,0)*DN_DX(1,0);
    K(4,3) = DN_DX(1,0) * nu * DN_DX(2,1)-(betta*nu-k)*DN_DX(2,0)*DN_DX(1,1);
    
    K(4,4) = 2.0 * pow(DN_DX(2,0), 2) * nu + nu * pow(DN_DX(2,1), 2)-(betta*nu-k)*DN_DX(2,0)*DN_DX(2,0);//k=K(4,4);
    K(4,5) = DN_DX(2,1) * nu * DN_DX(2,0)-(betta*nu-k)*DN_DX(2,0)*DN_DX(2,1);//k=K(4,5);

    K(5,0) = DN_DX(0,1) * nu * DN_DX(2,0)-(betta*nu-k)*DN_DX(2,1)*DN_DX(0,0);
    
    K(5,1) = 2.0 * DN_DX(0,1) * nu * DN_DX(2,1) + DN_DX(0,0) * nu * DN_DX(2,0)-(betta*nu-k)*DN_DX(2,1)*DN_DX(0,1);
    K(5,2) = DN_DX(1,1) * nu * DN_DX(2,0)-(betta*nu-k)*DN_DX(2,1)*DN_DX(1,0);
    K(5,3) = 2.0 * DN_DX(1,1) * nu * DN_DX(2,1) + DN_DX(1,0) * nu * DN_DX(2,0)-(betta*nu-k)*DN_DX(2,1)*DN_DX(1,1);
    K(5,4) = DN_DX(2,1) * nu * DN_DX(2,0)-(betta*nu-k)*DN_DX(2,1)*DN_DX(2,0);
    K(5,5) = 2.0 * nu * pow(DN_DX(2,1), 2) + pow(DN_DX(2,0), 2) * nu-(betta*nu-k)*DN_DX(2,1)*DN_DX(2,1);//k=K(5,5);*/
    
/********************************/
    /*K(0,0) *= (1/density);
    K(0,1) *= (1/density); 
    K(0,2) *= (1/density);
    K(0,3) *= (1/density);
    K(0,4) *= (1/density); 
    K(0,5) *= (1/density); 
    K(1,0) *= (1/density);
    K(1,1) *= (1/density);
    K(1,2) *= (1/density); 
    K(1,3) *= (1/density);
    K(1,4) *= (1/density);
    K(1,5) *= (1/density); 
    K(2,0) *= (1/density); 
    K(2,1) *= (1/density);
    K(2,2) *= (1/density);
    K(2,3) *= (1/density);
    K(2,4) *= (1/density); 
    K(2,5) *= (1/density);
    K(3,0) *= (1/density); 
    K(3,1) *= (1/density); 
    K(3,2) *= (1/density);
    K(3,3) *= (1/density);
    K(3,4) *= (1/density);
    K(3,5) *= (1/density); 
    K(4,0) *= (1/density); 
    K(4,1) *= (1/density);
    K(4,2) *= (1/density);
    K(4,3) *= (1/density); 
    K(4,4) *= (1/density);
    K(4,5) *= (1/density);
    K(5,0) *= (1/density); 
    K(5,1) *= (1/density);
    K(5,2) *= (1/density);
    K(5,3) *= (1/density); 
    K(5,4) *= (1/density); 
    K(5,5) *= (1/density); */


    //filling the symmetric part
    /*for(unsigned int i = 1; i<K.size1(); i++)
      for(unsigned int j = 0; j<i; j++)
      K(i,j) = K(j,i);*/
  }
  
  //***********************************************************************
  //***********************************************************************
  //performs the Kroneker product of the Reduced Matrix with the identity matrix of 
  //size "dimension" ADDING to the destination matrix
  inline void  QFluid2D::ExpandAndAddReducedMatrix( MatrixType& Destination,boost::numeric::ublas::bounded_matrix<double,3,3>& ReducedMatrix, const unsigned int dimension)
  {
    KRATOS_TRY
      for (unsigned int i=0;i<ReducedMatrix.size2();i++)
	{
	  int rowindex = i*dimension;
	  for (unsigned int j=0;j<ReducedMatrix.size2();j++)
	    {
	      unsigned int colindex = j*dimension;
	    double p=0;
	    for(unsigned int ii=0;ii<dimension;ii++){
	      
	      Destination(rowindex+ii,colindex+ii)+=ReducedMatrix(i,j);
	      p+=ReducedMatrix(i,j);
	    }
	  }
      }
  KRATOS_CATCH("")
    }

} // Namespace Kratos







