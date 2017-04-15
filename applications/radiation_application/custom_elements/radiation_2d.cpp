/*
==============================================================================
KratosConvectionDiffusionApplication 
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
the following condition:boost::numeric::ublas::bounded_matrix<double,6,2> msShapeFunc = ZeroMatrix(6,2);
	#pragma omp threadprivate(msShapeFunc)

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
//   Last modified by:    $Author: julio.marti $
//   Date:                $Date:  $
//   Revision:            $Revision: $
//
//
 

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/radiation_2d.h"
#include "radiation_application.h"
#include "includes/radiation_settings.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 

#define P1
namespace Kratos
{
  
  
  //************************************************************************************
  //************************************************************************************
  Rad2D::Rad2D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
  {		
    //DO NOT ADD DOFS HERE!!!
  }
  
  //************************************************************************************
  //************************************************************************************
  Rad2D::Rad2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
  {
    
  }
  
  Element::Pointer Rad2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
	  return Element::Pointer(new Rad2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}
  
  Rad2D::~Rad2D()
  {
  }
  
  //************************************************************************************
  Rad2D::IntegrationMethod Rad2D::GetIntegrationMethod1()
  {
    return mThisIntegrationMethod;
  }
  //************************************************************************************
  void Rad2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

      const unsigned int number_of_points = GetGeometry().size();
    const double lumping_factor = 1.00/double(number_of_points);
    //unsigned int TDim = 2;
    
    
    //KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");	
    mThisIntegrationMethod= GeometryData::GI_GAUSS_2;
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);
    
    
    
    if(rLeftHandSideMatrix.size1() != number_of_points)
      rLeftHandSideMatrix.resize(number_of_points,number_of_points,false);
    
    if(rRightHandSideVector.size() != number_of_points)
      rRightHandSideVector.resize(number_of_points,false);

    
    boost::numeric::ublas::bounded_matrix<double,3,3> msMassFactors = 1.0/3.0*IdentityMatrix(3,3);
    boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX = ZeroMatrix(3,2);
    array_1d<double,3> msN = ZeroVector(3); //dimension = number of nodes
    array_1d<double,2> ms_vel_gauss = ZeroVector(2); //dimesion coincides with space dimension
    array_1d<double,3> ms_temp_vec_np = ZeroVector(3); //dimension = number of nodes
    array_1d<double,3> ms_u_DN = ZeroVector(3); //dimension = number of nodes
    array_1d<double,2> grad_g = ZeroVector(2); //dimesion coincides with space dimension
    boost::numeric::ublas::bounded_matrix<double,2,2> Identity = 1.0*IdentityMatrix(2,2);
    boost::numeric::ublas::bounded_matrix<double,2,2> First = ZeroMatrix(2,2);
    boost::numeric::ublas::bounded_matrix<double,2,2> Second = ZeroMatrix(2,2);
    boost::numeric::ublas::bounded_matrix<double,2,3> Third = ZeroMatrix(2,3);
    boost::numeric::ublas::bounded_matrix<double,1,3> msConvOp = ZeroMatrix(1,3);
    boost::numeric::ublas::bounded_matrix<double,3,1> msShapeFunc = ZeroMatrix(3,1);
    boost::numeric::ublas::bounded_matrix<double,3,3> msAuxMat = ZeroMatrix(3,3);
    array_1d<double,3> msAuxVec = ZeroVector(3); //dimension = number of nodes
    boost::numeric::ublas::bounded_matrix<double,3,3> Tres = ZeroMatrix(3,3);
    
    
    
    //getting data for the given geometry
    double Area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);
    
    RadiationSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(RADIATION_SETTINGS);
    
    //const Variable<double>&  rDensityVar = my_settings->GetDensityVariable();
    
    const Variable<double>& rUnknownVar= my_settings->GetUnknownVariable();
    
    //const Variable<array_1d<double,3> >& rMeshVelocityVar =my_settings->GetMeshVelocityVariable();
    
    
    const double StefenBoltzmann =5.67e-8;  //1.0
    double absorptioncoefficient = 15.0;//0.01;// 0.0001; 0.01;// 960.5;//3000.0;//1000.0;//1000.0;//0.01;// 0.01;// 1000.0; //1000.0; //3000.0 //0.01
    
    
    
    
    double counti=0.0;
    for(unsigned int i = 0; i<number_of_points; i++) 
      {
	
	if(GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE)>0.5) counti++;
      }
    
    
    if(counti==3.00) 
	{
	  absorptioncoefficient = 0.0;	 
	  //absorptioncoefficient = 0.1;	
	}
    
    else
      {
	//KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
	absorptioncoefficient = 100.0; //100.0
      }
    
    
    double temperature ;
    ///*const*/ array_1d<double,3> w = GetGeometry()[0].FastGetSolutionStepValue(rMeshVelocityVar);
    
    //ms_vel_gauss=w;

    //absorptioncoefficient = 0.001;	
    double 	conductivity;
    conductivity=1.0/(3.0*absorptioncoefficient); 	
    

    double c1 = 4.00;
    //double c2 = 2.00;
    double h = sqrt(2.00*Area);
    //double norm_u =norm_2(ms_vel_gauss);
    //double tau1=( h)/( c2 * norm_u );
    //h=0.1;
    double tau1=1.00/( c1 * conductivity /(h * h ) + absorptioncoefficient );
    tau1=1.00/(  c1 * conductivity /(h * h ) + absorptioncoefficient);
    tau1 *= 0.0;	
    double pi=3.1416;
    //pi=1.0;

    //const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
    
    //CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
    //noalias(ms_u_DN) = prod(msDN_DX , ms_vel_gauss);
    
    //noalias(rLeftHandSideMatrix) =  outer_prod(msN,ms_u_DN) * Area;
    //////KRATOS_WATCH(rLeftHandSideMatrix);
    // tau1=1.0;
    //CONVECTION STABILIZING CONTRIBUTION (Suu)
    //noalias(rLeftHandSideMatrix) += tau1 * outer_prod(ms_u_DN,ms_u_DN) * Area;

    //filling the mass factors
    msMassFactors(0,0) = 1.00/3.00;     msMassFactors(0,1) = 0.00;		msMassFactors(0,2) = 0.00;
    msMassFactors(1,0) = 0.00;		msMassFactors(1,1) = 1.00/3.00;         msMassFactors(1,2) = 0.00;
    msMassFactors(2,0) = 0.00;		msMassFactors(2,1) = 0.00;		msMassFactors(2,2) = 1.00/3.00;
    
    double T0 = GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);
    double T1 = GetGeometry()[1].FastGetSolutionStepValue(TEMPERATURE);
    double T2 = GetGeometry()[2].FastGetSolutionStepValue(TEMPERATURE);
    
    //rLeftHandSideMatrix *=0.0;
    //rRightHandSideVector *=0.0;

    //noalias(rLeftHandSideMatrix) += absorptioncoefficient * msMassFactors * Area;

    //noalias(rLeftHandSideMatrix) = conductivity * prod(msDN_DX,trans(msDN_DX)) * Area;

    #if defined(P1)
    rLeftHandSideMatrix *=0.0;
    rRightHandSideVector *=0.0; 
    noalias(rLeftHandSideMatrix) = absorptioncoefficient * msMassFactors * Area;

    noalias(rLeftHandSideMatrix) = conductivity * prod(msDN_DX,trans(msDN_DX)) * Area;	 	

    #endif	

    //double e=2.7182;
    double t_gauss=0.0;
    //double tgauss=0.0;
    int nodes_number=3;
    for (unsigned int i=0;i<nodes_number;i++) t_gauss += msN[i]*GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE);	
    //tgauss=1e(tgauss); 
    //tgauss=pow(e,t_gauss);



    //if(T0>1600.0) T0=1600.0;
    //if(T1>1600.0) T1=1600.0;
    //if(T2>1600.0) T2=1600.0;
    msAuxVec[0]=pow(T0,4);  
    msAuxVec[1]=pow(T1,4);
    msAuxVec[2]=pow(T2,4);

    temperature = msAuxVec[0] + msAuxVec[1] + msAuxVec[2];
    temperature *=lumping_factor;
    
    double tempp= T0 + T1 + T2;
    tempp *=lumping_factor;
    //double tmedia=pow(tempp,4);


    //calculate convective term	
    //int nodes_number = 3;
    /*	
    int dof = 1;
    for (int ii = 0; ii< nodes_number; ii++)
    {
    int column = ii*dof;
    msConvOp(0,column) = msDN_DX(ii,0)*ms_vel_gauss[0] + msDN_DX(ii,1)*ms_vel_gauss[1];
    msShapeFunc(column,0) = msN[ii];
    }
    noalias(msAuxMat) = prod(msShapeFunc, msConvOp);*/
    
    array_1d<double, 3 > a_dot_grad;
    noalias(a_dot_grad) = prod(msDN_DX, ms_vel_gauss);
    array_1d<double, 3 > a_dot_grad_and_mass;
    a_dot_grad_and_mass = msN;
    noalias(msAuxMat)= outer_prod(a_dot_grad, a_dot_grad_and_mass);
    
    noalias(rLeftHandSideMatrix) += 1.0 * tau1 * absorptioncoefficient * msAuxMat * Area; 
    
    //    noalias(rRightHandSideVector)= tau1 * absorptioncoefficient * StefenBoltzmann * (1.0/pi) * prod(msAuxMat, msAuxVec) * Area * 0.0;//0
    
    //msShapeFunc = ZeroMatrix(3,1);
    //msAuxMat = ZeroMatrix(3,3);
    //noalias(rRightHandSideVector)=ZeroMatrix(1,3);
    rRightHandSideVector *=0.0;
    mInvJ0.resize(integration_points.size());
    mDetJ0.resize(integration_points.size(),false);
    
    GeometryType::JacobiansType J0;
    J0 = GetGeometry().Jacobian(J0, mThisIntegrationMethod);
    const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);
    
    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);
    
    
    for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
      {
    //getting informations for integration
    //double IntegrationWeight = integration_points[PointNumber].Weight();
    
    
    MathUtils<double>::InvertMatrix(J0[PointNumber],mInvJ0[PointNumber],mDetJ0[PointNumber]);				
    
    
    double Weight = integration_points[PointNumber].Weight()* mDetJ0[PointNumber];
    const Vector& N=row(Ncontainer,PointNumber);
    
    Tres = ZeroMatrix(3,3); 
    Tres(0,0)=N[0]*N[0];
    Tres(0,0)+=N[0]*N[1];
    Tres(0,0)+=N[0]*N[2];
    Tres(1,1)=N[1]*N[0];
    Tres(1,1)+=N[1]*N[1];
    Tres(1,1)+=N[1]*N[2];
    Tres(2,2)=N[2]*N[0];
    Tres(2,2)+=N[2]*N[1];
    Tres(2,2)+=N[2]*N[2];
    
    
    noalias(msDN_DX) = prod(DN_De[PointNumber],mInvJ0[PointNumber]);
    
    noalias(a_dot_grad) = prod(msDN_DX, ms_vel_gauss);
    a_dot_grad_and_mass = N;
    noalias(msAuxMat)= outer_prod(a_dot_grad, a_dot_grad_and_mass);
    
    double t_gauss=0.0;
    //double e=2.7182;
    for (unsigned int i=0;i<nodes_number;i++) t_gauss += N[i]*GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE);	
    //double tgauss = pow(t_gauss,4);
    
    
    //noalias(msAuxMat) = prod(msShapeFunc, msConvOp);
    
    noalias(rLeftHandSideMatrix) += 1.0 * Weight * absorptioncoefficient * Tres; 
    //rRightHandSideVector += 4.0 * Weight * absorptioncoefficient * StefenBoltzmann*  prod(Tres,msAuxVec ) ;
    double cc=1.0;
#if defined(P1)
    
    cc=4.0*pi;
#endif
    rRightHandSideVector += Weight * prod(Tres,msAuxVec ) *  absorptioncoefficient * StefenBoltzmann * (1.0/pi) * cc ; 
    
    
    rRightHandSideVector += 1.0 * 1.0 * (tau1 * Weight *  absorptioncoefficient * StefenBoltzmann * (1.0/pi)* prod(msAuxMat,msAuxVec )) * cc;
    
    
    
  }
    
    double constant=1.0/(4-2*1.0);
  array_1d<double, 2 > qrad=ZeroVector(2);
  array_1d<double,2>  interface_segment=ZeroVector(2);
  array_1d<double,2>  normaledge1=ZeroVector(2);
  
#if defined(P1)
  
  if((GetGeometry()[1].FastGetSolutionStepValue(IS_BOUNDARY)>0.5 && GetGeometry()[0].FastGetSolutionStepValue(IS_BOUNDARY)>0.5))  //IS_INTERFACE
    {
      double norm=0.0;
      
      interface_segment[0] = (GetGeometry()[0].X()-GetGeometry()[1].X());
      interface_segment[1] = (GetGeometry()[0].Y()-GetGeometry()[1].Y());
      norm = sqrt(  pow((interface_segment[0]),2) + pow((interface_segment[1]),2));
      double area1=norm;
      normaledge1(0)= -interface_segment[1]/norm;
      normaledge1(1)= interface_segment[0]/norm;
      
      
      T0=GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE); //if(T0>1600.0) T0=1600.0;
      T1=GetGeometry()[1].FastGetSolutionStepValue(TEMPERATURE); //if(T1>1600.0) T1=1600.0;
      
      rRightHandSideVector[0] +=  constant * 4.0 * StefenBoltzmann*(pow(T0,4) ) * 0.5 * area1;
      rRightHandSideVector[1] +=  constant * 4.0 * StefenBoltzmann*(pow(T1,4) ) * 0.5 * area1;	
      
      rLeftHandSideMatrix(0,0) += constant * 0.5 * area1;
      rLeftHandSideMatrix(1,1) += constant * 0.5 * area1;
      
    }
  
  if((GetGeometry()[1].FastGetSolutionStepValue(IS_BOUNDARY)>0.5 && GetGeometry()[2].FastGetSolutionStepValue(IS_BOUNDARY)>0.5))
    {
      
      double norm=0.0;
      
      interface_segment[0] = (GetGeometry()[1].X()-GetGeometry()[2].X());
      interface_segment[1] = (GetGeometry()[1].Y()-GetGeometry()[2].Y());
      norm = sqrt(  pow((interface_segment[0]),2) + pow((interface_segment[1]),2));
      double area1=norm;
      normaledge1(0)= -interface_segment[1]/norm;
      normaledge1(1)= interface_segment[0]/norm;
      
      T1=GetGeometry()[1].FastGetSolutionStepValue(TEMPERATURE); //if(T1>1600.0) T1=1600.0;
      T2=GetGeometry()[2].FastGetSolutionStepValue(TEMPERATURE); //if(T2>1600.0) T2=1600.0;

      rRightHandSideVector[1] +=  constant * 4.0 * StefenBoltzmann*(pow(T1,4) ) * 0.5 * area1;
      rRightHandSideVector[2] +=  constant * 4.0 * StefenBoltzmann*(pow(T2,4) ) * 0.5 * area1;	
      
      rLeftHandSideMatrix(1,1) += constant * 0.5 * area1;
      rLeftHandSideMatrix(2,2) += constant * 0.5 * area1;
      
    }
  
  
  if((GetGeometry()[2].FastGetSolutionStepValue(IS_BOUNDARY)>0.5 && GetGeometry()[0].FastGetSolutionStepValue(IS_BOUNDARY)>0.5))
    {
      double norm=0.0;
      
      interface_segment[0] = (GetGeometry()[2].X()-GetGeometry()[0].X());
      interface_segment[1] = (GetGeometry()[2].Y()-GetGeometry()[0].Y());
      norm = sqrt(  pow((interface_segment[0]),2) + pow((interface_segment[1]),2));
      double area1=norm;
      normaledge1(0)= -interface_segment[1]/norm;
      normaledge1(1)= interface_segment[0]/norm;
      
      
      T0=GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE); //if(T0>1600.0) T0=1600.0;
      T2=GetGeometry()[2].FastGetSolutionStepValue(TEMPERATURE); //if(T2>1600.0) T2=1600.0;
      
      rRightHandSideVector[0] +=  constant * 4.0 * StefenBoltzmann*(pow(T0,4) ) * 0.5 * area1;
      rRightHandSideVector[2] +=  constant * 4.0 * StefenBoltzmann*(pow(T2,4) ) * 0.5 * area1;	
      
      rLeftHandSideMatrix(0,0) += constant * 0.5 * area1;
      rLeftHandSideMatrix(2,2) += constant * 0.5 * area1;
      
      
    }
  
#endif
  
  
  for(unsigned int iii = 0; iii<number_of_points; iii++)
    ms_temp_vec_np[iii] = GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar); //TEMPERATURE
  noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp_vec_np);
  
  KRATOS_CATCH("");
  }
  
  //************************************************************************************
  //************************************************************************************
  void Rad2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
  }	 
  
  //************************************************************************************
  //************************************************************************************
  // this subroutine calculates the nodal contributions for the explicit steps of the 
  // fractional step procedure
  void Rad2D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
  {
    KRATOS_TRY
      /*int FractionalStpiepNumber = CurrentProcessInfo[FRACTIONAL_STEP];
	
      //getting data for the given geometry
      double Area;
      GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);		
      
      RadiationSettings::Pointer my_settings = CurrentProcessInfo.GetValue(RADIATION_SETTINGS);
      const Variable<double>& rUnknownVar= my_settings->GetUnknownVariable();
      const Variable<array_1d<double,3> >& rConvectionVar = my_settings->GetConvectionVariable();
      const Variable<array_1d<double,3> >& rMeshVelocityVar =my_settings->GetMeshVelocityVariable();
      
      if(FractionalStepNumber  == 2) //calculation of temperature convective projection
      {
      const unsigned int number_of_points = GetGeometry().size();
      const double lumping_factor = 1.00/double(number_of_points);
      unsigned int TDim = 2;
      
      //calculating viscosity
      ms_temp_vec_np[0] = GetGeometry()[0].FastGetSolutionStepValue(rUnknownVar);	//TEMP		
      const array_1d<double,3>& v = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY); //
      const array_1d<double,3>& w = GetGeometry()[0].FastGetSolutionStepValue(rMeshVelocityVar); //
      for(unsigned int j = 0; j<TDim; j++)
      ms_vel_gauss[j] = v[j] - w[j];
      
      
      for(unsigned int i = 1; i<number_of_points; i++)
      {
      ms_temp_vec_np[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);		//TEMP		
      const array_1d<double,3>& v = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
      const array_1d<double,3>& w = GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar);
      for(unsigned int j = 0; j<TDim; j++)
      ms_vel_gauss[j] += v[j] - w[j];
				
      }
      ms_vel_gauss *= lumping_factor;
      
      //calculating convective auxiliary vector
      noalias(ms_u_DN) = prod(msDN_DX , ms_vel_gauss);
      double temp_conv = inner_prod( ms_u_DN , ms_temp_vec_np);
      temp_conv *= Area;
      
      for(unsigned int i = 0; i<number_of_points; i++)
      {
      GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += lumping_factor*Area;
      GetGeometry()[i].FastGetSolutionStepValue(TEMP_CONV_PROJ) += lumping_factor*temp_conv; ;
      }
      }*/
      KRATOS_CATCH("");
  }

  
  //************************************************************************************
  //************************************************************************************
  void Rad2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
  {
    
    RadiationSettings::Pointer my_settings = CurrentProcessInfo.GetValue(RADIATION_SETTINGS);
    const Variable<double>& rUnknownVar= my_settings->GetUnknownVariable();
    
    
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if(rResult.size() != number_of_nodes)
      rResult.resize(number_of_nodes,false);	
    
    for (unsigned int i=0;i<number_of_nodes;i++)
      rResult[i] = GetGeometry()[i].GetDof(rUnknownVar).EquationId();//TEMP 
  }

//************************************************************************************
//************************************************************************************
void Rad2D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
{
  unsigned int number_of_nodes = GetGeometry().PointsNumber();
  
  RadiationSettings::Pointer my_settings = CurrentProcessInfo.GetValue(RADIATION_SETTINGS);
  const Variable<double>& rUnknownVar= my_settings->GetUnknownVariable();
  
  
  if(ElementalDofList.size() != number_of_nodes)
    ElementalDofList.resize(number_of_nodes);	
  
  for (unsigned int i=0;i<number_of_nodes;i++)
    ElementalDofList[i] = GetGeometry()[i].pGetDof(rUnknownVar);//TEMP
  
}

} // Namespace Kratos



