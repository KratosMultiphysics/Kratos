// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___ 
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _| 
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/conv_diff_3d.h"
#include "convection_diffusion_application.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "includes/deprecated_variables.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************

ConvDiff3D::ConvDiff3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

ConvDiff3D::ConvDiff3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

Element::Pointer ConvDiff3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<ConvDiff3D>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer ConvDiff3D::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<ConvDiff3D>(NewId, pGeom, pProperties);
}

ConvDiff3D::~ConvDiff3D()
{
}

//************************************************************************************
//************************************************************************************

void ConvDiff3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
   KRATOS_TRY

    const unsigned int number_of_points = GetGeometry().size();
  const double lumping_factor = 1.00 / double(number_of_points);
  unsigned int TDim = 3;
  mThisIntegrationMethod= GeometryData::GI_GAUSS_2;
    
  boost::numeric::ublas::bounded_matrix<double, 4, 4 > msMassFactors = 0.25 * IdentityMatrix(4, 4);
  boost::numeric::ublas::bounded_matrix<double, 4, 3 > msDN_DX;
  boost::numeric::ublas::bounded_matrix<double, 4, 3 > msDN_DX1;
  array_1d<double, 4 > msN;
  array_1d<double, 4 > N;
  array_1d<double, 3 > ms_vel_gauss;
  array_1d<double, 4 > ms_temp_vec_np;
  array_1d<double, 4 > ms_u_DN;
  double Area=0.0;	

  const unsigned int LocalSize = GetGeometry().size();
  unsigned int TNumNodes = GetGeometry().size();
  boost::numeric::ublas::bounded_matrix<double,4, 3 > coords;	

  
  if (rLeftHandSideMatrix.size1() != number_of_points)
    rLeftHandSideMatrix.resize(number_of_points, number_of_points, false);
    
  if (rRightHandSideVector.size() != number_of_points)
    rRightHandSideVector.resize(number_of_points, false);
    
    
  boost::numeric::ublas::bounded_matrix<double, 3, 3 > First = ZeroMatrix(3, 3);
  boost::numeric::ublas::bounded_matrix<double, 3, 3 > Second = ZeroMatrix(3, 3);
  boost::numeric::ublas::bounded_matrix<double, 3, 4 > Third = ZeroMatrix(3, 4);
  boost::numeric::ublas::bounded_matrix<double, 3, 3 > Identity = 1.0 * IdentityMatrix(3, 3);
    
  array_1d<double, 3 > grad_g = ZeroVector(3); //dimesion coincides with space dimension
    
  array_1d<double,3> aux_var= ZeroVector(3); //dimesion coincides with space dimension
  boost::numeric::ublas::bounded_matrix<double,4,1> msShapeFunc = ZeroMatrix(4,1);
  boost::numeric::ublas::bounded_matrix<double,1,4> msConvOp = ZeroMatrix(1,4);
  array_1d<double,12> msAuxVec = ZeroVector(4); //dimension = number of nodes
  array_1d<double,12> msAuxVec1 = ZeroVector(4); //dimension = number of nodes
  boost::numeric::ublas::bounded_matrix<double,4,4> msAuxMat = ZeroMatrix(4,4);
  boost::numeric::ublas::bounded_matrix<double,4,4> Tres = ZeroMatrix(4,4);
  boost::numeric::ublas::bounded_matrix<double,4,4> msResta = IdentityMatrix(4,4);
  boost::numeric::ublas::bounded_matrix<double, 4, 4 > Aux = ZeroMatrix(4, 4);
  boost::numeric::ublas::bounded_matrix<double, 4, 4 > other_matrix = ZeroMatrix(4,4);
    
  rLeftHandSideMatrix *= 0.0;
  rRightHandSideVector *= 0.0;
    
  double Dt = rCurrentProcessInfo[DELTA_TIME];

  GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);
  ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);

  const Variable<double>& rDensityVar = my_settings->GetDensityVariable();
  const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
  const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
  const Variable<double>& rSourceVar = my_settings->GetVolumeSourceVariable();
  const Variable<array_1d<double, 3 > >& rMeshVelocityVar = my_settings->GetMeshVelocityVariable();
  const Variable<double>& rProjectionVariable = my_settings->GetProjectionVariable();
  const Variable<array_1d<double, 3 > >& rVelocityVar = my_settings->GetVelocityVariable();
  const Variable<double>& rSpecificHeatVar = my_settings->GetSpecificHeatVariable();

  double conductivity = GetGeometry()[0].FastGetSolutionStepValue(rDiffusionVar);
  double density = GetGeometry()[0].FastGetSolutionStepValue(rDensityVar);
  double specific_heat = GetGeometry()[0].FastGetSolutionStepValue(rSpecificHeatVar);	
  double heat_flux = GetGeometry()[0].FastGetSolutionStepValue(rSourceVar);
  double proj = GetGeometry()[0].FastGetSolutionStepValue(rProjectionVariable);
  const array_1d<double, 3 > & v= GetGeometry()[0].FastGetSolutionStepValue(rVelocityVar); 
  const array_1d<double, 3 > & w = GetGeometry()[0].FastGetSolutionStepValue(rMeshVelocityVar);
  for (unsigned int j = 0; j < TDim; j++)
    ms_vel_gauss[j] = v[j] - w[j];

  for (unsigned int i = 1; i < number_of_points; i++)
    {
      conductivity += GetGeometry()[i].FastGetSolutionStepValue(rDiffusionVar);
      density += GetGeometry()[i].FastGetSolutionStepValue(rDensityVar);
      specific_heat += GetGeometry()[i].FastGetSolutionStepValue(rSpecificHeatVar);
      heat_flux += GetGeometry()[i].FastGetSolutionStepValue(rSourceVar);
      proj += GetGeometry()[i].FastGetSolutionStepValue(rProjectionVariable);
	
      const array_1d<double, 3 > & v = GetGeometry()[i].FastGetSolutionStepValue(rVelocityVar);
      const array_1d<double, 3 > & w = GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar);
      for (unsigned int j = 0; j < TDim; j++)
	ms_vel_gauss[j] += v[j] - w[j];
    }
  conductivity *= lumping_factor;
  density *= lumping_factor;
  specific_heat *= lumping_factor;
  heat_flux *= lumping_factor;
  proj *= lumping_factor;
  ms_vel_gauss *= lumping_factor;
  //ms_vel_gauss *= 0.0;
    
  const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
  //calculating parameter tau
  double c10 = 4.00;
  double c20 = 2.00;
  double h = pow(6.00 * Area, 0.3333333);
  double norm_uv = norm_2(ms_vel_gauss);
    
  if(rUnknownVar == TEMPERATURE or rUnknownVar == YCH4) 
    {
      array_1d<double,4> distances;
      bool has_negative_node=false;
      bool has_positive_node=false;
      bool has_distance_0=false;
      //for radiation
      double absorptioncoefficient2;
      const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);
	
      for(unsigned int iii = 0; iii<4; iii++)
	{
	  distances[iii] = this->GetGeometry()[iii].FastGetSolutionStepValue(DISTANCE);  
	}
	
      bool split_element=false;

      
      for(unsigned int iii = 0; iii<4; iii++)
	{
	  if (distances[iii]<0.0)
	    has_negative_node=true;
	  else
	    has_positive_node=true;		
	}
      
      
      if (has_positive_node && has_negative_node) split_element=true;
      double norm_grad;
      double res;
      double k_aux;
      
      split_element=false;
      if (split_element==true){
	
	
	KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "")
	
       
      }  //si el elemento estaba cortado!!!!
      else
	{
	  
	  double norm_grad;
      	  double res;
          double k_aux;
	  int coi=0;
	  double absorptioncoefficient=0.0;
	  ///NUEVVOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
	  
	  int suma=0;	 int sumai=0;
	  int s0, s1, s2, s3;
	  s0=0;s1=0;s2=0; s3=0; 	
	  int elem=0;	
	  int i0, i1,i2,i3;
	  int ii0, ii1,ii2,ii3;
	  has_negative_node=false;
	  has_positive_node=false;
          
	  
	  if(GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE)==1.0) s0=1; 	//IS_INTERFACE
	  if(GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE)==1.0) s1=1;
	  if(GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE)==1.0) s2=1;
	  if(GetGeometry()[3].FastGetSolutionStepValue(IS_INTERFACE)==1.0) s3=1;

	  if(GetGeometry()[0].FastGetSolutionStepValue(DISTANCE)==-1.0) s0=1; 	//IS_INTERFACE
	  if(GetGeometry()[1].FastGetSolutionStepValue(DISTANCE)==-1.0) s1=1;
	  if(GetGeometry()[2].FastGetSolutionStepValue(DISTANCE)==-1.0) s2=1;
	  if(GetGeometry()[3].FastGetSolutionStepValue(DISTANCE)==-1.0) s3=1;


	  elem=s0 + s1 + s2 + s3; 

	  int so_heat_flux=0;
          
	  if(elem==4) has_negative_node=true;

	  else has_negative_node=false;
	  
	  if (has_negative_node== true) 
	    {

	      //KRATOS_WATCH("eeeeeeeeeeeeestoy aqui");

	      //KRATOS_THROW_ERROR(std::logic_error, "method not implemented", ""); 
  
	      conductivity=0.16;  //0.16  //16
	      density=905.0;           
	      specific_heat=1900.0;
	      so_heat_flux=0;	
	      conductivity=0.12;  //0.16  //16
	
	      conductivity=0.06;
	      conductivity=0.12;
              density=0.0;
	      density=	GetGeometry()[0].FastGetSolutionStepValue(DENSITY) + GetGeometry()[1].FastGetSolutionStepValue(DENSITY)+ GetGeometry()[2].FastGetSolutionStepValue(DENSITY)+ GetGeometry()[3].FastGetSolutionStepValue(DENSITY); 
	      density *=0.25;

              conductivity=0.0;
              conductivity=GetGeometry()[0].FastGetSolutionStepValue(CONDUCTIVITY) + GetGeometry()[1].FastGetSolutionStepValue(CONDUCTIVITY)+ GetGeometry()[2].FastGetSolutionStepValue(CONDUCTIVITY)+ GetGeometry()[3].FastGetSolutionStepValue(CONDUCTIVITY); 
	      conductivity *=0.25;

              specific_heat=0.0;

              specific_heat=GetGeometry()[0].FastGetSolutionStepValue(SPECIFIC_HEAT) + GetGeometry()[1].FastGetSolutionStepValue(SPECIFIC_HEAT)+ GetGeometry()[2].FastGetSolutionStepValue(SPECIFIC_HEAT)+ GetGeometry()[3].FastGetSolutionStepValue(SPECIFIC_HEAT); 
	      specific_heat *=0.25;


	      
	    } 
	  else 
	    {
	      conductivity=0.0131;//0.0131;//0.05;//1.0; 0.05
	      density=1.0;//1000.0;
	      specific_heat=1310.0;//620.0;//1.0;
	      so_heat_flux=1;	
	      /*absorptioncoefficient=75.0;
              absorptioncoefficient=75.0;
              absorptioncoefficient=75.0;
              absorptioncoefficient=150.0;
              absorptioncoefficient=75.0;

              absorptioncoefficient=750.0;
              absorptioncoefficient=375.0;
              absorptioncoefficient=200.0;
              absorptioncoefficient=3000.0;
              absorptioncoefficient=1500.0;
              absorptioncoefficient=1000.0;*/
              absorptioncoefficient=1000.0;
              absorptioncoefficient=10.0;
              absorptioncoefficient=100.0;
	      //KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
	      density=0.0;
	      conductivity=0.0;
	      specific_heat=0.0;
  	      for (unsigned int i = 1; i < number_of_points; i++)
    	      {
	      //KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
              density += GetGeometry()[i].FastGetSolutionStepValue(DENSITY);
              specific_heat+= GetGeometry()[i].FastGetSolutionStepValue(SPECIFIC_HEAT);
              conductivity+= GetGeometry()[i].FastGetSolutionStepValue(CONDUCTIVITY);
    	      }
  	      density *= lumping_factor;
              specific_heat*= lumping_factor;
              conductivity*= lumping_factor;
	    }
	  
	 
	  
	  
	  const array_1d<double, 3 > & v= GetGeometry()[0].FastGetSolutionStepValue(rVelocityVar); 
	  const array_1d<double, 3 > & w = GetGeometry()[0].FastGetSolutionStepValue(rMeshVelocityVar); 
	  
	  for (unsigned int j = 0; j < TDim; j++)
            ms_vel_gauss[j] = v[j] - w[j];
	  
	  for (unsigned int i = 1; i < number_of_points; i++)
	    {
	      const array_1d<double, 3 > & v = GetGeometry()[i].FastGetSolutionStepValue(rVelocityVar);	
	      const array_1d<double, 3 > & w = GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar);
	      for (unsigned int j = 0; j < TDim; j++)
		ms_vel_gauss[j] += v[j] - w[j];
	    }
	  ms_vel_gauss *= lumping_factor;
	  
	  
	  const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
	  double c1 = 4.00;
	  double c2 = 2.00;
	  double h = sqrt(2.00 * Area);
	  double norm_u = norm_2(ms_vel_gauss);
	  double tau1;
	  
	  tau1 = (h * h) / (density * specific_heat * BDFcoeffs[0] * h * h + c1 * conductivity + c2 * density * specific_heat * (norm_u + 1e-6) * h);
	  
	  
	  for (unsigned int i = 0; i < TDim; i++)
	    {
	      for (unsigned int j = 0; j < number_of_points; j++)
		{
		  grad_g[i] += msDN_DX(j, i) * GetGeometry()[j].FastGetSolutionStepValue(rUnknownVar);
		}
	    }
	  
	  
	  h = pow(6.00 * Area, 0.3333333);
	  /*double */norm_u = norm_2(ms_vel_gauss);
	  res = density * specific_heat*(inner_prod(ms_vel_gauss,grad_g)) ;//+ 0.333333333333333 * (t0media+t1media+t2media)*(1/dt)*density*conductivity;
	  norm_grad=norm_2(grad_g);
	  k_aux=fabs(res) /(norm_grad + 0.000000000001);
	  
	  //noalias(First) = outer_prod(ms_vel_gauss, trans(ms_vel_gauss));
	  //First /= ((norm_u + 1e-6)*(norm_u + 1e-6));
	  //noalias(Second) = Identity - First;
	  //noalias(Third) = prod(Second, trans(msDN_DX));
	    
	  array_1d<double, 4 > a_dot_grad;
	  noalias(a_dot_grad) = prod(msDN_DX, ms_vel_gauss);
	  array_1d<double, 4 > a_dot_grad_and_mass;
	  a_dot_grad_and_mass = msN;
	  noalias(msAuxMat)= outer_prod(a_dot_grad, a_dot_grad_and_mass);
	  
          noalias(ms_u_DN) = prod(msDN_DX, ms_vel_gauss);


         
	  //KRATOS_WATCH(msAuxMat);
	  boost::numeric::ublas::bounded_matrix<double, 4, 4 > Aux = ZeroMatrix(4, 4);
	  
	  rLeftHandSideMatrix *=0.0;  
	  rRightHandSideVector *=0.0;
      
	  
	  noalias(rLeftHandSideMatrix) = (conductivity + so_heat_flux * k_aux * h) * prod(msDN_DX, trans(msDN_DX))* Area;  //0.0 * k_aux * h PARA PAPER RADIACION
	    
	  noalias(rLeftHandSideMatrix) += 1.0/Dt * (density * specific_heat) * msMassFactors * Area;	 
	  
	  noalias(ms_u_DN) = prod(msDN_DX, ms_vel_gauss);
	  
	  noalias(rLeftHandSideMatrix) += so_heat_flux * (density * specific_heat) * outer_prod(msN, ms_u_DN)* Area;
	  
	  noalias(rLeftHandSideMatrix) += so_heat_flux * density * specific_heat * tau1 * outer_prod(ms_u_DN, ms_u_DN)* Area;
	  
	  noalias(rLeftHandSideMatrix) += so_heat_flux * 1.0 / Dt * tau1 * msAuxMat * (density *specific_heat)* Area;
	  
 
	  
	  //  noalias(rRightHandSideVector) = (heat_flux * 1.0) * msN * Area *1.0 ;// * 0.0;	
	  
	  double yo=0.0;
	  double yf=0.0;
	  for (unsigned int j = 0; j < 4; j++) if (GetGeometry()[j].FastGetSolutionStepValue(YO)>0.0) yo += GetGeometry()[j].FastGetSolutionStepValue(YO);
	  yo *=0.25;
	  
	  
	  for (unsigned int j = 0; j < 4; j++) if (GetGeometry()[j].FastGetSolutionStepValue(YF)>0.0) yf += GetGeometry()[j].FastGetSolutionStepValue(YF);
	  yf *=0.25;
	  
	  
	  double temperature=  GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE)+GetGeometry()[1].FastGetSolutionStepValue(TEMPERATURE)+GetGeometry()[2].FastGetSolutionStepValue(TEMPERATURE)+GetGeometry()[3].FastGetSolutionStepValue(TEMPERATURE);
	  
	  temperature *=0.25;
	  
	  if(temperature>2600.0) temperature=2600.0;
	  //noalias(rLeftHandSideMatrix) +=(+1.0) *  26010.0 *1000.0 *  yf * yo * 6e9 * exp(-10700.0/(temperature)) * 10700.0 /(temperature * temperature) * msMassFactors * Area;
	  
	  //noalias(rRightHandSideVector) += (1.0) * msN * 26010.0 *1000.0 *  yf * yo * 6e9 * exp(-10700.0/(temperature)) * Area * so_heat_flux; //ULTIMO
	  //rRightHandSideVector *=0.0;
	  
	  //KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
	  //SOLO EN EL AIRE
	  //KRATOS_WATCH(elem);
	  if(elem<4)
	    {
	      double unknown = density * density * 26010.0 *1000.0 *  yf * yo * 6e9 * exp(-10700.0/(temperature))*0.70*1.0; // 1.0 //ULTIMO //250.0 //12desept double unknown = 26010.0 *1000.0 *  yf * yo * 6e9 * exp(-10700.0/(temperature))*0.7*1.2;  //last 
	      
	      if (unknown>1e9) //1e8
		{
		  
	  	  rRightHandSideVector[0] += 1.0 * msN[0] * 1e9 * Area * so_heat_flux;  ///1e8
		  rRightHandSideVector[1] += 1.0 * msN[1] * 1e9 * Area * so_heat_flux;
		  rRightHandSideVector[2] += 1.0 * msN[2] * 1e9 * Area * so_heat_flux;
		  rRightHandSideVector[3] += 1.0 * msN[3] * 1e9 * Area * so_heat_flux;

	  	  rRightHandSideVector[0] += 1.0 * msN[0] * Area * GetGeometry()[0].FastGetSolutionStepValue(HEAT_FLUX)* so_heat_flux * 1.0; /// a cero por el otro ejemplo
		  rRightHandSideVector[1] += 1.0 * msN[1] * Area * GetGeometry()[1].FastGetSolutionStepValue(HEAT_FLUX)* so_heat_flux * 1.0 ;
		  rRightHandSideVector[2] += 1.0 * msN[2] * Area * GetGeometry()[2].FastGetSolutionStepValue(HEAT_FLUX)* so_heat_flux * 1.0 ;
		  rRightHandSideVector[3] += 1.0 * msN[3] * Area * GetGeometry()[3].FastGetSolutionStepValue(HEAT_FLUX)* so_heat_flux * 1.0 ;



		  //KRATOS_WATCH("SE PASA");
		  //KRATOS_WATCH("SE PASA");
		 // KRATOS_WATCH(unknown);
		  //KRATOS_WATCH("temperatureeeeeeeeeeeeee");

		  //KRATOS_WATCH(temperature);
		}
	      else
		{
		  
		  rRightHandSideVector[0] += density * density * 1.0 * msN[0] * 0.70 * 1.0* 26010.0 *1000.0 * yf * yo * 6e9 * exp(-10700.0/(temperature)) * Area * so_heat_flux;
		  rRightHandSideVector[1] += density * density * 1.0 * msN[1] * 0.70 * 1.0* 26010.0 *1000.0 * yf * yo * 6e9 * exp(-10700.0/(temperature)) * Area * so_heat_flux;
		  rRightHandSideVector[2] += density * density * 1.0 * msN[2] * 0.70 * 1.0* 26010.0 *1000.0 * yf * yo * 6e9 * exp(-10700.0/(temperature)) * Area * so_heat_flux;
		  rRightHandSideVector[3] += density * density * 1.0 * msN[3] * 0.70 * 1.0* 26010.0 *1000.0 * yf * yo * 6e9 * exp(-10700.0/(temperature)) * Area * so_heat_flux; // rRightHandSideVector[3] += 1.0 * msN[3] * 0.7 * 1.2* 26010.0 *1000.0 * yf * yo * 6e9 * exp(-10700.0/(temperature)) * Area * so_heat_flux;


		  rRightHandSideVector[0] += 1.0 * msN[0] * Area * GetGeometry()[0].FastGetSolutionStepValue(HEAT_FLUX)* so_heat_flux * 1.0 ;
		  rRightHandSideVector[1] += 1.0 * msN[1] * Area * GetGeometry()[1].FastGetSolutionStepValue(HEAT_FLUX)* so_heat_flux * 1.0 ;
		  rRightHandSideVector[2] += 1.0 * msN[2] * Area * GetGeometry()[2].FastGetSolutionStepValue(HEAT_FLUX)* so_heat_flux * 1.0;
		  rRightHandSideVector[3] += 1.0 * msN[3] * Area * GetGeometry()[3].FastGetSolutionStepValue(HEAT_FLUX)* so_heat_flux * 1.0;

		}
	      
	    }	  
	  double StefenBoltzmann = 5.67e-8;
	  double emissivity = 1.0;
	  double aux = pow(298.0,4);
	  double convection_coefficient=0.0;
	  noalias(Aux)=rLeftHandSideMatrix;
	  
	  msAuxVec[0]=GetGeometry()[0].FastGetSolutionStepValue(rUnknownVar,1);
	  msAuxVec[1]=GetGeometry()[1].FastGetSolutionStepValue(rUnknownVar,1);
	  msAuxVec[2]=GetGeometry()[2].FastGetSolutionStepValue(rUnknownVar,1);
	  msAuxVec[3]=GetGeometry()[3].FastGetSolutionStepValue(rUnknownVar,1);
	  
	  noalias(rRightHandSideVector) += so_heat_flux * 1.0/Dt * tau1 * density* specific_heat * prod(msAuxMat, msAuxVec) * Area;

	  
	  for (unsigned int iii = 0; iii < number_of_points; iii++)
	    ms_temp_vec_np[iii] = -1.0/Dt * GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar, 1);
	  
	  noalias(rRightHandSideVector) -= 1.0 * prod(msMassFactors, ms_temp_vec_np * density * specific_heat* Area);
	  
	  //subtracting the dirichlet term
	  // RHS -= LHS*temperatures
	  
	  ///////////////////////////////////////
	  ///////////////////////////////////////
	  
	  double area=0.0;
	  double T0,T1,T2,T3;
	  suma=0.0;
	  //	  int i0,i1,i2,i3;
	  std::vector<array_1d<double,3> > PointsOfFSTriangle;
	  PointsOfFSTriangle.reserve(3);
	  const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);//GetGeometr
	  //double q55=  GetGeometry()[0].GSolutionStepValue(FACE_HEAT_FLUX);


	  if(has_negative_node== true){ 
	  

		
		double t1 = GetGeometry()[0].FastGetSolutionStepValue(YCH4);

		double t2 = GetGeometry()[1].FastGetSolutionStepValue(YCH4);

		double t3 = GetGeometry()[2].FastGetSolutionStepValue(YCH4);

		double t4 = GetGeometry()[3].FastGetSolutionStepValue(YCH4);

		double temp=t1 + t2 + t3 + t4;
		temp *= 0.25;
		//if(temp>1000.0) temp=1000.0;
	        if(temp>823.15) temp=823.15; //temp=900.0

		double E_over_R = 21877.25; // 24400.0;//28961.49;pp pure
	    	double C = 1.93e12;//2.18e12;//1.19e15; pp pure
                
                double heatofcombustion=800000.0;
		//TPU CABLES
		//E_over_R = 2411256.31;//28961.49;
                //heatofcombustion=26300000.0;
		double auxterm=exp(-E_over_R/(temp));
                //C=10000.0;

		/*TPU CABLES
		constexpr double E_over_R1 = 18626.00;//tpu07 17819.01;//tpu009  18626.0;//17819.01;//11256.31;//24400.0;//28961.49;
    		constexpr double C1 = 1.0e11;//tpu07 1.0e10;//tpu009 1.0e11;//1.0e10;// 10000.0;//2.18e12;//1.19e15;

		constexpr double E_over_R2 = 20720.81;//tpu07 25904.93;//tpu009 20720.81;//25904.93;//11256.31;//24400.0;//28961.49;
    		constexpr double C2 = 1.0e10;//tpu07 1.0e13;//tpu009 1.0e10;//1.0e13;// 10000.0;//2.18e12;//1.19e15;
		*/
		/* TPU CABLES
		double taux=0.0;
		if(temp>542.0) taux=542.0;
		else taux=temp;
	 
		double aux_var= C * exp(-E_over_R/(taux));
	
		double t1aux=0.0;
		if(temp>742.0) t1aux=742.0;
		else t1aux=temp;

        	double aux_var1= C1 * exp(-E_over_R1/(taux));
        	double aux_var2= C2 * exp(-E_over_R2/(t1aux));*/

		rRightHandSideVector[0] += 0.0 * (-1.0) * msN[0] * density * heatofcombustion *  C * auxterm * Area;
		rRightHandSideVector[1] += 0.0 * (-1.0) * msN[1] * density * heatofcombustion *  C * auxterm * Area;
		rRightHandSideVector[2] += 0.0 * (-1.0) * msN[2] * density * heatofcombustion *  C * auxterm * Area;
		rRightHandSideVector[3] += 0.0 * (-1.0) * msN[3] * density * heatofcombustion *  C * auxterm * Area;
		/* TPU CABLES
		rRightHandSideVector[0] += (-1.0) * msN[0] * density * heatofcombustion *  aux_var1 * Area;
		rRightHandSideVector[1] += (-1.0) * msN[1] * density * heatofcombustion *  aux_var1 * Area;
		rRightHandSideVector[2] += (-1.0) * msN[2] * density * heatofcombustion *  aux_var1 * Area;
		rRightHandSideVector[3] += (-1.0) * msN[3] * density * heatofcombustion *  aux_var1 * Area;

		rRightHandSideVector[0] += (-1.0) * msN[0] * density * heatofcombustion *  aux_var2 * Area;
		rRightHandSideVector[1] += (-1.0) * msN[1] * density * heatofcombustion *  aux_var2 * Area;
		rRightHandSideVector[2] += (-1.0) * msN[2] * density * heatofcombustion *  aux_var2 * Area;
		rRightHandSideVector[3] += (-1.0) * msN[3] * density * heatofcombustion *  aux_var2 * Area;*/
		//KRATOS_WATCH(800000.0 *  C * exp(-E_over_R/(temp)));
	        double Time = rCurrentProcessInfo[TIME];
	  
	      double qq=0.0;	      
	      i0=0;
	      i1=0;
	      i2=0;
	      i3=0;	
	      
	      suma=0;
	      if(GetGeometry()[1].FastGetSolutionStepValue(IS_FREE_SURFACE)==1.0) i1=1; //IS_IS_FREE_SURFACE //IS_BOUNDARY
	      if(GetGeometry()[2].FastGetSolutionStepValue(IS_FREE_SURFACE)==1.0) i2=1;
	      if(GetGeometry()[3].FastGetSolutionStepValue(IS_FREE_SURFACE)==1.0) i3=1;
	      ii0=0;
	      ii1=0;
	      ii2=0;
	      ii3=0;
	      sumai=0;
	      if(GetGeometry()[1].FastGetSolutionStepValue(IS_WATER)==1.0) ii1=1; //MATERIAL_VARIABLE
	      if(GetGeometry()[2].FastGetSolutionStepValue(IS_WATER)==1.0) ii2=1;
	      if(GetGeometry()[3].FastGetSolutionStepValue(IS_WATER)==1.0) ii3=1;
	
	      suma= i1 + i2 + i3;
	      sumai= ii1 + ii2 + ii3;	
	      
	      if(suma==3 or (suma==2 and sumai==1) or (suma==1 and sumai==2))// and sumai<3)
		{
		if(sumai==3) 
			{
							KRATOS_WATCH("aaaaaaaaaaaaaaaaaaaaaaaaaaa");				
			//KRATOS_THROW_ERROR(std::logic_error, "method not implemented", ""); 	
			}		 		  

	     // if(suma==3) 
		//{

		  PointsOfFSTriangle[0][0]=GetGeometry()[1].X();
		  PointsOfFSTriangle[0][1]=GetGeometry()[1].Y();
		  PointsOfFSTriangle[0][2]=GetGeometry()[1].Z();
		  
		  PointsOfFSTriangle[1][0]=GetGeometry()[2].X();
		  PointsOfFSTriangle[1][1]=GetGeometry()[2].Y();
		  PointsOfFSTriangle[1][2]=GetGeometry()[2].Z();
		  
		  PointsOfFSTriangle[2][0]=GetGeometry()[3].X();
		  PointsOfFSTriangle[2][1]=GetGeometry()[3].Y();
		  PointsOfFSTriangle[2][2]=GetGeometry()[3].Z();
		  
		  
		  area=CalculateTriangleArea3D(PointsOfFSTriangle[0], PointsOfFSTriangle[1], PointsOfFSTriangle[2]);
		  
		  T1=GetGeometry()[1].FastGetSolutionStepValue(YCH4);
		  T2=GetGeometry()[2].FastGetSolutionStepValue(YCH4);
		  T3=GetGeometry()[3].FastGetSolutionStepValue(YCH4);

		  int count=1;

		  if( GetGeometry()[1].FastGetSolutionStepValue(FACE_HEAT_FLUX)<=0.0) count=0;			  
		  if( GetGeometry()[2].FastGetSolutionStepValue(FACE_HEAT_FLUX)<=0.0) count=0;
		  if( GetGeometry()[3].FastGetSolutionStepValue(FACE_HEAT_FLUX)<=0.0) count=0;

 
		  rLeftHandSideMatrix(1,1) += count * 0.5 * ( 1.0*StefenBoltzmann*4.0*4.0*pow(T1,3))* 0.333333333333 * area;
		  rLeftHandSideMatrix(2,2) += count * 0.5 * ( 1.0*StefenBoltzmann*4.0*4.0*pow(T2,3))* 0.333333333333 * area;
		  rLeftHandSideMatrix(3,3) += count * 0.5 * ( 1.0*StefenBoltzmann*4.0*4.0*pow(T3,3))* 0.333333333333 * area;
		  
		  rRightHandSideVector[1] +=  GetGeometry()[1].FastGetSolutionStepValue(FACE_HEAT_FLUX) * 0.3333333333333*area - count * 0.5*1.0*StefenBoltzmann*4.0*(pow(T1,4))* 0.3333333333333*area ;
		  rRightHandSideVector[2] +=  GetGeometry()[2].FastGetSolutionStepValue(FACE_HEAT_FLUX) * 0.3333333333333*area - count * 0.5*1.0*StefenBoltzmann*4.0*(pow(T2,4))* 0.3333333333333*area ;
		  rRightHandSideVector[3] +=  GetGeometry()[3].FastGetSolutionStepValue(FACE_HEAT_FLUX) * 0.3333333333333*area - count * 0.5*1.0*StefenBoltzmann*4.0*(pow(T3,4))* 0.3333333333333*area ;


		  /* TPU CABLE	

		  
		  rRightHandSideVector[1] +=  GetGeometry()[1].FastGetSolutionStepValue(FACE_HEAT_FLUX) * 0.3333333333333*area -emissivity*StefenBoltzmann*(pow(T1,4) - aux)* 0.3333333333333*area -8.0* (T1 - 298)* 0.3333333333333*area ;
		  rRightHandSideVector[2] +=  GetGeometry()[2].FastGetSolutionStepValue(FACE_HEAT_FLUX) * 0.3333333333333*area - emissivity*StefenBoltzmann*(pow(T2,4) - aux)* 0.3333333333333*area -8.0* (T2 - 298)* 0.3333333333333*area ;
		  rRightHandSideVector[3] +=  GetGeometry()[3].FastGetSolutionStepValue(FACE_HEAT_FLUX) * 0.3333333333333*area - emissivity*StefenBoltzmann*(pow(T3,4) - aux)* 0.3333333333333*area -8.0* (T3 - 298)* 0.3333333333333*area ;

		  rLeftHandSideMatrix(1,1) += 1.0 * ( 8.0 + emissivity*StefenBoltzmann*4.0*pow(T1,3))* 0.333333333333 * area;
		  rLeftHandSideMatrix(2,2) += 1.0 * ( 8.0 + emissivity*StefenBoltzmann*4.0*pow(T2,3))* 0.333333333333 * area;
		  rLeftHandSideMatrix(3,3) += 1.0 * ( 8.0 + emissivity*StefenBoltzmann*4.0*pow(T3,3))* 0.333333333333 * area;
		 
		  GetGeometry()[1].FastGetSolutionStepValue(MIXTURE_FRACTION)+= 0.333333333333 * area;			  
		  GetGeometry()[2].FastGetSolutionStepValue(MIXTURE_FRACTION)+=0.333333333333 * area;
		  GetGeometry()[3].FastGetSolutionStepValue(MIXTURE_FRACTION)+= 0.333333333333 * area;

		 */

		}
	      i0=0;
	      i1=0;
	      i2=0;
	      i3=0;	
	      suma=0;
	      
	
	      if(GetGeometry()[0].FastGetSolutionStepValue(IS_FREE_SURFACE)==1.0) i0=1; 
	      if(GetGeometry()[2].FastGetSolutionStepValue(IS_FREE_SURFACE)==1.0) i2=1;
	      if(GetGeometry()[3].FastGetSolutionStepValue(IS_FREE_SURFACE)==1.0) i3=1;
	      suma= i0 + i2 + i3;
	      ii0=0;
	      ii1=0;
	      ii2=0;
	      ii3=0;	
	      sumai=0;
	      if(GetGeometry()[0].FastGetSolutionStepValue(IS_WATER)==1.0) ii0=1; 
	      if(GetGeometry()[2].FastGetSolutionStepValue(IS_WATER)==1.0) ii2=1;
	      if(GetGeometry()[3].FastGetSolutionStepValue(IS_WATER)==1.0) ii3=1;
	      sumai= ii0 + ii2 + ii3;	      		
	      if(suma==3 or (suma==2 and sumai==1) or (suma==1 and sumai==2))
		{	
		if(sumai==3) {
				KRATOS_WATCH("aaaaaaaaaaaaaaaaaaaaaaaaaaa");	
				//KRATOS_THROW_ERROR(std::logic_error, "method not implemented", ""); 
			}			
//	      if(suma==3) 
		//{

				//KRATOS_WATCH(sumai);	  
		  PointsOfFSTriangle[0][0]=GetGeometry()[0].X();
		  PointsOfFSTriangle[0][1]=GetGeometry()[0].Y();
		  PointsOfFSTriangle[0][2]=GetGeometry()[0].Z();
		  
		  PointsOfFSTriangle[1][0]=GetGeometry()[3].X();
		  PointsOfFSTriangle[1][1]=GetGeometry()[3].Y();
		  PointsOfFSTriangle[1][2]=GetGeometry()[3].Z();
		  
		  PointsOfFSTriangle[2][0]=GetGeometry()[2].X();
		  PointsOfFSTriangle[2][1]=GetGeometry()[2].Y();
		  PointsOfFSTriangle[2][2]=GetGeometry()[2].Z();
		  
		  area=CalculateTriangleArea3D(PointsOfFSTriangle[0], PointsOfFSTriangle[1], PointsOfFSTriangle[2]);
		  //KRATOS_WATCH(area); 
		  T0=GetGeometry()[0].FastGetSolutionStepValue(YCH4);
		  T2=GetGeometry()[2].FastGetSolutionStepValue(YCH4);
		  T3=GetGeometry()[3].FastGetSolutionStepValue(YCH4);
		  
		  int count=1;

		  if( GetGeometry()[0].FastGetSolutionStepValue(FACE_HEAT_FLUX)<=0.0) count=0;			  
		  if( GetGeometry()[3].FastGetSolutionStepValue(FACE_HEAT_FLUX)<=0.0) count=0;
		  if( GetGeometry()[2].FastGetSolutionStepValue(FACE_HEAT_FLUX)<=0.0) count=0;


 		  rLeftHandSideMatrix(0,0) += count * 0.5 * ( 1*StefenBoltzmann*4.0*4.0*pow(T0,3))* 0.333333333333 * area;
		  rLeftHandSideMatrix(2,2) += count * 0.5 * ( 1*StefenBoltzmann*4.0*4.0*pow(T2,3))* 0.333333333333 * area;
		  rLeftHandSideMatrix(3,3) += count * 0.5 * ( 1*StefenBoltzmann*4.0*4.0*pow(T3,3))* 0.333333333333 * area;
		  
		  rRightHandSideVector[0] +=  GetGeometry()[0].FastGetSolutionStepValue(FACE_HEAT_FLUX) * 0.3333333333333*area - count * 0.5*1*StefenBoltzmann*4.0*(pow(T0,4))* 0.3333333333333*area ;
		  rRightHandSideVector[2] +=  GetGeometry()[2].FastGetSolutionStepValue(FACE_HEAT_FLUX) * 0.3333333333333*area - count * 0.5*1*StefenBoltzmann*4.0*(pow(T2,4))* 0.3333333333333*area ;
		  rRightHandSideVector[3] +=  GetGeometry()[3].FastGetSolutionStepValue(FACE_HEAT_FLUX) * 0.3333333333333*area - count * 0.5*1*StefenBoltzmann*4.0*(pow(T3,4) )* 0.3333333333333*area ;
		   
		 
		/*TPU CABLE
		 rLeftHandSideMatrix(0,0) += 1.0*1.0 * ( 8.0 + emissivity*StefenBoltzmann*4.0*pow(T0,3))* 0.333333333333 * area;
		  rLeftHandSideMatrix(2,2) += 1.0*1.0 * ( 8.0 + emissivity*StefenBoltzmann*4.0*pow(T2,3))* 0.333333333333 * area;
		  rLeftHandSideMatrix(3,3) += 1.0*1.0 * ( 8.0 + emissivity*StefenBoltzmann*4.0*pow(T3,3))* 0.333333333333 * area;

		   rRightHandSideVector[0] +=  GetGeometry()[0].FastGetSolutionStepValue(FACE_HEAT_FLUX) * 0.3333333333333*area - emissivity*StefenBoltzmann*(pow(T0,4) - aux)* 0.3333333333333*area -8.0* (T0 - 298)* 0.3333333333333*area ;
		   rRightHandSideVector[2] +=  GetGeometry()[2].FastGetSolutionStepValue(FACE_HEAT_FLUX) * 0.3333333333333*area - emissivity*StefenBoltzmann*(pow(T2,4) - aux)* 0.3333333333333*area -8.0* (T2 - 298)* 0.3333333333333*area ;
		  rRightHandSideVector[3] +=  GetGeometry()[3].FastGetSolutionStepValue(FACE_HEAT_FLUX) * 0.3333333333333*area -  emissivity*StefenBoltzmann*(pow(T3,4) - aux)* 0.3333333333333*area -8.0* (T3 - 298)* 0.3333333333333*area ;

		  GetGeometry()[0].FastGetSolutionStepValue(MIXTURE_FRACTION)+= 0.333333333333 * area;			  
		  GetGeometry()[2].FastGetSolutionStepValue(MIXTURE_FRACTION)+=0.333333333333 * area;
		  GetGeometry()[3].FastGetSolutionStepValue(MIXTURE_FRACTION)+= 0.333333333333 * area;
		
		*/

		}
	      i0=0;
	      i1=0;
	      i2=0;
	      i3=0;	
	      suma=0;
	      if(GetGeometry()[0].FastGetSolutionStepValue(IS_FREE_SURFACE)==1.0) i0=1; 
	      if(GetGeometry()[1].FastGetSolutionStepValue(IS_FREE_SURFACE)==1.0) i1=1;
	      if(GetGeometry()[3].FastGetSolutionStepValue(IS_FREE_SURFACE)==1.0) i3=1;
	      suma= i0 + i1 + i3;

	      ii0=0;
	      ii1=0;
	      ii2=0;
	      ii3=0;	
	      sumai=0;
	      if(GetGeometry()[0].FastGetSolutionStepValue(IS_WATER)==1.0) ii0=1; 
	      if(GetGeometry()[1].FastGetSolutionStepValue(IS_WATER)==1.0) ii1=1;
	      if(GetGeometry()[3].FastGetSolutionStepValue(IS_WATER)==1.0) ii3=1;
	      sumai= ii0 + ii1 + ii3;	
	      
	      if(suma==3 or (suma==2 and sumai==1) or (suma==1 and sumai==2))//suma>=1 and sumai<3)
	      	{	 
	      	if(sumai==3) {
				KRATOS_WATCH("aaaaaaaaaaaaaaaaaaaaaaaaaaa");	
				//KRATOS_THROW_ERROR(std::logic_error, "method not implemented", ""); 
				}

	      //if(suma==3) 
		//{

		  PointsOfFSTriangle[0][0]=GetGeometry()[0].X();
		  PointsOfFSTriangle[0][1]=GetGeometry()[0].Y();
		  PointsOfFSTriangle[0][2]=GetGeometry()[0].Z();
		  
		  PointsOfFSTriangle[1][0]=GetGeometry()[1].X();
		  PointsOfFSTriangle[1][1]=GetGeometry()[1].Y();
		  PointsOfFSTriangle[1][2]=GetGeometry()[1].Z();
		  
		  PointsOfFSTriangle[2][0]=GetGeometry()[3].X();
		  PointsOfFSTriangle[2][1]=GetGeometry()[3].Y();
		  PointsOfFSTriangle[2][2]=GetGeometry()[3].Z();
		  
		  area=CalculateTriangleArea3D(PointsOfFSTriangle[0], PointsOfFSTriangle[1], PointsOfFSTriangle[2]);
		  
		  T0=GetGeometry()[0].FastGetSolutionStepValue(YCH4);
		  T1=GetGeometry()[1].FastGetSolutionStepValue(YCH4);
		  T3=GetGeometry()[3].FastGetSolutionStepValue(YCH4);
		  
		  int count=1;

		  if( GetGeometry()[0].FastGetSolutionStepValue(FACE_HEAT_FLUX)<=0.0) count=0;			  
		  if( GetGeometry()[1].FastGetSolutionStepValue(FACE_HEAT_FLUX)<=0.0) count=0;
		  if( GetGeometry()[3].FastGetSolutionStepValue(FACE_HEAT_FLUX)<=0.0) count=0;
   	          

		  rRightHandSideVector[0] +=  GetGeometry()[0].FastGetSolutionStepValue(FACE_HEAT_FLUX) * 0.3333333333333*area - count * 0.5*1*StefenBoltzmann*4.0*(pow(T0,4) )* 0.3333333333333*area ;
		  rRightHandSideVector[1] +=  GetGeometry()[1].FastGetSolutionStepValue(FACE_HEAT_FLUX) * 0.3333333333333*area - count * 0.5*1*StefenBoltzmann*4.0*(pow(T1,4) )* 0.3333333333333*area ;
		  rRightHandSideVector[3] +=  GetGeometry()[3].FastGetSolutionStepValue(FACE_HEAT_FLUX) * 0.3333333333333*area - count * 0.5*1*StefenBoltzmann*4.0*(pow(T3,4) )* 0.3333333333333*area ;


		  rLeftHandSideMatrix(0,0) += count * 0.5 * ( 1*StefenBoltzmann*4.0*4.0*pow(T0,3))* 0.333333333333 * area;
		  rLeftHandSideMatrix(1,1) += count * 0.5 * ( 1*StefenBoltzmann*4.0*4.0*pow(T1,3))* 0.333333333333 * area;
		  rLeftHandSideMatrix(3,3) += count * 0.5 * ( 1*StefenBoltzmann*4.0*4.0*pow(T3,3))* 0.333333333333 * area;
		
		  
		/*TPU CABLE
		rLeftHandSideMatrix(0,0) += 1.0 * ( 8.0 + emissivity*StefenBoltzmann*4.0*pow(T0,3))* 0.333333333333 * area;
		rLeftHandSideMatrix(1,1) += 1.0 * ( 8.0 + emissivity*StefenBoltzmann*4.0*pow(T1,3))* 0.333333333333 * area;
		rLeftHandSideMatrix(3,3) += 1.0 * ( 8.0 + emissivity*StefenBoltzmann*4.0*pow(T3,3))* 0.333333333333 * area;

		rRightHandSideVector[0] +=  GetGeometry()[0].FastGetSolutionStepValue(FACE_HEAT_FLUX) * 0.3333333333333*area -  emissivity*StefenBoltzmann*(pow(T0,4) - aux)* 0.3333333333333*area -8.0* (T0 - 298)* 0.3333333333333*area ;
		rRightHandSideVector[1] +=  GetGeometry()[1].FastGetSolutionStepValue(FACE_HEAT_FLUX) * 0.3333333333333*area - emissivity*StefenBoltzmann*(pow(T1,4) - aux)* 0.3333333333333*area -8.0* (T1 - 298)* 0.3333333333333*area ;
		rRightHandSideVector[3] +=  GetGeometry()[3].FastGetSolutionStepValue(FACE_HEAT_FLUX) * 0.3333333333333*area - emissivity*StefenBoltzmann*(pow(T3,4) - aux)* 0.3333333333333*area -8.0* (T3 - 298)* 0.3333333333333*area ;
			
	          GetGeometry()[0].FastGetSolutionStepValue(MIXTURE_FRACTION)+= 0.333333333333 * area;			  
		  GetGeometry()[1].FastGetSolutionStepValue(MIXTURE_FRACTION)+=0.333333333333 * area;
		  GetGeometry()[3].FastGetSolutionStepValue(MIXTURE_FRACTION)+= 0.333333333333 * area;

		*/	
		  
		}
	      i0=0;
	      i1=0;
	      i2=0;
	      i3=0;	
	      suma=0;
	      if(GetGeometry()[0].FastGetSolutionStepValue(IS_FREE_SURFACE)==1.0) i0=1; 
	      if(GetGeometry()[1].FastGetSolutionStepValue(IS_FREE_SURFACE)==1.0) i1=1;
	      if(GetGeometry()[2].FastGetSolutionStepValue(IS_FREE_SURFACE)==1.0) i2=1;
	      suma= i0 + i2 + i1;

	      ii0=0;
	      ii1=0;
	      ii2=0;
	      ii3=0;	
	      sumai=0;
	      if(GetGeometry()[0].FastGetSolutionStepValue(IS_WATER)==1.0) ii0=1; 
	      if(GetGeometry()[1].FastGetSolutionStepValue(IS_WATER)==1.0) ii1=1;
	      if(GetGeometry()[2].FastGetSolutionStepValue(IS_WATER)==1.0) ii2=1;
	      sumai= ii0 + ii2 + ii1;


	      
	      if(suma==3 or (suma==2 and sumai==1) or (suma==1 and sumai==2))//suma>=1 and sumai<3)
	      	{		
	      	  if(sumai==3) {
				KRATOS_WATCH("aaaaaaaaaaaaaaaaaaaaaaaaaaa");		
				//KRATOS_THROW_ERROR(std::logic_error, "method not implemented", ""); 
}
	 	//if(suma==3) 
		//{

		  //KRATOS_WATCH(sumai);
		  PointsOfFSTriangle[0][0]=GetGeometry()[0].X();
		  PointsOfFSTriangle[0][1]=GetGeometry()[0].Y();
		  PointsOfFSTriangle[0][2]=GetGeometry()[0].Z();
		  
		  PointsOfFSTriangle[1][0]=GetGeometry()[2].X();
		  PointsOfFSTriangle[1][1]=GetGeometry()[2].Y();
		  PointsOfFSTriangle[1][2]=GetGeometry()[2].Z();
		  
		  PointsOfFSTriangle[2][0]=GetGeometry()[1].X();
		  PointsOfFSTriangle[2][1]=GetGeometry()[1].Y();
		  PointsOfFSTriangle[2][2]=GetGeometry()[1].Z();
		  
		  area=CalculateTriangleArea3D(PointsOfFSTriangle[0], PointsOfFSTriangle[1], PointsOfFSTriangle[2]);
		  //KRATOS_WATCH(area); 
		  
		  double y=0.33333*(GetGeometry()[0].Y() + GetGeometry()[2].Y() + GetGeometry()[1].Y());
		  
		  T0=GetGeometry()[0].FastGetSolutionStepValue(YCH4);
		  T1=GetGeometry()[1].FastGetSolutionStepValue(YCH4);
		  T2=GetGeometry()[2].FastGetSolutionStepValue(YCH4);

                  int count=1;

		  if( GetGeometry()[0].FastGetSolutionStepValue(FACE_HEAT_FLUX)<=0.0) count=0;			  
		  if( GetGeometry()[1].FastGetSolutionStepValue(FACE_HEAT_FLUX)<=0.0) count=0;
		  if( GetGeometry()[2].FastGetSolutionStepValue(FACE_HEAT_FLUX)<=0.0) count=0;
		  
		  rRightHandSideVector[0] +=   GetGeometry()[0].FastGetSolutionStepValue(FACE_HEAT_FLUX) * 0.3333333333333*area - count * 0.5*1*StefenBoltzmann*4.0*(pow(T0,4) )* 0.3333333333333*area ; 
		  rRightHandSideVector[1] +=   GetGeometry()[1].FastGetSolutionStepValue(FACE_HEAT_FLUX) * 0.3333333333333*area - count * 0.5*1*StefenBoltzmann*4.0*(pow(T1,4) )* 0.3333333333333*area ;
		  rRightHandSideVector[2] +=   GetGeometry()[2].FastGetSolutionStepValue(FACE_HEAT_FLUX) * 0.3333333333333*area - count * 0.5*1*StefenBoltzmann*4.0*(pow(T2,4) )* 0.3333333333333*area ;

		  
		  rLeftHandSideMatrix(0,0) += count * 0.5 * ( 1*StefenBoltzmann*4.0*4.0*pow(T0,3))* 0.333333333333 * area;
		  rLeftHandSideMatrix(1,1) += count * 0.5 * ( 1*StefenBoltzmann*4.0*4.0*pow(T1,3))* 0.333333333333 * area;
		  rLeftHandSideMatrix(2,2) += count * 0.5 * ( 1*StefenBoltzmann*4.0*4.0*pow(T2,3))* 0.333333333333 * area;
		

		  /*TPU CABLE
		  rLeftHandSideMatrix(0,0)  += 1.0 * ( 8.0 + emissivity*StefenBoltzmann*4.0*pow(T0,3))* 0.333333333333 * area;
		  rLeftHandSideMatrix(1,1) += 1.0 * ( 8.0 + emissivity*StefenBoltzmann*4.0*pow(T1,3))* 0.333333333333 * area;
		  rLeftHandSideMatrix(2,2) += 1.0 * ( 8.0 + emissivity*StefenBoltzmann*4.0*pow(T2,3))* 0.333333333333 * area;
		
		   rRightHandSideVector[0] +=  GetGeometry()[0].FastGetSolutionStepValue(FACE_HEAT_FLUX) * 0.3333333333333*area - emissivity*StefenBoltzmann*(pow(T0,4) - aux)* 0.3333333333333*area -8.0* (T0 - 298)* 0.3333333333333*area ;
		  rRightHandSideVector[1] +=  GetGeometry()[1].FastGetSolutionStepValue(FACE_HEAT_FLUX) * 0.3333333333333*area - emissivity*StefenBoltzmann*(pow(T1,4) - aux)* 0.3333333333333*area -8.0* (T1 - 298)* 0.3333333333333*area ;
		  rRightHandSideVector[2] +=  GetGeometry()[2].FastGetSolutionStepValue(FACE_HEAT_FLUX) * 0.3333333333333*area - emissivity*StefenBoltzmann*(pow(T2,4) - aux)* 0.3333333333333*area -8.0* (T2 - 298)* 0.3333333333333*area ;
		
	          GetGeometry()[0].FastGetSolutionStepValue(MIXTURE_FRACTION)+= 0.333333333333 * area;			  
		  GetGeometry()[1].FastGetSolutionStepValue(MIXTURE_FRACTION)+=0.333333333333 * area;
		  GetGeometry()[2].FastGetSolutionStepValue(MIXTURE_FRACTION)+= 0.333333333333 * area;
		  */
		  
		}
//#endif
	    }

	  ///////////////////////////////////////////
	  /////////////////////////////////////////////
	  //double count=0.0;
	  double counti=0.0;
	  
//#if defined(RAD)

	   if(has_negative_node== false)
	     {	      
	       double b1=4.0;
	       double b2=3.0;
	       double R0 = GetGeometry()[0].FastGetSolutionStepValue(INCIDENT_RADIATION_FUNCTION);
	       double R1 = GetGeometry()[1].FastGetSolutionStepValue(INCIDENT_RADIATION_FUNCTION);
	       double R2 = GetGeometry()[2].FastGetSolutionStepValue(INCIDENT_RADIATION_FUNCTION);
	       double R3 = GetGeometry()[3].FastGetSolutionStepValue(INCIDENT_RADIATION_FUNCTION);
	       
	       T0=GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);
	       T1=GetGeometry()[1].FastGetSolutionStepValue(TEMPERATURE);
	       T2=GetGeometry()[2].FastGetSolutionStepValue(TEMPERATURE);
	       T3=GetGeometry()[3].FastGetSolutionStepValue(TEMPERATURE);
	      
	       msAuxVec[0]=pow(T0,4);
	       msAuxVec[1]=pow(T1,4);
	       msAuxVec[2]=pow(T2,4);
	       msAuxVec[3]=pow(T3,4);
	       
	       msAuxVec1[0]=R0;
	       msAuxVec1[1]=R1; 
	       msAuxVec1[2]=R2;
	       msAuxVec1[3]=R3;
	       
	       mInvJ0.resize(integration_points.size());
	       mDetJ0.resize(integration_points.size(),false);
	       
	       GeometryType::JacobiansType J0;
	       J0 = GetGeometry().Jacobian(J0, mThisIntegrationMethod);  
	       
	       for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
		 {
		   unsigned int number_of_nodes = GetGeometry().PointsNumber();
		   unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		   const Vector& N=row(Ncontainer,PointNumber);
		   
		   Tres = ZeroMatrix(4,4);	
		   Tres(0,0)=N[0]*N[0];
		   Tres(0,0)+=N[0]*N[1];
		   Tres(0,0)+=N[0]*N[2];
		   Tres(0,0)+=N[0]*N[3];
		   
		   
		   Tres(1,1)=N[1]*N[0];
		   Tres(1,1)+=N[1]*N[1];
		   Tres(1,1)+=N[1]*N[2];
		   Tres(1,1)+=N[1]*N[3];
		   
		   
		   Tres(2,2)=N[2]*N[0];
		   Tres(2,2)+=N[2]*N[1];
		   Tres(2,2)+=N[2]*N[2];
		   Tres(2,2)+=N[2]*N[3];
		   
		  
		   Tres(3,3)=N[3]*N[0];
		   Tres(3,3)+=N[3]*N[1];
		   Tres(3,3)+=N[3]*N[2];
		   Tres(3,3)+=N[3]*N[3];
		   
		   
		   MathUtils<double>::InvertMatrix(J0[PointNumber],mInvJ0[PointNumber],mDetJ0[PointNumber]);
		   double Weight = integration_points[PointNumber].Weight()* mDetJ0[PointNumber];
		   
		   double t_gauss=0.0;
		   
		   for (unsigned int i=0;i<number_of_nodes;i++) t_gauss += N[i]*GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE);	
		   double tgauss = pow(t_gauss,3);
		   //KRATOS_WATCH(absorptioncoefficient);
		   //KRATOS_WATCH(absorptioncoefficient2);
		   noalias(rLeftHandSideMatrix)+= 1.0 * Tres * Weight * tgauss * (4.0 * 4.0 * absorptioncoefficient * StefenBoltzmann) * so_heat_flux * 1.0;
	   
		   rRightHandSideVector += 1.0 * Weight * prod(Tres,msAuxVec ) *( -4.0  * absorptioncoefficient * StefenBoltzmann) * so_heat_flux * 1.0; //SIN ESTABILIX¡ZAR
		   rRightHandSideVector += 1.0 * Weight * absorptioncoefficient  *  prod(Tres,msAuxVec1 ) * so_heat_flux * 1.0; //SIN ESTABILIX¡ZAR
		    
		  }
	     }

//#endif   

	   //KRATOS_WATCH(Aux);
           //KRATOS_WATCH(rRightHandSideVector);
	   for (unsigned int iii = 0; iii < number_of_points; iii++)
	     ms_temp_vec_np[iii] = GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar) ;
	   noalias(rRightHandSideVector) -= prod(Aux, ms_temp_vec_np);    



	   	
 
	  }
	}
    else
      {
	
	const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
	array_1d<double, 3 > ms_vel_gauss;
	double conductivity = 0.0;
	double specific_heat = 0.0;
	double density = 1.0;
	double heat_flux = 0.0;
	const array_1d<double, 3 > & v= GetGeometry()[0].FastGetSolutionStepValue(rVelocityVar); 
	const array_1d<double, 3 > & w = GetGeometry()[0].FastGetSolutionStepValue(rMeshVelocityVar); 
	
	for (unsigned int j = 0; j < TDim; j++)
	  ms_vel_gauss[j] = v[j] - w[j];
	
	
	
	for (unsigned int i = 1; i < number_of_points; i++)
	  {
	
	const array_1d<double, 3 > & v = GetGeometry()[i].FastGetSolutionStepValue(rVelocityVar);	
	const array_1d<double, 3 > & w = GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar);
	for (unsigned int j = 0; j < TDim; j++)
	  ms_vel_gauss[j] += v[j] - w[j];
	
      }
    ms_vel_gauss *= lumping_factor;

    
    for (unsigned int i = 0; i < TDim; i++)
      {
	for (unsigned int j = 0; j < number_of_points; j++)
	{
		grad_g[i] += msDN_DX(j, i) * GetGeometry()[j].FastGetSolutionStepValue(rUnknownVar);
	}
     }
    conductivity=0.0131;//0.0454;//0.0131;//0.0454;//0.0131; //0.05;
    density=1.0;//1000.0;
    specific_heat=1310.0;//4540.0;//1310.0; //620.0;

    density=0.0;	
    for (unsigned int i = 1; i < number_of_points; i++)
    	      {
              density += GetGeometry()[i].FastGetSolutionStepValue(DENSITY);
    	      }
    density *= lumping_factor;

    double c1 = 4.00;
    double c2 = 2.00;
    double h = pow(6.00 * Area, 0.3333333);	
    //double h = sqrt(2.00 * Area);
    double norm_u = norm_2(ms_vel_gauss);
    double tau1;
    tau1 = (h * h) / (density * specific_heat * BDFcoeffs[0] * h * h + c1 * conductivity + c2 * density * specific_heat * (norm_u + 1e-6) * h);
    double res = density * specific_heat*(inner_prod(ms_vel_gauss,grad_g)) ;//+ 0.333333333333333 * (t0media+t1media+t2media)*(1/dt)*density*conductivity;
    double norm_grad=norm_2(grad_g);
    double k_aux=fabs(res) /(norm_grad + 0.000000000001);
    
    //calculating parameter tau
    //double c1 = 4.00;
    //double c2 = 2.00;
    //double h = pow(6.00 * Area, 0.3333333);
    //double norm_u = norm_2(ms_vel_gauss);
    //double tau1 = (h * h) / (density * specific_heat * BDFcoeffs[0] * h * h + c1 * conductivity + c2 * density * specific_heat * (norm_u + 1e-6) * h);

    //double pepe;
      
    noalias(First) = outer_prod(ms_vel_gauss, trans(ms_vel_gauss));
    First /= ((norm_u + 1e-6)*(norm_u + 1e-6));
    noalias(Second) = Identity - First;
    noalias(Third) = prod(Second, trans(msDN_DX));

    array_1d<double, 4 > a_dot_grad;
    noalias(a_dot_grad) = prod(msDN_DX, ms_vel_gauss);
    array_1d<double, 4 > a_dot_grad_and_mass;
    a_dot_grad_and_mass = msN;
    noalias(msAuxMat)= outer_prod(a_dot_grad, a_dot_grad_and_mass);
    
    //CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
    noalias(ms_u_DN) = prod(msDN_DX, ms_vel_gauss);

    noalias(rLeftHandSideMatrix) = (density * specific_heat) * outer_prod(msN, ms_u_DN);
    
    //CONVECTION STABILIZING CONTRIBUTION (Suu)
    
    noalias(rLeftHandSideMatrix) += density * specific_heat * tau1 * outer_prod(ms_u_DN, ms_u_DN);
    
    //VISCOUS CONTRIBUTION TO THE STIFFNESS MATRIX

    //KRATOS_WATCH(k_aux)	;
    //KRATOS_WATCH(h)	;
    noalias(rLeftHandSideMatrix) += ( (conductivity + 1.0 * k_aux * h) * prod(msDN_DX, trans(msDN_DX))); // + k_aux * h * prod(msDN_DX, Third));

    //rRightHandSideVector *= Area;
    // rLeftHandSideMatrix *= Area;
    
    //INERTIA CONTRIBUTION
    noalias(rLeftHandSideMatrix) += BDFcoeffs[0] * (density * specific_heat) * msMassFactors;
    
    // RHS = Fext
    //noalias(rRightHandSideVector) = 1.0 * (heat_flux * density) * msN;
    
    //RHS += Suy * proj[component]
    
    //rRightHandSideVector *= Area;
    rLeftHandSideMatrix *= Area;
    rRightHandSideVector *=0.0;  

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);
    mInvJ0.resize(integration_points.size());
    mDetJ0.resize(integration_points.size(),false);
    GeometryType::JacobiansType J0;
    J0 = GetGeometry().Jacobian(J0, mThisIntegrationMethod);  
    
    const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);
    
    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);
    
    for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
      {
        double IntegrationWeight = integration_points[PointNumber].Weight();
	MathUtils<double>::InvertMatrix(J0[PointNumber],mInvJ0[PointNumber],mDetJ0[PointNumber]);
	double Weight = integration_points[PointNumber].Weight()* mDetJ0[PointNumber];
	const Vector& N=row(Ncontainer,PointNumber);
	MathUtils<double>::InvertMatrix(J0[PointNumber],mInvJ0[PointNumber],mDetJ0[PointNumber]);
	noalias(msDN_DX) = prod(DN_De[PointNumber],mInvJ0[PointNumber]);
	
	Tres = ZeroMatrix(4,4);
	Tres(0,0)=N[0]*N[0];
	Tres(0,0)+=N[0]*N[1];
	Tres(0,0)+=N[0]*N[2];
	Tres(0,0)+=N[0]*N[3];
	
	Tres(1,1)=N[1]*N[0];
	Tres(1,1)+=N[1]*N[1];
	Tres(1,1)+=N[1]*N[2];
	Tres(1,1)+=N[1]*N[3];
	
	Tres(2,2)=N[2]*N[0];
	Tres(2,2)+=N[2]*N[1];
	Tres(2,2)+=N[2]*N[2];
	Tres(2,2)+=N[2]*N[3];
	
	Tres(3,3)=N[3]*N[0];
	Tres(3,3)+=N[3]*N[1];
	Tres(3,3)+=N[3]*N[2];
	Tres(3,3)+=N[3]*N[3];

	noalias(msDN_DX) = prod(DN_De[PointNumber],mInvJ0[PointNumber]);

	//KRATOS_WATCH(Tres);
	///

//	#if defined(IMPLICIT)





	if(rUnknownVar == YO) 
	  {    

	  double yo=0.0;
 	  double yf=0.0;

	  for (unsigned int j = 0; j < 4; j++) if (GetGeometry()[j].FastGetSolutionStepValue(YO)>0.0) yo += GetGeometry()[j].FastGetSolutionStepValue(YO);
          yo *=0.25;

	  for (unsigned int j = 0; j < 4; j++) if (GetGeometry()[j].FastGetSolutionStepValue(YF)>0.0) yf += GetGeometry()[j].FastGetSolutionStepValue(YF);
          yf *=0.25;
	

	    double temperature=  GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE)+GetGeometry()[1].FastGetSolutionStepValue(TEMPERATURE)+GetGeometry()[2].FastGetSolutionStepValue(TEMPERATURE)+GetGeometry()[3].FastGetSolutionStepValue(TEMPERATURE);
	    temperature *=0.25;

	    if(temperature>2600.0)temperature=2600.0;
 	
	    noalias(rLeftHandSideMatrix) += density* density* Weight * 1.92 * specific_heat *  yf * 6e9 * exp(-10700.0/(temperature)) * Tres;


	  }
	if(rUnknownVar == YF)
	  {


	  double yo=0.0;
 	  double yf=0.0;

	  for (unsigned int j = 0; j < 4; j++) if (GetGeometry()[j].FastGetSolutionStepValue(YO)>0.0) yo += GetGeometry()[j].FastGetSolutionStepValue(YO);
          yo *=0.25;

	  for (unsigned int j = 0; j < 4; j++) if (GetGeometry()[j].FastGetSolutionStepValue(YF)>0.0) yf += GetGeometry()[j].FastGetSolutionStepValue(YF);
          yf *=0.25;
	

	    double temperature=  GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE)+GetGeometry()[1].FastGetSolutionStepValue(TEMPERATURE)+GetGeometry()[2].FastGetSolutionStepValue(TEMPERATURE) + GetGeometry()[3].FastGetSolutionStepValue(TEMPERATURE);
	    temperature *=0.25;
        
	    if(temperature>2600.0)temperature=2600.0;

	
	    noalias(rLeftHandSideMatrix) +=  density * density * Weight * specific_heat *  yo * 6e9 * exp(-10700.0/(temperature)) * Tres;

	  }
//#endif


/*        double yo=0.0;
 	double yf=0.0;

	for (unsigned int j = 0; j < 4; j++) if (GetGeometry()[j].FastGetSolutionStepValue(YO)>0.0) yo += GetGeometry()[j].FastGetSolutionStepValue(YO);
        yo *=0.25;

	for (unsigned int j = 0; j < 4; j++) if (GetGeometry()[j].FastGetSolutionStepValue(YF)>0.0) yf += GetGeometry()[j].FastGetSolutionStepValue(YF);
        yf *=0.25;

        double temperature=  GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE)+GetGeometry()[1].FastGetSolutionStepValue(TEMPERATURE)+GetGeometry()[2].FastGetSolutionStepValue(TEMPERATURE)+GetGeometry()[3].FastGetSolutionStepValue(TEMPERATURE);
	temperature *=0.25;

	if(temperature>2600.0)temperature=2600.0;

	
	for (unsigned int iii = 0; iii < number_of_points; iii++)
      		ms_temp_vec_np[iii] = GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar);

	  
	if(rUnknownVar == YO) noalias(rRightHandSideVector) -=  1.92 * density * density * Weight * specific_heat *  yf  * 6e9 * exp(-10700.0/(temperature)) * prod(Tres,ms_temp_vec_np);

	if(rUnknownVar == YF) noalias(rRightHandSideVector) -=  density * density * Weight * specific_heat *  yo  * 6e9 * exp(-10700.0/(temperature)) * prod(Tres,ms_temp_vec_np);

*/
      }



	

	
	
    noalias(rLeftHandSideMatrix) +=  1.0/Dt * tau1 * msAuxMat * (density *specific_heat)* Area;  
    //adding the inertia terms
    // RHS += M*vhistory
    //calculating the historical velocity
    for (unsigned int iii = 0; iii < number_of_points; iii++)
      ms_temp_vec_np[iii] = GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar, 1);
    
    noalias(rRightHandSideVector) += 1.0/Dt * density * specific_heat* prod(msMassFactors, ms_temp_vec_np) * Area;
    
    noalias(rRightHandSideVector) += 1.0/Dt * tau1 * density * specific_heat * prod(msAuxMat, ms_temp_vec_np) * Area;
    
    //subtracting the dirichlet term
    // RHS -= LHS*temperatures
    for (unsigned int iii = 0; iii < number_of_points; iii++)
      ms_temp_vec_np[iii] = GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar);
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, ms_temp_vec_np);
  }
  KRATOS_CATCH("");
}   

//************************************************************************************
//************************************************************************************

void ConvDiff3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
}

//************************************************************************************
//************************************************************************************
// this subroutine calculates the nodal contributions for the explicit steps of the
// fractional step procedure

void ConvDiff3D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY
    int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

    BoundedMatrix<double, 4, 3 > msDN_DX;
    array_1d<double, 4 > msN;
    array_1d<double, 3 > ms_vel_gauss;
    array_1d<double, 4 > ms_temp_vec_np;
    array_1d<double, 4 > ms_u_DN;

    //getting data for the given geometry
    double Area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);
    ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
    const Variable<array_1d<double, 3 > >& rMeshVelocityVar = my_settings->GetMeshVelocityVariable();
    const Variable<array_1d<double, 3 > >& rVelocityVar = my_settings->GetVelocityVariable();
    const Variable<double>& rProjectionVariable = my_settings->GetProjectionVariable();

    if (FractionalStepNumber == 2) //calculation of temperature convective projection
    {
        const unsigned int number_of_points = GetGeometry().size();
        const double lumping_factor = 1.00 / double(number_of_points);
        unsigned int TDim = 3;

        //calculating viscosity
        ms_temp_vec_np[0] = GetGeometry()[0].FastGetSolutionStepValue(rUnknownVar);
        //const array_1d<double, 3 > & v = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
	const array_1d<double, 3 > & v = GetGeometry()[0].FastGetSolutionStepValue(rVelocityVar);
        const array_1d<double, 3 > & w = GetGeometry()[0].FastGetSolutionStepValue(rMeshVelocityVar);
        for (unsigned int j = 0; j < TDim; j++)
            ms_vel_gauss[j] = v[j] - w[j];

        for (unsigned int i = 1; i < number_of_points; i++)
        {
            ms_temp_vec_np[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);
	    const array_1d<double, 3 > & v = GetGeometry()[i].FastGetSolutionStepValue(rVelocityVar);
	    //const array_1d<double, 3 > & v = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3 > & w = GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar);
            for (unsigned int j = 0; j < TDim; j++)
                ms_vel_gauss[j] += v[j] - w[j];

        }
        ms_vel_gauss *= lumping_factor;

        //calculating convective auxiliary vector
        noalias(ms_u_DN) = prod(msDN_DX, ms_vel_gauss);
        double temp_conv = inner_prod(ms_u_DN, ms_temp_vec_np);
        temp_conv *= Area;

        for (unsigned int i = 0; i < number_of_points; i++)
        {
            GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += lumping_factor*Area;
            GetGeometry()[i].FastGetSolutionStepValue(rProjectionVariable) += lumping_factor*temp_conv;
            ;
        }
    }
    KRATOS_CATCH("");
}


//************************************************************************************
//************************************************************************************

void ConvDiff3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{

    ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if (rResult.size() != number_of_nodes)
        rResult.resize(number_of_nodes, false);

    for (unsigned int i = 0; i < number_of_nodes; i++)
        rResult[i] = GetGeometry()[i].GetDof(rUnknownVar).EquationId();
}

//************************************************************************************
//************************************************************************************

void ConvDiff3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

    if (ElementalDofList.size() != number_of_nodes)
        ElementalDofList.resize(number_of_nodes);

    for (unsigned int i = 0; i < number_of_nodes; i++)
        ElementalDofList[i] = GetGeometry()[i].pGetDof(rUnknownVar);

}

//************************************************************************************
//************************************************************************************

double ConvDiff3D::CalculateTriangleArea3D(	array_1d<double,3>& Point1, array_1d<double,3>& Point2, array_1d<double,3>& Point3	)
{
    //Heron's formula
    double a=Length(Point1, Point2);//sqrt((Point1[0]-Point2[0])*(Point1[0]-Point2[0]) + (Point1[1]-Point2[1])*(Point1[1]-Point2[1]) +(Point1[2]-Point2[2])*(Point1[2]-Point2[2]));
    double b=Length(Point1, Point3);//sqrt((Point3[0]-Point2[0])*(Point3[0]-Point2[0]) + (Point3[1]-Point2[1])*(Point3[1]-Point2[1]) +(Point3[2]-Point2[2])*(Point3[2]-Point2[2]));
    double c=Length(Point2, Point3);//sqrt((Point1[0]-Point3[0])*(Point1[0]-Point3[0]) + (Point1[1]-Point3[1])*(Point1[1]-Point3[1]) +(Point1[2]-Point3[2])*(Point1[2]-Point3[2]));
    double p=0.5*(a+b+c);
    return sqrt(p*(p-a)*(p-b)*(p-c));
}

//************************************************************************************
//************************************************************************************
  double ConvDiff3D::Length(array_1d<double,3>& Point1, array_1d<double,3>& Point2)
  {
    //KRATOS_WATCH("length calculation")
    return sqrt((Point1[0]-Point2[0])*(Point1[0]-Point2[0]) + (Point1[1]-Point2[1])*(Point1[1]-Point2[1]) +(Point1[2]-Point2[2])*(Point1[2]-Point2[2]));
  }
//************************************************************************************
//************************************************************************************

/*double ConvDiff3D::ComputeSmagorinskyViscosity(const BoundedMatrix<double, 4, 3 > & DN_DX, const double& h, const double& C, const double nu )
{
    BoundedMatrix<double, 3, 3 > dv_dx = ZeroMatrix(3, 3);

    const unsigned int nnodes = 4;

    ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<array_1d<double, 3 > >& rVelocityVar = my_settings->GetVelocityVariable();



    for (unsigned int k = 0; k < nnodes; ++k)
    {
        const array_1d< double, 3 > & rNodeVel = this->GetGeometry()[k].FastGetSolutionStepValue(rVelocityVar);

        for (unsigned int i = 0; i < 3; ++i)
        {
            for (unsigned int j = 0; j < i; ++j) // Off-diagonal
                dv_dx(i, j) += 0.5 * (DN_DX(k, j) * rNodeVel[i] + DN_DX(k, i) * rNodeVel[j]);
            dv_dx(i, i) += DN_DX(k, i) * rNodeVel[i]; // Diagonal
        }
    }

    // Norm[ Grad(u) ]
    double NormS(0.0);
    for (unsigned int i = 0; i < 3; ++i)
    {
        for (unsigned int j = 0; j < i; ++j)
            NormS += 2.0 * dv_dx(i, j) * dv_dx(i, j); // Using symmetry, lower half terms of the matrix are added twice
        NormS += dv_dx(i, i) * dv_dx(i, i); // Diagonal terms
    }

    NormS = sqrt(NormS);

    // Total Viscosity
    return 2.0 * C * C * h * h * NormS;
}*/

} // Namespace Kratos


