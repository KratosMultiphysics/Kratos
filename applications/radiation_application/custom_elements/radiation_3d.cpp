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
//   Last modified by:    $Author: julio.marti $
//   Date:                $Date:  $
//   Revision:            $Revision:  $
//
//
 

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/radiation_3d.h"
#include "radiation_application.h"
#include "includes/radiation_settings.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 
#include "utilities/enrichment_utilities.h"
#include "includes/kratos_flags.h"
#include "includes/deprecated_variables.h"
//#include "utilities/discont_utils.h"
#define P1
#define ENRICHMENT

namespace Kratos
{
  
  
  //************************************************************************************
  //************************************************************************************
  Rad3D::Rad3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
  {		
    //DO NOT ADD DOFS HERE!!!
  }
  
  //************************************************************************************
  //************************************************************************************
  Rad3D::Rad3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
  {
  }
  
  Element::Pointer Rad3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
  {
    return Element::Pointer(new Rad3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
  }
  
  Rad3D::~Rad3D()
  {
  }
  
  //************************************************************************************
  //************************************************************************************
  Rad3D::IntegrationMethod Rad3D::GetIntegrationMethod1()
  {
    return mThisIntegrationMethod;
  }
  //************************************************************************************
  //************************************************************************************
  void Rad3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY
      
      const unsigned int number_of_points = GetGeometry().size();
    //const double lumping_factor = 1.00/double(number_of_points);
    //unsigned int TDim = 3;
    //ddddddddddd
    KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
    boost::numeric::ublas::bounded_matrix<double,4,4> msMassFactors = 0.25*IdentityMatrix(4,4);
    boost::numeric::ublas::bounded_matrix<double,4,4> NN =  ZeroMatrix(4,4);
    boost::numeric::ublas::bounded_matrix<double,4,3> msDN_DX = ZeroMatrix(4,3);
    boost::numeric::ublas::bounded_matrix<double,4,3> msDN_DX_aux = ZeroMatrix(4,3);
    array_1d<double,4> msN = ZeroVector(4); //dimension = number of nodes
    array_1d<double,3> ms_vel_gauss = ZeroVector(3); //dimesion coincides with space dimension
    array_1d<double,4> ms_temp_vec_np = ZeroVector(4); //dimension = number of nodes
    array_1d<double,4> ms_u_DN = ZeroVector(4); //dimension = number of nodes
    array_1d<double,3> grad_g = ZeroVector(3); //dimesion coincides with space dimension
    boost::numeric::ublas::bounded_matrix<double,3,3> First = ZeroMatrix(3,3);
    boost::numeric::ublas::bounded_matrix<double,3,3> Second = ZeroMatrix(3,3);
    boost::numeric::ublas::bounded_matrix<double,3,4> Third = ZeroMatrix(3,4);
    boost::numeric::ublas::bounded_matrix<double,3,3> Identity = 1.0*IdentityMatrix(3,3);
    array_1d<double,3> aux_var= ZeroVector(3); //dimesion coincides with space dimension
    boost::numeric::ublas::bounded_matrix<double,4,1> msShapeFunc = ZeroMatrix(4,1);
    boost::numeric::ublas::bounded_matrix<double,1,4> msConvOp = ZeroMatrix(1,4);
    array_1d<double,4> msAuxVec; // = ZeroVector(4); //dimension = number of nodes
    boost::numeric::ublas::bounded_matrix<double,4,4> msAuxMat = ZeroMatrix(4,4);
    boost::numeric::ublas::bounded_matrix<double,4,4> Cuatro = ZeroMatrix(4,4);
    boost::numeric::ublas::bounded_matrix<double,4,4> Cuatros = ZeroMatrix(4,4);

    
    mThisIntegrationMethod= GeometryData::GI_GAUSS_2;
    // mThisIntegrationMethod= GeometryData::GI_GAUSS_3;
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);
    
    
    
    if(rLeftHandSideMatrix.size1() != number_of_points)
      rLeftHandSideMatrix.resize(number_of_points,number_of_points,false);
    
    if(rRightHandSideVector.size() != number_of_points)
      rRightHandSideVector.resize(number_of_points,false);
    
    //getting data for the given geometry
    double Area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);
    
    RadiationSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(RADIATION_SETTINGS);
    
    //Properties::Pointer properties = ThisModelPart.GetMesh().pGetProperties(1);
    
    //const Variable<double>&  rDensityVar = my_settings->GetDensityVariable();
    
    const Variable<double>& rUnknownVar= my_settings->GetUnknownVariable();
    
    //const Variable<array_1d<double,3> >& rMeshVelocityVar =my_settings->GetMeshVelocityVariable();
    
    const double StefenBoltzmann = 5.67e-8;
    double absorptioncoefficient = 10.0;
    
    const double T0 = GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);
    const double T1 = GetGeometry()[1].FastGetSolutionStepValue(TEMPERATURE);
    const double T2 = GetGeometry()[2].FastGetSolutionStepValue(TEMPERATURE);
    const double T3 = GetGeometry()[3].FastGetSolutionStepValue(TEMPERATURE);

    //double count=0.0;
    
    ///definitions
    bool has_negative_node=false;
    bool has_positive_node=false;
    
    boost::numeric::ublas::bounded_matrix<double,4, 3 > coords;
    array_1d<double,4> distances;
    array_1d<double,6>  volumes(6);
    boost::numeric::ublas::bounded_matrix<double,6, 4 > Ngauss;
    array_1d<double,6>  signs(6);
    std::vector< Matrix > gauss_gradients(6);
    //int number_interface_elements_1=0;
    
    unsigned int TNumNodes = GetGeometry().size();
    
    
    for(unsigned int iii = 0; iii<4; iii++)
      {
	distances[iii] = this->GetGeometry()[iii].FastGetSolutionStepValue(DISTANCE);  
      }
    
    ////////////////////////
    double unsigned_distance0=fabs(distances[0]);
    double unsigned_distance1=fabs(distances[1]);
    double unsigned_distance2=fabs(distances[2]);
    double unsigned_distance3=fabs(distances[3]);
    //we begin by finding the largest distance:
    double longest_distance=fabs(unsigned_distance0);
    if (unsigned_distance1>longest_distance) longest_distance=unsigned_distance1;
    if (unsigned_distance2>longest_distance) longest_distance=unsigned_distance2;
    if (unsigned_distance3>longest_distance) longest_distance=unsigned_distance3;
    //Now we set a maximum relative distance
    
    double tolerable_distance =longest_distance*1.0;   //0.000001
    
    if (unsigned_distance0<tolerable_distance)
      {
	if (unsigned_distance0>1.0e-10)
	  {
	    distances[0]=tolerable_distance*(distances[0]/fabs(distances[0]));
	  }
	else 
	  {	
	    distances[0]=tolerable_distance;
	  }
      }
    if (unsigned_distance1<tolerable_distance)
      {
	if (unsigned_distance1>1.0e-10)
	  {
	    distances[1]=tolerable_distance*(distances[1]/fabs(distances[1]));
	  }				
	else
	  {
	    distances[1]=tolerable_distance;
	    //////KRATOS_WATCH(rDistances[1]);	
	  }
	
      }
    if (unsigned_distance2<tolerable_distance)
      {
	if (unsigned_distance2>1.0e-10)
	  {
	    distances[2]=tolerable_distance*(distances[2]/fabs(distances[2]));
	  }
	else{
	  distances[2]=tolerable_distance;
	}
      }
    if (unsigned_distance3<tolerable_distance)
      {
	if (unsigned_distance3>1.0e-10)
	  {
	    distances[3]=tolerable_distance*(distances[3]/fabs(distances[3]));
	  }
	else
	  {
	    distances[3]=tolerable_distance;
	  }
      }
    
    
    for(unsigned int iii = 0; iii<4; iii++)
      {
	//distances[iii] = this->GetGeometry()[iii].FastGetSolutionStepValue(DISTANCE);
	if (distances[iii]<0.0)
	  has_negative_node=true;
	else
	  has_positive_node=true;		
      }
    
    bool split_element=false;
    
    if (has_positive_node && has_negative_node) split_element=true;
    
    split_element=false;
    
    if (split_element==true){
      array_1d<double,6> absorptioncoefficient;
      absorptioncoefficient *=0.0;
      for (unsigned int i = 0; i < TNumNodes; i++)
	{
	  const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
	  for (unsigned int j = 0; j < 3; j++) coords(i, j) = xyz[j];
	}
	//int number_interface_elements_1=0;
      std::vector< Matrix > gauss_gradients(6);
      std::vector< Matrix > gauss_gradients_aux(8);
      boost::numeric::ublas::bounded_matrix<double,2 , 3 >  coord_interface_nodes_1; //used to p	
      boost::numeric::ublas::bounded_matrix<double,6, 2> Nenriched;
      array_1d<double,6> area_interface_1; 
      array_1d<double,6> area_inter;
      array_1d<double,6>  N_Star(6);
      array_1d<double,6>  wei(6);
	//bool switch_off_e=false;
      std::vector< Matrix > edges_taux(8);
      std::vector< Matrix > nodes_aux(8);
      std::vector< Matrix > tot_edges(8);
      std::vector< Matrix > rGradientaux(6);
      int total_nodes=0;
      std::vector< Matrix > coord_interface(4);
      boost::numeric::ublas::bounded_matrix<double,6, 8 > Ngaussnew;
      std::vector< Matrix > mass_matrix(6);
      std::vector< Matrix > MASS_MATRIX2(6);
      std::vector< Matrix > Cuatros(6);
      boost::numeric::ublas::bounded_matrix<double, 4, 4 > Laplacian_matrix = ZeroMatrix(4,4);
      
      for (unsigned int i = 0; i < 6; i++)
      	{
	  /*mass_matrix[i].resize(4, 4, false);  
	    mass_matrix[i] *=0.0;*/
	  
	  mass_matrix[i].resize(8, 8, false);  
	  mass_matrix[i] *=0.0;
	  
	  MASS_MATRIX2[i].resize(4, 4, false); // resize(4, 4, false);
	  MASS_MATRIX2[i] *=0.0;
	  Cuatros[i].resize(4, 4, false);  
	  Cuatros[i] *=0.0;
      	}  
      
      
      
      for (unsigned int i = 0; i < 6; i++)
	{
	  gauss_gradients[i].resize(2, 3, false);  //2 values of the 2 shape functions, and derivates in (xyz) direction).
	  //	  gauss_gradients1[i].resize(4, 5, false);
	  gauss_gradients_aux[i].resize(8, 3, false);  //2 values of the 2 shape functions, and derivates in (xyz) direction).
	  
	}
      for (unsigned int i = 0; i < 8; i++)
	{
	  edges_taux[i].resize(5, 3, false);  //2 values of the 2 shape functions, and derivates in (xyz) direction).
	  nodes_aux[i].resize(1, 3, false);
	  tot_edges[i].resize(5, 3, false);
	}
      for (unsigned int i = 0; i < 6; i++)
	{
	  rGradientaux[i].resize(8,3,false);
	}
      for (unsigned int i = 0; i < 4; i++)
	{
	  coord_interface[i].resize(1,3,false);
	}
      unsigned int ndivisions = 2;//EnrichmentUtilities::CalculateEnrichedShapeFuncions_Simplified(this->GetGeometry(),coords, msDN_DX, distances, volumes, Ngauss, signs, gauss_gradients,gauss_gradients_aux,Nenriched, number_interface_elements_1, coord_interface_nodes_1, area_interface_1,area_inter, N_Star,switch_off_e,edges_taux,nodes_aux, tot_edges,rGradientaux, total_nodes, coord_interface,Ngaussnew,mass_matrix, MASS_MATRIX2, wei);//,mass_matrixone); 
      
      double coeffabsorptioncoefficient=0.0;
      double voltotal=0.0;
      //double alphav=0.0;
      double conductivity=0.0; //=1.0/(3.0*absorptioncoefficient);
      for (unsigned int i = 0; i!=ndivisions; i++) //the 3 partitions of the triangle
	{
	  if (signs(i)>0)
	    {
	      absorptioncoefficient[i] = 0.01;//0.01;  
	      coeffabsorptioncoefficient += absorptioncoefficient[i] * volumes[i];
	      voltotal +=volumes[i];
	      //alphav += absorptioncoefficient[i] * volumes[i];
		   conductivity += 1.0/(3.0 * absorptioncoefficient[i]);//*volumes[i];
	    }
	  else
	    {
	      absorptioncoefficient[i] = 3000.0;  		 
	      coeffabsorptioncoefficient += absorptioncoefficient[i] * volumes[i];
	      voltotal +=volumes[i];
	      //alphav += absorptioncoefficient[i] * volumes[i];
		   conductivity += 1.0/(3.0 * absorptioncoefficient[i]);//*volumes[i];
	    }
	}	
      coeffabsorptioncoefficient /=   voltotal;  
      //array_1d<double,3> w = GetGeometry()[0].FastGetSolutionStepValue(rMeshVelocityVar);
      //ms_vel_gauss=w;
	
      double c1 = 4.00;
      double c2 = 2.00;
      double h =  pow(6.00*Area,0.3333333);
      //double norm_u =norm_2(ms_vel_gauss);
      
      //	double tau1=1.0/(( c2 * norm_u /h )+ c1 * absorptioncoefficient);
      
      double tau1=1.0/(( 0.0/*c2 * norm_u /h*/ )+ c1 * coeffabsorptioncoefficient);
      //double pi=3.1416;
      //tau1*=0.0;
      //double h_aux;
      //double Volume_aux;
      //the coefficients INCLUDE the time step
      //const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
      
      //CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
      noalias(ms_u_DN) = prod(msDN_DX , ms_vel_gauss);
      
      noalias(rLeftHandSideMatrix) = 1.0 * outer_prod(msN,ms_u_DN) * Area;
      //double tau_total=0.0;
      noalias(rLeftHandSideMatrix) +=1.0 * tau1 * outer_prod(ms_u_DN,ms_u_DN) * Area;
      
      
      msMassFactors(0,0) = 1.00/4.00; msMassFactors(0,1) = 0.00;		msMassFactors(0,2) = 0.00;       msMassFactors(0,3) = 0.00;
      msMassFactors(1,0) = 0.00;		msMassFactors(1,1) = 1.00/4.00; msMassFactors(1,2) = 0.00;       msMassFactors(1,3) = 0.00;
      msMassFactors(2,0) = 0.00;		msMassFactors(2,1) = 0.00;		msMassFactors(2,2) = 1.00/4.00;  msMassFactors(2,3) = 0.00;
      msMassFactors(3,0) = 0.00;		msMassFactors(3,1) = 0.00;		msMassFactors(3,2) = 0.00;       msMassFactors(3,3) = 1.00/4.00;
      
      msAuxVec[0]=pow(T0,4); 
      msAuxVec[1]=pow(T1,4);
      msAuxVec[2]=pow(T2,4);
      msAuxVec[3]=pow(T3,4);
      
      //int nodes_number = 4;
      //int dof = 1;
      array_1d<double,4> N;
      noalias(rRightHandSideVector)=ZeroVector(4);
      
      
      //double xxx=0.0;
      for (unsigned int i = 0; i!=ndivisions; i++) //the 3 partitions of the triangle
	{
	  for (unsigned int kkk = 0; kkk<4; kkk++) //gauss points   kkk<4
	    {
	      N[0]=MASS_MATRIX2[i](kkk,0);
	      N[1]=MASS_MATRIX2[i](kkk,1);
	      N[2]=MASS_MATRIX2[i](kkk,2);
	      N[3]=MASS_MATRIX2[i](kkk,3);
	      
	      Cuatros[i](0,0)+=N[0]*N[0]*wei(i);
	      //Cuatros[i](0,1)+=N[0]*N[1]*wei(i);
	      //Cuatros[i](0,2)+=N[0]*N[2]*wei(i);
	      //Cuatros[i](0,3)+=N[0]*N[3]*wei(i);
	      
	      //Cuatros[i](0,0)+=N[0]*N[0]*wei(i);
	      Cuatros[i](0,0)+=N[0]*N[1]*wei(i);
	      Cuatros[i](0,0)+=N[0]*N[2]*wei(i);
	      Cuatros[i](0,0)+=N[0]*N[3]*wei(i);
	      
	      
	      
	      //		    Cuatros[i](1,0)+=N[1]*N[0]*wei(i);
	      Cuatros[i](1,1)+=N[1]*N[1]*wei(i);
	      //		    Cuatros[i](1,2)+=N[1]*N[2]*wei(i);
	      //		    Cuatros[i](1,3)+=N[1]*N[3]*wei(i);
	      
	      Cuatros[i](1,1)+=N[1]*N[0]*wei(i);
	      //		    Cuatros[i](1,1)+=N[1]*N[1]*wei(i);
	      Cuatros[i](1,1)+=N[1]*N[2]*wei(i);
	      Cuatros[i](1,1)+=N[1]*N[3]*wei(i);

	      
	      //		    Cuatros[i](2,0)+=N[2]*N[0]*wei(i);
	      //		    Cuatros[i](2,1)+=N[2]*N[1]*wei(i);
	      Cuatros[i](2,2)+=N[2]*N[2]*wei(i);
	      //		    Cuatros[i](2,3)+=N[2]*N[3]*wei(i);
	      
	      Cuatros[i](2,2)+=N[2]*N[0]*wei(i);
	      Cuatros[i](2,2)+=N[2]*N[1]*wei(i);
	      //Cuatros[i](2,2)+=N[2]*N[2]*wei(i);
	      Cuatros[i](2,2)+=N[2]*N[3]*wei(i);
	      
	      
	      //		    Cuatros[i](3,0)+=N[3]*N[0]*wei(i);
	      //		    Cuatros[i](3,1)+=N[3]*N[1]*wei(i);
	      //		    Cuatros[i](3,2)+=N[3]*N[2]*wei(i);
	      Cuatros[i](3,3)+=N[3]*N[3]*wei(i);
	      
	      Cuatros[i](3,3)+=N[3]*N[0]*wei(i);
	      Cuatros[i](3,3)+=N[3]*N[1]*wei(i);
	      Cuatros[i](3,3)+=N[3]*N[2]*wei(i);
	      //Cuatros[i](3,3)+=N[3]*N[3]*wei(i);
	      
	      
	      array_1d<double, 4 > a_dot_grad;
	      noalias(a_dot_grad) = prod(msDN_DX, ms_vel_gauss);
	      array_1d<double, 4 > a_dot_grad_and_mass;
	      a_dot_grad_and_mass = N;
	      noalias(msAuxMat)= outer_prod(a_dot_grad, a_dot_grad_and_mass);
	      
	      //Volume_aux=volumes[i];
	      msDN_DX_aux=rGradientaux[i];
	      //h_aux = CalculateH(msDN_DX_aux,Volume_aux);
	      
	      //noalias(rLeftHandSideMatrix) += 1.0 * tau1 * absorptioncoefficient[i] * msAuxMat * wei(i);
	      
	      //noalias(rRightHandSideVector) += 1.0 * tau1 * wei(i) *  absorptioncoefficient[i] * StefenBoltzmann * (1.0/pi)* prod(msAuxMat,msAuxVec );
	      
	    }
	}
      
      

      
      
      
      
#if defined(P1)
      rLeftHandSideMatrix *=0.0;
      rRightHandSideVector *=0.0;
      
	
      noalias(rLeftHandSideMatrix) = prod(msDN_DX,trans(msDN_DX)) * conductivity;
      rLeftHandSideMatrix *=0.0;
      rRightHandSideVector *=0.0;
      for (unsigned int i = 0; i!=ndivisions; i++) //the 3 partitions of the triangle
	{
	  noalias(rLeftHandSideMatrix) +=  1.0 * absorptioncoefficient[i] *   Cuatros[i];//MASS_MATRIX2[i];// * volumes[i];
	  noalias(rRightHandSideVector) +=  4.0 * StefenBoltzmann * absorptioncoefficient[i] *  prod(Cuatros[i],msAuxVec);// prod(MASS_MATRIX2[i
	}
      
#if defined(ENRICHMENT)
      //rLeftHandSideMatrix *=0.0;
      //rRightHandSideVector *=0.0;
      
      Matrix Laplacian_matrix_mod;
      Laplacian_matrix_mod.resize(total_nodes, total_nodes);
      Laplacian_matrix_mod = ZeroMatrix(total_nodes,total_nodes);
      
      Matrix Kaux1;
      Kaux1.resize(total_nodes,total_nodes);
      Kaux1 = ZeroMatrix(total_nodes,total_nodes);  
      
      Matrix Kaux;
      Kaux.resize(total_nodes - 4, 4);
      Kaux = ZeroMatrix(total_nodes - 4, 4);  
      
      Matrix Lenrichaux;
      Lenrichaux.resize(total_nodes - 4,total_nodes -4);
      Lenrichaux = ZeroMatrix(total_nodes -4,total_nodes-4);
      
      Matrix Lenrichaux1;
      
      Lenrichaux1.resize(total_nodes,total_nodes);
      Lenrichaux1 = ZeroMatrix(total_nodes,total_nodes); 
      
      Kaux1 *=0.0;
      Lenrichaux1 *=0.0;
      
      qi( ndivisions, edges_taux, nodes_aux, rGradientaux, absorptioncoefficient, Kaux1, Lenrichaux1);
      
      
      int max;
      if(ndivisions==2) 
	{
	  max=1; 
	}
      else if(ndivisions==3) 
	{
	  max=2;
	}
      else if(ndivisions==4) 
	{
	  max=3;
	}
      else {max=4;
      }
      max=total_nodes-4;
      
      for(unsigned int m = 0; m < max; m++)
	for (unsigned int s = 0; s < max; s++) 
	  Lenrichaux(m,s) = Lenrichaux1(4+m,s+4) ;
      
      for(unsigned int m = 0; m < max; m++)
	for (unsigned int s = 0; s < 4; s++) 
	  Kaux(m,s) = Kaux1(4+m,s);
      
      for (unsigned int i=0; i<ndivisions ; i++)
	{
	  Laplacian_matrix_mod +=prod (rGradientaux[i], trans(rGradientaux[i])) * 1.0/(3.0 * absorptioncoefficient[i]) * volumes[i];
	}
      
      Matrix mass_matrix_mod;
      mass_matrix_mod.resize(total_nodes, total_nodes);
      mass_matrix_mod = ZeroMatrix(total_nodes,total_nodes);
      
      
      Matrix Kbar;
      Matrix Abar;
      Matrix A;
      
      Kbar.resize(4,total_nodes -4);
      Kbar = ZeroMatrix(4,total_nodes-4);  
      Abar.resize(total_nodes-4,total_nodes-4);
      Abar = ZeroMatrix(total_nodes-4,total_nodes-4);
      A.resize(total_nodes-4,4);
      A = ZeroMatrix(total_nodes-4,4); 
      
      for(unsigned int m = 0; m < max; m++)  //3
	for (unsigned int j = 0; j < 4; j++) 
	  A(m,j)=Laplacian_matrix_mod(4 + m , j)  - Kaux(m,j);
      
      for(unsigned int m = 0; m < 4; m++)
	for (unsigned int j = 0; j < max; j++) 
	  Kbar(m,j)=Laplacian_matrix_mod(m,j + 4);
      
      for(unsigned int m = 0; m < max; m++)
	for (unsigned int j = 0; j < max; j++) 
	  Abar(m,j) = Laplacian_matrix_mod(4+m,4+j) - Lenrichaux(m,j);
      
      Matrix inverse_enrichments;
      Matrix temp_matrix;
      
      inverse_enrichments.resize(total_nodes-4,total_nodes-4);
      inverse_enrichments = ZeroMatrix(total_nodes-4,total_nodes-4);
      temp_matrix.resize(4,total_nodes-4);
      temp_matrix = ZeroMatrix(4,total_nodes-4);
      
      //int ninter=0;
      //////////////////////////////////
      double heat=1000.0;
      heat *=0.0;
      //Heat_Source(rRightHandSideVector, ndivisions, volumes, conductivities, Ngaussnew, heat);
      //double coefficientaux=0.0; 
      
      if(ndivisions==3 or ndivisions==2)  KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
      
      //////////////////////////////////
      //if(ndivisions==4) ninter=1;
      //if(ndivisions==6) ninter=2;
      //if(ndivisions==3) ninter=0;
      //if(ndivisions==2) ninter=0;
      //int number_of_nodes=0;
      //if(ndivisions==2) {number_of_nodes=1; /*KRATOS_ERROR(std::logic_error, "method not implemented", "");*/}
      //if(ndivisions==3) {number_of_nodes=2; /*KRATOS_ERROR(std::logic_error, "method not implemented", "");*/}
      //if(ndivisions==4) number_of_nodes=3;
      //if(ndivisions==6) number_of_nodes=4;
      
      Matrix Abar_aux;
      Abar_aux.resize(total_nodes-4,total_nodes-4);
      
      Abar_aux = ZeroMatrix(total_nodes-4,total_nodes-4);
      Abar_aux=Abar; ///NO TIENE EN CUENTA LA LINEARIZACIONNNNNNNNNNNNNNNNNNNNNNNNNNNN
		
      //bool prueba=false;
      //prueba=this->InvertMatrix(Abar,inverse_enrichments);
      
      temp_matrix= prod(Kbar, inverse_enrichments);
      
      for(unsigned int m = 0; m < 4; m++)
	for (unsigned int j = 0; j < 4; j++) 
	  Laplacian_matrix(m,j)=Laplacian_matrix_mod(m,j) + 0.0 * mass_matrix_mod(m,j);
      
      noalias(rLeftHandSideMatrix) += Laplacian_matrix;
      
      noalias(rLeftHandSideMatrix) -=  prod(temp_matrix,A);
      
      
      
#endif(ENRICHMENT)
#endif(P1)
      
      ///nuevooooooooooooooooooooooooooooooooooooo
      
      
      
    }
    
    else
      {
	
        //KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
	if (has_negative_node== true) //polymer
	  {
	    absorptioncoefficient=0.0; //10.0
	    KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
	  }
        else //air
	  {
	    absorptioncoefficient=100.0;//0.01;
	  }
	//absorptioncoefficient=0.01;//0.01;
	//array_1d<double,3> w = GetGeometry()[0].FastGetSolutionStepValue(rMeshVelocityVar);
	//ms_vel_gauss=w;
	double c1 = 4.00;
	double c2 = 2.00;
	double h =  pow(6.00*Area,0.3333333);
	//double norm_u =norm_2(ms_vel_gauss);
	
	double tau1=1.0/(( 0.0 )+ c1 * absorptioncoefficient);
	//double pi=3.1416;
	
	//the coefficients INCLUDE the time step
	//const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
	
	//CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
	noalias(ms_u_DN) = prod(msDN_DX , ms_vel_gauss);
	
	noalias(rLeftHandSideMatrix) = 1.0 * outer_prod(msN,ms_u_DN) * Area;
	//tau1*=0.0;
	noalias(rLeftHandSideMatrix) += 1.0 * tau1 * outer_prod(ms_u_DN,ms_u_DN) * Area;
	
	//noalias(rLeftHandSideMatrix) =ZeroMatrix(4,4);;

	//getting the BDF2 coefficients (not fixed to allow variable time step)
	//the coefficients INCLUDE the time step
	msMassFactors(0,0) = 1.00/4.00; msMassFactors(0,1) = 0.00;		msMassFactors(0,2) = 0.00;       msMassFactors(0,3) = 0.00;
	msMassFactors(1,0) = 0.00;		msMassFactors(1,1) = 1.00/4.00; msMassFactors(1,2) = 0.00;       msMassFactors(1,3) = 0.00;
	msMassFactors(2,0) = 0.00;		msMassFactors(2,1) = 0.00;		msMassFactors(2,2) = 1.00/4.00;  msMassFactors(2,3) = 0.00;
	msMassFactors(3,0) = 0.00;		msMassFactors(3,1) = 0.00;		msMassFactors(3,2) = 0.00;       msMassFactors(3,3) = 1.00/4.00;
    
	double 	conductivity=1.0/(3.0*absorptioncoefficient);
	
	msAuxVec[0]=pow(T0,4); 
	msAuxVec[1]=pow(T1,4);
	msAuxVec[2]=pow(T2,4);
	msAuxVec[3]=pow(T3,4);
	
	
	//int nodes_number = 4;
	//int dof = 1;
	
	msShapeFunc = ZeroMatrix(4,1);
	msAuxMat = ZeroMatrix(4,4);
	mInvJ0.resize(integration_points.size());
	mDetJ0.resize(integration_points.size(),false);
	
	GeometryType::JacobiansType J0;
	J0 = GetGeometry().Jacobian(J0, mThisIntegrationMethod);
	const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);
	
	const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);

	  noalias(rRightHandSideVector)=ZeroVector(4);
	  
	noalias(rLeftHandSideMatrix) = (conductivity) * prod(msDN_DX,trans(msDN_DX)) * Area; 
	noalias(rLeftHandSideMatrix) += absorptioncoefficient * msMassFactors * Area; //LO COMENTO AHORA

        double T0,T1,T2,T3;
	T0=GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);
  	T1=GetGeometry()[1].FastGetSolutionStepValue(TEMPERATURE);
	T2=GetGeometry()[2].FastGetSolutionStepValue(TEMPERATURE);
	T3=GetGeometry()[3].FastGetSolutionStepValue(TEMPERATURE);
	
	std::vector<array_1d<double,3> > PointsOfFSTriangle;
	PointsOfFSTriangle.reserve(3);
	int i0,i1,i2,i3;
	int suma;
	i0=0;
	i1=0;
	i2=0;
	i3=0;	
	suma=0;
	double area=0.0;	  
	if(GetGeometry()[1].FastGetSolutionStepValue(IS_BOUNDARY)==1.0) i1=1; 
	if(GetGeometry()[2].FastGetSolutionStepValue(IS_BOUNDARY)==1.0) i2=1;
	if(GetGeometry()[3].FastGetSolutionStepValue(IS_BOUNDARY)==1.0) i3=1;
	suma= i1 + i2 + i3;
	double constant;//=0;
	constant=1.0/(4-2*1.0);
	//constant=0.0;
	if(suma==3)
	{
	  PointsOfFSTriangle[0][0]=GetGeometry()[1].X();
	  PointsOfFSTriangle[0][1]=GetGeometry()[1].Y();
	  PointsOfFSTriangle[0][2]=GetGeometry()[1].Z();
	  
	  PointsOfFSTriangle[1][0]=GetGeometry()[2].X();
	  PointsOfFSTriangle[1][1]=GetGeometry()[2].Y();
	  PointsOfFSTriangle[1][2]=GetGeometry()[2].Z();
	  
	  PointsOfFSTriangle[2][0]=GetGeometry()[3].X();
	  PointsOfFSTriangle[2][1]=GetGeometry()[3].Y();
	  PointsOfFSTriangle[2][2]=GetGeometry()[3].Z();
	  
	  
	  area=0.5 * CalculateTriangleArea3D(PointsOfFSTriangle[0], PointsOfFSTriangle[1], PointsOfFSTriangle[2]);

	  T1=GetGeometry()[1].FastGetSolutionStepValue(TEMPERATURE);
	  T2=GetGeometry()[2].FastGetSolutionStepValue(TEMPERATURE);
	  T3=GetGeometry()[3].FastGetSolutionStepValue(TEMPERATURE);
	  
	  
	  rRightHandSideVector[1] +=  constant * 4.0 * StefenBoltzmann*(pow(T1,4) ) * 0.333333333333* area;
	  rRightHandSideVector[2] +=  constant * 4.0 * StefenBoltzmann*(pow(T2,4) ) * 0.333333333333* area; 
	  rRightHandSideVector[3] +=  constant * 4.0 * StefenBoltzmann*(pow(T3,4) ) * 0.333333333333*  area;
	  
	  
	  rLeftHandSideMatrix(1,1) +=  constant * 0.333333333333 * area;
	  rLeftHandSideMatrix(2,2) +=  constant * 0.333333333333 * area;
	  rLeftHandSideMatrix(3,3) +=  constant * 0.333333333333 * area;
	  
	  
	}
	i0=0;
	i1=0;
	i2=0;
	i3=0;	
	suma=0;
	if(GetGeometry()[0].FastGetSolutionStepValue(IS_BOUNDARY)==1.0) i0=1; 
	if(GetGeometry()[2].FastGetSolutionStepValue(IS_BOUNDARY)==1.0) i2=1;
	if(GetGeometry()[3].FastGetSolutionStepValue(IS_BOUNDARY)==1.0) i3=1;
	suma= i0 + i2 + i3;
	
	if(suma==3)
	  {
	    
	    PointsOfFSTriangle[0][0]=GetGeometry()[0].X();
	    PointsOfFSTriangle[0][1]=GetGeometry()[0].Y();
	    PointsOfFSTriangle[0][2]=GetGeometry()[0].Z();
	    
	    PointsOfFSTriangle[1][0]=GetGeometry()[3].X();
	    PointsOfFSTriangle[1][1]=GetGeometry()[3].Y();
	    PointsOfFSTriangle[1][2]=GetGeometry()[3].Z();
	    
	    PointsOfFSTriangle[2][0]=GetGeometry()[2].X();
	    PointsOfFSTriangle[2][1]=GetGeometry()[2].Y();
	    PointsOfFSTriangle[2][2]=GetGeometry()[2].Z();
	    
	    
	    area=0.5 * CalculateTriangleArea3D(PointsOfFSTriangle[0], PointsOfFSTriangle[1], PointsOfFSTriangle[2]);
	    //////KRATOS_WATCH(area); 
	    T0=GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);
	    T2=GetGeometry()[2].FastGetSolutionStepValue(TEMPERATURE);
	    T3=GetGeometry()[3].FastGetSolutionStepValue(TEMPERATURE);
	    
	    rRightHandSideVector[0] +=  constant * 4.0 * StefenBoltzmann*(pow(T0,4) )* 0.333333333333* area;
	    rRightHandSideVector[2] +=  constant * 4.0 * StefenBoltzmann*(pow(T2,4) )* 0.333333333333* area ;
	    rRightHandSideVector[3] +=  constant * 4.0 * StefenBoltzmann*(pow(T3,4) )* 0.333333333333* area; 
	    
	    rLeftHandSideMatrix(0,0) +=  constant * 0.333333333333 * area;
	    rLeftHandSideMatrix(2,2) +=  constant * 0.333333333333 * area;
	    rLeftHandSideMatrix(3,3) +=  constant * 0.333333333333 * area;
	    
	  }
	
	i0=0;
	i1=0;
	i2=0;
	i3=0;	
	suma=0;
	if(GetGeometry()[0].FastGetSolutionStepValue(IS_BOUNDARY)==1.0) i0=1; 
	if(GetGeometry()[1].FastGetSolutionStepValue(IS_BOUNDARY)==1.0) i1=1;
	if(GetGeometry()[3].FastGetSolutionStepValue(IS_BOUNDARY)==1.0) i3=1;
	suma= i0 + i1 + i3;
	
	if(suma==3)
	  {
	    
	    //KRATOS_ERROR(std::logic_error,  "method not implemented" , "");  
	    PointsOfFSTriangle[0][0]=GetGeometry()[0].X();
	    PointsOfFSTriangle[0][1]=GetGeometry()[0].Y();
	    PointsOfFSTriangle[0][2]=GetGeometry()[0].Z();
	    
	    PointsOfFSTriangle[1][0]=GetGeometry()[1].X();
	    PointsOfFSTriangle[1][1]=GetGeometry()[1].Y();
	    PointsOfFSTriangle[1][2]=GetGeometry()[1].Z();
	    
	    PointsOfFSTriangle[2][0]=GetGeometry()[3].X();
	    PointsOfFSTriangle[2][1]=GetGeometry()[3].Y();
	    PointsOfFSTriangle[2][2]=GetGeometry()[3].Z();
	    
	    area=0.5 * CalculateTriangleArea3D(PointsOfFSTriangle[0], PointsOfFSTriangle[1], PointsOfFSTriangle[2]);
	    //////KRATOS_WATCH(area); 
	    
	    T0=GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);
	    T1=GetGeometry()[1].FastGetSolutionStepValue(TEMPERATURE);
	    T3=GetGeometry()[3].FastGetSolutionStepValue(TEMPERATURE);
	    
	    
	    rRightHandSideVector[0] +=  constant * 4.0 * StefenBoltzmann*(pow(T0,4) )* 0.333333333333* area;  
	    rRightHandSideVector[1] +=  constant * 4.0 * StefenBoltzmann*(pow(T1,4) )* 0.333333333333* area;   
	    rRightHandSideVector[3] +=  constant * 4.0 * StefenBoltzmann*(pow(T3,4) )* 0.333333333333* area;  
	    
	    rLeftHandSideMatrix(0,0) += constant * 0.333333333333 * area;
	    rLeftHandSideMatrix(1,1) += constant * 0.333333333333 * area;
	    rLeftHandSideMatrix(3,3) += constant * 0.333333333333 * area;
	    
	    
	  }
	i0=0;
	i1=0;
	i2=0;
	i3=0;	
	suma=0;
	if(GetGeometry()[0].FastGetSolutionStepValue(IS_BOUNDARY)==1.0) i0=1; 
	if(GetGeometry()[1].FastGetSolutionStepValue(IS_BOUNDARY)==1.0) i1=1;
	if(GetGeometry()[2].FastGetSolutionStepValue(IS_BOUNDARY)==1.0) i2=1;
	suma= i0 + i2 + i1;
	
	if(suma==3)
	  {
	    
	    PointsOfFSTriangle[0][0]=GetGeometry()[0].X();
	    PointsOfFSTriangle[0][1]=GetGeometry()[0].Y();
	    PointsOfFSTriangle[0][2]=GetGeometry()[0].Z();

	    PointsOfFSTriangle[1][0]=GetGeometry()[2].X();
	    PointsOfFSTriangle[1][1]=GetGeometry()[2].Y();
	    PointsOfFSTriangle[1][2]=GetGeometry()[2].Z();
	    
	    PointsOfFSTriangle[2][0]=GetGeometry()[1].X();
	    PointsOfFSTriangle[2][1]=GetGeometry()[1].Y();
	    PointsOfFSTriangle[2][2]=GetGeometry()[1].Z();
	    
	    area=0.5*CalculateTriangleArea3D(PointsOfFSTriangle[0], PointsOfFSTriangle[1], PointsOfFSTriangle[2]);
	    //////KRATOS_WATCH(area); 
	    
	    
	    T0=GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);
	    T1=GetGeometry()[1].FastGetSolutionStepValue(TEMPERATURE);
	    T2=GetGeometry()[2].FastGetSolutionStepValue(TEMPERATURE);
	    
	    rRightHandSideVector[0] +=  constant * 4.0 * StefenBoltzmann*(pow(T0,4) )*0.333333333333*  area;
	    rRightHandSideVector[1] +=  constant * 4.0 * StefenBoltzmann*(pow(T1,4) )* 0.333333333333* area;	
	    rRightHandSideVector[2] +=  constant * 4.0 * StefenBoltzmann*(pow(T2,4) )*0.333333333333*  area;
	    
	    rLeftHandSideMatrix(0,0) += constant * 0.333333333333 * area;
	    rLeftHandSideMatrix(1,1) += constant * 0.333333333333 * area;
	    rLeftHandSideMatrix(2,2) += constant * 0.333333333333 * area;
	    
	    
	  }
	
	msAuxVec[0]=pow(T0,4); 
	msAuxVec[1]=pow(T1,4);
	msAuxVec[2]=pow(T2,4);
	msAuxVec[3]=pow(T3,4);
	
	for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
	  {
	    //getting informations for integration
	    //double IntegrationWeight = integration_points[PointNumber].Weight();
	    //
	    MathUtils<double>::InvertMatrix(J0[PointNumber],mInvJ0[PointNumber],mDetJ0[PointNumber]);
	    double Weight = integration_points[PointNumber].Weight()* mDetJ0[PointNumber];
	    const Vector& N=row(Ncontainer,PointNumber);
	    MathUtils<double>::InvertMatrix(J0[PointNumber],mInvJ0[PointNumber],mDetJ0[PointNumber]);
	    noalias(msDN_DX) = prod(DN_De[PointNumber],mInvJ0[PointNumber]);
	    
	    Cuatro=ZeroMatrix(4,4);
	    Cuatro(0,0)=N[0]*N[0]; //+ N[0]*N[1] + N[0]*N[2] + N[0]*N[3];
	    Cuatro(0,1)+=N[0]*N[1];
	    Cuatro(0,2)+=N[0]*N[2];
	    Cuatro(0,3)+=N[0]*N[3];
	    
	    Cuatro(1,1)=N[1]*N[0];
	    Cuatro(1,1)+=N[1]*N[1];// + N[1]*N[0] + N[1]*N[2] + N[1]*N[3];
	    Cuatro(1,1)+=N[1]*N[2];
	    Cuatro(1,1)+=N[1]*N[3];
	    
	    Cuatro(2,2)=N[2]*N[0];
	    Cuatro(2,2)+=N[2]*N[1];
	    Cuatro(2,2)+=N[2]*N[2];// + N[2]*N[0] + N[2]*N[1] + N[2]*N[3];
	    Cuatro(2,2)+=N[2]*N[3];
	    
	    Cuatro(3,3)=N[3]*N[0];
	    Cuatro(3,3)+=N[3]*N[1];
	    Cuatro(3,3)+=N[3]*N[2];
	    Cuatro(3,3)+=N[3]*N[3];// + N[3]*N[0] + N[3]*N[1] + N[3]*N[2];
	    
	    
	    //	    noalias(rLeftHandSideMatrix) += 1.0 * Weight * absorptioncoefficient * Cuatro; //LO COMENTO AHORA
	    noalias(rRightHandSideVector) += 4.0 * absorptioncoefficient * StefenBoltzmann * prod(Cuatro,msAuxVec)  * Weight;
	    
	  }
	
	
      }
    
    for(unsigned int iii = 0; iii<number_of_points; iii++)
      ms_temp_vec_np[iii] = GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar); 
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp_vec_np);
    
    //KRATOS_WATCH(rLeftHandSideMatrix);
    //KRATOS_WATCH(rRightHandSideVector);
    KRATOS_CATCH("");
  }
  
  //************************************************************************************
  //************************************************************************************
  void Rad3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
  }	 
  
  //************************************************************************************
  //************************************************************************************
  // this subroutine calculates the nodal contributions for the explicit steps of the 
  // fractional step procedure
  void Rad3D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
  {
    KRATOS_TRY
      
      KRATOS_CATCH("");
  }

  
  //************************************************************************************
  //************************************************************************************
  void Rad3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
  {
    
    RadiationSettings::Pointer my_settings = CurrentProcessInfo.GetValue(RADIATION_SETTINGS);
    
    const Variable<double>& rUnknownVar= my_settings->GetUnknownVariable();
    unsigned int number_of_nodes = GetGeometry().PointsNumber();

    
    if(rResult.size() != number_of_nodes)
      rResult.resize(number_of_nodes,false);	

    for (unsigned int i=0;i<number_of_nodes;i++)
      rResult[i] = GetGeometry()[i].GetDof(rUnknownVar).EquationId();
  }
  
  //************************************************************************************
  //************************************************************************************
  void Rad3D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
  {
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    RadiationSettings::Pointer my_settings = CurrentProcessInfo.GetValue(RADIATION_SETTINGS);
    
    const Variable<double>& rUnknownVar= my_settings->GetUnknownVariable();
    
    
    if(ElementalDofList.size() != number_of_nodes)
      ElementalDofList.resize(number_of_nodes);	
    
    for (unsigned int i=0;i<number_of_nodes;i++)
      ElementalDofList[i] = GetGeometry()[i].pGetDof(rUnknownVar);
    
  }
  
  inline double Rad3D::CalculateH(boost::numeric::ublas::bounded_matrix<double, 4, 3 > & DN_DX, double Volume)
  {
    
    double inv_h_max = 0.0;
    for(unsigned int i=0; i<4; i++)
      {
	double inv_h = 0.0;
	for(unsigned int k=0; k<3; k++)
	  inv_h += DN_DX(i,k)*DN_DX(i,k);
	
	if(inv_h > inv_h_max) inv_h_max = inv_h;
      }
    inv_h_max = sqrt(inv_h_max);
    double h = 1.0/inv_h_max;

    return h ;
  }
  
  double Rad3D::CalculateTriangleArea3D(	array_1d<double,3>& Point1, array_1d<double,3>& Point2, array_1d<double,3>& Point3	)
  {
    //Heron's formula
    double a=Length(Point1, Point2);//sqrt((Point1[0]-Point2[0])*(Point1[0]-Point2[0]) + (Point1[1]-Point2[1])*(Point1[1]-Point2[1]) +(Point1[2]-Point2[2])*(Point1[2]-Point2[2]));
    double b=Length(Point1, Point3);//sqrt((Point3[0]-Point2[0])*(Point3[0]-Point2[0]) + (Point3[1]-Point2[1])*(Point3[1]-Point2[1]) +(Point3[2]-Point2[2])*(Point3[2]-Point2[2]));
    double c=Length(Point2, Point3);//sqrt((Point1[0]-Point3[0])*(Point1[0]-Point3[0]) + (Point1[1]-Point3[1])*(Point1[1]-Point3[1]) +(Point1[2]-Point3[2])*(Point1[2]-Point3[2]));
    double p=0.5*(a+b+c);
    return sqrt(p*(p-a)*(p-b)*(p-c));
  }
  
  double Rad3D::Length(array_1d<double,3>& Point1, array_1d<double,3>& Point2)
  {
    //////KRATOS_WATCH("length calculation")
    return sqrt((Point1[0]-Point2[0])*(Point1[0]-Point2[0]) + (Point1[1]-Point2[1])*(Point1[1]-Point2[1]) +(Point1[2]-Point2[2])*(Point1[2]-Point2[2]));
  }
  void Rad3D::qi( const int ndivisionsp, std::vector< Matrix > edges_tauxp, std::vector< Matrix > nodes_auxp, std::vector< Matrix > rGradientauxp,  array_1d<double,6> conductivitiesp, Matrix& Kaux1p,Matrix& Lenrichaux1p)
  {
    array_1d<double,3> v3;
    array_1d<double,3> v4;
    array_1d<double,3> An;// =zero;
    
    array_1d<double,3>  Point1; 
    array_1d<double,3>  Point2; 
    array_1d<double,3>  Point3;
    array_1d<double,3>  normal=ZeroVector(3);
    array_1d<double,12> msAuxVec = ZeroVector(4); //dimension = number of nodes_aux
    
    //double area_normal=0.0;
    array_1d<double,3> area_normal;
    double c0, c1, c2;
    double norm_c;	      
    double norm_u; 
    array_1d<int, 3 > nodes; 
    nodes(0)=-1;nodes(1)=-1;nodes(2)=-1;
    
    for(unsigned int i=0; i<ndivisionsp;i++)
      for(unsigned int j=0; j<5;j++) //5 6 j<4
	if(edges_tauxp[i](j,0)>-1)
	  {
	    
	    for(unsigned int d=0; d<3;d++)
	      nodes(d)=edges_tauxp[i](j,d);
	    
	    for(unsigned int d=0; d<3;d++){
	     Point1(d)=nodes_auxp[nodes(0)](0,d);
	     Point2(d)=nodes_auxp[nodes(1)](0,d);
	     Point3(d)=nodes_auxp[nodes(2)](0,d);
	    }
	    
	    double trianglearea=CalculateTriangleArea3D( Point1, Point2, Point3);
	    
	    if(trianglearea<=0.0)
	      {
		//////KRATOS_WATCH("area_cero");
		////KRATOS_WATCH("area_cero");
		////KRATOS_WATCH("area_cero");
		////KRATOS_WATCH("area_cero");
		////KRATOS_WATCH(Point1);
		////KRATOS_WATCH(Point2);
		////KRATOS_WATCH(Point3);
		////KRATOS_WATCH(trianglearea);
	       KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
	      }
	    array_1d<double,3> v1;
	    array_1d<double,3> v2;
	   array_1d<double,3> An;// =zero;
	   //double area_normal=0.0;
	   array_1d<double,3> area_normal;
	   
	   v1 = Point2 - Point1;
	   
	   v2 = Point3 - Point2;
	   
	   v1=Point1-Point3;
	   v2=Point2-Point1;
	      
	   MathUtils<double>::CrossProduct(area_normal,v1,v2);
	   msAuxVec = ZeroVector(3);
	   c0 = fabs(area_normal[0]);
	   c1 = fabs(area_normal[1]);
	   c2 = fabs(area_normal[2]);
	   //////KRATOS_WATCH(area_normal);
	   //KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
	   msAuxVec[0]=c0;
	   msAuxVec[1]=c1;
	   msAuxVec[2]=c2;
	      
	   norm_u = msAuxVec[0]*msAuxVec[0] + msAuxVec[1]*msAuxVec[1] + msAuxVec[2]*msAuxVec[2];
	   norm_c =sqrt(norm_u);
	      
	   if(norm_c<=0.00000000000001)
	     {
	       
	       ////KRATOS_WATCH("norm_c_cero");
	       ////KRATOS_WATCH("norm_c_cero");
	       ////KRATOS_WATCH("norm_c_cero");
	       ////KRATOS_WATCH("norm_c_cero");
	       ////KRATOS_WATCH(Point1);
	       ////KRATOS_WATCH(Point2);
	       ////KRATOS_WATCH(Point3);
	       ////KRATOS_WATCH(ndivisionsp);
	       ////KRATOS_WATCH(area_normal);
	       ////KRATOS_WATCH(trianglearea);
	       ////KRATOS_WATCH(msAuxVec);
	       ////KRATOS_WATCH(norm_c);
	       KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
	       //sssss
	     }
	      
	   for(unsigned int d=0; d<3;d++){
	     normal(d) = area_normal[d]/ norm_c;
	   }   
	   /*normal(0) = area_normal[0]/ norm_c;
	     normal(1) = area_normal[1]/ norm_c;
	     normal(2) = area_normal[2]/ norm_c;*/
	   
	   
	   int max;
	   if(ndivisionsp==2) max=1;
	   else if(ndivisionsp==3) max=2;
	   else if(ndivisionsp==4) max=3;
	   else max=4;
	   /*////KRATOS_WATCH("max");
	     ////KRATOS_WATCH(max);
	     ////KRATOS_WATCH("normal");
	     ////KRATOS_WATCH(normal);
	     ////KRATOS_WATCH("trianglearea");
	     ////KRATOS_WATCH(trianglearea);*/
	   //////KRATOS_WATCH(Kaux);
	   //double z, xx;
	   ////////////////////////
	   /////////////////////////
	   for(unsigned int zz = 0; zz < 3; zz++)
	     {
	       if(edges_tauxp[i](j,zz)>3)
		 {
		   int position=edges_tauxp[i](j,zz);
		   for (unsigned int t = 0; t < 4; t++){ 
		     //(Kaux1(position,j));
		     //////KRATOS_WATCH(Kaux1);
		     for(unsigned int g=0; g<3; g++)
		       Kaux1p(position,t) +=(rGradientauxp[i](t,g) * normal(g)) *trianglearea * 0.33333333333333333333333333333 * 1.0/(3.0 * conductivitiesp[i]);
		   }
		   for (unsigned int t = 0; t < max; t++){
		     for(unsigned int g=0; g<3; g++)
		       Lenrichaux1p(position,t+4) +=(rGradientauxp[i](t+4,g) * normal(g) ) *trianglearea * 0.33333333333333333333333333333 * 1.0/(3.0 * conductivitiesp[i]);
		   }
		 }
	     }
	  }
  }
  /*void Rad3D::qi( const int ndivisionsp, std::vector< Matrix > edges_tauxp, std::vector< Matrix > nodes_auxp, std::vector< Matrix > rGradientauxp,  array_1d<double,6> conductivitiesp, Matrix& Kaux1p,Matrix& Lenrichaux1p)
    {
    array_1d<double,3> v3;
    array_1d<double,3> v4;
    array_1d<double,3> An;// =zero;
    
    array_1d<double,3>  Point1; 
    array_1d<double,3>  Point2; 
    array_1d<double,3>  Point3;
    array_1d<double,3>  normal=ZeroVector(3);
    array_1d<double,12> msAuxVec = ZeroVector(4); //dimension = number of nodes_aux
    
    //double area_normal=0.0;
    array_1d<double,3> area_normal;
    double c0, c1, c2;
    double norm_c;	      
    double norm_u; 
    array_1d<int, 3 > nodes; 
    nodes(0)=-1;nodes(1)=-1;nodes(2)=-1;
    
    for(unsigned int i=0; i<ndivisionsp;i++)
    for(unsigned int j=0; j<5;j++) //5 6 j<4
    if(edges_tauxp[i](j,0)>-1)
    {
    //int a=edges_taux[i](j,0);
    //int b=edges_taux[i](j,1);
    //int c=edges_taux[i](j,2);
    for(unsigned int d=0; d<3;d++)
    nodes(d)=edges_tauxp[i](j,d);
    
    for(unsigned int d=0; d<3;d++){
    Point1(d)=nodes_auxp[nodes(0)](0,d);
    Point2(d)=nodes_auxp[nodes(1)](0,d);
    Point3(d)=nodes_auxp[nodes(2)](0,d);
    }
    
    double trianglearea=CalculateTriangleArea3D( Point1, Point2, Point3);
    
    if(trianglearea<=0.0)
    {
    KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
    }
    array_1d<double,3> v1;
    array_1d<double,3> v2;
    array_1d<double,3> An;// =zero;
    //double area_normal=0.0;
    array_1d<double,3> area_normal;
    
    v1 = Point2 - Point1;
    
    v2 = Point3 - Point2;
    
    v1=Point1-Point3;
    v2=Point2-Point1;
    double soff=1.0;
    
    MathUtils<double>::CrossProduct(area_normal,v1,v2);
    msAuxVec = ZeroVector(3);
    c0 = fabs(area_normal[0]);
    c1 = fabs(area_normal[1]);
    c2 = fabs(area_normal[2]);
    msAuxVec[0]=c0;
    msAuxVec[1]=c1;
    msAuxVec[2]=c2;
    
    norm_u = msAuxVec[0]*msAuxVec[0] + msAuxVec[1]*msAuxVec[1] + msAuxVec[2]*msAuxVec[2];
    norm_c =sqrt(norm_u);
    
    if(norm_c<=0.00000000000001)
    {
    
    KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
    }
    
    for(unsigned int d=0; d<3;d++){
    normal(d) = area_normal[d]/ norm_c;
    }   
    
    
    int max;
    if(ndivisionsp==2) max=1;
    else if(ndivisionsp==3) max=2;
    else if(ndivisionsp==4) max=3;
    else max=4;
    double z, xx;
    ////////////////////////
    /////////////////////////
    for(unsigned int zz = 0; zz < 3; zz++)
    {
    if(edges_tauxp[i](j,zz)>3)
    {
    int position=edges_tauxp[i](j,zz);
    for (unsigned int t = 0; t < 4; t++){ 
    for(unsigned int g=0; g<3; g++)
    Kaux1p(position,t) +=soff*(rGradientauxp[i](t,g) * normal(g)) *trianglearea * 0.33333333333333333333333333333* 1.0/(3.0 * conductivitiesp[i]);
    }
    for (unsigned int t = 0; t < max; t++){
    for(unsigned int g=0; g<3; g++)
    Lenrichaux1p(position,t+4) +=soff*(rGradientauxp[i](t+4,g) * normal(g) ) *trianglearea * 0.33333333333333333333333333333 * 1.0/(3.0 * conductivitiesp[i]);
    }
    }
    }
    }
    }*/
  
  void Rad3D::Heat_Source(VectorType& rRightHandSideVector, const int ndivisionsp, array_1d<double,6>& volumesp, array_1d<double,6>& conductivitiesp, boost::numeric::ublas::bounded_matrix<double,6, 8 >& Ngaussnewp, const double heatp)
  {
    double coefficientaux=0.0;
    double heat_flux=0.0;
    
    for(unsigned int zz = 0; zz < 4; zz++) heat_flux += GetGeometry()[zz].FastGetSolutionStepValue(HEAT_FLUX);
    heat_flux *= 0.25;// *0.0;
    
    for (unsigned int i=0; i<ndivisionsp ; i++)
      {
        coefficientaux =  heat_flux * volumesp[i];// * conductivitiesp[i];
	rRightHandSideVector[0] += coefficientaux * Ngaussnewp(i,0); // * conductivities[i];
	rRightHandSideVector[1] += coefficientaux * Ngaussnewp(i,1); // * conductivities[i];
	rRightHandSideVector[2] += coefficientaux * Ngaussnewp(i,2); // * conductivities[i];
	rRightHandSideVector[3] += coefficientaux * Ngaussnewp(i,3); // * conductivities[i];  
      }
    
  }
  template<class T>
  
  bool Rad3D::InvertMatrix(const T& input, T& inverse)
    
  {
    
    typedef permutation_matrix<std::size_t> pmatrix;
    
    // create a working copy of the input
    
    T A(input);
    
    // create a permutation matrix for the LU-factorization
    
    pmatrix pm(A.size1());
    
    // perform LU-factorization
    
    int res = lu_factorize(A, pm);
    
    if (res != 0)
      
      return false;
    
    // create identity matrix of "inverse"
    
    inverse.assign(identity_matrix<double> (A.size1()));
    
    // backsubstitute to get the inverse
    
    lu_substitute(A, pm, inverse);
    
    return true;
    
  }
  
} // Namespace Kratos


