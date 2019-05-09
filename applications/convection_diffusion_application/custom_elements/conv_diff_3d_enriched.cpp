//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Julio Marti
//

#include "includes/define.h"
#include "custom_elements/conv_diff_3d_enriched.h"
#include "convection_diffusion_application.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "../../pfem_2_application/custom_utilities/enrichmentutilities.h"
#include "utilities/discont_utils.h"
#include "includes/kratos_flags.h"
#include "includes/deprecated_variables.h"

//#define FACE_HEAT_FLUX_ON_EDGES
//#define FACE_HEAT_FLUX
//#define NOT_STATIONARY
#define FACE_HEAT_FLUX_ON_FREE_SURFACE

#define RAD
//#define RAD_ON_CUT_ELEMENT
namespace Kratos
{

//************************************************************************************
//************************************************************************************

ConvDiff3Denriched::ConvDiff3Denriched(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
  ConvDiff3Denriched::IntegrationMethod ConvDiff3Denriched::GetIntegrationMethod()
  {
    return mThisIntegrationMethod;
  }
//************************************************************************************
//************************************************************************************

ConvDiff3Denriched::ConvDiff3Denriched(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

Element::Pointer ConvDiff3Denriched::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
  //return Kratos::make_shared<ConvDiff3Denriched>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    return Kratos::make_intrusive<ConvDiff3Denriched>(NewId, GetGeometry().Create(ThisNodes), pProperties);	
    //return Element::Pointer(new ConvDiff3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

Element::Pointer ConvDiff3Denriched::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
  //return Kratos::make_shared<ConvDiff3Denriched>(NewId, pGeom, pProperties);
    return Kratos::make_intrusive<ConvDiff3Denriched>(NewId, pGeom, pProperties);
}

ConvDiff3Denriched::~ConvDiff3Denriched()
{
}



//************************************************************************************
//************************************************************************************

void ConvDiff3Denriched::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
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

  //const unsigned int LocalSize = GetGeometry().size();
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
    
  //const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
  //calculating parameter tau
  //double c10 = 4.00;
  //double c20 = 2.00;
  //double h = pow(6.00 * Area, 0.3333333);
  //double norm_uv = norm_2(ms_vel_gauss);
    
  if(rUnknownVar == TEMPERATURE or rUnknownVar == YCH4) 
    {
      array_1d<double,4> distances;
      bool has_negative_node=false;
      bool has_positive_node=false;
      //bool has_distance_0=false;
      //for radiation
      //double absorptioncoefficient2;
      //const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);
	
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
    //double norm_grad;
    //double res;
    //double k_aux;
      
    //split_element=false;
		if (split_element==true){
		//KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "") 
		///NEW CHANGES
		
		double unsigned_distance0=fabs(distances[0]);
		double unsigned_distance1=fabs(distances[1]);
		double unsigned_distance2=fabs(distances[2]);
		double unsigned_distance3=fabs(distances[3]);
		
		double longest_distance=fabs(unsigned_distance0);
		if (unsigned_distance1>longest_distance) longest_distance=unsigned_distance1;
		if (unsigned_distance2>longest_distance) longest_distance=unsigned_distance2;
		if (unsigned_distance3>longest_distance) longest_distance=unsigned_distance3;
		
	
		double tolerable_distance =longest_distance*0.001;
		tolerable_distance=longest_distance * 0.001; //1.0;
					
    if(distances[0]==0) KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "")  
		if(distances[1]==0) KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "")  
		if(distances[2]==0) KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "")  
		if(distances[3]==0) KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "")  
	    
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
	    
	    else{
	      distances[3]=tolerable_distance;
	    }
	  }
	
	boost::numeric::ublas::bounded_matrix<double,6, 2> Nenriched;
	array_1d<double,6>  volumes(6);
	boost::numeric::ublas::bounded_matrix<double,4, 3 > coords;
	boost::numeric::ublas::bounded_matrix<double,6, 4 > Ngauss;
	boost::numeric::ublas::bounded_matrix<double,6, 8 > Ngaussnew;
	  
	array_1d<double,6>  signs(6);
	array_1d<double,6>  edge_areas(6);
	std::vector< Matrix > gauss_gradients(6);
	std::vector< Matrix > rGradientaux(6);
	std::vector< Matrix > edges_taux(8);
      
	std::vector< Matrix > gauss_gradients_aux(8);
	std::vector< Matrix > nodes_aux(8);
	std::vector< Matrix > tot_edges(8);
	boost::numeric::ublas::bounded_matrix<double, 4, 4 > Laplacian_matrix = ZeroMatrix(4,4);
	//boost::numeric::ublas::bounded_matrix<double, 4, 4 > Convective_matrix = ZeroMatrix(4,4);
	//boost::numeric::ublas::bounded_matrix<double, 4, 4 > Convective_matrix_tau = ZeroMatrix(4,4);
	boost::numeric::ublas::bounded_matrix<double, 4, 4 > Aux_matrix = ZeroMatrix(4,4);
	  
	boost::numeric::ublas::bounded_matrix<double, 4, 4 > Mass_Aux = ZeroMatrix(4,4);	
	array_1d<double,6> conductivities;
	array_1d<double,6> specific_heat;
	array_1d<double,6> volumetric_heat_capacities;
	array_1d<double,6> densities;
	  
	array_1d<double,4> distances1;
	array_1d<double,6>  volumes1(6);
	array_1d<double,6>  N_Star(6);
	boost::numeric::ublas::bounded_matrix<double,4, 3 > coords1;
	boost::numeric::ublas::bounded_matrix<double,6, 4 > Ngauss1;
	boost::numeric::ublas::bounded_matrix<double,6, 4> Nenriched1;
	array_1d<double,6>  signs1(6);
	array_1d<double,6>  edge_areas1(6);
	std::vector< Matrix > gauss_gradients1(6);
	std::vector< Matrix > mass_matrix(6);
	std::vector< Matrix > mass_matrixone(6);
	for (unsigned int i = 0; i < 6; i++)
	  {
	    mass_matrix[i].resize(4, 4, false);  
	    mass_matrix[i] *=0.0;
	    //mass_matrix[i].resize(8, 8, false);  
	    //mass_matrix[i] *=0.0;
	  }  
	
	bool switch_off_e;	
	switch_off_e=false;
	std::vector< Matrix > coord_interface(4);
	  
	int number_interface_elements_1=0;
	boost::numeric::ublas::bounded_matrix<double,2 , 3 >  coord_interface_nodes_1; //used to p
	array_1d<double,6> area_interface_1; 
	array_1d<double,6> area_inter;
	  
	//const array_1d<double, 3 > & x_c;
	array_1d<double, 3 >  x_g;
	//EXTRA
	  
	for (unsigned int i = 0; i < TNumNodes; i++)
	  {
	    const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
	      
	    for (unsigned int j = 0; j < 3; j++)
	      coords(i, j) = xyz[j];
	  }
	for (unsigned int i = 0; i < 6; i++)
	  {
	    gauss_gradients[i].resize(2, 3, false);  //2 values of the 2 shape functions, and derivates in (xyz) direction).
	    gauss_gradients1[i].resize(4, 5, false);
	    gauss_gradients_aux[i].resize(8, 3, false);  //2 values of the 2 shape functions, and derivates in (xyz) direction).
	  }
	  
	for (unsigned int i = 0; i < 6; i++)
	  {
	    rGradientaux[i].resize(8,3,false);
	  }
	  
	for (unsigned int i = 0; i < 4; i++)
	  {
	    coord_interface[i].resize(1,3,false);
	  }
	  
	for (unsigned int i = 0; i < 8; i++)
	  {
	    edges_taux[i].resize(5, 3, false);  //2 values of the 2 shape functions, and derivates in (xyz) direction).
	    nodes_aux[i].resize(1, 3, false);
	    //tot_edges[i].resize(1, 3, false);
	    tot_edges[i].resize(5, 3, false);
	  }
	  
	
	//Vector rEnrichedTemperature_oldit;// =  this->GetValue(ENRICHED_TEMPERATURE_OLDIT);
	//rEnrichedTemperature_oldit.resize(4);
	
	//double t1, t2,t3;
	double temperature=0.0;
	//int iteration_number = 0;//rCurrentProcessInfo[ITERATION_NUMBER];
	for(unsigned int j = 0; j < 4; j++) temperature += GetGeometry()[j].FastGetSolutionStepValue(TEMPERATURE);
	temperature *=0.25;
	  
	std::vector< Matrix > MASS_MATRIX2(6);
	std::vector< Matrix > Cuatros(6);
	
	array_1d<double,6>  wei(6);
	for (unsigned int i = 0; i < 6; i++)
	  {
	    MASS_MATRIX2[i].resize(4, 4, false);  
	    MASS_MATRIX2[i] *=0.0;
	    Cuatros[i].resize(4, 4, false);  
	    Cuatros[i] *=0.0;
	  }  
	int total_nodes=0;
	unsigned int ndivisions = EnrichmentUtilitiesforPFEM2::CalculateEnrichedShapeFuncions_Simplified(this->GetGeometry(),coords, msDN_DX, distances, volumes, Ngauss, signs, gauss_gradients,gauss_gradients_aux,Nenriched, number_interface_elements_1, coord_interface_nodes_1, area_interface_1,area_inter, N_Star,switch_off_e,edges_taux,nodes_aux, tot_edges,rGradientaux, total_nodes, coord_interface,Ngaussnew,mass_matrix, MASS_MATRIX2, wei);//,mass_matrixone);
	  
	double volume=0.0;
	double volumepositive=0.0;
	double volumenegative=0.0;
	for (unsigned int i = 0; i!=ndivisions; i++) //the 3 partitions of the triangle
	  {
	    if (signs(i)>0)
	      {
		volume +=volumes[i];
		volumepositive +=volumes[i];
	      }
	    else
	      {  
		volume +=volumes[i];
		volumenegative +=volumes[i];
	      }
	  }	
	  
      
	    Matrix Laplacian_matrix_mod;
	    //int size=ndivisions -1;
	    Laplacian_matrix_mod.resize(total_nodes, total_nodes);
	    Laplacian_matrix_mod = ZeroMatrix(total_nodes,total_nodes);
	      
	      
	    /****************************/
	    double density_elem=0.0;
	    double specific_heat_elem=0.0;
	    array_1d<double,6> absorptioncoefficient;
	    //double volume=0.0;
	    //double volumepositive=0.0;
	    //double volumenegative=0.0;
	    for (unsigned int i = 0; i!=ndivisions; i++) //the 3 partitions of the triangle
	      {
		if (signs(i)>0)
		  {
		    conductivities[i] = 0.0131; //1.0;
		    densities[i]=1.0;
		    specific_heat[i]=1310.0; //1.0
		    absorptioncoefficient[i]=0.01;
		    density_elem+=densities[i]* volumes[i];
		    specific_heat_elem+=specific_heat[i] * volumes[i];
		  }
		else
		  {
		    conductivities[i] = 100.0; //0.12;//0.16;//100.0;
		    densities[i]=905.0;
		    specific_heat[i]=1900.0; 
		    absorptioncoefficient[i]=3000.0;//1000.0;
		    density_elem+=densities[i] * volumes[i];
		    specific_heat_elem+=specific_heat[i] * volumes[i];
		  }
	      }	
	    
	    //Mass_matrix = ZeroMatrix(4, 4);  
	    boost::numeric::ublas::bounded_matrix<double, 4, 3 > msDN_DX_aux;
	    array_1d<double, 4 > msN_aux;
	      
	    //tau1*=0.0;
	    //boost::numeric::ublas::bounded_matrix<double, 4, 4 > mass_matrix_tau = ZeroMatrix(4,4);
	      
	    for (unsigned int i=0; i<ndivisions ; i++)
	      {
		Laplacian_matrix_mod +=prod (rGradientaux[i], trans(rGradientaux[i]))* conductivities[i] * volumes[i];
	      }
	    //KRATOS_WATCH(Laplacian_matrix_mod);
        
        
	    Matrix Kaux;
	    Matrix Kaux1;
	    
	    Kaux.resize(total_nodes - 4, 4);
	    Kaux = ZeroMatrix(total_nodes - 4, 4);  
	    Kaux1.resize(total_nodes,total_nodes);
	    Kaux1 = ZeroMatrix(total_nodes,total_nodes);  
	      
	    Matrix Lenrichaux;
	    Matrix Lenrichaux1;
	      
	    Lenrichaux.resize(total_nodes - 4,total_nodes -4);
	    Lenrichaux = ZeroMatrix(total_nodes -4,total_nodes-4);
	    Lenrichaux1.resize(total_nodes,total_nodes);
	    Lenrichaux1 = ZeroMatrix(total_nodes,total_nodes); 
	      
	    array_1d<double,3>  Point1; 
	    array_1d<double,3>  Point2; 
	    array_1d<double,3>  Point3;

	    array_1d<double,3>  Point4; 
	    array_1d<double,3>  Point5; 
	    array_1d<double,3>  Point6;
	      
	    array_1d<double,3>  normal=ZeroVector(3);
	      
	    //////////////////////////////////////////////////
	    //////////////////////////////////////////////////
	    ////////////////////////
	    //NORMAL
	    array_1d<double,3> v3;
	    array_1d<double,3> v4;
	    array_1d<double,3> An;// =zero;
	    //double area_normal=0.0;
	    array_1d<double,3> area_normal;
	    //double c0, c1, c2;
	    //double norm_c;	      
	    //double norm_u; 
	    array_1d<int, 3 > nodes; 
	    nodes(0)=-1;nodes(1)=-1;nodes(2)=-1;
	      
	    ///flujossssssssssssssssssssss
	    //KRATOS_WATCH("FLUJOS");
	    Kaux1 *=0.0;
	    Lenrichaux1 *=0.0;
	    qi( ndivisions, edges_taux, nodes_aux, rGradientaux, conductivities, Kaux1, Lenrichaux1);
	      
        //KRATOS_WATCH(Kaux1);
        //KRATOS_WATCH(Lenrichaux1);
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
	    // if(ndivisions==2 or ndivisions==3) KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
	    max=total_nodes-4;
	      
	      
	    for(unsigned int m = 0; m < max; m++)
	      for (unsigned int s = 0; s < max; s++) 
		Lenrichaux(m,s) = Lenrichaux1(4+m,s+4) ;
	   
	    for(unsigned int m = 0; m < max; m++)
	      for (unsigned int s = 0; s < 4; s++) 
		Kaux(m,s) = Kaux1(4+m,s);
	   
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
	   
	   
	    int ninter=0;
	   
	    double heat=1000.0;
	    heat *=0.0;
	    rRightHandSideVector *=0.0;
	    //Heat_Source(rRightHandSideVector, ndivisions, volumes, conductivities, Ngaussnew, heat);
	    //double coefficientaux=0.0; 

	    if(ndivisions==3 or ndivisions==2)  KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
      
	    //////////////////////////////////
	    if(ndivisions==4) ninter=1;
	    if(ndivisions==6) ninter=2;
	    if(ndivisions==3) ninter=0;
	    if(ndivisions==2) ninter=0;
	    int number_of_nodes=0;
	    if(ndivisions==2) {number_of_nodes=1; /*KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");*/}
	    if(ndivisions==3) {number_of_nodes=2; /*KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");*/}
	    if(ndivisions==4) number_of_nodes=3;
	    if(ndivisions==6) number_of_nodes=4;
      
	    //double qc=100000.0;

	    boost::numeric::ublas::bounded_matrix<double, 4, 4 > Mass_matrix = ZeroMatrix(4,4);
	    Mass_matrix = ZeroMatrix(4, 4); 
	    
	    for (unsigned int i=0; i<ndivisions ; i++)
	      {
		for (unsigned int j = 0; j < 4; j++) 
		  {		
		    Mass_matrix(j,j) += densities(i) * specific_heat(i) * Ngaussnew(i,j) * volumes(i) / Dt; //original
		    //msN_aux[j]= Ngaussnew(i,j);
		  }
	      }
	    
	    
	    std::vector< Matrix > area_interface(4);
	    for (unsigned int i = 0; i < 4; i++)
	      {
		area_interface[i].resize(1,1,false);
	      }    
	    
	    Vector Flux;
	    Flux.resize(number_of_nodes);
	    Flux=ZeroVector(number_of_nodes);
	    
	    //double aux_coeff=1.0;///(3.0*0.01);

   
	    //double area=0.0;
	    //double T0,T1,T2,T3;
	    //double suma1=0.0;
	    //int i00,i11,i22,i33;
	    std::vector<array_1d<double,3> > PointsOfFSTriangle1;
	    PointsOfFSTriangle1.reserve(3);
	    //double qcp1;
	    //double FHF=0.0;
	    
	    
	
		//array_1d<double,3> v1;
		//array_1d<double,3> v2;
		//array_1d<double,3> An;
		//array_1d<double,3> area_normal;
		
		//Point1=GetGeometry()[1];
		//Point2=GetGeometry()[2];
		//Point3=GetGeometry()[3];
		
		//v1=Point1-Point3;
		//v2=Point2-Point1;
		
		//MathUtils<double>::CrossProduct(area_normal,v1,v2);
		//msAuxVec = ZeroVector(3);
		//c0 = fabs(area_normal[0]);
		//c1 = fabs(area_normal[1]);
		//c2 = fabs(area_normal[2]);
		
		//msAuxVec[0]=c0;
		//msAuxVec[1]=c1;
		//msAuxVec[2]=c2;
		
		//norm_u = msAuxVec[0]*msAuxVec[0] + msAuxVec[1]*msAuxVec[1] + msAuxVec[2]*msAuxVec[2];
		//norm_c =sqrt(norm_u);
		
		
	  //   Vector Heat;
	  //  Heat.resize(number_of_nodes);
	  //  Heat=ZeroVector(number_of_nodes);
	    
	  //  double hp=0.0;
	  //  int np=0;
	  //  for (unsigned int i=0; i<4 ; i++)
	  //    {
		//if(this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE)<0.0) {hp += GetGeometry()[i].FastGetSolutionStepValue(HEAT_FLUX); np++;}
		//}     

     
	  //  hp *=(1.0/np); 
	    //hp *=0.25;      
                  
	  //  for (unsigned int i=0; i<ndivisions ; i++)
	  //    {
		//coefficientauxs = heat * volumes[i] * conductivities[i];
		//Heat[0] += hp *  Ngaussnew(i,4) * volumes[i];
		//Heat[1] += hp *  Ngaussnew(i,5) * volumes[i];
		//Heat[2] += hp *  Ngaussnew(i,6) * volumes[i];
		//if(number_of_nodes==4) Heat[3] += hp *  Ngaussnew(i,7) * volumes[i];
	  //    }


	    //double pruebas =  0.0;
	    //pruebas =fhf;

       
	  	   
	    Matrix Abar_aux;
	    Abar_aux.resize(total_nodes-4,total_nodes-4);
	   
	    Abar_aux = ZeroMatrix(total_nodes-4,total_nodes-4);
	    Abar_aux=Abar; ///NO TIENE EN CUENTA LA LINEARIZACIONNNNNNNNNNNNNNNNNNNNNNNNNNNN
	   
	    //int a,b,c;
	    double StefenBoltzmann = 5.67e-8;
	    double emissivity = 1.0;
	    //double temperature=0.0;
	    //double aux = pow(298.0,4);
	    //double t0;
	    boost::numeric::ublas::bounded_matrix<double, 4, 4 > Laplacian_matrix_auxiliar = ZeroMatrix(4,4);

	    bool prueba=false;
	    prueba=this->InvertMatrix(Abar,inverse_enrichments);
      

	    if(prueba==false) 
	      {
			KRATOS_WATCH("no se puede invertir");
			KRATOS_WATCH("no se puede invertir");
			KRATOS_WATCH("no se puede invertir");
			KRATOS_WATCH("no se puede invertir");
			KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
		
	      }
	    temp_matrix= prod(Kbar, inverse_enrichments);
      
	  
	    for(unsigned int m = 0; m < 4; m++)
	      for (unsigned int j = 0; j < 4; j++) 
		Laplacian_matrix(m,j)=Laplacian_matrix_mod(m,j);// + Mass_matrix_mod(m,j);
	   
	    noalias(rLeftHandSideMatrix) = Laplacian_matrix;
     
	   
	    for (unsigned int i = 0; i!=ndivisions; i++) //the 3 partitions of the triangle
	      {
		for (unsigned int kkk = 0; kkk<4; kkk++) //gauss points   kkk<4
		  {
		    N[0]=MASS_MATRIX2[i](kkk,0);
		    N[1]=MASS_MATRIX2[i](kkk,1);
		    N[2]=MASS_MATRIX2[i](kkk,2);
		    N[3]=MASS_MATRIX2[i](kkk,3);
		   
		    Cuatros[i](0,0)+=N[0]*N[0]*wei(i);
		    Cuatros[i](0,0)+=N[0]*N[1]*wei(i);
		    Cuatros[i](0,0)+=N[0]*N[2]*wei(i);
		    Cuatros[i](0,0)+=N[0]*N[3]*wei(i);
		   
		    Cuatros[i](1,1)+=N[1]*N[0]*wei(i);
		    Cuatros[i](1,1)+=N[1]*N[1]*wei(i);
		    Cuatros[i](1,1)+=N[1]*N[2]*wei(i);
		    Cuatros[i](1,1)+=N[1]*N[3]*wei(i);
		   
		    Cuatros[i](2,2)+=N[2]*N[0]*wei(i);
		    Cuatros[i](2,2)+=N[2]*N[1]*wei(i);
		    Cuatros[i](2,2)+=N[2]*N[2]*wei(i);
		    Cuatros[i](2,2)+=N[2]*N[3]*wei(i);

		    Cuatros[i](3,3)+=N[3]*N[0]*wei(i);
		    Cuatros[i](3,3)+=N[3]*N[1]*wei(i);
		    Cuatros[i](3,3)+=N[3]*N[2]*wei(i);
		    Cuatros[i](3,3)+=N[3]*N[3]*wei(i);
		   
		  }
	      }
	   
#if defined(NOT_STATIONARY)
	    for (unsigned int i = 0; i!=ndivisions; i++) 
	      {
		noalias(rLeftHandSideMatrix) +=  0.0 * densities[i] * specific_heat[i] * Cuatros[i] / Dt; 
		
	      }
#endif
	   
	    noalias(rLeftHandSideMatrix) += 0.0 * Mass_matrix;
	    //#if defined(FACE_HEAT_FLUX)	   
	    //	   noalias(rLeftHandSideMatrix) += Laplacian_matrix_auxiliar; ///POR SI LA INTERFASE PASA POR ALGUN NODO
	    //#endif    
	   
	   
	    ////////////////////////////////////////////////////////////////
	    ///////////////////////////////////////////////////////////////
	    ///nuevoooooooooooooooooooooooooooooo
	    array_1d<double, 4 > a_dot_grad;
	    noalias(a_dot_grad) = prod(msDN_DX, ms_vel_gauss);
	    array_1d<double, 4 > a_dot_grad_and_mass;
	    a_dot_grad_and_mass = msN;
	    noalias(msAuxMat)= outer_prod(a_dot_grad, a_dot_grad_and_mass);
      

	    noalias(rLeftHandSideMatrix) -=  prod(temp_matrix,A);
	   
	    rRightHandSideVector -=  0.0 * prod(temp_matrix,Flux);
	    //rRightHandSideVector -=  prod(temp_matrix,Heat);
	    
	   
	    Matrix inverse_enrichments_x;
	    inverse_enrichments_x.resize(total_nodes-4,total_nodes-4);
	    inverse_enrichments_x = ZeroMatrix(total_nodes-4,total_nodes-4);
      
	    Matrix temp_matrix_x;
	    temp_matrix_x.resize(4,total_nodes-4);
	    temp_matrix_x = ZeroMatrix(4,total_nodes-4);
	   
	    this->InvertMatrix(Abar_aux,inverse_enrichments_x);
	    temp_matrix_x= prod(Kbar, inverse_enrichments_x);
	   
	    other_matrix = Laplacian_matrix;
	   
#if defined(NOT_STATIONARY)

	    for (unsigned int i = 0; i!=ndivisions; i++) 
	      {
		other_matrix +=  0.0 * densities[i] * specific_heat[i] * Cuatros[i] / Dt; 
	      }
#endif
	   
	    other_matrix += 0.0 * Mass_matrix;//noalias(rLeftHandSideMatrix) +=
	    other_matrix -= 1.0 * prod(temp_matrix_x,A);
	   
#if defined(NOT_STATIONARY)
	   
	    for (unsigned int iii = 0; iii < number_of_points; iii++)
	      ms_temp_vec_np[iii] = GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar,1); 
	   
	    for (unsigned int i = 0; i!=ndivisions; i++) //the 3 partitions of the triangle
	      {
		noalias(rRightHandSideVector) +=  0.0 * densities[i] * specific_heat[i] *  prod(Cuatros[i],ms_temp_vec_np) / Dt ;
	      }
	      
	    noalias(rRightHandSideVector) +=  0.0 * prod(Mass_matrix,ms_temp_vec_np);
	      
#endif 	  
	   
	    rRightHandSideVector[0] += 0.0 * GetGeometry()[0].FastGetSolutionStepValue(HEAT_FLUX) * 0.25 * Area; 
	    rRightHandSideVector[1] += 0.0 * GetGeometry()[1].FastGetSolutionStepValue(HEAT_FLUX) * 0.25 * Area; 
	    rRightHandSideVector[2] += 0.0 * GetGeometry()[2].FastGetSolutionStepValue(HEAT_FLUX) * 0.25 * Area; 
	    rRightHandSideVector[3] += 0.0 * GetGeometry()[3].FastGetSolutionStepValue(HEAT_FLUX) * 0.25 * Area;

	    double pruebas1 =  0.0;//this->GetValue(heatsource);
	    pruebas1 =heat_flux * Area;

      
	 
	for (unsigned int iii = 0; iii < number_of_points; iii++)
	  ms_temp_vec_np[iii] = GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar);
	noalias(rRightHandSideVector) -= prod(other_matrix, ms_temp_vec_np);
       
      }  //si el elemento estaba cortado!!!!
      else
	{
	  
	  double norm_grad;
      	  double res;
          double k_aux;
	  //int coi=0;
	  double absorptioncoefficient=0.0;
	  
	  int suma=0;	 int sumai=0;
	  int s0, s1, s2, s3;
	  s0=0;s1=0;s2=0; s3=0; 	
	  int elem=0;	
	  int i0, i1,i2,i3;
	  int ii0, ii1,ii2,ii3;
	  //has_negative_node=false;
	  //has_positive_node=false;
          
	  /*
	  if(GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE)==1.0) s0=1; 	//IS_INTERFACE
	  if(GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE)==1.0) s1=1;
	  if(GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE)==1.0) s2=1;
	  if(GetGeometry()[3].FastGetSolutionStepValue(IS_INTERFACE)==1.0) s3=1;
	  */
	  //elem=s0 + s1 + s2 + s3; 

	  int so_heat_flux=0;
          
	  /*if(elem==4) has_negative_node=true;

	  else has_negative_node=false;
	  */
	  if (has_negative_node== true) 
	    {
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


	      density=905.0;           
	      specific_heat=1900.0;
	      conductivity=100.0;//0.12;

	      
	    } 
	  else 
	    {
	      conductivity=0.0131;//0.0131;//0.05;//1.0; 0.05
	      density=1.0;//1000.0;
	      specific_heat=1310.0;//620.0;//1.0;
	      so_heat_flux=1;	
	      absorptioncoefficient=75.0;
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

					  conductivity=0.0131;//0.0131;//0.05;//1.0; 0.05
	      	  density=1.0;//1000.0;
	      	  specific_heat=1310.0;//620.0;//1.0;


	    }
	  
	 
	  
	  so_heat_flux=0.0;
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
	  
	  	  
	  double c1 = 4.00;
	  double c2 = 2.00;
	  double h = sqrt(2.00 * Area);
	  double norm_u = norm_2(ms_vel_gauss);

	  double tau1=0.0;

		#if defined(NOT_STATIONARY)
	  const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
		tau1 = (h * h) / (density * specific_heat * BDFcoeffs[0] * h * h + c1 * conductivity + c2 * density * specific_heat * (norm_u + 1e-6) * h);
	  #endif
	  
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
      
	  
	  noalias(rLeftHandSideMatrix) = (conductivity + 0.0 * so_heat_flux * k_aux * h) * prod(msDN_DX, trans(msDN_DX))* Area;  //0.0 * k_aux * h PARA PAPER RADIACION

		#if defined(NOT_STATIONARY)
	  noalias(rLeftHandSideMatrix) += 1.0/Dt * (density * specific_heat) * msMassFactors * Area;	 
		#endif

	  noalias(ms_u_DN) = prod(msDN_DX, ms_vel_gauss);
	  
	  noalias(rLeftHandSideMatrix) += 0.0 * so_heat_flux * (density * specific_heat) * outer_prod(msN, ms_u_DN)* Area;
	  
	  noalias(rLeftHandSideMatrix) += 0.0 * so_heat_flux * density * specific_heat * tau1 * outer_prod(ms_u_DN, ms_u_DN)* Area;
	  
	  noalias(rLeftHandSideMatrix) += 0.0 * so_heat_flux * 1.0 / Dt * tau1 * msAuxMat * (density *specific_heat)* Area;
	  
 
	  
	//  noalias(rRightHandSideVector) = (heat_flux * 1.0) * msN * Area *1.0 ;// * 0.0;	
	  
	  
	//KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
	//SOLO EN EL AIRE
	//KRATOS_WATCH(elem);
	double StefenBoltzmann = 5.67e-8;
	double emissivity = 1.0;
	double aux = pow(298.0,4);
	//double convection_coefficient=0.0;
	noalias(Aux)=rLeftHandSideMatrix;
	  
	  msAuxVec[0]=GetGeometry()[0].FastGetSolutionStepValue(rUnknownVar,1);
	  msAuxVec[1]=GetGeometry()[1].FastGetSolutionStepValue(rUnknownVar,1);
	  msAuxVec[2]=GetGeometry()[2].FastGetSolutionStepValue(rUnknownVar,1);
	  msAuxVec[3]=GetGeometry()[3].FastGetSolutionStepValue(rUnknownVar,1);
	  
		#if defined(NOT_STATIONARY)
	  noalias(rRightHandSideVector) += so_heat_flux * 1.0/Dt * tau1 * density* specific_heat * prod(msAuxMat, msAuxVec) * Area;
		#endif

	  #if defined(NOT_STATIONARY)
	  for (unsigned int iii = 0; iii < number_of_points; iii++)
	    ms_temp_vec_np[iii] = -1.0/Dt * GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar, 1);
	  
	  noalias(rRightHandSideVector) -= 1.0 * prod(msMassFactors, ms_temp_vec_np * density * specific_heat* Area);
	  #endif
	  //subtracting the dirichlet term
	  // RHS -= LHS*temperatures
    for (unsigned int iii = 0; iii < number_of_points; iii++)
	     ms_temp_vec_np[iii] = GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar) ;
	   noalias(rRightHandSideVector) -= prod(Aux, ms_temp_vec_np);    

	  }
	}
    else
      {
	
	
  }
  KRATOS_CATCH("");
}
  
  //************************************************************************************
  //************************************************************************************

  void ConvDiff3Denriched::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
  }

  //************************************************************************************
  //************************************************************************************
  void ConvDiff3Denriched::Heat_Source(VectorType& rRightHandSideVector, const int ndivisionsp, array_1d<double,6>& volumesp, array_1d<double,6>& conductivitiesp, boost::numeric::ublas::bounded_matrix<double,6, 8 >& Ngaussnewp, const double heatp)
  {
    //double coefficientaux=0.0;
    double heat_flux=0.0;
      
    for(unsigned int zz = 0; zz < 4; zz++) heat_flux += GetGeometry()[zz].FastGetSolutionStepValue(HEAT_FLUX);
    heat_flux *= 0.25;// *0.0;
      
    for (unsigned int i=0; i<ndivisionsp ; i++)
      {
       // coefficientaux =  heat_flux * volumesp[i];// * conductivitiesp[i];
        
        if(this->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE)>0.0) 
        {
			rRightHandSideVector[0] += 0.0 * Ngaussnewp(i,0) * GetGeometry()[0].FastGetSolutionStepValue(HEAT_FLUX) * volumesp[i]; // * conductivities[i];
		}
		
        if(this->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE)>0.0) 
		{
		rRightHandSideVector[1] += 0.0 * Ngaussnewp(i,1) * GetGeometry()[1].FastGetSolutionStepValue(HEAT_FLUX) * volumesp[i]; // * conductivities[i];
		}
		
        if(this->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE)>0.0) 
		{
		rRightHandSideVector[2] += 0.0 * Ngaussnewp(i,2) * GetGeometry()[2].FastGetSolutionStepValue(HEAT_FLUX) * volumesp[i]; // * conductivities[i];
		}	
		
        if(this->GetGeometry()[3].FastGetSolutionStepValue(DISTANCE)>0.0) 
		{
		rRightHandSideVector[3] += 0.0 * Ngaussnewp(i,3) * GetGeometry()[3].FastGetSolutionStepValue(HEAT_FLUX) * volumesp[i]; // * conductivities[i];  
		}
      }
	
  }
//************************************************************************************
//************************************************************************************
void ConvDiff3Denriched::Face_Heat_Flux(ProcessInfo& rCurrentProcessInfo,VectorType& rRightHandSideVector, const int ninterp, const int number_of_nodesp, std::vector< Matrix >& coord_interfacep, std::vector< Matrix >& nodes_auxp, int qcp, Vector& Fluxp,std::vector< Matrix >& area_interfacep, double y_min , double q_x, double q_y, double q_z, const double fhf)
{

KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");	
}

  void ConvDiff3Denriched::qi( const int ndivisionsp, std::vector< Matrix > edges_tauxp, std::vector< Matrix > nodes_auxp, std::vector< Matrix > rGradientauxp,  array_1d<double,6> conductivitiesp, Matrix& Kaux1p,Matrix& Lenrichaux1p)
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


   //KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
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
	       KRATOS_WATCH("area_cero");
	       KRATOS_WATCH("area_cero");
	       KRATOS_WATCH("area_cero");
	       KRATOS_WATCH("area_cero");
	       KRATOS_WATCH(Point1);
	       KRATOS_WATCH(Point2);
	       KRATOS_WATCH(Point3);
	       KRATOS_WATCH(trianglearea);
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
	   //KRATOS_WATCH(area_normal);
	   //KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
	   msAuxVec[0]=c0;
	   msAuxVec[1]=c1;
	   msAuxVec[2]=c2;
	      
	   norm_u = msAuxVec[0]*msAuxVec[0] + msAuxVec[1]*msAuxVec[1] + msAuxVec[2]*msAuxVec[2];
	   norm_c =sqrt(norm_u);
	      
	   if(norm_c<=0.00000000000001)
	     {
		
	       KRATOS_WATCH("norm_c_cero");
	       KRATOS_WATCH("norm_c_cero");
	       KRATOS_WATCH("norm_c_cero");
	       KRATOS_WATCH("norm_c_cero");
	       KRATOS_WATCH(norm_c);
	       KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
	       //sssss
	     }
	      
	   for(unsigned int d=0; d<3;d++){
	     normal(d) = area_normal[d]/ norm_c;
	   }   
	
	
	     
	      
	   int max;
	   if(ndivisionsp==2) max=1;
	   else if(ndivisionsp==3) max=2;
	   else if(ndivisionsp==4) max=3;
	   else max=4;
	   
	   
	   //double z, xx;
	 
	 
	   for(unsigned int zz = 0; zz < 3; zz++)
	     {
	       if(edges_tauxp[i](j,zz)>3)
		 {
		   int position=edges_tauxp[i](j,zz);
		   for (unsigned int t = 0; t < 4; t++){ 
		  
		  
		     for(unsigned int g=0; g<3; g++)
		       Kaux1p(position,t) +=(rGradientauxp[i](t,g) * normal(g)) *trianglearea * 0.33333333333333333333333333333* conductivitiesp[i];
		   }
		   for (unsigned int t = 0; t < max; t++){
		     for(unsigned int g=0; g<3; g++)
		       Lenrichaux1p(position,t+4) +=(rGradientauxp[i](t+4,g) * normal(g) ) *trianglearea * 0.33333333333333333333333333333 * conductivitiesp[i];
		   }
		 }
	     }
	 }
 }
  // this subroutine calculates the nodal contributions for the explicit steps of the
  // fractional step procedure

void ConvDiff3Denriched::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
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
void ConvDiff3Denriched::Face_Heat_Flux_2(ProcessInfo& rCurrentProcessInfo,const int ninterp, const int number_of_nodesp, std::vector< Matrix >& coord_interfacep, std::vector< Matrix >& nodes_auxp, const int qcp, Vector& Fluxp,std::vector< Matrix >& area_interfacep )
    {

      //Vector Face_Heat_Flux;
      //Face_Heat_Flux.resize(number_of_nodes);
      //Face_Heat_Flux=ZeroVector(number_of_nodes);
      int a, b,c;
      array_1d<double,3>  Point1; 
      array_1d<double,3>  Point2; 
      array_1d<double,3>  Point3;
      double StefenBoltzmann = 5.67e-8;
      double emissivity = 1.0;
      //double temperature=0.0;
      double aux = pow(298.0,4);
      //for(unsigned int j = 0; j < 4; j++) temperature += GetGeometry()[j].FastGetSolutionStepValue(TEMPERATURE_OLDIT);
      //temperature *=0.25;
      //int iteration_number = 0;//rCurrentProcessInfo[ITERATION_NUMBER];
      
      Vector rEnrichedTemperature_oldit;
      rEnrichedTemperature_oldit.resize(4);
  /*    rEnrichedTemperature_oldit =  this->GetValue(ENRICHED_TEMPERATURE_OLDIT);
*/
		rEnrichedTemperature_oldit ;

      double t1, t2,t3;
      for(unsigned int m = 0; m < ninterp; m++)  
	{
	  a=coord_interfacep[m](0,0);
	  b=coord_interfacep[m](0,1);
	  c=coord_interfacep[m](0,2);
	  //if(iteration_number==1){
	  // t1= temperature;
	  // t2= temperature;
	  // t3= temperature;
	  // }
	  // else
	  // {
	  //double sum=0.3333333*( rEnrichedTemperature_oldit(a-4) + rEnrichedTemperature_oldit(b-4) +rEnrichedTemperature_oldit(c-4) );
	  t1= rEnrichedTemperature_oldit(a-4);
	  t2= rEnrichedTemperature_oldit(b-4);
	  t3= rEnrichedTemperature_oldit(c-4);
	  //}	  

	  for(unsigned int d=0; d<3;d++){
	    Point1(d)=nodes_auxp[a](0,d);
	    Point2(d)=nodes_auxp[b](0,d);
	    Point3(d)=nodes_auxp[c](0,d);
	  }
	  double trianglearea=CalculateTriangleArea3D( Point1, Point2, Point3);
	  
	  double y=0.33333*(Point1(1)+Point2(1)+Point3(1));
	  double qq=0.0;
	  if(y>-0.063) qq=100000.0 * -y/0.063;//pow(11,(-y))/(pow(11,0.063));
	  else qq=100000.0;
	  qq=100000.0;
	  qq=15000.0;
	  
	  Fluxp(a-4) += qq * 0.33333333333333333333333333333 * trianglearea - 1.0 * emissivity*StefenBoltzmann*(pow(t1,4) - aux) * 0.33333333333333333333333333333 * trianglearea;
	  Fluxp(b-4) += qq * 0.33333333333333333333333333333 * trianglearea - 1.0 * emissivity*StefenBoltzmann*(pow(t2,4) - aux) * 0.33333333333333333333333333333 * trianglearea; 
	  Fluxp(c-4) += qq * 0.33333333333333333333333333333 * trianglearea - 1.0 * emissivity*StefenBoltzmann*(pow(t3,4) - aux) * 0.33333333333333333333333333333 * trianglearea;
	  area_interfacep[m](0,0)=trianglearea;
	}
      /*
	if(number_of_nodesp==1)
	{
	array_1d<int, 2 > nodes; nodes(0)=-1; nodes(1)=-1;
	a=coord_interfacep[0](0,0);
	int pos=0;
	for(unsigned int iii = 0; iii<4; iii++)
	{
	if(this->GetGeometry()[iii].FastGetSolutionStepValue(DISTANCE)==0)
	{
	nodes(pos)=iii;
	pos++;
	}
	}
	
	Point1(0)=nodes_auxp[4](0,0);
	Point1(1)=nodes_auxp[4](0,1);
	Point1(2)=nodes_auxp[4](0,2);
	Point2(0)=nodes_auxp[nodes(0)](0,0);
	Point2(1)=nodes_auxp[nodes(0)](0,1);
	Point2(2)=nodes_auxp[nodes(0)](0,2);
	Point3(0)=nodes_auxp[nodes(1)](0,0);
	Point3(1)=nodes_auxp[nodes(1)](0,1);
	Point3(2)=nodes_auxp[nodes(1)](0,2);
	  
	double trianglearea=CalculateTriangleArea3D( Point1, Point2, Point3);
	//KRATOS_WATCH(  trianglearea);
	Fluxp(0) += qcp * 0.33333333333333333333333333333 * trianglearea; 
	 
	rRightHandSideVector[nodes(0)] += qcp * 0.33333333333333333333333333333 * trianglearea - emissivity*StefenBoltzmann*(pow(temperature,4) - aux) * 0.33333333333333333333333333333 * trianglearea; 
	rRightHandSideVector[nodes(1)] += qcp * 0.33333333333333333333333333333 * trianglearea - emissivity*StefenBoltzmann*(pow(temperature,4) - aux) * 0.33333333333333333333333333333 * trianglearea; 
	KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");  
	}
      */
      /*if(number_of_nodesp==2)
	{
	array_1d<int, 2 > nodes; nodes(0)=-1; nodes(1)=-1; 
	 
	int pos=0;
	for(unsigned int iii = 0; iii<4; iii++)
	{
	if(this->GetGeometry()[iii].FastGetSolutionStepValue(DISTANCE)==0)
	{
	nodes(pos)=iii;
	pos++;
	}
	}
	Point1(0)=nodes_auxp[4](0,0);
	Point1(1)=nodes_auxp[4](0,1);
	Point1(2)=nodes_auxp[4](0,2);
	Point2(0)=nodes_auxp[5](0,0);
	Point2(1)=nodes_auxp[5](0,0);
	Point2(2)=nodes_auxp[5](0,0);
	Point3(0)=nodes_auxp[nodes(0)](0,0);
	Point3(1)=nodes_auxp[nodes(0)](0,1);
	Point3(2)=nodes_auxp[nodes(0)](0,2);
	double trianglearea=CalculateTriangleArea3D( Point1, Point2, Point3);
	Fluxp(0) += qcp * 0.33333333333333333333333333333 * trianglearea - emissivity*StefenBoltzmann*(pow(temperature,4) - aux) * 0.33333333333333333333333333333 * trianglearea; 
	Fluxp(1) += qcp * 0.33333333333333333333333333333 * trianglearea - emissivity*StefenBoltzmann*(pow(temperature,4) - aux) * 0.33333333333333333333333333333 * trianglearea; 
	   
	rRightHandSideVector[nodes(0)] += qcp * 0.33333333333333333333333333333 * trianglearea - emissivity*StefenBoltzmann*(pow(temperature,4) - aux) * 0.33333333333333333333333333333 * trianglearea; 
	KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
	}
      */
 }
//************************************************************************************
//************************************************************************************

  void ConvDiff3Denriched::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
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

  void ConvDiff3Denriched::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
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
  double ConvDiff3Denriched::CalculateVol(const double x0, const double y0, const double z0,const double x1, const double y1, const double z1,const double x2, const double y2, const double z2,
				  const double x3, const double y3, const double z3 )
  {
    double x10 = x1 - x0;
    double y10 = y1 - y0;
    double z10 = z1 - z0;

    double x20 = x2 - x0;
    double y20 = y2 - y0;
    double z20 = z2 - z0;

    double x30 = x3 - x0;
    double y30 = y3 - y0;
    double z30 = z3 - z0;

    double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;
    return detJ * 0.1666666666666666666667;
   
  }
//************************************************************************************
//************************************************************************************
  void ConvDiff3Denriched::CalculateXG(boost::numeric::ublas::bounded_matrix<double,4, 3 >& rPoints, array_1d<double, 3 > & x_g){
 

    //intersection Points - at most we shall consider 4 points
    std::vector<array_1d<double,3> > IntersectionPoints;
    //std::vector<array_1d<double,3> > IntersectionVel;

    IntersectionPoints.reserve(8);
    //IntersectionVel.reserve(4);

    array_1d<double,3> Point;
    array_1d<double,3> xg;
    //identify
    int intersection_count=0;

    //iterating over edges (01 02 03 12 13 23)
    for (int i=0; i<3; i++)
      {
	for (int j=i+1; j<4; j++)
	  {
	    //std::cout<<"edge ij "<<i<<j<<std::endl;
	    double d0=this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
	    double d1=this->GetGeometry()[j].FastGetSolutionStepValue(DISTANCE); 

	    //double product=d0*d1;
	    //KRATOS_WATCH(product)

	    //if the product of distances of two nodes is negative - the edge is crossed
	    if (d0*d1<0.0)
	      {
		//std::cout<<"Intersected edge "<<i<<j<<std::endl;
		//KRATOS_WATCH(product)
		double x0=rPoints(i,0);//this->GetGeometry()[i].X();
		double x1=rPoints(j,0); //this->GetGeometry()[j].X();
		double k=0.0;
		double b=0.0;
		double x_inters=0.0;
		double y_inters=0.0;
		double z_inters=0.0;
		if (x1!=x0)
		  {
		    k=(d1-d0)/(x1-x0);
		    b=d0-k*x0;
		    x_inters=-b/k;
		  }
		else
		  x_inters=x0;

		double y0=rPoints(i,1); //this->GetGeometry()[i].Y();
		double y1=rPoints(j,1); //this->GetGeometry()[j].Y();
		if (y1!=y0)
		  {
		    k=(d1-d0)/(y1-y0);
		    b=d0-k*y0;
		    y_inters=-b/k;
		  }
		else
		  y_inters=y0;

		double z0= rPoints(i,2); //this->GetGeometry()[j].Z();
		double z1= rPoints(j,2); //this->GetGeometry()[j].Z();
		if (z1!=z0)
		  {
		    k=(d1-d0)/(z1-z0);
		    b=d0-k*z0;
		    z_inters=-b/k;
		  }
		else
		  z_inters=z0;

		Point[0]=x_inters;
		Point[1]=y_inters;
		Point[2]=z_inters;

		IntersectionPoints.push_back(Point);


		intersection_count++; 			
				
	      }

	  }
      }
    //KRATOS_WATCH(intersection_count)

    IntersectionPoints.resize(intersection_count);

	    
    //if the element is intersected by the embedded skin, create condition
    //if (intersection_count!=0)
    if (intersection_count==4 )
      {



	xg[0]= 0.25*(IntersectionPoints[0][0] + IntersectionPoints[1][0] +IntersectionPoints[2][0]+IntersectionPoints[3][0] );
	xg[1]= 0.25*(IntersectionPoints[0][1] + IntersectionPoints[1][1] +IntersectionPoints[2][1]+IntersectionPoints[3][1] );
	xg[2]= 0.25*(IntersectionPoints[0][2] + IntersectionPoints[1][2] +IntersectionPoints[2][2]+IntersectionPoints[3][2] );


	KRATOS_WATCH("PUNTOSSSSSSSSSS DE INTERSEPCION");
	KRATOS_WATCH("PUNTOSSSSSSSSSS DE INTERSEPCION");
	KRATOS_WATCH(xg);
	KRATOS_WATCH("PUNTOSSSSSSSSSS DE INTERSEPCION");
	KRATOS_WATCH("PUNTOSSSSSSSSSS DE INTERSEPCION");

      }

    if (intersection_count==3 )
      {


	xg[0]= 0.333*(IntersectionPoints[0][0] + IntersectionPoints[1][0] + IntersectionPoints[2][0]);
	xg[1]= 0.333*(IntersectionPoints[0][1] + IntersectionPoints[1][1] + IntersectionPoints[2][1]);
	xg[2]= 0.333*(IntersectionPoints[0][2] + IntersectionPoints[1][2] + IntersectionPoints[2][2]);
	     
	KRATOS_WATCH("PUNTOSSSSSSSSSS DE INTERSEPCION");
	KRATOS_WATCH("PUNTOSSSSSSSSSS DE INTERSEPCION");
	KRATOS_WATCH(xg);
	KRATOS_WATCH("PUNTOSSSSSSSSSS DE INTERSEPCION");
	KRATOS_WATCH("PUNTOSSSSSSSSSS DE INTERSEPCION");



      }
    if (intersection_count==2 )
      {
	KRATOS_WATCH("PUNTOSSSSSSSSSS DE INTERSEPCION");
	KRATOS_WATCH("PUNTOSSSSSSSSSS DE INTERSEPCION");
	KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
      }
    if (intersection_count==1 )
      {
	KRATOS_WATCH("PUNTOSSSSSSSSSS DE INTERSEPCION");
	KRATOS_WATCH("PUNTOSSSSSSSSSS DE INTERSEPCION");
	KRATOS_THROW_ERROR(std::logic_error, "method not implemented", "");
      }
  }

  double ConvDiff3Denriched::CalculateTriangleArea3D(	array_1d<double,3>& Point1, array_1d<double,3>& Point2, array_1d<double,3>& Point3	)
  {
    //Heron's formula
    double a=Length(Point1, Point2);//sqrt((Point1[0]-Point2[0])*(Point1[0]-Point2[0]) + (Point1[1]-Point2[1])*(Point1[1]-Point2[1]) +(Point1[2]-Point2[2])*(Point1[2]-Point2[2]));
    double b=Length(Point1, Point3);//sqrt((Point3[0]-Point2[0])*(Point3[0]-Point2[0]) + (Point3[1]-Point2[1])*(Point3[1]-Point2[1]) +(Point3[2]-Point2[2])*(Point3[2]-Point2[2]));
    double c=Length(Point2, Point3);//sqrt((Point1[0]-Point3[0])*(Point1[0]-Point3[0]) + (Point1[1]-Point3[1])*(Point1[1]-Point3[1]) +(Point1[2]-Point3[2])*(Point1[2]-Point3[2]));
    double p=0.5*(a+b+c);
    return sqrt(p*(p-a)*(p-b)*(p-c));
  }

  double ConvDiff3Denriched::Length(array_1d<double,3>& Point1, array_1d<double,3>& Point2)
  {
    //KRATOS_WATCH("length calculation")
    return sqrt((Point1[0]-Point2[0])*(Point1[0]-Point2[0]) + (Point1[1]-Point2[1])*(Point1[1]-Point2[1]) +(Point1[2]-Point2[2])*(Point1[2]-Point2[2]));
  }

  /*double ConvDiff3D::ComputeSmagorinskyViscosity(const boost::numeric::ublas::bounded_matrix<double, 4, 3 > & DN_DX, const double& h, const double& C, const double nu )
    {
    boost::numeric::ublas::bounded_matrix<double, 3, 3 > dv_dx = ZeroMatrix(3, 3);
  
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

  //************************************************************************************
  //************************************************************************************
  template<class T>
  
  bool ConvDiff3Denriched::InvertMatrix(const T& input, T& inverse)
    
  {
    
    typedef permutation_matrix<std::size_t> pmatrix;
    
    // create a working copy of the input
    
    T A(input);
    
    // create a permutation matrix for the LU-factorization
    
    pmatrix pm(A.size1());

    // perform LU-factorization
    
    int res = lu_factorize(A, pm);
    //KRATOS_WATCH(res);
    if (res != 0)
      
      return false;
    
    // create identity matrix of "inverse"
    
    inverse.assign(identity_matrix<double> (A.size1()));

    // backsubstitute to get the inverse
    
    lu_substitute(A, pm, inverse);
    
    return true;
    
  }
  

} // Namespace Kratos


