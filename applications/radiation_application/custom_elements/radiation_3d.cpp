//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Julio Marti
//


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
    array_1d<double,3> aux_var= ZeroVector(3); //dimesion coincides with space dimension
    boost::numeric::ublas::bounded_matrix<double,4,1> msShapeFunc = ZeroMatrix(4,1);
    array_1d<double,4> msAuxVec; // = ZeroVector(4); //dimension = number of nodes
    boost::numeric::ublas::bounded_matrix<double,4,4> msAuxMat = ZeroMatrix(4,4);
    boost::numeric::ublas::bounded_matrix<double,4,4> msAux = ZeroMatrix(4,4);
    
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
    
    const Variable<double>& rUnknownVar= my_settings->GetUnknownVariable();
    
    
    const double StefenBoltzmann = 5.67e-8;
    double absorptioncoefficient = 100.0;
    
    
    boost::numeric::ublas::bounded_matrix<double,4, 3 > coords;
    array_1d<double,4> distances;
    array_1d<double,6>  volumes(6);
    boost::numeric::ublas::bounded_matrix<double,6, 4 > Ngauss;
    array_1d<double,6>  signs(6);
    std::vector< Matrix > gauss_gradients(6);
    
    msMassFactors(0,0) = 1.00/4.00; msMassFactors(0,1) = 0.00;		msMassFactors(0,2) = 0.00;       msMassFactors(0,3) = 0.00;
    msMassFactors(1,0) = 0.00;		msMassFactors(1,1) = 1.00/4.00; msMassFactors(1,2) = 0.00;       msMassFactors(1,3) = 0.00;
    msMassFactors(2,0) = 0.00;		msMassFactors(2,1) = 0.00;		msMassFactors(2,2) = 1.00/4.00;  msMassFactors(2,3) = 0.00;
    msMassFactors(3,0) = 0.00;		msMassFactors(3,1) = 0.00;		msMassFactors(3,2) = 0.00;       msMassFactors(3,3) = 1.00/4.00;
    
    double 	conductivity=1.0/(3.0*absorptioncoefficient);
    
    
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
    noalias(rLeftHandSideMatrix) += absorptioncoefficient * msMassFactors * Area; 
    
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
    double constant;
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
	
	msAux=ZeroMatrix(4,4);
	msAux(0,0)=N[0]*N[0]; //+ N[0]*N[1] + N[0]*N[2] + N[0]*N[3];
	msAux(0,1)+=N[0]*N[1];
	msAux(0,2)+=N[0]*N[2];
	msAux(0,3)+=N[0]*N[3];
	    
	msAux(1,1)=N[1]*N[0];
	msAux(1,1)+=N[1]*N[1];// + N[1]*N[0] + N[1]*N[2] + N[1]*N[3];
	msAux(1,1)+=N[1]*N[2];
	msAux(1,1)+=N[1]*N[3];
	    
	msAux(2,2)=N[2]*N[0];
	msAux(2,2)+=N[2]*N[1];
	msAux(2,2)+=N[2]*N[2];// + N[2]*N[0] + N[2]*N[1] + N[2]*N[3];
	msAux(2,2)+=N[2]*N[3];
	    
	msAux(3,3)=N[3]*N[0];
	msAux(3,3)+=N[3]*N[1];
	msAux(3,3)+=N[3]*N[2];
	msAux(3,3)+=N[3]*N[3];// + N[3]*N[0] + N[3]*N[1] + N[3]*N[2];
	    
	    
	//	    noalias(rLeftHandSideMatrix) += 1.0 * Weight * absorptioncoefficient * msAux; //LO COMENTO AHORA
	noalias(rRightHandSideVector) += 4.0 * absorptioncoefficient * StefenBoltzmann * prod(msAux,msAuxVec)  * Weight;
	    
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
  void Rad3D::qi( const unsigned int ndivisionsp, std::vector< Matrix > edges_tauxp, std::vector< Matrix > nodes_auxp, std::vector< Matrix > rGradientauxp,  array_1d<double,6> conductivitiesp, Matrix& Kaux1p,Matrix& Lenrichaux1p)
  {
    array_1d<double,3> v3;
    array_1d<double,3> v4;
    array_1d<double,3> An;// =zero;
    //const unsigned int number_of_points = GetGeometry().size();
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


