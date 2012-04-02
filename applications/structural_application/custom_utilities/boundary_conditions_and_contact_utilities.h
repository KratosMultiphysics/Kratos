/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
//   Last Modified by:    $Author: Nelson Lafontaine $
//   Date:                $Date: 2007-03-06 10:30:34 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_BOUNDARY_CONDITIONS_AND_CONTACT_UTILITIES_INCLUDED )
#define  KRATOS_BOUNDARY_CONDITIONS_AND_CONTACT_UTILITIES_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <iomanip> 
#include <fstream> 
#include <algorithm>
#include <set>
#include <time.h>



#ifdef _OPENMP
#include <omp.h>
#endif

// External includes


// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/mesh.h"

#include "geometries/geometry.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"
#include "geometries/line_2d_2.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"


#include "spatial_containers/spatial_containers.h"
#include "spatial_containers/bounding_box.h"
#include "spatial_containers/cell.h"
#include "spatial_containers/bins_dynamic_objects.h"
#include "spatial_containers/bins_dynamic.h"
#include "spatial_containers/bins_static_objects.h"

#include "utilities/spatial_containers_configure.h"
#include "utilities/geometry_utilities.h"
#include "utilities/timer.h"
// #include "utilities/timer_CLabra.h"
#include "custom_conditions/slave_contact_point_2d.h"
#include "custom_conditions/slave_contact_point_3d.h"
#include "custom_conditions/master_contact_point_2d.h"
#include "custom_conditions/master_contact_face_2d.h"
#include "custom_conditions/master_contact_face_3D.h"
#include "custom_conditions/point_segment_contact_link.h"
#include "custom_conditions/point_point_contact_link.h"
#include "custom_conditions/contact_link_3D_explicit.h"


#include "geometries/plane.h"
#include "custom_utilities/segment_2d.h"
#include "custom_utilities/intersect_triangles_cases.h"
#include "processes/find_nodal_neighbours_process.h"
#include "processes/find_elements_neighbours_process.h"
#include "processes/find_conditions_neighbours_process.h"



namespace Kratos
{
  class BoundaryConditionsAndContactUtilities
	{
	public:
	   	   
	    #define EPSILON 1.0e-10
	    #define NEPSILON -1.0e-10
	    #define BEPSILON 1.0e+15
	    
	    static const int IT_POINT   = 0; 
	    static const int IT_SEGMENT = 1;
	    static const int IT_EMPTY   = 2;
	    
	    enum Exist_Node     {no_nodes   = 0, yes_nodes};
	    enum Near_Node      {no_near    = 0, yes_near};
	    enum Object         {is_node    = 0, is_object};

	    

	    KRATOS_CLASS_POINTER_DEFINITION(BoundaryConditionsAndContactUtilities);

	    /// Utilities
	    typedef IntersectionSegment2DToSegment2D  IntersectionSegments;

     
	    /// Elements
	    typedef ModelPart::ElementsContainerType                     ElementsArrayType;
	    /*
	    typedef ModelPart::ElementsContainerType::ContainerType      ContainerType; 
	    typedef ContainerType::value_type                            PointerType;
	    typedef ContainerType::iterator                              IteratorType; 
 	    typedef std::vector<PointerType>::iterator                   PointerTypeIterator;
 	    typedef ContactPair<PointerType>                             ContactPairType; 
 	    typedef std::vector<ContactPairType>                         ContainerContactPair;   
 	    typedef ContainerContactPair::iterator                       IteratorContainerContactPair;  
 	    typedef ContainerContactPair::value_type                     PointerContainerContactPair;  
	    */
	
	    /// Conditions General
	    typedef ModelPart::ConditionsContainerType                   ConditionsArrayType;
	    typedef ModelPart::ConditionsContainerType::ContainerType    ConditionsContainerType; 
	    typedef ConditionsContainerType::iterator                    ConditionsIteratorType; 
	    typedef ConditionsContainerType::value_type                  ConditionsPointerType;
	    typedef ContactPair<ConditionsPointerType>                   ConditionsContactPairType; 
	    typedef std::vector<ConditionsPointerType>::iterator         ConditionsPointerTypeIterator;
	    typedef std::vector<ConditionsContactPairType>               ConditionsContainerContactPair;   
	    typedef ConditionsContainerContactPair::iterator             ConditionsIteratorContainerContactPair;  
	    typedef ConditionsContainerContactPair::value_type           ConditionsPointerContainerContactPair;  
	    
	    /// Condition Especificas
	    typedef SlaveContactPoint2D                                  SlaveContactPointType; 
	    typedef MasterContactPoint2D                                 MasterContactPointType; 
	    typedef MasterContactFace2D                                  MasterContactFaceType;

	    
	    ///Nodes and properties
	    typedef Node<3>                                              NodeType;
	    typedef Node<3>::Pointer                                     NodePointerType;
	    typedef Geometry<NodeType>                                   GeometryType;
            typedef GeometryType::PointsArrayType                        PointsArrayType;
	    typedef ModelPart::NodesContainerType                        NodesArrayType;
	    typedef Element::GeometryType                                GeomType; 
	    typedef ModelPart::NodesContainerType::ContainerType         NodesContainerType;
	    typedef NodesContainerType::iterator                         NodesIteratorType;  
	    typedef Properties                                           PropertiesType;
	    
	    
	     
	    static const std::size_t space_dim                = 2;
	    typedef SpatialContainersConfigure<space_dim>     Configure;  
	    typedef Configure::PointType                      PointType; 
	    typedef PointType::CoordinatesArrayType           CoordinatesArrayType;
            typedef Configure::ContainerType                  ContainerType;   
            typedef Configure::PointerType                    PointerType;
            typedef Configure::IteratorType                   IteratorType; 
            typedef Configure::ResultContainerType            ResultContainerType;
	    typedef Configure::ResultPointerType              ResultPointerType;
            typedef Configure::ResultIteratorType             ResultIteratorType; 
            typedef Configure::ContactPairType                ContactPairType;
            typedef Configure::ContainerContactType           ContainerContactType; 
            typedef Configure::IteratorContactType            IteratorContactType; 
            typedef Configure::PointerContactType             PointerContactType; 
            typedef Configure::PointerTypeIterator            PointerTypeIterator;
	    typedef ContainerContactType                      ContainerContactPair;   
	    typedef IteratorContactType                       IteratorContainerContactPair;  
	    typedef PointerContactType                        PointerContainerContactPair;  
   
	    
	    
            BoundaryConditionsAndContactUtilities(){}
            BoundaryConditionsAndContactUtilities(ModelPart& model_part, const unsigned int& dimension, const double& penalty_factor) : mr_model_part(model_part), mrdimension(dimension) 
              {  
		mpenalty_factor = penalty_factor;
		mcompute_boundary_contour = true;
              }   
  
            virtual ~BoundaryConditionsAndContactUtilities(){}




//************************************************************************************
//************************************************************************************   
      // Crea las conciones de contacto valid for lagrage multiplier y setea los elementos
      // que son parte del contorno
      void CreateBoundaries(const unsigned int& initial_conditions_size)
      {
	KRATOS_TRY 
          Clear(initial_conditions_size);     
	  if(mcompute_boundary_contour){
	    std::cout<<"CREATING MASTER SURFACES"<< std::endl; 
	    if(mrdimension==2)
	       CalculateBoundaryContour2D(mMasterConditionsArray);  
	    else
	       CalculateBoundaryContour3D(mMasterConditionsArray);  
	  mcompute_boundary_contour = false;
	  }
	  return;
	KRATOS_CATCH("")
      }



//************************************************************************************
//************************************************************************************            


/// this function use the potencial contact force concept
void ComputeContactForce()
  {

    KRATOS_TRY
    
    if(mrdimension!= int(space_dim))
      KRATOS_ERROR(std::logic_error,  "The Dimension of Configure and ModelPart not iquals  "  , "");  
      
    IteratorType it_begin     =   mBoundaryElements.begin();
    IteratorType it_end       =   mBoundaryElements.end();  
    //BinsObjectDynamic<Configure>  rBinsObjectDynamic(it_begin, it_end);
    //BinsObjectDynamic<Configure>* rBins = &rBinsObjectDynamic;  

    ///Bins estatico
    BinsObjectStatic<Configure>  Bins(it_begin, it_end);
    BinsObjectStatic<Configure>* rBins    = &Bins; 
    const std::size_t MaxNumberOfResults  = 1000;
    std::size_t  NumberOfResults          = 0;
    ResultIteratorType  begin;
    
    #ifdef _OPENMP
    int number_of_threads = omp_get_max_threads();
    #else
    int number_of_threads = 1;
    #endif
 
    vector<unsigned int> partition;
    CreatePartition(number_of_threads, mBoundaryElements.size(), partition);
    //ContactPairType it_pair;
    ResultContainerType  Result(MaxNumberOfResults);
    
    
    std::cout<<"        PARTITION COMPUTING CONTACT CONDITIONS  = " << number_of_threads << std::endl; 
    if(mrdimension==2)
    {
    #pragma omp parallel for firstprivate(NumberOfResults,Result)  private(begin)  
    for(int k=0; k<number_of_threads; k++)          
    {
    IteratorType it_begin = mBoundaryElements.begin() + partition[k];
    IteratorType it_end   = mBoundaryElements.begin() + partition[k+1];
    for(IteratorType it =it_begin; it!=it_end; ++it)
      { 
	begin = Result.begin();
	NumberOfResults = rBins->SearchObjects(*it, begin, MaxNumberOfResults);         
	if(NumberOfResults!=0){
	    for(ResultIteratorType rthis = Result.begin();  rthis!= Result.begin() + NumberOfResults; rthis++){
	      if((*rthis)->GetValue(IS_TARGET)== false && (*it)->Id()!=(*rthis)->Id() && FiltratePairContacts(*it, *rthis)==true){
	            ComputeContactForce2D(*it, *rthis); 
		    ComputeContactDampingForces(*it, *rthis);
	       } } 
            (*it)->GetValue(IS_TARGET) = true; 
        } } } } 
        
    else
    { 
     #ifdef _OPENMP
     double start_prod = omp_get_wtime();
     #endif
     
     #pragma omp parallel for firstprivate(NumberOfResults,Result)  private(begin) 
     for(int k=0; k<number_of_threads; k++)          
     {
     IteratorType it_begin = mBoundaryElements.begin() + partition[k];
     IteratorType it_end   = mBoundaryElements.begin() + partition[k+1];
     for(IteratorType it =it_begin; it!=it_end; ++it)
      { 
	//Result[k].clear();
	//Result[k].reserve(100);
	//rBinsObjectDynamic.SearchObjectsInner(*it, Result[k]);
	begin = Result.begin();
	NumberOfResults = rBins->SearchObjects(*it, begin, MaxNumberOfResults);
	//if(Result[k].size()!=0){
	if(NumberOfResults!=0)
	   for(ResultIteratorType rthis = Result.begin(); rthis!= Result.begin() + NumberOfResults /*rthis!=Result[k].end()*/ ; rthis++){
	     if((*rthis)->GetValue(IS_TARGET)== false && (*it)->Id()!=(*rthis)->Id() && FiltratePairContacts(*it, *rthis)==true){
	         ComputeContactForce3D(*it, *rthis);
	  } }
	  (*it)->GetValue(IS_TARGET) = true; 
      } }
     
      #ifdef _OPENMP
      double stop_prod = omp_get_wtime();
      std::cout <<"          Time Calculating Forces Contact = " << stop_prod - start_prod << std::endl;
      #endif
    }     
	
    std::cout<<"        FINISHING COMPUTE CONTACT CONDITIONS  " << std::endl; 
    KRATOS_CATCH("")
 }


void ComputeContactDampingForces(const PointerType& Target, const PointerType& Contactor)
{
 KRATOS_TRY
 
 typedef Element::GeometryType::Pointer    GeometryPointer; 
 double dampT                              = (Target->GetProperties()[DAMPING_RATIO]);
 double dampC                              = (Contactor->GetProperties()[DAMPING_RATIO]);
 double damp                               = std::max(dampT, dampC);
 const GeometryPointer& GTarget            = (Target)->pGetGeometry();
 const GeometryPointer& GContactor         = (Contactor)->pGetGeometry();
 double  pen_vec_tar                       = mpenalty_factor * (Target->GetProperties()[YOUNG_MODULUS]);
 double  pen_vec_con                       = mpenalty_factor * (Contactor->GetProperties()[YOUNG_MODULUS]);
 double pen                                = std::min(pen_vec_tar, pen_vec_con); 
 double massC                              = 0.00;
 double massT                              = 0.00;
 
 array_1d<double ,3> vel_T = ZeroVector(3);
 array_1d<double ,3> vel_C = ZeroVector(3);
 array_1d<double ,3> vel   = ZeroVector(3);
 
 for(unsigned int i = 0; i<(*GTarget).size(); i++){
      massT+= (*GTarget)(i)->FastGetSolutionStepValue(NODAL_MASS);
      noalias(vel_T) += (*GTarget)(i)->FastGetSolutionStepValue(VELOCITY); 
 }
  for(unsigned int i = 0; i<(*GContactor).size(); i++){
      massC+=  (*GContactor)(i)->FastGetSolutionStepValue(NODAL_MASS);
      noalias(vel_C) += (*GContactor)(i)->FastGetSolutionStepValue(VELOCITY); 
  }
  
 
 double fact                  =  (massT * massC * pen)/(massT + massC);    
 double Ccr                   =  2.00 * std::sqrt(fact);
 array_1d<double ,3> normal_T = ZeroVector(3);
 array_1d<double ,3> normal_C = ZeroVector(3);
 array_1d<double ,3> Center_T = (GTarget)->Center();
 array_1d<double ,3> Center_C = (GContactor)->Center();
 noalias(normal_T)            = Center_C - Center_T;
 noalias(normal_C)            = Center_T - Center_C;
 noalias(normal_T)            = (1.00/norm_2(normal_T)) * normal_T; 
 noalias(vel)                 = vel_T - vel_C;  
 
 double vrn                   = inner_prod(vel, normal_T); 
 double fnd                   = std::fabs(damp * Ccr * vrn);
  
  for(unsigned int i = 0; i<(*GTarget).size(); i++){
      array_1d<double, 3>& rhs = (*GTarget)(i)->FastGetSolutionStepValue(RHS); 
      noalias(rhs)  = 0.33333333333333 * fnd * normal_T;
  } 
 
  for(unsigned int i = 0; i<(*GContactor).size(); i++){
      array_1d<double, 3>& rhs = (*GContactor)(i)->FastGetSolutionStepValue(RHS); 
      noalias(rhs) = 0.33333333333333 * fnd * normal_C; 
  }
    
 KRATOS_CATCH("")
}


/* Triangle to Triangle */
//Compute the normal contact force 
void ComputeContactForce2D(const PointerType& Target, const PointerType& Contactor)
{
  
  KRATOS_TRY
  
  const double R0   = 0.00;
  const double R1   = 1.00;
  const double R2   = 2.00;
  const double RP5  = 0.50;
  const double RP15 = 1.50;
  
  int icontact = 0; 
  int np = 0;
  double a0,a1,a2,b0,b1,b2,c0,c1,c2,n0,n1,n2,fn,fna,fnb;
  double pen,tmp,dmin2,smin,smax;
  double small =  EPSILON;
  double nsmall= -EPSILON;
  double big   =  BEPSILON;
  array_1d<double,10> p;
  array_1d<double,10> s;
  array_1d<double,3>  fx;
  array_1d<double,3>  fy;
  array_1d<double,2>  vol;
  array_1d<array_1d<double,3>,2> rx;
  array_1d<array_1d<double,3>,2> ry;
  array_1d<array_1d<double,3>,2> nx;
  array_1d<array_1d<double,3>,2> ny;
  array_1d<array_1d<array_1d<double,3>,3>,2> d;
 
  double  vol2        = 0.00;
  double  pen_vec_tar = mpenalty_factor * (Target->GetProperties()[YOUNG_MODULUS]);
  double  pen_vec_con = mpenalty_factor * (Contactor->GetProperties()[YOUNG_MODULUS]);
  
  pen = std::min(pen_vec_tar, pen_vec_con); ///penalty term 
  //Element::GeometryType& Tgeom = Target->GetGeometry();
  //Element::GeometryType& Cgeom = Contactor->GetGeometry();
  std::vector<Element::GeometryType::Pointer> Geom(2);
  Geom[0] = Contactor->pGetGeometry();
  Geom[1] = Target->pGetGeometry();
  
      for(unsigned int i=0;i<2;i++){ 
        for(unsigned int j=0;j<3; j++){
	  rx[i][j] = ((*Geom[i])(j))->X();  
	  ry[i][j] = ((*Geom[i])(j))->Y();
      } }
      
      for(unsigned int i=0; i<2; i++){ 
	vol[i]=(rx[i][1]-rx[i][0])*(ry[i][2]-ry[i][0])- (ry[i][1]-ry[i][0])*(rx[i][2]-rx[i][0]);
        //normales salientes  no unitarias  de las aristas de los elementos 
        unsigned int k = 0;
        for(unsigned int j=0; j<3; j++){
	  k= j+1; if(k>2) k=0;
	  nx[i][j]=ry[i][k]-ry[i][j];
	  ny[i][j]=rx[i][j]-rx[i][k];
      } }
      
      //computing the tranformation of nodal coordinate to the local coordinate
      for(unsigned int i=0; i<2; i++){
        unsigned int j = i+1; if(j>1) j=0;
	for(unsigned int k=0; k<3; k++){
	   for(unsigned int l=0; l<3; l++){
	     d[i][k][l]=((rx[j][l]-rx[i][k])*nx[j][l]+ (ry[j][l]-ry[i][k])*ny[j][l])/vol[j];
      } } }
      
      
      dmin2=big;

      /* main loop */
      for(unsigned int it=0; it<2; it++)
      { 
	unsigned int jt=it+1; if(jt>1)jt=0;
	noalias(fx) = ZeroVector(3);
	noalias(fy) = ZeroVector(3); 
	vol2        = vol[jt]*vol[jt];
	n0          = (nx[jt][0]*nx[jt][0]+ny[jt][0]*ny[jt][0])/(vol2);
	n1          = (nx[jt][1]*nx[jt][1]+ny[jt][1]*ny[jt][1])/(vol2);
	n2          = (nx[jt][2]*nx[jt][2]+ny[jt][2]*ny[jt][2])/(vol2);
	for(unsigned int in=0;in<3;in++)
	{ unsigned int jn=in+1; if(jn>2)jn=0;
	  a0=d[it][in][0];
	  a1=d[it][in][1];
	  a2=d[it][in][2];
	  b0=d[it][jn][0];
	  b1=d[it][jn][1];
	  b2=d[it][jn][2];
	  c0=d[jt][0][in];
	  c1=d[jt][1][in];
	  c2=d[jt][2][in];
	  
	  /* check if contact */
	  if((((c0>nsmall)&&(c1>nsmall)&&(c2>nsmall))||
	      ((c0<small)&&(c1<small)&&(c2<small)))||
	      (((a0<small)&&(b0<small))||((a1<small)&&(b1<small))||
	      ((a2<small)&&(b2<small))))
	  { if((a0<=a1)&&(a0<=a2))
	    { dmin2=std::min(dmin2,(a0*a0/n0));
	    }
	    else if((a1<=a0)&&(a1<=a2))
	    { dmin2=std::min(dmin2,(a1*a1/n1));
	    }
	    else
	    { dmin2=std::min(dmin2,(a2*a2/n2));
	    }
	   }
              else            
              { icontact=it;
                /* domain of contact  */
                smin=R0; smax=R1;
                if((a0<R0)&&(b0>small))smin=std::max(smin,(a0/(a0-b0)));
                if((a1<R0)&&(b1>small))smin=std::max(smin,(a1/(a1-b1)));
                if((a2<R0)&&(b2>small))smin=std::max(smin,(a2/(a2-b2)));
                if((a0>small)&&(b0<R0))smax=std::min(smax,(a0/(a0-b0)));
                if((a1>small)&&(b1<R0))smax=std::min(smax,(a1/(a1-b1)));
                if((a2>small)&&(b2<R0))smax=std::min(smax,(a2/(a2-b2)));
                if(smax>smin)
                { s[0]=smin;
                  p[0]=std::min((a0+smin*(b0-a0)),(a1+smin*(b1-a1)));
                  p[0]=std::min(p[0],(a2+smin*(b2-a2)));
                  np=1;
                  /* intermediate points */
                  tmp=b0-a0+a1-b1;
                  if((std::fabs(tmp))>small)
                  { tmp=(a1-a0)/tmp;
                    if((tmp>smin)&&(tmp<smax)&&
                       ((a0+tmp*(b0-a0))<(a2+tmp*(b2-a2))))
                    { s[np]=tmp;
                      p[np]=a0+tmp*(b0-a0);
                      np=np+1;
                  } }
                  tmp=b0-a0+a2-b2;
                  if((std::fabs(tmp))>small)
                  { tmp=(a2-a0)/tmp;
                    if((tmp>smin)&&(tmp<smax)&&
                       ((a0+tmp*(b0-a0))<(a1+tmp*(b1-a1))))
                    { s[np]=tmp;
                      p[np]=a0+tmp*(b0-a0);
                      np=np+1;
                  } }
                  tmp=b1-a1+a2-b2;
                  if((std::fabs(tmp))>small)
                  { tmp=(a2-a1)/tmp;
                    if((tmp>smin)&&(tmp<smax)&&
                       ((a1+tmp*(b1-a1))<(a0+tmp*(b0-a0))))
                    { s[np]=tmp;
                      p[np]=a1+tmp*(b1-a1);
                      np=np+1;
                  } }
                  s[np]=smax;
                  p[np]=std::min((a0+smax*(b0-a0)),(a1+smax*(b1-a1)));
                  p[np]=std::min(p[np],(a2+smax*(b2-a2)));
                  np=np+1;
                  /* order intermediate points */
                  for(unsigned ip=0;ip<(np-1);ip++)
                  { for(unsigned int jp=(ip+1);jp<np;jp++)
                    { if(s[ip]>s[jp])
                      { tmp=s[jp]; s[jp]=s[ip]; s[ip]=tmp;
                        tmp=p[jp]; p[jp]=p[ip]; p[ip]=tmp;
                  } } }
                  /* integrate normal force */
                  fn=p[0]*(s[1]-s[0])+p[np-1]*(s[np-1]-s[np-2]);
                  fnb=p[0]*(s[1]-s[0])*(s[1]+R2*s[0])+
                      p[np-1]*(s[np-1]-s[np-2])*(s[np-2]+R2*s[np-1]);
                  for(unsigned int ip=1;ip<(np-1);ip++)
                  { fn=fn+p[ip]*(s[ip+1]-s[ip-1]);
                    fnb=fnb+p[ip]*(
                    (s[ip]-s[ip-1])*(s[ip-1]+R2*s[ip])+
                    (s[ip+1]-s[ip])*(s[ip+1]+R2*s[ip]));
                  }
                  fnb=fnb*pen*RP5;
                  fn=fn*pen*RP15;
                  fna=fn-fnb;
                  /* update total force */
                  fx[in]=fx[in]-fna*nx[it][in];
                  fy[in]=fy[in]-fna*ny[it][in];
                  fx[jn]=fx[jn]-fnb*nx[it][in];
                  fy[jn]=fy[jn]-fnb*ny[it][in];               
            } } }
            
            //if(icontact==it) /* update nodal forces  */
            { 
	      Element::GeometryType& this_geom_1 = (*Geom[it]);
	      Element::GeometryType& this_geom_2 = (*Geom[jt]);
	      for(unsigned int in=0;in<3;in++)
              { 
		array_1d<double,3>& node_rhs_1 = this_geom_1(in)->FastGetSolutionStepValue(RHS);
		array_1d<double,3>& normal_1   = this_geom_1(in)->FastGetSolutionStepValue(NORMAL);
		this_geom_1[in].SetLock();
		node_rhs_1[0] += 0.50 * fx[in];
                node_rhs_1[1] += 0.50 * fy[in];
		node_rhs_1[2]  =  0.00;
		normal_1[0]   += 0.50 * fx[in]; 
		normal_1[1]   += 0.50 * fy[in]; 
		normal_1[2]    = 0.00; 
		this_geom_1[in].UnSetLock();
                unsigned int ie=in+1; if(ie>2)ie=0;
                for(unsigned int jn=0;jn<3;jn++)
                { 
		  array_1d<double,3>& node_rhs_2 = this_geom_2(in)->FastGetSolutionStepValue(RHS);
		  array_1d<double,3>& normal_2   = this_geom_2(in)->FastGetSolutionStepValue(NORMAL);
		  this_geom_2[jn].SetLock();
		  node_rhs_2[0] -= 0.50 * fx[jn]*d[it][jn][ie];
                  node_rhs_2[1] -= 0.50 * fy[jn]*d[it][jn][ie];
		  node_rhs_2[2] = 0.00;
		  normal_2[0]   -= 0.50 * fx[jn]*d[it][jn][ie];
	 	  normal_2[1]   -= 0.50 * fy[jn]*d[it][jn][ie];
		  normal_2[2]    = 0.00; 		  
		  this_geom_2[jn].UnSetLock();
                } } } }
    KRATOS_CATCH("")
}


/* tetrahedra to tetrahedra */
void ComputeContactForce3D(const PointerType& Target, const PointerType& Contactor)
{
  
	const double R0   = 0.00;
	const double R1   = 1.00;
	const double R2   = 2.00;
	const double R5   = 5.00;
	const double RP1  = 0.10;
	const double RP25 = 0.25; 
        const double RP5  = 0.50;
	
	
	double tmp,theigh,penetr,peneto,penetu,penetv,penalty;
	double force,forco,uforc,vforc,factor,fact0,facti,fact1;
	double xorig,yorig,zorig,xe[2],ye[2],ze[2],dct[4];
	double dsc[6][3],dcs[3][6],us[6],vs[6],ub[10],vb[10],anb[10],penetb[10];
	double xt[4],yt[4],zt[4],ut[4],vt[4],ft[4],xcent,ycent,zcent,xnt,ynt,znt; 
	double xc[4],yc[4],zc[4],uc[4],vc[4],fc[4],xcenc,ycenc,zcenc,xnc,ync,znc;
	double /*zone2,dmin2,*/factor1;

	long /*kprop,icontact,ielem,jelem,icoup,jcoup,*/fnonzero;
	long i,j,k,inext,jnext,itars,icons;
	long nspoin,ninerc,niners,nbpoin,innerc[3],inners[6];
	//long itarth,iconth;
	long iptn[4],ipcn[4];
	long iptn1[4],ipcn1[4],m;

	NodePointerType ipt[4], ipc[4];

	double  pen_vec_tar =  50.00 * (Target->GetProperties()[YOUNG_MODULUS]);
	double  pen_vec_con =  50.00 * (Contactor->GetProperties()[YOUNG_MODULUS]);
	Target->GetValue(IS_TARGET) = true;
	penalty = std::min(pen_vec_tar, pen_vec_con); ///penalty term
	Element::GeometryType& Tgeom = Target->GetGeometry();
	Element::GeometryType& Cgeom = Contactor->GetGeometry();
	//std::vector<Element::GeometryType::Pointer> Geom(2);
	//Geom[0] = Target->pGetGeometry();
	//Geom[1] = Contactor->pGetGeometry();
	  
        /*set centres of contactor and target object */
	xcent=R0; ycent=R0; zcent=R0; xcenc=R0; ycenc=R0; zcenc=R0;
	for(i=0;i<4;i++)
	{     xcenc=xcenc+RP25*(Cgeom(i)->X()); 
	      ycenc=ycenc+RP25*(Cgeom(i)->Y());
	      zcenc=zcenc+RP25*(Cgeom(i)->Z());
	      xcent=xcent+RP25*(Tgeom(i)->X());
	      ycent=ycent+RP25*(Tgeom(i)->Y());
	      zcent=zcent+RP25*(Tgeom(i)->Z());
	}
	
	/*********************************************************/     
	/*                loop over target surfaces              */
	/*********************************************************/
	//ipt->Guarda las conectividades del elemento
	for(itars=0;itars<4;itars++)
	{ 
	  ipt[0]  = Tgeom(itars);      //i2elto[itars][itarth];
	  iptn1[0]= itars;
	  ipt[1]  = Tgeom(1);          //i2elto[1][itarth];
	  iptn1[1]= 1;
	  ipt[2]  = Tgeom(2);          //i2elto[2][itarth];
	  iptn1[2]= 2;
	  if(itars>0)
	  { 
	    ipt[3]  = Tgeom(itars-1);   //i2elto[itars-1][itarth];
	    iptn1[3]= itars-1;
	  }
	  else
	  { 
	    ipt[3]  = Tgeom(3);          //i2elto[3][itarth];
	    iptn1[3]= 3;
	  }    
	  if((itars==1)||(itars==2))
	  {
	    ipt[1] = Tgeom(3);         //i2elto[3][itarth];
	    iptn1[1]=3;
	  }
	  if(itars>1)
	  {
	    ipt[2] = Tgeom(0);         //i2elto[0][itarth];
	    iptn1[2]= 0; 
	  }
	  
	/*****************************************************/ 
	/*           loop over contactor surfaces            */
	/*****************************************************/ 
	
	for(icons=0;icons<4;icons++)
	{ 
	 ipc[0]   = Cgeom(icons);            //i2elto[icons][iconth];
	 ipcn1[0] = icons;
	 ipc[1]   = Cgeom(1);                //i2elto[1][iconth];
	 ipcn1[1] = 1;
	 ipc[2]   = Cgeom(2);                //i2elto[2][iconth];
	 ipcn1[2] = 2;
	 if(icons>0)
	 { 
	  ipc[3]   = Cgeom(icons-1);        //i2elto[icons-1][iconth];
	  ipcn1[3] = icons-1;
	 }
	 else
	 { 
	  ipc[3]   = Cgeom(3);              //i2elto[3][iconth];
	  ipcn1[3] = 3;
	 }
	if((icons==1)||(icons==2))
	 {
	  ipc[1]  = Cgeom(3);               //i2elto[3][iconth];
	  ipcn1[1]= 3;
	 }
	if(icons>1)
	{
	  ipc[2]   = Cgeom(0);               //i2elto[0][iconth];
	  ipcn1[2] = 0;
	}

	for(m=0;m<4;m++)
	{
	  iptn[iptn1[m]]=m;
	  ipcn[ipcn1[m]]=m;
	}
	
	/* set nodal coordinates */
	for(i=0;i<3;i++)
	{ 
	  xt[i] = ipt[i]->X();
	  yt[i] = ipt[i]->Y();
	  zt[i] = ipt[i]->Z();
	  xc[i] = ipc[i]->X();
	  yc[i] = ipc[i]->Y();
	  zc[i] = ipc[i]->Z();
	}
	
	xt[3]=xcent; yt[3]=ycent; zt[3]=zcent;
	xc[3]=xcenc; yc[3]=ycenc; zc[3]=zcenc;
	xorig=xc[0]; yorig=yc[0]; zorig=zc[0];

	for(i=0;i<4;i++)
	{ 
	  xt[i]=xt[i]-xorig; yt[i]=yt[i]-yorig; zt[i]=zt[i]-zorig;
	  xc[i]=xc[i]-xorig; yc[i]=yc[i]-yorig; zc[i]=zc[i]-zorig; 
	} 
	
	/* contactor normal, e-base and target points in e-base */
	V3DCro(xnc,ync,znc,xc[1],yc[1],zc[1],xc[2],yc[2],zc[2]);
	V3DNor(xe[0],xnc,ync,znc);
	
	
	xe[0]=xc[1]; ye[0]=yc[1]; ze[0]=zc[1];
	V3DNor(xe[1],xe[0],ye[0],ze[0]); 
	V3DCro(xe[1],ye[1],ze[1],xnc,ync,znc,xe[0],ye[0],ze[0]);
	for(i=0;i<4;i++)
	{ 
	  V3DDot(dct[i],xnc,ync,znc,xt[i],yt[i],zt[i]);
	  V3DDot(ut[i],xt[i],yt[i],zt[i],xe[0],ye[0],ze[0]);
	  V3DDot(vt[i],xt[i],yt[i],zt[i],xe[1],ye[1],ze[1]);
	}
	
	/* u,v coordinates of S-points and C-points   */
	nspoin=0;
	for(i=0;i<3;i++)
	 { 
	   for(j=0;j<2;j++)
	    { 
	      inext=i+1; if(inext>2)inext=0; if(j==0)inext=3;
	      if(((dct[i]>EPSILON)&&(dct[inext]<NEPSILON))|| ((dct[i]<NEPSILON)&&(dct[inext]>EPSILON)))
	      //Modified by JXiang
	       { 
		 factor=std::fabs(dct[i]-dct[inext]);          
	         if(factor>EPSILON)
	         { 
		   factor=std::fabs(dct[i]/factor);
	           us[nspoin]=factor*ut[inext]+(R1-factor)*ut[i];
	           vs[nspoin]=factor*vt[inext]+(R1-factor)*vt[i];
	           inners[nspoin]=0;
	           nspoin=nspoin+1;
	
		 } 
	       } 
	    } 
	 }
	 
	if((nspoin<3)||(nspoin>4))continue;
	
	/* check odering of S-points  */
	if(((us[1]-us[0])*(vs[2]-vs[0])-(vs[1]-vs[0])*(us[2]-us[0]))<R0)
	{ 
	  i=0; j=nspoin-1;
	  while(i<j)
	  { 
	   k=inners[i]; inners[i]=inners[j]; inners[j]=k;
	   tmp=us[i];   us[i]=us[j];         us[j]=tmp;
	   tmp=vs[i];   vs[i]=vs[j];         vs[j]=tmp;
	   i++; j--;
	
	  } 
	}
	
	for(i=0;i<3;i++)
	{ 
	  V3DDot(uc[i],xc[i],yc[i],zc[i],xe[0],ye[0],ze[0]);
	  V3DDot(vc[i],xc[i],yc[i],zc[i],xe[1],ye[1],ze[1]);
	  innerc[i]=0;
	}
	
	/* distances of C-points from S edges */
	niners=0; ninerc=0;
	for(i=0;i<nspoin;i++)
	{ 
	  inext=i+1;
	  if(inext>=nspoin)inext=0;
	  for(j=0;j<3;j++) 
	  { 
	    jnext=j+1;
	    if(jnext>2)jnext=0;
	    dcs[j][i]=(uc[jnext]-uc[j])*(vs[i]-vc[j])-(vc[jnext]-vc[j])*(us[i]-uc[j]);
	    dsc[i][j]=(us[inext]-us[i])*(vc[j]-vs[i])-(vs[inext]-vs[i])*(uc[j]-us[i]);
	    if(dsc[i][j]>=R0)
	    { 
	      innerc[j]=innerc[j]+1;
  	      if(innerc[j]==nspoin) ninerc=ninerc+1;
	    }
	    if(dcs[j][i]>=R0)
	   { 
	     inners[i]=inners[i]+1;
	     if(inners[i]==3) niners = niners+1;
	
	   } 
	  } 
	}
	
	/* B-points   */
	if(ninerc==3)           /* triangle inside poligon      */
	{ 
	  nbpoin=3;
	  for(i=0;i<nbpoin;i++)
	  { 
	    ub[i]=uc[i]; vb[i]=vc[i];
	  } 
	}
	else if(niners==nspoin) /* poligon inside triangle */     
	{ 
	  nbpoin=nspoin;
	  for(i=0;i<nbpoin;i++)
	  { 
	    ub[i]=us[i]; vb[i]=vs[i];
	  } 
	}
	else            /* intersection points poligon triangle */
	{ 
	  nbpoin=0;
	  for(i=0;i<nspoin;i++)
	   { 
	     if(inners[i]==3)
	     { 
	       ub[nbpoin]=us[i]; vb[nbpoin]=vs[i]; nbpoin++; 
	
	     } 
	   }
	   for(i=0;i<3;i++)  /* grab inner C-points */ 
	   { 
	     if(innerc[i]==nspoin)  
	     { 
	     ub[nbpoin]=uc[i]; vb[nbpoin]=vc[i]; nbpoin++;       
	
	     } 
	   }       
	for(i=0;i<nspoin;i++)        /* intersection points */   
	{ 
	  inext=i+1; if(inext>=nspoin)inext=0;
	  for(j=0;j<3;j++)
	  {
	    jnext=j+1; if(jnext>2)jnext=0;
	   if((((dsc[i][j]>EPSILON)&&(dsc[i][jnext]<NEPSILON))||
	   ((dsc[i][j]<NEPSILON)&&(dsc[i][jnext]>EPSILON)))&&
	   (((dcs[j][i]>EPSILON)&&(dcs[j][inext]<NEPSILON))||
 	   ((dcs[j][i]<NEPSILON)&&(dcs[j][inext]>EPSILON))))
	    //modified by JXiang
	    { 
	     factor=std::fabs(dsc[i][j]-dsc[i][jnext]);
	     if(factor<EPSILON){ factor=RP5;}
	     else{factor=std::fabs(dsc[i][j]/factor);}
	     ub[nbpoin]=(R1-factor)*uc[j]+factor*uc[jnext];
	     vb[nbpoin]=(R1-factor)*vc[j]+factor*vc[jnext];
	     nbpoin++;
	 
	    } 
	  } 
	}
	
	for(i=1;i<nbpoin;i++)
	{ 
	  if(vb[i]<vb[0])
	   { 
	     tmp=vb[i]; vb[i]=vb[0]; vb[0]=tmp;
	     tmp=ub[i]; ub[i]=ub[0]; ub[0]=tmp;
	
	   } 
	}
	
	for(i=1;i<nbpoin;i++)
	{ 
	  tmp=ub[i]-ub[0];				 
	   if((tmp<R0)&&(tmp>(-EPSILON)))
	   { 
	     tmp=tmp-EPSILON;
	   }
	else if((tmp>=R0)&&(tmp<EPSILON))
	 { 
	   tmp=tmp+EPSILON;
	 }
	  anb[i]=(vb[i]-vb[0]+EPSILON)/tmp;
	}
	
	for(i=1;i<nbpoin;i++)  /* sort B-points */
	{ 
	  for(j=i+1;j<nbpoin;j++)
	  { 
	    if(((anb[i]>=R0)&&(anb[j]>=R0)&&(anb[j]<anb[i]))||
	    ((anb[i]<R0)&&((anb[j]>=R0)||(anb[j]<anb[i]))))
	     { 
	      tmp=vb[i];  vb[i]=vb[j];   vb[j]=tmp;
	      tmp=ub[i];  ub[i]=ub[j];   ub[j]=tmp;
	      tmp=anb[i]; anb[i]=anb[j]; anb[j]=tmp;
	
	      } 
	    } 
	  } 
	}
	
	if(nbpoin<3)continue; 
	/* Target-plain normal and penetration at B-points      */      
	V3DCro(xnt,ynt,znt,xt[1]-xt[0],yt[1]-yt[0],zt[1]-zt[0],xt[2]-xt[0],yt[2]-yt[0],zt[2]-zt[0]); 
	V3DDot(theigh,xt[3]-xt[0],yt[3]-yt[0],zt[3]-zt[0],xnt,ynt,znt); 
	/* penetration at origin of the e-base and dp/du dp/dv; */
	V3DDot(peneto,xc[0]-xt[0],yc[0]-yt[0],zc[0]-zt[0],xnt,ynt,znt); 
	V3DDot(penetu,xe[0],ye[0],ze[0],xnt,ynt,znt);
	V3DDot(penetv,xe[1],ye[1],ze[1],xnt,ynt,znt);
	peneto=peneto/theigh; 
	penetu=penetu/theigh; 
	penetv=penetv/theigh;
	for(i=0;i<nbpoin;i++)
	{ 
	  penetb[i]=peneto+ub[i]*penetu+vb[i]*penetv;
	}
	/* force and center of force */
	forco=R0; uforc=R0; vforc=R0;   
	for(i=1;i<(nbpoin-1);i++)
	{ 
	  penetr=penetb[0]+penetb[i]+penetb[i+1];
	  if(penetr>EPSILON){
	  force=((ub[i]-ub[0])*(vb[i+1]-vb[0])-(vb[i]-vb[0])*(ub[i+1]-ub[0]))*penetr*penalty;
	  fact0=(RP5*penetb[0]+RP25*(penetb[i]+penetb[i+1]))/penetr;
	  facti=(RP5*penetb[i]+RP25*(penetb[0]+
	  penetb[i+1]))/penetr;
	  fact1=R1-fact0-facti;
	  if(std::fabs(force+forco)>EPSILON)
	  { 
	    uforc=(forco*uforc+force*(fact0*ub[0]+
	    facti*ub[i]+fact1*ub[i+1]))/(forco+force); 
	    vforc=(forco*vforc+force*(fact0*vb[0]+
	    facti*vb[i]+fact1*vb[i+1]))/(forco+force);
	    forco=forco+force;
	
	   } 
	  } 
	} 
	
	/*resultant at C-points */
	for(i=0;i<4;i++)
	{ 
	  fc[i]=R0; ft[i]=R0;
	}
	
	tmp=((uc[1]-uc[0])*(vc[2]-vc[0])-
	(vc[1]-vc[0])*(uc[2]-uc[0]));
	for(i=0;i<3;i++)
	{ j=i+1; if(j>2)j=0; k=j+1; if(k>2)k=0;
	fc[k]=forco*(((uc[j]-uc[i])*(vforc-vc[i])-
	(vc[j]-vc[i])*(uforc-uc[i]))/tmp); 
	}
	
	/*resultant at T-points*/
	
	tmp=((ut[1]-ut[0])*(vt[2]-vt[0])-(vt[1]-vt[0])*(ut[2]-ut[0]));
	inext=-1;
	if(std::fabs(tmp)<RP1*theigh)
	{ 
	  inext=0; tmp=std::fabs(ut[1]-ut[0])+std::fabs(vt[1]-vt[0]); 
	  for(i=0;i<3;i++)
	  { 
	    j=i+1; if(j>2)j=0;
 	    if(tmp>(std::fabs(ut[j]-ut[i])+std::fabs(vt[j]-vt[i])))
	    { 
	      tmp=std::fabs(ut[j]-ut[i])+std::fabs(vt[j]-vt[i]);  inext=i;
	    }
	  }
	  j=inext+1; if(j>2)j=0;
	  if(std::fabs(zt[j])>std::fabs(zt[inext]))inext=j;
	  j=inext+1; if(j>2)j=0;  k=j+1;  if(k>2)k=0;  
	  tmp=(ut[k]-ut[j])*(vt[3]-vt[j])-
	  (vt[k]-vt[j])*(ut[3]-ut[j]);
	}
	
	for(jnext=0;jnext<3;jnext++)
	{ i=jnext; j=i+1; if(j>2)j=0; k=j+1; if(k>2)k=0;
	if(i==inext)i=3; if(j==inext)j=3; if(k==inext)k=3; 
	ft[k]=forco*(((ut[j]-ut[i])*(vforc-vt[i])-
	(vt[j]-vt[i])*(uforc-ut[i]))/tmp);                 
	}
	ft[3]=RP25*ft[3];
	for(i=0;i<3;i++)
	{ ft[i]=ft[i]+ft[3];
	}
	
	/* add forces into global vector    */
	factor1=R2/R5;
	fnonzero=1;
	
	for(i=0;i<4;i++)
	{ 
	  array_1d<double,3>& node_rhs_1 = (ipc[i])->FastGetSolutionStepValue(RHS);
	  array_1d<double,3>& node_rhs_2 = (ipt[i])->FastGetSolutionStepValue(RHS);
	  
	  node_rhs_1[0] += fc[i]*xnc*factor1;
	  node_rhs_1[1] += fc[i]*ync*factor1;
	  node_rhs_1[2] += fc[i]*znc*factor1;
	  node_rhs_2[0] -= ft[i]*xnc*factor1;
	  node_rhs_2[1] -= ft[i]*ync*factor1;
	  node_rhs_2[2] -= ft[i]*znc*factor1;
	}
     }
  }
}
  

//************************************************************************************
//************************************************************************************            

      // Funcion que se llama antes de la acualizacion de los desplazamientos
      void LocalSearch() 
      {
	  KRATOS_TRY
	  
	    IteratorType it_begin     =  mBoundaryElements.begin();
            IteratorType it_end       =  mBoundaryElements.end();  
	    BinsObjectDynamic<Configure>  rBinsObjectDynamic(it_begin, it_end ); 
	    if(mrdimension==2){
	      SearchNearNode2D(rBinsObjectDynamic, it_begin, it_end); 
	      LocalSearch2D(rBinsObjectDynamic, it_begin, it_end); 
	    }
	    else
	      LocalSearch3D(rBinsObjectDynamic, it_begin, it_end);
	    
	  KRATOS_CATCH("")
      }
       
       
   
//************************************************************************************
//************************************************************************************  

         bool SearchContactsPairs()
         {
	    KRATOS_TRY
	   
		std::cout<< std::endl;
		std::cout<<"  COMPUTING CONTACT CONDITIONS TO MODEL PART " << std::endl; 
		
		IteratorType it_begin     =  mBoundaryElements.begin();
		IteratorType it_end       =  mBoundaryElements.end(); 

		BinsObjectDynamic<Configure>  rBinsObjectDynamic(it_begin, it_end ); 
		rBinsObjectDynamic.SearchContact(mPairContacts);
		if(mrdimension==2){  
		    LocalSearch2D(rBinsObjectDynamic, it_begin, it_end); 
	            FiltratePairContacts2D(mPairContacts);
		}
		else
		{ 
		    //LocalSearch3D(rBinsObjectDynamic, it_begin, it_end); 
		    //std::cout<< "     NUMBER OF CONTACT PAIRS                 = " <<mPairContacts.size()<<std::endl; 
	            //FiltratePairContacts3D(mPairContacts);
		    std::cout<< "     NUMBER OF CONTACT PAIRS                 = " <<mPairContacts.size()<<std::endl; 
		}
		
		if(mPairContacts.size()!=0)
		{
		  std::cout<< "     NUMBER OF CONTACT PAIRS                 = " <<mPairContacts.size()<<std::endl; 
		  //KRATOS_ERROR(std::logic_error,  "GetValue", "");
		  return true;
		}
		
		std::cout<< "     NO CONTACTS PAIRS "<<std::endl; 
		return false;  
	    
	    KRATOS_CATCH("")
	 }
	 
	 
 
  //************************************************************************************
  //************************************************************************************ 

  void ResetValues()
  {
    
    	   KRATOS_TRY
           NodesArrayType& pNodes           =  mr_model_part.Nodes();
	   #ifdef _OPENMP
           int number_of_threads = omp_get_max_threads();
           #else
           int number_of_threads = 1;
           #endif

           vector<unsigned int> node_partition;
           CreatePartition(number_of_threads, pNodes.size(), node_partition);

           #pragma omp parallel for 
           for(int k=0; k<number_of_threads; k++)
             {
	        NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	        NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];
	        for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)      
	          {  
	 	      i->GetValue(IS_CONTACT_SLAVE)  = 0;
		      i->GetValue(IS_CONTACT_MASTER) = 0;
		      i->GetValue(NODAL_VALUES)      = 0;
		      i->GetValue(DISTANCE)          = DBL_MAX;
	          }
	     }
	     
	     KRATOS_CATCH("")  
  }


  //************************************************************************************
  //************************************************************************************ 

	 void Clear(const unsigned int&  initial_conditions_size)
	 {
	   
	   KRATOS_TRY
           NodesArrayType& pNodes           =  mr_model_part.Nodes();
	   ElementsArrayType& pElements     =  mr_model_part.Elements(); 
	   ConditionsArrayType& pConditions =  mr_model_part.Conditions(); 
	   //ProcessInfo& CurrentProcessInfo  =  mr_model_part.GetProcessInfo();
	   
	   
	   #ifdef _OPENMP
           int number_of_threads = omp_get_max_threads();
           #else
           int number_of_threads = 1;
           #endif

           vector<unsigned int> node_partition;
           CreatePartition(number_of_threads, pNodes.size(), node_partition);

           #pragma omp parallel for 
           for(int k=0; k<number_of_threads; k++)
             {
	        NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	        NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];
	        for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)      
	          {  
	 	      i->GetValue(IS_CONTACT_SLAVE)  = 0;
		      i->GetValue(IS_CONTACT_MASTER) = 0;
		      i->GetValue(NODAL_VALUES)      = 0;
		      i->GetValue(DISTANCE)          = DBL_MAX;
		      i->GetValue(NEAR_NODE)         = *(i.base());
	          }
	     }
	     
	     
	   //No se han producido nuevas condiciones de contorno
           if(mcompute_boundary_contour==true)
	   {
	      /// Borro las condiciones masters
	      if(pConditions.size()>initial_conditions_size){
	      ModelPart::ConditionIterator   end_previos  = pConditions.begin() + initial_conditions_size;
	      ModelPart::ConditionIterator   end_actual   = pConditions.end();
	      pConditions.erase(end_previos, end_actual);
	      }
	      
	      FindElementalNeighboursProcess ElementosVecinos(mr_model_part, mrdimension, 10);
	      FindNodalNeighboursProcess NodosVecinos(mr_model_part, mrdimension, 10);
	      FindConditionsNeighboursProcess CondicionesVecinas(mr_model_part, mrdimension, 10);

	      ElementosVecinos.ClearNeighbours();
	      NodosVecinos.ClearNeighbours();
	      CondicionesVecinas.ClearNeighbours();
	      ElementosVecinos.Execute();
	      NodosVecinos.Execute();
	      CondicionesVecinas.Execute();
	      
	      vector<unsigned int> condition_partition;
	      CreatePartition(number_of_threads, pConditions.size(), condition_partition);
	      #pragma omp parallel for
	      for(int k=0; k<number_of_threads; k++)
	      {
	         ConditionsArrayType::iterator it_begin=pConditions.ptr_begin()+condition_partition[k];
	         ConditionsArrayType::iterator it_end=pConditions.ptr_begin()+condition_partition[k+1];
	         for (ConditionsArrayType::iterator it= it_begin; it!=it_end; ++it)
		 {
		    WeakPointerVector<Element>& rC = it->GetValue(NEIGHBOUR_ELEMENTS);
		    rC.erase(rC.begin(),rC.end() );	
		 }
	      }
	      
	      vector<unsigned int> element_partition;
	      CreatePartition(number_of_threads, pElements.size(), element_partition);
	      #pragma omp parallel for
	      for(int k=0; k<number_of_threads; k++)
	      {
	         ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
	         ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];
	         for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
	          {
		    WeakPointerVector<Condition> & neighb_conds = it->GetValue(NEIGHBOUR_CONDITIONS);
		    for (WeakPointerVector< Condition >::iterator neighb_cond = neighb_conds.begin(); neighb_cond != neighb_conds.end(); neighb_cond++)
		        (neighb_cond->GetValue(NEIGHBOUR_ELEMENTS)).push_back( *(it.base()));
		    
		  }
	      }
	      
	      mPairContacts.clear(); 
	      mBoundaryElements.clear();
	      mMasterConditionsArray.clear();
	   }  
	   
	   KRATOS_CATCH("") 
	       
	 }
	  

//************************************************************************************
//************************************************************************************  

       void CreateLinkingConditionsBasedOnLocalSearch(const unsigned int&  initial_conditions_size) 
       {
	 KRATOS_TRY	  
	 if(mrdimension==2)
	    CreateLinkingConditionsBasedOnLocalSearch2D(initial_conditions_size);
	 else
	   CreateLinkingConditionsBasedOnLocalSearch3D(initial_conditions_size);
	 KRATOS_CATCH("")  
       }


//************************************************************************************
//************************************************************************************  

	 void CreateLinkingConditionsBasedOnLocalSearch2D(const unsigned int&  initial_conditions_size) 
	 {
	   
	    KRATOS_TRY
                  
            
            ConditionsArrayType& rConditions = mr_model_part.Conditions();
	    IntersectTriangleCases<Configure> IntersectTriangles(mr_model_part); 
  
	    array_1d<NodePointerType, 2>                Ids;      
	    std::vector<NodePointerType>                InsideNodes;
	    std::vector<array_1d<NodePointerType, 2 > > Ids_2;
	    std::vector<Near_Node>                      Is_Near;
	    std::vector<ConditionsArrayType>            LinkingConditions;
	    
	    unsigned int Id                        = rConditions.size() + 1;
	    unsigned int properties_index          = mr_model_part.NumberOfProperties();	    
	    PropertiesType::Pointer tempProperties = PropertiesType::Pointer(new PropertiesType(properties_index+1));
	    //mr_model_part.AddProperties(tempProperties);
 
	    int  Case                   =  0;    
	    unsigned int  master        =  0;
	    unsigned int  slave         =  1;    
	    bool is_repited             =  false;
	    bool corner                 =  false; 
	    bool Change                 =  true;    
            Near_Node  Near             =  no_near; 
            Exist_Node Exist            =  no_nodes;
	    
	    NodePointerType             Id_Node_Case_5; 
	    
	    #ifdef _OPENMP
            int number_of_threads = omp_get_max_threads();
            #else
            int number_of_threads = 1;
	    #endif
	    
            //std::vector<ResultContainerType> Result(number_of_threads);
            ResultContainerType Result;
 
	    
            /// creando bins de objetos
	    BinsObjectDynamic<Configure>  rBinsObjectDynamic(mBoundaryElements.begin(), mBoundaryElements.end()); 
	    
	    ///bins de puntos
	    //BinsDynamic<2, NodeType, NodesContainerType> BinsPoint(mBoundaryNodes.begin(), mBoundaryNodes.end());
	    
	    /// busqueda local de los segmentos mas cercanos a un nodo
	    //SearchNearNode2D();
	    SearchNearNode2D(rBinsObjectDynamic, mBoundaryElements.begin(), mBoundaryElements.end()); 
	    IdentifyMasterSegment2D();
	    //LocalSearch2D(rBinsObjectDynamic, mBoundaryElements.begin(), mBoundaryElements.end()); 
	      
	    
            LinkingConditions.resize(number_of_threads); 
            vector<unsigned int> partition;
            CreatePartition(number_of_threads, mBoundaryElements.size(), partition);
	    ContactPairType it_pair;  
	    
	    std::cout<<"     PARTITION COMPUTING CONTACT CONDITIONS  = " << number_of_threads << std::endl; 
	    std::cout<<"     NUMBER OF INITIAL CONDITIONS            = " << initial_conditions_size <<  std::endl;
            std::cout<<"     NUMBER OF MASTER SURFACES CONDITIONS    = " << rConditions.size()-initial_conditions_size <<  std::endl;
	    
	    #pragma omp parallel for  firstprivate (Case, Id_Node_Case_5, master, slave, is_repited, corner, Change, Near, Exist) private (Result, Id, it_pair, Ids, InsideNodes, Ids_2, Is_Near)
	    for(int k=0; k<number_of_threads; k++)          
	     {
	       IteratorType it_begin = mBoundaryElements.begin() + partition[k];
	       IteratorType it_end   = mBoundaryElements.begin() + partition[k+1];
	       for(IteratorType it =it_begin; it!=it_end; it++)
	         { 
		   Result.clear();
		   Result.reserve(100);
		   (*it)->GetValue(IS_TARGET)=true;
		   rBinsObjectDynamic.SearchObjectsInner(*it, Result);      
		   if(Result.size()!=0){    
		   for(ResultIteratorType rthis = Result.begin(); rthis!=Result.end(); ++rthis){ 
			 if(FiltratePairContacts(*it, *rthis)==true)
			 {    
			    Exist            =  no_nodes;
			    Case             =  0;       
			    master           =  0;
			    slave            =  1;
			    (it_pair)[master] = (*rthis); 
			    (it_pair)[slave]  = (*it); 
			    Ids_2.clear();
			    InsideNodes.clear(); 
			    NodeInside((it_pair)[0], (it_pair)[1], InsideNodes);  
			    if(InsideNodes.size()==0)
			    {
			      InsideNodes.clear();
			      ContactPairType it_pair_2;
			      (it_pair_2)[master] = (*it); 
			      (it_pair_2)[slave]  = (*rthis); 
			      NodeInside((it_pair_2)[0], (it_pair_2)[1], InsideNodes);  
			      if(InsideNodes.size()==0)
			         { Exist = no_nodes;} 
		              else{
		                 Exist = yes_nodes;
				 continue;
			        }
			     }
			   else
			      Exist = yes_nodes; 
			  
			   
			    switch(Exist)
			    {
			      case(yes_nodes):
			        {  
				Case = IntersectTriangles.LocateCaseItersection(Id_Node_Case_5, Change, InsideNodes, (it_pair)[master], (it_pair)[slave]);
				switch(Case)  
				{ 
				  /*
				  case 1: /// un solo nodo dentro
				  {  
				      Near = CheckNearNodes(master, slave, InsideNodes[0], (it_pair)[master], Ids);
				      if(CreateLinkingConditions(Id,  master, slave, Ids, it_pair, tempProperties, Exist, LinkingConditions[k])){
					Id++;
				      }
				      break;
				  }
				  */
				  case 1 : case 2:  /// dos nodos dentro
				  {
				  for(unsigned int in = 0; in<InsideNodes.size(); in++){
				          {
					   Near              = CheckNearNodes(master, slave, InsideNodes[in], (it_pair)[master], Ids);
					   if(CreateLinkingConditions(Id,  master, slave, Ids, it_pair, tempProperties, Exist, LinkingConditions[k])){ 
					     Id++;
				            }
				           }				           
				         }
				         break;
				  }
				  case 3:
				  {
				    break;
				  }  

				  case 5:
				   {
				    break;
				  }
				  
				}
				
				break;
			      }

			      case(no_nodes):
			      { 
				ComputeContactForce2D(((it_pair)[slave]), ((it_pair)[master]));
				/*
				///Penalty
				unsigned int size_master =   ((it_pair)[master])->GetValue(NEIGHBOUR_CONDITIONS).size();
				unsigned int size_slave  =   ((it_pair)[slave])->GetValue(NEIGHBOUR_CONDITIONS).size();
				if(size_master==1 && size_slave==1)
				 {
				   //std::cout<< "     MASTER OBJECT =  " <<  (it_pair)[master]->Id() <<"   SLAVE OBJECT = " << (it_pair)[slave]->Id() << std::endl;
				   CheckNearNodes(master, slave, (it_pair)[slave], (it_pair)[master], Ids_2, Is_Near);
				   //KRATOS_WATCH(Ids_2.size())
				   if(CreateLinkingConditions(Id,  master, slave, Ids_2[0], it_pair, tempProperties, Exist, LinkingConditions[k])){ 
				     Id++; 
				    }
				 }
				else if(size_master>1 || size_slave>1)
				 {
				   ComputeContactForce2D(((it_pair)[slave]), ((it_pair)[master]));
				 }
				 */
			      break;
			       }
			    }
		         }
		       }
		    }
	        }
             }
 
		    unsigned int size = 0;
		    //adding linking to model_part
		    for(int k=0; k<number_of_threads; k++){
			  size+=LinkingConditions[k].size(); 
			  for(ConditionsArrayType::ptr_iterator it=LinkingConditions[k].ptr_begin(); it!= LinkingConditions[k].ptr_end(); ++it )
			      mr_model_part.Conditions().push_back(*it);
		    }
		    
		    std::cout<<"     NUMBER OF LINKING CONTACT CONDITIONS    = " << size               <<  std::endl; 
		    std::cout<<"     TOTAL NUMBER CONDITIONS                 = " << rConditions.size() <<  std::endl; 
		    LinkingConditions.clear();
            KRATOS_CATCH("")
	 }
	 
//************************************************************************************
//************************************************************************************  


void CreateLinkingConditionsBasedOnLocalSearch3D(const unsigned int&  initial_conditions_size) 
	 {
	   
	    KRATOS_TRY
            #ifdef _OPENMP
            int number_of_threads = omp_get_max_threads();
            #else
            int number_of_threads = 1;
	    #endif
            
            ConditionsArrayType& rConditions = mr_model_part.Conditions();
  
	    array_1d<NodePointerType, 2 >            Ids;      
	    std::vector<NodePointerType>             InsideNodes;
	    std::vector<array_1d<unsigned int, 2 > > Ids_2;
	    std::vector<Near_Node>                   Is_Near;
	    std::vector<ConditionsArrayType>         LinkingConditions(number_of_threads);
	       
	    
	    
	    unsigned int Id                        = rConditions.size() + 1;
	    unsigned int properties_index          = mr_model_part.NumberOfProperties();	    
	    PropertiesType::Pointer tempProperties = PropertiesType::Pointer(new PropertiesType(properties_index+1) );
	    //mr_model_part.AddProperties(tempProperties);
 
	    int  Case                   =  0;    
	    unsigned int Id_Node_Case_5 =  0;
	    unsigned int  master        =  0;
	    unsigned int  slave         =  1;    
	    bool is_repited             =  false;
	    bool corner                 =  false; 
	    bool Change                 =  true; 	    
            Near_Node  Near             =  no_near; 
            Exist_Node Exist            =  no_nodes;
	    ResultContainerType         Result;
	    
	    #ifdef _OPENMP
	    double start_prod = omp_get_wtime();   
	    #endif
	    BinsObjectDynamic<Configure>  rBinsObjectDynamic(mBoundaryElements.begin(), mBoundaryElements.end()); 
	    #ifdef _OPENMP
	    double stop_prod  =  omp_get_wtime();
	    std::cout << "        Time creating bins                     = " << stop_prod - start_prod << "  seconds"  <<std::endl;
	    #endif
	    
	    LocalSearch3D(rBinsObjectDynamic, mBoundaryElements.begin(), mBoundaryElements.end()); 

            vector<unsigned int> partition;
            CreatePartition(number_of_threads, mBoundaryElements.size(), partition);
	    ContactPairType it_pair;
	    
	    
	    std::cout<<"        Number of threads used for contact     = " << number_of_threads << std::endl; 
	    std::cout<<"        Number of initial conditions           = " << initial_conditions_size <<  std::endl;
            std::cout<<"        Number of master surface conditions    = " << rConditions.size()-initial_conditions_size <<  std::endl;
	    
	    #ifdef _OPENMP
	    double start = omp_get_wtime();   
	    #endif
	    
	    #pragma omp parallel for  shared(LinkingConditions)  private(Id, it_pair, Ids, InsideNodes, Ids_2, Is_Near, Case, Id_Node_Case_5, master, slave, is_repited, corner, Change, Near,  Exist, Result) 
	    for(int k=0; k<number_of_threads; k++)          
	     {
	       IteratorType it_begin = mBoundaryElements.begin() + partition[k];
	       IteratorType it_end   = mBoundaryElements.begin() + partition[k+1];
	       for(IteratorType it =it_begin; it!=it_end; it++)
	         { 
		   Result.clear();
		   Result.reserve(100);
		   rBinsObjectDynamic.SearchObjectsInner(*it, Result);
		   if(Result.size()!=0){
		   for(ResultIteratorType rthis = Result.begin(); rthis!=Result.end(); rthis++)
		      {
			 if(FiltratePairContacts(*it, *rthis)==true)
			 {
			    Exist            =  no_nodes;
			    Case             =  0;       
			    master           =  0;
			    slave            =  1;
			    (it_pair)[master] = (*rthis); 
			    (it_pair)[slave]  = (*it); 
			    Ids_2.clear();
			    InsideNodes.clear(); 
			    NodeInside((it_pair)[0], (it_pair)[1], InsideNodes);  
			    if(InsideNodes.size()==0)
			     {
			       InsideNodes.clear();
			       ContactPairType it_pair_2;
			       (it_pair_2)[master] = (*it); 
			       (it_pair_2)[slave]  = (*rthis); 
			       NodeInside((it_pair_2)[0], (it_pair_2)[1], InsideNodes);  
			       if(InsideNodes.size()==0)
			        { Exist = no_nodes;} 
			        else{
			          Exist = yes_nodes;
			          continue;
			          }
			      }
			    else
			       Exist = yes_nodes;

			  switch(Exist)
			    {
			    case(yes_nodes):
			    {
				//std::cout<< "     Yes Nodes" << std::endl;
				//std::cout<< "     MASTER OBJECT =  " <<  (it_pair)[master]->Id() <<"   SLAVE OBJECT = " << (it_pair)[slave]->Id() << std::endl;
				Case = InsideNodes.size(); 
				
				switch(Case) 
				{
				  /*
				  case 1: // un solo nodo dentro
				  {  
				      (InsideNodes[0])->GetValue(IS_CONTACT_SLAVE) = 1;
				       Near              = CheckNearNodes(master, slave, InsideNodes[0], (it_pair)[master], Ids);
				       CreateLinkingConditions(Id,  master, slave, Ids, it_pair, tempProperties, Exist, LinkingConditions[k]);  
				       Id++;
				       break;
				       
				  }
				  */
				  case 1: case 2: case 3:
				  {
				  for(unsigned int in = 0; in<InsideNodes.size(); in++){  
				      if(InsideNodes[in]->GetValue(IS_CONTACT_SLAVE)==0)
				      {
					InsideNodes[in]->GetValue(IS_CONTACT_SLAVE) = 1;
					Near              = CheckNearNodes(master, slave, InsideNodes[in], (it_pair)[master], Ids);
					CreateLinkingConditions(Id,  master, slave, Ids, it_pair, tempProperties, Exist, LinkingConditions[k] );  
					Id++;
				      }
				    }
				  break;
				  }
				  
				  case 5:
				  {  
				    break;
				  }
				  
				}
				break;
			      }

			      case(no_nodes):
			      {   
			        break;
			      }
			  }
			  
			  Exist            =  no_nodes;
			  Near             =  no_near;
			  corner           =  false; 
			  Change           =  true; 
			  is_repited       =  false;
			  Case             =  0;       
			  Id_Node_Case_5   =  0;
			  master           =  0;
			  slave            =  1;
			  Ids_2.clear();
			  Is_Near.clear();
			  InsideNodes.clear();
		     }
		    }
	          }
                }
	       } 
 
                    #ifdef _OPENMP
	            double stop = omp_get_wtime();
	            std::cout << "        Time Creating Linking Conditions       = " << stop - start << "  seconds" << std::endl;
		    #endif
		    
		    int rId = rConditions.size() + 1;
		    for(int k=0; k<number_of_threads; k++){
			  for(ConditionsArrayType::ptr_iterator it=LinkingConditions[k].ptr_begin(); it!= LinkingConditions[k].ptr_end(); ++it ){
			      (*it.base())->SetId(rId);
			      rId++;
			  }
		    }
	           
		   
		    unsigned int size = 0;
		    for(int k=0; k<number_of_threads; k++){
			  size+=LinkingConditions[k].size(); 
			  for(ConditionsArrayType::ptr_iterator it=LinkingConditions[k].ptr_begin(); it!= LinkingConditions[k].ptr_end(); ++it )
			      mr_model_part.Conditions().push_back(*it);
		    }


		    std::cout<<"        Number of linking conditions           = " << size               <<  std::endl; 
		    std::cout<<"        Total number of conditions             = " << rConditions.size() <<  std::endl; 

		    LinkingConditions.clear();
		KRATOS_CATCH("")
	 }

			  
//*****************************************************************************************************
//*****************************************************************************************************

///Permite decidir si el nodo que no esta dentro de un elemento es el correcto para ser el slave
void VerifyCorrectSlaveNode(unsigned int&  master,
	                    unsigned int&  slave,   
			    const array_1d<NodePointerType, 2 >&  Ids)
	   
{
   const unsigned int master_aux          = master;
   const unsigned int slave_aux           = slave;
   
   array_1d<double, 2>               Distances; 
   array_1d<double, 2>               Points0;
   array_1d<double, 2>               Points1;
   array_1d<double, 2>               Points2;
   
   WeakPointerVector<Condition>& neighb_cond_slave    = (Ids[slave_aux])->GetValue(NEIGHBOUR_CONDITIONS);
   WeakPointerVector<Condition>& neighb_cond_master   = (Ids[master_aux])->GetValue(NEIGHBOUR_CONDITIONS);
   
   Segment2D rSegment; 
   Points0[0]              =  (Ids[slave_aux])->X(); 
   Points0[1]              =  (Ids[slave_aux])->Y();
   double distance         =  DBL_MAX;
   double compare_distance =  0.00;
   for(WeakPointerVector<Condition>::iterator neighb = neighb_cond_master.begin();  neighb!= neighb_cond_master.end(); neighb++)
    {
      Condition::GeometryType& geom_2 = (neighb)->GetGeometry(); 
      Points1[0]                      = geom_2[0].X(); 
      Points1[1]                      = geom_2[0].Y();
      Points2[0]                      = geom_2[1].X();
      Points2[1]                      = geom_2[1].Y(); 
    
      rSegment.AssignPointsAndComputeParameters(Points1, Points2);
      compare_distance = rSegment.DistPoint2Segment2D(Points0);
      if(compare_distance<distance)
       {
        distance = compare_distance;
       }
    } 
    
    Distances[0] = distance;
    
    
   Points0[0]       =  (Ids[master_aux])->X(); 
   Points0[1]       =  (Ids[master_aux])->Y();
   distance         =  DBL_MAX;
   compare_distance =  0.00;
   for(WeakPointerVector<Condition>::iterator neighb = neighb_cond_slave.begin();  neighb!= neighb_cond_slave.end(); neighb++)
    {
      Condition::GeometryType& geom_2 = (neighb)->GetGeometry(); 
      Points1[0]                      = geom_2[0].X(); 
      Points1[1]                      = geom_2[0].Y();
      Points2[0]                      = geom_2[1].X();
      Points2[1]                      = geom_2[1].Y(); 
    
      rSegment.AssignPointsAndComputeParameters(Points1, Points2);
      compare_distance = rSegment.DistPoint2Segment2D(Points0);
      if(compare_distance<distance)
       {
        distance = compare_distance;
       }
    } 
    
    Distances[1] = distance;

    if( Distances[1]< Distances[0])
    {
      master =  slave_aux;
      slave  =  master_aux;
    }
  
  
}



//*****************************************************************************************************
//*****************************************************************************************************
void CreatePointLinkingConditions(
	 const unsigned int& master,
         const unsigned int& slave, 
	 const array_1d<NodePointerType, 2 >& Ids,
	 const ContactPairType& it_pair,
	 const PropertiesType::Pointer& tempProperties, 
	 const unsigned int& Id,
	 ConditionsArrayType& LinkingConditions
	 )
	
	{
	  KRATOS_TRY	  
	  // Slave Node
	  Point2D<Node<3> >::Pointer point_geom_slave    =  Point2D<Node<3> >::Pointer( new Point2D<Node<3> >(Ids[slave]) );
	  Condition::Pointer SlaveNode                   =  Condition::Pointer(new SlaveContactPointType(Id, point_geom_slave) );  
	  // Master Node
	  Point2D<Node<3> >::Pointer point_geom_master   =  Point2D<Node<3> >::Pointer( new Point2D<Node<3> >(Ids[master]));
	  Condition::Pointer MasterNode                  =  Condition::Pointer(new MasterContactPointType(Id, point_geom_master) );

	  Condition::GeometryType& Mgeom                 =  MasterNode->GetGeometry();
	  Condition::GeometryType& Sgeom                 =  SlaveNode ->GetGeometry();  
	  Line2D2<Node<3> >::Pointer Lgeom               =  Line2D2<Node<3> >::Pointer( new Line2D2<Node<3> >(Sgeom(0), Mgeom(0) ) );
	  Condition::Pointer newLink                     =  Condition::Pointer( new PointPointContactLink(Id,
	  Lgeom,
	  tempProperties, 
	  SlaveNode,
	  MasterNode
	  ) );
	  
	  LinkingConditions.push_back(newLink);    
	  KRATOS_CATCH("")
	  
	}
	
	
//*****************************************************************************************************
//*****************************************************************************************************
	 
//      cuando dos objetos intersectan pero no sabes que nodo cae dentro
	void CheckNearNodes( 
	   const unsigned int&  master,
	   const unsigned int&  slave,   
	   const PointerType&   SlaveObject,
	   const PointerType&   MasterObject,
	   std::vector<array_1d<NodePointerType, 2 > >&  Ids,
	   std::vector<Near_Node>& Is_Near
	 )
	 {   
	   
	   KRATOS_TRY
	     
	    std::vector<double>           Distance;
	    std::vector<double>           Distance_aux;
	    std::vector<double>::iterator it;
	    std::vector<double>::iterator it_2;
	    array_1d<NodePointerType, 2 > Id;
	    array_1d<double, 3>        vector;
	    const Element::GeometryType& geom_0    =  MasterObject->GetGeometry();
	    const Element::GeometryType& geom_1    =  SlaveObject->GetGeometry();
	    double distance                        =  0.00;
	       
	    array_1d<unsigned int, 9 > M;
	    array_1d<unsigned int, 9 > S;
	    
	    M[0] = 0;  M[1] = 0;  M[2] = 0;
	    M[3] = 1;  M[4] = 1;  M[5] = 1;
	    M[6] = 2;  M[7] = 2;  M[8] = 2;
	    
            S[0] = 0;  S[1] = 1;  S[2] = 2;
	    S[3] = 0;  S[4] = 1;  S[5] = 2;
	    S[6] = 0;  S[7] = 1;  S[8] = 2;
	      
	    
	    // busco la distancia menor
	    for(unsigned int i = 0; i<geom_0.size(); i++){
	       for(unsigned int j = 0; j<geom_1.size(); j++){
	          noalias(vector) = ( geom_0[i]-geom_1[j]) ;
	          distance = norm_2(vector);
	          Distance.push_back(distance);
	       }
	    }
	    
	    const double min   = (*std::min_element(Distance.begin(), Distance.end() ) );
	    it                 = std::find(Distance.begin(), Distance.end(), min);
	    const int position = int(it-Distance.begin());
	    Id[master]         = geom_0(M[position]);   
	    Id[slave]          = geom_1(S[position]);
	    Ids.push_back(Id); 
	    Is_Near.push_back(no_near);
	    
	    /*
	    //Check si dos corner chocan
	    std::vector<NodePointerType> nodes;
	    const bool test_one = VerifyToCornerIntersect(nodes, SlaveObject, MasterObject);  
	    if( test_one==false &&  nodes.size()!=0) 
	    {   
	        KRATOS_WATCH("BBBBBBBBBBB")
	        if(nodes.size()==2)
		{
		Id[master]         = nodes[0]; 
		Id[slave]          = nodes[1]; 
		}
		else
		{
		 // si no se cumple lo anterior tomamos los nodos mas cercanos
		 // WARNING = Solo valido para un caso en que un solo nodo quede fuera
		 const double min   = (*std::min_element(Distance.begin(), Distance.end() ) );
		 it                 = std::find(Distance.begin(), Distance.end(), min);
		 const int position = int(it-Distance.begin());
		 Id[master]         = geom_0(M[position]);   
		 Id[slave]          = geom_1(S[position]);
		}
		
		Ids.push_back(Id); 
		Is_Near.push_back(no_near);
	     }
	     
	     
	    //NO VALIDO PARA ELEMTOS CON MAL RATIO
	    else
	    {

	        Distance_aux.resize(Distance.size());
		Distance_aux = Distance;
	        const double min   = (*std::min_element(Distance.begin(), Distance.end() ) );
		it                 = std::find(Distance.begin(), Distance.end(), min);
		const int position = int(it-Distance.begin()); 
		Id[master]         = geom_0(M[position]);   
		Id[slave]          = geom_1(S[position]);
		Ids.push_back(Id); 
		Is_Near.push_back(no_near);
		
		
		
		const double min_2   = (*min_element_2(Distance_aux.begin(), Distance_aux.end(), min ) );
		it_2                 =  std::find(Distance.begin(), Distance.end(), min_2);
		const int position_2 = int(it_2-Distance.begin()); 
		Id[master]           = geom_0(M[position_2]);   
		Id[slave]            = geom_1(S[position_2]);
		Ids.push_back(Id); 
		Is_Near.push_back(no_near);
	    }
	    
	    */
	    KRATOS_CATCH("")
	 }


//*****************************************************************************************************
//*****************************************************************************************************
 // Saca elsegundo min de un vector
 std::vector<double>::iterator  min_element_2( const std::vector<double>::iterator first,  const std::vector<double>::iterator last, const double& cond)
{
  std::vector<double>::iterator second_lowest = first;
  std::vector<double>::iterator first_1       = first;
  std::vector<double>::iterator first_2       = first;
  const int size  = int(last- first); 
  int count = 0;
  if (first==last) return last;
  for(first_1=first; first_1!=last; first_1++){
       if(*first_1!=cond) 
          for(first_2=first; first_2!=last; first_2++){   
	        if(*first_2>cond  && *first_2!=*first_1){
 	           if(*first_1<*first_2){
		     count++;
 		     continue;
		   }
 		   else
		     break;
		}
	}
	 if(count==size-2) {
	   *second_lowest = *first_1;
	   break;
	 }
	   else
	     count=0;
	 } 
	 
  return second_lowest;
}

//*****************************************************************************************************
//*****************************************************************************************************

    bool VerifyToCornerIntersect(  std::vector<NodePointerType>& Ids,
                                   const PointerType&    SlaveObject,
	                           const PointerType&    MasterObject
                                ) 

{

  KRATOS_TRY   
  
  WeakPointerVector<Condition>& neighb_cond_master  =  MasterObject->GetValue(NEIGHBOUR_CONDITIONS);
  WeakPointerVector<Condition>& neighb_cond_slave   =  SlaveObject->GetValue(NEIGHBOUR_CONDITIONS);
 
  std::vector<std::vector<unsigned int> >   segment;
  segment.resize(neighb_cond_slave.size());
  
  vector<array_1d<double, 2> >        Points0;
  vector<array_1d<double, 2> >        Points1;
  array_1d<double, 2>                 Point;

  unsigned int I   = 0;
  unsigned int II  = 1;
  unsigned int III = 1;
  unsigned int IV  = 0; 
  
  Points0.resize(2, false); 
  Points1.resize(2, false);
  
  for(WeakPointerVector<Condition>::iterator cond_slave  = neighb_cond_slave.begin(); cond_slave!= neighb_cond_slave.end(); ++cond_slave)
  {
      Condition::GeometryType& geom = cond_slave->GetGeometry();
      
      Point[0]      = 0.00;
      Point[1]      = 0.00;
      Points0(0)[0] = geom[0].X(); 
      Points0(0)[1] = geom[0].Y();
      Points0(1)[0] = geom[1].X();
      Points0(1)[1] = geom[1].Y(); 

      I    = 0;
      III  = 0; 
      for(WeakPointerVector< Condition >::iterator cond  = neighb_cond_master.begin(); cond!= neighb_cond_master.end(); ++cond){
        Condition::GeometryType& geom_3 = cond->GetGeometry();     
	Points1(0)[0] = geom_3[0].X(); 
	Points1(0)[1] = geom_3[0].Y();
	Points1(1)[0] = geom_3[1].X();
	Points1(1)[1] = geom_3[1].Y();

	if(IntersectionSegments::IntersectSegment(Point, Points0, Points1)!=IT_EMPTY)   
	      segment[IV].push_back(I);
	
        I++;
        III++;
        if(III>neighb_cond_master.size())
            break;
     }     

     II++;
     IV++;
     if(II>neighb_cond_slave.size())
         break;
  }  

        
     /// dos triangulos que se intersectan pero no tienen nodo dentro.
     /// Sus aristan chocan en dos partes del triangulo master
     if(segment.size()==3) 
	if(segment[0].size()== 2 && segment[1].size()== 2 && segment[2].size()== 2)
	  return true;    
        
    if(neighb_cond_master.size()==2 && neighb_cond_slave.size()==2)
    {
        
        Condition::GeometryType& geom_1 = (neighb_cond_master(0).lock())->GetGeometry();
	Condition::GeometryType& geom_2 = (neighb_cond_master(1).lock())->GetGeometry();
	
	if(geom_1[0].Id()==geom_2[0].Id())
	  Ids.push_back(geom_1(0));
	else if(geom_1[0].Id()==geom_2[1].Id())
	  Ids.push_back(geom_1(0));
	else if(geom_1[1].Id()==geom_2[0].Id())
	  Ids.push_back(geom_1(1));
        else if(geom_1[1].Id()==geom_2[1].Id())
	  Ids.push_back(geom_1(1));
	else
	  std::cout<< "No node A " << std::endl; 
	
	Condition::GeometryType& geom_3 = (neighb_cond_slave(0).lock())->GetGeometry();
	Condition::GeometryType& geom_4 = (neighb_cond_slave(1).lock())->GetGeometry();
	
	if(geom_3[0].Id()==geom_4[0].Id())
	  Ids.push_back(geom_3(0));
	else if(geom_3[0].Id()==geom_4[1].Id())
	  Ids.push_back(geom_3(0));
        else if(geom_3[1].Id()==geom_4[0].Id())
	  Ids.push_back(geom_3(1));
	else if(geom_3[1].Id()==geom_4[1].Id())
	  Ids.push_back(geom_3(1));
	else
	  std::cout<< "No node B " << std::endl; 

	if(Ids.size()==2)
	  return false;
      }
       
    return true;

   KRATOS_CATCH("")
       
     }


//*****************************************************************************************************
//*****************************************************************************************************
	 
	 Near_Node CheckNearNodes( 
	   const unsigned int&  master,
	   const unsigned int&  slave,   
	   const NodePointerType& SlaveNode,
	   const PointerType& MasterObject,
	   array_1d<NodePointerType, 2 >&  Ids
	 )
	 
	 {
	    //std::vector<double>        Distance;
	    array_1d<double, 3>        vector;
	    array_1d<double, 3>        coordinates =  SlaveNode->Coordinates();
	    const Element::GeometryType& geom_0    =  MasterObject->GetGeometry();
	    double distance                        = 0.00;
	    double distance2                       = 1E10;;

	    
	    // busco la distancia menor
	    for(unsigned int i = 0; i<geom_0.size(); i++){
	           noalias(vector) = ( geom_0(i)->Coordinates() -  coordinates);
	           distance  = norm_2(vector);
		   //Distance.push_back(distance);
	           if(distance<distance2){
	               distance2    =  distance;
	               Ids[master]  =  geom_0(i);
	            }
	       
	    }
	    
	    Ids[slave]   = SlaveNode;
	  
    
//             double max   = (*std::max_element(Distance.begin(), Distance.end() ) );  
// 	    double min   = (*std::min_element(Distance.begin(), Distance.end() ) );
// 	    double ratio = std::fabs(min/max); 
// 	    if(ratio < 1E-8) 
//                 return yes_near;
// 	     
 	    return no_near; 
	    
	 }
	 
	 
         bool  CreateLinkingConditions(
	 const unsigned int& Id,
	 const unsigned int& master,
         const unsigned int& slave,
         const array_1d<NodePointerType, 2 >& Ids,   
         const ContactPairType& it_pair,
	 const PropertiesType::Pointer& tempProperties, 
	 Exist_Node& Exist,
	 ConditionsArrayType& LinkingConditions
	 )
	 
	 {
	   if(mrdimension==2)
	     return CreateLinkingConditions2D(Id, master, slave, Ids, it_pair, tempProperties, Exist, LinkingConditions);
	   else
	     CreateLinkingConditions3D(Id, master, slave, Ids, it_pair, tempProperties, Exist, LinkingConditions);
	   
	   return false;
	   
	 }
	 
//*****************************************************************************************************
//*****************************************************************************************************

	 bool CreateLinkingConditions2D(
	 const unsigned int& Id,
	 const unsigned int& master,
         const unsigned int& slave,
         const array_1d<NodePointerType, 2 >& Ids, 
         const ContactPairType& it_pair,
	 const PropertiesType::Pointer& tempProperties, 
	 Exist_Node& Exist,
	 ConditionsArrayType& LinkingConditions
	 )
	 
	 {
	    KRATOS_TRY 
	   
	    Condition::Pointer MasterFace;
	    array_1d<double,3 > Normal_r;
	    array_1d<double,3 > GL;
	    ProcessInfo& CurrentProcessInfo  = mr_model_part.GetProcessInfo(); 
	    
	   if(Exist==yes_nodes){
	   const  bool exist_segment        =  LocateMasterSegment( Ids[slave], 
							    Ids[master], 
							    (it_pair)[master],
							    (it_pair)[slave],
							    MasterFace,
							    Exist);
	   if(exist_segment==true){
	   const double zero                   =  1.00E-6;  
	   Condition::GeometryType& geom       =  MasterFace->GetGeometry();                  
	   array_1d<double, 3>& point_slave    =  Ids[slave]->Coordinates();  
	   array_1d<double, 3>& point_left     =  geom.GetPoint(0);
	   array_1d<double,3>  seg             =  geom.GetPoint(0)-geom.GetPoint(1);
	   noalias(GL)                         =  point_slave - point_left; 
	   MasterFace->Calculate(NORMAL, Normal_r, CurrentProcessInfo);   
	   //const double distance = norm_2(seg);
	   //const double gat      = inner_prod(GL, (1.00/distance)*seg) ;
	   const double gap = inner_prod(GL, Normal_r); 
	   bool is_repited  = bool(Ids[slave]->GetValue(IS_CONTACT_SLAVE)==0);
	   if(gap<zero && is_repited==true) // && gat>zero && gat<distance)
	    {   
	        Ids[slave]->GetValue(IS_CONTACT_SLAVE)    =  1;
	        Point2D<Node<3> >::Pointer point_geom     =  Point2D<Node<3> >::Pointer( new Point2D<Node<3> >(Ids[slave]));
		Condition::Pointer SlaveNode              =  Condition::Pointer(new SlaveContactPointType(Id, point_geom) );  
		Condition::GeometryType& Mgeom            =  MasterFace->GetGeometry();
		Condition::GeometryType& Sgeom            =  SlaveNode->GetGeometry();
		//std::cout<<"Master = "<< (it_pair)[master]->Id() <<"   Slave = " << (it_pair)[slave]->Id() <<std::endl; 
		//std::cout<<"        Node (Y) =  " << Ids[slave]->Id() << " Master Face = " << MasterFace->Id() << std::endl;
		Triangle2D3<Node<3> >::Pointer Lgeom      =  Triangle2D3<Node<3> >::Pointer( new Triangle2D3<Node<3> >( Sgeom(0), Mgeom(0), Mgeom(1) ) );
		Condition::Pointer newLink                =  Condition::Pointer( new PointSegmentContactLink(Id,
													     Lgeom,
													     tempProperties,
													     MasterFace, 
													     SlaveNode));
	
		LinkingConditions.push_back( newLink );   
		return exist_segment;
	         }
	       }
	     }
	     
	     else if(Exist==no_nodes)
	      {
		//KRATOS_WATCH(Ids[slave]->Id()) 
	        Point2D<Node<3> >::Pointer point_geom     =  Point2D<Node<3> >::Pointer( new Point2D<Node<3> >(Ids[slave]));
		Condition::Pointer SlaveNode              =  Condition::Pointer(new SlaveContactPointType(Id, point_geom) ); 
		//KRATOS_WATCH((it_pair)[master]->Id())
		//KRATOS_WATCH((it_pair)[master]->GetValue(NEIGHBOUR_CONDITIONS).size())
		MasterFace                                = ((it_pair)[master]->GetValue(NEIGHBOUR_CONDITIONS)(0)).lock();
		//std::cout<<"        Master   =  " << (it_pair)[master]->Id() <<"   Slave = " << (it_pair)[slave]->Id() <<std::endl; 
		//std::cout<<"        Node (N) =  " << Ids[slave]->Id() << " Master Face = " << MasterFace->Id() << std::endl;
		Condition::GeometryType& Mgeom            =  MasterFace->GetGeometry();
		Condition::GeometryType& Sgeom            =  SlaveNode->GetGeometry();   
		Triangle2D3<Node<3> >::Pointer Lgeom      =  Triangle2D3<Node<3> >::Pointer( new Triangle2D3<Node<3> >( Sgeom(0), Mgeom(0), Mgeom(1) ) );
		Condition::Pointer newLink                =  Condition::Pointer( new PointSegmentContactLink(Id,
													     Lgeom,
													     tempProperties,
													     MasterFace, 
													     SlaveNode));
		LinkingConditions.push_back( newLink );   
		return true;
	      }
	      else
		return false;
	      
	      return false;
	      KRATOS_CATCH("")
	    }   
	 
 
//*****************************************************************************************************
//*****************************************************************************************************

	 void CreateLinkingConditions3D(
	 const unsigned int& Id,
	 const unsigned int& master,
         const unsigned int& slave,
         const array_1d<NodePointerType, 2 >& Ids, 
         const ContactPairType& it_pair,
	 const PropertiesType::Pointer& tempProperties, 
	 Exist_Node& Exist,
	 ConditionsArrayType& LinkingConditions
	 )
	 
	 {
	    KRATOS_TRY 
	    
	    //bool is_repited  = bool(Ids[slave]->GetValue(IS_CONTACT_SLAVE)==0)
	    //if(is_repited==true)
	    {
	    Point<3> MasterContactLocalPoint;
	    Point<3> SlaveContactLocalPoint;
	    int SlaveIntegrationPointIndex = 0;
	    Condition::Pointer MasterFace             =  (Ids[slave])->GetValue(CONTACT_LINK_MASTER);
	    Point3D<Node<3> >::Pointer point_geom     =  Point3D<Node<3> >::Pointer( new Point3D<Node<3> >(Ids[slave]));
	    Condition::Pointer SlaveNode              =  Condition::Pointer(new SlaveContactPoint3D(Id, point_geom) );
	    Condition::GeometryType& Mgeom            =  MasterFace->GetGeometry();
	    
	    
	    Condition::GeometryType& Sgeom            =  SlaveNode->GetGeometry();   
	    Tetrahedra3D4<Node<3> >::Pointer Lgeom    =  Tetrahedra3D4<Node<3> >::Pointer( new Tetrahedra3D4<Node<3> >( Sgeom(0), Mgeom(0), Mgeom(1), Mgeom(2) ) );
	   
	    Condition::Pointer newLink                =  Condition::Pointer( new ContactLink3DExplicit(
									    Id,
									    Lgeom,
									    tempProperties,
									    MasterFace, 
									    SlaveNode,
									    MasterContactLocalPoint,
									    SlaveContactLocalPoint,
									    SlaveIntegrationPointIndex
									    ));
		
    
	      LinkingConditions.push_back( newLink );    
	     }
	    
	     KRATOS_CATCH("")
	  }   


//*****************************************************************************************************
//*****************************************************************************************************

bool LocateMasterSegment( const NodePointerType& SlaveNode,
			  const NodePointerType& MasterNode,  //the most near
			  const PointerType& MasterObject,
			  const PointerType& SlaveObject,
			  Condition::Pointer& MasterFace,
			  Exist_Node& Exist
			 )
 	  {   
	    KRATOS_TRY
	   

	     WeakPointerVector<Condition>& neighb_cond = MasterObject->GetValue(NEIGHBOUR_CONDITIONS);     
	     if(neighb_cond.size()!=0){
	     switch(Exist)
	     {
	       case(yes_nodes):
	         {  
	          if(neighb_cond.size()==1){
		    
		     Condition::Pointer rCond_2;
		     Condition::Pointer& rCond_1 =  SlaveNode->GetValue(CONTACT_LINK_MASTER);
		     
		     //problema. Los Elementos internos no tienen master face
		     //KRATOS_WATCH(SlaveNode->Id())
		     //KRATOS_WATCH(rCond_1->Id())
		     const bool test_2 = Test_Four(SlaveNode, MasterObject, rCond_2);
		     //KRATOS_WATCH(MasterObject->Id())
		     //KRATOS_WATCH(SlaveNode->Id())
		     //KRATOS_WATCH(rCond_1->Id())
		     const bool test_3 = bool(neighb_cond(0).lock()->Id()==rCond_1->Id());
		     //KRATOS_WATCH("-----------------")
		     
		     
		     if(test_2==true){
		        if(rCond_1->Id()==rCond_2->Id())
			    MasterFace = rCond_1;
			else
			    MasterFace = rCond_2; ///WARNING 
		     }
		     else if( test_3==true)
		     {
		       MasterFace = rCond_1;
		     }
		     else
		       if(neighb_cond.size()!=0)
		          MasterFace = neighb_cond(0).lock();
		     
		     return true; 
		  }

	         if(neighb_cond.size()>=2){
                     Condition::Pointer& rCond_1 =  SlaveNode->GetValue(CONTACT_LINK_MASTER);
		    
		     const bool is_corner   = Is_Corner(SlaveNode, MasterObject);
		     if(is_corner==true)
		     { 
		       Test_One_B_Distances(SlaveNode, MasterObject, rCond_1);
		       MasterFace = rCond_1; 
		       return true;
		     }

		     const bool node_corner =   Is_Node_Corner(SlaveNode, SlaveObject);
		     if(node_corner==true)
		     { 
		        Condition::Pointer  rCond_2;
		        Condition::Pointer  rCond_3;
		        Condition::Pointer  rCond_4;
		        const bool test_1      = Test_One_C (SlaveNode, MasterNode,   rCond_2); //true==1
			if(test_1==true){
			    if(rCond_2->Id()==rCond_1->Id())
			      MasterFace = rCond_1;
			    else
			      MasterFace = rCond_2; ///WARNING
			  return true;
			}
		       //if(SlaveNode->Id()==29)
		       //std::cout<<"Test 1 " <<std::endl;
		       
			const bool test_2      = Test_Three (SlaveNode, MasterObject, rCond_3); 
			if(test_2==true){
				if(rCond_3->Id()==rCond_1->Id())
				      MasterFace = rCond_1;
				else
				  MasterFace = rCond_3; 
				return true;
				}
				
			//if(SlaveNode->Id()==29)
		        //std::cout<<"Test 2 " <<std::endl;
		       /*
		       if(SlaveNode->Id()==51)
		          std::cout<<"Test 3 " <<std::endl;
		       const bool test_3      = Test_Five  (SlaveNode, MasterObject, rCond_4);	
		       if(test_3==true){
				if(rCond_4->Id()==rCond_1->Id())
				      MasterFace = rCond_1;
				else
				  MasterFace = rCond_4; ///WARNING rCond_3;
				return true;
			    }
			    */
			}
		     
		     //std::cout<<" No test " <<std::endl;
		     Condition::Pointer  rCond_2;
		     const bool test_4 = Test_One_B_Distances(SlaveNode, MasterObject,  rCond_2); 
		     if(test_4==true){
		         if(rCond_2->Id()==rCond_1->Id())
			       MasterFace = rCond_1;
			     else
			        MasterFace = rCond_2; 
			return true;
			}
			
		     MasterFace = rCond_1;
		     return true; 
		     }
		     
		  break;
	       }
	       
	      case(no_nodes):
	       { 
		
		//Condition::Pointer& rCond_1 =  SlaveNode->GetValue(CONTACT_LINK_MASTER);
		//Condition::Pointer  rCond_2;
		//Condition::Pointer  rCond_3;
		//Condition::Pointer  rCond_4;
		
		//const bool test_1      = Test_One_C(SlaveNode, MasterNode,  rCond_2);
		//const bool test_2      = Test_Two(SlaveNode, MasterObject, MasterFace)
		//const bool test_3      = Test_Four(SlaveNode, MasterObject, MasterFace);
		//const bool test_4      = Test_One_A(SlaveNode, MasterNode, MasterFace);

// 		if(test_1==true){
// 		     if(rCond_2->Id()==rCond_1->Id())
// 		        MasterFace = rCond_1;
// 		     else
// 		        MasterFace = rCond_2;
// 		  return true;
// 		}
		
		if(Test_One_C(SlaveNode, MasterNode, MasterFace)){
		   return true;
		}
		  
		if(Test_Two(SlaveNode, MasterObject, MasterFace)){
		    return true;
		}
	         
	        if(Test_Four(SlaveNode, MasterObject, MasterFace) ){
	           return true;
		}
		
		if(Test_One_A(SlaveNode, MasterNode, MasterFace)){
		    return true;
		}
               
                MasterFace = SlaveNode->GetValue(CONTACT_LINK_MASTER);
                return true; 
		break;
	       }  
	     }
	     
              return false;
	     }
	     
	    MasterFace = SlaveNode->GetValue(CONTACT_LINK_MASTER);  // no estaba: Si no funciona comenatr y poner return false.
	    return true; // false
	    KRATOS_CATCH("")
   }
   
   
// Si el nodo esta dentro de elemento   
bool Is_Corner(
               const NodePointerType& SlaveNode,
	       const PointerType& MasterObject)

{

  KRATOS_TRY   
  
  WeakPointerVector<Condition>& neighb_cond_master  = MasterObject->GetValue(NEIGHBOUR_CONDITIONS);
  WeakPointerVector<Condition>& neighb_cond_slave   = SlaveNode->GetValue(NEIGHBOUR_CONDITIONS);
 
  std::vector<unsigned int>           segment;
  std::vector<unsigned int>::iterator it;
  vector<array_1d<double, 2> >        Points0;
  vector<array_1d<double, 2> >        Points1;
  array_1d<double, 2>                 Point;

  unsigned int I   = 0;
  unsigned int II  = 1;
  unsigned int III = 1;

    
   Points0.resize(2, false); 
   Points1.resize(2, false);
  
  for(WeakPointerVector<Condition>::iterator cond_slave  = neighb_cond_slave.begin(); cond_slave!= neighb_cond_slave.end(); ++cond_slave)
  {
      Condition::GeometryType& geom = cond_slave->GetGeometry();
      Point[0]      = 0.00;
      Point[1]      = 0.00;
      
      Points0(0)[0] = geom[0].X(); 
      Points0(0)[1] = geom[0].Y();
      Points0(1)[0] = geom[1].X();
      Points0(1)[1] = geom[1].Y(); 

      I   = 0; 
      III = 1;
      for(WeakPointerVector< Condition >::iterator cond  = neighb_cond_master.begin(); cond!= neighb_cond_master.end(); ++cond){
        Condition::GeometryType& geom_3 = cond->GetGeometry();   
	Points1(0)[0] = geom_3[0].X(); 
	Points1(0)[1] = geom_3[0].Y();
	Points1(1)[0] = geom_3[1].X();
	Points1(1)[1] = geom_3[1].Y();

	if(IntersectionSegments::IntersectSegment(Point, Points0, Points1)!=IT_EMPTY) 
	{  
	   if(segment.size()==0){
	         segment.push_back(I);}
	   else
	   {
	   it = std::find(segment.begin(), segment.end(), I);
	   if(it==segment.end())
	      segment.push_back(I); 
	   }
	}
        I++;
        III++;
        if(III>neighb_cond_master.size())
            break;
     }     
     
     II++;
     if(II>neighb_cond_slave.size())
         break;
  }  

   if(segment.size()>=2)
      return true;
   
   return false;

   KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************ 


///verifica si un nodo es corner
bool Is_Node_Corner(const NodePointerType& SlaveNode,
		    const PointerType& SlaveObject
		   )
 {
     KRATOS_TRY   
     WeakPointerVector<Condition>& neighb_cond_node    = SlaveNode->GetValue(NEIGHBOUR_CONDITIONS);
     WeakPointerVector<Condition>& neighb_cond_object  = SlaveObject->GetValue(NEIGHBOUR_CONDITIONS);
     
     unsigned int count = 0;
     ///WARNING = Verificar para tres lados 
     if(neighb_cond_object.size()>=2)
	  for(unsigned int i = 0; i<neighb_cond_node.size(); i++)
	     for(unsigned int j = 0; j<neighb_cond_object.size(); j++)
	       if((neighb_cond_node(i).lock())->Id()==(neighb_cond_object(j).lock())->Id())
		   count++;

     if(count==2)
        return true;
     
     
     return false;
     
     KRATOS_CATCH("")
   
 }

//************************************************************************************
//************************************************************************************ 



/// comparacion con desplazamientos
bool Test_One_A( const NodePointerType& SlaveNode,
	         const NodePointerType& MasterNode,
		 Condition::Pointer& rCond) 
{
  
  KRATOS_TRY   
  WeakPointerVector<Condition>& neighb_cond_master = MasterNode->GetValue(NEIGHBOUR_CONDITIONS);
  
  if(neighb_cond_master.size()!=0)
  {
  std::vector<unsigned int> segment;
  unsigned int I        = 0;
  unsigned int segmento = 0;

  std::vector<array_1d<double, 2> > Points;   // punto de interseccion del segmento 
  vector<array_1d<double, 2> >      Points0;
  vector<array_1d<double, 2> >      Points1;
  array_1d<double, 2>               Point;
     
  array_1d<double,3>& old_pos     = SlaveNode->FastGetSolutionStepValue(DISPLACEMENT,3);   


  Points0.resize(2, false); 
  Points1.resize(2, false);

  Points0(0)[0] = SlaveNode->X0() + old_pos[0];
  Points0(0)[1] = SlaveNode->Y0() + old_pos[1];
  Points0(1)[0] = SlaveNode->X(); 
  Points0(1)[1] = SlaveNode->Y(); 

  unsigned int JJ = 1; 
  for(WeakPointerVector< Condition >::iterator cond  = neighb_cond_master.begin(); cond!= neighb_cond_master.end(); cond++){
  Condition::GeometryType& geom_2 = cond->GetGeometry();
	    
    Points1(0)[0] = geom_2[0].X(); 
    Points1(0)[1] = geom_2[0].Y();
    Points1(1)[0] = geom_2[1].X();
    Points1(1)[1] = geom_2[1].Y();

    if(IntersectionSegments::IntersectSegment(Point, Points0, Points1)==IT_POINT){
        Points.push_back(Point); 
        segment.push_back(I);
    }
    I++;
    JJ++;
    if(JJ>neighb_cond_master.size())
       break;
  }

  if (Points.size()!=0)
  {
      if (Points.size()==1){
      segmento = segment[0];
      }

      else if (Points.size()>1)
      {
      
      double dist0 = 0;
      array_1d<double, 2> rect;
      std::vector<double> Distance;
      std::vector<double>::iterator it;
      int position  = 0; 
      
      for(unsigned int i = 0; i<Points.size(); i++){
         rect  = Points0[1] - Points[i];
         dist0 = std::sqrt(inner_prod(rect, rect ));
         Distance.push_back(dist0);
      }
     
      const double max   = (*std::max_element(Distance.begin(), Distance.end() ) );
      it                 = std::find(Distance.begin(), Distance.end(), max);
      position           = int(it-Distance.begin()); 
      segmento           = segment[position];
      }
      
     rCond = neighb_cond_master(segmento).lock(); 
     return true;
  }
  }
  return false;
  KRATOS_CATCH("")
}



bool Test_One_C( const NodePointerType& SlaveNode,
	         const NodePointerType& MasterNode,
		 Condition::Pointer& rCond) 
{
  
  KRATOS_TRY
  WeakPointerVector<Condition>& neighb_cond            = MasterNode->GetValue(NEIGHBOUR_CONDITIONS);
  WeakPointerVector<Condition>& neighb_cond_slave      = SlaveNode->GetValue(NEIGHBOUR_CONDITIONS);
  
  
  
  if(neighb_cond.size()!=0 && neighb_cond_slave.size()!=0)
  {  
  std::vector<unsigned int>         segment;
  std::vector<array_1d<double, 2> > Points; // punto de interseccion del segmento
  vector<array_1d<double, 2> >      Points0;
  vector<array_1d<double, 2> >      Points1;
  array_1d<double, 2>               Point;

  // test with edges
  unsigned int segmento = 0;
  unsigned int I        = 0;
  unsigned int II       = 1;
  unsigned int III      = 1;

   
  Points0.resize(2, false); 
  Points1.resize(2, false);
  
  for(WeakPointerVector<Condition>::iterator cond_slave  = neighb_cond_slave.begin(); cond_slave!= neighb_cond.end(); ++cond_slave)
  {
      Condition::GeometryType& geom = cond_slave->GetGeometry();
      Point[0]      = 0.00;
      Point[1]      = 0.00;
      
      Points0(0)[0] = geom[0].X(); 
      Points0(0)[1] = geom[0].Y();
      Points0(1)[0] = geom[1].X();
      Points0(1)[1] = geom[1].Y(); 

      I   = 0; 
      III = 1;
      for(WeakPointerVector< Condition >::iterator cond  = neighb_cond.begin(); cond!= neighb_cond.end(); ++cond){
        Condition::GeometryType& geom_3 = cond->GetGeometry();     
	Points1(0)[0] = geom_3[0].X(); 
	Points1(0)[1] = geom_3[0].Y();
	Points1(1)[0] = geom_3[1].X();
	Points1(1)[1] = geom_3[1].Y();
	
	if(IntersectionSegments::IntersectSegment(Point, Points0, Points1)==IT_POINT)
	{  
	  Points.push_back(Point); 
	  segment.push_back(I); 
	}
        I++;
        III++;
        if(III>neighb_cond.size())
            break;
     }     

     II++;
     if(II>neighb_cond_slave.size())
         break;
  }  

  if (Points.size()!=0)
  {
     if (Points.size()==1){
       segmento = segment[0];
     }
   // en caso de que el nodo quede fuera e intersecte con dos aristas 
  else if (Points.size()>1)
  {
      Points0(0)[0] = SlaveNode->X(); 
      Points0(0)[1] = SlaveNode->Y(); 
      Points0(1)[0] = SlaveNode->X();  
      Points0(1)[1] = SlaveNode->Y(); 

      double dist0 = 0.00;
      array_1d<double, 2> rect;
      std::vector<double> Distance;
      std::vector<double>::iterator it;
      int position  = 0; 
      
      for(unsigned int i = 0; i<Points.size(); i++){
         rect  = Points0[1] - Points[i];
         dist0 = std::sqrt(inner_prod(rect, rect ));
         Distance.push_back(dist0);
      }
     
      const double max   = (*std::max_element(Distance.begin(), Distance.end() ) );
      it                 = std::find(Distance.begin(), Distance.end(), max);
      position           = int(it-Distance.begin()); 
      segmento           = segment[position];
   
  }
  rCond = neighb_cond(segmento).lock();
  return true;
  }
  }
 
  return false;
  
  KRATOS_CATCH("")
}


//*****************************************************************************************************
//*****************************************************************************************************

/// Busca la interseccion  de la trayectoria del nodo esclavo
/// desde su tiempo actual hasta 3 paso atras con las aristas del
/// elemento master. La menor de todas es el segmento  de contacto.

bool Test_One_B(const NodePointerType& SlaveNode,
		const NodePointerType& MasterNode,
		Condition::Pointer& rCond) 
 {
  
  KRATOS_TRY   
  
  WeakPointerVector<Condition>& neighb_cond_master     = MasterNode->GetValue(NEIGHBOUR_CONDITIONS);
  if(neighb_cond_master.size()!=0 )
  {  
  std::vector<unsigned int> segment;
  unsigned int segmento = 0;
  unsigned int I        = 0;
  
     
  std::vector<array_1d<double, 2> > Points;   // punto de interseccion del segmento 
  vector<array_1d<double, 2> >      Points0;
  vector<array_1d<double, 2> >      Points1;
  array_1d<double, 2>               Point;
  
  array_1d<double,3>& old_pos     = SlaveNode->FastGetSolutionStepValue(DISPLACEMENT,3);   

  Points0.resize(2, false); 
  Points1.resize(2, false);

  Points0(0)[0] = SlaveNode->X0() + old_pos[0];
  Points0(0)[1] = SlaveNode->Y0() + old_pos[1];
  Points0(1)[0] = SlaveNode->X(); 
  Points0(1)[1] = SlaveNode->Y(); 
 /*
  array_1d<double, 2> Dir;
  Dir[0]       = Points0(1)[0] - Points0(0)[0];  
  Dir[1]       = Points0(1)[1] - Points0(0)[1];
  noalias(Dir) = (1.00/(sqrt(inner_prod(Dir, Dir)))) * Dir; 
  
  Points0(0)[0] -= Dir[0];
  Points0(0)[1] -= Dir[1];
  Points0(1)[0] += Dir[0]; 
  Points0(1)[1] += Dir[1]; 
 */ 
  
  unsigned int JJ = 1; 
  for(WeakPointerVector< Condition >::iterator cond  = neighb_cond_master.begin(); cond!= neighb_cond_master.end(); cond++){
  Condition::GeometryType& geom_2 = cond->GetGeometry();
	    
    Points1(0)[0] = geom_2[0].X(); 
    Points1(0)[1] = geom_2[0].Y();
    Points1(1)[0] = geom_2[1].X();
    Points1(1)[1] = geom_2[1].Y();
    
    if(IntersectionSegments::IntersectSegment(Point, Points0, Points1)!=IT_EMPTY){
        Points.push_back(Point); 
        segment.push_back(I);
    }
    I++;
    JJ++;
    if(JJ>neighb_cond_master.size())
       break;
  }

  if (Points.size()!=0)
  {
      if (Points.size()==1){
      segmento = segment[0];
      }

      else if (Points.size()>1)
      {

      double dist0 = 0.00;
      array_1d<double, 2> rect;
      std::vector<double> Distance;
      std::vector<double>::iterator it;
      int position  = 0; 
      
      for(unsigned int i = 0; i<Points.size(); i++){
         rect  = Points0[1] - Points[i];
         dist0 = std::sqrt(inner_prod(rect, rect ));
         Distance.push_back(dist0);
      }
     
      const double min   = (*std::min_element(Distance.begin(), Distance.end() ) );
      it                 = std::find(Distance.begin(), Distance.end(), min);
      position           = int(it-Distance.begin()); 
      segmento           = segment[position];
      }
      
     rCond = neighb_cond_master( segmento).lock(); 
     return true;
  }
  }
  return false;
  KRATOS_CATCH("")
}



//Calculando distancias de puntos a segmentos 
bool Test_One_B_Distances( const NodePointerType& SlaveNode,
	                   const PointerType& MasterObject,
	                   Condition::Pointer& rCond) 
{
  
  KRATOS_TRY   
  
  WeakPointerVector<Condition>& neighb_cond_master     = MasterObject->GetValue(NEIGHBOUR_CONDITIONS);

  std::vector<unsigned int> segment;
  unsigned int I        = 0;
  unsigned int segmento = 0;

     
  std::vector<double>               Distances;   // punto de interseccion del segmento 
  array_1d<double, 2>               Points0;
  vector<array_1d<double, 2> >      Points1;
    
  Points1.resize(2, false);

  Points0[0] = SlaveNode->X(); 
  Points0[1] = SlaveNode->Y();
  Segment2D Segment1;
  
  unsigned int JJ = 1; 
  for(WeakPointerVector< Condition >::iterator cond  = neighb_cond_master.begin(); cond!= neighb_cond_master.end(); cond++){
  Condition::GeometryType& geom_2 = cond->GetGeometry();
	    
    Points1(0)[0] = geom_2[0].X(); 
    Points1(0)[1] = geom_2[0].Y();
    Points1(1)[0] = geom_2[1].X();
    Points1(1)[1] = geom_2[1].Y();
    
    Segment1.AssignPointsAndComputeParameters(Points1[0], Points1[1]);
    Distances.push_back(Segment1.DistPoint2Segment2D(Points0));
    segment.push_back(I);
    
    I++;
    JJ++;
    if(JJ>neighb_cond_master.size())
       break;
  }

  
  if (Distances.size()!=0)
  {
      if (Distances.size()==1){
      segmento = segment[0];
      }

      else if (Distances.size()>1)
      {     
      std::vector<double>::iterator it;
      int position  = 0; 
      const double min   = (*std::min_element(Distances.begin(), Distances.end() ) );
      it                 = std::find(Distances.begin(), Distances.end(), min);
      position           = int(it-Distances.begin()); 
      segmento           = segment[position];
      }
      
     rCond = neighb_cond_master( segmento).lock(); 
     return true;
  }
  
  return false;
  KRATOS_CATCH("")
}

/// para cercanos
/// caso en que las aristas estan fuera
bool Test_Two( const NodePointerType& SlaveNode,
	       const PointerType& MasterObject,
	       Condition::Pointer& rCond) 
{
  
  KRATOS_TRY   
  
  WeakPointerVector<Condition>& neighb_cond     = MasterObject->GetValue(NEIGHBOUR_CONDITIONS);
 
  if(neighb_cond.size()!=0 )
  {
  std::vector<unsigned int> segment;
  unsigned int I         = 0;
  unsigned int segmento = 0;

  std::vector<array_1d<double, 2> > Points;   // punto de interseccion del segmento 
  vector<array_1d<double, 2> >      Points0;
  vector<array_1d<double, 2> >      Points1;
  array_1d<double, 2>               Point;
     
  array_1d<double,3>& old_pos     = SlaveNode->FastGetSolutionStepValue(DISPLACEMENT,3);   


  Points0.resize(2, false); 
  Points1.resize(2, false);

  Points0(0)[0] = SlaveNode->X0() + old_pos[0];
  Points0(0)[1] = SlaveNode->Y0() + old_pos[1];
  Points0(1)[0] = SlaveNode->X(); 
  Points0(1)[1] = SlaveNode->Y(); 
  
  
  unsigned int JJ = 1; 
  for(WeakPointerVector< Condition >::iterator cond  = neighb_cond.begin(); cond!= neighb_cond.end(); cond++){
  Condition::GeometryType& geom_2 = cond->GetGeometry();
	    
    Points1(0)[0] = geom_2[0].X(); 
    Points1(0)[1] = geom_2[0].Y();
    Points1(1)[0] = geom_2[1].X();
    Points1(1)[1] = geom_2[1].Y();

    if(IntersectionSegments::IntersectSegment(Point, Points0, Points1)==IT_POINT){
        Points.push_back(Point); 
        segment.push_back(I);
    }
    I++;
    JJ++;
    if(JJ>neighb_cond.size())
       break;
  }

  if (Points.size()!=0)
  {
      if (Points.size()==1){
      segmento = segment[0];
      }

      else if (Points.size()>1)
      {
      double dist0 = 0.00;
      array_1d<double, 2> rect;
      std::vector<double> Distance;
      std::vector<double>::iterator it;
      int position  = 0; 
      
      for(unsigned int i = 0; i<Points.size(); i++){
         rect  = Points0[1] - Points[i];
         dist0 = std::sqrt(inner_prod(rect, rect ));
         Distance.push_back(dist0);
      }
     
      const double max   = (*std::max_element(Distance.begin(), Distance.end() ) );
      it                 = std::find(Distance.begin(), Distance.end(), max);
      position           = int(it-Distance.begin()); 
      segmento           = segment[position];
      }
     
     rCond = neighb_cond(segmento).lock();
     return true;
  }
  }
  return false;
  KRATOS_CATCH("")
}



/// Busca la interseccion de la trayectoria de nodo esclavo
/// con los egbes de los elemntos. A Diferencia del test_One_B este lo hace con los elementos; no con su nodo master.
bool Test_Three(const NodePointerType& SlaveNode,
	        const PointerType& MasterObject,
	        Condition::Pointer& rCond) 
{
  
  KRATOS_TRY   
  
  WeakPointerVector<Condition>& neighb_cond     = MasterObject->GetValue(NEIGHBOUR_CONDITIONS);
  if(neighb_cond.size()!=0 )
  {
  std::vector<unsigned int> segment;
  unsigned int I         = 0;
  unsigned int segmento  = 0;

  std::vector<array_1d<double, 2> > Points;   // punto de interseccion del segmento 
  vector<array_1d<double, 2> >      Points0;
  vector<array_1d<double, 2> >      Points1;
  array_1d<double, 2>               Point;
     
  array_1d<double,3>& old_pos     = SlaveNode->FastGetSolutionStepValue(DISPLACEMENT,3);   


  Points0.resize(2, false); 
  Points1.resize(2, false);

  Points0(0)[0] = SlaveNode->X0() + old_pos[0];
  Points0(0)[1] = SlaveNode->Y0() + old_pos[1];
  Points0(1)[0] = SlaveNode->X(); 
  Points0(1)[1] = SlaveNode->Y(); 

  unsigned int JJ = 1; 
  for(WeakPointerVector< Condition >::iterator cond  = neighb_cond.begin(); cond!= neighb_cond.end(); cond++){
  Condition::GeometryType& geom_2 = cond->GetGeometry();
	    
    Points1(0)[0] = geom_2[0].X(); 
    Points1(0)[1] = geom_2[0].Y();
    Points1(1)[0] = geom_2[1].X();
    Points1(1)[1] = geom_2[1].Y();

    if(IntersectionSegments::IntersectSegment(Point, Points0, Points1)!=IT_EMPTY){
        Points.push_back(Point); 
        segment.push_back(I);
    }
    I++;
    JJ++;
    if(JJ>neighb_cond.size())
       break;
  }

  if (Points.size()!=0)
  {
      if (Points.size()==1){
      segmento = segment[0];
      }

      else if (Points.size()>1)
      {
      double dist0 = 0.00;
      array_1d<double, 2> rect;
      std::vector<double> Distance;
      std::vector<double>::iterator it;
      int position  = 0; 
      
      for(unsigned int i = 0; i<Points.size(); i++){
         rect  = Points0[1] - Points[i];
         dist0 = std::sqrt(inner_prod(rect, rect ));
         Distance.push_back(dist0);
      }
     
      const double min   = (*std::min_element(Distance.begin(), Distance.end() ) );
      it                 = std::find(Distance.begin(), Distance.end(), min);
      position           = int(it-Distance.begin()); 
      segmento           = segment[position];
      }
      
     rCond = neighb_cond(segmento).lock(); 
     return true;
  }
  }
  return false;
  KRATOS_CATCH("")
}


/// para cercanos
/// caso en que las aristas esten fuera de un elemento
bool Test_Four( const NodePointerType& SlaveNode,
	        const PointerType& MasterObject,
		Condition::Pointer& rCond)
{
 
  KRATOS_TRY   
  
  WeakPointerVector<Condition>& neighb_cond            = MasterObject->GetValue(NEIGHBOUR_CONDITIONS);
  WeakPointerVector<Condition>& neighb_cond_slave      = SlaveNode->GetValue(NEIGHBOUR_CONDITIONS);
    
  if(neighb_cond.size()!=0 && neighb_cond_slave.size()!=0)
  {
  array_1d<double,3>& old_pos                          = SlaveNode->FastGetSolutionStepValue(DISPLACEMENT,3);   
  std::vector<unsigned int>         segment;
  std::vector<array_1d<double, 2> > Points; // punto de interseccion del segmento
  vector<array_1d<double, 2> >      Points0;
  vector<array_1d<double, 2> >      Points1;
  array_1d<double, 2>               Point;

  // test with edges
  unsigned int segmento = 0;
  unsigned int I        = 0;
  unsigned int II       = 1;
  unsigned int III      = 1;

   
  Points0.resize(2, false); 
  Points1.resize(2, false);
  
  for(WeakPointerVector<Condition>::iterator cond_slave  = neighb_cond_slave.begin(); cond_slave!= neighb_cond.end(); ++cond_slave)
  {
      Condition::GeometryType& geom = cond_slave->GetGeometry();

      Point[0]      = 0.00;
      Point[1]      = 0.00;
      Points0(0)[0] = geom[0].X(); 
      Points0(0)[1] = geom[0].Y();
      Points0(1)[0] = geom[1].X();
      Points0(1)[1] = geom[1].Y(); 

      I   = 0; 
      III = 1; 
      for(WeakPointerVector< Condition >::iterator cond  = neighb_cond.begin(); cond!= neighb_cond.end(); ++cond){
        Condition::GeometryType& geom_3 = cond->GetGeometry();     
	Points1(0)[0] = geom_3[0].X(); 
	Points1(0)[1] = geom_3[0].Y();
	Points1(1)[0] = geom_3[1].X();
	Points1(1)[1] = geom_3[1].Y();

	if(IntersectionSegments::IntersectSegment(Point, Points0, Points1)==IT_POINT)
	{  
	  Points.push_back(Point); 
	  segment.push_back(I); 
	}
        I++;
        III++;
        if(III>neighb_cond.size())
            break;
     }     

     II++;
     if(II>neighb_cond_slave.size())
         break;
  }  

  if (Points.size()!=0)
  {
     if (Points.size()==1){
       segmento = segment[0];
     }
   // en caso de que el nodo quede fuera e intersecte con dos aristas 
  else if (Points.size()>1)
  {
      Points0(0)[0] = SlaveNode->X0() + old_pos[0];
      Points0(0)[1] = SlaveNode->Y0() + old_pos[1];
      Points0(1)[0] = SlaveNode->X();  
      Points0(1)[1] = SlaveNode->Y(); 

      double dist0 = 0.00;
      array_1d<double, 2> rect;
      std::vector<double> Distance;
      std::vector<double>::iterator it;
      int position  = 0; 
      
      for(unsigned int i = 0; i<Points.size(); i++){
         rect  = Points0[1] - Points[i];
         dist0 = std::sqrt(inner_prod(rect, rect ));
         Distance.push_back(dist0);
      }
     
      const double max   = (*std::max_element(Distance.begin(), Distance.end() ) );
      it                 = std::find(Distance.begin(), Distance.end(), max);
      position           = int(it-Distance.begin()); 
      segmento           = segment[position];     
  }
  rCond = neighb_cond(segmento).lock();
  return true;
  }
  }
  return false;
  
  KRATOS_CATCH("")
  
}




/// Buscla interseccion de las aristas del nodos slave con las aristas del nodo master
bool Test_Five( const NodePointerType& SlaveNode,
	        const PointerType& MasterObject,
	        Condition::Pointer& rCond)
{
 
  KRATOS_TRY   
  
  WeakPointerVector<Condition>& neighb_cond            = MasterObject->GetValue(NEIGHBOUR_CONDITIONS);
  WeakPointerVector<Condition>& neighb_cond_slave      = SlaveNode->GetValue(NEIGHBOUR_CONDITIONS);
  if(neighb_cond.size()!=0 && neighb_cond_slave.size()!=0)
  {
  array_1d<double,3>& old_pos                          = SlaveNode->FastGetSolutionStepValue(DISPLACEMENT,3);   
  std::vector<unsigned int>         segment;
  std::vector<array_1d<double, 2> > Points; // punto de interseccion del segmento
  vector<array_1d<double, 2> >      Points0;
  vector<array_1d<double, 2> >      Points1;
  array_1d<double, 2>               Point;

  // test with edges
  unsigned int segmento  = 0;
  unsigned int I         = 0;
  unsigned int II        = 1;
  unsigned int III       = 1;


  Points0.resize(2, false); 
  Points1.resize(2, false);
  
  for(WeakPointerVector<Condition>::iterator cond_slave  = neighb_cond_slave.begin(); cond_slave!= neighb_cond.end(); ++cond_slave)
  {
      Condition::GeometryType& geom = cond_slave->GetGeometry();
      Point[0] = 0.00; 
      Point[1] = 0.00;
      
      Points0(0)[0] = geom[0].X(); 
      Points0(0)[1] = geom[0].Y();
      Points0(1)[0] = geom[1].X();
      Points0(1)[1] = geom[1].Y(); 

      I   = 0; 
      III = 1;
      for(WeakPointerVector< Condition >::iterator cond  = neighb_cond.begin(); cond!= neighb_cond.end(); ++cond){
        Condition::GeometryType& geom_3 = cond->GetGeometry();     
	Points1(0)[0] = geom_3[0].X(); 
	Points1(0)[1] = geom_3[0].Y();
	Points1(1)[0] = geom_3[1].X();
	Points1(1)[1] = geom_3[1].Y();

	if(IntersectionSegments::IntersectSegment(Point, Points0, Points1)!=IT_EMPTY)
	{  
	  Points.push_back(Point); 
	  segment.push_back(I); 
	}
        I++;
        III++;
        if(III>neighb_cond.size())
            break;
     }     

     II++;
     if(II>neighb_cond_slave.size())
         break;
  }  

  if (Points.size()!=0)
  {
     if (Points.size()==1){
       segmento = segment[0];
     }
   // en caso de que el nodo quede fuera e intersecte con dos aristas 
  else if (Points.size()>1)
  {
      Points0(0)[0] = SlaveNode->X0() + old_pos[0];
      Points0(0)[1] = SlaveNode->Y0() + old_pos[1];
      Points0(1)[0] = SlaveNode->X();  
      Points0(1)[1] = SlaveNode->Y(); 

      double dist0 = 0.00;
      array_1d<double, 2> rect;
      std::vector<double> Distance;
      std::vector<double>::iterator it;
      int position  = 0; 
      
      for(unsigned int i = 0; i<Points.size(); i++){
         rect  = Points0[1] - Points[i];
         dist0 = std::sqrt(inner_prod(rect, rect ));
         Distance.push_back(dist0);
      }
     
      const double min   = (*std::min_element(Distance.begin(), Distance.end() ) );
      it                 = std::find(Distance.begin(), Distance.end(), min);
      position           = int(it-Distance.begin()); 
      segmento           = segment[position];      
  }
  
  rCond = neighb_cond(segmento).lock();
  return true;
  }
  }
  return false;
  
  KRATOS_CATCH("")
  
}




	
//************************************************************************************
//************************************************************************************ 
void CalculateBoundaryContour2D(ConditionsArrayType& MasterConditions)
	 {   
	      KRATOS_TRY
	      //std::cout<< std::endl; 
	      std::cout<<"  CALCULATING CONTOURS 2D" <<  std::endl; 
	     
	      typedef WeakPointerVector< Element >::iterator  ElementIteratorType; 
	      ContainerType& rElements           =  mr_model_part.ElementsArray();
	      ConditionsArrayType& rConditions   =  mr_model_part.Conditions();
	      IteratorType it_begin              =  rElements.begin();
	      IteratorType it_end                =  rElements.end(); 
	      
	      array_1d<NodePointerType,2>  Pair; 
	      
	      unsigned int face  = 0; 
	      unsigned int Id    = rConditions.size() + 1 ;
	      bool is_repited    =  false;  
	      for(IteratorType elem = it_begin; elem!=it_end; elem++)
	      {
		   Element::GeometryType& geom_1                 = (*elem)->GetGeometry();
		   WeakPointerVector< Element >& neighb_elems    = (*elem)->GetValue(NEIGHBOUR_ELEMENTS); 
		   //WeakPointerVector< Condition >& neighb_cond   = (*elem)->GetValue(NEIGHBOUR_CONDITIONS);
		   //neighb_cond.clear(); ///WARNING  
		   
		   //node_boundary.resize(neighb_elems.size(), false);
		   // Puede incluir como vecnino el mismo en caso de que hayan menos de 3 elemtos veninos.  
		   // ckeck si se repited elmento
		   // ElementIteratorType no necesita especificarse el * 
		   for( ElementIteratorType neighb_elem  = neighb_elems.begin(); neighb_elem!= neighb_elems.end(); neighb_elem++)
 		            { 
			       
			       if (neighb_elem->Id() ==  (*elem)->Id() )
			           {     
				     if(face == 0) // edge 1-2
				      {
					Pair[0]   =  geom_1(1);
					Pair[1]   =  geom_1(2);
					CreateMasterConditions2D(Pair, elem, Id, MasterConditions);
		                        geom_1[1].GetValue(IS_BOUNDARY) = 1; //FastGetSolutionStepValue(IS_BOUNDARY) = 1.00;
					geom_1[2].GetValue(IS_BOUNDARY) = 1; //FastGetSolutionStepValue(IS_BOUNDARY) = 1.00;
					Id++; 
			              }
			              
				      if (face==1)  // edge 2-0
				      {
					 Pair[0] =   geom_1(2);
					 Pair[1] =   geom_1(0);
					 CreateMasterConditions2D(Pair, elem, Id, MasterConditions);
					 geom_1[2].GetValue(IS_BOUNDARY) = 1; //FastGetSolutionStepValue(IS_BOUNDARY) = 1.00;
					 geom_1[0].GetValue(IS_BOUNDARY) = 1; //FastGetSolutionStepValue(IS_BOUNDARY) = 1.00;
					 Id++; 
					 
				      }
				      if (face==2) // edge 0-1
				      {
					 Pair[0] =   geom_1(0);
					 Pair[1] =   geom_1(1);
					 CreateMasterConditions2D(Pair, elem, Id, MasterConditions);
					 geom_1[0].GetValue(IS_BOUNDARY) = 1; //FastGetSolutionStepValue(IS_BOUNDARY) = 1.00;
					 geom_1[1].GetValue(IS_BOUNDARY) = 1; // FastGetSolutionStepValue(IS_BOUNDARY) = 1.00;
					 Id++; 
				      }
				      
				       if(is_repited==false)
				         {
					   (*elem)->GetValue(IS_BOUNDARY) = 1;
				           mBoundaryElements.push_back(*elem);
					   is_repited = true;
				         }
				   } 
			        face++;    
			      }    
		        face       = 0;  
			is_repited = false;
	           }
	           
	           unsigned int I = 0;
	           NodesContainerType& rNodes  =  mr_model_part.NodesArray();
	           for(NodesIteratorType inode =  rNodes.begin(); inode!=rNodes.end(); ++inode){
		       if((*inode)->GetValue(IS_BOUNDARY) == 1){
			   mBoundaryNodes.push_back(*inode);
		           WeakPointerVector<Element>& neighb_elems = (*inode)->GetValue(NEIGHBOUR_ELEMENTS);
			   I = 0;
		           for( ElementIteratorType neighb_elem  = neighb_elems.begin(); neighb_elem!= neighb_elems.end(); ++neighb_elem){
			       if(neighb_elem->GetValue(IS_BOUNDARY)!=1){
				     neighb_elem->GetValue(IS_BOUNDARY)=1;
                                     mBoundaryElements.push_back(neighb_elems(I).lock());
			       }
			      I++; 
			   }
				    
		       }
		   }
	                          
		   KRATOS_CATCH("")
	 }


//*****************************************************************************************************
//*****************************************************************************************************

void CalculateBoundaryContour3D(ConditionsArrayType& MasterConditions)
{
      KRATOS_TRY

      //std::cout<< std::endl; 
      std::cout<<"CALCULATING CONTOURS 3D"<<  std::endl; 

      typedef WeakPointerVector< Element >::iterator  ElementIteratorType; 
      ContainerType& rElements           =  mr_model_part.ElementsArray();
      ConditionsArrayType& rConditions   =  mr_model_part.Conditions();
      IteratorType it_begin              =  rElements.begin();
      IteratorType it_end                =  rElements.end(); 

      array_1d<NodePointerType,3>  Pair; 

      unsigned int face  = 0; 
      unsigned int Id    = rConditions.size() + 1 ;
      bool is_repited    =  false;  
      for(IteratorType elem = it_begin; elem!=it_end; elem++)
      {
        Element::GeometryType& geom_1                 = (*elem)->GetGeometry();
        WeakPointerVector< Element >& neighb_elems    = (*elem)->GetValue(NEIGHBOUR_ELEMENTS); 
        //WeakPointerVector< Condition >& neighb_cond   = (*elem)->GetValue(NEIGHBOUR_CONDITIONS);
        for( ElementIteratorType neighb_elem  = neighb_elems.begin(); neighb_elem!= neighb_elems.end(); neighb_elem++){
           if(neighb_elem->Id() ==  (*elem)->Id() )
	     {
		if(face == 0) 
		{
		  Pair[0]   =  geom_1(1);
		  Pair[1]   =  geom_1(2);
		  Pair[2]   =  geom_1(3);
		  CreateMasterConditions3D(Pair, elem, Id, MasterConditions);
		  geom_1[1].GetValue(IS_BOUNDARY) = 1; 
		  geom_1[2].GetValue(IS_BOUNDARY) = 1;
		  geom_1[3].GetValue(IS_BOUNDARY) = 1;
		  Id++; 
		}
	       if(face ==1) 
		{
		  Pair[0]   =  geom_1(0);
		  Pair[1]   =  geom_1(3);
		  Pair[2]   =  geom_1(2);
		  CreateMasterConditions3D(Pair, elem, Id, MasterConditions);
		  geom_1[0].GetValue(IS_BOUNDARY) = 1; 
		  geom_1[3].GetValue(IS_BOUNDARY) = 1;
		  geom_1[2].GetValue(IS_BOUNDARY) = 1;
		  Id++; 
		}
		if(face == 2) 
		{
		  Pair[0]   =  geom_1(0);
		  Pair[1]   =  geom_1(1);
		  Pair[2]   =  geom_1(3);
		  CreateMasterConditions3D(Pair, elem, Id, MasterConditions);
		  geom_1[0].GetValue(IS_BOUNDARY) = 1; 
		  geom_1[1].GetValue(IS_BOUNDARY) = 1;
		  geom_1[3].GetValue(IS_BOUNDARY) = 1;
		  Id++; 
		}
		if(face == 3) 
		{
		  Pair[0]   =  geom_1(0);
		  Pair[1]   =  geom_1(2);
		  Pair[2]   =  geom_1(1);
		  CreateMasterConditions3D(Pair, elem, Id, MasterConditions);
		  geom_1[0].GetValue(IS_BOUNDARY) = 1; 
		  geom_1[2].GetValue(IS_BOUNDARY) = 1;
		  geom_1[1].GetValue(IS_BOUNDARY) = 1;
		  Id++; 
		}
		  if(is_repited==false)
		  {
		  (*elem)->GetValue(IS_BOUNDARY) = 1;
		  mBoundaryElements.push_back(*elem);
		  is_repited = true;
		  }
		} 
		face++;    
		}    
		face       = 0;  
		is_repited = false;
             }
             
	      unsigned int I = 0;
	      NodesContainerType& rNodes  =  mr_model_part.NodesArray();
	      for(NodesIteratorType inode =  rNodes.begin(); inode!=rNodes.end(); ++inode){
	           if((*inode)->GetValue(IS_BOUNDARY) == 1){
	           WeakPointerVector<Element>& neighb_elems = (*inode)->GetValue(NEIGHBOUR_ELEMENTS);
	           I = 0;
	           for( ElementIteratorType neighb_elem  = neighb_elems.begin(); neighb_elem!= neighb_elems.end(); ++neighb_elem){
	              if(neighb_elem->GetValue(IS_BOUNDARY)!=1){
	                 neighb_elem->GetValue(IS_BOUNDARY)=1;
	                 mBoundaryElements.push_back(neighb_elems(I).lock());
	               }
	             I++; 
	           }
	         }
	      }
	      
             
  KRATOS_CATCH("")
}

//*****************************************************************************************************
//*****************************************************************************************************
void CreateMasterConditions2D(const array_1d<NodePointerType,2>&  rPair, const IteratorType& elem, const unsigned int& Id, ConditionsArrayType& MasterConditions)
	 {
	    KRATOS_TRY
	   
	    Line2D2<Node<3> >::Pointer pgeom =  Line2D2<Node<3> >::Pointer (new Line2D2<Node<3> >( rPair[0], rPair[1] ) ) ;  
	    Condition::Pointer MasterSegment = Condition::Pointer(new MasterContactFaceType(Id, pgeom ) ) ;
	    MasterSegment->GetValue(NEIGHBOUR_ELEMENTS).push_back(*(elem));
	    
	    ((rPair)[0])->GetValue(NEIGHBOUR_CONDITIONS).push_back(MasterSegment); 
	    ((rPair)[1])->GetValue(NEIGHBOUR_CONDITIONS).push_back(MasterSegment); 
	    
	    (*elem)->GetValue(NEIGHBOUR_CONDITIONS).push_back(MasterSegment); 
	    MasterConditions.push_back(MasterSegment);
    
	    ConditionsArrayType& pConditions = mr_model_part.Conditions();
	    pConditions.push_back(MasterSegment);
	    
	    KRATOS_CATCH("")
	 }


void CreateMasterConditions3D(const array_1d<NodePointerType,3>&  rPair, 
			      const IteratorType& elem, 
			      const unsigned int& Id, 
			      ConditionsArrayType& MasterConditions)
	 {
	    KRATOS_TRY
	   
	    Triangle2D3<Node<3> >::Pointer pgeom  =  Triangle2D3<Node<3> >::Pointer( new Triangle2D3<Node<3> > ( rPair[0], rPair[1], rPair[2]) );  
							
 	    Condition::Pointer MasterSurface      =  Condition::Pointer(new MasterContactFace3D(Id, pgeom) ); 
 	    MasterSurface->GetValue(NEIGHBOUR_ELEMENTS).push_back(*(elem));
 	    MasterSurface->GetValue(IS_BOUNDARY) = 1;
 	    rPair[0]->GetValue(NEIGHBOUR_CONDITIONS).push_back(MasterSurface); 
 	    rPair[1]->GetValue(NEIGHBOUR_CONDITIONS).push_back(MasterSurface); 
 	    rPair[2]->GetValue(NEIGHBOUR_CONDITIONS).push_back(MasterSurface); 
 	    
 	    (*elem)->GetValue(NEIGHBOUR_CONDITIONS).push_back(MasterSurface); 
 	    MasterConditions.push_back(MasterSurface);
     
 	    ConditionsArrayType& pConditions = mr_model_part.Conditions();
 	    pConditions.push_back(MasterSurface);
	    
	    KRATOS_CATCH("")
	 }



//*****************************************************************************************************
//*****************************************************************************************************
void FiltratePairContacts2D(ContainerContactPair& PairContacts)
	{
	  KRATOS_TRY
	 
	  ContainerContactPair temp;
	  
	  std::vector<unsigned int> id;  
	  for(IteratorContainerContactPair ipair = PairContacts.begin(); ipair!=PairContacts.end(); ipair++){
	     if(SearchCommonNode( (*ipair)[0], (*ipair)[1], id )) { 
	        /// Se localiza que comparte dos nodos en comun
                if( id.size()!=2 && (SearchInsideNode((*ipair)[0], (*ipair)[1], id[0]))==true) 
		 {
		   temp.push_back(*(ipair));  
		 }
	      }
	      else 
	      {
		temp.push_back(*(ipair)); 
	      }
	      
	     id.clear();
	  }
	  PairContacts.swap(temp);
	  
	  KRATOS_CATCH("")
	}
	
//*****************************************************************************************************
//*****************************************************************************************************
	void FiltratePairContacts3D(ContainerContactPair& PairContacts)
	{
	  KRATOS_TRY
	  
	  ContainerContactPair temp;
	  std::vector<unsigned int> id;  
	  for(IteratorContainerContactPair ipair = PairContacts.begin(); ipair!=PairContacts.end(); ipair++){
	     if(SearchCommonNode( (*ipair)[0], (*ipair)[1], id )==false) { 
	        /*
	       if( id.size()!=3){
		  
	          for(unsigned int i = 0; i<id.size(); i++)
                     if(SearchInsideNode((*ipair)[0], (*ipair)[1], id[i] )==true) 
		      {
			KRATOS_WATCH(id[i]) 
		        temp.push_back(*(ipair)); 
			break;
		     }
		 }
	      }
	      
	      else 
	      */	
	      //{
		temp.push_back(*(ipair)); 
	      //}
	      
	     id.clear();
	    }
	  }
	  
	  PairContacts.swap(temp);
	  
	  KRATOS_CATCH("")
	}
	
	
bool FiltratePairContacts(const PointerType& elem1, const PointerType& elem2)
{
   KRATOS_TRY
   std::vector<unsigned int> id;  
   if(mrdimension==2)
    { 
      const bool test_1 = SearchCommonNode(elem1, elem2, id);   
      if(test_1){ 
	/// Se localiza que comparte dos nodos en comun
	if(id.size()==2)
	  return false;
	else if(id.size()!=2 && (SearchInsideNode(elem1, elem2, id[0])==true)) 
	    return true; 
      }
      else 
	return true; 
      
    }
    else
    {
       if(SearchCommonNode(elem1, elem2, id)){
          return false;
       }
       else{
	 return true;
         }
    }
  return false;  
  KRATOS_CATCH("")
}


//*****************************************************************************************************
//*****************************************************************************************************	
/// Se busca el nodo comun entre los contactos
bool SearchCommonNode(const PointerType& elem1, const PointerType& elem2, std::vector<unsigned int>& id)
      {
	  KRATOS_TRY
	  
          Element::GeometryType& geom1  = (elem1)->GetGeometry();
	  Element::GeometryType& geom2  = (elem2)->GetGeometry();
	  /// buscando el nodo comun
	  for(unsigned int i = 0; i<geom1.size(); i++){
	     for(unsigned int j = 0; j<geom1.size(); j++) {
	         if(geom1[i].Id()==geom2[j].Id()) {id.push_back(geom1[i].Id());}
	     }
	  }
	  
	  return  id.size()!=0;
	  
	  //if( id.size()!=0) return true;
	  //return false;    
	  
	  KRATOS_CATCH("")
      }
	
/// Verifica si los nodos que no es el comun cae dentro del elemento
//*****************************************************************************************************
//*****************************************************************************************************
bool SearchInsideNode(const PointerType& elem1, const PointerType& elem2, const unsigned int& ide) 
      {  
	
	  KRATOS_TRY
	 
          Element::GeometryType& geom1  = (elem1)->GetGeometry();
	  Element::GeometryType& geom2  = (elem2)->GetGeometry();
	  array_1d<double, 3> result;
 
	  ///CoordinatesArrayType result;
	  /// buscando si uno de los nodos entra dentro del elemento
	  for(unsigned int i = 0; i<geom1.size(); i++){
	      if(geom2[i].Id()!=ide) {
	           if(geom1.IsInside(geom2[i], result)) 
		   {
		     return true;
		   }
	      }
	  }  
	  for(unsigned int i = 0; i<geom2.size(); i++){
	      if(geom1[i].Id()!=ide)
	      {
	         if(geom2.IsInside(geom1[i], result)) 
		 {
		   return true; 
		 }
	      }
	  }

	  return false;
	  
	  KRATOS_CATCH("")
      }


//*****************************************************************************************************
//*****************************************************************************************************

void NodeInside(const PointerType& MasterObject, const PointerType& SlaveObject,  std::vector<NodePointerType>& InsideNodes)
         {
	   
	   KRATOS_TRY
	   
	   Element::GeometryType& geom_master = MasterObject->GetGeometry();
	   Element::GeometryType& geom_slave  = SlaveObject->GetGeometry();
	   std::vector<unsigned> Nodes;
	   
	   /// buscando el nodo comun
	   bool commun = false;
	   for(unsigned int i = 0; i<geom_slave.size(); i++){
	     commun = false;
	     for(unsigned int j = 0; j<geom_master.size(); j++) {
	         if(geom_slave[i].Id()==geom_master[j].Id()) 
		 {
		   commun = true;
		 }
	      }
	          if(commun==false)
		     Nodes.push_back(i);
	     }
	  
	   array_1d<double, 3> result;
	   for (unsigned int i = 0; i<Nodes.size(); i++ ){
	       if(geom_master.IsInside(geom_slave[Nodes[i]], result)){
		     InsideNodes.push_back(geom_slave(Nodes[i])); }   
	        }
	   
	   
	   KRATOS_CATCH("")
	      
	 }

       void ResetFlagComputeBoundaryContour(const bool& rflag)
       {
	  mcompute_boundary_contour = rflag; 
       }
       
       
       private:
       ModelPart mr_model_part;    
       unsigned int mrdimension;
       double       mpenalty_factor;  
       bool mcompute_boundary_contour;
       
       NodesContainerType            mBoundaryNodes;
       ContainerType                 mBoundaryElements;
       ContainerContactPair          mPairContacts;
       ConditionsArrayType           mMasterConditionsArray;
       //WeakPointerVector<NodeType>   mBoundaryNodes;       

       

//*****************************************************************************************************
//*****************************************************************************************************

/// WARNING = To be parallel
void IdentifyMasterSegment2D()
{
  KRATOS_TRY
  
  std::cout<< "     IDENTIFYING THE MASTER 2D SEGMENT         " <<std::endl; 
  NodesArrayType& pNodes           = mr_model_part.Nodes();
  ProcessInfo& CurrentProcessInfo  = mr_model_part.GetProcessInfo(); 

  #ifdef _OPENMP
  int number_of_threads = 1; //omp_get_max_threads();
  #else
  int number_of_threads = 1;
  #endif

  vector<unsigned int> node_partition;
  CreatePartition(number_of_threads, pNodes.size(), node_partition);
  
  //int    I                     = 0;
  //double g                     = 0.00;
  //double g_old                 = 0.00;
  double gl                    = 0.00;
  double gr                    = 0.00;
  //double compare_distance      = 0.00;
  double pr                    = 0.00;
  double pl                    = 0.00;
  double wr                    = 0.00;
  double wl                    = 0.00;
  
  array_1d<double, 3> CS         = ZeroVector(3);
  array_1d<double, 3> CL         = ZeroVector(3);
  array_1d<double, 3> CR         = ZeroVector(3);
  array_1d<double, 3> Normal     = ZeroVector(3);
  array_1d<double, 3> Normal_s   = ZeroVector(3);
  array_1d<double, 3> Normal_r   = ZeroVector(3);
  array_1d<double, 3> Normal_l   = ZeroVector(3);
  array_1d<double, 3> GR         = ZeroVector(3);
  array_1d<double, 3> GL         = ZeroVector(3);
  double distance_r              = 0.00;
  double distance_l              = 0.00;
  double N                       = 0.00;
  array_1d<double, 3> e3;
  
  e3[0] = 0;
  e3[1] = 0;
  e3[2] = 1.00;
  Segment2D rSegment_r;
  Segment2D rSegment_l;
  
  //#pragma omp parallel for private(I, g, g_old, gl, gr, compare_distance, CS, CL, CR, Normal, rSegment)
  for(int k=0; k<number_of_threads; k++)
  {
    NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
    NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];
    for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i){
       if(i->GetValue(IS_BOUNDARY)==1){ 
	 WeakPointerVector<Condition>& neighb_cond       = ((i)->GetValue(NEAR_NODE))->GetValue(NEIGHBOUR_CONDITIONS);
	 WeakPointerVector<Condition>& neighb_cond_slave = (i)->GetValue(NEIGHBOUR_CONDITIONS);
	 
	 if( neighb_cond.size()!=0 &&  neighb_cond_slave.size()!=0 )
	 {
	   array_1d<int,2> RL;
	   //double& distance                    =  i->GetValue(DISTANCE); 
	   const unsigned int& ID              =  i->GetValue(NEAR_NODE)->Id();
	   const unsigned int& ID_slave        =  i->Id();
	   
	   /// the right and left
           RL[0] = 0;
	   RL[1] = 0;
	   if(neighb_cond(0).lock()->GetGeometry()(0)->Id()==ID)
	     RL[0]  = 1;
	   else
	     RL[1]  = 1;
	   const Condition::Pointer& first     =  neighb_cond(RL[0]).lock();
	   const Condition::Pointer& last      =  neighb_cond(RL[1]).lock();
	   
	   /// the right and left
           RL[0] = 0;
	   RL[1] = 0;
	   if(neighb_cond_slave(0).lock()->GetGeometry()(0)->Id()==ID_slave)
	     RL[0]  = 1;
	   else
	     RL[1]  = 1;
	   
	   const Condition::Pointer& first_s   =  neighb_cond_slave(RL[0]).lock();
	   const Condition::Pointer& last_s    =  neighb_cond_slave(RL[1]).lock();
	   
	   Condition::GeometryType& geom_first =  first->GetGeometry();
	   Condition::GeometryType& geom_last  =  last->GetGeometry();
	   array_1d<double, 3>& point_closest  =  i->GetValue(NEAR_NODE)->Coordinates();                     
	   array_1d<double, 3>& point_slave    =  i->Coordinates();  
	   array_1d<double, 3>& point_left     =  geom_last.GetPoint(1);
	   array_1d<double, 3>& point_right    =  geom_first.GetPoint(0);
	   noalias(CS)                         =  point_closest - point_slave;   
	   noalias(CL)                         =  point_left    - point_closest;  
	   noalias(CR)                         =  point_closest - point_right;
	   noalias(GR)                         =  point_slave   - point_right;
	   noalias(GL)                         =  point_slave   - point_closest; 
	   
	   //KRATOS_WATCH("CCCCCCCCCCC")
	   //KRATOS_WATCH(i->Id())
	   //KRATOS_WATCH(neighb_cond_slave.size())
	   //KRATOS_WATCH(first_s->Id())
	   //KRATOS_WATCH(last_s->Id())
	   first_s->Calculate(NORMAL, Normal_r, CurrentProcessInfo);
	   last_s->Calculate(NORMAL,  Normal_l, CurrentProcessInfo);
	   //KRATOS_WATCH("AKIIIIIIIIIIIIIIIII")
	   
	   noalias(Normal_s) = Normal_r + Normal_l;  
	   Normal_s          = (1.00/norm_2(Normal_s)) * Normal_s;
	      
	   first->Calculate(NORMAL, Normal_r, CurrentProcessInfo);
	   last->Calculate(NORMAL,  Normal_l, CurrentProcessInfo);
	   
	  
	   //const double& cs  =  norm_2(CS);
	   const double& cl  =  norm_2(CL);
	   const double& cr  =  norm_2(CR);
	   noalias(CL)       =  CL * (1.00/cl);
	   noalias(CR)       =  CR * (1.00/cr);

	   
           gr               =  inner_prod(GR,Normal_r);
	   gl               =  inner_prod(GL,Normal_l);
	   pr               =  -inner_prod(GL,CR);  if(pr<=0.00)  pr =  0.00;
	   pl               =  inner_prod(GL,CL);   if(pl<=0.00)  pl =  0.00;
	   
	   if(std::fabs(pr-pl)<=1E-14){ pr = 1.00; pl =1.00;}
	   wr               =  pr/(pr + pl);   
	   wl               =  1.00 - wr;   
	   
	   N                = norm_2( pr *  Normal_r + pl *  Normal_l );
	   noalias(Normal)  = (1.00/N) * ( pr * Normal_r + pl *  Normal_l);
	       
	   rSegment_r.AssignPointsAndComputeParameters(point_right, point_closest);
	   rSegment_l.AssignPointsAndComputeParameters(point_closest, point_left);
           distance_r = rSegment_r.DistPoint2Segment2D(point_slave);
	   distance_l = rSegment_l.DistPoint2Segment2D(point_slave);
           
	   /// decidiendo el mejor segmento
	   
	   const double& dot_l = inner_prod(Normal_s,Normal_l);
	   const double& dot_r = inner_prod(Normal_s,Normal_r);
	   
	   if(distance_r <= distance_l) { 
	       if(wr==wl)                                              /// cualquier segmento es valido
	          i->GetValue(CONTACT_LINK_MASTER) = first;
	       else if(wr > wl){                                        /// segmento derecho
		    if(dot_l < dot_r) 
	               i->GetValue(CONTACT_LINK_MASTER) = last;
		    else
		       i->GetValue(CONTACT_LINK_MASTER) = first;
	          }
	       else if(wr < wl){                                        /// segmento derecho
		  if(dot_l < dot_r) 
	               i->GetValue(CONTACT_LINK_MASTER) = last;
		  else
		       i->GetValue(CONTACT_LINK_MASTER) = first;
	        }  
	    }

	   if(distance_r > distance_l){ 
	       if(wr==wl)                                               /// cualquier segmento es valido
	          i->GetValue(CONTACT_LINK_MASTER) = first;
	       else if(wr > wl){                                        /// segmento derecho
		    if(dot_l < dot_r) 
	               i->GetValue(CONTACT_LINK_MASTER) = last;
		    else
		       i->GetValue(CONTACT_LINK_MASTER) = first;
	            }
	       else if(wr < wl){                                        /// segmento derecho
		    if(dot_l < dot_r) 
	               i->GetValue(CONTACT_LINK_MASTER) = last;
		    else
		       i->GetValue(CONTACT_LINK_MASTER) = first;
	       }
	   }
	   


//  	  if( i->Id()==44) //|| i->Id()==68 || i->Id()==49) //|| i->Id()==189) //(i->Id()==926 || i->Id()==927 || i->Id()==910 || i->Id()==919 || i->Id()==905 || i->Id()==904 )
//  	    {
//  	      KRATOS_WATCH(i->Id())
// // 	      KRATOS_WATCH(ID)
// // 	      KRATOS_WATCH(first->Id())
// // 	      KRATOS_WATCH(last->Id())
// // 	      KRATOS_WATCH(pr)
// //               KRATOS_WATCH(pl)
// //               KRATOS_WATCH(wr)
// //               KRATOS_WATCH(wl)
// //               KRATOS_WATCH(gr)
// //               KRATOS_WATCH(gl)
// //               KRATOS_WATCH(distance_r)
// //               KRATOS_WATCH(distance_l)
// //               KRATOS_WATCH(Normal_r)
// //               KRATOS_WATCH(Normal_l)
// //               KRATOS_WATCH(Normal_s)
// //               KRATOS_WATCH(dot_l)
// //               KRATOS_WATCH(dot_r)
//               KRATOS_WATCH(i->GetValue(NEAR_NODE)->Id());    
//  	      KRATOS_WATCH(i->GetValue(CONTACT_LINK_MASTER)->Id());
// // 	      const int& id = i->GetValue(CONTACT_LINK_MASTER)->Id();     
// //  	      if(id==58) { KRATOS_WATCH(id); KRATOS_ERROR(std::logic_error,  "" , "");}
// //  	      if(id==46) { KRATOS_WATCH(id); KRATOS_ERROR(std::logic_error,  "" , "");}
//  	      KRATOS_WATCH("---------------------")
//  	    }
	    
	  }
        }
      }
   }
	    
  KRATOS_CATCH("")
}   


//*****************************************************************************************************
//*****************************************************************************************************

       template<class TConfigure>
       void LocalSearch2D( BinsObjectDynamic<TConfigure>& rBins,
			   const IteratorType& it_begin,
			   const IteratorType& it_end)
       {
	  std::cout<< "     LOCAL SEARCH ALGORITHM         " <<std::endl; 
	  unsigned int I                    = 0; 
	  double compare_distance           = 0.00;
	  ResultContainerType               Result;
	  	  
	  array_1d<double, 3>               Normal    = ZeroVector(3);
	  array_1d<double, 3>               Mid_Point = ZeroVector(3);
	  array_1d<double, 3>               Vect      = ZeroVector(3);
	   
	  Segment2D rSegment;
	  for(IteratorType it = it_begin; it!=it_end; it++)
	  { 
	    std::size_t size = rBins.SearchObjectsInner(*it, Result);
	    
	    if(Result.size()!=0){
	    Element::GeometryType& geom = (*it)->GetGeometry();
	    for(unsigned int i = 0; i<geom.size(); i++){       
	      if(geom(i)->GetValue(IS_BOUNDARY) == 1){
		  array_1d<double, 3>& Points0 = geom.GetPoint(i);   
		  double& distance   = geom(i)->GetValue(DISTANCE); 
	          for(ResultIteratorType rthis = Result.begin(); rthis!=Result.end(); rthis++){            
		    I = 0;
		    WeakPointerVector<Condition>& neighb_cond = (*rthis)->GetValue(NEIGHBOUR_CONDITIONS);  
		    if(neighb_cond.size()!=0){ 
		    for(WeakPointerVector<Condition>::iterator neighb = neighb_cond.begin();  neighb!= neighb_cond.end(); neighb++){
		      //if(geom(i)->Id()==3)
		      //   KRATOS_WATCH(neighb->Id()) 
		      if(neighb->GetValue(NODAL_VALUES) == 0){
			neighb->GetValue(NODAL_VALUES)  = 1;
		        Condition::GeometryType& geom_2 = (neighb)->GetGeometry();
			//(neighb)->Calculate(NORMAL, Normal, CurrentProcessInfo);
			if( (geom_2(0)->Id() != geom(i)->Id()) && (geom_2(1)->Id() != geom(i)->Id()) ){
			  array_1d<double, 3>&  Points1  = geom_2.GetPoint(0);
			  array_1d<double, 3>&  Points2  = geom_2.GetPoint(1);
			  rSegment.AssignPointsAndComputeParameters(Points1, Points2);
			  compare_distance = rSegment.DistPoint2Segment2D(Points0);
			  if(compare_distance<distance){
				   distance = compare_distance;
				   geom(i)->GetValue(CONTACT_LINK_MASTER) = neighb_cond(I).lock();
			     } 
			   }
			 I++;
	                }
		       }
		      }
		     }
		    
		    //if(geom(i)->Id()==3){
		    //KRATOS_WATCH(geom(i)->GetValue(CONTACT_LINK_MASTER)->Id())
		    //KRATOS_WATCH("--------------------------")
		    //}
		    /// Reseting the values 
		    for(ResultIteratorType rthis = Result.begin(); rthis!=Result.end(); rthis++){            
		       WeakPointerVector<Condition>& neighb_cond = (*rthis)->GetValue(NEIGHBOUR_CONDITIONS);
		          if(neighb_cond.size()!=0)
		              for(WeakPointerVector<Condition>::iterator neighb = neighb_cond.begin();  neighb!= neighb_cond.end(); neighb++)
				  neighb->GetValue(NODAL_VALUES) = 0;
		       }
	            }
	           }
	         }
	         Result.clear(); 
	      }
	   }
      
      

//*****************************************************************************************************
//*****************************************************************************************************

 void SearchNearNode2D()
 {
  
   KRATOS_TRY
   
   std::cout<< "     SEARCHING NEAR NODE " <<std::endl; 
   #ifdef _OPENMP
   int number_of_threads = omp_get_max_threads();
   #else
   int number_of_threads = 1;
   #endif
   vector<unsigned int> node_partition;
   int distance   = std::distance(mBoundaryNodes.begin(), mBoundaryNodes.end());
   CreatePartition(number_of_threads, distance, node_partition);
   BinsDynamic<2, NodeType, NodesContainerType> BinsPoint(mBoundaryNodes.begin(), mBoundaryNodes.end());
   #pragma omp parallel for 
   for(int k=0; k<number_of_threads; k++)
      {
       NodesIteratorType it_begin = mBoundaryNodes.begin() + node_partition[k];
       NodesIteratorType it_end   = mBoundaryNodes.begin() + node_partition[k+1];
       for(NodesIteratorType inode = it_begin; inode!=it_end; inode++){
	  
          (*inode)->GetValue(NEAR_NODE) = BinsPoint.SearchNearestPoint(**inode);
	  KRATOS_WATCH((*inode)->Id()) 
	  KRATOS_WATCH((*inode)->GetValue(NEAR_NODE)->Id()) 
	  KRATOS_WATCH("----------------------------") 
       }
      }
   KRATOS_CATCH("")
 }
   
//*****************************************************************************************************
//*****************************************************************************************************

       template<class TConfigure>
       void SearchNearNode2D( BinsObjectDynamic<TConfigure>& rBins,
			 const IteratorType& it_begin,
			 const IteratorType& it_end)
       {
	 
	  KRATOS_TRY
	  std::cout<< "     SEARCHING NEAR NODE 2D  " <<std::endl; 
	  //ProcessInfo& CurrentProcessInfo  = mr_model_part.GetProcessInfo(); 
	  
	  unsigned int I                    = 0; 
	  double compare_distance           = 0.00;
	  ResultContainerType               Result;
          
	  array_1d<double, 3> Vect = ZeroVector(3);  
	  for(IteratorType it = it_begin; it!=it_end; it++)
	  { 
	    Result.clear();
	    rBins.SearchObjectsInner(*it, Result); ///SearchAroundObjectsInner(*it, Result); poner el cmentario para 2D
	    if(Result.size()!=0){
	    Element::GeometryType& geom = (*it)->GetGeometry();
	    for(unsigned int i = 0; i<geom.size(); i++){
	      if(geom(i)->GetValue(IS_BOUNDARY) == 1){ 
		  array_1d<double, 3>& Points0 = geom.GetPoint(i);   
		  double& distance             = geom(i)->GetValue(DISTANCE); 
	          for(ResultIteratorType rthis = Result.begin(); rthis!=Result.end(); rthis++){            
		    I = 0;
		    WeakPointerVector<Condition>& neighb_cond = (*rthis)->GetValue(NEIGHBOUR_CONDITIONS);      
		    if(neighb_cond.size()!=0){ 
		    for(WeakPointerVector<Condition>::iterator neighb = neighb_cond.begin();  neighb!= neighb_cond.end(); neighb++){
		       if(neighb->GetValue(NODAL_VALUES) == 0){
			neighb->GetValue(NODAL_VALUES)  = 1;
		        Condition::GeometryType& geom_2 = (neighb)->GetGeometry();
			if( (geom_2(0)->Id() != geom(i)->Id()) && (geom_2(1)->Id() != geom(i)->Id()) ){
			/// buscando el nodo mas cercano
			for(unsigned int k = 0; k<geom_2.size(); k++){
			    array_1d<double, 3>&  Points1  = geom_2.GetPoint(k);
			    noalias(Vect)                  = Points1 - Points0;
			    compare_distance               = norm_2(Vect); 
	                    if(compare_distance<distance)
			     { 
			       distance = compare_distance;
			       geom(i)->GetValue(NEAR_NODE) =  geom_2(k);
			     }
			   } 
			  }
		        }
		      }
		    }
		  }
		     
		  /// Reseting the values 
		    for(ResultIteratorType rthis = Result.begin(); rthis!=Result.end(); rthis++){            
		       WeakPointerVector<Condition>& neighb_cond = (*rthis)->GetValue(NEIGHBOUR_CONDITIONS);
		          if(neighb_cond.size()!=0)
		              for(WeakPointerVector<Condition>::iterator neighb = neighb_cond.begin();  neighb!= neighb_cond.end(); neighb++)
				  neighb->GetValue(NODAL_VALUES) = 0;
		       }
	         }
	       }
	     }
	   }
	   /// reseting the values of distance
	   NodesArrayType& pNodes           = mr_model_part.Nodes();
	   
	   #ifdef _OPENMP
           int number_of_threads = omp_get_max_threads();
           #else
           int number_of_threads = 1;
           #endif

           vector<unsigned int> node_partition;
           CreatePartition(number_of_threads, pNodes.size(), node_partition);

           #pragma omp parallel for 
           for(int k=0; k<number_of_threads; k++)
             {
	        NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	        NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];
	        for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i){   
		      i->GetValue(DISTANCE)          = DBL_MAX;
	          }
	     }
	   
	   
	  KRATOS_CATCH("")
       }
       
			  
			  
//*****************************************************************************************************
//*****************************************************************************************************

       
       template<class TConfigure>
       void LocalSearch3D( BinsObjectDynamic<TConfigure>& rBins,
			 const IteratorType& it_begin,
			 const IteratorType& it_end)
       {
	 
	 KRATOS_TRY 
	 
	  #ifdef _OPENMP
	  double start_prod = omp_get_wtime();   
	  #endif
	 
	  unsigned int I                    = 0; 
	  double distance                   = 0.00;
	  double compare_distance           = 0.00;

	  ResultContainerType               Result;
	  array_1d<double, 3>               Points0;
          array_1d<double, 3>               Points1;
	  array_1d<double, 3>               Points2;
	  array_1d<double, 3>               Points3;
	  //array_1d<double, 3>               Normal;
	  
	  Plane rPlane; 
	  
	  #ifdef _OPENMP
	  int number_of_threads = omp_get_max_threads();
	  #else
	  int number_of_threads = 1;
	  #endif

	  vector<unsigned int> partition;
	  int distance_2 = int(it_end-it_begin);
	  CreatePartition(number_of_threads, distance_2, partition);

	 #pragma omp parallel for private(I, distance, compare_distance, Result, Points0, Points1, Points2, Points3)
	 for(int k=0; k<number_of_threads; k++)          
	  {
	    IteratorType it_begin_1 = it_begin + partition[k];
	    IteratorType it_end_1   = it_begin + partition[k+1];
	    for(IteratorType it =it_begin_1; it!=it_end_1; it++) 
	    //for(IteratorType it = it_begin; it!=it_end; it++)
	   { 
	    rBins.SearchObjectsInner(*it, Result); 
	          
	    if(Result.size()!=0){
	    Element::GeometryType& geom = (*it)->GetGeometry();
	    for(unsigned int i = 0; i<geom.size(); i++)
	     {   
	       if(geom[i].GetValue(NODAL_VALUES)==0)
	       {
	        if(geom[i].GetValue(IS_BOUNDARY) == 1)
	        {
		  geom[i].SetLock();
		  geom[i].GetValue(NODAL_VALUES)  = 1;
		  geom[i].UnSetLock();
		  Points0[0]                      = geom[i].X(); 
		  Points0[1]                      = geom[i].Y();
		  Points0[2]                      = geom[i].Z();
		  distance                        = DBL_MAX;
		  
	          for(ResultIteratorType rthis = Result.begin(); rthis!=Result.end(); rthis++)
	          {
		    I = 0;
		    //KRATOS_WATCH((*rthis)->Id())
		    WeakPointerVector<Condition>& neighb_cond = (*rthis)->GetValue(NEIGHBOUR_CONDITIONS);
		    if(neighb_cond.size()!=0)
		    {  
		     for(WeakPointerVector<Condition>::iterator neighb = neighb_cond.begin();  neighb!= neighb_cond.end(); neighb++)
		      {   
		        Condition::GeometryType& geom_2 = (neighb)->GetGeometry();
			//KRATOS_WATCH((neighb)->Id())
			//KRATOS_WATCH(geom_2[0].Id())
			//KRATOS_WATCH(geom_2[1].Id())
			//KRATOS_WATCH(geom_2[2].Id())
			
			if((geom_2[0].Id() != geom[i].Id()) && (geom_2[1].Id() != geom[i].Id()) && (geom_2[2].Id()!= geom[i].Id()))
			{  
			  Points1[0]                      = geom_2[0].X(); 
			  Points1[1]                      = geom_2[0].Y();
			  Points1[2]                      = geom_2[0].Z();
			  
			  Points2[0]                      = geom_2[1].X();
			  Points2[1]                      = geom_2[1].Y(); 
			  Points2[2]                      = geom_2[1].Z(); 
			  
			  Points3[0]                      = geom_2[2].X();
			  Points3[1]                      = geom_2[2].Y(); 
			  Points3[2]                      = geom_2[2].Z(); 
			  
			  compare_distance = rPlane.DistPoint3Triangle3(Points0, Points1, Points2, Points3);
			  //KRATOS_WATCH(compare_distance)
			  if(compare_distance<distance) 
			    {  
				distance = compare_distance;
				geom[i].SetLock();
				geom[i].GetValue(CONTACT_LINK_MASTER) = neighb_cond(I).lock();
				geom[i].UnSetLock();
			    } 
			 }
			 I++;
	              } 
		     }
		    }
		   //KRATOS_WATCH(geom[i].Id())
		   //KRATOS_WATCH(Result.size())
	           //KRATOS_WATCH(geom[i].GetValue(CONTACT_LINK_MASTER))
	           //KRATOS_WATCH("-----------------") 
		     
	            }
	           }
	         }
	       }
	         Result.clear();   
	   }
	   
	  }
	   
	   //std::cout<< "     LOCAL SEARCH ALGORITHM   " <<std::endl; 
	   #ifdef _OPENMP
	   double stop_prod = omp_get_wtime();
	   std::cout << "        Time Searching  Masters Surfaces       = " << stop_prod - start_prod << "  seconds " << std::endl;
	   #endif
	   
	   
	   KRATOS_CATCH("")
       }
       
       
       
        inline void CreatePartition(unsigned int number_of_threads, const int number_of_rows, vector<unsigned int>& partitions)
        {
            partitions.resize(number_of_threads + 1);
            int partition_size = number_of_rows / number_of_threads;
            partitions[0] = 0;
            partitions[number_of_threads] = number_of_rows;
            for (unsigned int i = 1; i < number_of_threads; i++)
                partitions[i] = partitions[i - 1] + partition_size;
        }
       
       	  inline void V3DCro( double& x1,
			      double& y1,
			      double& z1,
			      const double& x2,
			      const double& y2,
			      const double& z2,
			      const double& x3,
			      const double& y3,
			      const double& z3)
              {
	        x1 =(y2)*(z3)-(z2)*(y3); 
                y1 =(z2)*(x3)-(x2)*(z3); 
                z1 =(x2)*(y3)-(y2)*(x3); 
	      }
	      
	                 
           inline void V3DNor( double& s,
			       double& x1,
			       double& y1,
			       double& z1)
              {
	         s= std::sqrt((x1)*(x1)+(y1)*(y1)+(z1)*(z1)); 
                 if((s)>EPSILON)(x1)=(x1)/(s);  
                 if((s)>EPSILON)(y1)=(y1)/(s);  
                 if((s)>EPSILON)(z1)=(z1)/(s);  
	      } 
	    
	    inline void V3DDot(double& s,
			       const double& x1,
			       const double& y1,
			       const double& z1,
			       const double& x2,
			       const double& y2,
			       const double& z2) 
	    {
	      s = ((x1)*(x2))+((y1)*(y2))+((z1)*(z2));
	    }

 
    };
	
	


}  // namespace Kratos.

#endif // KRATOS_GEOMETRY_UTILITIES_INCLUDED  defined 


