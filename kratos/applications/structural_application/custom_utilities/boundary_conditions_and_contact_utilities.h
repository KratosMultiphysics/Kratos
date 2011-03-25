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
//   Last Modified by:    $Author: nelson $
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
#include "geometries/line_2d_2.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/quadrilateral_2d_4.h"

#include "spatial_containers/spatial_containers.h"
#include "spatial_containers/bounding_box.h"
#include "spatial_containers/cell.h"
#include "spatial_containers/bins_dynamic_objects.h"

#include "utilities/spatial_containers_configure.h"
#include "utilities/geometry_utilities.h"
#include "utilities/timer.h"
// #include "utilities/timer_CLabra.h"
#include "custom_conditions/slave_contact_point_2d.h"
#include "custom_conditions/master_contact_point_2d.h"
#include "custom_conditions/master_contact_face_2d.h"
#include "custom_conditions/point_segment_contact_link.h"
#include "custom_conditions/point_point_contact_link.h"


namespace Kratos
{
  
///******************************************************************************************************************
///******************************************************************************************************************
 
 

  class BoundaryConditionsAndContactUtilities
	{
	public:
	  
	  
	enum Intersect{IT_POINT   = 0, IT_SEGMENT, IT_EMPTY}; 
        enum Exist_Node {no_nodes = 0, yes_nodes};
	enum Near_Node  {no_near  = 0, yes_near};
	enum Object     {is_node  = 0, is_object};
	
	class Segment2D
	{

	  public:
	  Segment2D(){}
	  Segment2D(array_1d<double,2 >& P0, array_1d<double,2 >& P1)
	  {
	  noalias(mP0) = P0;
	  noalias(mP1) = P1;
	  ComputeCenterDirectionExtent();
	  }

	  ~Segment2D(){}

	  double mExtent;
	  array_1d<double,2 > mP0;
	  array_1d<double,2 > mP1;
	  array_1d<double,2 > mCenter;
	  array_1d<double,2 > mDirection;

	  //--------------------------------------------------------------------------->
	  void ComputeCenterDirectionExtent ()
	  {
	  noalias(mCenter)    = (0.50) * (mP0 + mP1);
	  noalias(mDirection) = mP1 - mP0;
	  const double length = (std::sqrt(inner_prod(mDirection, mDirection) ) );
	  mExtent             =  0.500 * length;
	  mDirection          =  (1.00/(length ) )* mDirection;
	  }

	  //--------------------------------------------------------------------------->
	  inline double DotPerp (const array_1d<double,2 >& vec)
	  {
	     return mDirection[0]*vec[1] - mDirection[1]*vec[0];
	  }

	};
	  
	  
	    KRATOS_CLASS_POINTER_DEFINITION(BoundaryConditionsAndContactUtilities);
	    
	    
	     
	    /// Elements
	    typedef ModelPart::ElementsContainerType                     ElementsArrayType;
	    typedef ModelPart::ElementsContainerType::ContainerType      ContainerType; 
	    typedef ContainerType::value_type                            PointerType;
	    typedef ContainerType::iterator                              IteratorType; 
 	    typedef std::vector<PointerType>::iterator                   PointerTypeIterator;
 	    typedef ContactPair<PointerType>                             ContactPairType; 
 	    typedef std::vector<ContactPairType>                         ContainerContactPair;   
 	    typedef ContainerContactPair::iterator                       IteratorContainerContactPair;  
 	    typedef ContainerContactPair::value_type                     PointerContainerContactPair;  
	    
	
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
	    
	    
	    
            BoundaryConditionsAndContactUtilities(){}
            BoundaryConditionsAndContactUtilities(ModelPart& model_part, const unsigned int& dimension) : mr_model_part(model_part), mrdimension(dimension) 
              {  
		minitialize = false;
              }   
  
            virtual ~BoundaryConditionsAndContactUtilities(){}


            void  CreateBoundariesAndLinkingConditions()
                { 
		  KRATOS_TRY
		  
		  
		  Clear();
		  
		  /// Crea las conciones de contacto
		  //if (minitialize==false) 
		  CalculateBoundaryContour(mMasterConditionsArray);
  
		   
		  if(SearchContactsPairs()){ 
		     CreateLinkingConditions(); 
		  }
		  
	          KRATOS_CATCH("")
          }



       
	 private:

	   
         //************************************************************************************
	 //************************************************************************************   
         bool SearchContactsPairs()
         {
	  
	    KRATOS_TRY
	    
	    /// Configure 
	    std::cout<< std::endl;
	    std::cout<<"  COMPUTING CONTACT CONDITIONS TO MODEL PART " << std::endl; 
	    const std::size_t dim = 2;
	    typedef Point<dim , double>                   PointType;
	    typedef PointType::CoordinatesArrayType       CoordinatesArrayType;
	    typedef SpatialContainersConfigure<dim>       Configure;   

	    
	    IteratorType it_begin     =  mBoundaryElements.begin();
	    IteratorType it_end       =  mBoundaryElements.end(); 

	    
	    
	    BinsObjectDynamic<Configure>  rBinsObjectDynamic(it_begin, it_end ); 
	    rBinsObjectDynamic.SearchContact(mPairContacts);
	    FiltratePairContacts(mPairContacts);
	    
	    
	    if(mPairContacts.size()!=0)
	    {
	      std::cout<< "     NUMBER OF CONTACT PAIRS         = " <<mPairContacts.size()<<std::endl; 
	      return true;
	    }
	   
	    std::cout<< "     NO CONTACTS PAIRS "<<std::endl; 
	    return false;  
	    
	    KRATOS_CATCH("")
	 }
	    
	  	   
         //************************************************************************************
	 //************************************************************************************ 
	 
	 void Clear()
	 {
	   
	   KRATOS_TRY
     
	   NodesArrayType& pNodes = mr_model_part.Nodes(); 
	   for(ModelPart::NodeIterator i=pNodes.begin(); i!= pNodes.end(); i++){      
	       //if(i->FastGetSolutionStepValue(IS_BOUNDARY) == 1.00)
		   i->GetValue(NEIGHBOUR_CONDITIONS).clear();
	   }
	    
	    
	   ElementsArrayType& pElements = mr_model_part.Elements(); 
	   for(ModelPart::ElementIterator i=pElements.begin(); i!= pElements.end(); i++){      
	         i->GetValue(NEIGHBOUR_CONDITIONS).clear();
	   }
	    
	    
	    mBoundaryElements.clear();
            mPairContacts.clear();
            mpair.clear();
            mMasterConditionsArray.clear();     
	       
	   KRATOS_CATCH("") 
	       
	 }
	  
	  
	 //************************************************************************************
	 //************************************************************************************  
	    
	 void CreateLinkingConditions() 
	 {
	   
	    KRATOS_TRY
	    
	    ConditionsArrayType& rConditions = mr_model_part.Conditions();
	    ConditionsArrayType LinkingConditions;
	        
	    Exist_Node Exist =  no_nodes;
	    
	    unsigned int  master      = 0;
	    unsigned int  slave       = 1;
	    bool initialize           = false;
	    
	    std::vector<unsigned int>                InsideNodes;
	    std::vector<unsigned int>                TotalInsideNodes;
	    PointsArrayType                          Slaves;
	    std::vector<unsigned int>::iterator      repeated_object;  
	    
	    unsigned int Id               = rConditions.size() + 1;
	    unsigned int properties_index = mr_model_part.NumberOfProperties();
	    
	    PropertiesType::Pointer tempProperties = PropertiesType::Pointer(new PropertiesType(properties_index+1) );
	    mr_model_part.AddProperties( tempProperties );

	    
	    for(IteratorContainerContactPair it_pair = mPairContacts.begin(); it_pair!= mPairContacts.end(); it_pair++)
	        { 
		  
		   //Vericando los nodos que caen en el elemento.
		   master = 0;
		   slave  = 1;
		   NodeInside( (*it_pair)[master], (*it_pair)[slave], InsideNodes);
                   if(InsideNodes.size()==0) 
		     { 
		        master = 1;
		        slave  = 0;
		        InsideNodes.clear();
		        NodeInside( (*it_pair)[master], (*it_pair)[slave], InsideNodes);   
		     }
		   
		   
		   if(InsideNodes.size()!=0)
		      Exist = yes_nodes;
		   
		   std::cout<< "     MASTER OBJECT =  " <<  (*it_pair)[master]->Id() <<"   SLAVE OBJECT = " << (*it_pair)[slave]->Id() << std::endl; 
		   switch(Exist)
		   {
		    
		     case(yes_nodes):
		     {
		       std::cout<< "yes_nodes" << std::endl;  
		       array_1d<unsigned int, 2 >  Ids;       
		       Near_Node  Near  =  no_near;
		       bool corner      =  false; 
		       for(unsigned int in = 0; in<InsideNodes.size(); in++){
			 unsigned int& id  = InsideNodes[in]; 
		         Near   = CheckNearNodes(master, slave, mr_model_part.Nodes()(id), (*it_pair)[master], Ids);
			 corner = Is_Corner(mr_model_part.Nodes()(Ids[slave]), mr_model_part.Nodes()(Ids[master])); 
			 
			 if(Near==yes_near && corner==true){
			   std::cout<< "yes near " << std::endl;
			   // contact node - to - node  
			   CreatePointLinkingConditions(master, slave,  id, Ids, it_pair, tempProperties, initialize, Id, TotalInsideNodes, LinkingConditions); 
			 }
			 
			 else{
			     std::cout<< "no near " << std::endl;
			     //KRATOS_WATCH(Ids[master])
			     //KRATOS_WATCH(Ids[slave])
			     //KRATOS_WATCH(id)
			     // contact node - to - segment
			     CreateLinkingConditions(master, slave, Ids, Exist, it_pair, tempProperties,  initialize, Id, TotalInsideNodes, LinkingConditions );    
			    }
			 }
			 
			// intercambiando master por slave y viceversa 
			if(master==0){ 
			    std::cout<< "Changing " << std::endl;
			    InsideNodes.clear();
			    master = 1;
			    slave  = 0; 
			    NodeInside( (*it_pair)[master], (*it_pair)[slave], InsideNodes);
			    for(unsigned int in = 0; in<InsideNodes.size(); in++){
			    unsigned int& id  = InsideNodes[in];   
			    Near = CheckNearNodes(master, slave, mr_model_part.Nodes()(id), (*it_pair)[master], Ids);  
			    std::cout<< "no near " << std::endl;
			    //KRATOS_WATCH(Ids[master])
			    //KRATOS_WATCH(Ids[slave])
			    //unsigned int& id  = InsideNodes[in]; 
			    // contact node - to - segment
			    if(Near==no_near)
			        CreateLinkingConditions(master, slave, Ids, Exist,it_pair, tempProperties,  initialize, Id, TotalInsideNodes, LinkingConditions );    
			    }
			   }
		       
		       break;
		   }
		     
		     // Caso en que los triangulos se intersecten pero no hay nodo dentro de los elemetos
		     case(no_nodes):
		     {    
		        std::cout<< "no nodes " << std::endl;
			std::vector<array_1d<unsigned int, 2 > > Ids;
			std::vector<Near_Node> Is_Near;
			
		        CheckNearNodes(master, slave, (*it_pair)[slave], (*it_pair)[master], Ids, Is_Near);
		        for(unsigned int i = 0; i<Ids.size(); i++){   
			if(Is_Near[i]==yes_near){
			  std::cout<< "yes_near " << std::endl;
			  KRATOS_WATCH(Ids[i][master]) 
			  KRATOS_WATCH(Ids[i][slave])
			    // contact node - to - node  
			   CreatePointLinkingConditions(master, slave, Ids[i][slave], Ids[i], it_pair, tempProperties, initialize, Id, TotalInsideNodes, LinkingConditions); 
			 }
			 
			 
			 if(Is_Near[i]==no_near){  
			     std::cout<< "no_near " << std::endl;
			     KRATOS_WATCH(Ids[i][master]) 
			     KRATOS_WATCH(Ids[i][slave])
			     // contact node - to - segment
			     CreateLinkingConditions(master, slave, Ids[i], Exist, it_pair, tempProperties,  initialize, Id, TotalInsideNodes, LinkingConditions );    
			    }
			}
			
		        break;
		       }
		    }
		} 
		
		      std::cout<<"     NUMBER OF INITIAL CONDITIONS    = " << rConditions.size() <<  std::endl; 
		      //adding linking to model_part
		      for(ConditionsArrayType::ptr_iterator it=LinkingConditions.ptr_begin(); it != LinkingConditions.ptr_end(); ++it )
		            mr_model_part.Conditions().push_back( *it );
		     
		    
		      
		      std::cout<<"     NUMBER OF CALCULATED CONDITIONS = " << LinkingConditions.size() <<  std::endl; 
		      std::cout<<"     NUMBER OF FINAL CONDITIONS      = " << rConditions.size() <<  std::endl; 
		      
		      LinkingConditions.clear();
		      
		      
		      KRATOS_CATCH("")
	    }
	 
//*****************************************************************************************************
//*****************************************************************************************************

void CreatePointLinkingConditions(
	 const unsigned int& master,
         const unsigned int& slave, 
	 const unsigned int& InsideNodes,
	 const array_1d<unsigned int, 2 >& Ids,
	 const IteratorContainerContactPair& it_pair,
	 const PropertiesType::Pointer tempProperties, 
	 bool& initialize,
	 unsigned int& Id,
	 std::vector<unsigned int>& TotalInsideNodes,
	 ConditionsArrayType& LinkingConditions
	 )
	
	{
	  KRATOS_TRY
	  
	  const unsigned int& id              =  InsideNodes;
	  std::vector<unsigned int>::iterator    repeated_object;  

	  
	  
	  if(TotalInsideNodes.size()==0){   
	  TotalInsideNodes.push_back(id); 
	  repeated_object = TotalInsideNodes.end();
	  }

	  // verifcando que no se repitan los slaves nodes 
	  if(initialize==true)
	       repeated_object =  std::find(TotalInsideNodes.begin(), TotalInsideNodes.end(), id); 


	  if( repeated_object == (TotalInsideNodes.end())) 
	  {          
	     if(initialize==true)
	        TotalInsideNodes.push_back(id); 


	  // Slave Node
	  Point2D<Node<3> >::Pointer point_geom_slave    =  Point2D<Node<3> >::Pointer( new Point2D<Node<3> >(mr_model_part.Nodes()(Ids[slave]) ) );
	  Condition::Pointer SlaveNode                   =  Condition::Pointer(new SlaveContactPointType(Id, point_geom_slave) ); 
	 
	  // Master Node
	  Point2D<Node<3> >::Pointer point_geom_master   =  Point2D<Node<3> >::Pointer( new Point2D<Node<3> >(mr_model_part.Nodes()(Ids[master]) ) );
	  Condition::Pointer MasterNode                  =  Condition::Pointer(new MasterContactPointType(Id, point_geom_master) );

	  Condition::GeometryType& Mgeom   =  MasterNode->GetGeometry();
	  Condition::GeometryType& Sgeom   =  SlaveNode ->GetGeometry();  
	  Line2D2<Node<3> >::Pointer Lgeom =  Line2D2<Node<3> >::Pointer( new Line2D2<Node<3> >(Sgeom(0), Mgeom(0) ) );
	  Condition::Pointer newLink       =  Condition::Pointer( new PointPointContactLink(Id,
	  Lgeom,
	  tempProperties, 
	  SlaveNode,
	  MasterNode
	  ) );


	  LinkingConditions.push_back(newLink);    
	  Id++;  

	  if(initialize==false) 
	     initialize=true;
	  }
	  
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
	   std::vector<array_1d<unsigned int, 2 > >&  Ids,
	   std::vector<Near_Node>& Is_Near
	 )
	 {   
	   
	     
	    std::vector<double>           Distance;
	    std::vector<double>::iterator it;
	    std::vector<double>::iterator it2;
	    array_1d<unsigned int, 2 > Id;
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
	          distance = std::sqrt(inner_prod(vector, vector));
	          Distance.push_back(distance);
	       }
	    }

            // verificando si los nodos son muy cercanos  
            double ratio       = 0.00;
            const double max   = (*std::max_element(Distance.begin(), Distance.end() ) );  
	    for(unsigned int i = 0; i<Distance.size(); i++){
	       ratio =  std::fabs(Distance[i]/max);
	       if(ratio < 1E-6){
		  Id[master] = geom_0[M[i]].Id();   
		  Id[slave]  = geom_1[S[i]].Id();
		  Ids.push_back(Id); 
	          Is_Near.push_back(yes_near);
	       }
	    }
	    
	    
	    // si no se cumple lo anterior tomamos los nodos mas cercanos
	    // WARNING = Solo valido para un caso en que un solo nodo quede fuera
	    const double min   = (*std::min_element(Distance.begin(), Distance.end() ) );
	    it                 = std::find(Distance.begin(), Distance.end(), min);
	    const int position = int(it-Distance.begin()); 
	    Id[master]         = geom_0[M[position]].Id();   
	    Id[slave]          = geom_1[S[position]].Id();
	    bool is_repited    = false; 
	    for(unsigned int i = 0; i<Ids.size(); i++){
	      if(Ids[i][master]== Id[master] && Ids[i][slave]== Id[slave] )
	           is_repited = true;
	           break;
	    }
	    
	    if(is_repited==false ){
	      Ids.push_back(Id); 
	      Is_Near.push_back(no_near);
	    }
	    
	 }

	
//*****************************************************************************************************
//*****************************************************************************************************
	 
	 Near_Node CheckNearNodes( 
	   const unsigned int&  master,
	   const unsigned int&  slave,   
	   const NodePointerType& SlaveNode,
	   const PointerType& MasterObject,
	   array_1d<unsigned int, 2 >&  Ids
	 )
	 
	 {
	    unsigned int               Id_master;
	    std::vector<double>        Distance;
	    array_1d<double, 3>        vector;
	    array_1d<double, 3>        coordinates =  SlaveNode->Coordinates();
	    const Element::GeometryType& geom_0    =  MasterObject->GetGeometry();
	    double distance                        = 0.00;
	    double distance2                       = 1E10;;

	    
	    // busco la distancia menor
	    for(unsigned int i = 0; i<geom_0.size(); i++){
	           noalias(vector) = ( geom_0[i]-  coordinates);
	           distance = std::sqrt(inner_prod(vector, vector) );
		   Distance.push_back(distance);
	           if(distance<distance2){
	               distance2    =  distance;
	               Id_master    =  geom_0[i].Id();
	            }
	       
	    }
	    
	    
	    Ids[master]  = Id_master;
	    Ids[slave]   = SlaveNode->Id();
	  
    
            double max   = (*std::max_element(Distance.begin(), Distance.end() ) );  
	    double min   = (*std::min_element(Distance.begin(), Distance.end() ) );
	    double ratio = std::fabs(min/max); 
	    
	    if(ratio < 1E-3)
                return yes_near;
	     
	    return no_near; 
	    
	 }
	 
	 
//*****************************************************************************************************
//*****************************************************************************************************
	 
	 
	 void CreateLinkingConditions(
	 const unsigned int& master,
         const unsigned int& slave,
         const array_1d<unsigned int, 2 >& Ids,   
	 Exist_Node& Exist,
	 const IteratorContainerContactPair& it_pair,
	 const PropertiesType::Pointer tempProperties, 
	 bool& initialize,
	 unsigned int& Id,
	 std::vector<unsigned int>& TotalInsideNodes,
	 ConditionsArrayType& LinkingConditions
	 )
	 {
	    Object What_Is;
	    const unsigned int& InsideNodes = Ids[slave];
	    unsigned int segment = 0;
	    std::vector<unsigned int>::iterator    repeated_object;  
	    //PropertiesType::Pointer pProperties = mr_model_part.pGetProperties(1);   
	    

	    if(TotalInsideNodes.size()==0){   
	    TotalInsideNodes.push_back(InsideNodes); 
	    repeated_object = TotalInsideNodes.end();
	    }
            
	    // verifcando que no se repitan los slaves nodes 
	    if(initialize==true)
	        repeated_object =  std::find(TotalInsideNodes.begin(), TotalInsideNodes.end(), InsideNodes); 

	    //WARNING = No necesariamente para otros contenedores la comparacion se hace con end -1 
	    if( repeated_object == (TotalInsideNodes.end())) 
	    {          

	    if(initialize==true)
	        TotalInsideNodes.push_back(InsideNodes);  
	   
	      
	    bool exist_segment =  LocateMasterSegment(segment, Exist, What_Is, mr_model_part.Nodes()(Ids[slave]), mr_model_part.Nodes()(Ids[master]), (*it_pair)[master]);

	    if(exist_segment==true)
	    {
		// Slave Node
		Point2D<Node<3> >::Pointer point_geom     =  Point2D<Node<3> >::Pointer( new Point2D<Node<3> >(mr_model_part.Nodes()(InsideNodes) ) );
		Condition::Pointer SlaveNode              =  Condition::Pointer(new SlaveContactPointType(Id, point_geom) ); 
		 
		switch(What_Is)
		{
		  case ( is_node):
		  {
		   WeakPointerVector<Condition>& neighb_cond            =  mr_model_part.Nodes()(Ids[master])->GetValue(NEIGHBOUR_CONDITIONS);  
		   Condition::Pointer MasterFace                        =  (neighb_cond(segment).lock());  
		   Condition::GeometryType& Mgeom                       =  MasterFace->GetGeometry();
		   Condition::GeometryType& Sgeom                       =  SlaveNode->GetGeometry();   
		   Triangle2D3<Node<3> >::Pointer Lgeom                 =  Triangle2D3<Node<3> >::Pointer( new Triangle2D3<Node<3> >( Sgeom(0), Mgeom(0), Mgeom(1) ) );

		   
		   KRATOS_WATCH( "is node" )
		   KRATOS_WATCH(Sgeom[0].Id() )  
		   KRATOS_WATCH(Mgeom[0].Id() )
		   KRATOS_WATCH(Mgeom[1].Id() )  
		   KRATOS_WATCH("---------------------------")  
		   
		   Condition::Pointer newLink  = Condition::Pointer( new PointSegmentContactLink(Id,
		    Lgeom,
		    tempProperties,
		    MasterFace, 
		    SlaveNode) );
		   LinkingConditions.push_back( newLink );    
		   Id++;
		break;
		}
	        case(is_object):
		  {
	           WeakPointerVector<Condition>& neighb_cond =  (*it_pair)[master]->GetValue(NEIGHBOUR_CONDITIONS);
		   Condition::Pointer MasterFace             =  (neighb_cond(segment).lock());  
		   Condition::GeometryType& Mgeom            =  MasterFace->GetGeometry();
		   Condition::GeometryType& Sgeom            =  SlaveNode->GetGeometry();   
		   Triangle2D3<Node<3> >::Pointer Lgeom      =  Triangle2D3<Node<3> >::Pointer( new Triangle2D3<Node<3> >( Sgeom(0), Mgeom(0), Mgeom(1) ) );

		   KRATOS_WATCH( "is object" ) 
		   KRATOS_WATCH(Sgeom[0].Id() )  
		   KRATOS_WATCH(Mgeom[0].Id() )
		   KRATOS_WATCH(Mgeom[1].Id() )  
		   KRATOS_WATCH("---------------------------")
		   
		   Condition::Pointer newLink  = Condition::Pointer( new PointSegmentContactLink(Id,
		    Lgeom,
		    tempProperties,
		    MasterFace, 
		    SlaveNode) );
		    LinkingConditions.push_back( newLink );    
		    Id++;
		    break;
		}
	      }
	    }
	    
	    }

	    if(initialize==false) 
	        initialize=true;
	          
	    }
	 
 	
//*****************************************************************************************************
//*****************************************************************************************************

bool LocateMasterSegment(unsigned int& segmento,
			 Exist_Node& Exist,
			 Object& What_Is,
			 const NodePointerType& SlaveNode,
			 const NodePointerType& MasterNode, // the most near
			 const PointerType& MasterObject)
 	 {   
	    KRATOS_TRY
	   
	     WeakPointerVector<Condition>& neighb_cond = MasterObject->GetValue(NEIGHBOUR_CONDITIONS);
	     
	     segmento             = 0;               
	     switch(Exist)
	     { 
	       case(yes_nodes):
	       {
		 std::cout<< "TEST A " << std::endl; 
	       //el elemento tiene una sola condicion master
	      if(neighb_cond.size()==1)
	      {  What_Is = is_object;
	         segmento = 0;
	         return true;
	      } 
		 
	       std::cout<< "TEST B " << std::endl; 
	       if(Test_One_B(segmento, SlaveNode, MasterNode)) 
	          { 
		    What_Is = is_node;
		    return true;
		  }
	       std::cout<< "TEST C " << std::endl;	  
	       if(Test_Three(segmento, SlaveNode, MasterObject)) 
	          {
		    What_Is = is_object;
		    return true;
		  }
	       if(Test_Five(segmento, SlaveNode, MasterObject) )
	         {
		   What_Is = is_object;
		   return true;
	         }
	         
		  break;
	       }
	       
	      case(no_nodes):
	       {
		 
		//con desplazamientos 
		std::cout<< "TEST A " << std::endl; 
		if(Test_One_A(segmento, SlaveNode, MasterNode)) 
	          { 
		    What_Is = is_node;
		    return true;
		  }
		  
		 // con aristas
		std::cout<< "TEST B " << std::endl; 
		if(Test_One_C(segmento, SlaveNode, MasterNode)) 
	          { 
		    What_Is = is_node;
		    return true;
		  }
		  
	       std::cout<< "TEST C " << std::endl; 	  
	       if(Test_Two(segmento, SlaveNode, MasterObject)) 
	          {
		    What_Is = is_object;
		    return true;
		  }
		  
		 if(Test_Four(segmento, SlaveNode, MasterObject) )
	         {
		   What_Is = is_object;
		   return true;
	         }
	         
		  break;
	       }  
	    
	     }
	     
	     
	    segmento = 0; 
	    What_Is  = is_object;
            return false;
	    
	    
	    KRATOS_CATCH("")
   }
   
   
// Si el nodo esta dentro de elemento   
bool Is_Corner(
               const NodePointerType& SlaveNode,
	       const NodePointerType& MasterNode)

{

  KRATOS_TRY   
  
  WeakPointerVector<Condition>& neighb_cond_master  = MasterNode->GetValue(NEIGHBOUR_CONDITIONS);
  WeakPointerVector<Condition>& neighb_cond_slave   = SlaveNode->GetValue(NEIGHBOUR_CONDITIONS);
 
  std::vector<unsigned int>           segment;
  std::vector<unsigned int>::iterator it;
  vector<array_1d<double, 2> >        Points0;
  vector<array_1d<double, 2> >        Points1;
  array_1d<double, 2>                 Point;

  unsigned int I   = 0;
  unsigned int II  = 1;
  unsigned int III = 1;

  
  for(WeakPointerVector<Condition>::iterator cond_slave  = neighb_cond_slave.begin(); cond_slave!= neighb_cond_slave.end(); ++cond_slave)
  {
      Condition::GeometryType& geom = cond_slave->GetGeometry();
      Point.resize(2, false); 
      Points0.resize(2, false); 
      Points1.resize(2, false);

      Points0(0)[0] = geom[0].X(); 
      Points0(0)[1] = geom[0].Y();
      Points0(1)[0] = geom[1].X();
      Points0(1)[1] = geom[1].Y(); 

      I  = 0; 
      for(WeakPointerVector< Condition >::iterator cond  = neighb_cond_master.begin(); cond!= neighb_cond_master.end(); ++cond){
        Condition::GeometryType& geom_3 = cond->GetGeometry();     
	Points1(0)[0] = geom_3[0].X(); 
	Points1(0)[1] = geom_3[0].Y();
	Points1(1)[0] = geom_3[1].X();
	Points1(1)[1] = geom_3[1].Y();

	if(IntersectSegment(Point, Points0, Points1)!=IT_EMPTY) 
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

   if(segment.size()==1)
      return false;
   
   return true;

   
   
   KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************ 

/// comparacion con desplazamientos
bool Test_One_A(unsigned int& segmento,
	       const NodePointerType& SlaveNode,
	       const NodePointerType& MasterNode) 
{
  
  KRATOS_TRY   
  
  WeakPointerVector<Condition>& neighb_cond_master     = MasterNode->GetValue(NEIGHBOUR_CONDITIONS);

  std::vector<unsigned int> segment;
  unsigned int I       = 0;
  segmento             = 0;

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

    if(IntersectSegment(Point, Points0, Points1)==IT_POINT){
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
     return true;
  }
  
  return false;
  KRATOS_CATCH("")
}



/// comparacion con desplazamientos
bool Test_One_C(unsigned int& segmento,
	       const NodePointerType& SlaveNode,
	       const NodePointerType& MasterNode) 
{
  
  KRATOS_TRY
  WeakPointerVector<Condition>& neighb_cond            = MasterNode->GetValue(NEIGHBOUR_CONDITIONS);
  WeakPointerVector<Condition>& neighb_cond_slave      = SlaveNode->GetValue(NEIGHBOUR_CONDITIONS);
  
  std::vector<unsigned int>         segment;
  std::vector<array_1d<double, 2> > Points; // punto de interseccion del segmento
  vector<array_1d<double, 2> >      Points0;
  vector<array_1d<double, 2> >      Points1;
  array_1d<double, 2>               Point;

  // test with edges
  segmento         = 0;
  unsigned int I   = 0;
  unsigned int II  = 1;
  unsigned int III = 1;

  
  for(WeakPointerVector<Condition>::iterator cond_slave  = neighb_cond_slave.begin(); cond_slave!= neighb_cond.end(); ++cond_slave)
  {
      Condition::GeometryType& geom = cond_slave->GetGeometry();
      Point.resize(2, false); 
      Points0.resize(2, false); 
      Points1.resize(2, false);
      
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
	
	if(IntersectSegment(Point, Points0, Points1)==IT_POINT)
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
  return true;
  }
  
 
  return false;
  
  KRATOS_CATCH("")
}


/// comparacion con desplazamientos
bool Test_One_B(unsigned int& segmento,
	       const NodePointerType& SlaveNode,
	       const NodePointerType& MasterNode) 
{
  
  KRATOS_TRY   
  
  WeakPointerVector<Condition>& neighb_cond_master     = MasterNode->GetValue(NEIGHBOUR_CONDITIONS);

  std::vector<unsigned int> segment;
  unsigned int I       = 0;
  segmento             = 0;

     
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


  Points0.resize(2, false); 
  Points1.resize(2, false);

  Points0(0)[0] = SlaveNode->X();
  Points0(0)[1] = SlaveNode->Y(); 
  Points0(1)[0] = SlaveNode->X(); 
  Points0(1)[1] = SlaveNode->Y(); 

 
  unsigned int JJ = 1; 
  for(WeakPointerVector< Condition >::iterator cond  = neighb_cond_master.begin(); cond!= neighb_cond_master.end(); cond++){
  Condition::GeometryType& geom_2 = cond->GetGeometry();
	    
    Points1(0)[0] = geom_2[0].X(); 
    Points1(0)[1] = geom_2[0].Y();
    Points1(1)[0] = geom_2[1].X();
    Points1(1)[1] = geom_2[1].Y();

    if(IntersectSegment(Point, Points0, Points1)==IT_POINT){
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
     return true;
  }
  
  return false;
  KRATOS_CATCH("")
}


/// para cercanos
/// caso en que las aristas estan fuera

bool Test_Two( unsigned int& segmento, 
	       const NodePointerType& SlaveNode,
	       const PointerType& MasterObject) 
{
  
  KRATOS_TRY   
  
  WeakPointerVector<Condition>& neighb_cond     = MasterObject->GetValue(NEIGHBOUR_CONDITIONS);

  std::vector<unsigned int> segment;
  unsigned int I       = 0;
  segmento             = 0;

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

    if(IntersectSegment(Point, Points0, Points1)!=IT_EMPTY){
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
     return true;
  }
  
  return false;
  KRATOS_CATCH("")
}


/// para lejanos
/// caso en que las aristas esten dentro de un elemento
bool Test_Three( unsigned int& segmento, 
	       const NodePointerType& SlaveNode,
	       const PointerType& MasterObject) 
{
  
  KRATOS_TRY   
  
  WeakPointerVector<Condition>& neighb_cond     = MasterObject->GetValue(NEIGHBOUR_CONDITIONS);

  std::vector<unsigned int> segment;
  unsigned int I       = 0;
  segmento             = 0;

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

    if(IntersectSegment(Point, Points0, Points1)!=IT_EMPTY){
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
     return true;
  }
  
  return false;
  KRATOS_CATCH("")
}


/// para cercanos
/// caso en que las aristas esten fuera de un elemento
bool Test_Four(unsigned int& segmento,
	      const NodePointerType& SlaveNode,
	      const PointerType& MasterObject)
{
 
  KRATOS_TRY   
  
  WeakPointerVector<Condition>& neighb_cond            = MasterObject->GetValue(NEIGHBOUR_CONDITIONS);
  WeakPointerVector<Condition>& neighb_cond_slave      = SlaveNode->GetValue(NEIGHBOUR_CONDITIONS);
  array_1d<double,3>& old_pos                          = SlaveNode->FastGetSolutionStepValue(DISPLACEMENT,3);   
  std::vector<unsigned int>         segment;
  std::vector<array_1d<double, 2> > Points; // punto de interseccion del segmento
  vector<array_1d<double, 2> >      Points0;
  vector<array_1d<double, 2> >      Points1;
  array_1d<double, 2>               Point;

  // test with edges
  segmento         = 0;
  unsigned int I   = 0;
  unsigned int II  = 1;
  unsigned int III = 1;

  
  for(WeakPointerVector<Condition>::iterator cond_slave  = neighb_cond_slave.begin(); cond_slave!= neighb_cond.end(); ++cond_slave)
  {
      Condition::GeometryType& geom = cond_slave->GetGeometry();
      Point.resize(2, false); 
      Points0.resize(2, false); 
      Points1.resize(2, false);

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

	if(IntersectSegment(Point, Points0, Points1)==IT_POINT)
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
  return true;
  }
  
 
  return false;
  
  KRATOS_CATCH("")
  
}



/// para cercanos
/// caso en que las aristas esten fuera de un elemento
bool Test_Five(unsigned int& segmento,
	      const NodePointerType& SlaveNode,
	      const PointerType& MasterObject)
{
 
  KRATOS_TRY   
  
  WeakPointerVector<Condition>& neighb_cond            = MasterObject->GetValue(NEIGHBOUR_CONDITIONS);
  WeakPointerVector<Condition>& neighb_cond_slave      = SlaveNode->GetValue(NEIGHBOUR_CONDITIONS);
  array_1d<double,3>& old_pos                          = SlaveNode->FastGetSolutionStepValue(DISPLACEMENT,3);   
  std::vector<unsigned int>         segment;
  std::vector<array_1d<double, 2> > Points; // punto de interseccion del segmento
  vector<array_1d<double, 2> >      Points0;
  vector<array_1d<double, 2> >      Points1;
  array_1d<double, 2>               Point;

  // test with edges
  segmento         = 0;
  unsigned int I   = 0;
  unsigned int II  = 1;
  unsigned int III = 1;

  
  for(WeakPointerVector<Condition>::iterator cond_slave  = neighb_cond_slave.begin(); cond_slave!= neighb_cond.end(); ++cond_slave)
  {
      Condition::GeometryType& geom = cond_slave->GetGeometry();
      Point.resize(2, false); 
      Points0.resize(2, false); 
      Points1.resize(2, false);

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

	if(IntersectSegment(Point, Points0, Points1)==IT_POINT)
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
  return true;
  }
  
 
  return false;
  
  KRATOS_CATCH("")
  
}





	
//************************************************************************************
//************************************************************************************ 
void CalculateBoundaryContour(ConditionsArrayType& MasterConditions)
	 {   
	      
	   
	     KRATOS_TRY
	   
	      std::cout<<"     CALCULATING CONTOURS " <<  std::endl; 
	     
	      typedef WeakPointerVector< Element >::iterator  ElementIteratorType; 
	      ContainerType& rElements           =  mr_model_part.ElementsArray();
	      IteratorType it_begin              =  rElements.begin();
	      IteratorType it_end                =  rElements.end(); 
	      
	      array_1d<unsigned int,2>  Pair; 
	      bool is_boundary   = false;
	      
	      unsigned int face     = 0; 
	      unsigned int Id       = rElements.size() + 1 ;
	      for(IteratorType elem = it_begin; elem!=it_end; elem++)
	      {
		   Element::GeometryType& geom_1 = (*elem)->GetGeometry();
		   WeakPointerVector< Element >& neighb_elems    = (*elem)->GetValue(NEIGHBOUR_ELEMENTS); 
		   WeakPointerVector< Condition >& neighb_cond   = (*elem)->GetValue(NEIGHBOUR_CONDITIONS);
		   neighb_cond.clear();
		   
		   
		   //node_boundary.resize(neighb_elems.size(), false);
		   // Puede incluir como vecnino el mismo en caso de que hayan menos de 3 elemtos veninos.  
		   // ckeck si se repited elmento
		   // ElementIteratorType no necesita especificarse el * 
		   for( ElementIteratorType neighb_elem  = neighb_elems.begin(); neighb_elem!= neighb_elems.end(); neighb_elem++)
 		            {
			       
			       if ( neighb_elem->Id() ==  (*elem)->Id() )
			           {     
				     if(face == 0) // edge 1-2
				      {
					Pair[0]   =  geom_1[1].Id();
					Pair[1]   =  geom_1[2].Id();
					CreateMasterConditions(Pair, elem, Id, MasterConditions);
		                        geom_1[1].FastGetSolutionStepValue(IS_BOUNDARY) = 1.00;
					geom_1[2].FastGetSolutionStepValue(IS_BOUNDARY) = 1.00;
			              }
			              
				      if (face==1)  // edge 2-0
				      {
					 Pair[0] =   geom_1[2].Id();
					 Pair[1] =   geom_1[0].Id();
					 CreateMasterConditions(Pair, elem, Id, MasterConditions);
					 geom_1[2].FastGetSolutionStepValue(IS_BOUNDARY) = 1.00;
					 geom_1[0].FastGetSolutionStepValue(IS_BOUNDARY) = 1.00;
					 
				      }
				      if (face==2) // edge 0-1
				      {
					 Pair[0] =   geom_1[0].Id();
					 Pair[1] =   geom_1[1].Id();
					 CreateMasterConditions(Pair, elem, Id, MasterConditions);
					 geom_1[0].FastGetSolutionStepValue(IS_BOUNDARY) = 1.00;
					 geom_1[1].FastGetSolutionStepValue(IS_BOUNDARY) = 1.00;
				      }
				      
				       if (is_boundary==false)   
				       { 
					 is_boundary = true;
					 Id++;
					 mBoundaryElements.push_back(*elem);
				       }
				   }
				   
			       face++;    
			    }
			  
			is_boundary = false;  
		        face = 0;  
	           }
	           
	           minitialize = true;
		   
		   KRATOS_CATCH("")
	 }
	 

//*****************************************************************************************************
//*****************************************************************************************************
void CreateMasterConditions(array_1d<unsigned int,2>&  rPair, const IteratorType& elem, unsigned int& Id, ConditionsArrayType& MasterConditions)
	 {
	    KRATOS_TRY
	   
	    Line2D2<Node<3> >::Pointer pgeom =  Line2D2<Node<3> >::Pointer (new Line2D2<Node<3> >(mr_model_part.Nodes()((rPair)[0]), mr_model_part.Nodes()((rPair)[1]) ) ) ;  
	    Condition::Pointer MasterSegment = Condition::Pointer(new MasterContactFaceType(Id, pgeom ) ) ;
	    MasterSegment->GetValue(NEIGHBOUR_ELEMENTS).push_back(*elem);
	    
	    mr_model_part.Nodes()((rPair)[0])->GetValue(NEIGHBOUR_CONDITIONS).push_back(MasterSegment); 
	    mr_model_part.Nodes()((rPair)[1])->GetValue(NEIGHBOUR_CONDITIONS).push_back(MasterSegment); 
	    
	    (*elem)->GetValue(NEIGHBOUR_CONDITIONS).push_back(MasterSegment); 
	    MasterConditions.push_back(MasterSegment);
	    Id++;
	    
	    KRATOS_CATCH("")
	 }



//*****************************************************************************************************
//*****************************************************************************************************
void FiltratePairContacts(ContainerContactPair& PairContacts)
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
	  if( id.size()!=0) return true;
	  
	  return false;    
	  
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

void NodeInside(PointerType& MasterObject, PointerType& SlaveObject,  std::vector<unsigned int>& InsideNodes)
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
		      InsideNodes.push_back(geom_slave[Nodes[i]].Id()); }   
	   }
	   
	   
	   KRATOS_CATCH("")
	      
	 }


      
       private:
       ModelPart mr_model_part; 
       unsigned int mrdimension;
       bool minitialize;
       
       ContainerType             mBoundaryElements;
       ContainerContactPair      mPairContacts;
       std::vector<unsigned int> mpair;
       ConditionsArrayType       mMasterConditionsArray;
       
       
//*****************************************************************************************************
//*****************************************************************************************************  

       Intersect IntersectSegment(
       array_1d<double, 2>& Point, 
       vector<array_1d<double, 2> >& Points0,
       vector<array_1d<double, 2> >& Points1)
      {
	
	KRATOS_TRY
	
	double toler = 1E-14;
	Point        = ZeroVector(2);
	array_1d<double, 2> parameter = ZeroVector(2);
        
	
	Segment2D Segment0(Points0[0], Points0[1]);
	Segment2D Segment1(Points1[0], Points1[1]);
	Intersect IntersectionType = Classify(parameter, Segment0,  Segment1);
		
	if (IntersectionType == IT_POINT)
	{
	   // Test whether the line-line intersection is on the segments.
	   double a = std::fabs(parameter[0]) - Segment0.mExtent;
	   double b = std::fabs(parameter[1]) - Segment1.mExtent;
	   
	   if ( a<=toler && b<=toler)
	    {
	      Point    = Segment0.mCenter + parameter[0]*Segment0.mDirection;
	      //comprobando que el punto no sea el extremo del segmento
	      array_1d<double, 4 > aa; 
	      aa[0] = std::fabs(Points0(0)[0]- Point[0]);
	      aa[1] = std::fabs(Points0(0)[1]- Point[1]);
	      aa[2] = std::fabs(Points1(1)[0]- Point[0]);
	      aa[3] = std::fabs(Points1(1)[1]- Point[1]);

	      if( (aa[0]<toler && aa[1]<toler) || (aa[2]<toler && aa[3]<toler)){
		IntersectionType = IT_EMPTY;
	      }
	    }
	    else
	     IntersectionType = IT_EMPTY;
	 }
	 

	return IntersectionType;  // != IT_EMPTY;
	
	KRATOS_CATCH("")

      }
          

//*****************************************************************************************************
//*****************************************************************************************************  
          
      Intersect Classify(
      array_1d<double, 2>& s,
      Segment2D& Segment0,
      Segment2D& Segment1
      )
      
      {
	
	KRATOS_TRY
     
	// The intersection of two lines is a solution to P0+s0*D0 = P1+s1*D1.
	// Rewrite this as s0*D0 - s1*D1 = P1 - P0 = Q.  If D0.Dot(Perp(D1)) = 0,
	// the lines are parallel.  Additionally, if Q.Dot(Perp(D1)) = 0, the
	// lines are the same.  If D0.Dot(Perp(D1)) is not zero, then
	// s0 = Q.Dot(Perp(D1))/D0.Dot(Perp(D1))
	// produces the point of intersection.  Also,
	// s1 = Q.Dot(Perp(D0))/D0.Dot(Perp(D1))

        double toler                   = 1E-12;   
	array_1d<double, 2> originDiff = Segment1.mCenter - Segment0.mCenter;
	array_1d<double, 2> diff       = originDiff;
	double D0DotPerpD1             = Segment0.DotPerp(Segment1.mDirection);

	if ( std::fabs(D0DotPerpD1) > toler)
	{
	  // Lines intersect in a single point.
	  double invD0DotPerpD1 = 1.00/D0DotPerpD1;
	  double diffDotPerpD0  = originDiff[0]*Segment0.mDirection[1] - originDiff[1]*Segment0.mDirection[0]; 
	  double diffDotPerpD1  = originDiff[0]*Segment1.mDirection[1] - originDiff[1]*Segment1.mDirection[0]; 
	  s[0] = diffDotPerpD1*invD0DotPerpD1;
	  s[1] = diffDotPerpD0*invD0DotPerpD1;

	  return IT_POINT;
	}


	// Lines are parallel.
	originDiff = originDiff * (1.00 / ( std::sqrt(inner_prod(originDiff, originDiff) ) ) );    
	array_1d<double, 2> diffN = originDiff;
	

	double diffNDotPerpD1 = originDiff[0]*Segment1.mDirection[1] - originDiff[1]*Segment1.mDirection[0]; 
	if (std::fabs(diffNDotPerpD1) <= toler)
	{
	// Lines are colinear.	
	return IT_SEGMENT;
	}


	// Lines are parallel, but distinct.
	return IT_EMPTY;
	
	KRATOS_CATCH("")
	
      }

 

 
	};
	
	


}  // namespace Kratos.

#endif // KRATOS_GEOMETRY_UTILITIES_INCLUDED  defined 


