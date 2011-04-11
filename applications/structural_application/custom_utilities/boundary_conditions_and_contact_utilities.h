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
#include "geometries/line_2d_2.h"
#include "geometries/triangle_2d_3.h"


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


#include "custom_utilities/segment_2d.h"
#include "custom_utilities/intersect_triangles_cases.h"


namespace Kratos
{
   
  //template<std::size_t rdimension> 
  class BoundaryConditionsAndContactUtilities
	{
	public:
	   
	    
	    //enum Intersect{IT_POINT   = 0, IT_SEGMENT, IT_EMPTY};
	    static const int IT_POINT   = 0; 
	    static const int IT_SEGMENT = 1;
	    static const int IT_EMPTY   = 2;
	    
	    enum Exist_Node {no_nodes = 0, yes_nodes};
	    enum Near_Node  {no_near  = 0, yes_near};
	    enum Object     {is_node  = 0, is_object};

	    

	    KRATOS_CLASS_POINTER_DEFINITION(BoundaryConditionsAndContactUtilities);

	    /// Utilities
	    typedef IntersectionSegment2DToSegment2D  IntersectionSegments;
	    //enum IntersectionSegments::Intersect{IT_POINT = 0, IT_SEGMENT, IT_EMPTY};

     
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
	    
	    
	     
	    static const std::size_t space_dim = 2;
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
	    
	    typedef ContainerContactType  ContainerContactPair;   
 	    typedef IteratorContactType   IteratorContainerContactPair;  
 	    typedef PointerContactType    PointerContainerContactPair;  
	         
	    
	    
	    //typedef Point<space_dim , double>                             PointType;
	    //typedef PointType::CoordinatesArrayType                       CoordinatesArrayType;
	    

	    
	    
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
	    
	    IteratorType it_begin     =  mBoundaryElements.begin();
	    IteratorType it_end       =  mBoundaryElements.end(); 

	    
	    BinsObjectDynamic<Configure>  rBinsObjectDynamic(it_begin, it_end ); 
	    rBinsObjectDynamic.SearchContact(mPairContacts);
	    FiltratePairContacts(mPairContacts);
	    //mBinsObjectDynamic = rBinsObjectDynamic;
	    
	    
	    
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
		   i->GetValue(IS_CONTACT_SLAVE)  = 0;
		   i->GetValue(IS_CONTACT_MASTER) = 0;
		   
	   }
	    
	    
	   ElementsArrayType& pElements = mr_model_part.Elements(); 
	   for(ModelPart::ElementIterator i=pElements.begin(); i!= pElements.end(); i++){      
	         i->GetValue(NEIGHBOUR_CONDITIONS).clear();
	   }
	    
	    
	    mBoundaryElements.clear();
            mPairContacts.clear();
            //mpair.clear();
            mMasterConditionsArray.clear();     
	       
	   KRATOS_CATCH("") 
	       
	 }
	  
	  
	  //  caso en que un nodo este dentro de un elemento y al la vez fuera de otro
	  ///WARNING = NOT FINISHHHHHHHH YETTTTTTTt
	  bool CheckPenetrabilitySlaveNodeInOtherMasterElement(NodePointerType& SlaveNode, 
	                                                       PointerType& SlaveObject,
							       PointerType& MasterObject
	                                                       )
	  {
	      
	      std::vector<unsigned int> segment;
              unsigned int I                    = 0;
	      
              std::vector<array_1d<double, 2> > Points;   
              vector<array_1d<double, 2> >      Points0;
              vector<array_1d<double, 2> >      Points1;
              array_1d<double, 2>               Point;
	      ResultContainerType               Result;
	      
	      
              array_1d<double,3>& old_pos       = SlaveNode->FastGetSolutionStepValue(DISPLACEMENT,3);  
              //mBinsObjectDynamic.SearchObjectsInner(SlaveObject, Result);
           
	      // Solo contacta con uno 
	      if(Result.size()==1)
		 return true;
	      
	      // se busca la longitud de iterseccion mas larga.
	      // Si el master en cuestion es el correcto retunr true 
	      
	      else
	      {
	      Points0.resize(2, false); 
	      Points1.resize(2, false);

	      Points0(0)[0] = SlaveNode->X0() + old_pos[0];
	      Points0(0)[1] = SlaveNode->Y0() + old_pos[1];
	      Points0(1)[0] = SlaveNode->X(); 
	      Points0(1)[1] = SlaveNode->Y(); 
	      
	      
	      unsigned int JJ = 1;
	      for(ResultIteratorType it = Result.begin(); it!=Result.end(); it++)
	      {  
		JJ=1;
	        WeakPointerVector<Condition>& neighb_cond = (*it)->GetValue(NEIGHBOUR_CONDITIONS);
		
	        for(WeakPointerVector< Condition >::iterator cond  = neighb_cond.begin(); cond!= neighb_cond.end(); cond++){
	        Condition::GeometryType& geom_2 = cond->GetGeometry();

	         Points1(0)[0] = geom_2[0].X(); 
	         Points1(0)[1] = geom_2[0].Y();
	         Points1(1)[0] = geom_2[1].X();
	         Points1(1)[1] = geom_2[1].Y();

	         if(IntersectionSegments::IntersectSegment(Point, Points0, Points1)== IT_POINT){
	           Points.push_back(Point); 
	           segment.push_back(I);
	            }
	            
	         JJ++;
	         if(JJ>neighb_cond.size())
	            break;
	        }
	       }
	      }
	      return false;
	  }
	  
	 //************************************************************************************
	 //************************************************************************************  
	    
	 void CreateLinkingConditions() 
	 {
	   
	    KRATOS_TRY
	    
	    ConditionsArrayType& rConditions = mr_model_part.Conditions();
	    ConditionsArrayType LinkingConditions;
	        
	    Exist_Node Exist         =  no_nodes;
	    unsigned int  master      = 0;
	    unsigned int  slave       = 1;
	    bool initialize           = false;
	    
	    std::vector<unsigned int>                InsideNodes;
	    std::vector<unsigned int>                TotalInsideNodes;                
	    std::vector<unsigned int>::iterator      repeated_object;  
	    
	    unsigned int Id               = rConditions.size() + 1;
	    unsigned int properties_index = mr_model_part.NumberOfProperties();
	    
	    PropertiesType::Pointer tempProperties = PropertiesType::Pointer(new PropertiesType(properties_index+1) );
	    mr_model_part.AddProperties( tempProperties );

	    
	    array_1d<unsigned int, 2 >  Ids;       
	    Near_Node  Near             =  no_near;
	    bool is_repited             =  false;
	    bool corner                 =  false; 
	    bool Change                 =  true; 
	    int  Case                   =  0;     
	    unsigned int Id_Node_Case_5 =  0;
	    
	    //bool is_double_penetred =  false;
	    std::vector<array_1d<unsigned int, 2 > > Ids_2;
	    std::vector<Near_Node> Is_Near;
	    IntersectTriangleCases<Configure> IntersectTriangles(mr_model_part);  
	    
	    //#pragma omp parallel for 
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
		  {
		     unsigned int& id  = InsideNodes[0];  
		     Exist = yes_nodes;
		     if(InsideNodes.size()==1)
		         is_repited = (mr_model_part.Nodes()(id)->GetValue(IS_CONTACT_SLAVE)==1);
		  }
		  std::cout<< "     MASTER OBJECT =  " <<  (*it_pair)[master]->Id() <<"   SLAVE OBJECT = " << (*it_pair)[slave]->Id() << std::endl;
                  if(is_repited==false)
		  {
		    //std::cout<< "     MASTER OBJECT =  " <<  (*it_pair)[master]->Id() <<"   SLAVE OBJECT = " << (*it_pair)[slave]->Id() << std::endl;
//  		    for(unsigned int i = 0; i<InsideNodes.size(); i++)
//  		      std::cout<<InsideNodes[i]<<std::endl;
		    
		  switch(Exist)
		    {
		    case(yes_nodes):
		     {
		       //std::cout<< "Yes Nodes" << std::endl;
		       Case = IntersectTriangles.LocateCaseItersection(Id_Node_Case_5, Change, InsideNodes, (*it_pair)[master], (*it_pair)[slave]);
		       switch(Case) 
		       {
			 case 1: // un solo nodo dentro
			 {
			     unsigned int& id  = InsideNodes[0];  
			     Ids[master]       = 0;
			     Ids[slave]        = 0;
		             Near              = CheckNearNodes(master, slave, mr_model_part.Nodes()(id), (*it_pair)[master], Ids);
			     CreateLinkingConditions(master, slave, Ids, Exist, it_pair, tempProperties,  initialize, Id, TotalInsideNodes, LinkingConditions );  
			  break;
			 }
			 
			 case 2 :
			 {
			  for(unsigned int in = 0; in<InsideNodes.size(); in++){
			     unsigned int& id  = InsideNodes[in];  
			     Ids[master]       = 0;
			     Ids[slave]        = 0;
		             Near              = CheckNearNodes(master, slave, mr_model_part.Nodes()(id), (*it_pair)[master], Ids);
			     CreateLinkingConditions(master, slave, Ids, Exist, it_pair, tempProperties,  initialize, Id, TotalInsideNodes, LinkingConditions );  
			  }
			  break;
			 }
			 
			 case 3:
			 {
			   unsigned int& id  = InsideNodes[0];  
			   Ids[master]       = 0;
			   Ids[slave]        = 0;
		           Near              = CheckNearNodes(master, slave, mr_model_part.Nodes()(id), (*it_pair)[master], Ids);
			   CreatePointLinkingConditions(master, slave,  id, Ids, it_pair, tempProperties, initialize, Id, TotalInsideNodes, LinkingConditions);
			   break;
			 }  
			 
			 case 5:
			 {
			   unsigned int& id  = InsideNodes[0];  
			   Ids[master]       = 0;
			   Ids[slave]        = 0;
		           Near              = CheckNearNodes(master, slave, mr_model_part.Nodes()(id), (*it_pair)[master], Ids); 
			   CreateLinkingConditions(master, slave, Ids, Exist, it_pair, tempProperties,  initialize, Id, TotalInsideNodes, LinkingConditions );
			  
			   //WARNING = SPECIAL CASE
			   Exist             =  no_nodes;
			   unsigned int& id_2 =  Id_Node_Case_5;
			   Near              = CheckNearNodes(slave, master,  mr_model_part.Nodes()(id_2), (*it_pair)[slave], Ids);
			   CreateLinkingConditions(slave, master, Ids, Exist, it_pair, tempProperties,  initialize, Id, TotalInsideNodes, LinkingConditions );
			   break;
			 }
			 
		       }
			
		       break;
		     }
 
                     case(no_nodes):
		     {   
		        //std::cout<< "No Nodes" << std::endl;
		        CheckNearNodes(master, slave, (*it_pair)[slave], (*it_pair)[master], Ids_2, Is_Near);
		        for(unsigned int i = 0; i<Ids_2.size(); i++){   
			   if(Is_Near[i]==no_near){
			      CreateLinkingConditions(master, slave, Ids_2[i], Exist, it_pair, tempProperties,  initialize, Id, TotalInsideNodes, LinkingConditions );
			   }
		         }
		      break;
		    }
		    
		   }
		  }
		   
		// intercambiando master por slave y viceversa 
		  if(master==0 && Change==true)
		  { 
		    InsideNodes.clear();
		    master = 1;
		    slave  = 0; 
		    NodeInside( (*it_pair)[master], (*it_pair)[slave], InsideNodes);
		    
		   if(InsideNodes.size()!=0)
		   { 
		     unsigned int& id  = InsideNodes[0];  
		     Exist = yes_nodes;
		     if(InsideNodes.size()==1)
		         is_repited = (mr_model_part.Nodes()(id)->GetValue(IS_CONTACT_SLAVE)==1);
		   }
		   
		  if(is_repited==false)
		  {
		   //std::cout<< "Changing"<< std::endl;
		   //std::cout<< "     MASTER OBJECT =  " <<  (*it_pair)[master]->Id() <<"   SLAVE OBJECT = " << (*it_pair)[slave]->Id() << std::endl;
		  switch(Exist)
		    {
		    case(yes_nodes):
		     {
		       //identifico caso
		       Case = IntersectTriangles.LocateCaseItersection(Id_Node_Case_5, Change, InsideNodes, (*it_pair)[master], (*it_pair)[slave]);
		       switch(Case) 
		       {
			 case 1 : 
			 {
			    unsigned int& id  = InsideNodes[0];        
			    Ids[master]       = 0;
			    Ids[slave]        = 0;
			    Near              = CheckNearNodes(master, slave, mr_model_part.Nodes()(id), (*it_pair)[master], Ids);
			    CreateLinkingConditions(master, slave, Ids, Exist, it_pair, tempProperties,  initialize, Id, TotalInsideNodes, LinkingConditions );  
			    break;
			 }
			 
			 case 2 :
			 {
			  for(unsigned int in = 0; in<InsideNodes.size(); in++){
			     unsigned int& id  = InsideNodes[in];        
			     Ids[master]       = 0;
			     Ids[slave]        = 0;
		             Near              = CheckNearNodes(master, slave, mr_model_part.Nodes()(id), (*it_pair)[master], Ids);
			     CreateLinkingConditions(master, slave, Ids, Exist, it_pair, tempProperties,  initialize, Id, TotalInsideNodes, LinkingConditions );  
			     }
			    break;
			 }
			 
			 case 3:
			 {
			   unsigned int& id  = InsideNodes[0];  
			   Ids[master]       = 0;
			   Ids[slave]        = 0;
		           Near              = CheckNearNodes(master, slave, mr_model_part.Nodes()(id), (*it_pair)[master], Ids);
			   CreatePointLinkingConditions(master, slave,  id, Ids, it_pair, tempProperties, initialize, Id, TotalInsideNodes, LinkingConditions);
			   break;
			 }  
			 
			 case 5:
			 {
			   unsigned int& id  = InsideNodes[0];  
			   Ids[master]       = 0;
			   Ids[slave]        = 0;
		           Near              = CheckNearNodes(master, slave, mr_model_part.Nodes()(id), (*it_pair)[master], Ids); 
			   CreateLinkingConditions(master, slave, Ids, Exist, it_pair, tempProperties,  initialize, Id, TotalInsideNodes, LinkingConditions );
			   
			   ///WARNING = SPECIAL CASE
			   Exist              =  no_nodes;
			   unsigned int& id_2 =  Id_Node_Case_5;
			   Near               = CheckNearNodes(slave, master,  mr_model_part.Nodes()(id_2), (*it_pair)[slave], Ids);
			   CreateLinkingConditions(slave, master, Ids, Exist, it_pair, tempProperties,  initialize, Id, TotalInsideNodes, LinkingConditions );            
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
		  
		  //std::cout << "************************************"<< std::endl;
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
	 const PropertiesType::Pointer& tempProperties, 
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
	    std::vector<double>           Distance_aux;
	    std::vector<double>::iterator it;
	    std::vector<double>::iterator it_2;
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

            /*
            // verificando si los nodos son muy cercanos  
            double ratio       = 0.00;
            const double max   = (*std::max_element(Distance.begin(), Distance.end() ) );  
	    for(unsigned int i = 0; i<Distance.size(); i++){
	       ratio =  std::fabs(Distance[i]/max);
	       if(ratio < 1E-8){
		  Id[master] = geom_0[M[i]].Id();   
		  Id[slave]  = geom_1[S[i]].Id();
		  Ids.push_back(Id); 
	          Is_Near.push_back(yes_near);
	       }
	    }
	    */
	    
	    //Check si dos corner chocan
	    std::vector<unsigned int > nodes;
	    if(VerifyToCornerIntersect(nodes, SlaveObject, MasterObject)==false) 
	    { 
		if(nodes.size()==2)
		{
		Id[master]         = nodes[0]; //geom_0[M[position]].Id();   
		Id[slave]          = nodes[1]; //geom_1[S[position]].Id();
		}
		else
		{
		 // si no se cumple lo anterior tomamos los nodos mas cercanos
		 // WARNING = Solo valido para un caso en que un solo nodo quede fuera
		 const double min   = (*std::min_element(Distance.begin(), Distance.end() ) );
		 it                 = std::find(Distance.begin(), Distance.end(), min);
		 const int position = int(it-Distance.begin());
		 Id[master]         = geom_0[M[position]].Id();   
		 Id[slave]          = geom_1[S[position]].Id();
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
		Id[master]         = geom_0[M[position]].Id();   
		Id[slave]          = geom_1[S[position]].Id();
		Ids.push_back(Id); 
		Is_Near.push_back(no_near);
		
		
		
		const double min_2   = (*min_element_2(Distance_aux.begin(), Distance_aux.end(), min ) );
		it_2                 =  std::find(Distance.begin(), Distance.end(), min_2);
		const int position_2 = int(it_2-Distance.begin()); 
		Id[master]           = geom_0[M[position_2]].Id();   
		Id[slave]            = geom_1[S[position_2]].Id();
		Ids.push_back(Id); 
		Is_Near.push_back(no_near);
	    }
	    
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

    bool VerifyToCornerIntersect(  std::vector<unsigned int >& Ids,
                                   const PointerType&    SlaveObject,
	                           const PointerType&   MasterObject
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

        
    
    if(neighb_cond_master.size()==2 && neighb_cond_slave.size()==2)
    {

        Condition::GeometryType& geom_1 = (neighb_cond_master(0).lock())->GetGeometry();
	Condition::GeometryType& geom_2 = (neighb_cond_master(1).lock())->GetGeometry();
	
	if(geom_1[0].Id()==geom_2[0].Id())
	  Ids.push_back(geom_1[0].Id());
	else if(geom_1[0].Id()==geom_2[1].Id())
	  Ids.push_back(geom_1[0].Id());
	else if(geom_1[1].Id()==geom_2[0].Id())
	  Ids.push_back(geom_1[1].Id());
	else 
	  Ids.push_back(geom_1[1].Id());
	   
	
	Condition::GeometryType& geom_3 = (neighb_cond_slave(0).lock())->GetGeometry();
	Condition::GeometryType& geom_4 = (neighb_cond_slave(1).lock())->GetGeometry();
	
	if(geom_3[0].Id()==geom_4[0].Id())
	  Ids.push_back(geom_3[0].Id());
	else if(geom_3[0].Id()==geom_4[1].Id())
	  Ids.push_back(geom_3[0].Id());
        else if(geom_3[1].Id()==geom_4[0].Id())
	  Ids.push_back(geom_3[1].Id());
	else
	  Ids.push_back(geom_3[1].Id());
	
	return false;
      }
      
    
    if(segment.size()==3) 
	if(segment[0].size()== 2 && segment[1].size()== 2 && segment[2].size()== 2)
	  return true;
    
    return false;

   KRATOS_CATCH("")
       
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
	    if(ratio < 1E-8) 
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
	 const PropertiesType::Pointer& tempProperties, 
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
		//std::cout<< "TEST A " << std::endl; 
	       //el elemento tiene una sola condicion master
	        bool test_1 = bool(neighb_cond.size()==1);
		bool test_2 = Test_One_B(segmento, SlaveNode, MasterNode);
		
	        if(test_1==true )
		 {
		  if(test_2==true)
	           {  
		      What_Is = is_object;
	              segmento = 0;
	              return true;
	           } 
		 }
		 
	       //std::cout<< "TEST B " << std::endl; 
	       //if(test_1==false)
	       {
	         if(Test_One_B_Distances(segmento, SlaveNode, MasterObject)) 
	        {
		   What_Is = is_object;
		   return true; 
	        }
		 
	       //std::cout<< "TEST C " << std::endl; 
	       if(Test_One_B(segmento, SlaveNode, MasterNode)) 
	          { 
		    What_Is = is_node;
		    return true;
		  }
	       //std::cout<< "TEST D " << std::endl;	  
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
	       }
		  break;
	       }
	       
	      case(no_nodes):
	       {
		 
		//con desplazamientos 
		//std::cout<< "TEST A " << std::endl; 
		if(Test_One_A(segmento, SlaveNode, MasterNode)) 
	          { 
		    What_Is = is_node;
		    return true;
		  }
		  
		 // con aristas
		//std::cout<< "TEST B " << std::endl; 
		if(Test_One_C(segmento, SlaveNode, MasterNode)) 
	          { 
		    What_Is = is_node;
		    return true;
		  }
		  
	       //std::cout<< "TEST C " << std::endl; 	  
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

      I  = 0; 
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
     return true;
  }
  
  return false;
  KRATOS_CATCH("")
}



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

  //std::cout<< "FINISHHHHHHHH" << std::endl;
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



//Calculando distancias de puntos a segmentos 
bool Test_One_B_Distances(unsigned int& segmento,
	       const NodePointerType& SlaveNode,
	       const PointerType& MasterObject) 
{
  
  KRATOS_TRY   
  
  WeakPointerVector<Condition>& neighb_cond_master     = MasterObject->GetValue(NEIGHBOUR_CONDITIONS);

  std::vector<unsigned int> segment;
  unsigned int I       = 0;
  segmento             = 0;

     
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
		                        geom_1[1].GetValue(IS_BOUNDARY) = 1.00; //FastGetSolutionStepValue(IS_BOUNDARY) = 1.00;
					geom_1[2].GetValue(IS_BOUNDARY) = 1.00;
			              }
			              
				      if (face==1)  // edge 2-0
				      {
					 Pair[0] =   geom_1[2].Id();
					 Pair[1] =   geom_1[0].Id();
					 CreateMasterConditions(Pair, elem, Id, MasterConditions);
					 geom_1[2].GetValue(IS_BOUNDARY) = 1.00; //FastGetSolutionStepValue(IS_BOUNDARY) = 1.00;
					 geom_1[0].GetValue(IS_BOUNDARY) = 1.00; //FastGetSolutionStepValue(IS_BOUNDARY) = 1.00;
					 
				      }
				      if (face==2) // edge 0-1
				      {
					 Pair[0] =   geom_1[0].Id();
					 Pair[1] =   geom_1[1].Id();
					 CreateMasterConditions(Pair, elem, Id, MasterConditions);
					 geom_1[0].GetValue(IS_BOUNDARY) = 1.00; //FastGetSolutionStepValue(IS_BOUNDARY) = 1.00;
					 geom_1[1].GetValue(IS_BOUNDARY) = 1.00; // FastGetSolutionStepValue(IS_BOUNDARY) = 1.00;
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
       
       ContainerType                 mBoundaryElements;
       ContainerContactPair          mPairContacts;
       ConditionsArrayType           mMasterConditionsArray;
       //BinsObjectDynamic<Configure>  mBinsObjectDynamic;
       
       
 
    };
	
	


}  // namespace Kratos.

#endif // KRATOS_GEOMETRY_UTILITIES_INCLUDED  defined 


