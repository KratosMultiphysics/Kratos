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
#include "custom_conditions/master_contact_face_2d.h"
#include "custom_conditions/point_segment_contact_link.h"


namespace Kratos
{
  
///******************************************************************************************************************
///******************************************************************************************************************
 
 
 class BoundaryConditionsAndContactUtilities
	{
	public:
	  
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
		  
		  
		  /// Crea las conciones de contacto
		  //if (minitialize==false) 
		  CalculateBoundaryContour(mMasterConditionsArray);
  
		  if(SearchContactsPairs())  
		     CreateLinkingConditions(); 
		  
		  Clear();
		
	      
	          KRATOS_CATCH("")
          }



       
	 private:

	   
         //************************************************************************************
	 //************************************************************************************   
         bool SearchContactsPairs()
         {
	  
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
	 }
	    
	  	   
         //************************************************************************************
	 //************************************************************************************ 
	 
	 void Clear()
	 {
	   mPairContacts.clear();
	 }
	  
	 //************************************************************************************
	 //************************************************************************************  
	    
	 void CreateLinkingConditions() 
	 {
	   
	    ConditionsArrayType& rConditions = mr_model_part.Conditions();
	    ElementsArrayType& rElements   = mr_model_part.Elements();
	    ConditionsArrayType LinkingConditions;
	     
	    std::vector<unsigned int>              InsideNodes;
	    std::vector<unsigned int>              TotalInsideNodes;
	    PointsArrayType                        Slaves;
	    std::vector<unsigned int>::iterator    repeated_object;  
	       
	    PropertiesType::Pointer pProperties = mr_model_part.pGetProperties(1);
	    
	    unsigned int Id = rElements.size() + rConditions.size() + 1;
    
	    /// Aux variables for contact links
	    int properties_index                   = mr_model_part.NumberOfProperties();
	    PropertiesType::Pointer tempProperties = PropertiesType::Pointer(new PropertiesType(properties_index+1) );
	    mr_model_part.AddProperties( tempProperties );

	    int  master     = 0;
	    bool initialize = false;
	    for(IteratorContainerContactPair it_pair = mPairContacts.begin(); it_pair!= mPairContacts.end(); it_pair++)
	        {
		   /// Vericando los nodos que caen en el elemento.
		   master = 0;
		   NodeInside( (*it_pair)[0], (*it_pair)[1], InsideNodes);
                   if(InsideNodes.size()==0) 
		   { master = 1;
		     InsideNodes.clear();
		     NodeInside( (*it_pair)[1], (*it_pair)[0], InsideNodes);
		   }
		   
		   //std::cout<< "     MASTER  = "<< (*it_pair)[master]->Id() << "     MASTER  = " << master <<   std::endl; 
		   //std::cout<< "     PAIRS 1 = "<< (*it_pair)[0]->Id() <<  "  " <<  "PAIRS 2 = " <<  (*it_pair)[1]->Id() <<  std::endl; 
		   /// WARNING = Puede que en un tiempo determinado caigan dos nodos a la vez en un  elemento
		   /// Un nodo dentro del elemento
		   here:
		   for(unsigned int in = 0; in<InsideNodes.size(); in++)
		     {
		       
		       if(TotalInsideNodes.size()==0){   
		          TotalInsideNodes.push_back(InsideNodes[in] ); 
			  repeated_object = TotalInsideNodes.end();
		       }
		          
		          
			  /// verifcando que no se repitan los slaves nodes 
			  if(initialize==true)
			      repeated_object =  std::find(TotalInsideNodes.begin(), TotalInsideNodes.end(), InsideNodes[in]); 
			   
                                  
			  ///WARNING = No necesariamente para otros contenedores la comparacion se hace con end -1 
 			  if( repeated_object == (TotalInsideNodes.end())) 
 			  {          
			       
			      if(initialize==true)
                                  TotalInsideNodes.push_back(InsideNodes[in]);  
			     
 			      //std::cout<< "     Node Inside =  "<< InsideNodes[in] << std::endl;
 		              //std::cout<< "     MASTER OBJECT =  " <<  (*it_pair)[0]->Id() <<"   SLAVE OBJECT = " << (*it_pair)[1]->Id() << std::endl;
			     
			       /// Slave Node
			      Point2D<Node<3> >::Pointer point_geom    =  Point2D<Node<3> >::Pointer( new Point2D<Node<3> >(mr_model_part.Nodes()(InsideNodes[in]) ) );
			      Condition::Pointer SlaveNode             =  Condition::Pointer(new SlaveContactPointType(Id, point_geom) ); 
			     
			      WeakPointerVector<Condition> neighb_cond =  (*it_pair)[master]->GetValue(NEIGHBOUR_CONDITIONS); 			      
			      Condition::Pointer MasterFace            =  (neighb_cond(0).lock());   /// me devuelve el puntero
			     
			      /// creando geometria trinagular para el link
			      Condition::GeometryType& Mgeom = MasterFace->GetGeometry();
			      Condition::GeometryType& Sgeom = SlaveNode->GetGeometry();   
			      Triangle2D3<Node<3> >::Pointer Lgeom    =  Triangle2D3<Node<3> >::Pointer( new Triangle2D3<Node<3> >( Sgeom(0), Mgeom(0), Mgeom(1) ) );
			      
			      Condition::Pointer newLink    = Condition::Pointer( new PointSegmentContactLink(Id,
				  Lgeom,
				  tempProperties,
				  MasterFace, 
				  SlaveNode) );
				  
                              LinkingConditions.push_back( newLink );    
			      Id++;
			      
			  }
			  
			  if(initialize==false) {initialize=true;}
                              InsideNodes.clear();
		        }
		        
		        if(master==0){ 
                         InsideNodes.clear();
		            NodeInside( (*it_pair)[1], (*it_pair)[0], InsideNodes);
			    if(InsideNodes.size()!=0) {
			       master = 1;
			       //std::cout<< "     Second Searching "<<std::endl;  
			       goto here;
			    }
			}
		            
		        
		        
		       }
		       
		     
		      
		     std::cout<<"     NUMBER OF CALCULATED CONDITIONS = " << LinkingConditions.size() <<  std::endl;  
		     
		    ///adding linking to model_part
                    for(ConditionsArrayType::ptr_iterator it=LinkingConditions.ptr_begin();
                          it != LinkingConditions.ptr_end(); ++it )
                    {
                        mr_model_part.Conditions().push_back( *it );
                    }
                    LinkingConditions.clear();
		    
		           
          }
	 
 	
 	 void LocateMasterSegement(PointerType& MasterObject)
 	 {   
	 }
	
	
	  //************************************************************************************
	  //************************************************************************************ 
	 void CalculateBoundaryContour(ConditionsArrayType& MasterConditions)
	 {   
	      typedef WeakPointerVector< Element >::iterator  ElementIteratorType;  	   
	      ContainerType& rElements           =  mr_model_part.ElementsArray();
	      IteratorType it_begin              =  rElements.begin();
	      IteratorType it_end                =  rElements.end(); 
	      
	      array_1d<unsigned int,2>  Pair;
   
	      vector<unsigned int> node_boundary;
	      bool is_boundary   = false;
	      
	      unsigned int face     = 0; 
	      unsigned int Id       = rElements.size() + 1 ;
	      for(IteratorType elem = it_begin; elem!=it_end; elem++)
	      {
		   Element::GeometryType& geom_1 = (*elem)->GetGeometry();
		   WeakPointerVector< Element >& neighb_elems  = (*elem)->GetValue(NEIGHBOUR_ELEMENTS); 
		   node_boundary.resize(neighb_elems.size(), false);
		   /// Puede incluir como vecnino el mismo en caso de que hayan menos de 3 elemtos veninos.  
		   /// ckeck si se repited elmento
		   /// ElementIteratorType no necesita especificarse el * 
		   for( ElementIteratorType neighb_elem  = neighb_elems.begin(); neighb_elem!= neighb_elems.end(); neighb_elem++)
 		            {
			       
			       if ( neighb_elem->Id() ==  (*elem)->Id() )
			           {     
				     if(face == 0) /// edge 1-2
				      {
					Pair[0]   =  geom_1[1].Id();
					Pair[1]   =  geom_1[2].Id();
					CreateMasterConditions(Pair, elem, Id, MasterConditions);
			              }
			              
				      if (face==1)  /// edge 2-0
				      {
					 Pair[0] =   geom_1[2].Id();
					 Pair[1] =   geom_1[0].Id();
					 CreateMasterConditions(Pair, elem, Id, MasterConditions);
					 
				      }
				      if (face==2) /// edge 0-1
				      {
					 Pair[0] =   geom_1[0].Id();
					 Pair[1] =   geom_1[1].Id();
					 CreateMasterConditions(Pair, elem, Id, MasterConditions);
				      }
				      
				       if (is_boundary==false)   
				       { 
					 is_boundary = true;
					 mBoundaryElements.push_back(*elem);
				       }
				   }
				   
			       face++;    
			    }
			  
			is_boundary = false;  
		        face = 0;  
	           }
	           
	           minitialize = true;
	 }
	 
	 
	 void CreateMasterConditions(const array_1d<double, 3>& rPair, const IteratorType& elem, unsigned int& Id, ConditionsArrayType& MasterConditions)
	 {
	    Line2D2<Node<3> >::Pointer pgeom =  Line2D2<Node<3> >::Pointer (new Line2D2<Node<3> >(mr_model_part.Nodes()((rPair)[0]), mr_model_part.Nodes()((rPair)[1]) ) ) ;  
	    Condition::Pointer MasterSegment = Condition::Pointer(new MasterContactFaceType(Id, pgeom ) ) ;
	    MasterSegment->GetValue(NEIGHBOUR_ELEMENTS).push_back(*elem);
	    (*elem)->GetValue(NEIGHBOUR_CONDITIONS).push_back(MasterSegment); 
	    MasterConditions.push_back(MasterSegment);
	    Id++;
	 }



      void FiltratePairContacts(ContainerContactPair& PairContacts)
	{
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
	}
	
      /// Se busca el nodo comun entre los contactos
      bool SearchCommonNode(const PointerType& elem1, const PointerType& elem2, std::vector<unsigned int>& id)
      {
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
      }
	
      /// Verifica si los nodos que no es el comun cae dentro del elemento
      bool SearchInsideNode(const PointerType& elem1, const PointerType& elem2, const unsigned int& ide) 
      {  
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
      }



         void NodeInside(PointerType& MasterObject, PointerType& SlaveObject,  std::vector<unsigned int>& InsideNodes)
         {
	   
	   Element::GeometryType& geom_master = MasterObject->GetGeometry();
	   Element::GeometryType& geom_slave  = SlaveObject->GetGeometry();
	   array_1d<double, 3> result;
	   for (unsigned int i = 0; i<geom_master.size(); i++ )
	       if(geom_master.IsInside(geom_slave[i], result)) {InsideNodes.push_back(geom_slave[i].Id());}      
	  
	 }


      
       private:
       ModelPart mr_model_part; 
       unsigned int mrdimension;
       bool minitialize;
       
       ContainerType             mBoundaryElements;
       ContainerContactPair      mPairContacts;
       std::vector<unsigned int> mpair;
       ConditionsArrayType       mMasterConditionsArray;
       
       

       
       
	};
	
	


}  // namespace Kratos.

#endif // KRATOS_GEOMETRY_UTILITIES_INCLUDED  defined 


