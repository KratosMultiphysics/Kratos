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


#if !defined(KRATOS_BOUNDING_BOX_UTILITIES_INCLUDED )
#define  KRATOS_BOUNDING_BOX_UTILITIES_INCLUDED



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

#include "spatial_containers/spatial_containers.h"
#include "spatial_containers/bounding_box.h"
#include "spatial_containers/cell.h"
#include "spatial_containers/bins_dynamic_objects.h"

#include "utilities/spatial_containers_configure.h"
#include "utilities/geometry_utilities.h"
#include "utilities/timer.h"
#include "utilities/timer_CLabra.h"

namespace Kratos
{
  
///******************************************************************************************************************
///******************************************************************************************************************

template <std::size_t TDimension>
class BoxFunction
  {
   public:   
   template< class TPointType, class TPointerType> 
   void operator ()(TPointerType& rObject, TPointType& rLowPoint, TPointType& rHighPoint)
    { 
    rHighPoint = rObject.GetGeometry().GetPoint(0);
    rLowPoint  = rObject.GetGeometry().GetPoint(0);        
     for (unsigned int point = 0; point<rObject.GetGeometry().PointsNumber(); point++)      
       {
 	for(std::size_t i = 0; i<TDimension; i++)
	  {
 	      rLowPoint[i]  =  (rLowPoint[i]  >  rObject.GetGeometry().GetPoint(point)[i] ) ?  rObject.GetGeometry().GetPoint(point)[i] : rLowPoint[i]; 
 	      rHighPoint[i] =  (rHighPoint[i] <  rObject.GetGeometry().GetPoint(point)[i] ) ?  rObject.GetGeometry().GetPoint(point)[i] : rHighPoint[i];
 	  }
        }           
    }
  };



///******************************************************************************************************************
///******************************************************************************************************************
		
class TriBoxOverlapFunction
  {
    public:		
    template< class TPointType, class TPointerType> 
    bool operator ()(TPointerType& rObject,  const TPointType& rLowPoint, const TPointType& rHighPoint)
    { 
      return rObject.GetGeometry().HasIntersection(rLowPoint, rHighPoint); 
    }
 };
 

///******************************************************************************************************************
///******************************************************************************************************************

class TDistanceFunction
  {
    public:		
    template<class TPointerType> 
     bool operator ()(TPointerType& rObj_1, TPointerType& rObj_2)
      { 
	      Element::GeometryType& geom_1 = rObj_1.GetGeometry();
	      Element::GeometryType& geom_2 = rObj_2.GetGeometry();
	      return  geom_1.HasIntersection(geom_2); 
      }
   };

template<class TPointType, std::size_t TDimension>   
class Segment : public Point<TDimension, double> 
{
  enum {Dimension = TDimension };
  typedef Point<Dimension, double>    PointType; 
  typedef array_1d<double, Dimension> VectorType;
  public:
  
  Segment(const PointType& rPoint1, const PointType& rPoint2) :
             mPoint1(rPoint1), mPoint2(rPoint2)   
    {
    }
   
  PointType Center()
  {
    return 0.50 * (mPoint1 + mPoint2);
  }
   
  double Length()
  {
    return norm_2(Direction());
  }
  
  VectorType Direction()
  {
     return mPoint2 - mPoint1;
  }
  
  double Extent()
  {
    return    0.50 * Length();
  }
  
 double Normalize(const double epsilon = 1E-9)
   {
     const double length = Length();
     VectorType result   = Direction(); 
    if (length > epsilon)
    { 
        const double invLength = 1.00 / length;
        for(std::size_t i = 0; i<Dimension; i++)
           result[i] =  result * invLength;
    }
    else
    {
        for(std::size_t i = 0; i<Dimension; i++)
           result[i] =  0.00; 
    }
    
   }
   
   
  PointType  mPoint1;
  PointType  mPoint2;
  
  
    
};
   
   
///******************************************************************************************************************
///******************************************************************************************************************
 
 
 class BoundingBoxUtilities
	{
	public:
  
// 	    typedef Node<3> NodeType;  
// 	    typedef ModelPart::ElementsContainerType::ContainerType  ContainerType;
// 	    typedef ContainerType::value_type                        PointerType;
// 	    typedef ContainerType::iterator                          IteratorType; 
// 	    typedef BoundingBox<NodeType,  PointerType>              BoundingBoxType;
// 	    

            const static std::size_t dimension = 2;
	    typedef Point<dimension, double>                         PointType;
	    typedef Segment<PointType, 2>                            SegmentType;
	    typedef SegmentType*                                     SegmentPointer;
	    typedef std::vector<SegmentType>                         ContainerSegmentType;  
	    typedef ModelPart::ElementsContainerType::ContainerType  ContainerType; 
	    typedef ContainerType::value_type                        PointerType;
	    typedef ContainerType::iterator                          IteratorType;  
	    typedef SpatialContainersConfigure<dimension>            Configure2D;   
	    typedef Cell<Configure2D>                                CellType;
	    typedef std::vector<CellType>                            CellContainerType;
	    typedef CellContainerType::iterator                      CellContainerIterator;
	    typedef std::vector<PointerType>::iterator               PointerTypeIterator;
	    typedef ContactPair<PointerType>                         ContactPairType; 
	    typedef std::vector<ContactPairType>                     ContainerContactPair;   
	    typedef ContainerContactPair::iterator                   IteratorContainerContactPair;  
	    typedef ContainerContactPair::value_type                 PointerContainerContactPair;  
	    typedef Element::GeometryType                            GeomType; 
	    typedef Node<3>                                          NodeType;  
	    
	   
            BoundingBoxUtilities(){}
            BoundingBoxUtilities(ModelPart& model_part, const unsigned int& dimension) : mr_model_part(model_part), mrdimension(dimension) 
              {  
              }   

            virtual ~BoundingBoxUtilities(){}

            void Test() 
                { 
		  
		  timer Time;
		  Time.start("Starting Time");
		  
		  ContainerType& rElements  =  mr_model_part.ElementsArray();
		  IteratorType it_begin     =  rElements.begin();
		  IteratorType it_end       =  rElements.end(); 
		  
		  PointType MaxPoint, MinPoint;  
		  BinsObjectDynamic<Configure2D>  rBinsObjectDynamic(it_begin, it_end );
		  Time.stop("Stopping Time");
		  std::cout<< "Time = " << Time << std::endl; 
		  
		  /*std::size_t MaxNumberOfResults = 1E7; 
		  std::size_t NumberOfResults    = 0;  
		  ContainerType Results_1(MaxNumberOfResults); 
		  IteratorType ResultsIterator = Results_1.begin();  
		  
		  for(IteratorType elem = it_begin; elem!=it_end; elem++)      
                         NumberOfResults+=rBinsObjectDynamic.SearchObjects(*elem, ResultsIterator, MaxNumberOfResults);
		  

		  std::cout<< "NumberOfResults = "<< NumberOfResults << std::endl;
		  *///std::cout<< "Time = " << Time << std::endl; 
		  
 		  //std::size_t MaxNumberOfResults = 20; 
 		  //std::size_t NumberOfResults    = 0;     
		  //ElementsArrayType Results_1(MaxNumberOfResults); 
		  //ElementsArrayType Results_2;
		  //ElementsArrayTypeIterator ResultsIterator = Results_1.begin();  

		  
 		  //ContainerContactPair PairContacts;
		  //rBinsObjectDynamic.SearchContact(PairContacts);
		  //std::cout<< PairContacts.size()<<std::endl;
		  
		  
//   		  for (unsigned int i = 0; i<PairContacts.size(); i++)
//                         std::cout<< (PairContacts[i][0])->Id() << "   " << (PairContacts[i][1])->Id()<<std::endl;
		  
		  /*
		  std::vector<NodeType*> InsideNodes; InsideNodes.resize(0, false);
		  for(IteratorContainerContactPair it_pair = PairContacts.begin(); it_pair!=PairContacts.end(); it_pair++ )
		  {     
		     /// Computa el nodo slave, el primer par es el master el segundo es el slave
		     //NodeInside((*it_pair)[0], (*it_pair)[1], InsideNodes);
		     NodeInside((*it_pair)[1], (*it_pair)[0], InsideNodes);
		     
		     if(InsideNodes.size()!=0)
		     {
			    /// distancia mas cercana del punto al segmento del objeto
			  Element::GeometryType& geom_master = (*it_pair)[1]->GetGeometry();
			  Element::GeometryType& geom_slave  = (*it_pair)[0]->GetGeometry();
			  PointType sPoint, m1Point,  m2Point, m3Point;

			  sPoint[0] = (InsideNodes[0])->X();
			  sPoint[1] = (InsideNodes[0])->Y();
			  
			  
			  m1Point[0] = geom_master[0].X();
			  m1Point[1] = geom_master[0].Y();
			  
			  m2Point[0] = geom_master[1].X();
			  m2Point[1] = geom_master[1].Y();
			  
			  m3Point[0] = geom_master[2].X();
			  m3Point[1] = geom_master[2].Y();
			  
      		          //ContainerSegmentType MasterSegments(3); 
			  //const std::size_t size = geom_master.size();
			  //MasterSegments.resize(size); 
       		          SegmentType a(m1Point, m2Point);
       		          SegmentType b(m2Point, m3Point);
       		          SegmentType c(m3Point, m1Point);
			  
			  SegmentPointer aa = &a;
       		          DistPointSegment(sPoint, aa);
       		     
			  
			  //NodeInside((*it_pair)[1], (*it_pair)[0], InsideNodes);

				  
			  InsideNodes.resize(0, false);
		     }
		  }
		  
		  //ContainerContactPair PairContacts(MaxNumberOfResults);
 		  //IteratorContainerContactPair Pair = PairContacts.begin();
 		  //NumberOfResults = mybins.SearchContact(Pair, MaxNumberOfResults);
 
                        
		  
		  
		  
                     //AddBoundingBoxToMesh(mr_model_part);
	             //AddBigBoundingBoxToMesh(mr_model_part);
		     //std::cout<< "Begin AddCellsTomesh " << std::endl; 
                     //AddCellsTomesh(mr_model_part);
		     //std::cout<< "End AddCellsTomesh " << std::endl;
		     */
		}
		
	
	/*
	void GlobalSearchScheme(const PointType& rPoint, const SegmentType rSegment)
	{
	  PointType& P0  = rSegment[0];
	  PointType& P1  = rSegment[1];
	   
	}
		
        */
        /// In segmento es un vector de dos puntos
	/*
         double DistPointSegment(const PointType& rPoint, const SegmentPointer& rSegment)
         {
 	      const PointType diff =   rPoint - rSegment->Center();
              const double param   =   inner_prod(rSegment->Direction(), diff);
	      
	      if (-rSegment->Extent < param)
		{
		if (param < mSegment->Extent)
		{
		mClosestPoint1 = mSegment->Center + param*mSegment->Direction;
		}
		else
		{
		mClosestPoint1 = mSegment->P1;
		}
	      }
	      else
	      {
	        mClosestPoint1 = mSegment->P0;
	      }

	      mClosestPoint0 = *mPoint;
	      diff = mClosestPoint1 - mClosestPoint0;
	      return diff.SquaredLength();
	      
	      return 10.0;
         }
          
          
         void NodeInside(PointerType& MasterObject, PointerType& SlaveObject, std::vector<NodeType*>& InsideNodes)
         {
	   
	   Element::GeometryType& geom_master = MasterObject->GetGeometry();
	   Element::GeometryType& geom_slave  = SlaveObject->GetGeometry();
	   array_1d<double, 3> result;
	   for (unsigned int i = 0; i<geom_master.size(); i++ )
	       if(geom_master.IsInside(geom_slave[i], result)) {InsideNodes.push_back(&geom_slave[i]);}      
	  
	 }
	 	 
		
            */
            /*
            void AddBoundingBoxToMesh(
                    ModelPart& rThisModelPart)
                    {
                         

// 			typedef Node<3> NodeType;  
// 		        typedef ModelPart::ElementsContainerType::ContainerType  ContainerType;
//                         typedef ContainerType::value_type                        PointerType;
// 			typedef ContainerType::iterator                          IteratorType; 
// 			typedef BoundingBox<NodeType,  PointerType>              BoundingBoxType;

			ContainerType& rElements =    rThisModelPart.ElementsArray();
			IteratorType it_begin    =    rElements.begin();
			IteratorType it_end      =    rElements.end(); 
			
		        NodeType High, Low;

			
                        const unsigned int current_id =  (**(rElements.end()-1)).Id(); 			
                        unsigned int nofn =  0; //rThisMesh.NumberOfNodes();       
                        std::string name = rThisModelPart.Name(); 
                        name+="_boxes.msh";
                        std::ofstream output_file( name.c_str());
                         
                        if (mrdimension==3)
                        {   
                       
                        output_file << "MESH \"Boxes\" dimension 3 ElemType Hexahedra Nnode 8" << std::endl;
			output_file << "Coordinates" << std::endl;
			for (IteratorType it= it_begin; it!=it_end; ++it)
			 {	
			      
			      (**it).GetGeometry().BoundingBox(Low, High); 
                              BoundingBoxType rThisBoundingBox(Low, High); //(**it));
                              //Boxes.push_back( rThisBoundingBox );
                                

			      const double& xmin = (rThisBoundingBox.LowPoint()).X();
			      const double& ymin = (rThisBoundingBox.LowPoint()).Y();               
			      const double& zmin = (rThisBoundingBox.LowPoint()).Z();
			      const double& xmax = (rThisBoundingBox.HighPoint()).X();
			      const double& ymax = (rThisBoundingBox.HighPoint()).Y();               
			      const double& zmax = (rThisBoundingBox.HighPoint()).Z();
                                
			      
			      output_file << nofn+1 << "  " <<  xmin << "  " <<  ymin << "  " <<  zmin  << " " << std::endl;
                              output_file << nofn+2 << "  " <<  xmax << "  " <<  ymin << "  " <<  zmin  << " " << std::endl;
                              output_file << nofn+3 << "  " <<  xmax << "  " <<  ymax << "  " <<  zmin  << " " << std::endl;
                              output_file << nofn+4 << "  " <<  xmin << "  " <<  ymax << "  " <<  zmin  << " " << std::endl;
			      output_file << nofn+5 << "  " <<  xmin << "  " <<  ymin << "  " <<  zmax  << " " << std::endl;
			      output_file << nofn+6 << "  " <<  xmax << "  " <<  ymin << "  " <<  zmax  << " " << std::endl;
			      output_file << nofn+7 << "  " <<  xmax << "  " <<  ymax << "  " <<  zmax  << " " << std::endl;
			      output_file << nofn+8 << "  " <<  xmin << "  " <<  ymax << "  " <<  zmax  << " " << std::endl;
                              
                              nofn += 8;  
			      
			 }

                             output_file << "end coordinates" << std::endl;  
                             output_file << "Elements" << std::endl;
                             unsigned int nofe =  0;   
                             int material_id   = 0;   
                             for (unsigned int i = 1; i<=current_id; i++)
                             {
                               output_file << i << "  " << nofe + 1  << "  " <<  nofe + 2 << "  " <<  nofe + 3  << " " <<  nofe + 4 << "  " <<  nofe + 5 << "  " <<  nofe + 6 << "  " <<  nofe + 7 << "  " <<  nofe + 8 <<  "  " << material_id << std::endl;
                               nofe+=8; 
			       
                             }   
                                 
                              output_file << "end elements" << std::endl;    
			      
                            }

                           else
                            {
                              
                             output_file << "MESH \"Boxes\" dimension 2 ElemType Quadrilateral Nnode 4" << std::endl;
			     output_file << "Coordinates" << std::endl;
			     for (IteratorType it= it_begin; it!=it_end; ++it)
			     {
			      
                              (**it).GetGeometry().BoundingBox(Low, High); 
                              BoundingBoxType rThisBoundingBox(Low,High); // *it);
                              //Boxes.push_back(rThisBoundingBox ); 			     
 
			      const double& xmin = (rThisBoundingBox.LowPoint()).X();
			      const double& ymin = (rThisBoundingBox.LowPoint()).Y();               
			      //const double& zmin = 0.00;
			      const double& xmax = (rThisBoundingBox.HighPoint()).X();
			      const double& ymax = (rThisBoundingBox.HighPoint()).Y();               
			      //const double& zmax = 0.00;
                                
			      
			      output_file << nofn+1 << "  " <<  xmin << "  " <<  ymin <<  std::endl;
                              output_file << nofn+2 << "  " <<  xmax << "  " <<  ymin <<  std::endl;
                              output_file << nofn+3 << "  " <<  xmax << "  " <<  ymax <<  std::endl;
                              output_file << nofn+4 << "  " <<  xmin << "  " <<  ymax <<  std::endl;
			      
                              nofn += 4;  
			    }

                             output_file << "end coordinates" << std::endl;  
                             output_file << "Elements" << std::endl;
                             unsigned int nofe =  0;   
                             int material_id   = 0;   
                             for (unsigned int i = 1; i<=current_id; i++)
                             {
                               output_file << i << "  " << nofe + 1  << "  " <<  nofe + 2 << "  " <<  nofe + 3  << " " <<  nofe + 4 << "  " << material_id << std::endl;
                               nofe+=4;  
			      
                             }   
                                 
                              output_file << "end elements" << std::endl;  
                               
                            
			    }
		    }

                      
                     
                    
      void AddBigBoundingBoxToMesh(ModelPart& rThisModelPart)
      {
	
// 	  typedef Node<3> NodeType;  
// 	  typedef ModelPart::ElementsContainerType::ContainerType  ContainerType;
// 	  typedef ContainerType::value_type                        PointerType;
// 	  typedef ContainerType::iterator                          IteratorType; 
// 	  typedef BoundingBox<NodeType,  PointerType>              BoundingBoxType;

	  ContainerType& rElements =    rThisModelPart.ElementsArray();
	  IteratorType it_begin    =    rElements.begin();
	  IteratorType it_end      =    rElements.end(); 
	
          NodeType MaxPoint, MinPoint;
          std::vector<BoundingBoxType> Boxes;

          std::string name = rThisModelPart.Name(); 
          name+="_bigbox.msh";
          std::ofstream output_file( name.c_str());
	  unsigned int nofn = 0; 
	  const unsigned int current_id = (**(rElements.end()-1)).Id(); 
	  
	  for (IteratorType it= it_begin; it!=it_end; ++it)
	    {						                    
	      (**it).GetGeometry().BoundingBox(MinPoint, MaxPoint); 
	      BoundingBoxType rThisBoundingBox(MinPoint, MaxPoint ); // &it->GetGeometry());
	      Boxes.push_back( rThisBoundingBox );	      
	   }
	   
	  CalculateBoundingBox(Boxes, MinPoint, MaxPoint);
	  
	  if (mrdimension==3)
	  {   

	  output_file << "MESH \"Boxes\" dimension 3 ElemType Hexahedra Nnode 8" << std::endl;
	  output_file << "Coordinates" << std::endl;
	  const double& xmin = MinPoint.X();
	  const double& ymin = MinPoint.Y(); 
	  const double& zmin = MinPoint.Z();
	  const double& xmax = MaxPoint.X();
	  const double& ymax = MaxPoint.Y();   
	  const double& zmax = MaxPoint.Z();    

	  output_file << nofn+1 << "  " <<  xmin << "  " <<  ymin << "  " <<  zmin  << " " << std::endl;
	  output_file << nofn+2 << "  " <<  xmax << "  " <<  ymin << "  " <<  zmin  << " " << std::endl;
	  output_file << nofn+3 << "  " <<  xmax << "  " <<  ymax << "  " <<  zmin  << " " << std::endl;
	  output_file << nofn+4 << "  " <<  xmin << "  " <<  ymax << "  " <<  zmin  << " " << std::endl;
	  output_file << nofn+5 << "  " <<  xmin << "  " <<  ymin << "  " <<  zmax  << " " << std::endl;
	  output_file << nofn+6 << "  " <<  xmax << "  " <<  ymin << "  " <<  zmax  << " " << std::endl;
	  output_file << nofn+7 << "  " <<  xmax << "  " <<  ymax << "  " <<  zmax  << " " << std::endl;
	  output_file << nofn+8 << "  " <<  xmin << "  " <<  ymax << "  " <<  zmax  << " " << std::endl;

	  output_file << "end coordinates" << std::endl;  

	  output_file << "Elements" << std::endl;
	  unsigned int nofe =  0;   
	  int material_id   = 0;   output_file << 1 << "  " << nofe + 1  << "  " <<  nofe + 2 << "  " <<  nofe + 3  << " " <<  nofe + 4 << "  " <<  nofe + 5 << "  " <<  nofe + 6 << "  " <<  nofe + 7 << "  " <<  nofe + 8 <<  "  " << material_id << std::endl;output_file << 1 << "  " << nofe + 1  << "  " <<  nofe + 2 << "  " <<  nofe + 3  << " " <<  nofe + 4 << "  " << material_id << std::endl;  
	  output_file << "end elements" << std::endl; 
	  
	  
	  
	  }
	  else
	  { 
	  output_file << "MESH \"BigBox\" dimension 2 ElemType Quadrilateral Nnode 4" << std::endl;
	  output_file << "Coordinates" << std::endl;

	  const double& xmin = MinPoint.X();
	  const double& ymin = MinPoint.Y();               
	  const double& xmax = MaxPoint.X();
	  const double& ymax = MaxPoint.Y();               

	  output_file << nofn+1 << "  " <<  xmin << "  " <<  ymin <<  std::endl;
	  output_file << nofn+2 << "  " <<  xmax << "  " <<  ymin <<  std::endl;
	  output_file << nofn+3 << "  " <<  xmax << "  " <<  ymax <<  std::endl;
	  output_file << nofn+4 << "  " <<  xmin << "  " <<  ymax <<  std::endl;

	  output_file << "end coordinates" << std::endl;  

	  output_file << "Elements" << std::endl;
	  unsigned int nofe =  0;   
	  int material_id   = 0;   
	  output_file << 1 << "  " << nofe + 1  << "  " <<  nofe + 2 << "  " <<  nofe + 3  << " " <<  nofe + 4 << "  " << material_id << std::endl;  
	  output_file << "end elements" << std::endl;  
	  }
	  }
   
   
       void CalculateBoundingBox(std::vector<BoundingBoxType>& rBoxes, Node<3>& rMinPoint, Node<3>& rMaxPoint  )
           {
  
	       	  rMinPoint.X() = (rBoxes[0].LowPoint()).X();
	          rMinPoint.Y() = (rBoxes[0].LowPoint()).Y();               
	          rMinPoint.Z() = (rBoxes[0].LowPoint()).Z();
	          rMaxPoint.X() = (rBoxes[0].HighPoint()).X();
	          rMaxPoint.Y() = (rBoxes[0].HighPoint()).Y(); 
	          rMaxPoint.Z() = (rBoxes[0].HighPoint()).Z();    
	        
	    for(std::size_t k = 0 ; k <  rBoxes.size(); k++){	  

	            rMaxPoint.X()  =  (rMaxPoint.X()  < (rBoxes[k].HighPoint()).X()) ? (rBoxes[k].HighPoint()).X() : rMaxPoint.X(); 
		    rMaxPoint.Y()  =  (rMaxPoint.Y()  < (rBoxes[k].HighPoint()).Y()) ? (rBoxes[k].HighPoint()).Y() : rMaxPoint.Y(); 
		    rMaxPoint.Z()  =  (rMaxPoint.Y()  < (rBoxes[k].HighPoint()).Z()) ? (rBoxes[k].HighPoint()).Z() : rMaxPoint.Z(); 
		   
		    
		    rMinPoint.X()  =   (rMinPoint.X()  > (rBoxes[k].LowPoint()).X()) ? (rBoxes[k].LowPoint()).X() : rMinPoint.X(); 
		    rMinPoint.Y()  =   (rMinPoint.Y()  > (rBoxes[k].LowPoint()).Y()) ? (rBoxes[k].LowPoint()).Y() : rMinPoint.Y();
		    rMinPoint.Z()  =   (rMinPoint.Z()  > (rBoxes[k].LowPoint()).Z()) ? (rBoxes[k].LowPoint()).Z() : rMinPoint.Z();
		 }
	  }
	  
        
   
  
         
     void AddCellsTomesh(ModelPart& rThisModelPart)
      {
	
	  const std::size_t dimension = 2;
	  typedef Point<dimension, double>                         PointType;
          typedef ModelPart::ElementsContainerType::ContainerType  ElementsArrayType; 
	  typedef ElementsArrayType::value_type                    PointerType;
	  typedef ElementsArrayType::iterator                      ElementsArrayTypeIterator;  
	  typedef BoxFunction<dimension>                           BoxType;
	  typedef TriBoxOverlapFunction                            TriBoxType; 
	  typedef TDistanceFunction                                DistanceFunctionType;                    
	  
	  
           ElementsArrayType& rElements         =    rThisModelPart.ElementsArray();
           ElementsArrayTypeIterator it_begin   =    rElements.begin();
           ElementsArrayTypeIterator it_end     =    rElements.end(); 

	  std::string name = rThisModelPart.Name(); 
          name+="_Cells.msh";
          std::ofstream output_file( name.c_str());
	  
	   
	  typedef SpatialContainersConfigure<dimension>   Configure2D;   
	  typedef Cell<Configure2D>                        CellType;
	  typedef std::vector<CellType>                    CellContainerType;
          typedef CellContainerType::iterator              CellContainerIterator;
	  typedef std::vector<PointerType>::iterator       PointerTypeIterator;
	  
	  typedef ContactPair<PointerType>                 ContactPairType; 
	  //typedef array_1d<PointerType,2>                ContactPairType;
	  typedef std::vector<ContactPairType>             ContainerContactPair;   
	  typedef std::vector<ContactPairType>::iterator   IteratorContainerContactPair;   
	  
	  //typedef std::vector<double>::iterator          DistanceIteratorType;
	  
 	  PointType MaxPoint, MinPoint;  
//  	  std::cout<< "Making the Bins " << std::endl;
//  	  std::cout<<std::fixed<<std::setprecision(10);
          BinsObjectDynamic<Configure2D> mybins(it_begin, it_end );
// 
 	  std::size_t MaxNumberOfResults = 20; 
 	  std::size_t NumberOfResults    = 0;     
 	  //ElementsArrayType Results_1(MaxNumberOfResults); 
	  //ElementsArrayType Results_2;
 	  //ElementsArrayTypeIterator ResultsIterator = Results_1.begin();  

	  ContainerContactPair PairContacts(MaxNumberOfResults);
	  IteratorContainerContactPair Pair = PairContacts.begin();
	  
	   //mybins.SearchContact(PairContacts);
           NumberOfResults = mybins.SearchContact(Pair, MaxNumberOfResults);
	  
	  for (unsigned int i = 0; i<10; i++)
	      std::cout<< *(PairContacts[i][0]) << "   " << *(PairContacts[i][1]) <<std::endl;
          
//           for (std::size_t i = 0; i<NumberOfResults; i++)
//  	     {
// 	       if(PairContacts[i].first!=NULL) 
//  	          std::cout<< *PairContacts[i].first << "  " <<  *PairContacts[i].second << std::endl;
//  	     }
	  
	  
//                if(Results_1[i]!=NULL) 	     
//  	          std::cout<< *Results_1[i]<< std::endl;
// 	   }
// 	   
// 	   KRATOS_WATCH("******************************************************************* ")
// 
// 	   NumberOfResults = mybins.SearchObjects(*(it_begin+5), Results_2);
//            KRATOS_WATCH(NumberOfResults)  
// 	   for(unsigned int i=0; i<Results_2.size(); i++)
// 	   {
//                if(Results_2[i]!=NULL) 	     
//  	          std::cout<< *Results_2[i]<< std::endl;
// 	   }
	   
          //mybins.SearchObjects(Pair);
	  
// 	  unsigned int local = 0;
// 	  for( CellContainerType::iterator icell= mybins.GetCellContainer().begin(); icell != mybins.GetCellContainer().end(); icell++)
// 	  {
// 	    local++;
// 	    std::cout<<"cell = " << local<<std::endl;
// 	    icell->GetContainerSize();
// 	    for(PointerTypeIterator it = icell->GetContainer().begin(); it!= icell->GetContainer().end(); it++ )
// 	    {
// 
// 	      std::cout<<"elem Id = " <<  it->Id()<<std::endl;
// 	    }
// 	    std::cout<<"-------------------------------"<<std::endl;
// 	    
// 	  }

          
          /*
          if (mrdimension==3)
              {
		std::cout<< "No printing mesh yet " << std::endl;
	      }
	      
	  else
	  {
	  unsigned int sizecell =  1;	  
	  for(unsigned int i = 0; i< dimension; i++ )
	  {
	    sizecell *= mybins.GetDivisions()[i]; 
	  }
	 
	  std::vector< array_1d< array_1d<double,2 > ,2 > > Cell; 
	  Cell.resize(sizecell);
	  std::size_t& filas = mybins.GetDivisions()[0];  
	  
	  for (unsigned int y = 0; y< static_cast<unsigned int> (mybins.GetDivisions()[1]); y++)
	  {
	      for (unsigned int x = 0; x< static_cast<unsigned int>(mybins.GetDivisions()[0]); x++)
	      {
                double a =   static_cast<double>(mybins.GetCellSize()[0]);   
                double b =   static_cast<double>(mybins.GetCellSize()[1]);   
	        Cell[GetIndex(filas, x, y) ][0][0] = MinPoint[0] + x * a;
 		Cell[GetIndex(filas, x, y) ][0][1] = MinPoint[1] + y * b;
		Cell[GetIndex(filas, x, y) ][1][0] = MinPoint[0] + (x + 1.0) * a;
		Cell[GetIndex(filas, x, y) ][1][1] = MinPoint[1] + (y + 1.0) * b;
	      }
	  }
	  
	  output_file << "MESH \"Cells\" dimension 2 ElemType Quadrilateral Nnode 4" << std::endl;
	  output_file << "Coordinates" << std::endl;
	  
	  unsigned int nofn = 0;
	  for (unsigned int i = 0; i<Cell.size(); i++ )
	  {
	  const double& xmin = Cell[i][0][0];
	  const double& ymin = Cell[i][0][1];
	  const double& xmax = Cell[i][1][0];
	  const double& ymax = Cell[i][1][1];              
	  
	  output_file << nofn+1 << "  " <<  xmin << "  " <<  ymin <<  std::endl;
	  output_file << nofn+2 << "  " <<  xmax << "  " <<  ymin <<  std::endl;
	  output_file << nofn+3 << "  " <<  xmax << "  " <<  ymax <<  std::endl;
	  output_file << nofn+4 << "  " <<  xmin << "  " <<  ymax <<  std::endl;

	  nofn += 4;  
	  }
	  output_file << "end coordinates" << std::endl;  
	  output_file << "Elements" << std::endl;
	  unsigned int nofe =  0;   
	  int material_id   = 0;   
	  for (unsigned int i = 1; i<=sizecell; i++)
	  {
	  output_file << i << "  " << nofe + 1  << "  " <<  nofe + 2 << "  " <<  nofe + 3  << " " <<  nofe + 4 << "  " << material_id << std::endl;
	  nofe+=4;  
	  }   

	  output_file << "end elements" << std::endl;
	  
	  }
	  
           
// 	  int celda      = 1; 
// 	  int intersect = false; 
// 	  std::vector<std::pair<int,int> > Pair;
	  //clock_t init, final;
          //init=clock();
	  
//  	  for( CellContainerIterator icell= mybins.GetCellContainer().begin(); icell != mybins.GetCellContainer().end(); icell++)
// 	    {
// 	      std::cout<< " celda = " << celda++  << std::endl; 
// 	      for(ElementsArrayTypeIterator it_1 = icell->Begin(); it_1!= icell->End(); it_1++ )
// 	        { 
// 		  for(ElementsArrayTypeIterator it_2 = it_1 + 1 ; it_2!= icell->End(); it_2++ )
// 		  {
// 		    std::cout<< " Elem " << it_1->Id() << "  " << " with elem " << it_2->Id() << std::endl;     
// 		    Element::GeometryType& geom_1 = it_1->GetGeometry();
// 		    Element::GeometryType& geom_2 = it_2->GetGeometry();
// 
// 		        bool intersect = geom_1.HasIntersection( geom_2 ) ;
// 			if (intersect==true) {Pair.push_back(std::pair<int,int>(it_1->Id(),it_2->Id()));}
// 		    
// 		  }
// 		}
// 	    }
// 	    
// 	    for (unsigned int i = 0; i<Pair.size(); i++)
// 	     {
// 	       std::cout<< Pair[i].first << "  " <<  Pair[i].second << std::endl;
// 	     }
	      
	    
 	   //final=clock()-init;
           //std::cout << "Time Looping =" << (double)final / ((double)CLOCKS_PER_SEC) << std::endl;

      }
      
   
   unsigned int GetIndex( const std::size_t filas, const unsigned int& x, const unsigned int y)
   {
      return   filas * y + x;
   }

      */  

       private:
       ModelPart mr_model_part; 
       unsigned int mrdimension;


       
       
	};
	
	


}  // namespace Kratos.

#endif // KRATOS_GEOMETRY_UTILITIES_INCLUDED  defined 


