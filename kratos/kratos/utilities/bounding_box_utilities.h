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
#include "geometries/geometry.h"
#include "includes/mesh.h"
#include "spatial_containers/spatial_containers.h"
//#include "spatial_containers/bounding_box.h"
//#include "spatial_containers/cell.h"
//#include "spatial_containers/bins_dynamic_objects.h"

#include "utilities/geometry_utilities.h"

namespace Kratos
{
  
  
  
template <std::size_t TDimension>
class BoxFunction
  {
  public:
    
   template< class TPointType, class TPointerType> 
   void operator ()(const TPointerType& rObject, TPointType& rLowPoint, TPointType& rHighPoint)
    { 

    rHighPoint = rObject->GetGeometry().GetPoint(0);
    rLowPoint  = rObject->GetGeometry().GetPoint(0);           
      
    for (unsigned int point = 0; point<rObject->GetGeometry().PointsNumber(); point++)      
      {
	for(std::size_t i = 0; i<TDimension; i++)
	  {
	      rLowPoint[i]  =  (rLowPoint[i]  <  rObject->GetGeometry().GetPoint(point)[i] ) ?  rObject->GetGeometry().GetPoint(point)[i] : rLowPoint[i]; 
	      rHighPoint[i] =  (rHighPoint[i] >  rObject->GetGeometry().GetPoint(point)[i] ) ?  rObject->GetGeometry().GetPoint(point)[i] : rHighPoint[i];
	  }
      }       
    }
};
        
 class BoundingBoxUtilities
	{
	public:
 
	    typedef Node<3> NodeType;  
	    typedef Point<3, double> PointType;
	    typedef Geometry<NodeType> GeometryType;        
	    //typedef BoundingBox<NodeType,  GeometryType>  BoundingBoxType;
            //typedef std::vector<BoundingBoxType> BoundingBoxContainerType;
            //typedef BoundingBoxContainerType::iterator BoundingBoxIterator; 
	    
	    
	   
            BoundingBoxUtilities(){}
            BoundingBoxUtilities(ModelPart& model_part, const unsigned int& dimension) : mr_model_part(model_part), mrdimension(dimension) 
              {  
              }   

            virtual ~BoundingBoxUtilities(){}

            void Test() 
                {  
                     //AddBoundingBoxToMesh(mr_model_part);
		     //AddBigBoundingBoxToMesh(mr_model_part);
		     std::cout<< "Begin AddCellsTomesh " << std::endl; 
                     AddCellsTomesh(mr_model_part);
		     std::cout<< "End AddCellsTomesh " << std::endl;
		}  

/*
            void AddBoundingBoxToMesh(
                    ModelPart& rThisModelPart)
                    {
                         
                        NodeType High, Low;
                        std::vector<BoundingBoxType> Boxes;
                        //ModelPart::MeshType& rThisMesh;
                        typedef ModelPart::ElementsContainerType ElementsArrayType; 
			ElementsArrayType& rElements         =  rThisModelPart.Elements(); 
			ElementsArrayType::iterator it_begin =  rElements.ptr_begin();
			ElementsArrayType::iterator it_end   =  rElements.ptr_end(); 
                        const unsigned int current_id = (rElements.end()-1)->Id(); 			
                        unsigned int nofn =  0; //rThisMesh.NumberOfNodes();       
                        std::string name = rThisModelPart.Name(); 
                        name+="_boxes.msh";
                        std::ofstream output_file( name.c_str());
                         
                        if (mrdimension==3)
                        {   
                       
                        output_file << "MESH \"Boxes\" dimension 3 ElemType Hexahedra Nnode 8" << std::endl;
			output_file << "Coordinates" << std::endl;
			for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
			 {						
                              
                              
			      it->GetGeometry().Bounding_Box(High, Low); 
                              BoundingBoxType rThisBoundingBox(High,Low, &it->GetGeometry());
                              Boxes.push_back( rThisBoundingBox );
                                

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
			     for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
			     {
                              it->GetGeometry().Bounding_Box(High, Low); 
                              BoundingBoxType rThisBoundingBox(High, Low, &it->GetGeometry());
                              Boxes.push_back(rThisBoundingBox ); 			     
 
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
	
          NodeType MaxPoint, MinPoint;
          std::vector<BoundingBoxType> Boxes;
          typedef ModelPart::ElementsContainerType ElementsArrayType; 
          ElementsArrayType& rElements         =  rThisModelPart.Elements(); 
          ElementsArrayType::iterator it_begin =  rElements.ptr_begin();
          ElementsArrayType::iterator it_end   =  rElements.ptr_end(); 
          std::string name = rThisModelPart.Name(); 
          name+="_bigbox.msh";
          std::ofstream output_file( name.c_str());
	  unsigned int nofn = 0; 
	  const unsigned int current_id = (rElements.end()-1)->Id(); 
	  
	  for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
	    {						                    
	      it->GetGeometry().Bounding_Box(MaxPoint, MinPoint); 
	      BoundingBoxType rThisBoundingBox(MaxPoint, MinPoint, &it->GetGeometry());
	      Boxes.push_back( rThisBoundingBox );	      
	   }
	   
	  CalculateBoundingBox(Boxes, MinPoint, MaxPoint);
	  
	  
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
  */
      
     void AddCellsTomesh(ModelPart& rThisModelPart)
      {
	
	  const std::size_t dimension = 2;
          typedef ModelPart::ElementsContainerType    ElementsArrayType; 
	  typedef ElementsArrayType::value_type       PointerType;
	  typedef ElementsArrayType::iterator         ElementsArrayTypeIterator;  
	  typedef BoxFunction<dimension>              BoxType;
	  
	  
          ElementsArrayType& rElements         =    rThisModelPart.Elements(); 
          ElementsArrayTypeIterator it_begin   =    rElements.begin();
          ElementsArrayTypeIterator it_end     =    rElements.end(); 

	  //           std::string name = rThisModelPart.Name(); 
//           name+="_Cells.msh";
//           std::ofstream output_file( name.c_str());
	  
	  typedef Cell<PointerType>       CellType;
	  typedef std::vector<CellType>              CellContainerType;
          typedef CellContainerType::iterator        CellContainerIterator;
	  typedef std::vector<PointerType>::iterator PointerTypeIterator;
	  
	 
// 	  NodeType MaxPoint, MinPoint;
//           std::vector<BoundingBoxType> Boxes;
// 	  
// 	  for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
// 	    {						                    
// 	      it->GetGeometry().Bounding_Box(MaxPoint, MinPoint); 
// 	      BoundingBoxType rThisBoundingBox(MaxPoint, MinPoint, &it->GetGeometry());
// 	      Boxes.push_back( rThisBoundingBox );	      
// 	   }
// 	   
// 	  CalculateBoundingBox(Boxes, MinPoint, MaxPoint);
	  
	  
	  std::cout<< "Making the Bins " << std::endl;
	  std::cout<<std::fixed<<std::setprecision(10);
	  BinsObjectDynamic<dimension, PointType, ElementsArrayType,  BoxType> mybins(it_begin, it_end );
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
	  
// 	  if (mrdimension==3)
//               {
// 		std::cout<< "No printing mesh yet " << std::endl;
// 	      }
// 	      
// 	  else
// 	  {
// 	  unsigned int sizecell =  1;	  
// 	  for(unsigned int i = 0; i< dimension; i++ )
// 	  {
// 	    sizecell *= mybins.GetDivisions()[i]; 
// 	  }
// 	 
// 	  std::vector< array_1d< array_1d<double,2 > ,2 > > Cell; 
// 	  Cell.resize(sizecell);
// 	  std::size_t& filas = mybins.GetDivisions()[0];  
// 	  
// 	  for (unsigned int y = 0; y< static_cast<unsigned int> (mybins.GetDivisions()[1]); y++)
// 	  {
// 	      for (unsigned int x = 0; x< static_cast<unsigned int>(mybins.GetDivisions()[0]); x++)
// 	      {
//                 double a =   static_cast<double>(mybins.GetCellSize()[0]);   
//                 double b =   static_cast<double>(mybins.GetCellSize()[1]);   
// 	        Cell[GetIndex(filas, x, y) ][0][0] = MinPoint[0] + x * a;
//  		Cell[GetIndex(filas, x, y) ][0][1] = MinPoint[1] + y * b;
// 		Cell[GetIndex(filas, x, y) ][1][0] = MinPoint[0] + (x + 1.0) * a;
// 		Cell[GetIndex(filas, x, y) ][1][1] = MinPoint[1] + (y + 1.0) * b;
// 	      }
// 	  }
// 	  
// 	  output_file << "MESH \"Cells\" dimension 2 ElemType Quadrilateral Nnode 4" << std::endl;
// 	  output_file << "Coordinates" << std::endl;
// 	  
// 	  unsigned int nofn = 0;
// 	  for (unsigned int i = 0; i<Cell.size(); i++ )
// 	  {
// 	  const double& xmin = Cell[i][0][0];
// 	  const double& ymin = Cell[i][0][1];
// 	  const double& xmax = Cell[i][1][0];
// 	  const double& ymax = Cell[i][1][1];              
// 	  
// 	  output_file << nofn+1 << "  " <<  xmin << "  " <<  ymin <<  std::endl;
// 	  output_file << nofn+2 << "  " <<  xmax << "  " <<  ymin <<  std::endl;
// 	  output_file << nofn+3 << "  " <<  xmax << "  " <<  ymax <<  std::endl;
// 	  output_file << nofn+4 << "  " <<  xmin << "  " <<  ymax <<  std::endl;
// 
// 	  nofn += 4;  
// 	  }
// 	  output_file << "end coordinates" << std::endl;  
// 	  output_file << "Elements" << std::endl;
// 	  unsigned int nofe =  0;   
// 	  int material_id   = 0;   
// 	  for (unsigned int i = 1; i<=sizecell; i++)
// 	  {
// 	  output_file << i << "  " << nofe + 1  << "  " <<  nofe + 2 << "  " <<  nofe + 3  << " " <<  nofe + 4 << "  " << material_id << std::endl;
// 	  nofe+=4;  
// 	  }   
// 
// 	  output_file << "end elements" << std::endl;
// 	  
// 	  }
	  

	  int celda      = 1; 
	  int intersect = false; 
	  std::vector<std::pair<int,int> > results;
	  clock_t init, final;
          init=clock();
	  
	  for( CellContainerType::iterator icell= mybins.GetCellContainer().begin(); icell != mybins.GetCellContainer().end(); icell++)
	    {
	      //std::cout<< " celda = " << celda++  << std::endl; 
	      for(PointerTypeIterator it_1 = icell->GetContainer().begin(); it_1!= icell->GetContainer().end(); it_1++ )
	        { 
		  for(PointerTypeIterator it_2 = it_1 + 1 ; it_2!= icell->GetContainer().end(); it_2++ )
		  {
		    //std::cout<< " Elem " << it_1->Id() << "  " << "with elem  " << it_2->Id() << std::endl;     
		    Element::GeometryType& geom_1 = it_1->GetGeometry();
		    Element::GeometryType& geom_2 = it_2->GetGeometry();

		        bool intersect = geom_1.HasIntersection( geom_2 ) ;
			if (intersect==true) {results.push_back(std::pair<int,int>(it_1->Id(),it_2->Id()));}
		    
		  }
		}
	    }
	    
	    for (unsigned int i = 0; i<results.size(); i++)
	     {
	       std::cout<< results[i].first << "  " <<  results[i].second << std::endl;
	     }
	      
	    
	    final=clock()-init;
            std::cout << "Time Looping =" << (double)final / ((double)CLOCKS_PER_SEC) << std::endl;

      }
      
  /*    
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


