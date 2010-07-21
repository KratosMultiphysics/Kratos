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
#define  KRATOS_GEOMETRY_BOX_UTILITIES_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <fstream> 
#include <algorithm>

// External includes 


// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "geometries/geometry.h"
#include "includes/mesh.h"
#include "spatial_containers/bounding_box.h"
#include "geometries/hexahedra_3d_8.h"



namespace Kratos
{
        
        class BoundingBoxUtilities
	{
	public:
 
           typedef Node<3> NodeType;  
           typedef Geometry<NodeType> GeometryType;        
           typedef BoundingBox<NodeType,  GeometryType>  BoundingBoxType;    

            BoundingBoxUtilities(){}
            BoundingBoxUtilities(ModelPart& model_part, const unsigned int& dimension) : mr_model_part(model_part), mrdimension(dimension) 
              {  
              }   

            virtual ~BoundingBoxUtilities(){}

            void Test() 
                {  
                     AddBoundingBoxToMesh(mr_model_part);
                }  

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
                        name+=".msh";
                        std::ofstream output_file( name.c_str());
                         
                        if (mrdimension==3)
                        {   
                       
                        output_file << "MESH \"Boxes\" dimension 3 ElemType Hexahedra Nnode 8" << std::endl;
			output_file << "Coordinates" << std::endl;
			for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
			 {						
                              
                              
			      it->GetGeometry().Bounding_Box(Low, High); 
                              BoundingBoxType rThisBoundingBox(Low,High,&it->GetGeometry());
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
                              it->GetGeometry().Bounding_Box(Low, High); 
                              BoundingBoxType rThisBoundingBox(Low,High,&it->GetGeometry());
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

       private:
       ModelPart mr_model_part; 
       unsigned int mrdimension;  

	};

}  // namespace Kratos.

#endif // KRATOS_GEOMETRY_UTILITIES_INCLUDED  defined 


