/*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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
/* *********************************************************   
*          
*   Last Modified by:    $Author: Nelson Lafontaine $ 
*   Date:                $Date: 01-27-2010$
*   Revision:            $Revision: 1.00   $
*
* ***********************************************************/

#if !defined(SPLIT_TRIANGLE)
#define SPLIT_TRIANGLE


#ifdef _OPENMP
#include <omp.h>
#endif

#include "boost/smart_ptr.hpp"
#include <boost/timer.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/lu.hpp>


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "processes/find_nodal_neighbours_process.h"
#include "processes/find_elements_neighbours_process.h"
#include "containers/data_value_container.h"
#include "includes/mesh.h"
#include "utilities/math_utils.h"


#include <cmath>
#include <algorithm>



namespace Kratos
{
     class Split_Triangle_Elements
      {
        public:
          typedef ModelPart::NodesContainerType NodesArrayType;
	  typedef ModelPart::ElementsContainerType ElementsArrayType;
	  typedef ModelPart::ConditionsContainerType ConditionsArrayType;
          typedef boost::numeric::ublas::vector<Matrix> Matrix_Order_Tensor;
          typedef boost::numeric::ublas::vector<Vector> Vector_Order_Tensor;   
          typedef boost::numeric::ublas::vector<Vector_Order_Tensor> Node_Vector_Order_Tensor; 
          typedef Node<3> PointType;
          typedef Node<3>::Pointer PointPointerType;
          typedef std::vector<PointType::Pointer>  PointVector;
          typedef PointVector::iterator PointIterator;

          
          Split_Triangle_Elements(); //*ModelPart& model_part) : mr_model_part(model_part)
          ~Split_Triangle_Elements(){}

         
          static bool Split_Triangle
                (
                   array_1d<int,3>& triangle_ids, 
                   array_1d<int,3>& edge_ids,  
                   boost::numeric::ublas::matrix<int>& new_conectivity
                   //bool& create_element       
                )
            {     
             
             unsigned int edges_to_be_refined = 0;         

             array_1d<bool, 3>   topology; 
             topology[0] = false;
             topology[1] = false;
             topology[2] = false;
 
             if ( edge_ids[0] > 0 ){ edges_to_be_refined++;  topology[0] = true;} //create_element = true;  }   
             if ( edge_ids[1] > 0 ){ edges_to_be_refined++;  topology[1] = true;} //create_element = true;  }
             if ( edge_ids[2] > 0 ){ edges_to_be_refined++;  topology[2] = true;} //create_element = true;  }

             //KRATOS_WATCH(edge_ids)             

             if(edges_to_be_refined==1)
               {
                 new_conectivity.resize(2,3);
                 //create_element = true;
                      
                 /// caso 1
                 if(topology[0]== true)
                  {
                   new_conectivity(0,0) = edge_ids[0]; new_conectivity(0,1) = triangle_ids[2]; new_conectivity(0,2) = triangle_ids[0]; 
                   new_conectivity(1,0) = edge_ids[0]; new_conectivity(1,1) = triangle_ids[1]; new_conectivity(1,2) = triangle_ids[2];
                  } 

                /// caso 2
                else if( topology[1]== true)
                  {
                   new_conectivity(0,0) = edge_ids[1]; new_conectivity(0,1) = triangle_ids[0]; new_conectivity(0,2) = triangle_ids[1]; 
                   new_conectivity(1,0) = edge_ids[1], new_conectivity(1,1) = triangle_ids[2]; new_conectivity(1,2) = triangle_ids[0];
                  }

                /// caso 3 
                else if( topology[2]== true)
                    {
                      new_conectivity(0,0) = edge_ids[2];  new_conectivity(0,1) = triangle_ids[1]; new_conectivity(0,2) = triangle_ids[2]; 
                      new_conectivity(1,0) = edge_ids[2];  new_conectivity(1,1) = triangle_ids[0]; new_conectivity(1,2) = triangle_ids[1];
                    }

                return true;                    

               }
             
             
             else if(edges_to_be_refined==2)
               {
                      new_conectivity.resize(3,3); 
                      //create_element = true;
                     /// caso 4  
                     if( topology[0]== true &&  topology[1]== true )
                       {
                           new_conectivity(0,0) =  edge_ids[0];  new_conectivity(0,1) = triangle_ids[1];  new_conectivity(0,2) = edge_ids[1]; 
                           new_conectivity(1,0) =  edge_ids[0];  new_conectivity(1,1) = edge_ids[1];  new_conectivity(1,2) = triangle_ids[0];
                           new_conectivity(2,0) =  triangle_ids[0];  new_conectivity(2,1) = edge_ids[1];  new_conectivity(2,2) = triangle_ids[2];                                                                  
                       }
 
                     /// caso 5 
                     else if(topology[1]== true &&  topology[2]== true)
                       {
                          new_conectivity(0,0) =  edge_ids[1]; new_conectivity(0,1) = edge_ids[2]; new_conectivity(0,2) = triangle_ids[1]; 
                          new_conectivity(1,0) =  edge_ids[1]; new_conectivity(1,1) = triangle_ids[2]; new_conectivity(1,2) = edge_ids[2];
                          new_conectivity(2,0) =  edge_ids[2]; new_conectivity(2,1) = triangle_ids[0]; new_conectivity(2,2) = triangle_ids[1];   
                       } 

                     /// caso 3 
                     else if(topology[0]== true &&  topology[2]== true)
                       {                                 
                          new_conectivity(0,0) =  edge_ids[0]; new_conectivity(0,1) = edge_ids[2]; new_conectivity(0,2) = triangle_ids[0]; 
                          new_conectivity(1,0) =  edge_ids[0]; new_conectivity(1,1) = triangle_ids[1]; new_conectivity(1,2) = triangle_ids[2];
                          new_conectivity(2,0) =  edge_ids[2]; new_conectivity(2,1) = edge_ids[0]; new_conectivity(2,2) = triangle_ids[2]; 
                                                
                       }
                   return true;  
               }
 
             else if(edges_to_be_refined==3)
               {
                ///* case 8
                new_conectivity.resize(4,3); 
                //create_element = true;
		new_conectivity(0,0) =  edge_ids[0]; new_conectivity(0,1) = edge_ids[1]; new_conectivity(0,2) = edge_ids[2]; 
		new_conectivity(1,0) =  triangle_ids[1]; new_conectivity(1,1) = edge_ids[1]; new_conectivity(1,2) = edge_ids[0];
		new_conectivity(2,0) =  triangle_ids[2]; new_conectivity(2,1) = edge_ids[2]; new_conectivity(2,2) = edge_ids[1]; 
		new_conectivity(3,0) =  triangle_ids[0]; new_conectivity(3,1) = edge_ids[0]; new_conectivity(3,2) = edge_ids[2]; 
                return true; 
               }
             else 
               {      
                    //create_element = false;  
                     return false;  
                    //KRATOS_WATCH(" NO CONECTIVITY WAS CRAETED " )
               }            
           }    
    };
}
#endif 
