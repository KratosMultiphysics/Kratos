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

#if !defined(INTRA_FRACTURE_TRIANGLE_UTILITY)
#define INTRA_FRACTURE_TRIANGLE_UTILITY

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



// System includes
#include <string>
#include <iostream> 
#include <stdlib.h>
#include <cmath>
#include <algorithm>



/* Project includes */
#include "structural_application.h"

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
#include "utilities/split_triangle.h"
#include "custom_utilities/sd_math_utils.h"
#include "custom_utilities/smoothing_utility.h"
#include "geometries/triangle_2d_3.h"
#include "processes/node_erase_process.h" 
#include "spatial_containers/spatial_containers.h"


namespace Kratos
{
      class Intra_Fracture_Triangle
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

 


          Intra_Fracture_Triangle(ModelPart& model_part, int domain_size) : mr_model_part(model_part)
          {
             mdomain_size = domain_size;
             mInitialize  = false;  
          }
          
          ~Intra_Fracture_Triangle(){}
    

///************************************************************************************************          
///************************************************************************************************ 
 
///* Detecta los dos elementos en el pano que seran dividos.
///* En total son tres nodos que seran insertados uno de ellos se duplicara.
void Detect_And_Split_Elements(ModelPart& this_model_part, Node<3>::Pointer pNode, const array_1d<double,3>& Failure_Maps)

{
KRATOS_TRY  

    //ProcessInfo& CurrentProcessInfo    =  this_model_part.GetProcessInfo();  
    //ElementsArrayType& pElements       =  this_model_part.Elements(); 
    //NodesArrayType& pNodes             =  this_model_part.Nodes();

    compressed_matrix<int> Coord;
    boost::numeric::ublas::vector<int> List_New_Nodes;                         ///* the news nodes
    boost::numeric::ublas::vector<array_1d<int, 2 > > Position_Node;           ///* edges where are the news nodes
    boost::numeric::ublas::vector< array_1d<double, 3> > Coordinate_New_Node;  ///* the coordinate of the new nodes
    Node<3>::Pointer pduplicate_node;                                          ///* the node to be duplicated  
    WeakPointerVector< Element > Negative_Elements;                           
    WeakPointerVector< Element > Positive_Elements; 
    WeakPointerVector< Element > Splitted_Elements;

    ///* Nelson 
    Initialize(this_model_part);
    Detect_Elements(this_model_part, pNode, Failure_Maps, Splitted_Elements);
    CSR_Row_Matrix(this_model_part, Coord);  
    Calculate_Coordinate_New_Nodes(this_model_part, Failure_Maps, pNode, Splitted_Elements, List_New_Nodes, Coordinate_New_Node, Coord);
    Insert_New_Node_In_Mesh(this_model_part, Coordinate_New_Node, List_New_Nodes, pNode, pduplicate_node);
    Calculate_Negative_And_Positive_Elements(this_model_part, pNode, Failure_Maps,  Negative_Elements, Positive_Elements);    
    Change_Node_Of_NegativeElemensts(this_model_part,Negative_Elements, Splitted_Elements, pNode, pduplicate_node);
    Erase_Old_Element_And_Create_New_Triangle_Element(this_model_part, Coord, Failure_Maps, pNode, pduplicate_node, Splitted_Elements);   
    Renumering_Elements_And_Nodes(this_model_part);


  KRATOS_CATCH("")
}


///************************************************************************************************          
///************************************************************************************************  

inline void CreatePartition(unsigned int number_of_threads, const int number_of_rows, vector<unsigned int>& partitions)
    {
      partitions.resize(number_of_threads+1);
      int partition_size = number_of_rows / number_of_threads;
      partitions[0] = 0;
      partitions[number_of_threads] = number_of_rows;
      for(unsigned int i = 1; i<number_of_threads; i++)
      partitions[i] = partitions[i-1] + partition_size ;
    }


///************************************************************************************************          
///************************************************************************************************  

void Initialize(ModelPart& this_model_part)

{

        ElementsArrayType& pElements       =  this_model_part.Elements(); 

        #ifdef _OPENMP
	int number_of_threads = omp_get_max_threads();
	#else
	int number_of_threads = 1;
	#endif

	vector<unsigned int> element_partition;
	CreatePartition(number_of_threads, pElements.size(), element_partition);  
	#pragma omp parallel for    
	for(int k=0; k<number_of_threads; k++)
	{
	
	  ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
	  ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];
        
	  for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
	   {
                it->SetValue(REFINEMENT_LEVEL, 0);           
           }
       }
}


///************************************************************************************************          
///************************************************************************************************ 

void Detect_Elements(
ModelPart& this_model_part, Node<3>::Pointer pNode, 
const array_1d<double,3>& Failure_Maps,
WeakPointerVector< Element>& Splitted_Elements
)

{
    Splitted_Elements.reserve(1000);
    array_1d<double, 3> normal;
    array_1d<double, 3> unit;
    array_1d<double, 3> Coord_Point_1; 
    array_1d<double, 3> Coord_Point_2;  
    Calculate_Normal_Faliure_Maps(normal, Failure_Maps);

    //WeakPointerVector< Node<3> >& neighb_nodes  = pNode->GetValue(NEIGHBOUR_NODES);
    WeakPointerVector< Element >& neighb_elems  = pNode->GetValue(NEIGHBOUR_ELEMENTS);

    Coord_Point_1[0]  = pNode->X(); 
    Coord_Point_1[1]  = pNode->Y();   
    Coord_Point_1[2]  = pNode->Z();

    for(WeakPointerVector< Element >::iterator neighb_elem = neighb_elems.begin();
    neighb_elem != neighb_elems.end(); neighb_elem++)
    { 
       Element::GeometryType& geom = neighb_elem->GetGeometry(); // Nodos del elemento 
       std::vector<bool> dect(2);
       dect[0] = false; dect[1] = false;  
       unsigned int k = 0;
       for (unsigned int i = 0; i <geom.size(); i++)
 	   {
                if( geom[i].Id() != pNode->Id())
                 {
                        Coord_Point_2[0]  = geom[i].X(); 
                        Coord_Point_2[1]  = geom[i].Y();   
                        Coord_Point_2[2]  = geom[i].Z();
                        noalias(unit)     = Coord_Point_2 - Coord_Point_1;
                        noalias(unit)     = unit / norm_2(unit);
                        if (inner_prod(normal, unit) >= 0.00) { dect[k] = true;}
                        else {dect[k] = false;}   
                        k++;                      
                 }   
           }
         if(dect[0] != dect[1]) 
               { 
                 Splitted_Elements.push_back(*(neighb_elem.base() ));     
               }
          
         } 
       
   
  //std::cout<< "************ ELEMENT TO BE DIVIDED *************"<<std::endl;
  //KRATOS_WATCH(Splitted_Elements) 
  for(WeakPointerVector< Element >::iterator neighb_elem = Splitted_Elements.begin();
    neighb_elem != Splitted_Elements.end(); neighb_elem++)
    { 
      neighb_elem->SetValue(REFINEMENT_LEVEL, 1);
    }


} 

///************************************************************************************************          
///************************************************************************************************  

void Calculate_Coordinate_New_Nodes(
ModelPart& this_model_part,
const array_1d<double,3>&  failure_map, 
const Node<3>::Pointer pNode, ///* the node to be duplicated
WeakPointerVector< Element >& Splitted_Elements,
boost::numeric::ublas::vector<int>& List_New_Nodes,  
boost::numeric::ublas::vector< array_1d<double, 3> >& Coord_Point_New_Node,
compressed_matrix<int>& Coord
)
{

 NodesArrayType& pNodes =  this_model_part.Nodes();  
 int total_node =  pNodes.size();


 double m1 = 0.00;
 double m2 = 0.00;
 array_1d<array_1d<double, 3>, 2 >  Coord_Point; 
 //array_1d<array_1d<double, 3>, 2 >  Coord_Point_New_Node; 
 Coord_Point_New_Node.resize(3);
 
 ///* incluuimos la duplicacion del nodo como el primero
 List_New_Nodes.resize(3);
 List_New_Nodes[0] = total_node + 1;   
 Coord_Point_New_Node[0][0] =  pNode->X();
 Coord_Point_New_Node[0][1] =  pNode->Y();
 Coord_Point_New_Node[0][2] =  pNode->Z();

 ///* Bucle sobre los elementos fracturados 
 unsigned int el = 1;   
 for(WeakPointerVector< Element >::iterator frac_elem = Splitted_Elements.begin();
    frac_elem != Splitted_Elements.end(); frac_elem++)
   {
       List_New_Nodes[el] = total_node + 1 + el;            

         ///* Bucle sobre los nodos del elemento fracturado
       Element::GeometryType& geom = frac_elem->GetGeometry(); // Nodos del elemento
       unsigned int k = 0;
       for (unsigned int i = 0; i <geom.size(); i++ )
 	   {
              if( geom[i].Id() != pNode->Id())
                 {                    
                   Coord_Point[k][0] = geom[i].X();
                   Coord_Point[k][1] = geom[i].Y();
                   Coord_Point[k][2] = geom[i].Z();
                   k++;
                 }
           }

///* calculando los edges
if( geom[0].Id() == pNode->Id())   
      { int index_i = geom[1].Id()-1;
	int index_j = geom[2].Id()-1;  
        int index_k = geom[0].Id()-1;
        Coord(index_k, index_k) = -2;
        if( index_j > index_i) { Coord(index_i, index_j) = List_New_Nodes[el]; }
	else { Coord(index_j, index_i) = List_New_Nodes[el]; }  
      }

if( geom[1].Id() == pNode->Id())   
      { int index_i = geom[2].Id()-1;
	int index_j = geom[0].Id()-1;  
        int index_k = geom[1].Id()-1;
        Coord(index_k, index_k) = -2;              
        if( index_j > index_i) { Coord(index_i, index_j) = List_New_Nodes[el]; }
	else { Coord(index_j, index_i) = List_New_Nodes[el]; }  
      }

if( geom[2].Id() == pNode->Id())   
      { int index_i = geom[0].Id()-1;
	int index_j = geom[1].Id()-1;  
        int index_k = geom[2].Id()-1;
        Coord(index_k, index_k) = -2;          
        if( index_j > index_i) { Coord(index_i, index_j) = List_New_Nodes[el]; }
	else { Coord(index_j, index_i) = List_New_Nodes[el]; }  
      }


///* Calculando ecuacion de recta de elemento fracturado
  bool singular_1 = false;
  bool singular_2 = false;
 

  m1 = (Coord_Point[1][1] - Coord_Point[0][1]) / (Coord_Point[1][0] - Coord_Point[0][0]);
  m2 = tan(failure_map[1]/failure_map[0]); 

  if(fabs(Coord_Point[1][0] - Coord_Point[0][0]) < 1E-9 ) {singular_1 = true;}
  if(fabs(failure_map[0])<1E-9) { singular_2 = true;}
          
     

  double x1 =  Coord_Point[1][0];  
  double y1 =  Coord_Point[1][1];
  double x2 =  pNode->X(); 
  double y2 =  pNode->Y();     
  Coord_Point_New_Node[el][0] = (1.00 / (m1 - m2)) * ( y2 - y1 + m1 * x1 - m2 * x2 );
  Coord_Point_New_Node[el][1] = m1 * (Coord_Point_New_Node[el][0] - x1)  + y1;
  Coord_Point_New_Node[el][2] = 0.00;

  if(singular_1==true)
  {
    Coord_Point_New_Node[el][0] = x1;
    Coord_Point_New_Node[el][1] = m2 *(Coord_Point_New_Node[el][0] - x2)  + y2;
    Coord_Point_New_Node[el][2] = 0.00;
  }

  if(singular_2==true)
  {
    Coord_Point_New_Node[el][0] = x2;
    Coord_Point_New_Node[el][1] = m1 *(Coord_Point_New_Node[el][0] - x1)  + y1;
    Coord_Point_New_Node[el][2] = 0.00;
  }

  el++;

  }
                
}


///************************************************************************************************          
///************************************************************************************************  

void Calculate_Normal_Faliure_Maps(array_1d<double,3>&  normal, const array_1d<double,3>&  failure_map)
{
 
 // WARNING: SOLO 2D
 //double tetha = atan(failure_map[1]/failure_map[0]);
 //tetha += PI/2.00;
 normal[0] = -failure_map[1]; //cos(tetha);     
 normal[1] =  failure_map[0];
 normal[2] = 0.00;
 noalias(normal) = normal / norm_2(normal);

}

///************************************************************************************************          
///************************************************************************************************  

void CSR_Row_Matrix(ModelPart& this_model_part, compressed_matrix<int>& Coord )
{  
    NodesArrayType& pNodes =  this_model_part.Nodes();  
    Coord.resize(pNodes.size(),pNodes.size());    
    NodesArrayType::iterator i_begin=pNodes.ptr_begin();
    NodesArrayType::iterator i_end=pNodes.ptr_end();

    for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)     
	{    
          int index_i = i->Id()-1;
          WeakPointerVector< Node<3> >& neighb_nodes  = i->GetValue(NEIGHBOUR_NODES);   
	  std::vector<unsigned int> aux(neighb_nodes.size());
	  unsigned int active = 0;
	  for(WeakPointerVector< Node<3> >::iterator inode = neighb_nodes.begin();
                inode != neighb_nodes.end(); inode++)
          {    
	        int index_j = inode->Id()-1;
		if(index_j > index_i)
		{
		  aux[active] = index_j;
		   active++;
		 }
	  }
	  std::sort(aux.begin(), aux.begin()+active);

	  for(unsigned int k=0; k<active; k++)
	  {
		Coord.push_back(index_i, aux[k], -1);
			      
	  }

         }
  
}  




///************************************************************************************************          
///************************************************************************************************  

///* inserta los nuevos nodeos en el model part
void Insert_New_Node_In_Mesh(ModelPart& this_model_part,
boost::numeric::ublas::vector< array_1d<double, 3> >& Coordinate_New_Node,
boost::numeric::ublas::vector<int>& List_New_Nodes,
Node<3>::Pointer& pNode,
Node<3>::Pointer& pduplicate_node
)

{
//node to get the DOFs from	
unsigned int size = Coordinate_New_Node.size(); 
unsigned int step_data_size = this_model_part.GetNodalSolutionStepDataSize();
Node<3>::DofsContainerType& reference_dofs = (this_model_part.NodesBegin())->GetDofs();
Node<3>::Pointer pnode; 
 

bool duplicated_node = false;
for(unsigned int i = 0; i<size; i++)
 {
    pnode = this_model_part.CreateNewNode(List_New_Nodes[i],Coordinate_New_Node[i][0], Coordinate_New_Node[i][1], Coordinate_New_Node[i][2] );
    pnode->SetBufferSize(this_model_part.NodesBegin()->GetBufferSize() );
       
   //generating the dofs  
   for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin();    iii != reference_dofs.end(); iii++)  
     {
       Node<3>::DofType& rDof = *iii;
       Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );
       (p_new_dof)->FreeDof();
     }

     unsigned int buffer_size = pnode->GetBufferSize();
     for(unsigned int step = 0; step<buffer_size; step++)
      {
        double* step_data     = pnode->SolutionStepData().Data(step);			
        double* old_node_data = pNode->SolutionStepData().Data(step);
        //copying this data in the position of the vector we are interested in
        for(unsigned int j= 0; j<step_data_size; j++)
        {  
          step_data[j] = old_node_data[j];
        }						
     }
/*
    const array_1d<double,3>& disp = pnode->FastGetSolutionStepValue(DISPLACEMENT);
    pnode->X0() = pnode->X() -  disp[0];
    pnode->Y0() = pnode->Y() -  disp[1];
    pnode->Z0() = pnode->Z() -  disp[2];

    const array_1d<double,3>& vel_old = pNode->FastGetSolutionStepValue(VELOCITY);
    array_1d<double,3>& vel_new       = pnode->FastGetSolutionStepValue(VELOCITY);
    vel_new = vel_old ;

    const array_1d<double,3>& accel_old = pNode->FastGetSolutionStepValue(ACCELERATION);
    array_1d<double,3>& accel_new       = pnode->FastGetSolutionStepValue(ACCELERATION);
    accel_new = accel_old ;
*/


   if( duplicated_node == false)
     { pduplicate_node =  pnode;
       duplicated_node = true;
     }  

  }

}


///************************************************************************************************          
///************************************************************************************************  

void Erase_Old_Element_And_Create_New_Triangle_Element(
           ModelPart& this_model_part,
           const compressed_matrix<int>& Coord,
           const array_1d<double,3>&  Failure_Maps,
           Node<3>::Pointer& pNode, 
           Node<3>::Pointer& pduplicated_node, 
           WeakPointerVector< Element >& Splitted_Elements  
          )   
    {


       boost::numeric::ublas::matrix<int> new_conectivity;    
       ElementsArrayType& pElements         =  this_model_part.Elements(); 
       ElementsArrayType::iterator it_begin =  pElements.ptr_begin();
       ElementsArrayType::iterator it_end   =  pElements.ptr_end(); 
       Element const rReferenceElement;   
       unsigned int to_be_deleted=0;
       unsigned int large_id = (pElements.end()-1)->Id() * 7;
       array_1d<int,3> triangle_ids; 
       array_1d<int,3> edge_ids;
       bool create_element = false; 
       



	unsigned int current_id = (pElements.end()-1)->Id() + 1;
        for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
         { 
        
		Element::GeometryType& geom = it->GetGeometry();             
		triangle_ids[0] = geom[0].Id();     
		triangle_ids[1] = geom[1].Id();
		triangle_ids[2] = geom[2].Id();
		
		int index_i = geom[0].Id()-1;
		int index_j = geom[1].Id()-1;
		int index_k = geom[2].Id()-1;  

// 		edge_ids[0] = Coord(index_i, index_j );
// 		if(index_i > index_j ) {edge_ids[0] = Coord(index_j, index_i );}                 
// 		edge_ids[1] = Coord(index_j, index_k );
// 		if(index_j > index_k ) {edge_ids[1] = Coord(index_k, index_j );}                 
// 		edge_ids[2] = Coord(index_k, index_i);
// 		if(index_k > index_i ) {edge_ids[2] = Coord(index_i, index_k );}   

		if(index_i > index_j) edge_ids[0] = Coord(index_j, index_i );
		else edge_ids[0] = Coord(index_i, index_j );
		
		if(index_j > index_k) edge_ids[1] = Coord(index_k, index_j );
		else edge_ids[1] = Coord(index_j, index_k );
		
		if(index_k > index_i) edge_ids[2] = Coord(index_i, index_k );
		else edge_ids[2] = Coord(index_k, index_i ); 

			                                                      
                ///* WARNING = siemrpre tendremos los 3 primeros casos 
		create_element = Split_Triangle_Elements::Split_Triangle(triangle_ids, edge_ids, new_conectivity);           
		if(create_element==true)
		{
                  to_be_deleted++; 
		  for(unsigned int i=0; i<new_conectivity.size1();i++)
		  {	    
			  
			  Triangle2D3<Node<3> > geom(
						    this_model_part.Nodes()(new_conectivity(i,0))  ,
						    this_model_part.Nodes()(new_conectivity(i,1))  ,
						    this_model_part.Nodes()(new_conectivity(i,2)) 
						    );

			   
			  Element::Pointer p_element = it->Create(current_id, geom, it->pGetProperties());
			  pElements.push_back(p_element); 
                          current_id++;
			  
                         

                     ///* WARNING = Solo al elemento fracturado  
                    ///* cambiando Id para que el elemento creado vaya al nodo duplicado
                    ///* usamos el mismo criterio de elementos positivos(old node ) y negativos(new node) 
		    ///* Normal al plano de fractura
                    ///-------------------------------------------------------------------
                      if( Splitted_Elements[0].Id()==it->Id() or Splitted_Elements[1].Id()==it->Id()) 
		       {
                          array_1d<double,3>  normal;
		          Calculate_Normal_Faliure_Maps(normal, Failure_Maps);
			  //unsigned int Node_old = 0;
			  //unsigned int Node_new = 0;
			  array_1d<double, 3> Coord_Point_1; 
			  array_1d<double, 3> Coord_Point_2; 

			  Coord_Point_1[0]  = pNode->X(); 
			  Coord_Point_1[1]  = pNode->Y();   
			  Coord_Point_1[2]  = pNode->Z();


			  array_1d<double, 3> Unit;   
			  double prod  = 0.00; 
			  Element::GeometryType& geom_new = p_element->GetGeometry(); // Nodos del elemento nuevo 
			  Find_Coord_Gauss_Points(geom_new, Coord_Point_2);
			  noalias(Unit) = Coord_Point_2 - Coord_Point_1;
			  noalias(Unit) = Unit / norm_2(Unit);
			  prod = inner_prod(normal, Unit);
			  if( prod<0.00)
			      {
				for (unsigned int i = 0; i <geom.size(); i++)
				  {         
				      if (geom_new[i].Id()==pNode->Id())
				      { 
					geom_new(i) = pduplicated_node;
				      }      
				  }
			      }
      
                    }
                  
		  } 
                  it->SetId(large_id);     
                  large_id++;    
		}
	      }
          
 
          //interpolate internal variables

	  this_model_part.Elements().Sort(); //all of the elements to be erased are at the end
          
          ///*now remove all of the "old" elements
	  this_model_part.Elements().erase(this_model_part.Elements().end()-to_be_deleted, this_model_part.Elements().end());
					
  }



///************************************************************************************************          
///************************************************************************************************  


void Calculate_Negative_And_Positive_Elements(
ModelPart& this_model_part, 
Node<3>::Pointer pNode, 
const array_1d<double,3>& Failure_Maps,   
WeakPointerVector< Element >& Negative_Elements,
WeakPointerVector< Element >& Positive_Elements       
)   
{

KRATOS_TRY 

Positive_Elements.reserve(10); 
Negative_Elements.reserve(10); 

//WeakPointerVector< Node<3> >& neighb_nodes  = pNode->GetValue(NEIGHBOUR_NODES);
WeakPointerVector< Element >& neighb_elems  = pNode->GetValue(NEIGHBOUR_ELEMENTS);

// Normal al plano de fractura
array_1d<double,3>  normal;
Calculate_Normal_Faliure_Maps(normal, Failure_Maps);


array_1d<double, 3> Coord_Point_1; 
array_1d<double, 3> Coord_Point_2; 

Coord_Point_1[0]  = pNode->X(); 
Coord_Point_1[1]  = pNode->Y();   
Coord_Point_1[2]  = pNode->Z();


array_1d<double, 3> Unit;   
double prod  = 0.00;
for(WeakPointerVector< Element >::iterator neighb_elem = neighb_elems.begin();
neighb_elem != neighb_elems.end(); neighb_elem++)
 { 
   Element::GeometryType& geom = neighb_elem->GetGeometry(); // Nodos del elemento 
   Find_Coord_Gauss_Points(geom, Coord_Point_2);
   noalias(Unit) = Coord_Point_2 - Coord_Point_1;
   noalias(Unit) = Unit / norm_2(Unit);
   prod = inner_prod(normal, Unit);
   if( prod>=0.00)
     {
        KRATOS_WATCH("INSERTED_POSITIVE")
        Positive_Elements.push_back(*(neighb_elem.base() ));
     }
   else
     {
        KRATOS_WATCH("INSERTED_NEGATIVE")
        Negative_Elements.push_back(*(neighb_elem.base() ));
     }    

   Unit   = ZeroVector(3);
   prod   = 0.00;  
 
 }
KRATOS_CATCH("")
}

///************************************************************************************************          
///************************************************************************************************  

void Change_Node_Of_NegativeElemensts(
ModelPart& this_model_part,
WeakPointerVector< Element >& Negative_Elements,
WeakPointerVector< Element > Splitted_Elements,
Node<3>::Pointer& pNode,
Node<3>::Pointer& pduplicate_node
)
{

///* setting the node of the negative elements for the nw node
///* putting de new node to the negative element  

array_1d<unsigned int, 2 > Id_Split_Elem;
unsigned int k = 0; 
 
for(WeakPointerVector< Element >::iterator split_elem = Splitted_Elements.begin();
split_elem != Splitted_Elements.end(); split_elem++)
{
  Id_Split_Elem[k] = split_elem->Id();
  k++;  
}

for(WeakPointerVector< Element >::iterator neg_elem = Negative_Elements.begin();
neg_elem != Negative_Elements.end(); neg_elem++)
{
  if(neg_elem->Id()!= Id_Split_Elem[0] && neg_elem->Id()!= Id_Split_Elem[1]) {
    Element::GeometryType& geom = neg_elem->GetGeometry();
    for (unsigned int i = 0; i <geom.size(); i++)
 	{         
           if (geom[i].Id()==pNode->Id())
             { 
               geom(i) = pduplicate_node;
             }      
        }
   }
      
}


}




///************************************************************************************************          
///************************************************************************************************  


void Find_Coord_Gauss_Points(Element::GeometryType& geom, array_1d<double,3>&  Coord_Point)
{
        double x = 0.00;   
        double y = 0.00;
        double z = 0.00;
        double fact = 0.33333333333333333333333;

        if (geom.size()==4) {fact = 0.25;}  
        Coord_Point = ZeroVector(3);            
   	for (unsigned int i = 0; i <geom.size(); i++)
 	{
        
             x = geom[i].X();
             y = geom[i].Y();
             z = geom[i].Z();

             Coord_Point[0] += x;
             Coord_Point[1] += y;
             Coord_Point[2] += z;             
        }
       
        noalias(Coord_Point) = Coord_Point*fact;       

} 



void Renumering_Elements_And_Nodes( ModelPart& this_model_part)
{
    
   unsigned int id_node =  1;
   unsigned int id_elem =  1;
   NodesArrayType& pNodes =  this_model_part.Nodes();  
   NodesArrayType::iterator i_begin=pNodes.ptr_begin();
   NodesArrayType::iterator i_end=pNodes.ptr_end();
   //ProcessInfo& rCurrentProcessInfo  = this_model_part.GetProcessInfo();  


   for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)     
	{
           if(i->Id()!=id_node)
              {  //std::cout<< "Setting Id of Node  " << i->Id() << " by " <<  id_node << std::endl; 
                 i->SetId(id_node);
              }
            id_node++;     
        }

  ElementsArrayType& rElements         =  this_model_part.Elements(); 
  ElementsArrayType::iterator it_begin =  rElements.ptr_begin();
  ElementsArrayType::iterator it_end   =  rElements.ptr_end(); 

  for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
        {
          if(it->Id()!=id_elem)
          { 
            //std::cout<< "Setting Id of elelemnt  " << it->Id() << " by " <<  id_elem << std::endl;  
            it->SetId(id_elem);
          }
          id_elem++;
        }
}




bool mInitialize;
ModelPart& mr_model_part;
unsigned int mdomain_size;
WeakPointerVector< Element > mSplitted_Elements;


};
}
#endif 

