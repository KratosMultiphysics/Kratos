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

          typedef WeakPointerVector< Node<3> >::iterator Vec_Node_It; 
 


          Intra_Fracture_Triangle(ModelPart& model_part, int domain_size) :  mr_model_part(model_part)
          {
             mdomain_size = domain_size;
             mInitialize  = false;  
          }
          
          ~Intra_Fracture_Triangle(){}
 

///************************************************************************************************          
///************************************************************************************************ 
 
///* Detecta los dos elementos en el pano que seran dividos.
///* En total son tres nodos que seran insertados uno de ellos se duplicara.
void Detect_And_Split_Elements( bool refine_on_reference) //Node<3>::Pointer pNode, const array_1d<double,3>& Failure_Maps)

{
KRATOS_TRY  

    //ProcessInfo& CurrentProcessInfo    =  this_model_part.GetProcessInfo();  
    //ElementsArrayType& pElements       =  this_model_part.Elements(); 
    //NodesArrayType& pNodes             =  this_model_part.Nodes();
    // boost::numeric::ublas::vector<int> List_New_Nodes;


    ModelPart& this_model_part =  mr_model_part; 
    compressed_matrix<int> Coord;
    array_1d<int, 3>  List_New_Nodes;                                          ///* the news nodes
    array_1d<array_1d<int, 2 >, 2 > Position_Node;                             ///* edges where are the news nodes
    array_1d< array_1d<double, 3>, 3 > Coordinate_New_Node;                    ///* the coordinate of the new nodes
    Node<3>::Pointer pduplicate_node;                                          ///* the node to be duplicated  
    WeakPointerVector< Element > Negative_Elements;                           
    WeakPointerVector< Element > Positive_Elements; 
    WeakPointerVector< Element > Splitted_Elements;
    
   

    FindElementalNeighboursProcess ElementosVecinos(this_model_part, 2, 10);
    FindNodalNeighboursProcess     NodosVecinos(this_model_part, 2, 10);
     
    
    Initialize(this_model_part, refine_on_reference);
    Detect_Node_To_Be_Splitted(this_model_part);
     
      
    array_1d<double,3>  Failure_Maps;  
    for(WeakPointerVector< Node<3> >::iterator inode = mFail_Node.begin();
    inode != mFail_Node.end(); inode++)
      {
         
         Calculate_Map_Failure(inode, Failure_Maps);
         Detect_Elements(this_model_part, inode, Failure_Maps, Splitted_Elements); 
         Calculate_Negative_And_Positive_Elements(this_model_part, inode, Failure_Maps,  Negative_Elements, Positive_Elements);  
         

         CSR_Row_Matrix(this_model_part, Coord);  
         Calculate_Coordinate_New_Nodes(this_model_part, Failure_Maps, inode, Splitted_Elements, List_New_Nodes, Position_Node, Coordinate_New_Node, Coord);
         Insert_New_Node_In_Mesh(this_model_part, Coordinate_New_Node, Position_Node, List_New_Nodes, inode, pduplicate_node);
 
    
      
         Change_Node_Of_Negative_Elemensts(this_model_part,Negative_Elements, Splitted_Elements, inode, pduplicate_node);    
         Erase_Old_Element_And_Create_New_Triangle_Element(this_model_part, Coord, Failure_Maps, inode, pduplicate_node, Splitted_Elements);  
                   
      
         ElementosVecinos.ClearNeighbours();
         NodosVecinos.ClearNeighbours();
         ElementosVecinos.Execute(); 
         NodosVecinos.Execute();           
         Failure_Maps = ZeroVector(3);           
         Splitted_Elements.clear();         
      }
    
  
  
  Finalize(this_model_part, refine_on_reference);   

  //WeakPointerVector< Element >& neighb_elems  = pduplicate_node->GetValue(NEIGHBOUR_ELEMENTS);
  //KRATOS_WATCH( neighb_elems) 
  
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

void Initialize(ModelPart& this_model_part, bool refine_on_reference)

{

        ElementsArrayType& rElements         =  this_model_part.Elements(); 



        #ifdef _OPENMP
	int number_of_threads = omp_get_max_threads();
	#else
	int number_of_threads = 1;
	#endif

	vector<unsigned int> element_partition;
	CreatePartition(number_of_threads, rElements.size(), element_partition);  
	#pragma omp parallel for    
	for(int k=0; k<number_of_threads; k++)
	{
	
	  ElementsArrayType::iterator it_begin=rElements.ptr_begin()+element_partition[k];
	  ElementsArrayType::iterator it_end=rElements.ptr_begin()+element_partition[k+1];
        
	  for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
	   {
                it->SetValue(REFINEMENT_LEVEL, 0);           
           }
       }


      if(refine_on_reference==true)
        {
          for(ModelPart::NodesContainerType::iterator it= this_model_part.NodesBegin(); it!=this_model_part.NodesEnd(); it++)
              {
	       it->X()=it->X0();
	       it->Y()=it->Y0();
	       it->Z()=it->Z0();
             }
       }


}


void Finalize(ModelPart& this_model_part, bool refine_on_reference)
{
Renumering_Elements_And_Nodes(this_model_part);

if(refine_on_reference==true)
  {
     for(ModelPart::NodesContainerType::iterator it=this_model_part.NodesBegin(); it!=this_model_part.NodesEnd(); it++)
      {     
	     const array_1d<double,3>& disp = it->FastGetSolutionStepValue(DISPLACEMENT);
	     it->X()=it->X0() + disp[0];
	     it->Y()=it->Y0() + disp[1];
 	     it->Z()=it->Z0() + disp[2];
          }
  } 

}


///************************************************************************************************          
///************************************************************************************************      
void Detect_Node_To_Be_Splitted(ModelPart& this_model_part)
{
KRATOS_TRY  
/*
NodesArrayType& pNodes = this_model_part.Nodes();

#ifdef _OPENMP
int number_of_threads = omp_get_max_threads();
#else
int number_of_threads = 1;
#endif

mFail_Node.reserve(1000);

vector<unsigned int> node_partition;
CreatePartition(number_of_threads, pNodes.size(), node_partition);

#pragma omp parallel for  
for(int k=0; k<number_of_threads; k++)
{
NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

  for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)     
  {
    double& Condition = i->FastGetSolutionStepValue(NODAL_DAMAGE);
    if(Condition >= 1.00)
	  {  
	      i->FastGetSolutionStepValue(SPLIT_NODAL) = true;   
              mFail_Node.push_back(*(i.base()) );  
	  }
  }

}
  */



 for( 
     ModelPart::NodesContainerType::iterator inode = this_model_part.NodesBegin();
     inode != this_model_part.NodesEnd();
     inode++)	
      {
         bool & split = inode->GetValue(SPLIT_NODAL);  
	 if ( split == true) { 
              //KRATOS_WATCH(inode->Id())        
              mFail_Node.push_back(*(inode.base()) );  
	  }
     }

KRATOS_CATCH("")
}



///************************************************************************************************          
///************************************************************************************************      
///* Computa el vector, linea, o plano de falla.
void Calculate_Map_Failure(
Vec_Node_It& pNode, //Node<3>::Pointer pNode,  
array_1d<double,3>& Failure_Maps)
{
/*
Vector Eigen_Values  = ZeroVector(mdomain_size);               // Deformaciones o Tensiones principales
Matrix Eigen_Vector  = ZeroMatrix(mdomain_size,mdomain_size);  // Direcciones  principales  
Matrix Strain_Tensor = ZeroMatrix(mdomain_size,mdomain_size);
 
unsigned int size = 3;
unsigned int max_iterations    = 100;
double zero_tolerance = 1e-9;
if(mdomain_size==3){size = 6; }
Matrix Nodal_Values  = ZeroMatrix(1,size);
Failure_Maps         = ZeroVector(3);
Vector_Order_Tensor EigenVector;

EigenVector.resize(3);
EigenVector[0] = ZeroVector(3);
EigenVector[1] = ZeroVector(3);
EigenVector[2] = ZeroVector(3);


Nodal_Values  = pNode->GetValue(NODAL_STRAIN);
Vector temp   = ZeroVector(size);
for( unsigned int j=0; j<size; j++){temp[j] = Nodal_Values(0, j);} 

Strain_Tensor = SD_MathUtils<double>::StrainVectorToTensor(temp);
SD_MathUtils<double>::EigenVectors(Strain_Tensor, Eigen_Vector, Eigen_Values, zero_tolerance, max_iterations ); 

// traccion principal de traccion
double max =  (*std::max_element(Eigen_Values.begin(),Eigen_Values.end()));
int pos = 0;

for(unsigned int k=0; k<Eigen_Values.size(); k++ )     
{ if(max==Eigen_Values(k)) { pos = k; }
  for(unsigned int j=0; j<Eigen_Values.size(); j++ ) 
 {
  EigenVector[k][j] = Eigen_Vector(k,j);
 }
}

//Tomo el primer egenvector correspondiente al vector que es perpendicaular al vector de traccion principal mayor
if(pos==0) {pos=1;}
else
  pos = 0;

Failure_Maps[0] =  EigenVector[pos][0];
Failure_Maps[1] =  EigenVector[pos][1];
Failure_Maps[2] =  EigenVector[pos][2];
*/
  
Failure_Maps[0] =  1.00; //sin(PI/6.00);
Failure_Maps[1] =  0.00; //cos(PI/6.00);
Failure_Maps[2] =  0.00;
}

///************************************************************************************************          
///************************************************************************************************ 

void Detect_Elements(
ModelPart& this_model_part, 
Vec_Node_It& pNode, // Node<3>::Pointer pNode, 
const array_1d<double,3>& Failure_Maps,
WeakPointerVector< Element>& Splitted_Elements
)

{
    double toler  = 1E-2;  
    Splitted_Elements.reserve(10);
    array_1d<double, 3> normal;
    array_1d<double, 3> unit;
    array_1d<double, 3> Coord_Point_1; 
    array_1d<double, 3> Coord_Point_2;  
   
    //WeakPointerVector< Node<3> >& neighb_nodes  = pNode->GetValue(NEIGHBOUR_NODES);
    WeakPointerVector< Element >& neighb_elems  = pNode->GetValue(NEIGHBOUR_ELEMENTS);
    Calculate_Normal_Faliure_Maps(normal, Failure_Maps); 



    Coord_Point_1[0]  = pNode->X0(); 
    Coord_Point_1[1]  = pNode->Y0();   
    Coord_Point_1[2]  = pNode->Z0();

    std::vector<bool> dect(2);
    std::vector< double > prod(2);    
    for(WeakPointerVector< Element >::iterator neighb_elem = neighb_elems.begin();
    neighb_elem != neighb_elems.end(); neighb_elem++)
    { 

      
       Element::GeometryType& geom = neighb_elem->GetGeometry(); // Nodos del elemento 
       dect[0] = false; dect[1] = false;  
       prod[0] = 0.000; prod[1] = 0.000;     
       unsigned int k = 0;
       for (unsigned int i = 0; i <geom.size(); i++)
 	   {
                if( geom[i].Id() != pNode->Id())
                 {
                        Coord_Point_2[0]  = geom[i].X0(); 
                        Coord_Point_2[1]  = geom[i].Y0();   
                        Coord_Point_2[2]  = geom[i].Z0();
                        noalias(unit)     = Coord_Point_2 - Coord_Point_1;
                        noalias(unit)     = unit / norm_2(unit);
                        prod[k] = inner_prod(normal, unit);
                        if(prod[k] > 0.00) { dect[k] = true;}
                        else {dect[k] = false;}
                              
                        k++;                      
                 }   
           }
        
 
         if( fabs(prod[0])>=toler && fabs(prod[1])>=toler)
              {
                 if(dect[0] != dect[1]) 
                  { 
                   Splitted_Elements.push_back(*(neighb_elem.base() ));     
                  }
         } 
          
         } 
      
   
  //std::cout<< "************ ELEMENT TO BE DIVIDED *************"<<std::endl;
  for(WeakPointerVector< Element >::iterator it = Splitted_Elements.begin();
    it != Splitted_Elements.end(); it++)
    { 
      std::cout<<"Element To Be Splitted = "<<  it->Id() << std::endl;
      it->SetValue(REFINEMENT_LEVEL, 1);
    }

} 

///************************************************************************************************          
///************************************************************************************************  

void Calculate_Coordinate_New_Nodes(
ModelPart& this_model_part,
const array_1d<double,3>&  failure_map, 
const Vec_Node_It& pNode,
WeakPointerVector< Element >& Splitted_Elements,
array_1d<int, 3>&  List_New_Nodes,
array_1d<array_1d<int, 2 >, 2 >& Position_Node,
array_1d< array_1d<double, 3>, 3 >& Coordinate_New_Node,
compressed_matrix<int>& Coord
)

{


 List_New_Nodes[0] = 0;   
 List_New_Nodes[1] = 0;
 List_New_Nodes[2] = 0;
 
 Position_Node[0][0] = 0; Position_Node[0][1] = 0;
 Position_Node[1][0] = 0; Position_Node[1][1] = 0; 
 
 Coordinate_New_Node[0] = ZeroVector(3);
 Coordinate_New_Node[1] = ZeroVector(3);
 Coordinate_New_Node[2] = ZeroVector(3);


 NodesArrayType& pNodes =  this_model_part.Nodes();  
 int total_node =  pNodes.size();
 //unsigned int number_of_new_nodes = 2; // dos por cada elemento fracturado


 double m1 = 0.00;
 double m2 = 0.00;
 array_1d<array_1d<double, 3>, 2>  Coord_Point; 


 ///* incluuimos la duplicacion del nodo como el primero
 List_New_Nodes.resize(3);
 List_New_Nodes[0] = total_node + 1;   
 Coordinate_New_Node[0][0] =  pNode->X0();
 Coordinate_New_Node[0][1] =  pNode->Y0();
 Coordinate_New_Node[0][2] =  pNode->Z0();



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
                   Coord_Point[k][0] = geom[i].X0();
                   Coord_Point[k][1] = geom[i].Y0();
                   Coord_Point[k][2] = geom[i].Z0();
                   k++;
                 }
           }

///* calculando los edges
if( geom[0].Id() == pNode->Id())   
      { int index_i = geom[1].Id()-1;
	int index_j = geom[2].Id()-1;  
        //int index_k = geom[0].Id()-1;
        //Coord(index_k, index_k) = -2;
        if( index_j > index_i) 
            { 
              Coord(index_i, index_j) = List_New_Nodes[el];  
              Position_Node[el-1][0]  = index_i+1;      
              Position_Node[el-1][1]  = index_j+1;    
            }
	else {
              Coord(index_j, index_i) = List_New_Nodes[el];
              Position_Node[el-1][0]  = index_j+1;      
              Position_Node[el-1][1]  = index_i+1;    
             }  
      }

if( geom[1].Id() == pNode->Id())   
      { int index_i = geom[2].Id();
	int index_j = geom[0].Id();  
        //int index_k = geom[1].Id()-1;
        //Coord(index_k, index_k) = -2;              
        if( index_j > index_i) 
            { 
              Coord(index_i-1, index_j-1) = List_New_Nodes[el];  
              Position_Node[el-1][0]  = index_i;      
              Position_Node[el-1][1]  = index_j;    
            }
	else {
              Coord(index_j-1, index_i-1) = List_New_Nodes[el];
              Position_Node[el-1][0]  = index_j;      
              Position_Node[el-1][1]  = index_i;    
             }  
      }

if( geom[2].Id() == pNode->Id())   
      { int index_i = geom[0].Id()-1;
	int index_j = geom[1].Id()-1;  
        //int index_k = geom[2].Id()-1;
        //Coord(index_k, index_k) = -2;          
        if( index_j > index_i) 
            { 
              Coord(index_i, index_j) = List_New_Nodes[el];  
              Position_Node[el-1][0]  = index_i+1;      
              Position_Node[el-1][1]  = index_j+1;    
            }
	else {
              Coord(index_j, index_i) = List_New_Nodes[el];
              Position_Node[el-1][0]  = index_j+1;      
              Position_Node[el-1][1]  = index_i+1;    
             }   
      }


///* Calculando ecuacion de recta de elemento fracturado
  bool singular_1 = false;
  bool singular_2 = false;
 

  m1 = (Coord_Point[1][1] - Coord_Point[0][1]) / (Coord_Point[1][0] - Coord_Point[0][0]);
  m2 = (failure_map[1]/failure_map[0]); 
 

  if(fabs(Coord_Point[1][0] - Coord_Point[0][0]) < 1E-9 ) {singular_1 = true;}
  if(fabs(failure_map[0])<1E-9) { singular_2 = true;}
          
     

  double x1 =  Coord_Point[1][0];  
  double y1 =  Coord_Point[1][1];
  double x2 =  pNode->X0(); 
  double y2 =  pNode->Y0();     
  Coordinate_New_Node[el][0] = (1.00 / (m1 - m2)) * ( y2 - y1 + m1 * x1 - m2 * x2 );
  Coordinate_New_Node[el][1] = m1 * (Coordinate_New_Node[el][0] - x1)  + y1;
  Coordinate_New_Node[el][2] = 0.00;

  if(singular_1==true)
  {
    Coordinate_New_Node[el][0] = x1;
    Coordinate_New_Node[el][1] = m2 *(Coordinate_New_Node[el][0] - x2)  + y2;
    Coordinate_New_Node[el][2] = 0.00;
  }

  if(singular_2==true)
  {
    Coordinate_New_Node[el][0] = x2;
    Coordinate_New_Node[el][1] = m1 *(Coordinate_New_Node[el][0] - x1)  + y1;
    Coordinate_New_Node[el][2] = 0.00;
  }

  el++;

  }


//KRATOS_WATCH(pNode->Id())
//KRATOS_WATCH(List_New_Nodes)
//KRATOS_WATCH(Coordinate_New_Node)
//KRATOS_WATCH(Position_Node)
//KRATOS_WATCH(failure_map)


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
array_1d< array_1d<double, 3>, 3 >& Coordinate_New_Node,
array_1d<array_1d<int, 2 >, 2 >& Position_Node,
array_1d<int, 3>&  List_New_Nodes,
Vec_Node_It& pNode,
Node<3>::Pointer& pduplicate_node
)

{

//node to get the DOFs from	
unsigned int size = Position_Node.size();
unsigned int step_data_size = this_model_part.GetNodalSolutionStepDataSize();
Node<3>::DofsContainerType& reference_dofs = (this_model_part.NodesBegin())->GetDofs();
Node<3>::Pointer pnode; 
 
array_1d<double, 3 > Coord_Node_1; 
array_1d<double, 3 > Coord_Node_2;
array_1d<double, 3 > Coord_New_Node;
array_1d<double, 3 > dist;
double dist_1 = 0.00;
double dist_2 = 0.00;
double ratio  = 0.00;
const double zero = 0.00;


///* Para duplicated node
const int& node_i     = pNode->Id();
ModelPart::NodesContainerType::iterator it_node1 =  this_model_part.Nodes().find(node_i);
std::size_t pos1 = it_node1 - this_model_part.NodesBegin();
noalias(Coord_Node_1) = it_node1->Coordinates(); 

///* El primero es el nodo duplicado 
//noalias(Coord_New_Node) = Coordinate_New_Node[0];    
int    id   = List_New_Nodes[0]; 
double  x   = Coordinate_New_Node[0][0];
double  y   = Coordinate_New_Node[0][1];
double  z   = Coordinate_New_Node[0][2];

pduplicate_node = this_model_part.CreateNewNode(id, x, y, z);
pduplicate_node->SetBufferSize(this_model_part.NodesBegin()->GetBufferSize() );

//KRATOS_WATCH(*it_node1)
it_node1 = this_model_part.NodesBegin() + pos1;


///WARNING = si lo activo no me aparece el nodo duplicado consultar rossi
pduplicate_node->X0() = it_node1->X0(); 
pduplicate_node->Y0() = it_node1->Y0();
pduplicate_node->Z0() = it_node1->Z0();  



//generating the dofs  
for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin();    iii != reference_dofs.end(); iii++)  
 {
   Node<3>::DofType& rDof = *iii;
   Node<3>::DofType::Pointer p_new_dof =  pduplicate_node->pAddDof( rDof );
   if(it_node1->IsFixed(iii->GetVariable()) == true) 
     (p_new_dof)->FixDof();
      else
    { (p_new_dof)->FreeDof();
    }      
 }


///* intepolating the data
unsigned int buffer_size = pduplicate_node->GetBufferSize();
for(unsigned int step = 0; step<buffer_size; step++)
{
  double* new_step_data    = pduplicate_node->SolutionStepData().Data(step);
  double* step_data1       = it_node1->SolutionStepData().Data(step);			
 
  ///*copying this data in the position of the vector we are interested in 
   for(unsigned int j= 0; j<step_data_size; j++)
    {  
      new_step_data[j] = step_data1[j];
    }						
}


/// WARNING =  only for reactions;

for(Node<3>::DofsContainerType::iterator iii = pduplicate_node->GetDofs().begin();    iii != pduplicate_node->GetDofs().end(); iii++)  
{          
   if(pduplicate_node->IsFixed(iii->GetVariable())==false) 
    { 
       iii->GetSolutionStepReactionValue() = zero;
     }       
}  


Coord_New_Node = ZeroVector(3);
///* Para dos nodos
for(unsigned int i = 0; i<size; i++)
 {

   if(Position_Node[i][0]!=0 && Position_Node[i][1]!=0) 
   {

   const int& node_i     = Position_Node[i][0];
   const int& node_j     = Position_Node[i][1];
   KRATOS_WATCH(node_i)
   KRATOS_WATCH(node_j)

   ModelPart::NodesContainerType::iterator it_node1 =  this_model_part.Nodes().find(node_i);
   std::size_t pos1 = it_node1 - this_model_part.NodesBegin();
   noalias(Coord_Node_1) = it_node1->Coordinates(); 
   ModelPart::NodesContainerType::iterator it_node2 =  this_model_part.Nodes().find(node_j);
   std::size_t pos2 = it_node2 - this_model_part.NodesBegin();
   noalias(Coord_Node_2) = it_node2->Coordinates(); 

   ///* El primero es el nodo duplicado 
   noalias(Coord_New_Node) = Coordinate_New_Node[i+1];    
   int  id    = List_New_Nodes[i+1]; 
   double  x  = Coordinate_New_Node[i+1][0];
   double  y  = Coordinate_New_Node[i+1][1];
   double  z  = Coordinate_New_Node[i+1][2];
   double  toler = 1E-3; 


   noalias(dist) =   (Coord_Node_1 - Coord_Node_2);     
   dist_1        =   sqrt(inner_prod(dist, dist ));
   noalias(dist) =   Coord_Node_1 - Coord_New_Node; 
   dist_2        =   sqrt(inner_prod(dist, dist ) );   
   ratio         =   dist_2 / dist_1;

   if( (1.00 - ratio) <= toler )
   {  
   pnode = this_model_part.CreateNewNode(id, x, y, z);
   pnode->SetBufferSize(this_model_part.NodesBegin()->GetBufferSize() );


   
   it_node1 = this_model_part.NodesBegin() + pos1;
   it_node2 = this_model_part.NodesBegin() + pos2; 
  
   pnode->X0() = (1.00 - ratio ) * ( it_node1->X0()) + (ratio) * ( it_node2->X0()) ;
   pnode->Y0() = (1.00 - ratio ) * ( it_node1->Y0()) + (ratio) * ( it_node2->Y0()) ;
   pnode->Z0() = (1.00 - ratio ) * ( it_node1->Z0()) + (ratio) * ( it_node2->Z0()) ;  
   
      
   //generating the dofs  
   for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin();    iii != reference_dofs.end(); iii++)  
     {
       Node<3>::DofType& rDof = *iii;
       Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );
        if(it_node1->IsFixed(iii->GetVariable()) == true && it_node2->IsFixed(iii->GetVariable()) == true)
          (p_new_dof)->FixDof();
        else
          { (p_new_dof)->FreeDof();}    
     }


     ///* intepolating the data
     unsigned int buffer_size = pnode->GetBufferSize();
     for(unsigned int step = 0; step<buffer_size; step++)
      {
	double* new_step_data    = pnode->SolutionStepData().Data(step);
        double* step_data1       = it_node1->SolutionStepData().Data(step);			
        double* step_data2       = it_node2->SolutionStepData().Data(step); 
        ///*copying this data in the position of the vector we are interested in
        for(unsigned int j= 0; j<step_data_size; j++)
        {  
          new_step_data[j] = (1.00 -ratio) *(step_data1[j]) +  (ratio) * (step_data2[j]);
        }						
     }
    
      /// WARNING =  only for reactions;
      for(Node<3>::DofsContainerType::iterator iii = pnode->GetDofs().begin();    iii != pnode->GetDofs().end(); iii++)  
       {          
          if(pnode->IsFixed(iii->GetVariable())==false) 
                 { 
                    iii->GetSolutionStepReactionValue() = zero;
                  }       
       }  

        Coord_New_Node = ZeroVector(3); 
       
       }
        
    }
   }
   
 
}


///************************************************************************************************          
///************************************************************************************************  

void Erase_Old_Element_And_Create_New_Triangle_Element(
           ModelPart& this_model_part,
           const compressed_matrix<int>& Coord,
           const array_1d<double,3>&  Failure_Maps,
           Vec_Node_It& pNode, // Node<3>::Pointer& pNode, 
           Node<3>::Pointer& pduplicated_node, 
           WeakPointerVector< Element >& Splitted_Elements  
          )   
    {


       boost::numeric::ublas::matrix<int> new_conectivity;    
       ElementsArrayType& rElements         =  this_model_part.Elements(); 
       Element const rReferenceElement;   
       unsigned int to_be_deleted=0;
       unsigned int large_id = (rElements.end()-1)->Id() * 7;
       array_1d<int,3> triangle_ids; 
       array_1d<int,3> edge_ids;
       bool create_element = false; 
       PointerVector< Element > New_Elements; 
       ProcessInfo& rCurrentProcessInfo  = this_model_part.GetProcessInfo();      
       
	unsigned int current_id = (rElements.end()-1)->Id() + 1;

         for( 
         ModelPart::ElementsContainerType::iterator it = this_model_part.ElementsBegin();
         it != this_model_part.ElementsEnd();it++)
         { 
                  
                Element::GeometryType& geom = it->GetGeometry();  
		triangle_ids[0] = geom[0].Id();    
                triangle_ids[1] = geom[1].Id();
                triangle_ids[2] = geom[2].Id();
                
                
		int index_i = geom[0].Id()-1;
		int index_j = geom[1].Id()-1;
		int index_k = geom[2].Id()-1;  
                  
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

			   Element::Pointer p_element;  
 		           p_element = it->Create(current_id, geom, it->pGetProperties());
                           New_Elements.push_back(p_element);   
                           current_id++;
			   
                         

                     ///* WARNING = Solo al elemento fracturado  
                    ///* cambiando Id para que el elemento creado vaya al nodo duplicado
                    ///* usamos el mismo criterio de elementos positivos(old node ) y negativos(new node) 
		    ///* Normal al plano de fractura
                    ///-------------------------------------------------------------------
                      if( Splitted_Elements[0].Id()==it->Id() || Splitted_Elements[1].Id()==it->Id()) 
		       {
                          array_1d<double,3>  normal;
		          Calculate_Normal_Faliure_Maps(normal, Failure_Maps);
			  //unsigned int Node_old = 0;
			  //unsigned int Node_new = 0;
			  array_1d<double, 3> Coord_Point_1; 
			  array_1d<double, 3> Coord_Point_2; 

			  Coord_Point_1[0]  = pNode->X0(); 
			  Coord_Point_1[1]  = pNode->Y0();   
			  Coord_Point_1[2]  = pNode->Z0();


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
                  new_conectivity.resize(0,0);    
		}
	      }
          
          
          for(PointerVector< Element >::iterator it_new = New_Elements.begin(); it_new!=New_Elements.end(); it_new++)
            {                     
               //KRATOS_WATCH(it_new->Id());  
               it_new->Initialize();  
               it_new->InitializeSolutionStep(rCurrentProcessInfo);      
               it_new->FinalizeSolutionStep(rCurrentProcessInfo);           
               rElements.push_back(*(it_new.base()));
            } 


          ///* all of the elements to be erased are at the end
	  rElements.Sort(); 
          
          ///*now remove all of the "old" elements
	  rElements.erase(this_model_part.Elements().end()-to_be_deleted, this_model_part.Elements().end());
  	
  }



///************************************************************************************************          
///************************************************************************************************  


void Calculate_Negative_And_Positive_Elements(
ModelPart& this_model_part, 
Vec_Node_It& pNode, //Node<3>::Pointer pNode, 
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

Coord_Point_1[0]  = pNode->X0(); 
Coord_Point_1[1]  = pNode->Y0();   
Coord_Point_1[2]  = pNode->Z0();


array_1d<double, 3> Unit;   
double prod  = 0.00;
for(WeakPointerVector< Element >::iterator neighb_elem = neighb_elems.begin();
neighb_elem != neighb_elems.end(); neighb_elem++)
 { 
   KRATOS_WATCH(neighb_elem->Id()) 
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

void Change_Node_Of_Negative_Elemensts(
ModelPart& this_model_part,
WeakPointerVector< Element >& Negative_Elements,
WeakPointerVector< Element >& Splitted_Elements,
Vec_Node_It& pNode,  //Node<3>::Pointer& pNode,
Node<3>::Pointer& pduplicate_node
)
{

///* setting the node of the negative elements for the nw node
///* putting de new node to the negative element  

array_1d<unsigned int, 2 > Id_Split_Elem;
unsigned int k = 0; 
 
//KRATOS_WATCH(*pduplicate_node)
//KRATOS_WATCH(Splitted_Elements)

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
           //KRATOS_WATCH(*geom(i))    
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
        
             x = geom[i].X0();
             y = geom[i].Y0();
             z = geom[i].Z0();

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


private:

bool mInitialize;
ModelPart& mr_model_part;
unsigned int mdomain_size;
WeakPointerVector< Element > mSplitted_Elements;
WeakPointerVector< Node<3> > mFail_Node;


};
}
#endif 

