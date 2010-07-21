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

#if !defined(INTER_FRACTURE_TETRAHEDRA_UTILITY)
#define INTER_FRACTURE_TETRAHEDRA_UTILITY

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
#include "geometries/tetrahedra_3d_4.h"
#include "processes/node_erase_process.h" 
#include "spatial_containers/spatial_containers.h"


namespace Kratos
{
      class Inter_Fracture_Tetrahedra
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

 


          Inter_Fracture_Tetrahedra(ModelPart& model_part, int domain_size) : mr_model_part(model_part)
          {
             mdomain_size = domain_size;   
             //mInitialize  = false;  
          }
          
          ~Inter_Fracture_Tetrahedra(){}
    

///************************************************************************************************          
///************************************************************************************************ 
///* Detecta los dos elementos en el pano que seran dividos.
///* En total son tres nodos que seran insertados uno de ellos se duplicara.
void Detect_And_Split_Elements(ModelPart& this_model_part)
{ 
   KRATOS_TRY      
   NodesArrayType& pNodes = this_model_part.Nodes();  

   Detect_Node_To_Be_Splitted(this_model_part);
   NodesArrayType::iterator i_begin=pNodes.ptr_begin();
   NodesArrayType::iterator i_end=pNodes.ptr_end();  
   FindElementalNeighboursProcess ElementosVecinos(this_model_part, 2, 10);
   FindNodalNeighboursProcess     NodosVecinos(this_model_part, 2, 10);

   //ElementosVecinos.Execute(); 
   //NodosVecinos.Execute(); 

    for(ModelPart::NodeIterator inode=i_begin; inode!= i_end; ++inode)     
      {      
            bool& split = inode->FastGetSolutionStepValue(SPLIT_NODAL);   
            if(split == true) 
             {
              
             Node<3>::Pointer pNode = *(inode.base());   
             array_1d<double,3> Failure_Maps;  
             Calculate_Map_Failure(pNode,  Failure_Maps); 
             //KRATOS_WATCH(inode->Id()) 
             Split_Node(this_model_part, pNode, Failure_Maps);   
              
             ElementosVecinos.ClearNeighbours();
             NodosVecinos.ClearNeighbours();
             ElementosVecinos.Execute(); 
             NodosVecinos.Execute();  

            }                
    }

  
   Finalize(this_model_part); 

  KRATOS_CATCH("")
}

///************************************************************************************************          
///************************************************************************************************      
void Detect_Node_To_Be_Splitted(ModelPart& this_model_part)
{
KRATOS_TRY  

NodesArrayType& pNodes = this_model_part.Nodes();

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
    double& Condition = i->FastGetSolutionStepValue(NODAL_DAMAGE);
    if(Condition >= 1.00)
	  {  
	      i->FastGetSolutionStepValue(SPLIT_NODAL) = true;   
              //mfail_node.push_back(*i.base()); 
	  }
  }

}
  
KRATOS_CATCH("")
}

///************************************************************************************************          
///************************************************************************************************      
///* Computa el vector, linea, o plano de falla.
void Calculate_Map_Failure(Node<3>::Pointer pNode,  array_1d<double,3>& Failure_Maps)
{

Vector Eigen_Values  = ZeroVector(mdomain_size);               // Deformaciones o Tensiones principales
Matrix Eigen_Vector  = ZeroMatrix(mdomain_size,mdomain_size);  // Direcciones  principales  
Matrix Strain_Tensor = ZeroMatrix(mdomain_size,mdomain_size);
 
unsigned int size = 3;
unsigned int max_iterations    = 100;
double zero_tolerance = 1e-9;
if(mdomain_size==3){size = 6; }
Matrix Nodal_Values  = ZeroMatrix(1,size);

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
int pos    =  0;

for(unsigned int k=0; k<Eigen_Values.size(); k++ )     
{ if(max==Eigen_Values(k)) { pos = k; }
  for(unsigned int j=0; j<Eigen_Values.size(); j++ ) 
 {
  EigenVector[k][j] = Eigen_Vector(k,j);
 }
}

Failure_Maps[0] =  EigenVector[pos][0];
Failure_Maps[1] =  EigenVector[pos][1];
Failure_Maps[2] =  EigenVector[pos][2];
  
}


///************************************************************************************************          
///************************************************************************************************  

                 
void  Split_Node(ModelPart& this_model_part, Node<3>::Pointer pNode, const array_1d<double,3>&  failure_map)
{
    KRATOS_TRY 


    WeakPointerVector< Element > Negative_Elements;
    WeakPointerVector< Element > Positive_Elements;

    Positive_Elements.reserve(10); 
    Negative_Elements.reserve(10); 

    //WeakPointerVector< Node<3> >& neighb_nodes  = pNode->GetValue(NEIGHBOUR_NODES);
    WeakPointerVector< Element >& neighb_elems  = pNode->GetValue(NEIGHBOUR_ELEMENTS);

    //unsigned int Node_old = 0;
    //unsigned int Node_new = 0;
    array_1d<double, 3> Coord_Point_1; 
    array_1d<double, 3> Coord_Point_2; 
 
    Coord_Point_1[0]  = pNode->X(); 
    Coord_Point_1[1]  = pNode->Y();   
    Coord_Point_1[2]  = pNode->Z();

    array_1d<double, 3> Unit;   
    //unsigned int i = 0; 
    double prod  = 0.00;
    for(WeakPointerVector< Element >::iterator neighb_elem = neighb_elems.begin();
    neighb_elem != neighb_elems.end(); neighb_elem++)
    { 
      Element::GeometryType& geom = neighb_elem->GetGeometry();  
      Find_Coord_Gauss_Points(geom, Coord_Point_2);
      prod = Calculate(failure_map, Coord_Point_1, Coord_Point_2); 

      if( prod>=0.00)
	{
	    //KRATOS_WATCH("INSERTED_POSITIVE")
	    Positive_Elements.push_back(*(neighb_elem.base() ));
	}
      else
	{
	    //KRATOS_WATCH("INSERTED_NEGATIVE")
	    Negative_Elements.push_back(*(neighb_elem.base() ));
	}    

      Unit   = ZeroVector(3);
      prod   = 0.00;  
    
    }

    /*  
    ///Esquinas o superficie extrena donde no hay elementos negativos o positivos
    if (Positive_Elements.size()==0 or Negative_Elements.size()==0)
    {
	    KRATOS_WATCH("NOT_VALIOD") 
	    Positive_Elements.clear(); 
	    Negative_Elements.clear(); 

	    for(WeakPointerVector< Element >::iterator neighb_elem = neighb_elems.begin();
	    neighb_elem != neighb_elems.end(); neighb_elem++)
	    { 
		    if( inner_prod(failure_map, Unit)>=0.00)
		    {
		      Positive_Elements.push_back(*(neighb_elem.base() ) );
		    }
		    else
		    {
		      Negative_Elements.push_back(*(neighb_elem.base() ));
		    }    
	    }
    }

    */
    ///* Si el nodo no mas tiene un solo elemtno vecino
    bool& duplicated_pNode = pNode->FastGetSolutionStepValue(IS_DUPLICATED);
    if( (Positive_Elements.size()==0 && Negative_Elements.size()==0)  ||  (Positive_Elements.size()==1 && Negative_Elements.size()==0) || (Positive_Elements.size()==0 && Negative_Elements.size()==1) || duplicated_pNode==true)
    {
      std::cout<<"NO INSERTED NODE NO INSERTED NODE NO INSERTED NODE "<<std::endl;
    }
    else
    {

    ///* reseting internal varriables

    for(WeakPointerVector< Element >::iterator neighb_elem = neighb_elems.begin();
    neighb_elem != neighb_elems.end(); neighb_elem++)
    {
      neighb_elem->Initialize();
      }

    // creating new nodes
    duplicated_pNode=true;
    pNode->FastGetSolutionStepValue(IS_DUPLICATED)=true;
    unsigned int New_Id = this_model_part.Nodes().size() + 1; 
    Node<3>::Pointer pnode; // the new node   
    bool& duplicated_pnode = pNode->FastGetSolutionStepValue(IS_DUPLICATED);
    duplicated_pnode = true;


    Create_New_Node(this_model_part, pnode, New_Id, pNode);
    //std::cout<<"New_Node_Id = " << New_Id <<std::endl;

    ///* putting de new node to the negative element  
    for(WeakPointerVector< Element >::iterator neg_elem = Negative_Elements.begin();
    neg_elem != Negative_Elements.end(); neg_elem++)
    {
      //KRATOS_WATCH(neg_elem->Id())
      Element::GeometryType& geom = neg_elem->GetGeometry();
      Insert_New_Node_Id(geom, pnode,  pNode->Id());       
      //std::cout<<"---------------------- "<<std::endl;   
    }

    // Setting los vecinos
    pnode->GetValue(NEIGHBOUR_NODES).clear();
    pnode->GetValue(NEIGHBOUR_ELEMENTS).clear();
    pNode->GetValue(NEIGHBOUR_NODES).clear();
    pNode->GetValue(NEIGHBOUR_ELEMENTS).clear();
    //KRATOS_WATCH(*pnode)  
    }
    //KRATOS_WATCH(failure_map)
    KRATOS_CATCH("")

}


///************************************************************************************************
///************************************************************************************************

void Create_New_Node(ModelPart& this_model_part, Node<3>::Pointer& pnode, unsigned int New_Id, Node<3>::Pointer& pNode)
{
//node to get the DOFs from	
int step_data_size = this_model_part.GetNodalSolutionStepDataSize();
array_1d<double, 3>& Coord = pNode->Coordinates();				
Node<3>::DofsContainerType& reference_dofs = (this_model_part.NodesBegin())->GetDofs();

pnode = this_model_part.CreateNewNode(New_Id,Coord[0] ,Coord[1],Coord[2]);

pnode->SetBufferSize(this_model_part.NodesBegin()->GetBufferSize() );

//generating the dofs
for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin();    iii != reference_dofs.end(); iii++)
{
    Node<3>::DofType& rDof = *iii;
    Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );
    (p_new_dof)->FreeDof();
}


unsigned int buffer_size = pnode->GetBufferSize();
//KRATOS_WATCH(buffer_size)
//KRATOS_WATCH(step_data_size)

for(unsigned int step = 0; step<buffer_size; step++)
{	
double* step_data     = pnode->SolutionStepData().Data(step);			
double* old_node_data = pNode->SolutionStepData().Data(step);
//copying this data in the position of the vector we are interested in
for(signed int j= 0; j<step_data_size; j++)
{  
  step_data[j] = old_node_data[j];
}						
}


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

}


///************************************************************************************************          
///************************************************************************************************  

void Insert_New_Node_Id(Element::GeometryType& geom, Node<3>::Pointer pnode, unsigned int Node_Id_Old)
{
       for (unsigned int i = 0; i <geom.size(); i++)
 	{         
            if (geom[i].Id()==Node_Id_Old)
             { 
               geom(i) = pnode; 
             }
          //KRATOS_WATCH(geom[i].Id())      
        }
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

static void Find_Coord_Gauss_Points(Element::GeometryType& geom, array_1d<double,3>&  Coord_Point)
{
        double x    = 0.00;   
        double y    = 0.00;
        double z    = 0.00;
        double fact = 0.25;
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



void Finalize(ModelPart& this_model_part)
{
NodesArrayType& pNodes =  this_model_part.Nodes();  
#ifdef _OPENMP
int number_of_threads = omp_get_max_threads();
#else
int number_of_threads = 1;
#endif

//mfail_node.clear();
vector<unsigned int> node_partition;
CreatePartition(number_of_threads, pNodes.size(), node_partition);

#pragma omp parallel for 
for(int k=0; k<number_of_threads; k++)
 {
  NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
  NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

   for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)     
    {   
       //KRATOS_WATCH(i->Id()) 
       i->FastGetSolutionStepValue(SPLIT_NODAL) = false; 
    }
  }
}




double Calculate(const array_1d<double,3>&  normal_plane, 
const array_1d<double,3>&  Coord_Po,
const array_1d<double,3>&  Coord_P)
{

double d =   -( Coord_Po[0] * normal_plane[0]  +  Coord_Po[1]*normal_plane[1] + Coord_Po[2] * normal_plane[2]);
double fact = ( Coord_P[0] * normal_plane[0]  +  Coord_P[1]*normal_plane[1] + Coord_P[2] * normal_plane[2]);
return fact + d;
}




ModelPart& mr_model_part;
unsigned int mdomain_size;
//WeakPointerVector< Node<3> > mfail_node;



};
}
#endif 

