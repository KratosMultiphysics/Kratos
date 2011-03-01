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

#if !defined(INTER_FRACTURE_TRIANGLE_UTILITY)
#define INTER_FRACTURE_TRIANGLE_UTILITY

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
#include "processes/find_conditions_neighbours_process.h"
#include "containers/data_value_container.h"
#include "includes/mesh.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "custom_utilities/smoothing_utility.h"
#include "geometries/triangle_2d_3.h"
#include "processes/node_erase_process.h" 
#include "spatial_containers/spatial_containers.h"


namespace Kratos
{
      class Inter_Fracture_Triangle
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

 


          Inter_Fracture_Triangle(ModelPart& model_part, int domain_size) : mr_model_part(model_part)
          {
             mdomain_size = domain_size;
          }
          
          ~Inter_Fracture_Triangle(){}
    

///************************************************************************************************          
void Detect_And_Split_Elements(ModelPart& this_model_part)
{ 
   KRATOS_TRY      
   
   //NodesArrayType& pNodes = this_model_part.Nodes();  
   array_1d<double,3> Failure_Maps; 
   
   FindElementalNeighboursProcess    ElementosVecinos(this_model_part, 2, 10);
   FindNodalNeighboursProcess        NodosVecinos(this_model_part, 2, 10);
   FindConditionsNeighboursProcess   CondicionesVecinas(this_model_part, 2, 10);
   
   WeakPointerVector< Node<3> > Nodes_To_Be_Dupplicated; 
   unsigned int detect = Detect_Node_To_Be_Splitted(this_model_part, Nodes_To_Be_Dupplicated);
   KRATOS_WATCH(detect)
   
   
   if(detect!=0)
   {   
   WeakPointerVector< Node<3> >::iterator i_begin = Nodes_To_Be_Dupplicated.ptr_begin();
   WeakPointerVector< Node<3> >::iterator i_end   = Nodes_To_Be_Dupplicated.ptr_end();
      
   for(WeakPointerVector< Node<3> >::iterator inode=i_begin; inode!= i_end; ++inode)     
    {      
             Node<3>::Pointer pNode =  (*(inode.base())).lock();              
             Split_Node(this_model_part, pNode);       
	    
	     //ElementosVecinos.ClearNeighbours();
             //NodosVecinos.ClearNeighbours();
	     //CondicionesVecinas.ClearNeighbours();     
             //ElementosVecinos.Execute();
             //NodosVecinos.Execute();
	     //CondicionesVecinas.Execute();         
          }                
    
    }
    
         
  Finalize(this_model_part);
  KRATOS_CATCH("")
}


///************************************************************************************************          
///************************************************************************************************      
unsigned int Detect_Node_To_Be_Splitted(ModelPart& this_model_part,  WeakPointerVector< Node<3> >& Nodes_To_Be_Dupplicated)
{
KRATOS_TRY  

NodesArrayType& pNodes = this_model_part.Nodes();
NodesArrayType::iterator i_begin = pNodes.ptr_begin();
NodesArrayType::iterator i_end   = pNodes.ptr_end();

for(ModelPart::NodeIterator inode=i_begin; inode!= i_end; ++inode)  
   { 
      double& Condition = inode->GetValue(NODAL_DAMAGE);
       if(Condition >=0.01){
	       Nodes_To_Be_Dupplicated.push_back(*(inode.base()));     
            }
    }
    
return Nodes_To_Be_Dupplicated.size();
  
KRATOS_CATCH("")

}

///************************************************************************************************          
///************************************************************************************************      
///* Computa el vector, linea, o plano de falla.
void Calculate_Map_Failure(Node<3>::Pointer& pNode,  array_1d<double,3>& Failure_Maps)
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
  
}


///************************************************************************************************          
///************************************************************************************************  

                 
void  Split_Node(ModelPart& this_model_part, Node<3>::Pointer& pNode)
{
    KRATOS_TRY 
    
    array_1d<double,3>  failure_map;
    Calculate_Map_Failure(pNode, failure_map);
    
    WeakPointerVector< Element > Negative_Elements;
    WeakPointerVector< Element > Positive_Elements;
    
    
    Node<3>::Pointer child_node; 
    bool node_created = false;
    
    /// Crea el nuevo nodo y lo inserta al elemento
    //std::cout<<"INSERTING THE NEW NODE"<<std::endl;
    node_created = CalculateElements(this_model_part, pNode, child_node, failure_map, Negative_Elements, Positive_Elements );
    
    
    if (node_created==true)
    {
        ///UPDATING CONDITIONS
	/// En caso de que el nodo creado modifique la condicion de contorno
        CalculateConditions(this_model_part, pNode, child_node, failure_map );
   
        /// Updating vecinos    
        RecomputeLocalneighbourgs(pNode,      Positive_Elements);
        RecomputeLocalneighbourgs(child_node, Negative_Elements);
	
    }
    KRATOS_CATCH("")

}


void RecomputeLocalneighbourgs(Node<3>::Pointer& pNode, WeakPointerVector< Element>& Elements)
{
      WeakPointerVector< Element >& neighb_elems  = pNode->GetValue(NEIGHBOUR_ELEMENTS); 
      WeakPointerVector< Node<3> >& neighb_nodes  = pNode->GetValue(NEIGHBOUR_NODES); 
      neighb_elems.clear();
      neighb_nodes.clear();
      
      for(WeakPointerVector<Element>::iterator elem = Elements.begin(); elem!= Elements.end(); ++elem){
	    neighb_elems.push_back(*elem.base());
	    Element::GeometryType& geom = elem->GetGeometry(); 
	        for( unsigned int i = 0; i < geom.size(); i++ )
                   { 
		     WeakPointerVector< Node<3> >::iterator repeated_object =  std::find(neighb_nodes.begin(), neighb_nodes.end(), geom[i]);
		      if(repeated_object == (neighb_nodes.end())) 
 			  {
			    neighb_nodes.push_back(geom(i));
			  }
			  
		   }    
      }

}


///************************************************************************************************          
///************************************************************************************************

//pNode = the father node
bool CalculateElements(ModelPart& this_model_part, 
		       Node<3>::Pointer& pNode, // Nodo Padre 
		       Node<3>::Pointer& pnode, // Nodo Hijo
		       const array_1d<double,3>&  failure_map,
		       WeakPointerVector< Element >& Negative_Elements, 
		       WeakPointerVector< Element >& Positive_Elements 
		       )
{
  KRATOS_TRY
  
  
    double prod  = 0.00;
    array_1d<double, 3> Coord_Point_1; 
    array_1d<double, 3> Coord_Point_2;
    array_1d<double, 3> normal;
    array_1d<double, 3> Unit; 

      
    Positive_Elements.reserve(10); 
    Negative_Elements.reserve(10); 
     
    
    WeakPointerVector< Element >& neighb_elems  = pNode->GetValue(NEIGHBOUR_ELEMENTS);

    // Normal al plano de fractura
    Calculate_Normal_Faliure_Maps(normal,failure_map);

    noalias(Coord_Point_1)  = pNode->Coordinates(); 
    //Coord_Point_1[1]  = pNode->Y();   
    //Coord_Point_1[2]  = pNode->Z();
    
   
    for(WeakPointerVector< Element >::iterator neighb_elem = neighb_elems.begin();
    neighb_elem != neighb_elems.end(); neighb_elem++)
       { 
           Element::GeometryType& geom = neighb_elem->GetGeometry(); // Nodos del elemento 
	   noalias(Coord_Point_2) = geom.Center(); 
           //Find_Coord_Gauss_Points(geom, Coord_Point_2);
           noalias(Unit) = Coord_Point_2 - Coord_Point_1;
           noalias(Unit) = Unit / norm_2(Unit);
	   
           prod = inner_prod(normal, Unit);
           if( prod>=0.00)
              {
               //std::cout<<"INSERTED_POSITIVE"<<std::endl;
               Positive_Elements.push_back(*(neighb_elem.base() ));
              }
           else
             {
               //std::cout<<"INSERTED_NEGATIVE"<std::endl;
               Negative_Elements.push_back(*(neighb_elem.base() ));
             }    

            Unit   = ZeroVector(3);
            prod   = 0.00;  
        }

    ///* Esquinas o superficie extrena donde no hay elementos negativos o positivos
    if (Positive_Elements.size()==0  || Negative_Elements.size()==0)
        {
          Positive_Elements.clear(); 
          Negative_Elements.clear(); 
          for(WeakPointerVector< Element >::iterator neighb_elem = neighb_elems.begin();
          neighb_elem != neighb_elems.end(); neighb_elem++)
            { 
               if( inner_prod(failure_map, Unit)>=0.00){ Positive_Elements.push_back(*(neighb_elem.base() ) );}
               else { Negative_Elements.push_back(*(neighb_elem.base() )); }    
            }
         }


    ///* Si el nodo no mas tiene un solo elemtno vecino
    bool& duplicated_pNode = pNode->GetValue(IS_DUPLICATED);
    if( (Positive_Elements.size()==0 && Negative_Elements.size()==0)  ||  (Positive_Elements.size()==1 && Negative_Elements.size()==0) || (Positive_Elements.size()==0 && Negative_Elements.size()==1) || duplicated_pNode==true)
    {
     //std::cout<<"NO INSERTED NODE"<<std::endl;
        return false; 
    }
    else
    {

	///* reseting internal variables
	for(WeakPointerVector< Element >::iterator neighb_elem = neighb_elems.begin();
	neighb_elem != neighb_elems.end(); neighb_elem++)
	  {
	    neighb_elem->Initialize();
	  }

	// creating new nodes
	//pNode->FastGetSolutionStepValue(IS_DUPLICATED)=true;
	unsigned int New_Id = this_model_part.Nodes().size() + 1; 
	//Node<3>::Pointer pnode; // the new node   
	bool& duplicated_pnode = pNode->GetValue(IS_DUPLICATED);
	duplicated_pnode = true;
	Create_New_Node(this_model_part, pnode, New_Id, pNode);

	///* putting de new node to the negative element  
	for(WeakPointerVector< Element >::iterator neg_elem = Negative_Elements.begin();
	neg_elem != Negative_Elements.end(); neg_elem++)
	{
	    Element::GeometryType& geom = neg_elem->GetGeometry();
	    Insert_New_Node_In_Elements(geom, pnode,  pNode->Id());   
	}
	return true;
    }
    
  KRATOS_CATCH("")

}


///************************************************************************************************
///************************************************************************************************

void Create_New_Node(ModelPart& this_model_part, Node<3>::Pointer& pnode, unsigned int& New_Id, Node<3>::Pointer& pNode)
{
  
  
//node to get the DOFs from	
int step_data_size         = this_model_part.GetNodalSolutionStepDataSize();
array_1d<double, 3>& Coord = pNode->Coordinates();				
Node<3>::DofsContainerType& reference_dofs = (this_model_part.NodesBegin())->GetDofs();

pnode = this_model_part.CreateNewNode(New_Id,Coord[0],Coord[1],Coord[2]);
pnode->SetBufferSize(this_model_part.NodesBegin()->GetBufferSize() );

//generating the dofs
for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin();    iii != reference_dofs.end(); iii++)
{
    Node<3>::DofType& rDof = *iii;
    Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );
    if(pNode->IsFixed(iii->GetVariable()) == true )
      (p_new_dof)->FixDof();
    else
      { (p_new_dof)->FreeDof();}    
}


unsigned int buffer_size = pnode->GetBufferSize();
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

//pNode = the father node
void CalculateConditions(ModelPart& this_model_part, 
Node<3>::Pointer& pNode,  // father node
Node<3>::Pointer& pnode,  // child node			 
const array_1d<double,3>&  failure_map )
{
  KRATOS_TRY
  
    double prod  = 0.00;
    array_1d<double,3>  Coord_Point_1; 
    array_1d<double,3>  Coord_Point_2;  
    array_1d<double,3>  normal;
    array_1d<double,3>  Unit;
  
    WeakPointerVector< Condition > Negative_Conditions;
    Negative_Conditions.reserve(5); 
    
    WeakPointerVector< Condition >& neighb_conds  = pNode->GetValue(NEIGHBOUR_CONDITIONS);

    // Normal al plano de fractura
    Calculate_Normal_Faliure_Maps(normal,failure_map);

    Coord_Point_1  = pNode->Coordinates(); 
   

    
    for(WeakPointerVector< Condition >::iterator neighb_cond = neighb_conds.begin();
    neighb_cond != neighb_conds.end(); neighb_cond++)
       { 
           Condition::GeometryType& geom = neighb_cond->GetGeometry(); // Nodos de las condiciones 
           noalias(Coord_Point_2) = geom.Center(); 
           noalias(Unit) = Coord_Point_2 - Coord_Point_1;
           noalias(Unit) = Unit / norm_2(Unit);
           prod = inner_prod(normal, Unit);
           if( prod<0.00)
             {
               Negative_Conditions.push_back(*(neighb_cond.base() ));
             }    
            Unit   = ZeroVector(3);
            prod   = 0.00;  
	   
       }
        
    /// putting de new node to the negative conditions  
    for(WeakPointerVector< Condition >::iterator neg_cond = Negative_Conditions.begin();
    neg_cond != Negative_Conditions.end(); neg_cond++)
    {
             
         Condition::GeometryType& geom = neg_cond->GetGeometry();
         Insert_New_Node_In_Conditions(geom, pnode,  pNode->Id());
    }
    
  KRATOS_CATCH("")

}


///************************************************************************************************          
///************************************************************************************************  

void Insert_New_Node_In_Elements(Element::GeometryType& geom, Node<3>::Pointer& pnode, unsigned int Node_Id_Old)
{
       for (unsigned int i = 0; i <geom.size(); i++)
 	{         
            if (geom[i].Id()==Node_Id_Old)
             { 
               geom(i) = pnode; 
             }  
        }
}

///************************************************************************************************          
///************************************************************************************************  


void Insert_New_Node_In_Conditions(Condition::GeometryType& geom, Node<3>::Pointer& pnode, unsigned int Node_Id_Old)
{
       for (unsigned int i = 0; i <geom.size(); i++)
 	{         
            if (geom[i].Id()==Node_Id_Old)
             { 
               geom(i) = pnode; 
             }
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



void Finalize(ModelPart& this_model_part)
{
NodesArrayType& pNodes =  this_model_part.Nodes();  
#ifdef _OPENMP
int number_of_threads = omp_get_max_threads();
#else
int number_of_threads = 1;
#endif

vector<unsigned int> node_partition;
CreatePartition(number_of_threads, pNodes.size(), node_partition);
//mfail_node.clear();

#pragma omp parallel for 
for(int k=0; k<number_of_threads; k++)
 {
  NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
  NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

   for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)     
    {   
       //i->FastGetSolutionStepValue(NODAL_DAMAGE) = i->GetValue(NODAL_DAMAGE) ;
       //i->FastGetSolutionStepValue(NODAL_AREA)   = i->GetValue(NODAL_AREA) ;
       i->GetValue(SPLIT_NODAL)  = false; 
    }
  }
}




static void Calculate_Normal_Faliure_Maps(array_1d<double,3>&  normal, const array_1d<double,3>&  failure_map)
{
 
 // WARNING: SOLO 2D
 //double tetha = atan(failure_map[1]/failure_map[0]);
 //tetha += PI/2.00;
 normal[0] = -failure_map[1]; //cos(tetha);     
 normal[1] =  failure_map[0];
 normal[2] =   0.00;
 noalias(normal) = normal / norm_2(normal);

}


ModelPart& mr_model_part;
unsigned int mdomain_size;
//WeakPointerVector< Node<3> > mfail_node;


};
}
#endif 

