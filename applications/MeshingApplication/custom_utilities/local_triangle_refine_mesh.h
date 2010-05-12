//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: Nelson $
//   Date:                $Date: 2010-05-12 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_LOCAL_REFINE_TRIANGLE_MESH)
#define  KRATOS_LOCAL_REFINE_TRIANGLE_MESH


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
#include "geometries/triangle_2d_3.h"
#include "processes/node_erase_process.h" 
#include "spatial_containers/spatial_containers.h"


namespace Kratos
{
    class Local_Refine_Triangle_Mesh
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

              Local_Refine_Triangle_Mesh(ModelPart& model_part) : mr_model_part(model_part) {}

              ~Local_Refine_Triangle_Mesh(){}


///************************************************************************************************          
///************************************************************************************************  
void Local_Refine_Mesh()
{
      KRATOS_TRY  

      compressed_matrix<int> Coord;
      boost::numeric::ublas::vector<int> List_New_Nodes;                         ///* the news nodes
      boost::numeric::ublas::vector<array_1d<int, 2 > > Position_Node;           ///* edges where are the news nodes
      boost::numeric::ublas::vector< array_1d<double, 3> > Coordinate_New_Node;  ///* the coordinate of the new nodes

      WeakPointerVector< Element > Splitted_Elements;
      Splitted_Elements.reserve(1000);

      ElementsArrayType& pElements         =  mr_model_part.Elements(); 
      ElementsArrayType::iterator it_begin =  pElements.ptr_begin();
      ElementsArrayType::iterator it_end   =  pElements.ptr_end(); 
        
      for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
	{ 
            unsigned int level = it->GetValue(REFINEMENT_LEVEL);
            if( level > 0 )
	       {
                  Splitted_Elements.push_back(*(it.base() ));
               }
	}

	//std::cout<< "************ ELEMENT TO BE DIVIDED *************"<<std::endl;
        //KRATOS_WATCH(Splitted_Elements) 
	
	
	CSR_Row_Matrix(mr_model_part, Coord);  
	Search_Edge_To_Be_Refined(mr_model_part, Coord); 
	Create_List_Of_New_Nodes(mr_model_part, Coord, List_New_Nodes, Position_Node);
	Calculate_Coordinate_New_Nodes(mr_model_part, Position_Node,Coordinate_New_Node);
	Insert_New_Node_In_Mesh(mr_model_part, Coordinate_New_Node, List_New_Nodes);
	Erase_Old_Element_And_Create_New_Triangle_Element(mr_model_part, Coord ); 

	KRATOS_CATCH("")
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
          //KRATOS_WATCH(index_i)
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

void Search_Edge_To_Be_Refined(ModelPart& this_model_part, compressed_matrix<int>& Coord)
{

ElementsArrayType& pElements         =  this_model_part.Elements(); 

ElementsArrayType::iterator it_begin = pElements.ptr_begin(); //+element_partition[k];
ElementsArrayType::iterator it_end   = pElements.ptr_end();     //+element_partition[k+1];
for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
  {   
        unsigned int level = it->GetValue(REFINEMENT_LEVEL);
         if( level > 0 )
           {          
  	    Element::GeometryType& geom = it->GetGeometry(); // Nodos del elemento        
 	    for (unsigned int i = 0; i <geom.size(); i++)
 	      {  int index_i = geom[i].Id()-1;
                 for (unsigned int j = 0; j<geom.size(); j++)
 	          {  int index_j = geom[j].Id()-1;
                     if(index_j > index_i)
		      {                     
                        Coord(index_i, index_j) = -2;
                      }  
                  }
             }
        }

    }
}

///************************************************************************************************          
///************************************************************************************************ 

void Create_List_Of_New_Nodes(ModelPart& this_model_part, compressed_matrix<int>& Coord, boost::numeric::ublas::vector<int> &List_New_Nodes,
boost::numeric::ublas::vector<array_1d<int, 2 > >& Position_Node)
{
 
unsigned int numbrer_of_new_nodes = 0;
//boost::numeric::ublas::vector<array_1d<int, 2 > > Position_Node;   ///* Edges to be refided   

NodesArrayType& pNodes =  this_model_part.Nodes();  

typedef compressed_matrix<int>::iterator1 i1_t;
typedef compressed_matrix<int>::iterator2 i2_t;

for (i1_t i1 = Coord.begin1(); i1 != Coord.end1(); ++i1) {
    for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
      { if( Coord( i2.index1(), i2.index2() )==-2)
           { 
             numbrer_of_new_nodes++; 
           }
       }      
 }
///* New Id de los Nodos
List_New_Nodes.resize(numbrer_of_new_nodes);
int total_node =  pNodes.size();
for(unsigned int i = 0; i <numbrer_of_new_nodes; i++ )
{
    List_New_Nodes[i] = total_node + i + 1;
}


///* setting edges -2 to the new id of the new node
Position_Node.resize(numbrer_of_new_nodes);
unsigned int index = 0;
for (i1_t i1 = Coord.begin1(); i1 != Coord.end1(); ++i1) {
    for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
      { if( Coord( i2.index1(), i2.index2() )==-2)
           { 
             Coord(i2.index1(), i2.index2()) =  List_New_Nodes[index];
             Position_Node[index][0] = i2.index1()+1; 
             Position_Node[index][1] = i2.index2()+1;
             index++;             
           }
       }      
 }


}

///************************************************************************************************          
///************************************************************************************************  

///* WARNING = Computes the coordinate of the new nodes of the center of the edges s triangle.
void  Calculate_Coordinate_New_Nodes(ModelPart& this_model_part,
boost::numeric::ublas::vector<array_1d<int, 2 > >& Position_Node,
boost::numeric::ublas::vector< array_1d<double, 3> >& Coordinate_New_Node)
{

array_1d<double, 3 > Coord_Node_1; 
array_1d<double, 3 > Coord_Node_2;
//boost::numeric::ublas::vector< array_1d<double, 3> > Coordinate_New_Node; 
Coordinate_New_Node.resize(Position_Node.size());

for(unsigned int i = 0; i < Position_Node.size(); i++ )
  {
     int& node_i     = Position_Node[i][0];
     int& node_j     = Position_Node[i][1];
     
     Coord_Node_1[0] = this_model_part.Nodes()[node_i].X();
     Coord_Node_1[1] = this_model_part.Nodes()[node_i].Y();
     Coord_Node_1[2] = this_model_part.Nodes()[node_i].Z();

     Coord_Node_2[0] = this_model_part.Nodes()[node_j].X();
     Coord_Node_2[1] = this_model_part.Nodes()[node_j].Y();
     Coord_Node_2[2] = this_model_part.Nodes()[node_j].Z();
     noalias(Coordinate_New_Node[i]) =  0.50 * (Coord_Node_1 + Coord_Node_2);
      
  }

}


///************************************************************************************************          
///************************************************************************************************  

///* inserta los nuevos nodeos en el model part
void Insert_New_Node_In_Mesh(ModelPart& this_model_part,
boost::numeric::ublas::vector< array_1d<double, 3> >& Coordinate_New_Node,
boost::numeric::ublas::vector<int>& List_New_Nodes
)

{
//node to get the DOFs from	
unsigned int size = Coordinate_New_Node.size(); 
unsigned int step_data_size = this_model_part.GetNodalSolutionStepDataSize();
Node<3>::DofsContainerType& reference_dofs = (this_model_part.NodesBegin())->GetDofs();
 
for(unsigned int i = 0; i<size; i++)
 {
    Node<3>::Pointer  pnode = this_model_part.CreateNewNode(List_New_Nodes[i],Coordinate_New_Node[i][0], Coordinate_New_Node[i][1], Coordinate_New_Node[i][2] );
    pnode->SetBufferSize(this_model_part.NodesBegin()->GetBufferSize() );
       
   //generating the dofs  
   for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin();    iii != reference_dofs.end(); iii++)  
     {
       Node<3>::DofType& rDof = *iii;
       Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );
       (p_new_dof)->FreeDof();
     }
/*
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


    }  

 }

///************************************************************************************************          
///************************************************************************************************  

void Erase_Old_Element_And_Create_New_Triangle_Element(
           ModelPart& this_model_part,
           const compressed_matrix<int>& Coord          
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
       PointerVector< Element > New_Elements;  
       

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
      
             edge_ids[0] = Coord(index_i, index_j );
             if(index_i > index_j ) {edge_ids[0] = Coord(index_j, index_i );}                 
             edge_ids[1] = Coord(index_j, index_k );
             if(index_j > index_k ) {edge_ids[1] = Coord(index_k, index_j );}                 
             edge_ids[2] = Coord(index_k, index_i);
             if(index_k > index_i ) {edge_ids[2] = Coord(index_i, index_k );}   
                      

            ///* crea las nuevas conectividades
            create_element = Split_Triangle_Elements::Split_Triangle(triangle_ids, edge_ids, new_conectivity);   
            ///* crea los nuevos elementos
                         
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
                       //pElements.push_back(p_element); 
                       //KRATOS_WATCH(p_element)
                       current_id++;    
                                         
               }
                it->SetId(large_id);   		
                large_id++; 
	     }
             //interpolate
          }
     
          ///* adding news elements to the model part
        
         for(PointerVector< Element >::iterator it_new = New_Elements.begin(); it_new!=New_Elements.end(); it_new++)
          {      
              pElements.push_back(*it_new);
          }    
             
          ///* all of the elements to be erased are at the end
	  this_model_part.Elements().Sort(); 
          ///*now remove all of the "old" elements
	  this_model_part.Elements().erase(this_model_part.Elements().end()-to_be_deleted, this_model_part.Elements().end());
					
  }




protected:
ModelPart& mr_model_part;
         
 
      };

}  // namespace Kratos.

#endif // KRATOS_PROJECTION  defined 


