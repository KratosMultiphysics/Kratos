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
#include "includes/dof.h"
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
void Local_Refine_Mesh(bool refine_on_reference)
{
      KRATOS_TRY  

      NodesArrayType& pNodes = mr_model_part.Nodes();
      compressed_matrix<int> Coord;
      boost::numeric::ublas::vector<int> List_New_Nodes;                         ///* the news nodes
      boost::numeric::ublas::vector<array_1d<int, 2 > > Position_Node;           ///* edges where are the news nodes
      boost::numeric::ublas::vector< array_1d<double, 3> > Coordinate_New_Node;  ///* the coordinate of the new nodes
     
        
      PointerVector< Element > New_Elements; 
      
      if(refine_on_reference==true)
       {
          for(ModelPart::NodesContainerType::iterator it=mr_model_part.NodesBegin(); it!=mr_model_part.NodesEnd(); it++)
              {
	       it->X()=it->X0();
	       it->Y()=it->Y0();
	       it->Z()=it->Z0();
             }
       }
 	
      CSR_Row_Matrix(mr_model_part, Coord);  
      Search_Edge_To_Be_Refined(mr_model_part, Coord); 
      Create_List_Of_New_Nodes(mr_model_part, Coord, List_New_Nodes, Position_Node);
      Calculate_Coordinate_And_Insert_New_Nodes(mr_model_part, Position_Node, List_New_Nodes);
      //Insert_New_Node_In_Mesh(mr_model_part, Coordinate_New_Node, Position_Node, List_New_Nodes);
      Erase_Old_Element_And_Create_New_Triangle_Element(mr_model_part, Coord, New_Elements); 
      Renumering_Elements_And_Nodes(mr_model_part, New_Elements);


      
  if(refine_on_reference==true)
  {
     for(ModelPart::NodesContainerType::iterator it=mr_model_part.NodesBegin(); it!=mr_model_part.NodesEnd(); it++)
      {     
	     const array_1d<double,3>& disp = it->FastGetSolutionStepValue(DISPLACEMENT);
	     it->X()=it->X0() + disp[0];
	     it->Y()=it->Y0() + disp[1];
 	     it->Z()=it->Z0() + disp[2];
          }
  } 
     
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

ElementsArrayType& rElements         =  this_model_part.Elements(); 

ElementsArrayType::iterator it_begin = rElements.ptr_begin(); //+element_partition[k];
ElementsArrayType::iterator it_end   = rElements.ptr_end();     //+element_partition[k+1];
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
 
unsigned int number_of_new_nodes = 0;

NodesArrayType& pNodes =  this_model_part.Nodes();  

typedef compressed_matrix<int>::iterator1 i1_t;
typedef compressed_matrix<int>::iterator2 i2_t;

///*WARNING
for (i1_t i1 = Coord.begin1(); i1 != Coord.end1(); ++i1) {
    for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
      { if( Coord( i2.index1(), i2.index2() )==-2)
           { 
             number_of_new_nodes++; 
           }
       }      
 }


///* New Id de los Nodos
List_New_Nodes.resize(number_of_new_nodes);
int total_node =  pNodes.size();
for(unsigned int i = 0; i <number_of_new_nodes; i++ )
{
    List_New_Nodes[i] = total_node + i + 1;
}


///* setting edges -2 to the new id of the new node
///* WARNING
Position_Node.resize(number_of_new_nodes);
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
///* insert the news nodes in the model part and interopolate the variables

void  Calculate_Coordinate_And_Insert_New_Nodes(ModelPart& this_model_part,
const boost::numeric::ublas::vector<array_1d<int, 2 > >& Position_Node,
const boost::numeric::ublas::vector<int> &List_New_Nodes)
{

array_1d<double, 3 > Coord_Node_1; 
array_1d<double, 3 > Coord_Node_2;
boost::numeric::ublas::vector< array_1d<double, 3> > Coordinate_New_Node;
Coordinate_New_Node.resize(Position_Node.size());
unsigned int step_data_size = this_model_part.GetNodalSolutionStepDataSize();
Node<3>::DofsContainerType& reference_dofs = (this_model_part.NodesBegin())->GetDofs(); 



for(unsigned int i = 0; i < Position_Node.size(); i++ )
  {

     /// calculating the coordinate of the news nodes  
     const int& node_i     = Position_Node[i][0];
     const int& node_j     = Position_Node[i][1];
     ModelPart::NodesContainerType::iterator it_node1 =  this_model_part.Nodes().find(node_i);
     std::size_t pos1 = it_node1 - this_model_part.NodesBegin();
     noalias(Coord_Node_1) = it_node1->Coordinates(); 
     ModelPart::NodesContainerType::iterator it_node2 =  this_model_part.Nodes().find(node_j);
     std::size_t pos2 = it_node2 - this_model_part.NodesBegin();
     noalias(Coord_Node_2) = it_node2->Coordinates(); 
     noalias(Coordinate_New_Node[i]) =  0.50 * (Coord_Node_1 + Coord_Node_2);
     
     /// inserting the news node in the model part 
     Node<3>::Pointer  pnode = this_model_part.CreateNewNode(List_New_Nodes[i],Coordinate_New_Node[i][0], Coordinate_New_Node[i][1], Coordinate_New_Node[i][2] );
     pnode->SetBufferSize(this_model_part.NodesBegin()->GetBufferSize() );
 
     it_node1 = this_model_part.NodesBegin() + pos1;
     it_node2 = this_model_part.NodesBegin() + pos2;
//      it_node1 =  this_model_part.Nodes().find(node_i);
//      it_node2 =  this_model_part.Nodes().find(node_j);

    
     pnode->X0() = 0.5*(it_node1->X0() + it_node2->X0());
     pnode->Y0() = 0.5*(it_node1->Y0() + it_node2->Y0());
     pnode->Z0() = 0.5*(it_node1->Z0() + it_node2->Z0()); 
     
     for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin();    iii != reference_dofs.end(); iii++)  
      {
       Node<3>::DofType& rDof = *iii;
       Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );
       //KRATOS_WATCH(it_node1->Id())
       //KRATOS_WATCH(it_node2->Id())  
       //KRATOS_WATCH(iii->GetVariable())  
       //KRATOS_WATCH(it_node1->IsFixed(iii->GetVariable()))   
       //KRATOS_WATCH(it_node2->IsFixed(iii->GetVariable()))          
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
          new_step_data[j] = 0.5*(step_data1[j] + step_data2[j]);
        }						
     }
    
      /// WARNING =  only for reactions;
      const double zero = 0.00;
      for(Node<3>::DofsContainerType::iterator iii = pnode->GetDofs().begin();    iii != pnode->GetDofs().end(); iii++)  
       {          
          if(pnode->IsFixed(iii->GetVariable())==false) 
                 { 
                    //KRATOS_WATCH("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA") 
                    //KRATOS_WATCH( iii->GetSolutionStepReactionValue());
                    iii->GetSolutionStepReactionValue() = zero;
                    //KRATOS_WATCH( iii->GetSolutionStepReactionValue());
                  }       
       }          
  }
}


/*
/// ************************************************************************************************          
/// ************************************************************************************************  

///  inserta los nuevos nodeos en el model parts
void Insert_New_Node_In_Mesh(ModelPart& this_model_part,
boost::numeric::ublas::vector< array_1d<double, 3> >& Coordinate_New_Node,
boost::numeric::ublas::vector<array_1d<int, 2 > >& Position_Node,
boost::numeric::ublas::vector<int>& List_New_Nodes
)

{
//node to get the DOFs from	
unsigned int size = Coordinate_New_Node.size(); 
unsigned int step_data_size = this_model_part.GetNodalSolutionStepDataSize();
Node<3>::DofsContainerType& reference_dofs = (this_model_part.NodesBegin())->GetDofs();


 
for(unsigned int i = 0; i<size; i++)
 {
    int& node_i     = Position_Node[i][0];
    int& node_j     = Position_Node[i][1];
     
    ModelPart::NodesContainerType::iterator it_node1 =  this_model_part.Nodes().find(node_i);
    ModelPart::NodesContainerType::iterator it_node2 =  this_model_part.Nodes().find(node_j);

    Node<3>::Pointer  pnode = this_model_part.CreateNewNode(List_New_Nodes[i],Coordinate_New_Node[i][0], Coordinate_New_Node[i][1], Coordinate_New_Node[i][2] );
    pnode->SetBufferSize(this_model_part.NodesBegin()->GetBufferSize() );
       
   /// generating the dofs  
   for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin();    iii != reference_dofs.end(); iii++)  
     {
       Node<3>::DofType& rDof = *iii;
       Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );
       //if(it_node1->IsFixed(rDof) == true && it_node2->IsFixed(rDof) == true)
          (p_new_dof)->FixDof();
       //else
        //  { (p_new_dof)->FreeDof();}       

     }

     //interpolate
     unsigned int buffer_size = pnode->GetBufferSize();
     for(unsigned int step = 0; step<buffer_size; step++)
      {
	double* new_step_data    = pnode->SolutionStepData().Data(step);
        double* step_data1     = it_node1->SolutionStepData().Data(step);			
        double* step_data2     = it_node2->SolutionStepData().Data(step);
        //copying this data in the position of the vector we are interested in
        for(unsigned int j= 0; j<step_data_size; j++)
        {  
          new_step_data[j] = 0.5*(step_data1[j] + step_data2[j]);
        }						
     }

    }  

 }

*/


/// ************************************************************************************************          
/// ************************************************************************************************  


void Erase_Old_Element_And_Create_New_Triangle_Element(
           ModelPart& this_model_part,
           const compressed_matrix<int>& Coord,
           PointerVector< Element >& New_Elements            
          )   
    {

       boost::numeric::ublas::matrix<int> new_conectivity;    
       ElementsArrayType& rElements         =  this_model_part.Elements(); 
       ElementsArrayType::iterator it_begin =  rElements.ptr_begin();
       ElementsArrayType::iterator it_end   =  rElements.ptr_end(); 
       Element const rReferenceElement;   
       unsigned int to_be_deleted=0;
       unsigned int large_id = (rElements.end()-1)->Id() * 7;
       array_1d<int,3> triangle_ids; 
       array_1d<int,3> edge_ids;
       bool create_element = false; 
       //PointerVector< Element > New_Elements;  
       
       
       ProcessInfo& rCurrentProcessInfo  = this_model_part.GetProcessInfo();  


	unsigned int current_id = (rElements.end()-1)->Id() + 1;
        for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
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
                       //const std::string type = p_element->Info();  
                       //KRATOS_WATCH(p_element->Id())
                       //KRATOS_WATCH(type)
                       
                       New_Elements.push_back(p_element);                         
                       current_id++;    
                                         
               }
                it->SetId(large_id);   
                large_id++;
	     }
               
          }         
  
          ///* adding news elements to the model part
          //Vector temp;
          for(PointerVector< Element >::iterator it_new = New_Elements.begin(); it_new!=New_Elements.end(); it_new++)
            {      
               
               it_new->Initialize();  
               it_new->InitializeSolutionStep(rCurrentProcessInfo);      
               //it_new->CalculateRightHandSide(temp, rCurrentProcessInfo);
               it_new->FinalizeSolutionStep(rCurrentProcessInfo);           
               rElements.push_back(*(it_new.base()));
            } 
             
          ///* all of the elements to be erased are at the end
	  rElements.Sort(); 
          
          ///*now remove all of the "old" elements
	  rElements.erase(this_model_part.Elements().end()-to_be_deleted, this_model_part.Elements().end());
             
}


void Renumering_Elements_And_Nodes( ModelPart& this_model_part,
PointerVector< Element >& New_Elements)
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


  ///* adding news elements to the model part
//  Vector temp;  
/*
  for(PointerVector< Element >::iterator it_new = New_Elements.begin(); it_new!=New_Elements.end(); it_new++)
    {      
      std::cout<< "************************************************" <<std::endl;
      KRATOS_WATCH(it_new->Id())  
      it_new->Initialize();  
      it_new->InitializeSolutionStep(rCurrentProcessInfo);      
      it_new->CalculateRightHandSide(temp, rCurrentProcessInfo);
      it_new->FinalizeSolutionStep(rCurrentProcessInfo);           
      std::cout<< "************************************************" <<std::endl;
    }*/ 

}



inline void CreatePartition(unsigned int number_of_threads, const int number_of_rows, vector<unsigned int>& partitions)
    {
      partitions.resize(number_of_threads+1);
      int partition_size = number_of_rows / number_of_threads;
      partitions[0] = 0;
      partitions[number_of_threads] = number_of_rows;
      for(unsigned int i = 1; i<number_of_threads; i++)
      partitions[i] = partitions[i-1] + partition_size ;
    }


protected:
ModelPart& mr_model_part;
         
 
  };

}  // namespace Kratos.

#endif // KRATOS_PROJECTION  defined 


