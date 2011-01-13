
///WARNING
/// it is needed compute he neirgbourg

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
#include "includes/constitutive_law.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"

// #include "containers/array_1d.h"
// #include "processes/find_nodal_neighbours_process.h"
// #include "processes/find_elements_neighbours_process.h"
#include "containers/data_value_container.h"
#include "includes/mesh.h"
#include "utilities/math_utils.h"
//#include "utilities/split_triangle.h"
#include "utilities/split_triangle.c"
#include "geometries/triangle_2d_3.h"
#include "processes/node_erase_process.h" 
// #include "spatial_containers/spatial_containers.h"


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
void Local_Refine_Mesh(bool refine_on_reference, bool interpolate_internal_variables)
{
      KRATOS_TRY  
      
      if(refine_on_reference==true)
	 if(!(mr_model_part.NodesBegin()->SolutionStepsDataHas(DISPLACEMENT)) )
	     KRATOS_ERROR(std::logic_error,"DISPLACEMENT Variable is not in the model part -- needed if refine_on_reference = true","")

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
      Erase_Old_Element_And_Create_New_Triangle_Element(mr_model_part, Coord, New_Elements, interpolate_internal_variables);    
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
        
         //unsigned int level = it->GetValue(REFINEMENT_LEVEL);
         if (it->GetValue(SPLIT_ELEMENT) == true)
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

///*WARNING = it is can be optimizase
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
///* WARNING = it is can be optimizase
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

     pnode->GetValue(FATHER_NODES).resize(0);
     pnode->GetValue(FATHER_NODES).push_back( Node<3>::WeakPointer( *it_node1.base() ) );
     pnode->GetValue(FATHER_NODES).push_back( Node<3>::WeakPointer( *it_node2.base() ) );

     pnode->X0() = 0.5*(it_node1->X0() + it_node2->X0());
     pnode->Y0() = 0.5*(it_node1->Y0() + it_node2->Y0());
     pnode->Z0() = 0.5*(it_node1->Z0() + it_node2->Z0());

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
          new_step_data[j] = 0.5*(step_data1[j] + step_data2[j]);
        }						
     }

      /// WARNING =  only for reactions;
      const double zero = 0.00;
      for(Node<3>::DofsContainerType::iterator iii = pnode->GetDofs().begin();    iii != pnode->GetDofs().end(); iii++)  
       {
          if(pnode->IsFixed(iii->GetVariable())==false) 
                 {
                       KRATOS_TRY
                    iii->GetSolutionStepReactionValue() = zero;
                             KRATOS_CATCH("")
                  }       
       }
  }
}



/// ************************************************************************************************          
/// ************************************************************************************************  

/// WARNING = Computes the coordinate of the baricenter node of the element
/// insert the news nodes in the center of elements and interopolate the variables.
/// Create the element with the center node. No used yet.
void  Calculate_Coordinate_Center_Node_And_Insert_New_Nodes(ModelPart& this_model_part)	
{  
array_1d<double, 3 > Coord_Node_1; 
array_1d<double, 3 > Coord_Node_2;
array_1d<double, 3 > Coord_Node_3;
array_1d<double, 3 > Coordinate_center_node;
std::vector<int> node_center;
NodesArrayType& pNodes =  this_model_part.Nodes();  
int Id_Center = pNodes.size() + 1;
ElementsArrayType& rElements         =  this_model_part.Elements(); 
ElementsArrayType::iterator it_begin = rElements.ptr_begin(); //+element_partition[k];
ElementsArrayType::iterator it_end   = rElements.ptr_end();     //+element_partition[k+1];
for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
  {    
    Element::GeometryType& geom = it->GetGeometry();
    noalias(Coord_Node_1) = geom[0].Coordinates();
    noalias(Coord_Node_2) = geom[1].Coordinates();
    noalias(Coord_Node_3) = geom[2].Coordinates();

    unsigned int step_data_size = this_model_part.GetNodalSolutionStepDataSize();
    Node<3>::DofsContainerType& reference_dofs = (this_model_part.NodesBegin())->GetDofs(); 
    noalias(Coordinate_center_node) =  0.33333333333333333 * (Coord_Node_1 + Coord_Node_2 + Coord_Node_3);

    /// inserting the news node in the model part 
    Node<3>::Pointer  pnode = this_model_part.CreateNewNode(Id_Center,Coordinate_center_node[0], Coordinate_center_node[1], Coordinate_center_node[2] );
    pnode->SetBufferSize(this_model_part.NodesBegin()->GetBufferSize() );

    pnode->X0() = 0.33333333333333333*(geom[0].X0() + geom[1].X0() + geom[2].X0());
    pnode->Y0() = 0.33333333333333333*(geom[0].Y0() + geom[1].Y0() + geom[2].Y0());
    pnode->Z0() = 0.33333333333333333*(geom[0].Z0() + geom[1].Z0() + geom[2].Z0());

    for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin();    iii != reference_dofs.end(); iii++)  
	{
	  Node<3>::DofType& rDof = *iii;
	  Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );    
	  if(geom[0].IsFixed(iii->GetVariable()) == true && geom[1].IsFixed(iii->GetVariable()) == true &&  geom[2].IsFixed(iii->GetVariable()))
	    (p_new_dof)->FixDof();
	  else
	    { (p_new_dof)->FreeDof();}       
	  
	}

	///* intepolating the data
	unsigned int buffer_size = pnode->GetBufferSize();
	for(unsigned int step = 0; step<buffer_size; step++)
	{
	  double* new_step_data    = pnode->SolutionStepData().Data(step);
	  double* step_data1       = geom[0].SolutionStepData().Data(step);			
	  double* step_data2       = geom[1].SolutionStepData().Data(step); 
	  double* step_data3       = geom[2].SolutionStepData().Data(step); 
	  ///*copying this data in the position of the vector we are interested in
	  for(unsigned int j= 0; j<step_data_size; j++)
	  {  
	    new_step_data[j] = 0.333333333333333333*(step_data1[j] + step_data2[j] + step_data3[j]);
	  }						
	}
      
	/// WARNING =  only for reactions;
	const double zero = 0.00;
	for(Node<3>::DofsContainerType::iterator iii = pnode->GetDofs().begin();    iii != pnode->GetDofs().end(); iii++)  
	  {          
	    if(pnode->IsFixed(iii->GetVariable())==false) 
		    { 
		      iii->GetSolutionStepReactionValue() = zero;
		    }       
	  }
	  node_center.push_back(Id_Center);
	  Id_Center++;
  }
       
       
}
  

/// ************************************************************************************************          
/// ************************************************************************************************ 

void Erase_Old_Element_And_Create_New_Triangle_Element(
           ModelPart& this_model_part,
           const compressed_matrix<int>& Coord,
           PointerVector< Element >& New_Elements,
	   bool interpolate_internal_variables
          )   
    {

       boost::numeric::ublas::matrix<int> new_conectivity;    
       ElementsArrayType& rElements         =  this_model_part.Elements(); 
       ElementsArrayType::iterator it_begin =  rElements.ptr_begin();
       ElementsArrayType::iterator it_end   =  rElements.ptr_end(); 
       Element const rReferenceElement;   
       unsigned int to_be_deleted=0;
       unsigned int large_id = (rElements.end()-1)->Id() * 10;
       bool create_element = false; 
       int  edge_ids[3];       
       int  t[12];       
       int  nel             = 0;
       int  splitted_edges  = 0; 
       int  nint            = 0;
       array_1d<int,6> aux;
       
       
       ProcessInfo& rCurrentProcessInfo  = this_model_part.GetProcessInfo();  
       PointerVector< Element > Old_Elements;

	unsigned int current_id = (rElements.end()-1)->Id() + 1;
        for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
         { 
	     for (unsigned int i = 0; i<12; i++) {t[i]=-1; } 
  	     Element::GeometryType& geom = it->GetGeometry();
	     Calculate_Edges(geom, Coord, edge_ids, aux);		  

            ///* crea las nuevas conectividades
	    create_element =  Split_Triangle(edge_ids, t, &nel, &splitted_edges, &nint); 
	     
 
            ///* crea los nuevos elementos           
            if(create_element==true)
             {
              to_be_deleted++;  
	      for(int i=0; i<nel; i++)
	         {
                  
		  unsigned int base = i * 3;
		  unsigned int i0   = aux[t[base]];
		  unsigned int i1   = aux[t[base+1]];
		  unsigned int i2   = aux[t[base+2]];

		  	
   		  Triangle2D3<Node<3> > geom(
						this_model_part.Nodes()(i0),
						this_model_part.Nodes()(i1),
						this_model_part.Nodes()(i2) 
						);

                       
                       Element::Pointer p_element;  
 		       p_element = it->Create(current_id, geom, it->pGetProperties());
		       p_element->Initialize();  
                       p_element->InitializeSolutionStep(rCurrentProcessInfo);      
                       p_element->FinalizeSolutionStep(rCurrentProcessInfo);   
		       
		       /// setting the internal variables in the child elem
		       if(interpolate_internal_variables == true)
			  InterpolateInteralVariables(nel , *it.base(),p_element,rCurrentProcessInfo);

                       // Transfer elemental variables
                       p_element->Data() = it->Data();
		       //const unsigned int& level = it->GetValue(REFINEMENT_LEVEL);
                       p_element->GetValue(SPLIT_ELEMENT) = false;
		       //p_element->SetValue(REFINEMENT_LEVEL, 1);	       
                       New_Elements.push_back(p_element);                         
                       current_id++;    
                                         
               }
                it->SetId(large_id);   
                large_id++;
	     }
               
          }         
  
          ///* adding news elements to the model part
          for(PointerVector< Element >::iterator it_new = New_Elements.begin(); it_new!=New_Elements.end(); it_new++)
            {             
               rElements.push_back(*(it_new.base()));
            } 
             
          ///* all of the elements to be erased are at the end
	  rElements.Sort(); 
          
          ///*now remove all of the "old" elements
	  rElements.erase(this_model_part.Elements().end()-to_be_deleted, this_model_part.Elements().end());
             
}

/// ************************************************************************************************          
/// ************************************************************************************************ 
void  Calculate_Edges(Element::GeometryType& geom,
		      const compressed_matrix<int>& Coord,
		      int*  edge_ids,
		      array_1d<int,6>& aux
		      )
{           
             int index_0 = geom[0].Id()-1;
             int index_1 = geom[1].Id()-1;
             int index_2 = geom[2].Id()-1; 
	     
             aux[0] = geom[0].Id();
	     aux[1] = geom[1].Id();
	     aux[2] = geom[2].Id();
            //------------------------------------------------------------------------- 
            if(index_0 > index_1)
	      aux[3] = Coord(index_1, index_0);
            else
              aux[3] = Coord(index_0, index_1); 

  
            if(index_1 > index_2)
               aux[4] = Coord(index_2, index_1);
            else 
	       aux[4] = Coord(index_1, index_2 );

	    
            if(index_2 > index_0)
               aux[5] = Coord(index_0, index_2);
            else 
	       aux[5] = Coord(index_2, index_0 );           
	    //-------------------------------------------------------------------------
	    
	    //edge 01
	    if(aux[3] < 0)
		if(index_0 > index_1) edge_ids[0] = 0;                           
		else edge_ids[0] = 1;
	    else
		edge_ids[0] = 3;

	    //edge 12
	    if(aux[4] < 0)
		if(index_1 > index_2) edge_ids[1] = 1;
		else edge_ids[1] = 2;
	    else
		edge_ids[1] = 4;

	    //edge 20
	    if(aux[5] < 0)
		if(index_2 > index_0) edge_ids[2] = 2;
		else edge_ids[2] = 0;
	    else
		edge_ids[2] = 5;
  
}      


/// ************************************************************************************************          
/// ************************************************************************************************ 


void Renumering_Elements_And_Nodes( ModelPart& this_model_part,
PointerVector< Element >& New_Elements)
{
    
   unsigned int id_node =  1;
   unsigned int id_elem =  1;
   NodesArrayType& pNodes =  this_model_part.Nodes();  
   NodesArrayType::iterator i_begin=pNodes.ptr_begin();
   NodesArrayType::iterator i_end=pNodes.ptr_end();


   for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)     
	{
           if(i->Id()!=id_node){i->SetId(id_node);}
            id_node++;     
        }

  ElementsArrayType& rElements         =  this_model_part.Elements(); 
  ElementsArrayType::iterator it_begin =  rElements.ptr_begin();
  ElementsArrayType::iterator it_end   =  rElements.ptr_end(); 

  for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
        {
          if(it->Id()!=id_elem){it->SetId(id_elem);}
          id_elem++;
        }

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


void InterpolateInteralVariables(const int&  nel,
                                 const Element::Pointer father_elem,
				 Element::Pointer child_elem,
				 ProcessInfo& rCurrentProcessInfo)
    {
      std::vector<Vector> values;
      father_elem->GetValueOnIntegrationPoints(INTERNAL_VARIABLES, values, rCurrentProcessInfo); 
     /* /// WARNING =  Calculando la longitud ponderada de fractura del elemento. Solo valido para isotropic_damage 
      Element::GeometryType& geom_father = father_elem->GetGeometry();
      Element::GeometryType& geom_child  = child_elem->GetGeometry();
      double area_father = geom_father.Area();  
      double area_child  = geom_child.Area();
      values[0][4]       = (area_child/area_father) * values[0][4]; 
      */
      child_elem->SetValueOnIntegrationPoints (INTERNAL_VARIABLES, values, rCurrentProcessInfo);

    }

 };

}  // namespace Kratos.

#endif // KRATOS_PROJECTION  defined 


