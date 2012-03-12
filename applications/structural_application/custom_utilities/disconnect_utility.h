
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
*   Last Modified by:    $Author: Nelson $
*   Date:                $Date: 2011-03-29 11:41:31 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/

#if !defined(KRATOS_DISCONNECT_TRIANGLES_INCLUDED)
#define  KRATOS_DISCONNECT_TRIANGLES_INCLUDED
//System includes
#ifdef _OPENMP
#include <omp.h>
#endif
#include "utilities/openmp_utils.h"

//External includes
#include "boost/smart_ptr.hpp"
#include <cmath>

//Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "custom_utilities/sd_math_utils.h"
#include "custom_utilities/joint.h"

#include "includes/model_part.h"
#include "includes/mesh.h"
#include "geometries/geometry.h"
#include "includes/element.h"
#include "includes/variables.h"



namespace Kratos
{   
	
        class Disconnect_Triangle_Utilities
        {

	  public:
	    	  
	  typedef ModelPart::ConditionsContainerType ConditionsArrayType;
	  typedef ModelPart::ElementsContainerType   ElementsArrayType;
	  typedef ModelPart::NodesContainerType      NodesArrayType;
	  typedef Joint<4> Joint2D;
	  
	  KRATOS_CLASS_POINTER_DEFINITION(Disconnect_Triangle_Utilities);
	  
	  Disconnect_Triangle_Utilities(){}
	  
	  Disconnect_Triangle_Utilities(ModelPart& model_part) {}
	  
	 ~Disconnect_Triangle_Utilities(){}
	  
	  
	  
	  
	  //**************************************************************************************
	  //**************************************************************************************
	   
	  std::vector<Joint2D>::iterator Begin()
	  {
	    return mJointsArray.begin();
	  }

	  std::vector<Joint2D>::iterator End()
	  {
	    return mJointsArray.end();
	  }
	 
	 
	  //**************************************************************************************
	  //**************************************************************************************
	  
	  void CreateJoints(ModelPart& model_part, unsigned int& dimension)
	   {
              Disconnect_Elements(model_part, dimension); 
	      //Disconnect_Elements_DG(model_part, dimension); 
	   }
	  
	  
	  //**************************************************************************************
	  //**************************************************************************************
	  
	  
	  /// Desconecta los elementos y crea los joints
	  void Disconnect_Elements(ModelPart& model_part, unsigned int& dimension)
	  {
	     KRATOS_TRY
	     Disconnect(model_part);  
	     
	     if(dimension==2)
	         CreateJoints2D(model_part);
	     
	     KRATOS_CATCH("")
	  }
	 
	 
	 //**************************************************************************************
	 //**************************************************************************************
	 
	 /// Desconecta los elementos para permirtir hacer DG. Solo elementos Triangulares
	 void Disconnect_Elements_DG(ModelPart& model_part, unsigned int& dimension)
	  {
	     KRATOS_TRY
	     Disconnect(model_part);  
	     if(dimension==2)
	        CalculateNeighbourNodes2D(model_part);
	     else
	       CalculateNeighbourNodes3D(model_part);  
	     KRATOS_CATCH("")
	  }
	 //**************************************************************************************
	 //**************************************************************************************
	  
	 
	  private:
	    
	  /// Solo valido para triangulos
	  void CalculateNeighbourNodes2D(ModelPart& model_part)
	  {
	    KRATOS_TRY	    
	    ElementsArrayType& pElements  = model_part.Elements();
	    vector<unsigned int> element_partition;
	    int number_of_threads = OpenMPUtils::GetNumThreads();
	    OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition); 
	    int  count_2  = 0; 
	    int  count_1  = 0;
	    
	    #pragma omp parallel for private(count_2, count_1)  
	    for(int k=0; k<number_of_threads; k++){
	    ElementsArrayType::iterator it_begin = pElements.ptr_begin()+element_partition[k];
	    ElementsArrayType::iterator it_end   = pElements.ptr_begin()+element_partition[k+1];
	    for(ElementsArrayType::iterator it=it_begin; it!= it_end; it++){
	    if(it->GetProperties()[IS_DISCRETE]>=1.00){
		WeakPointerVector< Node<3> >& neighb_nodes  = it->GetValue(NEIGHBOUR_NODES);
		Element::GeometryType& geom_1 = it->GetGeometry();
		neighb_nodes.clear();
		neighb_nodes.resize(6);
		neighb_nodes(0) = geom_1(1);
		neighb_nodes(1) = geom_1(2);
		neighb_nodes(2) = geom_1(2);
		neighb_nodes(3) = geom_1(0);
		neighb_nodes(4) = geom_1(0);
		neighb_nodes(5) = geom_1(1);
		WeakPointerVector< Element >& neighb_elems  = it->GetValue(NEIGHBOUR_ELEMENTS);  
		count_1 = 0; 
		count_2 = 0;
		for(WeakPointerVector< Element >::iterator neighb = neighb_elems.begin(); neighb!=neighb_elems.end();  ++neighb){    
	          if(neighb->GetProperties()[IS_DISCRETE]>=1.00 && it->Id()!= neighb->Id()){
	              Element::GeometryType& geom_2 = neighb->GetGeometry();  
	              const int& Id = it->Id();
	              count_2 = SearchEdge(neighb, Id);
	              EdgesNodes(geom_2, count_2, count_1, neighb_nodes);
	            }
	          count_1++;
	          }         
	        }
	      }
	    }    
	    KRATOS_CATCH("")
	  }
	  
	  //**************************************************************************************
	  //**************************************************************************************
	  
	 
	  /// Solo valido para tetrahedros
	  void CalculateNeighbourNodes3D(ModelPart& model_part)
	  {
	    KRATOS_TRY	    
	    ElementsArrayType& pElements  = model_part.Elements();
	    vector<unsigned int> element_partition;
	    int number_of_threads = OpenMPUtils::GetNumThreads();
	    OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition); 
	    int  count_2  = 0; 
	    int  count_1  = 0;
	    
	    #pragma omp parallel for private(count_2, count_1)  
	    for(int k=0; k<number_of_threads; k++){
	    ElementsArrayType::iterator it_begin = pElements.ptr_begin()+element_partition[k];
	    ElementsArrayType::iterator it_end   = pElements.ptr_begin()+element_partition[k+1];
	    for(ElementsArrayType::iterator it=it_begin; it!= it_end; it++){
	    if(it->GetProperties()[IS_DISCRETE]>=1.00){
		WeakPointerVector< Node<3> >& neighb_nodes  = it->GetValue(NEIGHBOUR_NODES);
		Element::GeometryType& geom_1 = it->GetGeometry();
		neighb_nodes.clear();
		neighb_nodes.resize(12);
		
		/// usando regla de mano derecha
		//face 1 of neigh in 0
		neighb_nodes(0) = geom_1(1);
		neighb_nodes(1) = geom_1(3);
		neighb_nodes(2) = geom_1(2);
		
		//face 2 of neigh in 1
		neighb_nodes(3) = geom_1(0);
		neighb_nodes(4) = geom_1(2);
		neighb_nodes(5) = geom_1(3);
		
		//face 3 of neigh in 2
		neighb_nodes(6) = geom_1(0);
		neighb_nodes(7) = geom_1(3);
		neighb_nodes(8) = geom_1(1);
		
		//face 4 of neigh in 3
		neighb_nodes(9)  = geom_1(0);
		neighb_nodes(10) = geom_1(1);
		neighb_nodes(11) = geom_1(2);
		
		
		WeakPointerVector< Element >& neighb_elems  = it->GetValue(NEIGHBOUR_ELEMENTS);  
		count_1 = 0; 
		count_2 = 0;
		for(WeakPointerVector< Element >::iterator neighb = neighb_elems.begin(); neighb!=neighb_elems.end();  ++neighb){
	          if(neighb->GetProperties()[IS_DISCRETE]>=1.00 && it->Id()!= neighb->Id()){
	              Element::GeometryType& geom_2 = neighb->GetGeometry();  
	              const int& Id = it->Id();
	              count_2 = SearchEdge(neighb, Id);
	              EdgesNodes3D(geom_1, geom_2, count_1, count_2, neighb_nodes);
	            }
	          count_1++;
	          }          
	        }
	      }
	    }  
	    KRATOS_CATCH("")
	  }
	  
	  //**************************************************************************************
	  //**************************************************************************************
	  
	 /// Separate triangle and tetrahedra elements creating new nodes
	 void Disconnect(ModelPart& model_part)
	  {
	     KRATOS_TRY
	     NodesArrayType& pNodes        = model_part.Nodes(); 
	     unsigned int New_Id           = pNodes.size();
	     NodesArrayType New_pNodes; 
             NodesArrayType::iterator i_begin =  pNodes.begin();
	     NodesArrayType::iterator i_end   =  pNodes.end();
	     const int dis    =  pNodes.end() -  pNodes.begin();
	     std::size_t i    = 0;
	     ModelPart::NodeIterator inode = model_part.Nodes().begin() + i;
	     
	     std::cout<<std::endl;
	     std::cout<< "DUPLICATING NODES IN MODEL PART"<< std::endl;
	     while(inode!= model_part.Nodes().begin() + dis){  
	        inode = model_part.Nodes().begin() + i;
                WeakPointerVector< Element >& neighb_elems  = inode->GetValue(NEIGHBOUR_ELEMENTS); 
		if(neighb_elems.size()!=0){
		  for(WeakPointerVector<Element>::iterator ielem = neighb_elems.begin();  ielem!=neighb_elems.end()-1; ielem++){
		     if(ielem->GetProperties()[IS_DISCRETE]>=1.00){
		        inode = model_part.Nodes().begin() + i;
		        Element::GeometryType& geom = ielem->GetGeometry();
		          for(unsigned int j=0; j<geom.size(); j++){
		             if(geom(j)->Id()==inode->Id()){
			        New_Id++;
			        Create_New_Node(model_part,New_Id, geom(j));
			        break;
		              } } 
		       inode = model_part.Nodes().begin() + i;
		  } } }
		i++;
		inode = model_part.Nodes().begin() + i;
	     } 
	     
	     std::cout<< "DUPLICATING NODES IN MODEL PART FINISH" << std::endl; 
	     KRATOS_CATCH("")
	  }
	  
	  
	 
	 //**************************************************************************************
	 //**************************************************************************************
	 
	 void CreateJoints2D(ModelPart& model_part)
	 {
	   KRATOS_TRY
	      ElementsArrayType& pElements  = model_part.Elements(); 
	      unsigned int count   = 0;
	      unsigned int count_2 = 0;
              
	      Joint2D rJoint;
	      vector<unsigned int> element_partition;
	      
	      ElementsArrayType::iterator it_begin = pElements.ptr_begin();
	      ElementsArrayType::iterator it_end   = pElements.ptr_end();
	      int number_of_threads =  OpenMPUtils::GetNumThreads();
	      OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);
	      
	      
	      
	      #pragma omp parallel for
	      for(int k=0; k<number_of_threads; k++)
	        {
	            ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
	            ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];
	            for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
	                it->GetValue(IS_INACTIVE) = false;
		}
		
		
	      for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
	      {
		
		 WeakPointerVector< Condition >& neighb_cond  = it->GetValue(NEIGHBOUR_CONDITIONS); 
		 if(neighb_cond.size()!=0){
		   for(WeakPointerVector< Condition >::iterator rcond = neighb_cond.begin(); rcond!=neighb_cond.end(); ++rcond)
		   {
		     WeakPointerVector< Element >& neighb_elem_c = rcond->GetValue(NEIGHBOUR_ELEMENTS);
		     neighb_elem_c.push_back(*(it.base()));
		   }
		 }
		 
		 if(it->GetProperties()[IS_DISCRETE]>=1.00){
		 count = 0;
		 Element::GeometryType& geom_1 = it->GetGeometry();
		 WeakPointerVector< Element >& neighb_elems  = it->GetValue(NEIGHBOUR_ELEMENTS);  
	         for(WeakPointerVector< Element >::iterator neighb = neighb_elems.begin(); neighb!=neighb_elems.end();  ++neighb)
		 {    
		   if(neighb->GetProperties()[IS_DISCRETE]>=1.00 && neighb->GetValue(IS_INACTIVE)==false && it->Id()!= neighb->Id())
		   {
		     TheNodes(geom_1, count, 0, 1, rJoint);
		     Element::GeometryType& geom_2 = neighb->GetGeometry();
		     const unsigned int& Id        = it->Id();
		     count_2                       = SearchEdge(neighb, Id);     
		     TheNodes(geom_2, count_2, 2, 3, rJoint);
		     mJointsArray.push_back(rJoint); 
		     it->GetValue(IS_INACTIVE) = true;
		   }
		   count++;
		 }
	       }
	      }
	     
	     /// Check the conditions
	     ConditionsArrayType& pConditions = model_part.Conditions();
	     for(ConditionsArrayType::iterator rcond = pConditions.begin(); rcond!=pConditions.end(); ++rcond)
	      {
		 CheckConditions(rcond);
	      }
	   
	   KRATOS_CATCH("")
	 }
	 
	 
         //**************************************************************************************
	 //************************************************************************************** 
	 
	/// Sirve para actulizar los nodos de las condiciones cuando son duplicados solo en condiciones de linea
	void CheckConditions(ConditionsArrayType::iterator& rcond)
	{
	  int a              = 0; 
	  int b              = 1;
	  double check       = 0.00;
	  const double toler = 1E-8;
	  a = 0;
	  b = 1; 
	  Condition::GeometryType& geom_cond  = rcond->GetGeometry();
	  Element::Pointer relem              = (rcond->GetValue(NEIGHBOUR_ELEMENTS))(0).lock();
	  Element::GeometryType& geom_elem    = relem->GetGeometry();
	  array_1d<double,3> cond             = geom_cond.GetPoint(1) - geom_cond.GetPoint(0); 
	  array_1d<double,3> elem;              
	  noalias(cond) = (1.00/norm_2(cond)) * cond;
	  for(unsigned int i = 0; i<geom_elem.size(); i++)
	  {
	  if(i==2) {b = 0; a = 2;}
	  elem          = geom_elem.GetPoint(b) - geom_elem.GetPoint(a);
	  noalias(elem) = (1.00/norm_2(elem)) * elem;
	  check         = std::fabs(inner_prod(elem, cond));
	  if( std::fabs(check-1.00)<toler )
	  break;
	  b++; a++;
	  }
	  /// Las condiciones se nombran contrario a las manecillas del reloj 
	  geom_cond(0) = geom_elem(b); 
	  geom_cond(1) = geom_elem(a);
	}
	 
	 
        //**************************************************************************************
	//************************************************************************************** 
	 
	int SearchEdge(WeakPointerVector<Element>::iterator& this_elem, const unsigned int& Id)
	 {
	    int count = 0;
	    WeakPointerVector< Element >& neighb_elems = this_elem->GetValue(NEIGHBOUR_ELEMENTS); 
	    for(WeakPointerVector<Element>::iterator ielem = neighb_elems.begin();  ielem!=neighb_elems.end(); ++ielem)
	    {
	      if(ielem->Id()!=Id){count ++;}
	      else 
		break;
	    }
	    return count;
	 }
	  
	  
	 //**************************************************************************************
	 //**************************************************************************************
	  
	 ///Subrutina para tener los nodos correspondientes segun DG en 2D 
	 void EdgesNodes(const Element::GeometryType& geom, const int& count_2, const int& count_1, WeakPointerVector< Node<3> >& neighb_nodes)
	 {
	   
	      int i, j;
	      if(count_1==0){i = 0; j = 1;}
	      if(count_1==1){i = 2; j = 3;}
	      if(count_1==2){i = 4; j = 5;}

	      if(count_2==0)
	      {
	        neighb_nodes(i) =  geom(2);
	        neighb_nodes(j) =  geom(1); 
	      }
	      else if (count_2==1)
	      {
	        neighb_nodes(i)  =  geom(0);
	        neighb_nodes(j)  =  geom(2); 
	      }
	      else 
	      {
	        neighb_nodes(i)  =  geom(1);
	        neighb_nodes(j)  =  geom(0); 
	      }
	 }  
	 
	 //**************************************************************************************
	 //**************************************************************************************
	 
	 ///Subrutina para tener los nodos correspondientes segun DG en 2D
	 void EdgesNodes3D(const Element::GeometryType& geom_1, 
			   const Element::GeometryType& geom_2, 
			   const int& count_1, const int& count_2, 
			   WeakPointerVector< Node<3> >& neighb_nodes)
	 {
	      int i = 0, j = 0, k = 0;
	      array_1d<int,3> a;
	      array_1d<int,3> c;

	      array_1d<double,3> diff = ZeroVector(3);
	      if(count_1==0)
	      {
		i = 0; j = 1; k = 2;
		c[0] = 1; c[1] = 3; c[2] = 2;
	      }
	      else if(count_1==1){
		i = 3; j = 4; k = 5;
		c[0] = 0; c[1] = 2; c[2] = 3;
	      }
	      else if(count_1==2){
		i = 6; j = 7; k = 8;
		c[0] = 0; c[1] = 3; c[2] = 1;
	      }
              else {
		i = 9; j = 10; k = 11;
		c[0] = 0; c[1] = 1; c[2] = 2;
	      }
              
              array_1d<int,3> b;  b[0] = i; b[1] = j; b[2] = k;     
	      if(count_2==0)
	      {
		a[0] = 1; 
		a[1] = 3;
		a[2] = 2;
		for(int l = 0; l<3; l++){
		   for(int m=0; m<3; m++){
		    noalias(diff)   =  geom_1[c[l]] - geom_2[a[m]];
		    if(inner_prod(diff, diff)<1E-9){
		       neighb_nodes(b[l]) =  geom_2(a[m]);
		       break;
		    } } }
	      }
	      else if (count_2==1)
	      {
	        a[0] = 0; 
		a[1] = 2;
		a[2] = 3;
		for(int l = 0; l<3; l++){
		   for(int m=0; m<3; m++){
		    noalias(diff)   =  geom_1[c[l]] - geom_2[a[m]];
		    std::fabs(inner_prod(diff, diff));
		    if(inner_prod(diff, diff)<1E-9){
		       neighb_nodes(b[l]) =  geom_2(a[m]);
		       break;
		    } } }	      
	      }
	      else if(count_2==2) 
	      {
	        a[0] = 0; 
		a[1] = 3;
		a[2] = 1;
		for(int l = 0; l<3; l++){
		   for(int m=0; m<3; m++){
		    noalias(diff)   =  geom_1[c[l]] - geom_2[a[m]];
		    std::fabs(inner_prod(diff, diff));
		    if(inner_prod(diff, diff)<1E-9){
		       neighb_nodes(b[l]) =  geom_2(a[m]);
		       break;    
		    } } }
	      }
	      else
	      {
	        a[0] = 0; 
		a[1] = 1;
		a[2] = 2;
		for(int l = 0; l<3; l++){
		   for(int m=0; m<3; m++){
		    noalias(diff)   =  geom_1[c[l]] - geom_2[a[m]];
		    std::fabs(inner_prod(diff, diff));
		    if(inner_prod(diff, diff)<1E-9){
		       neighb_nodes(b[l]) =  geom_2(a[m]);
		       break;
		    } } }
	      }      
	 }  
	 
	 
	 //**************************************************************************************
	 //**************************************************************************************
	 
	 void TheNodes(const Element::GeometryType& geom, const int& count, const int& i, const int& j, Joint2D& rJoint)
	 {
	      if(count==0)
	      {
	        rJoint.InsertNode(i, geom(1));
	        rJoint.InsertNode(j, geom(2)); 
	      }
	      else if (count==1)
	      {
	        rJoint.InsertNode(i, geom(2));
	        rJoint.InsertNode(j, geom(0)); 
	      }

	      else if (count==2)
	      {
	        rJoint.InsertNode(i, geom(0));
	        rJoint.InsertNode(j, geom(1)); 
	      }
	 }
	  

	///************************************************************************************************
        ///************************************************************************************************
	
	void Create_New_Node(ModelPart& this_model_part, unsigned int& New_Id, Node<3>::Pointer& pNode)
        {
            //node to get the DOFs from
            int step_data_size           =  this_model_part.GetNodalSolutionStepDataSize();
            array_1d<double, 3 > & Coord =  pNode->Coordinates();
            Node < 3 > ::DofsContainerType& reference_dofs = (this_model_part.NodesBegin())->GetDofs();

            Node<3>::Pointer pnode = this_model_part.CreateNewNode(New_Id, Coord[0], Coord[1], Coord[2]);
            pnode->SetBufferSize(this_model_part.NodesBegin()->GetBufferSize());

            //generating the dofs
            for (Node < 3 > ::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
            {
                Node < 3 > ::DofType& rDof = *iii;
                Node < 3 > ::DofType::Pointer p_new_dof = pnode->pAddDof(rDof);
                if (pNode->IsFixed(iii->GetVariable()) == true)
                    (p_new_dof)->FixDof();
                else
                {
                    (p_new_dof)->FreeDof();
                }
            }

            unsigned int buffer_size = pnode->GetBufferSize();
            for (unsigned int step = 0; step < buffer_size; step++)
            {
                double* step_data = pnode->SolutionStepData().Data(step);
                double* old_node_data = pNode->SolutionStepData().Data(step);
                //copying this data in the position of the vector we are interested in
                for (signed int j = 0; j < step_data_size; j++)
                {
                    step_data[j] = old_node_data[j];
                }
            }


            const array_1d<double, 3 > & disp = pnode->FastGetSolutionStepValue(DISPLACEMENT);
            pnode->X0() = pnode->X() - disp[0];
            pnode->Y0() = pnode->Y() - disp[1];
            pnode->Z0() = pnode->Z() - disp[2];

            array_1d<double, 3 > & vel_old = pNode->FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3 > & vel_new = pnode->FastGetSolutionStepValue(VELOCITY);
            vel_new =  vel_old;
	    
            const array_1d<double, 3 > & accel_old = pNode->FastGetSolutionStepValue(ACCELERATION);
            array_1d<double, 3 > & accel_new = pnode->FastGetSolutionStepValue(ACCELERATION);
            accel_new = accel_old;
	    
	    pNode = pnode;
	    
        }
	  
	 
	 //**************************************************************************************
	 //**************************************************************************************  
	  
	  
	  private:
	  std::vector<Joint2D> mJointsArray;
	 

	  
       }; 
          
 
}//namespace Kratos.

#endif /*KRATOS_DISCONNECT_TRIANGLES*/

