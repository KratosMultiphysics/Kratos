
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
	    	  
	  typedef ModelPart::ElementsContainerType ElementsArrayType;
	  typedef ModelPart::NodesContainerType    NodesArrayType;
	  typedef Joint<4> Joint2D;
	  
	  KRATOS_CLASS_POINTER_DEFINITION(Disconnect_Triangle_Utilities);
	  
	  Disconnect_Triangle_Utilities(){}
	  
	  Disconnect_Triangle_Utilities(ModelPart& model_part) 
	  {
	    Disconnect_Elements(model_part);
	  }
	  
	 ~Disconnect_Triangle_Utilities(){}
	  
	  /// Desconecta todos los elementos creando nuevos nodos en el model part
	  void Disconnect_Elements(ModelPart& model_part)
	  {
	     KRATOS_TRY
	     
	     NodesArrayType& pNodes        = model_part.Nodes(); 
	     ElementsArrayType& pElements  = model_part.Elements();
	     unsigned int New_Id           = pNodes.size();
	     NodesArrayType New_pNodes; 
             NodesArrayType::iterator i_begin = pNodes.begin();
	     NodesArrayType::iterator i_end   = pNodes.end();
	     const int dis =  pNodes.end() -  pNodes.begin();
	     int i   = 0;
	     
	     for(ModelPart::NodeIterator inode=i_begin; inode!= pNodes.begin() + dis ; ++inode)     
	     {
	        inode = pNodes.begin() + i;  
                WeakPointerVector< Element >& neighb_elems  = inode->GetValue(NEIGHBOUR_ELEMENTS); 
		for(WeakPointerVector<Element>::iterator ielem = neighb_elems.begin();  ielem!=neighb_elems.end()-1; ielem++){
		  if(ielem->GetProperties()[IS_DISCRETE]>=1.00){
		       inode = pNodes.begin() + i;
		       Element::GeometryType& geom = ielem->GetGeometry();
		      if(geom(0)->Id()==inode->Id()){
			   New_Id++;
			   Create_New_Node(model_part,New_Id, geom(0));
		      }
	              if(geom(1)->Id()==inode->Id()){ 
			   New_Id++;
			   Create_New_Node(model_part, New_Id, geom(1));
		      }
		      if(geom(2)->Id()==inode->Id()){ 
			   New_Id++;
			   Create_New_Node(model_part, New_Id, geom(2));
		      }	
		  }  
		}
		i++;
	      }
	      
	     	      
	      int  count   = 0;
	      int  count_2 = 0;
//              bool test_1  = false;
//	      bool test_2  = false;   
	      Joint2D rJoint;
	      ElementsArrayType::iterator it_begin=pElements.ptr_begin();
	      ElementsArrayType::iterator it_end=pElements.ptr_end();
	      
	      for(ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
	        it->GetValue(IS_INACTIVE) = false;
	      
	      
	      for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
	      {
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
		     const int& Id = it->Id();
		     count_2 = SearchEdge(neighb, Id);     
		     TheNodes(geom_2, count_2, 2, 3, rJoint);
		     mJointsArray.push_back(rJoint); 
		     it->GetValue(IS_INACTIVE) = true;
		   }
		   count++;
		 }
	       }
	      }
	      
	     KRATOS_CATCH("")
	     
	  }
	  
	  
	int SearchEdge(WeakPointerVector<Element>::iterator& this_elem, const int& Id)
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

            const array_1d<double, 3 > & vel_old = pNode->FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3 > & vel_new = pnode->FastGetSolutionStepValue(VELOCITY);
            vel_new = vel_old;

            const array_1d<double, 3 > & accel_old = pNode->FastGetSolutionStepValue(ACCELERATION);
            array_1d<double, 3 > & accel_new = pnode->FastGetSolutionStepValue(ACCELERATION);
            accel_new = accel_old;

	    if(pnode->Id()==273)
	    {
	       array_1d<double, 3 > & force = pnode->FastGetSolutionStepValue(FORCE);
	       force =  ZeroVector(3);
	    }
	    
	    pNode = pnode;
	    
        }
	  
	  void CreateJoints(ModelPart& model_part)
	  {
	    Disconnect_Elements(model_part);
	  }
	  
	  std::vector<Joint2D>::iterator Begin()
	  {
	    return mJointsArray.begin();
	  }
	  
	  std::vector<Joint2D>::iterator End()
	  {
	    return mJointsArray.end();
	  }
	  
	  
	  private:
	  std::vector<Joint2D> mJointsArray;
	 

	  
       }; 
          
 
}//namespace Kratos.

#endif /*KRATOS_DISCONNECT_TRIANGLES*/

