/*
==============================================================================
KratosPFEMApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


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
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-11-19 15:38:01 $
//   Revision:            $Revision: 1.1 $
//
//
#if !defined(KRATOS_CHOOSE_ELEMENT_PROCESS_INCLUDED)
#define KRATOS_CHOOSE_ELEMENT_PROCESS_INCLUDED


#include <string>
#include <iostream>
#include <algorithm>

#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
#include "geometries/triangle_2d_3.h"
//#include "kratos/applications/MeshingApplication/meshing_application.h"

namespace Kratos
{ 
	class ChooseElementProcess 
		: public Process 
	 {
	   public:

	      ChooseElementProcess(ModelPart& model_part, unsigned int domain_size)
			:Process(), mr_model_part(model_part), mdomain_size(domain_size)
		{
		}

	      /// Destructor.
	      virtual ~ChooseElementProcess()
		{
		}
      

	      ///@}
	      ///@name Operators 
	      ///@{

	      void operator()()
		{
		  Execute();
		}
		
	   virtual void Execute()
		 {
			KRATOS_WATCH("++++++++++++++++++++BEGIN OF CHOOSE PROCESS ^^^^^^^^^^^^^^^^^^^^^^");
			ModelPart::ElementsContainerType ElemPart;
			//KRATOS_WATCH(mr_model_part.Elements().size());
			ElemPart.reserve(mr_model_part.Elements().size());

								
//			Element const& rEl1 = KratosComponents<Element>::Get("Fluid2DASGS");
//			Element const& rEl2 = KratosComponents<Element>::Get("ASGSPRDC");

			Element const& rEl1 = KratosComponents<Element>::Get("ASGSCOMPPRDC2D"); //water element
			Element const& rEl2 = KratosComponents<Element>::Get("ASGSCompressible2D"); // air element

			for(ModelPart::ElementsContainerType::iterator Belem = mr_model_part.ElementsBegin(); Belem != mr_model_part.ElementsEnd(); ++Belem)
			{
				Geometry< Node<3> >& geom = Belem->GetGeometry();

				//choose type:
				// everywhere El1 is chosen unless 3 nodes have IS_PROUS = 0.0 or when the node with IS_WATER = 1.0 is on the boundary
				double chooseflag = 0.0;
		/*	
				double first = geom[0].FastGetSolutionStepValue(IS_WATER);
				double second = geom[1].FastGetSolutionStepValue(IS_WATER);
				double third = geom[2].FastGetSolutionStepValue(IS_WATER);

				//three node of the same kind
			 	if(first == second && second==third)
					chooseflag = first;
				//IS_WATER = 1 is on the boundary
				else
				  {
					//KRATOS_WATCH("***********INSIDE NOT SIMILAR POINTS ******************");
				     for(int ii=0;ii<3;++ii)
					 if(geom[ii].GetSolutionStepValue(IS_WATER) == 1.0 //&& geom[ii].GetSolutionStepValue(IS_STRUCTURE) != 1.0//)
					    {
						chooseflag = 1.0;
						
				        	//KRATOS_WATCH("***********ASGSPRDC IS CHOSEN ******************");
					    }	
				    
						
				   }


				//take components to create element
				unsigned int ele_id = Belem->Id();
				
				
		*/		

				if(Belem->GetValue(IS_WATER_ELEMENT) == 1.0)
					chooseflag = 1.0;
				else
					chooseflag = 0.0;					


				if(chooseflag)
				  {
					Element::Pointer p_elem = rEl1.Create(Belem->Id(),geom, Belem->pGetProperties());
					ElemPart.push_back(p_elem);
					p_elem->GetValue(IS_WATER_ELEMENT) = 1.0;						
					//copy element of other type to consider two elements in divided element
						/*if(Belem->GetValue(IS_DIVIDED) == 1.0)
							{
								Element::Pointer p_elem_second = rEl2.Create(ele_id,geom, Belem->pGetProperties());
								p_elem_second->GetValue(IS_DIVIDED) = 1.0;
								ElemPart.push_back(p_elem_second);
								
				  			}*/
				  }
			        else
				  {
					Element::Pointer p_elem = rEl2.Create(Belem->Id(), geom, Belem->pGetProperties() );
					ElemPart.push_back(p_elem);
					p_elem->GetValue(IS_WATER_ELEMENT) = 0.0;

					//copy element of other type to consider two elements in divided element
						/*if(Belem->GetValue(IS_DIVIDED) == 1.0)
							{
								Element::Pointer p_elem_second = rEl1.Create(ele_id,geom, Belem->pGetProperties());
								p_elem_second->GetValue(IS_DIVIDED) = 1.0;
								ElemPart.push_back(p_elem_second);
	
				  			}*/
//KRATOS_WATCH("***********ASGSCompressible2D IS CHOSEN ******************");
				  }

				
			}
			//KRATOS_WATCH(ElemPart.size());
			//KRATOS_WATCH((mr_model_part.Elements()).size());

			mr_model_part.Elements() = ElemPart;

			//KRATOS_WATCH(mr_model_part.Elements().size());
			KRATOS_WATCH("++++++++++++++++++++END OF CHOOSE PROCESS ^^^^^^^^^^^^^^^^^^^^^^");
		 }

		private:
			ModelPart& mr_model_part;
			unsigned int mdomain_size;

	};

}//namespace kratos

#endif
