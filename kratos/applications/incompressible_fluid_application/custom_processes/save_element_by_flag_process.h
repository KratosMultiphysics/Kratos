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
#if !defined(KRATOS_SAVE_ELEMENT_BY_FLAG_PROCESS_INCLUDED)
#define KRATOS_SAVE_ELEMENT_BY_FLAG_PROCESS_INCLUDED


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
	class SaveElementByFlagProcess 
		: public Process 
	 {
	   public:

	      SaveElementByFlagProcess(ModelPart::ElementsContainerType& elements_model_part,ModelPart::ElementsContainerType& fluid_elements, Kratos::Variable<int>& flag, int value )
			:Process(), mr_elements(elements_model_part), mr_fluid_elements(fluid_elements), m_val(value), m_flag(flag)
		{
		}

	      /// Destructor.
	      virtual ~SaveElementByFlagProcess()
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
			//KRATOS_WATCH("++++++++++++++++++++BEGIN OF SaveElementByFlagProcess PROCESS ^^^^^^^^^^^^^^^^^^^^^^");
			//ModelPart::ElementsContainerType ElemPart;
			//KRATOS_WATCH(mr_model_part.Elements().size());
			mr_fluid_elements.clear();

								
//			Element const& rEl1 = KratosComponents<Element>::Get("Fluid2DASGS");
//			Element const& rEl2 = KratosComponents<Element>::Get("ASGSPRDC");

			//Element const& rEl1 = KratosComponents<Element>::Get("ASGSCOMPPRDC2D"); //water element
			//Element const& rEl2 = KratosComponents<Element>::Get("ASGSCompressible2D"); // air element

			for(ModelPart::ElementsContainerType::iterator Belem = mr_elements.begin(); Belem != mr_elements.end(); ++Belem)
			{
			      if(Belem->GetValue(m_flag) != m_val)
					mr_fluid_elements.push_back( *(Belem.base()) )  ;


				
			}
// 			KRATOS_WATCH("ALL ELEMENTS");
// 			KRATOS_WATCH(mr_elements.size());
// 
// 			KRATOS_WATCH("flag Elements");
// 			KRATOS_WATCH((mr_fluid_elements).size());

			KRATOS_WATCH("++++++++++++++++++++END OF SaveElementByFlagProcess PROCESS ^^^^^^^^^^^^^^^^^^^^^^");
		 }

		private:
			ModelPart::ElementsContainerType& mr_elements;
			ModelPart::ElementsContainerType& mr_fluid_elements;
		         int m_val;
			Kratos::Variable<int>& m_flag;

	};

}//namespace kratos

#endif
